from pathlib import Path
import subprocess


# SIMT-only GLASS sources vendored into every generated grid.cuh. The
# cuBLASDx-backed `glass::nvidia` namespace was removed in v2.0; see
# docs/source/user_guide/concepts/cublasdx_removal_design.rst for the
# rationale and the `archive/last-cublasdx` git tag for the historical
# vendoring list.
_GLASS_BASE_FILES = [
    "src/base/L1/reduce.cuh",
    "src/base/L1/dot.cuh",
    "src/base/L1/dot_strided.cuh",
    "src/base/L1/dot_strided_coalesced.cuh",
    "src/base/L2/gemv.cuh",
    "src/base/L2/gemv_strided.cuh",
    "src/base/L2/gemv_segmented.cuh",
    "src/base/L3/gemm.cuh",
    "src/base/L3/gemm_strided.cuh",
    "src/base/L3/gemm_batched_indexed.cuh",
]


def _grid_repo_root():
    return Path(__file__).resolve().parents[2]


def _glass_root():
    root = _grid_repo_root() / "GLASS"
    if not root.exists():
        raise FileNotFoundError(
            "GLASS submodule is missing. Run `git submodule update --init GLASS` "
            "from the GRiD-A2R repository root."
        )
    return root


def _glass_commit():
    try:
        return subprocess.check_output(
            ["git", "-C", str(_glass_root()), "rev-parse", "HEAD"],
            text=True,
        ).strip()
    except Exception:
        return "unknown"


def _emit_glass_source_file(self, relative_path):
    source_path = _glass_root() / relative_path
    if not source_path.exists():
        raise FileNotFoundError("Required GLASS source file is missing: " + str(source_path))
    self.gen_add_code_line("// BEGIN GLASS " + relative_path)
    for line in source_path.read_text().splitlines():
        stripped = line.strip()
        if stripped.startswith("#pragma once"):
            continue
        if stripped.startswith("#include"):
            continue
        self.gen_add_code_line(line)
    self.gen_add_code_line("// END GLASS " + relative_path)
    self.gen_add_code_line("")


def gen_linalg_smem_setup(self, temp_size):
    """Emit ``unsigned char *s_linalg_smem = nullptr;``.

    Historically set up an aligned smem pointer for cuBLASDx scratch.
    cuBLASDx was removed in v2.0 (see
    ``docs/source/user_guide/concepts/cublasdx_removal_design.rst``);
    the SIMT linalg path needs no scratch. The function (and the
    ``s_linalg_smem`` local) are kept so the ~100 callsites that pass
    it as the last ``grid_linalg_*`` argument continue to compile
    without per-callsite edits.
    """
    del temp_size  # unused; kept in signature for caller compatibility
    self.gen_add_code_line("unsigned char *s_linalg_smem = nullptr;")


def gen_grid_linalg_backend_helpers(self):
    """Emit the GRiD linear-algebra adapter.

    SIMT GLASS is vendored at codegen time, plus thin
    ``grid_linalg_*`` wrappers that delegate to it. The wrappers
    accept (and ignore) a trailing ``glass_nvidia_smem`` argument so
    pre-v2.0 callsites continue to compile without edits.
    """
    glass_commit = _glass_commit()
    self.gen_add_func_doc("Vendored GLASS linear algebra helpers (SIMT only)")
    self.gen_add_code_line("// Temporarily leave the generated namespace so GLASS keeps its public namespace.")
    self.gen_add_end_control_flow()
    self.gen_add_code_lines([
        "",
        "// Vendored from GLASS at codegen time.",
        "// Source repository: git@github.com:A2R-Lab/GLASS.git",
        "// Pinned commit: " + glass_commit,
        "namespace glass {",
        "",
    ])
    for relative_path in _GLASS_BASE_FILES:
        _emit_glass_source_file(self, relative_path)
    self.gen_add_code_line("} // namespace glass")
    self.gen_add_code_line("")

    self.gen_add_code_line("namespace " + self.file_namespace + " {", True)
    self.gen_add_func_doc("Linear algebra wrappers (SIMT GLASS)")
    self.gen_add_code_lines([
        "// SIMT-only linalg. cuBLASDx was removed in v2.0; the",
        "// `glass_nvidia_smem` parameter on each wrapper is retained for",
        "// caller compatibility and is always ignored. The stub",
        "// `GRID_LINALG_NVIDIA_MAX_HELPER_BYTES<T>()` below returns 0 so",
        "// shared-memory arena calculations continue to compile unchanged.",
        "",
        "template <typename T>",
        "__host__ __device__ constexpr size_t GRID_LINALG_NVIDIA_MAX_HELPER_BYTES() {",
        "    return static_cast<size_t>(0);",
        "}",
        "",
        "template <typename T, int M, int N, int K, bool TRANSPOSE_B = false, bool ROW_MAJOR_A = false, bool ROW_MAJOR_B = false, bool ROW_MAJOR_C = false>",
        "__device__ void grid_linalg_gemm(const T *A, const T *B, T *C, T alpha, T beta, unsigned char *glass_nvidia_smem = nullptr) {",
        "    (void)glass_nvidia_smem;",
        "    T *A_mut = const_cast<T *>(A);",
        "    T *B_mut = const_cast<T *>(B);",
        "    if (!TRANSPOSE_B && !ROW_MAJOR_A && !ROW_MAJOR_B && !ROW_MAJOR_C) {",
        "        glass::gemm<T, M, N, K>(alpha, A_mut, B_mut, beta, C);",
        "    }",
        "    else {",
        "        glass::gemm_ex<T, TRANSPOSE_B, ROW_MAJOR_A, ROW_MAJOR_B, ROW_MAJOR_C>(M, N, K, alpha, A_mut, B_mut, beta, C);",
        "    }",
        "    __syncthreads();",
        "}",
        "",
        "template <typename T, int M, int N, bool TRANSPOSE = false, bool ROW_MAJOR_A = false>",
        "__device__ void grid_linalg_gemv(const T *A, const T *x, T *y, T alpha, T beta, unsigned char *glass_nvidia_smem = nullptr) {",
        "    (void)glass_nvidia_smem;",
        "    T *A_mut = const_cast<T *>(A);",
        "    T *x_mut = const_cast<T *>(x);",
        "    if (!TRANSPOSE && !ROW_MAJOR_A) {",
        "        glass::gemv<T, M, N>(alpha, A_mut, x_mut, beta, y);",
        "    }",
        "    else {",
        "        glass::gemv_ex<T, TRANSPOSE, ROW_MAJOR_A>(M, N, alpha, A_mut, x_mut, beta, y);",
        "    }",
        "    __syncthreads();",
        "}",
        "",
        "template <typename T, int M, int N, int ROW_STRIDE>",
        "__device__ void grid_linalg_row_strided_gemv(const T *A, const T *x, T *y, T alpha, T beta, unsigned char *glass_nvidia_smem = nullptr) {",
        "    (void)glass_nvidia_smem;",
        "    ::glass::row_strided_gemv<T, M, N, ROW_STRIDE>(A, x, y, alpha, beta);",
        "    __syncthreads();",
        "}",
        "",
        "template <typename T, int M, int N, int K, int A_RS, int B_RS>",
        "__device__ void grid_linalg_row_strided_gemm(const T *A, const T *B, T *C, T alpha, T beta, unsigned char *glass_nvidia_smem = nullptr) {",
        "    (void)glass_nvidia_smem;",
        "    ::glass::row_strided_gemm<T, M, N, K, A_RS, B_RS>(A, B, C, alpha, beta);",
        "    __syncthreads();",
        "}",
        "",
        "template <typename T, int N, int S1, int S2>",
        "__device__ T grid_linalg_dot_strided(const T *vec1, const T *vec2) {",
        "    return ::glass::dot_strided<T, N, S1, S2>(vec1, vec2);",
        "}",
        "",
        "// Segmented (batched) row-strided GEMV: `segments` independent M x N GEMVs in one",
        "// block-cooperative pass, base offsets per segment via the descriptor arrays. With",
        "// FUSE_SCALED_ADD, folds a per-segment y += S*scalar add into the single y store.",
        "template <typename T, int M, int N, int ROW_STRIDE = M, bool FUSE_SCALED_ADD = false>",
        "__device__ void grid_linalg_segmented_row_strided_gemv(unsigned int segments, const int *seg_a_off, const int *seg_x_off, const int *seg_y_off, const T *A, const T *x, T *y, T alpha, T beta, const int *seg_s_off = nullptr, const T *S = nullptr, const T *scalar = nullptr, unsigned char *glass_nvidia_smem = nullptr) {",
        "    (void)glass_nvidia_smem;",
        "    ::glass::segmented_row_strided_gemv<T, M, N, ROW_STRIDE, FUSE_SCALED_ADD>(segments, seg_a_off, seg_x_off, seg_y_off, A, x, y, alpha, beta, seg_s_off, S, scalar);",
        "    __syncthreads();",
        "}",
        "",
        "// Indexed batched DIMxDIM (col-major) GEMM: C[c_idx[p]] = A[a_idx[p]] * B[b_idx[p]]",
        "// over a flat in-chain index list (e.g. compacted (ancestor,ee) pairs).",
        "template <typename T, int DIM = 4>",
        "__device__ void grid_linalg_indexed_batched_gemm(unsigned int pairs, const int *a_idx, const int *b_idx, const int *c_idx, const T *A_base, const T *B_base, T *C_base, unsigned char *glass_nvidia_smem = nullptr) {",
        "    (void)glass_nvidia_smem;",
        "    ::glass::indexed_batched_gemm<T, DIM>(pairs, a_idx, b_idx, c_idx, A_base, B_base, C_base);",
        "    __syncthreads();",
        "}",
        "",
        "// Coalesced block-cooperative strided dot: same value as grid_linalg_dot_strided but",
        "// the whole block cooperates on ONE dot with transposed/tiled iteration so consecutive",
        "// threads hit consecutive global addresses (fast when the operand is L2-pinned global).",
        "// Writes the scalar to *out (valid after the trailing barrier); needs ceil(blockDim/32)",
        "// T of s_scratch. NOT a drop-in for grid_linalg_dot_strided (that one is per-thread).",
        "template <typename T, int N, int SX = 1, int SY = 1>",
        "__device__ void grid_linalg_dot_strided_coalesced(const T *x, const T *y, T *out, T *s_scratch) {",
        "    ::glass::dot_strided_coalesced<T, N, SX, SY>(x, y, out, s_scratch);",
        "    __syncthreads();",
        "}",
        ""
    ])


def gen_invert_matrix(self, use_thread_group=False):
    """
    This function generates a matrix inversion function for cuda.
    The function employs Gaussian elimination.
    """

    self.gen_add_func_doc("Compute the inverse of a matrix", ["Uses gaussian elimination"], \
                          ['dimA is number of rows in A', \
                           'A is a pointer to the original invertible matrix. It is turned into an identity matrix', \
                           'Ainv is a pointer to an identity matrix that will be transformed into the inverse of A', \
                            's_temp is a pointer to temporary memory of size 4*dimA'])
    self.gen_add_code_line("template <typename T>")
    self.gen_add_code_line("__device__")
    self.gen_add_code_line("void invert_matrix(uint32_t dimA, T *A, T *Ainv, T *s_temp) {", True)
    self.gen_add_serial_ops(use_thread_group)
    self.gen_add_code_line("for (unsigned pivRC = 0; pivRC < dimA; pivRC++) {", True)   # iterate over diagonal
    self.gen_add_code_line("unsigned pivColOffset = pivRC*dimA;")
    self.gen_add_code_line("T pvInv = static_cast<T>(1)/A[pivRC + pivColOffset];")      # 1/pivot

    # save the pivot row and column values
    self.gen_add_code_line("for (unsigned ind = 0; ind < dimA; ind++) {", True)
    self.gen_add_code_line("s_temp[ind] = static_cast<T>(A[pivRC + dimA * ind]);")
    self.gen_add_code_line("s_temp[ind+dimA] = static_cast<T>(Ainv[pivRC + dimA * ind]);")
    self.gen_add_code_line("s_temp[ind+dimA*2] = static_cast<T>(A[ind + pivColOffset]);")
    self.gen_add_end_control_flow()

    # run gaussian elimination for the pivot row and column. Matrices are stored column-major.
    self.gen_add_code_line("for (unsigned ind = 0; ind < dimA*dimA; ind++) {", True)
    self.gen_add_code_line("unsigned row = ind % dimA, col = ind / dimA;")
    # apply to the pivot row
    self.gen_add_code_line("if (row == pivRC) {", True)
    self.gen_add_code_line("A[row + dimA * col] = s_temp[col] * pvInv;") # put 1 on the diagonal by multiplying row by inverse
    self.gen_add_code_line("Ainv[row + dimA * col] = s_temp[col+dimA] * pvInv;")
    self.gen_add_end_control_flow()
    # apply to other rows by reducing entries on the pivot column to 0s
    self.gen_add_code_line("else {", True)
    self.gen_add_code_line("T multiplier = s_temp[row+dimA*2] / s_temp[pivRC];")
    self.gen_add_code_line("A[row + dimA * col] -= multiplier * s_temp[col];")
    self.gen_add_code_line("Ainv[row + dimA * col] -= multiplier * s_temp[col+dimA];")
    self.gen_add_end_control_flow()

    self.gen_add_end_control_flow()
    self.gen_add_end_control_flow()
    self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)
    self.gen_add_end_function()
    return


def gen_matmul(self):
    """
    Generates the matrix multiplication helper function.
    This function allows for a transpose of B
    """
    self.gen_add_func_doc("Matrix multiplication helper function of AB", [], \
                          ['index - the index of the result vector', \
                           'A - pointer to the first matrix', \
                           'B - pointer to the second matrix', \
                           'dest - pointer to the destination matrix', \
                           'num - 36 or 6 depending on the indexing scheme', \
                           't - true => multiply with the transpose of B'])
    self.gen_add_code_line("template <typename T>")
    self.gen_add_code_line("__device__")
    self.gen_add_code_line("void matmul(int index, T *A, T *B, T *dest, int num, bool t) {", True)
    self.gen_add_code_line("int cur = 36*((index/num)%NUM_JOINTS);")
    self.gen_add_code_line("T *vec1 = &B[cur + (t*5+1)*(index%6)];")
    self.gen_add_code_line("T *vec2 = &A[6*(index/6)];")
    self.gen_add_code_line("dest[index] = dot_prod<T,6, 6, 1>(vec1, vec2);")
    self.gen_add_end_function()


def gen_matmul_trans(self):
    """
    Generates the matrix multiplication helper function where one of the
    matrices is transposed. Both A and B are 6x6 matrices.
    """
    self.gen_add_func_doc("Matrix multiplication helper function where one of the matrices is tranposed.", [], \
                          ['index - the index of the result vector', \
                           'A - pointer to the first 6x6 matrix', \
                           'B - pointer to the second 6x6 matrix', \
                           'dest - pointer to the destination matrix', \
                           'char trans_mat - a for A^TB, b for AB^T'])
    self.gen_add_code_line("template <typename T>")
    self.gen_add_code_line("__device__")
    self.gen_add_code_line("void matmul_trans(int index, T *A, T *B, T *dest, char trans_mat) {", True)
    self.gen_add_code_line("T *vec1;")
    self.gen_add_code_line("T *vec2;")
    self.gen_add_code_line("if (trans_mat == 'a'){", True)
    self.gen_add_code_line("vec1 = &A[6*(index%6)];")
    self.gen_add_code_line("vec2 = &B[6*(index/6)];")
    self.gen_add_code_line("dest[index] = dot_prod<T,6,1,1>(vec1, vec2);")
    self.gen_add_end_control_flow()
    self.gen_add_code_line("if (trans_mat == 'b'){", True)
    self.gen_add_code_line("vec1 = &A[index%6];")
    self.gen_add_code_line("vec2 = &B[index/6];")
    self.gen_add_code_line("dest[index] = dot_prod<T,6,6,6>(vec1, vec2);")
    self.gen_add_end_control_flow()
    self.gen_add_end_function()


def gen_outer_product(self):
    """
    This function generates the cuda for the outerProduct
    function.
    """
    self.gen_add_func_doc("Compute the outer product between two vectors: dest = ab^T", \
                          ["Function assumes it is called by a single thread."], \
                          ['a - first vector', \
                           'b - second vector', \
                           'dest - destination matrix', \
                           'aLength - length of a', \
                           'bLength - length of b', \
                           'idx - index of resulting matrix to be computed by this thread'])
    self.gen_add_code_line("template <typename T>")
    self.gen_add_code_line("__device__")
    self.gen_add_code_line("void outerProduct(T *a, T *b, T *dest, int aLength, int bLength, int idx) {", True)
    self.gen_add_code_line("int row = idx / bLength;")
    self.gen_add_code_line("int col = idx % bLength;")
    self.gen_add_code_line("if (row < aLength && col < bLength) dest[col * aLength + row] = a[row] * b[col];")
    self.gen_add_end_function()
