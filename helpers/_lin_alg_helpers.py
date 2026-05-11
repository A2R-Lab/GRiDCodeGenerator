from pathlib import Path
import subprocess


_GLASS_BASE_FILES = [
    "src/base/L1/reduce.cuh",
    "src/base/L1/dot.cuh",
    "src/base/L1/dot_strided.cuh",
    "src/base/L2/gemv.cuh",
    "src/base/L2/gemv_strided.cuh",
    "src/base/L3/gemm.cuh",
    "src/base/L3/gemm_strided.cuh",
]

_GLASS_NVIDIA_FILES = [
    "src/nvidia/l1.cuh",
    "src/nvidia/l2.cuh",
    "src/nvidia/l3.cuh",
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


def _ee_gradient_packed_gemm_k_values(self):
    if self.robot.is_serial_chain():
        return []
    n = self.robot.get_num_pos()
    all_ees = self.robot.get_leaf_nodes()
    num_ees = len(all_ees)
    n_bfs_levels = self.robot.get_max_bfs_level() + 1
    values = set()
    for bfs_level in range(1, n_bfs_levels):
        curr_parents = all_ees
        for _ in range(bfs_level):
            curr_parents = [(-1 if jid == -1 else self.robot.get_parent_id(jid)) for jid in curr_parents]
        run_start = 0
        while run_start < num_ees:
            parent_jid = curr_parents[run_start]
            run_end = run_start + 1
            while run_end < num_ees and curr_parents[run_end] == parent_jid:
                run_end += 1
            if parent_jid != -1:
                values.add(4 * n * (run_end - run_start))
            run_start = run_end
    if num_ees > 0:
        values.add(4 * n * num_ees)
    return sorted(values)


def _nvidia_gemm_sizes(self):
    sizes = {(4, 4, 4)}
    for k in _ee_gradient_packed_gemm_k_values(self):
        sizes.add((4, 4, k))
    for size in _nvidia_row_strided_gemm_sizes(self):
        sizes.add(size)
    return sorted(sizes)


def _nvidia_gemv_sizes(self):
    return sorted(set(_nvidia_row_strided_gemv_sizes(self)))


def _nvidia_row_strided_gemv_sizes(self):
    return [(6, 6)]


def _nvidia_row_strided_gemm_sizes(self):
    return [(6, 6, 6)]


def gen_grid_linalg_backend_helpers(self):
    """
    Generate GRiD-owned linear algebra backend wrappers.

    The generated header remains self-contained by copying the required GLASS
    source fragments from the GLASS submodule at codegen time.
    """
    glass_commit = _glass_commit()
    self.gen_add_func_doc("Vendored GLASS linear algebra helpers")
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
    self.gen_add_code_line("#if GRID_CUDA_USE_GLASS_NVIDIA")
    self.gen_add_code_line("#ifndef SMS")
    self.gen_add_code_line("#define SMS GRID_CUBLASDX_SM")
    self.gen_add_code_line("#endif")
    self.gen_add_code_line("namespace glass {")
    self.gen_add_code_line("namespace nvidia {")
    for relative_path in _GLASS_NVIDIA_FILES:
        _emit_glass_source_file(self, relative_path)
    self.gen_add_code_line("// Generated NVIDIA wrapper instantiations used by this robot.")
    for m, n in _nvidia_gemv_sizes(self):
        self.gen_add_code_line("DEFINE_NVIDIA_GEMV(" + str(m) + ", " + str(n) + ")")
    for m, n, k in _nvidia_gemm_sizes(self):
        self.gen_add_code_line("DEFINE_NVIDIA_GEMM(" + str(m) + ", " + str(n) + ", " + str(k) + ")")
    self.gen_add_code_line("} // namespace nvidia")
    self.gen_add_code_line("} // namespace glass")
    self.gen_add_code_line("#endif")
    self.gen_add_code_line("")
    self.gen_add_code_line("namespace " + self.file_namespace + " {", True)
    self.gen_add_func_doc("Compile-time linear algebra backend controls")
    self.gen_add_code_lines([
        "const int GRID_LINALG_GLASS_VALUE = GRID_LINALG_GLASS;",
        "const int GRID_LINALG_GLASS_NVIDIA_VALUE = GRID_LINALG_GLASS_NVIDIA;",
        "const int GRID_CUDA_LINALG_BACKEND_VALUE = GRID_CUDA_LINALG_BACKEND;",
        "const int GRID_CUDA_USE_GLASS_NVIDIA_VALUE = GRID_CUDA_USE_GLASS_NVIDIA;",
        "",
        "template <typename T, int M, int N, int K>",
        "__host__ __device__ constexpr size_t grid_linalg_nvidia_gemm_smem_bytes() {",
        "#if GRID_CUDA_USE_GLASS_NVIDIA",
        "    return glass::nvidia::gemm_smem_size<float, M, N, K>();",
        "#else",
        "    return static_cast<size_t>(0);",
        "#endif",
        "}",
        "",
        "template <typename T, int M, int N>",
        "__host__ __device__ constexpr size_t grid_linalg_nvidia_gemv_smem_bytes() {",
        "#if GRID_CUDA_USE_GLASS_NVIDIA",
        "    return glass::nvidia::gemv_smem_size<float, M, N>();",
        "#else",
        "    return static_cast<size_t>(0);",
        "#endif",
        "}",
        "",
        "template <typename T, int M, int N, int ROW_STRIDE>",
        "__host__ __device__ constexpr size_t grid_linalg_nvidia_row_strided_gemv_smem_bytes() {",
        "#if GRID_CUDA_USE_GLASS_NVIDIA",
        "    return glass::nvidia::row_strided_gemv_smem_size<float, M, N>();",
        "#else",
        "    return static_cast<size_t>(0);",
        "#endif",
        "}",
        "",
        "template <typename T, int M, int N, int K, int A_RS, int B_RS>",
        "__host__ __device__ constexpr size_t grid_linalg_nvidia_row_strided_gemm_smem_bytes() {",
        "#if GRID_CUDA_USE_GLASS_NVIDIA",
        "    return glass::nvidia::row_strided_gemm_smem_size<float, M, N, K>();",
        "#else",
        "    return static_cast<size_t>(0);",
        "#endif",
        "}",
        "",
        "template <typename T>",
        "__host__ __device__ constexpr size_t GRID_LINALG_NVIDIA_MAX_HELPER_BYTES() {",
        "#if GRID_CUDA_USE_GLASS_NVIDIA",
        "    size_t bytes = grid_linalg_nvidia_gemm_smem_bytes<T, 4, 4, 4>();",
        "    size_t b66 = grid_linalg_nvidia_gemm_smem_bytes<T, 6, 6, 6>();",
        "    size_t b661 = grid_linalg_nvidia_gemm_smem_bytes<T, 6, 6, 1>();",
        "    size_t bgemv66 = grid_linalg_nvidia_gemv_smem_bytes<T, 6, 6>();",
        "    size_t brsgemv666 = grid_linalg_nvidia_row_strided_gemv_smem_bytes<T, 6, 6, 6>();",
        "    bytes = bytes > b66 ? bytes : b66;",
        "    bytes = bytes > b661 ? bytes : b661;",
        "    bytes = bytes > bgemv66 ? bytes : bgemv66;",
        "    bytes = bytes > brsgemv666 ? bytes : brsgemv666;",
        "    return bytes;",
        "#else",
        "    return static_cast<size_t>(0);",
        "#endif",
        "}",
        "",
        "template <typename T, int M, int N, int K, bool TRANSPOSE_B = false, bool ROW_MAJOR_A = false, bool ROW_MAJOR_B = false, bool ROW_MAJOR_C = false>",
        "__device__ void grid_linalg_gemm_glass(const T *A, const T *B, T *C, T alpha, T beta) {",
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
        "template <typename T, int M, int N, int K, bool TRANSPOSE_B = false, bool ROW_MAJOR_A = false, bool ROW_MAJOR_B = false, bool ROW_MAJOR_C = false>",
        "__device__ void grid_linalg_gemm_default(const T *A, const T *B, T *C, T alpha, T beta) {",
        "    grid_linalg_gemm_glass<T, M, N, K, TRANSPOSE_B, ROW_MAJOR_A, ROW_MAJOR_B, ROW_MAJOR_C>(A, B, C, alpha, beta);",
        "}",
        "",
        "#if GRID_CUDA_USE_GLASS_NVIDIA",
        "template <typename T, int M, int N, int K>",
        "__device__ void grid_linalg_packed_gemm_nvidia_colmajor(const T *A, const T *B, T *C, T alpha, T beta, unsigned char *smem) {",
        "    static_assert(sizeof(T) == sizeof(float), \"glass-nvidia backend currently supports float only\");",
        "    glass::nvidia::gemm<float, M, N, K>(static_cast<float>(alpha), reinterpret_cast<float *>(const_cast<T *>(A)), reinterpret_cast<float *>(const_cast<T *>(B)), static_cast<float>(beta), reinterpret_cast<float *>(C), reinterpret_cast<char *>(smem));",
        "}",
        "#endif",
        "",
        "#if GRID_CUDA_USE_GLASS_NVIDIA",
        "template <typename T, int M, int N, int ROW_STRIDE>",
        "__device__ void grid_linalg_row_strided_gemv_nvidia(const T *A, const T *x, T *y, T alpha, T beta, unsigned char *smem) {",
        "    static_assert(sizeof(T) == sizeof(float), \"glass-nvidia row-strided GEMV currently supports float only\");",
        "    ::glass::nvidia::row_strided_gemv<float, M, N, ROW_STRIDE>(static_cast<float>(alpha), reinterpret_cast<float *>(const_cast<T *>(A)), reinterpret_cast<float *>(const_cast<T *>(x)), static_cast<float>(beta), reinterpret_cast<float *>(y), reinterpret_cast<char *>(smem));",
        "}",
        "",
        "template <typename T, int M, int N, int K, int A_RS, int B_RS>",
        "__device__ void grid_linalg_row_strided_gemm_nvidia(const T *A, const T *B, T *C, T alpha, T beta, unsigned char *smem) {",
        "    static_assert(sizeof(T) == sizeof(float), \"glass-nvidia row-strided GEMM currently supports float only\");",
        "    ::glass::nvidia::row_strided_gemm<float, M, N, K, A_RS, B_RS>(static_cast<float>(alpha), reinterpret_cast<float *>(const_cast<T *>(A)), reinterpret_cast<float *>(const_cast<T *>(B)), static_cast<float>(beta), reinterpret_cast<float *>(C), reinterpret_cast<char *>(smem));",
        "}",
        "#endif",
        "",
        "template <typename T, int M, int N, int K, bool TRANSPOSE_B = false, bool ROW_MAJOR_A = false, bool ROW_MAJOR_B = false, bool ROW_MAJOR_C = false>",
        "__device__ void grid_linalg_gemm(const T *A, const T *B, T *C, T alpha, T beta, unsigned char *glass_nvidia_smem = nullptr) {",
        "#if GRID_CUDA_USE_GLASS_NVIDIA",
        "    if (glass_nvidia_smem != nullptr && !TRANSPOSE_B && !ROW_MAJOR_A && !ROW_MAJOR_B && !ROW_MAJOR_C) {",
        "        grid_linalg_packed_gemm_nvidia_colmajor<T, M, N, K>(A, B, C, alpha, beta, glass_nvidia_smem);",
        "        __syncthreads();",
        "        return;",
        "    }",
        "#endif",
        "    grid_linalg_gemm_glass<T, M, N, K, TRANSPOSE_B, ROW_MAJOR_A, ROW_MAJOR_B, ROW_MAJOR_C>(A, B, C, alpha, beta);",
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
        "#if GRID_CUDA_USE_GLASS_NVIDIA",
        "    if (glass_nvidia_smem != nullptr) {",
        "        grid_linalg_row_strided_gemv_nvidia<T, M, N, ROW_STRIDE>(A, x, y, alpha, beta, glass_nvidia_smem);",
        "        __syncthreads();",
        "        return;",
        "    }",
        "#endif",
        "    ::glass::row_strided_gemv<T, M, N, ROW_STRIDE>(A, x, y, alpha, beta);",
        "    __syncthreads();",
        "}",
        "",
        "template <typename T, int M, int N, int K, int A_RS, int B_RS>",
        "__device__ void grid_linalg_row_strided_gemm(const T *A, const T *B, T *C, T alpha, T beta) {",
        "    ::glass::row_strided_gemm<T, M, N, K, A_RS, B_RS>(A, B, C, alpha, beta);",
        "    __syncthreads();",
        "}",
        "",
        "template <typename T, int N, int S1, int S2>",
        "__device__ T grid_linalg_dot_strided(const T *vec1, const T *vec2) {",
        "    return ::glass::dot_strided<T, N, S1, S2>(vec1, vec2);",
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
