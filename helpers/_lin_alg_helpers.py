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

# Files that carry their own `namespace glass { namespace nvidia { ... } }`
# wrapping (or other top-level scope) and must be vendored at global scope so
# we don't end up with nested glass::nvidia::glass::nvidia::... symbols.
_GLASS_NVIDIA_GLOBAL_SCOPE_FILES = [
    # tuning_table.cuh defines `namespace _glass_tuning { cublasdx_wins<>; ... }`
    # which the SIMT query (query_simt.cuh) looks up when deciding whether to
    # dispatch to cuBLASDx or fall through to the SIMT path. Must come before
    # any file (l3.cuh, l3_simt.cuh) that consumes `should_use_cublasdx<>`.
    "src/nvidia/tuning_table.cuh",
    # types.cuh defines glass::nvidia::layout (with its own namespace block)
    # plus the private helper macros (_GLASS_CUBLAS_LAYOUT, _GLASS_ASSERT_BLOCKDIM_GEQ).
    "src/nvidia/types.cuh",
]

# Files that contain bare contents intended to be included INSIDE the
# `namespace glass::nvidia { ... }` block (i.e., they don't open their own).
#
# Order matches glass-nvidia.cuh: query_simt.cuh declares the
# should_use_cublasdx_{gemv,row_strided_gemv,row_strided_gemm,batched}<>
# helpers that l2.cuh / l3.cuh / l3_simt.cuh consume, so it must come first.
_GLASS_NVIDIA_FILES = [
    "src/nvidia/sizes.cuh",
    "src/nvidia/l1.cuh",
    "src/nvidia/query_simt.cuh",
    # l3_simt.cuh — gemm_batched_1d / gemm_strided_batched_1d (1D-launch
    # batched APIs used by the eepose codegen for parent-shared chain GEMMs).
    "src/nvidia/l3_simt.cuh",
    "src/nvidia/l2.cuh",
    "src/nvidia/l3.cuh",
    # query.cuh — gemm_min_block_threads / gemm_block_threads_valid for
    # SUGGESTED_THREADS static_asserts.
    "src/nvidia/query.cuh",
]

# cuSOLVERDx-backed LAPACK wrappers. Vendored separately because they require
# -rdc=true -dlto -lcusolverdx -lcublas -lcusolver -lcudart at link time (the
# rest of glass-nvidia is header-only) and are guarded by a separate macro
# (GRID_CUDA_USE_GLASS_NVIDIA_LAPACK).
_GLASS_NVIDIA_LAPACK_FILES = [
    "src/nvidia/lapack.cuh",
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


def _cublasdx_wins_heuristic(m, n, k=None):
    """Mirror of `_glass_tuning::cublasdx_wins<>`'s primary-template
    heuristic: `max(M, N, K) >= 16 AND min(M, N, K) >= 4`.

    Used to gate `DEFINE_NVIDIA_GEMM_BLOCKDIM_SM` emission. Emitting the
    DEFINE creates an explicit cuBLASDx specialization that bypasses
    GLASS's compile-time auto-dispatch, so we only emit it for shapes
    where cuBLASDx is expected to win. Per-SM measurements in
    tuning_table.cuh may override the heuristic; we don't parse that
    file here.
    """
    dims = [d for d in (m, n, k) if d is not None]
    return max(dims) >= 16 and min(dims) >= 4


def _nvidia_gemm_sizes(self):
    """Shapes for which we emit DEFINE_NVIDIA_GEMM_BLOCKDIM_SM macros.

    Filtered by `_cublasdx_wins_heuristic(M,N,K)` so SIMT-winning shapes
    (e.g. (4,4,4), (6,6,6)) skip the DEFINE and let GLASS's primary
    template auto-dispatch to SIMT. Forcing cuBLASDx for a SIMT-winning
    shape both regresses perf and (historically, for (6,6,6)) caused
    cudaErrorIllegalAddress on glass-nvidia builds.
    """
    candidates = {(4, 4, 4)}
    for size in _nvidia_row_strided_gemm_sizes(self):
        candidates.add(size)
    return sorted(s for s in candidates if _cublasdx_wins_heuristic(*s))


def _nvidia_gemv_sizes(self):
    return sorted(
        s for s in set(_nvidia_row_strided_gemv_sizes(self))
        if _cublasdx_wins_heuristic(*s)
    )


def _nvidia_row_strided_gemv_sizes(self):
    return [(6, 6)]


def _nvidia_row_strided_gemm_sizes(self):
    return [(6, 6, 6)]


def gen_linalg_smem_setup(self, temp_size):
    """Emit a local `unsigned char *s_linalg_smem` set up to point past
    `s_temp[temp_size]` and re-aligned to 16 bytes (cuBLASDx's smem
    alignment requirement). Without the re-align, the offset can be up
    to 12 bytes off depending on temp_size + s_temp alignment, silently
    corrupting cuBLASDx's per-tile smem. On glass-only builds the
    pointer is nullptr (unused).
    """
    self.gen_add_code_line("#if GRID_CUDA_USE_GLASS_NVIDIA")
    self.gen_add_code_line(
        f"unsigned char *s_linalg_smem = reinterpret_cast<unsigned char *>("
        f"(reinterpret_cast<uintptr_t>(s_temp + {temp_size}) + 15u) & ~uintptr_t(15));"
    )
    self.gen_add_code_line("#else")
    self.gen_add_code_line("unsigned char *s_linalg_smem = nullptr;")
    self.gen_add_code_line("#endif")


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
    # Every DEFINE_NVIDIA_* macro and static_assert below passes
    # GRID_CUBLASDX_SM explicitly via the _SM variants, so SM is
    # per-instantiation rather than a single global #define — enables
    # multi-arch CUBIN builds (one header → many sm_xx).
    # Files that bring their own glass::nvidia namespace go at global scope.
    for relative_path in _GLASS_NVIDIA_GLOBAL_SCOPE_FILES:
        _emit_glass_source_file(self, relative_path)
    self.gen_add_code_line("namespace glass {")
    self.gen_add_code_line("namespace nvidia {")
    for relative_path in _GLASS_NVIDIA_FILES:
        _emit_glass_source_file(self, relative_path)
    # Pin cuBLASDx's BlockDim<TC,1,1> to this robot's SUGGESTED_THREADS so
    # callers can launch with SUGGESTED_THREADS without deadlocking inside
    # cuBLASDx's warp-level barriers.
    tc = self.suggested_threads
    self.gen_add_code_line("// Generated NVIDIA wrapper instantiations used by this robot.")
    self.gen_add_code_line(
        "// BLOCK_THREADS pinned to SUGGESTED_THREADS=" + str(tc)
        + " so kernels launched at that thread count can call cuBLASDx without"
    )
    self.gen_add_code_line(
        "// thread-count-mismatch deadlocks. Extra threads (if any) go idle"
        " inside execute() per cuBLASDx example 04_gemm_blockdim."
    )
    for m, n in _nvidia_gemv_sizes(self):
        self.gen_add_code_line(
            "static_assert(::glass::nvidia::gemv_block_threads_valid<float, "
            + str(m) + ", " + str(n) + ", " + str(tc) + ", GRID_CUBLASDX_SM>(),"
        )
        self.gen_add_code_line(
            '              "SUGGESTED_THREADS=' + str(tc)
            + ' is too small for cuBLASDx gemv<' + str(m) + ',' + str(n)
            + '> on this SM");'
        )
        self.gen_add_code_line(
            "DEFINE_NVIDIA_GEMV_BLOCKDIM_SM(" + str(m) + ", " + str(n)
            + ", " + str(tc) + ", GRID_CUBLASDX_SM)"
        )
    for m, n, k in _nvidia_gemm_sizes(self):
        self.gen_add_code_line(
            "static_assert(::glass::nvidia::gemm_block_threads_valid<float, "
            + str(m) + ", " + str(n) + ", " + str(k) + ", " + str(tc)
            + ", GRID_CUBLASDX_SM>(),"
        )
        self.gen_add_code_line(
            '              "SUGGESTED_THREADS=' + str(tc)
            + ' is too small for cuBLASDx gemm<' + str(m) + ',' + str(n) + ',' + str(k)
            + '> on this SM");'
        )
        # Plain col-major DEFINE only. Small TRANSB call sites (e.g. 6×6 in
        # _aba/_crba) auto-route to SIMT via glass::nvidia::gemm<>'s SIMT
        # branch, which honors (col, row, col) → TRANSPOSE_B=true natively
        # — no cuBLASDx specialization needed. If a larger cuBLASDx-routed
        # TRANSB shape gets added later, also emit
        # DEFINE_NVIDIA_GEMM_BLOCKDIM_TRANSB_SM for it, otherwise the
        # primary template will static_assert.
        self.gen_add_code_line(
            "DEFINE_NVIDIA_GEMM_BLOCKDIM_SM(" + str(m) + ", " + str(n) + ", " + str(k)
            + ", " + str(tc) + ", GRID_CUBLASDX_SM)"
        )
    self.gen_add_code_line("} // namespace nvidia")
    self.gen_add_code_line("} // namespace glass")
    self.gen_add_code_line("#endif")
    self.gen_add_code_line("")
    # cuSOLVERDx-backed LAPACK wrappers (chol_inplace, trsm). Vendored under a
    # separate macro guard so users who don't link cuSOLVERDx (which requires
    # -rdc=true -dlto -lcusolverdx -lcublas -lcusolver -lcudart) aren't forced
    # to. Algorithms don't wire posv/chol/trsm yet — Phase 4 territory. The
    # infrastructure is here so a .cu can `DEFINE_NVIDIA_CHOL_BLOCKDIM(N, TC)`
    # and call `glass::nvidia::chol_inplace<T, N, TC>(...)` directly.
    self.gen_add_code_line("#if GRID_CUDA_USE_GLASS_NVIDIA_LAPACK")
    self.gen_add_code_line("namespace glass {")
    self.gen_add_code_line("namespace nvidia {")
    for relative_path in _GLASS_NVIDIA_LAPACK_FILES:
        _emit_glass_source_file(self, relative_path)
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
        "// These smem-size queries must match the BLOCK_THREADS we pin via",
        "// DEFINE_NVIDIA_*_BLOCKDIM, otherwise GRID_LINALG_NVIDIA_MAX_HELPER_BYTES",
        "// returns the primary-template 0 and the kernel under-allocates smem.",
        "template <typename T, int M, int N, int K>",
        "__host__ __device__ constexpr size_t grid_linalg_nvidia_gemm_smem_bytes() {",
        "#if GRID_CUDA_USE_GLASS_NVIDIA",
        "    return glass::nvidia::gemm_smem_size<float, M, N, K, SUGGESTED_THREADS>();",
        "#else",
        "    return static_cast<size_t>(0);",
        "#endif",
        "}",
        "",
        "template <typename T, int M, int N>",
        "__host__ __device__ constexpr size_t grid_linalg_nvidia_gemv_smem_bytes() {",
        "#if GRID_CUDA_USE_GLASS_NVIDIA",
        "    return glass::nvidia::gemv_smem_size<float, M, N, SUGGESTED_THREADS>();",
        "#else",
        "    return static_cast<size_t>(0);",
        "#endif",
        "}",
        "",
        "// Returns 0 when the shape auto-dispatches to SIMT (no scratch needed).",
        "template <typename T, int M, int N, int ROW_STRIDE>",
        "__host__ __device__ constexpr size_t grid_linalg_nvidia_row_strided_gemv_smem_bytes() {",
        "#if GRID_CUDA_USE_GLASS_NVIDIA",
        "    return glass::nvidia::row_strided_gemv_smem_size<float, M, N, ROW_STRIDE, SUGGESTED_THREADS>();",
        "#else",
        "    return static_cast<size_t>(0);",
        "#endif",
        "}",
        "",
        "template <typename T, int M, int N, int K, int A_RS, int B_RS>",
        "__host__ __device__ constexpr size_t grid_linalg_nvidia_row_strided_gemm_smem_bytes() {",
        "#if GRID_CUDA_USE_GLASS_NVIDIA",
        "    return glass::nvidia::row_strided_gemm_smem_size<float, M, N, K, A_RS, B_RS, SUGGESTED_THREADS>();",
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
        "// Public linalg wrappers below route uniform col-major calls through",
        "// glass::nvidia::*, which auto-dispatches between cuBLASDx and SIMT at",
        "// compile time via should_use_cublasdx*<>. Mixed-layout / transposed",
        "// paths fall through to glass::*_ex (SIMT). The `glass_nvidia_smem`",
        "// param is the cuBLASDx scratch pointer; unused on glass-only builds.",
        "template <typename T, int M, int N, int K, bool TRANSPOSE_B = false, bool ROW_MAJOR_A = false, bool ROW_MAJOR_B = false, bool ROW_MAJOR_C = false>",
        "__device__ void grid_linalg_gemm(const T *A, const T *B, T *C, T alpha, T beta, unsigned char *glass_nvidia_smem = nullptr) {",
        "#if GRID_CUDA_USE_GLASS_NVIDIA",
        "    if constexpr (!ROW_MAJOR_A && !ROW_MAJOR_B && !ROW_MAJOR_C) {",
        "        static_assert(sizeof(T) == sizeof(float), \"glass-nvidia backend currently supports float only\");",
        "        // TRANSPOSE_B is expressed as LB=row_major; the SIMT branch in",
        "        // glass::nvidia::gemm<> maps (col, row, col) → TRANSPOSE_B=true.",
        "        constexpr auto _LB = TRANSPOSE_B ? ::glass::nvidia::layout::row_major : ::glass::nvidia::layout::col_major;",
        "        ::glass::nvidia::gemm<float, M, N, K, SUGGESTED_THREADS, ::glass::nvidia::layout::col_major, _LB, ::glass::nvidia::layout::col_major, GRID_CUBLASDX_SM>(static_cast<float>(alpha), reinterpret_cast<float *>(const_cast<T *>(A)), reinterpret_cast<float *>(const_cast<T *>(B)), static_cast<float>(beta), reinterpret_cast<float *>(C), reinterpret_cast<char *>(glass_nvidia_smem));",
        "        __syncthreads();",
        "        return;",
        "    }",
        "#endif",
        "    (void)glass_nvidia_smem;",
        "    grid_linalg_gemm_glass<T, M, N, K, TRANSPOSE_B, ROW_MAJOR_A, ROW_MAJOR_B, ROW_MAJOR_C>(A, B, C, alpha, beta);",
        "}",
        "",
        "template <typename T, int M, int N, bool TRANSPOSE = false, bool ROW_MAJOR_A = false>",
        "__device__ void grid_linalg_gemv(const T *A, const T *x, T *y, T alpha, T beta, unsigned char *glass_nvidia_smem = nullptr) {",
        "#if GRID_CUDA_USE_GLASS_NVIDIA",
        "    if constexpr (!TRANSPOSE && !ROW_MAJOR_A) {",
        "        static_assert(sizeof(T) == sizeof(float), \"glass-nvidia gemv currently supports float only\");",
        "        ::glass::nvidia::gemv<float, M, N, SUGGESTED_THREADS, ::glass::nvidia::layout::col_major, ::glass::nvidia::layout::col_major, ::glass::nvidia::layout::col_major, GRID_CUBLASDX_SM>(static_cast<float>(alpha), reinterpret_cast<float *>(const_cast<T *>(A)), reinterpret_cast<float *>(const_cast<T *>(x)), static_cast<float>(beta), reinterpret_cast<float *>(y), reinterpret_cast<char *>(glass_nvidia_smem));",
        "        __syncthreads();",
        "        return;",
        "    }",
        "#endif",
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
        "    static_assert(sizeof(T) == sizeof(float), \"glass-nvidia row-strided GEMV currently supports float only\");",
        "    ::glass::nvidia::row_strided_gemv<float, M, N, ROW_STRIDE, SUGGESTED_THREADS, ::glass::nvidia::layout::col_major, ::glass::nvidia::layout::col_major, ::glass::nvidia::layout::col_major, GRID_CUBLASDX_SM>(static_cast<float>(alpha), reinterpret_cast<float *>(const_cast<T *>(A)), reinterpret_cast<float *>(const_cast<T *>(x)), static_cast<float>(beta), reinterpret_cast<float *>(y), reinterpret_cast<char *>(glass_nvidia_smem));",
        "    __syncthreads();",
        "#else",
        "    (void)glass_nvidia_smem;",
        "    ::glass::row_strided_gemv<T, M, N, ROW_STRIDE>(A, x, y, alpha, beta);",
        "    __syncthreads();",
        "#endif",
        "}",
        "",
        "template <typename T, int M, int N, int K, int A_RS, int B_RS>",
        "__device__ void grid_linalg_row_strided_gemm(const T *A, const T *B, T *C, T alpha, T beta, unsigned char *glass_nvidia_smem = nullptr) {",
        "#if GRID_CUDA_USE_GLASS_NVIDIA",
        "    static_assert(sizeof(T) == sizeof(float), \"glass-nvidia row-strided GEMM currently supports float only\");",
        "    ::glass::nvidia::row_strided_gemm<float, M, N, K, A_RS, B_RS, SUGGESTED_THREADS, ::glass::nvidia::layout::col_major, ::glass::nvidia::layout::col_major, ::glass::nvidia::layout::col_major, GRID_CUBLASDX_SM>(static_cast<float>(alpha), reinterpret_cast<float *>(const_cast<T *>(A)), reinterpret_cast<float *>(const_cast<T *>(B)), static_cast<float>(beta), reinterpret_cast<float *>(C), reinterpret_cast<char *>(glass_nvidia_smem));",
        "    __syncthreads();",
        "#else",
        "    (void)glass_nvidia_smem;",
        "    ::glass::row_strided_gemm<T, M, N, K, A_RS, B_RS>(A, B, C, alpha, beta);",
        "    __syncthreads();",
        "#endif",
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
