MEMORY_THRESHOLD = 8 # Max num joints for shared mem allocation of result


def gen_fdsva_so_contract(self):
	# construct the boilerplate and function definition
    n = self.robot.get_num_vel()
    inner_arena_size = 4 * n**3
    func_params = ["s_df2 are the second derivatives of forward dynamics WRT q,qd,tau", \
                "s_idsva_so are the second derivative tensors of inverse dynamics", \
                "s_Minv is the inverse mass matrix", \
                "s_df_du is the gradient of the forward dynamics", \
                "s_temp is the (shared) scratch buffer; size FDSVA_SO_INNER_SMEM_BYTES<T, SCRATCH_IN_SMEM>() bytes (= " + str(inner_arena_size) + "*sizeof(T) when SCRATCH_IN_SMEM, else 0)", \
                "d_workspace is the global scratch buffer; size FDSVA_SO_INNER_WORKSPACE_BYTES<T, SCRATCH_IN_SMEM>() bytes (= " + str(inner_arena_size) + "*sizeof(T) when !SCRATCH_IN_SMEM, else 0). Pass nullptr when SCRATCH_IN_SMEM", \
                "gravity is the gravity constant"]
    func_def_start = "void fdsva_so_contract("
    func_def_middle = "T *s_df2, T *s_idsva_so, T *s_Minv, T *s_df_du, "
    func_def_end = "T *s_temp, T *d_workspace, const T gravity) {"
    func_notes = ["Assumes works with IDSVA",
                  "Inline-CUDA users: SCRATCH_IN_SMEM selects where the 4*NV^3 scratch arena lives.",
                  "  true  -> s_temp (shared memory; fastest, current default).",
                  "  false -> d_workspace (global; frees shared memory for the caller's outer kernel).",
                  "The placement is the INNER's choice (made at the top of this function) so the",
                  "kernel/device caller just sizes both arenas from the exposed *_BYTES constants",
                  "and hands both pointers in. Codegen maps each RESOURCE_TIER to a SCRATCH_IN_SMEM",
                  "value per robot (small robots keep all tiers in smem; large robots spill)."]
    func_def_middle, func_params = self.gen_insert_helpers_func_def_params(func_def_middle, func_params, -3)
    func_def = func_def_start + func_def_middle + func_def_end
    self.gen_add_func_doc("Second Order of Forward Dynamics with Spatial Vector Algebra", func_notes, func_params, None)
    self.gen_add_code_line("template <typename T, bool SCRATCH_IN_SMEM = true>")
    self.gen_add_code_line("__device__")
    self.gen_add_code_line(func_def, True)

    n = self.robot.get_num_vel()
    NV = self.robot.get_num_vel()

    self.gen_add_code_line('// Second Derivatives of Inverse Dynamics')
    self.gen_add_code_line("T *d2tau_dqdq = &s_idsva_so[" + str(n*n*n*0) + "];" )
    self.gen_add_code_line("T *d2tau_dvdv = &s_idsva_so[" + str(n*n*n*1) + "];" )
    self.gen_add_code_line("T *d2tau_dvdq = &s_idsva_so[" + str(n*n*n*2) + "];" )
    self.gen_add_code_line("T *dM_dq = &s_idsva_so[" + str(n*n*n*3) + "];" )
    self.gen_add_code_line('\n\n')
    self.gen_add_code_line('// First Derivatives of Forward Dynamics')
    self.gen_add_code_line("T *s_df_dq = s_df_du; T *s_df_dqd = &s_df_du[" + str(n*n) + "];")
    self.gen_add_code_line('\n\n')
    self.gen_add_code_line('// Second Derivatives of Forward Dynamics')
    self.gen_add_code_line('T *d2a_dqdq = s_df2;')
    self.gen_add_code_line('T *d2a_dvdq = &s_df2[' + str(n*n*n) + '];')
    self.gen_add_code_line('T *d2a_dvdv = &s_df2[' + str(2*n*n*n) + '];')
    self.gen_add_code_line('T *d2a_dtdq = &s_df2[' + str(3*n*n*n) + '];')
    self.gen_add_code_line('\n\n')
    self.gen_add_code_line('// Temporary Variables. The 4*n^3 scratch arena placement is the INNER\'s')
    self.gen_add_code_line('// choice, keyed on SCRATCH_IN_SMEM: s_temp (shared) when true, d_workspace')
    self.gen_add_code_line('// (global) when false. A surgical-spill change is local here: repoint a')
    self.gen_add_code_line('// sub-buffer + update FDSVA_SO_INNER_{SMEM,WORKSPACE}_BYTES.')
    self.gen_add_code_line('T *inner_arena;')
    self.gen_add_code_line('if constexpr (SCRATCH_IN_SMEM) {', True)
    self.gen_add_code_line('(void)d_workspace;')
    self.gen_add_code_line('inner_arena = s_temp;')
    self.gen_add_end_control_flow()
    self.gen_add_code_line('else {', True)
    self.gen_add_code_line('(void)s_temp;')
    self.gen_add_code_line('inner_arena = d_workspace;')
    self.gen_add_end_control_flow()
    self.gen_add_code_line(f'T *inner_dq = inner_arena; // Inner term for d2a_dqdq')
    self.gen_add_code_line(f'T *inner_cross = inner_dq + {n**3}; // Inner term for d2a_dvdq (dM_dq*Minv)')
    self.gen_add_code_line(f'T *inner_tau = inner_cross + {n**3}; // Inner term for d2a_dtdq (d2tau_dvdq + dM_dq*da_dv)')
    self.gen_add_code_line(f'T *rot_dq = inner_tau + {n**3}; // Rotated (dM_dq*da_dq)^R term used to compute inner_dq')
    self.gen_add_code_line(f'\n\n')

    # Start inner term for d2a_dqdq
    self.gen_add_code_line('// Start inner term for d2a_dqdq & Fill out Minv')
    self.gen_add_parallel_loop("ind",str(n**3 + n*n))
    self.gen_add_code_line(f'int i = ind / {n*n} % {n}; int j = ind / {n} % {n}; int k = ind % {n};')
    self.gen_add_code_line(f'if (ind < {n**3}) {{', True)
    self.gen_add_code_line(f'inner_dq[ind] = dot_prod<T, {n}, {n}, 1>(&dM_dq[{n*n}*i + k], &s_df_dq[{n}*j]);')
    self.gen_add_code_line(f'rot_dq[i*{n*n} + k*{n} + j] = inner_dq[ind];')
    self.gen_add_end_control_flow()
    self.gen_add_code_line(f'else if (k > j) s_Minv[j*{n} + k] = s_Minv[k*{n} + j];')
    self.gen_add_end_control_flow()
    self.gen_add_sync()
    
    # 3Dx2D Tensor Computation defined as iLk,Lj->ijk
    self.gen_add_code_line('// Compute relevant inner subterms in parallel')
    self.gen_add_parallel_loop("ind",str(3*n**3))
    self.gen_add_code_line(f'int i = ind / {n*n} % {n}; int j = ind / {n} % {n}; int k = ind % {n};')
    self.gen_add_code_line(f'if (ind < {n**3}) inner_dq[ind] += rot_dq[ind] + d2tau_dqdq[ind]; // Started with dM_dq*da_dq')
    self.gen_add_code_line(f'else if (ind < {2*n**3}) inner_cross[i*{n*n} + k*{n} + j] = dot_prod<T, {n}, {n}, 1>(&dM_dq[{n*n}*i + k], &s_df_dqd[{n}*j]) + d2tau_dvdq[i*{n*n} + j*{n} + k];')
    self.gen_add_code_line(f'else inner_tau[i*{n*n} + k*{n} + j] = dot_prod<T, {n}, {n}, 1>(&dM_dq[{n*n}*i + k], &s_Minv[{n}*j]);')
    self.gen_add_end_control_flow()
    self.gen_add_sync()

    # 2Dx3D tensor computation defined as iL,Ljk->ijk
    self.gen_add_code_line('// Multiply by -Minv to finish algorithm')
    # PERF EXPERIMENT CANDIDATE (2026-05-18): this 4*n^3 parallel_loop has the
    # highest FMA count in fdsva_so_contract (g1_floating: ~6M FMAs out of ~9M
    # total). Each thread does a serial n-element dot_prod, so total FMAs ~
    # 4*n^4. Structurally an iL,Ljk->ijk tensor contraction — could be
    # rewritten as 4*n batched n×n×n gemms or 4 large n × n² × n gemms.
    #
    # Standalone autotune on sm_120 says cuBLASDx wins 2.4-5.2× at 24-48
    # size range, BUT standalone-gemm wins do NOT directly translate into
    # hot-kernel wins. In-kernel concerns: register pressure with all SO
    # state live, shared-mem layout cost (current iL,Ljk strides aren't
    # gemm-friendly), tier pressure (g1_floating already in spill tier and
    # cuBLASDx wants more shared-mem), per-block-per-timestep sync overhead.
    #
    # First step before any refactor: profile this loop in-context (ptxas
    # occupancy / nsight compute trace) to confirm it really is the
    # bottleneck — could equally be memory-bandwidth or sync-bound, in
    # which case cuBLASDx won't help. Then prototype on a single robot
    # before generalizing.
    #
    # COALESCED-DOT EXPERIMENT (2026-05-26, EVALUATED, NOT TAKEN): the spilled
    # operand (inner_dq/inner_cross/inner_tau/d2tau_dvdv, in d_workspace when
    # CONTRACT_IN_SMEM=false) is gathered with stride n^2 over the contracted
    # axis L — inner_dq[j + k*n + L*n^2] reads the SAME buffer that was written
    # as [i][j][k], so L (the sum axis) maps to the buffer's slowest, stride-n^2
    # index. grid_linalg_dot_strided_coalesced coalesces ALONG the contraction
    # stride (consecutive thread ranks read x[rank*SX], y[rank*SY]), so applying
    # it here would issue stride-n^2 loads across the warp — WORSE than the
    # current stride-n warp gather, not better; the primitive only helps when the
    # summed axis is the contiguous (small-stride) one, which it is not in this
    # iL,Ljk->ijk layout. Worse still, it is BLOCK-COOPERATIVE (one scalar per
    # call): producing 4*n^3 outputs would mean 4*n^3 sequential block-reductions
    # (g1: ~171.5k dots, ~343k __syncthreads), each n-element dot using only
    # ~n/blockDim of the threads — collapsing the current embarrassingly-parallel
    # 4*n^3-way thread parallelism. Making the primitive coalesce would require
    # transposing inner_* so L is contiguous (an extra full global read+write of
    # the 4*n^3 tensor) AND still pay the per-output block-reduction serialization.
    # Both effects are large regressions, so the contraction is LEFT AS-IS (the
    # smem path is already fine; the spilled path's stride-n^2 gather is the cost
    # of spilling, not something this primitive fixes). Profiled/reasoned and
    # rejected — do not re-attempt with this primitive without a layout that puts
    # the contracted axis contiguous AND avoids one-dot-per-block serialization.
    self.gen_add_parallel_loop("ind",str(4*n**3))
    self.gen_add_code_line(f'int i = ind / {n*n} % {n}; int j = ind / {n} % {n}; int k = ind % {n};')
    self.gen_add_code_line(f'if (ind < {n**3}) d2a_dqdq[i*{n*n} + j*{n} + k] = -dot_prod<T, {n}, {n}, {n*n}>(&s_Minv[i], &inner_dq[j + k*{n}]);')
    self.gen_add_code_line(f'else if (ind < {2*n**3}) d2a_dvdq[i*{n*n} + j*{n} + k] = -dot_prod<T, {n}, {n}, {n*n}>(&s_Minv[i], &inner_cross[j + k*{n}]);')
    self.gen_add_code_line(f'else if (ind < {3*n**3}) d2a_dvdv[i*{n*n} + j*{n} + k] = -dot_prod<T, {n}, {n}, {n*n}>(&s_Minv[i], &d2tau_dvdv[j + k*{n}]);')
    self.gen_add_code_line(f'else d2a_dtdq[i*{n*n} + j*{n} + k] = -dot_prod<T, {n}, {n}, {n*n}>(&s_Minv[i], &inner_tau[j + k*{n}]);')
    self.gen_add_end_control_flow()
    self.gen_add_sync()

    self.gen_add_end_function()

def gen_fdsva_so_contract_temp_mem_size(self):
    n = self.robot.get_num_vel()
    return 4*n**3

def gen_fdsva_so_fd_gradient_inline_temp_mem_size(self):
    n = self.robot.get_num_vel()
    return 18*self.robot.get_num_joints() + 2*n*n + self.gen_inverse_dynamics_gradient_inner_temp_mem_size()


def gen_fdsva_so_fd_gradient_inline_temp_mem_size_spilled(self):
    """MEM1 variant: shared-mem footprint when id_du_gradient_inner uses spill.

    Same prefix layout (s_fd_vaf + s_fd_dc_du) but s_fd_temp is sized for
    `selective_shared_count` (= full minus the spilled da_dq..fxvi band)
    instead of full. inverse_dynamics_inner_vaf also needs s_fd_temp but
    its temp size (6*NUM_POS) is tiny — max() handles it.
    """
    n = self.robot.get_num_vel()
    layout = self.gen_inverse_dynamics_gradient_temp_layout()
    id_inner_vaf_temp = self.gen_inverse_dynamics_inner_temp_mem_size()
    s_fd_temp_size = max(id_inner_vaf_temp, layout["selective_shared_count"])
    return 18*self.robot.get_num_joints() + 2*n*n + s_fd_temp_size

def gen_fdsva_so_fd_gradient_inline(self, use_spill = False, spill_ptr_expr = "nullptr"):
    """Emit the inline FD-gradient computation.

    MEM1 (use_spill=True path): when the kernel is shared-mem-pressured (big
    floating-base robots like g1), have the inner id_du_gradient_inner spill
    its da_dq..fxvi block to ``spill_ptr_expr`` (typically a slice of
    d_workspace). This matches the well-tested spill pattern used by
    id_du_kernel — the WHOLE temp is NOT pushed to global; only the
    overflow band. s_fd_vaf/s_fd_dc_du/s_fd_temp still live in shared
    `s_temp`, but s_fd_temp is sized for ``selective_shared_count``, not
    ``full_count``. Saves the bulk of fd_grad_inline's ~40k-float footprint
    for NV=35 robots.
    """
    n = self.robot.get_num_vel()
    self.gen_add_code_line("// Compute FD gradient inline so nested device wrappers do not carve a second shared arena")
    self.gen_add_code_line("T *s_fd_vaf = s_temp;")
    self.gen_add_code_line(f"T *s_fd_dc_du = s_fd_vaf + {18*self.robot.get_num_joints()};")
    self.gen_add_code_line(f"T *s_fd_temp = s_fd_dc_du + {2*n*n};")
    self.gen_inverse_dynamics_inner_function_call(
        compute_c = False,
        use_qdd_input = True,
        updated_var_names = dict(
            s_vaf_name = "s_fd_vaf",
            s_q_name = "s_q",
            s_qd_name = "s_qd",
            s_qdd_name = "s_qdd",
            s_temp_name = "s_fd_temp",
            d_f_ext_name = "nullptr",  # fdsva_so does not support external forces
            gravity_name = "gravity",
        ),
    )
    self.gen_inverse_dynamics_gradient_inner_function_call(
        updated_var_names = dict(
            s_dc_du_name = "s_fd_dc_du",
            s_vaf_name = "s_fd_vaf",
            s_q_name = "s_q",
            s_qd_name = "s_qd",
            s_temp_name = "s_fd_temp",
            d_temp_spill_name = spill_ptr_expr if use_spill else "nullptr",
            temp_spill_flag_name = "true" if use_spill else "false",
            gravity_name = "gravity",
        ),
    )
    self.gen_add_parallel_loop("ind", str(2*n*n))
    self.gen_add_code_line(f"int row = ind % {n}; int dc_col_offset = ind - row;")
    self.gen_add_code_line("T val = static_cast<T>(0);")
    self.gen_add_code_line(f"for(int col = 0; col < {n}; col++) {{", True)
    self.gen_add_code_line(f"int index = (row <= col) * (col * {n} + row) + (row > col) * (row * {n} + col);")
    self.gen_add_code_line("val += s_Minv[index] * s_fd_dc_du[dc_col_offset + col];")
    self.gen_add_end_control_flow()
    self.gen_add_code_line("s_df_du[ind] = -val;")
    self.gen_add_end_control_flow()
    self.gen_add_sync()
    
def gen_fdsva_so_contract_function_call(self, updated_var_names = None,
                                     scratch_in_smem_expr = "true"):
    var_names = dict( \
        s_df2_name = "s_df2", \
        s_idsva_so_name = "s_idsva_so", \
        s_Minv_name = "s_Minv", \
        s_df_du_name = "s_df_du", \
        s_temp_name = "s_temp", \
        d_workspace_name = "nullptr", \
        gravity_name = "gravity"
    )
    if updated_var_names is not None:
        for key,value in updated_var_names.items():
            var_names[key] = value
    fdsva_so_code_start = "fdsva_so_contract<T, " + scratch_in_smem_expr + ">(" + var_names["s_df2_name"] + ", " + var_names["s_idsva_so_name"] + ", " + var_names["s_Minv_name"] + ", " + var_names["s_df_du_name"] + ", "
    fdsva_so_code_end = var_names["s_temp_name"] + ", " + var_names["d_workspace_name"] + ", " + var_names["gravity_name"] + ");"
    fdsva_so_code_middle = self.gen_insert_helpers_function_call()
    fdsva_so_code = fdsva_so_code_start + fdsva_so_code_middle + fdsva_so_code_end
    self.gen_add_code_line(fdsva_so_code)

def gen_fdsva_so_device_function_call(self,
                                          scratch_in_smem_expr = "true",
                                          fd_grad_use_spill_expr = "false",
                                          contract_in_smem_expr = "true",
                                          d_workspace_pool_name = "nullptr",
                                          d_fd_grad_spill_name = "nullptr",
                                          s_fdsva_temp_name = "nullptr"):
    """Emit the call to `fdsva_so_device`. Arg order MUST match the def in
    gen_fdsva_so_device. Pool/spill regions default to nullptr (unused under
    the matching if-constexpr); the kernel passes real pointers per tier."""
    tmpl = "<T, " + scratch_in_smem_expr + ", " + fd_grad_use_spill_expr + ", " + contract_in_smem_expr + ">"
    start = "fdsva_so_device" + tmpl + "(s_df2, s_idsva_so, s_Minv, s_df_du, s_qdd, s_q, s_qd, s_u, "
    middle = self.gen_insert_helpers_function_call()
    end = ("s_temp, " + d_workspace_pool_name + ", " + d_fd_grad_spill_name + ", "
           + s_fdsva_temp_name + ", d_robotModel, gravity);")
    self.gen_add_code_line(start + middle + end)

def gen_fdsva_so_device(self):
    """Emit `fdsva_so_device` — the whole fdsva_so orchestration as ONE
    inner that OWNS its scratch (s_temp) placement (inner-owns-placement; see
    docs/idsva_so_inner_refactor_notes.md). It wraps, in order:
      load_update_XImats -> direct_minv_inner -> forward_dynamics_inner ->
      fd-gradient-inline -> idsva_so_{world,body}_inner -> fdsva_so_contract.
    Because the s_temp repoint happens at the very top, EVERY consumer below —
    including the XImats helper's sincos scratch — follows the placement, so the
    kernel never repoints s_temp from the outside.

    Template flags:
      SCRATCH_IN_SMEM   : the shared s_temp pool lives in smem (true) or routes
                          to d_workspace (false). This is the dominant lever — on
                          big floating humanoids the fd-gradient pool (~180 KB)
                          is what overflows, so false makes it fit.
      FD_GRAD_USE_SPILL : the fd-gradient inline spills its da_dq..fxvi band to
                          d_fd_grad_spill (only meaningful when the pool is smem).
      CONTRACT_IN_SMEM  : the fdsva_so_contract 4*NV^3 contraction scratch placement.

    Pointer params are caller-supplied (the kernel decides where the OUTPUTS and
    df_du/Minv live — smem, device arrays, or workspace bands — and hands in the
    spill regions): only the s_temp POOL placement is the inner's call."""
    n = self.robot.get_num_vel()
    func_params = [
        "s_df2/s_idsva_so/s_Minv/s_df_du/s_qdd/s_q/s_qd/s_u are the in/out buffers (caller places)",
        "s_temp is the shared scratch pool (used when SCRATCH_IN_SMEM)",
        "d_workspace is the global scratch pool (used when !SCRATCH_IN_SMEM)",
        "d_fd_grad_spill is the fd-gradient da_df spill band (used when FD_GRAD_USE_SPILL)",
        "s_fdsva_temp is the 4*NV^3 contraction scratch (smem or workspace per CONTRACT_IN_SMEM)",
        "d_robotModel holds XImats/topology; gravity is the gravity constant",
    ]
    func_def_start = "void fdsva_so_device("
    func_def_middle = ("T *s_df2, T *s_idsva_so, T *s_Minv, T *s_df_du, T *s_qdd, "
                       "const T *s_q, const T *s_qd, const T *s_u, ")
    func_def_end = ("T *s_temp, T *d_workspace, T *d_fd_grad_spill, T *s_fdsva_temp, "
                    "const robotModel<T> *d_robotModel, const T gravity) {")
    func_def_middle, func_params = self.gen_insert_helpers_func_def_params(func_def_middle, func_params, -2)
    func_def = func_def_start + func_def_middle + func_def_end
    self.gen_add_func_doc("fdsva_so orchestration as a single inner-owns-placement device function",
                          [], func_params, None)
    self.gen_add_code_line("template <typename T, bool SCRATCH_IN_SMEM = true, bool FD_GRAD_USE_SPILL = false, bool CONTRACT_IN_SMEM = true>")
    # __forceinline__ so the whole orchestration inlines into the calling kernel.
    # Under -rdc (single-call/anti-LICM build) a separate __device__ wrapper keeps
    # its callees (e.g. direct_minv_inner, ~108 regs) as distinct functions whose
    # regcount must fit the kernel's launch_bounds budget (80 at LITE / 64 at
    # MINIMAL) -> ptxas regcount error. Inlining folds them into the kernel (as the
    # pre-refactor inline orchestration did). See HANDOFF.md / "Problem 1".
    self.gen_add_code_line("__device__ __forceinline__")
    self.gen_add_code_line(func_def, True)
    # Inner owns the pool placement; the repoint covers every consumer below
    # (incl. the XImats helper's sincos scratch), so no caller-side repoint.
    self.gen_add_code_line("if constexpr (!SCRATCH_IN_SMEM) { s_temp = d_workspace; } else { (void)d_workspace; }")
    self.gen_add_code_line("T *d_temp_spill = nullptr; (void)d_temp_spill;  // idsva uses the (placed) s_temp pool directly")
    self.gen_load_update_XImats_helpers_function_call()
    self.gen_direct_minv_inner_function_call(f_in_smem_expr = "true")
    self.gen_add_code_line("forward_dynamics_inner<T, true>(s_qdd, s_q, s_qd, s_u, " + self.gen_insert_helpers_function_call() + "s_temp, nullptr, nullptr, gravity);")
    self.gen_add_sync()
    # fd-gradient inline; the band-spill variant is a compile-time choice.
    self.gen_add_code_line("if constexpr (FD_GRAD_USE_SPILL) {", True)
    self.gen_fdsva_so_fd_gradient_inline(use_spill=True, spill_ptr_expr="d_fd_grad_spill")
    self.gen_add_end_control_flow()
    self.gen_add_code_line("else {", True)
    self.gen_add_code_line("(void)d_fd_grad_spill;")
    self.gen_fdsva_so_fd_gradient_inline(use_spill=False, spill_ptr_expr="nullptr")
    self.gen_add_end_control_flow()
    if self.robot.floating_base:
        self.gen_idsva_so_world_frame_inner_function_call()
    else:
        self.gen_idsva_so_body_frame_inner_function_call()
        self.gen_idsva_so_body_frame_public_dvdq_layout_repair()
    self.gen_fdsva_so_contract_function_call(
        updated_var_names = dict(d_workspace_name = "s_fdsva_temp"),
        scratch_in_smem_expr = "CONTRACT_IN_SMEM")
    self.gen_add_end_function()

_FDSVA_SO_PICK_FLAGS = [
    # (use_global_tensors, use_workspace_temp, fd_grad_use_spill, use_workspace_df_du, use_workspace_Minv, use_workspace_idsva_temp)
    (False, False, False, False, False, False),   # pick 0: full smem
    (True,  False, False, False, False, False),   # pick 1: outputs to global
    (True,  True,  False, False, False, False),   # pick 2: + inner temp to global
    (True,  True,  True,  False, False, False),   # pick 3: + fd_grad da_df band to global
    (True,  True,  True,  True,  False, False),   # pick 4 (Phase 3e): + s_df_du to global
    (True,  True,  True,  True,  True,  False),   # pick 5 (Phase 3e): + s_Minv to global
    (True,  True,  False, False, False, True),    # pick 6: pool->global (whole s_temp via full inner SCRATCH_IN_SMEM=false); df_du/Minv stay in smem (small). Works for BOTH bases because fdsva_so_device hands the placed pool to whichever idsva inner it composes (world for floating, body for fixed) and the inner does the repoint via its own SCRATCH_IN_SMEM=false.
]

def _emit_fdsva_so_kernel_body_for_flags(self, n, NUM_POS, use_global_tensors, use_workspace_temp,
                                         fd_grad_use_spill, use_workspace_df_du, use_workspace_Minv,
                                         single_call_timing, use_workspace_idsva_temp=False):
    """Emit fdsva_so kernel body for one tier's spill flags.
    use_workspace_idsva_temp: route the embedded idsva_so inner's scratch to
    d_workspace (via the inner's SCRATCH_IN_SMEM=false). The idsva inner is the
    dominant s_temp consumer on big robots; spilling it drops the smem arena to
    the next-largest sub-inner (fd_grad / minv / fd). Works for BOTH bases —
    both the world inner (floating) and the body inner (fixed) now own their
    placement and accept SCRATCH_IN_SMEM=false, so the whole-arena pool->global
    fallback composes uniformly. (Earlier revisions of this comment said
    floating-only because the body inner had not yet migrated; that has since
    landed — see docs/idsva_so_inner_refactor_notes.md.)"""
    inner_idsva_so_temp_size = (
        self.gen_idsva_so_world_frame_temp_mem_size() if self.robot.floating_base
        else self.gen_idsva_so_body_frame_inner_temp_mem_size()
    )
    fd_grad_temp_size = (
        self.gen_fdsva_so_fd_gradient_inline_temp_mem_size_spilled() if fd_grad_use_spill
        else self.gen_fdsva_so_fd_gradient_inline_temp_mem_size()
    )
    if use_workspace_idsva_temp:
        # Pool -> global: fdsva_so_device runs with SCRATCH_IN_SMEM=false, so the
        # WHOLE shared s_temp pool (helper sincos + minv + fd + fd_grad + idsva) lives
        # in d_workspace. The smem s_temp slot is unused -> size 0.
        shared_temp_size = 0
    else:
        shared_temp_size = max(inner_idsva_so_temp_size, fd_grad_temp_size)
        if not use_workspace_temp:
            shared_temp_size = max(shared_temp_size, self.gen_fdsva_so_contract_temp_mem_size())
    # Phase 3e: s_df_du and s_Minv can now be in workspace too. Drop them from
    # extra_t_buffers when spilled; declare workspace pointers in the body.
    extra_t_buffers = [("s_q_qd_u", NUM_POS + 2*n), ("s_qdd", n)]
    if not use_workspace_Minv:
        extra_t_buffers.append(("s_Minv", n*n))
    if not use_workspace_df_du:
        extra_t_buffers.append(("s_df_du", 2*n*n))
    if not use_global_tensors:
        extra_t_buffers.append(("s_idsva_so", n*n*n*4))
        extra_t_buffers.append(("s_df2", 4*n*n*n))
    self.gen_XImats_helpers_temp_shared_memory_code(shared_temp_size, extra_t_buffers = extra_t_buffers)
    needs_d_workspace = use_workspace_temp or fd_grad_use_spill or use_workspace_df_du or use_workspace_Minv
    if not needs_d_workspace:
        self.gen_add_code_line("(void)d_workspace;")
    if not use_global_tensors:
        self.gen_add_code_line("(void)d_idsva_so;")
    self.gen_add_code_line("T *s_q = s_q_qd_u; T *s_qd = &s_q_qd_u[" + str(NUM_POS) + "]; T *s_u = &s_q_qd_u[" + str(NUM_POS + n) + "];")
    self.gen_add_code_line("T *d_temp_spill = nullptr; (void)d_temp_spill;")
    # Inner-controlled placement: forward_dynamics_inner slices its own Minv-F
    # from s_temp (smem; this kernel does not surgically spill Minv/FD-F — its
    # tiers spill the SO outputs / df_du / Minv instead).
    # Canonical: build the shared helper ARGS via gen_insert_helpers_function_call
    # (was a bespoke "make a def-params string then .replace() the types out" hack,
    # which silently dropped/duplicated args when the helper signature changed).
    fd_start = "forward_dynamics_inner<T, true>(s_qdd, s_q, s_qd, s_u, " + self.gen_insert_helpers_function_call()
    fd_end = "s_temp, nullptr, nullptr, gravity);"  # trailing nullptr = no external forces
    if not single_call_timing:
        self.gen_add_parallel_loop("k","NUM_TIMESTEPS",block_level = True)
        self.gen_kernel_load_inputs("q_qd_u",str(NUM_POS + 2*n),stride="stride_q_qd_u")
        if use_global_tensors:
            self.gen_add_code_line(f'T *s_df2 = &d_df2[k*{4*n**3}];')
            self.gen_add_code_line(f'T *s_idsva_so = &d_idsva_so[k*{4*n**3}];')
        if use_workspace_temp:
            self.gen_add_code_line('T *s_fdsva_temp = reinterpret_cast<T *>(&d_workspace[k*GRID_WORKSPACE_BYTES_PER_TIMESTEP<T>() + GRID_SO_WORKSPACE_TEMP_OFFSET_BYTES<T>()]);')
        if fd_grad_use_spill:
            self.gen_add_code_line('T *d_fd_grad_spill = reinterpret_cast<T *>(&d_workspace[k*GRID_WORKSPACE_BYTES_PER_TIMESTEP<T>()]);')
        if use_workspace_df_du:
            # Phase 3e: s_df_du in L2-pinned workspace, in its own dedicated section
            # past grad + SO (avoids conflict with fd_grad_spill which is at offset 0).
            self.gen_add_code_line('T *s_df_du = reinterpret_cast<T *>(&d_workspace[k*GRID_WORKSPACE_BYTES_PER_TIMESTEP<T>() + GRID_FDSVA_SO_SPILL_OFFSET_BYTES<T>()]);')
        if use_workspace_Minv:
            # Phase 3e: s_Minv lives just past s_df_du in the FDSVA_SO spill section.
            self.gen_add_code_line('T *s_Minv = reinterpret_cast<T *>(&d_workspace[k*GRID_WORKSPACE_BYTES_PER_TIMESTEP<T>() + GRID_FDSVA_SO_SPILL_OFFSET_BYTES<T>() + ' + str(2*n*n) + '*sizeof(T)]);')
        self.gen_add_code_line("// compute — the orchestration inner owns its s_temp pool placement")
        # Pool->global reuses the (non-concurrent) fdsva SO-temp region; the
        # contraction uses the same region in its later phase. See full inner.
        self.gen_fdsva_so_device_function_call(
            scratch_in_smem_expr = "false" if use_workspace_idsva_temp else "true",
            fd_grad_use_spill_expr = "true" if fd_grad_use_spill else "false",
            contract_in_smem_expr = "false" if use_workspace_temp else "true",
            d_workspace_pool_name = "s_fdsva_temp" if use_workspace_idsva_temp else "nullptr",
            d_fd_grad_spill_name = "d_fd_grad_spill" if fd_grad_use_spill else "nullptr",
            s_fdsva_temp_name = "s_fdsva_temp" if use_workspace_temp else "nullptr")
        self.gen_add_sync()
        if not use_global_tensors: self.gen_kernel_save_result("df2",str(4*n*n*n),stride=f"{4*n**3}")
        self.gen_add_end_control_flow()
    else:
        self.gen_kernel_load_inputs("q_qd_u",str(NUM_POS + 2*n))
        self.gen_add_code_line("// compute with NUM_TIMESTEPS as NUM_REPS for timing")
        self.gen_add_code_line("for (int rep = 0; rep < NUM_TIMESTEPS; rep++){", True)
        self.gen_anti_licm_input_reload("q_qd_u",str(NUM_POS + 2*n))
        if use_global_tensors:
            self.gen_add_code_line('T *s_df2 = d_df2;')
            self.gen_add_code_line('T *s_idsva_so = d_idsva_so;')
        if use_workspace_temp:
            self.gen_add_code_line('T *s_fdsva_temp = reinterpret_cast<T *>(&d_workspace[GRID_SO_WORKSPACE_TEMP_OFFSET_BYTES<T>()]);')
        if fd_grad_use_spill:
            self.gen_add_code_line('T *d_fd_grad_spill = reinterpret_cast<T *>(d_workspace);')
        if use_workspace_df_du:
            self.gen_add_code_line('T *s_df_du = reinterpret_cast<T *>(&d_workspace[GRID_FDSVA_SO_SPILL_OFFSET_BYTES<T>()]);')
        if use_workspace_Minv:
            self.gen_add_code_line('T *s_Minv = reinterpret_cast<T *>(&d_workspace[GRID_FDSVA_SO_SPILL_OFFSET_BYTES<T>() + ' + str(2*n*n) + '*sizeof(T)]);')
        self.gen_fdsva_so_device_function_call(
            scratch_in_smem_expr = "false" if use_workspace_idsva_temp else "true",
            fd_grad_use_spill_expr = "true" if fd_grad_use_spill else "false",
            contract_in_smem_expr = "false" if use_workspace_temp else "true",
            d_workspace_pool_name = "s_fdsva_temp" if use_workspace_idsva_temp else "nullptr",
            d_fd_grad_spill_name = "d_fd_grad_spill" if fd_grad_use_spill else "nullptr",
            s_fdsva_temp_name = "s_fdsva_temp" if use_workspace_temp else "nullptr")
        self.gen_add_end_control_flow()
        if not use_global_tensors: self.gen_kernel_save_result("df2",str(4*n*n*n))


def gen_fdsva_so_kernel(self, single_call_timing = False):
    # NUM_VEL is the SO tensor dimension (rank-3 nv*nv*nv); NUM_POS is q-vector size.
    n = self.robot.get_num_vel()
    NUM_POS = self.robot.get_num_pos()
    func_params = ["d_df2 is the second derivatives of forward dynamics WRT q,qd,tau", \
                    "d_q_qd_u is the vector of joint positions, velocities, torques", \
                    "stride_q_qd_u is the stride between each q, qd, qdd", \
                    "d_workspace is the generated global spill workspace", \
                    "d_idsva_so is the pointer to the idsva_so output tensor in global memory", \
                    "d_robotModel is the pointer to the initialized model specific helpers on the GPU (XImats, topology_helpers, etc.)", \
                    "gravity is the gravity constant", \
                    "num_timesteps is the length of the trajectory points we need to compute over (or overloaded as test_iters for timing)"]
    func_notes = []
    func_def_start = "void fdsva_so_kernel(T *d_df2, const T *d_q_qd_u, const int stride_q_qd_u, unsigned char *d_workspace, T *d_idsva_so, "
    func_def_end = "const robotModel<T> *d_robotModel, const T gravity, const int NUM_TIMESTEPS) {"
    func_def = func_def_start + func_def_end
    if single_call_timing:
        func_def = func_def.replace("kernel(", "kernel_single_timing(")
    self.gen_add_func_doc("Compute the FDSVA_SO (Second Order of Forward Dynamics with Spacial Vector Algebra)", func_notes, func_params, None)
    self.gen_add_code_line("template <typename T, int RESOURCE_TIER = GRID_DEFAULT_RESOURCE_TIER>")
    self.gen_add_code_line("__global__")
    self.gen_add_code_line("__launch_bounds__(tier_max_threads<RESOURCE_TIER>())")
    self.gen_add_code_line(func_def, True)
    picks = getattr(self, "fdsva_so_spill_tier_3way", (5, 5, 5))
    if picks[0] == picks[1] == picks[2]:
        ugt, uwt, fgs, uwdfdu, uwminv, uwit = _FDSVA_SO_PICK_FLAGS[picks[0]]
        _emit_fdsva_so_kernel_body_for_flags(self, n, NUM_POS, ugt, uwt, fgs, uwdfdu, uwminv, single_call_timing, uwit)
    else:
        tier_names = ("TIER_PERF", "TIER_LITE", "TIER_MINIMAL")
        for tier_idx, (tier_name, pick) in enumerate(zip(tier_names, picks)):
            ugt, uwt, fgs, uwdfdu, uwminv, uwit = _FDSVA_SO_PICK_FLAGS[pick]
            head = "if constexpr (RESOURCE_TIER == " + tier_name + ") {" if tier_idx == 0 else \
                   "else if constexpr (RESOURCE_TIER == " + tier_name + ") {"
            self.gen_add_code_line(head, True)
            _emit_fdsva_so_kernel_body_for_flags(self, n, NUM_POS, ugt, uwt, fgs, uwdfdu, uwminv, single_call_timing, uwit)
            self.gen_add_end_control_flow()
    self.gen_add_end_function()

def gen_fdsva_so_host(self, mode = 0):
    # NUM_VEL is the SO tensor dimension. See gen_fdsva_so_kernel for details.
    n = self.robot.get_num_vel()
    # default is to do the full kernel call -- options are for single timing or compute only kernel wrapper
    single_call_timing = True if mode == 1 else False
    compute_only = True if mode == 2 else False

    # define function def and params
    func_params = ["hd_data is the packaged input and output pointers", \
                   "d_robotModel is the pointer to the initialized model specific helpers on the GPU (XImats, topology_helpers, etc.)", \
                   "gravity is the gravity constant,", \
                   "num_timesteps is the length of the trajectory points we need to compute over (or overloaded as test_iters for timing)", \
                   "streams are pointers to CUDA streams for async memory transfers (if needed)"]
    func_notes = []
    func_def_start = "void fdsva_so(gridData<T, KIND> *hd_data, const robotModel<T> *d_robotModel, const T gravity, const int num_timesteps,"
    func_def_end =   "                      const dim3 block_dimms, const dim3 thread_dimms, cudaStream_t *streams) {"
    if single_call_timing:
        func_def_start = func_def_start.replace("(", "_single_timing(")
        func_def_end = "              " + func_def_end
    if compute_only:
        func_def_start = func_def_start.replace("(", "_compute_only(")
        func_def_end = "             " + func_def_end.replace(", cudaStream_t *streams", "")
    # then generate the code
    self.gen_add_func_doc("Compute the FDSVA_SO (Second Order of Forward Dynamics with Spacial Vector Algebra)",\
                          func_notes,func_params,None)
    self.gen_add_code_line("template <typename T, gridDataKind KIND = GRID_DATA_ALL>")
    self.gen_add_code_line("__host__")
    self.gen_add_code_line(func_def_start)
    self.gen_add_code_line(func_def_end, True)
    self.gen_add_code_line("static_assert(KIND == GRID_DATA_ALL || KIND == GRID_DATA_DYNAMICS, \"fdsva_so requires all-data or dynamics gridData\");")

    func_call_start = "fdsva_so_kernel<T><<<block_dimms,thread_dimms,FDSVA_SO_DYNAMIC_SHARED_MEM_BYTES<T>()>>>(hd_data->d_df2,hd_data->d_q_qd_u,stride_q_qd_qdd,hd_data->d_workspace,hd_data->d_idsva_so,"
    func_call_end = "d_robotModel,gravity,num_timesteps);"
    self.gen_add_code_line("int stride_q_qd_qdd = Q_QD_U_STRIDE;")
    if single_call_timing:
        func_call_start = func_call_start.replace("kernel<T>","kernel_single_timing<T>")
    if not compute_only:
        # start code with memory transfer
        self.gen_add_code_lines(["// start code with memory transfer", \
                                 "gpuErrchk(cudaMemcpyAsync(hd_data->d_q_qd_u,hd_data->h_q_qd_u,stride_q_qd_qdd" + \
                                    ("*num_timesteps" if not single_call_timing else "") + "*sizeof(T),cudaMemcpyHostToDevice,streams[0]));", \
                                 "gpuErrchkKernel();"])    
    
    # then compute:
    self.gen_add_code_line("// call the kernel")
    func_call = func_call_start + func_call_end
    func_call_code = [func_call, "gpuErrchkKernel();"]
    # wrap function call in timing (if needed)
    if single_call_timing:
        func_call_code.insert(0,"struct timespec start, end; clock_gettime(CLOCK_MONOTONIC,&start);")
        func_call_code.append("clock_gettime(CLOCK_MONOTONIC,&end);")
    self.gen_add_code_line("gpuErrchk(grid_check_dynamic_shared_memory_bytes(\"fdsva_so\", FDSVA_SO_DYNAMIC_SHARED_MEM_BYTES<T>()));")
    workspace_bytes = "GRID_WORKSPACE_BYTES_PER_TIMESTEP<T>()*GRID_WORKSPACE_SLOTS*" + ("1" if single_call_timing else "num_timesteps")
    self.gen_add_code_line("if (GRID_FDSVA_SO_USES_WORKSPACE_ANY_TIER) {gpuErrchk(grid_begin_l2_persisting(0, hd_data->d_workspace, " + workspace_bytes + "));}")
    self.gen_add_code_lines(func_call_code)
    self.gen_add_code_line("if (GRID_FDSVA_SO_USES_WORKSPACE_ANY_TIER) {gpuErrchk(grid_end_l2_persisting(0));}")
    if not compute_only:
        # then transfer memory back
        self.gen_add_code_lines(["// finally transfer the result back", \
                                "gpuErrchk(cudaMemcpy(hd_data->h_df2,hd_data->d_df2," + \
                                ("num_timesteps*" if not single_call_timing else "") + str(4*n**3) + "*sizeof(T),cudaMemcpyDeviceToHost));",
                                "gpuErrchkKernel();"])
    # finally report out timing if requested
    if single_call_timing:
        from ..algo_registry import single_call_printf_line
        self.gen_add_code_line(single_call_printf_line("fdsva_so"))
    self.gen_add_end_function()

def gen_fdsva_so(self):
    # first the contraction sub-inner used by the orchestration _device
    self.gen_fdsva_so_contract()
    # then the canonical _device (orchestrator: owns s_temp placement; called from kernel)
    self.gen_fdsva_so_device()
    # then the kernels (call _device with caller-allocated smem)
    self.gen_fdsva_so_kernel(True)
    self.gen_fdsva_so_kernel(False)
    # then the host wrappers
    self.gen_fdsva_so_host(0)
    self.gen_fdsva_so_host(1)
    self.gen_fdsva_so_host(2)
    
