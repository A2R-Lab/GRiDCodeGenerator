MEMORY_THRESHOLD = 8 # Max num joints for shared mem allocation of result


def gen_fdsva_so_inner(self, use_thread_group = False):
	# construct the boilerplate and function definition
    n = self.robot.get_num_vel()
    inner_arena_size = 4 * n**3
    func_params = ["s_df2 are the second derivatives of forward dynamics WRT q,qd,tau", \
                "s_idsva_so are the second derivative tensors of inverse dynamics", \
                "s_Minv is the inverse mass matrix", \
                "s_df_du is the gradient of the forward dynamics", \
                "s_temp is the (shared) scratch buffer; size FDSVA_SO_INNER_SMEM_BYTES<T, SCRATCH_IN_SMEM>() bytes (= " + str(inner_arena_size) + "*sizeof(T) when SCRATCH_IN_SMEM, else 0)", \
                "s_workspace is the global scratch buffer; size FDSVA_SO_INNER_WORKSPACE_BYTES<T, SCRATCH_IN_SMEM>() bytes (= " + str(inner_arena_size) + "*sizeof(T) when !SCRATCH_IN_SMEM, else 0). Pass nullptr when SCRATCH_IN_SMEM", \
                "gravity is the gravity constant"]
    func_def_start = "void fdsva_so_inner("
    func_def_middle = "T *s_df2, T *s_idsva_so, T *s_Minv, T *s_df_du, "
    func_def_end = "T *s_temp, T *s_workspace, const T gravity) {"
    func_notes = ["Assumes works with IDSVA",
                  "Inline-CUDA users: SCRATCH_IN_SMEM selects where the 4*NV^3 scratch arena lives.",
                  "  true  -> s_temp (shared memory; fastest, current default).",
                  "  false -> s_workspace (global; frees shared memory for the caller's outer kernel).",
                  "The placement is the INNER's choice (made at the top of this function) so the",
                  "kernel/device caller just sizes both arenas from the exposed *_BYTES constants",
                  "and hands both pointers in. Codegen maps each RESOURCE_TIER to a SCRATCH_IN_SMEM",
                  "value per robot (small robots keep all tiers in smem; large robots spill)."]
    if use_thread_group:
        func_def_start = func_def_start.replace("(", "(cgrps::thread_group tgrp, ")
        func_params.insert(0,"tgrp is the handle to the thread_group running this function")
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
    self.gen_add_code_line('// choice, keyed on SCRATCH_IN_SMEM: s_temp (shared) when true, s_workspace')
    self.gen_add_code_line('// (global) when false. A surgical-spill change is local here: repoint a')
    self.gen_add_code_line('// sub-buffer + update FDSVA_SO_INNER_{SMEM,WORKSPACE}_BYTES.')
    self.gen_add_code_line('T *inner_arena;')
    self.gen_add_code_line('if constexpr (SCRATCH_IN_SMEM) {', True)
    self.gen_add_code_line('(void)s_workspace;')
    self.gen_add_code_line('inner_arena = s_temp;')
    self.gen_add_end_control_flow()
    self.gen_add_code_line('else {', True)
    self.gen_add_code_line('(void)s_temp;')
    self.gen_add_code_line('inner_arena = s_workspace;')
    self.gen_add_end_control_flow()
    self.gen_add_code_line(f'T *inner_dq = inner_arena; // Inner term for d2a_dqdq')
    self.gen_add_code_line(f'T *inner_cross = inner_dq + {n**3}; // Inner term for d2a_dvdq (dM_dq*Minv)')
    self.gen_add_code_line(f'T *inner_tau = inner_cross + {n**3}; // Inner term for d2a_dtdq (d2tau_dvdq + dM_dq*da_dv)')
    self.gen_add_code_line(f'T *rot_dq = inner_tau + {n**3}; // Rotated (dM_dq*da_dq)^R term used to compute inner_dq')
    self.gen_add_code_line(f'\n\n')

    # Start inner term for d2a_dqdq
    self.gen_add_code_line('// Start inner term for d2a_dqdq & Fill out Minv')
    self.gen_add_parallel_loop("ind",str(n**3 + n*n),use_thread_group)
    self.gen_add_code_line(f'int i = ind / {n*n} % {n}; int j = ind / {n} % {n}; int k = ind % {n};')
    self.gen_add_code_line(f'if (ind < {n**3}) {{', True)
    self.gen_add_code_line(f'inner_dq[ind] = dot_prod<T, {n}, {n}, 1>(&dM_dq[{n*n}*i + k], &s_df_dq[{n}*j]);')
    self.gen_add_code_line(f'rot_dq[i*{n*n} + k*{n} + j] = inner_dq[ind];')
    self.gen_add_end_control_flow()
    self.gen_add_code_line(f'else if (k > j) s_Minv[j*{n} + k] = s_Minv[k*{n} + j];')
    self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)
    
    # 3Dx2D Tensor Computation defined as iLk,Lj->ijk
    self.gen_add_code_line('// Compute relevant inner subterms in parallel')
    self.gen_add_parallel_loop("ind",str(3*n**3),use_thread_group)
    self.gen_add_code_line(f'int i = ind / {n*n} % {n}; int j = ind / {n} % {n}; int k = ind % {n};')
    self.gen_add_code_line(f'if (ind < {n**3}) inner_dq[ind] += rot_dq[ind] + d2tau_dqdq[ind]; // Started with dM_dq*da_dq')
    self.gen_add_code_line(f'else if (ind < {2*n**3}) inner_cross[i*{n*n} + k*{n} + j] = dot_prod<T, {n}, {n}, 1>(&dM_dq[{n*n}*i + k], &s_df_dqd[{n}*j]) + d2tau_dvdq[i*{n*n} + j*{n} + k];')
    self.gen_add_code_line(f'else inner_tau[i*{n*n} + k*{n} + j] = dot_prod<T, {n}, {n}, 1>(&dM_dq[{n*n}*i + k], &s_Minv[{n}*j]);')
    self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)

    # 2Dx3D tensor computation defined as iL,Ljk->ijk
    self.gen_add_code_line('// Multiply by -Minv to finish algorithm')
    # PERF EXPERIMENT CANDIDATE (2026-05-18): this 4*n^3 parallel_loop has the
    # highest FMA count in fdsva_so_inner (g1_floating: ~6M FMAs out of ~9M
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
    self.gen_add_parallel_loop("ind",str(4*n**3),use_thread_group)
    self.gen_add_code_line(f'int i = ind / {n*n} % {n}; int j = ind / {n} % {n}; int k = ind % {n};')
    self.gen_add_code_line(f'if (ind < {n**3}) d2a_dqdq[i*{n*n} + j*{n} + k] = -dot_prod<T, {n}, {n}, {n*n}>(&s_Minv[i], &inner_dq[j + k*{n}]);')
    self.gen_add_code_line(f'else if (ind < {2*n**3}) d2a_dvdq[i*{n*n} + j*{n} + k] = -dot_prod<T, {n}, {n}, {n*n}>(&s_Minv[i], &inner_cross[j + k*{n}]);')
    self.gen_add_code_line(f'else if (ind < {3*n**3}) d2a_dvdv[i*{n*n} + j*{n} + k] = -dot_prod<T, {n}, {n}, {n*n}>(&s_Minv[i], &d2tau_dvdv[j + k*{n}]);')
    self.gen_add_code_line(f'else d2a_dtdq[i*{n*n} + j*{n} + k] = -dot_prod<T, {n}, {n}, {n*n}>(&s_Minv[i], &inner_tau[j + k*{n}]);')
    self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)

    self.gen_add_end_function()

def gen_fdsva_so_device_temp_mem_size(self):
    # Same as idsva_so because idsva_so is called and takes more memory
    NV = self.robot.get_num_vel()
    jids_a, ancestors = self.robot.get_jid_ancestor_ids(include_joint=True)
    return int(36 * NV * 10 + 30 * NV + 6 + len(jids_a)*36)
    
def gen_fdsva_so_inner_temp_mem_size(self):
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

def gen_fdsva_so_fd_gradient_inline(self, use_thread_group = False, use_spill = False, spill_ptr_expr = "nullptr"):
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
        use_thread_group,
        compute_c = False,
        use_qdd_input = True,
        updated_var_names = dict(
            s_vaf_name = "s_fd_vaf",
            s_q_name = "s_q",
            s_qd_name = "s_qd",
            s_qdd_name = "s_qdd",
            s_temp_name = "s_fd_temp",
            gravity_name = "gravity",
        ),
    )
    self.gen_inverse_dynamics_gradient_inner_function_call(
        use_thread_group,
        updated_var_names = dict(
            s_dc_du_name = "s_fd_dc_du",
            s_vaf_name = "s_fd_vaf",
            s_q_name = "s_q",
            s_qd_name = "s_qd",
            s_temp_name = "s_fd_temp",
            s_temp_spill_name = spill_ptr_expr if use_spill else "nullptr",
            temp_spill_flag_name = "true" if use_spill else "false",
            gravity_name = "gravity",
        ),
    )
    self.gen_add_parallel_loop("ind", str(2*n*n), use_thread_group)
    self.gen_add_code_line(f"int row = ind % {n}; int dc_col_offset = ind - row;")
    self.gen_add_code_line("T val = static_cast<T>(0);")
    self.gen_add_code_line(f"for(int col = 0; col < {n}; col++) {{", True)
    self.gen_add_code_line(f"int index = (row <= col) * (col * {n} + row) + (row > col) * (row * {n} + col);")
    self.gen_add_code_line("val += s_Minv[index] * s_fd_dc_du[dc_col_offset + col];")
    self.gen_add_end_control_flow()
    self.gen_add_code_line("s_df_du[ind] = -val;")
    self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)
    
def gen_fdsva_so_inner_function_call(self, use_thread_group = False, updated_var_names = None,
                                     scratch_in_smem_expr = "true"):
    var_names = dict( \
        s_df2_name = "s_df2", \
        s_idsva_so_name = "s_idsva_so", \
        s_Minv_name = "s_Minv", \
        s_df_du_name = "s_df_du", \
        s_temp_name = "s_temp", \
        s_workspace_name = "nullptr", \
        gravity_name = "gravity"
    )
    if updated_var_names is not None:
        for key,value in updated_var_names.items():
            var_names[key] = value
    fdsva_so_code_start = "fdsva_so_inner<T, " + scratch_in_smem_expr + ">(" + var_names["s_df2_name"] + ", " + var_names["s_idsva_so_name"] + ", " + var_names["s_Minv_name"] + ", " + var_names["s_df_du_name"] + ", "
    fdsva_so_code_end = var_names["s_temp_name"] + ", " + var_names["s_workspace_name"] + ", " + var_names["gravity_name"] + ");"
    if use_thread_group:
        id_code_start = id_code_start.replace("(","(tgrp, ")
    fdsva_so_code_middle = self.gen_insert_helpers_function_call()
    fdsva_so_code = fdsva_so_code_start + fdsva_so_code_middle + fdsva_so_code_end
    self.gen_add_code_line(fdsva_so_code)

def gen_fdsva_so_device(self, use_thread_group = False):
    # NUM_VEL is the SO tensor dimension. See gen_fdsva_so_kernel for details.
    n = self.robot.get_num_vel()
    # construct the boilerplate and function definition
    func_params = ["s_df2 is the second derivatives of forward dynamics WRT q,qd,tau", \
                   "s_df_du is a pointer to memory for the derivative of forward dynamics WRT q,qd of size 2*NUM_JOINTS*NUM_JOINTS = " + str(2*n*n), \
                   "s_q is the vector of joint positions", \
                   "s_qd is the vector of joint velocities", \
                   "s_u is the vector of joint control inputs", \
                   "d_robotModel is the pointer to the initialized model specific helpers on the GPU (XImats, topology_helpers, etc.)", \
                   "gravity is the gravity constant"]
    func_notes = []
    func_def_start = "void fdsva_so_device("
    func_def_middle = "T *s_df2, T *s_df_du, const T *s_q, const T *s_qd, const T *s_u, "
    func_def_end = "const robotModel<T> *d_robotModel, const T gravity) {"
    if use_thread_group:
        func_def_start += "cgrps::thread_group tgrp, "
        func_params.insert(0,"tgrp is the handle to the thread_group running this function")
    func_def = func_def_start + func_def_middle + func_def_end

    # then generate the code
    self.gen_add_func_doc("Compute the FDSVA_SO (Second Order of Forward Dyamics with Spacial Vector Algebra)",\
                          func_notes,func_params,None)
    self.gen_add_code_line("template <typename T>")
    self.gen_add_code_line("__device__")
    self.gen_add_code_line(func_def, True)

    shared_mem_size = max(
        self.gen_fdsva_so_device_temp_mem_size(),
        self.gen_fdsva_so_inner_temp_mem_size(),
        self.gen_fdsva_so_fd_gradient_inline_temp_mem_size(),
    )
    self.gen_XImats_helpers_temp_shared_memory_code(shared_mem_size, extra_t_buffers = [("s_Minv", n*n), ("s_qdd", n), ("s_idsva_so", n*n*n*4)])
    # body_frame_inner takes a grav-shim spill on floating; world_frame_inner
    # doesn't need it. Floating-base now dispatches to world_frame, so
    # s_temp_spill is unused there. Fixed-base body_frame_inner also doesn't
    # dereference the spill, so nullptr is safe for both paths.
    self.gen_add_code_line("T *s_temp_spill = nullptr;")

    # then load/update XI and run the algo. Inner-controlled placement: Minv and
    # FD each slice their own F-region from s_temp (this device path keeps F in
    # smem — no surgical spill at this layer).
    self.gen_load_update_XImats_helpers_function_call(use_thread_group)
    self.gen_direct_minv_inner_function_call(use_thread_group, f_in_smem_expr = "true")
    self.gen_add_code_line(f"forward_dynamics_inner<T, true>(s_qdd, s_q, s_qd, s_u, s_XImats, s_temp, nullptr, gravity);")
    self.gen_add_sync(use_thread_group)
    self.gen_fdsva_so_fd_gradient_inline(use_thread_group)
    if self.robot.floating_base:
        self.gen_idsva_so_world_frame_inner_function_call(use_thread_group)
    else:
        self.gen_idsva_so_body_frame_inner_function_call(use_thread_group)
        self.gen_idsva_so_body_frame_public_dvdq_layout_repair(use_thread_group)
    self.gen_fdsva_so_inner_function_call(use_thread_group)
    self.gen_add_end_function()

_FDSVA_SO_PICK_FLAGS = [
    # (use_global_tensors, use_workspace_temp, fd_grad_use_spill, use_workspace_df_du, use_workspace_Minv)
    (False, False, False, False, False),   # pick 0: full smem
    (True,  False, False, False, False),   # pick 1: outputs to global
    (True,  True,  False, False, False),   # pick 2: + inner temp to global
    (True,  True,  True,  False, False),   # pick 3: + fd_grad da_df band to global
    (True,  True,  True,  True,  False),   # pick 4 (Phase 3e): + s_df_du to global
    (True,  True,  True,  True,  True),    # pick 5 (Phase 3e): + s_Minv to global
]

def _emit_fdsva_so_kernel_body_for_flags(self, n, NUM_POS, use_global_tensors, use_workspace_temp,
                                         fd_grad_use_spill, use_workspace_df_du, use_workspace_Minv,
                                         single_call_timing, use_thread_group):
    """Emit fdsva_so kernel body for one tier's spill flags."""
    inner_idsva_so_temp_size = (
        self.gen_idsva_so_world_frame_temp_mem_size() if self.robot.floating_base
        else self.gen_idsva_so_body_frame_inner_temp_mem_size()
    )
    fd_grad_temp_size = (
        self.gen_fdsva_so_fd_gradient_inline_temp_mem_size_spilled() if fd_grad_use_spill
        else self.gen_fdsva_so_fd_gradient_inline_temp_mem_size()
    )
    shared_temp_size = max(inner_idsva_so_temp_size, fd_grad_temp_size)
    if not use_workspace_temp:
        shared_temp_size = max(shared_temp_size, self.gen_fdsva_so_inner_temp_mem_size())
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
    self.gen_add_code_line("T *s_temp_spill = nullptr;")
    if use_thread_group:
        self.gen_add_code_line("cgrps::thread_group tgrp = TBD;")
    # Inner-controlled placement: forward_dynamics_inner slices its own Minv-F
    # from s_temp (smem; this kernel does not surgically spill Minv/FD-F — its
    # tiers spill the SO outputs / df_du / Minv instead).
    fd_start = "forward_dynamics_inner<T, true>(s_qdd, s_q, s_qd, s_u, "
    fd_end = "s_temp, nullptr, gravity);"
    fd_start, _ = self.gen_insert_helpers_func_def_params(fd_start, [], -2)
    if 'T *' in fd_start: fd_start = fd_start.replace("T *","")
    if 'int *' in fd_start: fd_start = fd_start.replace("int *","")
    if not single_call_timing:
        self.gen_add_parallel_loop("k","NUM_TIMESTEPS",use_thread_group,block_level = True)
        self.gen_kernel_load_inputs("q_qd_u","stride_q_qd_u",str(NUM_POS + 2*n),use_thread_group)
        if use_global_tensors:
            self.gen_add_code_line(f'T *s_df2 = &d_df2[k*{4*n**3}];')
            self.gen_add_code_line(f'T *s_idsva_so = &d_idsva_so[k*{4*n**3}];')
        if use_workspace_temp:
            self.gen_add_code_line('T *s_fdsva_temp = reinterpret_cast<T *>(&d_workspace[k*GRID_WORKSPACE_BYTES_PER_TIMESTEP<T>() + GRID_SO_WORKSPACE_TEMP_OFFSET_BYTES<T>()]);')
        if fd_grad_use_spill:
            self.gen_add_code_line('T *s_fd_grad_spill = reinterpret_cast<T *>(&d_workspace[k*GRID_WORKSPACE_BYTES_PER_TIMESTEP<T>()]);')
        if use_workspace_df_du:
            # Phase 3e: s_df_du in L2-pinned workspace, in its own dedicated section
            # past grad + SO (avoids conflict with fd_grad_spill which is at offset 0).
            self.gen_add_code_line('T *s_df_du = reinterpret_cast<T *>(&d_workspace[k*GRID_WORKSPACE_BYTES_PER_TIMESTEP<T>() + GRID_FDSVA_SO_SPILL_OFFSET_BYTES<T>()]);')
        if use_workspace_Minv:
            # Phase 3e: s_Minv lives just past s_df_du in the FDSVA_SO spill section.
            self.gen_add_code_line('T *s_Minv = reinterpret_cast<T *>(&d_workspace[k*GRID_WORKSPACE_BYTES_PER_TIMESTEP<T>() + GRID_FDSVA_SO_SPILL_OFFSET_BYTES<T>() + ' + str(2*n*n) + '*sizeof(T)]);')
        self.gen_add_code_line("// compute")
        self.gen_load_update_XImats_helpers_function_call(use_thread_group)
        # Phase 3a: Minv inner takes s_F + s_temp separately (pack F at offset 0 of s_temp).
        self.gen_direct_minv_inner_function_call(use_thread_group, f_in_smem_expr = "true")
        self.gen_add_code_line(fd_start + fd_end)
        self.gen_add_sync(use_thread_group)
        self.gen_fdsva_so_fd_gradient_inline(
            use_thread_group,
            use_spill=fd_grad_use_spill,
            spill_ptr_expr="s_fd_grad_spill" if fd_grad_use_spill else "nullptr",
        )
        if self.robot.floating_base:
            self.gen_idsva_so_world_frame_inner_function_call(use_thread_group)
        else:
            self.gen_idsva_so_body_frame_inner_function_call(use_thread_group)
            self.gen_idsva_so_body_frame_public_dvdq_layout_repair(use_thread_group)
        # Inner-controlled scratch placement: pass both arenas (s_temp = smem,
        # s_workspace = the L2-pinned spill slot when this tier spills) and let
        # fdsva_so_inner pick via SCRATCH_IN_SMEM. The kernel no longer aims a
        # single pointer — it just provides both, sized from the constants.
        fdsva_updates = dict(s_workspace_name = "s_fdsva_temp") if use_workspace_temp else None
        self.gen_fdsva_so_inner_function_call(use_thread_group, updated_var_names = fdsva_updates,
                                              scratch_in_smem_expr = "false" if use_workspace_temp else "true")
        self.gen_add_sync(use_thread_group)
        if not use_global_tensors: self.gen_kernel_save_result("df2",f"{4*n**3}",str(4*n*n*n),use_thread_group)
        self.gen_add_end_control_flow()
    else:
        self.gen_kernel_load_inputs_single_timing("q_qd_u",str(NUM_POS + 2*n),use_thread_group)
        self.gen_add_code_line("// compute with NUM_TIMESTEPS as NUM_REPS for timing")
        self.gen_add_code_line("for (int rep = 0; rep < NUM_TIMESTEPS; rep++){", True)
        self.gen_anti_licm_input_reload("q_qd_u",str(NUM_POS + 2*n),use_thread_group)
        if use_global_tensors:
            self.gen_add_code_line('T *s_df2 = d_df2;')
            self.gen_add_code_line('T *s_idsva_so = d_idsva_so;')
        if use_workspace_temp:
            self.gen_add_code_line('T *s_fdsva_temp = reinterpret_cast<T *>(&d_workspace[GRID_SO_WORKSPACE_TEMP_OFFSET_BYTES<T>()]);')
        if fd_grad_use_spill:
            self.gen_add_code_line('T *s_fd_grad_spill = reinterpret_cast<T *>(d_workspace);')
        if use_workspace_df_du:
            self.gen_add_code_line('T *s_df_du = reinterpret_cast<T *>(&d_workspace[GRID_FDSVA_SO_SPILL_OFFSET_BYTES<T>()]);')
        if use_workspace_Minv:
            self.gen_add_code_line('T *s_Minv = reinterpret_cast<T *>(&d_workspace[GRID_FDSVA_SO_SPILL_OFFSET_BYTES<T>() + ' + str(2*n*n) + '*sizeof(T)]);')
        self.gen_load_update_XImats_helpers_function_call(use_thread_group)
        # Phase 3a: Minv inner takes s_F + s_temp separately.
        self.gen_direct_minv_inner_function_call(use_thread_group, f_in_smem_expr = "true")
        self.gen_add_code_line(fd_start + fd_end)
        self.gen_add_sync(use_thread_group)
        self.gen_fdsva_so_fd_gradient_inline(
            use_thread_group,
            use_spill=fd_grad_use_spill,
            spill_ptr_expr="s_fd_grad_spill" if fd_grad_use_spill else "nullptr",
        )
        if self.robot.floating_base:
            self.gen_idsva_so_world_frame_inner_function_call(use_thread_group)
        else:
            self.gen_idsva_so_body_frame_inner_function_call(use_thread_group, updated_var_names = dict(s_mem_name = "s_temp"))
            self.gen_idsva_so_body_frame_public_dvdq_layout_repair(use_thread_group)
        fdsva_updates = dict(s_workspace_name = "s_fdsva_temp") if use_workspace_temp else None
        self.gen_fdsva_so_inner_function_call(use_thread_group, updated_var_names = fdsva_updates,
                                              scratch_in_smem_expr = "false" if use_workspace_temp else "true")
        self.gen_add_end_control_flow()
        if not use_global_tensors: self.gen_kernel_save_result_single_timing("df2",str(4*n*n*n),use_thread_group)


def gen_fdsva_so_kernel(self, use_thread_group = False, single_call_timing = False):
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
        ugt, uwt, fgs, uwdfdu, uwminv = _FDSVA_SO_PICK_FLAGS[picks[0]]
        _emit_fdsva_so_kernel_body_for_flags(self, n, NUM_POS, ugt, uwt, fgs, uwdfdu, uwminv, single_call_timing, use_thread_group)
    else:
        tier_names = ("TIER_PERF", "TIER_LITE", "TIER_MINIMAL")
        for tier_idx, (tier_name, pick) in enumerate(zip(tier_names, picks)):
            ugt, uwt, fgs, uwdfdu, uwminv = _FDSVA_SO_PICK_FLAGS[pick]
            head = "if constexpr (RESOURCE_TIER == " + tier_name + ") {" if tier_idx == 0 else \
                   "else if constexpr (RESOURCE_TIER == " + tier_name + ") {"
            self.gen_add_code_line(head, True)
            _emit_fdsva_so_kernel_body_for_flags(self, n, NUM_POS, ugt, uwt, fgs, uwdfdu, uwminv, single_call_timing, use_thread_group)
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
    self.gen_add_code_line("if (GRID_FDSVA_SO_USES_WORKSPACE_TEMP) {gpuErrchk(grid_begin_l2_persisting(0, hd_data->d_workspace, " + workspace_bytes + "));}")
    self.gen_add_code_lines(func_call_code)
    self.gen_add_code_line("if (GRID_FDSVA_SO_USES_WORKSPACE_TEMP) {gpuErrchk(grid_end_l2_persisting(0));}")
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

def gen_fdsva_so(self, use_thread_group = False):
    # first generate the inner helper
    self.gen_fdsva_so_inner(use_thread_group)
    # then generate the device wrapper
    self.gen_fdsva_so_device(use_thread_group)
    # then generate the kernels
    self.gen_fdsva_so_kernel(use_thread_group, True)
    self.gen_fdsva_so_kernel(use_thread_group, False)
    # then generate the host wrappers
    self.gen_fdsva_so_host(0)
    self.gen_fdsva_so_host(1)
    self.gen_fdsva_so_host(2)
    
