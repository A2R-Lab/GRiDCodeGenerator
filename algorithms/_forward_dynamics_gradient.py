def gen_forward_dynamics_gradient_inner_temp_mem_size(self, use_qdd_Minv_input = False):
    n = self.robot.get_num_vel()
    minv_temp = self.gen_minv_inner_temp_mem_size()
    id_du_temp = self.gen_inverse_dynamics_gradient_inner_temp_mem_size()
    return max(minv_temp,id_du_temp) if not use_qdd_Minv_input else id_du_temp

def gen_forward_dynamics_gradient_inner_python(self, use_qdd_Minv_input = False,
                                               s_df_du_name = "s_df_du",
                                               d_temp_spill_name = "nullptr",
                                               temp_spill_flag_name = "false",
                                               d_f_ext_name = "d_f_ext"):
    n = self.robot.get_num_vel()
    if not use_qdd_Minv_input:
        #
        # TODO: there is a slightly faster way as s_v does not change -- thus no recompute needed
        #       but that requires a custom function to be written
        #
        self.gen_add_code_line("//TODO: there is a slightly faster way as s_v does not change -- thus no recompute needed")
        # Inner-controlled placement: minv_inner slices its own F-region
        # from the tail of s_temp (FD_DU keeps Minv-F in smem; its surgical spill
        # is the id_du da_df band, handled separately). After Minv returns, the
        # c+vaf/ID code reuses these bytes (the steps run sequentially).
        self.gen_minv_inner_function_call(f_in_smem_expr = "true")
        # updated_var_names = dict(s_c_name = "s_temp", s_vaf_name = "&s_temp[" + str(n) + "]", s_temp_name = "&s_temp[" + str(19*n) + "]")
        updated_var_names = dict(s_c_name = "s_temp", s_temp_name = "&s_temp[" + str(n) + "]", d_f_ext_name = d_f_ext_name)
        self.gen_inverse_dynamics_inner_function_call(compute_c = True, use_qdd_input = False, updated_var_names = updated_var_names)
        self.gen_forward_dynamics_finish_function_call(updated_var_names)
        self.gen_add_sync()
        self.gen_inverse_dynamics_inner_function_call(compute_c = False, use_qdd_input = True, updated_var_names = dict(d_f_ext_name = d_f_ext_name))
    # else just compute vaf
    else:
        self.gen_inverse_dynamics_inner_function_call(compute_c = False, use_qdd_input = True, updated_var_names = dict(d_f_ext_name = d_f_ext_name))
    # then run the gradient code
    self.gen_inverse_dynamics_gradient_inner_function_call(
        dict(d_temp_spill_name = d_temp_spill_name, temp_spill_flag_name = temp_spill_flag_name)
    )

    if self.DEBUG_MODE:
        self.gen_add_sync()
        self.gen_add_serial_ops()
        self.gen_add_code_lines(["printf(\"Minv\\n\");", \
                                 "printMat<T," + str(n) + "," + str(n) + ">(s_Minv," + str(n) + ");", \
                                 "printf(\"qdd\\n\");", \
                                 "printMat<T,1," + str(n) + ">(s_qdd,1);", \
                                 "printf(\"v\\n\");", \
                                 "printMat<T,6," + str(n) + ">(s_vaf,6);", \
                                 "printf(\"a\\n\");", \
                                 "printMat<T,6," + str(n) + ">(&s_vaf[6*" + str(n) + "],6);", \
                                 "printf(\"f\\n\");", \
                                 "printMat<T,6," + str(n) + ">(&s_vaf[12*" + str(n) + "],6);", \
                                 "printf(\"dc/dq\\n\");", \
                                 "printMat<T," + str(n) + "," + str(n) + ">(&s_dc_du[0]," + str(n) + ");", \
                                 "printf(\"dc/dqd\\n\");", \
                                 "printMat<T," + str(n) + "," + str(n) + ">(&s_dc_du[" + str(n*n) + "]," + str(n) + ");"])
        self.gen_add_end_control_flow()
        self.gen_add_sync()

    # and finally finish with df/du = -Minv*dc/du
    self.gen_minv_apply(
        n, s_df_du_name + "[ind]", "s_dc_du[dc_col_offset + col]",
        loop_var = "ind", loop_max = str(n*2*n),
        pre_lines = ["int row = ind % " + str(n) + "; int dc_col_offset = ind - row;"],
        comment_in_loop = False, negate = True)

def gen_forward_dynamics_gradient_device_function_call(self,
                                                           use_qdd_Minv_input = False,
                                                           scratch_in_smem_expr = "true",
                                                           use_da_df_spill_expr = "false",
                                                           s_df_du_name = "s_df_du",
                                                           d_workspace_pool_name = "nullptr",
                                                           d_temp_spill_name = "nullptr",
                                                           d_f_ext_name = "d_f_ext"):
    """Emit the call to `forward_dynamics_gradient_device`. Arg order MUST
    match the def in gen_forward_dynamics_gradient_device. The caller decides
    where the OUTPUT s_df_du lives (smem buffer or the global d_df_du band) and
    hands in the pool/spill regions; these default to nullptr (unused under the
    matching if-constexpr). The _qdd C++ name variant additionally threads the
    caller-provided s_qdd / s_Minv inputs."""
    fname = "forward_dynamics_gradient_device_qdd" if use_qdd_Minv_input else "forward_dynamics_gradient_device"
    tmpl = "<T, " + scratch_in_smem_expr + ", " + use_da_df_spill_expr + ">"
    start = fname + tmpl + "(" + s_df_du_name + ", s_q, s_qd, "
    if use_qdd_Minv_input:
        start += "s_qdd, s_Minv, "
    else:
        start += "s_u, "
    start += "s_vaf, s_dc_du, s_qdd, s_Minv, "
    middle = self.gen_insert_helpers_function_call()
    end = ("s_temp, " + d_workspace_pool_name + ", " + d_temp_spill_name + ", "
           + "d_robotModel, " + d_f_ext_name + ", gravity);")
    self.gen_add_code_line(start + middle + end)

def gen_forward_dynamics_gradient_device(self, use_qdd_Minv_input = False):
    """Emit `forward_dynamics_gradient_device` — the whole fd_du orchestration
    as ONE inner that OWNS its scratch (s_temp) placement (inner-owns-placement;
    mirrors gen_inverse_dynamics_gradient_device / gen_fdsva_so_device). It
    wraps, in order:
      [repoint s_temp] -> load_update_XImats -> minv_inner (f_in_smem=true)
      -> inverse_dynamics_inner (c+vaf) -> forward_dynamics_finish -> id_inner (vaf)
      -> inverse_dynamics_gradient_inner (the id_du BAND sub-inner) -> df/du = -Minv*dc/du.
    Because the s_temp repoint happens at the very top, EVERY consumer below —
    including the XImats helper's sincos scratch and minv's own F-region —
    follows the placement, so the kernel never repoints s_temp from the outside.

    TWO independent template flags:
      SCRATCH_IN_SMEM  : the shared s_temp pool lives in smem (true) or routes the
                         WHOLE pool to d_workspace (false; the rung-2 global-temp
                         path).
      USE_DA_DF_SPILL  : the id_du band selectively spills its da_dq..fxvi band to
                         d_temp_spill (rung 1). Threaded through to the BAND
                         sub-inner's grid_id_du_temp_ptr<T, USE_DA_DF_SPILL> helper.
    The 3-rung menu (see _FD_DU_PICK_FLAGS): pick0=(SMEM=true, SPILL=false) full;
    pick1=(true, true) selective band; pick2=(false, false) whole-pool global.

    Pointer params are caller-supplied. The composed sub-inners (minv_inner,
    inverse_dynamics_inner, inverse_dynamics_gradient_inner) are FROZEN and
    placement-free: after the repoint, s_temp already points at the right pool, so
    passing it through is correct with no sub-inner change. The internal Minv/qdd
    smem buffers (s_Minv, s_qdd) are caller-placed too — in the qdd-input variant
    the caller supplies them as inputs; otherwise they are scratch outputs."""
    n = self.robot.get_num_vel()
    func_params = [
        "s_df_du is the output buffer (caller places); size 2*NUM_JOINTS*NUM_JOINTS = " + str(2*n*n),
        "s_q is the vector of joint positions",
        "s_qd is the vector of joint velocities",
    ]
    if use_qdd_Minv_input:
        func_params += [
            "s_qdd is the vector of joint accelerations (input)",
            "s_Minv is the mass matrix (input)",
        ]
    else:
        func_params.append("s_u is the vector of input torques")
    func_params += [
        "s_vaf is the id intermediate band (caller places); size 18*NUM_JOINTS = " + str(18*n),
        "s_dc_du is the id_du output band (caller places); size 2*NUM_JOINTS*NUM_JOINTS = " + str(2*n*n),
        "s_qdd is the joint-accel scratch (caller places); size NUM_JOINTS = " + str(n),
        "s_Minv is the mass-matrix scratch (caller places); size NUM_JOINTS*NUM_JOINTS = " + str(n*n),
        "s_temp is the shared scratch pool (used when SCRATCH_IN_SMEM)",
        "d_workspace is the global scratch pool (used when !SCRATCH_IN_SMEM)",
        "d_temp_spill is the id_du da_df band spill region (used when USE_DA_DF_SPILL)",
        "d_f_ext is the (optional) GLOBAL external forces, body-major 6*NUM_BODIES local-frame, or nullptr",
        "d_robotModel holds XImats/topology; gravity is the gravity constant",
    ]
    fname = "forward_dynamics_gradient_device_qdd" if use_qdd_Minv_input else "forward_dynamics_gradient_device"
    func_def_start = "void " + fname + "(T *s_df_du, const T *s_q, const T *s_qd, "
    if use_qdd_Minv_input:
        func_def_start += "const T *s_qdd, const T *s_Minv, "
        func_def_start += "T *s_vaf, T *s_dc_du, const T *s_qdd_unused, const T *s_Minv_unused, "
    else:
        func_def_start += "const T *s_u, "
        func_def_start += "T *s_vaf, T *s_dc_du, T *s_qdd, T *s_Minv, "
    func_def_end = ("T *s_temp, T *d_workspace, T *d_temp_spill, "
                    "const robotModel<T> *d_robotModel, T *d_f_ext, const T gravity) {")
    func_def_start, func_params = self.gen_insert_helpers_func_def_params(func_def_start, func_params, -2)
    func_def = func_def_start + func_def_end
    self.gen_add_func_doc("fd_du orchestration as a single inner-owns-placement device function",
                          ["Uses the fd/du = -Minv*id/du trick (Carpentier & Mansard 'Analytical Derivatives of Rigid Body Dynamics Algorithms')",
                           "Owns the s_temp pool placement; the repoint covers every consumer below (incl. the XImats helper's sincos scratch and minv's F-region)"],
                          func_params, None)
    self.gen_add_code_line("template <typename T, bool SCRATCH_IN_SMEM = true, bool USE_DA_DF_SPILL = false>")
    # __forceinline__ so the whole orchestration inlines into the calling kernel.
    # Under -rdc a separate __device__ wrapper keeps its callees as distinct
    # functions whose regcount must fit the kernel's launch_bounds budget -> ptxas
    # regcount error. Inlining folds them into the kernel. Mirrors id_du / fdsva_so.
    self.gen_add_code_line("__device__ __forceinline__")
    self.gen_add_code_line(func_def, True)
    if use_qdd_Minv_input:
        self.gen_add_code_line("(void)s_qdd_unused; (void)s_Minv_unused;")
    # Inner owns the pool placement; the repoint covers every consumer below
    # (incl. the XImats helper's sincos scratch + minv's F-region), so no
    # caller-side repoint. The XImats helper call goes AFTER this repoint so its
    # sincos scratch follows the placement (avoids a null-s_temp sincos crash).
    self.gen_add_code_line("if constexpr(!SCRATCH_IN_SMEM){ s_temp = d_workspace; } else { (void)d_workspace; }")
    self.gen_load_update_XImats_helpers_function_call()
    self.gen_forward_dynamics_gradient_inner_python(
        use_qdd_Minv_input,
        s_df_du_name = "s_df_du",
        d_temp_spill_name = "d_temp_spill",
        temp_spill_flag_name = "USE_DA_DF_SPILL")
    self.gen_add_end_function()

def gen_forward_dynamics_gradient_kernel_max_temp_mem_size(self):
    n = self.robot.get_num_vel()
    vaf_cnt = 18 * (self.robot.get_num_joints() if self.robot_has_mimic_joints() else n)
    base_size = 2*n + n*2*n + n*2*n + vaf_cnt + n + n*n + n
    temp_mem_size = self.gen_forward_dynamics_gradient_inner_temp_mem_size()
    return base_size + temp_mem_size

_FD_DU_PICK_FLAGS = [
    # (use_selective_spill, use_global_temp)
    (False, False),   # pick 0: full smem
    (True,  False),   # pick 1: selective spill
    (False, True),    # pick 2: global temp
]

def _emit_fd_du_kernel_body_for_flags(self, n, use_selective_spill, use_global_temp,
                                      use_qdd_Minv_input, single_call_timing):
    """Emit fd_du kernel body for one tier's spill flags."""
    # s_vaf is body-indexed (NB bodies, stride 6). For a MIMIC robot (fixed base)
    # NB > nv, so size 18*NB to keep the ID inner's writes from overflowing into
    # s_qdd/s_Minv. Non-mimic keeps 18*n (byte-identical; floating nv > NB).
    _vaf_cnt = 18 * (self.robot.get_num_joints() if self.robot_has_mimic_joints() else n)
    extra_t_buffers = [("s_q_qd", 2*n+self.robot.floating_base),
                       ("s_dc_du", n*2*n),
                       ("s_vaf", _vaf_cnt),
                       ("s_qdd", n),
                       ("s_Minv", n*n)]
    if not use_qdd_Minv_input:
        extra_t_buffers[0] = ("s_q_qd_u", 3*n+self.robot.floating_base)
    shared_mem_size = 0 if use_global_temp else (
        max(self.gen_minv_inner_temp_mem_size(), self.gen_inverse_dynamics_gradient_temp_layout()["selective_shared_count"])
        if use_selective_spill else self.gen_forward_dynamics_gradient_inner_temp_mem_size()
    )
    self.gen_XImats_helpers_temp_shared_memory_code(shared_mem_size, extra_t_buffers = extra_t_buffers, include_linalg_scratch=True)
    self.gen_add_code_line("T *d_temp_spill = nullptr; (void)d_temp_spill;")
    if use_qdd_Minv_input:
        self.gen_add_code_line(f"T *s_q = s_q_qd; T *s_qd = &s_q_qd[{n+self.robot.floating_base}];")
    else:
        self.gen_add_code_line(f"T *s_q = s_q_qd_u; T *s_qd = &s_q_qd_u[{n+self.robot.floating_base}]; T *s_u = &s_q_qd_u[{2*n+self.robot.floating_base}];")
    if not single_call_timing:
        self.gen_add_parallel_loop("k","NUM_TIMESTEPS",block_level = True)
        if use_qdd_Minv_input:
            self.gen_kernel_load_inputs("q_qd",str(2*n+self.robot.floating_base),"qdd",str(n),"Minv",str(n*n),stride="stride_q_qd",stride2=str(n),stride3=str(n*n))
        else:
            self.gen_kernel_load_inputs("q_qd_u",str(3*n+self.robot.floating_base),stride="stride_q_qd_u")
        # The kernel only SLICES the workspace band pointers; the device owns
        # the s_temp pool placement (the whole-pool global-temp repoint is its
        # SCRATCH_IN_SMEM=false path). Per-rung flags are passed as literals.
        if use_selective_spill or use_global_temp:
            self.gen_add_code_line("T *d_df_du_k = &d_df_du[k*" + str(n*2*n) + "];")
        if use_selective_spill:
            self.gen_add_code_line("d_temp_spill = reinterpret_cast<T *>(&d_workspace[k*GRID_WORKSPACE_BYTES_PER_TIMESTEP<T>()]);")
        self.gen_add_code_line("// compute — the orchestration inner owns its s_temp pool placement")
        self.gen_forward_dynamics_gradient_device_function_call(
            use_qdd_Minv_input,
            scratch_in_smem_expr = ("false" if use_global_temp else "true"),
            use_da_df_spill_expr = ("true" if use_selective_spill else "false"),
            s_df_du_name = ("d_df_du_k" if (use_global_temp or use_selective_spill) else "s_temp"),
            d_workspace_pool_name = ("reinterpret_cast<T *>(&d_workspace[k*GRID_WORKSPACE_BYTES_PER_TIMESTEP<T>()])" if use_global_temp else "nullptr"),
            d_temp_spill_name = ("d_temp_spill" if use_selective_spill else "nullptr"))
        if not (use_global_temp or use_selective_spill):
            self.gen_kernel_save_result("df_du",str(n*2*n),"s_temp",stride=str(n*2*n))
        self.gen_add_end_control_flow()
    else:
        if use_qdd_Minv_input:
            self.gen_kernel_load_inputs("q_qd",str(2*n+self.robot.floating_base),"qdd",str(n),"Minv",str(n*n))
        else:
            self.gen_kernel_load_inputs("q_qd_u",str(3*n+self.robot.floating_base))
        if use_selective_spill or use_global_temp:
            self.gen_add_code_line("T *d_df_du_k = d_df_du;")
        if use_selective_spill:
            self.gen_add_code_line("d_temp_spill = reinterpret_cast<T *>(d_workspace);")
        self.gen_add_code_line("// compute with NUM_TIMESTEPS as NUM_REPS for timing")
        self.gen_add_code_line("for (int rep = 0; rep < NUM_TIMESTEPS; rep++){", True)
        if use_qdd_Minv_input:
            self.gen_anti_licm_input_reload("q_qd", str(2*n + self.robot.floating_base), "qdd", str(n), "Minv", str(n*n))
        else:
            self.gen_anti_licm_input_reload("q_qd_u", str(3*n + self.robot.floating_base))
        # device owns s_temp placement (whole-pool global path = SCRATCH_IN_SMEM=false).
        self.gen_forward_dynamics_gradient_device_function_call(
            use_qdd_Minv_input,
            scratch_in_smem_expr = ("false" if use_global_temp else "true"),
            use_da_df_spill_expr = ("true" if use_selective_spill else "false"),
            s_df_du_name = ("d_df_du_k" if (use_global_temp or use_selective_spill) else "s_temp"),
            d_workspace_pool_name = ("reinterpret_cast<T *>(d_workspace)" if use_global_temp else "nullptr"),
            d_temp_spill_name = ("d_temp_spill" if use_selective_spill else "nullptr"))
        self.gen_add_code_line(
            "if ((threadIdx.x | threadIdx.y | threadIdx.z) == 0) { "
            "reinterpret_cast<volatile T *>(d_df_du)[rep & 63] = "
            + ("d_df_du_k[rep & 63];" if (use_global_temp or use_selective_spill) else "s_temp[rep & 63];")
            + " }"
        )
        self.gen_add_end_control_flow()
        if not (use_global_temp or use_selective_spill):
            self.gen_kernel_save_result("df_du",str(n*2*n),"s_temp")


def gen_forward_dynamics_gradient_kernel(self, use_qdd_Minv_input = False, single_call_timing = False):
    n = self.robot.get_num_vel()
    func_params = ["d_df_du is a pointer to memory for the final result of size 2*NUM_JOINTS*NUM_JOINTS = " + str(2*n*n), \
                   "d_q_dq is the vector of joint positions and velocities", \
                   "stride_q_qd is the stide between each q, qd", \
                   "d_robotModel is the pointer to the initialized model specific helpers on the GPU (XImats, topology_helpers, etc.)", \
                   "d_f_ext is the (optional) GLOBAL external forces, body-major 6*NUM_BODIES local-frame, or nullptr", \
                   "gravity is the gravity constant", \
                   "num_timesteps is the length of the trajectory points we need to compute over (or overloaded as test_iters for timing)"]
    func_notes = []
    func_def_start = "void forward_dynamics_gradient_kernel(T *d_df_du, unsigned char *d_workspace, const T *d_q_qd, const int stride_q_qd, "
    func_def_end = "T *d_f_ext, const robotModel<T> *d_robotModel, const T gravity, const int NUM_TIMESTEPS) {"
    if use_qdd_Minv_input:
        func_def_start += "const T *d_qdd, "
        func_params.insert(-2,"d_qdd is the vector of joint accelerations")
        func_def_start += "const T *d_Minv, "
        func_params.insert(-2,"d_Minv is the mass matrix")
    else:
        func_def_start = func_def_start.replace("_q_qd","_q_qd_u")
        func_params[1] = "d_q_dq is the vector of joint positions, velocities, and input torques"
        func_params[2] = "stride_q_qd_u is the stide between each q, qd, u"
    func_def = func_def_start + func_def_end
    if single_call_timing:
        func_def = func_def.replace("kernel(", "kernel_single_timing(")
    self.gen_add_func_doc("Computes the gradient of forward dynamics",func_notes,func_params,None)
    self.gen_add_code_line("template <typename T, int RESOURCE_TIER = GRID_DEFAULT_RESOURCE_TIER>")
    self.gen_add_code_line("__global__")
    self.gen_add_code_line("__launch_bounds__(tier_max_threads<RESOURCE_TIER>())")
    self.gen_add_code_line(func_def, True)
    picks = getattr(self, "fd_du_spill_tier_3way", (0, 0, 0))
    def _emit_fd_du_body(pick):
        uss, ugt = _FD_DU_PICK_FLAGS[pick]
        _emit_fd_du_kernel_body_for_flags(self, n, uss, ugt, use_qdd_Minv_input, single_call_timing)
    self.gen_tier_dispatch(picks, _emit_fd_du_body)
    self.gen_add_end_function()

def gen_forward_dynamics_gradient_host(self, mode = 0):
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
    func_def_start = "void forward_dynamics_gradient(gridData<T, KIND> *hd_data, const robotModel<T> *d_robotModel, const T gravity, const int num_timesteps,"
    func_def_end =   "                      const dim3 block_dimms, const dim3 thread_dimms, cudaStream_t *streams) {"
    if single_call_timing:
        func_def_start = func_def_start.replace("(", "_single_timing(")
        func_def_end = "              " + func_def_end
    if compute_only:
        func_def_start = func_def_start.replace("(", "_compute_only(")
        func_def_end = "             " + func_def_end.replace(", cudaStream_t *streams", "")
    # then generate the code
    self.gen_add_func_doc("Compute the RNEA (Recursive Newton-Euler Algorithm)",\
                          func_notes,func_params,None)
    self.gen_add_code_line("template <typename T, bool USE_QDD_MINV_FLAG = false, gridDataKind KIND = GRID_DATA_ALL>")
    self.gen_add_code_line("__host__")
    self.gen_add_code_line(func_def_start)
    self.gen_add_code_line(func_def_end, True)
    self.gen_add_code_line("static_assert(KIND == GRID_DATA_ALL || KIND == GRID_DATA_DYNAMICS, \"forward_dynamics_gradient requires all-data or dynamics gridData\");")
    func_call_start = "forward_dynamics_gradient_kernel<T><<<block_dimms,thread_dimms,FD_DU_DYNAMIC_SHARED_MEM_BYTES<T>()>>>(hd_data->d_df_du,hd_data->d_workspace,hd_data->d_q_qd_u,stride_q_qd,"
    func_call_end = "hd_data->d_f_ext,d_robotModel,gravity,num_timesteps);"
    if single_call_timing:
        func_call_start = func_call_start.replace("kernel<T>","kernel_single_timing<T>")
    self.gen_add_code_line("int stride_q_qd= 3*NUM_JOINTS;")
    if not compute_only:
        # start code with memory transfer
        self.gen_add_code_lines(["// start code with memory transfer", \
                                 "gpuErrchk(cudaMemcpyAsync(hd_data->d_q_qd_u,hd_data->h_q_qd_u,stride_q_qd*" + \
                                    ("num_timesteps*" if not single_call_timing else "") + "sizeof(T),cudaMemcpyHostToDevice,streams[0]));", \
                                 "if (USE_QDD_MINV_FLAG) {" ,\
                                 "    gpuErrchk(cudaMemcpyAsync(hd_data->d_qdd,hd_data->h_qdd,NUM_JOINTS*" + \
                                        ("num_timesteps*" if not single_call_timing else "") + "sizeof(T),cudaMemcpyHostToDevice,streams[1]));", \
                                 "    gpuErrchk(cudaMemcpyAsync(hd_data->d_Minv,hd_data->h_Minv,NUM_JOINTS*NUM_JOINTS*" + \
                                        ("num_timesteps*" if not single_call_timing else "") + "sizeof(T),cudaMemcpyHostToDevice,streams[2]));", \
                                 "}", \
                                 "gpuErrchkKernel();"])
    # then compute
    self.gen_add_code_line("// then call the kernel")
    func_call = func_call_start + func_call_end
    func_call_with_qdd_minv = func_call_start + "hd_data->d_qdd, hd_data->d_Minv, " + func_call_end
    func_call_code = ["if (USE_QDD_MINV_FLAG) {" + func_call_with_qdd_minv + "}", "else {" + func_call + "}", "gpuErrchkKernel();"]
    # wrap function call in timing (if needed)
    if single_call_timing:
        func_call_code.insert(0,"struct timespec start, end; clock_gettime(CLOCK_MONOTONIC,&start);")
        func_call_code.append("clock_gettime(CLOCK_MONOTONIC,&end);")
    self.gen_add_code_line("gpuErrchk(grid_check_dynamic_shared_memory_bytes(\"forward_dynamics_gradient\", FD_DU_DYNAMIC_SHARED_MEM_BYTES<T>()));")
    workspace_bytes = "GRID_WORKSPACE_BYTES_PER_TIMESTEP<T>()" if single_call_timing else "GRID_WORKSPACE_BYTES_PER_TIMESTEP<T>()*static_cast<size_t>(num_timesteps)"
    self.gen_add_code_line("if (GRID_FD_DU_USES_WORKSPACE_ANY_TIER) {gpuErrchk(grid_begin_l2_persisting(0, hd_data->d_workspace, " + workspace_bytes + "));}")
    self.gen_add_code_lines(func_call_code)
    self.gen_add_code_line("if (GRID_FD_DU_USES_WORKSPACE_ANY_TIER) {gpuErrchk(grid_end_l2_persisting(0));}")
    if not compute_only:
        # then transfer memory back
        self.gen_add_code_lines(["// finally transfer the result back", \
                                 "gpuErrchk(cudaMemcpy(hd_data->h_df_du,hd_data->d_df_du,NUM_JOINTS*2*NUM_JOINTS*" + \
                                    ("num_timesteps*" if not single_call_timing else "") + "sizeof(T),cudaMemcpyDeviceToHost));",
                                 "gpuErrchkKernel();"])
    # finally report out timing if requested
    if single_call_timing:
        from ..algo_registry import single_call_printf_line
        self.gen_add_code_line(single_call_printf_line("forward_dynamics_gradient"))
    self.gen_add_end_function()

def gen_forward_dynamics_gradient(self):
    # the canonical _device (orchestrator: owns s_temp placement; wraps XImats +
    # minv/id/finish/id + id_du band sub-inner + (-Minv*dc/du); called from kernel
    # and integrator_gradient). Both qdd-Minv-input variants.
    self.gen_forward_dynamics_gradient_device(False)
    self.gen_forward_dynamics_gradient_device(True)
    # then kernels
    self.gen_forward_dynamics_gradient_kernel(True,True)
    self.gen_forward_dynamics_gradient_kernel(True,False)
    self.gen_forward_dynamics_gradient_kernel(False,True)
    self.gen_forward_dynamics_gradient_kernel(False,False)
    # finally host wrappers
    self.gen_forward_dynamics_gradient_host(0)
    self.gen_forward_dynamics_gradient_host(1)
    self.gen_forward_dynamics_gradient_host(2)
