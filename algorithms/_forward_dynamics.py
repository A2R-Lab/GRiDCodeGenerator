def gen_forward_dynamics_inner_F_size(self):
    """Phase 3b: size of Minv's F-region used inside forward_dynamics_inner
    (6 * NV * NV floats). Now always a separate `s_minv_F` parameter so the
    caller can keep it in smem (extra t_buffer) or spill it to L2-pinned
    workspace on humanoid-scale robots."""
    n = self.robot.get_num_vel()
    return 6 * n * n

def gen_forward_dynamics_inner_temp_mem_size(self):
        """Phase 3b: s_temp arena = s_Minv (n*n, persistent) + max(Minv's no_F
        portion during the Minv call, c+vaf+ID-inner after Minv). Minv's F
        moved to a separate s_minv_F parameter (above)."""
        n = self.robot.get_num_pos()
        return n*n + max(self.gen_direct_minv_inner_no_F_size(),
                         19*n + self.gen_inverse_dynamics_inner_temp_mem_size())

def gen_forward_dynamics_finish_function_call(self, use_thread_group = False, updated_var_names = None):
    var_names = dict( \
        s_qdd_name = "s_qdd", \
        s_u_name = "s_u", \
        s_c_name = "s_c", \
        s_Minv_name = "s_Minv"
    )
    if updated_var_names is not None:
        for key,value in updated_var_names.items():
            var_names[key] = value
    code = "forward_dynamics_finish<T>(" + var_names["s_qdd_name"] + ", " + var_names["s_u_name"] + ", " + \
                                           var_names["s_c_name"] + ", " + var_names["s_Minv_name"] + ");"
    if use_thread_group:
        code = code.replace("(","(tgrp, ")
    self.gen_add_code_line(code)

def gen_forward_dynamics_finish(self, use_thread_group = False):
    n = self.robot.get_num_vel()
    # construct the boilerplate and function definition
    func_params = ["s_qdd is a pointer to memory for the final result", \
                   "s_u is the vector of joint input torques", \
                   "s_c is the bias vector", \
                   "s_Minv is the inverse mass matrix"]
    func_def = "void forward_dynamics_finish(T *s_qdd, const T *s_u, const T *s_c, const T *s_Minv) {"
    func_notes = ["Assumes s_Minv and s_c are already computed", 
                  "Does not internally sync the thread group, so it should be called after all threads have finished computing their values"]
    if use_thread_group:
        func_def = func_def.replace("(","(cgrps::thread_group tgrp, ")
        func_params.insert(0,"tgrp is the handle to the thread_group running this function")
    self.gen_add_func_doc("Finish the forward dynamics computation with qdd = Minv*(u-c)",func_notes,func_params,None)
    self.gen_add_code_line("template <typename T>")
    self.gen_add_code_line("__device__")
    self.gen_add_code_line(func_def, True)

    # compute the final answer qdd = Minv * (u - c)
    # remember that Minv is an SYMMETRIC_UPPER triangular matrix
    self.gen_add_parallel_loop("row",str(n),use_thread_group)
    self.gen_add_code_line("T val = static_cast<T>(0);")
    self.gen_add_code_line("for(int col = 0; col < " + str(n) + "; col++) {", True)
    self.gen_add_code_line("// account for the fact that Minv is an SYMMETRIC_UPPER triangular matrix")
    self.gen_add_code_line("int index = (row <= col) * (col * " + str(n) + " + row) + (row > col) * (row * " + str(n) + " + col);")
    self.gen_add_code_line("val += s_Minv[index] * (s_u[col] - s_c[col]);")
    self.gen_add_end_control_flow()
    self.gen_add_code_line("s_qdd[row] = val;")
    self.gen_add_end_control_flow()
    self.gen_add_end_function()

def gen_forward_dynamics_inner_function_call(self, use_thread_group = False, updated_var_names = None):
    var_names = dict( \
        s_q_name = "s_q", \
        s_qd_name = "s_qd", \
        s_qdd_name = "s_qdd", \
        s_u_name = "s_u", \
        s_minv_F_name = "s_minv_F", \
        s_temp_name = "s_temp", \
        gravity_name = "gravity"
    )
    if updated_var_names is not None:
        for key,value in updated_var_names.items():
            var_names[key] = value
    # Phase 3b: forward_dynamics_inner takes s_minv_F (6*NV*NV) as a separate
    # arg so the caller can choose smem or L2-pinned workspace for it.
    fd_code_start = "forward_dynamics_inner<T>(" + var_names["s_qdd_name"] + ", " + var_names["s_q_name"] + ", " + \
                                                   var_names["s_qd_name"] + ", " + var_names["s_u_name"] + ", " + \
                                                   var_names["s_minv_F_name"] + ", "
    fd_code_end = var_names["s_temp_name"] + ", " + var_names["gravity_name"] + ");"
    fd_code_middle = self.gen_insert_helpers_function_call()
    if use_thread_group:
        fd_code_start = fd_code_start.replace("(","(tgrp, ")
    fd_code = fd_code_start + fd_code_middle + fd_code_end
    self.gen_add_code_line(fd_code)

def gen_forward_dynamics_inner(self, use_thread_group = False):
    n = self.robot.get_num_vel()
    NJ = self.robot.get_num_joints()
    # construct the boilerplate and function definition
    func_params = ["s_qdd is a pointer to memory for the final result", \
                   "s_q is the vector of joint positions", \
                   "s_qd is the vector of joint velocities", \
                   "s_u is the vector of joint input torques", \
                   "s_minv_F is a pointer to the 6*NV*NV Minv-F scratch (Phase 3b: separate so caller can spill to L2-pinned workspace)", \
                   "s_temp is the pointer to the shared memory needed of size: " + \
                            str(self.gen_forward_dynamics_inner_temp_mem_size()), \
                   "gravity is the gravity constant"]
    func_def_start = "void forward_dynamics_inner(T *s_qdd, const T *s_q, const T *s_qd, const T *s_u, T *s_minv_F, "
    func_def_end = "T *s_temp, const T gravity) {"
    func_def_start, func_params = self.gen_insert_helpers_func_def_params(func_def_start, func_params, -2)
    func_notes = ["Assumes s_XImats is updated already for the current s_q",
                  "Does not internally sync the thread group, so it should be called after all threads have finished computing their values"]
    if use_thread_group:
        func_def_start = func_def_start.replace("(","(cgrps::thread_group tgrp, ")
        func_params.insert(0,"tgrp is the handle to the thread_group running this function")
    func_def = func_def_start + func_def_end
    # then generate the code
    self.gen_add_func_doc("Computes forward dynamics",func_notes,func_params,None)
    self.gen_add_code_line("template <typename T>")
    self.gen_add_code_line("__device__")
    self.gen_add_code_line(func_def, True)
    # Phase 3b layout: s_minv_F is a separate param. s_temp[0..n*n] = s_Minv
    # (persistent); s_temp[n*n..n*n+no_F] = Minv's no_F (during Minv only);
    # after Minv that region is reused for c+vaf+ID_inner.
    updated_var_names = dict(s_Minv_name = "s_temp",
                             s_F_name = "s_minv_F",
                             s_temp_name = "&s_temp[" + str(n*n) + "]")
    self.gen_direct_minv_inner_function_call(use_thread_group, updated_var_names)
    updated_var_names = dict(s_c_name = "&s_temp[" + str(n*n) + "]", s_vaf_name = "&s_temp[" + str(n*n + n) + "]", s_temp_name = "&s_temp[" + str(n*n + n + 18*NJ) + "]")
    self.gen_inverse_dynamics_inner_function_call(use_thread_group, compute_c = True, use_qdd_input = False, updated_var_names = updated_var_names)
    
    if self.DEBUG_MODE:
        self.gen_add_sync(use_thread_group)
        self.gen_add_serial_ops(use_thread_group)
        self.gen_add_code_lines(["printf(\"Minv\\n\"); printMat<T," + str(n) + "," + str(n) + ">(s_temp," + str(n) + ");",
                                 "printf(\"u\\n\"); printMat<T,1," + str(n) + ">(s_u,1);"
                                 "printf(\"c\\n\"); printMat<T,1," + str(n) + ">(&s_temp[" + str(n*n) + "],1);"])
        self.gen_add_end_control_flow()
        self.gen_add_sync(use_thread_group)

    # finally compute the final answer qdd = Minv * (u - c)
    updated_var_names = dict(s_Minv_name = "s_temp", s_c_name = "&s_temp[" + str(n*n) + "]")
    self.gen_forward_dynamics_finish_function_call(use_thread_group, updated_var_names)
    self.gen_add_end_function()

def gen_forward_dynamics_device(self, use_thread_group = False):
    n = self.robot.get_num_vel()
    # construct the boilerplate and function definition
    func_params = ["s_qdd is a pointer to memory for the final result", \
                   "s_q is the vector of joint positions", \
                   "s_qd is the vector of joint velocities", \
                   "s_u is the vector of joint input torques", \
                   "d_robotModel is the pointer to the initialized model specific helpers on the GPU (XImats, topology_helpers, etc.)", \
                   "gravity is the gravity constant"]
    func_def_start = "void forward_dynamics_device(T *s_qdd, const T *s_q, const T *s_qd, const T *s_u, "
    func_def_end = "const robotModel<T> *d_robotModel, const T gravity) {"
    func_notes = []
    if use_thread_group:
        func_def_start = func_def_start.replace("(","(cgrps::thread_group tgrp, ")
        func_params.insert(0,"tgrp is the handle to the thread_group running this function")
    func_def = func_def_start + func_def_end
    # then generate the code
    self.gen_add_func_doc("Computes forward dynamics",func_notes,func_params,None)
    self.gen_add_code_line("template <typename T>")
    self.gen_add_code_line("__device__")
    self.gen_add_code_line(func_def, True)
    # Phase 3b: allocate s_minv_F (6*NV*NV) alongside s_temp in smem. Inline-CUDA
    # device path keeps F in smem (no surgical spill at this layer).
    F_size = self.gen_forward_dynamics_inner_F_size()
    shared_mem_size = self.gen_forward_dynamics_inner_temp_mem_size() + F_size
    self.gen_XImats_helpers_temp_shared_memory_code(shared_mem_size, include_linalg_scratch=True)
    self.gen_add_code_line("T *s_minv_F = s_temp;")
    self.gen_add_code_line("T *fd_s_temp = &s_temp[" + str(F_size) + "];")
    # then load/update XI and run the algo
    self.gen_load_update_XImats_helpers_function_call(use_thread_group)
    self.gen_forward_dynamics_inner_function_call(use_thread_group,
        updated_var_names = dict(s_minv_F_name = "s_minv_F", s_temp_name = "fd_s_temp"))
    self.gen_add_end_function()

def _emit_fd_kernel_body_for_flags(self, n, spill_minv_F, single_call_timing, use_thread_group):
    """Emit forward_dynamics_kernel body for one tier's Minv-F spill flag.
    spill_minv_F=False: s_minv_F lives in extra smem (at start of s_temp);
    spill_minv_F=True:  s_minv_F lives in L2-pinned workspace."""
    F_size = self.gen_forward_dynamics_inner_F_size()
    if spill_minv_F:
        # F in workspace; s_temp arena holds only Minv_no_F + c+vaf+ID_inner
        shared_mem_size = self.gen_forward_dynamics_inner_temp_mem_size()
    else:
        # F embedded at start of s_temp; arena is bigger by F_size
        shared_mem_size = self.gen_forward_dynamics_inner_temp_mem_size() + F_size
    self.gen_XImats_helpers_temp_shared_memory_code(shared_mem_size, extra_t_buffers = [("s_q_qd_u", 3*n+self.robot.floating_base), ("s_qdd", n)], include_linalg_scratch=True)
    self.gen_add_code_line("T *s_q = s_q_qd_u; T *s_qd = &s_q_qd_u[" + str(n+self.robot.floating_base) + "]; T *s_u = &s_q_qd_u[" + str(2*n+self.robot.floating_base) + "];")
    if use_thread_group:
        self.gen_add_code_line("cgrps::thread_group tgrp = TBD;")
    if not single_call_timing:
        self.gen_add_parallel_loop("k","NUM_TIMESTEPS",use_thread_group,block_level = True)
        self.gen_kernel_load_inputs("q_qd_u","stride_q_qd_u",str(3*n),use_thread_group)
        if spill_minv_F:
            self.gen_add_code_line("T *s_minv_F = reinterpret_cast<T *>(&d_workspace[k*GRID_WORKSPACE_BYTES_PER_TIMESTEP<T>() + GRID_MINV_F_WORKSPACE_OFFSET_BYTES<T>()]);")
            self.gen_add_code_line("T *fd_s_temp = s_temp;")
        else:
            self.gen_add_code_line("(void)d_workspace;")
            self.gen_add_code_line("T *s_minv_F = s_temp;")
            self.gen_add_code_line("T *fd_s_temp = &s_temp[" + str(F_size) + "];")
        self.gen_add_code_line("// compute")
        self.gen_load_update_XImats_helpers_function_call(use_thread_group)
        self.gen_forward_dynamics_inner_function_call(use_thread_group,
            updated_var_names = dict(s_minv_F_name = "s_minv_F", s_temp_name = "fd_s_temp"))
        self.gen_add_sync(use_thread_group)
        self.gen_kernel_save_result("qdd",str(n),str(n),use_thread_group)
        self.gen_add_end_control_flow()
    else:
        input_count = 3*n + self.robot.floating_base
        self.gen_kernel_load_inputs_single_timing("q_qd_u",str(input_count))
        if spill_minv_F:
            self.gen_add_code_line("T *s_minv_F = reinterpret_cast<T *>(&d_workspace[GRID_MINV_F_WORKSPACE_OFFSET_BYTES<T>()]);")
            self.gen_add_code_line("T *fd_s_temp = s_temp;")
        else:
            self.gen_add_code_line("(void)d_workspace;")
            self.gen_add_code_line("T *s_minv_F = s_temp;")
            self.gen_add_code_line("T *fd_s_temp = &s_temp[" + str(F_size) + "];")
        self.gen_add_code_line("// compute with NUM_TIMESTEPS as NUM_REPS for timing")
        self.gen_add_code_line("for (int rep = 0; rep < NUM_TIMESTEPS; rep++){", True)
        self.gen_anti_licm_input_reload("q_qd_u",str(input_count),use_thread_group,feedback_from="qdd")
        self.gen_load_update_XImats_helpers_function_call(use_thread_group)
        self.gen_forward_dynamics_inner_function_call(use_thread_group,
            updated_var_names = dict(s_minv_F_name = "s_minv_F", s_temp_name = "fd_s_temp"))
        self.gen_anti_licm_output_write("qdd")
        self.gen_add_end_control_flow()
        # save to global
        self.gen_kernel_save_result_single_timing("qdd",str(n),use_thread_group)


def gen_forward_dynamics_kernel(self, use_thread_group = False, single_call_timing = False):
    n = self.robot.get_num_vel()
    func_params = ["d_qdd is a pointer to memory for the final result", \
                   "d_workspace is the L2-pinned global spill buffer (used when Minv-F overflows smem)", \
                   "d_q_qd_u is the vector of joint positions, velocities, and input torques", \
                   "d_robotModel is the pointer to the initialized model specific helpers on the GPU (XImats, topology_helpers, etc.)", \
                   "gravity is the gravity constant", \
                   "num_timesteps is the length of the trajectory points we need to compute over (or overloaded as test_iters for timing)"]
    func_def_start = "void forward_dynamics_kernel(T *d_qdd, unsigned char *d_workspace, const T *d_q_qd_u, const int stride_q_qd_u, "
    func_def_end = "const robotModel<T> *d_robotModel, const T gravity, const int NUM_TIMESTEPS) {"
    func_notes = []
    func_def = func_def_start + func_def_end
    if single_call_timing:
        func_def = func_def.replace("kernel(", "kernel_single_timing(")
    self.gen_add_func_doc("Computes forward dynamics",func_notes,func_params,None)
    self.gen_add_code_line("template <typename T, int RESOURCE_TIER = TIER_PERF>")
    self.gen_add_code_line("__global__")
    self.gen_add_code_line("__launch_bounds__(tier_max_threads<RESOURCE_TIER>())")
    self.gen_add_code_line(func_def, True)
    # Phase 3b: 3-way pick dispatch. Level 0 = Minv-F embedded in s_temp (extra
    # smem block); Level 1 = Minv-F in L2-pinned workspace.
    picks = getattr(self, "fd_spill_tier_3way", (0, 0, 0))
    if picks[0] == picks[1] == picks[2]:
        _emit_fd_kernel_body_for_flags(self, n, bool(picks[0]), single_call_timing, use_thread_group)
    else:
        tier_names = ("TIER_PERF", "TIER_LITE", "TIER_MINIMAL")
        for tier_idx, (tier_name, pick) in enumerate(zip(tier_names, picks)):
            head = "if constexpr (RESOURCE_TIER == " + tier_name + ") {" if tier_idx == 0 else \
                   "else if constexpr (RESOURCE_TIER == " + tier_name + ") {"
            self.gen_add_code_line(head, True)
            _emit_fd_kernel_body_for_flags(self, n, bool(pick), single_call_timing, use_thread_group)
            self.gen_add_end_control_flow()
    self.gen_add_end_function()


def gen_forward_dynamics_host(self, mode = 0):
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
    func_def_start = "void forward_dynamics(gridData<T, KIND> *hd_data, const robotModel<T> *d_robotModel, const T gravity, const int num_timesteps,"
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
    self.gen_add_code_line("template <typename T, gridDataKind KIND = GRID_DATA_ALL>")
    self.gen_add_code_line("__host__")
    self.gen_add_code_line(func_def_start)
    self.gen_add_code_line(func_def_end, True)
    self.gen_add_code_line("static_assert(KIND == GRID_DATA_ALL || KIND == GRID_DATA_DYNAMICS, \"forward_dynamics requires all-data or dynamics gridData\");")
    func_call_start = "forward_dynamics_kernel<T><<<block_dimms,thread_dimms,FD_DYNAMIC_SHARED_MEM_BYTES<T>()>>>(hd_data->d_qdd,hd_data->d_workspace,hd_data->d_q_qd_u,stride_q_qd_u,"
    func_call_end = "d_robotModel,gravity,num_timesteps);"
    if single_call_timing:
        func_call_start = func_call_start.replace("kernel<T>","kernel_single_timing<T>")
    self.gen_add_code_line("int stride_q_qd_u = 3*NUM_JOINTS;")
    if not compute_only:
        # start code with memory transfer
        self.gen_add_code_lines(["// start code with memory transfer", \
                                 "gpuErrchk(cudaMemcpyAsync(hd_data->d_q_qd_u,hd_data->h_q_qd_u,stride_q_qd_u*" + \
                                    ("num_timesteps*" if not single_call_timing else "") + "sizeof(T),cudaMemcpyHostToDevice,streams[0]));", \
                                 "gpuErrchkKernel();"])
    # then compute
    self.gen_add_code_line("// then call the kernel")
    func_call = func_call_start + func_call_end
    func_call_code = [func_call, "gpuErrchkKernel();"]
    # wrap function call in timing (if needed)
    if single_call_timing:
        func_call_code.insert(0,"struct timespec start, end; clock_gettime(CLOCK_MONOTONIC,&start);")
        func_call_code.append("clock_gettime(CLOCK_MONOTONIC,&end);")
    self.gen_add_code_line("gpuErrchk(grid_check_dynamic_shared_memory_bytes(\"forward_dynamics\", FD_DYNAMIC_SHARED_MEM_BYTES<T>()));")
    self.gen_add_code_lines(func_call_code)
    if not compute_only:
        # then transfer memory back
        self.gen_add_code_lines(["// finally transfer the result back", \
                                 "gpuErrchk(cudaMemcpy(hd_data->h_qdd,hd_data->d_qdd,NUM_JOINTS*" + \
                                    ("num_timesteps*" if not single_call_timing else "") + "sizeof(T),cudaMemcpyDeviceToHost));",
                                 "gpuErrchkKernel();"])
    # finally report out timing if requested
    if single_call_timing:
        from ..algo_registry import single_call_printf_line
        self.gen_add_code_line(single_call_printf_line("fd"))
    self.gen_add_end_function()

def gen_forward_dynamics(self, use_thread_group = False):
    # first helpers
    self.gen_forward_dynamics_finish(use_thread_group)
    self.gen_forward_dynamics_inner(use_thread_group)
    # then device wrapper
    self.gen_forward_dynamics_device(use_thread_group)
    # then kernels
    self.gen_forward_dynamics_kernel(use_thread_group,True)
    self.gen_forward_dynamics_kernel(use_thread_group,False)
    # then host launch
    self.gen_forward_dynamics_host(0)
    self.gen_forward_dynamics_host(1)
    self.gen_forward_dynamics_host(2)
