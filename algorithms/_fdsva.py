def gen_fdsva_inner(self, use_thread_group = False): 
	# construct the boilerplate and function definition
    func_params = ["s_fddq_fddqd_fddt are forward dynamics WRT q,qd,tau", \
                "s_q is the vector of joint positions", \
                "s_qd is the vector of joint velocities", \
                "s_qdd is the vector of joint accelerations", \
                "s_tau is the vector of joint torques", \
                "s_id_dq is the first order of inverse dyanmics WRT q", \
                "s_id_dqd is the first order of inverse dyanmics WRT qd", \
                "s_temp is the pointer to the shared memory needed of size: " + \
                            str(self.gen_fdsva_inner_temp_mem_size()), \
                "gravity is the gravity constant"]
    func_def_start = "void fdsva_inner("
    func_def_middle = "T *s_fddq_fddqd_fddt, T *s_q, T *s_qd, const T *s_qdd, const T *s_tau, const T *s_id_dq, const T *s_id_dqd, "
    func_def_end = "T *s_temp, const T gravity) {"
    func_notes = ["Assumes works with IDSVA"]
    if use_thread_group:
        func_def_start = func_def_start.replace("(", "(cgrps::thread_group tgrp, ")
        func_params.insert(0,"tgrp is the handle to the thread_group running this function")
    func_def_middle, func_params = self.gen_insert_helpers_func_def_params(func_def_middle, func_params, -2)
    func_def = func_def_start + func_def_middle + func_def_end
    self.gen_add_func_doc("Forward Dynamics with Spatial Vector Algebra", func_notes, func_params, None)
    self.gen_add_code_line("template <typename T>")
    self.gen_add_code_line("__device__")
    self.gen_add_code_line(func_def, True)

    #
    # Initial Debug Prints if Requested
    #
    if self.DEBUG_MODE:
        self.gen_add_sync(use_thread_group)
        self.gen_add_serial_ops(use_thread_group)
        self.gen_add_code_line("printf(\"q\\n\"); printMat<T,1," + str(n) + ">(s_q,1);")
        self.gen_add_code_line("printf(\"qd\\n\"); printMat<T,1," + str(n) + ">(s_qd,1);")
        self.gen_add_code_line("for (int i = 0; i < " + str(n) + "; i++){printf(\"X[%d]\\n\",i); printMat<T,6,6>(&s_XImats[36*i],6);}")
        self.gen_add_code_line("for (int i = 0; i < " + str(n) + "; i++){printf(\"I[%d]\\n\",i); printMat<T,6,6>(&s_XImats[36*(i+" + str(n) + ")],6);}")
        self.gen_add_end_control_flow()
        self.gen_add_sync(use_thread_group)

    # n = len(qd)
    n = self.robot.get_num_pos()
    n_bfs_levels = self.robot.get_max_bfs_level() + 1 # starts at 0

    # fddq = np.zeros((n,n))
    # fddqd = np.zeros((n,n))
    # fddt = np.zeros((n,n))

    # fddq_offset = 0
    # self.gen_add_code_line("T * fddq = &s_fddq_fddqd_fddt[" + str(fddq_offset) + "];")
    # fddqd_offset = fddq_offset + n*n
    # self.gen_add_code_line("T *fddqd = &s_fddq_fddqd_fddt[" + str(fddqd_offset) + "];")
    # fddt_offset = fddqd_offset + n*n
    # self.gen_add_code_line("T * fddt = &s_fddq_fddqd_fddt[" + str(fddt_offset) + "];")

    M_offset = 0
    self.gen_add_code_line("T * M = &s_temp[" + str(M_offset) + "];")

    # Minv = np.linalg.inv(id_M)
    #import M from crba
    # self.gen_add_code_line("M = &s_temp[" + str(M_offset) + "];")

    #inverse the matrix 

    # fddq = iddq
    self.gen_add_code_line("T * fddq = &s_id_dq[0];")
    # fddqd = iddqd
    self.gen_add_code_line("T * fddqd = &s_id_dqd[0];")


    # fddq = (-1)*Minv@fddq
    # fddqd = (-1)*Minv@fddqd
    # fddt = (-1)*Minv 
    self.gen_add_end_function()

    

    
def gen_fdsva_inner_temp_mem_size(self):
    n = self.robot.get_num_pos()
    return 140 * n
    
def gen_fdsva_inner_function_call(self, use_thread_group = False, updated_var_names = None):
    var_names = dict( \
        s_fddq_fddqd_fddt_name = "s_fddq_fddqd_fddt", \
        s_q_name = "s_q", \
        s_qd_name = "s_qd", \
        s_qdd_name = "s_qdd", \
        s_tau_name = "s_tau", \
        s_id_dq_name = "s_id_dq", \
        s_id_dqd_name = "s_id_dqd", \
        s_temp_name = "s_temp", \
        gravity_name = "gravity"
    )
    if updated_var_names is not None:
        for key,value in updated_var_names.items():
            var_names[key] = value
    fdsva_code_start = "fdsva_inner<T>(" + var_names["s_fddq_fddqd_fddt_name"] + ", " + var_names["s_q_name"] + ", " + var_names["s_qd_name"] + ", " + var_names["s_qdd_name"] + ", " + var_names["s_tau_name"] + ", " + var_names["s_id_dq_name"] + ", " + var_names["s_id_dqd_name"] + ", "
    fdsva_code_end = var_names["s_temp_name"] + ", " + var_names["gravity_name"] + ");"
    if use_thread_group:
        id_code_start = id_code_start.replace("(","(tgrp, ")
    fdsva_code_middle = self.gen_insert_helpers_function_call()
    fdsva_code = fdsva_code_start + fdsva_code_middle + fdsva_code_end
    self.gen_add_code_line(fdsva_code)

def gen_fdsva_device_temp_mem_size(self):
    n = self.robot.get_num_pos()
    wrapper_size = self.gen_topology_helpers_size() + 72*n # for XImats
    return self.gen_fdsva_inner_temp_mem_size() + wrapper_size

def gen_fdsva_device(self, use_thread_group = False):
    n = self.robot.get_num_pos()
    # construct the boilerplate and function definition
    func_params = ["s_q is the vector of joint positions", \
                   "s_qd is the vector of joint velocities", \
                   "s_qdd is the vector of joint accelerations", \
                   "s_tau is the vector of joint torques", \
                   "s_id_dq is the first order of inverse dyanmics WRT q", \
                   "s_id_dqd is the first order of inverse dyanmics WRT qd", \
                   "d_robotModel is the pointer to the initialized model specific helpers on the GPU (XImats, topology_helpers, etc.)", \
                   "gravity is the gravity constant"]
    func_notes = []
    func_def_start = "void fdsva_device("
    func_def_middle = "const T *s_q, const T *s_qd, const T *s_qdd, const T *s_tau, const T *s_id_dq, const T *s_id_dqd,"
    func_def_end = "const robotModel<T> *d_robotModel, const T gravity) {"
    if use_thread_group:
        func_def_start += "cgrps::thread_group tgrp, "
        func_params.insert(0,"tgrp is the handle to the thread_group running this function")
    func_def = func_def_start + func_def_middle + func_def_end

    # then generate the code
    self.gen_add_func_doc("Compute the FDSVA (Forward Dyamics with Spacial Vector Algebra)",\
                          func_notes,func_params,None)
    self.gen_add_code_line("template <typename T>")
    self.gen_add_code_line("__device__")
    self.gen_add_code_line(func_def, True)

    # add the shared memory variables
    shared_mem_size = self.gen_fdsva_device_temp_mem_size() if not self.use_dynamic_shared_mem_flag else None
    self.gen_XImats_helpers_temp_shared_memory_code(shared_mem_size)

    self.gen_add_code_line("extern __shared__ T s_fddq_fddqd_fddt[" + str(n*n*3) + "];")
    
    # then load/update XI and run the algo
    self.gen_load_update_XImats_helpers_function_call(use_thread_group)
    self.gen_fdsva_inner_function_call(use_thread_group)
    self.gen_add_end_function()

def gen_fdsva_kernel(self, use_thread_group = False, single_call_timing = False):
    n = self.robot.get_num_pos()
    # define function def and params
    func_params = ["d_q_qd_qdd_tau is the vector of joint positions, velocities, accelerations", \
                    "stride_q_qd_qdd is the stride between each q, qd, qdd", \
                    "d_id_dq_dqd is the first order of inverse dyanmics WRT q and qd", \
                    "d_robotModel is the pointer to the initialized model specific helpers on the GPU (XImats, topology_helpers, etc.)", \
                    "d_tau is the vector of joint torques", \
                    "gravity is the gravity constant", \
                    "num_timesteps is the length of the trajectory points we need to compute over (or overloaded as test_iters for timing)"]
    func_notes = []
    func_def_start = "void fdsva_kernel(T *d_fddq_fddqd_fddt, const T *d_q_qd_qdd_tau, const int stride_q_qd_qdd, const T *d_id_dq_dqd, "
    func_def_end = "const robotModel<T> *d_robotModel, const T gravity, const int NUM_TIMESTEPS) {"
    func_def = func_def_start + func_def_end
    if single_call_timing:
        func_def = func_def.replace("kernel(", "kernel_single_timing(")
    
    # then generate the code
    self.gen_add_func_doc("Compute the FDSVA (Forward Dynamics with Spacial Vector Algebra)", \
                            func_notes, func_params, None)
    self.gen_add_code_line("template <typename T>")
    self.gen_add_code_line("__global__")
    self.gen_add_code_line(func_def, True)

    # add shared memory variables
    shared_mem_vars = ["__shared__ T s_fddq_fddqd_fddt[3*" + str(n*n) + "];", \
                        "__shared__ T s_q_qd_qdd_tau[4*" + str(n) + "]; T *s_q = s_q_qd_qdd_tau; T *s_qd = &s_q_qd_qdd_tau[" + str(n) + "]; T *s_qdd = &s_q_qd_qdd_tau[2 * " + str(n) + "]; T *s_tau = &s_q_qd_qdd_tau[3 * " + str(n) + "];",\
                        "__shared__ T s_id_dq_dqd[2* " + str(n*n) + "]; T *s_id_dq = s_id_dq_dqd; T *s_id_dqd = &s_id_dq_dqd[" + str(n*n) + "];"]
    self.gen_add_code_lines(shared_mem_vars)
    shared_mem_size = self.gen_fdsva_inner_temp_mem_size() if not self.use_dynamic_shared_mem_flag else None
    self.gen_XImats_helpers_temp_shared_memory_code(shared_mem_size)
    if use_thread_group:
        self.gen_add_code_line("cgrps::thread_group tgrp = TBD;")
    if not single_call_timing:
        # load to shared mem and loop over blocks to compute all requested comps
        self.gen_add_parallel_loop("k","NUM_TIMESTEPS",use_thread_group,block_level = True)
        self.gen_kernel_load_inputs("q_qd_qdd_tau","stride_q_qd_qdd",str(3*n),use_thread_group)
        # compute
        self.gen_add_code_line("// compute")
        self.gen_load_update_XImats_helpers_function_call(use_thread_group)
        self.gen_fdsva_inner_function_call(use_thread_group)
        self.gen_add_sync(use_thread_group)
        # save to global
        self.gen_kernel_save_result("fddq_fddqd_fddt","1",str(n*n*3),use_thread_group)
        self.gen_add_end_control_flow()
    else:
        # repurpose NUM_TIMESTEPS for number of timing reps
        self.gen_kernel_load_inputs_single_timing("q_qd_qdd_tau",str(4*n),use_thread_group)
        # then compute in loop for timing
        self.gen_add_code_line("// compute with NUM_TIMESTEPS as NUM_REPS for timing")
        self.gen_add_code_line("for (int rep = 0; rep < NUM_TIMESTEPS; rep++){", True)
        self.gen_load_update_XImats_helpers_function_call(use_thread_group)
        self.gen_fdsva_inner_function_call(use_thread_group)
        self.gen_add_end_control_flow()
        # save to global
        self.gen_kernel_save_result_single_timing("fddq_fddqd_fddt",str(n*n*3),use_thread_group)
    self.gen_add_end_function()

def gen_fdsva_host(self, mode = 0):
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
    func_def_start = "void fdsva(gridData<T> *hd_data, const robotModel<T> *d_robotModel, const T gravity, const int num_timesteps,"
    func_def_end =   "                      const dim3 block_dimms, const dim3 thread_dimms, cudaStream_t *streams) {"
    if single_call_timing:
        func_def_start = func_def_start.replace("(", "_single_timing(")
        func_def_end = "              " + func_def_end
    if compute_only:
        func_def_start = func_def_start.replace("(", "_compute_only(")
        func_def_end = "             " + func_def_end.replace(", cudaStream_t *streams", "")
    # then generate the code
    self.gen_add_func_doc("Compute the FDSVA (Forward Dynamics with Spacial Vector Algebra)",\
                          func_notes,func_params,None)
    self.gen_add_code_line("template <typename T>")
    self.gen_add_code_line("__host__")
    self.gen_add_code_line(func_def_start)
    self.gen_add_code_line(func_def_end, True)

    func_call_start = "fdsva_kernel<T><<<block_dimms,thread_dimms,FDSVA_SHARED_MEM_COUNT*sizeof(T)>>>(hd_data->d_fddq_fddqd_fddt,hd_data->d_q_qd_u,stride_q_qd_qdd,hd_data->d_id_dq_dqd,"
    func_call_end = "d_robotModel,gravity,num_timesteps);"
    self.gen_add_code_line("int stride_q_qd_qdd = 3*NUM_JOINTS;")
    if single_call_timing:
        func_call_start = func_call_start.replace("kernel<T>","kernel_single_timing<T>")
    if not compute_only:
        # start code with memory transfer
        self.gen_add_code_lines(["// start code with memory transfer", \
                                 "gpuErrchk(cudaMemcpyAsync(hd_data->d_q_qd_u,hd_data->h_q_qd_u,stride_q_qd_qdd*" + \
                                    ("num_timesteps*" if not single_call_timing else "") + "sizeof(T),cudaMemcpyHostToDevice,streams[0]));", \
                                 "gpuErrchk(cudaDeviceSynchronize());"])
    # then compute:
    self.gen_add_code_line("// call the kernel")
    func_call = func_call_start + func_call_end
    func_call_code = [func_call, "gpuErrchk(cudaDeviceSynchronize());"]
    # wrap function call in timing (if needed)
    if single_call_timing:
        func_call_code.insert(0,"struct timespec start, end; clock_gettime(CLOCK_MONOTONIC,&start);")
        func_call_code.append("clock_gettime(CLOCK_MONOTONIC,&end);")
    self.gen_add_code_lines(func_call_code)
    if not compute_only:
        # then transfer memory back
        self.gen_add_code_lines(["// finally transfer the result back", \
                                "gpuErrchk(cudaMemcpy(hd_data->h_fddq_fddqd_fddt,hd_data->d_fddq_fddqd_fddt,NUM_JOINTS*" + \
                                ("num_timesteps*" if not single_call_timing else "") + "sizeof(T),cudaMemcpyDeviceToHost));",
                                "gpuErrchk(cudaDeviceSynchronize());"])
    # finally report out timing if requested
    if single_call_timing:
        self.gen_add_code_line("printf(\"Single Call FDSVA %fus\\n\",time_delta_us_timespec(start,end)/static_cast<double>(num_timesteps));")
    self.gen_add_end_function()

def gen_fdsva(self, use_thread_group = False):
    # first generate the inner helper
    self.gen_fdsva_inner(use_thread_group)
    # then generate the device wrapper
    self.gen_fdsva_device(use_thread_group)
    # then generate the kernels
    self.gen_fdsva_kernel(use_thread_group, True)
    self.gen_fdsva_kernel(use_thread_group, False)
    # then generate the host wrappers
    self.gen_fdsva_host(0)
    self.gen_fdsva_host(1)
    self.gen_fdsva_host(2)
    