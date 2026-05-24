"""
End Effector Posiitons

TODO: fix throughout this document for fixed_joint support for branched trees and for multiple fixed at once
"""
def gen_end_effector_pose_inner_temp_mem_size(self, fixed_target_name = ""):
    num_ees = self.robot.get_total_leaf_nodes() if fixed_target_name == "" else 1
    return 2*16*num_ees

def gen_end_effector_pose_inner_function_call(self, use_thread_group = False, updated_var_names = None, fixed_target_name = ""):
    var_names = dict( \
        s_Xhom_name = "s_XmatsHom", \
        s_eePos_name = "s_eePos", \
        s_q_name = "s_q", \
        s_topology_helpers_name = "s_topology_helpers", \
        s_temp_name = "s_temp", \
        s_linalg_smem_name = "s_linalg_smem", \
    )
    if updated_var_names is not None:
        for key,value in updated_var_names.items():
            var_names[key] = value
    code_start = "end_effector_pose_inner" + ("" if fixed_target_name == "" else "_" + fixed_target_name) + "<T>(" + var_names["s_eePos_name"] + ", " + var_names["s_q_name"] + ", "
    code_middle = var_names["s_Xhom_name"] + ", "
    code_end =  var_names["s_temp_name"] + ", " + var_names["s_linalg_smem_name"] + ");"
    # account for thread group and serial chains
    if use_thread_group:
        code_start = code_start.replace("(","(tgrp, ")
    # Canonical: append the shared topology-helper arg via the central helper
    # (NO_XI: the ee_pose family takes s_Xhom, not s_XImats). Mirrors the def's
    # gen_insert_helpers_func_def_params(NO_XI_FLAG=True) so def + call can't drift.
    code_middle += self.gen_insert_helpers_function_call(updated_var_names = var_names, NO_XI_FLAG = True)
    self.gen_add_code_line(code_start + code_middle + code_end)

def gen_end_effector_pose_inner(self, use_thread_group = False, fixed_target_name = ""):
    n = self.robot.get_num_pos()
    n_bfs_levels = self.robot.get_max_bfs_level() + 1 # starts at 0
    if fixed_target_name == "":
        all_ees = self.robot.get_leaf_nodes()
    else:
        all_ees = [self.robot.get_fixed_joint_by_name(fixed_target_name).get_id()]
    num_ees = len(all_ees)
    # construct the boilerplate and function definition
    func_params = ["s_eePos is a pointer to shared memory of size 6*NUM_EE where NUM_EE = " + str(num_ees), \
                   "s_q is the vector of joint positions", \
                   "s_Xhom is the pointer to the homogenous transformation matricies ", \
                   "s_temp is a pointer to helper shared memory of size " + \
                            str(self.gen_end_effector_pose_inner_temp_mem_size(fixed_target_name)), \
                   "s_linalg_smem is optional byte-addressed shared memory for cuBLASDx"]
    func_notes = ["Assumes the Xhom matricies have already been updated for the given q", "Defaults to all leave nodes if fixed_target_name is not provided"]
    func_def_start = "void end_effector_pose_inner" + ("" if fixed_target_name == "" else "_" + fixed_target_name) + "("
    func_def_middle = "T *s_eePos, const T *s_q, const T *s_Xhom, "
    func_def_end = "T *s_temp, unsigned char *s_linalg_smem) {"
    if use_thread_group:
        func_def_start += "cgrps::thread_group tgrp, "
        func_params.insert(0,"tgrp is the handle to the thread_group running this function")
    func_def_middle, func_params = self.gen_insert_helpers_func_def_params(func_def_middle, func_params, -1, NO_XI_FLAG = True)
    func_def = func_def_start + func_def_middle + func_def_end
    # now generate the code
    self.gen_add_func_doc("Computes the End Effector Position",\
                          func_notes,func_params,None)
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
        self.gen_add_code_line("for (int i = 0; i < " + str(n) + "; i++){printf(\"X[%d]\\n\",i); printMat<T,4,4>(&s_Xhom[16*i],4);}")
        self.gen_add_end_control_flow()
        self.gen_add_sync(use_thread_group)

    #
    # For each chain we need to (in parallel) multiply the Xmats
    # 
    self.gen_add_code_line("//")
    self.gen_add_code_line("// For each branch in parallel chain up the transform")
    self.gen_add_code_line("// Keep chaining until reaching the root (starting from the leaves)")
    self.gen_add_code_line("//")
    parent = -1
    for bfs_level in range(n_bfs_levels + (0 if fixed_target_name == "" else 1)): # at most bfs levels of parents to chain (unless with fixed target can be one larger)
        # if serial chain manipulator then this is easy
        if self.robot.is_serial_chain():
            self.gen_add_code_line("// Serial chain manipulator so optimize as parent is jid-1")
            if bfs_level == 0:
                self.gen_add_code_line("// First set to leaf (or fixed) transform")
                self.gen_add_parallel_loop("ind",str(16),use_thread_group)
                self.gen_add_code_line("s_temp[ind] = s_Xhom[16*" + str(all_ees[0]) + " + ind];")
                self.gen_add_end_control_flow()
                self.gen_add_sync(use_thread_group)
                if fixed_target_name == "":
                    parent = self.robot.get_parent_id(all_ees[0])
                else:
                    parent_name = self.robot.get_fixed_joint_by_id(all_ees[0]).get_parent()
                    parent = self.robot.get_joint_by_name(parent_name).get_id() if parent_name != "" else -1
            else:
                if parent == -1:
                    break # if no parent then we are done (this can happen if we have a fixed joint that is not at the end of the chain)
                self.gen_add_code_line("// Update with parent transform until you reach the base [level " + str(bfs_level) + "/" + str(n_bfs_levels-1) + "]")
                even = bfs_level % 2
                tempDstOffset = 16*(even)
                tempSrcOffset = 16*(not even)
                self.gen_add_parallel_loop("ind",str(16),use_thread_group)
                self.gen_add_code_line("int row = ind % 4; int col = ind / 4;")
                self.gen_add_code_line("s_temp[ind + " + str(tempDstOffset) + "] = dot_prod<T,4,4,1>" + \
                                       "(&s_Xhom[16*" + str(parent) + " + row], &s_temp[" + str(tempSrcOffset) + " + 4*col]);")
                self.gen_add_end_control_flow()
                self.gen_add_sync(use_thread_group)
                # update parent for next loop (if there is one)
                parent = self.robot.get_parent_id(parent)
        else:
            # if first loop then just set to transform at the leaf
            if bfs_level == 0:
                self.gen_add_code_line("// First set to leaf transform")
                self.gen_add_parallel_loop("ind",str(16*num_ees),use_thread_group)
                self.gen_add_code_line("int rc = ind % 16;")
                select_var_vals = [("int", "eeInd", [str(jid) for jid in all_ees])]
                self.gen_add_multi_threaded_select("ind", "<", [str(16*(i+1)) for i in range(num_ees)], select_var_vals)
                self.gen_add_code_line("s_temp[ind] = s_Xhom[16*eeInd + rc];")
                self.gen_add_end_control_flow()
                self.gen_add_sync(use_thread_group)
            else:
                self.gen_add_code_line("// Update with parent transform until you reach the base [level " + str(bfs_level) + "/" + str(n_bfs_levels-1) + "]")
                # get the parents we need at this level working backwards from all_ees
                curr_parents = all_ees
                for i in range(bfs_level):
                    curr_parents = [(-1 if jid == -1 else self.robot.get_parent_id(jid)) for jid in curr_parents]
                # need to swap dst and start each time
                even = bfs_level % 2
                tempDstOffset = 16*num_ees*(even)
                tempSrcOffset = 16*num_ees*(not even)
                self.gen_add_parallel_loop("ind",str(16*num_ees),use_thread_group)
                self.gen_add_code_line("int row = ind % 4; int col = (ind / 4) % 4; int eeOffset = ind - (ind % 16);")
                # get parents for this level
                select_var_vals = [("int", "parent_jid", [str(jid) for jid in curr_parents])]
                self.gen_add_multi_threaded_select("ind", "<", [str(16*(i+1)) for i in range(num_ees)], select_var_vals)
                if (-1 in curr_parents):
                    self.gen_add_code_line("if(parent_jid == -1){continue;}")
                self.gen_add_code_line("s_temp[ind + " + str(tempDstOffset) + "] = dot_prod<T,4,4,1>" + \
                                       "(&s_Xhom[16*parent_jid + row], &s_temp[" + str(tempSrcOffset) + " + eeOffset + 4*col]);")
                self.gen_add_end_control_flow()
                self.gen_add_sync(use_thread_group)
    
    self.gen_add_code_line("//")
    self.gen_add_code_line("// Now extract the eePos from the Tansforms")
    self.gen_add_code_line("// TODO: ADD OFFSETS")
    self.gen_add_code_line("//")
    tempOffset = 16*num_ees*(bfs_level % 2)
    # xyz position is easy (eePos_xyz1 = Xmat_hom * offset) where offset = [x,y,z,1]
    self.gen_add_parallel_loop("ind",str(3*num_ees),use_thread_group)
    self.gen_add_code_line("// xyz is easy")
    self.gen_add_code_line("int xyzInd = ind % 3; int eeInd = ind / 3; T *s_Xmat_hom = &s_temp[" + str(tempOffset) + " + 16*eeInd];")
    self.gen_add_code_line("s_eePos[6*eeInd + xyzInd] = s_Xmat_hom[12 + xyzInd];")
    # roll pitch yaw is a bit more difficult
    self.gen_add_code_line("// roll pitch yaw is a bit more difficult")
    self.gen_add_code_line("if(xyzInd > 0){continue;}")
    self.gen_add_code_line("s_eePos[6*eeInd + 3] = atan2(s_Xmat_hom[6],s_Xmat_hom[10]);")
    self.gen_add_code_line("s_eePos[6*eeInd + 4] = -atan2(s_Xmat_hom[2],sqrt(s_Xmat_hom[6]*s_Xmat_hom[6] + s_Xmat_hom[10]*s_Xmat_hom[10]));")
    self.gen_add_code_line("s_eePos[6*eeInd + 5] = atan2(s_Xmat_hom[1],s_Xmat_hom[0]);")
    self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)
    self.gen_add_end_function()

def gen_end_effector_pose_device_temp_mem_size(self, fixed_target_name = ""):
    n = self.robot.get_num_pos()
    XHom_size, dXhom_size, d2Xhom_size = self.gen_get_Xhom_size()
    wrapper_size = self.gen_topology_helpers_size() + XHom_size # for Xhom
    return self.gen_end_effector_pose_inner_temp_mem_size(fixed_target_name) + wrapper_size

def gen_end_effector_pose_device(self, use_thread_group = False, fixed_target_name = ""):
    n = self.robot.get_num_pos()
    num_ees = self.robot.get_total_leaf_nodes() if fixed_target_name == "" else 1
    # construct the boilerplate and function definition
    func_params = ["s_eePos is a pointer to shared memory of size 6*NUM_EE where NUM_EE = " + str(num_ees), \
                   "s_q is the vector of joint positions", \
                   "d_robotModel is the pointer to the initialized model specific helpers on the GPU (XImats, topology_helpers, etc.)"]
    func_notes = []
    func_def_start = "void end_effector_pose_device" + ("" if fixed_target_name == "" else "_" + fixed_target_name) + "("
    func_def_middle = "T *s_eePos, const T *s_q, "
    func_def_end = "const robotModel<T> *d_robotModel) {"
    if use_thread_group:
        func_def_start += "cgrps::thread_group tgrp, "
        func_params.insert(0,"tgrp is the handle to the thread_group running this function")
    func_def = func_def_start + func_def_middle + func_def_end
    # then generate the code
    self.gen_add_func_doc("Computes the End Effector Position",\
                          func_notes,func_params,None)
    self.gen_add_code_line("template <typename T>")
    self.gen_add_code_line("__device__")
    self.gen_add_code_line(func_def, True)
    # add the shared memory variables
    shared_mem_size = self.gen_end_effector_pose_inner_temp_mem_size(fixed_target_name)
    self.gen_XmatsHom_helpers_temp_shared_memory_code(shared_mem_size, include_linalg_scratch = True,
                                                      linalg_scratch_bytes = "GRID_EE_LINALG_SHARED_BYTES<T>()")
    # then load/update XI and run the algo
    self.gen_load_update_XmatsHom_helpers_function_call(use_thread_group)
    self.gen_end_effector_pose_inner_function_call(use_thread_group, fixed_target_name = fixed_target_name)
    self.gen_add_end_function()

def gen_end_effector_pose_kernel(self, use_thread_group = False, single_call_timing = False, fixed_target_name = ""):
    n = self.robot.get_num_pos()
    num_ees = self.robot.get_total_leaf_nodes() if fixed_target_name == "" else 1
    # define function def and params
    func_params = ["d_eePos is the vector of end effector positions", \
                   "d_q is the vector of joint positions", \
                   "stride_q is the stide between each q", \
                   "d_robotModel is the pointer to the initialized model specific helpers on the GPU (XImats, topology_helpers, etc.)", \
                   "num_timesteps is the length of the trajectory points we need to compute over (or overloaded as test_iters for timing)"]
    func_notes = []
    func_def_start = "void end_effector_pose_kernel" + ("" if fixed_target_name == "" else "_" + fixed_target_name) + "(T *d_eePos, const T *d_q, const int stride_q, "
    func_def_end = "const robotModel<T> *d_robotModel, const int NUM_TIMESTEPS) {"
    func_def = func_def_start + func_def_end
    if single_call_timing:
        func_def = func_def.replace("(", "_single_timing(")
    # then generate the code
    self.gen_add_func_doc("Compute the End Effector Position",\
                          func_notes,func_params,None)
    self.gen_add_code_line("template <typename T, int RESOURCE_TIER = GRID_DEFAULT_RESOURCE_TIER>")
    self.gen_add_code_line("__global__")
    self.gen_add_code_line("__launch_bounds__(tier_max_threads<RESOURCE_TIER>())")
    self.gen_add_code_line(func_def, True)
    # add shared memory variables
    shared_mem_size = self.gen_end_effector_pose_inner_temp_mem_size(fixed_target_name)
    self.gen_XmatsHom_helpers_temp_shared_memory_code(shared_mem_size, extra_t_buffers = [("s_q", n), ("s_eePos", 6*num_ees)],
                                                      include_linalg_scratch = True,
                                                      linalg_scratch_bytes = "GRID_EE_LINALG_SHARED_BYTES<T>()")
    if use_thread_group:
        self.gen_add_code_line("cgrps::thread_group tgrp = TBD;")
    if not single_call_timing:
        # load to shared mem and loop over blocks to compute all requested comps
        self.gen_add_parallel_loop("k","NUM_TIMESTEPS",use_thread_group,block_level = True)
        self.gen_kernel_load_inputs("q","stride_q",str(n),use_thread_group)
        # compute
        self.gen_add_code_line("// compute")
        # then load/update X and run the algo
        self.gen_load_update_XmatsHom_helpers_function_call(use_thread_group)
        self.gen_end_effector_pose_inner_function_call(use_thread_group, fixed_target_name = fixed_target_name)
        self.gen_add_sync(use_thread_group)
        # save to global
        self.gen_kernel_save_result("eePos",str(6*num_ees),str(6*num_ees),use_thread_group)
        self.gen_add_end_control_flow()
    else:
        #repurpose NUM_TIMESTEPS for number of timing reps
        self.gen_kernel_load_inputs_single_timing("q",str(n),use_thread_group)
        # then compute in loop for timing
        self.gen_add_code_line("// compute with NUM_TIMESTEPS as NUM_REPS for timing")
        self.gen_add_code_line("for (int rep = 0; rep < NUM_TIMESTEPS; rep++){", True)
        self.gen_anti_licm_input_reload("q",str(n),use_thread_group,feedback_from="eePos")
        # then load/update X and run the algo
        self.gen_load_update_XmatsHom_helpers_function_call(use_thread_group)
        self.gen_end_effector_pose_inner_function_call(use_thread_group, fixed_target_name = fixed_target_name)
        self.gen_anti_licm_output_write("eePos")
        self.gen_add_end_control_flow()
        # save to global
        self.gen_kernel_save_result_single_timing("eePos",str(6*num_ees),use_thread_group)
    self.gen_add_end_function()

def gen_end_effector_pose_host(self, mode = 0, fixed_target_name = ""):
    # default is to do the full kernel call -- options are for single timing or compute only kernel wrapper
    single_call_timing = True if mode == 1 else False
    compute_only = True if mode == 2 else False

    # define function def and params
    func_params = ["hd_data is the packaged input and output pointers", \
                   "d_robotModel is the pointer to the initialized model specific helpers on the GPU (XImats, topology_helpers, etc.)", \
                   "num_timesteps is the length of the trajectory points we need to compute over (or overloaded as test_iters for timing)", \
                   "streams are pointers to CUDA streams for async memory transfers (if needed)"]
    func_notes = []
    func_def_start = "void end_effector_pose" + ("" if fixed_target_name == "" else "_" + fixed_target_name) + \
                            "(gridData<T, KIND> *hd_data, const robotModel<T> *d_robotModel, const int num_timesteps,"
    func_def_end =   "                            const dim3 block_dimms, const dim3 thread_dimms, cudaStream_t *streams) {"
    if single_call_timing:
        func_def_start = func_def_start.replace("(", "_single_timing(")
        func_def_end = "              " + func_def_end
    if compute_only:
        func_def_start = func_def_start.replace("(", "_compute_only(")
        func_def_end = "             " + func_def_end.replace(", cudaStream_t *streams", "")
    # then generate the code
    self.gen_add_func_doc("Compute the End Effector Pose",\
                          func_notes,func_params,None)
    self.gen_add_code_line("template <typename T, bool USE_COMPRESSED_MEM = false, gridDataKind KIND = GRID_DATA_ALL>")
    self.gen_add_code_line("__host__")
    self.gen_add_code_line(func_def_start)
    self.gen_add_code_line(func_def_end, True)
    self.gen_add_code_line("static_assert(KIND == GRID_DATA_ALL || KIND == GRID_DATA_KINEMATICS, \"end_effector_pose requires all-data or kinematics gridData\");")
    func_call_start = "end_effector_pose_kernel" + ("" if fixed_target_name == "" else "_" + fixed_target_name) + \
                        "<T><<<block_dimms,thread_dimms,EE_POS_DYNAMIC_SHARED_MEM_BYTES<T>()>>>(hd_data->d_eePos,hd_data->d_q,stride_q,"
    func_call_end = "d_robotModel,num_timesteps);"
    if single_call_timing:
        func_call_start = func_call_start.replace("kernel<T>","kernel_single_timing<T>")
    if not compute_only:
        # start code with memory transfer
        self.gen_add_code_lines(["// start code with memory transfer", \
                                 "int stride_q;", \
                                 "if (USE_COMPRESSED_MEM) {stride_q = NUM_JOINTS; " + \
                                    "gpuErrchk(cudaMemcpyAsync(hd_data->d_q,hd_data->h_q,stride_q*" + \
                                    ("num_timesteps*" if not single_call_timing else "") + "sizeof(T),cudaMemcpyHostToDevice,streams[0]));}", \
                                 "else {stride_q = 3*NUM_JOINTS; " + \
                                    "gpuErrchk(cudaMemcpyAsync(hd_data->d_q_qd_u,hd_data->h_q_qd_u,stride_q*" + \
                                    ("num_timesteps*" if not single_call_timing else "") + "sizeof(T),cudaMemcpyHostToDevice,streams[0]));}", \
                                 "gpuErrchkKernel();"])
    else:
        self.gen_add_code_line("int stride_q = USE_COMPRESSED_MEM ? NUM_JOINTS: 3*NUM_JOINTS;")
    # then compute but adjust for compressed mem and qdd usage
    self.gen_add_code_line("// then call the kernel")
    func_call = func_call_start + func_call_end
    # add in compressed mem adjusts
    func_call_mem_adjust = "if (USE_COMPRESSED_MEM) {" + func_call + "}"
    func_call_mem_adjust2 = "else                    {" + func_call.replace("hd_data->d_q","hd_data->d_q_qd_u") + "}"
    # compule into a set of code
    func_call_code = [func_call_mem_adjust, func_call_mem_adjust2, "gpuErrchkKernel();"]
    # wrap function call in timing (if needed)
    if single_call_timing:
        func_call_code.insert(0,"struct timespec start, end; clock_gettime(CLOCK_MONOTONIC,&start);")
        func_call_code.append("clock_gettime(CLOCK_MONOTONIC,&end);")
    self.gen_add_code_line("gpuErrchk(grid_check_dynamic_shared_memory_bytes(\"end_effector_pose\", EE_POS_DYNAMIC_SHARED_MEM_BYTES<T>()));")
    self.gen_add_code_lines(func_call_code)
    if not compute_only:
        # then transfer memory back
        self.gen_add_code_lines(["// finally transfer the result back", \
                                 "gpuErrchk(cudaMemcpy(hd_data->h_eePos,hd_data->d_eePos,6*NUM_EES*" + \
                                    ("num_timesteps*" if not single_call_timing else "") + "sizeof(T),cudaMemcpyDeviceToHost));",
                                 "gpuErrchkKernel();"])
    # finally report out timing if requested
    if single_call_timing:
        from ..algo_registry import single_call_printf_line
        self.gen_add_code_line(single_call_printf_line("ee_pose"))
    self.gen_add_end_function()

def gen_end_effector_pose_gradient_inner_temp_mem_size(self, fixed_target_name = ""):
    n = self.robot.get_num_pos()
    num_ees = self.robot.get_total_leaf_nodes() if fixed_target_name == "" else 1
    return 2*2*16*num_ees*n

def gen_end_effector_pose_gradient_inner_function_call(self, use_thread_group = False, updated_var_names = None, fixed_target_name = "",
                                                       temp_in_smem_expr = "true"):
    var_names = dict( \
        s_Xhom_name = "s_XmatsHom", \
        s_dXhom_name = "s_dXmatsHom", \
        s_deePos_name = "s_deePos", \
        s_q_name = "s_q", \
        s_topology_helpers_name = "s_topology_helpers", \
        s_temp_name = "s_temp", \
        d_workspace_name = "nullptr", \
        s_linalg_smem_name = "s_linalg_smem", \
    )
    if updated_var_names is not None:
        for key,value in updated_var_names.items():
            var_names[key] = value
    code_start = "end_effector_pose_gradient_inner" + ("" if fixed_target_name == "" else "_" + fixed_target_name) + "<T, " + temp_in_smem_expr + ">(" + var_names["s_deePos_name"] + ", " + var_names["s_q_name"] + ", "
    code_middle = var_names["s_Xhom_name"] + ", " + var_names["s_dXhom_name"] + ", "
    code_end =  var_names["s_temp_name"] + ", " + var_names["d_workspace_name"] + ", " + var_names["s_linalg_smem_name"] + ");"
    # account for thread group
    if use_thread_group:
        code_start = code_start.replace("(","(tgrp, ")
    # Canonical: append the shared topology-helper arg via the central helper
    # (NO_XI: the ee_pose family takes s_Xhom, not s_XImats). Mirrors the def's
    # gen_insert_helpers_func_def_params(NO_XI_FLAG=True) so def + call can't drift.
    code_middle += self.gen_insert_helpers_function_call(updated_var_names = var_names, NO_XI_FLAG = True)
    self.gen_add_code_line(code_start + code_middle + code_end)

def gen_end_effector_pose_gradient_inner(self, use_thread_group = False, fixed_target_name = ""):
    n = self.robot.get_num_pos()
    n_bfs_levels = self.robot.get_max_bfs_level() + 1 # starts at 0
    if fixed_target_name == "":
        all_ees = self.robot.get_leaf_nodes()
    else:
        all_ees = [self.robot.get_fixed_joint_by_name(fixed_target_name).get_id()]
    num_ees = len(all_ees)
    # construct the boilerplate and function definition
    func_params = ["s_deePos is a pointer to shared memory of size 6*NUM_JOINTS*NUM_EE where NUM_JOINTS = " + str(n) + " and NUM_EE = " + str(num_ees), \
                   "s_q is the vector of joint positions", \
                   "s_Xhom is the pointer to the homogenous transformation matricies ", \
                   "s_dXhom is the pointer to the gradient of the homogenous transformation matricies ", \
                   "s_temp is a pointer to helper shared memory of size " + \
                            str(self.gen_end_effector_pose_gradient_inner_temp_mem_size()), \
                   "s_linalg_smem is optional byte-addressed shared memory for cuBLASDx"]
    func_notes = ["Assumes the Xhom and dXhom matricies have already been updated for the given q"]
    func_def_start = "void end_effector_pose_gradient_inner" + ("" if fixed_target_name == "" else "_" + fixed_target_name) + "("
    func_def_middle = "T *s_deePos, const T *s_q, const T *s_Xhom, const T *s_dXhom, "
    func_def_end = "T *s_temp, T *d_workspace, unsigned char *s_linalg_smem) {"
    if use_thread_group:
        func_def_start += "cgrps::thread_group tgrp, "
        func_params.insert(0,"tgrp is the handle to the thread_group running this function")
    func_def_middle, func_params = self.gen_insert_helpers_func_def_params(func_def_middle, func_params, -1, NO_XI_FLAG = True)
    func_def = func_def_start + func_def_middle + func_def_end
    # now generate the code
    self.gen_add_func_doc("Computes the Gradient of the End Effector Pose with respect to joint position",\
                          func_notes,func_params,None)
    self.gen_add_code_line("template <typename T, bool TEMP_IN_SMEM = true>")
    self.gen_add_code_line("__device__")
    self.gen_add_code_line(func_def, True)
    # Inner-controlled scratch placement: the double-buffered chain workspace
    # moves to d_workspace when !TEMP_IN_SMEM. Reassigning s_temp at the top keeps
    # every s_temp[...] reference below unchanged.
    self.gen_add_code_line("if constexpr (!TEMP_IN_SMEM) { s_temp = d_workspace; } else { (void)d_workspace; }")
    #
    # Initial Debug Prints if Requested
    #
    if self.DEBUG_MODE:
        self.gen_add_sync(use_thread_group)
        self.gen_add_serial_ops(use_thread_group)
        self.gen_add_code_line("printf(\"q\\n\"); printMat<T,1," + str(n) + ">(s_q,1);")
        self.gen_add_code_line("for (int i = 0; i < " + str(n) + "; i++){printf(\"X[%d]\\n\",i); printMat<T,4,4>(&s_Xhom[16*i],4);}")
        self.gen_add_code_line("for (int i = 0; i < " + str(n) + "; i++){printf(\"dX[%d]\\n\",i); printMat<T,4,4>(&s_dXhom[16*i],4);}")
        self.gen_add_end_control_flow()
        self.gen_add_sync(use_thread_group)

    #
    # For each chain we need to (in parallel) multiply the (d)Xmats
    # 
    self.gen_add_code_line("//")
    self.gen_add_code_line("// For each branch/gradient in parallel chain up the transform")
    self.gen_add_code_line("// Keep chaining until reaching the root (starting from the leaves)")
    self.gen_add_code_line("//")
    self.gen_add_code_line("T *s_eeTemp = &s_temp[0]; T *s_deeTemp = &s_temp[" + str(2*16*num_ees*n) + "];")
    parent = -1
    for bfs_level in range(n_bfs_levels + (0 if fixed_target_name == "" else 1)): # at most bfs levels of parents to chain (unless with fixed target can be one larger)
        # if serial chain manipulator then this is easy
        if self.robot.is_serial_chain():
            self.gen_add_code_line("// Serial chain manipulator so optimize as parent is jid-1")
            if bfs_level == 0:
                self.gen_add_code_line("// First set the leaf transforms for eePos and deePos")
                self.gen_add_parallel_loop("ind",str(16*n),use_thread_group)
                self.gen_add_code_line("int djid = ind / 16; int rc = ind % 16; int eeIndStart = 16*" + str(all_ees[0]) + ";")
                self.gen_add_code_line("s_eeTemp[ind] = s_Xhom[eeIndStart + rc];")
                self.gen_add_code_line("const T *s_Xhom_dXhom = grid_xhom_or_dxhom_ptr<T>(s_Xhom, s_dXhom, djid, " + str(all_ees[0]) + ");")
                self.gen_add_code_line("s_deeTemp[ind] = s_Xhom_dXhom[rc];")
                self.gen_add_end_control_flow()
                self.gen_add_sync(use_thread_group)
                # update parent for next loop (if there is one)
                if fixed_target_name == "":
                    parent = self.robot.get_parent_id(all_ees[0])
                else:
                    parent_name = self.robot.get_fixed_joint_by_id(all_ees[0]).get_parent()
                    parent = self.robot.get_joint_by_name(parent_name).get_id() if parent_name != "" else -1
                if self.DEBUG_MODE:
                    self.gen_add_sync(use_thread_group)
                    self.gen_add_serial_ops(use_thread_group)
                    self.gen_add_code_line("for (int i = 0; i < " + str(n) + "; i++){printf(\"X_chain0[%d]\\n\",i); printMat<T,4,4>(&s_eeTemp[16*i],4);}")
                    self.gen_add_code_line("for (int i = 0; i < " + str(n) + "; i++){printf(\"dX_chain0[%d]\\n\",i); printMat<T,4,4>(&s_deeTemp[16*i],4);}")
                    self.gen_add_end_control_flow()
                    self.gen_add_sync(use_thread_group)
            else:
                if parent == -1:
                    break # if no parent then we are done (this can happen if we have a fixed joint that is not at the end of the chain)
                self.gen_add_code_line("// Update with parent transform until you reach the base [level " + str(bfs_level) + "/" + str(n_bfs_levels-1) + "]")
                # need to swap dst and start each time
                even = bfs_level % 2
                tempDstOffset = 16*n*(even)
                tempSrcOffset = 16*n*(not even)
                self.gen_add_parallel_loop("ind",str(16*n),use_thread_group)
                self.gen_add_code_line("int djid = ind / 16; int rc = ind % 16; int row = rc % 4; int colInd = ind - row;")
                self.gen_add_code_line("s_eeTemp[ind + " + str(tempDstOffset) + "] = dot_prod<T,4,4,1>" + \
                                       "(&s_Xhom[16*" + str(parent) + " + row], &s_eeTemp[" + str(tempSrcOffset) + " + colInd]);")
                self.gen_add_code_line("const T *s_Xhom_dXhom = grid_xhom_or_dxhom_ptr<T>(s_Xhom, s_dXhom, djid, " + str(parent) + ");")
                self.gen_add_code_line("s_deeTemp[ind + " + str(tempDstOffset) + "] = dot_prod<T,4,4,1>" + \
                                       "(&s_Xhom_dXhom[row], &s_deeTemp[" + str(tempSrcOffset) + " + colInd]);")
                self.gen_add_end_control_flow()
                self.gen_add_sync(use_thread_group)
                # update parent for next loop (if there is one)
                parent = self.robot.get_parent_id(parent)
                if self.DEBUG_MODE:
                    self.gen_add_sync(use_thread_group)
                    self.gen_add_serial_ops(use_thread_group)
                    self.gen_add_code_line("for (int i = 0; i < " + str(n) + "; i++){printf(\"X_chain0[%d]\\n\",i); printMat<T,4,4>(&s_eeTemp[16*i + " + str(tempDstOffset) + "],4);}")
                    self.gen_add_code_line("for (int i = 0; i < " + str(n) + "; i++){printf(\"dX_chain0[%d]\\n\",i); printMat<T,4,4>(&s_deeTemp[16*i + " + str(tempDstOffset) + "],4);}")
                    self.gen_add_end_control_flow()
                    self.gen_add_sync(use_thread_group)
        else:
            # if first loop then just set to transform at the leaf
            if bfs_level == 0:
                self.gen_add_code_line("// First set the leaf transforms for eePos and deePos")
                self.gen_add_parallel_loop("ind",str(16*n*num_ees),use_thread_group)
                self.gen_add_code_line("int rc = ind % 16; int djid = (ind / 16) % " + str(n) + ";")
                select_var_vals = [("int", "eeInd", [str(jid) for jid in all_ees])]
                # make sure to zero out all things not in the chain
                jidChainCode = []
                for eejid in all_ees:
                    jidChain = sorted(self.robot.get_ancestors_by_id(eejid))
                    jidChain.append(eejid)
                    qinds = []
                    for jid in jidChain:
                        qind = self.robot.get_joint_index_q(jid)
                        qinds.extend(qind if isinstance(qind, list) else [qind])
                    code = self.gen_var_in_list("djid", [str(qind) for qind in qinds])
                    jidChainCode.append(code)
                select_var_vals.append(("bool", "inChain", jidChainCode))
                self.gen_add_multi_threaded_select("ind", "<", [str(16*n*(i+1)) for i in range(num_ees)], select_var_vals)
                self.gen_add_code_line("s_eeTemp[ind] = s_Xhom[16*eeInd + rc];")
                self.gen_add_code_line("const T *s_Xhom_dXhom = grid_xhom_or_dxhom_ptr<T>(s_Xhom, s_dXhom, djid, eeInd);")
                self.gen_add_code_line("s_deeTemp[ind] = inChain * s_Xhom_dXhom[rc];")
                self.gen_add_end_control_flow()
                self.gen_add_sync(use_thread_group)
                if self.DEBUG_MODE:
                    self.gen_add_sync(use_thread_group)
                    self.gen_add_serial_ops(use_thread_group)
                    self.gen_add_code_line("for (int i = 0; i < " + str(n*num_ees) + "; i++){printf(\"X_chain level[%d] with dj_ee_id [%d]\\n\"," + str(bfs_level) + ",i); printMat<T,4,4>(&s_eeTemp[16*i],4);}")
                    self.gen_add_code_line("for (int i = 0; i < " + str(n*num_ees) + "; i++){printf(\"dX_chain level[%d] with dj_ee_id [%d]\\n\"," + str(bfs_level) + ",i); printMat<T,4,4>(&s_deeTemp[16*i],4);}")
                    self.gen_add_end_control_flow()
                    self.gen_add_sync(use_thread_group)
            else:
                self.gen_add_code_line("// Update with parent transform until you reach the base [level " + str(bfs_level) + "/" + str(n_bfs_levels-1) + "]")
                # get the parents we need at this level working backwards from all_ees
                curr_parents = all_ees
                for i in range(bfs_level):
                    curr_parents = [(-1 if jid == -1 else self.robot.get_parent_id(jid)) for jid in curr_parents]
                # need to swap dst and start each time
                even = bfs_level % 2
                tempDstOffset = 16*n*num_ees*(even)
                tempSrcOffset = 16*n*num_ees*(not even)
                # SIMT parallel_loop computes s_eeTemp (the 4×4 chain product
                # against parent_jid's homogeneous transform) and s_deeTemp
                # (the gradient version) per (djid, ee, row, col) cell.
                # Pre-v2.0 had an opt-in `glass::nvidia::gemm_strided_batched_1d`
                # path for s_eeTemp on cuBLASDx builds; the 2026-05-18
                # autotune + sweep showed cuBLASDx loses to SIMT on 4×4×4 by
                # 2.6× so the dispatch was removed (see
                # docs/source/user_guide/concepts/cublasdx_removal_design.rst).
                self.gen_add_parallel_loop("ind",str(16*n*num_ees),use_thread_group)
                self.gen_add_code_line("int rc = ind % 16; int djid = (ind / 16) % " + str(n) + ";")
                self.gen_add_code_line("int row = rc % 4; int colInd = ind - row;")
                # get parents for this level
                select_var_vals = [("int", "parent_jid", [str(jid) for jid in curr_parents])]
                self.gen_add_multi_threaded_select("ind", "<", [str(16*n*(i+1)) for i in range(num_ees)], select_var_vals)
                if (-1 in curr_parents):
                    self.gen_add_code_line("if(parent_jid == -1){continue;}")
                self.gen_add_code_line("s_eeTemp[ind + " + str(tempDstOffset) + "] = dot_prod<T,4,4,1>" + \
                                       "(&s_Xhom[16*parent_jid + row], &s_eeTemp[" + str(tempSrcOffset) + " + colInd]);")
                self.gen_add_code_line("const T *s_Xhom_dXhom = grid_xhom_or_dxhom_ptr<T>(s_Xhom, s_dXhom, djid, parent_jid);")
                self.gen_add_code_line("s_deeTemp[ind + " + str(tempDstOffset) + "] = dot_prod<T,4,4,1>" + \
                                       "(&s_Xhom_dXhom[row], &s_deeTemp[" + str(tempSrcOffset) + " + colInd]);")
                self.gen_add_end_control_flow()
                self.gen_add_sync(use_thread_group)
                if self.DEBUG_MODE:
                    self.gen_add_sync(use_thread_group)
                    self.gen_add_serial_ops(use_thread_group)
                    self.gen_add_code_line("for (int i = 0; i < " + str(n*num_ees) + "; i++){printf(\"X_chain[%d]\\n\",i); printMat<T,4,4>(&s_eeTemp[16*i + " + str(tempDstOffset) + "],4);}")
                    self.gen_add_code_line("for (int i = 0; i < " + str(n*num_ees) + "; i++){printf(\"dX_chain[%d]\\n\",i); printMat<T,4,4>(&s_deeTemp[16*i + " + str(tempDstOffset) + "],4);}")
                    self.gen_add_end_control_flow()
                    self.gen_add_sync(use_thread_group)
    
    self.gen_add_code_line("//")
    self.gen_add_code_line("// Now extract the eePos from the Transforms")
    self.gen_add_code_line("// TODO: ADD OFFSETS")
    self.gen_add_code_line("//")
    tempOffset = 16*n*num_ees*(bfs_level % 2)
    self.gen_add_parallel_loop("ind",str(6*n*num_ees),use_thread_group)
    self.gen_add_code_line("int outputInd = ind % 6; int deeInd = ind / 6;")
    self.gen_add_code_line("T *s_Xmat_hom = &s_eeTemp[" + str(tempOffset) + " + 16*deeInd]; T *s_dXmat_hom = &s_deeTemp[" + str(tempOffset) + " + 16*deeInd];")
    # xyz position is easy (eePos_xyz1 = Xmat_hom * offset) where offset = [x,y,z,1]
    self.gen_add_code_line("// xyz is easy")
    self.gen_add_code_line("if (outputInd < 3){s_deePos[6*deeInd + outputInd] = s_dXmat_hom[12 + outputInd];}")
    # roll pitch yaw is a bit more difficult
    self.gen_add_code_line("// roll pitch yaw is a bit more difficult")
    self.gen_add_code_line("// note: d/dz of arctan2(y(z),x(z)) = [-x'(z)y(z)+x(z)y'(z)]/[(x(z)^2 + y(z)^2)]")
    self.gen_add_code_line("// Also note that d/dz of sqrt(f(z)) = f'(z)/2sqrt(f(z))")
    self.gen_add_code_line("else {", add_indent_after=True)
    self.gen_add_code_line("// simpler to recompute")
    self.gen_add_code_line("T sqrtTerm = sqrt(s_Xmat_hom[10]*s_Xmat_hom[10] + s_Xmat_hom[6]*s_Xmat_hom[6]);") 
    self.gen_add_code_line("T dsqrtTerm = (s_Xmat_hom[10]*s_dXmat_hom[10] + s_Xmat_hom[6]*s_dXmat_hom[6])/sqrtTerm;")
    select_var_vals = [("T", "y",       ["s_Xmat_hom[6]",  "-s_Xmat_hom[2]",  "s_Xmat_hom[1]"]), \
                       ("T", "x",       ["s_Xmat_hom[10]",  "sqrtTerm",        "s_Xmat_hom[0]"]), \
                       ("T", "y_prime", ["s_dXmat_hom[6]", "-s_dXmat_hom[2]", "s_dXmat_hom[1]"]), \
                       ("T", "x_prime", ["s_dXmat_hom[10]", "dsqrtTerm",       "s_dXmat_hom[0]"])]
    self.gen_add_multi_threaded_select("outputInd", "==", [str(i) for i in range(3,6)], select_var_vals)
    self.gen_add_code_line("s_deePos[6*deeInd + outputInd] = (-x_prime*y + x*y_prime)/(x*x + y*y);")
    self.gen_add_end_control_flow()
    self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)
    self.gen_add_end_function()

def gen_end_effector_pose_gradient_device_temp_mem_size(self, fixed_target_name = ""):
    n = self.robot.get_num_pos()
    XHom_size, dXhom_size, d2Xhom_size = self.gen_get_Xhom_size()
    wrapper_size = self.gen_topology_helpers_size() + XHom_size + dXhom_size # for Xhom and dXhom
    return self.gen_end_effector_pose_gradient_inner_temp_mem_size(fixed_target_name) + wrapper_size

def gen_end_effector_pose_gradient_device(self, use_thread_group = False, fixed_target_name = ""):
    n = self.robot.get_num_pos()
    num_ees = self.robot.get_total_leaf_nodes() if fixed_target_name == "" else 1
    # construct the boilerplate and function definition
    func_params = ["s_deePos is a pointer to shared memory of size 6*NUM_JOINTS*NUM_EE where NUM_JOINTS = " + str(n) + " and NUM_EE = " + str(num_ees), \
                   "s_q is the vector of joint positions", \
                   "d_robotModel is the pointer to the initialized model specific helpers on the GPU (XImats, topology_helpers, etc.)"]
    func_notes = []
    func_def_start = "void end_effector_pose_gradient_device" + ("" if fixed_target_name == "" else "_" + fixed_target_name) + "("
    func_def_middle = "T *s_deePos, const T *s_q, "
    func_def_end = "const robotModel<T> *d_robotModel) {"
    if use_thread_group:
        func_def_start += "cgrps::thread_group tgrp, "
        func_params.insert(0,"tgrp is the handle to the thread_group running this function")
    func_def = func_def_start + func_def_middle + func_def_end
    # then generate the code
    self.gen_add_func_doc("Computes the Gradient of the End Effector Pose with respect to joint position",\
                          func_notes,func_params,None)
    self.gen_add_code_line("template <typename T>")
    self.gen_add_code_line("__device__")
    self.gen_add_code_line(func_def, True)
    # add the shared memory variables
    shared_mem_size = self.gen_end_effector_pose_gradient_inner_temp_mem_size(fixed_target_name)
    self.gen_XmatsHom_helpers_temp_shared_memory_code(shared_mem_size, include_gradients = True, include_linalg_scratch = True,
                                                      linalg_scratch_bytes = "GRID_EE_LINALG_SHARED_BYTES<T>()")
    # then load/update XI and run the algo
    self.gen_load_update_XmatsHom_helpers_function_call(use_thread_group, include_gradients = True)
    self.gen_end_effector_pose_gradient_inner_function_call(use_thread_group, fixed_target_name = fixed_target_name)
    self.gen_add_end_function()

_EE_GRAD_PICK_FLAGS = [
    # (use_workspace_temp, use_workspace_dxhom)
    (False, False),   # pick 0: full smem (PERF)
    (True,  False),   # pick 1: inner_temp + s_deePos -> workspace/global (LITE)
    (True,  True),    # pick 2: also dXmatsHom -> workspace (MINIMAL)
]

def _emit_eepose_grad_kernel_body_for_flags(self, n, num_ees, fixed_target_name,
                                            use_workspace_temp, use_workspace_dxhom,
                                            single_call_timing, use_thread_group):
    """Emit the EE_POSE_GRAD kernel body specialized for one tier's spill flags.
    Wrapped in a brace pair (caller emits the `if constexpr (...)` head).
    Used by gen_end_effector_pose_gradient_kernel to emit either a single body
    (collapsed picks) or three branched bodies (divergent picks). Mirrors
    _emit_d2ee_kernel_body_for_flags."""
    shared_mem_size = 0 if use_workspace_temp else self.gen_end_effector_pose_gradient_inner_temp_mem_size(fixed_target_name)
    extra_t_buffers = [("s_q", n)] if use_workspace_temp else [("s_q", n), ("s_deePos", 6*n*num_ees)]
    self.gen_XmatsHom_helpers_temp_shared_memory_code(shared_mem_size, include_gradients = True,
                                                      extra_t_buffers = extra_t_buffers,
                                                      include_dxhom_shared = not use_workspace_dxhom,
                                                      include_linalg_scratch = True,
                                                      linalg_scratch_bytes = "GRID_EE_LINALG_SHARED_BYTES<T>()")
    if not use_workspace_temp:
        self.gen_add_code_line("(void)d_workspace;")
    if use_thread_group:
        self.gen_add_code_line("cgrps::thread_group tgrp = TBD;")
    if not single_call_timing:
        self.gen_add_parallel_loop("k","NUM_TIMESTEPS",use_thread_group,block_level = True)
        self.gen_kernel_load_inputs("q","stride_q",str(n),use_thread_group)
        if use_workspace_dxhom:
            self.gen_add_code_line("T *s_dXmatsHom = reinterpret_cast<T *>(&d_workspace[k*GRID_WORKSPACE_BYTES_PER_TIMESTEP<T>() + GRID_EE_GRAD_WORKSPACE_DXHOM_OFFSET_BYTES<T>()]);")
        if use_workspace_temp:
            self.gen_add_code_line("T *s_deePos = &d_deePos[k*" + str(6*n*num_ees) + "];")
            self.gen_add_code_line("T *s_eegrad_temp = reinterpret_cast<T *>(&d_workspace[k*GRID_WORKSPACE_BYTES_PER_TIMESTEP<T>() + GRID_EE_GRAD_WORKSPACE_TEMP_OFFSET_BYTES<T>()]);")
        self.gen_add_code_line("// compute")
        self.gen_load_update_XmatsHom_helpers_function_call(use_thread_group, include_gradients = True)
        # Inner-controlled: pass both arenas + placement; the inner picks where
        # the chain workspace lives via TEMP_IN_SMEM.
        updated = {"d_workspace_name": "s_eegrad_temp"} if use_workspace_temp else None
        self.gen_end_effector_pose_gradient_inner_function_call(use_thread_group, fixed_target_name = fixed_target_name,
            updated_var_names = updated, temp_in_smem_expr = ("false" if use_workspace_temp else "true"))
        self.gen_add_sync(use_thread_group)
        if not use_workspace_temp:
            self.gen_kernel_save_result("deePos",str(6*n*num_ees),str(6*n*num_ees),use_thread_group)
        self.gen_add_end_control_flow()
    else:
        self.gen_kernel_load_inputs_single_timing("q",str(n),use_thread_group)
        if use_workspace_dxhom:
            self.gen_add_code_line("T *s_dXmatsHom = reinterpret_cast<T *>(&d_workspace[GRID_EE_GRAD_WORKSPACE_DXHOM_OFFSET_BYTES<T>()]);")
        if use_workspace_temp:
            self.gen_add_code_line("T *s_deePos = d_deePos;")
            self.gen_add_code_line("T *s_eegrad_temp = reinterpret_cast<T *>(&d_workspace[GRID_EE_GRAD_WORKSPACE_TEMP_OFFSET_BYTES<T>()]);")
        self.gen_add_code_line("// compute with NUM_TIMESTEPS as NUM_REPS for timing")
        self.gen_add_code_line("for (int rep = 0; rep < NUM_TIMESTEPS; rep++){", True)
        # TODO(licm-eepose-grad): sm_86-specific, deprioritized. See pre-Phase-3d note in git history.
        self.gen_anti_licm_input_reload("q",str(n),use_thread_group,feedback_from="deePos")
        self.gen_load_update_XmatsHom_helpers_function_call(use_thread_group, include_gradients = True)
        updated = {"d_workspace_name": "s_eegrad_temp"} if use_workspace_temp else None
        self.gen_end_effector_pose_gradient_inner_function_call(use_thread_group, fixed_target_name = fixed_target_name,
            updated_var_names = updated, temp_in_smem_expr = ("false" if use_workspace_temp else "true"))
        self.gen_anti_licm_output_write("deePos")
        self.gen_add_end_control_flow()
        if not use_workspace_temp:
            self.gen_kernel_save_result_single_timing("deePos",str(6*n*num_ees),use_thread_group)


def gen_end_effector_pose_gradient_kernel(self, use_thread_group = False, single_call_timing = False, fixed_target_name = ""):
    n = self.robot.get_num_pos()
    num_ees = self.robot.get_total_leaf_nodes() if fixed_target_name == "" else 1
    func_params = ["d_deePos is the vector of end effector positions gradients", \
                   "d_workspace is the generated global spill workspace", \
                   "d_q is the vector of joint positions", \
                   "stride_q is the stide between each q", \
                   "d_robotModel is the pointer to the initialized model specific helpers on the GPU (XImats, topology_helpers, etc.)", \
                   "num_timesteps is the length of the trajectory points we need to compute over (or overloaded as test_iters for timing)"]
    func_notes = []
    func_def_start = "void end_effector_pose_gradient_kernel" + ("" if fixed_target_name == "" else "_" + fixed_target_name) + "(T *d_deePos, unsigned char *d_workspace, const T *d_q, const int stride_q, "
    func_def_end = "const robotModel<T> *d_robotModel, const int NUM_TIMESTEPS) {"
    func_def = func_def_start + func_def_end
    if single_call_timing:
        func_def = func_def.replace("(", "_single_timing(")
    self.gen_add_func_doc("Computes the Gradient of the End Effector Pose with respect to joint position",\
                          func_notes,func_params,None)
    self.gen_add_code_line("template <typename T, int RESOURCE_TIER = GRID_DEFAULT_RESOURCE_TIER>")
    self.gen_add_code_line("__global__")
    self.gen_add_code_line("__launch_bounds__(tier_max_threads<RESOURCE_TIER>())")
    self.gen_add_code_line(func_def, True)
    # Tier dispatch: when the 3 picks collapse, emit one body. When they
    # diverge, emit three if-constexpr branches — each specialized for that
    # tier's spill flags. DEE_POS_DYNAMIC_SHARED_MEM_BYTES<T, TIER>() is
    # tier-aware.
    picks = getattr(self, "ee_grad_spill_tier_3way", (0, 0, 0))
    if picks[0] == picks[1] == picks[2]:
        uwt, uwd = _EE_GRAD_PICK_FLAGS[picks[0]]
        _emit_eepose_grad_kernel_body_for_flags(self, n, num_ees, fixed_target_name, uwt, uwd, single_call_timing, use_thread_group)
    else:
        tier_names = ("TIER_PERF", "TIER_LITE", "TIER_MINIMAL")
        for tier_idx, (tier_name, pick) in enumerate(zip(tier_names, picks)):
            uwt, uwd = _EE_GRAD_PICK_FLAGS[pick]
            head = "if constexpr (RESOURCE_TIER == " + tier_name + ") {" if tier_idx == 0 else \
                   "else if constexpr (RESOURCE_TIER == " + tier_name + ") {"
            self.gen_add_code_line(head, True)
            _emit_eepose_grad_kernel_body_for_flags(self, n, num_ees, fixed_target_name, uwt, uwd, single_call_timing, use_thread_group)
            self.gen_add_end_control_flow()
    self.gen_add_end_function()

def gen_end_effector_pose_gradient_host(self, mode = 0, fixed_target_name = ""):
    # default is to do the full kernel call -- options are for single timing or compute only kernel wrapper
    single_call_timing = True if mode == 1 else False
    compute_only = True if mode == 2 else False

    # define function def and params
    func_params = ["hd_data is the packaged input and output pointers", \
                   "d_robotModel is the pointer to the initialized model specific helpers on the GPU (XImats, topology_helpers, etc.)", \
                   "num_timesteps is the length of the trajectory points we need to compute over (or overloaded as test_iters for timing)", \
                   "streams are pointers to CUDA streams for async memory transfers (if needed)"]
    func_notes = []
    func_def_start = "void end_effector_pose_gradient" + ("" if fixed_target_name == "" else "_" + fixed_target_name) + "(gridData<T, KIND> *hd_data, const robotModel<T> *d_robotModel, const int num_timesteps,"
    func_def_end =   "                            const dim3 block_dimms, const dim3 thread_dimms, cudaStream_t *streams) {"
    if single_call_timing:
        func_def_start = func_def_start.replace("(", "_single_timing(")
        func_def_end = "              " + func_def_end
    if compute_only:
        func_def_start = func_def_start.replace("(", "_compute_only(")
        func_def_end = "             " + func_def_end.replace(", cudaStream_t *streams", "")
    # then generate the code
    self.gen_add_func_doc("Computes the Gradient of the End Effector Pose with respect to joint position",\
                          func_notes,func_params,None)
    self.gen_add_code_line("template <typename T, bool USE_COMPRESSED_MEM = false, gridDataKind KIND = GRID_DATA_ALL>")
    self.gen_add_code_line("__host__")
    self.gen_add_code_line(func_def_start)
    self.gen_add_code_line(func_def_end, True)
    self.gen_add_code_line("static_assert(KIND == GRID_DATA_ALL || KIND == GRID_DATA_KINEMATICS, \"end_effector_pose_gradient requires all-data or kinematics gridData\");")
    func_call_start = "end_effector_pose_gradient_kernel" + ("" if fixed_target_name == "" else "_" + fixed_target_name) + "<T><<<block_dimms,thread_dimms,DEE_POS_DYNAMIC_SHARED_MEM_BYTES<T>()>>>(hd_data->d_deePos,hd_data->d_workspace,hd_data->d_q,stride_q,"
    func_call_end = "d_robotModel,num_timesteps);"
    if single_call_timing:
        func_call_start = func_call_start.replace("kernel<T>","kernel_single_timing<T>")
    if not compute_only:
        # start code with memory transfer
        self.gen_add_code_lines(["// start code with memory transfer", \
                                 "int stride_q;", \
                                 "if (USE_COMPRESSED_MEM) {stride_q = NUM_JOINTS; " + \
                                    "gpuErrchk(cudaMemcpyAsync(hd_data->d_q,hd_data->h_q,stride_q*" + \
                                    ("num_timesteps*" if not single_call_timing else "") + "sizeof(T),cudaMemcpyHostToDevice,streams[0]));}", \
                                 "else {stride_q = 3*NUM_JOINTS; " + \
                                    "gpuErrchk(cudaMemcpyAsync(hd_data->d_q_qd_u,hd_data->h_q_qd_u,stride_q*" + \
                                    ("num_timesteps*" if not single_call_timing else "") + "sizeof(T),cudaMemcpyHostToDevice,streams[0]));}", \
                                 "gpuErrchkKernel();"])
    else:
        self.gen_add_code_line("int stride_q = USE_COMPRESSED_MEM ? NUM_JOINTS: 3*NUM_JOINTS;")
    # then compute but adjust for compressed mem and qdd usage
    self.gen_add_code_line("// then call the kernel")
    func_call = func_call_start + func_call_end
    # add in compressed mem adjusts
    func_call_mem_adjust = "if (USE_COMPRESSED_MEM) {" + func_call + "}"
    func_call_mem_adjust2 = "else                    {" + func_call.replace("hd_data->d_q","hd_data->d_q_qd_u") + "}"
    # compule into a set of code
    func_call_code = [func_call_mem_adjust, func_call_mem_adjust2, "gpuErrchkKernel();"]
    # wrap function call in timing (if needed)
    if single_call_timing:
        func_call_code.insert(0,"struct timespec start, end; clock_gettime(CLOCK_MONOTONIC,&start);")
        func_call_code.append("clock_gettime(CLOCK_MONOTONIC,&end);")
    self.gen_add_code_line("gpuErrchk(grid_check_dynamic_shared_memory_bytes(\"end_effector_pose_gradient\", DEE_POS_DYNAMIC_SHARED_MEM_BYTES<T>()));")
    workspace_bytes = "GRID_WORKSPACE_BYTES_PER_TIMESTEP<T>()" if single_call_timing else "GRID_WORKSPACE_BYTES_PER_TIMESTEP<T>()*static_cast<size_t>(num_timesteps)"
    self.gen_add_code_line("if (GRID_EE_GRAD_USES_WORKSPACE_TEMP) {gpuErrchk(grid_begin_l2_persisting(0, hd_data->d_workspace, " + workspace_bytes + "));}")
    self.gen_add_code_lines(func_call_code)
    self.gen_add_code_line("if (GRID_EE_GRAD_USES_WORKSPACE_TEMP) {gpuErrchk(grid_end_l2_persisting(0));}")
    if not compute_only:
        # then transfer memory back
        self.gen_add_code_lines(["// finally transfer the result back", \
                                 "gpuErrchk(cudaMemcpy(hd_data->h_deePos,hd_data->d_deePos,6*NUM_EES*NUM_JOINTS*" + \
                                    ("num_timesteps*" if not single_call_timing else "") + "sizeof(T),cudaMemcpyDeviceToHost));",
                                 "gpuErrchkKernel();"])
    # finally report out timing if requested
    if single_call_timing:
        from ..algo_registry import single_call_printf_line
        self.gen_add_code_line(single_call_printf_line("ee_pose_gradient"))
    self.gen_add_end_function()

def gen_end_effector_pose_gradient_hessian_d2_temp_mem_size(self):
    n = self.robot.get_num_pos()
    num_ees = self.robot.get_total_leaf_nodes()
    return 2*16*num_ees*n*n

def gen_end_effector_pose_gradient_hessian_inner_temp_mem_size(self, include_d2_temp = True):
    n = self.robot.get_num_pos()
    num_ees = self.robot.get_total_leaf_nodes()
    d2_temp_size = self.gen_end_effector_pose_gradient_hessian_d2_temp_mem_size() if include_d2_temp else 0
    return 2*16*num_ees*(n+1) + d2_temp_size

def gen_end_effector_pose_gradient_hessian_inner_function_call(self, use_thread_group = False, updated_var_names = None):
    var_names = dict( \
        s_Xhom_name = "s_XmatsHom", \
        s_dXhom_name = "s_dXmatsHom", \
        s_d2Xhom_name = "s_d2XmatsHom", \
        s_deePos_name = "s_deePos", \
        s_d2eePos_name = "s_d2eePos", \
        s_q_name = "s_q", \
        s_topology_helpers_name = "s_topology_helpers", \
        s_temp_name = "s_temp", \
        s_d2eeTemp_name = "s_d2eeTemp", \
        s_linalg_smem_name = "s_linalg_smem", \
    )
    if updated_var_names is not None:
        for key,value in updated_var_names.items():
            var_names[key] = value
    code_start = "end_effector_pose_gradient_hessian_inner<T>(" + var_names["s_d2eePos_name"] + ", " + var_names["s_deePos_name"] + ", " + var_names["s_q_name"] + ", "
    code_middle = var_names["s_Xhom_name"] + ", " + var_names["s_dXhom_name"] + ", " + var_names["s_d2Xhom_name"] + ", "
    code_end =  var_names["s_temp_name"] + ", " + var_names["s_d2eeTemp_name"] + ", " + var_names["s_linalg_smem_name"] + ");"
    # account for thread group
    if use_thread_group:
        code_start = code_start.replace("(","(tgrp, ")
    # Canonical: append the shared topology-helper arg via the central helper
    # (NO_XI: the ee_pose family takes s_Xhom, not s_XImats). Mirrors the def's
    # gen_insert_helpers_func_def_params(NO_XI_FLAG=True) so def + call can't drift.
    code_middle += self.gen_insert_helpers_function_call(updated_var_names = var_names, NO_XI_FLAG = True)
    self.gen_add_code_line(code_start + code_middle + code_end)

def gen_end_effector_pose_gradient_hessian_inner(self, use_thread_group = False):
    n = self.robot.get_num_pos()
    n_bfs_levels = self.robot.get_max_bfs_level() + 1 # starts at 0
    all_ees = self.robot.get_leaf_nodes()
    num_ees = len(all_ees)
    # construct the boilerplate and function definition
    func_params = ["s_d2eePos is a pointer to shared memory of size 6*NUM_JOINTS*NUM_JOINTS*NUM_EE where NUM_JOINTS = " + str(n) + " and NUM_EE = " + str(num_ees), \
                   "s_deePos is a pointer to shared memory of size 6*NUM_JOINTS*NUM_EE where NUM_JOINTS = " + str(n) + " and NUM_EE = " + str(num_ees), \
                   "s_q is the vector of joint positions", \
                   "s_Xhom is the pointer to the homogenous transformation matricies ", \
                   "s_dXhom is the pointer to the 1st derivative of the homogenous transformation matricies ", \
                   "s_d2Xhom is the pointer to the 2nd derivative of the homogenous transformation matricies ", \
                   "s_temp is a pointer to helper shared memory of size " + \
                            str(self.gen_end_effector_pose_gradient_hessian_inner_temp_mem_size(include_d2_temp = False)), \
                   "s_d2eeTemp is a pointer to helper memory of size " + \
                            str(self.gen_end_effector_pose_gradient_hessian_d2_temp_mem_size()), \
                   "s_linalg_smem is optional byte-addressed shared memory for cuBLASDx"]
    func_notes = ["Assumes the Xhom and dXhom matricies have already been updated for the given q"]
    func_def_start = "void end_effector_pose_gradient_hessian_inner("
    func_def_middle = "T *s_d2eePos, T *s_deePos, const T *s_q, const T *s_Xhom, const T *s_dXhom, const T *s_d2Xhom, "
    func_def_end = "T *s_temp, T *s_d2eeTemp, unsigned char *s_linalg_smem) {"
    if use_thread_group:
        func_def_start += "cgrps::thread_group tgrp, "
        func_params.insert(0,"tgrp is the handle to the thread_group running this function")
    func_def_middle, func_params = self.gen_insert_helpers_func_def_params(func_def_middle, func_params, -1, NO_XI_FLAG = True)
    func_def = func_def_start + func_def_middle + func_def_end
    # now generate the code
    self.gen_add_func_doc("Computes the Gradient and Hessian of the End Effector Pose with respect to joint position",\
                          func_notes,func_params,None)
    self.gen_add_code_line("template <typename T>")
    self.gen_add_code_line("__device__")
    self.gen_add_code_line(func_def, True)
    #
    # Initial Debug Prints if Requested
    #
    if self.DEBUG_MODE:
        debug_code_lines = ["printf(\"q\\n\"); printMat<T,1," + str(n) + ">(s_q,1);", \
                            "for (int i = 0; i < " + str(n) + "; i++){printf(\"X[%d]\\n\",i); printMat<T,4,4>(&s_Xhom[16*i],4);}"
                            "for (int i = 0; i < " + str(n) + "; i++){printf(\"dX[%d]\\n\",i); printMat<T,4,4>(&s_dXhom[16*i],4);}"]
        self.gen_add_debug_print_code_lines(debug_code_lines,use_thread_group)
    #
    # For each chain we need to (in parallel) multiply the (d)Xmats
    # 
    self.gen_add_code_line("//")
    self.gen_add_code_line("// For each branch/gradient in parallel chain up the transform")
    self.gen_add_code_line("// Keep chaining until reaching the root (starting from the leaves)")
    self.gen_add_code_line("// Keep all gradients (and the value) as we need them for the hessian")
    self.gen_add_code_line("//")
    self.gen_add_code_line("T *s_eeTemp = &s_temp[0]; T *s_deeTemp = &s_temp[" + str(2*16*num_ees) + "];")
    for bfs_level in range(n_bfs_levels): # at most bfs levels of parents to chain
        # if serial chain manipulator then this is easy
        if self.robot.is_serial_chain():
            self.gen_add_code_line("// Serial chain manipulator so optimize as parent is jid-1")
            if bfs_level == 0:
                self.gen_add_code_line("// First set the leaf transforms for eePos and deePos")
                self.gen_add_parallel_loop("ind",str(16*n),use_thread_group)
                self.gen_add_code_line("int djid = ind / 16; int rc = ind % 16; int eeIndStart = 16*" + str(all_ees[0]) + ";")
                self.gen_add_code_line("if(djid == 0){s_eeTemp[ind] = s_Xhom[eeIndStart + rc];}")
                self.gen_add_code_line("const T *s_Xhom_dXhom = grid_xhom_or_dxhom_ptr<T>(s_Xhom, s_dXhom, djid, " + str(all_ees[0]) + ");")
                self.gen_add_code_line("s_deeTemp[ind] = s_Xhom_dXhom[rc];")
                self.gen_add_end_control_flow()
                self.gen_add_sync(use_thread_group)
                if self.DEBUG_MODE:
                    code_lines = ["printf(\"X_chain0\\n\"); printMat<T,4,4>(s_eeTemp,4);", \
                                  "for (int i = 0; i < " + str(n) + "; i++){printf(\"dX_chain0[%d]\\n\",i); printMat<T,4,4>(&s_deeTemp[16*i],4);}"]
                    self.gen_add_debug_print_code_lines(code_lines,use_thread_group)
            else:
                self.gen_add_code_line("// Update with parent transform until you reach the base [level " + str(bfs_level) + "/" + str(n_bfs_levels-1) + "]")
                # get the parents we need at this level working backwards from all_ees
                parent = all_ees[0]
                for i in range(bfs_level):
                    parent = self.robot.get_parent_id(parent)
                # need to swap dst and start each time
                even = bfs_level % 2
                tempDstOffset_ee = 16*num_ees*(even)
                tempSrcOffset_ee = 16*num_ees*(not even)
                tempDstOffset_dee = 16*n*num_ees*(even)
                tempSrcOffset_dee = 16*n*num_ees*(not even)
                self.gen_add_parallel_loop("ind",str(16*n),use_thread_group)
                self.gen_add_code_line("int djid = ind / 16; int rc = ind % 16; int row = rc % 4; int colInd = ind - row;")
                self.gen_add_code_line("if(djid == 0){s_eeTemp[ind + " + str(tempDstOffset_ee) + "] = dot_prod<T,4,4,1>" + \
                                       "(&s_Xhom[16*" + str(parent) + " + row], &s_eeTemp[" + str(tempSrcOffset_ee) + " + colInd]);}")
                self.gen_add_code_line("const T *s_Xhom_dXhom = grid_xhom_or_dxhom_ptr<T>(s_Xhom, s_dXhom, djid, " + str(parent) + ");")
                self.gen_add_code_line("s_deeTemp[ind + " + str(tempDstOffset_dee) + "] = dot_prod<T,4,4,1>" + \
                                       "(&s_Xhom_dXhom[row], &s_deeTemp[" + str(tempSrcOffset_dee) + " + colInd]);")
                self.gen_add_end_control_flow()
                self.gen_add_sync(use_thread_group)
                # make sure to save the offsets for next operations
                if bfs_level == n_bfs_levels - 1:
                    self.gen_add_code_line("int eeTemp_Offset = " + str(tempDstOffset_ee) + ";")
                    self.gen_add_code_line("int deeTemp_Offset = " + str(tempDstOffset_dee) + ";")
                if self.DEBUG_MODE:
                    code_lines = ["printf(\"X_chain_iter[%d]\\n\"," + str(bfs_level) + "); printMat<T,4,4>(&s_eeTemp[" + str(tempDstOffset_ee) + "],4);", \
                                  "for (int i = 0; i < " + str(n) + "; i++){printf(\"dX_chain_iter[%d][%d]\\n\"," + str(bfs_level) + ",i); printMat<T,4,4>(&s_deeTemp[16*i+" + str(tempDstOffset_dee) + "],4);}"]
                    self.gen_add_debug_print_code_lines(code_lines,use_thread_group)
        else:
            # if first loop then just set to transform at the leaf
            if bfs_level == 0:
                self.gen_add_code_line("// First set the leaf transforms for eePos and deePos")
                self.gen_add_parallel_loop("ind",str(16*n*num_ees),use_thread_group)
                self.gen_add_code_line("int rc = ind % 16; int curr_ee = ind / " + str(16*n) + "; int djid = (ind / 16) % " + str(n) + ";")
                select_var_vals = [("int", "eeInd", [str(jid) for jid in all_ees])]
                # make sure to zero out all things not in the chain
                jidChainCode = []
                for eejid in all_ees:
                    jidChain = sorted(self.robot.get_ancestors_by_id(eejid))
                    jidChain.append(eejid)
                    qinds = []
                    for jid in jidChain:
                        qind = self.robot.get_joint_index_q(jid)
                        qinds.extend(qind if isinstance(qind, list) else [qind])
                    code = self.gen_var_in_list("djid", [str(qind) for qind in qinds])
                    jidChainCode.append(code)
                select_var_vals.append(("bool", "inChain", jidChainCode))
                self.gen_add_multi_threaded_select("ind", "<", [str(16*n*(i+1)) for i in range(num_ees)], select_var_vals)
                self.gen_add_code_line("if(djid == 0){s_eeTemp[16*curr_ee + rc] = s_Xhom[16*eeInd + rc];}")
                self.gen_add_code_line("const T *s_Xhom_dXhom = grid_xhom_or_dxhom_ptr<T>(s_Xhom, s_dXhom, djid, eeInd);")
                self.gen_add_code_line("s_deeTemp[ind] = inChain * s_Xhom_dXhom[rc];")
                self.gen_add_end_control_flow()
                self.gen_add_sync(use_thread_group)
                if self.DEBUG_MODE:
                    code_lines = ["for (int i = 0; i < " + str(n*num_ees) + "; i++){printf(\"X_chain level[%d] with dj_ee_id [%d]\\n\"," + str(bfs_level) + ",i); printMat<T,4,4>(&s_eeTemp[16*i],4);}", \
                                  "for (int i = 0; i < " + str(n*num_ees) + "; i++){printf(\"dX_chain level[%d] with dj_ee_id [%d]\\n\"," + str(bfs_level) + ",i); printMat<T,4,4>(&s_deeTemp[16*i],4);}"]
            else:
                self.gen_add_code_line("// Update with parent transform until you reach the base [level " + str(bfs_level) + "/" + str(n_bfs_levels-1) + "]")
                # get the parents we need at this level working backwards from all_ees
                curr_parents = all_ees
                for i in range(bfs_level):
                    curr_parents = [(-1 if jid == -1 else self.robot.get_parent_id(jid)) for jid in curr_parents]
                # need to swap dst and start each time
                even = bfs_level % 2
                tempDstOffset_ee = 16*num_ees*(even)
                tempSrcOffset_ee = 16*num_ees*(not even)
                tempDstOffset_dee = 16*n*num_ees*(even)
                tempSrcOffset_dee = 16*n*num_ees*(not even)
                self.gen_add_parallel_loop("ind",str(16*n*num_ees),use_thread_group)
                self.gen_add_code_line("int rc = ind % 16; int curr_ee = ind / " + str(16*n) + "; int djid = (ind / 16) % " + str(n) + ";")
                self.gen_add_code_line("int row = rc % 4; int colInd = ind - row; int colInd_ee = 16*curr_ee + rc - row;")
                # get parents for this level
                select_var_vals = [("int", "parent_jid", [str(jid) for jid in curr_parents])]
                self.gen_add_multi_threaded_select("ind", "<", [str(16*n*(i+1)) for i in range(num_ees)], select_var_vals)
                if (-1 in curr_parents):
                    self.gen_add_code_line("if(parent_jid == -1){continue;}")
                self.gen_add_code_line("if(djid == 0){s_eeTemp[16*curr_ee + rc + " + str(tempDstOffset_ee) + "] = dot_prod<T,4,4,1>" + \
                                       "(&s_Xhom[16*parent_jid + row], &s_eeTemp[" + str(tempSrcOffset_ee) + " + colInd_ee]);}")
                self.gen_add_code_line("const T *s_Xhom_dXhom = grid_xhom_or_dxhom_ptr<T>(s_Xhom, s_dXhom, djid, parent_jid);")
                self.gen_add_code_line("s_deeTemp[ind + " + str(tempDstOffset_dee) + "] = dot_prod<T,4,4,1>" + \
                                       "(&s_Xhom_dXhom[row], &s_deeTemp[" + str(tempSrcOffset_dee) + " + colInd]);")
                self.gen_add_code_line("// if last loop then save the temp offsets")
                self.gen_add_end_control_flow()
                self.gen_add_sync(use_thread_group)
                # make sure to save the offsets for next operations
                if bfs_level == n_bfs_levels - 1:
                    self.gen_add_code_line("int eeTemp_Offset = " + str(tempDstOffset_ee) + ";")
                    self.gen_add_code_line("int deeTemp_Offset = " + str(tempDstOffset_dee) + ";")
                if self.DEBUG_MODE:
                    code_lines = ["printf(\"X_chain_iter[%d]\\n\"," + str(bfs_level) + "); printMat<T,4,4>(&s_eeTemp[" + str(tempDstOffset_ee) + "],4);", \
                                  "for (int i = 0; i < " + str(n) + "; i++){printf(\"dX_chain_iter[%d][%d]\\n\"," + str(bfs_level) + ",i); printMat<T,4,4>(&s_deeTemp[16*i+" + str(tempDstOffset_dee) + "],4);}"]
    #
    # For each chain we now need to (in parallel) form the d2Xmats
    # 
    self.gen_add_code_line("//")
    self.gen_add_code_line("// For each hessian term in parallel chain up the transform")
    self.gen_add_code_line("// Keep chaining until reaching the root (starting from the leaves)")
    self.gen_add_code_line("//")
    self.gen_add_code_line("// set eeTemp and deeTemp to the right offsets")
    self.gen_add_code_line("s_eeTemp = &s_eeTemp[eeTemp_Offset];")
    self.gen_add_code_line("s_deeTemp = &s_deeTemp[deeTemp_Offset];")
    for bfs_level in range(n_bfs_levels): # at most bfs levels of parents to chain
        # if serial chain manipulator then this is easy
        if self.robot.is_serial_chain():
            self.gen_add_code_line("// Serial chain manipulator so optimize as parent is jid-1")
            if bfs_level == 0:
                self.gen_add_code_line("// First set the leaf transforms")
                self.gen_add_parallel_loop("ind",str(16*n*n),use_thread_group)
                self.gen_add_code_line("int djid_ij = ind / 16; int rc = ind % 16; int djid_i = djid_ij / " + str(n) + "; int djid_j = djid_ij % " + str(n) + ";")
                self.gen_add_code_line("const T *s_Xhom_dXhom_d2Xhom = grid_xhom_or_dxhom_or_d2xhom_ptr<T>(s_Xhom, s_dXhom, s_d2Xhom, djid_i, djid_j, " + str(all_ees[0]) + ");")
                self.gen_add_code_line("s_d2eeTemp[ind] = s_Xhom_dXhom_d2Xhom[rc];")
                self.gen_add_end_control_flow()
                self.gen_add_sync(use_thread_group)
                if self.DEBUG_MODE:
                    code_lines = ["for (int ind = 0; ind < " + str(n*n) + "; ind++){int i = ind / " + str(n) + "; int j = ind % " + str(n) + "; printf(\"d2X_chain[0][%d][%d]\\n\",i,j); printMat<T,4,4>(&s_d2eeTemp[16*ind],4);}"]
                    self.gen_add_debug_print_code_lines(code_lines,use_thread_group)
            else:
                self.gen_add_code_line("// Update with parent transform until you reach the base [level " + str(bfs_level) + "/" + str(n_bfs_levels-1) + "]")
                # get the parents we need at this level working backwards from all_ees
                parent = all_ees[0]
                for i in range(bfs_level):
                    parent = self.robot.get_parent_id(parent)
                # need to swap dst and start each time
                even = bfs_level % 2
                tempDstOffset = 16*n*n*(even)
                tempSrcOffset = 16*n*n*(not even)
                self.gen_add_parallel_loop("ind",str(16*n*n),use_thread_group)
                self.gen_add_code_line("int djid_ij = ind / 16; int rc = ind % 16; int djid_i = djid_ij / " + str(n) + "; int djid_j = djid_ij % " + str(n) + "; int row = rc % 4; int colInd = ind - row;")
                self.gen_add_code_line("const T *s_Xhom_dXhom_d2Xhom = grid_xhom_or_dxhom_or_d2xhom_ptr<T>(s_Xhom, s_dXhom, s_d2Xhom, djid_i, djid_j, " + str(parent) + ");")
                self.gen_add_code_line("s_d2eeTemp[ind + " + str(tempDstOffset) + "] = dot_prod<T,4,4,1>" + \
                                       "(&s_Xhom_dXhom_d2Xhom[row], &s_d2eeTemp[" + str(tempSrcOffset) + " + colInd]);")
                self.gen_add_end_control_flow()
                self.gen_add_sync(use_thread_group)
                if bfs_level == n_bfs_levels - 1:
                    self.gen_add_code_line("int d2eeTemp_Offset = " + str(tempDstOffset) + ";")
                if self.DEBUG_MODE:
                    code_lines = ["for (int ind = 0; ind < " + str(n*n) + "; ind++){int i = ind / " + str(n) + "; int j = ind % " + str(n) + "; printf(\"d2X_chain_iter[%d][%d][%d]\\n\"," + str(bfs_level) + ",i,j); printMat<T,4,4>(&s_d2eeTemp[16*ind + " + str(tempDstOffset) + "],4);}"]
                    self.gen_add_debug_print_code_lines(code_lines,use_thread_group)
        else:
            if bfs_level == 0:
                self.gen_add_code_line("// First set the leaf transforms")
                self.gen_add_parallel_loop("ind",str(16*n*n*num_ees),use_thread_group)
                self.gen_add_code_line("int rc = ind % 16; int djid_ij = (ind / 16) % " + str(n*n) + "; int djid_i = djid_ij / " + str(n) + "; int djid_j = djid_ij % " + str(n) + ";")
                select_var_vals = [("int", "eeInd", [str(jid) for jid in all_ees])]
                # make sure to zero out all things not in the chain
                jidChainCode_i = []
                jidChainCode_j = []
                for eejid in all_ees:
                    jidChain = sorted(self.robot.get_ancestors_by_id(eejid))
                    jidChain.append(eejid)
                    qinds = []
                    for jid in jidChain:
                        qind = self.robot.get_joint_index_q(jid)
                        qinds.extend(qind if isinstance(qind, list) else [qind])
                    code_i = self.gen_var_in_list("djid_i", [str(qind) for qind in qinds])
                    code_j = self.gen_var_in_list("djid_j", [str(qind) for qind in qinds])
                    jidChainCode_i.append(code_i)
                    jidChainCode_j.append(code_j)
                select_var_vals.append(("bool", "inChain_i", jidChainCode_i))
                select_var_vals.append(("bool", "inChain_j", jidChainCode_j))
                self.gen_add_multi_threaded_select("ind", "<", [str(16*n*n*(i+1)) for i in range(num_ees)], select_var_vals)
                self.gen_add_code_line("bool inChain = inChain_i && inChain_j;")
                self.gen_add_code_line("const T *s_Xhom_dXhom_d2Xhom = grid_xhom_or_dxhom_or_d2xhom_ptr<T>(s_Xhom, s_dXhom, s_d2Xhom, djid_i, djid_j, eeInd);")
                self.gen_add_code_line("s_d2eeTemp[ind] = inChain * s_Xhom_dXhom_d2Xhom[rc];")
                self.gen_add_end_control_flow()
                self.gen_add_sync(use_thread_group)
                if self.DEBUG_MODE:
                    code_lines = ["for (int i = 0; i < " + str(n*n*num_ees) + "; i++){printf(\"d2X_chain0[%d]\\n\",i); printMat<T,4,4>(&s_eeTemp[16*i],4);}"]
                    self.gen_add_debug_print_code_lines(code_lines,use_thread_group)
            else:
                self.gen_add_code_line("// Update with parent transform until you reach the base [level " + str(bfs_level) + "/" + str(n_bfs_levels-1) + "]")
                # get the parents we need at this level working backwards from all_ees
                curr_parents = all_ees
                for i in range(bfs_level):
                    curr_parents = [(-1 if jid == -1 else self.robot.get_parent_id(jid)) for jid in curr_parents]
                # need to swap dst and start each time
                even = bfs_level % 2
                tempDstOffset = 16*n*n*num_ees*(even)
                tempSrcOffset = 16*n*n*num_ees*(not even)
                self.gen_add_parallel_loop("ind",str(16*n*n*num_ees),use_thread_group)
                self.gen_add_code_line("int djid_ij = (ind / 16) % " + str(n*n) + "; int djid_i = djid_ij / " + str(n) + "; int djid_j = djid_ij % " + str(n) + ";" + \
                                       "int rc = ind % 16; int row = rc % 4; int colInd = ind - row;")
                # get parents for this level
                select_var_vals = [("int", "parent_jid", [str(jid) for jid in curr_parents])]
                self.gen_add_multi_threaded_select("ind", "<", [str(16*n*n*(i+1)) for i in range(num_ees)], select_var_vals)
                if (-1 in curr_parents):
                    self.gen_add_code_line("if(parent_jid == -1){continue;}")
                self.gen_add_code_line("const T *s_Xhom_dXhom_d2Xhom = grid_xhom_or_dxhom_or_d2xhom_ptr<T>(s_Xhom, s_dXhom, s_d2Xhom, djid_i, djid_j, parent_jid);")
                self.gen_add_code_line("s_d2eeTemp[ind + " + str(tempDstOffset) + "] = dot_prod<T,4,4,1>" + \
                                       "(&s_Xhom_dXhom_d2Xhom[row], &s_d2eeTemp[" + str(tempSrcOffset) + " + colInd]);")
                self.gen_add_end_control_flow()
                self.gen_add_sync(use_thread_group)
                if bfs_level == n_bfs_levels - 1:
                    self.gen_add_code_line("int d2eeTemp_Offset = " + str(tempDstOffset) + ";")
                if self.DEBUG_MODE:
                    code_lines = ["for (int ind = 0; ind < " + str(n*n) + "; ind++){int i = ind / " + str(n) + "; int j = ind % " + str(n) + "; printf(\"d2X_chain_iter[%d][%d][%d]\\n\"," + str(bfs_level) + ",i,j); printMat<T,4,4>(&s_d2eeTemp[16*ind + " + str(tempDstOffset) + "],4);}"]
                    self.gen_add_debug_print_code_lines(code_lines,use_thread_group)

    # Then extract the end-effector position with the given offset(s)
    # TODO handle different offsets for different branches
    self.gen_add_code_line("//")
    self.gen_add_code_line("// Finally Extract the eePos from the Tansforms")
    self.gen_add_code_line("// TODO: ADD OFFSETS")
    self.gen_add_code_line("//")
    self.gen_add_code_line("// set d2eeTemp to the right offsets")
    self.gen_add_code_line("s_d2eeTemp = &s_d2eeTemp[d2eeTemp_Offset];")
    if self.DEBUG_MODE:
        code_lines = ["for (int ind = 0; ind < " + str(n*n) + "; ind++){int i = ind / " + str(n) + "; int j = ind % " + str(n) + "; printf(\"d2X_chain[%d][%d]\\n\",i,j); printMat<T,4,4>(&s_d2eeTemp[16*ind],4);}"]
        self.gen_add_debug_print_code_lines(code_lines,use_thread_group)
    self.gen_add_code_line("// For all n*n*num_ee in parallel")
    self.gen_add_parallel_loop("ind6",str(6*n*n*num_ees),use_thread_group)
    self.gen_add_code_line("int ind = ind6 / 6; int outputInd = ind6 % 6;")
    self.gen_add_code_line("int curr_ee = ind / " + str(n*n) + "; int djid_ij = ind % " + str(n*n) + "; int djid_i = djid_ij / " + str(n) + "; int djid_j = djid_ij % " + str(n) + ";")
    self.gen_add_code_line("int ee_offset = 16 * curr_ee; int dee_offset = ee_offset * " + str(n) + "; int d2ee_offset = dee_offset * " + str(n) + ";")
    self.gen_add_code_line("T *s_Xmat_hom = &s_eeTemp[ee_offset]; T *s_d2Xmat_hom_ij = &s_d2eeTemp[d2ee_offset + 16 * djid_i * " + str(n) + " + 16 * djid_j];")
    self.gen_add_code_line("T *s_dXmat_hom_i = &s_deeTemp[dee_offset + 16 * djid_i]; T *s_dXmat_hom_j = &s_deeTemp[dee_offset + 16 * djid_j];")
    self.gen_add_code_line("T *s_deePos_i = &s_deePos[6*(curr_ee*" + str(n) + " + djid_i)]; int d2eePosOffset_ij = djid_i * " + str(n) + " + djid_j; int d2eePosOffset_ind = " + str(n*n) + "; int d2eePosOffset_ee = curr_ee * 6 * " + str(n*n) + ";")
    # # make sure to zero out all things not in the chain
    # jidChainCode_i = []
    # jidChainCode_j = []
    # for eejid in all_ees:
    #     jidChain = sorted(self.robot.get_ancestors_by_id(eejid))
    #     jidChain.append(eejid)
    #     code_i = self.gen_var_in_list("djid_i", [str(jid) for jid in jidChain])
    #     code_j = self.gen_var_in_list("djid_j", [str(jid) for jid in jidChain])
    #     jidChainCode_i.append(code_i)
    #     jidChainCode_j.append(code_j)
    # select_var_vals.append(("bool", "inChain_i", jidChainCode_i))
    # select_var_vals.append(("bool", "inChain_j", jidChainCode_j))
    # self.gen_add_multi_threaded_select("ind", "<", [str(n*n*(i+1)) for i in range(num_ees)], select_var_vals)
    # self.gen_add_code_line("bool inChain = inChain_i || inChain_j;")

    # xyz is pretty straight forward
    self.gen_add_code_line("// Note: djid_j == 0 computes gradient too")
    self.gen_add_code_line("// xyz is easy")
    self.gen_add_code_line("if (outputInd < 3){", add_indent_after=True)
    self.gen_add_code_line("if (djid_j == 0){s_deePos[outputInd + 6*(curr_ee*" + str(n) + " + djid_i)] = s_dXmat_hom_i[12 + outputInd];}")
    self.gen_add_code_line("s_d2eePos[d2eePosOffset_ee + outputInd*d2eePosOffset_ind + d2eePosOffset_ij] = s_d2Xmat_hom_ij[12 + outputInd];")
    self.gen_add_end_control_flow()

    # roll pitch yaw is a bit more difficult
    self.gen_add_code_line("// roll pitch yaw is a bit more difficult")
    self.gen_add_code_line("// note: d/dz of arctan2(y(z),x(z)) = [-x'(z)y(z)+x(z)y'(z)]/[(x(z)^2 + y(z)^2)]")
    self.gen_add_code_line("//       d/dz of sqrt(f(z)) = f'(z)/2sqrt(f(z))")
    self.gen_add_code_line("//       d2/dz of arctan2(y(z),x(z)) is (bottom*dtop - top*dbottom) / (bottom*bottom) of:")
    self.gen_add_code_line("//          top = -x_prime_i*y + x*y_prime_i")
    self.gen_add_code_line("//          dtop = -x_prime_prime*y + x*y_prime_prime + (i != j)*(-x_prime_i*y_prime_j + x_prime_j*y_prime_i)")
    self.gen_add_code_line("//          bottom = x*x + y*y")
    self.gen_add_code_line("//          dbottom = 2*x*x_prime_j + 2*y*y_prime_j")
    self.gen_add_code_line("else {", add_indent_after=True)
    self.gen_add_code_line("T pitchSqrtTerm = sqrt(s_Xmat_hom[10]*s_Xmat_hom[10] + s_Xmat_hom[6]*s_Xmat_hom[6]);") 
    self.gen_add_code_line("T dpitchSqrtTerm_i = (s_Xmat_hom[10]*s_dXmat_hom_i[10] + s_Xmat_hom[6]*s_dXmat_hom_i[6])/pitchSqrtTerm;")
    self.gen_add_code_line("T dpitchSqrtTerm_j = (s_Xmat_hom[10]*s_dXmat_hom_j[10] + s_Xmat_hom[6]*s_dXmat_hom_j[6])/pitchSqrtTerm;")
    self.gen_add_code_line("T dpitchSqrtTerm_i_top_dj = s_dXmat_hom_j[10]*s_dXmat_hom_i[10] + s_Xmat_hom[10]*s_d2Xmat_hom_ij[10] + ")
    self.gen_add_code_line("                            s_dXmat_hom_j[6]*s_dXmat_hom_i[6] + s_Xmat_hom[6]*s_d2Xmat_hom_ij[6];")
    self.gen_add_code_line("T d2pitchSqrtTerm_ij = (pitchSqrtTerm*dpitchSqrtTerm_i_top_dj - dpitchSqrtTerm_i*dpitchSqrtTerm_j) / (pitchSqrtTerm*pitchSqrtTerm);")
    select_var_vals = [("T", "y",           ["s_Xmat_hom[6]",       "-s_Xmat_hom[2]",      "s_Xmat_hom[1]"]), \
                       ("T", "x",           ["s_Xmat_hom[10]",      "pitchSqrtTerm",       "s_Xmat_hom[0]"]), \
                       ("T", "y_prime_i",   ["s_dXmat_hom_i[6]",    "-s_dXmat_hom_i[2]",   "s_dXmat_hom_i[1]"]), \
                       ("T", "x_prime_i",   ["s_dXmat_hom_i[10]",   "dpitchSqrtTerm_i",    "s_dXmat_hom_i[0]"]), \
                       ("T", "y_prime_j",   ["s_dXmat_hom_j[6]",    "-s_dXmat_hom_j[2]",   "s_dXmat_hom_j[1]"]), \
                       ("T", "x_prime_j",   ["s_dXmat_hom_j[10]",   "dpitchSqrtTerm_j",    "s_dXmat_hom_j[0]"]), \
                       ("T", "y_dprime_ij", ["s_d2Xmat_hom_ij[6]",  "-s_d2Xmat_hom_ij[2]", "s_d2Xmat_hom_ij[1]"]), \
                       ("T", "x_dprime_ij", ["s_d2Xmat_hom_ij[10]", "d2pitchSqrtTerm_ij",  "s_d2Xmat_hom_ij[0]"])]
    self.gen_add_multi_threaded_select("outputInd", "==", [str(i) for i in range(3,6)], select_var_vals)
    self.gen_add_code_line("T top = -x_prime_i*y + x*y_prime_i;     T bottom = x*x + y*y;")
    self.gen_add_code_line("T dtop = -x_dprime_ij*y + x*y_dprime_ij + (djid_i != djid_j)*(-x_prime_i*y_prime_j + x_prime_j*y_prime_i);")
    self.gen_add_code_line("T dbottom = 2*x*x_prime_j + 2*y*y_prime_j;")
    self.gen_add_code_line("if (djid_j == 0){s_deePos[outputInd + 6*(curr_ee*" + str(n) + " + djid_i)] = top/bottom;}")
    self.gen_add_code_line("s_d2eePos[d2eePosOffset_ee + outputInd*d2eePosOffset_ind + d2eePosOffset_ij] = (bottom*dtop - top*dbottom) / (bottom*bottom);")
    self.gen_add_end_control_flow()
    self.gen_add_end_control_flow()
    self.gen_add_end_function()

def gen_end_effector_pose_gradient_hessian_device_temp_mem_size(self):
    n = self.robot.get_num_pos()
    XHom_size, dXhom_size, d2Xhom_size = self.gen_get_Xhom_size()
    wrapper_size = self.gen_topology_helpers_size() + XHom_size + dXhom_size + d2Xhom_size # for Xhom and dXhom and d2Xhom
    return self.gen_end_effector_pose_gradient_hessian_inner_temp_mem_size() + wrapper_size

def gen_end_effector_pose_gradient_hessian_device(self, use_thread_group = False):
    n = self.robot.get_num_pos()
    num_ees = self.robot.get_total_leaf_nodes()
    inner_no_d2_size = self.gen_end_effector_pose_gradient_hessian_inner_temp_mem_size(include_d2_temp = False)
    d2_temp_size = self.gen_end_effector_pose_gradient_hessian_d2_temp_mem_size()
    # construct the boilerplate and function definition
    func_params = ["s_d2eePos is a pointer to shared memory of size 6*NUM_JOINTS*NUM_JOINTS*NUM_EE where NUM_JOINTS = " + str(n) + " and NUM_EE = " + str(num_ees), \
                   "s_deePos is a pointer to shared memory of size 6*NUM_JOINTS*NUM_EE where NUM_JOINTS = " + str(n) + " and NUM_EE = " + str(num_ees), \
                   "s_q is the vector of joint positions", \
                   "d_robotModel is the pointer to the initialized model specific helpers on the GPU (XImats, topology_helpers, etc.)", \
                   "d_workspace is the global scratch buffer; size D2EE_DEVICE_INLINE_WORKSPACE_BYTES<T, RESOURCE_TIER>() bytes (= 0 at TIER_PERF, " + str(d2_temp_size) + "*sizeof(T) at TIER_LITE+). Pass nullptr at TIER_PERF"]
    func_notes = ["Inline-CUDA users: at TIER_LITE/TIER_MINIMAL the d2eeTemp scratch (~" + str(d2_temp_size) + "*sizeof(T) bytes) moves from shared memory to d_workspace, freeing smem for the caller's outer kernel"]
    func_def_start = "void end_effector_pose_gradient_hessian_device("
    func_def_middle = "T *s_d2eePos, T *s_deePos, const T *s_q, "
    func_def_end = "const robotModel<T> *d_robotModel, T *d_workspace = nullptr) {"
    if use_thread_group:
        func_def_start += "cgrps::thread_group tgrp, "
        func_params.insert(0,"tgrp is the handle to the thread_group running this function")
    func_def = func_def_start + func_def_middle + func_def_end
    # then generate the code
    self.gen_add_func_doc("Computes the Gradient and Hessian of the End Effector Pose with respect to joint position",\
                          func_notes,func_params,None)
    self.gen_add_code_line("template <typename T, int RESOURCE_TIER = TIER_PERF>")
    self.gen_add_code_line("__device__")
    self.gen_add_code_line(func_def, True)
    # Smem arena: s_temp is the inner-no-d2 size always (in smem at every tier).
    # s_d2eeTemp is tier-conditional: at TIER_PERF it extends the arena past
    # linalg/topology; at TIER_LITE+ it comes from d_workspace.
    self.gen_XmatsHom_helpers_temp_shared_memory_code(inner_no_d2_size, include_gradients = True, include_hessians = True,
                                                      include_linalg_scratch = True,
                                                      linalg_scratch_bytes = "GRID_EE_LINALG_SHARED_BYTES<T>()")
    self.gen_add_code_line("T *s_d2eeTemp;")
    self.gen_add_code_line("if constexpr (RESOURCE_TIER == TIER_PERF) {", True)
    self.gen_add_code_line("(void)d_workspace;")
    self.gen_add_code_line("s_arena_offset = grid_align_up(s_arena_offset, alignof(T));")
    self.gen_add_code_line("s_d2eeTemp = grid_arena_ptr<T>(s_arena, s_arena_offset);")
    self.gen_add_code_line("s_arena_offset += sizeof(T) * static_cast<size_t>(" + str(d2_temp_size) + ");")
    self.gen_add_end_control_flow()
    self.gen_add_code_line("else {", True)
    self.gen_add_code_line("s_d2eeTemp = d_workspace;")
    self.gen_add_end_control_flow()
    # then load/update XI and run the algo
    self.gen_load_update_XmatsHom_helpers_function_call(use_thread_group, include_gradients = True, include_hessians = True)
    self.gen_end_effector_pose_gradient_hessian_inner_function_call(use_thread_group)
    self.gen_add_end_function()

_D2EE_PICK_FLAGS = [
    # (use_workspace_temp, use_workspace_d2xhom)
    (False, False),   # pick 0: full smem
    (True,  False),   # pick 1: temp -> workspace
    (True,  True),    # pick 2: temp + d2xhom -> workspace
]

def _emit_d2ee_kernel_body_for_flags(self, n, num_ees, use_workspace_temp, use_workspace_d2xhom,
                                     single_call_timing, use_thread_group):
    """Emit the d2ee kernel body specialized for one tier's spill flags.
    Wrapped in a brace pair (caller emits the `if constexpr (...)` head).
    Used by gen_end_effector_pose_gradient_hessian_kernel to emit either a
    single body (collapsed picks) or three branched bodies (divergent picks)."""
    shared_mem_size = self.gen_end_effector_pose_gradient_hessian_inner_temp_mem_size(include_d2_temp = not use_workspace_temp)
    extra_t_buffers = [("s_q", n)] if use_workspace_temp else [("s_q", n), ("s_d2eePos", 6*n*n*num_ees), ("s_deePos", 6*n*num_ees)]
    self.gen_XmatsHom_helpers_temp_shared_memory_code(shared_mem_size, include_gradients = True, include_hessians = True,
                                                      extra_t_buffers = extra_t_buffers,
                                                      include_d2xhom_shared = not use_workspace_d2xhom,
                                                      include_linalg_scratch = True,
                                                      linalg_scratch_bytes = "GRID_EE_LINALG_SHARED_BYTES<T>()")
    if not use_workspace_temp:
        self.gen_add_code_line("(void)d_workspace;")
        self.gen_add_code_line("T *s_d2eeTemp = &s_temp[" + str(self.gen_end_effector_pose_gradient_hessian_inner_temp_mem_size(include_d2_temp = False)) + "];")
    if use_thread_group:
        self.gen_add_code_line("cgrps::thread_group tgrp = TBD;")
    if not single_call_timing:
        self.gen_add_parallel_loop("k","NUM_TIMESTEPS",use_thread_group,block_level = True)
        self.gen_kernel_load_inputs("q","stride_q",str(n),use_thread_group)
        if use_workspace_d2xhom:
            self.gen_add_code_line("T *s_d2XmatsHom = reinterpret_cast<T *>(&d_workspace[k*GRID_WORKSPACE_BYTES_PER_TIMESTEP<T>() + GRID_D2EE_WORKSPACE_D2XHOM_OFFSET_BYTES<T>()]);")
        if use_workspace_temp:
            self.gen_add_code_line("T *s_d2eePos = &d_d2eePos[k*" + str(6*n*n*num_ees) + "];")
            self.gen_add_code_line("T *s_deePos = &d_deePos[k*" + str(6*n*num_ees) + "];")
            self.gen_add_code_line("T *s_d2eeTemp = reinterpret_cast<T *>(&d_workspace[k*GRID_WORKSPACE_BYTES_PER_TIMESTEP<T>() + GRID_D2EE_WORKSPACE_D2EETEMP_OFFSET_BYTES<T>()]);")
        self.gen_add_code_line("// compute")
        self.gen_load_update_XmatsHom_helpers_function_call(use_thread_group, include_gradients = True, include_hessians = True)
        self.gen_end_effector_pose_gradient_hessian_inner_function_call(use_thread_group)
        self.gen_add_sync(use_thread_group)
        if not use_workspace_temp:
            self.gen_kernel_save_result("d2eePos",str(6*n*n*num_ees),str(6*n*n*num_ees),use_thread_group)
            self.gen_kernel_save_result("deePos",str(6*n*num_ees),str(6*n*num_ees),use_thread_group)
        self.gen_add_end_control_flow()
    else:
        self.gen_kernel_load_inputs_single_timing("q",str(n),use_thread_group)
        if use_workspace_d2xhom:
            self.gen_add_code_line("T *s_d2XmatsHom = reinterpret_cast<T *>(&d_workspace[GRID_D2EE_WORKSPACE_D2XHOM_OFFSET_BYTES<T>()]);")
        if use_workspace_temp:
            self.gen_add_code_line("T *s_d2eePos = d_d2eePos;")
            self.gen_add_code_line("T *s_deePos = d_deePos;")
            self.gen_add_code_line("T *s_d2eeTemp = reinterpret_cast<T *>(&d_workspace[GRID_D2EE_WORKSPACE_D2EETEMP_OFFSET_BYTES<T>()]);")
        self.gen_add_code_line("// compute with NUM_TIMESTEPS as NUM_REPS for timing")
        self.gen_add_code_line("for (int rep = 0; rep < NUM_TIMESTEPS; rep++){", True)
        self.gen_anti_licm_input_reload("q",str(n),use_thread_group,feedback_from="d2eePos")
        self.gen_load_update_XmatsHom_helpers_function_call(use_thread_group, include_gradients = True, include_hessians = True)
        self.gen_end_effector_pose_gradient_hessian_inner_function_call(use_thread_group)
        self.gen_anti_licm_output_write("d2eePos")
        self.gen_add_end_control_flow()
        if not use_workspace_temp:
            self.gen_kernel_save_result_single_timing("d2eePos",str(6*n*n*num_ees),use_thread_group)
            self.gen_kernel_save_result_single_timing("deePos",str(6*n*num_ees),use_thread_group)


def gen_end_effector_pose_gradient_hessian_kernel(self, use_thread_group = False, single_call_timing = False):
    n = self.robot.get_num_pos()
    num_ees = self.robot.get_total_leaf_nodes()
    func_params = ["d_d2eePos is the vector of end effector positions gradients", \
                   "d_deePos is the vector of end effector positions gradients", \
                   "d_workspace is the generated global spill workspace", \
                   "d_q is the vector of joint positions", \
                   "stride_q is the stide between each q", \
                   "d_robotModel is the pointer to the initialized model specific helpers on the GPU (XImats, topology_helpers, etc.)", \
                   "num_timesteps is the length of the trajectory points we need to compute over (or overloaded as test_iters for timing)"]
    func_notes = []
    func_def_start = "void end_effector_pose_gradient_hessian_kernel(T *d_d2eePos, T *d_deePos, unsigned char *d_workspace, const T *d_q, const int stride_q, "
    func_def_end = "const robotModel<T> *d_robotModel, const int NUM_TIMESTEPS) {"
    func_def = func_def_start + func_def_end
    if single_call_timing:
        func_def = func_def.replace("(", "_single_timing(")
    self.gen_add_func_doc("Computes the Gradient and Hessian of the End Effector Pose with respect to joint position",\
                          func_notes,func_params,None)
    self.gen_add_code_line("template <typename T, int RESOURCE_TIER = GRID_DEFAULT_RESOURCE_TIER>")
    self.gen_add_code_line("__global__")
    self.gen_add_code_line("__launch_bounds__(tier_max_threads<RESOURCE_TIER>())")
    self.gen_add_code_line(func_def, True)
    # Tier dispatch: when the 3 picks collapse, emit one body (current behavior).
    # When they diverge, emit three if-constexpr branches — each branch is a full
    # body specialized for that tier's spill flags. Smem-bytes constexpr
    # D2EE_POS_DYNAMIC_SHARED_MEM_BYTES<T,TIER>() is already tier-aware.
    picks = getattr(self, "d2ee_spill_tier_3way", (0, 0, 0))
    if picks[0] == picks[1] == picks[2]:
        uwt, uwd = _D2EE_PICK_FLAGS[picks[0]]
        _emit_d2ee_kernel_body_for_flags(self, n, num_ees, uwt, uwd, single_call_timing, use_thread_group)
    else:
        tier_names = ("TIER_PERF", "TIER_LITE", "TIER_MINIMAL")
        for tier_idx, (tier_name, pick) in enumerate(zip(tier_names, picks)):
            uwt, uwd = _D2EE_PICK_FLAGS[pick]
            head = "if constexpr (RESOURCE_TIER == " + tier_name + ") {" if tier_idx == 0 else \
                   "else if constexpr (RESOURCE_TIER == " + tier_name + ") {"
            self.gen_add_code_line(head, True)
            _emit_d2ee_kernel_body_for_flags(self, n, num_ees, uwt, uwd, single_call_timing, use_thread_group)
            self.gen_add_end_control_flow()
    self.gen_add_end_function()

def gen_end_effector_pose_gradient_hessian_host(self, mode = 0):
    # default is to do the full kernel call -- options are for single timing or compute only kernel wrapper
    single_call_timing = True if mode == 1 else False
    compute_only = True if mode == 2 else False

    # define function def and params
    func_params = ["hd_data is the packaged input and output pointers", \
                   "d_robotModel is the pointer to the initialized model specific helpers on the GPU (XImats, topology_helpers, etc.)", \
                   "num_timesteps is the length of the trajectory points we need to compute over (or overloaded as test_iters for timing)", \
                   "streams are pointers to CUDA streams for async memory transfers (if needed)"]
    func_notes = []
    func_def_start = "void end_effector_pose_gradient_hessian(gridData<T, KIND> *hd_data, const robotModel<T> *d_robotModel, const int num_timesteps,"
    func_def_end =   "                            const dim3 block_dimms, const dim3 thread_dimms, cudaStream_t *streams) {"
    if single_call_timing:
        func_def_start = func_def_start.replace("(", "_single_timing(")
        func_def_end = "              " + func_def_end
    if compute_only:
        func_def_start = func_def_start.replace("(", "_compute_only(")
        func_def_end = "             " + func_def_end.replace(", cudaStream_t *streams", "")
    # then generate the code
    self.gen_add_func_doc("Computes the Gradient and Hessian of the End Effector Pose with respect to joint position",\
                          func_notes,func_params,None)
    self.gen_add_code_line("template <typename T, bool USE_COMPRESSED_MEM = false, gridDataKind KIND = GRID_DATA_ALL>")
    self.gen_add_code_line("__host__")
    self.gen_add_code_line(func_def_start)
    self.gen_add_code_line(func_def_end, True)
    self.gen_add_code_line("static_assert(KIND == GRID_DATA_ALL || KIND == GRID_DATA_KINEMATICS, \"end_effector_pose_gradient_hessian requires all-data or kinematics gridData\");")
    func_call_start = "end_effector_pose_gradient_hessian_kernel<T><<<block_dimms,thread_dimms,D2EE_POS_DYNAMIC_SHARED_MEM_BYTES<T>()>>>(hd_data->d_d2eePos,hd_data->d_deePos,hd_data->d_workspace,hd_data->d_q,stride_q,"
    func_call_end = "d_robotModel,num_timesteps);"
    if single_call_timing:
        func_call_start = func_call_start.replace("kernel<T>","kernel_single_timing<T>")
    if not compute_only:
        # start code with memory transfer
        self.gen_add_code_lines(["// start code with memory transfer", \
                                 "int stride_q;", \
                                 "if (USE_COMPRESSED_MEM) {stride_q = NUM_JOINTS; " + \
                                    "gpuErrchk(cudaMemcpyAsync(hd_data->d_q,hd_data->h_q,stride_q*" + \
                                    ("num_timesteps*" if not single_call_timing else "") + "sizeof(T),cudaMemcpyHostToDevice,streams[0]));}", \
                                 "else {stride_q = 3*NUM_JOINTS; " + \
                                    "gpuErrchk(cudaMemcpyAsync(hd_data->d_q_qd_u,hd_data->h_q_qd_u,stride_q*" + \
                                    ("num_timesteps*" if not single_call_timing else "") + "sizeof(T),cudaMemcpyHostToDevice,streams[0]));}", \
                                 "gpuErrchkKernel();"])
    else:
        self.gen_add_code_line("int stride_q = USE_COMPRESSED_MEM ? NUM_JOINTS: 3*NUM_JOINTS;")
    # then compute but adjust for compressed mem and qdd usage
    self.gen_add_code_line("// then call the kernel")
    func_call = func_call_start + func_call_end
    # add in compressed mem adjusts
    func_call_mem_adjust = "if (USE_COMPRESSED_MEM) {" + func_call + "}"
    func_call_mem_adjust2 = "else                    {" + func_call.replace("hd_data->d_q","hd_data->d_q_qd_u") + "}"
    # compule into a set of code
    func_call_code = [func_call_mem_adjust, func_call_mem_adjust2, "gpuErrchkKernel();"]
    # wrap function call in timing (if needed)
    if single_call_timing:
        func_call_code.insert(0,"struct timespec start, end; clock_gettime(CLOCK_MONOTONIC,&start);")
        func_call_code.append("clock_gettime(CLOCK_MONOTONIC,&end);")
    self.gen_add_code_line("if (D2EE_POS_DYNAMIC_SHARED_MEM_BYTES<T>() > GRID_CUDA_TARGET_SHARED_MEM_BYTES) {fprintf(stderr,\"GRID end_effector_pose_gradient_hessian shared-memory request %zu exceeds compile target %d; regenerate with a deeper Hessian spill fallback or a higher GRID_CUDA_TARGET_SHARED_MEM_BYTES.\\n\", D2EE_POS_DYNAMIC_SHARED_MEM_BYTES<T>(), GRID_CUDA_TARGET_SHARED_MEM_BYTES); gpuErrchk(cudaErrorInvalidConfiguration);}")
    self.gen_add_code_line("gpuErrchk(grid_check_dynamic_shared_memory_bytes(\"end_effector_pose_gradient_hessian\", D2EE_POS_DYNAMIC_SHARED_MEM_BYTES<T>()));")
    workspace_bytes = "GRID_WORKSPACE_BYTES_PER_TIMESTEP<T>()" if single_call_timing else "GRID_WORKSPACE_BYTES_PER_TIMESTEP<T>()*static_cast<size_t>(num_timesteps)"
    self.gen_add_code_line("if (GRID_D2EE_USES_WORKSPACE_TEMP) {gpuErrchk(grid_begin_l2_persisting(0, hd_data->d_workspace, " + workspace_bytes + "));}")
    self.gen_add_code_lines(func_call_code)
    self.gen_add_code_line("if (GRID_D2EE_USES_WORKSPACE_TEMP) {gpuErrchk(grid_end_l2_persisting(0));}")
    if not compute_only:
        # then transfer memory back
        self.gen_add_code_lines(["// finally transfer the result back", \
                                 "gpuErrchk(cudaMemcpy(hd_data->h_deePos,hd_data->d_deePos,6*NUM_EES*NUM_JOINTS*" + \
                                    ("num_timesteps*" if not single_call_timing else "") + "sizeof(T),cudaMemcpyDeviceToHost));",
                                 "gpuErrchk(cudaMemcpy(hd_data->h_d2eePos,hd_data->d_d2eePos,6*NUM_EES*NUM_JOINTS*NUM_JOINTS*" + \
                                    ("num_timesteps*" if not single_call_timing else "") + "sizeof(T),cudaMemcpyDeviceToHost));",
                                 "gpuErrchkKernel();"])
    # finally report out timing if requested
    if single_call_timing:
        from ..algo_registry import single_call_printf_line
        self.gen_add_code_line(single_call_printf_line("ee_pose_gradient"))
    self.gen_add_end_function()

def gen_X_single_thread(self, fixed_target_name = ""):
    n = self.robot.get_num_pos()
    has_fixed_target = fixed_target_name != ""

    if has_fixed_target:
        fixed_joint = self.robot.get_fixed_joint_by_name(fixed_target_name)
        fixed_joints = self.robot.get_fixed_joints_ordered_by_id()
        fixed_offset = fixed_joints.index(fixed_joint)
        flange_idx = n + fixed_offset

    self.gen_add_func_doc(
        "Single thread joint transformation matrix accumulation up to joint (tid)",
        [],
        [
            "s_jointXforms is the pointer to the cumulative joint transfomration matrices",
            "s_XmatsHom is the pointer to the homogenous transformation matrices",
            "s_q is the vector of joint positions",
            "tid is the joint index up to compute",
        ],
        None
    )
    self.gen_add_code_line("template <typename T>")
    self.gen_add_code_line("__device__")
    self.gen_add_code_line("void X_single_thread(T *s_jointXforms, T *s_XmatsHom, T *s_q, int tid) {", True)

    if has_fixed_target:
        self.gen_add_code_line("constexpr int NJ = " + str(n) + ";")
        self.gen_add_code_line("constexpr int FLANGE_IDX = " + str(flange_idx) + ";")

    # Update s_XmatsHom from s_q
    # X_hom[0]
    self.gen_add_code_line("s_XmatsHom[0] = static_cast<T>(cos(s_q[0]));")
    self.gen_add_code_line("s_XmatsHom[1] = static_cast<T>(sin(s_q[0]));")
    self.gen_add_code_line("s_XmatsHom[4] = static_cast<T>(-sin(s_q[0]));")
    self.gen_add_code_line("s_XmatsHom[5] = static_cast<T>(cos(s_q[0]));")
    # X_hom[1]
    self.gen_add_code_line("s_XmatsHom[16] = static_cast<T>(cos(s_q[1]));")
    self.gen_add_code_line("s_XmatsHom[18] = static_cast<T>(-sin(s_q[1]));")
    self.gen_add_code_line("s_XmatsHom[20] = static_cast<T>(-sin(s_q[1]));")
    self.gen_add_code_line("s_XmatsHom[22] = static_cast<T>(-cos(s_q[1]));")
    # X_hom[2]
    self.gen_add_code_line("s_XmatsHom[32] = static_cast<T>(cos(s_q[2]));")
    self.gen_add_code_line("s_XmatsHom[34] = static_cast<T>(sin(s_q[2]));")
    self.gen_add_code_line("s_XmatsHom[36] = static_cast<T>(-sin(s_q[2]));")
    self.gen_add_code_line("s_XmatsHom[38] = static_cast<T>(cos(s_q[2]));")
    # X_hom[3]
    self.gen_add_code_line("s_XmatsHom[48] = static_cast<T>(cos(s_q[3]));")
    self.gen_add_code_line("s_XmatsHom[50] = static_cast<T>(sin(s_q[3]));")
    self.gen_add_code_line("s_XmatsHom[52] = static_cast<T>(-sin(s_q[3]));")
    self.gen_add_code_line("s_XmatsHom[54] = static_cast<T>(cos(s_q[3]));")
    # X_hom[4]
    self.gen_add_code_line("s_XmatsHom[64] = static_cast<T>(cos(s_q[4]));")
    self.gen_add_code_line("s_XmatsHom[66] = static_cast<T>(-sin(s_q[4]));")
    self.gen_add_code_line("s_XmatsHom[68] = static_cast<T>(-sin(s_q[4]));")
    self.gen_add_code_line("s_XmatsHom[70] = static_cast<T>(-cos(s_q[4]));")
    # X_hom[5]
    self.gen_add_code_line("s_XmatsHom[80] = static_cast<T>(cos(s_q[5]));")
    self.gen_add_code_line("s_XmatsHom[82] = static_cast<T>(sin(s_q[5]));")
    self.gen_add_code_line("s_XmatsHom[84] = static_cast<T>(-sin(s_q[5]));")
    self.gen_add_code_line("s_XmatsHom[86] = static_cast<T>(cos(s_q[5]));")
    # X_hom[6]
    self.gen_add_code_line("s_XmatsHom[96] = static_cast<T>(cos(s_q[6]));")
    self.gen_add_code_line("s_XmatsHom[98] = static_cast<T>(sin(s_q[6]));")
    self.gen_add_code_line("s_XmatsHom[100] = static_cast<T>(-sin(s_q[6]));")
    self.gen_add_code_line("s_XmatsHom[102] = static_cast<T>(cos(s_q[6]));")

    # Accumulate global transforms up to 'tid'
    if has_fixed_target:
        self.gen_add_code_line("if (tid >= NJ) tid = NJ - 1;")
    else:
        self.gen_add_code_line("if (tid >= " + str(n) + ") tid = " + str(n) + "-1;")
    self.gen_add_code_line("if (tid < 0) { return; }")

    self.gen_add_code_line("{", True)
    self.gen_add_code_line("const T* c0 = &s_XmatsHom[0];")
    self.gen_add_code_line("T* o0 = &s_jointXforms[0];")
    self.gen_add_code_line("o0[0]=c0[0];  o0[1]=c0[1];  o0[2]=c0[2];")
    self.gen_add_code_line("o0[4]=c0[4];  o0[5]=c0[5];  o0[6]=c0[6];")
    self.gen_add_code_line("o0[8]=c0[8];  o0[9]=c0[9];  o0[10]=c0[10];")
    self.gen_add_code_line("o0[12]=c0[12]; o0[13]=c0[13]; o0[14]=c0[14];")
    self.gen_add_code_line("o0[3]=(T)0; o0[7]=(T)0; o0[11]=(T)0; o0[15]=(T)1;")
    self.gen_add_end_control_flow()
    self.gen_add_code_line("if (tid == 0) { return; }")

    self.gen_add_code_line("#pragma unroll")
    self.gen_add_code_line("for (int j = 1; j <= tid; ++j) {", True)
    self.gen_add_code_line("const T* p = &s_jointXforms[(j - 1) * 16];")
    self.gen_add_code_line("const T* c = &s_XmatsHom[j * 16];")
    self.gen_add_code_line("T r0  = p[0]*c[0]   + p[4]*c[1]   + p[8]*c[2];")
    self.gen_add_code_line("T r1  = p[1]*c[0]   + p[5]*c[1]   + p[9]*c[2];")
    self.gen_add_code_line("T r2  = p[2]*c[0]   + p[6]*c[1]   + p[10]*c[2];")
    self.gen_add_code_line("T r4  = p[0]*c[4]   + p[4]*c[5]   + p[8]*c[6];")
    self.gen_add_code_line("T r5  = p[1]*c[4]   + p[5]*c[5]   + p[9]*c[6];")
    self.gen_add_code_line("T r6  = p[2]*c[4]   + p[6]*c[5]   + p[10]*c[6];")
    self.gen_add_code_line("T r8  = p[0]*c[8]   + p[4]*c[9]   + p[8]*c[10];")
    self.gen_add_code_line("T r9  = p[1]*c[8]   + p[5]*c[9]   + p[9]*c[10];")
    self.gen_add_code_line("T r10 = p[2]*c[8]   + p[6]*c[9]   + p[10]*c[10];")
    self.gen_add_code_line("T r12 = p[0]*c[12]  + p[4]*c[13]  + p[8]*c[14]  + p[12];")
    self.gen_add_code_line("T r13 = p[1]*c[12]  + p[5]*c[13]  + p[9]*c[14]  + p[13];")
    self.gen_add_code_line("T r14 = p[2]*c[12]  + p[6]*c[13]  + p[10]*c[14] + p[14];")
    self.gen_add_code_line("T* o = &s_jointXforms[j * 16];")
    self.gen_add_code_line("o[0]=r0;   o[1]=r1;   o[2]=r2;")
    self.gen_add_code_line("o[4]=r4;   o[5]=r5;   o[6]=r6;")
    self.gen_add_code_line("o[8]=r8;   o[9]=r9;   o[10]=r10;")
    self.gen_add_code_line("o[12]=r12; o[13]=r13; o[14]=r14;")
    self.gen_add_code_line("o[3]=(T)0; o[7]=(T)0; o[11]=(T)0; o[15]=(T)1;")
    self.gen_add_end_control_flow()

    if has_fixed_target:
        self.gen_add_code_line("if (tid == NJ - 1) {", True)
        self.gen_add_code_line("const T* T6 = &s_jointXforms[(NJ - 1) * 16];")
        self.gen_add_code_line("T* Tfl = &s_jointXforms[FLANGE_IDX * 16];")
        self.gen_add_code_line("#pragma unroll")
        self.gen_add_code_line("for (int k = 0; k < 16; ++k) Tfl[k] = T6[k];")
        self.gen_add_code_line("T* Tee = &s_jointXforms[NJ * 16];")
        self.gen_add_code_line("const T* Xfix = &s_XmatsHom[FLANGE_IDX * 16];")
        self.gen_add_code_line("mat4_mul(Tfl, Xfix, Tee);")
        self.gen_add_code_line("Tfl[3]=(T)0; Tfl[7]=(T)0; Tfl[11]=(T)0; Tfl[15]=(T)1;")
        self.gen_add_code_line("Tee[3]=(T)0; Tee[7]=(T)0; Tee[11]=(T)0; Tee[15]=(T)1;")
        self.gen_add_end_control_flow()

    self.gen_add_end_function()

def gen_X_warp(self, fixed_target_name = ""):
    n = self.robot.get_num_pos()
    has_fixed_target = fixed_target_name != ""

    if has_fixed_target:
        fixed_joint = self.robot.get_fixed_joint_by_name(fixed_target_name)
        fixed_joints = self.robot.get_fixed_joints_ordered_by_id()
        fixed_offset = fixed_joints.index(fixed_joint)
        flange_idx = n + fixed_offset

    self.gen_add_func_doc(
        "Warp-cooperative joint transformation matrix accumulation up to joint (tid)",
        [],
        [
            "s_jointXforms is the pointer to the cumulative joint transfomration matrices",
            "s_XmatsHom is the pointer to the homogenous transformation matrices",
            "s_q is the vector of joint positions",
            "tid is the joint index up to compute",
        ],
        None
    )
    self.gen_add_code_line("template <typename T>")
    self.gen_add_code_line("__device__ inline void X_warp(")
    self.gen_add_code_line("    T* __restrict__ s_jointXforms,")
    self.gen_add_code_line("    T* __restrict__ s_XmatsHom,")
    self.gen_add_code_line("    const T* __restrict__ s_q,")
    self.gen_add_code_line("    int tid)")
    self.gen_add_code_line("{", True)

    if has_fixed_target:
        self.gen_add_code_line("constexpr int NJ = " + str(n) + ";")
        self.gen_add_code_line("constexpr int FLANGE_IDX = " + str(flange_idx) + ";")

    self.gen_add_code_line("const int lane = threadIdx.x & 31;")
    self.gen_add_code_line("const unsigned mask = 0xFFFFFFFFu;")

    if has_fixed_target:
        self.gen_add_code_line("if (tid >= NJ) tid = NJ - 1;")
    else:
        self.gen_add_code_line("if (tid >= " + str(n) + ") tid = " + str(n) + "-1;")
    self.gen_add_code_line("if (tid < 0) return;")

    self.gen_add_code_line("if (lane <= 6) {", True)
    self.gen_add_code_line("const int j = lane;")
    self.gen_add_code_line("const T c = static_cast<T>(cos(s_q[j]));")
    self.gen_add_code_line("const T s = static_cast<T>(sin(s_q[j]));")
    self.gen_add_code_line("T* X = &s_XmatsHom[j * 16];")
    self.gen_add_code_line("if (j == 0) {", True)
    self.gen_add_code_line("X[0] = static_cast<T>(c);")
    self.gen_add_code_line("X[1] = static_cast<T>(s);")
    self.gen_add_code_line("X[4] = static_cast<T>(-s);")
    self.gen_add_code_line("X[5] = static_cast<T>(c);")
    self.gen_add_end_control_flow()
    self.gen_add_code_line("else if (j == 1) {", True)
    self.gen_add_code_line("X[0] = static_cast<T>(c);")
    self.gen_add_code_line("X[2] = static_cast<T>(-s);")
    self.gen_add_code_line("X[4] = static_cast<T>(-s);")
    self.gen_add_code_line("X[6] = static_cast<T>(-c);")
    self.gen_add_end_control_flow()
    self.gen_add_code_line("else if (j == 2) {", True)
    self.gen_add_code_line("X[0] = static_cast<T>(c);")
    self.gen_add_code_line("X[2] = static_cast<T>(s);")
    self.gen_add_code_line("X[4] = static_cast<T>(-s);")
    self.gen_add_code_line("X[6] = static_cast<T>(c);")
    self.gen_add_end_control_flow()
    self.gen_add_code_line("else if (j == 3) {", True)
    self.gen_add_code_line("X[0] = static_cast<T>(c);")
    self.gen_add_code_line("X[2] = static_cast<T>(s);")
    self.gen_add_code_line("X[4] = static_cast<T>(-s);")
    self.gen_add_code_line("X[6] = static_cast<T>(c);")
    self.gen_add_end_control_flow()
    self.gen_add_code_line("else if (j == 4) {", True)
    self.gen_add_code_line("X[0] = static_cast<T>(c);")
    self.gen_add_code_line("X[2] = static_cast<T>(-s);")
    self.gen_add_code_line("X[4] = static_cast<T>(-s);")
    self.gen_add_code_line("X[6] = static_cast<T>(-c);")
    self.gen_add_end_control_flow()
    self.gen_add_code_line("else if (j == 5) {", True)
    self.gen_add_code_line("X[0] = static_cast<T>(c);")
    self.gen_add_code_line("X[2] = static_cast<T>(s);")
    self.gen_add_code_line("X[4] = static_cast<T>(-s);")
    self.gen_add_code_line("X[6] = static_cast<T>(c);")
    self.gen_add_end_control_flow()
    self.gen_add_code_line("else if (j == 6) {", True)
    self.gen_add_code_line("X[0] = static_cast<T>(c);")
    self.gen_add_code_line("X[2] = static_cast<T>(s);")
    self.gen_add_code_line("X[4] = static_cast<T>(-s);")
    self.gen_add_code_line("X[6] = static_cast<T>(c);")
    self.gen_add_end_control_flow()
    self.gen_add_end_control_flow()
    self.gen_add_code_line("__syncwarp(mask);")

    self.gen_add_code_line("{", True)
    self.gen_add_code_line("const T* c = &s_XmatsHom[0];")
    self.gen_add_code_line("T* o = &s_jointXforms[0];")
    self.gen_add_code_line("if (lane == 0) {")
    self.gen_add_code_line("    o[0]  = c[0];  o[4]  = c[4];  o[8]  = c[8];  o[12] = c[12];")
    self.gen_add_code_line("    o[3]  = (T)0;  o[7]  = (T)0;  o[11] = (T)0;  o[15] = (T)1;")
    self.gen_add_code_line("}")
    self.gen_add_code_line("else if (lane == 1) {")
    self.gen_add_code_line("    o[1]  = c[1];  o[5]  = c[5];  o[9]  = c[9];  o[13] = c[13];")
    self.gen_add_code_line("}")
    self.gen_add_code_line("else if (lane == 2) {")
    self.gen_add_code_line("    o[2]  = c[2];  o[6]  = c[6];  o[10] = c[10]; o[14] = c[14];")
    self.gen_add_code_line("}")
    self.gen_add_end_control_flow()
    self.gen_add_code_line("__syncwarp(mask);")
    self.gen_add_code_line("if (tid == 0) return;")

    self.gen_add_code_line("#pragma unroll")
    self.gen_add_code_line("for (int j = 1; j <= tid; ++j) {", True)
    self.gen_add_code_line("const T* p = &s_jointXforms[(j - 1) * 16];")
    self.gen_add_code_line("const T* c = &s_XmatsHom[j * 16];")
    self.gen_add_code_line("T* o = &s_jointXforms[j * 16];")

    self.gen_add_code_line("if (lane == 0) {", True)
    self.gen_add_code_line("o[0]  = p[0]*c[0]  + p[4]*c[1]  + p[8]*c[2];")
    self.gen_add_code_line("o[4]  = p[0]*c[4]  + p[4]*c[5]  + p[8]*c[6];")
    self.gen_add_code_line("o[8]  = p[0]*c[8]  + p[4]*c[9]  + p[8]*c[10];")
    self.gen_add_code_line("o[12] = p[0]*c[12] + p[4]*c[13] + p[8]*c[14] + p[12];")
    self.gen_add_code_line("o[3]  = (T)0;      o[7]  = (T)0;      o[11] = (T)0;      o[15] = (T)1;")
    self.gen_add_end_control_flow()

    self.gen_add_code_line("else if (lane == 1) {", True)
    self.gen_add_code_line("o[1]  = p[1]*c[0]  + p[5]*c[1]  + p[9]*c[2];")
    self.gen_add_code_line("o[5]  = p[1]*c[4]  + p[5]*c[5]  + p[9]*c[6];")
    self.gen_add_code_line("o[9]  = p[1]*c[8]  + p[5]*c[9]  + p[9]*c[10];")
    self.gen_add_code_line("o[13] = p[1]*c[12] + p[5]*c[13] + p[9]*c[14] + p[13];")
    self.gen_add_end_control_flow()

    self.gen_add_code_line("else if (lane == 2) {", True)
    self.gen_add_code_line("o[2]  = p[2]*c[0]  + p[6]*c[1]  + p[10]*c[2];")
    self.gen_add_code_line("o[6]  = p[2]*c[4]  + p[6]*c[5]  + p[10]*c[6];")
    self.gen_add_code_line("o[10] = p[2]*c[8]  + p[6]*c[9]  + p[10]*c[10];")
    self.gen_add_code_line("o[14] = p[2]*c[12] + p[6]*c[13] + p[10]*c[14] + p[14];")
    self.gen_add_end_control_flow()

    self.gen_add_code_line("__syncwarp(mask);")
    self.gen_add_end_control_flow()

    if has_fixed_target:
        self.gen_add_code_line("if (tid == NJ - 1) {", True)
        self.gen_add_code_line("const T* T6 = &s_jointXforms[(NJ - 1) * 16];")
        self.gen_add_code_line("if (lane == 0) {", True)
        self.gen_add_code_line("T* Tfl = &s_jointXforms[FLANGE_IDX * 16];")
        self.gen_add_code_line("#pragma unroll")
        self.gen_add_code_line("for (int k = 0; k < 16; ++k) Tfl[k] = T6[k];")
        self.gen_add_code_line("T* Tee = &s_jointXforms[NJ * 16];")
        self.gen_add_code_line("const T* Xfix = &s_XmatsHom[FLANGE_IDX * 16];")
        self.gen_add_code_line("mat4_mul(Tfl, Xfix, Tee);")
        self.gen_add_code_line("Tfl[3]=(T)0; Tfl[7]=(T)0; Tfl[11]=(T)0; Tfl[15]=(T)1;")
        self.gen_add_code_line("Tee[3]=(T)0; Tee[7]=(T)0; Tee[11]=(T)0; Tee[15]=(T)1;")
        self.gen_add_end_control_flow()
        self.gen_add_code_line("__syncwarp(mask);")
        self.gen_add_end_control_flow()

    self.gen_add_end_function()

def gen_eepose_and_derivatives(self, use_thread_group = False, fixed_target_name = "",
                               include_pose = True, include_gradient = True, include_hessian = True):
    ee_target_names = [""]
    if fixed_target_name == "all":
        ee_target_names += [fj.name for fj in self.robot.fixed_joints]
    elif fixed_target_name != "":
        ee_target_names += [fixed_target_name]
    for target in ee_target_names:
        if include_pose:
            # first generate the inner helpers
            self.gen_end_effector_pose_inner(use_thread_group, fixed_target_name = target)
            # then generate the device wrappers
            self.gen_end_effector_pose_device(use_thread_group, fixed_target_name = target)
            # then generate the kernels
            self.gen_end_effector_pose_kernel(use_thread_group, single_call_timing = True, fixed_target_name = target)
            self.gen_end_effector_pose_kernel(use_thread_group, single_call_timing = False, fixed_target_name = target)
            # then the host launch wrappers
            self.gen_end_effector_pose_host(0, fixed_target_name = target)
            self.gen_end_effector_pose_host(1, fixed_target_name = target)
            self.gen_end_effector_pose_host(2, fixed_target_name = target)

        if include_gradient:
            # then for the gradient first generate the inner helpers
            self.gen_end_effector_pose_gradient_inner(use_thread_group, fixed_target_name = target)
            # then generate the device wrappers
            self.gen_end_effector_pose_gradient_device(use_thread_group, fixed_target_name = target)
            # then generate the kernels
            self.gen_end_effector_pose_gradient_kernel(use_thread_group,True, fixed_target_name = target)
            self.gen_end_effector_pose_gradient_kernel(use_thread_group,False, fixed_target_name = target)
            # then the host launch wrappers
            self.gen_end_effector_pose_gradient_host(0, fixed_target_name = target)
            self.gen_end_effector_pose_gradient_host(1, fixed_target_name = target)
            self.gen_end_effector_pose_gradient_host(2, fixed_target_name = target)

    if include_hessian:
        # then for the hessian first generate the inner helpers
        self.gen_end_effector_pose_gradient_hessian_inner(use_thread_group)
        # then generate the device wrappers
        self.gen_end_effector_pose_gradient_hessian_device(use_thread_group)
        # then generate the kernels
        self.gen_end_effector_pose_gradient_hessian_kernel(use_thread_group,True)
        self.gen_end_effector_pose_gradient_hessian_kernel(use_thread_group,False)
        # then the host launch wrappers
        self.gen_end_effector_pose_gradient_hessian_host(0)
        self.gen_end_effector_pose_gradient_hessian_host(1)
        self.gen_end_effector_pose_gradient_hessian_host(2)

    if include_pose or include_gradient or include_hessian:
        self.gen_X_single_thread(fixed_target_name = fixed_target_name)
        self.gen_X_warp(fixed_target_name = fixed_target_name)
