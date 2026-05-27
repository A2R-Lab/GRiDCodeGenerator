"""
End Effector Posiitons

TODO: fix throughout this document for fixed_joint support for branched trees and for multiple fixed at once
"""
def gen_end_effector_pose_inner_temp_mem_size(self, fixed_target_name = ""):
    num_ees = self.robot.get_total_leaf_nodes() if fixed_target_name == "" else 1
    return 2*16*num_ees

def gen_end_effector_pose_inner_function_call(self, use_thread_group = False, updated_var_names = None, fixed_target_name = "",
                                              temp_in_smem_expr = "true"):
    var_names = dict( \
        s_Xhom_name = "s_XmatsHom", \
        s_eePos_name = "s_eePos", \
        s_q_name = "s_q", \
        s_topology_helpers_name = "s_topology_helpers", \
        s_temp_name = "s_temp", \
        d_workspace_name = "nullptr", \
        s_linalg_smem_name = "s_linalg_smem", \
    )
    if updated_var_names is not None:
        for key,value in updated_var_names.items():
            var_names[key] = value
    code_start = "end_effector_pose_inner" + ("" if fixed_target_name == "" else "_" + fixed_target_name) + "<T, " + temp_in_smem_expr + ">(" + var_names["s_eePos_name"] + ", " + var_names["s_q_name"] + ", "
    code_middle = var_names["s_Xhom_name"] + ", "
    code_end =  var_names["s_temp_name"] + ", " + var_names["d_workspace_name"] + ", " + var_names["s_linalg_smem_name"] + ");"
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
                   "d_workspace is the global-memory chain workspace used in place of s_temp when !TEMP_IN_SMEM", \
                   "s_linalg_smem is optional byte-addressed shared memory (reserved; unused by this inner)"]
    func_notes = ["Assumes the Xhom matricies have already been updated for the given q", "Defaults to all leave nodes if fixed_target_name is not provided"]
    func_def_start = "void end_effector_pose_inner" + ("" if fixed_target_name == "" else "_" + fixed_target_name) + "("
    func_def_middle = "T *s_eePos, const T *s_q, const T *s_Xhom, "
    func_def_end = "T *s_temp, T *d_workspace, unsigned char *s_linalg_smem) {"
    if use_thread_group:
        func_def_start += "cgrps::thread_group tgrp, "
        func_params.insert(0,"tgrp is the handle to the thread_group running this function")
    func_def_middle, func_params = self.gen_insert_helpers_func_def_params(func_def_middle, func_params, -1, NO_XI_FLAG = True)
    func_def = func_def_start + func_def_middle + func_def_end
    # now generate the code
    self.gen_add_func_doc("Computes the End Effector Position",\
                          func_notes,func_params,None)
    self.gen_add_code_line("template <typename T, bool TEMP_IN_SMEM = true>")
    self.gen_add_code_line("__device__")
    self.gen_add_code_line(func_def, True)
    # Inner-controlled scratch placement: the (tiny, double-buffered) chain
    # workspace moves to d_workspace when !TEMP_IN_SMEM. Reassigning s_temp at the
    # top keeps every s_temp[...] reference below unchanged.
    self.gen_add_code_line("if constexpr (!TEMP_IN_SMEM) { s_temp = d_workspace; } else { (void)d_workspace; }")
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
            # NON-SERIAL / FLOATING-BASE: the legacy dense path computed ALL
            # (djid, ee) pairs at every BFS level and masked the out-of-chain
            # ones with `inChain` (most work computed then discarded). Instead we
            # emit ONE compacted chain-up over only the in-chain (ee, djid) pairs
            # via grid_linalg_indexed_batched_gemm, then break out of the
            # bfs_level loop (the compaction already walks every level). Output is
            # numerically identical: the masked entries contributed exactly 0.
            _emit_eepose_grad_compacted_nonserial(self, n, all_ees, num_ees, use_thread_group)
            break
    if self.robot.is_serial_chain():
        self.gen_add_code_line("//")
        self.gen_add_code_line("// Now extract the eePos from the Transforms")
        self.gen_add_code_line("// TODO: ADD OFFSETS")
        self.gen_add_code_line("//")
        tempOffset = 16*n*num_ees*(bfs_level % 2)
        _emit_eepose_grad_extraction(self, n, num_ees, tempOffset, tempOffset, use_thread_group, ee_compact = False)
    self.gen_add_end_function()

def _emit_eepose_grad_extraction(self, n, num_ees, ee_off, dee_off, use_thread_group, ee_compact = False):
    # Shared eePos extraction: reads the chained ee transform (s_eeTemp at ee_off)
    # and the gradient transform (s_deeTemp at dee_off) and writes s_deePos. When
    # ee_compact is True the ee transform is stored once per ee (slot = deeInd/n)
    # instead of redundantly per (ee, djid) pair (the serial/dense path used the
    # redundant layout and passes ee_off == dee_off, ee_compact == False).
    self.gen_add_parallel_loop("ind",str(6*n*num_ees),use_thread_group)
    self.gen_add_code_line("int outputInd = ind % 6; int deeInd = ind / 6;")
    if ee_compact:
        # ee transform stored once per ee (slot = deeInd / n); deeTemp still full
        # (ee*n + djid) layout (out-of-chain slots pre-zeroed -> 0 gradient out).
        self.gen_add_code_line("T *s_Xmat_hom = &s_eeTemp[" + str(ee_off) + " + 16*(deeInd / " + str(n) + ")]; T *s_dXmat_hom = &s_deeTemp[" + str(dee_off) + " + 16*deeInd];")
    else:
        self.gen_add_code_line("T *s_Xmat_hom = &s_eeTemp[" + str(ee_off) + " + 16*deeInd]; T *s_dXmat_hom = &s_deeTemp[" + str(dee_off) + " + 16*deeInd];")
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

def _emit_eepose_grad_compacted_nonserial(self, n, all_ees, num_ees, use_thread_group):
    # ---- Compacted (in-chain only) gradient chain-up for the non-serial /
    # floating-base case. Replaces the dense n*num_ees-per-level sweep that
    # masked out-of-chain pairs. Numerically identical: out-of-chain gradients
    # are exactly 0, and we keep the deeTemp pairs at the SAME (ee*n + djid)
    # slot layout (pre-zeroed) so the shared extraction reads them unchanged.
    #
    # Layout:
    #   s_eeTemp  : per-ee FK transform, double-buffered. ee `ei` -> slot
    #               parity*num_ees + ei  (16 floats / matrix).
    #   s_deeTemp : per in-chain (ei, djid) gradient transform, double-buffered
    #               at the full (ei*n + djid) layout. parity offset 16*n*num_ees.
    #
    # The chain for ee `ei` is J_e = [ee, parent(ee), ... , root]. At BFS level l
    # the incoming factor is joint J_e[l]; for the gradient pair (ei, djid) that
    # factor is dXhom[djid] iff djid affects J_e[l] (exactly one level per pair),
    # else Xhom[J_e[l]]. Each level is two grid_linalg_indexed_batched_gemm calls
    # (one A_base = s_dXhom for the dX-substituted pairs, one A_base = s_Xhom for
    # the rest) plus the per-ee FK multiply.
    import_q = self.robot.get_joint_index_q
    # `affects(q, j)` mirrors grid_q_index_affects_joint for BOTH base modes:
    # q affects joint j iff q is one of j's q-indices (floating root joint 0 owns
    # the 0..5/6 base coords; a revolute joint owns exactly its single q-index).
    def _affects_set(joint_id):
        q = import_q(joint_id)
        return set(q) if isinstance(q, list) else {q}
    affects = lambda q_index, joint_id: q_index in _affects_set(joint_id)
    # Build per-ee chains and the flat in-chain pair list.
    ee_chains = []          # ee_chains[ei] = [ee, p1, ..., root]
    max_len = 0
    for ee in all_ees:
        chain = [ee]
        cur = ee
        while True:
            par = self.robot.get_parent_id(cur)
            if par == -1:
                break
            chain.append(par)
            cur = par
        ee_chains.append(chain)
        max_len = max(max_len, len(chain))
    # in-chain q-indices per ee (sorted, matches dense djid ascending order)
    ee_qinds = []
    for ei, ee in enumerate(all_ees):
        qset = set()
        for j in ee_chains[ei]:
            q = import_q(j)
            for qi in (q if isinstance(q, list) else [q]):
                qset.add(qi)
        ee_qinds.append(sorted(qset))

    self.gen_add_code_line("// NON-SERIAL: compacted in-chain chain-up (GLASS indexed batched 4x4 GEMM)")
    # zero s_deeTemp (both buffers) so out-of-chain (ee,djid) slots read as 0 in
    # extraction -> exactly the dense `inChain * ...` masked result.
    self.gen_add_parallel_loop("ind", str(2*16*n*num_ees), use_thread_group)
    self.gen_add_code_line("s_deeTemp[ind] = static_cast<T>(0);")
    self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)

    # ---- Level 0: seed eeTemp (per ee) and deeTemp (per in-chain pair).
    # eeTemp[ei] = Xhom[ee]; deeTemp[ei,djid] = (djid affects ee ? dXhom[djid] : Xhom[ee]).
    self.gen_add_code_line("// level 0: seed per-ee FK transform")
    self.gen_add_parallel_loop("ind", str(16*num_ees), use_thread_group)
    self.gen_add_code_line("int rc = ind % 16; int ei = ind / 16;")
    select_var_vals = [("int", "eeInd", [str(jid) for jid in all_ees])]
    self.gen_add_multi_threaded_select("ind", "<", [str(16*(i+1)) for i in range(num_ees)], select_var_vals)
    self.gen_add_code_line("s_eeTemp[ind] = s_Xhom[16*eeInd + rc];")
    self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)

    # level-0 deeTemp seed: for each in-chain pair, the leaf factor.
    seed_pairs = []   # (dst_slot, src_matrix_slot_in_base, base) base in {"X","dX"}
    for ei, ee in enumerate(all_ees):
        for djid in ee_qinds[ei]:
            dst = ei*n + djid
            if affects(djid, ee):
                seed_pairs.append((dst, djid, "dX"))
            else:
                seed_pairs.append((dst, ee, "X"))
    self.gen_add_code_line("// level 0: seed per-(ee,djid) gradient transform (in-chain only)")
    self.gen_add_code_line("static const int grad_seed_dst[] = {" + ", ".join(str(p[0]) for p in seed_pairs) + "};")
    self.gen_add_code_line("static const int grad_seed_src[] = {" + ", ".join(str(p[1]) for p in seed_pairs) + "};")
    self.gen_add_code_line("static const int grad_seed_isdx[] = {" + ", ".join(("1" if p[2] == "dX" else "0") for p in seed_pairs) + "};")
    self.gen_add_parallel_loop("ind", str(16*len(seed_pairs)), use_thread_group)
    self.gen_add_code_line("int rc = ind % 16; int p = ind / 16;")
    self.gen_add_code_line("const T *s_src = grad_seed_isdx[p] ? &s_dXhom[16*grad_seed_src[p]] : &s_Xhom[16*grad_seed_src[p]];")
    self.gen_add_code_line("s_deeTemp[16*grad_seed_dst[p] + rc] = s_src[rc];")
    self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)

    # ---- Levels 1..max_len-1: chain up by the parent factor at that level.
    for level in range(1, max_len):
        even = level % 2
        ee_dst_off = 16*num_ees*even
        ee_src_off = 16*num_ees*(not even)
        dee_dst_off = 16*n*num_ees*even
        dee_src_off = 16*n*num_ees*(not even)
        # eeTemp: for each ee whose chain has a joint at this level: dst = parent * src.
        ee_a, ee_b, ee_c = [], [], []   # matrix slots
        for ei in range(num_ees):
            if level < len(ee_chains[ei]):
                par = ee_chains[ei][level]
                ee_a.append(par)                                  # s_Xhom slot
                ee_b.append((ee_src_off // 16) + ei)              # s_eeTemp src slot
                ee_c.append((ee_dst_off // 16) + ei)              # s_eeTemp dst slot
        # deeTemp: split into dX-substituted and X-only pair lists at this level.
        dee_dx_a, dee_dx_b, dee_dx_c = [], [], []
        dee_x_a,  dee_x_b,  dee_x_c  = [], [], []
        # carry-forward copies: ees whose chain already reached the root must keep
        # ping-ponging so EVERY ee lands in the same final-parity buffer (the dense
        # path assumed uniform final parity; copying makes us robust regardless).
        ee_carry_src, ee_carry_dst = [], []   # eeTemp matrix slots
        dee_carry_src, dee_carry_dst = [], []  # deeTemp matrix slots
        for ei in range(num_ees):
            if level >= len(ee_chains[ei]):
                ee_carry_src.append((ee_src_off // 16) + ei)
                ee_carry_dst.append((ee_dst_off // 16) + ei)
                for djid in ee_qinds[ei]:
                    dee_carry_src.append((dee_src_off // 16) + ei*n + djid)
                    dee_carry_dst.append((dee_dst_off // 16) + ei*n + djid)
        for ei in range(num_ees):
            if level >= len(ee_chains[ei]):
                continue
            par = ee_chains[ei][level]
            for djid in ee_qinds[ei]:
                src_slot = (dee_src_off // 16) + ei*n + djid
                dst_slot = (dee_dst_off // 16) + ei*n + djid
                if affects(djid, par):
                    dee_dx_a.append(djid)            # s_dXhom slot
                    dee_dx_b.append(src_slot)
                    dee_dx_c.append(dst_slot)
                else:
                    dee_x_a.append(par)              # s_Xhom slot
                    dee_x_b.append(src_slot)
                    dee_x_c.append(dst_slot)
        self.gen_add_code_line("// level " + str(level) + "/" + str(max_len-1) + ": chain up by parent factor")
        sfx = "_l" + str(level)
        # eeTemp FK multiply
        if ee_a:
            self.gen_add_code_line("static const int ee_a" + sfx + "[] = {" + ", ".join(map(str, ee_a)) + "};")
            self.gen_add_code_line("static const int ee_b" + sfx + "[] = {" + ", ".join(map(str, ee_b)) + "};")
            self.gen_add_code_line("static const int ee_c" + sfx + "[] = {" + ", ".join(map(str, ee_c)) + "};")
            self.gen_add_code_line("grid_linalg_indexed_batched_gemm<T, 4>(" + str(len(ee_a)) + ", ee_a" + sfx + ", ee_b" + sfx + ", ee_c" + sfx + ", s_Xhom, s_eeTemp, s_eeTemp);")
        # deeTemp X-only multiply
        if dee_x_a:
            self.gen_add_code_line("static const int dee_xa" + sfx + "[] = {" + ", ".join(map(str, dee_x_a)) + "};")
            self.gen_add_code_line("static const int dee_xb" + sfx + "[] = {" + ", ".join(map(str, dee_x_b)) + "};")
            self.gen_add_code_line("static const int dee_xc" + sfx + "[] = {" + ", ".join(map(str, dee_x_c)) + "};")
            self.gen_add_code_line("grid_linalg_indexed_batched_gemm<T, 4>(" + str(len(dee_x_a)) + ", dee_xa" + sfx + ", dee_xb" + sfx + ", dee_xc" + sfx + ", s_Xhom, s_deeTemp, s_deeTemp);")
        # deeTemp dX-substituted multiply
        if dee_dx_a:
            self.gen_add_code_line("static const int dee_da" + sfx + "[] = {" + ", ".join(map(str, dee_dx_a)) + "};")
            self.gen_add_code_line("static const int dee_db" + sfx + "[] = {" + ", ".join(map(str, dee_dx_b)) + "};")
            self.gen_add_code_line("static const int dee_dc" + sfx + "[] = {" + ", ".join(map(str, dee_dx_c)) + "};")
            self.gen_add_code_line("grid_linalg_indexed_batched_gemm<T, 4>(" + str(len(dee_dx_a)) + ", dee_da" + sfx + ", dee_db" + sfx + ", dee_dc" + sfx + ", s_dXhom, s_deeTemp, s_deeTemp);")
        # carry forward already-finished ees (copy src->dst parity, unchanged value)
        if ee_carry_src:
            self.gen_add_code_line("static const int ee_csrc" + sfx + "[] = {" + ", ".join(map(str, ee_carry_src)) + "};")
            self.gen_add_code_line("static const int ee_cdst" + sfx + "[] = {" + ", ".join(map(str, ee_carry_dst)) + "};")
            self.gen_add_parallel_loop("ind", str(16*len(ee_carry_src)), use_thread_group)
            self.gen_add_code_line("int rc = ind % 16; int c = ind / 16;")
            self.gen_add_code_line("s_eeTemp[16*ee_cdst" + sfx + "[c] + rc] = s_eeTemp[16*ee_csrc" + sfx + "[c] + rc];")
            self.gen_add_end_control_flow()
        if dee_carry_src:
            self.gen_add_code_line("static const int dee_csrc" + sfx + "[] = {" + ", ".join(map(str, dee_carry_src)) + "};")
            self.gen_add_code_line("static const int dee_cdst" + sfx + "[] = {" + ", ".join(map(str, dee_carry_dst)) + "};")
            self.gen_add_parallel_loop("ind", str(16*len(dee_carry_src)), use_thread_group)
            self.gen_add_code_line("int rc = ind % 16; int c = ind / 16;")
            self.gen_add_code_line("s_deeTemp[16*dee_cdst" + sfx + "[c] + rc] = s_deeTemp[16*dee_csrc" + sfx + "[c] + rc];")
            self.gen_add_end_control_flow()
        if ee_carry_src or dee_carry_src:
            self.gen_add_sync(use_thread_group)
        # NOTE: grid_linalg_indexed_batched_gemm already issues a trailing
        # __syncthreads(); the three calls in this level read distinct buffers /
        # disjoint c_idx slots so they are independent. The carry-forward copies
        # write the OTHER parity half (distinct slots) so they are independent too.
    # final parity after the last level
    final_level = max_len - 1
    final_even = final_level % 2
    ee_final_off = 16*num_ees*final_even
    dee_final_off = 16*n*num_ees*final_even
    self.gen_add_code_line("//")
    self.gen_add_code_line("// extract eePos from the compacted transforms (ee transform per-ee)")
    self.gen_add_code_line("//")
    # extraction reads eeTemp per-ee (slot deeInd/n) at ee_final_off and deeTemp
    # at the full (ee*n+djid) layout at dee_final_off.
    _emit_eepose_grad_extraction(self, n, num_ees, ee_final_off, dee_final_off, use_thread_group, ee_compact = True)

def _eepose_chain_metadata(self, all_ees):
    # Shared topology pre-compute for the compacted non-serial ee chains.
    # Returns (ee_chains, ee_qinds, max_len, affects) where ee_chains[ei] is the
    # leaf->root joint list and ee_qinds[ei] the sorted in-chain q-indices.
    # `affects(q, j)` mirrors grid_q_index_affects_joint for BOTH base modes.
    def _affects_set(joint_id):
        q = self.robot.get_joint_index_q(joint_id)
        return set(q) if isinstance(q, list) else {q}
    affects = lambda q_index, joint_id: q_index in _affects_set(joint_id)
    ee_chains, max_len = [], 0
    for ee in all_ees:
        chain, cur = [ee], ee
        while True:
            par = self.robot.get_parent_id(cur)
            if par == -1:
                break
            chain.append(par); cur = par
        ee_chains.append(chain); max_len = max(max_len, len(chain))
    ee_qinds = []
    for ei, ee in enumerate(all_ees):
        qset = set()
        for j in ee_chains[ei]:
            q = self.robot.get_joint_index_q(j)
            for qi in (q if isinstance(q, list) else [q]):
                qset.add(qi)
        ee_qinds.append(sorted(qset))
    return ee_chains, ee_qinds, max_len, affects

def _emit_eepose_hess_compacted_nonserial(self, n, all_ees, num_ees, use_thread_group):
    # ---- Compacted (in-chain only) gradient + hessian chain-up for the
    # non-serial / floating-base case. Replaces the dense n*n*num_ees-per-level
    # quadratic sweep that masked out-of-chain (i,j,ee) triples. Numerically
    # identical: out-of-chain entries are exactly 0, and we keep the SAME slot
    # layouts (pre-zeroed) so the shared extraction reads them unchanged.
    #
    # Layouts (matching the dense hessian extraction):
    #   s_eeTemp   : per-ee FK transform, slot ei                 (dbl-buf 16*num_ees)
    #   s_deeTemp  : per-(ei,djid) gradient, slot ei*n + djid      (dbl-buf 16*n*num_ees)
    #   s_d2eeTemp : per-(ei,i,j) hessian, slot ei*n*n + i*n + j   (dbl-buf 16*n*n*num_ees)
    ee_chains, ee_qinds, max_len, affects = _eepose_chain_metadata(self, all_ees)

    # ============ Phase 1: gradient chain (eeTemp + deeTemp) ============
    self.gen_add_code_line("// NON-SERIAL: compacted gradient chain (GLASS indexed batched 4x4 GEMM)")
    self.gen_add_parallel_loop("ind", str(2*16*n*num_ees), use_thread_group)
    self.gen_add_code_line("s_deeTemp[ind] = static_cast<T>(0);")
    self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)
    # level 0 eeTemp seed
    self.gen_add_code_line("// level 0: seed per-ee FK transform")
    self.gen_add_parallel_loop("ind", str(16*num_ees), use_thread_group)
    self.gen_add_code_line("int rc = ind % 16; int ei = ind / 16;")
    select_var_vals = [("int", "eeInd", [str(jid) for jid in all_ees])]
    self.gen_add_multi_threaded_select("ind", "<", [str(16*(i+1)) for i in range(num_ees)], select_var_vals)
    self.gen_add_code_line("s_eeTemp[ind] = s_Xhom[16*eeInd + rc];")
    self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)
    # level 0 deeTemp seed
    seed_pairs = []
    for ei, ee in enumerate(all_ees):
        for djid in ee_qinds[ei]:
            dst = ei*n + djid
            seed_pairs.append((dst, djid, "dX") if affects(djid, ee) else (dst, ee, "X"))
    self.gen_add_code_line("// level 0: seed per-(ee,djid) gradient transform (in-chain only)")
    self.gen_add_code_line("static const int hgrad_seed_dst[] = {" + ", ".join(str(p[0]) for p in seed_pairs) + "};")
    self.gen_add_code_line("static const int hgrad_seed_src[] = {" + ", ".join(str(p[1]) for p in seed_pairs) + "};")
    self.gen_add_code_line("static const int hgrad_seed_isdx[] = {" + ", ".join(("1" if p[2] == "dX" else "0") for p in seed_pairs) + "};")
    self.gen_add_parallel_loop("ind", str(16*len(seed_pairs)), use_thread_group)
    self.gen_add_code_line("int rc = ind % 16; int p = ind / 16;")
    self.gen_add_code_line("const T *s_src = hgrad_seed_isdx[p] ? &s_dXhom[16*hgrad_seed_src[p]] : &s_Xhom[16*hgrad_seed_src[p]];")
    self.gen_add_code_line("s_deeTemp[16*hgrad_seed_dst[p] + rc] = s_src[rc];")
    self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)
    # gradient chain-up levels
    for level in range(1, max_len):
        even = level % 2
        ee_dst, ee_src = 16*num_ees*even // 16, 16*num_ees*(not even) // 16
        dee_dst, dee_src = 16*n*num_ees*even // 16, 16*n*num_ees*(not even) // 16
        ee_a, ee_b, ee_c = [], [], []
        dee_dx_a, dee_dx_b, dee_dx_c = [], [], []
        dee_x_a, dee_x_b, dee_x_c = [], [], []
        ee_csrc, ee_cdst, dee_csrc, dee_cdst = [], [], [], []
        for ei in range(num_ees):
            if level >= len(ee_chains[ei]):
                ee_csrc.append(ee_src + ei); ee_cdst.append(ee_dst + ei)
                for djid in ee_qinds[ei]:
                    dee_csrc.append(dee_src + ei*n + djid); dee_cdst.append(dee_dst + ei*n + djid)
                continue
            par = ee_chains[ei][level]
            ee_a.append(par); ee_b.append(ee_src + ei); ee_c.append(ee_dst + ei)
            for djid in ee_qinds[ei]:
                s = dee_src + ei*n + djid; d = dee_dst + ei*n + djid
                if affects(djid, par):
                    dee_dx_a.append(djid); dee_dx_b.append(s); dee_dx_c.append(d)
                else:
                    dee_x_a.append(par); dee_x_b.append(s); dee_x_c.append(d)
        self.gen_add_code_line("// gradient level " + str(level) + "/" + str(max_len-1))
        sfx = "_hg" + str(level)
        _emit_idx_gemm(self, "ee" + sfx, ee_a, ee_b, ee_c, "s_Xhom", "s_eeTemp", "s_eeTemp")
        _emit_idx_gemm(self, "dx" + sfx, dee_x_a, dee_x_b, dee_x_c, "s_Xhom", "s_deeTemp", "s_deeTemp")
        _emit_idx_gemm(self, "dd" + sfx, dee_dx_a, dee_dx_b, dee_dx_c, "s_dXhom", "s_deeTemp", "s_deeTemp")
        _emit_carry_copy(self, "eec" + sfx, ee_csrc, ee_cdst, "s_eeTemp", use_thread_group)
        _emit_carry_copy(self, "dec" + sfx, dee_csrc, dee_cdst, "s_deeTemp", use_thread_group)
        if ee_csrc or dee_csrc:
            self.gen_add_sync(use_thread_group)
    final_even = (max_len - 1) % 2
    ee_final = 16*num_ees*final_even
    dee_final = 16*n*num_ees*final_even

    # ============ Phase 2: hessian chain (d2eeTemp) ============
    self.gen_add_code_line("// NON-SERIAL: compacted hessian chain (GLASS indexed batched 4x4 GEMM)")
    self.gen_add_parallel_loop("ind", str(2*16*n*n*num_ees), use_thread_group)
    self.gen_add_code_line("s_d2eeTemp[ind] = static_cast<T>(0);")
    self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)
    # in-chain (ei, i, j) triples. The d2Xhom matrix slot must match
    # grid_d2xhom_offset: floating base stores a dense (q_i,q_j) block (slot
    # i*n+j); fixed base stores only the diagonal (slot i, and i==j there since a
    # fixed-base joint is owned by a single q-index so "both affect" => i==j).
    floating = self.robot.floating_base
    def d2_factor(i, j, joint):
        ai, aj = affects(i, joint), affects(j, joint)
        if ai and aj:
            return ("d2", (i*n + j) if floating else i)   # s_d2Xhom slot grid_d2xhom_offset
        if ai:
            return ("dxi", i)               # s_dXhom slot i
        if aj:
            return ("dxj", j)               # s_dXhom slot j
        return ("x", joint)                 # s_Xhom slot
    # level 0 seed
    seed_d2, seed_dxi, seed_dxj, seed_x = [], [], [], []  # (dst, src)
    for ei, ee in enumerate(all_ees):
        for i in ee_qinds[ei]:
            for j in ee_qinds[ei]:
                dst = ei*n*n + i*n + j
                kind, src = d2_factor(i, j, ee)
                if kind == "d2": seed_d2.append((dst, src))
                elif kind == "dxi": seed_dxi.append((dst, src))
                elif kind == "dxj": seed_dxj.append((dst, src))
                else: seed_x.append((dst, src))
    self.gen_add_code_line("// hessian level 0: seed per-(ee,i,j) transform (in-chain only)")
    _emit_seed_copy(self, "hs_d2", seed_d2, "s_d2Xhom", "s_d2eeTemp", use_thread_group)
    _emit_seed_copy(self, "hs_di", seed_dxi, "s_dXhom", "s_d2eeTemp", use_thread_group)
    _emit_seed_copy(self, "hs_dj", seed_dxj, "s_dXhom", "s_d2eeTemp", use_thread_group)
    _emit_seed_copy(self, "hs_x", seed_x, "s_Xhom", "s_d2eeTemp", use_thread_group)
    self.gen_add_sync(use_thread_group)
    # hessian chain-up levels
    for level in range(1, max_len):
        even = level % 2
        d2_dst, d2_src = 16*n*n*num_ees*even // 16, 16*n*n*num_ees*(not even) // 16
        g = {"d2": ([], [], []), "dxi": ([], [], []), "dxj": ([], [], []), "x": ([], [], [])}
        carry_src, carry_dst = [], []
        for ei in range(num_ees):
            if level >= len(ee_chains[ei]):
                for i in ee_qinds[ei]:
                    for j in ee_qinds[ei]:
                        carry_src.append(d2_src + ei*n*n + i*n + j); carry_dst.append(d2_dst + ei*n*n + i*n + j)
                continue
            par = ee_chains[ei][level]
            for i in ee_qinds[ei]:
                for j in ee_qinds[ei]:
                    s = d2_src + ei*n*n + i*n + j; d = d2_dst + ei*n*n + i*n + j
                    kind, src = d2_factor(i, j, par)
                    g[kind][0].append(src); g[kind][1].append(s); g[kind][2].append(d)
        self.gen_add_code_line("// hessian level " + str(level) + "/" + str(max_len-1))
        sfx = "_hh" + str(level)
        _emit_idx_gemm(self, "d2" + sfx, g["d2"][0], g["d2"][1], g["d2"][2], "s_d2Xhom", "s_d2eeTemp", "s_d2eeTemp")
        _emit_idx_gemm(self, "di" + sfx, g["dxi"][0], g["dxi"][1], g["dxi"][2], "s_dXhom", "s_d2eeTemp", "s_d2eeTemp")
        _emit_idx_gemm(self, "dj" + sfx, g["dxj"][0], g["dxj"][1], g["dxj"][2], "s_dXhom", "s_d2eeTemp", "s_d2eeTemp")
        _emit_idx_gemm(self, "xx" + sfx, g["x"][0], g["x"][1], g["x"][2], "s_Xhom", "s_d2eeTemp", "s_d2eeTemp")
        _emit_carry_copy(self, "d2c" + sfx, carry_src, carry_dst, "s_d2eeTemp", use_thread_group)
        if carry_src:
            self.gen_add_sync(use_thread_group)
    d2_final = 16*n*n*num_ees*final_even

    # ============ Phase 3: extraction (rebase pointers to final parity) ============
    self.gen_add_code_line("// rebase to the final-parity buffers, then extract")
    self.gen_add_code_line("s_eeTemp = &s_eeTemp[" + str(ee_final) + "];")
    self.gen_add_code_line("s_deeTemp = &s_deeTemp[" + str(dee_final) + "];")
    self.gen_add_code_line("s_d2eeTemp = &s_d2eeTemp[" + str(d2_final) + "];")
    _emit_eepose_hess_extraction(self, n, num_ees, use_thread_group)

def _emit_idx_gemm(self, name, a, b, c, A_base, B_base, C_base):
    if not a:
        return
    self.gen_add_code_line("static const int " + name + "_a[] = {" + ", ".join(map(str, a)) + "};")
    self.gen_add_code_line("static const int " + name + "_b[] = {" + ", ".join(map(str, b)) + "};")
    self.gen_add_code_line("static const int " + name + "_c[] = {" + ", ".join(map(str, c)) + "};")
    self.gen_add_code_line("grid_linalg_indexed_batched_gemm<T, 4>(" + str(len(a)) + ", " + name + "_a, " + name + "_b, " + name + "_c, " + A_base + ", " + B_base + ", " + C_base + ");")

def _emit_seed_copy(self, name, pairs, src_base, dst_base, use_thread_group):
    # pairs: list of (dst_slot, src_slot). Copies 4x4 from src_base[src] to dst_base[dst].
    if not pairs:
        return
    self.gen_add_code_line("static const int " + name + "_dst[] = {" + ", ".join(str(p[0]) for p in pairs) + "};")
    self.gen_add_code_line("static const int " + name + "_src[] = {" + ", ".join(str(p[1]) for p in pairs) + "};")
    self.gen_add_parallel_loop("ind", str(16*len(pairs)), use_thread_group)
    self.gen_add_code_line("int rc = ind % 16; int p = ind / 16;")
    self.gen_add_code_line(dst_base + "[16*" + name + "_dst[p] + rc] = " + src_base + "[16*" + name + "_src[p] + rc];")
    self.gen_add_end_control_flow()

def _emit_carry_copy(self, name, src, dst, base, use_thread_group):
    # carry-forward: copy already-finished slots into the other parity half.
    if not src:
        return
    self.gen_add_code_line("static const int " + name + "_src[] = {" + ", ".join(map(str, src)) + "};")
    self.gen_add_code_line("static const int " + name + "_dst[] = {" + ", ".join(map(str, dst)) + "};")
    self.gen_add_parallel_loop("ind", str(16*len(src)), use_thread_group)
    self.gen_add_code_line("int rc = ind % 16; int c = ind / 16;")
    self.gen_add_code_line(base + "[16*" + name + "_dst[c] + rc] = " + base + "[16*" + name + "_src[c] + rc];")
    self.gen_add_end_control_flow()

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
    # Per-tier eegrad_temp byte offset. The shared GRID_EE_GRAD_WORKSPACE_TEMP_OFFSET_BYTES
    # macro keys off the single-valued PERF-pick GRID_EE_GRAD_USES_WORKSPACE_DXHOM, so it
    # would collide with the spilled dXhom region at tiers whose pick spills dXhom but whose
    # PERF pick does not (e.g. go2 ee_grad = (0,0,2)). Compute the offset locally from THIS
    # tier's use_workspace_dxhom so the temp arena always lands past the spilled dXhom region.
    eegrad_temp_off = "GRID_EE_GRAD_WORKSPACE_DXHOM_OFFSET_BYTES<T>()"
    if use_workspace_dxhom:
        eegrad_temp_off += " + sizeof(T) * static_cast<size_t>(DXHOM_T_COUNT)"
    if use_thread_group:
        self.gen_add_code_line("cgrps::thread_group tgrp = TBD;")
    if not single_call_timing:
        self.gen_add_parallel_loop("k","NUM_TIMESTEPS",use_thread_group,block_level = True)
        self.gen_kernel_load_inputs("q","stride_q",str(n),use_thread_group)
        if use_workspace_dxhom:
            self.gen_add_code_line("T *s_dXmatsHom = reinterpret_cast<T *>(&d_workspace[k*GRID_WORKSPACE_BYTES_PER_TIMESTEP<T>() + GRID_EE_GRAD_WORKSPACE_DXHOM_OFFSET_BYTES<T>()]);")
        if use_workspace_temp:
            self.gen_add_code_line("T *s_deePos = &d_deePos[k*" + str(6*n*num_ees) + "];")
            self.gen_add_code_line("T *s_eegrad_temp = reinterpret_cast<T *>(&d_workspace[k*GRID_WORKSPACE_BYTES_PER_TIMESTEP<T>() + " + eegrad_temp_off + "]);")
            # Whole inner arena spilled -> smem s_temp is null. Repoint it at the
            # spilled workspace so the XmatsHom helper's sincos scratch is backed.
            self.gen_add_code_line("s_temp = s_eegrad_temp;")
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
            self.gen_add_code_line("T *s_eegrad_temp = reinterpret_cast<T *>(&d_workspace[" + eegrad_temp_off + "]);")
            # See note above: repoint the null smem s_temp at the spilled workspace.
            self.gen_add_code_line("s_temp = s_eegrad_temp;")
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
    # Per-tier gate: arm L2 persistence if ANY tier routes the chain workspace
    # through d_workspace (runtime RESOURCE_TIER may differ from the PERF pick).
    self.gen_add_code_line("if (GRID_EE_GRAD_USES_WORKSPACE_TEMP_ANY) {gpuErrchk(grid_begin_l2_persisting(0, hd_data->d_workspace, " + workspace_bytes + "));}")
    self.gen_add_code_lines(func_call_code)
    self.gen_add_code_line("if (GRID_EE_GRAD_USES_WORKSPACE_TEMP_ANY) {gpuErrchk(grid_end_l2_persisting(0));}")
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

def gen_end_effector_pose_gradient_hessian_inner_function_call(self, use_thread_group = False, updated_var_names = None,
                                                               d2temp_in_smem_expr = "true"):
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
        d_workspace_name = "nullptr", \
        s_linalg_smem_name = "s_linalg_smem", \
    )
    if updated_var_names is not None:
        for key,value in updated_var_names.items():
            var_names[key] = value
    code_start = "end_effector_pose_gradient_hessian_inner<T, " + d2temp_in_smem_expr + ">(" + var_names["s_d2eePos_name"] + ", " + var_names["s_deePos_name"] + ", " + var_names["s_q_name"] + ", "
    code_middle = var_names["s_Xhom_name"] + ", " + var_names["s_dXhom_name"] + ", " + var_names["s_d2Xhom_name"] + ", "
    code_end =  var_names["s_temp_name"] + ", " + var_names["s_d2eeTemp_name"] + ", " + var_names["d_workspace_name"] + ", " + var_names["s_linalg_smem_name"] + ");"
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
                            str(self.gen_end_effector_pose_gradient_hessian_inner_temp_mem_size(include_d2_temp = False)) + \
                            " (the hot eeTemp/deeTemp chain; kept in smem at every tier)", \
                   "s_d2eeTemp is a pointer to helper memory of size " + \
                            str(self.gen_end_effector_pose_gradient_hessian_d2_temp_mem_size()) + \
                            " (the large n^2 Hessian arena). When !D2TEMP_IN_SMEM the inner repoints it at d_workspace", \
                   "d_workspace is the global spill arena s_d2eeTemp is repointed at when !D2TEMP_IN_SMEM (else unused)", \
                   "s_linalg_smem is optional byte-addressed shared memory for cuBLASDx"]
    func_notes = ["Assumes the Xhom and dXhom matricies have already been updated for the given q",
                  "Inner-owns scratch placement: the large n^2 s_d2eeTemp arena moves to d_workspace when !D2TEMP_IN_SMEM. The smaller eeTemp/deeTemp chain (s_temp) stays in smem always (it is hot)."]
    func_def_start = "void end_effector_pose_gradient_hessian_inner("
    func_def_middle = "T *s_d2eePos, T *s_deePos, const T *s_q, const T *s_Xhom, const T *s_dXhom, const T *s_d2Xhom, "
    func_def_end = "T *s_temp, T *s_d2eeTemp, T *d_workspace, unsigned char *s_linalg_smem) {"
    if use_thread_group:
        func_def_start += "cgrps::thread_group tgrp, "
        func_params.insert(0,"tgrp is the handle to the thread_group running this function")
    func_def_middle, func_params = self.gen_insert_helpers_func_def_params(func_def_middle, func_params, -1, NO_XI_FLAG = True)
    func_def = func_def_start + func_def_middle + func_def_end
    # now generate the code
    self.gen_add_func_doc("Computes the Gradient and Hessian of the End Effector Pose with respect to joint position",\
                          func_notes,func_params,None)
    self.gen_add_code_line("template <typename T, bool D2TEMP_IN_SMEM = true>")
    self.gen_add_code_line("__device__")
    self.gen_add_code_line(func_def, True)
    # Inner-controlled scratch placement: only the large n^2 s_d2eeTemp arena
    # moves to d_workspace when !D2TEMP_IN_SMEM. Reassigning it at the top keeps
    # every s_d2eeTemp[...] reference below unchanged. The hot eeTemp/deeTemp
    # chain (s_temp) stays in smem at every tier.
    self.gen_add_code_line("if constexpr (!D2TEMP_IN_SMEM) { s_d2eeTemp = d_workspace; } else { (void)d_workspace; }")
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
    if not self.robot.is_serial_chain():
        # NON-SERIAL / FLOATING-BASE: the legacy dense path computed ALL
        # (djid_i, djid_j, ee) triples per BFS level (quadratically wasteful) and
        # masked the out-of-chain ones. Emit ONE compacted in-chain chain-up
        # (gradient + hessian) via grid_linalg_indexed_batched_gemm instead.
        # Numerically identical: masked entries contributed exactly 0.
        _emit_eepose_hess_compacted_nonserial(self, n, all_ees, num_ees, use_thread_group)
        self.gen_add_end_function()
        return
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
    _emit_eepose_hess_extraction(self, n, num_ees, use_thread_group)
    self.gen_add_end_function()

def _emit_eepose_hess_extraction(self, n, num_ees, use_thread_group):
    # Shared d2eePos/deePos extraction. Assumes s_eeTemp (per-ee), s_deeTemp
    # (per-(ee,djid)) and s_d2eeTemp (per-(ee,i,j)) already point at their final
    # base (caller rebased). Out-of-chain deeTemp/d2eeTemp slots must be 0 so the
    # arctan2 derivative formulas evaluate to 0 there (matches the dense mask).
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
    self.gen_add_code_line("T dpitchSqrtTerm_i_top = s_Xmat_hom[10]*s_dXmat_hom_i[10] + s_Xmat_hom[6]*s_dXmat_hom_i[6];")
    self.gen_add_code_line("T dpitchSqrtTerm_i = dpitchSqrtTerm_i_top/pitchSqrtTerm;")
    self.gen_add_code_line("T dpitchSqrtTerm_j = (s_Xmat_hom[10]*s_dXmat_hom_j[10] + s_Xmat_hom[6]*s_dXmat_hom_j[6])/pitchSqrtTerm;")
    self.gen_add_code_line("T dpitchSqrtTerm_i_top_dj = s_dXmat_hom_j[10]*s_dXmat_hom_i[10] + s_Xmat_hom[10]*s_d2Xmat_hom_ij[10] + ")
    self.gen_add_code_line("                            s_dXmat_hom_j[6]*s_dXmat_hom_i[6] + s_Xmat_hom[6]*s_d2Xmat_hom_ij[6];")
    self.gen_add_code_line("// s'_i = T_i/s; s''_ij = d/dq_j(T_i/s) uses T_i (=_top), not s'_i, in the quotient numerator.")
    self.gen_add_code_line("T d2pitchSqrtTerm_ij = (pitchSqrtTerm*dpitchSqrtTerm_i_top_dj - dpitchSqrtTerm_i_top*dpitchSqrtTerm_j) / (pitchSqrtTerm*pitchSqrtTerm);")
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
    # s_d2eeTemp placement is now INNER-OWNED: the device only carves the smem
    # arena slot when this tier keeps it in smem (D2EE_D2TEMP_IN_SMEM<TIER>());
    # otherwise it hands the inner d_workspace and the per-tier D2TEMP_IN_SMEM
    # flag, and the inner repoints s_d2eeTemp at d_workspace itself.
    self.gen_XmatsHom_helpers_temp_shared_memory_code(inner_no_d2_size, include_gradients = True, include_hessians = True,
                                                      include_linalg_scratch = True,
                                                      linalg_scratch_bytes = "GRID_EE_LINALG_SHARED_BYTES<T>()")
    self.gen_add_code_line("T *s_d2eeTemp = nullptr;")
    self.gen_add_code_line("if constexpr (D2EE_D2TEMP_IN_SMEM<RESOURCE_TIER>()) {", True)
    self.gen_add_code_line("s_arena_offset = grid_align_up(s_arena_offset, alignof(T));")
    self.gen_add_code_line("s_d2eeTemp = grid_arena_ptr<T>(s_arena, s_arena_offset);")
    self.gen_add_code_line("s_arena_offset += sizeof(T) * static_cast<size_t>(" + str(d2_temp_size) + ");")
    self.gen_add_end_control_flow()
    # then load/update XI and run the algo
    self.gen_load_update_XmatsHom_helpers_function_call(use_thread_group, include_gradients = True, include_hessians = True)
    # Inner-owns placement: pass d_workspace + the per-tier flag. When the flag is
    # false the inner repoints s_d2eeTemp (currently nullptr) at d_workspace.
    self.gen_end_effector_pose_gradient_hessian_inner_function_call(use_thread_group,
        updated_var_names = {"d_workspace_name": "d_workspace"},
        d2temp_in_smem_expr = "D2EE_D2TEMP_IN_SMEM<RESOURCE_TIER>()")
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
    # Inner-owns s_d2eeTemp placement. When this tier keeps it in smem, carve it
    # from the tail of the (d2-inclusive) s_temp arena and pass the inner
    # d2temp_in_smem='true' with a null workspace. When this tier spills it, leave
    # s_d2eeTemp null and hand the inner the per-timestep workspace slice + 'false';
    # the inner repoints s_d2eeTemp at d_workspace itself.
    if not use_workspace_temp:
        self.gen_add_code_line("(void)d_workspace;")
        self.gen_add_code_line("T *s_d2eeTemp = &s_temp[" + str(self.gen_end_effector_pose_gradient_hessian_inner_temp_mem_size(include_d2_temp = False)) + "];")
    d2temp_in_smem_expr = "false" if use_workspace_temp else "true"
    # Per-tier d2eeTemp byte offset into the workspace slice. The shared
    # GRID_D2EE_WORKSPACE_D2EETEMP_OFFSET_BYTES macro keys off the single-valued
    # PERF-pick GRID_D2EE_USES_WORKSPACE_D2XHOM, so it would collide with the
    # spilled d2Xhom region at tiers whose pick spills d2Xhom but whose PERF pick
    # does not. Compute the offset locally from THIS tier's use_workspace_d2xhom
    # so d2eeTemp always lands past the (tier-local) spilled d2Xhom region.
    d2eetemp_off = "GRID_D2EE_WORKSPACE_TEMP_OFFSET_BYTES<T>()"
    if use_workspace_d2xhom:
        d2eetemp_off += " + sizeof(T) * static_cast<size_t>(D2XHOM_T_COUNT)"
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
            # Inner-owns: s_d2eeTemp stays null; the inner repoints it at this slice.
            self.gen_add_code_line("T *s_d2eeTemp = nullptr;")
            self.gen_add_code_line("T *s_d2eeTemp_ws = reinterpret_cast<T *>(&d_workspace[k*GRID_WORKSPACE_BYTES_PER_TIMESTEP<T>() + " + d2eetemp_off + "]);")
        self.gen_add_code_line("// compute")
        self.gen_load_update_XmatsHom_helpers_function_call(use_thread_group, include_gradients = True, include_hessians = True)
        updated = {"d_workspace_name": "s_d2eeTemp_ws"} if use_workspace_temp else None
        self.gen_end_effector_pose_gradient_hessian_inner_function_call(use_thread_group,
            updated_var_names = updated, d2temp_in_smem_expr = d2temp_in_smem_expr)
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
            # Inner-owns: s_d2eeTemp stays null; the inner repoints it at this slice.
            self.gen_add_code_line("T *s_d2eeTemp = nullptr;")
            self.gen_add_code_line("T *s_d2eeTemp_ws = reinterpret_cast<T *>(&d_workspace[" + d2eetemp_off + "]);")
        self.gen_add_code_line("// compute with NUM_TIMESTEPS as NUM_REPS for timing")
        self.gen_add_code_line("for (int rep = 0; rep < NUM_TIMESTEPS; rep++){", True)
        self.gen_anti_licm_input_reload("q",str(n),use_thread_group,feedback_from="d2eePos")
        self.gen_load_update_XmatsHom_helpers_function_call(use_thread_group, include_gradients = True, include_hessians = True)
        updated = {"d_workspace_name": "s_d2eeTemp_ws"} if use_workspace_temp else None
        self.gen_end_effector_pose_gradient_hessian_inner_function_call(use_thread_group,
            updated_var_names = updated, d2temp_in_smem_expr = d2temp_in_smem_expr)
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
    # Per-tier gate: arm L2 persistence if ANY tier routes the d2eeTemp arena
    # through d_workspace (runtime RESOURCE_TIER may differ from the PERF pick).
    self.gen_add_code_line("if (GRID_D2EE_USES_WORKSPACE_TEMP_ANY) {gpuErrchk(grid_begin_l2_persisting(0, hd_data->d_workspace, " + workspace_bytes + "));}")
    self.gen_add_code_lines(func_call_code)
    self.gen_add_code_line("if (GRID_D2EE_USES_WORKSPACE_TEMP_ANY) {gpuErrchk(grid_end_l2_persisting(0));}")
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
        self.gen_add_code_line(single_call_printf_line("ee_pose_hessian"))
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
