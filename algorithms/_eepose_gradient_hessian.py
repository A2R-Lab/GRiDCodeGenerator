"""
End Effector Posiitons

TODO: fix throughout this document for fixed_joint support for branched trees and for multiple fixed at once
"""
def gen_end_effector_pose_inner_temp_mem_size(self, fixed_target_name = ""):
    num_ees = self.robot.get_total_leaf_nodes() if fixed_target_name == "" else 1
    return 2*16*num_ees

def gen_end_effector_pose_inner_function_call(self, updated_var_names = None, fixed_target_name = "",
                                              temp_in_smem_expr = "true"):
    var_names = dict( \
        s_Xhom_name = "s_XmatsHom", \
        s_end_effector_pose_name = "s_end_effector_pose", \
        s_q_name = "s_q", \
        s_topology_helpers_name = "s_topology_helpers", \
        s_temp_name = "s_temp", \
        d_workspace_name = "nullptr", \
        s_linalg_smem_name = "s_linalg_smem", \
    )
    if updated_var_names is not None:
        for key,value in updated_var_names.items():
            var_names[key] = value
    code_start = "end_effector_pose_inner" + ("" if fixed_target_name == "" else "_" + fixed_target_name) + "<T, " + temp_in_smem_expr + ">(" + var_names["s_end_effector_pose_name"] + ", " + var_names["s_q_name"] + ", "
    code_middle = var_names["s_Xhom_name"] + ", "
    code_end =  var_names["s_temp_name"] + ", " + var_names["d_workspace_name"] + ", " + var_names["s_linalg_smem_name"] + ");"
    # account for thread group and serial chains
    # Canonical: append the shared topology-helper arg via the central helper
    # (NO_XI: the ee_pose family takes s_Xhom, not s_XImats). Mirrors the def's
    # gen_insert_helpers_func_def_params(NO_XI_FLAG=True) so def + call can't drift.
    code_middle += self.gen_insert_helpers_function_call(updated_var_names = var_names, NO_XI_FLAG = True)
    self.gen_add_code_line(code_start + code_middle + code_end)

def gen_end_effector_pose_inner(self, fixed_target_name = ""):
    n = self.robot.get_num_pos()
    n_bfs_levels = self.robot.get_max_bfs_level() + 1 # starts at 0
    if fixed_target_name == "":
        all_ees = self.robot.get_leaf_nodes()
    else:
        all_ees = [self.robot.get_fixed_joint_by_name(fixed_target_name).get_id()]
    num_ees = len(all_ees)
    # construct the boilerplate and function definition
    func_params = ["s_end_effector_pose is a pointer to shared memory of size 6*NUM_EE where NUM_EE = " + str(num_ees), \
                   "s_q is the vector of joint positions", \
                   "s_Xhom is the pointer to the homogenous transformation matricies ", \
                   "s_temp is a pointer to helper shared memory of size " + \
                            str(self.gen_end_effector_pose_inner_temp_mem_size(fixed_target_name)), \
                   "d_workspace is the global-memory chain workspace used in place of s_temp when !TEMP_IN_SMEM", \
                   "s_linalg_smem is optional byte-addressed shared memory (reserved; unused by this inner)"]
    func_notes = ["Assumes the Xhom matricies have already been updated for the given q", "Defaults to all leave nodes if fixed_target_name is not provided"]
    func_def_start = "void end_effector_pose_inner" + ("" if fixed_target_name == "" else "_" + fixed_target_name) + "("
    func_def_middle = "T *s_end_effector_pose, const T *s_q, const T *s_Xhom, "
    func_def_end = "T *s_temp, T *d_workspace, unsigned char *s_linalg_smem) {"
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
        self.gen_add_sync()
        self.gen_add_serial_ops()
        self.gen_add_code_line("printf(\"q\\n\"); printMat<T,1," + str(n) + ">(s_q,1);")
        self.gen_add_code_line("for (int i = 0; i < " + str(n) + "; i++){printf(\"X[%d]\\n\",i); printMat<T,4,4>(&s_Xhom[16*i],4);}")
        self.gen_add_end_control_flow()
        self.gen_add_sync()

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
                self.gen_add_parallel_loop("ind",str(16))
                self.gen_add_code_line("s_temp[ind] = s_Xhom[16*" + str(all_ees[0]) + " + ind];")
                self.gen_add_end_control_flow()
                self.gen_add_sync()
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
                self.gen_add_parallel_loop("ind",str(16))
                self.gen_add_code_line("int row = ind % 4; int col = ind / 4;")
                self.gen_add_code_line("s_temp[ind + " + str(tempDstOffset) + "] = dot_prod<T,4,4,1>" + \
                                       "(&s_Xhom[16*" + str(parent) + " + row], &s_temp[" + str(tempSrcOffset) + " + 4*col]);")
                self.gen_add_end_control_flow()
                self.gen_add_sync()
                # update parent for next loop (if there is one)
                parent = self.robot.get_parent_id(parent)
        else:
            # if first loop then just set to transform at the leaf
            if bfs_level == 0:
                self.gen_add_code_line("// First set to leaf transform")
                self.gen_add_parallel_loop("ind",str(16*num_ees))
                self.gen_add_code_line("int rc = ind % 16;")
                select_var_vals = [("int", "eeInd", [str(jid) for jid in all_ees])]
                self.gen_add_multi_threaded_select("ind", "<", [str(16*(i+1)) for i in range(num_ees)], select_var_vals)
                self.gen_add_code_line("s_temp[ind] = s_Xhom[16*eeInd + rc];")
                self.gen_add_end_control_flow()
                self.gen_add_sync()
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
                self.gen_add_parallel_loop("ind",str(16*num_ees))
                self.gen_add_code_line("int row = ind % 4; int col = (ind / 4) % 4; int eeOffset = ind - (ind % 16);")
                # get parents for this level
                select_var_vals = [("int", "parent_jid", [str(jid) for jid in curr_parents])]
                self.gen_add_multi_threaded_select("ind", "<", [str(16*(i+1)) for i in range(num_ees)], select_var_vals)
                if (-1 in curr_parents):
                    self.gen_add_code_line("if(parent_jid == -1){continue;}")
                self.gen_add_code_line("s_temp[ind + " + str(tempDstOffset) + "] = dot_prod<T,4,4,1>" + \
                                       "(&s_Xhom[16*parent_jid + row], &s_temp[" + str(tempSrcOffset) + " + eeOffset + 4*col]);")
                self.gen_add_end_control_flow()
                self.gen_add_sync()
    
    self.gen_add_code_line("//")
    self.gen_add_code_line("// Now extract the end_effector_pose from the Tansforms")
    self.gen_add_code_line("// TODO: ADD OFFSETS")
    self.gen_add_code_line("//")
    tempOffset = 16*num_ees*(bfs_level % 2)
    # xyz position is easy (end_effector_pose_xyz1 = Xmat_hom * offset) where offset = [x,y,z,1]
    self.gen_add_parallel_loop("ind",str(3*num_ees))
    self.gen_add_code_line("// xyz is easy")
    self.gen_add_code_line("int xyzInd = ind % 3; int eeInd = ind / 3; T *s_Xmat_hom = &s_temp[" + str(tempOffset) + " + 16*eeInd];")
    self.gen_add_code_line("s_end_effector_pose[6*eeInd + xyzInd] = s_Xmat_hom[12 + xyzInd];")
    # roll pitch yaw is a bit more difficult
    self.gen_add_code_line("// roll pitch yaw is a bit more difficult")
    self.gen_add_code_line("if(xyzInd > 0){continue;}")
    self.gen_add_code_line("s_end_effector_pose[6*eeInd + 3] = atan2(s_Xmat_hom[6],s_Xmat_hom[10]);")
    self.gen_add_code_line("s_end_effector_pose[6*eeInd + 4] = -atan2(s_Xmat_hom[2],sqrt(s_Xmat_hom[6]*s_Xmat_hom[6] + s_Xmat_hom[10]*s_Xmat_hom[10]));")
    self.gen_add_code_line("s_end_effector_pose[6*eeInd + 5] = atan2(s_Xmat_hom[1],s_Xmat_hom[0]);")
    self.gen_add_end_control_flow()
    self.gen_add_sync()
    self.gen_add_end_function()

def gen_end_effector_pose_device_temp_mem_size(self, fixed_target_name = ""):
    n = self.robot.get_num_pos()
    XHom_size, dXhom_size, d2Xhom_size = self.gen_get_Xhom_size()
    wrapper_size = self.gen_topology_helpers_size() + XHom_size # for Xhom
    return self.gen_end_effector_pose_inner_temp_mem_size(fixed_target_name) + wrapper_size

def gen_end_effector_pose_device(self, fixed_target_name = ""):
    n = self.robot.get_num_pos()
    num_ees = self.robot.get_total_leaf_nodes() if fixed_target_name == "" else 1
    # construct the boilerplate and function definition
    func_params = ["s_end_effector_pose is a pointer to shared memory of size 6*NUM_EE where NUM_EE = " + str(num_ees), \
                   "s_q is the vector of joint positions", \
                   "d_robotModel is the pointer to the initialized model specific helpers on the GPU (XImats, topology_helpers, etc.)"]
    func_notes = []
    func_def_start = "void end_effector_pose_device" + ("" if fixed_target_name == "" else "_" + fixed_target_name) + "("
    func_def_middle = "T *s_end_effector_pose, const T *s_q, "
    func_def_end = "const robotModel<T> *d_robotModel) {"
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
    self.gen_load_update_XmatsHom_helpers_function_call()
    self.gen_end_effector_pose_inner_function_call(fixed_target_name = fixed_target_name)
    self.gen_add_end_function()

def gen_end_effector_pose_kernel(self, single_call_timing = False, fixed_target_name = ""):
    n = self.robot.get_num_pos()
    num_ees = self.robot.get_total_leaf_nodes() if fixed_target_name == "" else 1
    # define function def and params
    func_params = ["d_end_effector_pose is the vector of end effector positions", \
                   "d_q is the vector of joint positions", \
                   "stride_q is the stide between each q", \
                   "d_robotModel is the pointer to the initialized model specific helpers on the GPU (XImats, topology_helpers, etc.)", \
                   "num_timesteps is the length of the trajectory points we need to compute over (or overloaded as test_iters for timing)"]
    func_notes = []
    func_def_start = "void end_effector_pose_kernel" + ("" if fixed_target_name == "" else "_" + fixed_target_name) + "(T *d_end_effector_pose, const T *d_q, const int stride_q, "
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
    self.gen_XmatsHom_helpers_temp_shared_memory_code(shared_mem_size, extra_t_buffers = [("s_q", n), ("s_end_effector_pose", 6*num_ees)],
                                                      include_linalg_scratch = True,
                                                      linalg_scratch_bytes = "GRID_EE_LINALG_SHARED_BYTES<T>()")
    if not single_call_timing:
        # load to shared mem and loop over blocks to compute all requested comps
        self.gen_add_parallel_loop("k","NUM_TIMESTEPS",block_level = True)
        self.gen_kernel_load_inputs("q",str(n),stride="stride_q")
        # compute
        self.gen_add_code_line("// compute")
        # then load/update X and run the algo
        self.gen_load_update_XmatsHom_helpers_function_call()
        self.gen_end_effector_pose_inner_function_call(fixed_target_name = fixed_target_name)
        self.gen_add_sync()
        # save to global
        self.gen_kernel_save_result("end_effector_pose",str(6*num_ees),stride=str(6*num_ees))
        self.gen_add_end_control_flow()
    else:
        #repurpose NUM_TIMESTEPS for number of timing reps
        self.gen_kernel_load_inputs("q",str(n))
        # then compute in loop for timing
        self.gen_add_code_line("// compute with NUM_TIMESTEPS as NUM_REPS for timing")
        self.gen_add_code_line("for (int rep = 0; rep < NUM_TIMESTEPS; rep++){", True)
        self.gen_anti_licm_input_reload("q",str(n),feedback_from="end_effector_pose")
        # then load/update X and run the algo
        self.gen_load_update_XmatsHom_helpers_function_call()
        self.gen_end_effector_pose_inner_function_call(fixed_target_name = fixed_target_name)
        self.gen_anti_licm_output_write("end_effector_pose")
        self.gen_add_end_control_flow()
        # save to global
        self.gen_kernel_save_result("end_effector_pose",str(6*num_ees))
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
                        "<T><<<block_dimms,thread_dimms,END_EFFECTOR_POSE_DYNAMIC_SHARED_MEM_BYTES<T>()>>>(hd_data->d_end_effector_pose,hd_data->d_q,stride_q,"
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
    self.gen_add_code_line("gpuErrchk(grid_check_dynamic_shared_memory_bytes(\"end_effector_pose\", END_EFFECTOR_POSE_DYNAMIC_SHARED_MEM_BYTES<T>()));")
    self.gen_add_code_lines(func_call_code)
    if not compute_only:
        # then transfer memory back
        self.gen_add_code_lines(["// finally transfer the result back", \
                                 "gpuErrchk(cudaMemcpy(hd_data->h_end_effector_pose,hd_data->d_end_effector_pose,6*NUM_EES*" + \
                                    ("num_timesteps*" if not single_call_timing else "") + "sizeof(T),cudaMemcpyDeviceToHost));",
                                 "gpuErrchkKernel();"])
    # finally report out timing if requested
    if single_call_timing:
        from ..algo_registry import single_call_printf_line
        self.gen_add_code_line(single_call_printf_line("end_effector_pose"))
    self.gen_add_end_function()

def gen_end_effector_pose_gradient_inner_temp_mem_size(self, fixed_target_name = ""):
    # Scratch for the shared-chain geometric Jacobian:
    #   s_Xworld  : 16 * NUM_JOINTS   (world transforms of every joint)
    #   s_Jv,s_Jw : 2 * (3 * nv * num_ees)
    #   s_E       : 4 * num_ees       (cy, sy, cp, sp per ee for E(rpy) inversion)
    n_joints = self.robot.get_num_joints()
    nv = self.robot.get_num_vel()
    num_ees = self.robot.get_total_leaf_nodes() if fixed_target_name == "" else 1
    return 16*n_joints + 2*3*nv*num_ees + 4*num_ees

def gen_end_effector_pose_gradient_inner_function_call(self, updated_var_names = None, fixed_target_name = "",
                                                       temp_in_smem_expr = "true"):
    var_names = dict( \
        s_Xhom_name = "s_XmatsHom", \
        s_dXhom_name = "s_dXmatsHom", \
        s_end_effector_pose_gradient_name = "s_end_effector_pose_gradient", \
        s_q_name = "s_q", \
        s_topology_helpers_name = "s_topology_helpers", \
        s_temp_name = "s_temp", \
        d_workspace_name = "nullptr", \
        s_linalg_smem_name = "s_linalg_smem", \
    )
    if updated_var_names is not None:
        for key,value in updated_var_names.items():
            var_names[key] = value
    code_start = "end_effector_pose_gradient_inner" + ("" if fixed_target_name == "" else "_" + fixed_target_name) + "<T, " + temp_in_smem_expr + ">(" + var_names["s_end_effector_pose_gradient_name"] + ", " + var_names["s_q_name"] + ", "
    code_middle = var_names["s_Xhom_name"] + ", " + var_names["s_dXhom_name"] + ", "
    code_end =  var_names["s_temp_name"] + ", " + var_names["d_workspace_name"] + ", " + var_names["s_linalg_smem_name"] + ");"
    # account for thread group
    # Canonical: append the shared topology-helper arg via the central helper
    # (NO_XI: the ee_pose family takes s_Xhom, not s_XImats). Mirrors the def's
    # gen_insert_helpers_func_def_params(NO_XI_FLAG=True) so def + call can't drift.
    code_middle += self.gen_insert_helpers_function_call(updated_var_names = var_names, NO_XI_FLAG = True)
    self.gen_add_code_line(code_start + code_middle + code_end)

def _eepose_grad_chain_metadata(self, all_ees, fixed_target_name):
    """Bake out per-ee chain-joint fill jobs for the geometric-Jacobian rewrite.

    For each end-effector returns:
      (chain_jids, ee_anchor_jid, ee_uses_fixed_offset_chain)
      jobs: list of dicts { j, vi, ang_local (len-3), lin_local (len-3), revolute }

    `ee_anchor_jid` is the joint whose world transform is the EE's world frame
    (for a leaf-joint EE this is the leaf itself; for a fixed-joint EE this is
    the parent joint and the fixed transform is composed in C++ separately —
    not implemented in this first cut and is asserted out)."""
    import numpy as _np
    chains, anchors, jobs_all = [], [], []
    for ee in all_ees:
        chain = sorted(self.robot.get_ancestors_by_id(ee)) + [ee]
        chains.append(chain)
        anchors.append(ee)
        jobs = []
        for j in chain:
            S = _np.asarray(self.robot.get_S_by_id(j), dtype=_np.float64)
            if S.ndim == 1:
                S = S.reshape(-1, 1)
            try:
                vinds = self.robot.get_joint_index_v(j)
            except Exception:
                vinds = self.robot.get_joint_index_q(j)
            if not isinstance(vinds, (list, tuple, _np.ndarray)):
                vinds = [vinds]
            vinds = list(vinds)
            for c in range(S.shape[1]):
                vi = vinds[c] if c < len(vinds) else vinds[-1]
                ang_local = [float(x) for x in S[:3, c]]
                lin_local = [float(x) for x in S[3:6, c]]
                revolute = max(abs(x) for x in ang_local) > 0.5
                jobs.append({
                    "j": int(j), "vi": int(vi),
                    "ang": ang_local, "lin": lin_local,
                    "revolute": bool(revolute),
                })
        jobs_all.append(jobs)
    return chains, anchors, jobs_all

def gen_end_effector_pose_gradient_inner(self, fixed_target_name = ""):
    """Shared-chain geometric (spatial) Jacobian for d(pose)/dv (tangent).

    Output `s_end_effector_pose_gradient` is sized 6 * nv * NUM_EE (NOT 6 * nq) so the floating-base
    base block is the spatial Jacobian (omega; v_world) rather than the older
    non-standard quaternion-component derivs. Convention matches pinocchio's
    LOCAL_WORLD_ALIGNED frame Jacobian (mapped through E(rpy)^{-1} for the rpy
    rows). Algorithm:
      1. One forward-kinematics pass builds the world transform of every joint
         via BFS-level chain-up (s_Xworld).
      2. Per ee, per chain joint j, per S column c: J_w[:, ee, vi] = R_j_world *
         ang_local (revolute) or 0 (prismatic); J_v[:, ee, vi] = J_w x (p_ee -
         p_j) (revolute) or R_j_world * lin_local (prismatic). vi is the
         joint's velocity index from get_joint_index_v.
      3. Per ee, extract (cy, sy, cp, sp) from R_ee_world.
      4. Write s_end_effector_pose_gradient: rows 0..2 = J_v columns; rows 3..5 = E(rpy)^{-1} * J_w
         columns. Closed-form E^{-1} avoids an explicit matrix inverse:
            row 3 (droll/dv): (cy*Jw[0] + sy*Jw[1]) / cp
            row 4 (dpitch/dv): -sy*Jw[0] + cy*Jw[1]
            row 5 (dyaw/dv):  sp/cp * (cy*Jw[0] - sy*Jw[1]) + Jw[2]
    """
    n = self.robot.get_num_pos()
    nv = self.robot.get_num_vel()
    n_joints = self.robot.get_num_joints()
    n_bfs_levels = self.robot.get_max_bfs_level() + 1

    if fixed_target_name == "":
        all_ees = self.robot.get_leaf_nodes()
    else:
        # The fixed-target gradient path predates this rewrite and is not
        # exercised by the current bench/equivalence harness; flag if it
        # ever surfaces so we know to extend the shared-chain emission.
        raise NotImplementedError(
            "gen_end_effector_pose_gradient_inner: fixed_target_name='" + fixed_target_name +
            "' not yet supported by the shared-chain geometric-Jacobian rewrite."
        )
    num_ees = len(all_ees)
    chains, anchors, fill_jobs = _eepose_grad_chain_metadata(self, all_ees, fixed_target_name)

    # function header
    func_params = ["s_end_effector_pose_gradient is a pointer to shared memory of size 6*NUM_VEL*NUM_EE where NUM_VEL = " + str(nv) + " and NUM_EE = " + str(num_ees), \
                   "s_q is the vector of joint positions (unused; kept for signature compatibility)", \
                   "s_Xhom is the pointer to the LOCAL homogeneous transformation matrices (per-joint Xhom_local)", \
                   "s_dXhom is the pointer to the LOCAL d-transforms (unused by the geometric-Jacobian path; kept for signature compatibility)", \
                   "s_temp is a pointer to helper shared memory of size " + \
                            str(self.gen_end_effector_pose_gradient_inner_temp_mem_size()), \
                   "d_workspace is the global-memory chain workspace used in place of s_temp when !TEMP_IN_SMEM", \
                   "s_linalg_smem is optional byte-addressed shared memory (reserved; unused)"]
    func_notes = ["Assumes s_Xhom has been populated with the per-joint LOCAL transforms for the given q.",
                  "Output d/dv (TANGENT) is 6 x nv per ee (was 6 x nq for d/dq) -- matches pinocchio."]
    func_def_start = "void end_effector_pose_gradient_inner" + ("" if fixed_target_name == "" else "_" + fixed_target_name) + "("
    func_def_middle = "T *s_end_effector_pose_gradient, const T *s_q, const T *s_Xhom, const T *s_dXhom, "
    func_def_end = "T *s_temp, T *d_workspace, unsigned char *s_linalg_smem) {"
    func_def_middle, func_params = self.gen_insert_helpers_func_def_params(func_def_middle, func_params, -1, NO_XI_FLAG = True)
    func_def = func_def_start + func_def_middle + func_def_end
    self.gen_add_func_doc("Computes the Gradient of the End Effector Pose with respect to generalized velocity (d/dv tangent, pinocchio convention)",
                          func_notes,func_params,None)
    self.gen_add_code_line("template <typename T, bool TEMP_IN_SMEM = true>")
    self.gen_add_code_line("__device__")
    self.gen_add_code_line(func_def, True)
    self.gen_add_code_line("if constexpr (!TEMP_IN_SMEM) { s_temp = d_workspace; } else { (void)d_workspace; }")
    self.gen_add_code_line("(void)s_q; (void)s_dXhom; (void)s_linalg_smem;")

    # scratch layout (matches gen_end_effector_pose_gradient_inner_temp_mem_size)
    off_Xworld = 0
    off_Jv = off_Xworld + 16 * n_joints
    off_Jw = off_Jv + 3 * nv * num_ees
    off_E  = off_Jw + 3 * nv * num_ees   # 4 * num_ees floats: cy, sy, cp, sp per ee
    self.gen_add_code_line("// scratch layout: Xworld | Jv (3 x nv x ee) | Jw (3 x nv x ee) | E_sincos (4 x ee)")
    self.gen_add_code_line("T *s_Xworld = &s_temp[" + str(off_Xworld) + "];")
    self.gen_add_code_line("T *s_Jv     = &s_temp[" + str(off_Jv)     + "];")
    self.gen_add_code_line("T *s_Jw     = &s_temp[" + str(off_Jw)     + "];")
    self.gen_add_code_line("T *s_E_sc   = &s_temp[" + str(off_E)      + "];   // cy,sy,cp,sp per ee")

    # ============ Step 1: world transforms by BFS level ============
    self.gen_add_code_line("//")
    self.gen_add_code_line("// Step 1: build world transforms for every joint via BFS-level chain-up")
    self.gen_add_code_line("//")
    for level in range(n_bfs_levels):
        ids_at_level = self.robot.get_ids_by_bfs_level(level)
        if not ids_at_level:
            continue
        njs = len(ids_at_level)
        self.gen_add_code_line("// BFS level " + str(level) + " -> joints " + str(ids_at_level))
        self.gen_add_parallel_loop("ind", str(16 * njs))
        self.gen_add_code_line("int slot = ind / 16; int ele = ind % 16;")
        self.gen_add_code_line("int row = ele & 3; int col = ele >> 2;")
        # bake the joint id and parent id per slot
        jid_list = [str(j) for j in ids_at_level]
        par_list = [str(self.robot.get_parent_id(j)) for j in ids_at_level]
        select_var_vals = [("int", "jid", jid_list), ("int", "par", par_list)]
        self.gen_add_multi_threaded_select("slot", "<", [str(i+1) for i in range(njs)], select_var_vals)
        # If par == -1 (root), world := local; else world[jid] = world[par] @ local[jid].
        # local[jid] is s_Xhom[16*jid]; world[jid] is s_Xworld[16*jid].
        self.gen_add_code_line("if (par == -1) {", True)
        self.gen_add_code_line("s_Xworld[16*jid + ele] = s_Xhom[16*jid + ele];")
        self.gen_add_end_control_flow()
        self.gen_add_code_line("else {", True)
        # dot_prod<T,4,4,1>(&s_Xworld[16*par + row], &s_Xhom[16*jid + 4*col])
        self.gen_add_code_line("s_Xworld[16*jid + ele] = dot_prod<T,4,4,1>(&s_Xworld[16*par + row], &s_Xhom[16*jid + 4*col]);")
        self.gen_add_end_control_flow()
        self.gen_add_end_control_flow()
        self.gen_add_sync()

    # ============ Step 2: zero Jv, Jw ============
    self.gen_add_code_line("//")
    self.gen_add_code_line("// Step 2: zero the J_v and J_w scratch (out-of-chain columns stay zero)")
    self.gen_add_code_line("//")
    self.gen_add_parallel_loop("ind", str(2 * 3 * nv * num_ees))
    self.gen_add_code_line("s_Jv[ind] = static_cast<T>(0);")
    self.gen_add_end_control_flow()
    self.gen_add_sync()

    # ============ Step 3: per-ee, per-chain-joint, per-S-col column fills ============
    # Each (ee, vi) pair emits one block that:
    #   - reads R_j_world (3x3 in column-major from s_Xworld[16*j+0..10])
    #   - reads p_j_world (s_Xworld[16*j + 12..14])
    #   - reads p_ee_world (s_Xworld[16*ee_anchor + 12..14])
    #   - computes axis_world (3-vector via R_j @ S_local)
    #   - writes Jv, Jw columns
    # Flatten all (ee, job) pairs across ees so each work item is one column-fill.
    flat_jobs = []
    for ee_idx, jobs in enumerate(fill_jobs):
        ee_anchor = anchors[ee_idx]
        for job in jobs:
            flat_jobs.append((ee_idx, ee_anchor, job))
    # Collapse jobs that target the same destination column (same (ee_idx, vi)) —
    # these occur when several joint ids in a chain map to one velocity coordinate.
    # For a NON-mimic robot this only happens degenerately and the legacy emission
    # let the LAST writer win. For a MIMIC robot a mimic joint SHARES its target's
    # v-slot and BOTH contribute to that geometric-Jacobian column, each scaled by
    # its mimic multiplier alpha (the mimic body moves alpha * the target's joint
    # rate). So group jobs by (ee_idx, vi): single-job groups go to the disjoint
    # block-parallel fill (byte-identical to the legacy path for non-mimic, where
    # every group is a singleton with alpha == 1); multi-job groups (mimic) are
    # emitted as a serial alpha-accumulate fold so the shared column sums all
    # contributions. Mirrors the inverse_dynamics_gradient / crba mimic v-slot accumulate.
    HAS_MIMIC = self.robot_has_mimic_joints()
    _groups = {}
    for entry in flat_jobs:
        ee_idx, _anc, job = entry
        _groups.setdefault((ee_idx, job["vi"]), []).append(entry)
    single_jobs = [grp[0] for grp in _groups.values() if len(grp) == 1]
    multi_groups = [grp for grp in _groups.values() if len(grp) > 1]
    flat_jobs = single_jobs
    n_flat = len(flat_jobs)
    if n_flat > 0:
        self.gen_add_code_line("//")
        self.gen_add_code_line("// Step 3: per-chain-joint columns of J_v, J_w (one block-parallel work-item per (ee, S-column))")
        self.gen_add_code_line("//")
        # Each job fills a DISJOINT column of J_v / J_w (the destination base
        # offset 3*nv*ee_idx + 3*vi is unique per (ee_idx, vi)), so the jobs are
        # fully independent and distribute across the block with no atomics. The
        # per-job compile-time constants (joint id j, S-column vi, ee anchor, the
        # 3-component local axis, the revolute flag) are baked into const arrays
        # indexed by job_idx; the existing Step-2 zero-fill + sync above leaves
        # out-of-chain columns at zero.
        job_j    = [job["j"] for (_ee, _anc, job) in flat_jobs]
        job_anc  = [ee_anchor for (_ee, ee_anchor, _job) in flat_jobs]
        job_rev  = [1 if job["revolute"] else 0 for (_ee, _anc, job) in flat_jobs]
        # Destination column base into s_Jv / s_Jw is unique per (ee_idx, vi);
        # bake it (with nv folded in) at codegen time so the device code needs no
        # runtime nv symbol and the disjointness is manifest.
        job_base = [3*nv*ee_idx + 3*job["vi"] for (ee_idx, _anc, job) in flat_jobs]
        job_ax   = []  # the local axis (angular for revolute, linear for prismatic)
        for (_ee, _anc, job) in flat_jobs:
            ax = job["ang"] if job["revolute"] else job["lin"]
            # Snap sub-threshold components to exact 0 so axw = a*x + 0*0 + 0*0 is
            # bit-identical to the original "drop near-zero terms" emission (adding
            # exact 0.0 never perturbs a finite float).
            job_ax.append([float(ax[c]) if abs(ax[c]) >= 1e-15 else 0.0 for c in range(3)])

        def _int_arr(vals):
            return "{ " + ", ".join(str(v) for v in vals) + " }"
        def _ax_arr(vals):
            return "{ " + ", ".join("static_cast<T>({:.17g})".format(v) for v in vals) + " }"

        self.gen_add_code_line("static const int eeg_job_j[]    = " + _int_arr(job_j) + ";")
        self.gen_add_code_line("static const int eeg_job_anc[]  = " + _int_arr(job_anc) + ";")
        self.gen_add_code_line("static const int eeg_job_rev[]  = " + _int_arr(job_rev) + ";")
        self.gen_add_code_line("static const int eeg_job_base[] = " + _int_arr(job_base) + ";")
        self.gen_add_code_line("const T eeg_job_ax[] = " + _ax_arr([a for ax in job_ax for a in ax]) + ";")
        self.gen_add_parallel_loop("job_idx", str(n_flat))
        self.gen_add_code_line("int j   = eeg_job_j[job_idx];")
        self.gen_add_code_line("int ee_anchor = eeg_job_anc[job_idx];")
        self.gen_add_code_line("int col_base = eeg_job_base[job_idx];")
        self.gen_add_code_line("T ax0 = eeg_job_ax[3*job_idx + 0]; T ax1 = eeg_job_ax[3*job_idx + 1]; T ax2 = eeg_job_ax[3*job_idx + 2];")
        # column-major rotation: R_j[r,c] = s_Xworld[16*j + r + 4*c], r,c in 0..2
        # axis_world[r] = sum_c R_j[r,c] * ax[c]
        self.gen_add_code_line("T axw_0 = s_Xworld[16*j + 0]*ax0 + s_Xworld[16*j + 4]*ax1 + s_Xworld[16*j + 8]*ax2;")
        self.gen_add_code_line("T axw_1 = s_Xworld[16*j + 1]*ax0 + s_Xworld[16*j + 5]*ax1 + s_Xworld[16*j + 9]*ax2;")
        self.gen_add_code_line("T axw_2 = s_Xworld[16*j + 2]*ax0 + s_Xworld[16*j + 6]*ax1 + s_Xworld[16*j + 10]*ax2;")
        self.gen_add_code_line("if (eeg_job_rev[job_idx]) {", True)
        # J_w[ee, vi, r] = axw_r
        self.gen_add_code_line("s_Jw[col_base + 0] = axw_0; s_Jw[col_base + 1] = axw_1; s_Jw[col_base + 2] = axw_2;")
        # arm = p_ee - p_j -> dx, dy, dz ; J_v = axw cross (p_ee - p_j)
        self.gen_add_code_line("T dx = s_Xworld[16*ee_anchor + 12] - s_Xworld[16*j + 12];")
        self.gen_add_code_line("T dy = s_Xworld[16*ee_anchor + 13] - s_Xworld[16*j + 13];")
        self.gen_add_code_line("T dz = s_Xworld[16*ee_anchor + 14] - s_Xworld[16*j + 14];")
        self.gen_add_code_line("s_Jv[col_base + 0] = axw_1*dz - axw_2*dy;")
        self.gen_add_code_line("s_Jv[col_base + 1] = axw_2*dx - axw_0*dz;")
        self.gen_add_code_line("s_Jv[col_base + 2] = axw_0*dy - axw_1*dx;")
        self.gen_add_end_control_flow()
        self.gen_add_code_line("else {", True)
        # J_v[ee, vi, r] = axw_r; J_w already zero from init
        self.gen_add_code_line("s_Jv[col_base + 0] = axw_0; s_Jv[col_base + 1] = axw_1; s_Jv[col_base + 2] = axw_2;")
        self.gen_add_end_control_flow()
        self.gen_add_end_control_flow()
        self.gen_add_sync()

    # ============ Step 3b: MIMIC shared-column alpha-accumulate (serial) ===========
    # A mimic joint and its target share one velocity coordinate; both contribute
    # to that geometric-Jacobian column scaled by their mimic multiplier alpha.
    # Emit these shared columns serially (thread 0), accumulating alpha_j * J_col(j)
    # over every chain joint j in the group. Only reached for mimic robots (non-mimic
    # has no multi-job groups), so non-mimic output is byte-identical.
    if HAS_MIMIC and multi_groups:
        self.gen_add_code_line("//")
        self.gen_add_code_line("// Step 3b: mimic shared-v-slot columns (serial alpha-accumulate)")
        self.gen_add_code_line("//")
        self.gen_add_serial_ops()
        for grp in multi_groups:
            ee_idx0, _anc0, job0 = grp[0]
            col_base = 3 * nv * ee_idx0 + 3 * job0["vi"]
            self.gen_add_code_line("// (ee " + str(ee_idx0) + ", v-slot " + str(job0["vi"]) +
                                   ") <- " + str(len(grp)) + " chain joints")
            self.gen_add_code_line("s_Jv[" + str(col_base) + " + 0] = static_cast<T>(0); s_Jv[" + str(col_base) + " + 1] = static_cast<T>(0); s_Jv[" + str(col_base) + " + 2] = static_cast<T>(0);")
            self.gen_add_code_line("s_Jw[" + str(col_base) + " + 0] = static_cast<T>(0); s_Jw[" + str(col_base) + " + 1] = static_cast<T>(0); s_Jw[" + str(col_base) + " + 2] = static_cast<T>(0);")
            for (ee_idx, ee_anchor, job) in grp:
                j = job["j"]
                alpha = self._alpha_for_jid(j)
                ax = job["ang"] if job["revolute"] else job["lin"]
                ax = [float(ax[c]) if abs(ax[c]) >= 1e-15 else 0.0 for c in range(3)]
                self.gen_add_code_line("{")
                self.gen_add_code_line("  T ax0 = static_cast<T>({:.17g}); T ax1 = static_cast<T>({:.17g}); T ax2 = static_cast<T>({:.17g});".format(ax[0], ax[1], ax[2]))
                self.gen_add_code_line("  T axw_0 = s_Xworld[16*" + str(j) + " + 0]*ax0 + s_Xworld[16*" + str(j) + " + 4]*ax1 + s_Xworld[16*" + str(j) + " + 8]*ax2;")
                self.gen_add_code_line("  T axw_1 = s_Xworld[16*" + str(j) + " + 1]*ax0 + s_Xworld[16*" + str(j) + " + 5]*ax1 + s_Xworld[16*" + str(j) + " + 9]*ax2;")
                self.gen_add_code_line("  T axw_2 = s_Xworld[16*" + str(j) + " + 2]*ax0 + s_Xworld[16*" + str(j) + " + 6]*ax1 + s_Xworld[16*" + str(j) + " + 10]*ax2;")
                a = repr(float(alpha))
                if job["revolute"]:
                    self.gen_add_code_line("  s_Jw[" + str(col_base) + " + 0] += static_cast<T>(" + a + ") * axw_0; s_Jw[" + str(col_base) + " + 1] += static_cast<T>(" + a + ") * axw_1; s_Jw[" + str(col_base) + " + 2] += static_cast<T>(" + a + ") * axw_2;")
                    self.gen_add_code_line("  T dx = s_Xworld[16*" + str(ee_anchor) + " + 12] - s_Xworld[16*" + str(j) + " + 12];")
                    self.gen_add_code_line("  T dy = s_Xworld[16*" + str(ee_anchor) + " + 13] - s_Xworld[16*" + str(j) + " + 13];")
                    self.gen_add_code_line("  T dz = s_Xworld[16*" + str(ee_anchor) + " + 14] - s_Xworld[16*" + str(j) + " + 14];")
                    self.gen_add_code_line("  s_Jv[" + str(col_base) + " + 0] += static_cast<T>(" + a + ") * (axw_1*dz - axw_2*dy);")
                    self.gen_add_code_line("  s_Jv[" + str(col_base) + " + 1] += static_cast<T>(" + a + ") * (axw_2*dx - axw_0*dz);")
                    self.gen_add_code_line("  s_Jv[" + str(col_base) + " + 2] += static_cast<T>(" + a + ") * (axw_0*dy - axw_1*dx);")
                else:
                    self.gen_add_code_line("  s_Jv[" + str(col_base) + " + 0] += static_cast<T>(" + a + ") * axw_0; s_Jv[" + str(col_base) + " + 1] += static_cast<T>(" + a + ") * axw_1; s_Jv[" + str(col_base) + " + 2] += static_cast<T>(" + a + ") * axw_2;")
                self.gen_add_code_line("}")
        self.gen_add_end_control_flow()
        self.gen_add_sync()

    # ============ Step 4: per-ee rpy sincos ============
    self.gen_add_code_line("//")
    self.gen_add_code_line("// Step 4: extract (cy, sy, cp, sp) from each ee's world rotation for E(rpy)^{-1}")
    self.gen_add_code_line("//")
    self.gen_add_parallel_loop("ee", str(num_ees))
    # bake the ee anchor jid via select
    if num_ees > 1:
        select_var_vals = [("int", "ee_jid", [str(a) for a in anchors])]
        self.gen_add_multi_threaded_select("ee", "<", [str(i+1) for i in range(num_ees)], select_var_vals)
    else:
        self.gen_add_code_line("const int ee_jid = " + str(anchors[0]) + ";")
    # R world is column-major: R[r,c] = s_Xworld[16*ee_jid + r + 4*c], r,c in 0..2
    # roll  = atan2(R[2,1], R[2,2]) -> ind 2 + 4*1 = 6 ; 2 + 4*2 = 10
    # pitch = atan2(-R[2,0], sqrt(R[2,2]^2 + R[2,1]^2)) -> ind 2 + 4*0 = 2 ; 10, 6
    # yaw   = atan2(R[1,0], R[0,0]) -> ind 1, 0
    self.gen_add_code_line("T R20 = s_Xworld[16*ee_jid + 2];")
    self.gen_add_code_line("T R21 = s_Xworld[16*ee_jid + 6];")
    self.gen_add_code_line("T R22 = s_Xworld[16*ee_jid + 10];")
    self.gen_add_code_line("T R10 = s_Xworld[16*ee_jid + 1];")
    self.gen_add_code_line("T R00 = s_Xworld[16*ee_jid + 0];")
    self.gen_add_code_line("T cp_term = sqrt(R22*R22 + R21*R21);")
    self.gen_add_code_line("T yaw = atan2(R10, R00);")
    self.gen_add_code_line("T pitch = atan2(-R20, cp_term);")
    self.gen_add_code_line("s_E_sc[4*ee + 0] = cos(yaw);")
    self.gen_add_code_line("s_E_sc[4*ee + 1] = sin(yaw);")
    self.gen_add_code_line("s_E_sc[4*ee + 2] = cos(pitch);")
    self.gen_add_code_line("s_E_sc[4*ee + 3] = sin(pitch);")
    self.gen_add_end_control_flow()
    self.gen_add_sync()

    # ============ Step 5: write s_end_effector_pose_gradient = [J_v ; E^{-1} J_w] ============
    self.gen_add_code_line("//")
    self.gen_add_code_line("// Step 5: write s_end_effector_pose_gradient (rows 0..2 = J_v, rows 3..5 = E(rpy)^{-1} J_w)")
    self.gen_add_code_line("//")
    self.gen_add_parallel_loop("ind", str(6 * nv * num_ees))
    self.gen_add_code_line("int row = ind % 6; int rem = ind / 6; int vi = rem % " + str(nv) + "; int ee = rem / " + str(nv) + ";")
    self.gen_add_code_line("int jv_base = 3 * (" + str(nv) + " * ee + vi);")
    self.gen_add_code_line("if (row < 3) {", True)
    self.gen_add_code_line("s_end_effector_pose_gradient[ind] = s_Jv[jv_base + row];")
    self.gen_add_end_control_flow()
    self.gen_add_code_line("else {", True)
    self.gen_add_code_line("T cy = s_E_sc[4*ee + 0]; T sy = s_E_sc[4*ee + 1]; T cp = s_E_sc[4*ee + 2]; T sp = s_E_sc[4*ee + 3];")
    self.gen_add_code_line("T Jw0 = s_Jw[jv_base + 0]; T Jw1 = s_Jw[jv_base + 1]; T Jw2 = s_Jw[jv_base + 2];")
    self.gen_add_code_line("T outv;")
    self.gen_add_code_line("if (row == 3) { outv = (cy*Jw0 + sy*Jw1) / cp; }")
    self.gen_add_code_line("else if (row == 4) { outv = -sy*Jw0 + cy*Jw1; }")
    self.gen_add_code_line("else { outv = (sp / cp) * (cy*Jw0 + sy*Jw1) + Jw2; }")
    self.gen_add_code_line("s_end_effector_pose_gradient[ind] = outv;")
    self.gen_add_end_control_flow()
    self.gen_add_end_control_flow()
    self.gen_add_sync()
    self.gen_add_end_function()

def gen_end_effector_pose_gradient_device(self, fixed_target_name = ""):
    n = self.robot.get_num_pos()
    nv = self.robot.get_num_vel()
    num_ees = self.robot.get_total_leaf_nodes() if fixed_target_name == "" else 1
    # construct the boilerplate and function definition
    func_params = ["s_end_effector_pose_gradient is a pointer to shared memory of size 6*NUM_VEL*NUM_EE where NUM_VEL = " + str(nv) + " and NUM_EE = " + str(num_ees), \
                   "s_q is the vector of joint positions", \
                   "d_robotModel is the pointer to the initialized model specific helpers on the GPU (XImats, topology_helpers, etc.)"]
    func_notes = []
    func_def_start = "void end_effector_pose_gradient_device" + ("" if fixed_target_name == "" else "_" + fixed_target_name) + "("
    func_def_middle = "T *s_end_effector_pose_gradient, const T *s_q, "
    func_def_end = "const robotModel<T> *d_robotModel) {"
    func_def = func_def_start + func_def_middle + func_def_end
    # then generate the code
    self.gen_add_func_doc("Computes the Gradient of the End Effector Pose with respect to joint position",\
                          func_notes,func_params,None)
    self.gen_add_code_line("template <typename T>")
    self.gen_add_code_line("__device__")
    self.gen_add_code_line(func_def, True)
    # add the shared memory variables. The shared-chain geometric-Jacobian inner
    # uses ONLY local Xhom (s_dXhom is unused, marked `(void)`) so skip the
    # per-joint local d-transform allocation + computation entirely. On floating
    # base this also skips the (expensive) quaternion derivative of the base
    # transform, which dominated the old per-(djid, ee) re-chain cost.
    shared_mem_size = self.gen_end_effector_pose_gradient_inner_temp_mem_size(fixed_target_name)
    self.gen_XmatsHom_helpers_temp_shared_memory_code(shared_mem_size, include_gradients = False, include_linalg_scratch = True,
                                                      linalg_scratch_bytes = "GRID_EE_LINALG_SHARED_BYTES<T>()")
    # then load/update XI and run the algo
    self.gen_load_update_XmatsHom_helpers_function_call(include_gradients = False)
    self.gen_end_effector_pose_gradient_inner_function_call(fixed_target_name = fixed_target_name,
        updated_var_names = {"s_dXhom_name": "nullptr"})
    self.gen_add_end_function()

_EE_GRAD_PICK_FLAGS = [
    # (use_workspace_temp, use_workspace_dxhom)
    (False, False),   # pick 0: full smem (PERF)
    (True,  False),   # pick 1: inner_temp + s_end_effector_pose_gradient -> workspace/global (LITE)
    (True,  True),    # pick 2: also dXmatsHom -> workspace (MINIMAL)
]

def _emit_eepose_grad_kernel_body_for_flags(self, n, num_ees, fixed_target_name,
                                            use_workspace_temp, use_workspace_dxhom,
                                            single_call_timing):
    """Emit the EE_POSE_GRAD kernel body specialized for one tier's spill flags.
    Wrapped in a brace pair (caller emits the `if constexpr (...)` head).
    Used by gen_end_effector_pose_gradient_kernel to emit either a single body
    (collapsed picks) or three branched bodies (divergent picks). Mirrors
    _emit_d2ee_kernel_body_for_flags."""
    nv = self.robot.get_num_vel()
    shared_mem_size = 0 if use_workspace_temp else self.gen_end_effector_pose_gradient_inner_temp_mem_size(fixed_target_name)
    extra_t_buffers = [("s_q", n)] if use_workspace_temp else [("s_q", n), ("s_end_effector_pose_gradient", 6*nv*num_ees)]
    # Geometric-Jacobian inner doesn't use s_dXhom -> skip its allocation/computation entirely.
    self.gen_XmatsHom_helpers_temp_shared_memory_code(shared_mem_size, include_gradients = False,
                                                      extra_t_buffers = extra_t_buffers,
                                                      include_linalg_scratch = True,
                                                      linalg_scratch_bytes = "GRID_EE_LINALG_SHARED_BYTES<T>()")
    if not use_workspace_temp:
        self.gen_add_code_line("(void)d_workspace;")
    # Per-tier eegrad_temp byte offset. The shared GRID_END_EFFECTOR_POSE_GRADIENT_WORKSPACE_TEMP_OFFSET_BYTES
    # macro keys off the single-valued PERF-pick GRID_END_EFFECTOR_POSE_GRADIENT_USES_WORKSPACE_DXHOM, so it
    # would collide with the spilled dXhom region at tiers whose pick spills dXhom but whose
    # PERF pick does not (e.g. go2 end_effector_pose_gradient = (0,0,2)). Compute the offset locally from THIS
    # tier's use_workspace_dxhom so the temp arena always lands past the spilled dXhom region.
    eegrad_temp_off = "GRID_END_EFFECTOR_POSE_GRADIENT_WORKSPACE_DXHOM_OFFSET_BYTES<T>()"
    if use_workspace_dxhom:
        eegrad_temp_off += " + sizeof(T) * static_cast<size_t>(DXHOM_T_COUNT)"
    if not single_call_timing:
        self.gen_add_parallel_loop("k","NUM_TIMESTEPS",block_level = True)
        self.gen_kernel_load_inputs("q",str(n),stride="stride_q")
        if use_workspace_dxhom:
            self.gen_add_code_line("T *s_dXmatsHom = reinterpret_cast<T *>(&d_workspace[k*GRID_WORKSPACE_BYTES_PER_TIMESTEP<T>() + GRID_END_EFFECTOR_POSE_GRADIENT_WORKSPACE_DXHOM_OFFSET_BYTES<T>()]);")
        if use_workspace_temp:
            self.gen_add_code_line("T *s_end_effector_pose_gradient = &d_end_effector_pose_gradient[k*" + str(6*nv*num_ees) + "];")
            self.gen_add_code_line("T *s_eegrad_temp = reinterpret_cast<T *>(&d_workspace[k*GRID_WORKSPACE_BYTES_PER_TIMESTEP<T>() + " + eegrad_temp_off + "]);")
            # Whole inner arena spilled -> smem s_temp is null. Repoint it at the
            # spilled workspace so the XmatsHom helper's sincos scratch is backed.
            self.gen_add_code_line("s_temp = s_eegrad_temp;")
        self.gen_add_code_line("// compute")
        self.gen_load_update_XmatsHom_helpers_function_call(include_gradients = False)
        # Inner-controlled: pass both arenas + placement; the inner picks where
        # the chain workspace lives via TEMP_IN_SMEM. s_dXhom is unused by the
        # shared-chain geometric-Jacobian inner -> pass nullptr.
        updated = {"d_workspace_name": "s_eegrad_temp"} if use_workspace_temp else {}
        updated["s_dXhom_name"] = "nullptr"
        self.gen_end_effector_pose_gradient_inner_function_call(fixed_target_name = fixed_target_name,
            updated_var_names = updated, temp_in_smem_expr = ("false" if use_workspace_temp else "true"))
        self.gen_add_sync()
        if not use_workspace_temp:
            self.gen_kernel_save_result("end_effector_pose_gradient",str(6*nv*num_ees),stride=str(6*nv*num_ees))
        self.gen_add_end_control_flow()
    else:
        self.gen_kernel_load_inputs("q",str(n))
        if use_workspace_dxhom:
            self.gen_add_code_line("T *s_dXmatsHom = reinterpret_cast<T *>(&d_workspace[GRID_END_EFFECTOR_POSE_GRADIENT_WORKSPACE_DXHOM_OFFSET_BYTES<T>()]);")
        if use_workspace_temp:
            self.gen_add_code_line("T *s_end_effector_pose_gradient = d_end_effector_pose_gradient;")
            self.gen_add_code_line("T *s_eegrad_temp = reinterpret_cast<T *>(&d_workspace[" + eegrad_temp_off + "]);")
            # See note above: repoint the null smem s_temp at the spilled workspace.
            self.gen_add_code_line("s_temp = s_eegrad_temp;")
        self.gen_add_code_line("// compute with NUM_TIMESTEPS as NUM_REPS for timing")
        self.gen_add_code_line("for (int rep = 0; rep < NUM_TIMESTEPS; rep++){", True)
        # TODO(licm-eepose-grad): sm_86-specific, deprioritized. See pre-Phase-3d note in git history.
        self.gen_anti_licm_input_reload("q",str(n),feedback_from="end_effector_pose_gradient")
        self.gen_load_update_XmatsHom_helpers_function_call(include_gradients = False)
        updated = {"d_workspace_name": "s_eegrad_temp"} if use_workspace_temp else {}
        updated["s_dXhom_name"] = "nullptr"
        self.gen_end_effector_pose_gradient_inner_function_call(fixed_target_name = fixed_target_name,
            updated_var_names = updated, temp_in_smem_expr = ("false" if use_workspace_temp else "true"))
        self.gen_anti_licm_output_write("end_effector_pose_gradient")
        self.gen_add_end_control_flow()
        if not use_workspace_temp:
            self.gen_kernel_save_result("end_effector_pose_gradient",str(6*nv*num_ees))


def gen_end_effector_pose_gradient_kernel(self, single_call_timing = False, fixed_target_name = ""):
    n = self.robot.get_num_pos()
    num_ees = self.robot.get_total_leaf_nodes() if fixed_target_name == "" else 1
    func_params = ["d_end_effector_pose_gradient is the vector of end effector positions gradients", \
                   "d_workspace is the generated global spill workspace", \
                   "d_q is the vector of joint positions", \
                   "stride_q is the stide between each q", \
                   "d_robotModel is the pointer to the initialized model specific helpers on the GPU (XImats, topology_helpers, etc.)", \
                   "num_timesteps is the length of the trajectory points we need to compute over (or overloaded as test_iters for timing)"]
    func_notes = []
    func_def_start = "void end_effector_pose_gradient_kernel" + ("" if fixed_target_name == "" else "_" + fixed_target_name) + "(T *d_end_effector_pose_gradient, unsigned char *d_workspace, const T *d_q, const int stride_q, "
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
    # tier's spill flags. END_EFFECTOR_POSE_GRADIENT_DYNAMIC_SHARED_MEM_BYTES<T, TIER>() is
    # tier-aware.
    picks = getattr(self, "end_effector_pose_gradient_spill_tier_3way", (0, 0, 0))
    def _emit_end_effector_pose_gradient_body(pick):
        uwt, uwd = _EE_GRAD_PICK_FLAGS[pick]
        _emit_eepose_grad_kernel_body_for_flags(self, n, num_ees, fixed_target_name, uwt, uwd, single_call_timing)
    self.gen_tier_dispatch(picks, _emit_end_effector_pose_gradient_body)
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
    func_call_start = "end_effector_pose_gradient_kernel" + ("" if fixed_target_name == "" else "_" + fixed_target_name) + "<T><<<block_dimms,thread_dimms,END_EFFECTOR_POSE_GRADIENT_DYNAMIC_SHARED_MEM_BYTES<T>()>>>(hd_data->d_end_effector_pose_gradient,hd_data->d_workspace,hd_data->d_q,stride_q,"
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
    self.gen_add_code_line("gpuErrchk(grid_check_dynamic_shared_memory_bytes(\"end_effector_pose_gradient\", END_EFFECTOR_POSE_GRADIENT_DYNAMIC_SHARED_MEM_BYTES<T>()));")
    workspace_bytes = "GRID_WORKSPACE_BYTES_PER_TIMESTEP<T>()" if single_call_timing else "GRID_WORKSPACE_BYTES_PER_TIMESTEP<T>()*static_cast<size_t>(num_timesteps)"
    # Per-tier gate: arm L2 persistence if ANY tier routes the chain workspace
    # through d_workspace (runtime RESOURCE_TIER may differ from the PERF pick).
    self.gen_add_code_line("if (GRID_END_EFFECTOR_POSE_GRADIENT_USES_WORKSPACE_TEMP_ANY) {gpuErrchk(grid_begin_l2_persisting(0, hd_data->d_workspace, " + workspace_bytes + "));}")
    self.gen_add_code_lines(func_call_code)
    self.gen_add_code_line("if (GRID_END_EFFECTOR_POSE_GRADIENT_USES_WORKSPACE_TEMP_ANY) {gpuErrchk(grid_end_l2_persisting(0));}")
    if not compute_only:
        # then transfer memory back
        self.gen_add_code_lines(["// finally transfer the result back", \
                                 "gpuErrchk(cudaMemcpy(hd_data->h_end_effector_pose_gradient,hd_data->d_end_effector_pose_gradient,6*NUM_EES*NUM_VEL*" + \
                                    ("num_timesteps*" if not single_call_timing else "") + "sizeof(T),cudaMemcpyDeviceToHost));",
                                 "gpuErrchkKernel();"])
    # finally report out timing if requested
    if single_call_timing:
        from ..algo_registry import single_call_printf_line
        self.gen_add_code_line(single_call_printf_line("end_effector_pose_gradient"))
    self.gen_add_end_function()

def gen_end_effector_pose_hessian_output_count(self):
    """Number of T elements in the end_effector_pose_hessian output: 6 * nv * nv * num_ees.

    Output is now d^2(pose)/dv^2 (TANGENT, pinocchio convention). For fixed-base
    nv == nq so the size is unchanged; for floating-base the (nv x nv) block now
    indexes spatial twist components rather than the older non-standard
    quaternion derivatives.
    """
    nv = self.robot.get_num_vel()
    num_ees = self.robot.get_total_leaf_nodes()
    return 6 * nv * nv * num_ees

def gen_end_effector_pose_hessian_inner_temp_mem_size(self):
    """Size (in T elements) of the analytic d2ee inner's s_temp.

    The closed-form per-chain second-order Taylor algorithm (see
    docs/d2ee_analytic_derivation.md) needs:

      [0 .. 16*n_joints)               s_Xworld     world transform of every joint
                                                    (shared FK pass, identical to
                                                    end_effector_pose_gradient_inner)
      [+ 16*nv*num_ees)                s_Sworld     per-DOF world-frame 4x4 generator
                                                    L_a * A_i_local * L_a^{-1}
                                                    (top-left 3x3 = skew for revolute /
                                                    zeros for prismatic; column 3 = the
                                                    "twist origin offset" piece)
      [+ 4*num_ees)                    s_E_sc       cy, sy, cp, sp per ee for E(rpy)^-1
    """
    n_joints = self.robot.get_num_joints()
    nv = self.robot.get_num_vel()
    num_ees = self.robot.get_total_leaf_nodes()
    return 16*n_joints + 16*nv*num_ees + 4*num_ees

def gen_end_effector_pose_hessian_inner_function_call(self, updated_var_names = None,
                                                               out_in_smem_expr = "true"):
    var_names = dict( \
        s_Xhom_name = "s_XmatsHom", \
        s_end_effector_pose_gradient_name = "s_end_effector_pose_gradient", \
        s_end_effector_pose_hessian_name = "s_end_effector_pose_hessian", \
        s_q_name = "s_q", \
        s_topology_helpers_name = "s_topology_helpers", \
        s_temp_name = "s_temp", \
        d_workspace_name = "nullptr", \
        d_robotModel_name = "d_robotModel", \
        s_linalg_smem_name = "s_linalg_smem", \
    )
    if updated_var_names is not None:
        for key,value in updated_var_names.items():
            var_names[key] = value
    code_start = "end_effector_pose_hessian_inner<T, " + out_in_smem_expr + ">(" + var_names["s_end_effector_pose_hessian_name"] + ", " + var_names["s_end_effector_pose_gradient_name"] + ", " + var_names["s_q_name"] + ", "
    code_middle = var_names["s_Xhom_name"] + ", "
    code_end = var_names["s_temp_name"] + ", " + var_names["d_workspace_name"] + ", " + var_names["d_robotModel_name"] + ", " + var_names["s_linalg_smem_name"] + ");"
    # account for thread group
    # Canonical: append the shared topology-helper arg via the central helper
    # (NO_XI: the ee_pose family takes s_Xhom, not s_XImats). Mirrors the def's
    # gen_insert_helpers_func_def_params(NO_XI_FLAG=True) so def + call can't drift.
    code_middle += self.gen_insert_helpers_function_call(updated_var_names = var_names, NO_XI_FLAG = True)
    self.gen_add_code_line(code_start + code_middle + code_end)

def _eepose_hessian_chain_metadata(self, all_ees):
    """Per-ee chain bookkeeping for the analytic d2ee inner.

    Returns (chains, anchors, per_ee_dof_info, intra_joint_pairs_per_ee):
      chains[ee_idx]                 = sorted list of joint ids on the chain root..ee
      anchors[ee_idx]                = the joint id whose Xworld is the EE world transform
      per_ee_dof_info[ee_idx]        = list of dicts:
        {vi, chain_pos, S_col, joint_jid, ang (3), lin (3), revolute (bool)}
        one entry per chain DOF; vi is the v-space index, S_col is the column of
        the joint's S matrix (0 for single-DOF joints, 0..5 for the floating base).
      intra_joint_pairs_per_ee[ee_idx] = list of (vi_a, vi_b, joint_chain_pos, c_a, c_b, joint_jid)
        for every UNORDERED pair (a, b) of DOFs that live in the same multi-DOF
        joint on the chain. For typical revolute/prismatic 1-DOF joints there are
        no such pairs; only the floating-base jid=0 contributes (15 unique pairs
        for nv >= 6 of the 6 base DOFs).
    """
    import numpy as _np
    chains, anchors, per_ee_dof_info, intra_joint_pairs_per_ee = [], [], [], []
    for ee in all_ees:
        chain = sorted(self.robot.get_ancestors_by_id(ee)) + [ee]
        chains.append(chain)
        anchors.append(ee)
        dof_info = []
        intra_pairs = []
        for chain_pos, j in enumerate(chain):
            S = _np.asarray(self.robot.get_S_by_id(j), dtype=_np.float64)
            if S.ndim == 1:
                S = S.reshape(-1, 1)
            try:
                vinds = self.robot.get_joint_index_v(j)
            except Exception:
                vinds = self.robot.get_joint_index_q(j)
            if not isinstance(vinds, (list, tuple, _np.ndarray)):
                vinds = [vinds]
            vinds = list(vinds)
            this_joint_dofs = []
            for c in range(S.shape[1]):
                vi = vinds[c] if c < len(vinds) else vinds[-1]
                ang_local = [float(x) for x in S[:3, c]]
                lin_local = [float(x) for x in S[3:6, c]]
                revolute = max(abs(x) for x in ang_local) > 0.5
                dof_info.append({
                    "vi": int(vi),
                    "chain_pos": chain_pos,
                    "S_col": c,
                    "joint_jid": int(j),
                    "ang": ang_local,
                    "lin": lin_local,
                    "revolute": bool(revolute),
                })
                this_joint_dofs.append((int(vi), c))
            # collect intra-joint UNORDERED pairs (a <= b in S column order)
            if len(this_joint_dofs) > 1:
                for ai, (vi_a, c_a) in enumerate(this_joint_dofs):
                    for bi, (vi_b, c_b) in enumerate(this_joint_dofs):
                        if ai > bi:
                            continue
                        intra_pairs.append((vi_a, vi_b, chain_pos, c_a, c_b, int(j)))
        per_ee_dof_info.append(dof_info)
        intra_joint_pairs_per_ee.append(intra_pairs)
    return chains, anchors, per_ee_dof_info, intra_joint_pairs_per_ee


def gen_end_effector_pose_hessian_inner(self):
    """Analytic d^2(pose)/dv^2 of the end-effector pose via per-chain second-
    order Taylor expansion (see docs/d2ee_analytic_derivation.md).

    Replaces the previous FD-on-d/dv-Jacobian implementation (2*nv + 1 gradient
    calls) with a single closed-form pass. Mirrors
    `RBDReference.end_effector_pose_hessian_analytic`, which agrees with
    pinocchio's analytic `getJointKinematicHessian(LOCAL_WORLD_ALIGNED)` to the
    FD-noise floor (~1e-5 rel) across the full manifest fleet (iiwa14/go2/g1/
    h1_2/fr3/rizon4/gen3/fetch/baxter, fixed + floating). The CUDA path is
    confirmed vs the pinocchio oracle on iiwa14-fixed + go2-floating (the
    floating orientation-hessian block was the old B1 bug; the chain-composition
    d^2/dv^2 derivation is correct fleet-wide on both surfaces).

    Output convention: d^2(pose)/dv^2 (TANGENT, pinocchio convention), shape
    (num_ees, 6, nv, nv) in row-major (C-order): linear index
    e*6*nv*nv + c*nv*nv + j*nv + i.

    Algorithm summary:
      1. Forward kinematics: world transform of every joint (s_Xworld).
      2. Build per-DOF world-frame 4x4 generator S_i_world = L_a*A_i_local*L_a^{-1}.
         For revolute axis a_local (chain joint a with world transform Xw_a):
           S_world = [[ [Rw_a a_local]_x, -[Rw_a a_local]_x * pw_a ],
                      [ 0,                0                       ]]
         For prismatic (linear) axis a_local:
           S_world = [[ 0, Rw_a a_local ], [ 0, 0 ]]
         The per-DOF angular axis is then skew_inv(S_world[:3,:3]) = Rw_a*ang_local
         (revolute) or 0 (prismatic). J_v = S_world[:3,3] + S_world[:3,:3]*p_ee.
      3. Emit s_end_effector_pose_gradient = [J_v; E^{-1}*J_w] using the same closed-form E^{-1}
         the gradient inner uses.
      4. For each DOF pair (i, j) with proximal/distal joints (a<=b in chain):
           if a < b: d2M = S_prox_world * S_dist_world * X_ee
           if a == b (intra-joint, only floating base): d2M = B_world * X_ee
             where B_world = L_a*B_local*L_a^{-1} (closed form -- see comments).
         Then:
           H_xyz[:, i, j]      = (d2M * ee_offset)[:3]  with ee_offset = [0,0,0,1]
           d2R_R^T            = d2M[:3,:3]_top_of_factor  (the X_ee factor cancels
                                with R_chain^T since they are equal for joint-EEs)
           H_w[:, i, j]       = skew_inv(d2R_R^T - [J_w_i]_x * [J_w_j]_x)
           H_rpy[:, i, j]     = (dEinv/dv_j) * J_w[:, i] + Einv * H_w[:, i, j]
                                (closed-form dE/drpy chain rule)
      5. Symmetrize H_rpy over the (i, j) Hessian axes; H_xyz is symmetric by
         construction (d2M[i,j] == d2M[j,i] in the formula above).

    No FD step; no perturbed q recomputes; no integrate(). Pure closed-form,
    O(nv^2 * const) per ee.
    """
    nv = self.robot.get_num_vel()
    nq = self.robot.get_num_pos()
    n_joints = self.robot.get_num_joints()
    n_bfs_levels = self.robot.get_max_bfs_level() + 1
    all_ees = self.robot.get_leaf_nodes()
    num_ees = len(all_ees)
    chains, anchors, per_ee_dof_info, intra_joint_pairs_per_ee = \
        _eepose_hessian_chain_metadata(self, all_ees)

    # scratch offsets
    off_Xworld = 0
    off_Sworld = off_Xworld + 16 * n_joints
    off_Esc    = off_Sworld + 16 * nv * num_ees

    func_params = [
        "s_end_effector_pose_hessian is a pointer to memory of size 6*NUM_VEL*NUM_VEL*NUM_EE where NUM_VEL = " + str(nv) + " and NUM_EE = " + str(num_ees) +
            " (d^2(pose)/dv^2 tangent-space Hessian, pinocchio convention)",
        "s_end_effector_pose_gradient is a pointer to memory of size 6*NUM_VEL*NUM_EE (the d/dv tangent Jacobian at q)",
        "s_q is the vector of joint positions (size NUM_POS = " + str(nq) + "; kept for signature compatibility, unused by the analytic path)",
        "s_Xhom is the per-joint LOCAL homogeneous-transform buffer (read-only)",
        "s_temp is helper shared memory of size " + str(self.gen_end_effector_pose_hessian_inner_temp_mem_size()) +
            " (s_Xworld | s_Sworld | s_E_sc; always kept in smem)",
        "d_workspace is the global spill arena s_end_effector_pose_hessian is repointed at when !OUT_IN_SMEM (else unused)",
        "d_robotModel is the model-specific helper struct (kept for signature compatibility, unused by the analytic path)",
        "s_linalg_smem is optional byte-addressed shared memory (reserved; unused by this inner)",
    ]
    func_notes = [
        "Closed-form analytic d2(pose)/dv2; matches RBDReference.end_effector_pose_hessian_analytic (which agrees with pinocchio getJointKinematicHessian(LOCAL_WORLD_ALIGNED) to the FD floor fleet-wide; CUDA confirmed on iiwa14-fixed + go2-floating).",
        "Inner-owns scratch placement: the large nv^2 output s_end_effector_pose_hessian moves to d_workspace when !OUT_IN_SMEM. The s_Xworld+s_Sworld+s_E_sc scratch in s_temp stays in smem at every tier.",
    ]
    func_def_start = "void end_effector_pose_hessian_inner("
    func_def_middle = "T *s_end_effector_pose_hessian, T *s_end_effector_pose_gradient, const T *s_q, T *s_Xhom, "
    func_def_end = "T *s_temp, T *d_workspace, const robotModel<T> *d_robotModel, unsigned char *s_linalg_smem) {"
    func_def_middle, func_params = self.gen_insert_helpers_func_def_params(func_def_middle, func_params, -1, NO_XI_FLAG = True)
    func_def = func_def_start + func_def_middle + func_def_end
    self.gen_add_func_doc(
        "Computes the Hessian (and Jacobian) of the End Effector Pose with respect to generalized velocity (d^2/dv^2 tangent, pinocchio convention)",
        func_notes, func_params, None)
    self.gen_add_code_line("template <typename T, bool OUT_IN_SMEM = true>")
    self.gen_add_code_line("__device__")
    self.gen_add_code_line(func_def, True)
    # Inner-controlled scratch placement: the large nv^2 output moves to
    # d_workspace when !OUT_IN_SMEM. Reassigning s_end_effector_pose_hessian here keeps every
    # s_end_effector_pose_hessian[...] reference below unchanged.
    self.gen_add_code_line("if constexpr (!OUT_IN_SMEM) { s_end_effector_pose_hessian = d_workspace; } else { (void)d_workspace; }")
    self.gen_add_code_line("(void)s_q; (void)d_robotModel; (void)s_linalg_smem;")
    self.gen_add_code_line("// scratch in s_temp: s_Xworld (16*n_joints) | s_Sworld (16*nv*num_ees) | s_E_sc (4*num_ees)")
    self.gen_add_code_line("T *s_Xworld = &s_temp[" + str(off_Xworld) + "];")
    self.gen_add_code_line("T *s_Sworld = &s_temp[" + str(off_Sworld) + "];  // per-DOF world-frame 4x4 generator (S_i_world)")
    self.gen_add_code_line("T *s_E_sc   = &s_temp[" + str(off_Esc)    + "];  // cy, sy, cp, sp per ee")

    # ===== Step 1: forward kinematics — world transforms by BFS level =====
    # (identical to the gradient inner's Step 1; populates s_Xworld[16*j] for every j)
    self.gen_add_code_line("//")
    self.gen_add_code_line("// Step 1: forward kinematics -- build s_Xworld[16*j] for every joint via BFS-level chain-up")
    self.gen_add_code_line("//")
    for level in range(n_bfs_levels):
        ids_at_level = self.robot.get_ids_by_bfs_level(level)
        if not ids_at_level:
            continue
        njs = len(ids_at_level)
        self.gen_add_code_line("// BFS level " + str(level) + " -> joints " + str(ids_at_level))
        self.gen_add_parallel_loop("ind", str(16 * njs))
        self.gen_add_code_line("int slot = ind / 16; int ele = ind % 16;")
        self.gen_add_code_line("int row = ele & 3; int col = ele >> 2;")
        jid_list = [str(j) for j in ids_at_level]
        par_list = [str(self.robot.get_parent_id(j)) for j in ids_at_level]
        select_var_vals = [("int", "jid", jid_list), ("int", "par", par_list)]
        self.gen_add_multi_threaded_select("slot", "<", [str(i+1) for i in range(njs)], select_var_vals)
        self.gen_add_code_line("if (par == -1) {", True)
        self.gen_add_code_line("s_Xworld[16*jid + ele] = s_Xhom[16*jid + ele];")
        self.gen_add_end_control_flow()
        self.gen_add_code_line("else {", True)
        self.gen_add_code_line("s_Xworld[16*jid + ele] = dot_prod<T,4,4,1>(&s_Xworld[16*par + row], &s_Xhom[16*jid + 4*col]);")
        self.gen_add_end_control_flow()
        self.gen_add_end_control_flow()
        self.gen_add_sync()

    # ===== Step 2: per-DOF world-frame generator s_Sworld =====
    # Layout: s_Sworld[16 * (ee*nv + vi) + ele] (4x4 per DOF per ee, column-major)
    # For revolute axis a_local in chain joint a (world Xw_a = (Rw_a, pw_a)):
    #   ω_w = Rw_a @ a_local
    #   S[:3,:3] = [ω_w]_×; S[:3,3] = -[ω_w]_× * pw_a = pw_a × ω_w; bottom row = 0
    # For prismatic axis a_local:
    #   S[:3,:3] = 0; S[:3,3] = Rw_a @ a_local; bottom row = 0
    self.gen_add_code_line("//")
    self.gen_add_code_line("// Step 2: build per-DOF world-frame 4x4 generator S_i_world")
    self.gen_add_code_line("//")
    # First zero all of s_Sworld (out-of-chain DOFs stay zero — they contribute nothing).
    self.gen_add_parallel_loop("ind", str(16 * nv * num_ees))
    self.gen_add_code_line("s_Sworld[ind] = static_cast<T>(0);")
    self.gen_add_end_control_flow()
    self.gen_add_sync()
    # Then emit per (ee, chain-joint, S-col) the explicit 4x4 fill. Serial-ops
    # per slot — total work is small (chain_depth * dofs_per_joint * num_ees blocks).
    #
    # MIMIC fold: a mimic joint and its target share one velocity slot vi, so
    # several chain joints map to the SAME s_Sworld[16*(ee*nv+vi)] generator. S_world
    # is LINEAR in the joint rate, so the shared column is the alpha-weighted SUM of
    # each contributing joint's generator (the mimic body moves alpha*target_rate).
    # Every downstream step (d2M products, J_w/J_v readout, rpy chain rule) reads only
    # s_Sworld[vi] / s_end_effector_pose_gradient[vi], so folding the generator here is sufficient. The
    # signed S column is baked into `ax`, so (unlike the unit-axis inverse_dynamics_gradient path) the only
    # scalar fold is alpha — no separate s_sign. Non-mimic: each (ee,vi) slot written
    # once with alpha==1.0 => byte-identical to the legacy "=" assignment.
    #
    # PARALLELIZATION: each (ee,vi) s_Sworld slot is a DISJOINT 16-float region,
    # zero-filled above with a sync, so slots are independent. Group chain DOFs by
    # (ee_idx, vi) and dispatch one-thread-per-slot over a flat index — every writer
    # to a slot lives in the SAME guard (first writer "=", later mimic writers "+=")
    # so the mimic reduction is never split across threads.
    HAS_MIMIC = self.robot_has_mimic_joints()
    slot_groups = []     # list of (ee_idx, vi, [dof, ...]) in stable first-seen order
    slot_index = {}      # (ee_idx, vi) -> position in slot_groups
    for ee_idx in range(num_ees):
        for dof in per_ee_dof_info[ee_idx]:
            key = (ee_idx, dof["vi"])
            if key not in slot_index:
                slot_index[key] = len(slot_groups)
                slot_groups.append((ee_idx, dof["vi"], []))
            slot_groups[slot_index[key]][2].append(dof)

    def _emit_sworld_slot(ee_idx, slot_vi, dofs):
        # Emit every chain-DOF writer that folds into this (ee_idx, slot_vi) slot,
        # in chain order, with first-writer "=" and later (mimic) writers "+=".
        for w, dof in enumerate(dofs):
            vi = dof["vi"]
            j = dof["joint_jid"]
            ang = dof["ang"]
            lin = dof["lin"]
            rev = dof["revolute"]
            base = 16 * (ee_idx * nv + vi)
            ax = ang if rev else lin
            alpha = self._alpha_for_jid(j) if HAS_MIMIC else 1.0
            first_writer = (w == 0)
            # First writer to a fresh slot assigns ("="); a later writer (only the
            # mimic case) accumulates ("+="). With alpha == 1.0 and first_writer the
            # emitted text matches the legacy path exactly.
            assign = "=" if first_writer else "+="
            def _val(rhs):
                # Wrap rhs by the mimic scale when alpha != 1.0; otherwise leave it
                # untouched so non-mimic output is byte-identical.
                if alpha == 1.0:
                    return rhs
                return "static_cast<T>(" + repr(float(alpha)) + ") * (" + rhs + ")"
            def _set(slot, rhs):
                # When first_writer the zero-fill above already cleared the slot, so
                # "=" and "+=" are numerically identical; we keep "=" so the legacy
                # (non-mimic) text is preserved exactly. A scaled accumulate from a
                # later mimic joint folds in additively.
                self.gen_add_code_line("s_Sworld[" + str(slot) + "] " + assign + " " + _val(rhs) + ";")
            self.gen_add_code_line(
                "// ee=" + str(ee_idx) + " vi=" + str(vi) + " jid=" + str(j) +
                (" rev" if rev else " prism") + " ax_local=" + str(ax) +
                ("" if alpha == 1.0 else " alpha=" + repr(float(alpha)) +
                 (" (fold)" if not first_writer else " (mimic-target)")))
            self.gen_add_code_line("{", True)
            # axis_world = R_j_world @ ax_local
            # R_j_world is column-major in s_Xworld[16*j]: R[r,c] = s_Xworld[16*j + r + 4*c]
            for r in range(3):
                terms = []
                for c in range(3):
                    if abs(ax[c]) < 1e-15:
                        continue
                    coef = "static_cast<T>(" + "{:.17g}".format(ax[c]) + ")"
                    terms.append("s_Xworld[" + str(16*j + r + 4*c) + "] * " + coef)
                expr = " + ".join(terms) if terms else "static_cast<T>(0)"
                self.gen_add_code_line("T axw_" + str(r) + " = " + expr + ";")
            # p_j_world (last column, rows 0..2)
            self.gen_add_code_line("T pjx = s_Xworld[" + str(16*j + 12) + "];")
            self.gen_add_code_line("T pjy = s_Xworld[" + str(16*j + 13) + "];")
            self.gen_add_code_line("T pjz = s_Xworld[" + str(16*j + 14) + "];")
            if rev:
                # S[:3, :3] = [axw]_x, S[:3, 3] = p_j x axw (= -[axw]_x p_j)
                # Column-major: S[r + 4*c] = S[r, c]
                # [axw]_x  =  [[0, -wz,  wy],
                #              [wz, 0,  -wx],
                #              [-wy, wx, 0]]
                if first_writer and alpha == 1.0:
                    # Legacy fast path: reproduce the original emission CHARACTER
                    # FOR CHARACTER (incl. the aligned double-space before bare
                    # axw_* terms) so non-mimic grid.cuh stays byte-identical.
                    self.gen_add_code_line("s_Sworld[" + str(base +  0) + "] = static_cast<T>(0);")  # S[0,0]
                    self.gen_add_code_line("s_Sworld[" + str(base +  1) + "] =  axw_2;")             # S[1,0] =  wz
                    self.gen_add_code_line("s_Sworld[" + str(base +  2) + "] = -axw_1;")             # S[2,0] = -wy
                    self.gen_add_code_line("s_Sworld[" + str(base +  4) + "] = -axw_2;")             # S[0,1] = -wz
                    self.gen_add_code_line("s_Sworld[" + str(base +  5) + "] = static_cast<T>(0);")  # S[1,1]
                    self.gen_add_code_line("s_Sworld[" + str(base +  6) + "] =  axw_0;")             # S[2,1] =  wx
                    self.gen_add_code_line("s_Sworld[" + str(base +  8) + "] =  axw_1;")             # S[0,2] =  wy
                    self.gen_add_code_line("s_Sworld[" + str(base +  9) + "] = -axw_0;")             # S[1,2] = -wx
                    self.gen_add_code_line("s_Sworld[" + str(base + 10) + "] = static_cast<T>(0);")  # S[2,2]
                    # Column 3: p_j x axw  (= -[axw]_x p_j)
                    self.gen_add_code_line("s_Sworld[" + str(base + 12) + "] = pjy*axw_2 - pjz*axw_1;")
                    self.gen_add_code_line("s_Sworld[" + str(base + 13) + "] = pjz*axw_0 - pjx*axw_2;")
                    self.gen_add_code_line("s_Sworld[" + str(base + 14) + "] = pjx*axw_1 - pjy*axw_0;")
                else:
                    # Mimic fold (later writer or scaled): accumulate the skew axis
                    # entries; the diagonal zeros need no accumulate (the slot was
                    # zero-filled and no contributing generator touches the diagonal).
                    _set(base +  1, "axw_2")
                    _set(base +  2, "-axw_1")
                    _set(base +  4, "-axw_2")
                    _set(base +  6, "axw_0")
                    _set(base +  8, "axw_1")
                    _set(base +  9, "-axw_0")
                    # Column 3: p_j x axw  (= -[axw]_x p_j)
                    _set(base + 12, "pjy*axw_2 - pjz*axw_1")
                    _set(base + 13, "pjz*axw_0 - pjx*axw_2")
                    _set(base + 14, "pjx*axw_1 - pjy*axw_0")
            else:
                if first_writer and alpha == 1.0:
                    # Prismatic legacy fast path (byte-identical).
                    self.gen_add_code_line("s_Sworld[" + str(base + 12) + "] = axw_0;")
                    self.gen_add_code_line("s_Sworld[" + str(base + 13) + "] = axw_1;")
                    self.gen_add_code_line("s_Sworld[" + str(base + 14) + "] = axw_2;")
                else:
                    # Prismatic mimic fold: S[:3, 3] = alpha * axis_world.
                    _set(base + 12, "axw_0")
                    _set(base + 13, "axw_1")
                    _set(base + 14, "axw_2")
            self.gen_add_end_control_flow()

    # Dispatch all (ee, vi) slots one-thread-per-slot via a single block-parallel loop.
    n_slots = len(slot_groups)
    self.gen_add_parallel_loop("sworld_slot", str(n_slots))
    for k, (ee_idx, slot_vi, dofs) in enumerate(slot_groups):
        self.gen_add_code_line("if (sworld_slot == " + str(k) + ") {", True)
        _emit_sworld_slot(ee_idx, slot_vi, dofs)
        self.gen_add_end_control_flow()
    self.gen_add_end_control_flow()
    self.gen_add_sync()

    # ===== Step 3: extract (cy, sy, cp, sp) from each ee's world rotation =====
    self.gen_add_code_line("//")
    self.gen_add_code_line("// Step 3: extract (cy, sy, cp, sp) for E(rpy)^-1 / dE/drpy")
    self.gen_add_code_line("//")
    self.gen_add_parallel_loop("ee", str(num_ees))
    if num_ees > 1:
        select_var_vals = [("int", "ee_jid", [str(a) for a in anchors])]
        self.gen_add_multi_threaded_select("ee", "<", [str(i+1) for i in range(num_ees)], select_var_vals)
    else:
        self.gen_add_code_line("const int ee_jid = " + str(anchors[0]) + ";")
    self.gen_add_code_line("T R20 = s_Xworld[16*ee_jid + 2];")
    self.gen_add_code_line("T R21 = s_Xworld[16*ee_jid + 6];")
    self.gen_add_code_line("T R22 = s_Xworld[16*ee_jid + 10];")
    self.gen_add_code_line("T R10 = s_Xworld[16*ee_jid + 1];")
    self.gen_add_code_line("T R00 = s_Xworld[16*ee_jid + 0];")
    self.gen_add_code_line("T cp_term = sqrt(R22*R22 + R21*R21);")
    self.gen_add_code_line("T yaw = atan2(R10, R00);")
    self.gen_add_code_line("T pitch = atan2(-R20, cp_term);")
    self.gen_add_code_line("s_E_sc[4*ee + 0] = cos(yaw);")
    self.gen_add_code_line("s_E_sc[4*ee + 1] = sin(yaw);")
    self.gen_add_code_line("s_E_sc[4*ee + 2] = cos(pitch);")
    self.gen_add_code_line("s_E_sc[4*ee + 3] = sin(pitch);")
    self.gen_add_end_control_flow()
    self.gen_add_sync()

    # ===== Step 4: emit s_end_effector_pose_gradient = [J_v; E^-1 * J_w] from s_Sworld and s_Xworld =====
    # J_w[:, vi] = skew_inv(S_world[:3, :3]) = (axw_x, axw_y, axw_z) (the angular axis)
    #   In our column-major layout: skew[2,1] = wx -> s_Sworld[base+6]; skew[0,2] = wy -> s_Sworld[base+8]; skew[1,0] = wz -> s_Sworld[base+1].
    # J_v[:, vi] = (S_world @ p_ee_world)[:3] - hmm actually J_v[:, vi] = dM[vi][:3, 3] = (S_world @ X_ee)[:3, 3]
    #   = S_world[:3, :3] @ X_ee[:3, 3] + S_world[:3, 3]
    #   = [w_world]_x @ p_ee_world + (pj_world x w_world)     (revolute)
    #   = ω × p_ee_world - ω × pj_world = ω × (p_ee - pj)     ✓ matches the gradient inner
    #   = 0                       + Rw_a @ ax_local           (prismatic)
    # Use the latter expansion directly.
    self.gen_add_code_line("//")
    self.gen_add_code_line("// Step 4: write s_end_effector_pose_gradient = [J_v ; E(rpy)^-1 * J_w] from S_world")
    self.gen_add_code_line("//")
    self.gen_add_parallel_loop("ind", str(6 * nv * num_ees))
    self.gen_add_code_line("int row = ind % 6; int rem = ind / 6; int vi = rem % " + str(nv) + "; int ee = rem / " + str(nv) + ";")
    self.gen_add_code_line("int s_base = 16 * (ee * " + str(nv) + " + vi);")
    # angular axis components (top-left skew of S_world, read out)
    self.gen_add_code_line("T wx = s_Sworld[s_base + 6];   // S[2,1]")
    self.gen_add_code_line("T wy = s_Sworld[s_base + 8];   // S[0,2]")
    self.gen_add_code_line("T wz = s_Sworld[s_base + 1];   // S[1,0]")
    # Column 3 of S_world (the "translation" piece in the world-frame generator)
    self.gen_add_code_line("T s03 = s_Sworld[s_base + 12]; // S[0,3]")
    self.gen_add_code_line("T s13 = s_Sworld[s_base + 13]; // S[1,3]")
    self.gen_add_code_line("T s23 = s_Sworld[s_base + 14]; // S[2,3]")
    self.gen_add_code_line("if (row < 3) {", True)
    # J_v = S[:3,:3] @ p_ee + S[:3,3]
    # = ([w]_x @ p_ee) + s03/13/23
    # Compose ee_anchor index for current ee via select (one of `anchors`)
    if num_ees > 1:
        sel_vals = [("int", "ee_jid", [str(a) for a in anchors])]
        self.gen_add_multi_threaded_select("ee", "<", [str(i+1) for i in range(num_ees)], sel_vals)
    else:
        self.gen_add_code_line("const int ee_jid = " + str(anchors[0]) + ";")
    self.gen_add_code_line("T pex = s_Xworld[16*ee_jid + 12]; T pey = s_Xworld[16*ee_jid + 13]; T pez = s_Xworld[16*ee_jid + 14];")
    # [w]_x @ pe + s_col3
    self.gen_add_code_line("T Jv0 = (wy*pez - wz*pey) + s03;")
    self.gen_add_code_line("T Jv1 = (wz*pex - wx*pez) + s13;")
    self.gen_add_code_line("T Jv2 = (wx*pey - wy*pex) + s23;")
    self.gen_add_code_line("T outv;")
    self.gen_add_code_line("if (row == 0) outv = Jv0; else if (row == 1) outv = Jv1; else outv = Jv2;")
    self.gen_add_code_line("s_end_effector_pose_gradient[ind] = outv;")
    self.gen_add_end_control_flow()
    self.gen_add_code_line("else {", True)
    # rows 3..5 = E^-1 @ J_w with closed form (same as gradient inner)
    self.gen_add_code_line("T cy = s_E_sc[4*ee + 0]; T sy = s_E_sc[4*ee + 1]; T cp = s_E_sc[4*ee + 2]; T sp = s_E_sc[4*ee + 3];")
    self.gen_add_code_line("T outv;")
    self.gen_add_code_line("if (row == 3) { outv = (cy*wx + sy*wy) / cp; }")
    self.gen_add_code_line("else if (row == 4) { outv = -sy*wx + cy*wy; }")
    self.gen_add_code_line("else { outv = (sp / cp) * (cy*wx + sy*wy) + wz; }")
    self.gen_add_code_line("s_end_effector_pose_gradient[ind] = outv;")
    self.gen_add_end_control_flow()
    self.gen_add_end_control_flow()
    self.gen_add_sync()

    # ===== Step 5: per (ee, i, j) pair compute d2M, extract H_xyz + d2R_R^T =====
    # Strategy:
    #  - For each ee, iterate over all (i, j) DOF pairs with i, j on the chain.
    #  - Determine chain ordering: a (i's chain pos) vs b (j's chain pos).
    #  - If a < b: d2M[i,j] = S_i_world * S_j_world * X_ee
    #    H_xyz: column 3 of d2M.
    #    d2R_R^T: top-left 3x3 of (S_i_world * S_j_world).
    #  - If a > b: swap (proximal/distal).
    #  - If a == b (intra-joint, only floating base): handle separately.
    #  - We write d2M values into s_end_effector_pose_hessian. We'll do the rpy rows in a second
    #    pass once H_w / E^-1 are known.
    #
    # We pre-zero s_end_effector_pose_hessian (covers out-of-chain entries and is also needed
    # because the intra-joint case only writes the unique (i, j) ordered pair).
    self.gen_add_code_line("//")
    self.gen_add_code_line("// Step 5a: zero the full end_effector_pose_hessian output (out-of-chain pairs stay zero)")
    self.gen_add_code_line("//")
    self.gen_add_parallel_loop("ind", str(6 * nv * nv * num_ees))
    self.gen_add_code_line("s_end_effector_pose_hessian[ind] = static_cast<T>(0);")
    self.gen_add_end_control_flow()
    self.gen_add_sync()

    # Step 5b: per-pair d2M -> H_xyz + (temporarily, into end_effector_pose_hessian rpy rows) d2R_R^T
    # We use the rpy rows (c=3,4,5 of s_end_effector_pose_hessian) as a SCRATCH BUFFER for d2R_R^T's
    # skew axis (a 3-vector per pair). Specifically we write the WORLD-ANGULAR
    # kinematic Hessian H_w[:, i, j] into rows 3,4,5 here, and overwrite with
    # rpy in Step 6 via the chain rule.
    self.gen_add_code_line("//")
    self.gen_add_code_line("// Step 5b: per (ee, i, j) pair: d2M = S_i_world * S_j_world * X_ee (for a < b)")
    self.gen_add_code_line("//   or B_world * X_ee (intra-joint). Writes H_xyz to rows 0..2 and (temporarily)")
    self.gen_add_code_line("//   the world-angular Hessian H_w to rows 3..5; the rpy chain rule in Step 6")
    self.gen_add_code_line("//   then overwrites rows 3..5 with the proper H_rpy.")
    self.gen_add_code_line("//")
    # Emit per-ee per-pair code. For each ee, we have len(chain_dofs)^2 pairs.
    # Each pair fires once with explicit constants (chain_pos, S_col, joint_jid).
    #
    # PARALLELIZATION: each output cell (ee, vi, vj) writes a DISJOINT 6-element
    # region s_end_effector_pose_hessian[ee*6nv2 + c*nv2 + vi*nv + vj] (pre-zeroed in Step 5a, with
    # a sync, and no read-after-write between cells), so the cells are fully
    # independent. We collect one deferred-emit closure per cell, then dispatch
    # them one-thread-per-cell via a single block-parallel loop over a flat cell
    # index (`d2m_cell`). Each closure emits the SAME per-cell scalar arithmetic
    # as the old thread-0 serial path -> bit-identical output. For a MIMIC v-slot
    # pair the whole block-pair SUM (the reduction into that one cell) lives in a
    # SINGLE closure (one thread owns the cell and sums its block-pairs serially),
    # so the reduction is never split across threads.
    HAS_MIMIC = self.robot_has_mimic_joints()
    # Cross-joint cells (the overwhelming majority for big robots: every (vi, vj)
    # pair whose proximal/distal joints are DISTINCT chain joints) all run the
    # SAME per-cell scalar arithmetic, differing ONLY in six integer offsets
    # (prox/dist/si/sj s_Sworld bases, the ee world-transform base, and the output
    # base). The legacy emit inlined that arithmetic once per cell behind an
    # `if (d2m_cell == k)` ladder -> O(n_cross) copies of a ~90-line body, which is
    # the nvcc compile-time / code-size blow-up on h1_2/big-floating. We instead
    # BAKE the six offsets into a per-cell `int` table and emit the cross body
    # ONCE, reading its offsets from the table indexed by `d2m_cell`. Output is
    # byte-identical (same scalar ops, just offsets sourced from a table row).
    #
    # The few same-joint (intra-multi-DoF, floating-base only) and mimic
    # block-pair cells keep their explicit per-cell bodies (their arithmetic SHAPE
    # differs per cell: rev/prism axis literals, variable-length alpha sums), but
    # they are a small minority, so the `if==k` ladder over THEM stays cheap.
    cross_table = []     # list of (prox_base, dist_base, si_base, sj_base, pee_base, out_base)
    cell_emitters = []   # list of (comment_str, emit_callable) — one per NON-cross cell
    for ee_idx in range(num_ees):
        ee_jid = anchors[ee_idx]
        chain_dofs = per_ee_dof_info[ee_idx]
        chain_jids = chains[ee_idx]
        intra_pairs = intra_joint_pairs_per_ee[ee_idx]
        # Index by (vi_a, vi_b) for quick lookup
        intra_pair_lookup = {(p[0], p[1]): p for p in intra_pairs}
        intra_pair_lookup.update({(p[1], p[0]): p for p in intra_pairs})

        # Group chain DOFs by v-slot. A MIMIC joint folds into its target's slot,
        # so a v-slot can collect MULTIPLE chain blocks (the target + each mimic,
        # e.g. the h1_2 thumb: proximal + 2 mimics -> 3 blocks on one slot). The
        # Hessian cell for any (vi, vj) where either slot is multi-block must SUM
        # over every block-pair (RBDReference convention); a single-writer emit
        # would last-writer-win and drop all but one term (the exact-zero bug).
        vi_to_blocks = {}
        for d in chain_dofs:
            vi_to_blocks.setdefault(d["vi"], []).append({
                "chain_pos": d["chain_pos"], "jid": d["joint_jid"],
                "alpha": (self._alpha_for_jid(d["joint_jid"]) if HAS_MIMIC else 1.0),
                "ang": d["ang"], "lin": d["lin"], "revolute": d["revolute"],
            })
        for vlist in vi_to_blocks.values():
            vlist.sort(key=lambda b: b["chain_pos"])
        multi_slots = {vi for vi, blks in vi_to_blocks.items() if len(blks) > 1}

        # Pair iteration: per (vi_i, vi_j), with i,j enumerated over chain DOFs.
        for di in chain_dofs:
            vi = di["vi"]
            for dj in chain_dofs:
                vj = dj["vi"]
                # MIMIC multi-block routing: if either slot has >1 chain block,
                # emit the alpha-weighted block-pair SUM exactly ONCE per ordered
                # (vi, vj) pair (skip the duplicate chain-DOF iterations that map
                # to the same slot pair). Non-mimic slots are always singletons
                # so this branch is never taken -> output byte-identical.
                if vi in multi_slots or vj in multi_slots:
                    # Emit the slot pair ONCE: only when this (di, dj) is the
                    # first chain block of slot vi AND the first of slot vj.
                    if (di["chain_pos"] != vi_to_blocks[vi][0]["chain_pos"]
                            or dj["chain_pos"] != vi_to_blocks[vj][0]["chain_pos"]):
                        continue  # already emitted for this (vi, vj) slot pair
                    si_base = 16 * (ee_idx * nv + vi)
                    sj_base = 16 * (ee_idx * nv + vj)
                    cell_emitters.append((
                        "// ee=" + str(ee_idx) + " MIMIC pair (vi=" + str(vi) +
                        ", vj=" + str(vj) + ")",
                        (lambda ee_idx=ee_idx, ee_jid=ee_jid, vi=vi, vj=vj,
                                bi=vi_to_blocks[vi], bj=vi_to_blocks[vj],
                                sib=si_base, sjb=sj_base:
                            _emit_d2M_mimic_vslot_pair_block(
                                self, ee_idx, ee_jid, vi, vj, nv, num_ees,
                                bi, bj, sib, sjb))))
                    continue
                # Determine ordering: a = di['chain_pos'], b = dj['chain_pos']
                a = di["chain_pos"]; b = dj["chain_pos"]
                # H index: c*nv*nv + i*nv + j  with i = "outer" (Hessian row j_v),
                # j = "inner" (Hessian col i_v).  Our chosen layout from the FD
                # path was: idx = e*6*nv*nv + c*nv*nv + j_outer*nv + i_inner.
                # The C-order (6, nv, nv) tensor here is H[c, i_h, j_h] with the
                # convention that mid axis = i, last axis = j. Match the FD path
                # (which writes h_idx = ... + j*nv + i_fd) so the public layout
                # is identical: index = e*6*nv*nv + c*nv*nv + vi*nv + vj.
                # Per-pair scoped block — local names don't collide across pairs.
                si_base = 16 * (ee_idx * nv + vi)
                sj_base = 16 * (ee_idx * nv + vj)
                comment = ("// ee=" + str(ee_idx) + " pair (vi=" + str(vi) + ", vj=" + str(vj) +
                           ", a=" + str(a) + ", b=" + str(b) + ")")
                # We need the top-left 3x3 product (S_prox @ S_dist)[:3,:3] and
                # the column-3 expansion (S_prox @ S_dist @ X_ee)[:3, 3].
                # Compute it for the proximal/distal ordering.
                if a == b:
                    # Same chain joint. For single-DOF joints, B_local = 0 (revolute or
                    # prismatic intra-pair doesn't exist except for the diagonal where
                    # B = A_x^2 for revolute, 0 for prismatic). For multi-DOF (floating
                    # base) joints, use the closed form B_world below.
                    # Diagonal vi == vj case (always present); off-diagonal same-joint
                    # pairs only exist for multi-DOF (floating base) joints.
                    if vi == vj:
                        sj_block = di  # same as dj on the diagonal
                    elif (vi, vj) in intra_pair_lookup:
                        sj_block = dj
                    else:
                        # Shouldn't happen (a == b but DOFs not in same joint).
                        # Out-of-chain cell stays at its Step-5a zero -> emit nothing.
                        continue
                    emit = (lambda di=di, dj=sj_block, ee_idx=ee_idx, ee_jid=ee_jid, vi=vi, vj=vj,
                                   jid=chain_jids[a], sib=si_base, sjb=sj_base:
                                _emit_d2M_same_joint_block(self, di, dj,
                                                           ee_idx, ee_jid, vi, vj, nv, num_ees,
                                                           jid, sib, sjb))
                else:
                    # Different chain joints: proximal = smaller chain_pos.
                    # This is the data-driven cross-joint path: bake the six
                    # offsets into cross_table and let the single shared body read
                    # them by `d2m_cell`. (Same final values as the inlined body.)
                    prox_base = si_base if a < b else sj_base
                    dist_base = sj_base if a < b else si_base
                    out_base = ee_idx * 6 * nv * nv + (vi * nv + vj)
                    cross_table.append((prox_base, dist_base, si_base, sj_base,
                                        16 * ee_jid, out_base))
                    continue
                def _emit_nonmimic_cell(ee_jid=ee_jid, _emit=emit):
                    # Read X_ee column 3 (p_ee); helpers consume pex/pey/pez.
                    self.gen_add_code_line("T pex = s_Xworld[" + str(16*ee_jid + 12) + "];")
                    self.gen_add_code_line("T pey = s_Xworld[" + str(16*ee_jid + 13) + "];")
                    self.gen_add_code_line("T pez = s_Xworld[" + str(16*ee_jid + 14) + "];")
                    _emit()
                cell_emitters.append((comment, _emit_nonmimic_cell))

    # Dispatch all cells one-thread-per-cell via a single block-parallel loop.
    # Flat cell index layout: [0, n_cross) cross-joint cells (one shared
    # table-driven body), then [n_cross, n_cells) the explicit non-cross bodies
    # (same-joint + mimic) behind the residual `if==k` ladder. Same total cell
    # count and launch geometry as before; the cross bulk is now a single body.
    n_cross = len(cross_table)
    n_noncross = len(cell_emitters)
    n_cells = n_cross + n_noncross
    self.gen_add_parallel_loop("d2m_cell", str(n_cells))
    if n_cross > 0:
        # Baked offset table for the cross-joint cells. Row layout (6 ints):
        #   [0]=prox_base [1]=dist_base [2]=si_base [3]=sj_base
        #   [4]=pee_base (s_Xworld base of the ee transform; +12/13/14 = p_ee)
        #   [5]=out_base (ee*6*nv*nv + vi*nv + vj; +c*nv*nv selects the 6 rows)
        # Flattened row-major so a single static array drives every cross cell.
        flat = []
        for row in cross_table:
            flat.extend(int(x) for x in row)
        table_literal = ", ".join(str(x) for x in flat)
        self.gen_add_code_line("// Cross-joint cells (vast majority): one shared body driven by a baked")
        self.gen_add_code_line("// per-cell offset table; collapses the old O(n_cross) if==k ladder.")
        self.gen_add_code_line("static const int s_d2ee_cross_tab[" + str(len(flat)) + "] = {" + table_literal + "};")
        self.gen_add_code_line("if (d2m_cell < " + str(n_cross) + ") {", True)
        _emit_d2M_cross_joint_table_body(self, nv)
        self.gen_add_end_control_flow()
    for k, (comment, emit) in enumerate(cell_emitters):
        self.gen_add_code_line(comment)
        self.gen_add_code_line("if (d2m_cell == " + str(n_cross + k) + ") {", True)
        emit()
        self.gen_add_end_control_flow()
    self.gen_add_end_control_flow()
    self.gen_add_sync()

    # ===== Step 6: rpy chain rule on rows 3..5 =====
    # Currently rows 3..5 hold H_w[:, i, j]. We want H_rpy[:, i, j] =
    # (dEinv/dv_j) * J_w[:, i] + Einv * H_w[:, i, j], with closed-form Einv
    # and dE/drpy. Recall the gradient inner already wrote drpy/dv = Einv*J_w
    # to s_end_effector_pose_gradient rows 3..5: but we need that for each (ee, vj).
    #
    # NOTE on race-safety: each thread owns a single (ee, vi, vj) cell and reads
    # all 3 components of H_w[:, vi, vj] from rows 3..5 of s_end_effector_pose_hessian BEFORE
    # writing any rpy. We then write all 3 rpy components. There's no read-after-
    # write hazard within the parallel pass because each (ee, vi, vj) is owned by
    # exactly one thread (no thread reads another thread's writes).
    self.gen_add_code_line("//")
    self.gen_add_code_line("// Step 6: rpy chain rule -- replace rows 3..5 with H_rpy = dEinv/dv_j @ J_w_i + Einv @ H_w")
    self.gen_add_code_line("//")
    self.gen_add_parallel_loop("ind", str(nv * nv * num_ees))
    self.gen_add_code_line("int vj = ind % " + str(nv) + ";")
    self.gen_add_code_line("int vi = (ind / " + str(nv) + ") % " + str(nv) + ";")
    self.gen_add_code_line("int ee = ind / " + str(nv * nv) + ";")
    self.gen_add_code_line("T cy = s_E_sc[4*ee + 0]; T sy = s_E_sc[4*ee + 1]; T cp = s_E_sc[4*ee + 2]; T sp = s_E_sc[4*ee + 3];")
    # Einv (closed form):
    # E = [[cy*cp, -sy, 0], [sy*cp, cy, 0], [-sp, 0, 1]]
    # det(E) = cp; Einv = (1/cp) * [[cy, sy, 0], [-sy*cp, cy*cp, 0], [cy*sp, sy*sp, cp]]
    self.gen_add_code_line("T inv_cp = static_cast<T>(1) / cp;")
    self.gen_add_code_line("// Einv (3x3); row 0 = roll, row 1 = pitch, row 2 = yaw")
    self.gen_add_code_line("T Einv00 = cy * inv_cp;       T Einv01 = sy * inv_cp;       T Einv02 = static_cast<T>(0);")
    self.gen_add_code_line("T Einv10 = -sy;               T Einv11 = cy;                T Einv12 = static_cast<T>(0);")
    self.gen_add_code_line("T Einv20 = cy * sp * inv_cp;  T Einv21 = sy * sp * inv_cp;  T Einv22 = static_cast<T>(1);")
    # drpy_j (from s_end_effector_pose_gradient rows 3,4,5)
    self.gen_add_code_line("int dee_base_j = ee * " + str(6 * nv) + " + 6 * vj;")
    self.gen_add_code_line("T drpy_j0 = s_end_effector_pose_gradient[dee_base_j + 3];")
    self.gen_add_code_line("T drpy_j1 = s_end_effector_pose_gradient[dee_base_j + 4];")
    self.gen_add_code_line("T drpy_j2 = s_end_effector_pose_gradient[dee_base_j + 5];")
    # drpy_i (also from s_end_effector_pose_gradient)
    self.gen_add_code_line("int dee_base_i = ee * " + str(6 * nv) + " + 6 * vi;")
    self.gen_add_code_line("T drpy_i0 = s_end_effector_pose_gradient[dee_base_i + 3];")
    self.gen_add_code_line("T drpy_i1 = s_end_effector_pose_gradient[dee_base_i + 4];")
    self.gen_add_code_line("T drpy_i2 = s_end_effector_pose_gradient[dee_base_i + 5];")
    # READ ALL 3 H_w components BEFORE any rpy writes (critical for correctness:
    # we will overwrite rows 3..5 below; if we read after writing the race would
    # silently corrupt the other two components in this thread's row).
    self.gen_add_code_line("int hw_base = ee * " + str(6 * nv * nv) + " + 3 * " + str(nv * nv) + " + vi * " + str(nv) + " + vj;")
    self.gen_add_code_line("T Hw_x = s_end_effector_pose_hessian[hw_base + 0 * " + str(nv * nv) + "];")
    self.gen_add_code_line("T Hw_y = s_end_effector_pose_hessian[hw_base + 1 * " + str(nv * nv) + "];")
    self.gen_add_code_line("T Hw_z = s_end_effector_pose_hessian[hw_base + 2 * " + str(nv * nv) + "];")
    # dE_total = sum_k dE/drpy_k * drpy_j[k]:
    # dE/droll = 0 (irrelevant; drops out)
    # dE/dpitch = [[-cy*sp, 0, 0], [-sy*sp, 0, 0], [-cp, 0, 0]]
    # dE/dyaw   = [[-sy*cp, -cy, 0], [cy*cp, -sy, 0], [0, 0, 0]]
    # dE_total[r,c] = drpy_j1 * dE/dpitch[r,c] + drpy_j2 * dE/dyaw[r,c]
    self.gen_add_code_line("T dE00 = -cy * sp * drpy_j1 - sy * cp * drpy_j2;")
    self.gen_add_code_line("T dE01 = -cy * drpy_j2;")
    self.gen_add_code_line("// dE02 = 0")
    self.gen_add_code_line("T dE10 = -sy * sp * drpy_j1 + cy * cp * drpy_j2;")
    self.gen_add_code_line("T dE11 = -sy * drpy_j2;")
    self.gen_add_code_line("// dE12 = 0")
    self.gen_add_code_line("T dE20 = -cp * drpy_j1;")
    self.gen_add_code_line("// dE21 = dE22 = 0")
    # For each row r of the rpy Hessian: H_rpy_r = -(U_r @ drpy_i) + (Einv_r @ Hw)
    # where U_r = Einv[r,:] @ dE_total (a 3-vector; only U_r[0] and U_r[1] are nonzero
    # because dE_total's columns 2, and entries with r==2,c=1, etc., are zero).
    # We just unroll all three rows.
    self.gen_add_code_line("// row 0 (roll): U_0 = Einv[0,:] @ dE_total")
    self.gen_add_code_line("T U0_0 = Einv00 * dE00 + Einv01 * dE10 + Einv02 * dE20;")
    self.gen_add_code_line("T U0_1 = Einv00 * dE01 + Einv01 * dE11;")
    self.gen_add_code_line("// row 1 (pitch)")
    self.gen_add_code_line("T U1_0 = Einv10 * dE00 + Einv11 * dE10 + Einv12 * dE20;")
    self.gen_add_code_line("T U1_1 = Einv10 * dE01 + Einv11 * dE11;")
    self.gen_add_code_line("// row 2 (yaw)")
    self.gen_add_code_line("T U2_0 = Einv20 * dE00 + Einv21 * dE10 + Einv22 * dE20;")
    self.gen_add_code_line("T U2_1 = Einv20 * dE01 + Einv21 * dE11;")
    # H_rpy[r] = -(U_r[0]*drpy_i0 + U_r[1]*drpy_i1 + U_r[2]*drpy_i2) + Einv[r,:] @ Hw
    # U_r[2] = Einv[r,0]*dE[0,2] + Einv[r,1]*dE[1,2] + Einv[r,2]*dE[2,2] = 0 (all zero)
    self.gen_add_code_line("T H_rpy_0 = -(U0_0 * drpy_i0 + U0_1 * drpy_i1) + (Einv00 * Hw_x + Einv01 * Hw_y + Einv02 * Hw_z);")
    self.gen_add_code_line("T H_rpy_1 = -(U1_0 * drpy_i0 + U1_1 * drpy_i1) + (Einv10 * Hw_x + Einv11 * Hw_y + Einv12 * Hw_z);")
    self.gen_add_code_line("T H_rpy_2 = -(U2_0 * drpy_i0 + U2_1 * drpy_i1) + (Einv20 * Hw_x + Einv21 * Hw_y + Einv22 * Hw_z);")
    # Write all three rpy components (rows 3, 4, 5)
    self.gen_add_code_line("int out_base = ee * " + str(6 * nv * nv) + " + 3 * " + str(nv * nv) + " + vi * " + str(nv) + " + vj;")
    self.gen_add_code_line("s_end_effector_pose_hessian[out_base + 0 * " + str(nv * nv) + "] = H_rpy_0;")
    self.gen_add_code_line("s_end_effector_pose_hessian[out_base + 1 * " + str(nv * nv) + "] = H_rpy_1;")
    self.gen_add_code_line("s_end_effector_pose_hessian[out_base + 2 * " + str(nv * nv) + "] = H_rpy_2;")
    self.gen_add_end_control_flow()
    self.gen_add_sync()

    # ===== Step 7: symmetrize rows 3..5 over (i, j) =====
    # H_xyz is symmetric by construction. H_rpy is computed asymmetrically (the
    # dEinv_dvj branch only sees the j-direction), so we average with the
    # transpose. The Python reference does the same final symmetrization.
    self.gen_add_code_line("//")
    self.gen_add_code_line("// Step 7: symmetrize H_rpy (rows 3..5) over (i, j)")
    self.gen_add_code_line("//")
    self.gen_add_parallel_loop("ind", str(num_ees * 3 * nv * nv))
    self.gen_add_code_line("int e = ind / " + str(3 * nv * nv) + ";")
    self.gen_add_code_line("int cji = ind % " + str(3 * nv * nv) + ";")
    self.gen_add_code_line("int rrow = cji / " + str(nv * nv) + ";   // 0..2 -> c = 3 + rrow")
    self.gen_add_code_line("int ji = cji % " + str(nv * nv) + ";")
    self.gen_add_code_line("int i = ji / " + str(nv) + "; int j = ji % " + str(nv) + ";")
    self.gen_add_code_line("if (i <= j) {", True)
    self.gen_add_code_line("int c = 3 + rrow;")
    self.gen_add_code_line("int idx_ij = e * " + str(6 * nv * nv) + " + c * " + str(nv * nv) + " + i * " + str(nv) + " + j;")
    self.gen_add_code_line("int idx_ji = e * " + str(6 * nv * nv) + " + c * " + str(nv * nv) + " + j * " + str(nv) + " + i;")
    self.gen_add_code_line("T avg = static_cast<T>(0.5) * (s_end_effector_pose_hessian[idx_ij] + s_end_effector_pose_hessian[idx_ji]);")
    self.gen_add_code_line("s_end_effector_pose_hessian[idx_ij] = avg;")
    self.gen_add_code_line("if (i != j) { s_end_effector_pose_hessian[idx_ji] = avg; }")
    self.gen_add_end_control_flow()
    self.gen_add_end_control_flow()
    self.gen_add_sync()

    self.gen_add_end_function()


def _emit_d2M_cross_joint_block(self, prox_base, dist_base,
                                 ee_idx, vi, vj, nv, num_ees, si_base, sj_base):
    """Emit explicit per-pair scalar code for d2M = S_prox * S_dist * X_ee in a
    cross-joint pair (chain ordering a < b after prox/dist resolution).

    Writes:
      - rows 0..2 of s_end_effector_pose_hessian[(ee, vi, vj)]   = (d2M @ ee_offset)[:3] with
        ee_offset = [0, 0, 0, 1] -> just column 3 of d2M.
      - rows 3..5 of s_end_effector_pose_hessian[(ee, vi, vj)]   = H_w[:, vi, vj] =
        skew_inv(d2R @ R_chain^T - [Jw_i]_x @ [Jw_j]_x)
        with d2R_R^T = (S_prox * S_dist)[:3, :3]
        and [Jw_i]_x = S_i_world[:3, :3]
        and [Jw_j]_x = S_j_world[:3, :3]

    Caller has already declared and set: pex, pey, pez (X_ee column 3 in world).
    """
    # Read S_prox and S_dist top 3 rows (column-major: S[r + 4*c]).
    # S has zero bottom row so we only need rows 0..2.
    self.gen_add_code_line("// Read S_prox (chain proximal) rows 0..2")
    for c in range(4):
        for r in range(3):
            self.gen_add_code_line("T P" + str(r) + str(c) + " = s_Sworld[" + str(prox_base + r + 4*c) + "];")
    self.gen_add_code_line("// Read S_dist (chain distal) rows 0..2")
    for c in range(4):
        for r in range(3):
            self.gen_add_code_line("T D" + str(r) + str(c) + " = s_Sworld[" + str(dist_base + r + 4*c) + "];")
    # Compute M = S_prox * S_dist  (4x4, but bottom row of result is 0).
    # M[r, c] = sum_k P[r, k] * D[k, c]; since P[3, :] = 0 and D[3, :] = 0,
    # we only need top 3 rows of M, and for each (r, c) we sum k = 0..2.
    # M[r, c] = P[r, 0]*D[0, c] + P[r, 1]*D[1, c] + P[r, 2]*D[2, c]
    self.gen_add_code_line("// M = S_prox * S_dist (top 3 rows, all 4 cols)")
    for r in range(3):
        for c in range(4):
            self.gen_add_code_line(
                "T M" + str(r) + str(c) + " = P" + str(r) + "0*D0" + str(c) +
                " + P" + str(r) + "1*D1" + str(c) +
                " + P" + str(r) + "2*D2" + str(c) + ";")
    # d2M = M * X_ee. Top-left 3x3 of d2M = M[:3, :3] * X_ee[:3, :3] (since
    # M[:3, 3] only enters column 3 of d2M times X_ee[3, :3] which is 0).
    # Column 3 of d2M[:3] = M[:3, :3] * X_ee[:3, 3] + M[:3, 3] * X_ee[3, 3]
    #                     = M[:3, :3] * p_ee + M[:3, 3].
    # We only need column 3 (for H_xyz) — the top-left 3x3 of d2R_R^T cancels
    # to M[:3, :3] anyway (since R_chain^T = X_ee[:3,:3]^T cancels with the
    # X_ee[:3, :3] factor on the right).
    self.gen_add_code_line("// H_xyz[:, vi, vj] = (S_prox*S_dist*X_ee)[:3, 3] = M[:3,:3] * p_ee + M[:3, 3]")
    self.gen_add_code_line("T Hxyz_x = M00*pex + M01*pey + M02*pez + M03;")
    self.gen_add_code_line("T Hxyz_y = M10*pex + M11*pey + M12*pez + M13;")
    self.gen_add_code_line("T Hxyz_z = M20*pex + M21*pey + M22*pez + M23;")
    # d2R @ R_chain^T = M[:3, :3]. Need [Jw_i]_x and [Jw_j]_x = S_i_world[:3,:3] and
    # S_j_world[:3,:3]. Note: the "prox"/"dist" assignment may have swapped i↔j,
    # but H_w is the same value either way (the formula skew_inv(d2R_R^T - [Jwi]_x [Jwj]_x)
    # uses the original i,j indexing). So we read S_i_world, S_j_world directly via si_base/sj_base.
    self.gen_add_code_line("// Read [Jw_i]_x (S_i_world top-left)")
    for c in range(3):
        for r in range(3):
            self.gen_add_code_line("T Si" + str(r) + str(c) + " = s_Sworld[" + str(si_base + r + 4*c) + "];")
    self.gen_add_code_line("// Read [Jw_j]_x (S_j_world top-left)")
    for c in range(3):
        for r in range(3):
            self.gen_add_code_line("T Sj" + str(r) + str(c) + " = s_Sworld[" + str(sj_base + r + 4*c) + "];")
    # Compute SiSj = S_i_world[:3,:3] * S_j_world[:3,:3]
    self.gen_add_code_line("// SiSj = [Jw_i]_x @ [Jw_j]_x")
    for r in range(3):
        for c in range(3):
            self.gen_add_code_line(
                "T SiSj" + str(r) + str(c) + " = Si" + str(r) + "0*Sj0" + str(c) +
                " + Si" + str(r) + "1*Sj1" + str(c) +
                " + Si" + str(r) + "2*Sj2" + str(c) + ";")
    # H_w_skew = M[:3,:3] - SiSj  (note: M[:3,:3] is d2R @ R_chain^T = top-left 3x3 of S_prox*S_dist).
    # skew_inv(A) = 0.5 * (A[2,1] - A[1,2], A[0,2] - A[2,0], A[1,0] - A[0,1])
    self.gen_add_code_line("// H_w[:, vi, vj] = skew_inv(M[:3,:3] - SiSj)")
    self.gen_add_code_line("T HW_x = static_cast<T>(0.5) * ((M21 - SiSj21) - (M12 - SiSj12));")
    self.gen_add_code_line("T HW_y = static_cast<T>(0.5) * ((M02 - SiSj02) - (M20 - SiSj20));")
    self.gen_add_code_line("T HW_z = static_cast<T>(0.5) * ((M10 - SiSj10) - (M01 - SiSj01));")
    # Write into s_end_effector_pose_hessian: idx = ee*6*nv*nv + c*nv*nv + vi*nv + vj
    base = "(" + str(ee_idx * 6 * nv * nv) + " + " + str(vi * nv + vj) + ")"
    self.gen_add_code_line("s_end_effector_pose_hessian[" + base + " + 0 * " + str(nv*nv) + "] = Hxyz_x;")
    self.gen_add_code_line("s_end_effector_pose_hessian[" + base + " + 1 * " + str(nv*nv) + "] = Hxyz_y;")
    self.gen_add_code_line("s_end_effector_pose_hessian[" + base + " + 2 * " + str(nv*nv) + "] = Hxyz_z;")
    self.gen_add_code_line("s_end_effector_pose_hessian[" + base + " + 3 * " + str(nv*nv) + "] = HW_x;")
    self.gen_add_code_line("s_end_effector_pose_hessian[" + base + " + 4 * " + str(nv*nv) + "] = HW_y;")
    self.gen_add_code_line("s_end_effector_pose_hessian[" + base + " + 5 * " + str(nv*nv) + "] = HW_z;")


def _emit_d2M_cross_joint_table_body(self, nv):
    """Data-driven cross-joint d2M body shared by EVERY cross-joint cell.

    Numerically identical to `_emit_d2M_cross_joint_block` (same scalar ops, same
    float association), except the six per-cell offsets are read at runtime from
    the baked `s_d2ee_cross_tab` table indexed by the loop counter `d2m_cell`,
    instead of being inlined as compile-time constants. This replaces the
    O(n_cross) `if==k` ladder of inlined bodies with ONE body -> the nvcc
    compile-time + code-size win at identical Hessian values.

    Table row (6 ints, see emitter): prox_base, dist_base, si_base, sj_base,
    pee_base, out_base. s_Sworld is column-major 4x4 per DOF (S[r + 4*c]); the ee
    world transform's p_ee is s_Xworld[pee_base + 12/13/14]; the output cell base
    is out_base with row stride nv*nv (rows 0..2 = H_xyz, rows 3..5 = H_w temp).
    """
    nn = nv * nv
    # Load the six per-cell offsets for this thread's cross cell.
    self.gen_add_code_line("const int *row = &s_d2ee_cross_tab[6 * d2m_cell];")
    self.gen_add_code_line("int prox_base = row[0]; int dist_base = row[1];")
    self.gen_add_code_line("int si_base = row[2]; int sj_base = row[3];")
    self.gen_add_code_line("int pee_base = row[4]; int out_base = row[5];")
    # p_ee from the ee world transform (column 3).
    self.gen_add_code_line("T pex = s_Xworld[pee_base + 12];")
    self.gen_add_code_line("T pey = s_Xworld[pee_base + 13];")
    self.gen_add_code_line("T pez = s_Xworld[pee_base + 14];")
    # Read S_prox / S_dist top 3 rows (column-major). Bottom row of S is zero.
    self.gen_add_code_line("// Read S_prox (chain proximal) rows 0..2")
    for c in range(4):
        for r in range(3):
            self.gen_add_code_line("T P" + str(r) + str(c) + " = s_Sworld[prox_base + " + str(r + 4*c) + "];")
    self.gen_add_code_line("// Read S_dist (chain distal) rows 0..2")
    for c in range(4):
        for r in range(3):
            self.gen_add_code_line("T D" + str(r) + str(c) + " = s_Sworld[dist_base + " + str(r + 4*c) + "];")
    # M = S_prox * S_dist (top 3 rows, all 4 cols).
    self.gen_add_code_line("// M = S_prox * S_dist (top 3 rows, all 4 cols)")
    for r in range(3):
        for c in range(4):
            self.gen_add_code_line(
                "T M" + str(r) + str(c) + " = P" + str(r) + "0*D0" + str(c) +
                " + P" + str(r) + "1*D1" + str(c) +
                " + P" + str(r) + "2*D2" + str(c) + ";")
    self.gen_add_code_line("// H_xyz[:, vi, vj] = (S_prox*S_dist*X_ee)[:3, 3] = M[:3,:3] * p_ee + M[:3, 3]")
    self.gen_add_code_line("T Hxyz_x = M00*pex + M01*pey + M02*pez + M03;")
    self.gen_add_code_line("T Hxyz_y = M10*pex + M11*pey + M12*pez + M13;")
    self.gen_add_code_line("T Hxyz_z = M20*pex + M21*pey + M22*pez + M23;")
    self.gen_add_code_line("// Read [Jw_i]_x (S_i_world top-left)")
    for c in range(3):
        for r in range(3):
            self.gen_add_code_line("T Si" + str(r) + str(c) + " = s_Sworld[si_base + " + str(r + 4*c) + "];")
    self.gen_add_code_line("// Read [Jw_j]_x (S_j_world top-left)")
    for c in range(3):
        for r in range(3):
            self.gen_add_code_line("T Sj" + str(r) + str(c) + " = s_Sworld[sj_base + " + str(r + 4*c) + "];")
    self.gen_add_code_line("// SiSj = [Jw_i]_x @ [Jw_j]_x")
    for r in range(3):
        for c in range(3):
            self.gen_add_code_line(
                "T SiSj" + str(r) + str(c) + " = Si" + str(r) + "0*Sj0" + str(c) +
                " + Si" + str(r) + "1*Sj1" + str(c) +
                " + Si" + str(r) + "2*Sj2" + str(c) + ";")
    self.gen_add_code_line("// H_w[:, vi, vj] = skew_inv(M[:3,:3] - SiSj)")
    self.gen_add_code_line("T HW_x = static_cast<T>(0.5) * ((M21 - SiSj21) - (M12 - SiSj12));")
    self.gen_add_code_line("T HW_y = static_cast<T>(0.5) * ((M02 - SiSj02) - (M20 - SiSj20));")
    self.gen_add_code_line("T HW_z = static_cast<T>(0.5) * ((M10 - SiSj10) - (M01 - SiSj01));")
    self.gen_add_code_line("s_end_effector_pose_hessian[out_base + 0 * " + str(nn) + "] = Hxyz_x;")
    self.gen_add_code_line("s_end_effector_pose_hessian[out_base + 1 * " + str(nn) + "] = Hxyz_y;")
    self.gen_add_code_line("s_end_effector_pose_hessian[out_base + 2 * " + str(nn) + "] = Hxyz_z;")
    self.gen_add_code_line("s_end_effector_pose_hessian[out_base + 3 * " + str(nn) + "] = HW_x;")
    self.gen_add_code_line("s_end_effector_pose_hessian[out_base + 4 * " + str(nn) + "] = HW_y;")
    self.gen_add_code_line("s_end_effector_pose_hessian[out_base + 5 * " + str(nn) + "] = HW_z;")


def _emit_d2M_same_joint_block(self, di, dj, ee_idx, ee_jid, vi, vj, nv, num_ees,
                                joint_jid, si_base, sj_base):
    """Emit per-pair code for the SAME-JOINT (a == b) case: d2M = L_a * B_local * L_a^-1 * X_ee
    where the joint at chain position a contributes its intrinsic
    second-order Lie-group term B_local (a 4x4).

    For 1-DOF revolute joints with axis a_local and vi == vj:
       B_local = [[ [a]_x^2, 0 ], [ 0, 0 ]]
    For 1-DOF prismatic: B_local = 0.
    For multi-DOF (floating base, jid=0) intra-joint pairs (c_a, c_b):
       - lin-lin: B_local = 0
       - lin-ang or ang-lin: B_local[:3, 3] = 0.5 * (ang_local x lin_local) (third col only)
       - ang-ang: B_local[:3, :3] = 0.5 * ([a]_x [b]_x + [b]_x [a]_x); col 3 = 0.

    We resolve to (c_a = di['S_col'], c_b = dj['S_col']), grab the body-frame
    ang/lin axes from the metadata, and inline emit the world-frame B (after
    L_a conjugation) -> d2M = B_world * X_ee.

    Closed-form world-frame B_world (= L_a B_local L_a^{-1}):
      Let L_a = [[Ra, pa], [0, 1]]; L_a^-1 = [[Ra^T, -Ra^T pa], [0, 1]].
      Let B_local = [[Br, Bt], [0, 0]] (top-left rotation 3x3 Br, top-right col Bt).
      Then L_a B_local = [[Ra Br, Ra Bt], [0, 0]]
           L_a B_local L_a^-1 = [[(Ra Br) Ra^T, -(Ra Br)(Ra^T pa) + (Ra Bt)], [0, 0]]
                              = [[ Ra Br Ra^T, Ra Bt - (Ra Br Ra^T) pa ], [ 0, 0 ]]
      So B_world[:3, :3] = Ra @ Br @ Ra^T
         B_world[:3, 3]  = Ra @ Bt - B_world[:3, :3] @ pa

    Then d2M = B_world * X_ee, exactly the same final step as the cross-joint
    case (apart from the B_world matrix sourcing).

    Writes:
      - rows 0..2 of s_end_effector_pose_hessian[(ee, vi, vj)] = (d2M)[:3, 3]
      - rows 3..5                                = H_w[:, vi, vj] =
        skew_inv(B_world[:3, :3] - [Jwi]_x @ [Jwj]_x)
    """
    c_a = di["S_col"]
    c_b = dj["S_col"]
    ang_a = di["ang"]; lin_a = di["lin"]; rev_a = di["revolute"]
    ang_b = dj["ang"]; lin_b = dj["lin"]; rev_b = dj["revolute"]

    # We need Ra, pa for chain joint a. Since a == b == di['chain_pos'] (and the
    # joint is joint_jid), L_a is s_Xworld[16*joint_jid].
    self.gen_add_code_line("// same-joint pair (intra-joint), joint_jid=" + str(joint_jid) +
                           " c_a=" + str(c_a) + " c_b=" + str(c_b))
    # Compute B_world[:3, :3] (Br_w) and B_world[:3, 3] (Bt_w) symbolically.
    # We branch on (rev_a, rev_b) configurations:
    #   rev-rev:   Br_local = 0.5 * ([a]_x [b]_x + [b]_x [a]_x); Bt_local = 0
    #     Br_world = Ra @ Br_local @ Ra^T;  Bt_world = -Br_world @ pa
    #     But Ra @ [a]_x @ Ra^T = [Ra a]_x = [a_w]_x (the world axis). So:
    #     Br_world = 0.5 * ([a_w]_x [b_w]_x + [b_w]_x [a_w]_x)  -- can be expressed
    #     using S_world top-left for both DOFs (already stored).
    #   rev-pris (a rev, b pris): Br_local = 0; Bt_local = 0.5 * (a × b_lin)
    #     Br_world = 0; Bt_world = 0.5 * Ra @ (a × b_lin) = 0.5 * (a_w × b_lin_w)
    #     where a_w = Ra @ a (rotational axis), b_lin_w = Ra @ b_lin (lin axis).
    #     Both a_w and b_lin_w can be read from the corresponding S_world entries.
    #   pris-rev: symmetric to rev-pris (we treat B as symmetric in c_a, c_b).
    #   pris-pris: B_local = 0 → d2M = 0 (no contribution).
    a_w = "Ra @ ang_a"  # placeholder; we emit via S_world entries
    # Get world axes from stored S entries:
    # If revolute: skew block in S_world has axis (wx, wy, wz) = (S[2,1], S[0,2], S[1,0]).
    # If prismatic: S_world[:3, 3] = R_a @ ax_local = axis_world.
    # For the "ang_a" axis we need it whether a is revolute or whether it
    # contributes only via the cross-product term. We'll compute axis_world
    # FROM s_Xworld[16*joint_jid] @ ax_local for each axis we need.
    # That keeps the code uniform and doesn't depend on extra interim variables.
    def _emit_world_axis(name, ax_local):
        # ax_world = Ra @ ax_local where Ra = top-left 3x3 of s_Xworld[16*joint_jid] (column-major)
        for r in range(3):
            terms = []
            for c in range(3):
                if abs(ax_local[c]) < 1e-15:
                    continue
                coef = "static_cast<T>(" + "{:.17g}".format(ax_local[c]) + ")"
                terms.append("s_Xworld[" + str(16*joint_jid + r + 4*c) + "] * " + coef)
            expr = " + ".join(terms) if terms else "static_cast<T>(0)"
            self.gen_add_code_line("T " + name + "_" + str(r) + " = " + expr + ";")
    self.gen_add_code_line("// Compute joint-a world frame axes of c_a and c_b body axes")
    # For Br computation we need rotational world axes for revolute DOFs.
    # For Bt (lin-ang) we need lin world axis and ang world axis. Compute both
    # always (cheap) so the branching logic below is straightforward.
    _emit_world_axis("aw", ang_a)  # ang_a in world
    _emit_world_axis("bw", ang_b)
    _emit_world_axis("alw", lin_a)  # lin_a in world
    _emit_world_axis("blw", lin_b)
    # pa = column 3 of s_Xworld[16*joint_jid]
    self.gen_add_code_line("T pax = s_Xworld[" + str(16*joint_jid + 12) + "];")
    self.gen_add_code_line("T pay = s_Xworld[" + str(16*joint_jid + 13) + "];")
    self.gen_add_code_line("T paz = s_Xworld[" + str(16*joint_jid + 14) + "];")

    if rev_a and rev_b:
        # Br_world = 0.5 * ([a_w]_x @ [b_w]_x + [b_w]_x @ [a_w]_x)
        # Using axw cross product identity:  [a]_x @ [b]_x = b @ a^T - (a . b) I
        # So 0.5 * ([a]_x [b]_x + [b]_x [a]_x) = 0.5 * (a b^T + b a^T) - (a . b) I
        # (the symmetric symmetric product of skews equals the symmetrized outer minus dot*I)
        self.gen_add_code_line("// Br_world = 0.5 * ([aw]_x [bw]_x + [bw]_x [aw]_x)")
        self.gen_add_code_line("//   = 0.5 * (aw bw^T + bw aw^T) - (aw . bw) I")
        self.gen_add_code_line("T adotb = aw_0*bw_0 + aw_1*bw_1 + aw_2*bw_2;")
        for r in range(3):
            for c in range(3):
                # Br[r, c] = 0.5 * (aw[r]*bw[c] + bw[r]*aw[c]) - adotb * (r == c)
                diag = " - adotb" if r == c else ""
                self.gen_add_code_line("T Br_" + str(r) + str(c) +
                                       " = static_cast<T>(0.5) * (aw_" + str(r) + "*bw_" + str(c) +
                                       " + bw_" + str(r) + "*aw_" + str(c) + ")" + diag + ";")
        # Bt_world = -Br_world @ pa
        for r in range(3):
            self.gen_add_code_line(
                "T Bt_" + str(r) + " = -(Br_" + str(r) + "0*pax + Br_" + str(r) + "1*pay + Br_" + str(r) + "2*paz);")
    elif (rev_a and not rev_b) or ((not rev_a) and rev_b):
        # Mixed lin-ang. Choose: a is rot if rev_a else b is rot.
        if rev_a:
            ang_var_x, ang_var_y, ang_var_z = "aw_0", "aw_1", "aw_2"
            lin_var_x, lin_var_y, lin_var_z = "blw_0", "blw_1", "blw_2"
        else:
            ang_var_x, ang_var_y, ang_var_z = "bw_0", "bw_1", "bw_2"
            lin_var_x, lin_var_y, lin_var_z = "alw_0", "alw_1", "alw_2"
        # Br_world = 0; Bt_world = 0.5 * (ang_world x lin_world)
        for r in range(3):
            for c in range(3):
                self.gen_add_code_line("T Br_" + str(r) + str(c) + " = static_cast<T>(0);")
        self.gen_add_code_line(
            "T Bt_0 = static_cast<T>(0.5) * (" + ang_var_y + "*" + lin_var_z + " - " + ang_var_z + "*" + lin_var_y + ");")
        self.gen_add_code_line(
            "T Bt_1 = static_cast<T>(0.5) * (" + ang_var_z + "*" + lin_var_x + " - " + ang_var_x + "*" + lin_var_z + ");")
        self.gen_add_code_line(
            "T Bt_2 = static_cast<T>(0.5) * (" + ang_var_x + "*" + lin_var_y + " - " + ang_var_y + "*" + lin_var_x + ");")
    else:
        # pris-pris: B_local = 0
        for r in range(3):
            for c in range(3):
                self.gen_add_code_line("T Br_" + str(r) + str(c) + " = static_cast<T>(0);")
        for r in range(3):
            self.gen_add_code_line("T Bt_" + str(r) + " = static_cast<T>(0);")
    # Now d2M = B_world * X_ee. We need d2M[:3, 3] (for H_xyz) and (d2R @ R_chain^T) = B_world[:3,:3].
    # d2M[:3, 3] = B_world[:3, :3] @ p_ee + B_world[:3, 3] (since X_ee[3, 3] = 1).
    self.gen_add_code_line("T Hxyz_x = Br_00*pex + Br_01*pey + Br_02*pez + Bt_0;")
    self.gen_add_code_line("T Hxyz_y = Br_10*pex + Br_11*pey + Br_12*pez + Bt_1;")
    self.gen_add_code_line("T Hxyz_z = Br_20*pex + Br_21*pey + Br_22*pez + Bt_2;")
    # [Jwi]_x and [Jwj]_x from S_i_world / S_j_world top-left
    self.gen_add_code_line("// Read [Jw_i]_x and [Jw_j]_x for the H_w correction")
    for c in range(3):
        for r in range(3):
            self.gen_add_code_line("T Si" + str(r) + str(c) + " = s_Sworld[" + str(si_base + r + 4*c) + "];")
    for c in range(3):
        for r in range(3):
            self.gen_add_code_line("T Sj" + str(r) + str(c) + " = s_Sworld[" + str(sj_base + r + 4*c) + "];")
    for r in range(3):
        for c in range(3):
            self.gen_add_code_line(
                "T SiSj" + str(r) + str(c) + " = Si" + str(r) + "0*Sj0" + str(c) +
                " + Si" + str(r) + "1*Sj1" + str(c) +
                " + Si" + str(r) + "2*Sj2" + str(c) + ";")
    # H_w = skew_inv(Br_world - SiSj)
    self.gen_add_code_line("T HW_x = static_cast<T>(0.5) * ((Br_21 - SiSj21) - (Br_12 - SiSj12));")
    self.gen_add_code_line("T HW_y = static_cast<T>(0.5) * ((Br_02 - SiSj02) - (Br_20 - SiSj20));")
    self.gen_add_code_line("T HW_z = static_cast<T>(0.5) * ((Br_10 - SiSj10) - (Br_01 - SiSj01));")
    base = "(" + str(ee_idx * 6 * nv * nv) + " + " + str(vi * nv + vj) + ")"
    self.gen_add_code_line("s_end_effector_pose_hessian[" + base + " + 0 * " + str(nv*nv) + "] = Hxyz_x;")
    self.gen_add_code_line("s_end_effector_pose_hessian[" + base + " + 1 * " + str(nv*nv) + "] = Hxyz_y;")
    self.gen_add_code_line("s_end_effector_pose_hessian[" + base + " + 2 * " + str(nv*nv) + "] = Hxyz_z;")
    self.gen_add_code_line("s_end_effector_pose_hessian[" + base + " + 3 * " + str(nv*nv) + "] = HW_x;")
    self.gen_add_code_line("s_end_effector_pose_hessian[" + base + " + 4 * " + str(nv*nv) + "] = HW_y;")
    self.gen_add_code_line("s_end_effector_pose_hessian[" + base + " + 5 * " + str(nv*nv) + "] = HW_z;")

def _emit_d2M_mimic_vslot_pair_block(self, ee_idx, ee_jid, vi, vj, nv, num_ees,
                                     blocks_i, blocks_j, si_base, sj_base):
    """Emit the d2(pose)/dv2 cell for a MIMIC-shared v-slot pair (vi, vj).

    When several chain joints fold into one velocity coordinate (a mimic joint
    and its target, or a multi-mimic finger like the h1_2 thumb: target + two
    mimics), the Hessian column for that v-slot is the SUM over EVERY block-pair
    (a in v-slot vi) x (b in v-slot vj), each scaled by alpha_a * alpha_b, with
    per-block-pair chain ordering -- exactly mirroring the RBDReference analytic
    double-block accumulate (the `vi_to_blocks` nested loop). The legacy
    per-chain-DOF emission overwrote the cell once per block-pair
    (last-writer-wins), which dropped every contribution but one (e.g. the h1_2
    thumb diagonal collapsed to a single mimic's term and the proximal/cross
    terms vanished -> exact-zero output columns 58798/101906/104948).

    blocks_i / blocks_j are lists of dicts:
        {chain_pos, jid, alpha, ang (3), lin (3), revolute (bool)}
    sorted by chain_pos. si_base / sj_base are the FOLDED s_Sworld bases used for
    the H_w (Jw_i x Jw_j) correction, which is built from the full folded angular
    columns (s_Sworld already holds the alpha-folded generator from Step 2).
    """
    self.gen_add_code_line("// MIMIC v-slot pair (vi=" + str(vi) + ", vj=" + str(vj) +
                           "): sum over " + str(len(blocks_i)) + "x" + str(len(blocks_j)) +
                           " block-pairs (alpha-weighted)")
    self.gen_add_code_line("T pex = s_Xworld[" + str(16*ee_jid + 12) + "];")
    self.gen_add_code_line("T pey = s_Xworld[" + str(16*ee_jid + 13) + "];")
    self.gen_add_code_line("T pez = s_Xworld[" + str(16*ee_jid + 14) + "];")
    # d2M accumulators: top 3 rows x 4 cols (column-major Mrc), zeroed.
    for r in range(3):
        for c in range(4):
            self.gen_add_code_line("T M" + str(r) + str(c) + " = static_cast<T>(0);")

    def _world_axis_lines(prefix, jid, ax_local):
        # axw_r = R_jid_world @ ax_local; R column-major in s_Xworld[16*jid].
        for r in range(3):
            terms = []
            for c in range(3):
                if abs(ax_local[c]) < 1e-15:
                    continue
                coef = "static_cast<T>(" + "{:.17g}".format(ax_local[c]) + ")"
                terms.append("s_Xworld[" + str(16*jid + r + 4*c) + "] * " + coef)
            expr = " + ".join(terms) if terms else "static_cast<T>(0)"
            self.gen_add_code_line("T " + prefix + "_" + str(r) + " = " + expr + ";")

    def _emit_block_generator(name, blk):
        # Per-block world generator G (top 3 rows, 4 cols) as scalars name_rc.
        # Revolute: G[:3,:3]=[axw]_x, G[:3,3]=pj x axw. Prismatic: G[:3,3]=axw.
        jid = blk["jid"]
        ax = blk["ang"] if blk["revolute"] else blk["lin"]
        ax = [float(ax[c]) if abs(ax[c]) >= 1e-15 else 0.0 for c in range(3)]
        _world_axis_lines(name + "w", jid, ax)
        self.gen_add_code_line("T " + name + "pjx = s_Xworld[" + str(16*jid + 12) + "];")
        self.gen_add_code_line("T " + name + "pjy = s_Xworld[" + str(16*jid + 13) + "];")
        self.gen_add_code_line("T " + name + "pjz = s_Xworld[" + str(16*jid + 14) + "];")
        if blk["revolute"]:
            self.gen_add_code_line("T " + name + "00 = static_cast<T>(0); T " + name + "11 = static_cast<T>(0); T " + name + "22 = static_cast<T>(0);")
            self.gen_add_code_line("T " + name + "10 =  " + name + "w_2; T " + name + "20 = -" + name + "w_1;")
            self.gen_add_code_line("T " + name + "01 = -" + name + "w_2; T " + name + "21 =  " + name + "w_0;")
            self.gen_add_code_line("T " + name + "02 =  " + name + "w_1; T " + name + "12 = -" + name + "w_0;")
            self.gen_add_code_line("T " + name + "03 = " + name + "pjy*" + name + "w_2 - " + name + "pjz*" + name + "w_1;")
            self.gen_add_code_line("T " + name + "13 = " + name + "pjz*" + name + "w_0 - " + name + "pjx*" + name + "w_2;")
            self.gen_add_code_line("T " + name + "23 = " + name + "pjx*" + name + "w_1 - " + name + "pjy*" + name + "w_0;")
        else:
            for r in range(3):
                for c in range(3):
                    self.gen_add_code_line("T " + name + str(r) + str(c) + " = static_cast<T>(0);")
            self.gen_add_code_line("T " + name + "03 = " + name + "w_0;")
            self.gen_add_code_line("T " + name + "13 = " + name + "w_1;")
            self.gen_add_code_line("T " + name + "23 = " + name + "w_2;")

    for blk_a in blocks_i:
        for blk_b in blocks_j:
            a = blk_a["chain_pos"]; b = blk_b["chain_pos"]
            scale = float(blk_a["alpha"]) * float(blk_b["alpha"])
            s_lit = "static_cast<T>(" + repr(scale) + ")"
            self.gen_add_code_line("{  // block-pair a_cp=" + str(a) + " b_cp=" + str(b) +
                                   " alpha_i*alpha_j=" + repr(scale))
            if a == b:
                # Same chain joint: scale * (L_a @ B_local @ Linv_a @ X_ee).
                jid = blk_a["jid"]
                ang_a = blk_a["ang"]; lin_a = blk_a["lin"]; rev_a = blk_a["revolute"]
                ang_b = blk_b["ang"]; lin_b = blk_b["lin"]; rev_b = blk_b["revolute"]
                _world_axis_lines("baw", jid, ang_a)
                _world_axis_lines("bbw", jid, ang_b)
                _world_axis_lines("balw", jid, lin_a)
                _world_axis_lines("bblw", jid, lin_b)
                if rev_a and rev_b:
                    self.gen_add_code_line("T badotb = baw_0*bbw_0 + baw_1*bbw_1 + baw_2*bbw_2;")
                    for r in range(3):
                        for c in range(3):
                            diag = " - badotb" if r == c else ""
                            self.gen_add_code_line("T Br" + str(r) + str(c) +
                                                   " = static_cast<T>(0.5)*(baw_" + str(r) + "*bbw_" + str(c) +
                                                   " + bbw_" + str(r) + "*baw_" + str(c) + ")" + diag + ";")
                    self.gen_add_code_line("T bpax = s_Xworld[" + str(16*jid + 12) + "]; T bpay = s_Xworld[" + str(16*jid + 13) + "]; T bpaz = s_Xworld[" + str(16*jid + 14) + "];")
                    for r in range(3):
                        self.gen_add_code_line("T Bt" + str(r) + " = -(Br" + str(r) + "0*bpax + Br" + str(r) + "1*bpay + Br" + str(r) + "2*bpaz);")
                elif (rev_a and not rev_b) or ((not rev_a) and rev_b):
                    if rev_a:
                        avx, avy, avz = "baw_0", "baw_1", "baw_2"
                        lvx, lvy, lvz = "bblw_0", "bblw_1", "bblw_2"
                    else:
                        avx, avy, avz = "bbw_0", "bbw_1", "bbw_2"
                        lvx, lvy, lvz = "balw_0", "balw_1", "balw_2"
                    for r in range(3):
                        for c in range(3):
                            self.gen_add_code_line("T Br" + str(r) + str(c) + " = static_cast<T>(0);")
                    self.gen_add_code_line("T Bt0 = static_cast<T>(0.5)*(" + avy + "*" + lvz + " - " + avz + "*" + lvy + ");")
                    self.gen_add_code_line("T Bt1 = static_cast<T>(0.5)*(" + avz + "*" + lvx + " - " + avx + "*" + lvz + ");")
                    self.gen_add_code_line("T Bt2 = static_cast<T>(0.5)*(" + avx + "*" + lvy + " - " + avy + "*" + lvx + ");")
                else:
                    for r in range(3):
                        for c in range(3):
                            self.gen_add_code_line("T Br" + str(r) + str(c) + " = static_cast<T>(0);")
                    for r in range(3):
                        self.gen_add_code_line("T Bt" + str(r) + " = static_cast<T>(0);")
                for r in range(3):
                    for c in range(3):
                        self.gen_add_code_line("M" + str(r) + str(c) + " += " + s_lit + " * Br" + str(r) + str(c) + ";")
                self.gen_add_code_line("M03 += " + s_lit + " * (Br00*pex + Br01*pey + Br02*pez + Bt0);")
                self.gen_add_code_line("M13 += " + s_lit + " * (Br10*pex + Br11*pey + Br12*pez + Bt1);")
                self.gen_add_code_line("M23 += " + s_lit + " * (Br20*pex + Br21*pey + Br22*pez + Bt2);")
            else:
                if a < b:
                    prox, dist = blk_a, blk_b
                else:
                    prox, dist = blk_b, blk_a
                _emit_block_generator("P", prox)
                _emit_block_generator("D", dist)
                for r in range(3):
                    for c in range(4):
                        self.gen_add_code_line("T Mb" + str(r) + str(c) + " = P" + str(r) + "0*D0" + str(c) +
                                               " + P" + str(r) + "1*D1" + str(c) + " + P" + str(r) + "2*D2" + str(c) + ";")
                for r in range(3):
                    for c in range(3):
                        self.gen_add_code_line("M" + str(r) + str(c) + " += " + s_lit + " * Mb" + str(r) + str(c) + ";")
                self.gen_add_code_line("M03 += " + s_lit + " * (Mb00*pex + Mb01*pey + Mb02*pez + Mb03);")
                self.gen_add_code_line("M13 += " + s_lit + " * (Mb10*pex + Mb11*pey + Mb12*pez + Mb13);")
                self.gen_add_code_line("M23 += " + s_lit + " * (Mb20*pex + Mb21*pey + Mb22*pez + Mb23);")
            self.gen_add_code_line("}")

    self.gen_add_code_line("T Hxyz_x = M03; T Hxyz_y = M13; T Hxyz_z = M23;")
    for c in range(3):
        for r in range(3):
            self.gen_add_code_line("T Si" + str(r) + str(c) + " = s_Sworld[" + str(si_base + r + 4*c) + "];")
    for c in range(3):
        for r in range(3):
            self.gen_add_code_line("T Sj" + str(r) + str(c) + " = s_Sworld[" + str(sj_base + r + 4*c) + "];")
    for r in range(3):
        for c in range(3):
            self.gen_add_code_line("T SiSj" + str(r) + str(c) + " = Si" + str(r) + "0*Sj0" + str(c) +
                                   " + Si" + str(r) + "1*Sj1" + str(c) + " + Si" + str(r) + "2*Sj2" + str(c) + ";")
    self.gen_add_code_line("T HW_x = static_cast<T>(0.5) * ((M21 - SiSj21) - (M12 - SiSj12));")
    self.gen_add_code_line("T HW_y = static_cast<T>(0.5) * ((M02 - SiSj02) - (M20 - SiSj20));")
    self.gen_add_code_line("T HW_z = static_cast<T>(0.5) * ((M10 - SiSj10) - (M01 - SiSj01));")
    base = "(" + str(ee_idx * 6 * nv * nv) + " + " + str(vi * nv + vj) + ")"
    self.gen_add_code_line("s_end_effector_pose_hessian[" + base + " + 0 * " + str(nv*nv) + "] = Hxyz_x;")
    self.gen_add_code_line("s_end_effector_pose_hessian[" + base + " + 1 * " + str(nv*nv) + "] = Hxyz_y;")
    self.gen_add_code_line("s_end_effector_pose_hessian[" + base + " + 2 * " + str(nv*nv) + "] = Hxyz_z;")
    self.gen_add_code_line("s_end_effector_pose_hessian[" + base + " + 3 * " + str(nv*nv) + "] = HW_x;")
    self.gen_add_code_line("s_end_effector_pose_hessian[" + base + " + 4 * " + str(nv*nv) + "] = HW_y;")
    self.gen_add_code_line("s_end_effector_pose_hessian[" + base + " + 5 * " + str(nv*nv) + "] = HW_z;")


def gen_end_effector_pose_hessian_device(self):
    n = self.robot.get_num_pos()
    nv = self.robot.get_num_vel()
    num_ees = self.robot.get_total_leaf_nodes()
    inner_temp_size = self.gen_end_effector_pose_hessian_inner_temp_mem_size()
    output_count = self.gen_end_effector_pose_hessian_output_count()
    # construct the boilerplate and function definition
    func_params = ["s_end_effector_pose_hessian is a pointer to shared memory of size 6*NUM_VEL*NUM_VEL*NUM_EE where NUM_VEL = " + str(nv) + " and NUM_EE = " + str(num_ees) + " (d^2/dv^2 tangent, pinocchio convention)", \
                   "s_end_effector_pose_gradient is a pointer to shared memory of size 6*NUM_VEL*NUM_EE (d/dv tangent Jacobian)", \
                   "s_q is the vector of joint positions", \
                   "d_robotModel is the pointer to the initialized model specific helpers on the GPU (XImats, topology_helpers, etc.)", \
                   "d_workspace is the global scratch buffer; size END_EFFECTOR_POSE_HESSIAN_DEVICE_INLINE_WORKSPACE_BYTES<T, RESOURCE_TIER>() bytes (= 0 at TIER_SHARED, " + str(output_count) + "*sizeof(T) at TIER_LITE+). Pass nullptr at TIER_SHARED"]
    func_notes = ["Inline-CUDA users: at TIER_LITE/TIER_MINIMAL the large s_end_effector_pose_hessian output (~" + str(output_count) + "*sizeof(T) bytes) moves from shared memory to d_workspace, freeing smem for the caller's outer kernel"]
    func_def_start = "void end_effector_pose_hessian_device("
    func_def_middle = "T *s_end_effector_pose_hessian, T *s_end_effector_pose_gradient, const T *s_q, "
    func_def_end = "const robotModel<T> *d_robotModel, T *d_workspace = nullptr) {"
    func_def = func_def_start + func_def_middle + func_def_end
    # then generate the code
    self.gen_add_func_doc("Computes the Hessian (and Jacobian) of the End Effector Pose with respect to generalized velocity (d^2/dv^2 tangent, pinocchio convention)",\
                          func_notes,func_params,None)
    self.gen_add_code_line("template <typename T, int RESOURCE_TIER = TIER_SHARED>")
    self.gen_add_code_line("__device__")
    self.gen_add_code_line(func_def, True)
    # Smem arena: s_temp is always the inner-temp size. The output s_end_effector_pose_hessian lives
    # in smem at TIER_SHARED (carved from the arena tail) and in d_workspace at
    # TIER_LITE/MINIMAL (inner repoints internally). Note: the geometric-Jacobian
    # path uses ONLY s_Xhom (LOCAL transforms); s_dXmatsHom and s_d2XmatsHom are
    # no longer needed (saves substantial smem on big robots).
    self.gen_XmatsHom_helpers_temp_shared_memory_code(inner_temp_size, include_gradients = False, include_hessians = False,
                                                      include_linalg_scratch = True,
                                                      linalg_scratch_bytes = "GRID_EE_LINALG_SHARED_BYTES<T>()")
    # At TIER_SHARED s_end_effector_pose_hessian is allocated by the caller; at LITE/MINIMAL it's
    # the inner's job to repoint via OUT_IN_SMEM=false + d_workspace.
    # then load Xhom (Jacobian only needs local transforms) and run the algo
    self.gen_load_update_XmatsHom_helpers_function_call(include_gradients = False, include_hessians = False)
    # Inner-owns placement: pass d_workspace + the per-tier flag. When the flag is
    # false the inner repoints s_end_effector_pose_hessian at d_workspace.
    self.gen_end_effector_pose_hessian_inner_function_call(
        updated_var_names = {"d_workspace_name": "d_workspace", "s_Xhom_name": "s_XmatsHom", "d_robotModel_name": "d_robotModel"},
        out_in_smem_expr = "D2EE_OUT_IN_SMEM<RESOURCE_TIER>()")
    self.gen_add_end_function()

_D2EE_PICK_FLAGS = [
    # use_workspace_output (s_end_effector_pose_hessian lives in d_workspace?)
    # Only one spill bit: the big nv^2 output. dXhom/d2Xhom are no longer
    # used by the FD-on-Jacobian inner, so the prior dXhom/d2Xhom spill bits
    # are gone.
    False,   # pick 0: full smem (PERF)
    True,    # pick 1: output in workspace (LITE)
    True,    # pick 2: same as pick 1 (MINIMAL -- no remaining smem to spill)
]

def _emit_d2ee_kernel_body_for_flags(self, n, num_ees, use_workspace_output,
                                     single_call_timing):
    """Emit the d2ee kernel body specialized for one tier's spill flags.
    Wrapped in a brace pair (caller emits the `if constexpr (...)` head).
    Used by gen_end_effector_pose_hessian_kernel to emit either a
    single body (collapsed picks) or three branched bodies (divergent picks)."""
    nv = self.robot.get_num_vel()
    output_count = self.gen_end_effector_pose_hessian_output_count()
    inner_temp_size = self.gen_end_effector_pose_hessian_inner_temp_mem_size()
    extra_t_buffers = [("s_q", n)] if use_workspace_output else [("s_q", n), ("s_end_effector_pose_hessian", output_count), ("s_end_effector_pose_gradient", 6*nv*num_ees)]
    self.gen_XmatsHom_helpers_temp_shared_memory_code(inner_temp_size, include_gradients = False, include_hessians = False,
                                                      extra_t_buffers = extra_t_buffers,
                                                      include_linalg_scratch = True,
                                                      linalg_scratch_bytes = "GRID_EE_LINALG_SHARED_BYTES<T>()")
    out_in_smem_expr = "false" if use_workspace_output else "true"
    if not single_call_timing:
        self.gen_add_parallel_loop("k","NUM_TIMESTEPS",block_level = True)
        self.gen_kernel_load_inputs("q",str(n),stride="stride_q")
        if use_workspace_output:
            self.gen_add_code_line("T *s_end_effector_pose_hessian = nullptr;  // inner repoints at d_workspace slice")
            self.gen_add_code_line("T *s_end_effector_pose_gradient = &d_end_effector_pose_gradient[k*" + str(6*nv*num_ees) + "];")
            self.gen_add_code_line("// Use d_end_effector_pose_hessian directly as the spill target so the inner writes into the persistent output buffer (one allocation, no extra copy).")
            self.gen_add_code_line("T *s_end_effector_pose_hessian_ws = &d_end_effector_pose_hessian[k*" + str(output_count) + "];")
        else:
            self.gen_add_code_line("(void)d_workspace;")
        self.gen_add_code_line("// compute")
        self.gen_load_update_XmatsHom_helpers_function_call(include_gradients = False, include_hessians = False)
        updated = {"s_Xhom_name": "s_XmatsHom", "d_robotModel_name": "d_robotModel"}
        if use_workspace_output:
            updated["d_workspace_name"] = "s_end_effector_pose_hessian_ws"
        self.gen_end_effector_pose_hessian_inner_function_call(
            updated_var_names = updated, out_in_smem_expr = out_in_smem_expr)
        self.gen_add_sync()
        if not use_workspace_output:
            self.gen_kernel_save_result("end_effector_pose_hessian",str(output_count),stride=str(output_count))
            self.gen_kernel_save_result("end_effector_pose_gradient",str(6*nv*num_ees),stride=str(6*nv*num_ees))
        else:
            # gradient still needs the smem -> global copy; the Hessian was already written to d_end_effector_pose_hessian directly via s_end_effector_pose_hessian_ws.
            self.gen_kernel_save_result("end_effector_pose_gradient",str(6*nv*num_ees),stride=str(6*nv*num_ees))
        self.gen_add_end_control_flow()
    else:
        self.gen_kernel_load_inputs("q",str(n))
        if use_workspace_output:
            self.gen_add_code_line("T *s_end_effector_pose_hessian = nullptr;  // inner repoints at d_workspace")
            self.gen_add_code_line("T *s_end_effector_pose_gradient = d_end_effector_pose_gradient;")
            self.gen_add_code_line("T *s_end_effector_pose_hessian_ws = d_end_effector_pose_hessian;  // use output buffer directly as spill target")
        else:
            self.gen_add_code_line("(void)d_workspace;")
        self.gen_add_code_line("// compute with NUM_TIMESTEPS as NUM_REPS for timing")
        self.gen_add_code_line("for (int rep = 0; rep < NUM_TIMESTEPS; rep++){", True)
        self.gen_anti_licm_input_reload("q",str(n),feedback_from="end_effector_pose_hessian")
        self.gen_load_update_XmatsHom_helpers_function_call(include_gradients = False, include_hessians = False)
        updated = {"s_Xhom_name": "s_XmatsHom", "d_robotModel_name": "d_robotModel"}
        if use_workspace_output:
            updated["d_workspace_name"] = "s_end_effector_pose_hessian_ws"
        self.gen_end_effector_pose_hessian_inner_function_call(
            updated_var_names = updated, out_in_smem_expr = out_in_smem_expr)
        self.gen_anti_licm_output_write("end_effector_pose_hessian")
        self.gen_add_end_control_flow()
        if not use_workspace_output:
            self.gen_kernel_save_result("end_effector_pose_hessian",str(output_count))
            self.gen_kernel_save_result("end_effector_pose_gradient",str(6*nv*num_ees))
        else:
            self.gen_kernel_save_result("end_effector_pose_gradient",str(6*nv*num_ees))


def gen_end_effector_pose_hessian_kernel(self, single_call_timing = False):
    n = self.robot.get_num_pos()
    num_ees = self.robot.get_total_leaf_nodes()
    func_params = ["d_end_effector_pose_hessian is the vector of end effector pose Hessians (6 x nv x nv per ee)", \
                   "d_end_effector_pose_gradient is the vector of end effector pose Jacobians (6 x nv per ee)", \
                   "d_workspace is the generated global spill workspace", \
                   "d_q is the vector of joint positions", \
                   "stride_q is the stide between each q", \
                   "d_robotModel is the pointer to the initialized model specific helpers on the GPU (XImats, topology_helpers, etc.)", \
                   "num_timesteps is the length of the trajectory points we need to compute over (or overloaded as test_iters for timing)"]
    func_notes = ["Output d^2(pose)/dv^2 is in tangent-space convention (d/dv), shape 6 x nv x nv per ee, C-order. Matches pinocchio."]
    func_def_start = "void end_effector_pose_hessian_kernel(T *d_end_effector_pose_hessian, T *d_end_effector_pose_gradient, unsigned char *d_workspace, const T *d_q, const int stride_q, "
    func_def_end = "const robotModel<T> *d_robotModel, const int NUM_TIMESTEPS) {"
    func_def = func_def_start + func_def_end
    if single_call_timing:
        func_def = func_def.replace("(", "_single_timing(")
    self.gen_add_func_doc("Computes the Hessian (and Jacobian) of the End Effector Pose with respect to generalized velocity (d^2/dv^2 tangent, pinocchio convention)",\
                          func_notes,func_params,None)
    self.gen_add_code_line("template <typename T, int RESOURCE_TIER = GRID_DEFAULT_RESOURCE_TIER>")
    self.gen_add_code_line("__global__")
    self.gen_add_code_line("__launch_bounds__(tier_max_threads<RESOURCE_TIER>())")
    self.gen_add_code_line(func_def, True)
    # Tier dispatch: when the 3 picks collapse, emit one body. When they
    # diverge, emit three if-constexpr branches -- each specialized for that
    # tier's spill flag. Smem-bytes constexpr END_EFFECTOR_POSE_HESSIAN_DYNAMIC_SHARED_MEM_BYTES<T,TIER>()
    # is already tier-aware.
    picks = getattr(self, "d2ee_spill_tier_3way", (0, 0, 0))
    def _emit_d2ee_body(pick):
        uwo = _D2EE_PICK_FLAGS[pick]
        _emit_d2ee_kernel_body_for_flags(self, n, num_ees, uwo, single_call_timing)
    self.gen_tier_dispatch(picks, _emit_d2ee_body)
    self.gen_add_end_function()

def gen_end_effector_pose_hessian_host(self, mode = 0):
    # default is to do the full kernel call -- options are for single timing or compute only kernel wrapper
    single_call_timing = True if mode == 1 else False
    compute_only = True if mode == 2 else False

    # define function def and params
    func_params = ["hd_data is the packaged input and output pointers", \
                   "d_robotModel is the pointer to the initialized model specific helpers on the GPU (XImats, topology_helpers, etc.)", \
                   "num_timesteps is the length of the trajectory points we need to compute over (or overloaded as test_iters for timing)", \
                   "streams are pointers to CUDA streams for async memory transfers (if needed)"]
    func_notes = []
    func_def_start = "void end_effector_pose_hessian(gridData<T, KIND> *hd_data, const robotModel<T> *d_robotModel, const int num_timesteps,"
    func_def_end =   "                            const dim3 block_dimms, const dim3 thread_dimms, cudaStream_t *streams) {"
    if single_call_timing:
        func_def_start = func_def_start.replace("(", "_single_timing(")
        func_def_end = "              " + func_def_end
    if compute_only:
        func_def_start = func_def_start.replace("(", "_compute_only(")
        func_def_end = "             " + func_def_end.replace(", cudaStream_t *streams", "")
    # then generate the code
    self.gen_add_func_doc("Computes the Hessian (and Jacobian) of the End Effector Pose with respect to generalized velocity (d^2/dv^2 tangent, pinocchio convention)",\
                          func_notes,func_params,None)
    self.gen_add_code_line("template <typename T, bool USE_COMPRESSED_MEM = false, gridDataKind KIND = GRID_DATA_ALL>")
    self.gen_add_code_line("__host__")
    self.gen_add_code_line(func_def_start)
    self.gen_add_code_line(func_def_end, True)
    self.gen_add_code_line("static_assert(KIND == GRID_DATA_ALL || KIND == GRID_DATA_KINEMATICS, \"end_effector_pose_hessian requires all-data or kinematics gridData\");")
    func_call_start = "end_effector_pose_hessian_kernel<T><<<block_dimms,thread_dimms,END_EFFECTOR_POSE_HESSIAN_DYNAMIC_SHARED_MEM_BYTES<T>()>>>(hd_data->d_end_effector_pose_hessian,hd_data->d_end_effector_pose_gradient,hd_data->d_workspace,hd_data->d_q,stride_q,"
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
    self.gen_add_code_line("if (END_EFFECTOR_POSE_HESSIAN_DYNAMIC_SHARED_MEM_BYTES<T>() > GRID_CUDA_TARGET_SHARED_MEM_BYTES) {fprintf(stderr,\"GRID end_effector_pose_hessian shared-memory request %zu exceeds compile target %d; regenerate with a deeper Hessian spill fallback or a higher GRID_CUDA_TARGET_SHARED_MEM_BYTES.\\n\", END_EFFECTOR_POSE_HESSIAN_DYNAMIC_SHARED_MEM_BYTES<T>(), GRID_CUDA_TARGET_SHARED_MEM_BYTES); gpuErrchk(cudaErrorInvalidConfiguration);}")
    self.gen_add_code_line("gpuErrchk(grid_check_dynamic_shared_memory_bytes(\"end_effector_pose_hessian\", END_EFFECTOR_POSE_HESSIAN_DYNAMIC_SHARED_MEM_BYTES<T>()));")
    # No L2 persistence: at LITE/MINIMAL the end_effector_pose_hessian spill target IS the output
    # buffer (d_end_effector_pose_hessian), which is written once and read once -- no benefit from
    # L2 pinning.
    self.gen_add_code_lines(func_call_code)
    if not compute_only:
        # then transfer memory back
        self.gen_add_code_lines(["// finally transfer the result back", \
                                 "gpuErrchk(cudaMemcpy(hd_data->h_end_effector_pose_gradient,hd_data->d_end_effector_pose_gradient,6*NUM_EES*NUM_VEL*" + \
                                    ("num_timesteps*" if not single_call_timing else "") + "sizeof(T),cudaMemcpyDeviceToHost));",
                                 "gpuErrchk(cudaMemcpy(hd_data->h_end_effector_pose_hessian,hd_data->d_end_effector_pose_hessian,6*NUM_EES*NUM_VEL*NUM_VEL*" + \
                                    ("num_timesteps*" if not single_call_timing else "") + "sizeof(T),cudaMemcpyDeviceToHost));",
                                 "gpuErrchkKernel();"])
    # finally report out timing if requested
    if single_call_timing:
        from ..algo_registry import single_call_printf_line
        self.gen_add_code_line(single_call_printf_line("end_effector_pose_hessian"))
    self.gen_add_end_function()

def gen_ee_pose_inner_xform_from_q_lines(self, lane_guarded = False):
    # Robot-general per-joint homogeneous-transform refresh from s_q.
    #
    # Emits the q-DEPENDENT cells of every joint's 4x4 s_XmatsHom block,
    # computing sin/cos of the joint angle inline. The constant cells (fixed
    # rotation pattern + translation offsets) are assumed to already be present
    # in s_XmatsHom -- the caller pre-loads them once (e.g. from
    # d_robotModel->d_XImats) and only the trig cells change with q. This
    # mirrors the serial section of gen_load_update_XmatsHom_helpers but as a
    # standalone, d_robotModel-free device snippet.
    #
    # Generality: driven by the symbolic per-joint Xmats (any DoF / axis /
    # joint type the parser emits), NOT a hardcoded iiwa14 7R pattern. Returns
    # the list of "joint j touches its block" booleans so warp lane-partitioning
    # can mirror it.
    import sympy as sp
    NJ = self.robot.get_num_joints()
    Xmats_hom = self.robot.get_Xmats_hom_ordered_by_id(include_fixed_joints = False)
    fb = self.robot.floating_base
    has_mimic = self.robot_has_mimic_joints()
    if fb or has_mimic:
        # The standalone (d_robotModel-free) inner only supports the plain
        # fixed-base, non-mimic angle->sincos substitution. Floating-base roots
        # and mimic q-folding need the s_temp/s_q_eff scratch that the general
        # load_update path provides; route those through end_effector_pose_inner
        # instead. Flag loudly so a bad emit fails at codegen, not silently.
        raise NotImplementedError(
            "ee_pose_inner_{thread,warp}: floating-base / mimic robots are not "
            "supported by the standalone FK inner; use end_effector_pose for those.")
    # which joints actually have q-dependent cells (so warp can skip the rest)
    joint_has_q = []
    for jid in range(NJ):
        M = Xmats_hom[jid]
        joint_has_q.append(any(not self.custom_is_constant(M[r, c])
                               for r in range(4) for c in range(4)))
    return Xmats_hom, joint_has_q

def gen_ee_pose_inner_thread(self, fixed_target_name = ""):
    import sympy as sp
    NJ = self.robot.get_num_joints()
    parents = [self.robot.get_parent_id(jid) for jid in range(NJ)]
    Xmats_hom, joint_has_q = self.gen_ee_pose_inner_xform_from_q_lines()

    self.gen_add_func_doc(
        "Thread-per-sample forward kinematics: serial chain walk that fills the "
        "full cumulative (world-frame) joint transforms s_jointXforms from s_q. "
        "Robot-general (any DoF / joint types). Reads s_jointXforms[16*target_idx] "
        "for the world pose of frame target_idx.",
        ["Assumes the constant (q-independent) cells of s_XmatsHom are pre-loaded; "
         "only the q-dependent cells are refreshed here."],
        [
            "s_jointXforms is the pointer to the cumulative (world) joint transforms (16 per joint)",
            "s_XmatsHom is the pointer to the per-joint homogeneous transforms (16 per joint)",
            "s_q is the vector of joint positions",
            "target_idx is the joint index whose world transform is desired (full array is filled)",
        ],
        None
    )
    self.gen_add_code_line("template <typename T>")
    self.gen_add_code_line("__device__")
    self.gen_add_code_line("void ee_pose_inner_thread(T *s_jointXforms, T *s_XmatsHom, T *s_q, int target_idx) {", True)
    self.gen_add_code_line("(void)target_idx;")

    # --- refresh q-dependent cells of s_XmatsHom from s_q (general) ---
    for jid in range(NJ):
        if not joint_has_q[jid]:
            continue
        qslot = self.robot.get_joint_index_q(jid)
        if isinstance(qslot, (list, tuple)):
            qslot = qslot[0]
        self.gen_add_code_line("// X_hom[" + str(jid) + "] q-dependent cells")
        self.gen_add_code_line("{", True)
        self.gen_add_code_line("const T s = static_cast<T>(sin(s_q[" + str(qslot) + "]));")
        self.gen_add_code_line("const T c = static_cast<T>(cos(s_q[" + str(qslot) + "]));")
        self.gen_add_code_line("(void)s; (void)c;")
        M = Xmats_hom[jid]
        for col in range(4):
            for row in range(4):
                val = M[row, col]
                if self.custom_is_constant(val):
                    continue
                str_val = sp.ccode(val)
                str_val = str_val.replace("sin(theta)", "s").replace("cos(theta)", "c")
                cell = self.gen_static_array_ind_3d(jid, col, row, ind_stride=16, col_stride=4)
                self.gen_add_code_line("s_XmatsHom[16*" + str(jid) + " + " + str(cell - 16*jid) +
                                       "] = static_cast<T>(" + str_val + ");")
        self.gen_add_end_control_flow()

    # --- chain walk: world transform of every joint by id order (parent < child) ---
    self.gen_add_code_line("// chain up cumulative world transforms (parent always < child)")
    self.gen_add_code_line("#pragma unroll")
    self.gen_add_code_line("for (int j = 0; j < " + str(NJ) + "; ++j) {", True)
    self.gen_add_code_line("const T* c = &s_XmatsHom[j * 16];")
    self.gen_add_code_line("T* o = &s_jointXforms[j * 16];")
    self.gen_add_code_line("int par = " + self.gen_ee_pose_inner_parent_lookup(parents) + ";")
    self.gen_add_code_line("if (par < 0) {", True)
    self.gen_add_code_line("o[0]=c[0];   o[1]=c[1];   o[2]=c[2];")
    self.gen_add_code_line("o[4]=c[4];   o[5]=c[5];   o[6]=c[6];")
    self.gen_add_code_line("o[8]=c[8];   o[9]=c[9];   o[10]=c[10];")
    self.gen_add_code_line("o[12]=c[12]; o[13]=c[13]; o[14]=c[14];")
    self.gen_add_end_control_flow()
    self.gen_add_code_line("else {", True)
    self.gen_add_code_line("const T* p = &s_jointXforms[par * 16];")
    self.gen_add_code_line("o[0]  = p[0]*c[0]   + p[4]*c[1]   + p[8]*c[2];")
    self.gen_add_code_line("o[1]  = p[1]*c[0]   + p[5]*c[1]   + p[9]*c[2];")
    self.gen_add_code_line("o[2]  = p[2]*c[0]   + p[6]*c[1]   + p[10]*c[2];")
    self.gen_add_code_line("o[4]  = p[0]*c[4]   + p[4]*c[5]   + p[8]*c[6];")
    self.gen_add_code_line("o[5]  = p[1]*c[4]   + p[5]*c[5]   + p[9]*c[6];")
    self.gen_add_code_line("o[6]  = p[2]*c[4]   + p[6]*c[5]   + p[10]*c[6];")
    self.gen_add_code_line("o[8]  = p[0]*c[8]   + p[4]*c[9]   + p[8]*c[10];")
    self.gen_add_code_line("o[9]  = p[1]*c[8]   + p[5]*c[9]   + p[9]*c[10];")
    self.gen_add_code_line("o[10] = p[2]*c[8]   + p[6]*c[9]   + p[10]*c[10];")
    self.gen_add_code_line("o[12] = p[0]*c[12]  + p[4]*c[13]  + p[8]*c[14]  + p[12];")
    self.gen_add_code_line("o[13] = p[1]*c[12]  + p[5]*c[13]  + p[9]*c[14]  + p[13];")
    self.gen_add_code_line("o[14] = p[2]*c[12]  + p[6]*c[13]  + p[10]*c[14] + p[14];")
    self.gen_add_end_control_flow()
    self.gen_add_code_line("o[3]=(T)0; o[7]=(T)0; o[11]=(T)0; o[15]=(T)1;")
    self.gen_add_end_control_flow()

    self.gen_add_end_function()

def gen_ee_pose_inner_parent_lookup(self, parents):
    # Emit a constant lookup expression mapping the loop index j -> parent id.
    # For a serial chain this is just (j-1); otherwise emit a small static table.
    if all(parents[j] == j - 1 for j in range(len(parents))):
        return "j - 1"
    arr = "{" + ",".join(str(p) for p in parents) + "}"
    return "((const int[]) " + arr + ")[j]"

def gen_ee_pose_inner_warp(self, fixed_target_name = ""):
    import sympy as sp
    NJ = self.robot.get_num_joints()
    parents = [self.robot.get_parent_id(jid) for jid in range(NJ)]
    Xmats_hom, joint_has_q = self.gen_ee_pose_inner_xform_from_q_lines()

    self.gen_add_func_doc(
        "Warp-per-sample forward kinematics: warp-cooperative chain walk that fills "
        "the full cumulative (world-frame) joint transforms s_jointXforms from s_q. "
        "Robot-general (any DoF / joint types). 3 lanes own the matrix rows; "
        "__syncwarp between levels. Reads s_jointXforms[16*target_idx] for the "
        "world pose of frame target_idx.",
        ["Assumes the constant (q-independent) cells of s_XmatsHom are pre-loaded; "
         "only the q-dependent cells are refreshed here."],
        [
            "s_jointXforms is the pointer to the cumulative (world) joint transforms (16 per joint)",
            "s_XmatsHom is the pointer to the per-joint homogeneous transforms (16 per joint)",
            "s_q is the vector of joint positions",
            "target_idx is the joint index whose world transform is desired (full array is filled)",
        ],
        None
    )
    self.gen_add_code_line("template <typename T>")
    self.gen_add_code_line("__device__ inline void ee_pose_inner_warp(")
    self.gen_add_code_line("    T* __restrict__ s_jointXforms,")
    self.gen_add_code_line("    T* __restrict__ s_XmatsHom,")
    self.gen_add_code_line("    const T* __restrict__ s_q,")
    self.gen_add_code_line("    int target_idx)")
    self.gen_add_code_line("{", True)
    self.gen_add_code_line("(void)target_idx;")
    self.gen_add_code_line("const int lane = threadIdx.x & 31;")
    self.gen_add_code_line("const unsigned mask = 0xFFFFFFFFu;")

    # --- refresh q-dependent cells: one lane per joint (lane == jid) ---
    # joints with q-dependent cells; <=31 lanes cover them, else fall back to a
    # lane-strided loop so any NJ works.
    q_joints = [jid for jid in range(NJ) if joint_has_q[jid]]
    self.gen_add_code_line("// refresh q-dependent X_hom cells: lane j owns joint j")
    for jid in q_joints:
        qslot = self.robot.get_joint_index_q(jid)
        if isinstance(qslot, (list, tuple)):
            qslot = qslot[0]
        self.gen_add_code_line("if (lane == " + str(jid) + ") {", True)
        self.gen_add_code_line("const T s = static_cast<T>(sin(s_q[" + str(qslot) + "]));")
        self.gen_add_code_line("const T c = static_cast<T>(cos(s_q[" + str(qslot) + "]));")
        self.gen_add_code_line("(void)s; (void)c;")
        M = Xmats_hom[jid]
        for col in range(4):
            for row in range(4):
                val = M[row, col]
                if self.custom_is_constant(val):
                    continue
                str_val = sp.ccode(val)
                str_val = str_val.replace("sin(theta)", "s").replace("cos(theta)", "c")
                cell = self.gen_static_array_ind_3d(jid, col, row, ind_stride=16, col_stride=4)
                self.gen_add_code_line("s_XmatsHom[16*" + str(jid) + " + " + str(cell - 16*jid) +
                                       "] = static_cast<T>(" + str_val + ");")
        self.gen_add_end_control_flow()
    self.gen_add_code_line("__syncwarp(mask);")

    # --- chain walk: 3 row-lanes per joint, parent always < child ---
    self.gen_add_code_line("// chain up cumulative world transforms (parent always < child)")
    self.gen_add_code_line("#pragma unroll")
    self.gen_add_code_line("for (int j = 0; j < " + str(NJ) + "; ++j) {", True)
    self.gen_add_code_line("const T* c = &s_XmatsHom[j * 16];")
    self.gen_add_code_line("T* o = &s_jointXforms[j * 16];")
    self.gen_add_code_line("int par = " + self.gen_ee_pose_inner_parent_lookup(parents) + ";")

    self.gen_add_code_line("if (par < 0) {", True)
    self.gen_add_code_line("if (lane == 0) { o[0]=c[0]; o[4]=c[4]; o[8]=c[8];  o[12]=c[12]; o[3]=(T)0; o[7]=(T)0; o[11]=(T)0; o[15]=(T)1; }")
    self.gen_add_code_line("else if (lane == 1) { o[1]=c[1]; o[5]=c[5]; o[9]=c[9];  o[13]=c[13]; }")
    self.gen_add_code_line("else if (lane == 2) { o[2]=c[2]; o[6]=c[6]; o[10]=c[10]; o[14]=c[14]; }")
    self.gen_add_end_control_flow()
    self.gen_add_code_line("else {", True)
    self.gen_add_code_line("const T* p = &s_jointXforms[par * 16];")
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
    self.gen_add_end_control_flow()
    self.gen_add_code_line("__syncwarp(mask);")
    self.gen_add_end_control_flow()

    self.gen_add_end_function()

def gen_ee_pose_fk_batched_kernel(self):
    # Large-batch FK kernel: ONE BLOCK PER SAMPLE (b = blockIdx.x), mirroring
    # the HJCD-IK launch <<<B, threads>>>. Each block walks its sample's whole
    # chain via the ee_pose_inner_{thread,warp} device inner and writes a
    # 7-element pose (position + quaternion). USE_WARP selects the warp- vs
    # thread-cooperative inner; both produce identical poses.
    #   d_q layout:    q[b*stride_q + j]      (batch-major, stride_q == NUM_POS)
    #   d_pose7 layout: pose7[b*7 + 0..2] = translation, [3..6] = quaternion (w,x,y,z)
    n = self.robot.get_num_pos()
    NJ = self.robot.get_num_joints()
    Xhom_size, _, _ = self.gen_get_Xhom_size()
    temp_size = self.gen_load_update_XImats_helpers_temp_mem_size()
    default_ee = self.robot.get_leaf_nodes()[0]
    self.gen_add_func_doc(
        "Batched forward kinematics: one block per sample, pos+quat output.",
        ["USE_WARP picks the warp-cooperative inner (warp 0) vs the thread inner (thread 0).",
         "target_idx selects the output frame (defaults to the leaf EE joint id)."],
        ["d_pose7 is the (B x 7) output: [tx,ty,tz, qw,qx,qy,qz] per sample",
         "d_q is the (B x NUM_POS) joint-position input (batch-major, stride stride_q)",
         "stride_q is the stride between samples in d_q (== NUM_POS)",
         "d_robotModel holds the per-robot constants (XImats, topology helpers)",
         "B is the batch size (== gridDim.x)",
         "target_idx is the joint frame whose world pose is written"],
        None)
    self.gen_add_code_line("template <typename T, bool USE_WARP = false>")
    self.gen_add_code_line("__global__")
    self.gen_add_code_line("void ee_pose_fk_batched_kernel(T *d_pose7, const T *d_q, const int stride_q, "
                           "const robotModel<T> *d_robotModel, const int B, const int target_idx = " + str(default_ee) + ") {", True)
    # static shared per-block scratch (sizes are compile-time constants)
    self.gen_add_code_line("__shared__ T s_q[" + str(n) + "];")
    self.gen_add_code_line("__shared__ T s_XmatsHom[" + str(Xhom_size) + "];")
    self.gen_add_code_line("__shared__ T s_jointXforms[" + str(16*NJ) + "];")
    self.gen_add_code_line("__shared__ T s_temp[" + str(max(temp_size,1)) + "];")
    if not self.robot.is_serial_chain() or not self.robot.are_Ss_identical(list(range(n))):
        self.gen_add_code_line("__shared__ int s_topology_helpers[" + str(self.gen_topology_helpers_size()) + "];")
    self.gen_add_code_line("for (int b = blockIdx.x; b < B; b += gridDim.x) {", True)
    # cooperative load of this sample's q
    self.gen_add_code_line("for (int j = threadIdx.x; j < " + str(n) + "; j += blockDim.x) { s_q[j] = d_q[b*stride_q + j]; }")
    self.gen_add_sync()
    # fill s_XmatsHom (constant + q-dependent cells) via the canonical path
    self.gen_load_update_XmatsHom_helpers_function_call()
    self.gen_add_sync()
    # walk the chain with the requested cooperative inner
    self.gen_add_code_line("if (USE_WARP) {", True)
    self.gen_add_code_line("if ((threadIdx.x >> 5) == 0) { ee_pose_inner_warp<T>(s_jointXforms, s_XmatsHom, s_q, target_idx); }")
    self.gen_add_end_control_flow()
    self.gen_add_code_line("else {", True)
    self.gen_add_code_line("if (threadIdx.x == 0) { ee_pose_inner_thread<T>(s_jointXforms, s_XmatsHom, s_q, target_idx); }")
    self.gen_add_end_control_flow()
    self.gen_add_sync()
    # extract pos + quaternion from the target's 4x4 (column-major homogeneous)
    self.gen_add_code_line("if (threadIdx.x == 0) {", True)
    self.gen_add_code_line("const T* X = &s_jointXforms[16*target_idx];")
    self.gen_add_code_line("// rotation block (column-major): R[r][c] = X[4*c + r]")
    self.gen_add_code_line("const T r00=X[0], r10=X[1], r20=X[2];")
    self.gen_add_code_line("const T r01=X[4], r11=X[5], r21=X[6];")
    self.gen_add_code_line("const T r02=X[8], r12=X[9], r22=X[10];")
    self.gen_add_code_line("T* o = &d_pose7[b*7];")
    self.gen_add_code_line("o[0]=X[12]; o[1]=X[13]; o[2]=X[14];")
    self.gen_add_code_line("// quaternion (w,x,y,z) from rotation (Shepperd's method)")
    self.gen_add_code_line("const T tr = r00 + r11 + r22;")
    self.gen_add_code_line("T qw,qx,qy,qz;")
    self.gen_add_code_line("if (tr > (T)0) {", True)
    self.gen_add_code_line("T S = sqrt(tr + (T)1) * (T)2; qw = (T)0.25*S; qx=(r21-r12)/S; qy=(r02-r20)/S; qz=(r10-r01)/S;")
    self.gen_add_end_control_flow()
    self.gen_add_code_line("else if (r00 > r11 && r00 > r22) {", True)
    self.gen_add_code_line("T S = sqrt((T)1 + r00 - r11 - r22) * (T)2; qw=(r21-r12)/S; qx=(T)0.25*S; qy=(r01+r10)/S; qz=(r02+r20)/S;")
    self.gen_add_end_control_flow()
    self.gen_add_code_line("else if (r11 > r22) {", True)
    self.gen_add_code_line("T S = sqrt((T)1 + r11 - r00 - r22) * (T)2; qw=(r02-r20)/S; qx=(r01+r10)/S; qy=(T)0.25*S; qz=(r12+r21)/S;")
    self.gen_add_end_control_flow()
    self.gen_add_code_line("else {", True)
    self.gen_add_code_line("T S = sqrt((T)1 + r22 - r00 - r11) * (T)2; qw=(r10-r01)/S; qx=(r02+r20)/S; qy=(r12+r21)/S; qz=(T)0.25*S;")
    self.gen_add_end_control_flow()
    self.gen_add_code_line("o[3]=qw; o[4]=qx; o[5]=qy; o[6]=qz;")
    self.gen_add_end_control_flow()
    self.gen_add_sync()
    self.gen_add_end_control_flow()
    self.gen_add_end_function()

def gen_ee_pose_fk_batched_host(self):
    # Host launcher for the batched FK kernel. Mirrors the existing batched
    # end_effector_pose host convention (device buffers + streams) but takes
    # the (B x NUM_POS) input -> (B x 7) pos+quat output directly.
    n = self.robot.get_num_pos()
    default_ee = self.robot.get_leaf_nodes()[0]
    self.gen_add_func_doc(
        "Host launcher for batched FK (<<<B, threads>>>, one block per sample).",
        ["USE_WARP selects the warp- vs thread-cooperative per-sample inner."],
        ["d_pose7 is the device (B x 7) output buffer",
         "d_q is the device (B x NUM_POS) input buffer (batch-major)",
         "B is the batch size",
         "d_robotModel holds the per-robot constants",
         "threads is the per-block thread count (>=32 for the warp variant)",
         "target_idx selects the output frame (defaults to the leaf EE)",
         "stream is the CUDA stream"],
        None)
    self.gen_add_code_line("template <typename T, bool USE_WARP = false>")
    self.gen_add_code_line("__host__")
    self.gen_add_code_line("void ee_pose_fk_batched(T *d_pose7, const T *d_q, const int B, "
                           "const robotModel<T> *d_robotModel, const int threads = 32, "
                           "const int target_idx = " + str(default_ee) + ", cudaStream_t stream = (cudaStream_t)0) {", True)
    self.gen_add_code_line("const int stride_q = " + str(n) + ";")
    self.gen_add_code_line("ee_pose_fk_batched_kernel<T, USE_WARP><<<B, threads, 0, stream>>>("
                           "d_pose7, d_q, stride_q, d_robotModel, B, target_idx);")
    self.gen_add_end_function()
    self.gen_add_code_line("#define GRID_HAS_FK_BATCHED 1")

def gen_eepose_and_derivatives(self, fixed_target_name = "",
                               include_pose = True, include_gradient = True, include_hessian = True):
    ee_target_names = [""]
    if fixed_target_name == "all":
        ee_target_names += [fj.name for fj in self.robot.fixed_joints]
    elif fixed_target_name != "":
        ee_target_names += [fixed_target_name]
    for target in ee_target_names:
        if include_pose:
            # first generate the inner helpers
            self.gen_end_effector_pose_inner(fixed_target_name = target)
            # then generate the device wrappers
            self.gen_end_effector_pose_device(fixed_target_name = target)
            # then generate the kernels
            self.gen_end_effector_pose_kernel(single_call_timing = True, fixed_target_name = target)
            self.gen_end_effector_pose_kernel(single_call_timing = False, fixed_target_name = target)
            # then the host launch wrappers
            self.gen_end_effector_pose_host(0, fixed_target_name = target)
            self.gen_end_effector_pose_host(1, fixed_target_name = target)
            self.gen_end_effector_pose_host(2, fixed_target_name = target)

        if include_gradient:
            # then for the gradient first generate the inner helpers
            self.gen_end_effector_pose_gradient_inner(fixed_target_name = target)
            # then generate the device wrappers
            self.gen_end_effector_pose_gradient_device(fixed_target_name = target)
            # then generate the kernels
            self.gen_end_effector_pose_gradient_kernel(True, fixed_target_name = target)
            self.gen_end_effector_pose_gradient_kernel(False, fixed_target_name = target)
            # then the host launch wrappers
            self.gen_end_effector_pose_gradient_host(0, fixed_target_name = target)
            self.gen_end_effector_pose_gradient_host(1, fixed_target_name = target)
            self.gen_end_effector_pose_gradient_host(2, fixed_target_name = target)

    if include_hessian:
        # then for the hessian first generate the inner helpers
        self.gen_end_effector_pose_hessian_inner()
        # then generate the device wrappers
        self.gen_end_effector_pose_hessian_device()
        # then generate the kernels
        self.gen_end_effector_pose_hessian_kernel(True)
        self.gen_end_effector_pose_hessian_kernel(False)
        # then the host launch wrappers
        self.gen_end_effector_pose_hessian_host(0)
        self.gen_end_effector_pose_hessian_host(1)
        self.gen_end_effector_pose_hessian_host(2)

    if include_pose or include_gradient or include_hessian:
        # standalone warp/thread FK inners + batched convenience path.
        # Skip ENTIRELY for floating-base / mimic robots: the standalone inner
        # does not support those (it raises) — they route through
        # end_effector_pose instead.
        if not self.robot.floating_base and not self.robot_has_mimic_joints():
            self.gen_ee_pose_inner_thread(fixed_target_name = fixed_target_name)
            self.gen_ee_pose_inner_warp(fixed_target_name = fixed_target_name)
            self.gen_ee_pose_fk_batched_kernel()
            self.gen_ee_pose_fk_batched_host()
