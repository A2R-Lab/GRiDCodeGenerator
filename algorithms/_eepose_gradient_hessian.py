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
    self.gen_add_code_line("// Now extract the eePos from the Tansforms")
    self.gen_add_code_line("// TODO: ADD OFFSETS")
    self.gen_add_code_line("//")
    tempOffset = 16*num_ees*(bfs_level % 2)
    # xyz position is easy (eePos_xyz1 = Xmat_hom * offset) where offset = [x,y,z,1]
    self.gen_add_parallel_loop("ind",str(3*num_ees))
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
    func_params = ["s_eePos is a pointer to shared memory of size 6*NUM_EE where NUM_EE = " + str(num_ees), \
                   "s_q is the vector of joint positions", \
                   "d_robotModel is the pointer to the initialized model specific helpers on the GPU (XImats, topology_helpers, etc.)"]
    func_notes = []
    func_def_start = "void end_effector_pose_device" + ("" if fixed_target_name == "" else "_" + fixed_target_name) + "("
    func_def_middle = "T *s_eePos, const T *s_q, "
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
        self.gen_kernel_save_result("eePos",str(6*num_ees),stride=str(6*num_ees))
        self.gen_add_end_control_flow()
    else:
        #repurpose NUM_TIMESTEPS for number of timing reps
        self.gen_kernel_load_inputs("q",str(n))
        # then compute in loop for timing
        self.gen_add_code_line("// compute with NUM_TIMESTEPS as NUM_REPS for timing")
        self.gen_add_code_line("for (int rep = 0; rep < NUM_TIMESTEPS; rep++){", True)
        self.gen_anti_licm_input_reload("q",str(n),feedback_from="eePos")
        # then load/update X and run the algo
        self.gen_load_update_XmatsHom_helpers_function_call()
        self.gen_end_effector_pose_inner_function_call(fixed_target_name = fixed_target_name)
        self.gen_anti_licm_output_write("eePos")
        self.gen_add_end_control_flow()
        # save to global
        self.gen_kernel_save_result("eePos",str(6*num_ees))
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

    Output `s_deePos` is sized 6 * nv * NUM_EE (NOT 6 * nq) so the floating-base
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
      4. Write s_deePos: rows 0..2 = J_v columns; rows 3..5 = E(rpy)^{-1} * J_w
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
    func_params = ["s_deePos is a pointer to shared memory of size 6*NUM_VEL*NUM_EE where NUM_VEL = " + str(nv) + " and NUM_EE = " + str(num_ees), \
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
    func_def_middle = "T *s_deePos, const T *s_q, const T *s_Xhom, const T *s_dXhom, "
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
    # these occur when several joint ids in a chain map to one velocity coordinate
    # (e.g. h1_2's finger joints sharing a v-index). The original serial emission
    # wrote them in order and let the LAST writer win; keeping only the last
    # occurrence per (ee_idx, vi) reproduces that last-writer-wins result exactly
    # while making every remaining job's destination column disjoint (a hard
    # requirement for the block-parallel fill — no two work-items touch one cell).
    _dedup = {}
    for entry in flat_jobs:
        ee_idx, _anc, job = entry
        _dedup[(ee_idx, job["vi"])] = entry  # later entries overwrite earlier
    flat_jobs = list(_dedup.values())
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

    # ============ Step 5: write s_deePos = [J_v ; E^{-1} J_w] ============
    self.gen_add_code_line("//")
    self.gen_add_code_line("// Step 5: write s_deePos (rows 0..2 = J_v, rows 3..5 = E(rpy)^{-1} J_w)")
    self.gen_add_code_line("//")
    self.gen_add_parallel_loop("ind", str(6 * nv * num_ees))
    self.gen_add_code_line("int row = ind % 6; int rem = ind / 6; int vi = rem % " + str(nv) + "; int ee = rem / " + str(nv) + ";")
    self.gen_add_code_line("int jv_base = 3 * (" + str(nv) + " * ee + vi);")
    self.gen_add_code_line("if (row < 3) {", True)
    self.gen_add_code_line("s_deePos[ind] = s_Jv[jv_base + row];")
    self.gen_add_end_control_flow()
    self.gen_add_code_line("else {", True)
    self.gen_add_code_line("T cy = s_E_sc[4*ee + 0]; T sy = s_E_sc[4*ee + 1]; T cp = s_E_sc[4*ee + 2]; T sp = s_E_sc[4*ee + 3];")
    self.gen_add_code_line("T Jw0 = s_Jw[jv_base + 0]; T Jw1 = s_Jw[jv_base + 1]; T Jw2 = s_Jw[jv_base + 2];")
    self.gen_add_code_line("T outv;")
    self.gen_add_code_line("if (row == 3) { outv = (cy*Jw0 + sy*Jw1) / cp; }")
    self.gen_add_code_line("else if (row == 4) { outv = -sy*Jw0 + cy*Jw1; }")
    self.gen_add_code_line("else { outv = (sp / cp) * (cy*Jw0 + sy*Jw1) + Jw2; }")
    self.gen_add_code_line("s_deePos[ind] = outv;")
    self.gen_add_end_control_flow()
    self.gen_add_end_control_flow()
    self.gen_add_sync()
    self.gen_add_end_function()

def _emit_eepose_grad_extraction(self, n, num_ees, ee_off, dee_off, ee_compact = False):
    # Shared eePos extraction: reads the chained ee transform (s_eeTemp at ee_off)
    # and the gradient transform (s_deeTemp at dee_off) and writes s_deePos. When
    # ee_compact is True the ee transform is stored once per ee (slot = deeInd/n)
    # instead of redundantly per (ee, djid) pair (the serial/dense path used the
    # redundant layout and passes ee_off == dee_off, ee_compact == False).
    self.gen_add_parallel_loop("ind",str(6*n*num_ees))
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
    self.gen_add_sync()

def _emit_eepose_grad_compacted_nonserial(self, n, all_ees, num_ees):
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
    self.gen_add_parallel_loop("ind", str(2*16*n*num_ees))
    self.gen_add_code_line("s_deeTemp[ind] = static_cast<T>(0);")
    self.gen_add_end_control_flow()
    self.gen_add_sync()

    # ---- Level 0: seed eeTemp (per ee) and deeTemp (per in-chain pair).
    # eeTemp[ei] = Xhom[ee]; deeTemp[ei,djid] = (djid affects ee ? dXhom[djid] : Xhom[ee]).
    self.gen_add_code_line("// level 0: seed per-ee FK transform")
    self.gen_add_parallel_loop("ind", str(16*num_ees))
    self.gen_add_code_line("int rc = ind % 16; int ei = ind / 16;")
    select_var_vals = [("int", "eeInd", [str(jid) for jid in all_ees])]
    self.gen_add_multi_threaded_select("ind", "<", [str(16*(i+1)) for i in range(num_ees)], select_var_vals)
    self.gen_add_code_line("s_eeTemp[ind] = s_Xhom[16*eeInd + rc];")
    self.gen_add_end_control_flow()
    self.gen_add_sync()

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
    self.gen_add_parallel_loop("ind", str(16*len(seed_pairs)))
    self.gen_add_code_line("int rc = ind % 16; int p = ind / 16;")
    self.gen_add_code_line("const T *s_src = grad_seed_isdx[p] ? &s_dXhom[16*grad_seed_src[p]] : &s_Xhom[16*grad_seed_src[p]];")
    self.gen_add_code_line("s_deeTemp[16*grad_seed_dst[p] + rc] = s_src[rc];")
    self.gen_add_end_control_flow()
    self.gen_add_sync()

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
            self.gen_add_parallel_loop("ind", str(16*len(ee_carry_src)))
            self.gen_add_code_line("int rc = ind % 16; int c = ind / 16;")
            self.gen_add_code_line("s_eeTemp[16*ee_cdst" + sfx + "[c] + rc] = s_eeTemp[16*ee_csrc" + sfx + "[c] + rc];")
            self.gen_add_end_control_flow()
        if dee_carry_src:
            self.gen_add_code_line("static const int dee_csrc" + sfx + "[] = {" + ", ".join(map(str, dee_carry_src)) + "};")
            self.gen_add_code_line("static const int dee_cdst" + sfx + "[] = {" + ", ".join(map(str, dee_carry_dst)) + "};")
            self.gen_add_parallel_loop("ind", str(16*len(dee_carry_src)))
            self.gen_add_code_line("int rc = ind % 16; int c = ind / 16;")
            self.gen_add_code_line("s_deeTemp[16*dee_cdst" + sfx + "[c] + rc] = s_deeTemp[16*dee_csrc" + sfx + "[c] + rc];")
            self.gen_add_end_control_flow()
        if ee_carry_src or dee_carry_src:
            self.gen_add_sync()
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
    _emit_eepose_grad_extraction(self, n, num_ees, ee_final_off, dee_final_off, ee_compact = True)

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

def _emit_eepose_hess_compacted_nonserial(self, n, all_ees, num_ees):
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
    self.gen_add_parallel_loop("ind", str(2*16*n*num_ees))
    self.gen_add_code_line("s_deeTemp[ind] = static_cast<T>(0);")
    self.gen_add_end_control_flow()
    self.gen_add_sync()
    # level 0 eeTemp seed
    self.gen_add_code_line("// level 0: seed per-ee FK transform")
    self.gen_add_parallel_loop("ind", str(16*num_ees))
    self.gen_add_code_line("int rc = ind % 16; int ei = ind / 16;")
    select_var_vals = [("int", "eeInd", [str(jid) for jid in all_ees])]
    self.gen_add_multi_threaded_select("ind", "<", [str(16*(i+1)) for i in range(num_ees)], select_var_vals)
    self.gen_add_code_line("s_eeTemp[ind] = s_Xhom[16*eeInd + rc];")
    self.gen_add_end_control_flow()
    self.gen_add_sync()
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
    self.gen_add_parallel_loop("ind", str(16*len(seed_pairs)))
    self.gen_add_code_line("int rc = ind % 16; int p = ind / 16;")
    self.gen_add_code_line("const T *s_src = hgrad_seed_isdx[p] ? &s_dXhom[16*hgrad_seed_src[p]] : &s_Xhom[16*hgrad_seed_src[p]];")
    self.gen_add_code_line("s_deeTemp[16*hgrad_seed_dst[p] + rc] = s_src[rc];")
    self.gen_add_end_control_flow()
    self.gen_add_sync()
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
        _emit_carry_copy(self, "eec" + sfx, ee_csrc, ee_cdst, "s_eeTemp")
        _emit_carry_copy(self, "dec" + sfx, dee_csrc, dee_cdst, "s_deeTemp")
        if ee_csrc or dee_csrc:
            self.gen_add_sync()
    final_even = (max_len - 1) % 2
    ee_final = 16*num_ees*final_even
    dee_final = 16*n*num_ees*final_even

    # ============ Phase 2: hessian chain (d2eeTemp) ============
    self.gen_add_code_line("// NON-SERIAL: compacted hessian chain (GLASS indexed batched 4x4 GEMM)")
    self.gen_add_parallel_loop("ind", str(2*16*n*n*num_ees))
    self.gen_add_code_line("s_d2eeTemp[ind] = static_cast<T>(0);")
    self.gen_add_end_control_flow()
    self.gen_add_sync()
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
    _emit_seed_copy(self, "hs_d2", seed_d2, "s_d2Xhom", "s_d2eeTemp")
    _emit_seed_copy(self, "hs_di", seed_dxi, "s_dXhom", "s_d2eeTemp")
    _emit_seed_copy(self, "hs_dj", seed_dxj, "s_dXhom", "s_d2eeTemp")
    _emit_seed_copy(self, "hs_x", seed_x, "s_Xhom", "s_d2eeTemp")
    self.gen_add_sync()
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
        _emit_carry_copy(self, "d2c" + sfx, carry_src, carry_dst, "s_d2eeTemp")
        if carry_src:
            self.gen_add_sync()
    d2_final = 16*n*n*num_ees*final_even

    # ============ Phase 3: extraction (rebase pointers to final parity) ============
    self.gen_add_code_line("// rebase to the final-parity buffers, then extract")
    self.gen_add_code_line("s_eeTemp = &s_eeTemp[" + str(ee_final) + "];")
    self.gen_add_code_line("s_deeTemp = &s_deeTemp[" + str(dee_final) + "];")
    self.gen_add_code_line("s_d2eeTemp = &s_d2eeTemp[" + str(d2_final) + "];")
    _emit_eepose_hess_extraction(self, n, num_ees)

def _emit_idx_gemm(self, name, a, b, c, A_base, B_base, C_base):
    if not a:
        return
    self.gen_add_code_line("static const int " + name + "_a[] = {" + ", ".join(map(str, a)) + "};")
    self.gen_add_code_line("static const int " + name + "_b[] = {" + ", ".join(map(str, b)) + "};")
    self.gen_add_code_line("static const int " + name + "_c[] = {" + ", ".join(map(str, c)) + "};")
    self.gen_add_code_line("grid_linalg_indexed_batched_gemm<T, 4>(" + str(len(a)) + ", " + name + "_a, " + name + "_b, " + name + "_c, " + A_base + ", " + B_base + ", " + C_base + ");")

def _emit_seed_copy(self, name, pairs, src_base, dst_base):
    # pairs: list of (dst_slot, src_slot). Copies 4x4 from src_base[src] to dst_base[dst].
    if not pairs:
        return
    self.gen_add_code_line("static const int " + name + "_dst[] = {" + ", ".join(str(p[0]) for p in pairs) + "};")
    self.gen_add_code_line("static const int " + name + "_src[] = {" + ", ".join(str(p[1]) for p in pairs) + "};")
    self.gen_add_parallel_loop("ind", str(16*len(pairs)))
    self.gen_add_code_line("int rc = ind % 16; int p = ind / 16;")
    self.gen_add_code_line(dst_base + "[16*" + name + "_dst[p] + rc] = " + src_base + "[16*" + name + "_src[p] + rc];")
    self.gen_add_end_control_flow()

def _emit_carry_copy(self, name, src, dst, base):
    # carry-forward: copy already-finished slots into the other parity half.
    if not src:
        return
    self.gen_add_code_line("static const int " + name + "_src[] = {" + ", ".join(map(str, src)) + "};")
    self.gen_add_code_line("static const int " + name + "_dst[] = {" + ", ".join(map(str, dst)) + "};")
    self.gen_add_parallel_loop("ind", str(16*len(src)))
    self.gen_add_code_line("int rc = ind % 16; int c = ind / 16;")
    self.gen_add_code_line(base + "[16*" + name + "_dst[c] + rc] = " + base + "[16*" + name + "_src[c] + rc];")
    self.gen_add_end_control_flow()

def gen_end_effector_pose_gradient_device(self, fixed_target_name = ""):
    n = self.robot.get_num_pos()
    nv = self.robot.get_num_vel()
    num_ees = self.robot.get_total_leaf_nodes() if fixed_target_name == "" else 1
    # construct the boilerplate and function definition
    func_params = ["s_deePos is a pointer to shared memory of size 6*NUM_VEL*NUM_EE where NUM_VEL = " + str(nv) + " and NUM_EE = " + str(num_ees), \
                   "s_q is the vector of joint positions", \
                   "d_robotModel is the pointer to the initialized model specific helpers on the GPU (XImats, topology_helpers, etc.)"]
    func_notes = []
    func_def_start = "void end_effector_pose_gradient_device" + ("" if fixed_target_name == "" else "_" + fixed_target_name) + "("
    func_def_middle = "T *s_deePos, const T *s_q, "
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
    (True,  False),   # pick 1: inner_temp + s_deePos -> workspace/global (LITE)
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
    extra_t_buffers = [("s_q", n)] if use_workspace_temp else [("s_q", n), ("s_deePos", 6*nv*num_ees)]
    # Geometric-Jacobian inner doesn't use s_dXhom -> skip its allocation/computation entirely.
    self.gen_XmatsHom_helpers_temp_shared_memory_code(shared_mem_size, include_gradients = False,
                                                      extra_t_buffers = extra_t_buffers,
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
    if not single_call_timing:
        self.gen_add_parallel_loop("k","NUM_TIMESTEPS",block_level = True)
        self.gen_kernel_load_inputs("q",str(n),stride="stride_q")
        if use_workspace_dxhom:
            self.gen_add_code_line("T *s_dXmatsHom = reinterpret_cast<T *>(&d_workspace[k*GRID_WORKSPACE_BYTES_PER_TIMESTEP<T>() + GRID_EE_GRAD_WORKSPACE_DXHOM_OFFSET_BYTES<T>()]);")
        if use_workspace_temp:
            self.gen_add_code_line("T *s_deePos = &d_deePos[k*" + str(6*nv*num_ees) + "];")
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
            self.gen_kernel_save_result("deePos",str(6*nv*num_ees),stride=str(6*nv*num_ees))
        self.gen_add_end_control_flow()
    else:
        self.gen_kernel_load_inputs("q",str(n))
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
        self.gen_anti_licm_input_reload("q",str(n),feedback_from="deePos")
        self.gen_load_update_XmatsHom_helpers_function_call(include_gradients = False)
        updated = {"d_workspace_name": "s_eegrad_temp"} if use_workspace_temp else {}
        updated["s_dXhom_name"] = "nullptr"
        self.gen_end_effector_pose_gradient_inner_function_call(fixed_target_name = fixed_target_name,
            updated_var_names = updated, temp_in_smem_expr = ("false" if use_workspace_temp else "true"))
        self.gen_anti_licm_output_write("deePos")
        self.gen_add_end_control_flow()
        if not use_workspace_temp:
            self.gen_kernel_save_result("deePos",str(6*nv*num_ees))


def gen_end_effector_pose_gradient_kernel(self, single_call_timing = False, fixed_target_name = ""):
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
        _emit_eepose_grad_kernel_body_for_flags(self, n, num_ees, fixed_target_name, uwt, uwd, single_call_timing)
    else:
        tier_names = ("TIER_PERF", "TIER_LITE", "TIER_MINIMAL")
        for tier_idx, (tier_name, pick) in enumerate(zip(tier_names, picks)):
            uwt, uwd = _EE_GRAD_PICK_FLAGS[pick]
            head = "if constexpr (RESOURCE_TIER == " + tier_name + ") {" if tier_idx == 0 else \
                   "else if constexpr (RESOURCE_TIER == " + tier_name + ") {"
            self.gen_add_code_line(head, True)
            _emit_eepose_grad_kernel_body_for_flags(self, n, num_ees, fixed_target_name, uwt, uwd, single_call_timing)
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
                                 "gpuErrchk(cudaMemcpy(hd_data->h_deePos,hd_data->d_deePos,6*NUM_EES*NUM_VEL*" + \
                                    ("num_timesteps*" if not single_call_timing else "") + "sizeof(T),cudaMemcpyDeviceToHost));",
                                 "gpuErrchkKernel();"])
    # finally report out timing if requested
    if single_call_timing:
        from ..algo_registry import single_call_printf_line
        self.gen_add_code_line(single_call_printf_line("ee_pose_gradient"))
    self.gen_add_end_function()

def gen_end_effector_pose_gradient_hessian_output_count(self):
    """Number of T elements in the d2eePos output: 6 * nv * nv * num_ees.

    Output is now d^2(pose)/dv^2 (TANGENT, pinocchio convention). For fixed-base
    nv == nq so the size is unchanged; for floating-base the (nv x nv) block now
    indexes spatial twist components rather than the older non-standard
    quaternion derivatives.
    """
    nv = self.robot.get_num_vel()
    num_ees = self.robot.get_total_leaf_nodes()
    return 6 * nv * nv * num_ees

def gen_end_effector_pose_gradient_hessian_inner_temp_mem_size(self):
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

def gen_end_effector_pose_gradient_hessian_inner_function_call(self, updated_var_names = None,
                                                               out_in_smem_expr = "true"):
    var_names = dict( \
        s_Xhom_name = "s_XmatsHom", \
        s_deePos_name = "s_deePos", \
        s_d2eePos_name = "s_d2eePos", \
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
    code_start = "end_effector_pose_gradient_hessian_inner<T, " + out_in_smem_expr + ">(" + var_names["s_d2eePos_name"] + ", " + var_names["s_deePos_name"] + ", " + var_names["s_q_name"] + ", "
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


def gen_end_effector_pose_gradient_hessian_inner(self):
    """Analytic d^2(pose)/dv^2 of the end-effector pose via per-chain second-
    order Taylor expansion (see docs/d2ee_analytic_derivation.md).

    Replaces the previous FD-on-d/dv-Jacobian implementation (2*nv + 1 gradient
    calls) with a single closed-form pass. Mirrors
    `RBDReference.end_effector_pose_hessian_analytic` (validated to ~1e-9 vs
    the FD oracle on iiwa14-fixed / floating + go2-floating).

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
      3. Emit s_deePos = [J_v; E^{-1}*J_w] using the same closed-form E^{-1}
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
        "s_d2eePos is a pointer to memory of size 6*NUM_VEL*NUM_VEL*NUM_EE where NUM_VEL = " + str(nv) + " and NUM_EE = " + str(num_ees) +
            " (d^2(pose)/dv^2 tangent-space Hessian, pinocchio convention)",
        "s_deePos is a pointer to memory of size 6*NUM_VEL*NUM_EE (the d/dv tangent Jacobian at q)",
        "s_q is the vector of joint positions (size NUM_POS = " + str(nq) + "; kept for signature compatibility, unused by the analytic path)",
        "s_Xhom is the per-joint LOCAL homogeneous-transform buffer (read-only)",
        "s_temp is helper shared memory of size " + str(self.gen_end_effector_pose_gradient_hessian_inner_temp_mem_size()) +
            " (s_Xworld | s_Sworld | s_E_sc; always kept in smem)",
        "d_workspace is the global spill arena s_d2eePos is repointed at when !OUT_IN_SMEM (else unused)",
        "d_robotModel is the model-specific helper struct (kept for signature compatibility, unused by the analytic path)",
        "s_linalg_smem is optional byte-addressed shared memory (reserved; unused by this inner)",
    ]
    func_notes = [
        "Closed-form analytic d2(pose)/dv2; matches RBDReference.end_effector_pose_hessian_analytic (validated ~1e-9 vs the FD oracle on iiwa14 fixed/floating + go2 floating).",
        "Inner-owns scratch placement: the large nv^2 output s_d2eePos moves to d_workspace when !OUT_IN_SMEM. The s_Xworld+s_Sworld+s_E_sc scratch in s_temp stays in smem at every tier.",
    ]
    func_def_start = "void end_effector_pose_gradient_hessian_inner("
    func_def_middle = "T *s_d2eePos, T *s_deePos, const T *s_q, T *s_Xhom, "
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
    # d_workspace when !OUT_IN_SMEM. Reassigning s_d2eePos here keeps every
    # s_d2eePos[...] reference below unchanged.
    self.gen_add_code_line("if constexpr (!OUT_IN_SMEM) { s_d2eePos = d_workspace; } else { (void)d_workspace; }")
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
    self.gen_add_serial_ops()
    for ee_idx in range(num_ees):
        for dof in per_ee_dof_info[ee_idx]:
            vi = dof["vi"]
            j = dof["joint_jid"]
            ang = dof["ang"]
            lin = dof["lin"]
            rev = dof["revolute"]
            base = 16 * (ee_idx * nv + vi)
            ax = ang if rev else lin
            self.gen_add_code_line(
                "// ee=" + str(ee_idx) + " vi=" + str(vi) + " jid=" + str(j) +
                (" rev" if rev else " prism") + " ax_local=" + str(ax))
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
                # Prismatic: S[:3, 3] = axis_world, rest zero
                self.gen_add_code_line("s_Sworld[" + str(base + 12) + "] = axw_0;")
                self.gen_add_code_line("s_Sworld[" + str(base + 13) + "] = axw_1;")
                self.gen_add_code_line("s_Sworld[" + str(base + 14) + "] = axw_2;")
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

    # ===== Step 4: emit s_deePos = [J_v; E^-1 * J_w] from s_Sworld and s_Xworld =====
    # J_w[:, vi] = skew_inv(S_world[:3, :3]) = (axw_x, axw_y, axw_z) (the angular axis)
    #   In our column-major layout: skew[2,1] = wx -> s_Sworld[base+6]; skew[0,2] = wy -> s_Sworld[base+8]; skew[1,0] = wz -> s_Sworld[base+1].
    # J_v[:, vi] = (S_world @ p_ee_world)[:3] - hmm actually J_v[:, vi] = dM[vi][:3, 3] = (S_world @ X_ee)[:3, 3]
    #   = S_world[:3, :3] @ X_ee[:3, 3] + S_world[:3, 3]
    #   = [w_world]_x @ p_ee_world + (pj_world x w_world)     (revolute)
    #   = ω × p_ee_world - ω × pj_world = ω × (p_ee - pj)     ✓ matches the gradient inner
    #   = 0                       + Rw_a @ ax_local           (prismatic)
    # Use the latter expansion directly.
    self.gen_add_code_line("//")
    self.gen_add_code_line("// Step 4: write s_deePos = [J_v ; E(rpy)^-1 * J_w] from S_world")
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
    self.gen_add_code_line("s_deePos[ind] = outv;")
    self.gen_add_end_control_flow()
    self.gen_add_code_line("else {", True)
    # rows 3..5 = E^-1 @ J_w with closed form (same as gradient inner)
    self.gen_add_code_line("T cy = s_E_sc[4*ee + 0]; T sy = s_E_sc[4*ee + 1]; T cp = s_E_sc[4*ee + 2]; T sp = s_E_sc[4*ee + 3];")
    self.gen_add_code_line("T outv;")
    self.gen_add_code_line("if (row == 3) { outv = (cy*wx + sy*wy) / cp; }")
    self.gen_add_code_line("else if (row == 4) { outv = -sy*wx + cy*wy; }")
    self.gen_add_code_line("else { outv = (sp / cp) * (cy*wx + sy*wy) + wz; }")
    self.gen_add_code_line("s_deePos[ind] = outv;")
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
    #  - We write d2M values into s_d2eePos. We'll do the rpy rows in a second
    #    pass once H_w / E^-1 are known.
    #
    # We pre-zero s_d2eePos (covers out-of-chain entries and is also needed
    # because the intra-joint case only writes the unique (i, j) ordered pair).
    self.gen_add_code_line("//")
    self.gen_add_code_line("// Step 5a: zero the full d2eePos output (out-of-chain pairs stay zero)")
    self.gen_add_code_line("//")
    self.gen_add_parallel_loop("ind", str(6 * nv * nv * num_ees))
    self.gen_add_code_line("s_d2eePos[ind] = static_cast<T>(0);")
    self.gen_add_end_control_flow()
    self.gen_add_sync()

    # Step 5b: per-pair d2M -> H_xyz + (temporarily, into d2eePos rpy rows) d2R_R^T
    # We use the rpy rows (c=3,4,5 of s_d2eePos) as a SCRATCH BUFFER for d2R_R^T's
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
    self.gen_add_serial_ops()
    for ee_idx in range(num_ees):
        ee_jid = anchors[ee_idx]
        chain_dofs = per_ee_dof_info[ee_idx]
        chain_jids = chains[ee_idx]
        intra_pairs = intra_joint_pairs_per_ee[ee_idx]
        # Index by (vi_a, vi_b) for quick lookup
        intra_pair_lookup = {(p[0], p[1]): p for p in intra_pairs}
        intra_pair_lookup.update({(p[1], p[0]): p for p in intra_pairs})

        # Pair iteration: per (vi_i, vi_j), with i,j enumerated over chain DOFs.
        for di in chain_dofs:
            vi = di["vi"]
            for dj in chain_dofs:
                vj = dj["vi"]
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
                self.gen_add_code_line(
                    "// ee=" + str(ee_idx) + " pair (vi=" + str(vi) + ", vj=" + str(vj) +
                    ", a=" + str(a) + ", b=" + str(b) + ")")
                self.gen_add_code_line("{", True)
                # Read X_ee column 3 (p_ee) and X_ee R block (used in skew/inv)
                self.gen_add_code_line("T pex = s_Xworld[" + str(16*ee_jid + 12) + "];")
                self.gen_add_code_line("T pey = s_Xworld[" + str(16*ee_jid + 13) + "];")
                self.gen_add_code_line("T pez = s_Xworld[" + str(16*ee_jid + 14) + "];")
                # Read S_i_world and S_j_world skew axis components and column 3
                si_base = 16 * (ee_idx * nv + vi)
                sj_base = 16 * (ee_idx * nv + vj)
                # We need the top-left 3x3 product (S_prox @ S_dist)[:3,:3] and
                # the column-3 expansion (S_prox @ S_dist @ X_ee)[:3, 3].
                # Compute it for the proximal/distal ordering.
                if a == b:
                    # Same chain joint. For single-DOF joints, B_local = 0 (revolute or
                    # prismatic intra-pair doesn't exist except for the diagonal where
                    # B = A_x^2 for revolute, 0 for prismatic). For multi-DOF (floating
                    # base) joints, use the closed form B_world below.
                    # Diagonal vi == vj case (always present): B_local at v=0 = A_x^2.
                    # For revolute, A = [[ω_×, 0]; 0] so A^2 = [[ω_×^2, 0]; 0]: this is
                    # the "centripetal" term.
                    # We handle this via the intra-pair lookup (which includes c_a == c_b).
                    if vi == vj:
                        # Diagonal: B_local for column c_a alone.
                        dof = di  # same as dj
                        c_a = dof["S_col"]
                        c_b = dof["S_col"]
                        is_intra_multi = (vi, vj) in intra_pair_lookup
                        # Even for single-DOF joints we hit this code path on the
                        # diagonal. Handle revolute / prismatic / floating-base
                        # uniformly via the on-the-fly B formula.
                        _emit_d2M_same_joint_block(self, dof, dof,
                                                   ee_idx, ee_jid, vi, vj, nv, num_ees,
                                                   chain_jids[a],
                                                   si_base, sj_base)
                    else:
                        # Off-diagonal same-joint pair: only exists for multi-DOF
                        # (floating base) joints.
                        if (vi, vj) in intra_pair_lookup:
                            # Look up which is c_a, c_b by S_col
                            _emit_d2M_same_joint_block(self, di, dj,
                                                       ee_idx, ee_jid, vi, vj, nv, num_ees,
                                                       chain_jids[a],
                                                       si_base, sj_base)
                        else:
                            # Shouldn't happen (a == b but DOFs not in same joint)
                            # Emit a zero-write defensively (rows 0..5 already zero
                            # from the bulk zero in Step 5a).
                            pass
                else:
                    # Different chain joints: pick the proximal/distal.
                    # Proximal = the one with smaller chain_pos.
                    if a < b:
                        prox_base = si_base
                        dist_base = sj_base
                    else:
                        prox_base = sj_base
                        dist_base = si_base
                    # Read all 4x4 entries of S_prox and S_dist (column-major)
                    # (Skip the bottom row — known zero)
                    _emit_d2M_cross_joint_block(self, prox_base, dist_base,
                                                ee_idx, vi, vj, nv, num_ees, si_base, sj_base)
                self.gen_add_end_control_flow()
    self.gen_add_end_control_flow()
    self.gen_add_sync()

    # ===== Step 6: rpy chain rule on rows 3..5 =====
    # Currently rows 3..5 hold H_w[:, i, j]. We want H_rpy[:, i, j] =
    # (dEinv/dv_j) * J_w[:, i] + Einv * H_w[:, i, j], with closed-form Einv
    # and dE/drpy. Recall the gradient inner already wrote drpy/dv = Einv*J_w
    # to s_deePos rows 3..5: but we need that for each (ee, vj).
    #
    # NOTE on race-safety: each thread owns a single (ee, vi, vj) cell and reads
    # all 3 components of H_w[:, vi, vj] from rows 3..5 of s_d2eePos BEFORE
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
    # drpy_j (from s_deePos rows 3,4,5)
    self.gen_add_code_line("int dee_base_j = ee * " + str(6 * nv) + " + 6 * vj;")
    self.gen_add_code_line("T drpy_j0 = s_deePos[dee_base_j + 3];")
    self.gen_add_code_line("T drpy_j1 = s_deePos[dee_base_j + 4];")
    self.gen_add_code_line("T drpy_j2 = s_deePos[dee_base_j + 5];")
    # drpy_i (also from s_deePos)
    self.gen_add_code_line("int dee_base_i = ee * " + str(6 * nv) + " + 6 * vi;")
    self.gen_add_code_line("T drpy_i0 = s_deePos[dee_base_i + 3];")
    self.gen_add_code_line("T drpy_i1 = s_deePos[dee_base_i + 4];")
    self.gen_add_code_line("T drpy_i2 = s_deePos[dee_base_i + 5];")
    # READ ALL 3 H_w components BEFORE any rpy writes (critical for correctness:
    # we will overwrite rows 3..5 below; if we read after writing the race would
    # silently corrupt the other two components in this thread's row).
    self.gen_add_code_line("int hw_base = ee * " + str(6 * nv * nv) + " + 3 * " + str(nv * nv) + " + vi * " + str(nv) + " + vj;")
    self.gen_add_code_line("T Hw_x = s_d2eePos[hw_base + 0 * " + str(nv * nv) + "];")
    self.gen_add_code_line("T Hw_y = s_d2eePos[hw_base + 1 * " + str(nv * nv) + "];")
    self.gen_add_code_line("T Hw_z = s_d2eePos[hw_base + 2 * " + str(nv * nv) + "];")
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
    self.gen_add_code_line("s_d2eePos[out_base + 0 * " + str(nv * nv) + "] = H_rpy_0;")
    self.gen_add_code_line("s_d2eePos[out_base + 1 * " + str(nv * nv) + "] = H_rpy_1;")
    self.gen_add_code_line("s_d2eePos[out_base + 2 * " + str(nv * nv) + "] = H_rpy_2;")
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
    self.gen_add_code_line("T avg = static_cast<T>(0.5) * (s_d2eePos[idx_ij] + s_d2eePos[idx_ji]);")
    self.gen_add_code_line("s_d2eePos[idx_ij] = avg;")
    self.gen_add_code_line("if (i != j) { s_d2eePos[idx_ji] = avg; }")
    self.gen_add_end_control_flow()
    self.gen_add_end_control_flow()
    self.gen_add_sync()

    self.gen_add_end_function()


def _emit_d2M_cross_joint_block(self, prox_base, dist_base,
                                 ee_idx, vi, vj, nv, num_ees, si_base, sj_base):
    """Emit explicit per-pair scalar code for d2M = S_prox * S_dist * X_ee in a
    cross-joint pair (chain ordering a < b after prox/dist resolution).

    Writes:
      - rows 0..2 of s_d2eePos[(ee, vi, vj)]   = (d2M @ ee_offset)[:3] with
        ee_offset = [0, 0, 0, 1] -> just column 3 of d2M.
      - rows 3..5 of s_d2eePos[(ee, vi, vj)]   = H_w[:, vi, vj] =
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
    # Write into s_d2eePos: idx = ee*6*nv*nv + c*nv*nv + vi*nv + vj
    base = "(" + str(ee_idx * 6 * nv * nv) + " + " + str(vi * nv + vj) + ")"
    self.gen_add_code_line("s_d2eePos[" + base + " + 0 * " + str(nv*nv) + "] = Hxyz_x;")
    self.gen_add_code_line("s_d2eePos[" + base + " + 1 * " + str(nv*nv) + "] = Hxyz_y;")
    self.gen_add_code_line("s_d2eePos[" + base + " + 2 * " + str(nv*nv) + "] = Hxyz_z;")
    self.gen_add_code_line("s_d2eePos[" + base + " + 3 * " + str(nv*nv) + "] = HW_x;")
    self.gen_add_code_line("s_d2eePos[" + base + " + 4 * " + str(nv*nv) + "] = HW_y;")
    self.gen_add_code_line("s_d2eePos[" + base + " + 5 * " + str(nv*nv) + "] = HW_z;")


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
      - rows 0..2 of s_d2eePos[(ee, vi, vj)] = (d2M)[:3, 3]
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
    self.gen_add_code_line("s_d2eePos[" + base + " + 0 * " + str(nv*nv) + "] = Hxyz_x;")
    self.gen_add_code_line("s_d2eePos[" + base + " + 1 * " + str(nv*nv) + "] = Hxyz_y;")
    self.gen_add_code_line("s_d2eePos[" + base + " + 2 * " + str(nv*nv) + "] = Hxyz_z;")
    self.gen_add_code_line("s_d2eePos[" + base + " + 3 * " + str(nv*nv) + "] = HW_x;")
    self.gen_add_code_line("s_d2eePos[" + base + " + 4 * " + str(nv*nv) + "] = HW_y;")
    self.gen_add_code_line("s_d2eePos[" + base + " + 5 * " + str(nv*nv) + "] = HW_z;")

def gen_end_effector_pose_gradient_hessian_device(self):
    n = self.robot.get_num_pos()
    nv = self.robot.get_num_vel()
    num_ees = self.robot.get_total_leaf_nodes()
    inner_temp_size = self.gen_end_effector_pose_gradient_hessian_inner_temp_mem_size()
    output_count = self.gen_end_effector_pose_gradient_hessian_output_count()
    # construct the boilerplate and function definition
    func_params = ["s_d2eePos is a pointer to shared memory of size 6*NUM_VEL*NUM_VEL*NUM_EE where NUM_VEL = " + str(nv) + " and NUM_EE = " + str(num_ees) + " (d^2/dv^2 tangent, pinocchio convention)", \
                   "s_deePos is a pointer to shared memory of size 6*NUM_VEL*NUM_EE (d/dv tangent Jacobian)", \
                   "s_q is the vector of joint positions", \
                   "d_robotModel is the pointer to the initialized model specific helpers on the GPU (XImats, topology_helpers, etc.)", \
                   "d_workspace is the global scratch buffer; size D2EE_DEVICE_INLINE_WORKSPACE_BYTES<T, RESOURCE_TIER>() bytes (= 0 at TIER_PERF, " + str(output_count) + "*sizeof(T) at TIER_LITE+). Pass nullptr at TIER_PERF"]
    func_notes = ["Inline-CUDA users: at TIER_LITE/TIER_MINIMAL the large s_d2eePos output (~" + str(output_count) + "*sizeof(T) bytes) moves from shared memory to d_workspace, freeing smem for the caller's outer kernel"]
    func_def_start = "void end_effector_pose_gradient_hessian_device("
    func_def_middle = "T *s_d2eePos, T *s_deePos, const T *s_q, "
    func_def_end = "const robotModel<T> *d_robotModel, T *d_workspace = nullptr) {"
    func_def = func_def_start + func_def_middle + func_def_end
    # then generate the code
    self.gen_add_func_doc("Computes the Hessian (and Jacobian) of the End Effector Pose with respect to generalized velocity (d^2/dv^2 tangent, pinocchio convention)",\
                          func_notes,func_params,None)
    self.gen_add_code_line("template <typename T, int RESOURCE_TIER = TIER_PERF>")
    self.gen_add_code_line("__device__")
    self.gen_add_code_line(func_def, True)
    # Smem arena: s_temp is always the inner-temp size. The output s_d2eePos lives
    # in smem at TIER_PERF (carved from the arena tail) and in d_workspace at
    # TIER_LITE/MINIMAL (inner repoints internally). Note: the geometric-Jacobian
    # path uses ONLY s_Xhom (LOCAL transforms); s_dXmatsHom and s_d2XmatsHom are
    # no longer needed (saves substantial smem on big robots).
    self.gen_XmatsHom_helpers_temp_shared_memory_code(inner_temp_size, include_gradients = False, include_hessians = False,
                                                      include_linalg_scratch = True,
                                                      linalg_scratch_bytes = "GRID_EE_LINALG_SHARED_BYTES<T>()")
    # At TIER_PERF s_d2eePos is allocated by the caller; at LITE/MINIMAL it's
    # the inner's job to repoint via OUT_IN_SMEM=false + d_workspace.
    # then load Xhom (Jacobian only needs local transforms) and run the algo
    self.gen_load_update_XmatsHom_helpers_function_call(include_gradients = False, include_hessians = False)
    # Inner-owns placement: pass d_workspace + the per-tier flag. When the flag is
    # false the inner repoints s_d2eePos at d_workspace.
    self.gen_end_effector_pose_gradient_hessian_inner_function_call(
        updated_var_names = {"d_workspace_name": "d_workspace", "s_Xhom_name": "s_XmatsHom", "d_robotModel_name": "d_robotModel"},
        out_in_smem_expr = "D2EE_OUT_IN_SMEM<RESOURCE_TIER>()")
    self.gen_add_end_function()

_D2EE_PICK_FLAGS = [
    # use_workspace_output (s_d2eePos lives in d_workspace?)
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
    Used by gen_end_effector_pose_gradient_hessian_kernel to emit either a
    single body (collapsed picks) or three branched bodies (divergent picks)."""
    nv = self.robot.get_num_vel()
    output_count = self.gen_end_effector_pose_gradient_hessian_output_count()
    inner_temp_size = self.gen_end_effector_pose_gradient_hessian_inner_temp_mem_size()
    extra_t_buffers = [("s_q", n)] if use_workspace_output else [("s_q", n), ("s_d2eePos", output_count), ("s_deePos", 6*nv*num_ees)]
    self.gen_XmatsHom_helpers_temp_shared_memory_code(inner_temp_size, include_gradients = False, include_hessians = False,
                                                      extra_t_buffers = extra_t_buffers,
                                                      include_linalg_scratch = True,
                                                      linalg_scratch_bytes = "GRID_EE_LINALG_SHARED_BYTES<T>()")
    out_in_smem_expr = "false" if use_workspace_output else "true"
    if not single_call_timing:
        self.gen_add_parallel_loop("k","NUM_TIMESTEPS",block_level = True)
        self.gen_kernel_load_inputs("q",str(n),stride="stride_q")
        if use_workspace_output:
            self.gen_add_code_line("T *s_d2eePos = nullptr;  // inner repoints at d_workspace slice")
            self.gen_add_code_line("T *s_deePos = &d_deePos[k*" + str(6*nv*num_ees) + "];")
            self.gen_add_code_line("// Use d_d2eePos directly as the spill target so the inner writes into the persistent output buffer (one allocation, no extra copy).")
            self.gen_add_code_line("T *s_d2eePos_ws = &d_d2eePos[k*" + str(output_count) + "];")
        else:
            self.gen_add_code_line("(void)d_workspace;")
        self.gen_add_code_line("// compute")
        self.gen_load_update_XmatsHom_helpers_function_call(include_gradients = False, include_hessians = False)
        updated = {"s_Xhom_name": "s_XmatsHom", "d_robotModel_name": "d_robotModel"}
        if use_workspace_output:
            updated["d_workspace_name"] = "s_d2eePos_ws"
        self.gen_end_effector_pose_gradient_hessian_inner_function_call(
            updated_var_names = updated, out_in_smem_expr = out_in_smem_expr)
        self.gen_add_sync()
        if not use_workspace_output:
            self.gen_kernel_save_result("d2eePos",str(output_count),stride=str(output_count))
            self.gen_kernel_save_result("deePos",str(6*nv*num_ees),stride=str(6*nv*num_ees))
        else:
            # gradient still needs the smem -> global copy; the Hessian was already written to d_d2eePos directly via s_d2eePos_ws.
            self.gen_kernel_save_result("deePos",str(6*nv*num_ees),stride=str(6*nv*num_ees))
        self.gen_add_end_control_flow()
    else:
        self.gen_kernel_load_inputs("q",str(n))
        if use_workspace_output:
            self.gen_add_code_line("T *s_d2eePos = nullptr;  // inner repoints at d_workspace")
            self.gen_add_code_line("T *s_deePos = d_deePos;")
            self.gen_add_code_line("T *s_d2eePos_ws = d_d2eePos;  // use output buffer directly as spill target")
        else:
            self.gen_add_code_line("(void)d_workspace;")
        self.gen_add_code_line("// compute with NUM_TIMESTEPS as NUM_REPS for timing")
        self.gen_add_code_line("for (int rep = 0; rep < NUM_TIMESTEPS; rep++){", True)
        self.gen_anti_licm_input_reload("q",str(n),feedback_from="d2eePos")
        self.gen_load_update_XmatsHom_helpers_function_call(include_gradients = False, include_hessians = False)
        updated = {"s_Xhom_name": "s_XmatsHom", "d_robotModel_name": "d_robotModel"}
        if use_workspace_output:
            updated["d_workspace_name"] = "s_d2eePos_ws"
        self.gen_end_effector_pose_gradient_hessian_inner_function_call(
            updated_var_names = updated, out_in_smem_expr = out_in_smem_expr)
        self.gen_anti_licm_output_write("d2eePos")
        self.gen_add_end_control_flow()
        if not use_workspace_output:
            self.gen_kernel_save_result("d2eePos",str(output_count))
            self.gen_kernel_save_result("deePos",str(6*nv*num_ees))
        else:
            self.gen_kernel_save_result("deePos",str(6*nv*num_ees))


def gen_end_effector_pose_gradient_hessian_kernel(self, single_call_timing = False):
    n = self.robot.get_num_pos()
    num_ees = self.robot.get_total_leaf_nodes()
    func_params = ["d_d2eePos is the vector of end effector pose Hessians (6 x nv x nv per ee)", \
                   "d_deePos is the vector of end effector pose Jacobians (6 x nv per ee)", \
                   "d_workspace is the generated global spill workspace", \
                   "d_q is the vector of joint positions", \
                   "stride_q is the stide between each q", \
                   "d_robotModel is the pointer to the initialized model specific helpers on the GPU (XImats, topology_helpers, etc.)", \
                   "num_timesteps is the length of the trajectory points we need to compute over (or overloaded as test_iters for timing)"]
    func_notes = ["Output d^2(pose)/dv^2 is in tangent-space convention (d/dv), shape 6 x nv x nv per ee, C-order. Matches pinocchio."]
    func_def_start = "void end_effector_pose_gradient_hessian_kernel(T *d_d2eePos, T *d_deePos, unsigned char *d_workspace, const T *d_q, const int stride_q, "
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
    # tier's spill flag. Smem-bytes constexpr D2EE_POS_DYNAMIC_SHARED_MEM_BYTES<T,TIER>()
    # is already tier-aware.
    picks = getattr(self, "d2ee_spill_tier_3way", (0, 0, 0))
    if picks[0] == picks[1] == picks[2]:
        uwo = _D2EE_PICK_FLAGS[picks[0]]
        _emit_d2ee_kernel_body_for_flags(self, n, num_ees, uwo, single_call_timing)
    else:
        tier_names = ("TIER_PERF", "TIER_LITE", "TIER_MINIMAL")
        for tier_idx, (tier_name, pick) in enumerate(zip(tier_names, picks)):
            uwo = _D2EE_PICK_FLAGS[pick]
            head = "if constexpr (RESOURCE_TIER == " + tier_name + ") {" if tier_idx == 0 else \
                   "else if constexpr (RESOURCE_TIER == " + tier_name + ") {"
            self.gen_add_code_line(head, True)
            _emit_d2ee_kernel_body_for_flags(self, n, num_ees, uwo, single_call_timing)
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
    self.gen_add_func_doc("Computes the Hessian (and Jacobian) of the End Effector Pose with respect to generalized velocity (d^2/dv^2 tangent, pinocchio convention)",\
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
    # No L2 persistence: at LITE/MINIMAL the d2eePos spill target IS the output
    # buffer (d_d2eePos), which is written once and read once -- no benefit from
    # L2 pinning.
    self.gen_add_code_lines(func_call_code)
    if not compute_only:
        # then transfer memory back
        self.gen_add_code_lines(["// finally transfer the result back", \
                                 "gpuErrchk(cudaMemcpy(hd_data->h_deePos,hd_data->d_deePos,6*NUM_EES*NUM_VEL*" + \
                                    ("num_timesteps*" if not single_call_timing else "") + "sizeof(T),cudaMemcpyDeviceToHost));",
                                 "gpuErrchk(cudaMemcpy(hd_data->h_d2eePos,hd_data->d_d2eePos,6*NUM_EES*NUM_VEL*NUM_VEL*" + \
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
        self.gen_end_effector_pose_gradient_hessian_inner()
        # then generate the device wrappers
        self.gen_end_effector_pose_gradient_hessian_device()
        # then generate the kernels
        self.gen_end_effector_pose_gradient_hessian_kernel(True)
        self.gen_end_effector_pose_gradient_hessian_kernel(False)
        # then the host launch wrappers
        self.gen_end_effector_pose_gradient_hessian_host(0)
        self.gen_end_effector_pose_gradient_hessian_host(1)
        self.gen_end_effector_pose_gradient_hessian_host(2)

    if include_pose or include_gradient or include_hessian:
        self.gen_X_single_thread(fixed_target_name = fixed_target_name)
        self.gen_X_warp(fixed_target_name = fixed_target_name)
