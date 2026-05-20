import numpy as np
import copy
#np.set_printoptions(precision=4, suppress=True, linewidth = 100)

def gen_crba_inner_temp_mem_size(self):
    if self.robot.floating_base:
        NJ = self.robot.get_num_joints()
        return 36*NJ + 36 + 36 + 6
    n = self.robot.get_num_pos()
    return 140*n

def gen_crba_inner_function_call(self, use_thread_group = False, updated_var_names = None):
    var_names = dict( \
        s_M_name = "s_M", \
        s_q_name = "s_q", \
        s_qd_name = "s_qd", \
        s_temp_name = "s_temp", \
        #s_XI = "s_XImats", \
        gravity_name = "gravity"
    )
    if updated_var_names is not None:
        for key,value in updated_var_names.items():
            var_names[key] = value
    crba_code_start = "crba_inner<T>(" + var_names["s_M_name"] + ", " +  var_names["s_q_name"] + ", " + var_names["s_qd_name"] + ", " 
    crba_code_end = var_names["s_temp_name"] + ", " + var_names["gravity_name"] + ");"
    if use_thread_group:
        id_code_start = id_code_start.replace("(","(tgrp, ")
    crba_code_middle = self.gen_insert_helpers_function_call()
    crba_code = crba_code_start + crba_code_middle + crba_code_end
    self.gen_add_code_line(crba_code)


def gen_crba_inner(self, use_thread_group = False):
    if self.robot.floating_base:
        return gen_crba_inner_floating(self, use_thread_group)
    
    n = self.robot.get_num_joints()
    n_bfs_levels = self.robot.get_max_bfs_level() + 1
    has_linear_axis = any(self.robot.get_S_index_by_id(jid) >= 3 for jid in range(n))
    imat_offset = n if has_linear_axis else 7

    #construct the boilerplate and function definition
    func_params = [ "s_q is the vector of joint positions", \
                    "s_qd is the vector of joint velocities", \
                    "s_M is a pointer to the matrix of inertia" \
                    "s_XI is the pointer to the transformation and inertia matricies ", \
                    "s_temp is a pointer to helper shared memory of size 6*NUM_JOINTS = " + \
                            str(self.gen_crba_inner_temp_mem_size()), \
                    "gravity is the gravity constant"]
    func_notes = [] #insert notes abt function 
    func_def_start = "void crba_inner("
    func_def_middle = "T *s_M, const T *s_q, const T *s_qd, "
    func_def_end = "T *s_temp, const T gravity) {"
    if use_thread_group:
        func_def_start = func_def_start.replace("(", "(cgrps::thread_group tgrp, ")
        func_params.insert(0,"tgrp is the handle to the thread_group running this function")

    func_def_middle, func_params = self.gen_insert_helpers_func_def_params(func_def_middle, func_params, -2)
    func_def = func_def_start + func_def_middle + func_def_end
    self.gen_add_func_doc("Compute the Composite Rigid Body Algorithm", func_notes, func_params, None)
    self.gen_add_code_line("template <typename T>")
    self.gen_add_code_line("__device__")
    self.gen_add_code_line(func_def, True)
    temp_size = self.gen_crba_inner_temp_mem_size()
    self.gen_linalg_smem_setup(temp_size)


    # first clear the matrix
    self.gen_add_parallel_loop("i",str(n*n),use_thread_group)
    self.gen_add_code_line('s_M[i] = static_cast<T>(0);')
    self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)
    #deal with like memory for variables --> memory is taken care of in device and host 
    alpha_offset = 0
    beta_offset = alpha_offset + 36*n
    fh_offset = beta_offset + 36*n 
    parent_offset = fh_offset + 6*n #bc j is 1 int 
    jid_offset = parent_offset + n
    self.gen_add_code_line("T *alpha = &s_temp[" + str(alpha_offset) + "];")
    self.gen_add_code_line("T *beta = &s_temp[" + str(beta_offset) + "];")
    self.gen_add_code_line("T *s_fh = &s_temp[" + str(fh_offset) + "];")
    # self.gen_add_code_line("T *s_parent_inds = &s_temp[" + str(parent_offset) + "];")
    self.gen_add_code_line("T *s_jid_list = &s_temp[" + str(jid_offset) + "];")

    self.gen_add_code_line("//")
    self.gen_add_code_line("// first loop (split into 2 parallel loops in bfs loop)")
    self.gen_add_code_line("// each bfs level runs in parallel")
    self.gen_add_code_line("//")
 
    for bfs_level in range(n_bfs_levels-1,0,-1):
        inds = self.robot.get_ids_by_bfs_level(bfs_level)
 
        joint_names = [self.robot.get_joint_by_id(indj).get_name() for indj in inds]
        link_names = [self.robot.get_link_by_id(indl).get_name() for indl in inds]

        self.gen_add_code_line("// pass updates where bfs_level is " + str(bfs_level))
        self.gen_add_code_line("//     joints are: " + ", ".join(joint_names))
        self.gen_add_code_line("//     links are: " + ", ".join(link_names))

        parent_ind_cpp, S_ind_cpp = self.gen_topology_helpers_pointers_for_cpp(inds, NO_GRAD_FLAG = True)

        if len(inds) > 1 and has_linear_axis:
            for jid in inds:
                parent_ind = self.robot.get_parent_id(jid)
                self.gen_add_code_line(f"grid_linalg_gemm<T,6,6,6,false,true>(&s_XImats[{36*jid}], &s_XImats[{36*(jid+n)}], &alpha[{36*jid}], static_cast<T>(1), static_cast<T>(0), s_linalg_smem);")
                self.gen_add_code_line(f"grid_linalg_gemm<T,6,6,6>(&alpha[{36*jid}], &s_XImats[{36*jid}], &s_XImats[{36*(parent_ind+n)}], static_cast<T>(1), static_cast<T>(1), s_linalg_smem);")

        elif len(inds) > 1:
            for jid in inds:
                parent_ind = self.robot.get_parent_id(jid)
                self.gen_add_code_line(f"grid_linalg_gemm<T,6,6,6,false,true>(&s_XImats[{36*jid}], &s_XImats[{36*(jid+n)}], &alpha[{36*jid}], static_cast<T>(1), static_cast<T>(0), s_linalg_smem);")
                self.gen_add_code_line(f"grid_linalg_gemm<T,6,6,6>(&alpha[{36*jid}], &s_XImats[{36*jid}], &s_XImats[{36*(parent_ind+n)}], static_cast<T>(1), static_cast<T>(1), s_linalg_smem);")

        else:
            jid = inds[0]
            parent_ind = self.robot.get_parent_id(jid)
            self.gen_add_code_line(f"grid_linalg_gemm<T,6,6,6,false,true>(&s_XImats[{36*jid}], &s_XImats[{36*(jid+n)}], &alpha[{36*jid}], static_cast<T>(1), static_cast<T>(0), s_linalg_smem);")
            self.gen_add_code_line(f"grid_linalg_gemm<T,6,6,6>(&alpha[{36*jid}], &s_XImats[{36*jid}], &s_XImats[{36*(parent_ind+n)}], static_cast<T>(1), static_cast<T>(1), s_linalg_smem);")

    # Calculation of M[ind,ind]
    self.gen_add_code_line("//")
    self.gen_add_code_line("// Calculation of M[ind, ind] ")
    self.gen_add_code_line("//")
    _, S_ind_cpp = self.gen_topology_helpers_pointers_for_cpp(NO_GRAD_FLAG = True)
    S_sign_cpp = self.gen_topology_S_sign_for_cpp()
    
    self.gen_add_parallel_loop("jid",str(n),use_thread_group)

    ImatOffset = 36*n   # Offset in XImats to Imats
    self.gen_add_code_line(f"s_M[jid+jid*{n}] = s_XImats[{ImatOffset} + 36*jid + 6*{S_ind_cpp} + {S_ind_cpp}];") # take the S_ind row and S_ind column of appropriate Imat
    self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)
    

    self.gen_add_code_line("//")
    self.gen_add_code_line("// Calculation of M[ind, parent]")
    self.gen_add_code_line("//")

    # initialize fh as (XS)^T
    self.gen_add_parallel_loop('i',str(n*6),use_thread_group)
    self.gen_add_code_line('int jid = i / 6; int ind = i % 6;')
    self.gen_add_code_line(f's_fh[i] = ({S_sign_cpp}) * s_XImats[{ImatOffset} + 36*jid + 6*{S_ind_cpp} + ind];')
    self.gen_add_end_control_flow()

    # M[jid, parent_id] = dot product between (S_lambda^T * X_lambda^t) * IS
    # IS = I[6*S_ind]
    # (S_lambda^T * X_lambda^t) = X_lambda[S_ind row]

    # loop over all joints and parent chains in parallel
    max_ancestors = self.robot.get_max_num_ancestors()
    self.gen_add_parallel_loop("jid",str(n),use_thread_group)
    parent_chain = '{' + "-1, " * (max_ancestors-1) + "-1}"
    self.gen_add_code_line(f"int jid_parents[] = {parent_chain};")
    self.gen_add_code_line('int num_parents = 0;')


    self.gen_add_code_line('switch (jid) {', True)
    for jid in range(n):
        self.gen_add_code_line(f"case {jid}:", True)
        parent_chain = self.robot.get_ancestors_by_id(jid)
        length = len(parent_chain)

        for i, parent_ind in enumerate(parent_chain): 
            self.gen_add_code_line(f"jid_parents[{i}] = {parent_ind};")
        self.gen_add_code_line(f"num_parents += {length};")
        self.gen_add_code_line("break;")
        self.indent_level -= 1
    self.gen_add_end_control_flow()

    self.gen_add_code_line("T s_alpha[6];")
    self.gen_add_code_line('for (int i = 0; i < num_parents; i++) {', True)
    
    self.gen_add_code_line('int X_ind = i==0 ? jid : jid_parents[i-1];')
    self.gen_add_code_line('for (int k = 0; k < 6; k++) s_alpha[k] = s_fh[jid*6+k];')
    self.gen_add_code_line('for (int k = 0; k < 6; k++) s_fh[jid*6 + k] = dot_prod<T,6,1,1>(&s_XImats[36*X_ind+k*6], &s_alpha[0]);')

    self.gen_add_code_line('int parent_ind = jid_parents[i];')
    if has_linear_axis:
        _, parent_S_ind_cpp = self.gen_topology_helpers_pointers_for_cpp(
            list(range(n)),
            updated_var_names=dict(jid_name="parent_ind"),
            NO_GRAD_FLAG=True,
        )
        parent_S_sign_cpp = self.gen_topology_S_sign_for_cpp(
            list(range(n)),
            updated_var_names=dict(jid_name="parent_ind"),
        )
        self.gen_add_code_line(f"s_M[jid*{n} + parent_ind] = ({parent_S_sign_cpp}) * s_fh[jid*6 + {parent_S_ind_cpp}];")
    else:
        self.gen_add_code_line(f"s_M[jid*{n} + parent_ind] = s_fh[jid*6 + {S_ind_cpp}];")
    self.gen_add_code_line(f"s_M[parent_ind*{n} + jid] = s_M[jid*{n} + parent_ind];") # M symmetric

    self.gen_add_end_control_flow()
    self.gen_add_end_control_flow()

    self.gen_add_sync(use_thread_group)
    self.gen_add_end_function()


def gen_crba_inner_floating(self, use_thread_group = False):
    NJ = self.robot.get_num_joints()
    nv = self.robot.get_num_vel()
    ICOffset = 0
    alphaOffset = 36 * NJ
    betaOffset = alphaOffset + 36
    fhOffset = betaOffset + 36

    func_params = [ "s_q is the vector of joint positions", \
                    "s_qd is the vector of joint velocities", \
                    "s_M is a pointer to the matrix of inertia", \
                    "s_XI is the pointer to the transformation and inertia matricies", \
                    "s_temp is a pointer to helper shared memory of size " + \
                            str(self.gen_crba_inner_temp_mem_size()), \
                    "gravity is the gravity constant"]
    func_notes = ["Floating-base CRBA keeps composite inertias in body order and writes a public-order NUM_VEL x NUM_VEL mass matrix."]
    func_def_start = "void crba_inner("
    func_def_middle = "T *s_M, const T *s_q, const T *s_qd, "
    func_def_end = "T *s_temp, const T gravity) {"
    if use_thread_group:
        func_def_start = func_def_start.replace("(", "(cgrps::thread_group tgrp, ")
        func_params.insert(0,"tgrp is the handle to the thread_group running this function")
    func_def_middle, func_params = self.gen_insert_helpers_func_def_params(func_def_middle, func_params, -2)
    func_def = func_def_start + func_def_middle + func_def_end
    self.gen_add_func_doc("Compute the Floating-Base Composite Rigid Body Algorithm", func_notes, func_params, None)
    self.gen_add_code_line("template <typename T>")
    self.gen_add_code_line("__device__")
    self.gen_add_code_line(func_def, True)
    temp_size = self.gen_crba_inner_temp_mem_size()
    self.gen_linalg_smem_setup(temp_size)

    self.gen_add_code_line("// Initialize IC = I and clear H")
    self.gen_add_parallel_loop("ind", str(36 * NJ + nv * nv), use_thread_group)
    self.gen_add_code_line("if (ind < " + str(36 * NJ) + ") { s_temp[" + str(ICOffset) + " + ind] = s_XImats[" + str(36 * NJ) + " + ind]; }")
    self.gen_add_code_line("else { s_M[ind - " + str(36 * NJ) + "] = static_cast<T>(0); }")
    self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)

    for jid in range(NJ - 1, 0, -1):
        parent = self.robot.get_parent_id(jid)
        S_ind = self.robot.get_S_index_by_id(jid)
        S_sign = self.robot.get_S_sign_by_id(jid)
        dof = jid + 5
        self.gen_add_code_line("// CRBA body " + str(jid))
        self.gen_add_code_line(f"grid_linalg_gemm<T,6,6,6,false,true>(&s_XImats[{36*jid}], &s_temp[{ICOffset + 36*jid}], &s_temp[{alphaOffset}], static_cast<T>(1), static_cast<T>(0), s_linalg_smem);")
        self.gen_add_code_line(f"grid_linalg_gemm<T,6,6,6>(&s_temp[{alphaOffset}], &s_XImats[{36*jid}], &s_temp[{ICOffset + 36*parent}], static_cast<T>(1), static_cast<T>(1), s_linalg_smem);")

        self.gen_add_parallel_loop("row", "6", use_thread_group)
        self.gen_add_code_line("s_temp[" + str(fhOffset) + " + row] = static_cast<T>(" + str(S_sign) + ") * s_temp[" + str(ICOffset + 36 * jid + 6 * S_ind) + " + row];")
        self.gen_add_end_control_flow()
        self.gen_add_sync(use_thread_group)
        self.gen_add_serial_ops(use_thread_group)
        self.gen_add_code_line("s_M[" + str(dof + nv * dof) + "] = static_cast<T>(" + str(S_sign) + ") * s_temp[" + str(fhOffset + S_ind) + "];")
        self.gen_add_end_control_flow()
        self.gen_add_sync(use_thread_group)

        chain = []
        curr = jid
        while self.robot.get_parent_id(curr) > 0:
            chain.append((curr, self.robot.get_parent_id(curr)))
            curr = self.robot.get_parent_id(curr)
        chain.append((curr, 0))
        for x_jid, ancestor in chain:
            self.gen_add_parallel_loop("row", "6", use_thread_group)
            self.gen_add_code_line("s_temp[" + str(alphaOffset) + " + row] = dot_prod<T,6,1,1>(&s_XImats[" + str(36 * x_jid) + " + 6*row], &s_temp[" + str(fhOffset) + "]);")
            self.gen_add_end_control_flow()
            self.gen_add_sync(use_thread_group)
            self.gen_add_parallel_loop("row", "6", use_thread_group)
            self.gen_add_code_line("s_temp[" + str(fhOffset) + " + row] = s_temp[" + str(alphaOffset) + " + row];")
            self.gen_add_end_control_flow()
            self.gen_add_sync(use_thread_group)
            if ancestor > 0:
                a_S_ind = self.robot.get_S_index_by_id(ancestor)
                a_S_sign = self.robot.get_S_sign_by_id(ancestor)
                a_dof = ancestor + 5
                self.gen_add_serial_ops(use_thread_group)
                self.gen_add_code_line("s_M[" + str(dof + nv * a_dof) + "] = static_cast<T>(" + str(a_S_sign) + ") * s_temp[" + str(fhOffset + a_S_ind) + "];")
                self.gen_add_code_line("s_M[" + str(a_dof + nv * dof) + "] = s_M[" + str(dof + nv * a_dof) + "];")
                self.gen_add_end_control_flow()
                self.gen_add_sync(use_thread_group)
            else:
                self.gen_add_parallel_loop("col", "6", use_thread_group)
                self.gen_add_code_line("int S_col = col < 3 ? col + 3 : col - 3;")
                self.gen_add_code_line("s_M[" + str(dof) + " + " + str(nv) + "*col] = s_temp[" + str(fhOffset) + " + S_col];")
                self.gen_add_code_line("s_M[col + " + str(nv * dof) + "] = s_M[" + str(dof) + " + " + str(nv) + "*col];")
                self.gen_add_end_control_flow()
                self.gen_add_sync(use_thread_group)

    self.gen_add_code_line("// floating-base root block H[:6,:6] = S^T * IC[0] * S")
    self.gen_add_parallel_loop("ind", "36", use_thread_group)
    self.gen_add_code_line("int row = ind % 6; int col = ind / 6;")
    self.gen_add_code_line("int S_row = row < 3 ? row + 3 : row - 3;")
    self.gen_add_code_line("int S_col = col < 3 ? col + 3 : col - 3;")
    self.gen_add_code_line("s_M[row + " + str(nv) + "*col] = s_temp[" + str(ICOffset) + " + S_row + 6*S_col];")
    self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)

    self.gen_add_end_function()




def gen_crba_device_temp_mem_size(self):
    n = self.robot.get_num_joints()
    wrapper_size = self.gen_topology_helpers_size() + 72*n # for XImats
    return self.gen_crba_inner_temp_mem_size() + wrapper_size

def gen_crba_device(self, use_thread_group = False):
    n = self.robot.get_num_joints()

    # construct the boilerplate and function definition
    func_params = ["s_M is a pointer to the matrix of inertia", \
                   "s_q is the vector of joint positions", \
                   "s_qd is the vector of joint velocities", \
                   "d_robotModel is the pointer to the initialized model specific helpers on the GPU (XImats, topology_helpers, etc.)", \
                   "gravity is the gravity constant"]
    func_notes = []
    func_def_start = "void crba_device("
    func_def_middle = "T *s_M, const T *s_q, const T *s_qd,"
    func_def_end = "const robotModel<T> *d_robotModel, const T gravity) {"
    if use_thread_group:
        func_def_start += "cgrps::thread_group tgrp, "
        func_params.insert(0,"tgrp is the handle to the thread_group running this function")

    func_def = func_def_start + func_def_middle + func_def_end

    # then generate the code
    self.gen_add_func_doc("Compute the CRBA (Composite Rigid Body Algorithm)",\
                          func_notes,func_params,None)
    self.gen_add_code_line("template <typename T>")
    self.gen_add_code_line("__device__")
    self.gen_add_code_line(func_def, True)

    # add the shared memory variables
    shared_mem_size = self.gen_crba_device_temp_mem_size()
    self.gen_XImats_helpers_temp_shared_memory_code(shared_mem_size, include_linalg_scratch=True)

    # then load/update XI and run the algo
    self.gen_load_update_XImats_helpers_function_call(use_thread_group)
    self.gen_crba_inner_function_call(use_thread_group)
    self.gen_add_end_function()

def gen_crba_kernel(self, use_thread_group = False, single_call_timing = False):
    nq = self.robot.get_num_pos()
    nv = self.robot.get_num_vel()
    n = self.robot.get_num_joints()
    input_count = nq + nv
    # define function def and params
    func_params = ["d_M is the pointer to the matrix of inertia", \
                    "d_q_qd is the vector of joint positions and velocities", \
                    "stride_q_qd is the stride between each q, qd", \
                    "d_robotModel is the pointer to the initialized model specific helpers on the GPU (XImats, topology_helpers, etc.)", \
                    "gravity is the gravity constant", \
                    "num_timesteps is the length of the trajectory points we need to compute over (or overloaded as test_iters for timing)"]
    func_notes = []
    func_def_start = "void crba_kernel(T *d_M, const T *d_q_qd, const int stride_q_qd, "
    func_def_end = "const robotModel<T> *d_robotModel, const T gravity, const int NUM_TIMESTEPS) {"
    func_def = func_def_start + func_def_end
    if single_call_timing:
        func_def = func_def.replace("kernel(", "kernel_single_timing(")
    
    # then generate the code
    self.gen_add_func_doc("Compute the CRBA (Composite Rigid Body Algorithm)", \
                            func_notes, func_params, None)
    self.gen_add_code_line("template <typename T, int RESOURCE_TIER = GRID_DEFAULT_RESOURCE_TIER>")
    self.gen_add_code_line("__global__")
    self.gen_add_code_line("__launch_bounds__(tier_max_threads<RESOURCE_TIER>())")
    self.gen_add_code_line(func_def, True)

    # add shared memory variables
    shared_mem_size = self.gen_crba_inner_temp_mem_size()
    self.gen_XImats_helpers_temp_shared_memory_code(shared_mem_size, extra_t_buffers = [("s_M", nv*nv), ("s_q_qd", input_count)], include_linalg_scratch=True)
    self.gen_add_code_line("T *s_q = s_q_qd; T *s_qd = &s_q_qd[" + str(nq) + "];")
    if use_thread_group:
        self.gen_add_code_line("cgrps::thread_group tgrp = TBD;")
    if not single_call_timing:
        # load to shared mem and loop over blocks to compute all requested comps
        self.gen_add_parallel_loop("k","NUM_TIMESTEPS",use_thread_group,block_level = True)
        self.gen_kernel_load_inputs("q_qd","stride_q_qd",str(input_count),use_thread_group)
        # compute
        self.gen_add_code_line("// compute")
        self.gen_load_update_XImats_helpers_function_call(use_thread_group)
        self.gen_crba_inner_function_call(use_thread_group)
        self.gen_add_sync(use_thread_group)
        # save to global  (stride = nv*nv per timestep — without this, batches
        # overlap since each block writes nv*nv elements starting at offset k*1)
        self.gen_kernel_save_result("M",str(nv*nv),str(nv*nv),use_thread_group)
        self.gen_add_end_control_flow()
    else:
        # repurpose NUM_TIMESTEPS for number of timing reps
        self.gen_kernel_load_inputs_single_timing("q_qd",str(input_count),use_thread_group)
        # then compute in loop for timing
        self.gen_add_code_line("// compute with NUM_TIMESTEPS as NUM_REPS for timing")
        self.gen_add_code_line("for (int rep = 0; rep < NUM_TIMESTEPS; rep++){", True)
        self.gen_anti_licm_input_reload("q_qd",str(input_count),use_thread_group,feedback_from="M")
        self.gen_load_update_XImats_helpers_function_call(use_thread_group)
        self.gen_crba_inner_function_call(use_thread_group)
        self.gen_anti_licm_output_write("M")
        self.gen_add_end_control_flow()
        # save to global
        self.gen_kernel_save_result_single_timing("M",str(nv*nv),use_thread_group)
    self.gen_add_end_function()

def gen_crba_host(self, mode = 0):


    #old version that works for iiwa but not for hyq
    single_call_timing = True if mode == 1 else False
    compute_only = True if mode == 2 else False
    # define function def and params
    func_params = ["hd_data is the packaged input and output pointers", \
                   "d_robotModel is the pointer to the initialized model specific helpers on the GPU (XImats, topology_helpers, etc.)", \
                   "gravity is the gravity constant,", \
                   "num_timesteps is the length of the trajectory points we need to compute over (or overloaded as test_iters for timing)", \
                   "streams are pointers to CUDA streams for async memory transfers (if needed)"]
    func_notes = []
    func_def_start = "void crba(gridData<T, KIND> *hd_data, const robotModel<T> *d_robotModel, const T gravity, const int num_timesteps,"
    func_def_end =   "                      const dim3 block_dimms, const dim3 thread_dimms, cudaStream_t *streams) {"
    if single_call_timing:
        func_def_start = func_def_start.replace("(", "_single_timing(")
        func_def_end = "              " + func_def_end
    if compute_only:
        func_def_start = func_def_start.replace("(", "_compute_only(")
        func_def_end = "             " + func_def_end.replace(", cudaStream_t *streams", "")
    # then generate the code
    self.gen_add_func_doc("Compute the CRBA (Composite Rigid Body Algorithm)",\
                          func_notes,func_params,None)
    self.gen_add_code_line("template <typename T, bool USE_COMPRESSED_MEM = false, gridDataKind KIND = GRID_DATA_ALL>")
    self.gen_add_code_line("__host__")
    self.gen_add_code_line(func_def_start)
    self.gen_add_code_line(func_def_end, True)
    self.gen_add_code_line("static_assert(KIND == GRID_DATA_ALL || KIND == GRID_DATA_DYNAMICS, \"crba requires all-data or dynamics gridData\");")
    func_call_start = "crba_kernel<T><<<block_dimms,thread_dimms,CRBA_DYNAMIC_SHARED_MEM_BYTES<T>()>>>(hd_data->d_M,hd_data->d_q_qd,stride_q_qd,"
    func_call_end = "d_robotModel,gravity,num_timesteps);"
    if single_call_timing:
        func_call_start = func_call_start.replace("kernel<T>","kernel_single_timing<T>")
    if not compute_only:
        # start code with memory transfer
        self.gen_add_code_lines(["// start code with memory transfer", \
                                 "int stride_q_qd;", \
                                 "if (USE_COMPRESSED_MEM) {stride_q_qd = NUM_JOINTS + NUM_VEL; " + \
                                    "gpuErrchk(cudaMemcpyAsync(hd_data->d_q_qd,hd_data->h_q_qd,stride_q_qd*" + \
                                    ("num_timesteps*" if not single_call_timing else "") + "sizeof(T),cudaMemcpyHostToDevice,streams[0]));}", \
                                 "else {stride_q_qd = NUM_JOINTS + 2*NUM_VEL; " + \
                                    "gpuErrchk(cudaMemcpyAsync(hd_data->d_q_qd_u,hd_data->h_q_qd_u,stride_q_qd*" + \
                                    ("num_timesteps*" if not single_call_timing else "") + "sizeof(T),cudaMemcpyHostToDevice,streams[0]));}"])
    else:
        self.gen_add_code_line("int stride_q_qd = USE_COMPRESSED_MEM ? NUM_JOINTS + NUM_VEL : NUM_JOINTS + 2*NUM_VEL;")
    self.gen_add_code_line("// then call the kernel")
    func_call = func_call_start + func_call_end
    # add in compressed mem adjusts
    func_call_mem_adjust = "if (USE_COMPRESSED_MEM) {" + func_call + "}"
    func_call_mem_adjust2 = "else                    {" + func_call.replace("hd_data->d_q_qd","hd_data->d_q_qd_u") + "}"
    # compule into a set of code
    func_call_code = [func_call_mem_adjust, func_call_mem_adjust2, "gpuErrchkKernel();"]
    # wrap function call in timing (if needed)
    if single_call_timing:
        func_call_code.insert(0,"struct timespec start, end; clock_gettime(CLOCK_MONOTONIC,&start);")
        func_call_code.append("clock_gettime(CLOCK_MONOTONIC,&end);")
    self.gen_add_code_line("gpuErrchk(grid_check_dynamic_shared_memory_bytes(\"crba\", CRBA_DYNAMIC_SHARED_MEM_BYTES<T>()));")
    self.gen_add_code_lines(func_call_code)
    if not compute_only:
        # then transfer memory back
        self.gen_add_code_lines(["// finally transfer the result back", \
                                 "gpuErrchk(cudaMemcpy(hd_data->h_M,hd_data->d_M,NUM_VEL*NUM_VEL*" + \
                                    ("num_timesteps*" if not single_call_timing else "") + "sizeof(T),cudaMemcpyDeviceToHost));",
                                 "gpuErrchkKernel();"])
    # finally report out timing if requested
    if single_call_timing:
        from ..algo_registry import single_call_printf_line
        self.gen_add_code_line(single_call_printf_line("crba"))
    self.gen_add_end_function()

def gen_crba(self, use_thread_group = False):
    # first generate the inner helpers
    self.gen_crba_inner(use_thread_group)
    # then generate the device wrappers
    self.gen_crba_device(use_thread_group)
    # then generate the kernels
    self.gen_crba_kernel(use_thread_group,True)
    self.gen_crba_kernel(use_thread_group,False)
    # then the host launch wrappers
    self.gen_crba_host(0)
    self.gen_crba_host(1)
    self.gen_crba_host(2)
