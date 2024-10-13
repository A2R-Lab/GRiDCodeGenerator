def gen_aba_inner(self, use_thread_group = False): 
    n = self.robot.get_num_pos()
    n_bfs_levels = self.robot.get_max_bfs_level() + 1 # starts at 0
	# construct the boilerplate and function definition
    func_params = ["s_qdd is the vector of joint accelerations", \
                "s_va is a pointer to shared memory of size 2*6*NUM_JOINTS = " + str(12*n), \
                "s_q is the vector of joint positions", \
                "s_qd is the vector of joint velocities", \
                "s_tau is the vector of joint torques", \
                "s_temp is the pointer to the shared memory needed of size: " + \
                            str(self.gen_forward_dynamics_inner_temp_mem_size()), \
                "gravity is the gravity constant"]
    func_def_start = "void aba_inner("
    func_def_middle = "T *s_qdd, T *s_va, const T *s_q, const T *s_qd, const T *s_tau, "
    func_def_end = "T *s_temp, const T gravity) {"
    func_notes = ["Assumes the XI matricies have already been updated for the given q"]
    if use_thread_group:
        func_def_start = func_def_start.replace("(", "(cgrps::thread_group tgrp, ")
        func_params.insert(0,"tgrp is the handle to the thread_group running this function")
    func_def_middle, func_params = self.gen_insert_helpers_func_def_params(func_def_middle, func_params, -2)
    func_def = func_def_start + func_def_middle + func_def_end
    self.gen_add_func_doc("Computes the Articulated Body Algorithm", func_notes, func_params, None)
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

    #
    # Forward Pass we are going to go in bfs_level waves
    # 
    self.gen_add_code_line("//")
    self.gen_add_code_line("// Forward Pass")
    self.gen_add_code_line("//")
    for bfs_level in range(n_bfs_levels):
        inds = self.robot.get_ids_by_bfs_level(bfs_level)
        joint_names = [self.robot.get_joint_by_id(ind).get_name() for ind in inds]
        link_names = [self.robot.get_link_by_id(ind).get_name() for ind in inds]
        parent_ind_cpp, S_ind_cpp = self.gen_topology_helpers_pointers_for_cpp(inds, NO_GRAD_FLAG = True)

        if bfs_level == 0:
            self.gen_add_code_line("// s_v where parent is base")
            self.gen_add_code_line("//     joints are: " + ", ".join(joint_names))
            self.gen_add_code_line("//     links are: " + ", ".join(link_names))
            # compute the initial v which is just S*qd
            self.gen_add_code_line("// s_v[k] = S[k]*qd[k]")
            if len(inds) > 1:
                self.gen_add_parallel_loop("ind",str(6*len(inds)),use_thread_group)
                self.gen_add_code_line("int row = ind % 6;")
                select_var_vals = [("int", "jid", [str(jid) for jid in inds])]
                self.gen_add_multi_threaded_select("ind", "<", [str(6*(i+1)) for i in range(len(inds))], select_var_vals)
                jid = "jid"
            else:
                self.gen_add_parallel_loop("row",str(6),use_thread_group)
                jid = str(inds[0])
                self.gen_add_code_line("int jid = " + jid + ";")
            # load in 0 to v 
            self.gen_add_code_lines(["int jid6 = 6*jid;", \
                                     "s_va[jid6 + row] = static_cast<T>(0);"])
            # add in qd
            self.gen_add_code_line("if (row == " + S_ind_cpp + "){s_va[jid6 + " + S_ind_cpp + "] += s_qd[" + jid + "];}")
            self.gen_add_end_control_flow()
            self.gen_add_sync(use_thread_group)

            # add debug if requested
            if self.DEBUG_MODE:
                self.gen_add_sync(use_thread_group)
                self.gen_add_serial_ops(use_thread_group)
                for ind in inds:
                    self.gen_add_code_line("printf(\"s_v[" + str(ind) + "]\\n\"); printMat<T,1,6>(&s_va[6*" + str(ind) + "],1);")
                self.gen_add_end_control_flow()
                self.gen_add_sync(use_thread_group)
        
        else:
            self.gen_add_code_line("// s_v where bfs_level is " + str(bfs_level))
            self.gen_add_code_line("//     joints are: " + ", ".join(joint_names))
            self.gen_add_code_line("//     links are: " + ", ".join(link_names))
            self.gen_add_code_line("// s_v[k] = X[k]*v[parent_k] + S[k]*qd[k]")

            self.gen_add_parallel_loop("ind",str(6*len(inds)),use_thread_group)
            if len(inds) > 1:
                self.gen_add_code_line("int row = ind % 6; int comp = ind / 6; int comp_mod = comp % " + str(len(inds)) + ";")
                select_var_vals = [("int", "jid", [str(jid) for jid in inds])]
                jid = "jid"
                self.gen_add_multi_threaded_select("comp_mod", "==", [str(i) for i in range(len(inds))], select_var_vals)
            else:
                self.gen_add_code_line("int row = ind % 6;")
                jid = str(inds[0])
                self.gen_add_code_line("int jid = " + jid + ";")
            self.gen_add_code_line("int jid6 = 6 * jid;")
            self.gen_add_code_line("T qd_val = (row == " + S_ind_cpp + ") * (s_qd[" + jid + "]);")
            self.gen_add_code_line("s_va[jid6 + row] = dot_prod<T,6,6,1>(&s_XImats[6*jid6 + row], &s_va[6*" + parent_ind_cpp + "]) + qd_val;")
            

            self.gen_add_end_control_flow()
            self.gen_add_sync(use_thread_group)            
        
            # add debug if requested
            if self.DEBUG_MODE:
                self.gen_add_sync(use_thread_group)
                self.gen_add_serial_ops(use_thread_group)
                for ind in inds:
                    self.gen_add_code_line("printf(\"s_v[" + str(ind) + "] = X*s_v[" + parent_ind_cpp + "] + S*qd[" + str(ind) + "]\\n\"); printMat<T,1,6>(&s_va[6*" + str(ind) + "],1);")
                self.gen_add_end_control_flow()
                self.gen_add_sync(use_thread_group)

    # calculate c
    self.gen_add_code_line("// c[k] = mxS(v[k])*qd[k]")
    self.gen_add_parallel_loop("ind", str(n), use_thread_group)
    self.gen_add_code_line("int jid = ind;")
    self.gen_add_code_line("int jid6 = 6 * jid;")
    self.gen_add_code_line("mx2_scaled<T>(&s_temp[72 * " + str(n) + "+jid6], &s_va[jid6], s_qd[jid]);")
    self.gen_add_end_control_flow()
    
    # add debug if requested
    if self.DEBUG_MODE:
        self.gen_add_sync(use_thread_group)
        self.gen_add_serial_ops(use_thread_group)
        self.gen_add_code_line("printf(\"c\\n\"); printMat<T,6,"+str(n)+">(&s_temp[72 * "+str(n)+"], 6);")
        self.gen_add_end_control_flow()
        self.gen_add_sync(use_thread_group)

    # set IA = I
    self.gen_add_code_line("// Initialize IA = I")
    self.gen_add_parallel_loop("ind",str(36*n),use_thread_group)
    self.gen_add_code_line("s_temp[ind] = s_XImats[" + str(36*n) + " + ind];")
    self.gen_add_end_control_flow()
    
    # initialize vcross from v
    self.gen_add_code_line("// Initialize vcross[k]")
    self.gen_add_parallel_loop("ind", str(n), use_thread_group)
    self.gen_add_code_line("int jid = ind;")
    self.gen_add_code_line("int jid6 = 6 * jid;")

    self.gen_add_code_line("vcross<T>(&s_temp[36*("+str(n)+"+jid)], &s_va[jid6]);")
    self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)

    self.gen_add_code_line("// temp[k] = -vcross.T*I[k]")
    self.gen_add_parallel_loop("ind", str(36*n), use_thread_group)
    self.gen_add_code_line("int row = ind % 6; int col = (ind / 6) %6; int jid = ind / 36;")
    self.gen_add_code_line("int jid6 = 6 * jid;")
    self.gen_add_code_line("s_temp[98 * " + str(n) + " + jid6*6 + row+col*6] = -1 * dot_prod<T,6,1,1>(&s_temp[36*("+str(n)+"+jid)+row*6], &s_XImats[36 * ("+str(n)+"+jid) + col*6]);")
    self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)
    # calculate pA
    self.gen_add_code_line("// pA[k] = temp[k]*v[k][0]")
    self.gen_add_parallel_loop("ind", str(6*n), use_thread_group)
    self.gen_add_code_line("int row = ind % 6; int comp = ind / 6; int jid = comp % " + str(n) + ";")
    self.gen_add_code_line("int jid6 = 6 * jid;")
    self.gen_add_code_line("s_temp[78 * " + str(n) + " + jid6 + row] = dot_prod<T,6,6,1>(&s_temp[98 * " + str(n) + " + 6*jid6+row], &s_va[jid6]);")

    self.gen_add_end_control_flow()

    if self.DEBUG_MODE:
        self.gen_add_sync(use_thread_group)
        self.gen_add_serial_ops(use_thread_group)
        self.gen_add_code_line("for (int i = 0; i < " + str(n) + "; i++){printf(\"IA[%d]\\n\",i); printMat<T,6,6>(&s_temp[36*(i)],6);}")
        self.gen_add_code_line("printf(\"pA\\n\"); printMat<T,6,"+str(n)+">(&s_temp[78 * "+str(n)+"], 6);")
        self.gen_add_end_control_flow()
        self.gen_add_sync(use_thread_group)

    #
    # Then compute the Backward Pass again in bfs waves
    #
    self.gen_add_code_line("//")
    self.gen_add_code_line("// Backward Pass")
    self.gen_add_code_line("//")
    for bfs_level in range(n_bfs_levels - 1, -1, -1): 
        inds = self.robot.get_ids_by_bfs_level(bfs_level)
        joint_names = [self.robot.get_joint_by_id(ind).get_name() for ind in inds]
        link_names = [self.robot.get_link_by_id(ind).get_name() for ind in inds]
        parent_ind_cpp, S_ind_cpp = self.gen_topology_helpers_pointers_for_cpp(inds, NO_GRAD_FLAG = True)
        self.gen_add_code_line("// Backward pass where bfs_level is " + str(bfs_level))
        self.gen_add_code_line("//     joints are: " + ", ".join(joint_names))
        self.gen_add_code_line("//     links are: " + ", ".join(link_names))
        # caclulate U, which is just IA*S
        self.gen_add_code_line("// U[k] = IA[k]*S[k]")
        self.gen_add_parallel_loop("ind", str(6*len(inds)), use_thread_group)
        self.gen_add_code_line("int row = ind % 6;")
        if len(inds) > 1:
            select_var_vals = [("int", "jid", [str(jid) for jid in inds])]
            self.gen_add_multi_threaded_select("ind", "<", [str(6*(i+1)) for i in range(len(inds))], select_var_vals)
            jid = "jid"
        else:
            jid = str(inds[0])
            self.gen_add_code_line("int jid = " + jid + ";")
        self.gen_add_code_line("int jid6 = 6 * " + jid + ";")

        self.gen_add_code_line("s_temp[84*"+str(n)+"+jid6+row] = s_temp[36*jid+row+6*("+ S_ind_cpp+")];")
        self.gen_add_end_control_flow()
        self.gen_add_sync(use_thread_group)

        # caclulate d which is S*U and u which is tau - S*pA
        self.gen_add_code_line("// d[k] = S[k]*U[k], u[k] = tau[k] - S[k].T*pA[k]")
        self.gen_add_parallel_loop("ind", str(len(inds)), use_thread_group)
        if len(inds) > 1:
            select_var_vals = [("int", "jid", [str(jid) for jid in inds])]
            self.gen_add_multi_threaded_select("ind", "<", [str((i+1)) for i in range(len(inds))], select_var_vals)
            jid = "jid"
        else:
            jid = str(inds[0])
            self.gen_add_code_line("int jid = " + jid + ";")
        self.gen_add_code_line("int jid6 = 6 * " + jid + ";")
        self.gen_add_code_line("s_temp[96 * "+ str(n) +" + jid] = s_temp[84 * " + str(n) + " + jid6 + " + S_ind_cpp + "];")
        
        self.gen_add_code_line("T tempval = s_temp[78 * " + str(n) + " + jid6 + " + S_ind_cpp +"];") 
        self.gen_add_code_line("s_temp[97 * " + str(n) + " + jid] = s_tau[jid] - tempval;")
        self.gen_add_end_control_flow()
        self.gen_add_sync(use_thread_group)
        
        # calculate Ia from IA, U, and d
        self.gen_add_code_line("// Ia[k] = IA[k] - U[k]*U[k].T/d[k]")
        self.gen_add_parallel_loop("ind", str(36 * len(inds)), use_thread_group)
        self.gen_add_code_line("int row = ind % 6; int col = (ind / 6) %6;")
        if len(inds) > 1:
            select_var_vals = [("int", "jid", [str(jid) for jid in inds])]
            self.gen_add_multi_threaded_select("ind", "<", [str(36*(i+1)) for i in range(len(inds))], select_var_vals)
            jid = "jid"
        else:
            jid = str(inds[0])
            self.gen_add_code_line("int jid = " + jid + ";")
        self.gen_add_code_line("int jid6 = 6 * " + jid + ";")

        self.gen_add_code_line("s_temp[36 * "+str(n)+"+6*jid6+row+6*col] = s_temp[84*"+str(n)+"+jid6+row]*s_temp[84*"+str(n)+"+jid6+col]/s_temp[96 *"+str(n)+"+jid];")

        self.gen_add_code_line("s_temp[36 * "+str(n)+"+6*jid6+row+6*col] = s_temp[6*jid6+row+6*col] - s_temp[36 * "+str(n)+"+6*jid6+row+6*col];")
        self.gen_add_end_control_flow()
        self.gen_add_sync(use_thread_group)

        # caclulate pa
        self.gen_add_code_line("// pa[k] = pA[k] + Ia[k]*c[k]+U[k]*u[k]/d[k]")
        self.gen_add_parallel_loop("ind", str(6*len(inds)), use_thread_group)
        self.gen_add_code_line("int row = ind % 6;")
        if len(inds) > 1:
            select_var_vals = [("int", "jid", [str(jid) for jid in inds])]
            self.gen_add_multi_threaded_select("ind", "<", [str(6*(i+1)) for i in range(len(inds))], select_var_vals)
            jid = "jid"
        else:
            jid = str(inds[0])
            self.gen_add_code_line("int jid = " + jid + ";")
        self.gen_add_code_line("int jid6 = 6 * " + jid + ";")

        self.gen_add_code_line("T Uval = s_temp[84 * "+str(n)+"+jid6+row]*s_temp[97*"+str(n)+"+jid]/s_temp[96*"+str(n)+"+jid];")
        self.gen_add_code_line("s_temp[90 * "+str(n)+" + jid6 + row] = s_temp[78 * "+str(n)+" + jid6+row] + dot_prod<T,6,6,1>(&s_temp[36*("+str(n)+"+jid)+row], &s_temp[72*"+str(n)+"+jid6]) + Uval;")
        self.gen_add_end_control_flow()
        
        if bfs_level != 0:
            self.gen_add_code_line("// temp[k] = X[k].T*Ia[k]")
            self.gen_add_parallel_loop("ind", str(36 * len(inds)), use_thread_group)
            self.gen_add_code_line("int row = ind % 6; int col = (ind / 6) %6;")
            if len(inds) > 1:
                select_var_vals = [("int", "jid", [str(jid) for jid in inds])]
                self.gen_add_multi_threaded_select("ind", "<", [str(36*(i+1)) for i in range(len(inds))], select_var_vals)
                jid = "jid"
            else:
                jid = str(inds[0])
                self.gen_add_code_line("int jid = " + jid + ";")
            self.gen_add_code_line("int jid6 = 6 * jid;")

            self.gen_add_code_line("s_temp[98 * " + str(n) + " + 6 * jid6 + row + 6*col] = dot_prod<T,6,1,1>(&s_XImats[6*jid6+6*row], &s_temp[36 * "+str(n)+"+jid6*6+6*col]);")
            self.gen_add_end_control_flow()
            self.gen_add_sync(use_thread_group)
            
            # update IA of the parent
            self.gen_add_code_line("// IA[parent] += temp[k]*X[k]")
            self.gen_add_parallel_loop("ind", str(36 * len(inds)), use_thread_group)
            self.gen_add_code_line("int row = ind % 6; int col = (ind / 6) %6;")
            if len(inds) > 1:
                select_var_vals = [("int", "jid", [str(jid) for jid in inds])]
                self.gen_add_multi_threaded_select("ind", "<", [str(36*(i+1)) for i in range(len(inds))], select_var_vals)
                jid = "jid"
            else:
                jid = str(inds[0])
                self.gen_add_code_line("int jid = " + jid + ";")
            self.gen_add_code_line("int jid6 = 6 * jid;")
            self.gen_add_code_line("T prodtemp = static_cast<T>(0);")
            self.gen_add_code_line("prodtemp =  dot_prod<T,6,6,1>(&s_temp[98 * " + str(n) + " + 6 * jid6 + row], &s_XImats[6*jid6+6*col]);")
            self.gen_add_code_line("atomicAdd(&s_temp[36 * " + parent_ind_cpp +" + row + 6*col], prodtemp);")
            self.gen_add_end_control_flow()
            self.gen_add_sync(use_thread_group)

            # update pA of the parent
            self.gen_add_code_line("// pA[parent] += X[k].T*pa[k]")
            self.gen_add_parallel_loop("ind", str(6*len(inds)), use_thread_group)
            self.gen_add_code_line("int row = ind % 6;")
            if len(inds) > 1:
                select_var_vals = [("int", "jid", [str(jid) for jid in inds])]
                self.gen_add_multi_threaded_select("ind", "<", [str(6*(i+1)) for i in range(len(inds))], select_var_vals)
                jid = "jid"
            else:
                jid = str(inds[0])
                self.gen_add_code_line("int jid = " + jid + ";")
            self.gen_add_code_line("int jid6 = 6 * " + jid + ";")
            self.gen_add_code_line("s_temp[134 * " + str(n) + " + jid6 + row] = dot_prod<T,6,1,1>(&s_XImats[36*jid+6*row],&s_temp[90*"+str(n)+"+jid6]);")
            self.gen_add_code_line("atomicAdd(&s_temp[78 * " + str(n) + " + 6 * "+ parent_ind_cpp +" + row], s_temp[134 * " + str(n) + " + jid6 + row]);")
            self.gen_add_end_control_flow()
            self.gen_add_sync(use_thread_group)

    # add debug if requested
    if self.DEBUG_MODE:
        self.gen_add_sync(use_thread_group)
        self.gen_add_serial_ops(use_thread_group)
        self.gen_add_code_line("printf(\"U \\n\"); printMat<T,6,"+str(n)+">(&s_temp[84 * "+str(n)+"], 6);")
        self.gen_add_code_line("printf(\"d \\n\"); printMat<T,1,"+str(n)+">(&s_temp[96 * "+str(n)+"], 1);")
        self.gen_add_code_line("printf(\"u \\n\"); printMat<T,1,"+str(n)+">(&s_temp[97 * "+str(n)+"], 1);")
        self.gen_add_code_line("for (int i = 0; i < " + str(n) + "; i++){printf(\"Ia[%d]\\n\",i); printMat<T,6,6>(&s_temp[36*("+str(n)+"+i)],6);}")
        self.gen_add_code_line("for (int i = 0; i < " + str(n) + "; i++){printf(\"IA[%d]\\n\",i); printMat<T,6,6>(&s_temp[36*(i)],6);}")
        self.gen_add_code_line("printf(\"pA\\n\"); printMat<T,6,"+str(n)+">(&s_temp[78 * "+str(n)+"], 6);")
        self.gen_add_end_control_flow()
        self.gen_add_sync(use_thread_group)

    self.gen_add_code_line("//")
    self.gen_add_code_line("// Second Forward Pass")
    self.gen_add_code_line("//")
    for bfs_level in range(n_bfs_levels):
        inds = self.robot.get_ids_by_bfs_level(bfs_level)
        parent_ind_cpp, S_ind_cpp = self.gen_topology_helpers_pointers_for_cpp(inds, NO_GRAD_FLAG = True)
        joint_names = [self.robot.get_joint_by_id(ind).get_name() for ind in inds]
        link_names = [self.robot.get_link_by_id(ind).get_name() for ind in inds]
        # calculate a where parent is base
        if bfs_level == 0:
            self.gen_add_code_line("// s_a, qdd where parent is base")
            self.gen_add_code_line("//     joints are: " + ", ".join(joint_names))
            self.gen_add_code_line("//     links are: " + ", ".join(link_names))
            self.gen_add_code_line("// a[k] = X[k]*gravity_vec + c[k]")
            if len(inds) > 1:
                self.gen_add_parallel_loop("ind",str(6*len(inds)),use_thread_group)
                self.gen_add_code_line("int row = ind % 6;")
                select_var_vals = [("int", "jid", [str(jid) for jid in inds])]
                self.gen_add_multi_threaded_select("ind", "<", [str(6*(i+1)) for i in range(len(inds))], select_var_vals)
                jid = "jid"
            else:
                self.gen_add_parallel_loop("row",str(6),use_thread_group)
                jid = str(inds[0])
                self.gen_add_code_line("int jid = " + jid + ";")
            self.gen_add_code_line("int jid6 = 6*" + jid + ";")
            self.gen_add_code_line("T gravity_vec[] = {0,0,0,0,0,gravity};")
            self.gen_add_code_line("s_va[6*"+str(n)+"+jid6+row] = dot_prod<T,6,6,1>(&s_XImats[36 * jid + row], &gravity_vec[0]) + s_temp[72*"+str(n)+"+jid6+row];")
        # calculate a where parent is not base
        else:
            self.gen_add_code_line("// s_a, s_qdd where bfs_level is " + str(bfs_level))
            self.gen_add_code_line("//     joints are: " + ", ".join(joint_names))
            self.gen_add_code_line("//     links are: " + ", ".join(link_names))
            self.gen_add_code_line("// a[k] = X[k]*a[parent] + c[k]")
            self.gen_add_parallel_loop("ind",str(6*len(inds)),use_thread_group)
            if len(inds) > 1:
                self.gen_add_code_line("int row = ind % 6; int comp = ind / 6; int comp_mod = comp % " + str(len(inds)) + ";")
                select_var_vals = [("int", "jid", [str(jid) for jid in inds])]
                jid = "jid"
                self.gen_add_multi_threaded_select("comp_mod", "==", [str(i) for i in range(len(inds))], select_var_vals)
            else:
                self.gen_add_code_line("int row = ind % 6;")
                jid = str(inds[0])
                self.gen_add_code_line("int jid = " + jid + ";")
            self.gen_add_code_line("int jid6 = 6 * " + jid + ";")
            self.gen_add_code_line("s_va[6*"+str(n)+"+jid6+row] = dot_prod<T,6,6,1>(&s_XImats[36 * jid + row], &s_va[6*"+str(n)+"+(6 * "+ parent_ind_cpp + ")]) + s_temp[72*"+str(n)+"+jid6+row];")
        self.gen_add_end_control_flow()
        self.gen_add_sync(use_thread_group)
        
        # calculate qdd which is (u - U*a)/d 
        self.gen_add_code_line("// qdd[k] = (u[k] - U[k].T*a[k])/d[k]")
        self.gen_add_parallel_loop("ind",str(len(inds)),use_thread_group)
        if len(inds) > 1:
            self.gen_add_code_line("int comp_mod = ind % "+ str(len(inds)) + ";")
            select_var_vals = [("int", "jid", [str(jid) for jid in inds])]
            jid = "jid"
            self.gen_add_multi_threaded_select("comp_mod", "==", [str(i) for i in range(len(inds))], select_var_vals)
        else:
            jid = str(inds[0])
            self.gen_add_code_line("int jid = " + jid + ";")
        self.gen_add_code_line("int jid6 = 6 * " + jid + ";")
        self.gen_add_code_line("T tempval = s_temp[97 * "+str(n)+"+jid] - dot_prod<T,6,1,1>(&s_temp[84*"+str(n)+"+jid6], &s_va[6*"+str(n)+"+jid6]);")
        self.gen_add_code_line("s_qdd[jid] = tempval / s_temp[96*"+str(n)+"+jid];")
        self.gen_add_end_control_flow()
        self.gen_add_sync(use_thread_group)

        # update a by adding qdd*S
        self.gen_add_code_line("// a[k] += qdd[k]*S[k]")
        self.gen_add_parallel_loop("ind",str(6*len(inds)),use_thread_group)
        
        if len(inds) > 1:
            self.gen_add_code_line("int row = ind % 6; int comp = ind / 6; int comp_mod = comp % " + str(len(inds)) + ";")
            select_var_vals = [("int", "jid", [str(jid) for jid in inds])]
            jid = "jid"
            self.gen_add_multi_threaded_select("comp_mod", "==", [str(i) for i in range(len(inds))], select_var_vals)
        else:
            self.gen_add_code_line("int row = ind % 6;")
            jid = str(inds[0])
            self.gen_add_code_line("int jid = " + jid + ";")
        self.gen_add_code_line("int jid6 = 6 * " + jid + ";")
        self.gen_add_code_line("T qdd_val = (row == " + S_ind_cpp + ") * (s_qdd[jid]);")
        self.gen_add_code_line("s_va[6*"+str(n)+"+jid6+row] += qdd_val;")

        self.gen_add_end_control_flow()
        self.gen_add_sync(use_thread_group)
    
    # add debug if requested
    if self.DEBUG_MODE:
        self.gen_add_sync(use_thread_group)
        self.gen_add_serial_ops(use_thread_group)
        self.gen_add_code_line("printf(\"a\\n\"); printMat<T,6,"+str(n)+">(&s_va[6 * "+str(n)+"], 6);")
        self.gen_add_code_line("printf(\"qdd\\n\"); printMat<T,1," + str(n) + ">(s_qdd,1);")
        self.gen_add_end_control_flow()
        self.gen_add_sync(use_thread_group)
        

    self.gen_add_end_function()

def gen_aba_inner_temp_mem_size(self):
    n = self.robot.get_num_pos()
    return 140 * n

def gen_aba_inner_function_call(self, use_thread_group = False, updated_var_names = None):
    var_names = dict( \
        s_va_name = "s_va", \
        s_q_name = "s_q", \
        s_qd_name = "s_qd", \
        s_qdd_name = "s_qdd", \
        s_tau_name = "s_tau", \
        s_temp_name = "s_temp", \
        gravity_name = "gravity"
    )
    if updated_var_names is not None:
        for key,value in updated_var_names.items():
            var_names[key] = value
    aba_code_start = "aba_inner<T>(" + var_names["s_qdd_name"] + ", " + var_names["s_va_name"] + ", " + var_names["s_q_name"] + ", " + var_names["s_qd_name"] + ", " + var_names["s_tau_name"] + ", "
    aba_code_end = var_names["s_temp_name"] + ", " + var_names["gravity_name"] + ");"
    if use_thread_group:
        id_code_start = id_code_start.replace("(","(tgrp, ")
    aba_code_middle = self.gen_insert_helpers_function_call()
    aba_code = aba_code_start + aba_code_middle + aba_code_end
    self.gen_add_code_line(aba_code)

def gen_aba_device(self, use_thread_group = False):
    n = self.robot.get_num_pos()
    # construct the boilerplate and function definition
    func_params = ["s_q is the vector of joint positions", \
                   "s_qd is the vector of joint velocities", \
                    "s_tau is the vector of joint torques", \
                   "d_robotModel is the pointer to the initialized model specific helpers on the GPU (XImats, topology_helpers, etc.)", \
                   "gravity is the gravity constant"]
    func_notes = []
    func_def_start = "void aba_device("
    func_def_middle = "const T *s_q, const T *s_qd, const T *s_tau, "
    func_def_end = "const robotModel<T> *d_robotModel, const T gravity) {"
    if use_thread_group:
        func_def_start += "cgrps::thread_group tgrp, "
        func_params.insert(0,"tgrp is the handle to the thread_group running this function")
    func_def = func_def_start + func_def_middle + func_def_end

    # then generate the code
    self.gen_add_func_doc("Compute the ABA (Articulated Body Algorithm)",\
                          func_notes,func_params,None)
    self.gen_add_code_line("template <typename T>")
    self.gen_add_code_line("__device__")
    self.gen_add_code_line(func_def, True)

    # add the shared memory variables
    shared_mem_size = self.gen_aba_inner_temp_mem_size() if not self.use_dynamic_shared_mem_flag else None
    self.gen_XImats_helpers_temp_shared_memory_code(shared_mem_size)

    self.gen_add_code_line("extern __shared__ T s_va[2*6*"+ str(n) +"];")
    self.gen_add_code_line("extern __shared__ T s_qdd[" + str(n) + "];")
    
    # then load/update XI and run the algo
    self.gen_load_update_XImats_helpers_function_call(use_thread_group)
    self.gen_aba_inner_function_call(use_thread_group)
    self.gen_add_end_function()

def gen_aba_kernel(self, use_thread_group = False, single_call_timing = False):
    n = self.robot.get_num_pos()
    # define function def and params
    func_params = ["d_q_qd_tau is the vector of joint positions and velocities", \
                    "stride_q_qd is the stride between each q, qd", \
                    "d_robotModel is the pointer to the initialized model specific helpers on the GPU (XImats, topology_helpers, etc.)", \
                    "d_tau is the vector of joint torques", \
                    "gravity is the gravity constant", \
                    "num_timesteps is the length of the trajectory points we need to compute over (or overloaded as test_iters for timing)"]
    func_notes = []
    func_def_start = "void aba_kernel(T *d_qdd, const T *d_q_qd_tau, const int stride_q_qd, "
    func_def_end = "const robotModel<T> *d_robotModel, const T gravity, const int NUM_TIMESTEPS) {"
    func_def = func_def_start + func_def_end
    if single_call_timing:
        func_def = func_def.replace("kernel(", "kernel_single_timing(")
    
    # then generate the code
    self.gen_add_func_doc("Compute the ABA (Articulated Body Algorithm)", \
                            func_notes, func_params, None)
    self.gen_add_code_line("template <typename T>")
    self.gen_add_code_line("__global__")
    self.gen_add_code_line(func_def, True)

    # add shared memory variables
    shared_mem_vars = ["__shared__ T s_qdd[" + str(n) + "];", \
                        "__shared__ T s_q_qd_tau[3*" + str(n) + "]; T *s_q = s_q_qd_tau; T *s_qd = &s_q_qd_tau[" + str(n) + "]; T *s_tau = &s_q_qd_tau[2 * " + str(n) + "];", \
                        "__shared__ T s_va[" + str(12*n) + "];"]
    self.gen_add_code_lines(shared_mem_vars)
    shared_mem_size = self.gen_aba_inner_temp_mem_size() if not self.use_dynamic_shared_mem_flag else None
    self.gen_XImats_helpers_temp_shared_memory_code(shared_mem_size)
    if use_thread_group:
        self.gen_add_code_line("cgrps::thread_group tgrp = TBD;")
    if not single_call_timing:
        # load to shared mem and loop over blocks to compute all requested comps
        self.gen_add_parallel_loop("k","NUM_TIMESTEPS",use_thread_group,block_level = True)
        self.gen_kernel_load_inputs("q_qd_tau","stride_q_qd",str(3*n),use_thread_group)
        # compute
        self.gen_add_code_line("// compute")
        self.gen_load_update_XImats_helpers_function_call(use_thread_group)
        self.gen_aba_inner_function_call(use_thread_group)
        self.gen_add_sync(use_thread_group)
        # save to global
        self.gen_kernel_save_result("qdd","1",str(n),use_thread_group)
        self.gen_add_end_control_flow()
    else:
        # repurpose NUM_TIMESTEPS for number of timing reps
        self.gen_kernel_load_inputs_single_timing("q_qd_tau",str(3*n),use_thread_group)
        # then compute in loop for timing
        self.gen_add_code_line("// compute with NUM_TIMESTEPS as NUM_REPS for timing")
        self.gen_add_code_line("for (int rep = 0; rep < NUM_TIMESTEPS; rep++){", True)
        self.gen_load_update_XImats_helpers_function_call(use_thread_group)
        self.gen_aba_inner_function_call(use_thread_group)
        self.gen_add_end_control_flow()
        # save to global
        self.gen_kernel_save_result_single_timing("qdd",str(n),use_thread_group)
    self.gen_add_end_function()

def gen_aba_host(self, mode = 0):
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
    func_def_start = "void aba(gridData<T> *hd_data, const robotModel<T> *d_robotModel, const T gravity, const int num_timesteps,"
    func_def_end =   "                      const dim3 block_dimms, const dim3 thread_dimms, cudaStream_t *streams) {"
    if single_call_timing:
        func_def_start = func_def_start.replace("(", "_single_timing(")
        func_def_end = "              " + func_def_end
    if compute_only:
        func_def_start = func_def_start.replace("(", "_compute_only(")
        func_def_end = "             " + func_def_end.replace(", cudaStream_t *streams", "")
    # then generate the code
    self.gen_add_func_doc("Compute the ABA (Articulated Body Algorithm)",\
                          func_notes,func_params,None)
    self.gen_add_code_line("template <typename T>")
    self.gen_add_code_line("__host__")
    self.gen_add_code_line(func_def_start)
    self.gen_add_code_line(func_def_end, True)

    func_call_start = "aba_kernel<T><<<block_dimms,thread_dimms,ABA_DYNAMIC_SHARED_MEM_COUNT*sizeof(T)>>>(hd_data->d_qdd,hd_data->d_q_qd_u,stride_q_qd,"
    func_call_end = "d_robotModel,gravity,num_timesteps);"
    self.gen_add_code_line("int stride_q_qd = 3*NUM_JOINTS;")
    if single_call_timing:
        func_call_start = func_call_start.replace("kernel<T>","kernel_single_timing<T>")
    if not compute_only:
        # start code with memory transfer
        self.gen_add_code_lines(["// start code with memory transfer", \
                                 "gpuErrchk(cudaMemcpyAsync(hd_data->d_q_qd_u,hd_data->h_q_qd_u,stride_q_qd*" + \
                                    ("num_timesteps*" if not single_call_timing else "") + "sizeof(T),cudaMemcpyHostToDevice,streams[0]));", \
                                 "gpuErrchk(cudaDeviceSynchronize());"])
    # then compute:
    self.gen_add_code_line("// then call the kernel")
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
                                "gpuErrchk(cudaMemcpy(hd_data->h_qdd,hd_data->d_qdd,NUM_JOINTS*" + \
                                ("num_timesteps*" if not single_call_timing else "") + "sizeof(T),cudaMemcpyDeviceToHost));",
                                "gpuErrchk(cudaDeviceSynchronize());"])
    # finally report out timing if requested
    if single_call_timing:
        self.gen_add_code_line("printf(\"Single Call ABA %fus\\n\",time_delta_us_timespec(start,end)/static_cast<double>(num_timesteps));")
    self.gen_add_end_function()

def gen_aba(self, use_thread_group = False):
    # first generate the inner helper
    self.gen_aba_inner(use_thread_group)
    # then generate the device wrapper
    self.gen_aba_device(use_thread_group)
    # then generate the kernels
    self.gen_aba_kernel(use_thread_group, True)
    self.gen_aba_kernel(use_thread_group, False)
    # then generate the host wrappers
    self.gen_aba_host(0)
    self.gen_aba_host(1)
    self.gen_aba_host(2)
    
