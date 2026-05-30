def gen_aba_inner_floating(self):
    NJ = self.robot.get_num_joints()
    nv = self.robot.get_num_vel()
    n_bfs_levels = self.robot.get_max_bfs_level() + 1

    IAOffset = 0
    vcrossOffset = 36 * NJ
    cOffset = 72 * NJ
    pAOffset = 78 * NJ
    UOffset = 84 * NJ
    paOffset = 90 * NJ
    dOffset = 96 * NJ
    uOffset = 97 * NJ
    tempMatOffset = 98 * NJ
    tempVecOffset = 134 * NJ
    fbUOffset = 140 * NJ
    fbDOffset = fbUOffset + 36
    fbDinvOffset = fbDOffset + 36
    fbRhsOffset = fbDinvOffset + 36
    fbInvTempOffset = fbRhsOffset + 6

    func_params = ["s_qdd is the vector of joint accelerations", \
                "s_va is a pointer to shared memory of size 2*6*NUM_BODIES = " + str(12*NJ), \
                "s_q is the vector of joint positions", \
                "s_qd is the vector of joint velocities", \
                "s_tau is the vector of generalized forces", \
                "s_temp is the (shared) scratch; size ABA_INNER_SMEM_BYTES<T, TEMP_IN_SMEM>() (the band when TEMP_IN_SMEM, else 0)", \
                "d_workspace is the global scratch. !TEMP_IN_SMEM: the whole band (ABA_INNER_WORKSPACE_BYTES). TEMP_IN_SMEM && !COLD_IN_SMEM: the cold slab d_cold (ABA_INNER_COLD_BYTES = vcross 36*NJ + fb* root tail 138, packed back-to-back). Pass nullptr at PERF (TEMP_IN_SMEM && COLD_IN_SMEM).", \
                "gravity is the gravity constant"]
    func_def_start = "void aba_inner("
    func_def_middle = "T *s_qdd, T *s_va, const T *s_q, const T *s_qd, const T *s_tau, "
    func_def_end = "T *s_temp, T *d_workspace, const T gravity) {"
    func_notes = ["Assumes the XI matricies have already been updated for the given q",
                  "Floating-base implementation keeps the scalar-joint ABA recursion and solves the 6x6 root block explicitly.",
                  "Inner-controlled placement, two orthogonal levers decided at the top:",
                  "  TEMP_IN_SMEM=false                : whole scratch band -> d_workspace (blunt MINIMAL fallback).",
                  "  TEMP_IN_SMEM=true, COLD_IN_SMEM=true  : PERF, everything in s_temp (byte-identical to the original).",
                  "  TEMP_IN_SMEM=true, COLD_IN_SMEM=false : SURGICAL -- hot recursion stays in s_temp, only the cold vcross slab [36*NJ,72*NJ) and the fb* root tail [140*NJ,140*NJ+138) spill to d_cold (=d_workspace sub-offset), packed back-to-back.",
                  "Caller sizes the arenas from ABA_INNER_{SMEM,WORKSPACE,COLD}_BYTES."]
    func_def_middle, func_params = self.gen_insert_helpers_func_def_params(func_def_middle, func_params, -2)
    func_def = func_def_start + func_def_middle + func_def_end
    self.gen_add_func_doc("Computes the Floating-Base Articulated Body Algorithm", func_notes, func_params, None)
    self.gen_add_code_line("template <typename T, bool TEMP_IN_SMEM = true, bool COLD_IN_SMEM = true>")
    self.gen_add_code_line("__device__")
    self.gen_add_code_line(func_def, True)
    # Inner-controlled scratch-band placement (two levers):
    #   !TEMP_IN_SMEM           : whole band -> d_workspace (s_temp reassigned).
    #   TEMP_IN_SMEM,!COLD_IN_SMEM: hot band stays in s_temp; the cold vcross
    #     slab [36*NJ,72*NJ) and the fb* root tail [140*NJ,140*NJ+138) are
    #     repointed through s_vcross_cold / s_fb_cold to d_cold (=d_workspace),
    #     packed back-to-back (vcross at d_cold[0..36*NJ), fb tail after it).
    # The biases line up the absolute offsets onto d_cold; when COLD_IN_SMEM
    # both pointers are just s_temp -> byte-identical to the original.
    self.gen_add_code_line("if constexpr (!TEMP_IN_SMEM) { s_temp = d_workspace; } else { (void)d_workspace; }")
    self.gen_add_code_line("T *s_vcross_cold = s_temp;")
    self.gen_add_code_line("T *s_fb_cold = s_temp;")
    self.gen_add_code_line("if constexpr (TEMP_IN_SMEM && !COLD_IN_SMEM) { s_vcross_cold = d_workspace - " + str(vcrossOffset) + "; s_fb_cold = d_workspace + " + str(36 * NJ - fbUOffset) + "; }")
    temp_size = self.gen_aba_inner_temp_mem_size()
    self.gen_linalg_smem_setup(temp_size)
    self.gen_add_code_line("// Recursive floating ABA root-port.")

    self.gen_add_code_line("// Initialize IA = I and clear c")
    self.gen_add_parallel_loop("ind", str(36 * NJ + 6 * NJ))
    self.gen_add_code_line("if (ind < " + str(36 * NJ) + ") { s_temp[" + str(IAOffset) + " + ind] = s_XImats[" + str(36 * NJ) + " + ind]; }")
    self.gen_add_code_line("else { s_temp[" + str(cOffset) + " + ind - " + str(36 * NJ) + "] = static_cast<T>(0); }")
    self.gen_add_end_control_flow()
    self.gen_add_sync()

    self.gen_add_code_line("//")
    self.gen_add_code_line("// Forward Pass")
    self.gen_add_code_line("//")
    for bfs_level in range(n_bfs_levels):
        inds = self.robot.get_ids_by_bfs_level(bfs_level)
        joint_names = [self.robot.get_joint_by_id(ind).get_name() for ind in inds]
        link_names = [self.robot.get_link_by_id(ind).get_name() for ind in inds]
        self.gen_add_code_line("// forward pass where bfs_level is " + str(bfs_level))
        self.gen_add_code_line("//     joints are: " + ", ".join(joint_names))
        self.gen_add_code_line("//     links are: " + ", ".join(link_names))

        if bfs_level == 0:
            self.gen_add_parallel_loop("row", "6")
            self.gen_add_code_line("int fb_col = row < 3 ? row + 3 : row - 3;")
            self.gen_add_code_line("s_va[row] = s_qd[fb_col];")
            self.gen_add_code_line("s_temp[" + str(cOffset) + " + row] = static_cast<T>(0);")
            self.gen_add_end_control_flow()
            self.gen_add_sync()
            continue

        for jid in inds:
            parent = self.robot.get_parent_id(jid)
            S_ind = self.robot.get_S_index_by_id(jid)
            S_sign = self.robot.get_S_sign_by_id(jid)
            dof = jid + 5
            jid6 = 6 * jid
            parent6 = 6 * parent
            self.gen_add_code_line("// v[" + str(jid) + "] = X[" + str(jid) + "]*v[" + str(parent) + "] + S*qdot")
            self.gen_add_code_line(f"grid_linalg_row_strided_gemv<T,6,6,6>(&s_XImats[{36*jid}], &s_va[{parent6}], &s_va[{jid6}], static_cast<T>(1), static_cast<T>(0), s_linalg_smem);")
            self.gen_add_serial_ops()
            self.gen_add_code_line(f"s_va[{jid6 + S_ind}] += static_cast<T>({S_sign}) * s_qd[{dof}];")
            self.gen_add_end_control_flow()
            self.gen_add_sync()
            self.gen_add_serial_ops()
            self.gen_mx_func_call_for_cpp([jid], updated_var_names = dict(S_ind_name = str(S_ind),
                                                                          s_dst_name = "&s_temp[" + str(cOffset + jid6) + "]",
                                                                          s_src_name = "&s_va[" + str(jid6) + "]",
                                                                          s_scale_name = "static_cast<T>(" + str(S_sign) + ") * s_qd[" + str(dof) + "]"),
                                              PEQ_FLAG = False, SCALE_FLAG = True)
            self.gen_add_end_control_flow()
            self.gen_add_sync()

    self.gen_add_code_line("// Initialize vcross[k]")
    self.gen_add_parallel_loop("jid", str(NJ))
    self.gen_add_code_line("int jid6 = 6 * jid;")
    self.gen_add_code_line("vcross<T>(&s_vcross_cold[" + str(vcrossOffset) + " + 36*jid], &s_va[jid6]);")
    self.gen_add_end_control_flow()
    self.gen_add_sync()

    self.gen_add_code_line("// temp[k] = -vcross.T*I[k]")
    self.gen_add_parallel_loop("ind", str(36 * NJ))
    self.gen_add_code_line("int row = ind % 6; int col = (ind / 6) % 6; int jid = ind / 36;")
    self.gen_add_code_line("int jid6 = 6 * jid;")
    self.gen_add_code_line("s_temp[" + str(tempMatOffset) + " + jid6*6 + row + col*6] = -dot_prod<T,6,1,1>(&s_vcross_cold[" + str(vcrossOffset) + " + 36*jid + row*6], &s_XImats[" + str(36 * NJ) + " + 36*jid + col*6]);")
    self.gen_add_end_control_flow()
    self.gen_add_sync()

    self.gen_add_code_line("// pA[k] = temp[k]*v[k]")
    self.gen_add_parallel_loop("ind", str(6 * NJ))
    self.gen_add_code_line("int row = ind % 6; int jid = ind / 6;")
    self.gen_add_code_line("int jid6 = 6 * jid;")
    self.gen_add_code_line("s_temp[" + str(pAOffset) + " + jid6 + row] = dot_prod<T,6,6,1>(&s_temp[" + str(tempMatOffset) + " + 6*jid6 + row], &s_va[jid6]);")
    self.gen_add_end_control_flow()
    self.gen_add_sync()

    self.gen_add_code_line("//")
    self.gen_add_code_line("// Backward Pass")
    self.gen_add_code_line("//")
    for jid in range(NJ - 1, -1, -1):
        jid6 = 6 * jid
        parent = self.robot.get_parent_id(jid)
        if jid == 0:
            self.gen_add_code_line("// floating-base root: U = IA*S, D = S^T*U")
            self.gen_add_parallel_loop("ind", "36")
            self.gen_add_code_line("int row = ind % 6; int col = ind / 6;")
            self.gen_add_code_line("int S_col = col < 3 ? col + 3 : col - 3;")
            self.gen_add_code_line("int S_row = row < 3 ? row + 3 : row - 3;")
            self.gen_add_code_line("s_fb_cold[" + str(fbUOffset) + " + ind] = s_temp[" + str(IAOffset) + " + row + 6*S_col];")
            self.gen_add_code_line("s_fb_cold[" + str(fbDOffset) + " + ind] = s_temp[" + str(IAOffset) + " + S_row + 6*S_col];")
            # Ainv=I pre-init dropped 2026-05-29: glass::invertMatrix_dense seeds Ainv internally.
            self.gen_add_end_control_flow()
            self.gen_add_sync()
            self.gen_add_code_line("invert_matrix(6, &s_fb_cold[" + str(fbDOffset) + "], &s_fb_cold[" + str(fbDinvOffset) + "], &s_fb_cold[" + str(fbInvTempOffset) + "]);")
            self.gen_add_parallel_loop("col", "6")
            self.gen_add_code_line("int S_col = col < 3 ? col + 3 : col - 3;")
            self.gen_add_code_line("s_fb_cold[" + str(fbRhsOffset) + " + col] = s_tau[col] - s_temp[" + str(pAOffset) + " + S_col];")
            self.gen_add_end_control_flow()
            self.gen_add_sync()
            continue

        S_ind = self.robot.get_S_index_by_id(jid)
        S_sign = self.robot.get_S_sign_by_id(jid)
        dof = jid + 5
        self.gen_add_code_line("// scalar joint " + str(jid) + ": U, d, u")
        self.gen_add_parallel_loop("row", "6")
        self.gen_add_code_line("s_temp[" + str(UOffset + jid6) + " + row] = static_cast<T>(" + str(S_sign) + ") * s_temp[" + str(IAOffset + 36 * jid + 6 * S_ind) + " + row];")
        self.gen_add_end_control_flow()
        self.gen_add_sync()
        self.gen_add_serial_ops()
        self.gen_add_code_line("s_temp[" + str(dOffset + jid) + "] = static_cast<T>(" + str(S_sign) + ") * s_temp[" + str(UOffset + jid6 + S_ind) + "];")
        self.gen_add_code_line("s_temp[" + str(uOffset + jid) + "] = s_tau[" + str(dof) + "] - static_cast<T>(" + str(S_sign) + ") * s_temp[" + str(pAOffset + jid6 + S_ind) + "] - dot_prod<T,6,1,1>(&s_temp[" + str(UOffset + jid6) + "], &s_temp[" + str(cOffset + jid6) + "]);")
        self.gen_add_end_control_flow()
        self.gen_add_sync()

        if parent != -1:
            parent6 = 6 * parent
            self.gen_add_code_line("// transform U into the parent frame for joint " + str(jid))
            self.gen_add_code_line(f"grid_linalg_gemv<T,6,6,true>(&s_XImats[{36*jid}], &s_temp[{UOffset + jid6}], &s_temp[{tempVecOffset}], static_cast<T>(1), static_cast<T>(0));")
            self.gen_add_parallel_loop("row", "6")
            self.gen_add_code_line("s_temp[" + str(UOffset + jid6) + " + row] = s_temp[" + str(tempVecOffset) + " + row];")
            self.gen_add_end_control_flow()
            self.gen_add_sync()

            self.gen_add_code_line("// temp = X.T*IA*X")
            self.gen_add_code_line(f"grid_linalg_gemm<T,6,6,6,false,true>(&s_XImats[{36*jid}], &s_temp[{IAOffset + 36*jid}], &s_temp[{tempVecOffset}], static_cast<T>(1), static_cast<T>(0), s_linalg_smem);")
            self.gen_add_code_line(f"grid_linalg_gemm<T,6,6,6>(&s_temp[{tempVecOffset}], &s_XImats[{36*jid}], &s_temp[{tempMatOffset}], static_cast<T>(1), static_cast<T>(0), s_linalg_smem);")
            self.gen_add_parallel_loop("ind", "36")
            self.gen_add_code_line("int row = ind % 6; int col = ind / 6;")
            self.gen_add_code_line("s_temp[" + str(tempMatOffset) + " + row + 6*col] -= s_temp[" + str(UOffset + jid6) + " + row] * s_temp[" + str(UOffset + jid6) + " + col] / s_temp[" + str(dOffset + jid) + "];")
            self.gen_add_end_control_flow()
            self.gen_add_sync()

            self.gen_add_code_line("// pa = X.T*(pA + IA*c) + U*u/d")
            self.gen_add_parallel_loop("row", "6")
            self.gen_add_code_line("s_temp[" + str(paOffset + jid6) + " + row] = s_temp[" + str(pAOffset + jid6) + " + row] + dot_prod<T,6,6,1>(&s_temp[" + str(IAOffset + 36 * jid) + " + row], &s_temp[" + str(cOffset + jid6) + "]);")
            self.gen_add_end_control_flow()
            self.gen_add_sync()
            self.gen_add_code_line(f"grid_linalg_gemv<T,6,6,true>(&s_XImats[{36*jid}], &s_temp[{paOffset + jid6}], &s_temp[{tempVecOffset}], static_cast<T>(1), static_cast<T>(0));")
            self.gen_add_serial_ops()
            for _ind in range(6):
                self.gen_add_code_line(f"s_temp[{tempVecOffset + _ind}] += s_temp[{UOffset + jid6 + _ind}] * s_temp[{uOffset + jid}] / s_temp[{dOffset + jid}];")
            self.gen_add_end_control_flow()
            self.gen_add_sync()
            self.gen_add_parallel_loop("ind", "42")
            self.gen_add_code_line("int row = ind % 6; int col = ind / 6;")
            self.gen_add_code_line("if (ind < 36) { s_temp[" + str(IAOffset + 36 * parent) + " + row + 6*col] += s_temp[" + str(tempMatOffset) + " + row + 6*col]; }")
            self.gen_add_code_line("else { s_temp[" + str(pAOffset + parent6) + " + row] += s_temp[" + str(tempVecOffset) + " + row]; }")
            self.gen_add_end_control_flow()
            self.gen_add_sync()

    self.gen_add_code_line("//")
    self.gen_add_code_line("// Second Forward Pass")
    self.gen_add_code_line("//")
    for bfs_level in range(n_bfs_levels):
        inds = self.robot.get_ids_by_bfs_level(bfs_level)
        if bfs_level == 0:
            self.gen_add_code_line("// root acceleration from gravity, then solve root qdd")
            self.gen_add_parallel_loop("ind", "36")
            self.gen_add_code_line("s_temp[" + str(tempMatOffset) + " + ind] = s_XImats[ind];")
            # Ainv=I pre-init (and row/col) dropped 2026-05-29: glass::invertMatrix_dense seeds Ainv internally.
            self.gen_add_end_control_flow()
            self.gen_add_sync()
            self.gen_add_code_line("invert_matrix(6, &s_temp[" + str(tempMatOffset) + "], &s_temp[" + str(tempVecOffset) + "], &s_fb_cold[" + str(fbInvTempOffset) + "]);")
            self.gen_add_parallel_loop("row", "6")
            self.gen_add_code_line("s_va[" + str(6 * NJ) + " + row] = s_temp[" + str(tempVecOffset) + " + row + 6*5] * gravity;")
            self.gen_add_end_control_flow()
            self.gen_add_sync()
            self.gen_add_parallel_loop("row", "6")
            self.gen_add_code_line("s_fb_cold[" + str(fbRhsOffset) + " + row] -= dot_prod<T,6,1,1>(&s_fb_cold[" + str(fbUOffset) + " + 6*row], &s_va[" + str(6 * NJ) + "]);")
            self.gen_add_end_control_flow()
            self.gen_add_sync()
            self.gen_add_parallel_loop("row", "6")
            self.gen_add_code_line("s_qdd[row] = dot_prod<T,6,6,1>(&s_fb_cold[" + str(fbDinvOffset) + " + row], &s_fb_cold[" + str(fbRhsOffset) + "]);")
            self.gen_add_end_control_flow()
            self.gen_add_sync()
            self.gen_add_parallel_loop("row", "6")
            self.gen_add_code_line("int fb_col = row < 3 ? row + 3 : row - 3;")
            self.gen_add_code_line("s_va[" + str(6 * NJ) + " + row] += s_qdd[fb_col];")
            self.gen_add_end_control_flow()
            self.gen_add_sync()
            continue

        for jid in inds:
            parent = self.robot.get_parent_id(jid)
            S_ind = self.robot.get_S_index_by_id(jid)
            S_sign = self.robot.get_S_sign_by_id(jid)
            dof = jid + 5
            jid6 = 6 * jid
            parent6 = 6 * parent
            self.gen_add_serial_ops()
            self.gen_add_code_line("T tempval = s_temp[" + str(uOffset + jid) + "] - dot_prod<T,6,1,1>(&s_temp[" + str(UOffset + jid6) + "], &s_va[" + str(6 * NJ + parent6) + "]);")
            self.gen_add_code_line("s_qdd[" + str(dof) + "] = tempval / s_temp[" + str(dOffset + jid) + "];")
            self.gen_add_end_control_flow()
            self.gen_add_sync()
            self.gen_add_code_line(f"grid_linalg_row_strided_gemv<T,6,6,6>(&s_XImats[{36*jid}], &s_va[{6*NJ + parent6}], &s_va[{6*NJ + jid6}], static_cast<T>(1), static_cast<T>(0), s_linalg_smem);")
            self.gen_add_parallel_loop("row", "6")
            self.gen_add_code_line("s_va[" + str(6 * NJ + jid6) + " + row] += s_temp[" + str(cOffset + jid6) + " + row];")
            self.gen_add_end_control_flow()
            self.gen_add_sync()
            self.gen_add_serial_ops()
            self.gen_add_code_line("s_va[" + str(6 * NJ + jid6 + S_ind) + "] += static_cast<T>(" + str(S_sign) + ") * s_qdd[" + str(dof) + "];")
            self.gen_add_end_control_flow()
            self.gen_add_sync()

    self.gen_add_end_function()


def gen_aba_inner(self):
    if self.robot.floating_base:
        return gen_aba_inner_floating(self)
    n = self.robot.get_num_joints()
    n_bfs_levels = self.robot.get_max_bfs_level() + 1 # starts at 0
	# construct the boilerplate and function definition
    func_params = ["s_qdd is the vector of joint accelerations", \
                "s_va is a pointer to shared memory of size 2*6*NUM_JOINTS = " + str(12*n), \
                "s_q is the vector of joint positions", \
                "s_qd is the vector of joint velocities", \
                "s_tau is the vector of joint torques", \
                "s_temp is the (shared) scratch; size ABA_INNER_SMEM_BYTES<T, TEMP_IN_SMEM>() (the 140*NJ+ band when TEMP_IN_SMEM, else 0)", \
                "d_workspace is the global scratch. !TEMP_IN_SMEM: the whole band (ABA_INNER_WORKSPACE_BYTES). TEMP_IN_SMEM && !COLD_IN_SMEM: the cold slab d_cold (ABA_INNER_COLD_BYTES = the [98*NJ,140*NJ) tempMat slab). Pass nullptr at PERF (TEMP_IN_SMEM && COLD_IN_SMEM).", \
                "gravity is the gravity constant"]
    func_def_start = "void aba_inner("
    func_def_middle = "T *s_qdd, T *s_va, const T *s_q, const T *s_qd, const T *s_tau, "
    func_def_end = "T *s_temp, T *d_workspace, const T gravity) {"
    func_notes = ["Assumes the XI matricies have already been updated for the given q",
                  "Inner-controlled placement, two orthogonal levers decided at the top:",
                  "  TEMP_IN_SMEM=false                : whole scratch band -> d_workspace (blunt MINIMAL fallback).",
                  "  TEMP_IN_SMEM=true, COLD_IN_SMEM=true  : PERF, everything in s_temp (byte-identical to the original).",
                  "  TEMP_IN_SMEM=true, COLD_IN_SMEM=false : SURGICAL -- hot recursion stays in s_temp, only the cold tempMat slab [98*NJ,140*NJ) spills to d_cold (=d_workspace sub-offset).",
                  "Caller sizes the arenas from ABA_INNER_{SMEM,WORKSPACE,COLD}_BYTES."]
    func_def_middle, func_params = self.gen_insert_helpers_func_def_params(func_def_middle, func_params, -2)
    func_def = func_def_start + func_def_middle + func_def_end
    self.gen_add_func_doc("Computes the Articulated Body Algorithm", func_notes, func_params, None)
    self.gen_add_code_line("template <typename T, bool TEMP_IN_SMEM = true, bool COLD_IN_SMEM = true>")
    self.gen_add_code_line("__device__")
    self.gen_add_code_line(func_def, True)
    # Inner-controlled scratch-band placement (two levers):
    #   !TEMP_IN_SMEM           : whole band -> d_workspace (s_temp reassigned).
    #   TEMP_IN_SMEM,!COLD_IN_SMEM: hot band stays in s_temp, the cold tempMat
    #     slab [98*n,140*n) is repointed through s_cold to d_cold (=d_workspace).
    # s_cold is biased by the cold base so the cold s_cold[98*n+...] references
    # land at d_cold[0...]; when COLD_IN_SMEM it is just s_temp -> byte-identical.
    self.gen_add_code_line("if constexpr (!TEMP_IN_SMEM) { s_temp = d_workspace; } else { (void)d_workspace; }")
    self.gen_add_code_line("T *s_cold = s_temp;")
    self.gen_add_code_line("if constexpr (TEMP_IN_SMEM && !COLD_IN_SMEM) { s_cold = d_workspace - " + str(98 * n) + "; }")
    temp_size = self.gen_aba_inner_temp_mem_size()
    self.gen_linalg_smem_setup(temp_size)

    #
    # Initial Debug Prints if Requested
    #
    if self.DEBUG_MODE:
        self.gen_add_sync()
        self.gen_add_serial_ops()
        self.gen_add_code_line("printf(\"q\\n\"); printMat<T,1," + str(n) + ">(s_q,1);")
        self.gen_add_code_line("printf(\"qd\\n\"); printMat<T,1," + str(n) + ">(s_qd,1);")
        self.gen_add_code_line("for (int i = 0; i < " + str(n) + "; i++){printf(\"X[%d]\\n\",i); printMat<T,6,6>(&s_XImats[36*i],6);}")
        self.gen_add_code_line("for (int i = 0; i < " + str(n) + "; i++){printf(\"I[%d]\\n\",i); printMat<T,6,6>(&s_XImats[36*(i+" + str(n) + ")],6);}")
        self.gen_add_end_control_flow()
        self.gen_add_sync()

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
        S_sign_cpp = self.gen_topology_S_sign_for_cpp(inds)

        if bfs_level == 0:
            self.gen_add_code_line("// s_v where parent is base")
            self.gen_add_code_line("//     joints are: " + ", ".join(joint_names))
            self.gen_add_code_line("//     links are: " + ", ".join(link_names))
            # compute the initial v which is just S*qd
            self.gen_add_code_line("// s_v[k] = S[k]*qd[k]")
            if len(inds) > 1:
                self.gen_add_parallel_loop("ind",str(6*len(inds)))
                self.gen_add_code_line("int row = ind % 6;")
                select_var_vals = [("int", "jid", [str(jid) for jid in inds])]
                self.gen_add_multi_threaded_select("ind", "<", [str(6*(i+1)) for i in range(len(inds))], select_var_vals)
                jid = "jid"
            else:
                self.gen_add_parallel_loop("row",str(6))
                jid = str(inds[0])
                self.gen_add_code_line("int jid = " + jid + ";")
            # load in 0 to v 
            self.gen_add_code_lines(["int jid6 = 6*jid;", \
                                     "s_va[jid6 + row] = static_cast<T>(0);"])
            # add in qd
            self.gen_add_code_line("if (row == " + S_ind_cpp + "){s_va[jid6 + " + S_ind_cpp + "] += (" + S_sign_cpp + ") * s_qd[" + jid + "];}")
            self.gen_add_end_control_flow()
            self.gen_add_sync()

            # add debug if requested
            if self.DEBUG_MODE:
                self.gen_add_sync()
                self.gen_add_serial_ops()
                for ind in inds:
                    self.gen_add_code_line("printf(\"s_v[" + str(ind) + "]\\n\"); printMat<T,1,6>(&s_va[6*" + str(ind) + "],1);")
                self.gen_add_end_control_flow()
                self.gen_add_sync()
        
        else:
            self.gen_add_code_line("// s_v where bfs_level is " + str(bfs_level))
            self.gen_add_code_line("//     joints are: " + ", ".join(joint_names))
            self.gen_add_code_line("//     links are: " + ", ".join(link_names))
            self.gen_add_code_line("// s_v[k] = X[k]*v[parent_k] + S[k]*qd[k]")
            # per-jid row_strided_gemv for v
            for jid_val in inds:
                parent_val = self.robot.get_parent_id(jid_val)
                s_ind_val = self.robot.get_S_index_by_id(jid_val)
                s_sign_val = self.robot.get_S_sign_by_id(jid_val)
                self.gen_add_code_line(f"grid_linalg_row_strided_gemv<T,6,6,6>(&s_XImats[{36*jid_val}], &s_va[{6*parent_val}], &s_va[{6*jid_val}], static_cast<T>(1), static_cast<T>(0), s_linalg_smem);")
                self.gen_add_serial_ops()
                self.gen_add_code_line(f"s_va[{6*jid_val + s_ind_val}] += ({s_sign_val}) * s_qd[{jid_val}];")
                self.gen_add_end_control_flow()
                self.gen_add_sync()

        
            # add debug if requested
            if self.DEBUG_MODE:
                self.gen_add_sync()
                self.gen_add_serial_ops()
                for ind in inds:
                    self.gen_add_code_line("printf(\"s_v[" + str(ind) + "] = X*s_v[" + parent_ind_cpp + "] + S*qd[" + str(ind) + "]\\n\"); printMat<T,1,6>(&s_va[6*" + str(ind) + "],1);")
                self.gen_add_end_control_flow()
                self.gen_add_sync()

    # calculate c
    self.gen_add_code_line("// c[k] = mxS(v[k])*qd[k]")
    self.gen_add_parallel_loop("ind", str(n))
    self.gen_add_code_line("int jid = ind;")
    self.gen_add_code_line("int jid6 = 6 * jid;")
    _, S_ind_cpp = self.gen_topology_helpers_pointers_for_cpp(NO_GRAD_FLAG = True)
    S_sign_cpp = self.gen_topology_S_sign_for_cpp()
    self.gen_mx_func_call_for_cpp(list(range(n)), updated_var_names = dict(S_ind_name = S_ind_cpp, s_dst_name = "&s_temp[72 * " + str(n) + " + jid6]", s_src_name = "&s_va[jid6]", s_scale_name = "(" + S_sign_cpp + ") * s_qd[jid]"), PEQ_FLAG = False, SCALE_FLAG = True)
    self.gen_add_end_control_flow()
    
    # add debug if requested
    if self.DEBUG_MODE:
        self.gen_add_sync()
        self.gen_add_serial_ops()
        self.gen_add_code_line("printf(\"c\\n\"); printMat<T,6,"+str(n)+">(&s_temp[72 * "+str(n)+"], 6);")
        self.gen_add_end_control_flow()
        self.gen_add_sync()

    # set IA = I
    self.gen_add_code_line("// Initialize IA = I")
    self.gen_add_parallel_loop("ind",str(36*n))
    self.gen_add_code_line("s_temp[ind] = s_XImats[" + str(36*n) + " + ind];")
    self.gen_add_end_control_flow()
    
    # initialize vcross from v
    self.gen_add_code_line("// Initialize vcross[k]")
    self.gen_add_parallel_loop("ind", str(n))
    self.gen_add_code_line("int jid = ind;")
    self.gen_add_code_line("int jid6 = 6 * jid;")

    self.gen_add_code_line("vcross<T>(&s_temp[36*("+str(n)+"+jid)], &s_va[jid6]);")
    self.gen_add_end_control_flow()
    self.gen_add_sync()

    self.gen_add_code_line("// temp[k] = -vcross.T*I[k]")
    self.gen_add_parallel_loop("ind", str(36*n))
    self.gen_add_code_line("int row = ind % 6; int col = (ind / 6) %6; int jid = ind / 36;")
    self.gen_add_code_line("int jid6 = 6 * jid;")
    self.gen_add_code_line("s_cold[98 * " + str(n) + " + jid6*6 + row+col*6] = -1 * dot_prod<T,6,1,1>(&s_temp[36*("+str(n)+"+jid)+row*6], &s_XImats[36 * ("+str(n)+"+jid) + col*6]);")
    self.gen_add_end_control_flow()
    self.gen_add_sync()
    # calculate pA
    self.gen_add_code_line("// pA[k] = temp[k]*v[k][0]")
    self.gen_add_parallel_loop("ind", str(6*n))
    self.gen_add_code_line("int row = ind % 6; int comp = ind / 6; int jid = comp % " + str(n) + ";")
    self.gen_add_code_line("int jid6 = 6 * jid;")
    self.gen_add_code_line("s_temp[78 * " + str(n) + " + jid6 + row] = dot_prod<T,6,6,1>(&s_cold[98 * " + str(n) + " + 6*jid6+row], &s_va[jid6]);")

    self.gen_add_end_control_flow()

    if self.DEBUG_MODE:
        self.gen_add_sync()
        self.gen_add_serial_ops()
        self.gen_add_code_line("for (int i = 0; i < " + str(n) + "; i++){printf(\"IA[%d]\\n\",i); printMat<T,6,6>(&s_temp[36*(i)],6);}")
        self.gen_add_code_line("printf(\"pA\\n\"); printMat<T,6,"+str(n)+">(&s_temp[78 * "+str(n)+"], 6);")
        self.gen_add_end_control_flow()
        self.gen_add_sync()

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
        S_sign_cpp = self.gen_topology_S_sign_for_cpp(inds)
        self.gen_add_code_line("// Backward pass where bfs_level is " + str(bfs_level))
        self.gen_add_code_line("//     joints are: " + ", ".join(joint_names))
        self.gen_add_code_line("//     links are: " + ", ".join(link_names))
        # caclulate U, which is just IA*S
        self.gen_add_code_line("// U[k] = IA[k]*S[k]")
        self.gen_add_parallel_loop("ind", str(6*len(inds)))
        self.gen_add_code_line("int row = ind % 6;")
        if len(inds) > 1:
            select_var_vals = [("int", "jid", [str(jid) for jid in inds])]
            self.gen_add_multi_threaded_select("ind", "<", [str(6*(i+1)) for i in range(len(inds))], select_var_vals)
            jid = "jid"
        else:
            jid = str(inds[0])
            self.gen_add_code_line("int jid = " + jid + ";")
        self.gen_add_code_line("int jid6 = 6 * " + jid + ";")

        self.gen_add_code_line("s_temp[84*"+str(n)+"+jid6+row] = (" + S_sign_cpp + ") * s_temp[36*jid+row+6*("+ S_ind_cpp+")];")
        self.gen_add_end_control_flow()
        self.gen_add_sync()

        # caclulate d which is S*U and u which is tau - S*pA
        self.gen_add_code_line("// d[k] = S[k]*U[k], u[k] = tau[k] - S[k].T*pA[k]")
        self.gen_add_parallel_loop("ind", str(len(inds)))
        if len(inds) > 1:
            select_var_vals = [("int", "jid", [str(jid) for jid in inds])]
            self.gen_add_multi_threaded_select("ind", "<", [str((i+1)) for i in range(len(inds))], select_var_vals)
            jid = "jid"
        else:
            jid = str(inds[0])
            self.gen_add_code_line("int jid = " + jid + ";")
        self.gen_add_code_line("int jid6 = 6 * " + jid + ";")
        self.gen_add_code_line("s_temp[96 * "+ str(n) +" + jid] = (" + S_sign_cpp + ") * s_temp[84 * " + str(n) + " + jid6 + " + S_ind_cpp + "];")
        
        self.gen_add_code_line("T tempval = (" + S_sign_cpp + ") * s_temp[78 * " + str(n) + " + jid6 + " + S_ind_cpp +"];") 
        self.gen_add_code_line("s_temp[97 * " + str(n) + " + jid] = s_tau[jid] - tempval;")
        self.gen_add_end_control_flow()
        self.gen_add_sync()
        
        # calculate Ia from IA, U, and d
        self.gen_add_code_line("// Ia[k] = IA[k] - U[k]*U[k].T/d[k]")
        self.gen_add_parallel_loop("ind", str(36 * len(inds)))
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
        self.gen_add_sync()

        # caclulate pa
        self.gen_add_code_line("// pa[k] = pA[k] + Ia[k]*c[k]+U[k]*u[k]/d[k]")
        self.gen_add_parallel_loop("ind", str(6*len(inds)))
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
            if len(inds) > 1 and self.robot.has_repeated_parents(inds):
                # Repeated parents: keep atomic dot_prod loops to avoid write conflicts
                self.gen_add_code_line("// temp[k] = X[k].T*Ia[k]")
                self.gen_add_parallel_loop("ind", str(36 * len(inds)))
                self.gen_add_code_line("int row = ind % 6; int col = (ind / 6) %6;")
                select_var_vals = [("int", "jid", [str(jid) for jid in inds])]
                self.gen_add_multi_threaded_select("ind", "<", [str(36*(i+1)) for i in range(len(inds))], select_var_vals)
                self.gen_add_code_line("int jid6 = 6 * jid;")
                self.gen_add_code_line("s_cold[98 * " + str(n) + " + 6 * jid6 + row + 6*col] = dot_prod<T,6,1,1>(&s_XImats[6*jid6+6*row], &s_temp[36 * "+str(n)+"+jid6*6+6*col]);")
                self.gen_add_end_control_flow()
                self.gen_add_sync()
                # update IA of the parent
                self.gen_add_code_line("// IA[parent] += temp[k]*X[k]")
                self.gen_add_parallel_loop("ind", str(36 * len(inds)))
                self.gen_add_code_line("int row = ind % 6; int col = (ind / 6) %6;")
                select_var_vals = [("int", "jid", [str(jid) for jid in inds])]
                self.gen_add_multi_threaded_select("ind", "<", [str(36*(i+1)) for i in range(len(inds))], select_var_vals)
                self.gen_add_code_line("int jid6 = 6 * jid;")
                self.gen_add_code_line("T prodtemp = static_cast<T>(0);")
                self.gen_add_code_line("prodtemp =  dot_prod<T,6,6,1>(&s_cold[98 * " + str(n) + " + 6 * jid6 + row], &s_XImats[6*jid6+6*col]);")
                self.gen_add_code_line("atomicAdd(&s_temp[36 * " + parent_ind_cpp +" + row + 6*col], prodtemp);")
                self.gen_add_end_control_flow()
                self.gen_add_sync()
            else:
                # GEMM path: X^T*Ia*X per jid (single jid or all-distinct parents)
                for jid_val in inds:
                    parent_val = self.robot.get_parent_id(jid_val)
                    self.gen_add_code_line("// X[" + str(jid_val) + "].T*Ia[" + str(jid_val) + "]*X[" + str(jid_val) + "] -> IA[" + str(parent_val) + "]")
                    self.gen_add_code_line(f"grid_linalg_gemm<T,6,6,6,false,true>(&s_XImats[{36*jid_val}], &s_temp[{36*(n+jid_val)}], &s_cold[{98*n + 36*jid_val}], static_cast<T>(1), static_cast<T>(0), s_linalg_smem);")
                    self.gen_add_code_line(f"grid_linalg_gemm<T,6,6,6>(&s_cold[{98*n + 36*jid_val}], &s_XImats[{36*jid_val}], &s_temp[{36*parent_val}], static_cast<T>(1), static_cast<T>(1), s_linalg_smem);")

            # update pA of the parent (sequential GEMVs safe even for repeated parents)
            self.gen_add_code_line("// pA[parent] += X[k].T*pa[k]")
            for jid_val in inds:
                parent_val = self.robot.get_parent_id(jid_val)
                self.gen_add_code_line(f"grid_linalg_gemv<T,6,6,true>(&s_XImats[{36*jid_val}], &s_temp[{90*n + 6*jid_val}], &s_temp[{78*n + 6*parent_val}], static_cast<T>(1), static_cast<T>(1));")

    # add debug if requested
    if self.DEBUG_MODE:
        self.gen_add_sync()
        self.gen_add_serial_ops()
        self.gen_add_code_line("printf(\"U \\n\"); printMat<T,6,"+str(n)+">(&s_temp[84 * "+str(n)+"], 6);")
        self.gen_add_code_line("printf(\"d \\n\"); printMat<T,1,"+str(n)+">(&s_temp[96 * "+str(n)+"], 1);")
        self.gen_add_code_line("printf(\"u \\n\"); printMat<T,1,"+str(n)+">(&s_temp[97 * "+str(n)+"], 1);")
        self.gen_add_code_line("for (int i = 0; i < " + str(n) + "; i++){printf(\"Ia[%d]\\n\",i); printMat<T,6,6>(&s_temp[36*("+str(n)+"+i)],6);}")
        self.gen_add_code_line("for (int i = 0; i < " + str(n) + "; i++){printf(\"IA[%d]\\n\",i); printMat<T,6,6>(&s_temp[36*(i)],6);}")
        self.gen_add_code_line("printf(\"pA\\n\"); printMat<T,6,"+str(n)+">(&s_temp[78 * "+str(n)+"], 6);")
        self.gen_add_end_control_flow()
        self.gen_add_sync()

    self.gen_add_code_line("//")
    self.gen_add_code_line("// Second Forward Pass")
    self.gen_add_code_line("//")
    for bfs_level in range(n_bfs_levels):
        inds = self.robot.get_ids_by_bfs_level(bfs_level)
        parent_ind_cpp, S_ind_cpp = self.gen_topology_helpers_pointers_for_cpp(inds, NO_GRAD_FLAG = True)
        S_sign_cpp = self.gen_topology_S_sign_for_cpp(inds)
        joint_names = [self.robot.get_joint_by_id(ind).get_name() for ind in inds]
        link_names = [self.robot.get_link_by_id(ind).get_name() for ind in inds]
        # calculate a where parent is base
        if bfs_level == 0:
            self.gen_add_code_line("// s_a, qdd where parent is base")
            self.gen_add_code_line("//     joints are: " + ", ".join(joint_names))
            self.gen_add_code_line("//     links are: " + ", ".join(link_names))
            self.gen_add_code_line("// a[k] = X[k]*gravity_vec + c[k]")
            if len(inds) > 1:
                self.gen_add_parallel_loop("ind",str(6*len(inds)))
                self.gen_add_code_line("int row = ind % 6;")
                select_var_vals = [("int", "jid", [str(jid) for jid in inds])]
                self.gen_add_multi_threaded_select("ind", "<", [str(6*(i+1)) for i in range(len(inds))], select_var_vals)
                jid = "jid"
            else:
                self.gen_add_parallel_loop("row",str(6))
                jid = str(inds[0])
                self.gen_add_code_line("int jid = " + jid + ";")
            self.gen_add_code_line("int jid6 = 6*" + jid + ";")
            self.gen_add_code_line("T gravity_vec[] = {0,0,0,0,0,gravity};")
            self.gen_add_code_line("s_va[6*"+str(n)+"+jid6+row] = dot_prod<T,6,6,1>(&s_XImats[36 * jid + row], &gravity_vec[0]) + s_temp[72*"+str(n)+"+jid6+row];")
            self.gen_add_end_control_flow()
            self.gen_add_sync()
        # calculate a where parent is not base
        else:
            self.gen_add_code_line("// s_a, s_qdd where bfs_level is " + str(bfs_level))
            self.gen_add_code_line("//     joints are: " + ", ".join(joint_names))
            self.gen_add_code_line("//     links are: " + ", ".join(link_names))
            self.gen_add_code_line("// a[k] = X[k]*a[parent] + c[k]")
            # per-jid GEMV: a[jid] = X[jid]*a[parent] + c[jid]
            for jid_val in inds:
                parent_val = self.robot.get_parent_id(jid_val)
                self.gen_add_code_line(f"grid_linalg_row_strided_gemv<T,6,6,6>(&s_XImats[{36*jid_val}], &s_va[{6*n + 6*parent_val}], &s_va[{6*n + 6*jid_val}], static_cast<T>(1), static_cast<T>(0), s_linalg_smem);")
                self.gen_add_parallel_loop("row", "6")
                self.gen_add_code_line(f"s_va[{6*n + 6*jid_val} + row] += s_temp[{72*n + 6*jid_val} + row];")
                self.gen_add_end_control_flow()
                self.gen_add_sync()
        
        # calculate qdd which is (u - U*a)/d 
        self.gen_add_code_line("// qdd[k] = (u[k] - U[k].T*a[k])/d[k]")
        self.gen_add_parallel_loop("ind",str(len(inds)))
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
        self.gen_add_sync()

        # update a by adding qdd*S
        self.gen_add_code_line("// a[k] += qdd[k]*S[k]")
        self.gen_add_parallel_loop("ind",str(6*len(inds)))
        
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
        self.gen_add_code_line("T qdd_val = (row == " + S_ind_cpp + ") * (" + S_sign_cpp + ") * (s_qdd[jid]);")
        self.gen_add_code_line("s_va[6*"+str(n)+"+jid6+row] += qdd_val;")

        self.gen_add_end_control_flow()
        self.gen_add_sync()
    
    # add debug if requested
    if self.DEBUG_MODE:
        self.gen_add_sync()
        self.gen_add_serial_ops()
        self.gen_add_code_line("printf(\"a\\n\"); printMat<T,6,"+str(n)+">(&s_va[6 * "+str(n)+"], 6);")
        self.gen_add_code_line("printf(\"qdd\\n\"); printMat<T,1," + str(n) + ">(s_qdd,1);")
        self.gen_add_end_control_flow()
        self.gen_add_sync()
        

    self.gen_add_end_function()

def gen_aba_inner_temp_mem_size(self):
    n = self.robot.get_num_joints()
    if self.robot.floating_base:
        return max(140 * n + 138, self.gen_forward_dynamics_inner_temp_mem_size(minv_f_in_smem=False))
    return 140 * n

def gen_aba_inner_cold_mem_size(self):
    """Float count of the COLD sub-band aba_inner spills under the surgical
    rung (TEMP_IN_SMEM=true, COLD_IN_SMEM=false). The hot recursion stays in
    s_temp; only this cold slab moves to d_cold (a sub-offset of d_workspace).

    FIXED-base : the tempMat/vcross build-scratch slab [98*n, 140*n) -> 42*n.
    FLOATING   : vcross [36*NJ, 72*NJ) (36*NJ) packed contiguously ahead of the
                 fb* root block tail [140*NJ, 140*NJ+138) (138) -> 36*NJ + 138.
    The two floating regions are laid out back-to-back in d_cold so a single
    d_cold pointer covers them (vcross at d_cold[0..36*NJ), fb tail at
    d_cold[36*NJ..36*NJ+138))."""
    n = self.robot.get_num_joints()
    if self.robot.floating_base:
        return 36 * n + 138
    return 42 * n

def gen_aba_inner_function_call(self, updated_var_names = None,
                                temp_in_smem_expr = "true", cold_in_smem_expr = "true"):
    var_names = dict( \
        s_va_name = "s_va", \
        s_q_name = "s_q", \
        s_qd_name = "s_qd", \
        s_qdd_name = "s_qdd", \
        s_tau_name = "s_tau", \
        s_temp_name = "s_temp", \
        d_workspace_name = "nullptr", \
        gravity_name = "gravity"
    )
    if updated_var_names is not None:
        for key,value in updated_var_names.items():
            var_names[key] = value
    aba_code_start = "aba_inner<T, " + temp_in_smem_expr + ", " + cold_in_smem_expr + ">(" + var_names["s_qdd_name"] + ", " + var_names["s_va_name"] + ", " + var_names["s_q_name"] + ", " + var_names["s_qd_name"] + ", " + var_names["s_tau_name"] + ", "
    aba_code_end = var_names["s_temp_name"] + ", " + var_names["d_workspace_name"] + ", " + var_names["gravity_name"] + ");"
    aba_code_middle = self.gen_insert_helpers_function_call()
    aba_code = aba_code_start + aba_code_middle + aba_code_end
    self.gen_add_code_line(aba_code)

def gen_aba_device(self):
    n = self.robot.get_num_joints()
    nv = self.robot.get_num_vel()
    # construct the boilerplate and function definition
    func_params = ["s_qdd is the vector of joint accelerations", \
                   "s_q is the vector of joint positions", \
                   "s_qd is the vector of joint velocities", \
                    "s_tau is the vector of joint torques", \
                   "d_robotModel is the pointer to the initialized model specific helpers on the GPU (XImats, topology_helpers, etc.)", \
                   "gravity is the gravity constant"]
    func_notes = []
    func_def_start = "void aba_device("
    func_def_middle = "T *s_qdd, const T *s_q, const T *s_qd, const T *s_tau, "
    func_def_end = "const robotModel<T> *d_robotModel, const T gravity) {"
    func_def = func_def_start + func_def_middle + func_def_end

    # then generate the code
    self.gen_add_func_doc("Compute the ABA (Articulated Body Algorithm)",\
                          func_notes,func_params,None)
    self.gen_add_code_line("template <typename T>")
    self.gen_add_code_line("__device__")
    self.gen_add_code_line(func_def, True)

    # add the shared memory variables
    shared_mem_size = self.gen_aba_inner_temp_mem_size()
    self.gen_XImats_helpers_temp_shared_memory_code(shared_mem_size, extra_t_buffers = [("s_va", 12*n)], include_linalg_scratch=True)

    # then load/update XI and run the algo
    self.gen_load_update_XImats_helpers_function_call()
    self.gen_aba_inner_function_call()
    self.gen_add_end_function()

def _aba_surgical_inner_smem_size(self):
    """Float count the smem s_temp arena needs at the SURGICAL rung. The cold
    band relocates to d_cold, but the remaining hot references run up to a
    fixed top offset, so the contiguous arena must reach that offset.
      FIXED   : hot ends at the cold base 98*n  -> reclaims the whole 42*n tail.
      FLOATING: hot tempVec ends at 140*NJ; only the 138-float fb* tail above it
                is reclaimed from smem (the interior vcross slot still spills to
                d_cold but its smem hole cannot be compacted byte-identically)."""
    n = self.robot.get_num_joints()
    if self.robot.floating_base:
        return self.gen_aba_inner_temp_mem_size() - 138
    return 98 * n

def _emit_aba_kernel_body_for_flags(self, nq, nv, n, input_count, level, single_call_timing):
    """Emit aba_kernel body for one tier's spill level.
    level 0 (full)     : s_temp in smem, whole inner arena in smem (PERF; byte-identical to original).
    level 1 (surgical) : hot recursion stays in smem; only the cold band spills to d_cold
                         (= d_workspace sub-offset GRID_ABA_COLD_OFFSET_BYTES). smem holds the hot arena.
    level 2 (workspace): whole inner arena redirected to L2-pinned workspace; smem holds only extra_t_buffers."""
    use_workspace_temp = (level == 2)
    use_cold_spill     = (level == 1)
    if use_workspace_temp:
        shared_mem_size = 0
    elif use_cold_spill:
        shared_mem_size = _aba_surgical_inner_smem_size(self)
    else:
        shared_mem_size = self.gen_aba_inner_temp_mem_size()
    self.gen_XImats_helpers_temp_shared_memory_code(shared_mem_size, extra_t_buffers = [("s_qdd", nv), ("s_q_qd_tau", input_count), ("s_va", 12*n)], include_linalg_scratch=True)
    self.gen_add_code_line("T *s_q = s_q_qd_tau; T *s_qd = &s_q_qd_tau[" + str(nq) + "]; T *s_tau = &s_q_qd_tau[" + str(nq + nv) + "];")
    # per-timestep workspace base expr (k-indexed in the batched kernel, slot 0 for single-timing)
    ws_base = "&d_workspace[k*GRID_WORKSPACE_BYTES_PER_TIMESTEP<T>()]" if not single_call_timing else "d_workspace"
    if not single_call_timing:
        self.gen_add_parallel_loop("k","NUM_TIMESTEPS",block_level = True)
        self.gen_kernel_load_inputs("q_qd_tau",str(input_count),stride="stride_q_qd")
    else:
        self.gen_kernel_load_inputs("q_qd_tau",str(input_count))
    if use_workspace_temp:
        self.gen_add_code_line("T *aba_d_workspace = reinterpret_cast<T *>(" + ws_base + ");")
        # The whole inner arena spilled to global, so the smem s_temp slot is
        # null. Repoint s_temp at the workspace so the XImats helper's sincos
        # scratch (and the inner) have a valid backing store, not nullptr.
        self.gen_add_code_line("s_temp = aba_d_workspace;")
    elif use_cold_spill:
        # Surgical rung: hot band stays in smem s_temp; only the cold sub-band
        # lives in d_cold, a sub-offset of the per-timestep workspace. Reuse the
        # SO/grad band base (ABA never runs concurrently with SO/grad).
        self.gen_add_code_line("T *aba_d_cold = reinterpret_cast<T *>(" + ws_base + " + GRID_ABA_COLD_OFFSET_BYTES<T>());")
    else:
        self.gen_add_code_line("(void)d_workspace;")
    if not single_call_timing:
        self.gen_add_code_line("// compute")
        self.gen_load_update_XImats_helpers_function_call()
        self.gen_aba_inner_function_call(
            updated_var_names = (dict(d_workspace_name = "aba_d_workspace") if use_workspace_temp else
                                 (dict(d_workspace_name = "aba_d_cold") if use_cold_spill else None)),
            temp_in_smem_expr = ("false" if use_workspace_temp else "true"),
            cold_in_smem_expr = ("false" if use_cold_spill else "true"))
        self.gen_add_sync()
        self.gen_kernel_save_result("qdd",str(nv),stride=str(nv))
        self.gen_add_end_control_flow()
    else:
        self.gen_add_code_line("// compute with NUM_TIMESTEPS as NUM_REPS for timing")
        self.gen_add_code_line("for (int rep = 0; rep < NUM_TIMESTEPS; rep++){", True)
        self.gen_anti_licm_input_reload("q_qd_tau",str(input_count),feedback_from="qdd")
        self.gen_load_update_XImats_helpers_function_call()
        self.gen_aba_inner_function_call(
            updated_var_names = (dict(d_workspace_name = "aba_d_workspace") if use_workspace_temp else
                                 (dict(d_workspace_name = "aba_d_cold") if use_cold_spill else None)),
            temp_in_smem_expr = ("false" if use_workspace_temp else "true"),
            cold_in_smem_expr = ("false" if use_cold_spill else "true"))
        self.gen_anti_licm_output_write("qdd")
        self.gen_add_end_control_flow()
        self.gen_kernel_save_result("qdd",str(nv))


def gen_aba_kernel(self, single_call_timing = False):
    nq = self.robot.get_num_pos()
    nv = self.robot.get_num_vel()
    n = self.robot.get_num_joints()
    input_count = nq + 2 * nv
    func_params = ["d_qdd is the vector of joint accelerations (output)", \
                    "d_workspace is the L2-pinned global spill buffer (used at LITE/MINIMAL on h1_2-scale)", \
                    "d_q_qd_tau is the vector of joint positions, velocities, torques", \
                    "stride_q_qd is the stride between each q, qd", \
                    "d_robotModel is the pointer to the initialized model specific helpers on the GPU (XImats, topology_helpers, etc.)", \
                    "gravity is the gravity constant", \
                    "num_timesteps is the length of the trajectory points we need to compute over (or overloaded as test_iters for timing)"]
    func_notes = []
    func_def_start = "void aba_kernel(T *d_qdd, unsigned char *d_workspace, const T *d_q_qd_tau, const int stride_q_qd, "
    func_def_end = "const robotModel<T> *d_robotModel, const T gravity, const int NUM_TIMESTEPS) {"
    func_def = func_def_start + func_def_end
    if single_call_timing:
        func_def = func_def.replace("kernel(", "kernel_single_timing(")
    self.gen_add_func_doc("Compute the ABA (Articulated Body Algorithm)", func_notes, func_params, None)
    self.gen_add_code_line("template <typename T, int RESOURCE_TIER = GRID_DEFAULT_RESOURCE_TIER>")
    self.gen_add_code_line("__global__")
    self.gen_add_code_line("__launch_bounds__(tier_max_threads<RESOURCE_TIER>())")
    self.gen_add_code_line(func_def, True)
    # Surgical-spill ladder, 3 rungs. ABA's inner scratch (140*NJ + 138) keeps
    # its hot recursion in smem and spills only the cold sub-band when it can.
    #   level 0 (full)     : whole arena in smem (PERF; byte-identical to original).
    #   level 1 (surgical) : hot band in smem, cold sub-band -> d_cold.
    #   level 2 (workspace): whole arena -> L2-pinned workspace (blunt fallback).
    # picks[tier] IS the level for that tier (see aba_spill_tier_3way).
    picks = getattr(self, "aba_spill_tier_3way", (0, 0, 0))
    if picks[0] == picks[1] == picks[2]:
        _emit_aba_kernel_body_for_flags(self, nq, nv, n, input_count, picks[0], single_call_timing)
    else:
        tier_names = ("TIER_SHARED", "TIER_LITE", "TIER_MINIMAL")
        for tier_idx, (tier_name, pick) in enumerate(zip(tier_names, picks)):
            head = "if constexpr (RESOURCE_TIER == " + tier_name + ") {" if tier_idx == 0 else \
                   "else if constexpr (RESOURCE_TIER == " + tier_name + ") {"
            self.gen_add_code_line(head, True)
            _emit_aba_kernel_body_for_flags(self, nq, nv, n, input_count, pick, single_call_timing)
            self.gen_add_end_control_flow()
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
    func_def_start = "void aba(gridData<T, KIND> *hd_data, const robotModel<T> *d_robotModel, const T gravity, const int num_timesteps,"
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
    self.gen_add_code_line("template <typename T, gridDataKind KIND = GRID_DATA_ALL>")
    self.gen_add_code_line("__host__")
    self.gen_add_code_line(func_def_start)
    self.gen_add_code_line(func_def_end, True)
    self.gen_add_code_line("static_assert(KIND == GRID_DATA_ALL || KIND == GRID_DATA_DYNAMICS, \"aba requires all-data or dynamics gridData\");")

    func_call_start = "aba_kernel<T><<<block_dimms,thread_dimms,ABA_DYNAMIC_SHARED_MEM_BYTES<T>()>>>(hd_data->d_qdd,hd_data->d_workspace,hd_data->d_q_qd_u,stride_q_qd,"
    func_call_end = "d_robotModel,gravity,num_timesteps);"
    self.gen_add_code_line("int stride_q_qd = NUM_JOINTS + 2*NUM_VEL;")
    if single_call_timing:
        func_call_start = func_call_start.replace("kernel<T>","kernel_single_timing<T>")
    if not compute_only:
        # start code with memory transfer
        self.gen_add_code_lines(["// start code with memory transfer", \
                                 "gpuErrchk(cudaMemcpyAsync(hd_data->d_q_qd_u,hd_data->h_q_qd_u,stride_q_qd*" + \
                                    ("num_timesteps*" if not single_call_timing else "") + "sizeof(T),cudaMemcpyHostToDevice,streams[0]));", \
                                 "gpuErrchkKernel();"])
    # then compute:
    self.gen_add_code_line("// then call the kernel")
    func_call = func_call_start + func_call_end
    func_call_code = [func_call, "gpuErrchkKernel();"]
    # wrap function call in timing (if needed)
    if single_call_timing:
        func_call_code.insert(0,"struct timespec start, end; clock_gettime(CLOCK_MONOTONIC,&start);")
        func_call_code.append("clock_gettime(CLOCK_MONOTONIC,&end);")
    self.gen_add_code_line("gpuErrchk(grid_check_dynamic_shared_memory_bytes(\"aba\", ABA_DYNAMIC_SHARED_MEM_BYTES<T>()));")
    self.gen_add_code_lines(func_call_code)
    if not compute_only:
        # then transfer memory back
        self.gen_add_code_lines(["// finally transfer the result back", \
                                "gpuErrchk(cudaMemcpy(hd_data->h_qdd,hd_data->d_qdd,NUM_VEL*" + \
                                ("num_timesteps*" if not single_call_timing else "") + "sizeof(T),cudaMemcpyDeviceToHost));",
                                "gpuErrchkKernel();"])
    # finally report out timing if requested
    if single_call_timing:
        from ..algo_registry import single_call_printf_line
        self.gen_add_code_line(single_call_printf_line("aba"))
    self.gen_add_end_function()

def gen_aba(self):
    # first generate the inner helper
    self.gen_aba_inner()
    # then generate the device wrapper
    self.gen_aba_device()
    # then generate the kernels
    self.gen_aba_kernel(True)
    self.gen_aba_kernel(False)
    # then generate the host wrappers
    self.gen_aba_host(0)
    self.gen_aba_host(1)
    self.gen_aba_host(2)
    
