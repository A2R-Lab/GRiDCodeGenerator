def gen_inverse_dynamics_gradient_inner_temp_mem_size(self):
        if self.robot_has_mimic_joints():
            # The mimic path emits a DENSE serial fold (6 dense per-body buffers
            # + Iv) rather than the sparse-compressed band, so it needs its own
            # (larger) scratch. Big-NB humanoids route this whole pool to
            # d_workspace at the global-temp tier (SCRATCH_IN_SMEM=false).
            return _id_du_mimic_temp_count(self)
        return self.gen_inverse_dynamics_gradient_temp_layout()["full_count"]

def gen_inverse_dynamics_gradient_temp_layout(self):
    n = self.robot.get_num_vel()
    (dva_cols_per_partial, _, _, df_cols_per_partial, _, _, _) = self.gen_topology_sparsity_helpers_python()
    offset_dv_dq = 0
    offset_dv_dqd = offset_dv_dq + 6*dva_cols_per_partial
    offset_da_dq = offset_dv_dqd + 6*dva_cols_per_partial
    offset_da_dqd = offset_da_dq + 6*dva_cols_per_partial
    offset_df_dq = offset_da_dqd + 6*dva_cols_per_partial
    offset_df_dqd = offset_df_dq + 6*df_cols_per_partial
    offset_fxvi = offset_df_dqd + 6*df_cols_per_partial
    offset_mxxv = offset_fxvi + 36*n
    offset_mxxa = offset_mxxv + 6*n
    offset_mxv = offset_mxxa + 6*n
    offset_mxf = offset_mxv + 6*n
    offset_iv = offset_mxf + 6*n
    full_count = offset_iv + 6*n
    spill_start = offset_da_dq
    spill_end = offset_fxvi
    spill_count = spill_end - spill_start
    return {
        "dva_cols_per_partial": dva_cols_per_partial,
        "df_cols_per_partial": df_cols_per_partial,
        "offset_dv_dq": offset_dv_dq,
        "offset_dv_dqd": offset_dv_dqd,
        "offset_da_dq": offset_da_dq,
        "offset_da_dqd": offset_da_dqd,
        "offset_df_dq": offset_df_dq,
        "offset_df_dqd": offset_df_dqd,
        "offset_fxvi": offset_fxvi,
        "offset_mxxv": offset_mxxv,
        "offset_mxxa": offset_mxxa,
        "offset_mxv": offset_mxv,
        "offset_mxf": offset_mxf,
        "offset_iv": offset_iv,
        "full_count": full_count,
        "spill_start": spill_start,
        "spill_end": spill_end,
        "spill_count": spill_count,
        "selective_shared_count": full_count - spill_count,
    }

def _rewrite_id_du_temp_accesses_for_spill(code):
    def replace_accesses(text, address_of):
        token = "&s_temp[" if address_of else "s_temp["
        out = []
        i = 0
        while i < len(text):
            start = text.find(token, i)
            if start < 0:
                out.append(text[i:])
                break
            out.append(text[i:start])
            idx_start = start + len(token)
            depth = 1
            j = idx_start
            while j < len(text) and depth > 0:
                if text[j] == "[":
                    depth += 1
                elif text[j] == "]":
                    depth -= 1
                j += 1
            idx_expr = text[idx_start:j-1]
            ptr_expr = "grid_id_du_temp_ptr<T, USE_DA_DF_SPILL>(s_temp, d_temp_spill, " + idx_expr + ")"
            out.append(ptr_expr if address_of else "(*" + ptr_expr + ")")
            i = j
        return "".join(out)

    return replace_accesses(replace_accesses(code, True), False)

def gen_inverse_dynamics_gradient_inner_function_call(self, updated_var_names = None):
    var_names = dict( \
        s_dc_du_name = "s_dc_du", \
        s_vaf_name = "s_vaf", \
        s_q_name = "s_q", \
        s_qd_name = "s_qd", \
        s_qdd_name = "s_qdd", \
        s_temp_name = "s_temp", \
        d_temp_spill_name = "nullptr", \
        temp_spill_flag_name = "false", \
        gravity_name = "gravity"
    )
    if updated_var_names is not None:
        for key,value in updated_var_names.items():
            var_names[key] = value
    id_du_code_start = "inverse_dynamics_gradient_inner<T, " + var_names["temp_spill_flag_name"] + ">(" + var_names["s_dc_du_name"] + ", " + var_names["s_q_name"] + ", " + var_names["s_qd_name"] + ", "
    id_du_code_middle = var_names["s_vaf_name"] + ", " + self.gen_insert_helpers_function_call()
    id_du_code_end = var_names["s_temp_name"] + ", " + var_names["d_temp_spill_name"] + ", " + var_names["gravity_name"] + ");"
    id_du_code = id_du_code_start + id_du_code_middle + id_du_code_end
    self.gen_add_code_line(id_du_code)

def gen_inverse_dynamics_gradient_inner(self):
    function_start = len(self.code_str)
    n = self.robot.get_num_vel()
    NJ = self.robot.get_num_joints()
    max_bfs_levels = self.robot.get_max_bfs_level()
    n_bfs_levels = max_bfs_levels + 1 # starts at 0

    # construct the boilerplate and function definition
    func_params = ["s_dc_du is a pointer to memory for the final result of size 2*NUM_JOINTS*NUM_JOINTS = " + str(2*n*n), \
                   "s_q is the vector of joint positions", \
                   "s_qd is the vector of joint velocities", \
                   "s_vaf are the helper intermediate variables computed by inverse_dynamics", \
                   "s_temp is a pointer to helper shared memory of size 66*NUM_JOINTS + 6*sparse_dv,da,df_col_needs = " + \
                            str(self.gen_inverse_dynamics_gradient_inner_temp_mem_size()), \
                   "gravity is the gravity constant"]
    func_def_start = "void inverse_dynamics_gradient_inner(T *s_dc_du, const T *s_q, const T *s_qd, const T *s_vaf, "
    func_def_end = "T *s_temp, T *d_temp_spill, const T gravity) {"
    func_def_start, func_params = self.gen_insert_helpers_func_def_params(func_def_start, func_params, -2)
    func_notes = ["Assumes s_XImats is updated already for the current s_q",
                  "This is the id_du band sub-inner (the stable surface composed by fd_du / integrator_gradient). It does NOT own s_temp placement; the USE_DA_DF_SPILL band selectively spills its da_dq..fxvi band to d_temp_spill via grid_id_du_temp_ptr<T, USE_DA_DF_SPILL>. The whole-pool placement is owned by the wrapping inverse_dynamics_gradient_device."]
    func_def = func_def_start + func_def_end
    # then generate the code
    self.gen_add_func_doc("Computes the gradient of inverse dynamics",func_notes,func_params,None)
    self.gen_add_code_line("template <typename T, bool USE_DA_DF_SPILL = false>")
    self.gen_add_code_line("__device__")
    self.gen_add_code_line(func_def, True)

    if self.robot_has_mimic_joints():
        # MIMIC (T3-finisher P3): the sparse NJ-indexed gradient assembly below
        # writes s_dc_du by raw body id and assumes NJ == NV, so it can't fold a
        # mimic model (NB > NV, shared v-slots). Emit instead a DENSE serial
        # reduced-space fold that mirrors RBDReference.rnea_grad exactly:
        # alpha-scaled velocity/accel reads in the forward pass + per-body
        # v-slot accumulate (+=) in the backward pass. Correctness, not perf, is
        # the goal here (mimic robots are the gripper/hand class). The large
        # dense per-body buffers spill to d_temp_spill / d_workspace via the
        # inner's mimic temp-size; see gen_inverse_dynamics_gradient_inner_temp_mem_size.
        _gen_id_du_mimic_inner(self, n, NJ)
        self.gen_add_end_function()
        return

    #
    # Optimize memory requirements due to sparsity induced by branching
    # requires more complex pointer math but/and saves a lot of space
    #
    (dva_cols_per_partial, dva_cols_per_jid, running_sum_dva_cols_per_jid, \
      df_cols_per_partial, df_cols_per_jid, running_sum_df_cols_per_jid, df_col_that_is_jid) = self.gen_topology_sparsity_helpers_python()
    if not self.robot.floating_base:
        self.gen_add_code_line("//")
        self.gen_add_code_line("// dv and da need " + str(dva_cols_per_partial) + " cols per dq,dqd")
        self.gen_add_code_line("// df needs " + str(df_cols_per_partial) + " cols per dq,dqd")
        self.gen_add_code_line("//    out of a possible " + str(n*n) + " cols per dq,dqd")
        self.gen_add_code_line("// Gradients are stored compactly as dv_i/dq_[0...a], dv_i+1/dq_[0...b], etc")
        self.gen_add_code_line("//    where a and b are the needed number of columns")
        self.gen_add_code_line("//")
    # gen som aditional helpers
    running_sum_delta_df_dva_cols_per_jid = [running_sum_df_cols_per_jid[jid] - running_sum_dva_cols_per_jid[jid] for jid in range(NJ)]

    # add shared memory note
    Offset_dv_dq = 0
    Offset_dv_dqd = Offset_dv_dq + 6*dva_cols_per_partial
    Offset_da_dq = Offset_dv_dqd + 6*dva_cols_per_partial
    Offset_da_dqd = Offset_da_dq + 6*dva_cols_per_partial
    Offset_df_dq = Offset_da_dqd + 6*dva_cols_per_partial
    Offset_df_dqd = Offset_df_dq + 6*df_cols_per_partial
    Offset_FxvI = Offset_df_dqd + 6*df_cols_per_partial
    Offset_MxXv = Offset_FxvI + 36*n
    Offset_MxXa = Offset_MxXv + 6*n
    Offset_Mxv = Offset_MxXa + 6*n
    Offset_Mxf = Offset_Mxv + 6*n
    Offset_Iv = Offset_Mxf + 6*n
    # Offset_dva_cols = Offset_Iv + 6*n
    # Offset_df_cols = Offset_dva_cols + n

    self.gen_add_code_line("// Temp memory offsets are as follows:")
    self.gen_add_code_line("// T *s_dv_dq = &s_temp[" + str(Offset_dv_dq) + "]; " + \
                              "T *s_dv_dqd = &s_temp[" + str(Offset_dv_dqd) + "]; " + \
                              "T *s_da_dq = &s_temp[" + str(Offset_da_dq) + "];")
    self.gen_add_code_line("// T *s_da_dqd = &s_temp[" + str(Offset_da_dqd) + "]; " + \
                              "T *s_df_dq = &s_temp[" + str(Offset_df_dq) + "]; " + \
                              "T *s_df_dqd = &s_temp[" + str(Offset_df_dqd) + "];")
    self.gen_add_code_line("// T *s_FxvI = &s_temp[" + str(Offset_FxvI) + "]; T *s_MxXv = &s_temp[" + str(Offset_MxXv) + "]; " + \
                              "T *s_MxXa = &s_temp[" + str(Offset_MxXa) + "];")
    self.gen_add_code_line("// T *s_Mxv = &s_temp[" + str(Offset_Mxv) + "]; T *s_Mxf = &s_temp[" + str(Offset_Mxf) + "]; " + \
                              "T *s_Iv = &s_temp[" + str(Offset_Iv) + "];")

    if self.DEBUG_MODE:
        self.gen_add_sync()
        self.gen_add_serial_ops()
        self.gen_add_code_lines(["printf(\"-------------------------\\n\");", \
                                 "printf(\"Validating Function Inputs\\n\");", \
                                 "printf(\"-------------------------\\n\");", \
                                 "printf(\"q\\n\"); printMat<T,1," + str(n) + ">(s_q,1);", \
                                 "printf(\"qd\\n\"); printMat<T,1," + str(n) + ">(s_qd,1);", \
                                 "printf(\"vaf-v\\n\"); printMat<T,6," + str(n) + ">(s_vaf,6);", \
                                 "printf(\"vaf-a\\n\"); printMat<T,6," + str(n) + ">(&s_vaf[6*" + str(n) + "],6);", \
                                 "printf(\"vaf-f\\n\"); printMat<T,6," + str(n) + ">(&s_vaf[12*" + str(n) + "],6);"])
        self.gen_add_code_line("for (int i = 0; i < " + str(n) + "; i++){printf(\"X[%d]\\n\",i); printMat<T,6,6>(&s_XImats[36*i],6);}")
        self.gen_add_code_line("for (int i = 0; i < " + str(n) + "; i++){printf(\"I[%d]\\n\",i); printMat<T,6,6>(&s_XImats[36*(i+" + str(n) + ")],6);}")
        self.gen_add_code_line("printf(\"-------------------------\\n\");")
        self.gen_add_end_control_flow()
        self.gen_add_sync()
    
    #
    # Initial temp comps
    #
    self.gen_add_code_line("//")
    self.gen_add_code_line("// Initial Temp Comps")
    self.gen_add_code_line("//")
    # first compute temporary values by type of operation
    # we can use part of FxvI temp mem for Xv and Xa initial comps also compute Iv
    self.gen_add_code_line("// First compute Imat*v and Xmat*v_parent, Xmat*a_parent (store in FxvI for now)")
    self.gen_add_code_line("// Note that if jid_parent == -1 then v_parent = 0 and a_parent = gravity")
    self.gen_add_parallel_loop("ind",str(6*3*NJ))
    self.gen_add_code_line("int row = ind % 6; int col = ind / 6; int jid = col % " + str(NJ) + "; int jid6 = 6*jid;")
    # get the parent (note that in some cases we have more efficient ways of computing this so add some special cases)
    parent_ind_cpp, S_ind_cpp = self.gen_topology_helpers_pointers_for_cpp(NO_GRAD_FLAG = True, OFFSET=False)
    S_sign_cpp = self.gen_topology_S_sign_for_cpp(OFFSET=False)
    self.gen_add_code_line("bool parentIsBase = " + parent_ind_cpp + " == -1;")
    # then get the offsets
    self.gen_add_code_lines(["bool comp1 = col < " + str(NJ) + "; bool comp3 = col >= " + str(2*NJ) + ";",
                             "int XIOffset  =  comp1 * " + str(36*NJ) + " + 6*jid6 + row; // rowCol of I (comp1) or X (comp 2 and 3)",
                             "int vaOffset  = comp1 * jid6 + !comp1 * 6*" + parent_ind_cpp + " + comp3 * " + str(6*NJ) + "; // v_i (comp1) or va_parent (comp 2 and 3)",
                             "int dstOffset = comp1 * " + str(Offset_Iv) + " + !comp1 * " + str(Offset_FxvI) + " + comp3 * " + str(6*NJ) + " + jid6 + row; // rowCol of dst"])
    if self.robot.floating_base:
        self.gen_add_code_lines(["s_temp[dstOffset] = (parentIsBase && !comp1) ?",
                                 "                           (comp3 ? (row < 3 ? static_cast<T>(0) : s_XImats[6*jid6 + 6*row + 5] * gravity) : static_cast<T>(0)) :",
                                 "                           dot_prod<T,6,6,1>(&s_XImats[XIOffset],&s_vaf[vaOffset]);"])
    else:
        self.gen_add_code_lines(["s_temp[dstOffset] = (parentIsBase && !comp1) ? comp3 * s_XImats[XIOffset + 30] * gravity : ",
                                 "                                               dot_prod<T,6,6,1>(&s_XImats[XIOffset],&s_vaf[vaOffset]);"])
    self.gen_add_end_control_flow()
    self.gen_add_sync()

    if self.DEBUG_MODE:
        self.gen_add_sync()
        self.gen_add_serial_ops()
        self.gen_add_code_lines(["printf(\"-------------------------\\n\");", \
                                 "printf(\"Temp Comps Part 1\\n\");", \
                                 "printf(\"-------------------------\\n\");", \
                                 "printf(\"Iv\\n\"); printMat<T,6," + str(n) + ">(&s_temp[" + str(Offset_Iv) + "],6);", \
                                 "printf(\"Xv\\n\"); printMat<T,6," + str(n) + ">(&s_temp[" + str(Offset_FxvI) + "],6);", \
                                 "printf(\"Xa\\n\"); printMat<T,6," + str(n) + ">(&s_temp[" + str(Offset_FxvI + 6*n) + "],6);"])
        self.gen_add_code_line("printf(\"-------------------------\\n\");")
        self.gen_add_end_control_flow()
        self.gen_add_sync()

    # then do the mx comps
    self.gen_add_code_line("// Then compute Mx(Xv), Mx(Xa), Mx(v), Mx(f)")
    self.gen_add_parallel_loop("col",str(4*n))
    self.gen_add_code_line("int dof_id = col / 4; int selector = col % 4; int dof_id6 = 6*dof_id;")
    if self.robot.floating_base: self.gen_add_code_line("int jid = dof_id < 6 ? 0 : dof_id - 5; int jid6 = jid*6;") # First 6 dof belong to fb
    else: self.gen_add_code_line("int jid6 = dof_id6;")
    select_var_vals = [("int", "dstOffset", [str(Offset_MxXv), str(Offset_MxXa), str(Offset_Mxv), str(Offset_Mxf)]), \
                       ("const T *", "src", ["&s_temp[" + str(Offset_FxvI) + "]", "&s_temp[" + str(Offset_FxvI + 6*NJ) + "]", \
                                       "&s_vaf[0]", "&s_vaf[" + str(12*NJ) + "]"])]
    self.gen_add_multi_threaded_select("selector", "==", [str(i) for i in range(4)], select_var_vals)
    if 'jid' in S_ind_cpp: S_ind_cpp = S_ind_cpp.replace('jid', 'dof_id')
    if 'jid' in S_sign_cpp: S_sign_cpp = S_sign_cpp.replace('jid', 'dof_id')
    updated_var_names = dict(S_ind_name = S_ind_cpp, s_dst_name = "&s_temp[dstOffset + dof_id6]", s_src_name = "&src[jid6]", s_scale_name = S_sign_cpp)
    self.gen_mx_func_call_for_cpp(PEQ_FLAG = False, SCALE_FLAG = True, updated_var_names = updated_var_names)
    self.gen_add_end_control_flow()
    self.gen_add_sync()

    has_linear_axis = any(
        self.robot.get_S_index_by_id(jid) >= 3 for jid in range(n)
    )
    if has_linear_axis:
        # For prismatic axes, the force derivative term needs the force cross-product
        # column. The motion and force columns coincide for the revolute axes covered
        # by the original path, but differ for linear axes.
        self.gen_add_parallel_loop("dof_id", str(n))
        _, S_ind_dof_cpp = self.gen_topology_helpers_pointers_for_cpp(
            list(range(n)),
            updated_var_names=dict(jid_name="dof_id"),
            NO_GRAD_FLAG=True,
            OFFSET=False,
        )
        S_sign_dof_cpp = self.gen_topology_S_sign_for_cpp(
            list(range(n)),
            updated_var_names=dict(jid_name="dof_id"),
            OFFSET=False,
        )
        self.gen_add_code_line("int S_ind = " + S_ind_dof_cpp + ";")
        self.gen_add_code_line("T S_sign = static_cast<T>(" + S_sign_dof_cpp + ");")
        self.gen_add_code_line("if (S_ind >= 3) {", True)
        self.gen_add_code_line(f"T *dst = &s_temp[{Offset_Mxf} + 6*dof_id];")
        if self.robot.floating_base:
            self.gen_add_code_line("int jid = dof_id < 6 ? 0 : dof_id - 5;")
            self.gen_add_code_line(f"const T *src = &s_vaf[{12*NJ} + 6*jid];")
        else:
            self.gen_add_code_line(f"const T *src = &s_vaf[{12*NJ} + 6*dof_id];")
        self.gen_add_code_line("for (int row = 0; row < 6; ++row) dst[row] = static_cast<T>(0);")
        self.gen_add_code_line("if (S_ind == 3) {", True)
        self.gen_add_code_line("dst[1] = S_sign * src[5];")
        self.gen_add_code_line("dst[2] = -S_sign * src[4];")
        self.gen_add_end_control_flow()
        self.gen_add_code_line("else if (S_ind == 4) {", True)
        self.gen_add_code_line("dst[0] = -S_sign * src[5];")
        self.gen_add_code_line("dst[2] = S_sign * src[3];")
        self.gen_add_end_control_flow()
        self.gen_add_code_line("else {", True)
        self.gen_add_code_line("dst[0] = S_sign * src[4];")
        self.gen_add_code_line("dst[1] = -S_sign * src[3];")
        self.gen_add_end_control_flow()
        self.gen_add_end_control_flow()
        self.gen_add_end_control_flow()
        self.gen_add_sync()

    if self.DEBUG_MODE:
        self.gen_add_sync()
        self.gen_add_serial_ops()
        self.gen_add_code_lines(["printf(\"-------------------------\\n\");", \
                                 "printf(\"Temp Comps Part 2\\n\");", \
                                 "printf(\"-------------------------\\n\");", \
                                 "printf(\"Mx(Xv)\\n\"); printMat<T,6," + str(n) + ">(&s_temp[" + str(Offset_MxXv) + "],6);", \
                                 "printf(\"Mx(Xa)\\n\"); printMat<T,6," + str(n) + ">(&s_temp[" + str(Offset_MxXa) + "],6);", \
                                 "printf(\"Mx(v)\\n\"); printMat<T,6," + str(n) + ">(&s_temp[" + str(Offset_Mxv) + "],6);",\
                                 "printf(\"Mx(f)\\n\"); printMat<T,6," + str(n) + ">(&s_temp[" + str(Offset_Mxf) + "],6);"])
        self.gen_add_code_line("printf(\"-------------------------\\n\");")
        self.gen_add_end_control_flow()
        self.gen_add_sync()

    #
    # FORWARD PASS
    #
    self.gen_add_code_line("//")
    self.gen_add_code_line("// Forward Pass")
    self.gen_add_code_line("//")
    self.gen_add_code_line("// We start with dv/du noting that we only have values")
    self.gen_add_code_line("//    for ancestors and for the current index else 0")
    # then serial dv/du in bfs waves
    for bfs_level in range(n_bfs_levels):
        inds = self.robot.get_ids_by_bfs_level(bfs_level)
        joint_names = [self.robot.get_joint_by_id(ind).get_name() for ind in inds]
        link_names = [self.robot.get_link_by_id(ind).get_name() for ind in inds]
        _, S_ind_cpp, dva_col_offset_for_jid_cpp, _, dva_col_offset_for_parent_cpp, _, _, _ = self.gen_topology_helpers_pointers_for_cpp(inds, OFFSET=False)
        S_sign_cpp = self.gen_topology_S_sign_for_cpp(inds, OFFSET=False)
        self.gen_add_code_line("// dv/du where bfs_level is " + str(bfs_level))
        self.gen_add_code_line("//     joints are: " + ", ".join(joint_names))
        self.gen_add_code_line("//     links are: " + ", ".join(link_names))

        # when parent is base dv_dq = 0, dv_dqd = S
        if bfs_level == 0:
            self.gen_add_code_line("// when parent is base dv_dq = 0, dv_dqd = S")
            if self.robot.floating_base:
                self.gen_add_parallel_loop("ind",str(6*2*n))
                self.gen_add_code_line("bool dq_flag = ind < " + str(6*n) + ";")
                self.gen_add_code_line("int row = ind % 6; int col = (!dq_flag * " + str(-n) + ") + (ind / 6);")
                self.gen_add_code_line("int du_offset = dq_flag ? " + str(Offset_dv_dq) + " : " + str(Offset_dv_dqd) + ";")
                self.gen_add_code_line("int fb_col = row < 3 ? row + 3 : row - 3;")
                self.gen_add_code_line("s_temp[du_offset + 6*col + row] = !dq_flag * (fb_col == col);")
            else:
                self.gen_add_parallel_loop("ind",str(6*2*len(inds)))
                if len(inds) > 1:
                    self.gen_add_code_line("int row = ind % 6; int col = ind / 6; int col_du = col % " + str(len(inds)) + "; bool dq_flag = col == col_du;")
                    select_var_vals = [("int", "jid", [str(jid) for jid in inds])]
                    self.gen_add_multi_threaded_select("col_du", "<", [str((i+1)) for i in range(len(inds))], select_var_vals)
                else:
                    self.gen_add_code_line("int row = ind % 6; int dq_flag = (ind / 6) == 0;")
                self.gen_add_code_line("int du_offset = dq_flag ? " + str(Offset_dv_dq) + " : " + str(Offset_dv_dqd) + ";")
                self.gen_add_code_line("s_temp[du_offset + 6*" + dva_col_offset_for_jid_cpp + " + row] = " + \
                                    "(!dq_flag && row == " + S_ind_cpp + ") * static_cast<T>(" + S_sign_cpp + ");") 
            self.gen_add_end_control_flow()

        # dv/du = X dv_parent/du + {MxXv or S for col ind}
        # there are 2*(bfs_level + 1) columns per du with 2*bfs mults with X and then the addition in the last col
        else:
            self.gen_add_code_line("// dv/du = Xmat*dv_parent/du + {Mx(Xv) or S for col ind}")
            self.gen_add_code_line("// first compute dv/du = Xmat*dv_parent/du")
            if self.robot.floating_base:
                self.gen_add_parallel_loop("ind",str(6*2*n*len(inds)))
                self.gen_add_code_line(f"bool dq_flag = ind < {6*n*len(inds)};")
                self.gen_add_code_line(f"int row = ind % 6; int col = (ind / 6) % {n};")
                if len(inds) > 1: 
                    self.gen_add_code_line(f'int ind_du = ind % {6*n*len(inds)};')
                    select_var_vals = [("int", "jid", [str(jid) for jid in inds])]
                    jid = "jid"
                    self.gen_add_multi_threaded_select("(ind_du)", "<", [str((idx+1)*n*6) for idx, jid in enumerate(inds)], select_var_vals)
                    select_var_vals = [("int", "parent_jid", [str(self.robot.get_parent_id(jid)) for jid in inds])]
                    parent_ind_cpp = "parent_jid"
                    self.gen_add_multi_threaded_select("(ind_du)", "<", [str((idx+1)*n*6) for idx, jid in enumerate(inds)], select_var_vals)
                else:
                    jid = inds[0]
                    parent_ind_cpp = self.robot.get_parent_id(jid)
                self.gen_add_code_line(f"int dof_id = {jid}+5; (void)dof_id;")
                self.gen_add_code_line(f"int du_offset = dq_flag ? {Offset_dv_dq} + {6*n}*{jid} : {Offset_dv_dqd} + {6*n}*{jid};")
                self.gen_add_code_line(f"int parent_du_offset = dq_flag ? {Offset_dv_dq} + {6*n}*{parent_ind_cpp} : {Offset_dv_dqd} + {6*n}*{parent_ind_cpp};")
                self.gen_add_code_line("s_temp[du_offset + 6*col + row] = dot_prod<T,6,6,1>(&s_XImats[36 * " + str(jid) + " + row]," + \
                                        " &s_temp[6*col + parent_du_offset]);")
                # then add in S or Mx(Xv); dof_id is only referenced when the S_ind/S_sign cpp expressions actually substitute it.
                self.gen_add_code_line(f"if (col == {jid} + 5)" + ' {', True)
                if 'jid' in S_ind_cpp: S_ind_cpp = S_ind_cpp.replace('jid', 'dof_id')
                if 'jid' in S_sign_cpp: S_sign_cpp = S_sign_cpp.replace('jid', 'dof_id')
                self.gen_add_code_line(f"s_temp[du_offset + 6*col + row] += !dq_flag * (row == {S_ind_cpp})" + \
                                        f" * ({S_sign_cpp}) + dq_flag * s_temp[{Offset_MxXv} + ({jid}+5)*6 + row];")
                self.gen_add_end_control_flow()
            else:
                self.gen_add_parallel_loop("ind",str(6*2*(bfs_level)*len(inds)))
                self.gen_add_code_line("int row = ind % 6; int col = ind / 6; int col_du = col % " + str(bfs_level*len(inds)) + "; " + \
                                                                            "int col_jid = col_du % " + str(bfs_level) + ";")
                if bfs_level > 1 or len(inds) > 1:
                    self.gen_add_code_line("int dq_flag = col == col_du;")
                else:
                    self.gen_add_code_line("int dq_flag = col < 1;")
                if len(inds) > 1:
                    select_var_vals = [("int", "jid", [str(jid) for jid in inds])]
                    jid = "jid"
                    self.gen_add_multi_threaded_select("col_du", "<", [str((i+1)*bfs_level) for i in range(len(inds))], select_var_vals)
                else:
                    jid = str(inds[0])
                self.gen_add_code_line("int du_col_offset = dq_flag * " + str(Offset_dv_dq) + " + !dq_flag * " + str(Offset_dv_dqd) + " + 6 * col_jid;")
                self.gen_add_code_line("s_temp[du_col_offset + 6*" + dva_col_offset_for_jid_cpp + " + row] = ")
                self.gen_add_code_line("    dot_prod<T,6,6,1>(&s_XImats[36*" + jid + " + row]," + \
                                                            "&s_temp[du_col_offset + 6*" + dva_col_offset_for_parent_cpp + "]);")
                self.gen_add_code_line("// then add {Mx(Xv) or S for col ind}")
                # all cols add if bfs_level is 1 so skip the if statement
                if bfs_level > 1:
                    self.gen_add_code_line("if (col_jid == " + str(bfs_level-1) + ") {", True)
                # do the non-branching if/else
                self.gen_add_code_line("s_temp[du_col_offset + 6*" + dva_col_offset_for_jid_cpp + " + 6 + row] = ")
                self.gen_add_code_line("    dq_flag * s_temp[" + str(Offset_MxXv) + " + 6*" + jid + " + row] + " + \
                                        "(!dq_flag && row == " + S_ind_cpp + ") * static_cast<T>(" + S_sign_cpp + ");")
                # all cols add if bfs_level is 1 so skip the if statement
                if bfs_level > 1:
                    self.gen_add_end_control_flow()
            self.gen_add_end_control_flow()
        self.gen_add_sync()

        if self.DEBUG_MODE:
            self.gen_add_sync()
            self.gen_add_serial_ops()
            if bfs_level == 0:
                self.gen_add_code_lines(["printf(\"-------------------------\\n\");", \
                                         "printf(\"dv/du in bfs waves\\n\");", \
                                         "printf(\"-------------------------\\n\");"])
            self.gen_add_code_line("printf(\"dv/du in for bfs wave[%d]\\n\"," + str(bfs_level) + ");")
            for ind in inds:
                self.gen_add_code_lines(["printf(\"dv[%d]/dq\\n\"," + str(ind) + ");", \
                                         "printMat<T,6," + str(bfs_level+1) + ">(&s_temp[" + \
                                                str(Offset_dv_dq + 6*running_sum_dva_cols_per_jid[ind]) + "],6);", \
                                         "printf(\"dv[%d]/dqd\\n\"," + str(ind) + ");", \
                                         "printMat<T,6," + str(bfs_level+1) + ">(&s_temp[" + \
                                                str(Offset_dv_dqd + 6*running_sum_dva_cols_per_jid[ind]) + "],6);"])
            self.gen_add_end_control_flow()
            self.gen_add_sync()

    # Start the da/du comp with da/du = MxS(dv/du)*qd + {MxXa, Mxv}
    self.gen_add_code_line("// start da/du by setting = MxS(dv/du)*qd + {MxXa, Mxv} for all n in parallel")
    self.gen_add_code_line("// start with da/du = MxS(dv/du)*qd")
    _ , S_ind_cpp , _ , _ , _ , _ , dva_col_offset_for_jidp1_cpp, _ = self.gen_topology_helpers_pointers_for_cpp(list(range(n)), OFFSET=False)
    S_sign_cpp = self.gen_topology_S_sign_for_cpp(OFFSET=False)
    add_col_for_jid = "(" + dva_col_offset_for_jidp1_cpp + " - 1)"
    if self.robot.floating_base:
        # set da/du = 0
        self.gen_add_code_line("// First zero da/du")
        self.gen_add_parallel_loop('ind',str(2*n*NJ*6))
        self.gen_add_code_line(f"s_temp[{Offset_da_dq} + ind] = static_cast<T>(0);")
        self.gen_add_end_control_flow()
        # Sync before the += accumulation below: the zeroing loop and the
        # MxS(dv/du)*qd accumulation write the same s_temp[Offset_da_dq] region
        # from different threads. Without this barrier the result is correct
        # only within a single warp (<=32 threads) and races at larger blocks.
        self.gen_add_sync()
        if 'jid' in S_ind_cpp: S_ind_cpp = S_ind_cpp.replace('jid', 'dof_id')
        if 'jid' in S_sign_cpp: S_sign_cpp = S_sign_cpp.replace('jid', 'dof_id')
        # Axis-indexed S helpers for the serialized root accumulation below.
        S_ind_ax = S_ind_cpp.replace('dof_id', 'ax')
        S_sign_ax = S_sign_cpp.replace('dof_id', 'ax')
        self.gen_add_parallel_loop("col",str(2*n*n))
        self.gen_add_code_line(f"int dof = col % {n};") # column within each joint that is being focused
        self.gen_add_code_line(f"int dof_id = (col / {n}) % {n}; int jid = dof_id < 6 ? 0 : dof_id - 5;") # dof_id being applied with S, to jid
        self.gen_add_code_line(f"bool dq_flag = col < {n*n}; int dqd_offset = !dq_flag * {6*dva_cols_per_partial};")
        # Floating root (jid==0): all 6 root axes (dof_id 0..5) accumulate into
        # the SAME da/du column, but mxX_peq_scaled assumes a single writer per
        # destination. Run the whole accumulation for the root on one lane
        # (dof_id==0), summing over the 6 axes, so there is no multi-thread +=
        # race. (Correct only within one warp otherwise -> wrong J_qv at
        # MAX_PERF_LEVEL_THREADS.) Non-root joints have a unique axis per lane.
        self.gen_add_code_line("if (jid == 0 && dof_id == 0) {", True)
        self.gen_add_code_line(f"T *root_dst = &s_temp[{Offset_da_dq} + dof*6 + dqd_offset];")
        self.gen_add_code_line(f"const T *root_src = &s_temp[{Offset_dv_dq} + dof*6 + dqd_offset];")
        self.gen_add_code_line("for (int ax = 0; ax < 6; ax++) {", True)
        self.gen_mx_func_call_for_cpp(PEQ_FLAG = True, SCALE_FLAG = True, updated_var_names = dict(
            S_ind_name = S_ind_ax, s_dst_name = "root_dst", s_src_name = "root_src",
            s_scale_name = "(" + S_sign_ax + ") * s_qd[ax]"))
        self.gen_add_end_control_flow()  # for ax
        # The {MxXa, Mxv} add applies to root columns dof in [0,6) (axis == column).
        self.gen_add_code_line("if (dof < 6) {", True)
        self.gen_add_code_line(f"int src_offset = dq_flag * {Offset_MxXa} + !dq_flag * {Offset_Mxv} + 6*dof;")
        self.gen_add_code_line("for (int row = 0; row < 6; row++) { root_dst[row] += s_temp[src_offset + row]; }")
        self.gen_add_end_control_flow()  # if dof < 6
        self.gen_add_end_control_flow()  # if jid == 0 && dof_id == 0
        # Non-root joints: one lane per (jid, dof), no collision.
        self.gen_add_code_line("if (jid != 0) {", True)
        updated_var_names = dict(S_ind_name = S_ind_cpp, s_dst_name = f"&s_temp[{Offset_da_dq} + jid*{6*n} + dof*6 + dqd_offset]", \
                                 s_src_name = f"&s_temp[{Offset_dv_dq} + jid*{6*n} + dof*6 + dqd_offset]", s_scale_name = "(" + S_sign_cpp + ") * s_qd[dof_id]")
        self.gen_mx_func_call_for_cpp(PEQ_FLAG = True, SCALE_FLAG = True, updated_var_names = updated_var_names)
        self.gen_add_code_line("// then add {MxXa, Mxv} to the appropriate column")
        self.gen_add_code_line("if (dof == dof_id) {", True)
        self.gen_add_code_line(f"int src_offset = dq_flag * {Offset_MxXa} + !dq_flag * {Offset_Mxv} + 6*dof_id;")
        self.gen_add_code_line(f"for (int row = 0; row < 6; row++) {{ s_temp[{Offset_da_dq} + 6*dof + row + jid*{6*n} + dqd_offset] += s_temp[src_offset + row]; }}")
        self.gen_add_end_control_flow()  # if dof == dof_id
        self.gen_add_end_control_flow()  # if jid != 0
        self.gen_add_end_control_flow()  # parallel loop
        self.gen_add_sync()
    else:
        self.gen_add_parallel_loop("col",str(2*dva_cols_per_partial))
        self.gen_add_code_line("int col_du = col % " + str(dva_cols_per_partial) + ";") # signifies col of corresponding du
        select_var_vals = [("int", "jid", [str(jid) for jid in range(NJ)])]
        self.gen_add_multi_threaded_select("col_du", "<", [str(running_sum_dva_cols_per_jid[jid+1]) for jid in range(NJ)], select_var_vals)
        updated_var_names = dict(S_ind_name = S_ind_cpp, s_dst_name = "&s_temp[" + str(Offset_da_dq) + " + 6*col]", \
                                s_src_name = "&s_temp[" + str(Offset_dv_dq) + " + 6*col]", s_scale_name = "(" + S_sign_cpp + ") * s_qd[jid]")
        # call the mx func
        self.gen_mx_func_call_for_cpp(PEQ_FLAG = False, SCALE_FLAG = True, updated_var_names = updated_var_names)
        # then add to the add col
        self.gen_add_code_lines(["// then add {MxXa, Mxv} to the appropriate column", \
                                "int dq_flag = col == col_du; int src_offset = dq_flag * " + str(Offset_MxXa) + " + !dq_flag * " + str(Offset_Mxv) + " + 6*jid;"])
        self.gen_add_code_line("if(col_du == " + add_col_for_jid + "){", True)
        self.gen_add_code_line("for(int row = 0; row < 6; row++){", True)
        self.gen_add_code_line("s_temp[" + str(Offset_da_dq) + " + 6*col + row] += s_temp[src_offset + row];")
        self.gen_add_end_control_flow()
        self.gen_add_end_control_flow()
        self.gen_add_end_control_flow()
        self.gen_add_sync()

    if self.DEBUG_MODE:
        self.gen_add_sync()
        self.gen_add_serial_ops()
        self.gen_add_code_lines(["printf(\"-------------------------\\n\");", \
                                 "printf(\"da/du part 1 = MxS(dv/du)*qd + {MxXa, Mxf}\\n\");", \
                                 "printf(\"-------------------------\\n\");"])
        for ind in range(n):
            num_cols = self.robot.get_bfs_level_by_id(ind) + 1
            self.gen_add_code_lines(["printf(\"da[%d]/dq\\n\"," + str(ind) + ");", \
                                     "printMat<T,6," + str(num_cols) + ">(&s_temp[" + \
                                            str(Offset_da_dq + 6*running_sum_dva_cols_per_jid[ind]) + "],6);", \
                                     "printf(\"da[%d]/dqd\\n\"," + str(ind) + ");", \
                                     "printMat<T,6," + str(num_cols) + ">(&s_temp[" + \
                                            str(Offset_da_dqd + 6*running_sum_dva_cols_per_jid[ind]) + "],6);"])
        self.gen_add_end_control_flow()
        self.gen_add_sync()

    # then serial da/du in bfs waves
    self.gen_add_code_line("// Finish da/du with parent updates noting that we only have values")
    self.gen_add_code_line("//    for ancestors and for the current index and nothing for bfs 0")
    for bfs_level in range(1,n_bfs_levels):
        inds = self.robot.get_ids_by_bfs_level(bfs_level)
        joint_names = [self.robot.get_joint_by_id(ind).get_name() for ind in inds]
        link_names = [self.robot.get_link_by_id(ind).get_name() for ind in inds]
        _, _, dva_col_offset_for_jid_cpp, _, dva_col_offset_for_parent_cpp, _, _, _ = self.gen_topology_helpers_pointers_for_cpp(inds, OFFSET=False)
        parent_inds = [self.robot.get_parent_id(ind) for ind in inds]
        self.gen_add_code_line("// da/du where bfs_level is " + str(bfs_level))
        self.gen_add_code_line("//     joints are: " + ", ".join(joint_names))
        self.gen_add_code_line("//     links are: " + ", ".join(link_names))
        
        # da/du += X da_parent/du
        # there are 2*(bfs_level + 1) columns per du with 2*bfs mults with X and then the addition in the last col
        self.gen_add_code_line("// da/du += Xmat*da_parent/du")    
        if self.robot.floating_base: 
            self.gen_add_parallel_loop("ind",str(6*2*n*len(inds)))
            self.gen_add_code_line(f"bool dq_flag = ind < {6*n*len(inds)};")
            self.gen_add_code_line(f"int row = ind % 6; int col = (ind / 6) % {n};")
            if len(inds) > 1: 
                self.gen_add_code_line(f'int ind_du = ind % {6*n*len(inds)};')
                select_var_vals = [("int", "jid", [str(jid) for jid in inds])]
                jid = "jid"
                self.gen_add_multi_threaded_select("(ind_du)", "<", [str((idx+1)*n*6) for idx, jid in enumerate(inds)], select_var_vals)
                select_var_vals = [("int", "parent_jid", [str(self.robot.get_parent_id(jid)) for jid in inds])]
                parent_ind_cpp = "parent_jid"
                self.gen_add_multi_threaded_select("(ind_du)", "<", [str((idx+1)*n*6) for idx, jid in enumerate(inds)], select_var_vals)
            else:
                jid = inds[0]
                parent_ind_cpp = self.robot.get_parent_id(jid)
            self.gen_add_code_line(f"int du_offset = dq_flag ? {Offset_da_dq} + {6*n}*{jid} : {Offset_da_dqd} + {6*n}*{jid};")
            self.gen_add_code_line(f"int parent_du_offset = dq_flag ? {Offset_da_dq} + {6*n}*{parent_ind_cpp} : {Offset_da_dqd} + {6*n}*{parent_ind_cpp};")
            self.gen_add_code_line("s_temp[du_offset + 6*col + row] += dot_prod<T,6,6,1>(&s_XImats[36 * " + str(jid) + " + row]," + \
                                        " &s_temp[6*col + parent_du_offset]);")
        else:
            self.gen_add_parallel_loop("ind",str(6*2*bfs_level*len(inds)))
            self.gen_add_code_lines(["int row = ind % 6; int col = ind / 6; int col_du = col % " + str(bfs_level*len(inds)) + ";", \
                                    "int dq_flag = col == col_du; int col_jid = col_du % " + str(bfs_level) + ";"])
            if len(inds) > 1:
                select_var_vals = [("int", "jid", [str(jid) for jid in inds])]
                jid = "jid"
                self.gen_add_multi_threaded_select("col_du", "<", [str((i+1)*bfs_level) for i in range(len(inds))], select_var_vals)
            else:
                jid = str(inds[0])
            self.gen_add_code_line("int du_col_offset = dq_flag * " + str(Offset_da_dq) + " + !dq_flag * " + str(Offset_da_dqd) + " + 6 * col_jid;")
            self.gen_add_code_lines(["s_temp[du_col_offset + 6*" + dva_col_offset_for_jid_cpp + " + row] += ", \
                                    "    dot_prod<T,6,6,1>(&s_XImats[36*" + jid + " + row]," + \
                                                            "&s_temp[du_col_offset + 6*" + dva_col_offset_for_parent_cpp + "]);"])
        self.gen_add_end_control_flow()
        self.gen_add_sync()

        if self.DEBUG_MODE:
            self.gen_add_sync()
            self.gen_add_serial_ops()
            if bfs_level == 1:
                self.gen_add_code_lines(["printf(\"-------------------------\\n\");", \
                                         "printf(\"da/du in bfs waves\\n\");", \
                                         "printf(\"-------------------------\\n\");"])
            self.gen_add_code_line("printf(\"da/du for bfs wave[%d]\\n\"," + str(bfs_level) + ");")
            for ind in inds:
                self.gen_add_code_lines(["printf(\"da[%d]/dq\\n\"," + str(ind) + ");", \
                                         "printMat<T,6," + str(bfs_level+1) + ">(&s_temp[" + \
                                                str(Offset_da_dq + 6*running_sum_dva_cols_per_jid[ind]) + "],6);", \
                                         "printf(\"da[%d]/dqd\\n\"," + str(ind) + ");", \
                                         "printMat<T,6," + str(bfs_level+1) + ">(&s_temp[" + \
                                                str(Offset_da_dqd + 6*running_sum_dva_cols_per_jid[ind]) + "],6);"])
            self.gen_add_end_control_flow()
            self.gen_add_sync()

    
    # Intiialize df/du to 0 to make sure we don't have issues with remaining values later when we do +=
    self.gen_add_code_line("// Init df/du to 0")
    self.gen_add_parallel_loop("ind",str(6*2*df_cols_per_partial))
    self.gen_add_code_line("s_temp[" + str(Offset_df_dq) + " + ind] = static_cast<T>(0);")
    self.gen_add_end_control_flow()
    self.gen_add_sync()

    # Start the df/du by setting = fx(dv/du)*Iv and also compute the temp = Fx(v)*I 
    # aka do all of the Fx comps in parallel
    self.gen_add_code_lines(["// Start the df/du by setting = fx(dv/du)*Iv and also compute the temp = Fx(v)*I ", \
                             "//    aka do all of the Fx comps in parallel", \
                             "// note that while df has more cols than dva the dva cols are the first few df cols"])
    _, _, dva_col_offset_for_jid_cpp, df_col_offset_for_jid_cpp, _, _, _, _ = self.gen_topology_helpers_pointers_for_cpp(list(range(n)), OFFSET=False)
    self.gen_add_parallel_loop("col",str(2*dva_cols_per_partial + 6*NJ))
    self.gen_add_code_line("int col_du = col % " + str(dva_cols_per_partial) + ";")
    if self.robot.floating_base: 
        self.gen_add_code_line(f'int jid = col_du / {n};')
        self.gen_add_code_lines(["// Compute Offsets and Pointers", \
                                "int dq_flag = col == col_du;", \
                                "int Offset_col_du_src = dq_flag * " + str(Offset_dv_dq) + " + !dq_flag * " + str(Offset_dv_dqd) + " + 6*col_du;", \
                                "int Offset_col_du_dst = dq_flag * " + str(Offset_df_dq) + " + !dq_flag * " + str(Offset_df_dqd) + " + 6*col_du;"])
    else:
        select_var_vals = [("int", "jid", [str(jid) for jid in range(NJ)])]
        self.gen_add_multi_threaded_select("col_du", "<", [str(running_sum_dva_cols_per_jid[jid+1]) for jid in range(NJ)], select_var_vals)
        self.gen_add_code_lines(["// Compute Offsets and Pointers", \
                                "int dq_flag = col == col_du; int dva_to_df_adjust = " + df_col_offset_for_jid_cpp + " - " + dva_col_offset_for_jid_cpp + ";", \
                                "int Offset_col_du_src = dq_flag * " + str(Offset_dv_dq) + " + !dq_flag * " + str(Offset_dv_dqd) + " + 6*col_du;", \
                                "int Offset_col_du_dst = dq_flag * " + str(Offset_df_dq) + " + !dq_flag * " + str(Offset_df_dqd) + " + 6*(col_du + dva_to_df_adjust);"])
    self.gen_add_code_line("T *dst = &s_temp[Offset_col_du_dst]; " + \
                        "const T *fx_src = &s_temp[Offset_col_du_src]; " + \
                        "const T *mult_src = &s_temp[" + str(Offset_Iv) + " + 6*jid];")
    # do the adjust for the temp comps
    self.gen_add_code_line("// Adjust pointers for temp comps (if applicable)")
    self.gen_add_code_line("if (col >= " + str(2*dva_cols_per_partial) + ") {", True)
    self.gen_add_code_lines(["int comp = col - " + str(2*dva_cols_per_partial) + "; int comp_col = comp % 6; // int jid = comp / 6;", \
                            "int jid6 = comp - comp_col; int jid36_col6 = 6*jid6 + 6*comp_col;"])
    self.gen_add_code_line("dst = &s_temp[" + str(Offset_FxvI) + " + jid36_col6]; " + \
                        "fx_src = &s_vaf[jid6]; " + \
                        "mult_src = &s_XImats[" + str(36*NJ) + " + jid36_col6];")
    self.gen_add_end_control_flow()
    self.gen_add_code_line("fx_times_v<T>(dst, fx_src, mult_src);")
    self.gen_add_end_control_flow()
    self.gen_add_sync()

    if self.DEBUG_MODE:
        self.gen_add_sync()
        self.gen_add_serial_ops()
        self.gen_add_code_lines(["printf(\"-------------------------\\n\");", \
                                 "printf(\"df/du part 1 = fx(dv/du)*Iv\\n\");", \
                                 "printf(\"     and Temp = Fx(v)*I\\n\");", \
                                 "printf(\"-------------------------\\n\");"])
        for ind in range(n):
            num_cols = df_cols_per_jid[ind]
            self.gen_add_code_lines(["printf(\"df[%d]/dq\\n\"," + str(ind) + ");", \
                                     "printMat<T,6," + str(num_cols) + ">(&s_temp[" + str(Offset_df_dq + 6*running_sum_df_cols_per_jid[ind]) + "],6);", \
                                     "printf(\"df[%d]/dqd\\n\"," + str(ind) + ");", \
                                     "printMat<T,6," + str(num_cols) + ">(&s_temp[" + str(Offset_df_dqd + 6*running_sum_df_cols_per_jid[ind]) + "],6);"])
            self.gen_add_code_lines(["printf(\"Fx(v)*I[%d]\\n\"," + str(ind) + ");", \
                                     "printMat<T,6,6>(&s_temp[" + str(Offset_FxvI) + " + 36*" + str(ind) + "],6);"])
        self.gen_add_end_control_flow()
        self.gen_add_sync()

    # then in parallel finish df/du += I*da/du + FxvI*dv/du
    self.gen_add_code_line("// Then in parallel finish df/du += I*da/du + (Fx(v)I)*dv/du")
    self.gen_add_parallel_loop("ind",str(6*2*dva_cols_per_partial))
    if self.robot.floating_base: 
        self.gen_add_code_line(f"int row = ind % 6; int col = ind / 6; int jid = (ind / {6*n}) % {NJ};")
        self.gen_add_code_lines(["T *df_row_col = &s_temp[" + str(Offset_df_dq) + " + 6*col + row];",
                                "const T *dv_col = &s_temp[" + str(Offset_dv_dq) + " + 6*col]; " + \
                                        "const T *da_col = &s_temp[" + str(Offset_da_dq) + " + 6*col];",
                                "int jid36 = 36*jid; const T *I_row = &s_XImats[" + str(36*NJ) + " + jid36 + row]; " + \
                                                "const T *FxvI_row = &s_temp[" + str(Offset_FxvI) + " + jid36 + row];"])
    else:
        self.gen_add_code_line("int row = ind % 6; int col = ind / 6; int col6 = ind - row; int col_du = (col % " + str(dva_cols_per_partial) + ");")
        select_var_vals = [("int", "jid", [str(jid) for jid in range(NJ)])]
        self.gen_add_multi_threaded_select("col_du", "<", [str(running_sum_dva_cols_per_jid[jid+1]) for jid in range(NJ)], select_var_vals)
        self.gen_add_code_lines(["// Compute Offsets and Pointers", \
                                "int dva_to_df_adjust = " + df_col_offset_for_jid_cpp + " - " + dva_col_offset_for_jid_cpp + ";", \
                                "if (col >= " + str(dva_cols_per_partial) + "){dva_to_df_adjust += " + str(df_cols_per_partial - dva_cols_per_partial) + ";}", \
                                "T *df_row_col = &s_temp[" + str(Offset_df_dq) + " + 6*dva_to_df_adjust + ind];",
                                "const T *dv_col = &s_temp[" + str(Offset_dv_dq) + " + col6]; " + \
                                        "const T *da_col = &s_temp[" + str(Offset_da_dq) + " + col6];",
                                "int jid36 = 36*jid; const T *I_row = &s_XImats[" + str(36*n) + " + jid36 + row]; " + \
                                                "const T *FxvI_row = &s_temp[" + str(Offset_FxvI) + " + jid36 + row];"])
    self.gen_add_code_lines(["// Compute the values", \
                             "*df_row_col += dot_prod<T,6,6,1>(I_row,da_col) + dot_prod<T,6,6,1>(FxvI_row,dv_col);"])
    self.gen_add_end_control_flow()

    if self.DEBUG_MODE:
        self.gen_add_sync()
        self.gen_add_serial_ops()
        self.gen_add_code_lines(["printf(\"-------------------------\\n\");", \
                                 "printf(\"df/du += I*da/du + FxvI*dv/du\\n\");", \
                                 "printf(\"-------------------------\\n\");"])
        for ind in range(n):
            num_cols = df_cols_per_jid[ind]
            self.gen_add_code_lines(["printf(\"df[%d]/dq\\n\"," + str(ind) + ");", \
                                     "printMat<T,6," + str(num_cols) + ">(&s_temp[" + str(Offset_df_dq + 6*running_sum_df_cols_per_jid[ind]) + "],6);", \
                                     "printf(\"df[%d]/dqd\\n\"," + str(ind) + ");", \
                                     "printMat<T,6," + str(num_cols) + ">(&s_temp[" + str(Offset_df_dqd + 6*running_sum_df_cols_per_jid[ind]) + "],6);"])
        self.gen_add_end_control_flow()
        self.gen_add_sync()

    # and also at the same time compute the temp var -X^T * mxf
    # since all temps are done re-use one in practice
    self.gen_add_code_line("// At the same time compute the last temp var: -X^T * mx(f)")
    self.gen_add_code_line("// use Mx(Xv) temp memory as those values are no longer needed")
    self.gen_add_parallel_loop("ind",str(6*n))
    if self.robot.floating_base: 
        self.gen_add_code_line("int XTcol = ind % 6; int jid = ind / 6; int dof_id6 = (jid+5)*6; int jid6 = jid*6;")
        self.gen_add_code_line("s_temp[" + str(Offset_MxXv) + " + ind] = -dot_prod<T,6,1,1>(" + \
                                        "&s_XImats[6*(jid6 + XTcol)], &s_temp[" + str(Offset_Mxf) + " + dof_id6]);")
    else:
        self.gen_add_code_line("int XTcol = ind % 6; int jid6 = ind - XTcol;")
        self.gen_add_code_line("s_temp[" + str(Offset_MxXv) + " + ind] = -dot_prod<T,6,1,1>(" + \
                                        "&s_XImats[6*(jid6 + XTcol)], &s_temp[" + str(Offset_Mxf) + " + jid6]);")
    self.gen_add_end_control_flow()
    self.gen_add_sync()

    if self.DEBUG_MODE:
        self.gen_add_sync()
        self.gen_add_serial_ops()
        self.gen_add_code_lines(["printf(\"-------------------------\\n\");", \
                                 "printf(\"Temp = -X^T * mx(f)\\n\");", \
                                 "printf(\"-------------------------\\n\");"])
        for ind in range(n):
            self.gen_add_code_lines(["printf(\"-X^T*mx(f)[%d]\\n\"," + str(ind) + ");", \
                                     "printMat<T,1,6>(&s_temp[" + str(Offset_MxXv) + " + 6*" + str(ind) + "],1);"])
        self.gen_add_end_control_flow()
        self.gen_add_sync()

    #
    # BACKWARD PASS
    #
    self.gen_add_code_line("//")
    self.gen_add_code_line("// BACKWARD Pass")
    self.gen_add_code_line("//")
    # update df serially (df_lambda/du = X^T * df/du + {Xmx(f), 0})
    for bfs_level in range(max_bfs_levels,0,-1): # STOP AT 1 because updating parent and last is 0 ---- !!!!!
        inds = self.robot.get_ids_by_bfs_level(bfs_level)
        _, _, dva_col_offset_for_jid_cpp, df_col_offset_for_jid_cpp, _, df_col_offset_for_parent_cpp, _, df_col_that_is_jid_cpp = self.gen_topology_helpers_pointers_for_cpp(inds, OFFSET=False)
        joint_names = [self.robot.get_joint_by_id(ind).get_name() for ind in inds]
        link_names = [self.robot.get_link_by_id(ind).get_name() for ind in inds]
        self.gen_add_code_line("// df/du update where bfs_level is " + str(bfs_level))
        self.gen_add_code_line("//     joints are: " + ", ".join(joint_names))
        self.gen_add_code_line("//     links are: " + ", ".join(link_names))

        if self.DEBUG_MODE:
            self.gen_add_sync()
            self.gen_add_serial_ops()
            for ind in self.robot.get_unique_parent_ids(inds):
                self.gen_add_code_lines(["printf(\"df[%d]/dq (parent update) BEFORE UPDATE\\n\"," + str(ind) + ");", \
                                         "printMat<T,6," + str(df_cols_per_jid[ind]) + ">(&s_temp[" + str(Offset_df_dq + \
                                                    6*running_sum_df_cols_per_jid[ind]) + "],6);", \
                                         "printf(\"df[%d]/dqd (parent update) BEFORE UPDATE\\n\"," + str(ind) + ");", \
                                         "printMat<T,6," + str(df_cols_per_jid[ind]) + ">(&s_temp[" + str(Offset_df_dqd + \
                                                    6*running_sum_df_cols_per_jid[ind]) + "],6);"])
            self.gen_add_end_control_flow()
            self.gen_add_sync()

        # df_lambda/du = X^T * df/du + {Xmx(f), 0}
        # there are 2*(bfs_level + 1) columns per du
        df_cols_per_this_bfs = [df_cols_per_jid[ind] for ind in inds]
        curr_cols_per_du = sum(df_cols_per_this_bfs)
        breakpoints = [sum(df_cols_per_this_bfs[0:i+1]) for i in range(len(inds))]
        col_adjusts = [sum(df_cols_per_this_bfs[0:i]) for i in range(len(inds))]
        sparsity_branch_corrector_vals = [str(jid - self.robot.get_parent_id(jid) - 1) for jid in inds]
        sparsity_branch_corrector_needed = any([int(i) for i in sparsity_branch_corrector_vals])
        if not sparsity_branch_corrector_needed:
            sparsity_branch_corrector = str(0)
        self.gen_add_code_line("// df_lambda/du += X^T * df/du + {Xmx(f), 0}")
        if self.robot.floating_base: 
            self.gen_add_parallel_loop("ind",str(6*2*n*len(inds)))
            self.gen_add_code_line(f"int row = ind % 6; int col = (ind / 6) % {n};")
            if len(inds) > 1:
                self.gen_add_code_line(f'int ind_du = ind % {6*n*len(inds)};')
                select_var_vals = [("int", "jid", [str(jid) for jid in inds])]
                jid = "jid"
                self.gen_add_multi_threaded_select("(ind_du)", "<", [str((idx+1)*n*6) for idx, jid in enumerate(inds)], select_var_vals)
                select_var_vals = [("int", "parent_jid", [str(self.robot.get_parent_id(jid)) for jid in inds])]
                
                self.gen_add_multi_threaded_select("(ind_du)", "<", [str((idx+1)*n*6) for idx, jid in enumerate(inds)], select_var_vals)
                parent_jid = f'parent_jid*{6*n}'
        else:
            self.gen_add_parallel_loop("ind",str(6*2*curr_cols_per_du))
            self.gen_add_code_line(f"int row = ind % 6; int col = ind / 6; int col_du = col % {curr_cols_per_du};")
            if len(inds) > 1:
                select_var_vals = [("int", "jid", [str(jid) for jid in inds])]
                select_var_vals.append(("int", "col_adjust", [str(val) for val in col_adjusts]))
                jid = "jid"
                if sparsity_branch_corrector_needed:
                    select_var_vals.append(("int", "sparsity_branch_corrector", sparsity_branch_corrector_vals))
                    sparsity_branch_corrector = "sparsity_branch_corrector"
                self.gen_add_multi_threaded_select("col_du", "<", [str(val) for val in breakpoints], select_var_vals, True)
                adjustments = "col_du -= col_adjust; // adjust for variable number of columns"
            else: adjustments = ''
        if len(inds) <= 1:
            jid = str(inds[0])
            if self.robot.floating_base: parent_jid = self.robot.get_parent_id(inds[0])*6*n
            elif sparsity_branch_corrector_needed:
                sparsity_branch_corrector = str(self.robot.get_parent_id(inds[0]) - inds[0])
                adjustments = ""
        if self.robot.floating_base: self.gen_add_code_line(f"bool dq_flag = ind < {6*n*len(inds)};")
        else: self.gen_add_code_line(f"int dq_flag = col == col_du;")
        if not self.robot.floating_base and adjustments:
            self.gen_add_code_line(adjustments)
        if self.robot.floating_base: 
            self.gen_add_code_line("int du_col_offset = dq_flag * " + str(Offset_df_dq) + " + !dq_flag * " + str(Offset_df_dqd) + " + 6*col;")
            self.gen_add_code_line(f"T *dst = &s_temp[du_col_offset + {parent_jid} + row];")
            self.gen_add_code_lines(["T update_val = dot_prod<T,6,1,1>(&s_XImats[36*" + jid + f" + 6*row],&s_temp[du_col_offset + {jid}*{6*n}])",
                                    f"              + dq_flag * (col == {jid}+5) * s_temp[" + str(Offset_MxXv) + " + 6*" + jid + " + row];"])
        else:
            self.gen_add_code_line("int du_col_offset = dq_flag * " + str(Offset_df_dq) + " + !dq_flag * " + str(Offset_df_dqd) + " + 6*col_du;")
            self.gen_add_code_line("int dst_adjust = (col_du >= " + df_col_that_is_jid_cpp + ") * 6 * " + sparsity_branch_corrector + "; // adjust for sparsity compression offsets")
            self.gen_add_code_line("T *dst = &s_temp[du_col_offset + 6*" + df_col_offset_for_parent_cpp + " + dst_adjust + row];")
            self.gen_add_code_lines(["T update_val = dot_prod<T,6,1,1>(&s_XImats[36*" + jid + " + 6*row],&s_temp[du_col_offset + 6*" + df_col_offset_for_jid_cpp + "])",
                                    "              + dq_flag * (col_du == " + df_col_that_is_jid_cpp + ") * s_temp[" + str(Offset_MxXv) + " + 6*" + jid + " + row];"])
        # check for repeated parent and add atomics
        if self.robot.has_repeated_parents(inds):
            self.gen_add_code_line("// Atomics required for shared parent")
            self.gen_add_code_line("atomicAdd(dst,update_val);")
        else:
            self.gen_add_code_line("*dst += update_val;")
        self.gen_add_end_control_flow()
        self.gen_add_sync()

        if self.DEBUG_MODE:
            self.gen_add_sync()
            self.gen_add_serial_ops()
            if bfs_level == 0:
                self.gen_add_code_lines(["printf(\"-------------------------\\n\");", \
                                         "printf(\"df/du in bfs waves\\n\");", \
                                         "printf(\"-------------------------\\n\");"])
            self.gen_add_code_line("printf(\"df/du for bfs wave[%d]\\n\"," + str(bfs_level) + ");")
            for ind in self.robot.get_unique_parent_ids(inds):
                self.gen_add_code_lines(["printf(\"df[%d]/dq (parent update)\\n\"," + str(ind) + ");", \
                                         "printMat<T,6," + str(df_cols_per_jid[ind]) + ">(&s_temp[" + str(Offset_df_dq + 6*running_sum_df_cols_per_jid[ind]) + "],6);", \
                                         "printf(\"df[%d]/dqd (parent update)\\n\"," + str(ind) + ");", \
                                         "printMat<T,6," + str(df_cols_per_jid[ind]) + ">(&s_temp[" + str(Offset_df_dqd + 6*running_sum_df_cols_per_jid[ind]) + "],6);"])
            self.gen_add_end_control_flow()
            self.gen_add_sync()

    if self.DEBUG_MODE:
        self.gen_add_sync()
        self.gen_add_serial_ops()
        self.gen_add_code_lines(["printf(\"-------------------------\\n\");", \
                                 "printf(\"Final dvaf/du\\n\");", \
                                 "printf(\"-------------------------\\n\");"])
        for ind in range(n):
            num_cols = self.robot.get_bfs_level_by_id(ind) + 1
            self.gen_add_code_lines(["printf(\"dv[%d]/dq\\n\"," + str(ind) + ");", \
                                     "printMat<T,6," + str(num_cols) + ">(&s_temp[" + str(Offset_dv_dq + 6*running_sum_dva_cols_per_jid[ind]) + "],6);"])
        for ind in range(n):
            num_cols = self.robot.get_bfs_level_by_id(ind) + 1
            self.gen_add_code_lines(["printf(\"dv[%d]/dqd\\n\"," + str(ind) + ");", \
                                     "printMat<T,6," + str(num_cols) + ">(&s_temp[" + str(Offset_dv_dqd + 6*running_sum_dva_cols_per_jid[ind]) + "],6);"])
        for ind in range(n):
            num_cols = self.robot.get_bfs_level_by_id(ind) + 1
            self.gen_add_code_lines(["printf(\"da[%d]/dq\\n\"," + str(ind) + ");", \
                                     "printMat<T,6," + str(num_cols) + ">(&s_temp[" + str(Offset_da_dq + 6*running_sum_dva_cols_per_jid[ind]) + "],6);"])
        for ind in range(n):
            num_cols = self.robot.get_bfs_level_by_id(ind) + 1
            self.gen_add_code_lines(["printf(\"da[%d]/dqd\\n\"," + str(ind) + ");", \
                                     "printMat<T,6," + str(num_cols) + ">(&s_temp[" + str(Offset_da_dqd + 6*running_sum_dva_cols_per_jid[ind]) + "],6);"])
        for ind in range(n):
            num_cols = self.robot.get_bfs_level_by_id(ind) + 1
            self.gen_add_code_lines(["printf(\"df[%d]/dq\\n\"," + str(ind) + ");", \
                                     "printMat<T,6," + str(num_cols) + ">(&s_temp[" + str(Offset_df_dq + 6*running_sum_df_cols_per_jid[ind]) + "],6);"])
        for ind in range(n):
            num_cols = self.robot.get_bfs_level_by_id(ind) + 1
            self.gen_add_code_lines(["printf(\"df[%d]/dqd\\n\"," + str(ind) + ");", \
                                     "printMat<T,6," + str(num_cols) + ">(&s_temp[" + str(Offset_df_dqd + 6*running_sum_df_cols_per_jid[ind]) + "],6);"])
        self.gen_add_end_control_flow()

    # extract dc/du
    self.gen_add_code_line("// Finally dc[i]/du = S[i]^T*df[i]/du")
    _, S_ind_cpp, _, df_col_offset_for_jid_cpp, _, _, _, _ = self.gen_topology_helpers_pointers_for_cpp(list(range(n)), OFFSET=False)
    S_sign_cpp = self.gen_topology_S_sign_for_cpp(OFFSET=False)
    # Note that for a serial chain this is straightforward (all df are size n) but otherwise gets complicated
    if self.robot.is_serial_chain() or self.robot.floating_base:
        self.gen_add_parallel_loop("ind",str(2*n*n))
        if self.robot.floating_base: 
            self.gen_add_code_line(f"bool dq_flag = ind < {n*n}; int row = ind % {n}; int col = (ind / {n}) % {n};")
            self.gen_add_code_line("int jid = row < 6 ? 0 : row - 5;")
            if 'jid' in S_ind_cpp: S_ind_cpp = S_ind_cpp.replace('jid', 'row')
            if 'jid' in S_sign_cpp: S_sign_cpp = S_sign_cpp.replace('jid', 'row')
            self.gen_add_code_line(f"int srcOffset = dq_flag * {Offset_df_dq} + !dq_flag * {Offset_df_dqd} + {6*n}*jid + 6*col + {S_ind_cpp};")
            self.gen_add_code_line(f"s_dc_du[!dq_flag * {n*n} + {n}*col + row] = (" + S_sign_cpp + ") * s_temp[srcOffset];")
        else:
            self.gen_add_code_line("int jid = ind % " + str(n) + "; int jid_dq_qd = ind / " + str(n) + "; " + 
                                "int jid_du = jid_dq_qd % " + str(n) + "; int dq_flag = jid_du == jid_dq_qd;")
            self.gen_add_code_lines(["int Offset_src = dq_flag * " + str(Offset_df_dq) + " + !dq_flag * " + str(Offset_df_dqd) + \
                                        " + 6 * " + str(n) + " * jid + 6 * jid_du + " + S_ind_cpp + ";",
                                    "int Offset_dst = !dq_flag * " + str(n*n) + " + " + str(n) + " * jid_du + jid;"])
            self.gen_add_code_line("s_dc_du[Offset_dst] = (" + S_sign_cpp + ") * s_temp[Offset_src];")
        self.gen_add_end_control_flow()
        self.gen_add_sync()

    else:
        self.gen_add_parallel_loop("jid_dq_qd",str(2*NJ))
        self.gen_add_code_line("int jid = jid_dq_qd % " + str(NJ) + "; int dq_flag = jid == jid_dq_qd;")
        # now we need to get a local pointer and loop over all n filling in 0 or col data based on the local pointer
        # and the specific topology of the robot
        self.gen_add_code_line("// Note that this gets a tad complicated due to memory compression and variable column length")
        self.gen_add_code_line("//    so we need to fully unroll the loop -- this will not be the most efficient for a serial")
        self.gen_add_code_line("//    chain manipulator but will generalize to branched robots")
        self.gen_add_code_lines(["int Offset_src = dq_flag * " + str(Offset_df_dq) + " + !dq_flag * " + str(Offset_df_dqd) + \
                                                 " + 6*" + df_col_offset_for_jid_cpp + " + " + S_ind_cpp + ";",
                                 "int Offset_dst = !dq_flag * " + str(n*n) + " + jid; bool flag = 0;"])
        for djid in range(NJ):
            self.gen_add_code_line("// dc[jid]/du[" + str(djid) + "]")
            # extract all the inds we care about for this du
            is_in_subtree_or_ancestor = [self.robot.get_is_in_subtree_of(djid,ind) or self.robot.get_is_ancestor_of(djid,ind) for ind in range(NJ)]
            non_zero_inds = [i for (i, x) in enumerate(is_in_subtree_or_ancestor) if x == True]
            # then set the if statement (if applicable)
            if len(non_zero_inds) != n:
                zero_inds = list(set(list(range(n))).difference(set(non_zero_inds)))
                if len(non_zero_inds) == 0:
                    # No body in [0,n) has djid in its subtree/ancestor set: this
                    # du-column contributes nothing. gen_var_in_list([]) would emit
                    # an empty "()" expression, so hardcode flag=false.
                    jid_du_check_for_jid = "false"
                elif len(zero_inds) == 0:
                    # Every body in [0,n) couples to djid (non_zero_inds covers the
                    # whole reduced range, e.g. when a mimic body extends the raw
                    # body set past NV): gen_var_not_in_list([]) would emit "()";
                    # the flag is unconditionally true.
                    jid_du_check_for_jid = "true"
                elif len(non_zero_inds) > n/2:
                    jid_du_check_for_jid = self.gen_var_not_in_list("jid",[str(i) for i in zero_inds])
                else:
                    jid_du_check_for_jid = self.gen_var_in_list("jid",[str(i) for i in non_zero_inds])
                self.gen_add_code_line("flag = " + jid_du_check_for_jid + ";")
                # compute the val and pointer updates accordingly
                self.gen_add_code_line("s_dc_du[Offset_dst] = flag * (" + S_sign_cpp + ") * s_temp[Offset_src]; Offset_src += flag*6; Offset_dst += " + str(n) + ";")
            else: # else everyone updates and updates their pointer
                self.gen_add_code_line("s_dc_du[Offset_dst] = (" + S_sign_cpp + ") * s_temp[Offset_src]; Offset_src += 6; Offset_dst += " + str(n) + ";")
        self.gen_add_end_control_flow()
        self.gen_add_sync()

    if self.DEBUG_MODE:
        self.gen_add_sync()
        self.gen_add_serial_ops()
        self.gen_add_code_lines(["printf(\"-------------------------\\n\");", \
                                 "printf(\"Final dc/du\\n\");", \
                                 "printf(\"-------------------------\\n\");"])
        self.gen_add_code_lines(["printf(\"dc/dq\\n\");", \
                                 "printMat<T," + str(n) + "," + str(n) + ">(&s_dc_du[0]," + str(n) + ");", \
                                 "printf(\"dc/dqd\\n\");", \
                                 "printMat<T," + str(n) + "," + str(n) + ">(&s_dc_du[" + str(n*n) + "]," + str(n) + ");"])
        self.gen_add_end_control_flow()
        self.gen_add_sync()

    self.gen_add_end_function()
    function_code = self.code_str[function_start:]
    self.code_str = self.code_str[:function_start] + _rewrite_id_du_temp_accesses_for_spill(function_code)

def gen_inverse_dynamics_gradient_device_function_call(self,
                                                           use_qdd_input = False,
                                                           scratch_in_smem_expr = "true",
                                                           use_da_df_spill_expr = "false",
                                                           d_workspace_pool_name = "nullptr",
                                                           d_temp_spill_name = "nullptr",
                                                           d_f_ext_name = "d_f_ext"):
    """Emit the call to `inverse_dynamics_gradient_device`. Arg order MUST
    match the def in gen_inverse_dynamics_gradient_device. Pool/spill regions
    default to nullptr (unused under the matching if-constexpr); the kernel passes
    real pointers per tier. The _qdd C++ name variant additionally threads s_qdd."""
    fname = "inverse_dynamics_gradient_device_qdd" if use_qdd_input else "inverse_dynamics_gradient_device"
    tmpl = "<T, " + scratch_in_smem_expr + ", " + use_da_df_spill_expr + ">"
    start = fname + tmpl + "(s_dc_du, s_q, s_qd, s_vaf, "
    if use_qdd_input:
        start += "s_qdd, "
    middle = self.gen_insert_helpers_function_call()
    end = ("s_temp, " + d_workspace_pool_name + ", " + d_temp_spill_name + ", "
           + "d_robotModel, " + d_f_ext_name + ", gravity);")
    self.gen_add_code_line(start + middle + end)

def gen_inverse_dynamics_gradient_device(self, use_qdd_input = False):
    """Emit `inverse_dynamics_gradient_device` — the whole id_du orchestration
    as ONE inner that OWNS its scratch (s_temp) placement (inner-owns-placement;
    mirrors gen_fdsva_so_device). It wraps, in order:
      [repoint s_temp] -> load_update_XImats -> inverse_dynamics_inner (vaf) ->
      inverse_dynamics_gradient_inner (the id_du band sub-inner).
    Because the s_temp repoint happens at the very top, EVERY consumer below —
    including the XImats helper's sincos scratch — follows the placement, so the
    kernel never repoints s_temp from the outside.

    TWO independent template flags:
      SCRATCH_IN_SMEM  : the shared s_temp pool lives in smem (true) or routes the
                         WHOLE pool to d_workspace (false; the rung-2 global-temp
                         path). Dominant lever on big floating humanoids.
      USE_DA_DF_SPILL  : the id_du band selectively spills its da_dq..fxvi band to
                         d_temp_spill (rung 1). Threaded through to the band
                         sub-inner's grid_id_du_temp_ptr<T, USE_DA_DF_SPILL> helper.
    The 3-rung menu (see _ID_DU_PICK_FLAGS): pick0=(SMEM=true, SPILL=false) full;
    pick1=(true, true) selective band; pick2=(false, false) whole-pool global.

    Pointer params are caller-supplied (the kernel decides where the OUTPUT s_dc_du
    lives and hands in the spill regions); only the s_temp POOL placement is the
    inner's call. The id inner `inverse_dynamics_inner_vaf` is FROZEN and
    placement-free: after the repoint, s_temp already points at the right pool, so
    passing it through is correct with no id-side change."""
    n = self.robot.get_num_vel()
    func_params = [
        "s_dc_du is the output buffer (caller places); size 2*NUM_JOINTS*NUM_JOINTS = " + str(2*n*n),
        "s_q is the vector of joint positions",
        "s_qd is the vector of joint velocities",
        "s_vaf is the id intermediate band (caller places); size 18*NUM_JOINTS = " + str(18*n),
        "s_temp is the shared scratch pool (used when SCRATCH_IN_SMEM)",
        "d_workspace is the global scratch pool (used when !SCRATCH_IN_SMEM)",
        "d_temp_spill is the id_du da_df band spill region (used when USE_DA_DF_SPILL)",
        "d_robotModel holds XImats/topology; gravity is the gravity constant",
    ]
    fname = "inverse_dynamics_gradient_device_qdd" if use_qdd_input else "inverse_dynamics_gradient_device"
    func_def_start = "void " + fname + "(T *s_dc_du, const T *s_q, const T *s_qd, T *s_vaf, "
    if use_qdd_input:
        func_def_start += "const T *s_qdd, "
    func_def_end = ("T *s_temp, T *d_workspace, T *d_temp_spill, "
                    "const robotModel<T> *d_robotModel, T *d_f_ext, const T gravity) {")
    func_params.append("d_f_ext is the (optional) GLOBAL external forces, body-major 6*NUM_BODIES local-frame, or nullptr")
    if use_qdd_input:
        func_params.insert(4, "s_qdd is the vector of joint accelerations")
    func_def_start, func_params = self.gen_insert_helpers_func_def_params(func_def_start, func_params, -2)
    func_def = func_def_start + func_def_end
    self.gen_add_func_doc("id_du orchestration as a single inner-owns-placement device function",
                          ["Owns the s_temp pool placement; the repoint covers every consumer below (incl. the XImats helper's sincos scratch)"],
                          func_params, None)
    self.gen_add_code_line("template <typename T, bool SCRATCH_IN_SMEM = true, bool USE_DA_DF_SPILL = false>")
    # __forceinline__ so the whole orchestration inlines into the calling kernel.
    # Under -rdc a separate __device__ wrapper keeps its callees as distinct
    # functions whose regcount must fit the kernel's launch_bounds budget
    # (80 at LITE / 64 at MINIMAL) -> ptxas regcount error. Inlining folds them
    # into the kernel. See _fdsva_so.py:295-300 / HANDOFF.md "Problem 1".
    self.gen_add_code_line("__device__ __forceinline__")
    self.gen_add_code_line(func_def, True)
    # Inner owns the pool placement; the repoint covers every consumer below
    # (incl. the XImats helper's sincos scratch), so no caller-side repoint.
    self.gen_add_code_line("if constexpr (!SCRATCH_IN_SMEM) { s_temp = d_workspace; } else { (void)d_workspace; }")
    self.gen_load_update_XImats_helpers_function_call()
    self.gen_inverse_dynamics_inner_function_call(False, use_qdd_input)
    self.gen_inverse_dynamics_gradient_inner_function_call(
        dict(d_temp_spill_name = "d_temp_spill", temp_spill_flag_name = "USE_DA_DF_SPILL")
    )
    self.gen_add_end_function()

def gen_inverse_dynamics_gradient_kernel_max_temp_mem_size(self):
    n = self.robot.get_num_vel()
    # s_vaf is 18*NB (body-indexed) for mimic robots; 18*n otherwise (non-mimic
    # floating has nv > NB so 18*n is the safe/byte-identical size).
    vaf_cnt = 18 * (self.robot.get_num_joints() if self.robot_has_mimic_joints() else n)
    base_size = 2*n + n*2*n + vaf_cnt + n
    temp_mem_size = self.gen_inverse_dynamics_gradient_inner_temp_mem_size()
    return base_size + temp_mem_size

_ID_DU_PICK_FLAGS = [
    # (use_selective_spill, use_global_temp)
    (False, False),   # pick 0: full smem
    (True,  False),   # pick 1: selective spill (da_df band to workspace)
    (False, True),    # pick 2: global temp (entire s_temp to workspace)
]

def _emit_id_du_kernel_body_for_flags(self, NUM_POS, n, use_selective_spill, use_global_temp,
                                      use_qdd_input, single_call_timing):
    """Emit the id_du kernel body for one tier's spill flags."""
    # s_vaf is body-indexed (the ID inner writes NB bodies, stride 6). For a
    # MIMIC robot (fixed base) NB > nv, so size it 18*NB to avoid overflowing
    # into the adjacent arena buffers. Non-mimic keeps 18*n (byte-identical;
    # for floating non-mimic nv > NB so 18*n already covers the body writes).
    _vaf_cnt = 18 * (self.robot.get_num_joints() if self.robot_has_mimic_joints() else n)
    extra_t_buffers = [("s_q_qd", n + NUM_POS), ("s_dc_du", n*2*n), ("s_vaf", _vaf_cnt)]
    if use_qdd_input:
        extra_t_buffers.append(("s_qdd", n))
    # Mimic robots use a dense inner with NO sparse-band selective spill, so the
    # selective-spill smem size collapses to the dense full size (the inner
    # ignores d_temp_spill and reads the whole dense pool from s_temp/workspace).
    _selective_shared = (
        self.gen_inverse_dynamics_gradient_inner_temp_mem_size()
        if self.robot_has_mimic_joints()
        else self.gen_inverse_dynamics_gradient_temp_layout()["selective_shared_count"]
    )
    shared_mem_size = 0 if use_global_temp else (
        _selective_shared
        if use_selective_spill else self.gen_inverse_dynamics_gradient_inner_temp_mem_size()
    )
    self.gen_XImats_helpers_temp_shared_memory_code(shared_mem_size, extra_t_buffers = extra_t_buffers, include_linalg_scratch=True)
    self.gen_add_code_line("T *d_temp_spill = nullptr; (void)d_temp_spill;")
    self.gen_add_code_line("T *s_q = s_q_qd; T *s_qd = &s_q_qd[" + str(NUM_POS) + "];")
    if not single_call_timing:
        self.gen_add_parallel_loop("k","NUM_TIMESTEPS",block_level = True)
        if use_qdd_input:
            self.gen_kernel_load_inputs("q_qd",str(n + NUM_POS),"qdd",str(n),stride="stride_q_qd",stride2=str(n))
        else:
            self.gen_kernel_load_inputs("q_qd",str(n + NUM_POS),stride="stride_q_qd")
        # The kernel only SLICES the workspace band pointers; the device owns
        # the s_temp pool placement (the whole-pool global-temp repoint is its
        # SCRATCH_IN_SMEM=false path). Per-rung flags are passed as literals.
        if use_selective_spill:
            self.gen_add_code_line("d_temp_spill = reinterpret_cast<T *>(&d_workspace[k*GRID_WORKSPACE_BYTES_PER_TIMESTEP<T>()]);")
        self.gen_add_code_line("// compute — the orchestration inner owns its s_temp pool placement")
        self.gen_inverse_dynamics_gradient_device_function_call(
            use_qdd_input,
            scratch_in_smem_expr = ("false" if use_global_temp else "true"),
            use_da_df_spill_expr = ("true" if use_selective_spill else "false"),
            d_workspace_pool_name = ("reinterpret_cast<T *>(&d_workspace[k*GRID_WORKSPACE_BYTES_PER_TIMESTEP<T>()])" if use_global_temp else "nullptr"),
            d_temp_spill_name = ("d_temp_spill" if use_selective_spill else "nullptr"))
        self.gen_add_sync()
        self.gen_kernel_save_result("dc_du",str(n*2*n),stride=str(n*2*n))
        self.gen_add_end_control_flow()
    else:
        if use_qdd_input:
            self.gen_kernel_load_inputs("q_qd",str(n + NUM_POS),"qdd",str(n))
        else:
            self.gen_kernel_load_inputs("q_qd",str(n + NUM_POS))
        if use_selective_spill:
            self.gen_add_code_line("d_temp_spill = reinterpret_cast<T *>(d_workspace);")
        self.gen_add_code_line("// compute with NUM_TIMESTEPS as NUM_REPS for timing")
        self.gen_add_code_line("for (int rep = 0; rep < NUM_TIMESTEPS; rep++){", True)
        if use_qdd_input:
            self.gen_anti_licm_input_reload("q_qd",str(n + NUM_POS),"qdd",str(n),feedback_from="dc_du")
        else:
            self.gen_anti_licm_input_reload("q_qd",str(n + NUM_POS),feedback_from="dc_du")
        # device owns s_temp placement (whole-pool global path = SCRATCH_IN_SMEM=false).
        self.gen_inverse_dynamics_gradient_device_function_call(
            use_qdd_input,
            scratch_in_smem_expr = ("false" if use_global_temp else "true"),
            use_da_df_spill_expr = ("true" if use_selective_spill else "false"),
            d_workspace_pool_name = ("reinterpret_cast<T *>(d_workspace)" if use_global_temp else "nullptr"),
            d_temp_spill_name = ("d_temp_spill" if use_selective_spill else "nullptr"))
        self.gen_anti_licm_output_write("dc_du")
        self.gen_add_end_control_flow()
        self.gen_kernel_save_result("dc_du",str(n*2*n))


def gen_inverse_dynamics_gradient_kernel(self, use_qdd_input = False, single_call_timing = False):
    NUM_POS = self.robot.get_num_pos()
    n = self.robot.get_num_vel()
    func_params = ["d_dc_du is a pointer to memory for the final result of size 2*NUM_JOINTS*NUM_JOINTS = " + str(2*n*n), \
                   "d_q_dq is the vector of joint positions and velocities", \
                   "stride_q_qd is the stide between each q, qd", \
                   "d_robotModel is the pointer to the initialized model specific helpers on the GPU (XImats, topology_helpers, etc.)", \
                   "d_f_ext is the (optional) GLOBAL external forces, body-major 6*NUM_BODIES local-frame, or nullptr", \
                   "gravity is the gravity constant", \
                   "num_timesteps is the length of the trajectory points we need to compute over (or overloaded as test_iters for timing)"]
    func_notes = []
    func_def_start = "void inverse_dynamics_gradient_kernel(T *d_dc_du, unsigned char *d_workspace, const T *d_q_qd, const int stride_q_qd, "
    func_def_end = "T *d_f_ext, const robotModel<T> *d_robotModel, const T gravity, const int NUM_TIMESTEPS) {"
    if use_qdd_input:
        func_def_start += "const T *d_qdd, "
        func_params.insert(-2,"d_qdd is the vector of joint accelerations")
    else:
        func_notes.append("optimized for qdd = 0")
    func_def = func_def_start + func_def_end
    if single_call_timing:
        func_def = func_def.replace("kernel(", "kernel_single_timing(")
    self.gen_add_func_doc("Computes the gradient of inverse dynamics",func_notes,func_params,None)
    self.gen_add_code_line("template <typename T, int RESOURCE_TIER = GRID_DEFAULT_RESOURCE_TIER>")
    self.gen_add_code_line("__global__")
    self.gen_add_code_line("__launch_bounds__(tier_max_threads<RESOURCE_TIER>())")
    self.gen_add_code_line(func_def, True)
    # Tier dispatch: collapsed picks → single body (current behavior);
    # divergent picks → 3 if-constexpr branches.
    picks = getattr(self, "id_du_spill_tier_3way", (0, 0, 0))
    def _emit_id_du_body(pick):
        uss, ugt = _ID_DU_PICK_FLAGS[pick]
        _emit_id_du_kernel_body_for_flags(self, NUM_POS, n, uss, ugt, use_qdd_input, single_call_timing)
    self.gen_tier_dispatch(picks, _emit_id_du_body)
    self.gen_add_end_function()

def gen_inverse_dynamics_gradient_host(self, mode = 0):
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
    func_def_start = "void inverse_dynamics_gradient(gridData<T, KIND> *hd_data, const robotModel<T> *d_robotModel, const T gravity, const int num_timesteps,"
    func_def_end =   "                               const dim3 block_dimms, const dim3 thread_dimms, cudaStream_t *streams) {"
    if single_call_timing:
        func_def_start = func_def_start.replace("(", "_single_timing(")
        func_def_end = "              " + func_def_end
    if compute_only:
        func_def_start = func_def_start.replace("(", "_compute_only(")
        func_def_end = "             " + func_def_end.replace(", cudaStream_t *streams", "")
    # then generate the code
    self.gen_add_func_doc("Compute the RNEA (Recursive Newton-Euler Algorithm)",\
                          func_notes,func_params,None)
    self.gen_add_code_line("template <typename T, bool USE_QDD_FLAG = false, bool USE_COMPRESSED_MEM = false, gridDataKind KIND = GRID_DATA_ALL>")
    self.gen_add_code_line("__host__")
    self.gen_add_code_line(func_def_start)
    self.gen_add_code_line(func_def_end, True)
    self.gen_add_code_line("static_assert(KIND == GRID_DATA_ALL || KIND == GRID_DATA_DYNAMICS, \"inverse_dynamics_gradient requires all-data or dynamics gridData\");")
    func_call_start = "inverse_dynamics_gradient_kernel<T><<<block_dimms,thread_dimms,ID_DU_DYNAMIC_SHARED_MEM_BYTES<T>()>>>(hd_data->d_dc_du,hd_data->d_workspace,hd_data->d_q_qd,stride_q_qd,"
    func_call_end = "hd_data->d_f_ext,d_robotModel,gravity,num_timesteps);"
    if single_call_timing:
        func_call_start = func_call_start.replace("kernel<T>","kernel_single_timing<T>")
    if not compute_only:
        # start code with memory transfer
        self.gen_add_code_lines(["// start code with memory transfer", \
                                 "int stride_q_qd;", \
                                 "if (USE_COMPRESSED_MEM) {stride_q_qd = 2*NUM_JOINTS; " + \
                                    "gpuErrchk(cudaMemcpyAsync(hd_data->d_q_qd,hd_data->h_q_qd,stride_q_qd*" + \
                                    ("num_timesteps*" if not single_call_timing else "") + "sizeof(T),cudaMemcpyHostToDevice,streams[0]));}", \
                                 "else {stride_q_qd = 3*NUM_JOINTS; " + \
                                    "gpuErrchk(cudaMemcpyAsync(hd_data->d_q_qd_u,hd_data->h_q_qd_u,stride_q_qd*" + \
                                    ("num_timesteps*" if not single_call_timing else "") + "sizeof(T),cudaMemcpyHostToDevice,streams[0]));}", \
                                 "if (USE_QDD_FLAG) {gpuErrchk(cudaMemcpyAsync(hd_data->d_qdd,hd_data->h_qdd,NUM_JOINTS*" + \
                                    ("num_timesteps*" if not single_call_timing else "") + "sizeof(T),cudaMemcpyHostToDevice,streams[1]));}", \
                                 "gpuErrchkKernel();"])
    else:
        self.gen_add_code_line("int stride_q_qd = USE_COMPRESSED_MEM ? 2*NUM_JOINTS: 3*NUM_JOINTS;")
    # then compute but adjust for compressed mem and qdd usage
    self.gen_add_code_line("// then call the kernel")
    func_call = func_call_start + func_call_end
    func_call_with_qdd = func_call_start + "hd_data->d_qdd, " + func_call_end
    # add in compressed mem adjusts
    func_call_mem_adjust = "    if (USE_COMPRESSED_MEM) {" + func_call + "}"
    func_call_mem_adjust2 = "    else                    {" + func_call.replace("hd_data->d_q_qd","hd_data->d_q_qd_u") + "}"
    func_call_with_qdd_mem_adjust = "    if (USE_COMPRESSED_MEM) {" + func_call_with_qdd + "}"
    func_call_with_qdd_mem_adjust2 = "    else                    {" + func_call_with_qdd.replace("hd_data->d_q_qd","hd_data->d_q_qd_u") + "}"
    # compule into a set of code
    func_call_code = ["if (USE_QDD_FLAG) {", func_call_with_qdd_mem_adjust, func_call_with_qdd_mem_adjust2, "}", \
                      "else {", func_call_mem_adjust, func_call_mem_adjust2, "}", "gpuErrchkKernel();"]
    # wrap function call in timing (if needed)
    if single_call_timing:
        func_call_code.insert(0,"struct timespec start, end; clock_gettime(CLOCK_MONOTONIC,&start);")
        func_call_code.append("clock_gettime(CLOCK_MONOTONIC,&end);")
    self.gen_add_code_line("gpuErrchk(grid_check_dynamic_shared_memory_bytes(\"inverse_dynamics_gradient\", ID_DU_DYNAMIC_SHARED_MEM_BYTES<T>()));")
    workspace_bytes = "GRID_WORKSPACE_BYTES_PER_TIMESTEP<T>()" if single_call_timing else "GRID_WORKSPACE_BYTES_PER_TIMESTEP<T>()*static_cast<size_t>(num_timesteps)"
    self.gen_add_code_line("if (GRID_ID_DU_USES_WORKSPACE_ANY_TIER) {gpuErrchk(grid_begin_l2_persisting(0, hd_data->d_workspace, " + workspace_bytes + "));}")
    self.gen_add_code_lines(func_call_code)
    self.gen_add_code_line("if (GRID_ID_DU_USES_WORKSPACE_ANY_TIER) {gpuErrchk(grid_end_l2_persisting(0));}")
    if not compute_only:
        # then transfer memory back
        self.gen_add_code_lines(["// finally transfer the result back", \
                                 "gpuErrchk(cudaMemcpy(hd_data->h_dc_du,hd_data->d_dc_du,NUM_JOINTS*2*NUM_JOINTS*" + \
                                    ("num_timesteps*" if not single_call_timing else "") + "sizeof(T),cudaMemcpyDeviceToHost));",
                                 "gpuErrchkKernel();"])
    # finally report out timing if requested
    if single_call_timing:
        from ..algo_registry import single_call_printf_line
        self.gen_add_code_line(single_call_printf_line("id_du"))
    self.gen_add_end_function()

def gen_inverse_dynamics_gradient(self):
    # first the id_du band sub-inner (internal helper composed by the orchestration
    # _device; also called by fd_du / integrator_gradient orchestrators).
    self.gen_inverse_dynamics_gradient_inner()
    # then the canonical _device (orchestrator: owns s_temp placement; wraps XImats
    # + id-inner + id_du band sub-inner; called from kernel and fd_du / integrator).
    # Both qdd variants (qdd-input vs qdd=0 specialization).
    self.gen_inverse_dynamics_gradient_device(True)
    self.gen_inverse_dynamics_gradient_device(False)
    # and the kernels
    self.gen_inverse_dynamics_gradient_kernel(True,True)
    self.gen_inverse_dynamics_gradient_kernel(True,False)
    self.gen_inverse_dynamics_gradient_kernel(False,True)
    self.gen_inverse_dynamics_gradient_kernel(False,False)
    # and host wrapeprs
    self.gen_inverse_dynamics_gradient_host(0)
    self.gen_inverse_dynamics_gradient_host(1)
    self.gen_inverse_dynamics_gradient_host(2)


def _gen_id_du_mimic_inner(self, nv, NB):
    """Dense serial mimic ID-gradient (T3-finisher P3, fixed-base).

    Mirrors RBDReference.rnea_grad exactly, in REDUCED v-space:
      forward pass dq/dqd -> per-body dv_du, da_du, df_du (6 x nv x NB),
      backward pass       -> dc_dq, dc_dqd (nv x nv).
    Every joint-velocity read scales by the body's mimic multiplier alpha and
    every dc_du / df_du write accumulates (+=) into the body's reduced v-slot,
    so a mimic body and its target fold together. Single-DoF bodies only
    (fixed base); the multi-DoF floating root mimic-gradient is a separate
    follow-on and is refused upstream.

    Output s_dc_du is 2*nv*nv, column-major nv x nv per half:
      dc_dq  at [0, nv*nv)       element [v_i, c] -> s_dc_du[c*nv + v_i]
      dc_dqd at [nv*nv, 2*nv*nv)  element [v_i, c] -> s_dc_du[nv*nv + c*nv + v_i]
    s_vaf is body-indexed (stride NB): v @ s_vaf[6*ind], a @ s_vaf[6*NB+6*ind],
    f @ s_vaf[12*NB+6*ind]. The big dense buffers live in s_temp (routed to
    workspace at the global-temp tier for humanoid-scale NB)."""
    import numpy as _np
    assert not self.robot.floating_base, \
        "_gen_id_du_mimic_inner is fixed-base only (floating mimic gradient refused upstream)"
    GRAV_NEG = "gravity"  # s_a base row 5 holds X*gravity already via s_vaf

    bw = 6 * nv * NB  # one dense buffer (6 rows x nv cols x NB bodies)
    off_dv_dq  = 0
    off_da_dq  = off_dv_dq  + bw
    off_df_dq  = off_da_dq  + bw
    off_dv_dqd = off_df_dq  + bw
    off_da_dqd = off_dv_dqd + bw
    off_df_dqd = off_da_dqd + bw
    off_iv     = off_df_dqd + bw          # Iv per body: 6*NB
    off_scr    = off_iv + 6 * NB          # scratch 6-vectors (a few)

    def cell(base, ind, c):
        # &buffer[base] element column c of body ind (6-vector)
        return base + ind * (6 * nv) + 6 * c

    self.gen_add_code_line("// === mimic ID-gradient (dense serial reduced-space fold) ===")
    self.gen_add_code_line("// zero the dense fwd buffers + output")
    self.gen_add_parallel_loop("i", str(off_iv))
    self.gen_add_code_line("s_temp[i] = static_cast<T>(0);")
    self.gen_add_end_control_flow()
    self.gen_add_parallel_loop("i", str(2 * nv * nv))
    self.gen_add_code_line("s_dc_du[i] = static_cast<T>(0);")
    self.gen_add_end_control_flow()
    self.gen_add_sync()

    # ---- forward pass (serial over bodies, root first) ----
    self.gen_add_serial_ops()
    self.gen_add_code_line("T s_iv6[6];")
    self.gen_add_code_line("T s_mtmp[6];")
    self.gen_add_code_line("T s_ftmp[6];")
    self.gen_add_code_line("T s_Svec[6];")
    for ind in range(NB):
        parent = self.robot.get_parent_id(ind)
        idx = self._v_slot_cpp(ind)
        alpha = self._alpha_for_jid(ind)
        s_ind = self.robot.get_S_index_by_id(ind)
        s_sign = float(self.robot.get_S_sign_by_id(ind))
        v_ind = 6 * ind                 # s_vaf v
        a_ind = 6 * NB + 6 * ind        # s_vaf a
        Xoff = 36 * ind                 # X[ind] col-major in s_XImats
        Ioff = 36 * NB + 36 * ind       # I[ind]
        self.gen_add_code_line("// --- body " + str(ind) + " (v-slot " + str(idx) +
                               ", alpha=" + repr(alpha) + ", S_ind=" + str(s_ind) +
                               ", S_sign=" + repr(s_sign) + ") ---")
        self.gen_add_code_line("{", True)  # per-body scope (avoid local redeclare)
        # Build the S 6-vector (s_sign at row s_ind) for fxS-style products.
        self.gen_add_code_line("for (int r = 0; r < 6; r++) s_Svec[r] = static_cast<T>(0);")
        self.gen_add_code_line("s_Svec[" + str(s_ind) + "] = static_cast<T>(" + repr(s_sign) + ");")

        # Iv = I[ind] * v[ind]
        self.gen_add_code_line("for (int r = 0; r < 6; r++) { s_iv6[r] = static_cast<T>(0);")
        self.gen_add_code_line("  for (int p = 0; p < 6; p++) s_iv6[r] += s_XImats[" + str(Ioff) + " + r + 6*p] * s_vaf[" + str(v_ind) + " + p]; }")
        self.gen_add_code_line("for (int r = 0; r < 6; r++) s_temp[" + str(off_iv + 6*ind) + " + r] = s_iv6[r];")

        if parent != -1:
            p_dv_dq  = cell(off_dv_dq,  parent, 0)
            p_da_dq  = cell(off_da_dq,  parent, 0)
            p_dv_dqd = cell(off_dv_dqd, parent, 0)
            p_da_dqd = cell(off_da_dqd, parent, 0)
            c_dv_dq  = cell(off_dv_dq,  ind, 0)
            c_da_dq  = cell(off_da_dq,  ind, 0)
            c_dv_dqd = cell(off_dv_dqd, ind, 0)
            c_da_dqd = cell(off_da_dqd, ind, 0)
            # dv_du[ind] = X[ind] * dv_du[parent] ; da_du[ind] = X[ind]*da_du[parent]
            self.gen_add_code_line("for (int c = 0; c < " + str(nv) + "; c++) {", True)
            self.gen_add_code_line("for (int r = 0; r < 6; r++) {", True)
            self.gen_add_code_line("T acc_vq=static_cast<T>(0), acc_aq=static_cast<T>(0), acc_vqd=static_cast<T>(0), acc_aqd=static_cast<T>(0);")
            self.gen_add_code_line("for (int p = 0; p < 6; p++) {", True)
            self.gen_add_code_line("T xrp = s_XImats[" + str(Xoff) + " + r + 6*p];")
            self.gen_add_code_line("acc_vq  += xrp * s_temp[" + str(p_dv_dq)  + " + 6*c + p];")
            self.gen_add_code_line("acc_aq  += xrp * s_temp[" + str(p_da_dq)  + " + 6*c + p];")
            self.gen_add_code_line("acc_vqd += xrp * s_temp[" + str(p_dv_dqd) + " + 6*c + p];")
            self.gen_add_code_line("acc_aqd += xrp * s_temp[" + str(p_da_dqd) + " + 6*c + p];")
            self.gen_add_end_control_flow()
            self.gen_add_code_line("s_temp[" + str(c_dv_dq)  + " + 6*c + r] = acc_vq;")
            self.gen_add_code_line("s_temp[" + str(c_da_dq)  + " + 6*c + r] = acc_aq;")
            self.gen_add_code_line("s_temp[" + str(c_dv_dqd) + " + 6*c + r] = acc_vqd;")
            self.gen_add_code_line("s_temp[" + str(c_da_dqd) + " + 6*c + r] = acc_aqd;")
            self.gen_add_end_control_flow()
            self.gen_add_end_control_flow()

            # dv_dq[:,idx,ind]  += alpha * mxS(S, X*v_parent)  = alpha*s_sign*mx_Sind(X v_parent)
            # X*v_parent: dv contribution uses v[parent]
            self.gen_add_code_line("// dv_dq[:,idx] += alpha*mxS(S, X*v_parent); dv_dqd[:,idx] += alpha*S")
            self.gen_add_code_line("for (int r = 0; r < 6; r++) { s_mtmp[r] = static_cast<T>(0);")
            self.gen_add_code_line("  for (int p = 0; p < 6; p++) s_mtmp[r] += s_XImats[" + str(Xoff) + " + r + 6*p] * s_vaf[" + str(6*parent) + " + p]; }")
            # mx<s_ind>_peq_scaled into dv_dq[:,idx,ind]
            self.gen_add_code_line("mx" + str(s_ind) + "_peq_scaled<T>(&s_temp[" + str(cell(off_dv_dq, ind, 0)) + " + 6*" + str(idx) + "], s_mtmp, static_cast<T>(" + repr(alpha) + "));")

        # dv_dqd[:,idx,ind] += alpha*S  (S = s_sign*e_{s_ind}). NOTE: the oracle
        # adds this for EVERY body including the root (it sits OUTSIDE the
        # parent!=-1 guard in rnea_grad_fpass_dqd) — the joint's own velocity
        # subspace contributes to dv/dqd regardless of having a parent.
        self.gen_add_code_line("// dv_dqd[:,idx] += alpha*S (all bodies incl. root)")
        self.gen_add_code_line("s_temp[" + str(cell(off_dv_dqd, ind, 0)) + " + 6*" + str(idx) + " + " + str(s_ind) + "] += static_cast<T>(" + repr(alpha * s_sign) + ");")

        # da_du[:,c,ind] += mxS(S, dv_du[:,c,ind], alpha*qd[idx])   for every column c
        self.gen_add_code_line("// da_du[:,c] += mxS(S, dv_du[:,c], alpha*qd[idx])")
        self.gen_add_code_line("T qd_a = static_cast<T>(" + repr(alpha) + ") * s_qd[" + str(idx) + "];")
        self.gen_add_code_line("for (int c = 0; c < " + str(nv) + "; c++) {", True)
        self.gen_add_code_line("mx" + str(s_ind) + "_peq_scaled<T>(&s_temp[" + str(cell(off_da_dq, ind, 0)) + " + 6*c], &s_temp[" + str(cell(off_dv_dq, ind, 0)) + " + 6*c], static_cast<T>(" + repr(s_sign) + ") * qd_a);")
        self.gen_add_code_line("mx" + str(s_ind) + "_peq_scaled<T>(&s_temp[" + str(cell(off_da_dqd, ind, 0)) + " + 6*c], &s_temp[" + str(cell(off_dv_dqd, ind, 0)) + " + 6*c], static_cast<T>(" + repr(s_sign) + ") * qd_a);")
        self.gen_add_end_control_flow()

        # da_dq[:,idx,ind] += alpha*mxS(S, X*a_parent or root_gravity)
        # da_dqd[:,idx,ind] += alpha*mxS(S, v[ind])
        self.gen_add_code_line("// da_dq[:,idx] += alpha*mxS(S, X*a_parent); da_dqd[:,idx] += alpha*mxS(S, v[ind])")
        if parent != -1:
            self.gen_add_code_line("for (int r = 0; r < 6; r++) { s_mtmp[r] = static_cast<T>(0);")
            self.gen_add_code_line("  for (int p = 0; p < 6; p++) s_mtmp[r] += s_XImats[" + str(Xoff) + " + r + 6*p] * s_vaf[" + str(6*NB + 6*parent) + " + p]; }")
            self.gen_add_code_line("mx" + str(s_ind) + "_peq_scaled<T>(&s_temp[" + str(cell(off_da_dq, ind, 0)) + " + 6*" + str(idx) + "], s_mtmp, static_cast<T>(" + repr(alpha) + "));")
        else:
            # root: the base's accel is PURE gravity (NOT the body's own a, which
            # also carries S*qdd when use_qdd_input — that would corrupt fd_du).
            # X*gravity is column 5 of X scaled by `gravity`:
            #   (X*gravity)[r] = s_XImats[36*root + 30 + r] * gravity   (col 5 = +30).
            self.gen_add_code_line("for (int r = 0; r < 6; r++) s_mtmp[r] = s_XImats[" + str(Xoff) + " + 30 + r] * gravity;")
            self.gen_add_code_line("mx" + str(s_ind) + "_peq_scaled<T>(&s_temp[" + str(cell(off_da_dq, ind, 0)) + " + 6*" + str(idx) + "], s_mtmp, static_cast<T>(" + repr(alpha) + "));")
        # da_dqd[:,idx,ind] += alpha*mxS(S, v[ind])
        self.gen_add_code_line("mx" + str(s_ind) + "_peq_scaled<T>(&s_temp[" + str(cell(off_da_dqd, ind, 0)) + " + 6*" + str(idx) + "], &s_vaf[" + str(v_ind) + "], static_cast<T>(" + repr(alpha) + "));")

        # df_du[:,:,ind] = I*da_du + fxv(dv_du, Iv) + fxv(v, I*dv_du)
        self.gen_add_code_line("// df_du[:,c] = I*da_du[:,c] + fx(dv_du[:,c])*Iv + fx(v)*I*dv_du[:,c]")
        self.gen_add_code_line("for (int c = 0; c < " + str(nv) + "; c++) {", True)
        for (dabuf, dfbuf, dvbuf) in [(off_da_dq, off_df_dq, off_dv_dq), (off_da_dqd, off_df_dqd, off_dv_dqd)]:
            # I*da_du
            self.gen_add_code_line("for (int r = 0; r < 6; r++) { T acc=static_cast<T>(0);")
            self.gen_add_code_line("  for (int p = 0; p < 6; p++) acc += s_XImats[" + str(Ioff) + " + r + 6*p] * s_temp[" + str(cell(dabuf, ind, 0)) + " + 6*c + p];")
            self.gen_add_code_line("  s_temp[" + str(cell(dfbuf, ind, 0)) + " + 6*c + r] = acc; }")
            # fxv(dv_du[:,c], Iv)
            self.gen_add_code_line("fx_times_v<T>(s_ftmp, &s_temp[" + str(cell(dvbuf, ind, 0)) + " + 6*c], &s_temp[" + str(off_iv + 6*ind) + "]);")
            self.gen_add_code_line("for (int r = 0; r < 6; r++) s_temp[" + str(cell(dfbuf, ind, 0)) + " + 6*c + r] += s_ftmp[r];")
            # I*dv_du[:,c]
            self.gen_add_code_line("for (int r = 0; r < 6; r++) { T acc=static_cast<T>(0);")
            self.gen_add_code_line("  for (int p = 0; p < 6; p++) acc += s_XImats[" + str(Ioff) + " + r + 6*p] * s_temp[" + str(cell(dvbuf, ind, 0)) + " + 6*c + p];")
            self.gen_add_code_line("  s_mtmp[r] = acc; }")
            # fxv(v[ind], I*dv_du[:,c])
            self.gen_add_code_line("fx_times_v<T>(s_ftmp, &s_vaf[" + str(v_ind) + "], s_mtmp);")
            self.gen_add_code_line("for (int r = 0; r < 6; r++) s_temp[" + str(cell(dfbuf, ind, 0)) + " + 6*c + r] += s_ftmp[r];")
        self.gen_add_end_control_flow()  # end for c (df_du)
        self.gen_add_end_control_flow()  # end per-body scope
    self.gen_add_end_control_flow()  # end serial fwd
    self.gen_add_sync()

    # ---- backward pass (serial, deepest first) ----
    # dc_du[idx,:] += alpha * S^T * df_du[:,:,ind]
    # df_du[:,idx,parent] += alpha * (X^T * fxS(S, f[ind]))   [dq only]
    # df_du[:,:,parent]  += X^T * df_du[:,:,ind]
    self.gen_add_serial_ops()
    self.gen_add_code_line("T s_fxs[6];")
    self.gen_add_code_line("T s_xtfxs[6];")
    self.gen_add_code_line("T s_Svec[6];")
    for ind in range(NB - 1, -1, -1):
        parent = self.robot.get_parent_id(ind)
        idx = self._v_slot_cpp(ind)
        alpha = self._alpha_for_jid(ind)
        s_ind = self.robot.get_S_index_by_id(ind)
        s_sign = float(self.robot.get_S_sign_by_id(ind))
        Xoff = 36 * ind
        f_ind = 12 * NB + 6 * ind
        self.gen_add_code_line("// --- bpass body " + str(ind) + " (v-slot " + str(idx) + ") ---")
        # dc_dq[idx, c]  += alpha * s_sign * df_dq[s_ind, c, ind]
        # dc_dqd[idx, c] += alpha * s_sign * df_dqd[s_ind, c, ind]
        coeff = alpha * s_sign
        self.gen_add_code_line("for (int c = 0; c < " + str(nv) + "; c++) {", True)
        self.gen_add_code_line("s_dc_du[c*" + str(nv) + " + " + str(idx) + "] += static_cast<T>(" + repr(coeff) + ") * s_temp[" + str(cell(off_df_dq, ind, 0)) + " + 6*c + " + str(s_ind) + "];")
        self.gen_add_code_line("s_dc_du[" + str(nv*nv) + " + c*" + str(nv) + " + " + str(idx) + "] += static_cast<T>(" + repr(coeff) + ") * s_temp[" + str(cell(off_df_dqd, ind, 0)) + " + 6*c + " + str(s_ind) + "];")
        self.gen_add_end_control_flow()
        if parent != -1:
            # df_dq[:,idx,parent] += alpha * X^T * fxS(S, f[ind])
            # fxS(S, f) = Fx(S)*f = fx_times_v(S, f); S = s_sign*e_{s_ind}
            self.gen_add_code_line("for (int r = 0; r < 6; r++) s_Svec[r] = static_cast<T>(0);")
            self.gen_add_code_line("s_Svec[" + str(s_ind) + "] = static_cast<T>(" + repr(s_sign) + ");")
            self.gen_add_code_line("fx_times_v<T>(s_fxs, s_Svec, &s_vaf[" + str(f_ind) + "]);")
            # X^T * s_fxs : (X^T)[r,p] = X[p,r] = s_XImats[Xoff + p + 6*r]
            self.gen_add_code_line("for (int r = 0; r < 6; r++) { s_xtfxs[r] = static_cast<T>(0);")
            self.gen_add_code_line("  for (int p = 0; p < 6; p++) s_xtfxs[r] += s_XImats[" + str(Xoff) + " + p + 6*r] * s_fxs[p]; }")
            self.gen_add_code_line("for (int r = 0; r < 6; r++) s_temp[" + str(cell(off_df_dq, parent, 0)) + " + 6*" + str(idx) + " + r] += static_cast<T>(" + repr(alpha) + ") * s_xtfxs[r];")
            # df_du[:,:,parent] += X^T * df_du[:,:,ind]   (both dq and dqd)
            self.gen_add_code_line("for (int c = 0; c < " + str(nv) + "; c++) {", True)
            self.gen_add_code_line("for (int r = 0; r < 6; r++) {", True)
            self.gen_add_code_line("T acc_q=static_cast<T>(0), acc_qd=static_cast<T>(0);")
            self.gen_add_code_line("for (int p = 0; p < 6; p++) {", True)
            self.gen_add_code_line("T xtr = s_XImats[" + str(Xoff) + " + p + 6*r];")
            self.gen_add_code_line("acc_q  += xtr * s_temp[" + str(cell(off_df_dq,  ind, 0)) + " + 6*c + p];")
            self.gen_add_code_line("acc_qd += xtr * s_temp[" + str(cell(off_df_dqd, ind, 0)) + " + 6*c + p];")
            self.gen_add_end_control_flow()
            self.gen_add_code_line("s_temp[" + str(cell(off_df_dq,  parent, 0)) + " + 6*c + r] += acc_q;")
            self.gen_add_code_line("s_temp[" + str(cell(off_df_dqd, parent, 0)) + " + 6*c + r] += acc_qd;")
            self.gen_add_end_control_flow()
            self.gen_add_end_control_flow()
    self.gen_add_end_control_flow()  # end serial bpass
    self.gen_add_sync()


def _id_du_mimic_temp_count(self):
    """Dense mimic ID-gradient scratch size: 6 buffers of 6*nv*NB + Iv(6*NB)."""
    nv = self.robot.get_num_vel()
    NB = self.robot.get_num_joints()
    return 6 * (6 * nv * NB) + 6 * NB
