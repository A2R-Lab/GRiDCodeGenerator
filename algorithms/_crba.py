import numpy as np
import copy
#np.set_printoptions(precision=4, suppress=True, linewidth = 100)

def gen_crba_inner_temp_mem_size(self):
    if self.robot.floating_base:
        NJ = self.robot.get_num_joints()
        # IC (36*NJ) + alpha slab (36 * max BFS width) — the alpha slab is
        # per-sibling so all gemms at a BFS level can run with one shared
        # forward/backward sync (BFS-parallel body recursion). On a 1-wide
        # level the slab is exactly 36 (same as the old single-alpha buffer).
        max_bfs_width = max(1, self.robot.get_max_bfs_width())
        return 36*NJ + 36*max_bfs_width
    n = self.robot.get_num_pos()
    return 140*n

def gen_crba_inner_function_call(self, updated_var_names = None,
                                 temp_in_smem_expr = "true"):
    var_names = dict( \
        s_M_name = "s_M", \
        s_q_name = "s_q", \
        s_qd_name = "s_qd", \
        s_temp_name = "s_temp", \
        d_workspace_name = "nullptr", \
        #s_XI = "s_XImats", \
        gravity_name = "gravity"
    )
    if updated_var_names is not None:
        for key,value in updated_var_names.items():
            var_names[key] = value
    crba_code_start = "crba_inner<T, " + temp_in_smem_expr + ">(" + var_names["s_M_name"] + ", " +  var_names["s_q_name"] + ", " + var_names["s_qd_name"] + ", "
    crba_code_end = var_names["s_temp_name"] + ", " + var_names["d_workspace_name"] + ", " + var_names["gravity_name"] + ");"
    crba_code_middle = self.gen_insert_helpers_function_call()
    crba_code = crba_code_start + crba_code_middle + crba_code_end
    self.gen_add_code_line(crba_code)


def gen_crba_inner(self):
    if self.robot.floating_base:
        return gen_crba_inner_floating(self)
    
    n = self.robot.get_num_joints()
    n_bfs_levels = self.robot.get_max_bfs_level() + 1
    has_linear_axis = any(self.robot.get_S_index_by_id(jid) >= 3 for jid in range(n))
    imat_offset = n if has_linear_axis else 7

    #construct the boilerplate and function definition
    func_params = [ "s_q is the vector of joint positions", \
                    "s_qd is the vector of joint velocities", \
                    "s_M is a pointer to the matrix of inertia" \
                    "s_XI is the pointer to the transformation and inertia matricies ", \
                    "s_temp is the (shared) scratch; size CRBA_INNER_SMEM_BYTES<T, TEMP_IN_SMEM>() (the 140*NJ band when TEMP_IN_SMEM, else 0)", \
                    "d_workspace is the global scratch; size CRBA_INNER_WORKSPACE_BYTES<T, TEMP_IN_SMEM>() (the band when !TEMP_IN_SMEM, else 0). Pass nullptr when TEMP_IN_SMEM", \
                    "gravity is the gravity constant"]
    func_notes = ["Inner-controlled placement: TEMP_IN_SMEM selects where the scratch band lives (s_temp vs d_workspace). Decided at the top; caller sizes both arenas from CRBA_INNER_*_BYTES."]
    func_def_start = "void crba_inner("
    func_def_middle = "T *s_M, const T *s_q, const T *s_qd, "
    func_def_end = "T *s_temp, T *d_workspace, const T gravity) {"

    func_def_middle, func_params = self.gen_insert_helpers_func_def_params(func_def_middle, func_params, -2)
    func_def = func_def_start + func_def_middle + func_def_end
    self.gen_add_func_doc("Compute the Composite Rigid Body Algorithm", func_notes, func_params, None)
    self.gen_add_code_line("template <typename T, bool TEMP_IN_SMEM = true>")
    self.gen_add_code_line("__device__")
    self.gen_add_code_line(func_def, True)
    # Inner-controlled scratch-band placement: the whole band moves to
    # d_workspace when !TEMP_IN_SMEM. Reassigning s_temp at the top keeps every
    # s_temp[...] reference below unchanged.
    self.gen_add_code_line("if constexpr (!TEMP_IN_SMEM) { s_temp = d_workspace; } else { (void)d_workspace; }")
    temp_size = self.gen_crba_inner_temp_mem_size()
    self.gen_linalg_smem_setup(temp_size)


    # first clear the matrix
    self.gen_add_parallel_loop("i",str(n*n))
    self.gen_add_code_line('s_M[i] = static_cast<T>(0);')
    self.gen_add_end_control_flow()
    self.gen_add_sync()
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
    
    self.gen_add_parallel_loop("jid",str(n))

    ImatOffset = 36*n   # Offset in XImats to Imats
    self.gen_add_code_line(f"s_M[jid+jid*{n}] = s_XImats[{ImatOffset} + 36*jid + 6*{S_ind_cpp} + {S_ind_cpp}];") # take the S_ind row and S_ind column of appropriate Imat
    self.gen_add_end_control_flow()
    self.gen_add_sync()
    

    self.gen_add_code_line("//")
    self.gen_add_code_line("// Calculation of M[ind, parent]")
    self.gen_add_code_line("//")

    # initialize fh as (XS)^T
    self.gen_add_parallel_loop('i',str(n*6))
    self.gen_add_code_line('int jid = i / 6; int ind = i % 6;')
    self.gen_add_code_line(f's_fh[i] = ({S_sign_cpp}) * s_XImats[{ImatOffset} + 36*jid + 6*{S_ind_cpp} + ind];')
    self.gen_add_end_control_flow()
    self.gen_add_sync()

    # M[jid, parent] = S_parent^T * (X_lambda^T chain) * IS.
    #
    # Each thread owns ONE jid and walks ITS ancestor chain serially, advancing
    # s_fh[jid] by one Xmat^T per step (the loop-carried dependence is ALONG the
    # chain). Threads are independent — thread `jid` only touches s_fh[jid*6..]
    # and the M cells (jid,parent)/(parent,jid) for its own ancestors — so the
    # whole fill needs NO inter-thread syncs. (A depth-stepped variant that split
    # each chain step into separate parallel loops added 3 __syncthreads PER
    # DEPTH — ~21 for a 7-deep chain — and regressed crba 2.5-6x; reverted.)
    #
    # CORRECTNESS: the M entry indexes s_fh by the PARENT's S index/sign, not the
    # owning jid's (they differ on branched robots). s_Sidx/s_Ssgn_by_jid are
    # compile-time per-joint tables looked up at runtime by parent id.
    max_ancestors = self.robot.get_max_num_ancestors()
    S_idx_arr = "{" + ", ".join(str(self.robot.get_S_index_by_id(j)) for j in range(n)) + "}"
    S_sgn_arr = "{" + ", ".join(str(self.robot.get_S_sign_by_id(j)) for j in range(n)) + "}"
    self.gen_add_parallel_loop("jid", str(n))
    self.gen_add_code_line(f"const int s_Sidx_by_jid[{n}] = {S_idx_arr};")
    self.gen_add_code_line(f"const T s_Ssgn_by_jid[{n}] = {S_sgn_arr};")
    parent_chain_init = "{" + "-1, " * (max_ancestors - 1) + "-1}" if max_ancestors >= 1 else "{-1}"
    self.gen_add_code_line(f"int jid_parents[] = {parent_chain_init};")
    self.gen_add_code_line("int num_parents = 0;")
    self.gen_add_code_line("switch (jid) {", True)
    for jid in range(n):
        self.gen_add_code_line(f"case {jid}:", True)
        parent_chain = self.robot.get_ancestors_by_id(jid)
        for i, parent_ind in enumerate(parent_chain):
            self.gen_add_code_line(f"jid_parents[{i}] = {parent_ind};")
        self.gen_add_code_line(f"num_parents += {len(parent_chain)};")
        self.gen_add_code_line("break;")
        self.indent_level -= 1
    self.gen_add_end_control_flow()
    self.gen_add_code_line("T s_alpha[6];")
    self.gen_add_code_line("for (int i = 0; i < num_parents; i++) {", True)
    self.gen_add_code_line("int X_ind = i==0 ? jid : jid_parents[i-1];")
    self.gen_add_code_line("for (int k = 0; k < 6; k++) s_alpha[k] = s_fh[jid*6+k];")
    self.gen_add_code_line("for (int k = 0; k < 6; k++) s_fh[jid*6 + k] = dot_prod<T,6,1,1>(&s_XImats[36*X_ind+k*6], &s_alpha[0]);")
    self.gen_add_code_line("int parent_ind = jid_parents[i];")
    self.gen_add_code_line(f"s_M[jid*{n} + parent_ind] = s_Ssgn_by_jid[parent_ind] * s_fh[jid*6 + s_Sidx_by_jid[parent_ind]];")
    self.gen_add_code_line(f"s_M[parent_ind*{n} + jid] = s_M[jid*{n} + parent_ind];") # M symmetric
    self.gen_add_end_control_flow()
    self.gen_add_end_control_flow()

    self.gen_add_end_function()


def gen_crba_inner_floating(self):
    NJ = self.robot.get_num_joints()
    nv = self.robot.get_num_vel()
    n_bfs_levels = self.robot.get_max_bfs_level() + 1
    ICOffset = 0
    # Per-sibling alpha slab (one 6x6 block per joint at the current BFS level).
    # Slot i in [0, len(inds)) lives at alphaOffset + 36*i.
    alphaOffset = 36 * NJ

    func_params = [ "s_q is the vector of joint positions", \
                    "s_qd is the vector of joint velocities", \
                    "s_M is a pointer to the matrix of inertia", \
                    "s_XI is the pointer to the transformation and inertia matricies", \
                    "s_temp is the (shared) scratch; size CRBA_INNER_SMEM_BYTES<T, TEMP_IN_SMEM>() (the band when TEMP_IN_SMEM, else 0)", \
                    "d_workspace is the global scratch; size CRBA_INNER_WORKSPACE_BYTES<T, TEMP_IN_SMEM>() (the band when !TEMP_IN_SMEM, else 0). Pass nullptr when TEMP_IN_SMEM", \
                    "gravity is the gravity constant"]
    func_notes = ["Floating-base CRBA keeps composite inertias in body order and writes a public-order NUM_VEL x NUM_VEL mass matrix.",
                  "Inner-controlled placement: TEMP_IN_SMEM selects where the scratch band lives (s_temp vs d_workspace). Decided at the top; caller sizes both arenas from CRBA_INNER_*_BYTES."]
    func_def_start = "void crba_inner("
    func_def_middle = "T *s_M, const T *s_q, const T *s_qd, "
    func_def_end = "T *s_temp, T *d_workspace, const T gravity) {"
    func_def_middle, func_params = self.gen_insert_helpers_func_def_params(func_def_middle, func_params, -2)
    func_def = func_def_start + func_def_middle + func_def_end
    self.gen_add_func_doc("Compute the Floating-Base Composite Rigid Body Algorithm", func_notes, func_params, None)
    self.gen_add_code_line("template <typename T, bool TEMP_IN_SMEM = true>")
    self.gen_add_code_line("__device__")
    self.gen_add_code_line(func_def, True)
    # Inner-controlled scratch-band placement: the whole band moves to
    # d_workspace when !TEMP_IN_SMEM. Reassigning s_temp at the top keeps every
    # s_temp[...] reference below unchanged.
    self.gen_add_code_line("if constexpr (!TEMP_IN_SMEM) { s_temp = d_workspace; } else { (void)d_workspace; }")
    temp_size = self.gen_crba_inner_temp_mem_size()
    self.gen_linalg_smem_setup(temp_size)

    self.gen_add_code_line("// Initialize IC = I and clear H")
    self.gen_add_parallel_loop("ind", str(36 * NJ + nv * nv))
    self.gen_add_code_line("if (ind < " + str(36 * NJ) + ") { s_temp[" + str(ICOffset) + " + ind] = s_XImats[" + str(36 * NJ) + " + ind]; }")
    self.gen_add_code_line("else { s_M[ind - " + str(36 * NJ) + "] = static_cast<T>(0); }")
    self.gen_add_end_control_flow()
    self.gen_add_sync()

    # Phase 1 — BFS-level body recursions.
    # For each BFS level (deepest to root, level 0 = floating root is skipped):
    #   alpha[i] = X[jid_i]^T * IC[jid_i]               (forward, per sibling slot i)
    #   IC[parent[jid_i]] += alpha[i] * X[jid_i]        (backward)
    # Levels are data-dependent (each writes its parents' IC), so the LEVEL loop is
    # sequential. WITHIN a level, sibling joints are independent (their alpha slots
    # are disjoint, and their parents are typically disjoint too). Fusing the
    # per-sibling 6x6 gemms into one block-cooperative pass per phase collapses
    # 2 * len(inds) syncs down to 2 per BFS level on branched robots
    # (go2-floating bfs widths 4,4,4 → 24 syncs → 6; g1-floating → 58 → 20).
    # Single-sibling levels (entire iiwa14 chain) emit the same byte-identical
    # GLASS gemm pair as before, so chain robots are unchanged.
    for bfs_level in range(n_bfs_levels - 1, 0, -1):
        inds = self.robot.get_ids_by_bfs_level(bfs_level)
        k = len(inds)
        joint_names = [self.robot.get_joint_by_id(j).get_name() for j in inds]
        self.gen_add_code_line(f"// CRBA Phase 1 BFS level {bfs_level} (jids {inds})")
        self.gen_add_code_line(f"//     joints: {', '.join(joint_names)}")

        if k == 1:
            # Serial-chain fast path: unchanged from prior emit. iiwa14 floating
            # falls entirely here (BFS levels 1..7 each have exactly one jid).
            jid = inds[0]
            parent = self.robot.get_parent_id(jid)
            self.gen_add_code_line(f"grid_linalg_gemm<T,6,6,6,false,true>(&s_XImats[{36*jid}], &s_temp[{ICOffset + 36*jid}], &s_temp[{alphaOffset}], static_cast<T>(1), static_cast<T>(0), s_linalg_smem);")
            self.gen_add_code_line(f"grid_linalg_gemm<T,6,6,6>(&s_temp[{alphaOffset}], &s_XImats[{36*jid}], &s_temp[{ICOffset + 36*parent}], static_cast<T>(1), static_cast<T>(1), s_linalg_smem);")
            continue

        # k siblings ≥ 2 — fused per-level forward + backward. Wrap in a block
        # scope so per-level compile-time tables (s_jid_lvl, s_par_lvl) don't
        # collide across BFS levels emitted into the same function body.
        jids_csv = ", ".join(str(j) for j in inds)
        parents_csv = ", ".join(str(self.robot.get_parent_id(j)) for j in inds)

        self.gen_add_code_line("{", True)
        # Forward: alpha[slot] = X[jid]^T * IC[jid] for slot in [0, k).
        # Thread `el` owns one (slot, row, col) triple → one output scalar.
        # 6x6 gemm with transposed A (= X^T): C[r,c] = sum_p X[p,r] * IC[p,c]
        # X is column-major in s_XImats: X[p,r] = s_XImats[36*jid + p + 6*r],
        # so X^T[r,p] = X[p,r] is read via Xj[p + 6*r] in the column-major slab.
        # IC is also column-major: IC[p,c] = s_temp[ICOffset + 36*jid + p + 6*c].
        self.gen_add_code_line(f"// fused forward: {k} siblings, each computes 36 outputs (alpha[slot] = X[jid_slot]^T * IC[jid_slot])")
        self.gen_add_code_line(f"const int s_jid_lvl[{k}] = {{{jids_csv}}};")
        self.gen_add_code_line(f"const int s_par_lvl[{k}] = {{{parents_csv}}};")
        self.gen_add_parallel_loop("el", str(36 * k))
        self.gen_add_code_line("int slot = el / 36;")
        self.gen_add_code_line("int rc = el % 36;")
        self.gen_add_code_line("int row = rc % 6;")
        self.gen_add_code_line("int col = rc / 6;")
        self.gen_add_code_line("int jid_l = s_jid_lvl[slot];")
        self.gen_add_code_line("const T *Xj = &s_XImats[36*jid_l];")
        self.gen_add_code_line(f"const T *ICj = &s_temp[{ICOffset} + 36*jid_l];")
        # X^T row `row` = X column `row` reading; dot with IC column `col`.
        self.gen_add_code_line("T acc = static_cast<T>(0);")
        self.gen_add_code_line("for (int p = 0; p < 6; p++) { acc += Xj[p + 6*row] * ICj[p + 6*col]; }")
        self.gen_add_code_line(f"s_temp[{alphaOffset} + 36*slot + row + 6*col] = acc;")
        self.gen_add_end_control_flow()
        self.gen_add_sync()

        # Backward: IC[parent] += alpha[slot] * X[jid].
        # When siblings share a parent (e.g. quadruped legs all feeding the
        # floating-base body), the accumulation cells COLLIDE and we must
        # atomicAdd into IC[parent]. When parents are disjoint (the common
        # case below the root), each (slot, r, c) writes a unique cell.
        if self.robot.has_repeated_parents(inds):
            self.gen_add_code_line(f"// fused backward (shared parents → atomicAdd): IC[parent] += alpha[slot] * X[jid_slot]")
            self.gen_add_parallel_loop("el", str(36 * k))
            self.gen_add_code_line("int slot = el / 36;")
            self.gen_add_code_line("int rc = el % 36;")
            self.gen_add_code_line("int row = rc % 6;")
            self.gen_add_code_line("int col = rc / 6;")
            self.gen_add_code_line("int jid_l = s_jid_lvl[slot];")
            self.gen_add_code_line("int par_l = s_par_lvl[slot];")
            self.gen_add_code_line(f"const T *alphaSlot = &s_temp[{alphaOffset} + 36*slot];")
            self.gen_add_code_line("const T *Xj = &s_XImats[36*jid_l];")
            # alpha is row-major effectively but stored col-major as a 6x6: alpha[r, p] at offset r + 6*p.
            # X column-major: X[p, c] at p + 6*c. Result[r, c] = sum_p alpha[r,p] * X[p,c].
            self.gen_add_code_line("T acc = static_cast<T>(0);")
            self.gen_add_code_line("for (int p = 0; p < 6; p++) { acc += alphaSlot[row + 6*p] * Xj[p + 6*col]; }")
            self.gen_add_code_line(f"atomicAdd(&s_temp[{ICOffset} + 36*par_l + row + 6*col], acc);")
            self.gen_add_end_control_flow()
            self.gen_add_sync()
            self.gen_add_end_control_flow()  # close the per-level scope
        else:
            # Disjoint parents → no collisions; one fused parallel loop over
            # (slot, r, c) writes IC[par[slot]] += alpha[slot] * X[jid_slot]
            # in registers, then commits to its unique IC cell.
            self.gen_add_code_line(f"// fused backward (disjoint parents): IC[parent] += alpha[slot] * X[jid_slot]")
            self.gen_add_parallel_loop("el", str(36 * k))
            self.gen_add_code_line("int slot = el / 36;")
            self.gen_add_code_line("int rc = el % 36;")
            self.gen_add_code_line("int row = rc % 6;")
            self.gen_add_code_line("int col = rc / 6;")
            self.gen_add_code_line("int jid_l = s_jid_lvl[slot];")
            self.gen_add_code_line("int par_l = s_par_lvl[slot];")
            self.gen_add_code_line(f"const T *alphaSlot = &s_temp[{alphaOffset} + 36*slot];")
            self.gen_add_code_line("const T *Xj = &s_XImats[36*jid_l];")
            self.gen_add_code_line("T acc = static_cast<T>(0);")
            self.gen_add_code_line("for (int p = 0; p < 6; p++) { acc += alphaSlot[row + 6*p] * Xj[p + 6*col]; }")
            self.gen_add_code_line(f"s_temp[{ICOffset} + 36*par_l + row + 6*col] += acc;")
            self.gen_add_end_control_flow()
            self.gen_add_sync()
            self.gen_add_end_control_flow()  # close the per-level scope

    # Phase 2 — per-jid thread-parallel walk for M's diagonal + scalar-joint
    # off-diagonal + floating-root coupling cells. ONE __syncthreads at the
    # end of the loop, vs the ~3-per-(jid, ancestor) of the prior impl
    # (~42 syncs/call on iiwa14-floating). See
    # docs/a3_core_dynamics_floating_loss_audit.md for the full refactor plan.
    self.gen_add_code_line("//")
    self.gen_add_code_line("// Phase 2: per-jid thread-parallel chain walk filling M's scalar-joint")
    self.gen_add_code_line("// diagonal + scalar/scalar off-diagonals + scalar/floating-root coupling.")
    self.gen_add_code_line("// Each thread owns ONE jid in [1, NJ) and walks its ancestor chain in")
    self.gen_add_code_line("// thread-local registers (s_fh, s_alpha); writes are race-free because")
    self.gen_add_code_line("// each thread only touches M cells indexed by its own dof = jid+5.")
    self.gen_add_code_line("//")
    max_ancestors = max(1, self.robot.get_max_num_ancestors())
    S_idx_arr = "{" + ", ".join(str(self.robot.get_S_index_by_id(j)) for j in range(NJ)) + "}"
    S_sgn_arr = "{" + ", ".join(str(self.robot.get_S_sign_by_id(j)) for j in range(NJ)) + "}"
    self.gen_add_parallel_loop("jid_off", str(NJ - 1))
    self.gen_add_code_line("int jid = jid_off + 1;          // jid in [1, NJ)")
    self.gen_add_code_line(f"int dof = jid + 5;")
    self.gen_add_code_line(f"const int s_Sidx_by_jid[{NJ}] = {S_idx_arr};")
    self.gen_add_code_line(f"const T s_Ssgn_by_jid[{NJ}] = {S_sgn_arr};")
    # Per-jid compile-time ancestor chain (matches fixed-base's pattern, _crba.py:165-183).
    parent_chain_init = "{" + ", ".join(["-1"] * max_ancestors) + "}"
    self.gen_add_code_line(f"int jid_parents[{max_ancestors}] = {parent_chain_init};")
    self.gen_add_code_line("int num_parents = 0;")
    self.gen_add_code_line("switch (jid) {", True)
    for jid in range(1, NJ):
        self.gen_add_code_line(f"case {jid}:", True)
        parent_chain = self.robot.get_ancestors_by_id(jid)
        for i, parent_ind in enumerate(parent_chain):
            self.gen_add_code_line(f"jid_parents[{i}] = {parent_ind};")
        self.gen_add_code_line(f"num_parents = {len(parent_chain)};")
        self.gen_add_code_line("break;")
        self.indent_level -= 1
    self.gen_add_end_control_flow()
    # Initialize fh = S_sgn[jid] * IC[jid][:, S_ind[jid]] in thread-local regs.
    self.gen_add_code_line("int sidx = s_Sidx_by_jid[jid];")
    self.gen_add_code_line("T   ssgn = s_Ssgn_by_jid[jid];")
    self.gen_add_code_line("T s_fh[6];")
    self.gen_add_code_line(f"for (int k = 0; k < 6; k++) s_fh[k] = ssgn * s_temp[{ICOffset} + 36*jid + 6*sidx + k];")
    # Diagonal M[dof, dof] = S_sgn * fh[sidx] = ssgn^2 * IC[jid][sidx, sidx] = IC[jid][sidx, sidx].
    self.gen_add_code_line(f"s_M[dof + {nv}*dof] = ssgn * s_fh[sidx];")
    # Chain walk: at step i, transform fh via X[X_id]^T and write M[jid, jid_parents[i]].
    self.gen_add_code_line("T s_alpha[6];")
    self.gen_add_code_line("for (int i = 0; i < num_parents; i++) {", True)
    self.gen_add_code_line("int X_id = (i == 0) ? jid : jid_parents[i-1];")
    self.gen_add_code_line("for (int k = 0; k < 6; k++) s_alpha[k] = s_fh[k];")
    self.gen_add_code_line("for (int k = 0; k < 6; k++) s_fh[k] = dot_prod<T,6,1,1>(&s_XImats[36*X_id + 6*k], &s_alpha[0]);")
    self.gen_add_code_line("int anc = jid_parents[i];")
    self.gen_add_code_line("if (anc > 0) {", True)
    self.gen_add_code_line("// scalar-joint ancestor: single M cell + its symmetric partner")
    self.gen_add_code_line("int a_dof = anc + 5;")
    self.gen_add_code_line("int a_sidx = s_Sidx_by_jid[anc];")
    self.gen_add_code_line("T   a_ssgn = s_Ssgn_by_jid[anc];")
    self.gen_add_code_line(f"s_M[dof + {nv}*a_dof] = a_ssgn * s_fh[a_sidx];")
    self.gen_add_code_line(f"s_M[a_dof + {nv}*dof] = s_M[dof + {nv}*a_dof];")
    self.gen_add_end_control_flow()
    self.gen_add_code_line("else {", True)
    self.gen_add_code_line("// floating-base root coupling: 6-wide M[dof, 0..5] row + symmetric col")
    self.gen_add_code_line("for (int col = 0; col < 6; col++) {", True)
    self.gen_add_code_line("int S_col = col < 3 ? col + 3 : col - 3;")
    self.gen_add_code_line(f"s_M[dof + {nv}*col] = s_fh[S_col];")
    self.gen_add_code_line(f"s_M[col + {nv}*dof] = s_M[dof + {nv}*col];")
    self.gen_add_end_control_flow()
    self.gen_add_end_control_flow()
    self.gen_add_end_control_flow()
    self.gen_add_end_control_flow()
    self.gen_add_sync()

    self.gen_add_code_line("// floating-base root block H[:6,:6] = S^T * IC[0] * S")
    self.gen_add_parallel_loop("ind", "36")
    self.gen_add_code_line("int row = ind % 6; int col = ind / 6;")
    self.gen_add_code_line("int S_row = row < 3 ? row + 3 : row - 3;")
    self.gen_add_code_line("int S_col = col < 3 ? col + 3 : col - 3;")
    self.gen_add_code_line("s_M[row + " + str(nv) + "*col] = s_temp[" + str(ICOffset) + " + S_row + 6*S_col];")
    self.gen_add_end_control_flow()
    self.gen_add_sync()

    self.gen_add_end_function()




def gen_crba_device_temp_mem_size(self):
    n = self.robot.get_num_joints()
    wrapper_size = self.gen_topology_helpers_size() + 72*n # for XImats
    return self.gen_crba_inner_temp_mem_size() + wrapper_size

def gen_crba_device(self):
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

    # then load/update XI and run the algo.
    # A.3 surgical lever: CRBA never reads s_XImats[0..35] (the floating root
    # spatial transform X[0]) — Phase 1 BFS starts at level 1, and Phase 2's
    # chain walk only dereferences X[X_id] for X_id in {jid, ancestors[:-1]};
    # the root id 0 only appears as `anc`, never as `X_id`. So we skip
    # recomputing X[0] from the floating-base quaternion on every call.
    self.gen_load_update_XImats_helpers_function_call(skip_floating_base_X=True)
    self.gen_crba_inner_function_call()
    self.gen_add_end_function()

def _emit_crba_kernel_body_for_flags(self, nq, nv, n, input_count, use_workspace_temp, single_call_timing):
    """Emit crba_kernel body for one tier's spill flag.
    use_workspace_temp=False: s_temp in smem (full arena); Level 0 / current.
    use_workspace_temp=True:  s_temp redirected to L2-pinned workspace; smem arena holds only extra_t_buffers."""
    shared_mem_size = 0 if use_workspace_temp else self.gen_crba_inner_temp_mem_size()
    self.gen_XImats_helpers_temp_shared_memory_code(shared_mem_size, extra_t_buffers = [("s_M", nv*nv), ("s_q_qd", input_count)], include_linalg_scratch=True)
    self.gen_add_code_line("T *s_q = s_q_qd; T *s_qd = &s_q_qd[" + str(nq) + "];")
    if not single_call_timing:
        # load to shared mem and loop over blocks to compute all requested comps
        self.gen_add_parallel_loop("k","NUM_TIMESTEPS",block_level = True)
        self.gen_kernel_load_inputs("q_qd",str(input_count),stride="stride_q_qd")
        if use_workspace_temp:
            self.gen_add_code_line("T *crba_d_workspace = reinterpret_cast<T *>(&d_workspace[k*GRID_WORKSPACE_BYTES_PER_TIMESTEP<T>()]);")
            # The whole inner arena spilled to global, so the smem s_temp slot is
            # null. Repoint s_temp at the workspace BEFORE the XImats helper call
            # so its sincos scratch (and the inner) have a valid backing store.
            self.gen_add_code_line("s_temp = crba_d_workspace;")
        else:
            self.gen_add_code_line("(void)d_workspace;")
        # compute
        # A.3 surgical lever: skip the floating-base X[0] recompute on every
        # CRBA call (see gen_crba_device for the trace argument).
        self.gen_add_code_line("// compute")
        self.gen_load_update_XImats_helpers_function_call(skip_floating_base_X=True)
        self.gen_crba_inner_function_call(
            updated_var_names = (dict(d_workspace_name = "crba_d_workspace") if use_workspace_temp else None),
            temp_in_smem_expr = ("false" if use_workspace_temp else "true"))
        self.gen_add_sync()
        # save to global  (stride = nv*nv per timestep — without this, batches
        # overlap since each block writes nv*nv elements starting at offset k*1)
        self.gen_kernel_save_result("M",str(nv*nv),stride=str(nv*nv))
        self.gen_add_end_control_flow()
    else:
        # repurpose NUM_TIMESTEPS for number of timing reps
        self.gen_kernel_load_inputs("q_qd",str(input_count))
        if use_workspace_temp:
            self.gen_add_code_line("T *crba_d_workspace = reinterpret_cast<T *>(d_workspace);")
            # See note above: repoint the null smem s_temp at the spilled workspace
            # before the XImats helper call so its scratch is valid (global) memory.
            self.gen_add_code_line("s_temp = crba_d_workspace;")
        else:
            self.gen_add_code_line("(void)d_workspace;")
        # then compute in loop for timing
        self.gen_add_code_line("// compute with NUM_TIMESTEPS as NUM_REPS for timing")
        self.gen_add_code_line("for (int rep = 0; rep < NUM_TIMESTEPS; rep++){", True)
        self.gen_anti_licm_input_reload("q_qd",str(input_count),feedback_from="M")
        # A.3 surgical lever: skip floating-base X[0] recompute (CRBA-safe).
        self.gen_load_update_XImats_helpers_function_call(skip_floating_base_X=True)
        self.gen_crba_inner_function_call(
            updated_var_names = (dict(d_workspace_name = "crba_d_workspace") if use_workspace_temp else None),
            temp_in_smem_expr = ("false" if use_workspace_temp else "true"))
        self.gen_anti_licm_output_write("M")
        self.gen_add_end_control_flow()
        # save to global
        self.gen_kernel_save_result("M",str(nv*nv))


def gen_crba_kernel(self, single_call_timing = False):
    nq = self.robot.get_num_pos()
    nv = self.robot.get_num_vel()
    n = self.robot.get_num_joints()
    input_count = nq + nv
    # define function def and params
    func_params = ["d_M is the pointer to the matrix of inertia", \
                    "d_workspace is the L2-pinned global spill buffer (used at LITE/MINIMAL on large robots)", \
                    "d_q_qd is the vector of joint positions and velocities", \
                    "stride_q_qd is the stride between each q, qd", \
                    "d_robotModel is the pointer to the initialized model specific helpers on the GPU (XImats, topology_helpers, etc.)", \
                    "gravity is the gravity constant", \
                    "num_timesteps is the length of the trajectory points we need to compute over (or overloaded as test_iters for timing)"]
    func_notes = []
    func_def_start = "void crba_kernel(T *d_M, unsigned char *d_workspace, const T *d_q_qd, const int stride_q_qd, "
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

    # Whole-arena spill lever: the CRBA inner scratch band is spilled as one
    # band to L2-pinned workspace at LITE/MINIMAL. Level 0 = arena in smem
    # (current); Level 1 = redirected to workspace.
    picks = getattr(self, "crba_spill_tier_3way", (0, 0, 0))
    if picks[0] == picks[1] == picks[2]:
        _emit_crba_kernel_body_for_flags(self, nq, nv, n, input_count, bool(picks[0]), single_call_timing)
    else:
        tier_names = ("TIER_SHARED", "TIER_LITE", "TIER_MINIMAL")
        for tier_idx, (tier_name, pick) in enumerate(zip(tier_names, picks)):
            head = "if constexpr (RESOURCE_TIER == " + tier_name + ") {" if tier_idx == 0 else \
                   "else if constexpr (RESOURCE_TIER == " + tier_name + ") {"
            self.gen_add_code_line(head, True)
            _emit_crba_kernel_body_for_flags(self, nq, nv, n, input_count, bool(pick), single_call_timing)
            self.gen_add_end_control_flow()
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
    func_call_start = "crba_kernel<T><<<block_dimms,thread_dimms,CRBA_DYNAMIC_SHARED_MEM_BYTES<T>()>>>(hd_data->d_M,hd_data->d_workspace,hd_data->d_q_qd,stride_q_qd,"
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

def gen_crba(self):
    # first generate the inner helpers
    self.gen_crba_inner()
    # then generate the device wrappers
    self.gen_crba_device()
    # then generate the kernels
    self.gen_crba_kernel(True)
    self.gen_crba_kernel(False)
    # then the host launch wrappers
    self.gen_crba_host(0)
    self.gen_crba_host(1)
    self.gen_crba_host(2)
