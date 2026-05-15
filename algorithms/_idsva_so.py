SHARED_MEMORY_JOINT_THRESHOLD = 10 # Max shared memory threshold => Write directly to RAM

def _idsva_so_as_index_list(index):
    if isinstance(index, list):
        return [int(v) for v in index]
    if isinstance(index, tuple):
        return [int(v) for v in index]
    if hasattr(index, "flatten"):
        return [int(v) for v in index.flatten()]
    return [int(index)]

def _idsva_so_int_array(values):
    return ", ".join(map(str, values)) if values else "0"

def _idsva_so_unit_axis(column):
    values = column.reshape(-1).tolist() if hasattr(column, "reshape") else list(column)
    for row, value in enumerate(values):
        value = float(value)
        if abs(value) == 1.0:
            return row, (1 if value > 0.0 else -1)
    raise ValueError("Floating IDSVA-SO expected unit joint subspace columns.")

def _idsva_so_floating_velocity_metadata(robot):
    num_bodies = robot.get_num_bodies()
    body_v_start = [0]
    body_v_index = []
    vel_to_body = [0] * robot.get_num_vel()
    vel_to_local_col = [0] * robot.get_num_vel()
    vel_s_index = [0] * robot.get_num_vel()
    vel_s_sign = [1] * robot.get_num_vel()

    for body_id in range(num_bodies):
        v_inds = _idsva_so_as_index_list(robot.get_joint_index_v(body_id))
        S = robot.get_S_by_id(body_id)
        if len(S.shape) == 1:
            S_cols = [S]
        else:
            S_cols = [S[:, col] for col in range(S.shape[1])]
        if len(v_inds) != len(S_cols):
            raise ValueError(
                "Floating IDSVA-SO velocity metadata expected one S column per velocity index."
            )
        for local_col, vel_index in enumerate(v_inds):
            body_v_index.append(vel_index)
            vel_to_body[vel_index] = body_id
            vel_to_local_col[vel_index] = local_col
            s_index, s_sign = _idsva_so_unit_axis(S_cols[local_col])
            vel_s_index[vel_index] = s_index
            vel_s_sign[vel_index] = s_sign
        body_v_start.append(len(body_v_index))

    subtree_v_start = [0]
    subtree_v_index = []
    successor_v_start = [0]
    successor_v_index = []
    ancestor_body_start = [0]
    ancestor_body_index = []
    for body_id in range(num_bodies):
        subtree = list(robot.get_subtree_by_id(body_id))
        successors = [subtree_body for subtree_body in subtree if subtree_body != body_id]
        ancestors = list(robot.get_ancestors_by_id(body_id))
        ancestors.insert(0, body_id)
        ancestors = ancestors[::-1]

        for subtree_body in subtree:
            subtree_v_index.extend(
                body_v_index[body_v_start[subtree_body]:body_v_start[subtree_body + 1]]
            )
        subtree_v_start.append(len(subtree_v_index))

        for successor_body in successors:
            successor_v_index.extend(
                body_v_index[body_v_start[successor_body]:body_v_start[successor_body + 1]]
            )
        successor_v_start.append(len(successor_v_index))

        ancestor_body_index.extend(ancestors)
        ancestor_body_start.append(len(ancestor_body_index))

    return {
        "body_v_start": body_v_start,
        "body_v_index": body_v_index,
        "vel_to_body": vel_to_body,
        "vel_to_local_col": vel_to_local_col,
        "vel_s_index": vel_s_index,
        "vel_s_sign": vel_s_sign,
        "subtree_v_start": subtree_v_start,
        "subtree_v_index": subtree_v_index,
        "successor_v_start": successor_v_start,
        "successor_v_index": successor_v_index,
        "ancestor_body_start": ancestor_body_start,
        "ancestor_body_index": ancestor_body_index,
    }

def gen_idsva_so_inner_temp_mem_size(self):
    """
    Returns the total size of the temporary memory required for the
    second order idsva inner function.

    Returns:
        int: The total size of the temporary memory required for the
        second order idsva inner function.
    """
    NV = self.robot.get_num_vel()
    num_bodies = self.robot.get_num_bodies()
    if self.robot.floating_base:
        body_mat_count = 8 * 36 * num_bodies
        body_vec_count = 7 * 6 * num_bodies
        vel_vec_count = 12 * 6 * NV
        vel_mat_count = 9 * 36 * NV
        base_count = body_mat_count + body_vec_count + vel_vec_count + vel_mat_count + 6 + 64
        # Gravity-shim shared portion only — d2X / d2a / d2f spill to `d_workspace`
        # (see `gen_floating_gravity_d2tau_dq_spill_count`).
        return int(base_count + gen_floating_gravity_d2tau_dq_shared_count(self))
    jids_a, ancestors = self.robot.get_jid_ancestor_ids(include_joint=True)
    return int(36 * NV * 10 + 30 * NV + 6 + len(jids_a)*36)

def _floating_gravity_lie_metadata(robot):
    """Per-velocity Lie-generator metadata for the floating-base gravity Hessian.

    For each velocity coordinate `vi` belonging to body `jid`:
      - `body[vi] = jid`
      - `s_index[vi]` is the row in S_local that the unit-axis column points along.
      - `s_sign[vi]` is the unit-axis sign in S_local.
      - `lie_sign_flip[vi] = 1` iff the right-Lie generator `B` such that
        `dXmat/dq = Xmat @ B` is `-crm(S_local)` rather than `+crm(S_local)` in this
        codebase's Xmat convention.

    Empirically (from `_spatial_xmat_derivative_func` for non-root joints), the GRiD
    codebase's `Xmat(q)` is defined such that `dXmat/dq = -Xmat @ crm(S_local)` for
    all non-root single-DoF joints. For the floating-base root, the Python helper at
    `_floating_gravity_d2tau_dq_lie_direct` uses `B = +crm` for rotation columns and
    `B = -crm` for translation columns (Featherstone xlt sign convention).

    Combined rule for `lie_sign_flip`:
      - Root rotation columns (local_col >= 3 for the floating-base joint): 0.
      - Root translation columns (local_col < 3): 1.
      - All non-root joints: 1.
    """
    num_bodies = robot.get_num_bodies()
    num_vel = robot.get_num_vel()
    body = [0] * num_vel
    s_index = [0] * num_vel
    s_sign = [1] * num_vel
    lie_sign_flip = [0] * num_vel

    for body_id in range(num_bodies):
        v_inds = _idsva_so_as_index_list(robot.get_joint_index_v(body_id))
        S = robot.get_S_by_id(body_id)
        if len(S.shape) == 1:
            S_cols = [S]
        else:
            S_cols = [S[:, col] for col in range(S.shape[1])]
        for local_col, vel_index in enumerate(v_inds):
            row_idx, sign = _idsva_so_unit_axis(S_cols[local_col])
            body[vel_index] = body_id
            s_index[vel_index] = row_idx
            s_sign[vel_index] = sign
            # Flip sign for non-root joints (sympy Xmat convention) and for root
            # translation columns (Featherstone xlt sign convention).
            if body_id != 0:
                lie_sign_flip[vel_index] = 1
            elif local_col < 3:
                lie_sign_flip[vel_index] = 1

    return {
        "body": body,
        "s_index": s_index,
        "s_sign": s_sign,
        # Keep the field name `is_root_translation` for backwards-compat with
        # the existing emission code; it now means "needs Lie-sign flip".
        "is_root_translation": lie_sign_flip,
    }

def gen_floating_gravity_d2tau_dq_spill_count(self):
    """Floats of the gravity-Hessian helper that live in the global `d_workspace`
    spill region (per timestep). The spill set is {grav_d2X, grav_d2a, grav_d2f} —
    the three O(NV²·NB) tensors that dominate the memory budget for larger robots.
    """
    NV = self.robot.get_num_vel()
    NB = self.robot.get_num_bodies()
    d2X_count = 36 * NV * NV
    d2a_count = 6 * NV * NV * NB
    d2f_count = 6 * NV * NV * NB
    return int(d2X_count + d2a_count + d2f_count)


def gen_floating_gravity_d2tau_dq_shared_count(self):
    """Floats of the gravity-Hessian helper that stay in shared memory.
    Includes dX (sparse-but-stored-dense), a/da, f/df, and the 6x6 scratch buffers.
    """
    NV = self.robot.get_num_vel()
    NB = self.robot.get_num_bodies()
    dX_count = 36 * NV
    a_count = 6 * NB
    da_count = 6 * NV * NB
    f_count = a_count
    df_count = da_count
    scratch_count = 4 * 36
    return int(dX_count + a_count + da_count + f_count + df_count + scratch_count)


def gen_floating_gravity_d2tau_dq_temp_mem_size(self):
    """Total floats the gravity-Hessian helper needs (shared + spill).

    The kernel emitter splits this into a shared-memory portion (the smaller arrays
    plus scratch) and a `d_workspace` spill portion (d2X / d2a / d2f) so larger
    robots (g1, etc.) still fit. See `gen_floating_gravity_d2tau_dq_shared_count`
    and `gen_floating_gravity_d2tau_dq_spill_count`.
    """
    return int(gen_floating_gravity_d2tau_dq_shared_count(self)
               + gen_floating_gravity_d2tau_dq_spill_count(self))

def gen_floating_gravity_d2tau_dq_lie_inline(self, use_thread_group=False):
    """Emit (inline) the floating-base gravity-Hessian addition into `d2tau_dq2`.

    Translates the Python helper `_floating_gravity_d2tau_dq_lie_direct` (see
    `RBDReference/RBDReference.py`) into CUDA emission for use inside
    `gen_idsva_so_floating_reference_inner`.

    Preconditions (set up by the caller):
      - `Xup[NB*36]` contains cumulative world-frame joint transforms.
      - `I[NB*36]` is body-frame inertia (the constant inertia tensor for each body).
      - `s_q` holds joint position parameters.
      - `d2tau_dq2[NV*NV*NV]` is the output buffer; we add to it here.

    Memory allocated from the caller's scratch (carved by the caller). Phase A
    assumes the carve fits; Phase D moves the large tensors to `d_workspace`.

    Algorithm (one-pass body-frame propagation):
      1. dX[vi] = X[jid(vi)] @ B_vi   (Lie generator, with translation sign flip
         for the floating-base root).
      2. d2X[vi][vj] = X[jid] @ B_vj @ B_vi when vi, vj share a body; else 0.
      3. Forward: a[jid] = X[jid] @ a[parent] with `a[parent_of_root] = -gravity`;
         carry first/second derivatives via the standard chain rule.
      4. f = I @ a per body (all bodies, batched conceptually).
      5. Backward: project onto each joint's S, propagate f-derivatives to parent.
      6. Output: d2tau_dq2[v_index_of_jid_dof, :, :] += S^T @ d2f[jid, :, :].

    See `RBDReference.py` line ~2371 for the body-major Python reference this mirrors.
    """
    NV = self.robot.get_num_vel()
    NB = self.robot.get_num_bodies()
    lie_meta = _floating_gravity_lie_metadata(self.robot)
    parent_ids = [self.robot.get_parent_id(b) for b in range(NB)]

    # Wrap the entire emission in a single-thread block so the helper is self-contained
    # and safe to call regardless of whether the caller is already in a thread-zero scope.
    self.gen_add_code_line("if (threadIdx.x == 0 && threadIdx.y == 0 && threadIdx.z == 0) {", True)
    # Compute the shared-memory base offset: existing main-sweep temp size MINUS the
    # gravity-shim's shared portion (we want grav_scratch to point at where the helper's
    # shared arrays live, which is right after the main-sweep allocations).
    main_sweep_count = self.gen_idsva_so_inner_temp_mem_size() - gen_floating_gravity_d2tau_dq_shared_count(self)
    self.gen_add_code_lines([
        "// ===== Gravity-Hessian (Lie-tangent) addition into d2tau_dq2 =====",
        "// Shared portion (dX / a / da / f / df / 4x36 scratch) lives in s_temp; the",
        "// three O(NV*NV*NB) tensors (d2X / d2a / d2f) spill to s_temp_spill which the",
        "// kernel emitter points into d_workspace per-timestep.",
        f"static const int grav_lie_body[] = {{ {_idsva_so_int_array(lie_meta['body'])} }};",
        f"static const int grav_lie_s_index[] = {{ {_idsva_so_int_array(lie_meta['s_index'])} }};",
        f"static const int grav_lie_s_sign[] = {{ {_idsva_so_int_array(lie_meta['s_sign'])} }};",
        f"static const int grav_lie_is_root_translation[] = {{ {_idsva_so_int_array(lie_meta['is_root_translation'])} }};",
        f"static const int grav_lie_parent[] = {{ {_idsva_so_int_array(parent_ids)} }};",
        "",
        "// Shared-memory carve (small arrays).",
        f"T *grav_scratch = s_temp + {main_sweep_count};",
        "T *grav_dX      = grav_scratch;",
        "T *grav_a       = grav_dX      + 36*NUM_VEL;",
        "T *grav_da      = grav_a       + 6*NUM_BODIES;",
        "T *grav_f       = grav_da      + 6*NUM_VEL*NUM_BODIES;",
        "T *grav_df      = grav_f       + 6*NUM_BODIES;",
        "T *grav_invX    = grav_df      + 6*NUM_VEL*NUM_BODIES;",
        "T *grav_tmpA    = grav_invX    + 36;",
        "T *grav_tmpB    = grav_tmpA    + 36;",
        "T *grav_tmpC    = grav_tmpB    + 36;",
        "// Global-memory spill carve (large per-timestep tensors).",
        "T *grav_d2X     = s_temp_spill;",
        "T *grav_d2a     = grav_d2X     + 36*NUM_VEL*NUM_VEL;",
        "T *grav_d2f     = grav_d2a     + 6*NUM_VEL*NUM_VEL*NUM_BODIES;",
        "",
        "// gravity_vec mirrors Python's `gravity_vec[5] = -GRAVITY`. In CUDA the `gravity`",
        "// parameter is the *positive magnitude* of gravitational acceleration (= 9.81),",
        "// matching the convention used by the rest of the GRiD CUDA functions; the Python",
        "// helper takes signed GRAVITY = -9.81. So Python's `gravity_vec[5] = -GRAVITY = +9.81`",
        "// equals CUDA's `gravity_vec[5] = +gravity`.",
        "T grav_gravity_vec[6] = { static_cast<T>(0), static_cast<T>(0), static_cast<T>(0),",
        "                          static_cast<T>(0), static_cast<T>(0), gravity };",
        "",
        "// Zero scratch (single-thread; could be parallelised across threadIdx for speed).",
        "for (int idx = 0; idx < 36*NUM_VEL; ++idx) grav_dX[idx] = static_cast<T>(0);",
        "for (int idx = 0; idx < 36*NUM_VEL*NUM_VEL; ++idx) grav_d2X[idx] = static_cast<T>(0);",
        "for (int idx = 0; idx < 6*NUM_BODIES; ++idx) { grav_a[idx] = static_cast<T>(0); grav_f[idx] = static_cast<T>(0); }",
        "for (int idx = 0; idx < 6*NUM_VEL*NUM_BODIES; ++idx) { grav_da[idx] = static_cast<T>(0); grav_df[idx] = static_cast<T>(0); }",
        "for (int idx = 0; idx < 6*NUM_VEL*NUM_VEL*NUM_BODIES; ++idx) { grav_d2a[idx] = static_cast<T>(0); grav_d2f[idx] = static_cast<T>(0); }",
        "",
        "// ---- (1) Build Lie generators B[vi] = +/- crm(unit_axis_si).",
        "//        Root (jid==0):  dX[vi] = X_local[0] @ B[vi]   (Featherstone xlt convention,",
        "//                        sign flip on translation columns baked into B).",
        "//        Non-root joints: dX[vi] = B[vi] @ X_local[jid] (codebase's sympy convention:",
        "//                         dXmat/dq = -crm(S_local) @ Xmat, with B = -crm(S_local)).",
        "for (int vi = 0; vi < NUM_VEL; ++vi) {", True,
        "int jid = grav_lie_body[vi];",
        "int s_row = grav_lie_s_index[vi];",
        "T sign = static_cast<T>(grav_lie_s_sign[vi]);",
        "if (grav_lie_is_root_translation[vi]) sign = -sign;",
        "T B[36];",
        "T e_vec[6] = { static_cast<T>(0), static_cast<T>(0), static_cast<T>(0),",
        "               static_cast<T>(0), static_cast<T>(0), static_cast<T>(0) };",
        "e_vec[s_row] = sign;",
        "for (int idx = 0; idx < 36; ++idx) B[idx] = crm<T>(idx, e_vec);",
        "for (int idx = 0; idx < 36; ++idx) {", True,
        "int row = idx % 6; int col = idx / 6;",
        "T acc = static_cast<T>(0);",
        "if (jid == 0) {",
        "    for (int kk = 0; kk < 6; ++kk) acc += s_XImats[jid*36 + row + 6*kk] * B[kk + 6*col];",
        "} else {",
        "    for (int kk = 0; kk < 6; ++kk) acc += B[row + 6*kk] * s_XImats[jid*36 + kk + 6*col];",
        "}",
        "grav_dX[vi*36 + idx] = acc;",
        "",  # close inner idx loop
        ])
    self.gen_add_end_control_flow()
    self.gen_add_end_control_flow()

    self.gen_add_code_lines([
        "",
        "// ---- (2) Build d2X[vi][vj] when vi, vj share a body.",
        "//        Root (jid_i==0):  d2X[vi,vj] = X_local[0] @ Bj @ Bi   (right multiplication).",
        "//        Non-root:         d2X[vi,vj] = Bj @ Bi @ X_local[jid] (left multiplication;",
        "//                          for single-DoF non-root vi==vj this is B^2 @ X_local).",
        "for (int vi = 0; vi < NUM_VEL; ++vi) {", True,
        "int jid_i = grav_lie_body[vi];",
        "for (int vj = 0; vj < NUM_VEL; ++vj) {", True,
        "if (grav_lie_body[vj] != jid_i) continue;  // d2X is zero for cross-body pairs.",
        "T sign_i = static_cast<T>(grav_lie_s_sign[vi]);",
        "if (grav_lie_is_root_translation[vi]) sign_i = -sign_i;",
        "T sign_j = static_cast<T>(grav_lie_s_sign[vj]);",
        "if (grav_lie_is_root_translation[vj]) sign_j = -sign_j;",
        "T Bi[36], Bj[36];",
        "T ei[6] = { static_cast<T>(0), static_cast<T>(0), static_cast<T>(0),",
        "            static_cast<T>(0), static_cast<T>(0), static_cast<T>(0) };",
        "T ej[6] = { static_cast<T>(0), static_cast<T>(0), static_cast<T>(0),",
        "            static_cast<T>(0), static_cast<T>(0), static_cast<T>(0) };",
        "ei[grav_lie_s_index[vi]] = sign_i;",
        "ej[grav_lie_s_index[vj]] = sign_j;",
        "for (int idx = 0; idx < 36; ++idx) { Bi[idx] = crm<T>(idx, ei); Bj[idx] = crm<T>(idx, ej); }",
        "// tmpA = Bj @ Bi  (6x6 @ 6x6).",
        "for (int idx = 0; idx < 36; ++idx) {", True,
        "int row = idx % 6; int col = idx / 6;",
        "T acc = static_cast<T>(0);",
        "for (int kk = 0; kk < 6; ++kk) acc += Bj[row + 6*kk] * Bi[kk + 6*col];",
        "grav_tmpA[idx] = acc;",
        ])
    self.gen_add_end_control_flow()
    self.gen_add_code_line("// Multiply tmpA by X_local on the correct side (right for root, left for non-root).")
    self.gen_add_code_line("for (int idx = 0; idx < 36; ++idx) {", True)
    self.gen_add_code_line("int row = idx % 6; int col = idx / 6;")
    self.gen_add_code_line("T acc = static_cast<T>(0);")
    self.gen_add_code_line("if (jid_i == 0) {")
    self.gen_add_code_line("    for (int kk = 0; kk < 6; ++kk) acc += s_XImats[jid_i*36 + row + 6*kk] * grav_tmpA[kk + 6*col];")
    self.gen_add_code_line("} else {")
    self.gen_add_code_line("    for (int kk = 0; kk < 6; ++kk) acc += grav_tmpA[row + 6*kk] * s_XImats[jid_i*36 + kk + 6*col];")
    self.gen_add_code_line("}")
    self.gen_add_code_line("grav_d2X[(vi*NUM_VEL + vj)*36 + idx] = acc;")
    self.gen_add_end_control_flow()
    self.gen_add_end_control_flow()
    self.gen_add_end_control_flow()

    self.gen_add_code_lines([
        "",
        "// ---- (3) Forward sweep: propagate a, da, d2a through the tree (body frame).",
        "//        Root (parent < 0): a[0] = inv(X_local[0]) @ gravity_vec. Since X_local[0]",
        "//        equals the cumulative Xup[0] (no parent in the chain), Xdown[0] is its inverse.",
        "//        Non-root: a[jid] = X_local[jid] @ a[parent]; carry first/second derivatives.",
        "for (int jid = 0; jid < NUM_BODIES; ++jid) {", True,
        "int parent = grav_lie_parent[jid];",
        "if (parent < 0) {", True,
        "    // ---- Root: a[0] = inv(X[0]) @ gravity_vec.",
        "    for (int row = 0; row < 6; ++row) {", True,
        "T acc = static_cast<T>(0);",
        "for (int kk = 0; kk < 6; ++kk) acc += Xdown[jid*36 + row + 6*kk] * grav_gravity_vec[kk];",
        "grav_a[jid*6 + row] = acc;",
        ])
    self.gen_add_end_control_flow()  # close root row loop

    # da[0, vi] = -inv_X @ dX[vi] @ a[0], only for vi whose body is root (else stays 0).
    self.gen_add_code_lines([
        "    // da[0, vi] = -inv_X @ dX[vi] @ a[0] for root coords; zero otherwise.",
        "    for (int vi = 0; vi < NUM_VEL; ++vi) {", True,
        "if (grav_lie_body[vi] != jid) continue;",
        "// tmp = dX[vi] @ a[0]   (6-vector).",
        "T tmp_v[6];",
        "for (int row = 0; row < 6; ++row) {", True,
        "T acc = static_cast<T>(0);",
        "for (int kk = 0; kk < 6; ++kk) acc += grav_dX[vi*36 + row + 6*kk] * grav_a[jid*6 + kk];",
        "tmp_v[row] = acc;",
        ])
    self.gen_add_end_control_flow()  # close tmp_v row loop
    self.gen_add_code_lines([
        "// grav_da[jid, vi] = -inv_X @ tmp_v = -Xdown[jid] @ tmp_v.",
        "for (int row = 0; row < 6; ++row) {", True,
        "T acc = static_cast<T>(0);",
        "for (int kk = 0; kk < 6; ++kk) acc += Xdown[jid*36 + row + 6*kk] * tmp_v[kk];",
        "grav_da[(jid*NUM_VEL + vi)*6 + row] = -acc;",
        ])
    self.gen_add_end_control_flow()  # close da row loop
    self.gen_add_end_control_flow()  # close vi loop

    # d2a[0, vi, vj] for root coords only.
    self.gen_add_code_lines([
        "    // d2a[0, vi, vj] = -inv_X @ (dX[vi]@da[0,vj] + dX[vj]@da[0,vi] + d2X[vi,vj]@a[0]).",
        "    // Nonzero only when both vi and vj are coords of the root body.",
        "    for (int vi = 0; vi < NUM_VEL; ++vi) {", True,
        "if (grav_lie_body[vi] != jid) continue;",
        "for (int vj = 0; vj < NUM_VEL; ++vj) {", True,
        "if (grav_lie_body[vj] != jid) continue;",
        "T inner_v[6] = { static_cast<T>(0), static_cast<T>(0), static_cast<T>(0),",
        "                 static_cast<T>(0), static_cast<T>(0), static_cast<T>(0) };",
        "// dX[vi] @ da[0, vj].",
        "for (int row = 0; row < 6; ++row) {", True,
        "for (int kk = 0; kk < 6; ++kk) inner_v[row] += grav_dX[vi*36 + row + 6*kk] * grav_da[(jid*NUM_VEL + vj)*6 + kk];",
        ])
    self.gen_add_end_control_flow()  # close row loop term1
    self.gen_add_code_lines([
        "// dX[vj] @ da[0, vi].",
        "for (int row = 0; row < 6; ++row) {", True,
        "for (int kk = 0; kk < 6; ++kk) inner_v[row] += grav_dX[vj*36 + row + 6*kk] * grav_da[(jid*NUM_VEL + vi)*6 + kk];",
        ])
    self.gen_add_end_control_flow()  # close row loop term2
    self.gen_add_code_lines([
        "// d2X[vi, vj] @ a[0].",
        "for (int row = 0; row < 6; ++row) {", True,
        "for (int kk = 0; kk < 6; ++kk) inner_v[row] += grav_d2X[(vi*NUM_VEL + vj)*36 + row + 6*kk] * grav_a[jid*6 + kk];",
        ])
    self.gen_add_end_control_flow()  # close row loop term3
    self.gen_add_code_lines([
        "// grav_d2a[0, vi, vj] = -Xdown[0] @ inner_v.",
        "for (int row = 0; row < 6; ++row) {", True,
        "T acc = static_cast<T>(0);",
        "for (int kk = 0; kk < 6; ++kk) acc += Xdown[jid*36 + row + 6*kk] * inner_v[kk];",
        "grav_d2a[((jid*NUM_VEL + vi)*NUM_VEL + vj)*6 + row] = -acc;",
        ])
    self.gen_add_end_control_flow()  # close d2a row loop
    self.gen_add_end_control_flow()  # close vj loop
    self.gen_add_end_control_flow()  # close vi loop

    self.gen_add_end_control_flow()  # close parent < 0 branch
    # Non-root branch.
    self.gen_add_code_lines([
        "else {", True,
        "    // ---- Non-root: a[jid] = X_local[jid] @ a[parent].",
        "    for (int row = 0; row < 6; ++row) {", True,
        "T acc = static_cast<T>(0);",
        "for (int kk = 0; kk < 6; ++kk) acc += s_XImats[jid*36 + row + 6*kk] * grav_a[parent*6 + kk];",
        "grav_a[jid*6 + row] = acc;",
        ])
    self.gen_add_end_control_flow()  # close a row loop
    self.gen_add_code_lines([
        "    // da[jid, vi] = dX[jid, vi] @ a[parent]   (only when vi is a coord of jid)",
        "    //            + X_local[jid] @ da[parent, vi]   (always for any vi that affects parent).",
        "    for (int vi = 0; vi < NUM_VEL; ++vi) {", True,
        "T row_vals[6] = { static_cast<T>(0), static_cast<T>(0), static_cast<T>(0),",
        "                  static_cast<T>(0), static_cast<T>(0), static_cast<T>(0) };",
        "// First term: dX[jid, vi] @ a[parent], only when vi is a coord of jid.",
        "if (grav_lie_body[vi] == jid) {", True,
        "for (int row = 0; row < 6; ++row) {", True,
        "for (int kk = 0; kk < 6; ++kk) row_vals[row] += grav_dX[vi*36 + row + 6*kk] * grav_a[parent*6 + kk];",
        ])
    self.gen_add_end_control_flow()  # row loop
    self.gen_add_end_control_flow()  # if dX nonzero
    self.gen_add_code_lines([
        "// Second term: X_local[jid] @ da[parent, vi].",
        "for (int row = 0; row < 6; ++row) {", True,
        "T acc = static_cast<T>(0);",
        "for (int kk = 0; kk < 6; ++kk) acc += s_XImats[jid*36 + row + 6*kk] * grav_da[(parent*NUM_VEL + vi)*6 + kk];",
        "row_vals[row] += acc;",
        ])
    self.gen_add_end_control_flow()  # row loop second term
    self.gen_add_code_lines([
        "for (int row = 0; row < 6; ++row) grav_da[(jid*NUM_VEL + vi)*6 + row] = row_vals[row];",
        ])
    self.gen_add_end_control_flow()  # vi loop
    # d2a non-root: 4 terms.
    self.gen_add_code_lines([
        "    // d2a[jid, vi, vj] = d2X[jid, vi, vj] @ a[parent]   (only when both vi, vj coords of jid)",
        "    //               + dX[jid, vi] @ da[parent, vj]   (only when vi is a coord of jid)",
        "    //               + dX[jid, vj] @ da[parent, vi]   (only when vj is a coord of jid)",
        "    //               + X_local[jid] @ d2a[parent, vi, vj]   (always).",
        "    for (int vi = 0; vi < NUM_VEL; ++vi) {", True,
        "for (int vj = 0; vj < NUM_VEL; ++vj) {", True,
        "T row_vals[6] = { static_cast<T>(0), static_cast<T>(0), static_cast<T>(0),",
        "                  static_cast<T>(0), static_cast<T>(0), static_cast<T>(0) };",
        "// Term 1: d2X[jid, vi, vj] @ a[parent].",
        "if (grav_lie_body[vi] == jid && grav_lie_body[vj] == jid) {", True,
        "for (int row = 0; row < 6; ++row) {", True,
        "for (int kk = 0; kk < 6; ++kk) row_vals[row] += grav_d2X[(vi*NUM_VEL + vj)*36 + row + 6*kk] * grav_a[parent*6 + kk];",
        ])
    self.gen_add_end_control_flow()  # row loop
    self.gen_add_end_control_flow()  # if d2X nonzero
    self.gen_add_code_lines([
        "// Term 2: dX[jid, vi] @ da[parent, vj].",
        "if (grav_lie_body[vi] == jid) {", True,
        "for (int row = 0; row < 6; ++row) {", True,
        "for (int kk = 0; kk < 6; ++kk) row_vals[row] += grav_dX[vi*36 + row + 6*kk] * grav_da[(parent*NUM_VEL + vj)*6 + kk];",
        ])
    self.gen_add_end_control_flow()  # row loop
    self.gen_add_end_control_flow()  # if dX[vi] nonzero
    self.gen_add_code_lines([
        "// Term 3: dX[jid, vj] @ da[parent, vi].",
        "if (grav_lie_body[vj] == jid) {", True,
        "for (int row = 0; row < 6; ++row) {", True,
        "for (int kk = 0; kk < 6; ++kk) row_vals[row] += grav_dX[vj*36 + row + 6*kk] * grav_da[(parent*NUM_VEL + vi)*6 + kk];",
        ])
    self.gen_add_end_control_flow()  # row loop
    self.gen_add_end_control_flow()  # if dX[vj] nonzero
    self.gen_add_code_lines([
        "// Term 4: X_local[jid] @ d2a[parent, vi, vj].",
        "for (int row = 0; row < 6; ++row) {", True,
        "T acc = static_cast<T>(0);",
        "for (int kk = 0; kk < 6; ++kk) acc += s_XImats[jid*36 + row + 6*kk] * grav_d2a[((parent*NUM_VEL + vi)*NUM_VEL + vj)*6 + kk];",
        "row_vals[row] += acc;",
        ])
    self.gen_add_end_control_flow()  # row loop term 4
    self.gen_add_code_lines([
        "for (int row = 0; row < 6; ++row) grav_d2a[((jid*NUM_VEL + vi)*NUM_VEL + vj)*6 + row] = row_vals[row];",
        ])
    self.gen_add_end_control_flow()  # vj loop
    self.gen_add_end_control_flow()  # vi loop
    self.gen_add_end_control_flow()  # else branch
    self.gen_add_end_control_flow()  # jid loop

    # ---- (4) f = I @ a, df = I @ da, d2f = I @ d2a per body. Body-frame inertia
    # `I` is the per-body Imat from the layout (set up by the caller, see line ~324).
    self.gen_add_code_lines([
        "",
        "// ---- (4) f = I @ a (and derivatives) per body, body-frame inertia.",
        "for (int jid = 0; jid < NUM_BODIES; ++jid) {", True,
        "for (int row = 0; row < 6; ++row) {", True,
        "T acc = static_cast<T>(0);",
        "for (int kk = 0; kk < 6; ++kk) acc += I[jid*36 + row + 6*kk] * grav_a[jid*6 + kk];",
        "grav_f[jid*6 + row] = acc;",
        ])
    self.gen_add_end_control_flow()  # f row loop
    self.gen_add_code_lines([
        "for (int vi = 0; vi < NUM_VEL; ++vi) {", True,
        "for (int row = 0; row < 6; ++row) {", True,
        "T acc = static_cast<T>(0);",
        "for (int kk = 0; kk < 6; ++kk) acc += I[jid*36 + row + 6*kk] * grav_da[(jid*NUM_VEL + vi)*6 + kk];",
        "grav_df[(jid*NUM_VEL + vi)*6 + row] = acc;",
        ])
    self.gen_add_end_control_flow()  # df row loop
    self.gen_add_end_control_flow()  # vi loop for df
    self.gen_add_code_lines([
        "for (int vi = 0; vi < NUM_VEL; ++vi) {", True,
        "for (int vj = 0; vj < NUM_VEL; ++vj) {", True,
        "for (int row = 0; row < 6; ++row) {", True,
        "T acc = static_cast<T>(0);",
        "for (int kk = 0; kk < 6; ++kk) acc += I[jid*36 + row + 6*kk] * grav_d2a[((jid*NUM_VEL + vi)*NUM_VEL + vj)*6 + kk];",
        "grav_d2f[((jid*NUM_VEL + vi)*NUM_VEL + vj)*6 + row] = acc;",
        ])
    self.gen_add_end_control_flow()  # d2f row loop
    self.gen_add_end_control_flow()  # vj loop
    self.gen_add_end_control_flow()  # vi loop for d2f
    self.gen_add_end_control_flow()  # jid loop

    self.gen_add_code_lines([
        "",
        "// ---- (5) Backward sweep: project d2f onto each joint's S (body-frame motion subspace),",
        "//        ADD into d2tau_dq2 (caller's main sweep is expected to have run with gravity=0).",
        "//        Then bubble f, df, d2f to the parent via X_local transposes.",
        "for (int jid = NUM_BODIES - 1; jid >= 0; --jid) {", True,
        "// Projection: for each velocity coord finds_c of jid, d2tau_dq2[finds_c, k, l] += S_body[:, c] . d2f[jid, k, l, :].",
        "// Body-frame S has a single nonzero entry per column: row=grav_lie_s_index[finds_c], value=grav_lie_s_sign[finds_c].",
        "for (int pos = body_v_start[jid]; pos < body_v_start[jid + 1]; ++pos) {", True,
        "int finds_c = body_v_index[pos];",
        "int s_row = grav_lie_s_index[finds_c];",
        "T s_sign = static_cast<T>(grav_lie_s_sign[finds_c]);",
        "for (int k = 0; k < NUM_VEL; ++k) {", True,
        "for (int l = 0; l < NUM_VEL; ++l) {", True,
        "T proj = s_sign * grav_d2f[((jid*NUM_VEL + k)*NUM_VEL + l)*6 + s_row];",
        "d2tau_dq2[(finds_c*NUM_VEL + k)*NUM_VEL + l] += proj;",
        ])
    self.gen_add_end_control_flow()  # l loop
    self.gen_add_end_control_flow()  # k loop
    self.gen_add_end_control_flow()  # pos loop (velocity coords of jid)

    # Parent update: skip when root.
    self.gen_add_code_lines([
        "int parent = grav_lie_parent[jid];",
        "if (parent < 0) continue;",
        "// Build X_local transpose once per body (used in all three parent updates).",
        "T Xt[36];",
        "for (int idx = 0; idx < 36; ++idx) {", True,
        "int row = idx % 6; int col = idx / 6;",
        "Xt[idx] = s_XImats[jid*36 + col + 6*row];",
        ])
    self.gen_add_end_control_flow()  # Xt loop
    self.gen_add_code_lines([
        "// f[parent] += Xt @ f[jid].",
        "for (int row = 0; row < 6; ++row) {", True,
        "T acc = static_cast<T>(0);",
        "for (int kk = 0; kk < 6; ++kk) acc += Xt[row + 6*kk] * grav_f[jid*6 + kk];",
        "grav_f[parent*6 + row] += acc;",
        ])
    self.gen_add_end_control_flow()  # f update row loop

    # df parent update.
    self.gen_add_code_lines([
        "// df[parent, vi] += dXt[vi] @ f[jid]   (only when vi is a coord of jid)",
        "//                + Xt @ df[jid, vi]   (always).",
        "for (int vi = 0; vi < NUM_VEL; ++vi) {", True,
        "T row_vals[6] = { static_cast<T>(0), static_cast<T>(0), static_cast<T>(0),",
        "                  static_cast<T>(0), static_cast<T>(0), static_cast<T>(0) };",
        "if (grav_lie_body[vi] == jid) {", True,
        "// dXt[vi] uses dX[vi].T, i.e., swap row<->col indices.",
        "for (int row = 0; row < 6; ++row) {", True,
        "for (int kk = 0; kk < 6; ++kk) row_vals[row] += grav_dX[vi*36 + kk + 6*row] * grav_f[jid*6 + kk];",
        ])
    self.gen_add_end_control_flow()  # dXt row loop
    self.gen_add_end_control_flow()  # if dX nonzero
    self.gen_add_code_lines([
        "for (int row = 0; row < 6; ++row) {", True,
        "T acc = static_cast<T>(0);",
        "for (int kk = 0; kk < 6; ++kk) acc += Xt[row + 6*kk] * grav_df[(jid*NUM_VEL + vi)*6 + kk];",
        "row_vals[row] += acc;",
        ])
    self.gen_add_end_control_flow()  # Xt @ df row loop
    self.gen_add_code_lines([
        "for (int row = 0; row < 6; ++row) grav_df[(parent*NUM_VEL + vi)*6 + row] += row_vals[row];",
        ])
    self.gen_add_end_control_flow()  # vi loop for df

    # d2f parent update.
    self.gen_add_code_lines([
        "// d2f[parent, vi, vj] += d2Xt[vi,vj]@f[jid] + dXt[vi]@df[jid,vj] + dXt[vj]@df[jid,vi] + Xt@d2f[jid,vi,vj].",
        "for (int vi = 0; vi < NUM_VEL; ++vi) {", True,
        "for (int vj = 0; vj < NUM_VEL; ++vj) {", True,
        "T row_vals[6] = { static_cast<T>(0), static_cast<T>(0), static_cast<T>(0),",
        "                  static_cast<T>(0), static_cast<T>(0), static_cast<T>(0) };",
        "// Term 1: d2Xt[vi, vj] @ f[jid] — nonzero only when both vi, vj are coords of jid.",
        "if (grav_lie_body[vi] == jid && grav_lie_body[vj] == jid) {", True,
        "for (int row = 0; row < 6; ++row) {", True,
        "for (int kk = 0; kk < 6; ++kk) row_vals[row] += grav_d2X[(vi*NUM_VEL + vj)*36 + kk + 6*row] * grav_f[jid*6 + kk];",
        ])
    self.gen_add_end_control_flow()  # term 1 row loop
    self.gen_add_end_control_flow()  # if d2X nonzero
    self.gen_add_code_lines([
        "// Term 2: dXt[vi] @ df[jid, vj].",
        "if (grav_lie_body[vi] == jid) {", True,
        "for (int row = 0; row < 6; ++row) {", True,
        "for (int kk = 0; kk < 6; ++kk) row_vals[row] += grav_dX[vi*36 + kk + 6*row] * grav_df[(jid*NUM_VEL + vj)*6 + kk];",
        ])
    self.gen_add_end_control_flow()  # term 2 row loop
    self.gen_add_end_control_flow()  # if dX[vi] nonzero
    self.gen_add_code_lines([
        "// Term 3: dXt[vj] @ df[jid, vi].",
        "if (grav_lie_body[vj] == jid) {", True,
        "for (int row = 0; row < 6; ++row) {", True,
        "for (int kk = 0; kk < 6; ++kk) row_vals[row] += grav_dX[vj*36 + kk + 6*row] * grav_df[(jid*NUM_VEL + vi)*6 + kk];",
        ])
    self.gen_add_end_control_flow()  # term 3 row loop
    self.gen_add_end_control_flow()  # if dX[vj] nonzero
    self.gen_add_code_lines([
        "// Term 4: Xt @ d2f[jid, vi, vj].",
        "for (int row = 0; row < 6; ++row) {", True,
        "T acc = static_cast<T>(0);",
        "for (int kk = 0; kk < 6; ++kk) acc += Xt[row + 6*kk] * grav_d2f[((jid*NUM_VEL + vi)*NUM_VEL + vj)*6 + kk];",
        "row_vals[row] += acc;",
        ])
    self.gen_add_end_control_flow()  # term 4 row loop
    self.gen_add_code_lines([
        "for (int row = 0; row < 6; ++row) grav_d2f[((parent*NUM_VEL + vi)*NUM_VEL + vj)*6 + row] += row_vals[row];",
        ])
    self.gen_add_end_control_flow()  # vj loop
    self.gen_add_end_control_flow()  # vi loop for d2f
    self.gen_add_end_control_flow()  # jid backward loop
    self.gen_add_end_control_flow()  # close thread-zero wrap
    self.gen_add_sync(use_thread_group)

def gen_idsva_so_inner_function_call(self, use_thread_group = False, use_qdd_input = False, updated_var_names = None):
    var_names = dict( \
        s_idsva_so_name = "s_idsva_so", \
        s_q_name = "s_q", \
        s_qd_name = "s_qd", \
        s_qdd_name = "s_qdd", \
        s_temp_name = "s_temp", \
        s_temp_spill_name = "s_temp_spill", \
        gravity_name = "gravity"
    )
    if updated_var_names is not None:
        for key,value in updated_var_names.items():
            var_names[key] = value
    id_so_code_start = "idsva_so_inner<T>(" + var_names["s_idsva_so_name"] + ", " + var_names["s_q_name"] + ", " + var_names["s_qd_name"] + ", " + var_names["s_qdd_name"] + ", "
    id_so_code_middle = self.gen_insert_helpers_function_call()
    if self.robot.floating_base:
        # Floating-base inner signature includes the gravity-shim spill pointer.
        id_so_code_end = var_names["s_temp_name"] + ", " + var_names["s_temp_spill_name"] + ", " + var_names["gravity_name"] + ");"
    else:
        id_so_code_end = var_names["s_temp_name"] + ", " + var_names["gravity_name"] + ");"
    if use_thread_group:
        id_so_code_start = id_so_code_start.replace("(","(tgrp, ")
    id_so_code = id_so_code_start + id_so_code_middle + id_so_code_end
    self.gen_add_code_line(id_so_code)

def idsva_so_parent_topology_needs_reference_order_output_repair(parent_ids):
    """
    Return true when moving-joint fanout requires the reference-order repair.

    Fanout directly from the fixed world/base parent is treated as independent
    root chains and can keep the optimized assembly path. Fanout below a moving
    joint can create order-sensitive duplicate tensor writes and needs repair.
    """
    direct_child_counts = {}
    for parent_jid in parent_ids:
        if parent_jid < 0:
            continue
        direct_child_counts[parent_jid] = direct_child_counts.get(parent_jid, 0) + 1
    return any(count > 1 for count in direct_child_counts.values())

def idsva_so_needs_reference_order_output_repair(self):
    """
    Returns true for fixed-base topologies where final IDSVA-SO tensor writes
    are order-sensitive because at least one moving joint has multiple direct
    children. Serial chains and base-rooted independent chain forests keep the
    optimized parallel tensor assembly.
    """
    if self.robot.is_serial_chain():
        return False
    parent_ids = [self.robot.get_parent_id(jid) for jid in range(self.robot.get_num_joints())]
    return idsva_so_parent_topology_needs_reference_order_output_repair(parent_ids)

def gen_idsva_so_reference_order_output_repair(self, use_thread_group = False):
    """
    Emits a serial final tensor assembly pass that mirrors RBDReference.idsva_so.

    The preceding generated code computes all reusable intermediates in
    parallel. The final second-order tensors, however, have many symmetry and
    duplicate-write relationships. Those writes are order-dependent for
    branched topologies, so replay the reference loop order with one thread.
    """
    num_bodies = self.robot.get_num_bodies()
    st_start = [0]
    st_values = []
    succ_start = [0]
    succ_values = []
    anc_start = [0]
    anc_values = []
    for jid in range(num_bodies):
        subtree = list(self.robot.get_subtree_by_id(jid))
        successors = [st_j for st_j in subtree if st_j != jid]
        ancestors = list(self.robot.get_ancestors_by_id(jid))
        ancestors.insert(0, jid)
        ancestors = ancestors[::-1]

        st_values.extend(subtree)
        st_start.append(len(st_values))
        succ_values.extend(successors)
        succ_start.append(len(succ_values))
        anc_values.extend(ancestors)
        anc_start.append(len(anc_values))

    def int_array(values):
        if values:
            return ", ".join(map(str, values))
        return "0"

    self.gen_add_sync(use_thread_group)
    self.gen_add_code_line("\n\n")
    self.gen_add_code_line("// Reference-order final IDSVA-SO tensor assembly")
    self.gen_add_code_line(f"static const int idsva_ref_st_start[] = {{ {int_array(st_start)} }};")
    self.gen_add_code_line(f"static const int idsva_ref_st_values[] = {{ {int_array(st_values)} }};")
    self.gen_add_code_line(f"static const int idsva_ref_succ_start[] = {{ {int_array(succ_start)} }};")
    self.gen_add_code_line(f"static const int idsva_ref_succ_values[] = {{ {int_array(succ_values)} }};")
    self.gen_add_code_line(f"static const int idsva_ref_anc_start[] = {{ {int_array(anc_start)} }};")
    self.gen_add_code_line(f"static const int idsva_ref_anc_values[] = {{ {int_array(anc_values)} }};")
    self.gen_add_code_line("if (threadIdx.x == 0 && threadIdx.y == 0 && threadIdx.z == 0) {", True)
    self.gen_add_code_line("for (int out_idx = 0; out_idx < SECOND_ORDER_TENSOR_SIZE; ++out_idx) s_idsva_so[out_idx] = static_cast<T>(0);")
    self.gen_add_code_line("T rt1[36], rt2[36], rt3[36], rt4[36], rt5[36], rt6[36], rt7[36], rt8[36], rt9[36];")
    self.gen_add_code_line("T rp1[6], rp2[6], rp3[6], rp4[6], rp5[6], rp6[6];")
    self.gen_add_code_line("for (int jid = NUM_BODIES - 1; jid >= 0; --jid) {", True)
    self.gen_add_code_line("int st_begin = idsva_ref_st_start[jid];")
    self.gen_add_code_line("int st_end = idsva_ref_st_start[jid + 1];")
    self.gen_add_code_line("int succ_begin = idsva_ref_succ_start[jid];")
    self.gen_add_code_line("int succ_end = idsva_ref_succ_start[jid + 1];")
    self.gen_add_code_line("int anc_begin = idsva_ref_anc_start[jid];")
    self.gen_add_code_line("int anc_end = idsva_ref_anc_start[jid + 1];")
    self.gen_add_code_line("for (int anc_pos = anc_begin; anc_pos < anc_end; ++anc_pos) {", True)
    self.gen_add_code_line("int ancestor_j = idsva_ref_anc_values[anc_pos];")
    self.gen_add_code_line("for (int idx = 0; idx < 36; ++idx) {", True)
    self.gen_add_code_line("int row = idx % 6;")
    self.gen_add_code_line("int col = idx / 6;")
    self.gen_add_code_line("rt1[idx] = S[jid*6 + row] * psid[ancestor_j*6 + col];")
    self.gen_add_code_line("rt2[idx] = S[jid*6 + row] * S[ancestor_j*6 + col];")
    self.gen_add_code_line("rt3[idx] = psid[jid*6 + row] * psid[ancestor_j*6 + col];")
    self.gen_add_code_line("rt4[idx] = S[jid*6 + row] * psidd[ancestor_j*6 + col];")
    self.gen_add_code_line("rt5[idx] = S[jid*6 + row] * psid_Sd[ancestor_j*6 + col];")
    self.gen_add_code_line("rt6[idx] = S[ancestor_j*6 + row] * psid[jid*6 + col];")
    self.gen_add_code_line("rt7[idx] = S[ancestor_j*6 + row] * psidd[jid*6 + col];")
    self.gen_add_code_line("rt8[idx] = S[ancestor_j*6 + row] * S[jid*6 + col];")
    self.gen_add_code_line("rt9[idx] = S[ancestor_j*6 + row] * psid_Sd[jid*6 + col];")
    self.gen_add_end_control_flow()
    self.gen_add_code_line("for (int row = 0; row < 6; ++row) {", True)
    self.gen_add_code_line("rp1[row] = crm_mul<T>(row, &psid[ancestor_j*6], &S[jid*6]);")
    self.gen_add_code_line("rp2[row] = crm_mul<T>(row, &psidd[ancestor_j*6], &S[jid*6]);")
    self.gen_add_code_line("rp3[row] = crm_mul<T>(row, &S[ancestor_j*6], &S[jid*6]);")
    self.gen_add_code_line("rp4[row] = crm_mul<T>(row, &psid_Sd[ancestor_j*6], &S[jid*6]) - static_cast<T>(2) * crm_mul<T>(row, &psid[jid*6], &S[ancestor_j*6]);")
    self.gen_add_code_line("rp5[row] = crm_mul<T>(row, &S[jid*6], &S[ancestor_j*6]);")
    self.gen_add_code_line("rp6[row] = dot_prod<T, 6, 1, 1>(&IC_S[jid*6], &crm_S[ancestor_j*36 + row*6]) + dot_prod<T, 6, 1, 1>(&S[ancestor_j*6], &crf_S_IC[jid*36 + row*6]);")
    self.gen_add_end_control_flow()
    self.gen_add_code_line("for (int st_pos = st_begin; st_pos < st_end; ++st_pos) {", True)
    self.gen_add_code_line("int st_j = idsva_ref_st_values[st_pos];")
    self.gen_add_code_line("d2tau_dq2[st_j*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + jid*SECOND_ORDER_COORDS + ancestor_j] = -dot_prod<T, 36, 1, 1>(rt3, &D3[st_j*36]) - dot_prod<T, 6, 1, 1>(rp1, &T2[st_j*6]) + dot_prod<T, 6, 1, 1>(rp2, &T1[st_j*6]);")
    self.gen_add_code_line("d2tau_dvdq[st_j*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + jid*SECOND_ORDER_COORDS + ancestor_j] = -dot_prod<T, 36, 1, 1>(rt1, &D3[st_j*36]);")
    self.gen_add_end_control_flow()
    self.gen_add_code_line("if (ancestor_j < jid) {", True)
    self.gen_add_code_line("for (int st_pos = st_begin; st_pos < st_end; ++st_pos) {", True)
    self.gen_add_code_line("int st_j = idsva_ref_st_values[st_pos];")
    self.gen_add_code_line("d2tau_dq2[st_j*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + ancestor_j*SECOND_ORDER_COORDS + jid] = d2tau_dq2[st_j*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + jid*SECOND_ORDER_COORDS + ancestor_j];")
    self.gen_add_code_line("d2tau_dqd2[st_j*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + ancestor_j*SECOND_ORDER_COORDS + jid] = -dot_prod<T, 36, 1, 1>(rt2, &D3[st_j*36]);")
    self.gen_add_code_line("d2tau_dqd2[st_j*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + jid*SECOND_ORDER_COORDS + ancestor_j] = d2tau_dqd2[st_j*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + ancestor_j*SECOND_ORDER_COORDS + jid];")
    self.gen_add_code_line("d2tau_dvdq[st_j*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + ancestor_j*SECOND_ORDER_COORDS + jid] = -dot_prod<T, 36, 1, 1>(rt6, &D3[st_j*36]) - dot_prod<T, 6, 1, 1>(rp3, &T2[st_j*6]) + dot_prod<T, 6, 1, 1>(rp4, &T1[st_j*6]);")
    self.gen_add_code_line("d2tau_dq2[ancestor_j*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + st_j*SECOND_ORDER_COORDS + jid] = dot_prod<T, 36, 1, 1>(rt6, &D2[st_j*36]) + dot_prod<T, 36, 1, 1>(rt7, &D1[st_j*36]) - dot_prod<T, 6, 1, 1>(rp5, &T3[st_j*6]);")
    self.gen_add_code_line("d2tau_dvdq[ancestor_j*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + st_j*SECOND_ORDER_COORDS + jid] = dot_prod<T, 36, 1, 1>(rt6, &D3[st_j*36]) - dot_prod<T, 6, 1, 1>(rp5, &T4[st_j*6]);")
    self.gen_add_code_line("dM_dq[ancestor_j*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + st_j*SECOND_ORDER_COORDS + jid] = dot_prod<T, 36, 1, 1>(rt8, &D4[st_j*36]);")
    self.gen_add_code_line("dM_dq[st_j*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + ancestor_j*SECOND_ORDER_COORDS + jid] = dM_dq[ancestor_j*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + st_j*SECOND_ORDER_COORDS + jid];")
    self.gen_add_end_control_flow()
    self.gen_add_code_line("d2tau_dqd2[ancestor_j*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + jid*SECOND_ORDER_COORDS + jid] = dot_prod<T, 6, 1, 1>(rp6, &S[jid*6]);")
    self.gen_add_code_line("for (int succ_pos = succ_begin; succ_pos < succ_end; ++succ_pos) {", True)
    self.gen_add_code_line("int succ_j = idsva_ref_succ_values[succ_pos];")
    self.gen_add_code_line("d2tau_dqd2[ancestor_j*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + succ_j*SECOND_ORDER_COORDS + jid] = dot_prod<T, 36, 1, 1>(rt8, &D3[succ_j*36]);")
    self.gen_add_code_line("d2tau_dqd2[ancestor_j*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + jid*SECOND_ORDER_COORDS + succ_j] = d2tau_dqd2[ancestor_j*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + succ_j*SECOND_ORDER_COORDS + jid];")
    self.gen_add_code_line("d2tau_dvdq[ancestor_j*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + jid*SECOND_ORDER_COORDS + succ_j] = dot_prod<T, 36, 1, 1>(rt8, &D2[succ_j*36]) + dot_prod<T, 36, 1, 1>(rt9, &D1[succ_j*36]);")
    self.gen_add_code_line("d2tau_dq2[ancestor_j*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + jid*SECOND_ORDER_COORDS + succ_j] = d2tau_dq2[ancestor_j*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + succ_j*SECOND_ORDER_COORDS + jid];")
    self.gen_add_end_control_flow()
    self.gen_add_end_control_flow()
    self.gen_add_code_line("for (int succ_pos = succ_begin; succ_pos < succ_end; ++succ_pos) {", True)
    self.gen_add_code_line("int succ_j = idsva_ref_succ_values[succ_pos];")
    self.gen_add_code_line("d2tau_dq2[jid*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + ancestor_j*SECOND_ORDER_COORDS + succ_j] = dot_prod<T, 36, 1, 1>(rt1, &D2[succ_j*36]) + dot_prod<T, 36, 1, 1>(rt4, &D1[succ_j*36]);")
    self.gen_add_code_line("d2tau_dqd2[jid*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + ancestor_j*SECOND_ORDER_COORDS + succ_j] = dot_prod<T, 36, 1, 1>(rt2, &D3[succ_j*36]);")
    self.gen_add_code_line("d2tau_dqd2[jid*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + succ_j*SECOND_ORDER_COORDS + ancestor_j] = d2tau_dqd2[jid*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + ancestor_j*SECOND_ORDER_COORDS + succ_j];")
    self.gen_add_code_line("d2tau_dvdq[jid*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + succ_j*SECOND_ORDER_COORDS + ancestor_j] = dot_prod<T, 36, 1, 1>(rt1, &D3[succ_j*36]);")
    self.gen_add_code_line("d2tau_dq2[jid*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + succ_j*SECOND_ORDER_COORDS + ancestor_j] = d2tau_dq2[jid*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + ancestor_j*SECOND_ORDER_COORDS + succ_j];")
    self.gen_add_code_line("d2tau_dvdq[jid*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + ancestor_j*SECOND_ORDER_COORDS + succ_j] = dot_prod<T, 36, 1, 1>(rt2, &D2[succ_j*36]) + dot_prod<T, 36, 1, 1>(rt5, &D1[succ_j*36]);")
    self.gen_add_code_line("dM_dq[ancestor_j*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + jid*SECOND_ORDER_COORDS + succ_j] = dot_prod<T, 36, 1, 1>(rt8, &D1[succ_j*36]);")
    self.gen_add_code_line("dM_dq[jid*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + ancestor_j*SECOND_ORDER_COORDS + succ_j] = dM_dq[ancestor_j*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + jid*SECOND_ORDER_COORDS + succ_j];")
    self.gen_add_end_control_flow()
    self.gen_add_code_line("if (ancestor_j == jid) {", True)
    self.gen_add_code_line("for (int st_pos = st_begin; st_pos < st_end; ++st_pos) {", True)
    self.gen_add_code_line("int st_j = idsva_ref_st_values[st_pos];")
    self.gen_add_code_line("d2tau_dqd2[st_j*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + jid*SECOND_ORDER_COORDS + ancestor_j] = -dot_prod<T, 36, 1, 1>(rt2, &D1[st_j*36]);")
    self.gen_add_end_control_flow()
    self.gen_add_end_control_flow()
    self.gen_add_end_control_flow()
    self.gen_add_end_control_flow()
    self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)

def gen_idsva_so_floating_reference_inner(self, use_thread_group = False, use_qdd_input = False):
    """
    Emits a floating-base diagnostic IDSVA-SO path with explicit body/velocity
    split memory. Fixed-base keeps the optimized generator path below.
    """
    NV = self.robot.get_num_vel()
    num_bodies = self.robot.get_num_bodies()
    metadata = _idsva_so_floating_velocity_metadata(self.robot)
    parent_ids = [self.robot.get_parent_id(body_id) for body_id in range(num_bodies)]

    func_params = ["s_idsva_so is a pointer to memory for the final result of size 4*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS = " + str(4*NV**3), \
                   "s_q is the vector of joint positions", \
                   "s_qd is the vector of joint velocities", \
                   "s_qdd is the vector of joint accelerations", \
                   "s_temp is a pointer to helper shared memory of size  = " + \
                            str(self.gen_idsva_so_inner_temp_mem_size()), \
                   "gravity is the gravity constant"]
    func_def_start = "void idsva_so_inner(T *s_idsva_so, const T *s_q, const T *s_qd, T *s_qdd, "
    func_def_end = "T *s_temp, const T gravity) {"
    func_params.insert(-1, "s_temp_spill is a pointer to global-memory scratch (per-timestep) of size = " + \
                            str(gen_floating_gravity_d2tau_dq_spill_count(self)) + " floats")
    func_def_end = "T *s_temp, T *s_temp_spill, const T gravity) {"
    func_def_start, func_params = self.gen_insert_helpers_func_def_params(func_def_start, func_params, -2)
    func_notes = [
        "Floating diagnostic path: body-indexed spatial state plus packed velocity-indexed derivative columns.",
        "d2tau_dq is assembled analytically in velocity-coordinate tensor space.",
    ]
    if use_thread_group:
        func_def_start = func_def_start.replace("(", "(cgrps::thread_group tgrp, ")
        func_params.insert(0,"tgrp is the handle to the thread_group running this function")
    func_def = func_def_start + func_def_end

    self.gen_add_func_doc("Computes floating-base second-order inverse dynamics diagnostics",func_notes,func_params,None)
    self.gen_add_code_line("template <typename T>")
    self.gen_add_code_line("__device__")
    self.gen_add_code_line(func_def, True)

    self.gen_add_code_lines([
        "// Floating IDSVA-SO split-memory layout.",
        "T *I = s_XImats + XIMAT_SIZE*NUM_BODIES;",
        "T *Xdown = s_temp;",
        "T *Xup = Xdown + 36*NUM_BODIES;",
        "T *IC = Xup + 36*NUM_BODIES;",
        "T *I_Xup = IC + 36*NUM_BODIES;",
        "T *BC = I_Xup + 36*NUM_BODIES;",
        "T *crm_v = BC + 36*NUM_BODIES;",
        "T *crf_v = crm_v + 36*NUM_BODIES;",
        "T *vJ = crf_v + 36*NUM_BODIES;",
        "T *v = vJ + 6*NUM_BODIES;",
        "T *aJ = v + 6*NUM_BODIES;",
        "T *a = aJ + 6*NUM_BODIES;",
        "T *f = a + 6*NUM_BODIES;",
        "T *IC_v = f + 6*NUM_BODIES;",
        "T *a_world = IC_v + 6*NUM_BODIES;",
        "T *S_vel = a_world + 6;",
        "T *Sd_vel = S_vel + 6*NUM_VEL;",
        "T *psid_vel = Sd_vel + 6*NUM_VEL;",
        "T *psidd_vel = psid_vel + 6*NUM_VEL;",
        "T *psid_Sd_vel = psidd_vel + 6*NUM_VEL;",
        "T *IC_S = psid_Sd_vel + 6*NUM_VEL;",
        "T *IC_psid = IC_S + 6*NUM_VEL;",
        "T *ICT_S = IC_psid + 6*NUM_VEL;",
        "T *T1 = IC_S;",
        "T *T2 = ICT_S + 6*NUM_VEL;",
        "T *T3 = T2 + 6*NUM_VEL;",
        "T *T4 = T3 + 6*NUM_VEL;",
        "T *crm_S = T4 + 6*NUM_VEL;",
        "T *crf_S = crm_S + 36*NUM_VEL;",
        "T *crm_psid = crf_S + 36*NUM_VEL;",
        "T *crf_psid = crm_psid + 36*NUM_VEL;",
        "T *icrf_f = crf_psid + 36*NUM_VEL;",
        "T *B_IC_S = icrf_f + 36*NUM_BODIES;",
        "T *D1 = B_IC_S + 36*NUM_VEL;",
        "T *D2 = D1 + 36*NUM_VEL;",
        "T *D3 = B_IC_S;",
        "T *D4 = D2 + 36*NUM_VEL;",
        "T *crf_S_IC = D4 + 36*NUM_VEL;",
        "T *d2tau_dq2 = s_idsva_so;",
        "T *d2tau_dqd2 = d2tau_dq2 + SECOND_ORDER_COORDS*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS;",
        "T *d2tau_dvdq = d2tau_dqd2 + SECOND_ORDER_COORDS*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS;",
        "T *dM_dq = d2tau_dvdq + SECOND_ORDER_COORDS*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS;",
        "",
        f"static const int idsva_float_parent[] = {{ {_idsva_so_int_array(parent_ids)} }};",
        f"static const int body_v_start[] = {{ {_idsva_so_int_array(metadata['body_v_start'])} }};",
        f"static const int body_v_index[] = {{ {_idsva_so_int_array(metadata['body_v_index'])} }};",
        f"static const int vel_to_body[] = {{ {_idsva_so_int_array(metadata['vel_to_body'])} }};",
        f"static const int vel_s_index[] = {{ {_idsva_so_int_array(metadata['vel_s_index'])} }};",
        f"static const int vel_s_sign[] = {{ {_idsva_so_int_array(metadata['vel_s_sign'])} }};",
        f"static const int subtree_v_start[] = {{ {_idsva_so_int_array(metadata['subtree_v_start'])} }};",
        f"static const int subtree_v_index[] = {{ {_idsva_so_int_array(metadata['subtree_v_index'])} }};",
        f"static const int successor_v_start[] = {{ {_idsva_so_int_array(metadata['successor_v_start'])} }};",
        f"static const int successor_v_index[] = {{ {_idsva_so_int_array(metadata['successor_v_index'])} }};",
        f"static const int ancestor_body_start[] = {{ {_idsva_so_int_array(metadata['ancestor_body_start'])} }};",
        f"static const int ancestor_body_index[] = {{ {_idsva_so_int_array(metadata['ancestor_body_index'])} }};",
        "",
    ])

    self.gen_add_code_line("if (threadIdx.x == 0 && threadIdx.y == 0 && threadIdx.z == 0) {", True)
    self.gen_add_code_line("for (int out_idx = 0; out_idx < SECOND_ORDER_TENSOR_SIZE; ++out_idx) s_idsva_so[out_idx] = static_cast<T>(0);")
    self.gen_add_code_line("// Floating-base: run main sweep with gravity = 0; gravity Hessian added below by")
    self.gen_add_code_line("// `gen_floating_gravity_d2tau_dq_lie_inline` (mirrors Python idsva_gravity = 0.0 + shim).")
    self.gen_add_code_line("a_world[0] = static_cast<T>(0); a_world[1] = static_cast<T>(0); a_world[2] = static_cast<T>(0);")
    self.gen_add_code_line("a_world[3] = static_cast<T>(0); a_world[4] = static_cast<T>(0); a_world[5] = static_cast<T>(0);")

    self.gen_add_code_line("// Compute accumulated Xup transforms.")
    self.gen_add_code_line("for (int jid = 0; jid < NUM_BODIES; ++jid) {", True)
    self.gen_add_code_line("int parent = idsva_float_parent[jid];")
    self.gen_add_code_line("for (int idx = 0; idx < 36; ++idx) {", True)
    self.gen_add_code_line("int row = idx % 6;")
    self.gen_add_code_line("int col = idx / 6;")
    self.gen_add_code_line("if (parent < 0) {", True)
    self.gen_add_code_line("Xup[jid*36 + idx] = s_XImats[jid*36 + idx];")
    self.gen_add_end_control_flow()
    self.gen_add_code_line("else {", True)
    self.gen_add_code_line("T acc = static_cast<T>(0);")
    self.gen_add_code_line("for (int kk = 0; kk < 6; ++kk) acc += s_XImats[jid*36 + row + 6*kk] * Xup[parent*36 + kk + 6*col];")
    self.gen_add_code_line("Xup[jid*36 + idx] = acc;")
    self.gen_add_end_control_flow()
    self.gen_add_end_control_flow()
    self.gen_add_end_control_flow()

    self.gen_add_code_line("// Compute IC = Xup^T * I * Xup.")
    self.gen_add_code_line("for (int jid = 0; jid < NUM_BODIES; ++jid) {", True)
    self.gen_add_code_line("for (int idx = 0; idx < 36; ++idx) {", True)
    self.gen_add_code_line("int row = idx % 6; int col = idx / 6;")
    self.gen_add_code_line("T acc = static_cast<T>(0);")
    self.gen_add_code_line("for (int kk = 0; kk < 6; ++kk) acc += I[jid*36 + row + 6*kk] * Xup[jid*36 + kk + 6*col];")
    self.gen_add_code_line("I_Xup[jid*36 + idx] = acc;")
    self.gen_add_end_control_flow()
    self.gen_add_code_line("for (int idx = 0; idx < 36; ++idx) {", True)
    self.gen_add_code_line("int row = idx % 6; int col = idx / 6;")
    self.gen_add_code_line("T acc = static_cast<T>(0);")
    self.gen_add_code_line("for (int kk = 0; kk < 6; ++kk) acc += Xup[jid*36 + kk + 6*row] * I_Xup[jid*36 + kk + 6*col];")
    self.gen_add_code_line("IC[jid*36 + idx] = acc;")
    self.gen_add_end_control_flow()
    self.gen_add_end_control_flow()

    self.gen_add_code_line("// Compute Xdown using the same spatial-transform inverse pattern as the fixed path.")
    self.gen_add_code_line("for (int flat = 0; flat < 36*NUM_BODIES; ++flat) Xdown[flat] = static_cast<T>(0);")
    self.gen_add_code_line("for (int flat = 0; flat < 36*NUM_BODIES; ++flat) {", True)
    self.gen_add_code_line("int idx = flat % 36;")
    self.gen_add_code_line("int sub_idx = idx % 18;")
    self.gen_add_code_line("if (idx % 18 == 1 || idx % 18 == 4 || idx % 18 == 8 || idx % 18 == 11) {", True)
    self.gen_add_code_line("Xdown[flat] = Xup[flat + 5];")
    self.gen_add_code_line("Xdown[flat + 5] = Xup[flat];")
    self.gen_add_end_control_flow()
    self.gen_add_code_line("else if (idx % 18 == 2 || idx % 18 == 5) {", True)
    self.gen_add_code_line("Xdown[flat] = Xup[flat + 10];")
    self.gen_add_code_line("Xdown[flat + 10] = Xup[flat];")
    self.gen_add_end_control_flow()
    self.gen_add_code_line("else if (sub_idx != 6 && sub_idx != 9 && sub_idx != 13 && sub_idx != 16 && sub_idx != 12 && sub_idx != 15) {", True)
    self.gen_add_code_line("Xdown[flat] = Xup[flat];")
    self.gen_add_end_control_flow()
    self.gen_add_end_control_flow()

    self.gen_add_code_line("// Transform each velocity-coordinate S column.")
    self.gen_add_code_line("for (int vel = 0; vel < NUM_VEL; ++vel) {", True)
    self.gen_add_code_line("int jid = vel_to_body[vel];")
    self.gen_add_code_line("int s_col = vel_s_index[vel];")
    self.gen_add_code_line("T s_sign = static_cast<T>(vel_s_sign[vel]);")
    self.gen_add_code_line("for (int row = 0; row < 6; ++row) S_vel[vel*6 + row] = s_sign * Xdown[jid*36 + s_col*6 + row];")
    self.gen_add_end_control_flow()

    self.gen_add_code_line("// Forward pass in body order.")
    self.gen_add_code_line("for (int jid = 0; jid < NUM_BODIES; ++jid) {", True)
    self.gen_add_code_line("int parent = idsva_float_parent[jid];")
    self.gen_add_code_line("for (int row = 0; row < 6; ++row) {", True)
    self.gen_add_code_line("vJ[jid*6 + row] = static_cast<T>(0);")
    self.gen_add_code_line("aJ[jid*6 + row] = static_cast<T>(0);")
    self.gen_add_code_line("for (int pos = body_v_start[jid]; pos < body_v_start[jid + 1]; ++pos) {", True)
    self.gen_add_code_line("int vel = body_v_index[pos];")
    self.gen_add_code_line("vJ[jid*6 + row] += S_vel[vel*6 + row] * s_qd[vel];")
    self.gen_add_code_line("aJ[jid*6 + row] += S_vel[vel*6 + row] * s_qdd[vel];")
    self.gen_add_end_control_flow()
    self.gen_add_code_line("if (parent < 0) {", True)
    self.gen_add_code_line("v[jid*6 + row] = static_cast<T>(0);")
    self.gen_add_code_line("a[jid*6 + row] = dot_prod<T, 6, 6, 1>(&Xdown[jid*36 + row], a_world);")
    self.gen_add_end_control_flow()
    self.gen_add_code_line("else {", True)
    self.gen_add_code_line("v[jid*6 + row] = v[parent*6 + row];")
    self.gen_add_code_line("a[jid*6 + row] = a[parent*6 + row];")
    self.gen_add_end_control_flow()
    self.gen_add_end_control_flow()
    self.gen_add_code_line("for (int row = 0; row < 6; ++row) aJ[jid*6 + row] += crm_mul<T>(row, &v[jid*6], &vJ[jid*6]);")
    self.gen_add_code_line("for (int pos = body_v_start[jid]; pos < body_v_start[jid + 1]; ++pos) {", True)
    self.gen_add_code_line("int vel = body_v_index[pos];")
    self.gen_add_code_line("for (int row = 0; row < 6; ++row) {", True)
    self.gen_add_code_line("psid_vel[vel*6 + row] = crm_mul<T>(row, &v[jid*6], &S_vel[vel*6]);")
    self.gen_add_code_line("psidd_vel[vel*6 + row] = crm_mul<T>(row, &a[jid*6], &S_vel[vel*6]) + crm_mul<T>(row, &v[jid*6], &psid_vel[vel*6]);")
    self.gen_add_end_control_flow()
    self.gen_add_end_control_flow()
    self.gen_add_code_line("for (int row = 0; row < 6; ++row) {", True)
    self.gen_add_code_line("v[jid*6 + row] += vJ[jid*6 + row];")
    self.gen_add_code_line("a[jid*6 + row] += aJ[jid*6 + row];")
    self.gen_add_end_control_flow()
    self.gen_add_code_line("for (int pos = body_v_start[jid]; pos < body_v_start[jid + 1]; ++pos) {", True)
    self.gen_add_code_line("int vel = body_v_index[pos];")
    self.gen_add_code_line("for (int row = 0; row < 6; ++row) Sd_vel[vel*6 + row] = crm_mul<T>(row, &v[jid*6], &S_vel[vel*6]);")
    self.gen_add_end_control_flow()
    self.gen_add_code_line("for (int idx = 0; idx < 36; ++idx) {", True)
    self.gen_add_code_line("int row = idx % 6; int col = idx / 6;")
    self.gen_add_code_line("crm_v[jid*36 + idx] = crm<T>(idx, &v[jid*6]);")
    self.gen_add_code_line("crf_v[jid*36 + row*6 + col] = -crm<T>(idx, &v[jid*6]);")
    self.gen_add_end_control_flow()
    self.gen_add_code_line("for (int row = 0; row < 6; ++row) IC_v[jid*6 + row] = dot_prod<T, 6, 6, 1>(&IC[jid*36 + row], &v[jid*6]);")
    self.gen_add_code_line("for (int idx = 0; idx < 36; ++idx) {", True)
    self.gen_add_code_line("int row = idx % 6; int col = idx / 6;")
    self.gen_add_code_line("BC[jid*36 + idx] = dot_prod<T, 6, 6, 1>(&crf_v[jid*36 + row], &IC[jid*36 + col*6]) + icrf<T>(idx, &IC_v[jid*6]) - dot_prod<T, 6, 6, 1>(&IC[jid*36 + row], &crm_v[jid*36 + col*6]);")
    self.gen_add_end_control_flow()
    self.gen_add_code_line("for (int row = 0; row < 6; ++row) f[jid*6 + row] = dot_prod<T, 6, 6, 1>(&IC[jid*36 + row], &a[jid*6]) + dot_prod<T, 6, 6, 1>(&crf_v[jid*36 + row], &IC_v[jid*6]);")
    self.gen_add_end_control_flow()

    self.gen_add_code_line("// Backward accumulation of body-indexed composite quantities.")
    self.gen_add_code_line("for (int jid = NUM_BODIES - 1; jid >= 0; --jid) {", True)
    self.gen_add_code_line("int parent = idsva_float_parent[jid];")
    self.gen_add_code_line("if (parent >= 0) {", True)
    self.gen_add_code_line("for (int idx = 0; idx < 36; ++idx) { IC[parent*36 + idx] += IC[jid*36 + idx]; BC[parent*36 + idx] += BC[jid*36 + idx]; }")
    self.gen_add_code_line("for (int row = 0; row < 6; ++row) f[parent*6 + row] += f[jid*6 + row];")
    self.gen_add_end_control_flow()
    self.gen_add_end_control_flow()

    self.gen_add_code_line("// Build velocity-indexed T and D intermediates.")
    self.gen_add_code_line("for (int jid = NUM_BODIES - 1; jid >= 0; --jid) {", True)
    self.gen_add_code_line("for (int pos = body_v_start[jid]; pos < body_v_start[jid + 1]; ++pos) {", True)
    self.gen_add_code_line("int vel = body_v_index[pos];")
    self.gen_add_code_line("for (int idx = 0; idx < 36; ++idx) {", True)
    self.gen_add_code_line("int row = idx % 6; int col = idx / 6;")
    self.gen_add_code_line("crm_S[vel*36 + idx] = crm<T>(idx, &S_vel[vel*6]);")
    self.gen_add_code_line("crf_S[vel*36 + row*6 + col] = -crm<T>(idx, &S_vel[vel*6]);")
    self.gen_add_code_line("crm_psid[vel*36 + idx] = crm<T>(idx, &psid_vel[vel*6]);")
    self.gen_add_code_line("crf_psid[vel*36 + row*6 + col] = -crm<T>(idx, &psid_vel[vel*6]);")
    self.gen_add_end_control_flow()
    self.gen_add_code_line("for (int row = 0; row < 6; ++row) {", True)
    self.gen_add_code_line("IC_S[vel*6 + row] = dot_prod<T, 6, 6, 1>(&IC[jid*36 + row], &S_vel[vel*6]);")
    self.gen_add_code_line("IC_psid[vel*6 + row] = dot_prod<T, 6, 6, 1>(&IC[jid*36 + row], &psid_vel[vel*6]);")
    self.gen_add_code_line("psid_Sd_vel[vel*6 + row] = psid_vel[vel*6 + row] + Sd_vel[vel*6 + row];")
    self.gen_add_code_line("ICT_S[vel*6 + row] = dot_prod<T, 6, 1, 1>(&IC[jid*36 + row*6], &S_vel[vel*6]);")
    self.gen_add_end_control_flow()
    self.gen_add_code_line("for (int idx = 0; idx < 36; ++idx) {", True)
    self.gen_add_code_line("int row = idx % 6; int col = idx / 6;")
    self.gen_add_code_line("icrf_f[jid*36 + idx] = icrf<T>(idx, &f[jid*6]);")
    self.gen_add_code_line("B_IC_S[vel*36 + idx] = dot_prod<T, 6, 6, 1>(&crf_S[vel*36 + row], &IC[jid*36 + col*6]) + icrf<T>(idx, &IC_S[vel*6]) - dot_prod<T, 6, 6, 1>(&IC[jid*36 + row], &crm_S[vel*36 + col*6]);")
    self.gen_add_code_line("D2[vel*36 + idx] = dot_prod<T, 6, 6, 1>(&crf_psid[vel*36 + row], &IC[jid*36 + col*6]) + icrf<T>(idx, &IC_psid[vel*6]) - dot_prod<T, 6, 6, 1>(&IC[jid*36 + row], &crm_psid[vel*36 + col*6]);")
    self.gen_add_code_line("// RBDReference stores D1 with NumPy default flatten() order, unlike D2/D3/D4.")
    self.gen_add_code_line("D1[vel*36 + row*6 + col] = dot_prod<T, 6, 6, 1>(&crf_S[vel*36 + row], &IC[jid*36 + col*6]) - dot_prod<T, 6, 6, 1>(&IC[jid*36 + row], &crm_S[vel*36 + col*6]);")
    self.gen_add_code_line("D2[vel*36 + idx] += dot_prod<T, 6, 6, 1>(&crf_S[vel*36 + row], &BC[jid*36 + col*6]) - dot_prod<T, 6, 6, 1>(&BC[jid*36 + row], &crm_S[vel*36 + col*6]);")
    self.gen_add_code_line("D4[vel*36 + idx] = icrf<T>(idx, &ICT_S[vel*6]);")
    self.gen_add_code_line("crf_S_IC[vel*36 + idx] = dot_prod<T, 6, 6, 1>(&crf_S[vel*36 + row], &IC[jid*36 + col*6]);")
    self.gen_add_end_control_flow()
    self.gen_add_code_line("for (int row = 0; row < 6; ++row) {", True)
    self.gen_add_code_line("T2[vel*6 + row] = -dot_prod<T, 6, 1, 1>(&BC[jid*36 + row*6], &S_vel[vel*6]);")
    self.gen_add_code_line("T3[vel*6 + row] = dot_prod<T, 6, 6, 1>(&BC[jid*36 + row], &psid_vel[vel*6]) + dot_prod<T, 6, 6, 1>(&IC[jid*36 + row], &psidd_vel[vel*6]) + dot_prod<T, 6, 6, 1>(&icrf_f[jid*36 + row], &S_vel[vel*6]);")
    self.gen_add_code_line("T4[vel*6 + row] = dot_prod<T, 6, 6, 1>(&BC[jid*36 + row], &S_vel[vel*6]) + dot_prod<T, 6, 6, 1>(&IC[jid*36 + row], &psid_Sd_vel[vel*6]);")
    self.gen_add_end_control_flow()
    self.gen_add_end_control_flow()
    self.gen_add_end_control_flow()

    self.gen_add_code_line("// Reference-order velocity-indexed tensor assembly.")
    self.gen_add_code_line("T rt1[36], rt2[36], rt3[36], rt4[36], rt5[36], rt6[36], rt7[36], rt8[36], rt9[36];")
    self.gen_add_code_line("T rp1[6], rp2[6], rp3[6], rp4[6], rp5[6], rp6[6];")
    self.gen_add_code_line("for (int jid = NUM_BODIES - 1; jid >= 0; --jid) {", True)
    self.gen_add_code_line("int st_begin = subtree_v_start[jid]; int st_end = subtree_v_start[jid + 1];")
    self.gen_add_code_line("int succ_begin = successor_v_start[jid]; int succ_end = successor_v_start[jid + 1];")
    self.gen_add_code_line("int anc_begin = ancestor_body_start[jid]; int anc_end = ancestor_body_start[jid + 1];")
    self.gen_add_code_line("for (int dpos = body_v_start[jid]; dpos < body_v_start[jid + 1]; ++dpos) {", True)
    self.gen_add_code_line("int dd = body_v_index[dpos];")
    self.gen_add_code_line("for (int anc_pos = anc_begin; anc_pos < anc_end; ++anc_pos) {", True)
    self.gen_add_code_line("int ancestor_body = ancestor_body_index[anc_pos];")
    self.gen_add_code_line("for (int cpos = body_v_start[ancestor_body]; cpos < body_v_start[ancestor_body + 1]; ++cpos) {", True)
    self.gen_add_code_line("int cc = body_v_index[cpos];")
    self.gen_add_code_line("for (int idx = 0; idx < 36; ++idx) {", True)
    self.gen_add_code_line("int row = idx % 6; int col = idx / 6;")
    self.gen_add_code_line("rt1[idx] = S_vel[dd*6 + row] * psid_vel[cc*6 + col];")
    self.gen_add_code_line("rt2[idx] = S_vel[dd*6 + row] * S_vel[cc*6 + col];")
    self.gen_add_code_line("rt3[idx] = psid_vel[dd*6 + row] * psid_vel[cc*6 + col];")
    self.gen_add_code_line("rt4[idx] = S_vel[dd*6 + row] * psidd_vel[cc*6 + col];")
    self.gen_add_code_line("rt5[idx] = S_vel[dd*6 + row] * psid_Sd_vel[cc*6 + col];")
    self.gen_add_code_line("rt6[idx] = S_vel[cc*6 + row] * psid_vel[dd*6 + col];")
    self.gen_add_code_line("rt7[idx] = S_vel[cc*6 + row] * psidd_vel[dd*6 + col];")
    self.gen_add_code_line("rt8[idx] = S_vel[cc*6 + row] * S_vel[dd*6 + col];")
    self.gen_add_code_line("rt9[idx] = S_vel[cc*6 + row] * psid_Sd_vel[dd*6 + col];")
    self.gen_add_end_control_flow()
    self.gen_add_code_line("for (int row = 0; row < 6; ++row) {", True)
    self.gen_add_code_line("rp1[row] = crm_mul<T>(row, &psid_vel[cc*6], &S_vel[dd*6]);")
    self.gen_add_code_line("rp2[row] = crm_mul<T>(row, &psidd_vel[cc*6], &S_vel[dd*6]);")
    self.gen_add_code_line("rp3[row] = crm_mul<T>(row, &S_vel[cc*6], &S_vel[dd*6]);")
    self.gen_add_code_line("rp4[row] = crm_mul<T>(row, &psid_Sd_vel[cc*6], &S_vel[dd*6]) - static_cast<T>(2) * crm_mul<T>(row, &psid_vel[dd*6], &S_vel[cc*6]);")
    self.gen_add_code_line("rp5[row] = crm_mul<T>(row, &S_vel[dd*6], &S_vel[cc*6]);")
    self.gen_add_code_line("rp6[row] = dot_prod<T, 6, 1, 1>(&IC_S[dd*6], &crm_S[cc*36 + row*6]) + dot_prod<T, 6, 1, 1>(&S_vel[cc*6], &crf_S_IC[dd*36 + row*6]);")
    self.gen_add_end_control_flow()
    self.gen_add_code_line("for (int st_pos = st_begin; st_pos < st_end; ++st_pos) {", True)
    self.gen_add_code_line("int st_vel = subtree_v_index[st_pos];")
    self.gen_add_code_line("d2tau_dq2[st_vel*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + dd*SECOND_ORDER_COORDS + cc] = -dot_prod<T, 36, 1, 1>(rt3, &D3[st_vel*36]) - dot_prod<T, 6, 1, 1>(rp1, &T2[st_vel*6]) + dot_prod<T, 6, 1, 1>(rp2, &T1[st_vel*6]);")
    self.gen_add_code_line("d2tau_dvdq[st_vel*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + dd*SECOND_ORDER_COORDS + cc] = -dot_prod<T, 36, 1, 1>(rt1, &D3[st_vel*36]);")
    self.gen_add_end_control_flow()
    self.gen_add_code_line("if (ancestor_body < jid) {", True)
    self.gen_add_code_line("for (int st_pos = st_begin; st_pos < st_end; ++st_pos) {", True)
    self.gen_add_code_line("int st_vel = subtree_v_index[st_pos];")
    self.gen_add_code_line("d2tau_dq2[st_vel*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + cc*SECOND_ORDER_COORDS + dd] = d2tau_dq2[st_vel*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + dd*SECOND_ORDER_COORDS + cc];")
    self.gen_add_code_line("T dqd_val = -dot_prod<T, 36, 1, 1>(rt2, &D3[st_vel*36]);")
    self.gen_add_code_line("d2tau_dqd2[st_vel*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + cc*SECOND_ORDER_COORDS + dd] = dqd_val;")
    self.gen_add_code_line("d2tau_dqd2[st_vel*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + dd*SECOND_ORDER_COORDS + cc] = dqd_val;")
    self.gen_add_code_line("d2tau_dvdq[st_vel*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + cc*SECOND_ORDER_COORDS + dd] = -dot_prod<T, 36, 1, 1>(rt6, &D3[st_vel*36]) - dot_prod<T, 6, 1, 1>(rp3, &T2[st_vel*6]) + dot_prod<T, 6, 1, 1>(rp4, &T1[st_vel*6]);")
    self.gen_add_code_line("d2tau_dq2[cc*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + st_vel*SECOND_ORDER_COORDS + dd] = dot_prod<T, 36, 1, 1>(rt6, &D2[st_vel*36]) + dot_prod<T, 36, 1, 1>(rt7, &D1[st_vel*36]) - dot_prod<T, 6, 1, 1>(rp5, &T3[st_vel*6]);")
    self.gen_add_code_line("d2tau_dvdq[cc*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + st_vel*SECOND_ORDER_COORDS + dd] = dot_prod<T, 36, 1, 1>(rt6, &D3[st_vel*36]) - dot_prod<T, 6, 1, 1>(rp5, &T4[st_vel*6]);")
    self.gen_add_code_line("T dm_val = dot_prod<T, 36, 1, 1>(rt8, &D4[st_vel*36]);")
    self.gen_add_code_line("dM_dq[cc*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + st_vel*SECOND_ORDER_COORDS + dd] = dm_val;")
    self.gen_add_code_line("dM_dq[st_vel*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + cc*SECOND_ORDER_COORDS + dd] = dm_val;")
    self.gen_add_end_control_flow()
    self.gen_add_code_line("T dqd_diag = static_cast<T>(0);")
    self.gen_add_code_line("for (int row = 0; row < 6; ++row) {", True)
    self.gen_add_code_line("dqd_diag += S_vel[dd*6 + row] * dot_prod<T, 6, 6, 1>(&IC[jid*36 + row], rp3);")
    self.gen_add_code_line("dqd_diag += S_vel[cc*6 + row] * dot_prod<T, 6, 6, 1>(&crf_S[dd*36 + row], &IC_S[dd*6]);")
    self.gen_add_end_control_flow()
    self.gen_add_code_line("d2tau_dqd2[cc*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + dd*SECOND_ORDER_COORDS + dd] = dqd_diag;")
    self.gen_add_code_line("for (int succ_pos = succ_begin; succ_pos < succ_end; ++succ_pos) {", True)
    self.gen_add_code_line("int succ_vel = successor_v_index[succ_pos];")
    self.gen_add_code_line("T dqd_succ = dot_prod<T, 36, 1, 1>(rt8, &D3[succ_vel*36]);")
    self.gen_add_code_line("d2tau_dqd2[cc*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + succ_vel*SECOND_ORDER_COORDS + dd] = dqd_succ;")
    self.gen_add_code_line("d2tau_dqd2[cc*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + dd*SECOND_ORDER_COORDS + succ_vel] = dqd_succ;")
    self.gen_add_code_line("d2tau_dvdq[cc*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + dd*SECOND_ORDER_COORDS + succ_vel] = dot_prod<T, 36, 1, 1>(rt8, &D2[succ_vel*36]) + dot_prod<T, 36, 1, 1>(rt9, &D1[succ_vel*36]);")
    self.gen_add_code_line("d2tau_dq2[cc*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + dd*SECOND_ORDER_COORDS + succ_vel] = d2tau_dq2[cc*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + succ_vel*SECOND_ORDER_COORDS + dd];")
    self.gen_add_end_control_flow()
    self.gen_add_end_control_flow()
    self.gen_add_code_line("for (int succ_pos = succ_begin; succ_pos < succ_end; ++succ_pos) {", True)
    self.gen_add_code_line("int succ_vel = successor_v_index[succ_pos];")
    self.gen_add_code_line("d2tau_dq2[dd*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + cc*SECOND_ORDER_COORDS + succ_vel] = dot_prod<T, 36, 1, 1>(rt1, &D2[succ_vel*36]) + dot_prod<T, 36, 1, 1>(rt4, &D1[succ_vel*36]);")
    self.gen_add_code_line("T dqd_child = dot_prod<T, 36, 1, 1>(rt2, &D3[succ_vel*36]);")
    self.gen_add_code_line("d2tau_dqd2[dd*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + cc*SECOND_ORDER_COORDS + succ_vel] = dqd_child;")
    self.gen_add_code_line("d2tau_dqd2[dd*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + succ_vel*SECOND_ORDER_COORDS + cc] = dqd_child;")
    self.gen_add_code_line("d2tau_dvdq[dd*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + succ_vel*SECOND_ORDER_COORDS + cc] = dot_prod<T, 36, 1, 1>(rt1, &D3[succ_vel*36]);")
    self.gen_add_code_line("d2tau_dq2[dd*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + succ_vel*SECOND_ORDER_COORDS + cc] = d2tau_dq2[dd*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + cc*SECOND_ORDER_COORDS + succ_vel];")
    self.gen_add_code_line("d2tau_dvdq[dd*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + cc*SECOND_ORDER_COORDS + succ_vel] = dot_prod<T, 36, 1, 1>(rt2, &D2[succ_vel*36]) + dot_prod<T, 36, 1, 1>(rt5, &D1[succ_vel*36]);")
    self.gen_add_code_line("T dm_child = dot_prod<T, 36, 1, 1>(rt8, &D1[succ_vel*36]);")
    self.gen_add_code_line("dM_dq[cc*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + dd*SECOND_ORDER_COORDS + succ_vel] = dm_child;")
    self.gen_add_code_line("dM_dq[dd*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + cc*SECOND_ORDER_COORDS + succ_vel] = dm_child;")
    self.gen_add_end_control_flow()
    self.gen_add_code_line("if (ancestor_body == jid) {", True)
    self.gen_add_code_line("for (int st_pos = st_begin; st_pos < st_end; ++st_pos) {", True)
    self.gen_add_code_line("int st_vel = subtree_v_index[st_pos];")
    self.gen_add_code_line("d2tau_dqd2[st_vel*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + dd*SECOND_ORDER_COORDS + cc] = -dot_prod<T, 36, 1, 1>(rt2, &D1[st_vel*36]);")
    self.gen_add_end_control_flow()
    self.gen_add_end_control_flow()
    self.gen_add_end_control_flow()
    self.gen_add_end_control_flow()
    self.gen_add_end_control_flow()
    self.gen_add_end_control_flow()

    self.gen_add_end_control_flow()
    if self.robot.floating_base and self.robot.using_quaternion:
        self.gen_add_code_line("// Reduced quaternion-vector q columns differentiate rotation with a factor of two at the identity convention used by GRiD/RBDReference.")
        self.gen_add_code_line("for (int q_col = 3; q_col < 6 && q_col < NUM_VEL; ++q_col) {", True)
        self.gen_add_code_line("for (int rc = 0; rc < NUM_VEL*NUM_VEL; ++rc) {", True)
        self.gen_add_code_line("d2tau_dq2[rc*NUM_VEL + q_col] *= static_cast<T>(2);")
        self.gen_add_end_control_flow()
        self.gen_add_end_control_flow()
    # Phase B: emit the Lie-tangent gravity-Hessian addition. Main sweep above ran with
    # a_world[5] = 0, so this call provides the missing gravity contribution to d2tau_dq2
    # (mirrors the Python `idsva_so` + `_floating_gravity_d2tau_dq_lie_direct` pattern).
    self.gen_floating_gravity_d2tau_dq_lie_inline(use_thread_group)
    self.gen_add_sync(use_thread_group)
    self.gen_add_end_function()

def gen_idsva_so_inner(self, use_thread_group = False, use_qdd_input = False):
    """
    Generates the inner device function to compute the second order
    idsva.
    """
    if self.robot.floating_base:
        self.gen_idsva_so_floating_reference_inner(use_thread_group, use_qdd_input)
        return

    NV = self.robot.get_num_vel()
    num_bodies = self.robot.get_num_bodies()
    max_bfs_levels = self.robot.get_max_bfs_level()
    n_bfs_levels = max_bfs_levels + 1 # starts at 0

    # construct the boilerplate and function definition
    func_params = ["s_idsva_so is a pointer to memory for the final result of size 4*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS = " + str(4*NV**3), \
                   "s_q is the vector of joint positions", \
                   "s_qd is the vector of joint velocities", \
                   "s_qdd is the vector of joint accelerations", \
                   "s_temp is a pointer to helper shared memory of size  = " + \
                            str(self.gen_idsva_so_inner_temp_mem_size()), \
                   "gravity is the gravity constant"]
    func_def_start = "void idsva_so_inner(T *s_idsva_so, const T *s_q, const T *s_qd, T *s_qdd, "
    func_def_end = "T *s_temp, const T gravity) {"
    func_def_start, func_params = self.gen_insert_helpers_func_def_params(func_def_start, func_params, -2)
    func_notes = ["Assumes s_XImats is updated already for the current s_q"]
    if use_thread_group:
        func_def_start = func_def_start.replace("(", "(cgrps::thread_group tgrp, ")
        func_params.insert(0,"tgrp is the handle to the thread_group running this function")
    func_def = func_def_start + func_def_end
    # then generate the code
    self.gen_add_func_doc("Computes the second order derivatives of inverse dynamics",func_notes,func_params,None)
    self.gen_add_code_line("template <typename T>")
    self.gen_add_code_line("__device__")
    self.gen_add_code_line(func_def, True)


    # MEMORY LAYOUT (s_temp):
        # Xdown(36*NJ)/
            # vJ(6 * NJ)/f(6*NJ)/ICT_S(6*NJ) & 
            # v(6 * NJ)/psid_Sd(6*NJ) & 
            # Sd(6 * NJ) & 
            # aJ(6 * NJ)/IC_v (6 * NJ)/IC_S (6*NJ)/T1 (6*NJ) & 
            # a(6 * NJ)/IC_psid (6 * NJ) & 
            # psid(6 * NJ)
        # IC (36 * NJ)
        # I_Xup (36 * NJ)/
            # S (6*NJ) & 
            # psidd (6 * NJ) & 
            # a_world (6)
            # T2 (6*NJ)
            # T3 (6*NJ)
            # T4 (6*NJ)
        # BC (36 * NJ)
        # B_IC_S (36*NJ)/D3 (36*NJ)
        # crm_v (36 * NJ)/crm_S (36*NJ)
        # crf_v (36 * NJ)/crf_S (36*NJ)
        # crm_psid (36 * NJ)/crf_S_IC (36*NJ)
        # crf_psid (36 * NJ)/D4 (36*NJ)
        # icrf_f (36 * NJ)/D1 (36*NJ)
        # D2 (36*NJ)
        # Xup(36*NJ)/t - t1/t2/t3/t4/t5/t6/t7/t8/t9 [(len(jids_a) * 36]/p1 & p2 & p3 & p4 & p5 & p6 ([len(jids_a) * 6]*6)


    jids_a, ancestors = self.robot.get_jid_ancestor_ids(include_joint=True)
    var_offset = len(jids_a)
    vars = [
        '// Relevant Tensors in the order they appear',
        'T *I = s_XImats + XIMAT_SIZE*NUM_BODIES;', # Inertia Matrices (6x6 for each joint)
        f'T *Xup = s_temp + 11*XIMAT_SIZE*NUM_BODIES;', # Spatial Transforms from parent to child (6x6 for each joint)
        'T *IC = s_temp + XIMAT_SIZE*NUM_BODIES;', # Centroidal Inertia (6x6 for each joint)
        # Xup, I_Xup done being used
        'T *Xdown = s_temp;\n', # Spatial Transforms from child to parent (6x6 for each joint)
        'T *S = IC + XIMAT_SIZE*NUM_BODIES;', # Transformed Joint Subspace Tensors (6x1 for each joint)
        # Xdown done being used
        'T *vJ = Xdown;', # Non-propogated Joint Spatial velocities (6x1 for each joint),
        'T *v = vJ + 6*NUM_BODIES;', # Joint Spatial velocities (6x1 for each joint),
        'T *Sd = v + 6*NUM_BODIES;', # Time derivative of Joint subspace tensor due to each joint moving (6x1 for each joint),
        'T *aJ = Sd + 6*NUM_BODIES;', # Non-propogated Joint Spatial accelerations (6x1 for each joint),
        'T *a = aJ + 6*NUM_BODIES;', # Joint Spatial accelerations (6x1 for each joint),
        'T *psid = a + 6*NUM_BODIES;', # Time derivative of joint subspace tensor due to each joint's parent moving (6x1 for each joint),
        'T *psidd = S + 6*NUM_BODIES;', # 2nd Time derivative of joint subspace tensor due to each joint's parent moving (6x1 for each joint),
        'T *a_world = psidd + 6*NUM_BODIES;', # Acceleration of the world frame (6x1)',
        'T *BC = S + 30*NUM_BODIES + 6;', # Composite body-Coriolis Bias tensor (6x6 for each joint)',
        'T *f = vJ;', # Joint Spatial forces (6x1 for each joint),
        'T *B_IC_S = BC + 36*NUM_BODIES;', # Body coriolis tensor wrt joint subspace (6x6 for each joint)',

        '\n\n',
        '// Temporary Variables for Computations',
        'T *I_Xup = S;', # Temporary to compute IC for I * Xup (6x6 for each joint)
        'T *crm_v = B_IC_S + 36*NUM_BODIES;', # Motion cross product of v (6x6 for each joint)
        'T *crf_v = crm_v + 36*NUM_BODIES;', # Force cross product of v (6x6 for each joint)',
        'T *IC_v = aJ;', # IC @ v (6x1 for each joint)',
        'T *crm_S = crm_v;', # Motion cross product of S (6x6 for each joint),
        'T *crf_S = crf_v;', # Force cross product of S (6x6 for each joint)',
        'T *IC_S = IC_v;', # IC @ S (6x1 for each joint)',
        'T *crm_psid = crf_v + 36*NUM_BODIES;', # Motion cross product of psid (6x6 for each joint)',
        'T *crf_psid = crm_psid + 36*NUM_BODIES;', # Force cross product of psid (6x6 for each joint)',
        'T *IC_psid = a;', # IC @ psid (6x6 for each joint)',
        'T *icrf_f = crf_psid + 36*NUM_BODIES;', # icrf(f) (6x6 for each joint)',
        'T *psid_Sd = v;', # psid + Sd (6x1 for each joint)',
        'T *ICT_S = f;', # IC^T @ S (6x1 for each joint)',

        '\n\n',
        '// Main Temporary Tensors For Backward Pass',
        'T *T1 = IC_S;', # Temporary for IC @ S (6x1 for each joint)',
        'T *T2 = a_world + 6;', # Temporary for -BC.T @ S (6x1 for each joint)',
        'T *T3 = T2 + 6*NUM_BODIES;', # Temporary matrix (6x1 for each joint)',
        'T *T4 = T3 + 6*NUM_BODIES;', # Temporary matrix (6x1 for each joint)',
        'T *D1 = icrf_f;', # Temporary D1 tensor (6x6 for each joint)',
        'T *D2 = D1 + 36*NUM_BODIES;', # Temporary D2 tensor (6x6 for each joint)',
        'T *D3 = B_IC_S;', # Temporary D3 tensor - same as B(IC, S) (6x6 for each joint)',
        'T *D4 = crf_psid;', # Temporary D4 tensor (6x6 for each joint)',
        f'T *t = D2 + 36*NUM_BODIES;', # Temporary outer product tensor for t1-t9 (6x6 for each joint and its ancestors)',
        'T *p1 = t;', # Temporary cross product vector for p1 (6x1 for each joint and its ancestors)',
        f'T *p2 = p1 + 6*{var_offset};', # Temporary cross product vector for p2 (6x1 for each joint and its ancestors)',
        f'T *p3 = p2 + 6*{var_offset};', # Temporary cross product vector for p3 (6x1 for each joint and its ancestors)',
        f'T *p4 = p3 + 6*{var_offset};', # Temporary cross product vector for p4 (6x1 for each joint and its ancestors)',
        f'T *p5 = p4 + 6*{var_offset};', # Temporary cross product vector for p5 (6x1 for each joint and its ancestors)',
        f'T *p6 = p5 + 6*{var_offset};', # Temporary cross product vector used in computation of d2tau_dqd2[ancestor, joint, joint] (6x1 for each joint and its ancestors)',
        'T *crf_S_IC = crm_psid;', # Cross product of S and IC (6x6 for each joint)',
        

        '\n\n',
        '// Final Tensors for Output',
        'T *d2tau_dq2 = s_idsva_so;', # Second positional derivative of the joint torques (NJxNJXNJ)',
        'T *d2tau_dqd2 = d2tau_dq2+ SECOND_ORDER_COORDS*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS;', # Second velocity derivative of the joint torques (NJxNJXNJ)',
        'T *d2tau_dvdq = d2tau_dqd2 + SECOND_ORDER_COORDS*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS;', # Cross velocity/position derivative of the joint torques (NJxNJXNJ)',
        'T *dM_dq = d2tau_dvdq + SECOND_ORDER_COORDS*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS;', # Positional Derivative of the mass matrix (NJxNJXNJ)',
    ]
    
    self.gen_add_code_lines(vars)

    self.gen_add_code_line("// Initialize output tensor; optimized assembly paths only write structurally nonzero entries.")
    self.gen_add_parallel_loop('i', '4*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS', use_thread_group)
    self.gen_add_code_line("s_idsva_so[i] = static_cast<T>(0);")
    self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)

    parent_ind_cpp, S_ind_cpp = self.gen_topology_helpers_pointers_for_cpp([i for i in range(num_bodies)], NO_GRAD_FLAG = True)
    S_sign_cpp = self.gen_topology_S_sign_for_cpp([i for i in range(num_bodies)])
    parent_ind_cpp_for_jid = parent_ind_cpp


    # Compute Xup transformations
    self.gen_add_code_line("\n")
    self.gen_add_code_line("// Compute Xup - parent to child transformation matrices")
    # If parent is base, Copy X to Xup - X matrices always 6x6
    if self.robot.is_serial_chain():
        self.gen_add_code_line('#pragma unroll')
        self.gen_add_code_line('for (int jid = 0; jid < NUM_BODIES; ++jid) {', 1)
        self.gen_add_code_line('// Compute Xup[joint]')
        self.gen_add_code_line('int X_idx = jid*XIMAT_SIZE;')
        self.gen_add_parallel_loop('i','XIMAT_SIZE',use_thread_group)
        self.gen_add_code_line(f'if ({parent_ind_cpp } == -1) Xup[X_idx + i] = s_XImats[X_idx + i]; // Parent is base')
        self.gen_add_code_line(f'else matmul<T>(i, &Xup[{parent_ind_cpp} * XIMAT_SIZE], &s_XImats[X_idx], &Xup[X_idx], XIMAT_SIZE, 0);')
        self.gen_add_end_control_flow()
        self.gen_add_sync(use_thread_group)
        self.gen_add_end_control_flow()
    else:
        for bfs_level in range(n_bfs_levels):
            inds = self.robot.get_ids_by_bfs_level(bfs_level)
            self.gen_add_code_line(f'// Compute Xup for bfs_level {bfs_level}')
            self.gen_add_parallel_loop('i', str(36*len(inds)), use_thread_group)
            if len(inds) > 1: 
                    select_var_vals = [("int", "jid", [str(jid) for jid in inds])]
                    jid_cpp = "jid"
                    level_parent_ind_cpp = parent_ind_cpp_for_jid
                    self.gen_add_multi_threaded_select("(i)", "<", [str((idx+1)*36) for idx, jid in enumerate(inds)], select_var_vals)
            else:
                jid_cpp = str(inds[0])
                level_parent_ind_cpp = str(self.robot.get_parent_id(inds[0]))
            self.gen_add_code_line(f'int X_idx = {jid_cpp}*XIMAT_SIZE;')
            if bfs_level == 0: self.gen_add_code_line(f'Xup[X_idx + i % XIMAT_SIZE] = s_XImats[X_idx + i % XIMAT_SIZE]; // Parent is base')
            else: self.gen_add_code_line(f'matmul<T>(i % 36, &Xup[{level_parent_ind_cpp} * XIMAT_SIZE], &s_XImats[X_idx], &Xup[X_idx], XIMAT_SIZE, 0);')
            self.gen_add_end_control_flow()
            self.gen_add_sync(use_thread_group)
            

    # Next compute IC - Centroidal Rigid Body Inertia
    self.gen_add_code_line("\n\n")
    self.gen_add_code_line("// Compute IC - Centroidal Rigid Body Inertia")
    # First I @ Xup
    self.gen_add_code_line('// First I @ Xup')
    self.gen_add_parallel_loop('i','XIMAT_SIZE*NUM_BODIES',use_thread_group)
    self.gen_add_code_line('// All involved matrices are 6x6')
    self.gen_add_code_line('matmul<T>(i, Xup, I, I_Xup, 36, false);')
    self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)
    # Next Xup.T @ I
    self.gen_add_code_line('// Next Xup.T @ I')
    self.gen_add_parallel_loop('i','XIMAT_SIZE*NUM_BODIES',use_thread_group)
    self.gen_add_code_line('// All involved matrices are 6x6')
    self.gen_add_code_line('int mat_idx = (i / 36) * 36;')
    self.gen_add_code_line("matmul_trans<T>(i % 36, &Xup[mat_idx], &I_Xup[mat_idx], &IC[mat_idx], 'a');")
    self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)

    # Next compute Xdown transformations
    # Just the transpose of internal 3x3 submatrices
    self.gen_add_code_line("\n\n")
    self.gen_add_code_line("// Compute Xdown - child to parent transformation matrices")
    self.gen_add_parallel_loop('i','XIMAT_SIZE*NUM_BODIES',use_thread_group)
    self.gen_add_code_line('size_t idx = i % XIMAT_SIZE;')
    self.gen_add_code_line('size_t sub_idx = idx % 18;')
    # TODO fix magic numbers
    self.gen_add_code_line('if (idx % 18 == 1 || idx % 18 == 4 || idx % 18 == 8 || idx % 18 == 11) {', True)
    self.gen_add_code_line(f'Xdown[i] = Xup[i+5];')
    self.gen_add_code_line(f'Xdown[i+5] = Xup[i];')
    self.gen_add_end_control_flow()
    self.gen_add_code_line('else if (idx % 18 == 2 || idx % 18 == 5) {', True)
    self.gen_add_code_line(f'Xdown[i] = Xup[i+10];')
    self.gen_add_code_line('Xdown[i+10] = Xup[i];')
    self.gen_add_end_control_flow()
    self.gen_add_code_line('else if (sub_idx != 6 && sub_idx != 9 && sub_idx != 13 && sub_idx != 16 &&')
    self.gen_add_code_line('            sub_idx != 12 && sub_idx != 15)', True)
    self.gen_add_code_line(f'Xdown[i] = Xup[i];')
    self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)

    # Transform S
    self.gen_add_code_line("\n\n")
    self.gen_add_code_line('// Transform S')
    self.gen_add_parallel_loop('i','6*NUM_BODIES',use_thread_group)
    self.gen_add_code_line('int jid = i / 6;')
    self.gen_add_code_line(f'S[i] = ({S_sign_cpp}) * Xdown[jid*XIMAT_SIZE + {S_ind_cpp}*6 + (i % 6)];')
    self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)

    # Compute vJ = S @ qd & aJ = S @ qdd in parallel
    self.gen_add_code_line("\n\n")
    self.gen_add_code_line('// Compute vJ = S @ qd & aJ = S @ qdd')
    self.gen_add_parallel_loop('i','2*6*NUM_BODIES',use_thread_group)
    self.gen_add_code_line('int joint = i / 6;')
    self.gen_add_code_line('if (joint < NUM_BODIES) vJ[i] = S[i] * s_qd[joint];')
    self.gen_add_code_line('else aJ[i - 6*NUM_BODIES] = S[i - 6*NUM_BODIES] * s_qdd[joint - NUM_BODIES];')
    self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)

    # Compute v = v[parent] + vJ
    self.gen_add_code_line("\n\n")
    self.gen_add_code_line('// Compute v = v[parent] + vJ')
    if self.robot.is_serial_chain():
        self.gen_add_code_line('#pragma unroll')
        self.gen_add_code_line('for (int jid = 0; jid < NUM_BODIES; ++jid) {', 1)
        self.gen_add_parallel_loop('i','6',use_thread_group)
        self.gen_add_code_line(f'if ({parent_ind_cpp} == -1) v[jid*6 + i] = vJ[jid*6 + i];')
        self.gen_add_code_line(f'else v[jid*6 + i] = v[{parent_ind_cpp}*6 + i] + vJ[jid*6 + i];')
        self.gen_add_end_control_flow()
        self.gen_add_sync(use_thread_group)
        self.gen_add_end_control_flow()
    else:
        for bfs_level in range(n_bfs_levels):
            inds = self.robot.get_ids_by_bfs_level(bfs_level)
            self.gen_add_code_line(f'// Compute v for bfs_level {bfs_level}')
            self.gen_add_parallel_loop('i', str(6*len(inds)), use_thread_group)
            if len(inds) > 1: 
                    select_var_vals = [("int", "jid", [str(jid) for jid in inds])]
                    jid_cpp = "jid"
                    level_parent_ind_cpp = parent_ind_cpp_for_jid
                    self.gen_add_multi_threaded_select("(i)", "<", [str((idx+1)*6) for idx, jid in enumerate(inds)], select_var_vals)
            else:
                jid_cpp = str(inds[0])
                level_parent_ind_cpp = str(self.robot.get_parent_id(inds[0]))
            self.gen_add_code_line(f'int idx = i % 6;')
            if bfs_level == 0: self.gen_add_code_line(f'v[{jid_cpp}*6 + idx] = vJ[{jid_cpp}*6 + idx]; // Parent is base')
            else: self.gen_add_code_line(f'v[{jid_cpp}*6 + idx] = v[{level_parent_ind_cpp}*6 + idx] + vJ[{jid_cpp}*6 + idx];')
            self.gen_add_end_control_flow()
            self.gen_add_sync(use_thread_group)

    # Finish aJ += crm(v[parent])@vJ
    self.gen_add_code_line("\n\n")
    self.gen_add_code_line('// Finish aJ += crm(v[parent])@vJ')
    self.gen_add_code_line('// For base, v[parent] = 0')
    self.gen_add_parallel_loop('i','6*NUM_BODIES',use_thread_group)
    self.gen_add_code_line('int jid = i / 6;')
    self.gen_add_code_line('int index = i % 6;')
    self.gen_add_code_line(f'if ({parent_ind_cpp_for_jid} != -1) aJ[i] += crm_mul<T>(index, &v[{parent_ind_cpp_for_jid}*6], &vJ[jid*6]);')
    self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)

    # Compute Sd = crm(v) @ S & psid = crm(v[parent]) @ S
    self.gen_add_code_line("\n\n")
    self.gen_add_code_line('// Compute Sd = crm(v) @ S & psid = crm(v[parent]) @ S')
    self.gen_add_code_line('// For base, v[parent] = 0')
    self.gen_add_parallel_loop('i','2*6*NUM_BODIES',use_thread_group)
    self.gen_add_code_line('int jid = (i / 6) % NUM_BODIES;')
    self.gen_add_code_line('int index = i % 6;')
    self.gen_add_code_line('if (i < 6*NUM_BODIES) Sd[i] = crm_mul<T>(index, &v[jid*6], &S[jid*6]);')
    self.gen_add_code_line('else {', True)
    self.gen_add_code_line(f'if ({parent_ind_cpp_for_jid} == -1) psid[jid*6 + index] = 0;')
    self.gen_add_code_line(f'else psid[i - 6 * NUM_BODIES] = crm_mul<T>(index, &v[{parent_ind_cpp_for_jid}*6], &S[jid*6]);')   
    self.gen_add_end_control_flow()
    self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)

    # Compute a = a[parent] + aJ
    self.gen_add_code_line("\n\n")
    self.gen_add_code_line('// Compute a = a[parent] + aJ')
    if self.robot.is_serial_chain():
        self.gen_add_code_line('#pragma unroll')
        self.gen_add_code_line('for (int jid = 0; jid < NUM_BODIES; ++jid) {', 1)
        self.gen_add_parallel_loop('i','6',use_thread_group)
        self.gen_add_code_line(f"if ({parent_ind_cpp} == -1) a[jid*6+ i] = aJ[jid*6 + i] + gravity * (i == 5); // Base joint's parent is the world")
        self.gen_add_code_line(f'else a[jid*6 + i] = a[{parent_ind_cpp}*6 + i] + aJ[jid*6 + i];')
        self.gen_add_end_control_flow()
        self.gen_add_sync(use_thread_group)
        self.gen_add_end_control_flow()
    else:
        for bfs_level in range(n_bfs_levels):
            inds = self.robot.get_ids_by_bfs_level(bfs_level)
            self.gen_add_code_line(f'// Compute a for bfs_level {bfs_level}')
            self.gen_add_parallel_loop('i', str(6*len(inds)), use_thread_group)
            if len(inds) > 1: 
                    select_var_vals = [("int", "jid", [str(jid) for jid in inds])]
                    jid_cpp = "jid"
                    level_parent_ind_cpp = parent_ind_cpp_for_jid
                    self.gen_add_multi_threaded_select("(i)", "<", [str((idx+1)*6) for idx, jid in enumerate(inds)], select_var_vals)
            else:
                jid_cpp = str(inds[0])
                level_parent_ind_cpp = str(self.robot.get_parent_id(inds[0]))
            self.gen_add_code_line(f'int idx = i % 6;')
            if bfs_level == 0: self.gen_add_code_line(f"a[{jid_cpp}*6+ idx] = aJ[{jid_cpp}*6 + idx] + gravity * (idx == 5); // Base joint's parent is the world")
            else: self.gen_add_code_line(f'a[{jid_cpp}*6 + idx] = a[{level_parent_ind_cpp}*6 + idx] + aJ[{jid_cpp}*6 + idx];')
            self.gen_add_end_control_flow()
            self.gen_add_sync(use_thread_group)
        

    # Initialize a_world
    self.gen_add_code_line("\n\n")
    self.gen_add_code_line('// Initialize a_world')
    self.gen_add_parallel_loop('i','6',use_thread_group)
    self.gen_add_code_line('if (i < 5) a_world[i] = 0;')
    self.gen_add_code_line('else a_world[5] = gravity;')
    self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)
    
    # Compute psidd = crm(a[parent])@S + crm(v[parent])@psid & IC_v
    self.gen_add_code_line("\n\n")
    self.gen_add_code_line('// Compute psidd = crm(a[parent])@S + crm(v[:,i])@psid[:,i] & IC @ v (for BC) in parallel')
    self.gen_add_parallel_loop('i','2*6*NUM_BODIES',use_thread_group)
    self.gen_add_code_line('int jid = (i / 6) % NUM_BODIES;')
    self.gen_add_code_line('int index = i % 6;')
    self.gen_add_code_line('if (i < 6*NUM_BODIES) {', True)
    self.gen_add_code_line(f'if ({parent_ind_cpp_for_jid} == -1) psidd[i] = crm_mul<T>(index, a_world, &S[jid*6]);')
    self.gen_add_code_line(f'else psidd[i] = crm_mul<T>(index, &a[{parent_ind_cpp_for_jid}*6], &S[jid*6]) + crm_mul<T>(index, &v[{parent_ind_cpp_for_jid}*6], &psid[jid*6]);')
    self.gen_add_end_control_flow()
    self.gen_add_code_line(f'else IC_v[i - 6*NUM_BODIES] = dot_prod<T, 6, 6, 1>(&IC[index + jid*36], &v[jid*6]);')
    self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)

    # Begin BC Computation
    # First Compute crm(v) & crf(v)
    self.gen_add_code_line("\n\n")
    self.gen_add_code_line('// Need crm(v), crf(v) for BC computation')
    self.gen_add_parallel_loop('i','2*36*NUM_BODIES',use_thread_group)
    self.gen_add_code_line('int jid = (i / 36) % NUM_BODIES;')
    self.gen_add_code_line('int col = (i / 6) % 6;')
    self.gen_add_code_line('int row = i % 6;')
    self.gen_add_code_line('if (i < 36*NUM_BODIES) crm_v[i] = crm<T>(i % 36, &v[jid*6]);')
    self.gen_add_code_line('else crf_v[(jid*36) + row*6 + col] = -crm<T>(i % 36, &v[jid*6]); // crf is negative tranpose of crm')
    self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)


    # Finish BC = crf(v) @ IC + icrf(IC @ v) - IC @ crm(v)
    self.gen_add_code_line("\n\n")
    self.gen_add_code_line('// Finish BC = crf(v) @ IC + icrf(IC @ v) - IC @ crm(v)')
    self.gen_add_parallel_loop('i','36*NUM_BODIES',use_thread_group)
    self.gen_add_code_line('int jid = i / 36;')
    self.gen_add_code_line('int row = i % 6;')
    self.gen_add_code_line('int col_idx = (i / 6) * 6;')
    self.gen_add_code_line('BC[i] = dot_prod<T, 6, 6, 1>(&crf_v[jid*36 + row], &IC[col_idx]) +')
    self.gen_add_code_line('        icrf<T>(i % 36, &IC_v[jid*6]) -')
    self.gen_add_code_line('        dot_prod<T, 6, 6, 1>(&IC[jid*36 + row], &crm_v[col_idx]);')
    self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)

    # Next f = IC @ a + crf(v) @ IC @ v
    self.gen_add_code_line("\n\n")
    self.gen_add_code_line('// Compute f = IC @ a + crf(v) @ IC @ v')
    self.gen_add_parallel_loop('i','6*NUM_BODIES',use_thread_group)
    self.gen_add_code_line('int jid = i / 6;')
    self.gen_add_code_line('int row = i % 6;')
    self.gen_add_code_line('f[i] = dot_prod<T, 6, 6, 1>(&IC[jid*36 + row], &a[jid*6]) +')
    self.gen_add_code_line('        dot_prod<T, 6, 6, 1>(&crf_v[jid*36 + row], &IC_v[jid*6]);')
    self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)

    # Forward Pass Completed
    self.gen_add_code_line("\n\n")
    self.gen_add_code_line('// Forward Pass Completed')
    self.gen_add_code_line('// Now compute the backward pass')


    # Compute IC[parent] += IC[i], BC[parent] += BC[i], f[parent] += f[i]
    self.gen_add_code_line("\n\n")
    self.gen_add_code_line('// Compute IC[parent] += IC[i], BC[parent] += BC[i], f[parent] += f[i]')
    if self.robot.is_serial_chain():
        self.gen_add_code_line('#pragma unroll')
        self.gen_add_code_line('for (int jid = NUM_BODIES-1; jid > 0; --jid) {', 1)
        self.gen_add_parallel_loop('i','36*2 + 6',use_thread_group)
        self.gen_add_code_line(f'if ({parent_ind_cpp} != -1) {{', True)
        self.gen_add_code_line(f'if (i < 36) IC[{parent_ind_cpp}*36 + i] += IC[jid*36 + i];')
        self.gen_add_code_line(f'else if (i < 36*2) BC[{parent_ind_cpp}*36 + i - 36] += BC[jid*36 + i - 36];')
        self.gen_add_code_line(f'else f[{parent_ind_cpp}*6 + i - 36*2] += f[jid*6 + i - 36*2];')
        self.gen_add_end_control_flow()
        self.gen_add_end_control_flow()
        self.gen_add_sync(use_thread_group)
        self.gen_add_end_control_flow()
    else:
        for bfs_level in range(n_bfs_levels-1, 0, -1):
            inds = self.robot.get_ids_by_bfs_level(bfs_level)
            self.gen_add_code_line(f'// Compute propogations for bfs_level {bfs_level}')
            for jid in inds:
                parent_ind = self.robot.get_parent_id(jid)
                if parent_ind == -1:
                    continue
                self.gen_add_code_line(
                    f'// Accumulate joint {jid} into parent {parent_ind}'
                )
                self.gen_add_parallel_loop('i','36*2 + 6', use_thread_group)
                self.gen_add_code_line('int idx = i;')
                self.gen_add_code_line(f'if (idx < 36) IC[{parent_ind}*36 + idx] += IC[{jid}*36 + idx];')
                self.gen_add_code_line(f'else if (idx < 36*2) BC[{parent_ind}*36 + idx - 36] += BC[{jid}*36 + idx - 36];')
                self.gen_add_code_line(f'else f[{parent_ind}*6 + idx - 36*2] += f[{jid}*6 + idx - 36*2];')
                self.gen_add_end_control_flow()
                self.gen_add_sync(use_thread_group)

    # Begin B(IC, S) & B(IC, psid) computation
    # First compute crm(S), crf(S), IC @ S && crm(psid), crf(psid), IC @ psid, icrf(f), psid+Sd
    self.gen_add_code_line("\n\n")
    self.gen_add_code_line('// Need crm(S), crf(S), IC@S, crm(psid), crf(psid), IC@psid for B computations & icrf(f), psid+Sd for T3,T4')
    self.gen_add_parallel_loop('i','5*36*NUM_BODIES + 3*6*NUM_BODIES',use_thread_group)
    self.gen_add_code_line('int jid = (i / 36) % NUM_BODIES;')
    self.gen_add_code_line('int jidMatmul = (i / 6) % NUM_BODIES;')
    self.gen_add_code_line('int col = (i / 6) % 6;')
    self.gen_add_code_line('int row = i % 6;')
    self.gen_add_code_line('if (i < 36*NUM_BODIES) crm_S[i] = crm<T>(i % 36, &S[jid*6]);')
    self.gen_add_code_line('else if (i < 2*36*NUM_BODIES) crf_S[(jid*36) + row*6 + col] = -crm<T>(i % 36, &S[jid*6]); // crf is negative tranpose of crm')
    self.gen_add_code_line('else if (i < 3*36*NUM_BODIES) crm_psid[jid*36 + col*6 + row] = crm<T>(i % 36, &psid[jid*6]);')
    self.gen_add_code_line('else if (i < 4*36*NUM_BODIES) crf_psid[(jid*36) + row*6 + col] = -crm<T>(i % 36, &psid[jid*6]); // crf is negative tranpose of crm')
    self.gen_add_code_line('else if (i < 5*36*NUM_BODIES) icrf_f[i - 4*36*NUM_BODIES] = icrf<T>(i % 36, &f[jid*6]);')
    self.gen_add_code_line('else if (i < 5*36*NUM_BODIES + 6*NUM_BODIES) IC_S[i - 5*36*NUM_BODIES] = dot_prod<T, 6, 6, 1>(&IC[row + jidMatmul*36], &S[jidMatmul*6]);')
    self.gen_add_code_line('else if (i < 5*36*NUM_BODIES + 2*6*NUM_BODIES) psid_Sd[i - 5*36*NUM_BODIES - 6*NUM_BODIES] = psid[i - 5*36*NUM_BODIES - 6*NUM_BODIES] + Sd[i - 5*36*NUM_BODIES - 6*NUM_BODIES];')
    self.gen_add_code_line('else IC_psid[i - 5*36*NUM_BODIES - 2*6*NUM_BODIES] = dot_prod<T, 6, 6, 1>(&IC[row + jidMatmul*36], &psid[jidMatmul*6]);')
    self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)

    # Finish B_IC_S, Start D2
    self.gen_add_code_line('\n\n')
    self.gen_add_code_line('// Finish B_IC_S, Start D2')
    self.gen_add_code_line('// B_IC_S = crf(S) @ IC + icrf(IC @ S) - IC @ crm(S)')
    self.gen_add_code_line('// D2 = crf(psid) @ IC + icrf(IC @ psid) - IC @ crm(psid)')
    self.gen_add_parallel_loop('i','2*36*NUM_BODIES',use_thread_group)
    self.gen_add_code_line('int jid = (i / 36) % NUM_BODIES;')
    self.gen_add_code_line('int row = i % 6;')
    self.gen_add_code_line('int col = (i / 6) % 6;')
    self.gen_add_code_line('if (i < 36*NUM_BODIES) {', True)
    self.gen_add_code_line('B_IC_S[i] = dot_prod<T, 6, 6, 1>(&crf_S[jid*36 + row], &IC[jid*36 + col*6]) + ')
    self.gen_add_code_line('            icrf<T>(i % 36, &IC_S[jid*6]) -') 
    self.gen_add_code_line('            dot_prod<T, 6, 6, 1>(&IC[jid*36 + row], &crm_S[jid*36 + col*6]);')
    self.gen_add_end_control_flow()
    self.gen_add_code_line('else {', True)
    self.gen_add_code_line('D2[i - 36*NUM_BODIES] = dot_prod<T, 6, 6, 1>(&crf_psid[jid*36 + row], &IC[jid*36 + col*6]) + ')
    self.gen_add_code_line('                                icrf<T>(i % 36, &IC_psid[jid*6]) -') 
    self.gen_add_code_line('                                dot_prod<T, 6, 6, 1>(&IC[jid*36 + row], &crm_psid[jid*36 + col*6]);')
    self.gen_add_end_control_flow()
    self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)

    # Compute T2 = -BC.T @ S & T3 = BC @ psid + IC @ psidd + icrf(f) @ S, & T4 = BC @ S + IC @ (psid + Sd), & IC.T @ S for D4
    self.gen_add_code_line('\n\n')
    self.gen_add_code_line('// Compute T2 = -BC.T @ S')
    self.gen_add_code_line('// Compute T3 = BC @ psid + IC @ psidd + icrf(f) @ S')
    self.gen_add_code_line('// Compute T4 = BC @ S + IC @ (psid + Sd)')
    self.gen_add_code_line('// Compute IC.T @ S for D4')
    self.gen_add_parallel_loop('i','4*6*NUM_BODIES',use_thread_group)
    self.gen_add_code_line('int jid = (i / 6) % NUM_BODIES;')
    self.gen_add_code_line('int row = i % 6;')
    self.gen_add_code_line('if (i < 6*NUM_BODIES) T2[i] = -dot_prod<T, 6, 1, 1>(&BC[jid*36 + row*6], &S[jid*6]);')
    self.gen_add_code_line('else if (i < 2*6*NUM_BODIES) {', True)
    self.gen_add_code_line('T3[i - 6*NUM_BODIES] = dot_prod<T, 6, 6, 1>(&BC[jid*36 + row], &psid[jid*6]) +')
    self.gen_add_code_line('                    dot_prod<T, 6, 6, 1>(&IC[jid*36 + row], &psidd[jid*6]) +')
    self.gen_add_code_line('                    dot_prod<T, 6, 6, 1>(&icrf_f[jid*36 + row], &S[jid*6]);')
    self.gen_add_end_control_flow()
    self.gen_add_code_line('else if (i < 3*6*NUM_BODIES) {', True)
    self.gen_add_code_line('T4[i - 2*6*NUM_BODIES] = dot_prod<T, 6, 6, 1>(&BC[jid*36 + row], &S[jid*6]) +')
    self.gen_add_code_line('                    dot_prod<T, 6, 6, 1>(&IC[jid*36 + row], &psid_Sd[jid*6]);')
    self.gen_add_end_control_flow()
    self.gen_add_code_line('else ICT_S[i - 3*6*NUM_BODIES] = dot_prod<T, 6, 6, 1>(&IC[jid*36 + row], &S[jid*6]);')
    self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)

    # Compute D1..D4
    self.gen_add_code_line('\n\n')
    self.gen_add_code_line('// Compute D1, D2, D4, crf_S_IC')
    self.gen_add_parallel_loop('i','4*36*NUM_BODIES',use_thread_group)
    self.gen_add_code_line('int jid = (i / 36) % NUM_BODIES;')
    self.gen_add_code_line('int row = i % 6;')
    self.gen_add_code_line('int col = (i / 6) % 6;')
    self.gen_add_code_line('if (i < 36*NUM_BODIES) {', True)
    self.gen_add_code_line('D1[i] = dot_prod<T, 6, 6, 1>(&crf_S[jid*36 + row], &IC[jid*36 + col*6]) -')
    self.gen_add_code_line('        dot_prod<T, 6, 6, 1>(&IC[jid*36 + row], &crm_S[jid*36 + col*6]);')
    self.gen_add_end_control_flow()
    self.gen_add_code_line('else if (i < 2*36*NUM_BODIES) {', True)
    self.gen_add_code_line('D2[i - 36*NUM_BODIES] += dot_prod<T, 6, 6, 1>(&crf_S[jid*36 + row], &BC[jid*36 + col*6]) -')
    self.gen_add_code_line('                        dot_prod<T, 6, 6, 1>(&BC[jid*36 + row], &crm_S[jid*36 + col*6]);')
    self.gen_add_end_control_flow()
    self.gen_add_code_line('else if (i < 3*36*NUM_BODIES) D4[i - 2*36*NUM_BODIES] = icrf<T>(i % 36, &ICT_S[jid*6]);')
    self.gen_add_code_line('else crf_S_IC[i - 3*36*NUM_BODIES] = dot_prod<T, 6, 6, 1>(&crf_S[jid*36 + row], &IC[jid*36 + col*6]);')
    self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)
    


    # Compute t1
    self.gen_add_code_line('\n\n')
    jids_a, ancestors = self.robot.get_jid_ancestor_ids(include_joint=True)
    self.gen_add_code_line('// Compute t1 = outer(S[j], psid[ancestor])')
    self.gen_add_code_line('// t1[j][k] is stored at t[((j*(j+1)/2) + k)*36]')
    self.gen_add_code_line(f'static const int jids[] = {{ {", ".join(map(str, jids_a))} }}; // Joints with ancestor at equivalent index of ancestors_j') 
    self.gen_add_code_line(f'static const int ancestors_j[] = {{ {", ".join(map(str, ancestors))} }}; // Joint or ancestor of joint at equivalent index of jids_a')
    
    # Create t indexing map
    # Initialize the matrix with -1
    t_index_map = [[-1 for _ in range(NV)] for _ in range(NV)]

    # Fill in the map with t_idx
    for t_idx, (j, a) in enumerate(zip(jids_a, ancestors)):
        t_index_map[j][a] = t_idx

    # Emit CUDA code
    self.gen_add_code_line("const int t_index_map[{}][{}] = {{".format(NV, NV))
    for row in t_index_map:
        self.gen_add_code_line("    { " + ", ".join("{:2}".format(x) for x in row) + " },")
    self.gen_add_code_line("};")
    
    self.gen_add_parallel_loop('i',f'{len(jids_a)}*36',use_thread_group)
    self.gen_add_code_line('int jid = jids[i / 36];')
    self.gen_add_code_line('int ancestor_j = ancestors_j[i / 36];')
    self.gen_add_code_line(f'int t_idx = t_index_map[jid][ancestor_j]*36;')
    self.gen_add_code_line('outerProduct<T>(&S[jid*6], &psid[ancestor_j*6], &t[t_idx], 6, 6, i%36);')
    self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)

    # Perform all computations with t1
    self.gen_add_code_line('\n\n')
    self.gen_add_code_line('// Perform all computations with t1')
    jids, ancestors, st = self.robot.get_jid_ancestor_st_ids(True) # Generate indices for the joint, ancestor, and subtree
    self.gen_add_code_line(f'static const int jids_compute[] = {{ {", ".join(map(str, jids))} }}; // Joints with ancestor at equivalent index of ancestors_j') 
    self.gen_add_code_line(f'static const int ancestors_j_compute[] = {{ {", ".join(map(str, ancestors))} }}; // Joint or ancestor of joint at equivalent index of jids')
    self.gen_add_code_line(f'static const int st[] = {{ {", ".join(map(str, st))} }}; // Subtree of joint at equivalent index of jids')
    self.gen_add_code_lines(['// d2tau_dvdq[child, joint, ancestor] = -np.dot(t1, D3[:, child])', \
                             '// d2tau_dq[joint, ancestor, child] = np.dot(t1, D2[:, child])', \
                             '// d2tau_dq[joint, child, ancestor] = -np.dot(t1, D2[:, child])', \
                             '// d2tau_dvdq[joint, child, ancestor] = np.dot(t1, D3[:, child])'])
    self.gen_add_parallel_loop('i',f'{4*len(jids)}',use_thread_group)
    self.gen_add_code_line(f'int index = i % {len(jids)};')
    self.gen_add_code_line(f'int jid = jids_compute[index];')
    self.gen_add_code_line(f'int ancestor_j = ancestors_j_compute[index];')
    self.gen_add_code_line(f'int st_j = st[index];')
    self.gen_add_code_line(f'int t_idx = t_index_map[jid][ancestor_j]*36;')
    self.gen_add_code_line(f'if (i < {len(jids)}) d2tau_dvdq[st_j*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + ancestor_j * SECOND_ORDER_COORDS + jid] = -dot_prod<T, 36, 1, 1>(&t[t_idx], &D3[st_j*36]);')
    self.gen_add_code_line(f'else if (i < {len(jids)*2} && jid != st_j) d2tau_dq2[jid*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + st_j * SECOND_ORDER_COORDS + ancestor_j] = dot_prod<T, 36, 1, 1>(&t[t_idx], &D2[st_j*36]);')
    self.gen_add_code_line(f'else if (i < {len(jids)*3} && jid != st_j) d2tau_dq2[jid*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + ancestor_j * SECOND_ORDER_COORDS + st_j] = dot_prod<T, 36, 1, 1>(&t[t_idx], &D2[st_j*36]);')
    self.gen_add_code_line(f'else if (jid != st_j) d2tau_dvdq[jid*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + ancestor_j * SECOND_ORDER_COORDS + st_j] = dot_prod<T, 36, 1, 1>(&t[t_idx], &D3[st_j*36]);')
    self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)

    # Compute t2
    self.gen_add_code_line('\n\n')
    jids_a, ancestors = self.robot.get_jid_ancestor_ids(include_joint=True)
    self.gen_add_code_line('// Compute t2 = outer(S[j], S[ancestor])')
    self.gen_add_code_line('// t2[j][k] is stored at t[((j*(j+1)/2) + k)*36]')
    self.gen_add_parallel_loop('i',f'{len(jids_a)}*36',use_thread_group)
    self.gen_add_code_line('int jid = jids[i / 36];')
    self.gen_add_code_line('int ancestor_j = ancestors_j[i / 36];')
    self.gen_add_code_line(f'int t_idx = t_index_map[jid][ancestor_j]*36;')
    self.gen_add_code_line('outerProduct<T>(&S[jid*6], &S[ancestor_j*6], &t[t_idx], 6, 6, i%36);')
    self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)

    # Perform all computations with t2
    self.gen_add_code_line('\n\n')
    self.gen_add_code_line('// Perform all computations with t2')
    self.gen_add_code_lines(['// for ancestor d2tau_dqd[child, ancestor, joint] = -np.dot(t2, D3[child])', \
                             '// for joint d2tau_dqd[child, joint, joint] = -np.dot(t2, D1[child])', \
                             '// for child d2tau_dqd[joint, ancestor, child] = np.dot(t2, D3[child])', \
                             '// for ancestor d2tau_dqd[child, joint, ancestor] = -np.dot(t2, D3[child])', \
                             '// for child d2tau_dqd[joint, child, ancestor] = np.dot(t2, D3[child])', \
                             '// for child d2tau_dvdq[joint, ancestor, child] = np.dot(t2, D2[child])'])
    self.gen_add_parallel_loop('i',f'{5*len(jids)}',use_thread_group)
    self.gen_add_code_line(f'int index = i % {len(jids)};')
    self.gen_add_code_line(f'int jid = jids_compute[index];')
    self.gen_add_code_line(f'int ancestor_j = ancestors_j_compute[index];')
    self.gen_add_code_line(f'int st_j = st[index];')
    self.gen_add_code_line(f'int t_idx = t_index_map[jid][ancestor_j]*36;')
    self.gen_add_code_line(f'if (i < {len(jids)} && ancestor_j < jid) d2tau_dqd2[st_j*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + jid * SECOND_ORDER_COORDS + ancestor_j] = -dot_prod<T, 36, 1, 1>(&t[t_idx], &D3[st_j*36]);')
    self.gen_add_code_line(f'else if (i < {len(jids)} && jid == ancestor_j) d2tau_dqd2[st_j*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + ancestor_j * SECOND_ORDER_COORDS + jid] = -dot_prod<T, 36, 1, 1>(&t[t_idx], &D1[st_j*36]);')
    self.gen_add_code_line(f'else if (i < {2*len(jids)} && jid != st_j) d2tau_dqd2[jid*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + st_j * SECOND_ORDER_COORDS + ancestor_j] = dot_prod<T, 36, 1, 1>(&t[t_idx], &D3[st_j*36]);')
    self.gen_add_code_line(f'else if (i < {3*len(jids)} && ancestor_j < jid) d2tau_dqd2[st_j*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + ancestor_j * SECOND_ORDER_COORDS + jid] = -dot_prod<T, 36, 1, 1>(&t[t_idx], &D3[st_j*36]);')
    self.gen_add_code_line(f'else if (i < {4*len(jids)} && jid != st_j) d2tau_dqd2[jid*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + ancestor_j * SECOND_ORDER_COORDS + st_j] = dot_prod<T, 36, 1, 1>(&t[t_idx], &D3[st_j*36]);')
    self.gen_add_code_line(f'else if (i >= {4*len(jids)} && jid != st_j) d2tau_dvdq[jid*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + st_j * SECOND_ORDER_COORDS + ancestor_j] = dot_prod<T, 36, 1, 1>(&t[t_idx], &D2[st_j*36]);')
    self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)



    # Compute t3
    self.gen_add_code_line('\n\n')
    jids_a, ancestors = self.robot.get_jid_ancestor_ids(include_joint=True)
    self.gen_add_code_line('// Compute t3 = outer(psid[j], psid[ancestor])')
    self.gen_add_code_line('// t3[j][k] is stored at t[((j*(j+1)/2) + k)*36]')
    self.gen_add_parallel_loop('i',f'{len(jids_a)}*36',use_thread_group)
    self.gen_add_code_line('int jid = jids[i / 36];')
    self.gen_add_code_line('int ancestor_j = ancestors_j[i / 36];')
    self.gen_add_code_line(f'int t_idx = t_index_map[jid][ancestor_j]*36;')
    self.gen_add_code_line('outerProduct<T>(&psid[jid*6], &psid[ancestor_j*6], &t[t_idx], 6, 6, i%36);')
    self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)

    # Perform all computations with t3
    self.gen_add_code_line('\n\n')
    self.gen_add_code_line('// Perform all computations with t3')
    self.gen_add_code_lines(['// for joint d2tau_dqd[child, joint, ancestor] = -np.dot(t3, D3[:, st_j])', \
                             '// for ancestor d2tau_dqd[child, ancestor, joint] = -np.dot(t3, D3[:, st_j])'])
    self.gen_add_parallel_loop('i',f'{2*len(jids)}',use_thread_group)
    self.gen_add_code_line(f'int index = i % {len(jids)};')
    self.gen_add_code_line(f'int jid = jids_compute[index];')
    self.gen_add_code_line(f'int ancestor_j = ancestors_j_compute[index];')
    self.gen_add_code_line(f'int st_j = st[index];')
    self.gen_add_code_line(f'int t_idx = t_index_map[jid][ancestor_j]*36;')
    self.gen_add_code_line(f'if (i < {len(jids)}) d2tau_dq2[st_j*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + ancestor_j * SECOND_ORDER_COORDS + jid] = -dot_prod<T, 36, 1, 1>(&t[t_idx], &D3[st_j*36]);')
    self.gen_add_code_line(f'else if (ancestor_j < jid) d2tau_dq2[st_j*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + jid * SECOND_ORDER_COORDS + ancestor_j] = -dot_prod<T, 36, 1, 1>(&t[t_idx], &D3[st_j*36]);')
    self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)


    # Compute t4
    self.gen_add_code_line('\n\n')
    jids_a, ancestors = self.robot.get_jid_ancestor_ids(include_joint=True)
    self.gen_add_code_line('// Compute t4 = outer(S[j], psidd[ancestor])')
    self.gen_add_code_line('// t4[j][k] is stored at t[((j*(j+1)/2) + k)*36]')
    self.gen_add_parallel_loop('i',f'{len(jids_a)}*36',use_thread_group)
    self.gen_add_code_line('int jid = jids[i / 36];')
    self.gen_add_code_line('int ancestor_j = ancestors_j[i / 36];')
    self.gen_add_code_line(f'int t_idx = t_index_map[jid][ancestor_j]*36;')
    self.gen_add_code_line('outerProduct<T>(&S[jid*6], &psidd[ancestor_j*6], &t[t_idx], 6, 6, i%36);')
    self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)

    # Perform all computations with t4
    self.gen_add_code_line('\n\n')
    self.gen_add_code_line('// Perform all computations with t4')
    self.gen_add_code_lines(['// for child d2tau_dq[dd, cc, succ_j] += np.dot(t4, D1[:, succ_j])', \
                             '// for child d2tau_dq[dd, succ_j, cc] += np.dot(t4, D1[:, succ_j])'])
    self.gen_add_parallel_loop('i',f'{2*len(jids)}',use_thread_group)
    self.gen_add_code_line(f'int index = i % {len(jids)};')
    self.gen_add_code_line(f'int jid = jids_compute[index];')
    self.gen_add_code_line(f'int ancestor_j = ancestors_j_compute[index];')
    self.gen_add_code_line(f'int st_j = st[index];')
    self.gen_add_code_line(f'int t_idx = t_index_map[jid][ancestor_j]*36;')
    self.gen_add_code_line(f'if (i < {len(jids)} && jid != st_j) d2tau_dq2[jid*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + st_j * SECOND_ORDER_COORDS + ancestor_j] += dot_prod<T, 36, 1, 1>(&t[t_idx], &D1[st_j*36]);')
    self.gen_add_code_line(f'else if (jid != st_j) d2tau_dq2[jid*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + ancestor_j * SECOND_ORDER_COORDS + st_j] += dot_prod<T, 36, 1, 1>(&t[t_idx], &D1[st_j*36]);')
    self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)


    # Compute t5
    self.gen_add_code_line('\n\n')
    jids_a, ancestors = self.robot.get_jid_ancestor_ids(include_joint=True)
    self.gen_add_code_line('// Compute t5 = outer(S[j], (Sd+psid)[ancestor])')
    self.gen_add_code_line('// t5[j][k] is stored at t[((j*(j+1)/2) + k)*36]')
    self.gen_add_parallel_loop('i',f'{len(jids_a)}*36',use_thread_group)
    self.gen_add_code_line('int jid = jids[i / 36];')
    self.gen_add_code_line('int ancestor_j = ancestors_j[i / 36];')
    self.gen_add_code_line(f'int t_idx = t_index_map[jid][ancestor_j]*36;')
    self.gen_add_code_line('outerProduct<T>(&S[jid*6], &psid_Sd[ancestor_j*6], &t[t_idx], 6, 6, i%36);')
    self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)

    # Perform all computations with t5
    self.gen_add_code_line('\n\n')
    self.gen_add_code_line('// Perform all computations with t5')
    self.gen_add_code_lines(['// for child d2tau_dvdq[dd, cc, succ_j] += np.dot(t5, D1[:, succ_j])'])
    self.gen_add_parallel_loop('i',f'{len(jids)}',use_thread_group)
    self.gen_add_code_line(f'int index = i % {len(jids)};')
    self.gen_add_code_line(f'int jid = jids_compute[index];')
    self.gen_add_code_line(f'int ancestor_j = ancestors_j_compute[index];')
    self.gen_add_code_line(f'int st_j = st[index];')
    self.gen_add_code_line(f'int t_idx = t_index_map[jid][ancestor_j]*36;')
    self.gen_add_code_line(f'if (st_j != jid) d2tau_dvdq[jid*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + st_j * SECOND_ORDER_COORDS + ancestor_j] += dot_prod<T, 36, 1, 1>(&t[t_idx], &D1[st_j*36]);')
    self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)


    # Compute t6
    self.gen_add_code_line('\n\n')
    jids_a, ancestors = self.robot.get_jid_ancestor_ids(include_joint=True)
    self.gen_add_code_line('// Compute t6 = outer(S[ancestor], psid[joint])')
    self.gen_add_code_line('// t6[j][k] is stored at t[((j*(j+1)/2) + k)*36]')
    self.gen_add_parallel_loop('i',f'{len(jids_a)}*36',use_thread_group)
    self.gen_add_code_line('int jid = jids[i / 36];')
    self.gen_add_code_line('int ancestor_j = ancestors_j[i / 36];')
    self.gen_add_code_line(f'int t_idx = t_index_map[jid][ancestor_j]*36;')
    self.gen_add_code_line('outerProduct<T>(&S[ancestor_j*6], &psid[jid*6], &t[t_idx], 6, 6, i%36);')
    self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)

    # Perform all computations with t6
    self.gen_add_code_line('\n\n')
    self.gen_add_code_line('// Perform all computations with t6')
    self.gen_add_code_lines(['// for ancestor d2tau_dvdq[st_j, cc, dd] = -np.dot(t6, D3[:, st_j])', \
                             '// for ancestor d2tau_dq[cc, st_j, dd] = np.dot(t6, D2[:, st_j])', \
                             '// for ancestor d2tau_dvdq[cc, st_j, dd] = np.dot(t6, D3[:, st_j])'])
    self.gen_add_parallel_loop('i',f'{3*len(jids)}',use_thread_group)
    self.gen_add_code_line(f'int index = i % {len(jids)};')
    self.gen_add_code_line(f'int jid = jids_compute[index];')
    self.gen_add_code_line(f'int ancestor_j = ancestors_j_compute[index];')
    self.gen_add_code_line(f'int st_j = st[index];')
    self.gen_add_code_line(f'int t_idx = t_index_map[jid][ancestor_j]*36;')
    self.gen_add_code_line('if (ancestor_j < jid) {', True)
    self.gen_add_code_line(f'if (i < {len(jids)}) d2tau_dvdq[st_j*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + jid * SECOND_ORDER_COORDS + ancestor_j] = -dot_prod<T, 36, 1, 1>(&t[t_idx], &D3[st_j*36]);')
    self.gen_add_code_line(f'else if (i < {2*len(jids)}) d2tau_dq2[ancestor_j*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + jid * SECOND_ORDER_COORDS + st_j] = dot_prod<T, 36, 1, 1>(&t[t_idx], &D2[st_j*36]);')
    self.gen_add_code_line('else d2tau_dvdq[ancestor_j*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + jid * SECOND_ORDER_COORDS + st_j] = dot_prod<T, 36, 1, 1>(&t[t_idx], &D3[st_j*36]);')
    self.gen_add_end_control_flow()
    self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)


    # Compute t7
    self.gen_add_code_line('\n\n')
    jids_a, ancestors = self.robot.get_jid_ancestor_ids(include_joint=True)
    self.gen_add_code_line('// Compute t7 = outer(S[ancestor], psidd[joint])')
    self.gen_add_code_line('// t7[j][k] is stored at t[((j*(j+1)/2) + k)*36]')
    self.gen_add_parallel_loop('i',f'{len(jids_a)}*36',use_thread_group)
    self.gen_add_code_line('int jid = jids[i / 36];')
    self.gen_add_code_line('int ancestor_j = ancestors_j[i / 36];')
    self.gen_add_code_line(f'int t_idx = t_index_map[jid][ancestor_j]*36;')
    self.gen_add_code_line('outerProduct<T>(&S[ancestor_j*6], &psidd[jid*6], &t[t_idx], 6, 6, i%36);')
    self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)

    # Perform all computations with t7
    self.gen_add_code_line('\n\n')
    self.gen_add_code_line('// Perform all computations with t7')
    self.gen_add_code_lines(['// for ancestor d2tau_dq[cc, st_j, dd] += np.dot(t7, D1[:, st_j])'])
    self.gen_add_parallel_loop('i',f'{len(jids)}',use_thread_group)
    self.gen_add_code_line(f'int index = i % {len(jids)};')
    self.gen_add_code_line(f'int jid = jids_compute[index];')
    self.gen_add_code_line(f'int ancestor_j = ancestors_j_compute[index];')
    self.gen_add_code_line(f'int st_j = st[index];')
    self.gen_add_code_line(f'int t_idx = t_index_map[jid][ancestor_j]*36;')
    self.gen_add_code_line(f'if (ancestor_j < jid) d2tau_dq2[ancestor_j*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + jid * SECOND_ORDER_COORDS + st_j] += dot_prod<T, 36, 1, 1>(&t[t_idx], &D1[st_j*36]);')
    self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)


    # Compute t8
    self.gen_add_code_line('\n\n')
    jids_a, ancestors = self.robot.get_jid_ancestor_ids(include_joint=True)
    self.gen_add_code_line('// Compute t8 = outer(S[ancestor], S[joint])')
    self.gen_add_code_line('// t8[j][k] is stored at t[((j*(j+1)/2) + k)*36]')
    self.gen_add_parallel_loop('i',f'{len(jids_a)}*36',use_thread_group)
    self.gen_add_code_line('int jid = jids[i / 36];')
    self.gen_add_code_line('int ancestor_j = ancestors_j[i / 36];')
    self.gen_add_code_line(f'int t_idx = t_index_map[jid][ancestor_j]*36;')
    self.gen_add_code_line('outerProduct<T>(&S[ancestor_j*6], &S[jid*6], &t[t_idx], 6, 6, i%36);')
    self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)

    # Perform all computations with t8
    self.gen_add_code_line('\n\n')
    self.gen_add_code_line('// Perform all computations with t8')
    self.gen_add_code_lines(['// for ancestor dM_dq[cc,st_j,dd] = t8.T @ D4[:, st_j]', \
                             '// for ancestor dM_dq[st_j,cc,dd] = t8.T @ D4[:, st_j]', \
                             '// for child dM_dq[cc, dd, succ_j] = np.dot(t8, D1[:, succ_j])', \
                             '// for child dM_dq[dd, cc, succ_j] = np.dot(t8, D1[:, succ_j])'
                             '// for child & ancestor d2tau_dqd[cc, succ_j, dd] = np.dot(t8, D3[:, succ_j])', \
                             '// for child & ancestor d2tau_dqd[cc, dd, succ_j] = np.dot(t8, D3[:, succ_j])', \
                             '// for child & ancestor d2tau_dvdq[cc, dd, succ_j] = np.dot(t8, D2[:, succ_j])'])
    self.gen_add_parallel_loop('i',f'{7*len(jids)}',use_thread_group)
    self.gen_add_code_line(f'int index = i % {len(jids)};')
    self.gen_add_code_line(f'int jid = jids_compute[index];')
    self.gen_add_code_line(f'int ancestor_j = ancestors_j_compute[index];')
    self.gen_add_code_line(f'int st_j = st[index];')
    self.gen_add_code_line(f'int t_idx = t_index_map[jid][ancestor_j]*36;')
    self.gen_add_code_line('if (ancestor_j < jid) {', True)
    self.gen_add_code_line(f'if (i < {len(jids)}) dM_dq[ancestor_j*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + st_j * SECOND_ORDER_COORDS + jid] = dot_prod<T, 36, 1, 1>(&t[t_idx], &D4[st_j*36]);')
    self.gen_add_code_line(f'else if (i < {2*len(jids)}) dM_dq[st_j*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + ancestor_j * SECOND_ORDER_COORDS + jid] = dot_prod<T, 36, 1, 1>(&t[t_idx], &D4[st_j*36]);')
    self.gen_add_code_line('if (st_j != jid) {', True)
    self.gen_add_code_line(f'if (i < {3*len(jids)}) d2tau_dqd2[ancestor_j*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + jid * SECOND_ORDER_COORDS + st_j] = dot_prod<T, 36, 1, 1>(&t[t_idx], &D3[st_j*36]);')
    self.gen_add_code_line(f'else if (i < {4*len(jids)}) d2tau_dqd2[ancestor_j*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + st_j * SECOND_ORDER_COORDS + jid] = dot_prod<T, 36, 1, 1>(&t[t_idx], &D3[st_j*36]);')
    self.gen_add_code_line(f'else if (i < {5*len(jids)}) d2tau_dvdq[ancestor_j*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + st_j * SECOND_ORDER_COORDS + jid] = dot_prod<T, 36, 1, 1>(&t[t_idx], &D2[st_j*36]);')
    self.gen_add_end_control_flow()
    self.gen_add_end_control_flow()
    self.gen_add_code_line(f'if (jid != st_j && i < {6*len(jids)}) dM_dq[ancestor_j*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + jid * SECOND_ORDER_COORDS + st_j] = dot_prod<T, 36, 1, 1>(&t[t_idx], &D1[st_j*36]);')
    self.gen_add_code_line(f'else if (jid != st_j) dM_dq[jid*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + ancestor_j * SECOND_ORDER_COORDS + st_j] = dot_prod<T, 36, 1, 1>(&t[t_idx], &D1[st_j*36]);')
    self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)


    # Compute t9
    self.gen_add_code_line('\n\n')
    jids_a, ancestors = self.robot.get_jid_ancestor_ids(include_joint=True)
    self.gen_add_code_line('// Compute t9 = outer(S[ancestor], (Sd+psid)[joint])')
    self.gen_add_code_line('// t9[j][k] is stored at t[((j*(j+1)/2) + k)*36]')
    self.gen_add_parallel_loop('i',f'{len(jids_a)}*36',use_thread_group)
    self.gen_add_code_line('int jid = jids[i / 36];')
    self.gen_add_code_line('int ancestor_j = ancestors_j[i / 36];')
    self.gen_add_code_line(f'int t_idx = t_index_map[jid][ancestor_j]*36;')
    self.gen_add_code_line('outerProduct<T>(&S[ancestor_j*6], &psid_Sd[jid*6], &t[t_idx], 6, 6, i%36);')
    self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)

    # Perform all computations with t9
    self.gen_add_code_line('\n\n')
    self.gen_add_code_line('// Perform all computations with t9')
    self.gen_add_code_lines(['// for ancestor & child d2tau_dvdq[cc, dd, succ_j] += np.dot(t9, D1[:, succ_j])', \
                             '// for ancestor & child d2tau_dq[cc, dd, succ_j] = d2tau_dq[cc, succ_j, dd]'])
    self.gen_add_parallel_loop('i',f'{2*len(jids)}',use_thread_group)
    self.gen_add_code_line(f'int index = i % {len(jids)};')
    self.gen_add_code_line(f'int jid = jids_compute[index];')
    self.gen_add_code_line(f'int ancestor_j = ancestors_j_compute[index];')
    self.gen_add_code_line(f'int st_j = st[index];')
    self.gen_add_code_line(f'int t_idx = t_index_map[jid][ancestor_j]*36;')
    self.gen_add_code_line(f'if (i < {len(jids)} && ancestor_j < jid && st_j != jid) d2tau_dvdq[ancestor_j*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + st_j * SECOND_ORDER_COORDS + jid] += dot_prod<T, 36, 1, 1>(&t[t_idx], &D1[st_j*36]);')
    self.gen_add_code_line(f'else if (ancestor_j < jid & st_j != jid) d2tau_dq2[ancestor_j*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + st_j * SECOND_ORDER_COORDS + jid] = d2tau_dq2[ancestor_j*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + jid * SECOND_ORDER_COORDS + st_j];')
    self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)
    
    # Compute p1..p6 in parallel
    jids_a, ancestors = self.robot.get_jid_ancestor_ids(include_joint=True)
    self.gen_add_code_line('\n\n')
    self.gen_add_code_line('// Compute p1..p6 in parallel')
    self.gen_add_code_lines(['// p1 = self.crm(psid_c) @ S_d', \
                             '// p2 = self.crm(psidd[:, k]) @ S_d', \
                             '// p3 = self.crm(S_c) @ S_d', \
                             '// p4 = self.crm(Sd_c + psid_c) @ S_d - 2 * self.crm(psid_d) @ S_c', \
                             '// p5 = self.crm(S_d) @ S_c', \
                             '// p6 = IC_S[joint] @ crm(S[ancestor]) + S[ancestor] @ crf_S_IC[joint]'])
    self.gen_add_parallel_loop('i',f'{6*6*len(jids_a)}',use_thread_group)
    self.gen_add_code_line(f'int index = i % {6*len(jids_a)};')
    self.gen_add_code_line(f'int jid = jids[index / 6];')
    self.gen_add_code_line(f'int ancestor_j = ancestors_j[index / 6];')
    self.gen_add_code_line(f'int p_idx = t_index_map[jid][ancestor_j]*6;')
    self.gen_add_code_line(f'if (i < {len(jids_a)*6}) p1[p_idx + i % 6] = crm_mul<T>(i % 6, &psid[ancestor_j*6], &S[jid*6]);')
    self.gen_add_code_line(f'else if (i < {2*len(jids_a)*6}) p2[p_idx + i % 6] = crm_mul<T>(i % 6, &psidd[ancestor_j*6], &S[jid*6]);')
    self.gen_add_code_line(f'else if (i < {3*len(jids_a)*6}) p3[p_idx + i % 6] = crm_mul<T>(i % 6, &S[ancestor_j*6], &S[jid*6]);')
    self.gen_add_code_line(f'else if (i < {4*len(jids_a)*6}) p4[p_idx + i % 6] = crm_mul<T>(i % 6, &psid_Sd[ancestor_j*6], &S[jid*6]) - 2 * crm_mul<T>(i % 6, &psid[jid*6], &S[ancestor_j*6]);')
    self.gen_add_code_line(f'else if (i < {5*len(jids_a)*6}) p5[p_idx + i % 6] = crm_mul<T>(i % 6, &S[jid*6], &S[ancestor_j*6]);')
    self.gen_add_code_line(f'else p6[p_idx + i % 6] = dot_prod<T, 6, 1, 1>(&IC_S[jid*6], &crm_S[ancestor_j*36 + (i % 6)*6]) + dot_prod<T, 6, 1, 1>(&S[ancestor_j*6], &crf_S_IC[jid*36 + (i % 6)*6]);')
    self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)

    # Finish all computations with p1..p6
    self.gen_add_code_line('\n\n')
    self.gen_add_code_line('// Finish all computations with p1..p5')
    self.gen_add_code_lines(['// for joint d2tau_dq[st_j, dd, cc] += -np.dot(p1, T2[:, st_j]) + np.dot(p2, T1[:, st_j])', \
                             '// for ancestor d2tau_dq[st_j, cc, dd] += -np.dot(p1, T2[:, st_j]) + np.dot(p2, T1[:, st_j])', \
                             '// for ancestor d2tau_dvdq[st_j, cc, dd] += -np.dot(p3, T2[:, st_j]) + np.dot(p4, T1[:, st_j])', \
                             '// for ancestor d2tau_dq[cc, st_j, dd] -= np.dot(p5, T3[:, st_j])', \
                             '// for ancestor && child d2tau_dq[cc, dd, succ_j] -= np.dot(p5, T3[:, st_j])', \
                             '// for ancestor d2tau_dvdq[cc, st_j, dd] -= np.dot(p5, T4[:, st_j])'])
    self.gen_add_parallel_loop('i',f'{6*len(jids)}',use_thread_group)
    self.gen_add_code_line(f'int index = i % {len(jids)};')
    self.gen_add_code_line(f'int jid = jids_compute[index];')
    self.gen_add_code_line(f'int ancestor_j = ancestors_j_compute[index];')
    self.gen_add_code_line(f'int st_j = st[index];')
    self.gen_add_code_line(f'int p_idx = t_index_map[jid][ancestor_j]*6;')
    self.gen_add_code_line(f'if (i < {len(jids)}) d2tau_dq2[st_j*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + ancestor_j * SECOND_ORDER_COORDS + jid] += -dot_prod<T, 6, 1, 1>(&p1[p_idx], &T2[st_j*6]) + dot_prod<T, 6, 1, 1>(&p2[p_idx], &T1[st_j*6]);')
    self.gen_add_code_line('else if (ancestor_j < jid) {', True)
    self.gen_add_code_line(f'if (i < {2*len(jids)}) d2tau_dq2[st_j*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + jid * SECOND_ORDER_COORDS + ancestor_j] += -dot_prod<T, 6, 1, 1>(&p1[p_idx], &T2[st_j*6]) + dot_prod<T, 6, 1, 1>(&p2[p_idx], &T1[st_j*6]);')
    self.gen_add_code_line(f'else if (i < {3*len(jids)}) d2tau_dvdq[st_j*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + jid * SECOND_ORDER_COORDS + ancestor_j] += -dot_prod<T, 6, 1, 1>(&p3[p_idx], &T2[st_j*6]) + dot_prod<T, 6, 1, 1>(&p4[p_idx], &T1[st_j*6]);')
    self.gen_add_code_line(f'else if (i < {4*len(jids)}) d2tau_dq2[ancestor_j*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + jid * SECOND_ORDER_COORDS + st_j] -= dot_prod<T, 6, 1, 1>(&p5[p_idx], &T3[st_j*6]);')
    self.gen_add_code_line(f'else if (i < {5*len(jids)} && st_j != jid) d2tau_dq2[ancestor_j*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + st_j * SECOND_ORDER_COORDS + jid] -= dot_prod<T, 6, 1, 1>(&p5[p_idx], &T3[st_j*6]);')
    self.gen_add_code_line(f'else if (i >= {5*len(jids)}) d2tau_dvdq[ancestor_j*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + jid * SECOND_ORDER_COORDS + st_j] -= dot_prod<T, 6, 1, 1>(&p5[p_idx], &T4[st_j*6]);')
    self.gen_add_end_control_flow()
    self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)

    # Finish computation with p6
    self.gen_add_code_line('\n\n')
    self.gen_add_code_line('// Finish computation with p6')
    self.gen_add_code_line('// d2tau_dqd[ancestor, joint, joint] = p6[joint][ancestor] @ S[joint]')
    self.gen_add_parallel_loop('i',f'{len(jids_a)}',use_thread_group)
    self.gen_add_code_line(f'int jid = jids[i];')
    self.gen_add_code_line(f'int ancestor_j = ancestors_j[i];')
    self.gen_add_code_line(f'int p_idx = t_index_map[jid][ancestor_j]*6;')
    self.gen_add_code_line(f'if (ancestor_j < jid) d2tau_dqd2[ancestor_j*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + jid * SECOND_ORDER_COORDS + jid] = dot_prod<T, 6, 1, 1>(&p6[p_idx], &S[jid*6]);')
    self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)

    if self.idsva_so_needs_reference_order_output_repair():
        self.gen_idsva_so_reference_order_output_repair(use_thread_group)

    self.gen_add_end_function()

        
def gen_idsva_so_public_dvdq_layout_repair(self, use_thread_group = False):
    """
    Emit a final public-output repair for the optimized IDSVA-SO assembly path.

    The optimized path stores the d2tau_dvdq block with the last two axes
    transposed relative to RBDReference/public CUDA output. FDSVA consumes the
    inner tensor directly, so callers that consume public-layout tensors should
    run this repair after inner assembly instead of changing the optimized
    assembly order in the first corrective pass.
    """
    if self.robot.floating_base or self.idsva_so_needs_reference_order_output_repair():
        return

    NV = self.robot.get_num_vel()
    block_offset = 2 * NV**3
    self.gen_add_sync(use_thread_group)
    self.gen_add_code_line("// Repair public d2tau_dvdq layout for optimized IDSVA-SO output")
    self.gen_add_parallel_loop("dvdq_swap_idx", "SECOND_ORDER_COORDS*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS", use_thread_group)
    self.gen_add_code_line("int dvdq_i = dvdq_swap_idx / (SECOND_ORDER_COORDS*SECOND_ORDER_COORDS);")
    self.gen_add_code_line("int dvdq_j = (dvdq_swap_idx / SECOND_ORDER_COORDS) % SECOND_ORDER_COORDS;")
    self.gen_add_code_line("int dvdq_k = dvdq_swap_idx % SECOND_ORDER_COORDS;")
    self.gen_add_code_line("if (dvdq_j < dvdq_k) {", True)
    self.gen_add_code_line(f"T dvdq_tmp = s_idsva_so[{block_offset} + dvdq_i*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + dvdq_j*SECOND_ORDER_COORDS + dvdq_k];")
    self.gen_add_code_line(f"s_idsva_so[{block_offset} + dvdq_i*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + dvdq_j*SECOND_ORDER_COORDS + dvdq_k] = s_idsva_so[{block_offset} + dvdq_i*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + dvdq_k*SECOND_ORDER_COORDS + dvdq_j];")
    self.gen_add_code_line(f"s_idsva_so[{block_offset} + dvdq_i*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS + dvdq_k*SECOND_ORDER_COORDS + dvdq_j] = dvdq_tmp;")
    self.gen_add_end_control_flow()
    self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)


def gen_idsva_so_device_temp_mem_size(self):
    return self.gen_idsva_so_inner_temp_mem_size()
    

def gen_idsva_so_device(self, use_thread_group = False, use_qdd_input = False, single_call_timing=False):
    # TODO --- this is all wrong
    NV = self.robot.get_num_vel()
    # construct the boilerplate and function definition
    func_params = ["s_idsva_so is a pointer to memory for the final result of size 4*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS = " + str(4*NV**3), \
                   "s_q is the vector of joint positions", \
                   "s_qd is the vector of joint velocities", \
                   "d_robotModel is the pointer to the initialized model specific helpers on the GPU (XImats, topology_helpers, etc.)", \
                   "gravity is the gravity constant"]
    func_def_start = "void idsva_so_device(T *s_idsva_so, const T *s_q, const T *s_qd, "
    func_def_end = "const robotModel<T> *d_robotModel, const T gravity) {"
    func_notes = []
    if use_thread_group:
        func_def_start += "cgrps::thread_group tgrp, "
        func_params.insert(0,"tgrp is the handle to the thread_group running this function")
    if use_qdd_input:
        func_def_start += "const T *s_qdd, "
        func_params.insert(-2,"s_qdd is the vector of joint accelerations")
    else:
        func_notes.append("optimized for qdd = 0")
    func_def = func_def_start + func_def_end
    self.gen_add_func_doc("Computes the second order derivates of idsva",func_notes,func_params,None)
    self.gen_add_code_line("template <typename T>")
    self.gen_add_code_line("__device__")
    self.gen_add_code_line(func_def, True)
    # add the shared memory variables
    shared_mem_size = self.gen_idsva_so_inner_temp_mem_size()
    self.gen_XImats_helpers_temp_shared_memory_code(shared_mem_size)
    # then load/update XI and run the algo
    self.gen_load_update_XImats_helpers_function_call(use_thread_group)
    self.gen_idsva_so_inner_function_call(use_thread_group)
    self.gen_idsva_so_public_dvdq_layout_repair(use_thread_group)
    self.gen_add_end_function()

def gen_idsva_so_kernel(self, use_thread_group = False, use_qdd_input = False, single_call_timing = False):
    NUM_POS = self.robot.get_num_pos()
    n = self.robot.get_num_vel()
    NJ = self.robot.get_num_joints()
    use_global_output = getattr(self, "idsva_so_use_global_output", NJ > SHARED_MEMORY_JOINT_THRESHOLD)
    # define function def and params
    func_params = ["d_idsva_so is a pointer to memory for the final result of size 4*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS = " + str(4*n**3), \
                   "d_q_dq_u is the vector of joint positions, velocities, and accelerations", \
                   "stride_q_qd_u is the stide between each q, qd, u", \
                   "d_robotModel is the pointer to the initialized model specific helpers on the GPU (XImats, topology_helpers, etc.)", \
                   "gravity is the gravity constant", \
                   "num_timesteps is the length of the trajectory points we need to compute over (or overloaded as test_iters for timing)"]
    func_notes = []
    # Floating-base SO kernels take a global-memory workspace pointer for the
    # gravity-shim's d2X/d2a/d2f spill (Phase D). Fixed-base SO doesn't need it
    # but we still emit the parameter so the kernel signature is uniform across
    # base modes and host launchers don't fork.
    func_def_start = "void idsva_so_kernel(T *d_idsva_so, unsigned char *d_workspace, const T *d_q_qd_u, const int stride_q_qd_u, "
    func_params.insert(1, "d_workspace is a per-timestep global-memory scratch buffer of " + \
                          str(gen_floating_gravity_d2tau_dq_spill_count(self)) + " floats per timestep " + \
                          "(only used when robot.floating_base is True)")
    func_def_end = "const robotModel<T> *d_robotModel, const T gravity, const int NUM_TIMESTEPS) {"
    if use_qdd_input: # TODO
        func_def_start += "const T *d_qdd, "
        func_params.insert(-2,"d_qdd is the vector of joint accelerations")
    func_def = func_def_start + func_def_end
    if single_call_timing:
        func_def = func_def.replace("kernel(", "kernel_single_timing(")
    # then generate the code
    self.gen_add_func_doc("Computes the second order derivatives of inverse dynamics",func_notes,func_params,None)
    self.gen_add_code_line("template <typename T>")
    self.gen_add_code_line("__global__")
    self.gen_add_code_line(func_def, True)
    # add shared memory variables
    extra_t_buffers = [("s_q_qd_u", n*2+NUM_POS)]
    if not use_global_output:
        extra_t_buffers.append(("s_idsva_so", 4*n**3))
    if use_qdd_input:
        extra_t_buffers.append(("s_qdd", n))
    shared_mem_size = self.gen_idsva_so_inner_temp_mem_size()
    self.gen_XImats_helpers_temp_shared_memory_code(shared_mem_size, extra_t_buffers = extra_t_buffers)
    spill_floats = gen_floating_gravity_d2tau_dq_spill_count(self) if self.robot.floating_base else 0
    # `s_temp_spill` is always declared so the inner-function call site has a uniform
    # shape; for fixed-base SO it stays nullptr (inner doesn't dereference it).
    self.gen_add_code_line("T *s_temp_spill = nullptr;")
    if use_qdd_input:
        self.gen_add_code_line(f"T *s_q = s_q_qd_u; T *s_qd = &s_q_qd_u[{NUM_POS}];")
    else:
        self.gen_add_code_line(f"T *s_q = s_q_qd_u; T *s_qd = &s_q_qd_u[{NUM_POS}]; T *s_qdd = &s_q_qd_u[{NUM_POS + n}];")
    if not single_call_timing:
        # load to shared mem and loop over blocks to compute all requested comps
        self.gen_add_parallel_loop("k","NUM_TIMESTEPS",use_thread_group,block_level = True)
        if use_qdd_input: # TODO
            self.gen_kernel_load_inputs("q_qd","stride_q_qd",str(n + NUM_POS),use_thread_group,"qdd",str(n),str(n))
        else:
            self.gen_kernel_load_inputs("q_qd_u","stride_q_qd_u",str(2*n + NUM_POS),use_thread_group)
        if self.robot.floating_base:
            # Per-timestep spill region for the gravity-Hessian helper. The workspace
            # is shared with the gradient kernels, so SO lives at the SO offset within
            # each timestep slot (see `GRID_SO_WORKSPACE_TEMP_OFFSET_BYTES` in grid.cuh).
            self.gen_add_code_line(
                "s_temp_spill = reinterpret_cast<T *>(&d_workspace[k*GRID_WORKSPACE_BYTES_PER_TIMESTEP<T>() + GRID_SO_WORKSPACE_TEMP_OFFSET_BYTES<T>()]);"
            )
        # compute
        self.gen_add_code_line("// compute")
        self.gen_load_update_XImats_helpers_function_call(use_thread_group)
        if use_global_output:
            self.gen_add_code_line("// Write directly to RAM due to output tensor size")
            self.gen_add_code_line(f"T *s_idsva_so = &d_idsva_so[k*{4*n**3}];")
        self.gen_idsva_so_inner_function_call(use_thread_group)
        self.gen_idsva_so_public_dvdq_layout_repair(use_thread_group)
        if not use_global_output: self.gen_kernel_save_result("idsva_so",str(4*n**3),str(4*n**3),use_thread_group)
        self.gen_add_end_control_flow()
    else:
        #repurpose NUM_TIMESTEPS for number of timing reps
        if use_qdd_input: # TODO
            self.gen_kernel_load_inputs_single_timing("q_qd",str(2*n),use_thread_group,"qdd",str(n))
        else:
            self.gen_kernel_load_inputs_single_timing("q_qd_u",str(NUM_POS + 2*n),use_thread_group)
        if self.robot.floating_base:
            # Single-timing path reuses one slot in the workspace across all reps.
            self.gen_add_code_line("s_temp_spill = reinterpret_cast<T *>(&d_workspace[GRID_SO_WORKSPACE_TEMP_OFFSET_BYTES<T>()]);")
        # then compute in loop for timing
        self.gen_add_code_line("// compute with NUM_TIMESTEPS as NUM_REPS for timing")
        self.gen_add_code_line("for (int rep = 0; rep < NUM_TIMESTEPS; rep++){", True)
        if use_qdd_input:
            self.gen_anti_licm_input_reload("q_qd",str(2*n),use_thread_group,"qdd",str(n))
        else:
            self.gen_anti_licm_input_reload("q_qd_u",str(NUM_POS + 2*n),use_thread_group)
        self.gen_load_update_XImats_helpers_function_call(use_thread_group)
        if use_global_output:
            self.gen_add_code_line("// Write directly to RAM due to output tensor size")
            self.gen_add_code_line("T *s_idsva_so = d_idsva_so;")
        self.gen_idsva_so_inner_function_call(use_thread_group)
        self.gen_idsva_so_public_dvdq_layout_repair(use_thread_group)
        self.gen_add_end_control_flow()
        # save to global
        if not use_global_output: self.gen_kernel_save_result_single_timing("idsva_so",str(4*n**3),use_thread_group)
    self.gen_add_end_function()

def gen_idsva_so_host(self, mode = 0):
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
    func_def_start = "void idsva_so_host(gridData<T, KIND> *hd_data, const robotModel<T> *d_robotModel, const T gravity, const int num_timesteps,"
    func_def_end =   "                      const dim3 block_dimms, const dim3 thread_dimms, cudaStream_t *streams) {"
    if single_call_timing:
        func_def_start = func_def_start.replace("(", "_single_timing(")
        func_def_end = "              " + func_def_end
    if compute_only:
        func_def_start = func_def_start.replace("(", "_compute_only(")
        func_def_end = "             " + func_def_end.replace(", cudaStream_t *streams", "")
    # then generate the code
    self.gen_add_func_doc("Compute IDSVA-SO (Inverse Dynamics - Spatial Vector Algebra - Second Order)",\
                          func_notes,func_params,None)
    self.gen_add_code_line("template <typename T, gridDataKind KIND = GRID_DATA_ALL>")
    self.gen_add_code_line("__host__")
    self.gen_add_code_line(func_def_start)
    self.gen_add_code_line(func_def_end, True)
    self.gen_add_code_line("static_assert(KIND == GRID_DATA_ALL || KIND == GRID_DATA_DYNAMICS, \"idsva_so_host requires all-data or dynamics gridData\");")
    func_call_start = "idsva_so_kernel<T><<<block_dimms,thread_dimms,IDSVA_SO_DYNAMIC_SHARED_MEM_BYTES<T>()>>>(hd_data->d_idsva_so," + \
        "hd_data->d_workspace,hd_data->d_q_qd_u,stride_q_qd,"
    func_call_end = "d_robotModel,gravity,num_timesteps);"
    if single_call_timing:
        func_call_start = func_call_start.replace("kernel<T>","kernel_single_timing<T>")

    self.gen_add_code_line("int stride_q_qd = Q_QD_U_STRIDE;")
    if not compute_only:
        # start code with memory transfer
        self.gen_add_code_lines(["// start code with memory transfer", \
                                "gpuErrchk(cudaMemcpyAsync(hd_data->d_q_qd_u,hd_data->h_q_qd_u,stride_q_qd*" + \
                                ("num_timesteps*" if not single_call_timing else "") + "sizeof(T),cudaMemcpyHostToDevice,streams[0]));", \
                                 "gpuErrchk(cudaDeviceSynchronize());"])
    # TODO then compute but adjust for compressed mem and qdd usage
    self.gen_add_code_line("// then call the kernel")
    # TODO - qdd=0 optimization
    
    func_call_code = [f'{func_call_start}{func_call_end}']
    # wrap function call in timing (if needed)
    if single_call_timing:
        func_call_code.insert(0,"struct timespec start, end; clock_gettime(CLOCK_MONOTONIC,&start);")
        func_call_code.append("clock_gettime(CLOCK_MONOTONIC,&end);")
    self.gen_add_code_line("gpuErrchk(grid_check_dynamic_shared_memory_bytes(\"idsva_so\", IDSVA_SO_DYNAMIC_SHARED_MEM_BYTES<T>()));")
    self.gen_add_code_lines(func_call_code)
    if not compute_only:
        # then transfer memory back
        self.gen_add_code_lines(["// finally transfer the result back", \
                                 "gpuErrchk(cudaMemcpy(hd_data->h_idsva_so,hd_data->d_idsva_so,SECOND_ORDER_TENSOR_SIZE*" + \
                                    ("num_timesteps*" if not single_call_timing else "") + "sizeof(T),cudaMemcpyDeviceToHost));",
                                 "gpuErrchk(cudaDeviceSynchronize());"])
        # (Removed) Historical floating-base runtime FD diagnostic that compared the
        # analytic d2tau_dq path against a `inverse_dynamics_gradient`-based finite
        # difference. Gated by `GRID_FLOATING_SO_DQ_MODE`. Deleted alongside the macros
        # — Phase A+B make the analytic floating-base d2tau_dq correct in one pass.

    # finally report out timing if requested
    if single_call_timing:
        self.gen_add_code_line("printf(\"Single Call ID-SO %fus\\n\",time_delta_us_timespec(start,end)/static_cast<double>(num_timesteps));")
    self.gen_add_end_function()

def gen_idsva_so(self, use_thread_group = False):
    # gen the inner code
    self.gen_idsva_so_inner(use_thread_group)
    # gen the wrapper code for device fn
    # self.gen_idsva_so_device(use_thread_group,False) TODO
    # and the kernels
    self.gen_idsva_so_kernel(use_thread_group,False,True)
    self.gen_idsva_so_kernel(use_thread_group,False,False)
    # and host wrapeprs
    self.gen_idsva_so_host(0)
    self.gen_idsva_so_host(1)
    self.gen_idsva_so_host(2)
