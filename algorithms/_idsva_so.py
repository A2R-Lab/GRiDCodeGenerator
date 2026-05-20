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

def gen_idsva_so_body_frame_inner_temp_mem_size(self):
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
        # icrf_f used to be in shared (one of the 8 36*NB body matrices); it's now a
        # 36-float kernel-local array per body, so subtract it from the count.
        body_mat_count = 7 * 36 * num_bodies
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

def _gravity_shim_full_spill_count(robot):
    """Floats of the gravity-Hessian helper's shared portion (when fully spilled).

    Includes dX, a, da, f, df, and 4*6x6 scratch. Pure-shape function; no `self`.
    """
    NV = robot.get_num_vel()
    NB = robot.get_num_bodies()
    return int(36 * NV + 6 * NB + 6 * NV * NB + 6 * NB + 6 * NV * NB + 4 * 36)


def _gravity_shim_use_full_spill(self):
    """Decide whether to spill the gravity-shim's shared portion to d_workspace.

    Triggered by robot size: when leaving the shared portion in `s_temp` would push
    `idsva_so` total shared bytes over the target, we move dX/a/da/f/df/scratch to
    `d_workspace` (in addition to the d2X/d2a/d2f that always spill). Saves 50-60 KB
    for large floating-base robots. The 4*36 scratch buffers become kernel-local arrays.
    """
    if not self.robot.floating_base:
        return False
    return bool(getattr(self, "idsva_so_body_frame_grav_full_spill", False))


def gen_floating_gravity_d2tau_dq_spill_count(self):
    """Floats of the gravity-Hessian helper that live in the global `d_workspace`
    spill region (per timestep). Always includes the three O(NV²·NB) tensors
    {d2X, d2a, d2f}; when `idsva_so_body_frame_grav_full_spill` is set (large robots), also
    includes the previously-shared {dX, a, da, f, df} arrays.
    """
    NV = self.robot.get_num_vel()
    NB = self.robot.get_num_bodies()
    d2X_count = 36 * NV * NV
    d2a_count = 6 * NV * NV * NB
    d2f_count = 6 * NV * NV * NB
    total = d2X_count + d2a_count + d2f_count
    if _gravity_shim_use_full_spill(self):
        # When fully spilled, dX/a/da/f/df move to d_workspace too. The 4*36
        # scratch buffers become kernel-local arrays (not in workspace).
        NV_ = NV; NB_ = NB
        total += 36 * NV_ + 6 * NB_ + 6 * NV_ * NB_ + 6 * NB_ + 6 * NV_ * NB_
    return int(total)


def gen_floating_gravity_d2tau_dq_shared_count(self):
    """Floats of the gravity-Hessian helper that stay in shared memory.

    Default: dX (sparse-but-stored-dense), a/da, f/df, and the 4*36 scratch
    buffers. When `idsva_so_body_frame_grav_full_spill` is set, returns 0 (everything moves
    to d_workspace except the 4*36 scratch which becomes kernel-local).
    """
    if _gravity_shim_use_full_spill(self):
        return 0
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
    floating-base robots still fit. See `gen_floating_gravity_d2tau_dq_shared_count`
    and `gen_floating_gravity_d2tau_dq_spill_count`.
    """
    return int(gen_floating_gravity_d2tau_dq_shared_count(self)
               + gen_floating_gravity_d2tau_dq_spill_count(self))

def gen_floating_gravity_d2tau_dq_lie_inline(self, use_thread_group=False):
    """Emit (inline) the floating-base gravity-Hessian addition into `d2tau_dq2`.

    Translates the Python helper `_floating_gravity_d2tau_dq_lie_direct` (see
    `RBDReference/RBDReference.py`) into CUDA emission for use inside
    `gen_idsva_so_body_frame_floating_reference_inner`.

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
    main_sweep_count = self.gen_idsva_so_body_frame_inner_temp_mem_size() - gen_floating_gravity_d2tau_dq_shared_count(self)
    full_spill = _gravity_shim_use_full_spill(self)
    layout_comment = (
        "// Full-spill layout (size-triggered): dX/a/da/f/df spill to d_workspace; 4*36 scratch is kernel-local."
        if full_spill else
        "// Shared portion (dX / a / da / f / df / 4x36 scratch) lives in s_temp; the\n"
        "// three O(NV*NV*NB) tensors (d2X / d2a / d2f) spill to s_temp_spill which the\n"
        "// kernel emitter points into d_workspace per-timestep."
    )
    self.gen_add_code_lines([
        "// ===== Gravity-Hessian (Lie-tangent) addition into d2tau_dq2 =====",
        *layout_comment.split("\n"),
        f"static const int grav_lie_body[] = {{ {_idsva_so_int_array(lie_meta['body'])} }};",
        f"static const int grav_lie_s_index[] = {{ {_idsva_so_int_array(lie_meta['s_index'])} }};",
        f"static const int grav_lie_s_sign[] = {{ {_idsva_so_int_array(lie_meta['s_sign'])} }};",
        f"static const int grav_lie_is_root_translation[] = {{ {_idsva_so_int_array(lie_meta['is_root_translation'])} }};",
        f"static const int grav_lie_parent[] = {{ {_idsva_so_int_array(parent_ids)} }};",
        "",
    ])
    if full_spill:
        # Everything that used to be in s_temp moves into s_temp_spill AFTER d2X/d2a/d2f.
        # The 4*36 scratch buffers become kernel-local stack arrays.
        self.gen_add_code_lines([
            "// Global-memory spill carve (large per-timestep tensors + the previously-shared dX/a/da/f/df).",
            "T *grav_d2X     = s_temp_spill;",
            "T *grav_d2a     = grav_d2X     + 36*NUM_VEL*NUM_VEL;",
            "T *grav_d2f     = grav_d2a     + 6*NUM_VEL*NUM_VEL*NUM_BODIES;",
            "T *grav_dX      = grav_d2f     + 6*NUM_VEL*NUM_VEL*NUM_BODIES;",
            "T *grav_a       = grav_dX      + 36*NUM_VEL;",
            "T *grav_da      = grav_a       + 6*NUM_BODIES;",
            "T *grav_f       = grav_da      + 6*NUM_VEL*NUM_BODIES;",
            "T *grav_df      = grav_f       + 6*NUM_BODIES;",
            "T grav_invX_buf[36];",
            "T grav_tmpA_buf[36];",
            "T grav_tmpB_buf[36];",
            "T grav_tmpC_buf[36];",
            "T *grav_invX    = grav_invX_buf;",
            "T *grav_tmpA    = grav_tmpA_buf;",
            "T *grav_tmpB    = grav_tmpB_buf;",
            "T *grav_tmpC    = grav_tmpC_buf;",
        ])
    else:
        self.gen_add_code_lines([
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
        ])
    self.gen_add_code_lines([
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

def gen_idsva_so_body_frame_inner_function_call(self, use_thread_group = False, use_qdd_input = False, updated_var_names = None):
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
    id_so_code_start = "idsva_so_body_frame_inner<T>(" + var_names["s_idsva_so_name"] + ", " + var_names["s_q_name"] + ", " + var_names["s_qd_name"] + ", " + var_names["s_qdd_name"] + ", "
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

def gen_idsva_so_body_frame_reference_order_output_repair(self, use_thread_group = False):
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

def gen_idsva_so_body_frame_floating_reference_inner(self, use_thread_group = False, use_qdd_input = False):
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
                            str(self.gen_idsva_so_body_frame_inner_temp_mem_size()), \
                   "gravity is the gravity constant"]
    func_def_start = "void idsva_so_body_frame_inner(T *s_idsva_so, const T *s_q, const T *s_qd, T *s_qdd, "
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
        "// icrf_f used to live here (size 36*NUM_BODIES) but is now a kernel-local",
        "// per-body 36-float array; saves NUM_BODIES * 36 floats of shared memory.",
        "T *B_IC_S = crf_psid + 36*NUM_VEL;",
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
    # BUG FIX (2026-05-15): the previous emission fused the psid_vel and psidd_vel rows
    # into one inner loop, which caused `crm_mul(row, &v, &psid_vel[vel*6])` to read
    # rows of `psid_vel[vel*6 + 1..5]` that had not yet been written in the current
    # outer iteration. That left a stale (or zero) value in those slots, producing
    # wrong psidd_vel rows 0, 1, 3, 4 (rows 2 and 5 happened to be correct because
    # crm_mul only reads psid_vel components <= row for those indices). The fix is to
    # fully populate psid_vel for the current vel BEFORE consuming it in psidd_vel.
    self.gen_add_code_line("for (int pos = body_v_start[jid]; pos < body_v_start[jid + 1]; ++pos) {", True)
    self.gen_add_code_line("int vel = body_v_index[pos];")
    self.gen_add_code_line("for (int row = 0; row < 6; ++row) psid_vel[vel*6 + row] = crm_mul<T>(row, &v[jid*6], &S_vel[vel*6]);")
    self.gen_add_code_line("for (int row = 0; row < 6; ++row) psidd_vel[vel*6 + row] = crm_mul<T>(row, &a[jid*6], &S_vel[vel*6]) + crm_mul<T>(row, &v[jid*6], &psid_vel[vel*6]);")
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
    self.gen_add_code_line("// icrf(f[jid]) is consumed only within this per-jid loop's T3 build, so")
    self.gen_add_code_line("// keep it as a 36-float kernel-local array instead of NUM_BODIES*36 in shared.")
    self.gen_add_code_line("T icrf_f_local[36];")
    self.gen_add_code_line("for (int idx = 0; idx < 36; ++idx) icrf_f_local[idx] = icrf<T>(idx, &f[jid*6]);")
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
    self.gen_add_code_line("T3[vel*6 + row] = dot_prod<T, 6, 6, 1>(&BC[jid*36 + row], &psid_vel[vel*6]) + dot_prod<T, 6, 6, 1>(&IC[jid*36 + row], &psidd_vel[vel*6]) + dot_prod<T, 6, 6, 1>(&icrf_f_local[row], &S_vel[vel*6]);")
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

def gen_idsva_so_body_frame_inner(self, use_thread_group = False, use_qdd_input = False):
    """
    Generates the inner device function to compute the second order
    idsva.
    """
    if self.robot.floating_base:
        self.gen_idsva_so_body_frame_floating_reference_inner(use_thread_group, use_qdd_input)
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
                            str(self.gen_idsva_so_body_frame_inner_temp_mem_size()), \
                   "gravity is the gravity constant"]
    func_def_start = "void idsva_so_body_frame_inner(T *s_idsva_so, const T *s_q, const T *s_qd, T *s_qdd, "
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
        self.gen_idsva_so_body_frame_reference_order_output_repair(use_thread_group)

    self.gen_add_end_function()

        
def gen_idsva_so_body_frame_public_dvdq_layout_repair(self, use_thread_group = False):
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


def gen_idsva_so_body_frame_device_temp_mem_size(self):
    return self.gen_idsva_so_body_frame_inner_temp_mem_size()
    

def gen_idsva_so_body_frame_device(self, use_thread_group = False, use_qdd_input = False, single_call_timing=False):
    # TODO --- this is all wrong
    NV = self.robot.get_num_vel()
    # construct the boilerplate and function definition
    func_params = ["s_idsva_so is a pointer to memory for the final result of size 4*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS = " + str(4*NV**3), \
                   "s_q is the vector of joint positions", \
                   "s_qd is the vector of joint velocities", \
                   "d_robotModel is the pointer to the initialized model specific helpers on the GPU (XImats, topology_helpers, etc.)", \
                   "gravity is the gravity constant"]
    func_def_start = "void idsva_so_body_frame_device(T *s_idsva_so, const T *s_q, const T *s_qd, "
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
    shared_mem_size = self.gen_idsva_so_body_frame_inner_temp_mem_size()
    self.gen_XImats_helpers_temp_shared_memory_code(shared_mem_size)
    # then load/update XI and run the algo
    self.gen_load_update_XImats_helpers_function_call(use_thread_group)
    self.gen_idsva_so_body_frame_inner_function_call(use_thread_group)
    self.gen_idsva_so_body_frame_public_dvdq_layout_repair(use_thread_group)
    self.gen_add_end_function()

def gen_idsva_so_body_frame_kernel(self, use_thread_group = False, use_qdd_input = False, single_call_timing = False):
    NUM_POS = self.robot.get_num_pos()
    n = self.robot.get_num_vel()
    NJ = self.robot.get_num_joints()
    use_global_output = getattr(self, "idsva_so_body_frame_use_global_output", NJ > SHARED_MEMORY_JOINT_THRESHOLD)
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
    func_def_start = "void idsva_so_body_frame_kernel(T *d_idsva_so, unsigned char *d_workspace, const T *d_q_qd_u, const int stride_q_qd_u, "
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
    self.gen_add_code_line("template <typename T, int RESOURCE_TIER = TIER_PERF>")
    self.gen_add_code_line("__global__")
    self.gen_add_code_line("__launch_bounds__(tier_max_threads<RESOURCE_TIER>())")
    self.gen_add_code_line(func_def, True)
    # add shared memory variables
    extra_t_buffers = [("s_q_qd_u", n*2+NUM_POS)]
    if not use_global_output:
        extra_t_buffers.append(("s_idsva_so", 4*n**3))
    if use_qdd_input:
        extra_t_buffers.append(("s_qdd", n))
    shared_mem_size = self.gen_idsva_so_body_frame_inner_temp_mem_size()
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
        self.gen_idsva_so_body_frame_inner_function_call(use_thread_group)
        self.gen_idsva_so_body_frame_public_dvdq_layout_repair(use_thread_group)
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
        self.gen_idsva_so_body_frame_inner_function_call(use_thread_group)
        self.gen_idsva_so_body_frame_public_dvdq_layout_repair(use_thread_group)
        self.gen_add_end_control_flow()
        # save to global
        if not use_global_output: self.gen_kernel_save_result_single_timing("idsva_so",str(4*n**3),use_thread_group)
    self.gen_add_end_function()

def gen_idsva_so_body_frame_host(self, mode = 0):
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
    func_def_start = "void idsva_so_body_frame_host(gridData<T, KIND> *hd_data, const robotModel<T> *d_robotModel, const T gravity, const int num_timesteps,"
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
    self.gen_add_code_line("static_assert(KIND == GRID_DATA_ALL || KIND == GRID_DATA_DYNAMICS, \"idsva_so_body_frame_host requires all-data or dynamics gridData\");")
    func_call_start = "idsva_so_body_frame_kernel<T><<<block_dimms,thread_dimms,IDSVA_SO_BODY_FRAME_DYNAMIC_SHARED_MEM_BYTES<T>()>>>(hd_data->d_idsva_so," + \
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
                                 "gpuErrchkKernel();"])
    # TODO then compute but adjust for compressed mem and qdd usage
    self.gen_add_code_line("// then call the kernel")
    # TODO - qdd=0 optimization
    
    func_call_code = [f'{func_call_start}{func_call_end}']
    # wrap function call in timing (if needed). The sync between the launch
    # and `clock_gettime(end)` is REQUIRED: kernel launch is async, so
    # without it the timer captures only host-side launch overhead, not
    # actual kernel work. Other algorithms (FD/ABA/etc.) already do this.
    if single_call_timing:
        func_call_code.insert(0,"struct timespec start, end; clock_gettime(CLOCK_MONOTONIC,&start);")
        func_call_code.append("gpuErrchkKernel();")
        func_call_code.append("clock_gettime(CLOCK_MONOTONIC,&end);")
    self.gen_add_code_line("gpuErrchk(grid_check_dynamic_shared_memory_bytes(\"idsva_so\", IDSVA_SO_BODY_FRAME_DYNAMIC_SHARED_MEM_BYTES<T>()));")
    self.gen_add_code_lines(func_call_code)
    if not compute_only:
        # then transfer memory back
        self.gen_add_code_lines(["// finally transfer the result back", \
                                 "gpuErrchk(cudaMemcpy(hd_data->h_idsva_so,hd_data->d_idsva_so,SECOND_ORDER_TENSOR_SIZE*" + \
                                    ("num_timesteps*" if not single_call_timing else "") + "sizeof(T),cudaMemcpyDeviceToHost));",
                                 "gpuErrchkKernel();"])
    else:
        # compute_only path needs an explicit sync after the kernel launch
        # so the caller's batch-timing loop captures actual kernel completion
        # time (otherwise the async launch returns ~immediately and we
        # measure ~launch-overhead per call regardless of N).
        self.gen_add_code_line("gpuErrchkKernel();")

    # finally report out timing if requested. Label format matches the
    # bench's `parse_grid_output` parser, which keys on "single call idsva_so_body_frame"
    # (lowercase): emit "IDSVA_SO_BODY_FRAME" so the parser picks it up. The old
    # "ID-SO" label was silently dropped by the parser → null timings.
    if single_call_timing:
        from ..algo_registry import single_call_printf_line
        self.gen_add_code_line(single_call_printf_line("idsva_so_body_frame"))
    self.gen_add_end_function()

def gen_idsva_so_body_frame(self, use_thread_group = False):
    # gen the inner code
    self.gen_idsva_so_body_frame_inner(use_thread_group)
    # gen the wrapper code for device fn
    # self.gen_idsva_so_body_frame_device(use_thread_group,False) TODO
    # and the kernels
    self.gen_idsva_so_body_frame_kernel(use_thread_group,False,True)
    self.gen_idsva_so_body_frame_kernel(use_thread_group,False,False)
    # and host wrapeprs
    self.gen_idsva_so_body_frame_host(0)
    self.gen_idsva_so_body_frame_host(1)
    self.gen_idsva_so_body_frame_host(2)


# =============================================================================
# world-frame IDSVA-SO: a separate, single-pass CUDA emission that mirrors
# `RBDReference.idsva_so_world_frame` (a faithful port of spatial_v2_extended's
# `ID_SO_derivatives.m`). World-frame propagation; gravity baked into the main
# sweep at the floating-base root; no separate gravity-shim. Co-exists with the
# existing shim-based `gen_idsva_so_body_frame_floating_reference_inner` path.
# =============================================================================

def gen_idsva_so_world_frame_temp_mem_size(self):
    """Shared-memory float count for the world-frame inner.

    Layout:
      - Xup, Xdown, IC, BC: 4 * 36 * NB
      - v, a, f:           3 *  6 * NB
      - S, Sd, psid, psidd: 4 *  6 * NV
      - Per-(i, p) scratch (A0..A7, Bic_phi, Bic_psid): 10 * 36
      - Per-(j, t) scratch (u1..u12): 12 * 6
      - a_grav scratch: 6
      - Per-(i, p) helper vectors (ICi_S, ICi_psid, ICi_psidd, BCi_S, BCi_psid,
        BCiT_S, crf_S_f_i, A5_vec, A7_vec): 9 * 6 — moved from thread-0 stack
        to shared so the idx-over-36 parallel loop in phase 5a can read them.
    """
    NV = self.robot.get_num_vel()
    NB = self.robot.get_num_bodies()
    return int(4 * 36 * NB + 3 * 6 * NB + 4 * 6 * NV + 10 * 36 + 12 * 6 + 6 + 9 * 6)


def gen_idsva_so_world_frame_inner(self, use_thread_group = False, use_qdd_input = False):
    """Emit `idsva_so_world_frame_inner` — a clean world-frame IDSVA-SO.

    Mirrors `RBDReference.idsva_so_world_frame`:
      - World-frame quantities: `S[i] = Xdown[i] @ S_local`, `IC[i] = Xup[i].T @ I @ Xup[i]`.
      - Root acceleration `a[:, 0] = -a_grav` (world frame, gravity baked in).
      - Floating-base root has `Xup[0] = inv(X_local[0])` (Featherstone xlt-inverse pattern).
      - Triple ancestor walk `(i over bodies reverse, p over body i's velocity columns,
        j over ancestors-or-self of i, t over body j's velocity columns, k over ancestors-of-j,
        r over body k's velocity columns)` produces d2tau_dq, d2tau_dqd, d2tau_dvdq, dM_dq.
      - NO gravity-shim. NO `*= 2` quaternion scaling.

    SIMT-parallel: Phases that are inherently parent-dependent (Xup forward pass,
    forward sweep, per-(i,p) and per-(j,t) intermediate builds, IC/BC/f aggregate
    bubble-up) run under a thread-0 guard; phases that are pleasingly parallel
    (Xdown per-body, S_vel per-velocity, output-init, final transpose) use a
    parallel_loop. The dominant triple ancestor walk's innermost (k, rr)
    iteration set is flattened and distributed across threads (each thread owns
    a disjoint vel_k, so output-cell writes are race-free).
    """
    NV = self.robot.get_num_vel()
    NB = self.robot.get_num_bodies()
    metadata = _idsva_so_floating_velocity_metadata(self.robot)
    parent_ids = [self.robot.get_parent_id(body_id) for body_id in range(NB)]

    func_params = [
        "s_idsva_so is a pointer to memory for the final result of size 4*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS = " + str(4*NV**3),
        "s_q is the vector of joint positions",
        "s_qd is the vector of joint velocities",
        "s_qdd is the vector of joint accelerations",
        "s_temp is a pointer to helper shared memory of size = " + str(self.gen_idsva_so_world_frame_temp_mem_size()),
        "gravity is the gravity constant",
    ]
    func_def_start = "void idsva_so_world_frame_inner(T *s_idsva_so, const T *s_q, const T *s_qd, T *s_qdd, "
    func_def_end = "T *s_temp, const T gravity) {"
    func_def_start, func_params = self.gen_insert_helpers_func_def_params(func_def_start, func_params, -2)
    func_notes = [
        "world-frame propagation, gravity baked into main sweep.",
        "Mirrors RBDReference.idsva_so_world_frame (port of spatial_v2_extended ID_SO_derivatives.m).",
        "SIMT-parallel: pleasingly parallel phases distribute across threads; sequential phases (Xup, forward sweep, per-(i,p)/(j,t) intermediates) run under a thread-0 guard.",
    ]
    if use_thread_group:
        func_def_start = func_def_start.replace("(", "(cgrps::thread_group tgrp, ")
        func_params.insert(0, "tgrp is the handle to the thread_group running this function")
    func_def = func_def_start + func_def_end

    self.gen_add_func_doc(
        "Computes IDSVA second-order derivatives via the world-frame single-pass formulation",
        func_notes, func_params, None,
    )
    self.gen_add_code_line("template <typename T>")
    self.gen_add_code_line("__device__")
    self.gen_add_code_line(func_def, True)

    self.gen_add_code_lines([
        "// world-frame IDSVA-SO shared-memory layout.",
        "T *Ipool   = s_XImats + XIMAT_SIZE*NUM_BODIES;",
        "T *Xup     = s_temp;",
        "T *Xdown   = Xup     + 36*NUM_BODIES;",
        "T *IC      = Xdown   + 36*NUM_BODIES;",
        "T *BC      = IC      + 36*NUM_BODIES;",
        "T *v_w     = BC      + 36*NUM_BODIES;",
        "T *a_w     = v_w     +  6*NUM_BODIES;",
        "T *f_w     = a_w     +  6*NUM_BODIES;",
        "T *S_vel   = f_w     +  6*NUM_BODIES;",
        "T *Sd_vel  = S_vel   +  6*NUM_VEL;",
        "T *psid_v  = Sd_vel  +  6*NUM_VEL;",
        "T *psidd_v = psid_v  +  6*NUM_VEL;",
        "T *scratch = psidd_v +  6*NUM_VEL;",
        "// Per-(i,p) scratch blocks (each 6x6 column-major, total 10).",
        "T *S_Bphi  = scratch;            // Bic_phi  (Bic(IC[i], S_p))",
        "T *S_Bpsid = S_Bphi    + 36;     // Bic_psid (Bic(IC[i], psid_p))",
        "T *S_A0    = S_Bpsid   + 36;",
        "T *S_A1    = S_A0      + 36;",
        "T *S_A2    = S_A1      + 36;",
        "T *S_A3    = S_A2      + 36;",
        "T *S_A4    = S_A3      + 36;",
        "T *S_A5    = S_A4      + 36;",
        "T *S_A6    = S_A5      + 36;",
        "T *S_A7    = S_A6      + 36;",
        "// Per-(j,t) scratch: u1..u12 (each 6-vector).",
        "T *S_u1    = S_A7      + 36;",
        "T *S_u2    = S_u1      +  6;",
        "T *S_u3    = S_u2      +  6;",
        "T *S_u4    = S_u3      +  6;",
        "T *S_u5    = S_u4      +  6;",
        "T *S_u6    = S_u5      +  6;",
        "T *S_u7    = S_u6      +  6;",
        "T *S_u8    = S_u7      +  6;",
        "T *S_u9    = S_u8      +  6;",
        "T *S_u10   = S_u9      +  6;",
        "T *S_u11   = S_u10     +  6;",
        "T *S_u12   = S_u11     +  6;",
        "T *S_agrav = S_u12     +  6;",
        "// Per-(i,p) vector intermediates (used by phase 5a parallel idx-over-36 build of A0..A7).",
        "T *S_ICi_S      = S_agrav      +  6;",
        "T *S_ICi_psid   = S_ICi_S      +  6;",
        "T *S_ICi_psidd  = S_ICi_psid   +  6;",
        "T *S_BCi_S      = S_ICi_psidd  +  6;",
        "T *S_BCi_psid   = S_BCi_S      +  6;",
        "T *S_BCiT_S     = S_BCi_psid   +  6;",
        "T *S_crf_S_f_i  = S_BCiT_S     +  6;",
        "T *S_A5_vec     = S_crf_S_f_i  +  6;",
        "T *S_A7_vec     = S_A5_vec     +  6;",
        "T *d2tau_dq2  = s_idsva_so;",
        "T *d2tau_dqd2 = d2tau_dq2  + SECOND_ORDER_COORDS*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS;",
        "T *d2tau_dvdq = d2tau_dqd2 + SECOND_ORDER_COORDS*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS;",
        "T *dM_dq      = d2tau_dvdq + SECOND_ORDER_COORDS*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS;",
        "",
        # constexpr (rather than static const) so the values get baked into use
        # sites at compile time and no per-TU linker symbol is emitted. nvcc's
        # nvlink with -rdc=true otherwise complains "Size doesn't match" when
        # the same template instantiation lives in multiple per-algo TUs.
        f"constexpr int wf_parent[] = {{ {_idsva_so_int_array(parent_ids)} }};",
        f"constexpr int wf_body_v_start[] = {{ {_idsva_so_int_array(metadata['body_v_start'])} }};",
        f"constexpr int wf_body_v_index[] = {{ {_idsva_so_int_array(metadata['body_v_index'])} }};",
        f"constexpr int wf_vel_s_index[]  = {{ {_idsva_so_int_array(metadata['vel_s_index'])} }};",
        f"constexpr int wf_vel_s_sign[]   = {{ {_idsva_so_int_array(metadata['vel_s_sign'])} }};",
        "",
    ])

    # ---- Init: zero output tensor in parallel + init S_agrav with thread 0.
    self.gen_add_parallel_loop("out_idx", "SECOND_ORDER_TENSOR_SIZE", use_thread_group)
    self.gen_add_code_line("s_idsva_so[out_idx] = static_cast<T>(0);")
    self.gen_add_end_control_flow()
    self.gen_add_serial_ops(use_thread_group)
    self.gen_add_code_line("// MATLAB convention: a_grav vector with a_grav[5] = GRAVITY (signed, e.g. -9.81).")
    self.gen_add_code_line("// The CUDA `gravity` parameter is the positive magnitude (+9.81) by GRiD convention,")
    self.gen_add_code_line("// so use -gravity here to match RBDReference.idsva_so_world_frame's `a_grav[5] = GRAVITY`.")
    self.gen_add_code_line("S_agrav[0] = static_cast<T>(0); S_agrav[1] = static_cast<T>(0); S_agrav[2] = static_cast<T>(0);")
    self.gen_add_code_line("S_agrav[3] = static_cast<T>(0); S_agrav[4] = static_cast<T>(0); S_agrav[5] = -gravity;")
    self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)

    # ---- Step 1: Build cumulative Xup. For floating-base root, Xup[0] = inv(X_local[0]).
    # BFS parent-dependent — keep sequential under thread-0 guard.
    self.gen_add_code_line("// Build cumulative Xup. Floating-base root: Xup[0] = inv(X_local[0]).")
    floating_base = self.robot.floating_base
    self.gen_add_serial_ops(use_thread_group)
    self.gen_add_code_line("for (int jid = 0; jid < NUM_BODIES; ++jid) {", True)
    self.gen_add_code_line("int parent = wf_parent[jid];")
    self.gen_add_code_line("if (parent < 0) {", True)
    if floating_base:
        # Spatial Plücker `X = [E 0; B E]` with B = -E*r̂. Inverse:
        # `X^{-1} = [E^T 0; -E^T*B*E^T E^T]`. Both blocks of E transpose; the bottom-left
        # block computes -E^T * B * E^T.
        self.gen_add_code_lines([
            "// Floating-base root: Xup[0] = inv(X_local[0]).",
            "// Plücker form X = [E 0; B E] (column-major) with E orthogonal.",
            "// X^{-1} = [E^T 0; -E^T*B*E^T  E^T].",
            "// Step A: write the four blocks of inv into Xup[jid].",
            "// Top-right block (cols 3..5, rows 0..2) of inv is 0.",
            "for (int idx = 0; idx < 36; ++idx) Xup[jid*36 + idx] = static_cast<T>(0);",
            "// Top-left = E^T:  inv[a, b] = X[b, a] for a,b < 3.",
            "for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) Xup[jid*36 + a + 6*b] = s_XImats[jid*36 + b + 6*a];",
            "// Bottom-right = E^T: inv[a+3, b+3] = X[b+3, a+3] for a,b < 3.",
            "for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) Xup[jid*36 + (a + 3) + 6*(b + 3)] = s_XImats[jid*36 + (b + 3) + 6*(a + 3)];",
            "// Bottom-left = -E^T * B * E^T where B = X[3..6, 0..3].",
            "// Compute tmp1 = E^T * B (3x3 @ 3x3).",
            "T wf_tmp_invX_1[9];",
            "for (int a = 0; a < 3; ++a) {",
            "    for (int b = 0; b < 3; ++b) {",
            "        T acc = static_cast<T>(0);",
            "        for (int kk = 0; kk < 3; ++kk) acc += s_XImats[jid*36 + kk + 6*a] * s_XImats[jid*36 + (kk + 3) + 6*b];",
            "        wf_tmp_invX_1[a + 3*b] = acc;",
            "    }",
            "}",
            "// inv[a+3, b] = -(tmp1 @ E^T)[a, b] = -sum_kk tmp1[a, kk] * E^T[kk, b] = -sum_kk tmp1[a, kk] * X[b, kk].",
            "for (int a = 0; a < 3; ++a) {",
            "    for (int b = 0; b < 3; ++b) {",
            "        T acc = static_cast<T>(0);",
            "        for (int kk = 0; kk < 3; ++kk) acc += wf_tmp_invX_1[a + 3*kk] * s_XImats[jid*36 + b + 6*kk];",
            "        Xup[jid*36 + (a + 3) + 6*b] = -acc;",
            "    }",
            "}",
        ])
    else:
        self.gen_add_code_lines([
            "// Fixed-base root: Xup[0] = X_local[0] (no inversion).",
            "for (int idx = 0; idx < 36; ++idx) Xup[jid*36 + idx] = s_XImats[jid*36 + idx];",
        ])
    self.gen_add_end_control_flow()
    self.gen_add_code_line("else {", True)
    self.gen_add_code_line("// Xup[jid] = X_local[jid] @ Xup[parent] (column-major matmul).")
    self.gen_add_code_line("for (int idx = 0; idx < 36; ++idx) {", True)
    self.gen_add_code_line("int row = idx % 6; int col = idx / 6;")
    self.gen_add_code_line("T acc = static_cast<T>(0);")
    self.gen_add_code_line("for (int kk = 0; kk < 6; ++kk) acc += s_XImats[jid*36 + row + 6*kk] * Xup[parent*36 + kk + 6*col];")
    self.gen_add_code_line("Xup[jid*36 + idx] = acc;")
    self.gen_add_end_control_flow()
    self.gen_add_end_control_flow()
    self.gen_add_end_control_flow()
    self.gen_add_end_control_flow()  # close thread-0 guard for Xup
    self.gen_add_sync(use_thread_group)

    # ---- Step 2: Xdown[i] = inv(Xup[i]) — parallel over jid.
    self.gen_add_code_line("// Build Xdown[i] = inv(Xup[i]) using Plücker block inverse:")
    self.gen_add_code_line("// Xup = [E 0; B E]  =>  Xdown = [E^T 0; -E^T*B*E^T  E^T].")
    self.gen_add_parallel_loop("jid", "NUM_BODIES", use_thread_group)
    self.gen_add_code_line("for (int idx = 0; idx < 36; ++idx) Xdown[jid*36 + idx] = static_cast<T>(0);")
    self.gen_add_code_line("// Top-left = E^T and Bottom-right = E^T.")
    self.gen_add_code_line("for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) {", True)
    self.gen_add_code_line("Xdown[jid*36 + a + 6*b] = Xup[jid*36 + b + 6*a];")
    self.gen_add_code_line("Xdown[jid*36 + (a + 3) + 6*(b + 3)] = Xup[jid*36 + (b + 3) + 6*(a + 3)];")
    self.gen_add_end_control_flow()
    self.gen_add_code_line("// Bottom-left = -E^T * B * E^T where B = Xup[3..6, 0..3] (column-major).")
    self.gen_add_code_line("T wf_tmpEt_B[9];")
    self.gen_add_code_line("for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) {", True)
    self.gen_add_code_line("T acc = static_cast<T>(0);")
    self.gen_add_code_line("for (int kk = 0; kk < 3; ++kk) acc += Xup[jid*36 + kk + 6*a] * Xup[jid*36 + (kk + 3) + 6*b];")
    self.gen_add_code_line("wf_tmpEt_B[a + 3*b] = acc;")
    self.gen_add_end_control_flow()
    self.gen_add_code_line("for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) {", True)
    self.gen_add_code_line("T acc = static_cast<T>(0);")
    self.gen_add_code_line("for (int kk = 0; kk < 3; ++kk) acc += wf_tmpEt_B[a + 3*kk] * Xup[jid*36 + b + 6*kk];")
    self.gen_add_code_line("Xdown[jid*36 + (a + 3) + 6*b] = -acc;")
    self.gen_add_end_control_flow()
    self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)

    # ---- Step 3: S_vel = Xdown @ S_local — parallel over vel.
    self.gen_add_code_line("// S_vel[vel] = Xdown[body(vel)] @ S_local[vel] = sign * Xdown[body][:, s_index].")
    self.gen_add_parallel_loop("vel", "NUM_VEL", use_thread_group)
    # vel_to_body deduce from wf_body_v_index. Easier: pre-compute a vel_to_body table.
    self.gen_add_code_line("// Find body containing this vel.")
    self.gen_add_code_line("int jid = -1;")
    self.gen_add_code_line("for (int b = 0; b < NUM_BODIES && jid < 0; ++b) {", True)
    self.gen_add_code_line("for (int pos = wf_body_v_start[b]; pos < wf_body_v_start[b + 1]; ++pos) {", True)
    self.gen_add_code_line("if (wf_body_v_index[pos] == vel) { jid = b; break; }")
    self.gen_add_end_control_flow()
    self.gen_add_end_control_flow()
    self.gen_add_code_line("int s_col = wf_vel_s_index[vel];")
    self.gen_add_code_line("T s_sign = static_cast<T>(wf_vel_s_sign[vel]);")
    self.gen_add_code_line("for (int row = 0; row < 6; ++row) S_vel[vel*6 + row] = s_sign * Xdown[jid*36 + s_col*6 + row];")
    self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)

    # ---- Step 4: Forward sweep — v, a, f, IC, BC, psid, psidd, Sd.
    # Parent-dependent (sequential by jid); under thread-0 guard, sync after.
    self.gen_add_code_line("// Forward sweep: build v, a, f, IC, BC, psid, psidd, Sd.")
    self.gen_add_serial_ops(use_thread_group)
    self.gen_add_code_line("for (int jid = 0; jid < NUM_BODIES; ++jid) {", True)
    self.gen_add_code_line("int parent = wf_parent[jid];")
    # Initialize v[jid], a[jid]
    self.gen_add_code_line("if (parent < 0) {", True)
    self.gen_add_code_line("for (int row = 0; row < 6; ++row) { v_w[jid*6 + row] = static_cast<T>(0); a_w[jid*6 + row] = -S_agrav[row]; }")
    self.gen_add_end_control_flow()
    self.gen_add_code_line("else {", True)
    self.gen_add_code_line("for (int row = 0; row < 6; ++row) { v_w[jid*6 + row] = v_w[parent*6 + row]; a_w[jid*6 + row] = a_w[parent*6 + row]; }")
    self.gen_add_end_control_flow()

    # vJ, aJ, psid, psidd (referring to v[jid], a[jid] which haven't been updated yet).
    self.gen_add_code_line("// vJ = sum_p S_vel[p] * qd[p]; aJ = sum_p S_vel[p] * qdd[p].")
    self.gen_add_code_line("T vJ[6] = {0,0,0,0,0,0}; T aJ[6] = {0,0,0,0,0,0};")
    self.gen_add_code_line("for (int pos = wf_body_v_start[jid]; pos < wf_body_v_start[jid + 1]; ++pos) {", True)
    self.gen_add_code_line("int vel = wf_body_v_index[pos];")
    self.gen_add_code_line("for (int row = 0; row < 6; ++row) { vJ[row] += S_vel[vel*6 + row] * s_qd[vel]; aJ[row] += S_vel[vel*6 + row] * s_qdd[vel]; }")
    self.gen_add_end_control_flow()
    self.gen_add_code_line("// aJ += crm(v[jid]) @ vJ.")
    self.gen_add_code_line("for (int row = 0; row < 6; ++row) aJ[row] += crm_mul<T>(row, &v_w[jid*6], vJ);")

    # psid[vel] = crm(v[jid]) @ S_vel[vel], psidd[vel] = crm(a[jid]) @ S + crm(v) @ psid
    self.gen_add_code_line("// psid[vel] = crm(v[jid]) @ S; psidd[vel] = crm(a[jid]) @ S + crm(v[jid]) @ psid.")
    self.gen_add_code_line("for (int pos = wf_body_v_start[jid]; pos < wf_body_v_start[jid + 1]; ++pos) {", True)
    self.gen_add_code_line("int vel = wf_body_v_index[pos];")
    self.gen_add_code_line("for (int row = 0; row < 6; ++row) psid_v[vel*6 + row] = crm_mul<T>(row, &v_w[jid*6], &S_vel[vel*6]);")
    self.gen_add_code_line("for (int row = 0; row < 6; ++row) psidd_v[vel*6 + row] = crm_mul<T>(row, &a_w[jid*6], &S_vel[vel*6]) + crm_mul<T>(row, &v_w[jid*6], &psid_v[vel*6]);")
    self.gen_add_end_control_flow()

    # Update v[jid] += vJ, a[jid] += aJ
    self.gen_add_code_line("for (int row = 0; row < 6; ++row) { v_w[jid*6 + row] += vJ[row]; a_w[jid*6 + row] += aJ[row]; }")

    # Sd[vel] = crm(v[jid]_new) @ S
    self.gen_add_code_line("for (int pos = wf_body_v_start[jid]; pos < wf_body_v_start[jid + 1]; ++pos) {", True)
    self.gen_add_code_line("int vel = wf_body_v_index[pos];")
    self.gen_add_code_line("for (int row = 0; row < 6; ++row) Sd_vel[vel*6 + row] = crm_mul<T>(row, &v_w[jid*6], &S_vel[vel*6]);")
    self.gen_add_end_control_flow()

    # IC[jid] = Xup[jid]^T @ I_body @ Xup[jid].
    self.gen_add_code_line("// IC[jid] = Xup[jid]^T @ I_body @ Xup[jid].")
    self.gen_add_code_line("T I_Xup[36];")
    self.gen_add_code_line("for (int idx = 0; idx < 36; ++idx) {", True)
    self.gen_add_code_line("int row = idx % 6; int col = idx / 6;")
    self.gen_add_code_line("T acc = static_cast<T>(0);")
    self.gen_add_code_line("for (int kk = 0; kk < 6; ++kk) acc += Ipool[jid*36 + row + 6*kk] * Xup[jid*36 + kk + 6*col];")
    self.gen_add_code_line("I_Xup[idx] = acc;")
    self.gen_add_end_control_flow()
    self.gen_add_code_line("for (int idx = 0; idx < 36; ++idx) {", True)
    self.gen_add_code_line("int row = idx % 6; int col = idx / 6;")
    self.gen_add_code_line("T acc = static_cast<T>(0);")
    self.gen_add_code_line("for (int kk = 0; kk < 6; ++kk) acc += Xup[jid*36 + kk + 6*row] * I_Xup[kk + 6*col];")
    self.gen_add_code_line("IC[jid*36 + idx] = acc;")
    self.gen_add_end_control_flow()

    # BC[jid] = crf(v) @ IC + icrf(IC @ v) - IC @ crm(v).
    self.gen_add_code_line("// BC[jid] = crf(v) @ IC + icrf(IC @ v) - IC @ crm(v).  (crf(v) = -crm(v)^T.)")
    self.gen_add_code_line("T IC_v[6];")
    self.gen_add_code_line("for (int row = 0; row < 6; ++row) IC_v[row] = dot_prod<T, 6, 6, 1>(&IC[jid*36 + row], &v_w[jid*6]);")
    self.gen_add_code_line("for (int idx = 0; idx < 36; ++idx) {", True)
    self.gen_add_code_line("int row = idx % 6; int col = idx / 6;")
    # crf(v) @ IC: row r of crf(v) is row r of -crm(v)^T = -(crm(v)[:, r]). So (crf(v) @ IC)[r, c]
    # = sum_k crf(v)[r,k] * IC[k,c] = -sum_k crm(v)[k,r] * IC[k,c]. In column-major,
    # crm(v)[k,r] = crm<T>(k + 6*r, v).
    # Easier: form crf_v[r,c] = -crm<T>(c + 6*r, v) and dot it with IC[jid][:, c].
    # The existing shim does this — copy pattern.
    self.gen_add_code_line("T crf_v_row[6];  T crm_v_col[6];")
    self.gen_add_code_line("for (int kk = 0; kk < 6; ++kk) crf_v_row[kk] = -crm<T>(kk + 6*row, &v_w[jid*6]);")
    self.gen_add_code_line("for (int kk = 0; kk < 6; ++kk) crm_v_col[kk] = crm<T>(kk + 6*col, &v_w[jid*6]);")
    self.gen_add_code_line("T t_crfv_IC = dot_prod<T, 6, 1, 1>(crf_v_row, &IC[jid*36 + 6*col]);")
    self.gen_add_code_line("T t_IC_crmv = dot_prod<T, 6, 6, 1>(&IC[jid*36 + row], crm_v_col);")
    self.gen_add_code_line("BC[jid*36 + idx] = t_crfv_IC + icrf<T>(idx, IC_v) - t_IC_crmv;")
    self.gen_add_end_control_flow()

    # f[jid] = IC @ a + crf(v) @ IC @ v.
    self.gen_add_code_line("// f[jid] = IC[jid] @ a[jid] + crf(v[jid]) @ (IC[jid] @ v[jid]).")
    self.gen_add_code_line("for (int row = 0; row < 6; ++row) {", True)
    self.gen_add_code_line("T crf_v_row2[6];")
    self.gen_add_code_line("for (int kk = 0; kk < 6; ++kk) crf_v_row2[kk] = -crm<T>(kk + 6*row, &v_w[jid*6]);")
    self.gen_add_code_line("f_w[jid*6 + row] = dot_prod<T, 6, 6, 1>(&IC[jid*36 + row], &a_w[jid*6]) + dot_prod<T, 6, 1, 1>(crf_v_row2, IC_v);")
    self.gen_add_end_control_flow()
    self.gen_add_end_control_flow()  # end forward jid loop
    self.gen_add_end_control_flow()  # close thread-0 guard for Phase 4
    self.gen_add_sync(use_thread_group)

    # ---- Step 5: Triple ancestor walk (reverse over bodies).
    # All threads run the outer (i, pp, j, tt) sequencing; inside, thread-0 builds
    # the (i, p) and (j, t) intermediates in shared mem, then a parallel loop
    # distributes the inner (k, rr) work across threads (each thread handles a
    # disjoint vel_k, so output writes are race-free).
    self.gen_add_code_line("// Triple ancestor walk: i over bodies (reverse), p over body i columns,")
    self.gen_add_code_line("// j over ancestors-or-self of i, t over body j columns, k over ancestors-of-j, r over body k columns.")
    self.gen_add_code_line("for (int i = NUM_BODIES - 1; i >= 0; --i) {", True)
    self.gen_add_code_line("for (int pp = wf_body_v_start[i]; pp < wf_body_v_start[i + 1]; ++pp) {", True)
    self.gen_add_code_line("int vel_i = wf_body_v_index[pp];")
    # Per-(i, p) intermediates: build A0..A7, Bphi, Bpsid in shared mem.
    # Pattern: thread-0 builds the small "global per-(i,p)" helper vectors
    # (ICi_S, ICi_psid, ..., A5_vec, A7_vec) into shared, then a parallel
    # loop over idx ∈ [0, 36) builds each A-matrix element using the shared
    # helpers. This unlocks ~10× wall-time speedup on the per-(i,p) work.
    self.gen_add_code_line("// === Per-(i, p) intermediates: thread-0 helpers, then parallel idx-over-36 build of A0..A7 ===")
    self.gen_add_code_line("T *S_p     = &S_vel[vel_i*6];")
    self.gen_add_code_line("T *Sd_p    = &Sd_vel[vel_i*6];")
    self.gen_add_code_line("T *psid_p  = &psid_v[vel_i*6];")
    self.gen_add_code_line("T *psidd_p = &psidd_v[vel_i*6];")
    # Helpers — small (6-vec each). Parallel over 7 "helper_id" × 6 "r" = 42 elements.
    # Each thread computes one element of one helper. A5_vec/A7_vec depend on
    # other helpers so they're emitted in a second parallel_loop after a sync.
    self.gen_add_parallel_loop("h_idx", "42", use_thread_group)
    self.gen_add_code_lines([
        "int helper_id = h_idx / 6;",
        "int r = h_idx % 6;",
        "switch (helper_id) {",
        "  case 0: S_ICi_S[r]   = dot_prod<T, 6, 6, 1>(&IC[i*36 + r], S_p); break;",
        "  case 1: S_ICi_psid[r] = dot_prod<T, 6, 6, 1>(&IC[i*36 + r], psid_p); break;",
        "  case 2: S_ICi_psidd[r] = dot_prod<T, 6, 6, 1>(&IC[i*36 + r], psidd_p); break;",
        "  case 3: S_BCi_S[r]   = dot_prod<T, 6, 6, 1>(&BC[i*36 + r], S_p); break;",
        "  case 4: S_BCi_psid[r] = dot_prod<T, 6, 6, 1>(&BC[i*36 + r], psid_p); break;",
        "  case 5: S_BCiT_S[r]  = dot_prod<T, 6, 1, 1>(&BC[i*36 + 6*r], S_p); break;",
        "  case 6: {",
        "    T crf_S_row[6]; for (int kk = 0; kk < 6; ++kk) crf_S_row[kk] = -crm<T>(kk + 6*r, S_p);",
        "    S_crf_S_f_i[r] = dot_prod<T, 6, 1, 1>(crf_S_row, &f_w[i*6]);",
        "    break;",
        "  }",
        "}",
    ])
    self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)
    # A5_vec depends on S_BCi_psid + S_ICi_psidd + S_crf_S_f_i; A7_vec depends on S_BCi_S + IC@(psid+Sd).
    # Parallel over 12 = 6 (A5_vec) + 6 (A7_vec).
    self.gen_add_parallel_loop("v_idx", "12", use_thread_group)
    self.gen_add_code_lines([
        "int r = v_idx % 6;",
        "if (v_idx < 6) {",
        "    S_A5_vec[r] = S_BCi_psid[r] + S_ICi_psidd[r] + S_crf_S_f_i[r];",
        "} else {",
        "    T s_sum = static_cast<T>(0);",
        "    for (int kk = 0; kk < 6; ++kk) s_sum += IC[i*36 + r + 6*kk] * (psid_p[kk] + Sd_p[kk]);",
        "    S_A7_vec[r] = S_BCi_S[r] + s_sum;",
        "}",
    ])
    self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)
    # Parallel build of A0/A1/Bphi/Bpsid (each over idx ∈ [0, 36)).
    self.gen_add_parallel_loop("idx", "36", use_thread_group)
    self.gen_add_code_lines([
        "int row = idx % 6; int col = idx / 6;",
        "T crf_Sp_row[6]; for (int kk = 0; kk < 6; ++kk) crf_Sp_row[kk] = -crm<T>(kk + 6*row, S_p);",
        "T crm_Sp_col[6]; for (int kk = 0; kk < 6; ++kk) crm_Sp_col[kk] = crm<T>(kk + 6*col, S_p);",
        "T crf_psid_row[6]; for (int kk = 0; kk < 6; ++kk) crf_psid_row[kk] = -crm<T>(kk + 6*row, psid_p);",
        "T crm_psid_col[6]; for (int kk = 0; kk < 6; ++kk) crm_psid_col[kk] = crm<T>(kk + 6*col, psid_p);",
        "T t_crfSp_IC = dot_prod<T, 6, 1, 1>(crf_Sp_row, &IC[i*36 + 6*col]);",
        "T t_IC_crmSp = dot_prod<T, 6, 6, 1>(&IC[i*36 + row], crm_Sp_col);",
        "T t_crfpsid_IC = dot_prod<T, 6, 1, 1>(crf_psid_row, &IC[i*36 + 6*col]);",
        "T t_IC_crmpsid = dot_prod<T, 6, 6, 1>(&IC[i*36 + row], crm_psid_col);",
        "S_Bphi[idx]  = t_crfSp_IC  + icrf<T>(idx, S_ICi_S)    - t_IC_crmSp;",
        "S_Bpsid[idx] = t_crfpsid_IC + icrf<T>(idx, S_ICi_psid) - t_IC_crmpsid;",
        "S_A0[idx] = icrf<T>(idx, S_ICi_S);",
        "S_A1[idx] = t_crfSp_IC - t_IC_crmSp;",
    ])
    self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)
    # Parallel build of A2/A3/A4/A5/A6/A7 (depends on prior A0/A1/Bphi/Bpsid).
    self.gen_add_parallel_loop("idx", "36", use_thread_group)
    self.gen_add_code_lines([
        "int row = idx % 6; int col = idx / 6;",
        "// A2 = 2*A0 - Bphi",
        "S_A2[idx] = static_cast<T>(2) * S_A0[idx] - S_Bphi[idx];",
        "// A3 = Bpsid + dot_matrix(BC[i], S_p)",
        "T crf_Sp_row[6]; for (int kk = 0; kk < 6; ++kk) crf_Sp_row[kk] = -crm<T>(kk + 6*row, S_p);",
        "T crm_Sp_col[6]; for (int kk = 0; kk < 6; ++kk) crm_Sp_col[kk] = crm<T>(kk + 6*col, S_p);",
        "T t_crfSp_BC = dot_prod<T, 6, 1, 1>(crf_Sp_row, &BC[i*36 + 6*col]);",
        "T t_BC_crmSp = dot_prod<T, 6, 6, 1>(&BC[i*36 + row], crm_Sp_col);",
        "S_A3[idx] = S_Bpsid[idx] + t_crfSp_BC - t_BC_crmSp;",
        "// A4 = icrf(BCiT_S)",
        "S_A4[idx] = icrf<T>(idx, S_BCiT_S);",
        "// A5 = icrf(A5_vec)",
        "S_A5[idx] = icrf<T>(idx, S_A5_vec);",
        "// A6 = crf(S_p) @ IC[i][:, col] + A0[idx]",
        "T t_crfSp_IC = dot_prod<T, 6, 1, 1>(crf_Sp_row, &IC[i*36 + 6*col]);",
        "S_A6[idx] = t_crfSp_IC + S_A0[idx];",
        "// A7 = icrf(A7_vec)",
        "S_A7[idx] = icrf<T>(idx, S_A7_vec);",
    ])
    self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)

    # j-loop (all threads run sequentially through ancestor chain).
    self.gen_add_code_line("int j = i;")
    self.gen_add_code_line("while (j >= 0) {", True)
    self.gen_add_code_line("for (int tt = wf_body_v_start[j]; tt < wf_body_v_start[j + 1]; ++tt) {", True)
    self.gen_add_code_line("int vel_j = wf_body_v_index[tt];")
    # Per-(j, t) intermediates: parallel build of u1..u12 over 72 element-slots
    # (12 u-vectors × 6 elements each). Each thread builds one element of one u.
    self.gen_add_code_line("T *S_t     = &S_vel[vel_j*6];")
    self.gen_add_code_line("T *Sd_t    = &Sd_vel[vel_j*6];")
    self.gen_add_code_line("T *psid_t  = &psid_v[vel_j*6];")
    self.gen_add_code_line("T *psidd_t = &psidd_v[vel_j*6];")
    self.gen_add_parallel_loop("u_idx", "72", use_thread_group)
    self.gen_add_code_lines([
        "int which_u = u_idx / 6;",
        "int r = u_idx % 6;",
        "switch (which_u) {",
        "  case 0:  S_u1[r]  = dot_prod<T, 6, 1, 1>(&S_A3[6*r], S_t); break;",
        "  case 1:  S_u2[r]  = dot_prod<T, 6, 1, 1>(&S_A1[6*r], S_t); break;",
        "  case 2:  S_u3[r]  = dot_prod<T, 6, 6, 1>(&S_A3[r], psid_t) + dot_prod<T, 6, 6, 1>(&S_A1[r], psidd_t) + dot_prod<T, 6, 6, 1>(&S_A5[r], S_t); break;",
        "  case 3:  S_u4[r]  = dot_prod<T, 6, 6, 1>(&S_A6[r], S_t); break;",
        "  case 4:  S_u5[r]  = dot_prod<T, 6, 6, 1>(&S_A2[r], psid_t) + dot_prod<T, 6, 6, 1>(&S_A4[r], S_t); break;",
        "  case 5:  S_u6[r]  = dot_prod<T, 6, 6, 1>(&S_Bphi[r], psid_t) + dot_prod<T, 6, 6, 1>(&S_A7[r], S_t); break;",
        "  case 6: {",
        "    T psd_Sd[6]; for (int kk = 0; kk < 6; ++kk) psd_Sd[kk] = psid_t[kk] + Sd_t[kk];",
        "    S_u7[r] = dot_prod<T, 6, 6, 1>(&S_A3[r], S_t) + dot_prod<T, 6, 6, 1>(&S_A1[r], psd_Sd);",
        "    break;",
        "  }",
        "  case 7:  S_u8[r]  = dot_prod<T, 6, 6, 1>(&S_A4[r], S_t) - dot_prod<T, 6, 1, 1>(&S_Bphi[6*r], psid_t); break;",
        "  case 8:  S_u9[r]  = dot_prod<T, 6, 6, 1>(&S_A0[r], S_t); break;",
        "  case 9:  S_u10[r] = dot_prod<T, 6, 6, 1>(&S_Bphi[r], S_t); break;",
        "  case 10: S_u11[r] = dot_prod<T, 6, 1, 1>(&S_Bphi[6*r], S_t); break;",
        "  case 11: S_u12[r] = dot_prod<T, 6, 6, 1>(&S_A1[r], S_t); break;",
        "}",
    ])
    self.gen_add_end_control_flow()  # end parallel_loop u_idx
    self.gen_add_sync(use_thread_group)

    # ---- Parallel inner (k, rr) walk ----
    # Each thread is assigned a unique kr_idx in [0, total_kr). The thread chases
    # the ancestor chain of j to map kr_idx -> (k_local, rr_local, vel_k). All
    # writes to output tensors contain `vel_k` in the index, and vel_k values are
    # disjoint across kr_idx (each velocity belongs to exactly one body in the
    # ancestor chain), so threads never race on output cells.
    self.gen_add_code_line("// Flatten (k, rr) iterations of the ancestor chain of j into a single index space.")
    self.gen_add_code_line("int wf_total_kr = 0;")
    self.gen_add_code_line("for (int _kk = j; _kk >= 0; _kk = wf_parent[_kk]) wf_total_kr += wf_body_v_start[_kk + 1] - wf_body_v_start[_kk];")
    self.gen_add_parallel_loop("kr_idx", "wf_total_kr", use_thread_group)
    self.gen_add_code_lines([
        "// Map kr_idx -> (k, vel_k) by walking the chain.",
        "int k = -1; int vel_k = -1; int rr = -1;",
        "{",
        "int _seen = 0;",
        "for (int _kk = j; _kk >= 0; _kk = wf_parent[_kk]) {",
        "    int _nk = wf_body_v_start[_kk + 1] - wf_body_v_start[_kk];",
        "    if (kr_idx < _seen + _nk) {",
        "        rr = wf_body_v_start[_kk] + (kr_idx - _seen);",
        "        vel_k = wf_body_v_index[rr];",
        "        k = _kk;",
        "        break;",
        "    }",
        "    _seen += _nk;",
        "}",
        "}",
        "T *S_r     = &S_vel[vel_k*6];",
        "T *Sd_r    = &Sd_vel[vel_k*6];",
        "T *psid_r  = &psid_v[vel_k*6];",
        "T *psidd_r = &psidd_v[vel_k*6];",
        "T p1 = dot_prod<T, 6, 1, 1>(S_u11, psid_r);",
        "T p2 = dot_prod<T, 6, 1, 1>(S_u8, psid_r) + dot_prod<T, 6, 1, 1>(S_u9, psidd_r);",
        "d2tau_dq2[(vel_i*NUM_VEL + vel_j)*NUM_VEL + vel_k] = p2;",
        "d2tau_dvdq[(vel_i*NUM_VEL + vel_k)*NUM_VEL + vel_j] = -p1;",
        "",
        "T u1_psid_r  = dot_prod<T, 6, 1, 1>(S_u1, psid_r);",
        "T u2_psidd_r = dot_prod<T, 6, 1, 1>(S_u2, psidd_r);",
        "T u11_S_r    = dot_prod<T, 6, 1, 1>(S_u11, S_r);",
        "T u1_S_r     = dot_prod<T, 6, 1, 1>(S_u1, S_r);",
        "T u2_psd_Sd  = static_cast<T>(0);",
        "for (int kk = 0; kk < 6; ++kk) u2_psd_Sd += S_u2[kk] * (psid_r[kk] + Sd_r[kk]);",
        "T S_r_u3 = dot_prod<T, 6, 1, 1>(S_r, S_u3);",
        "T S_r_u4 = dot_prod<T, 6, 1, 1>(S_r, S_u4);",
        "T S_r_u5 = dot_prod<T, 6, 1, 1>(S_r, S_u5);",
        "T S_r_u6 = dot_prod<T, 6, 1, 1>(S_r, S_u6);",
        "T S_r_u7 = dot_prod<T, 6, 1, 1>(S_r, S_u7);",
        "T S_r_u9 = dot_prod<T, 6, 1, 1>(S_r, S_u9);",
        "T S_r_u10 = dot_prod<T, 6, 1, 1>(S_r, S_u10);",
        "T S_r_u12 = dot_prod<T, 6, 1, 1>(S_r, S_u12);",
        "T u9_psd_Sd = static_cast<T>(0);",
        "for (int kk = 0; kk < 6; ++kk) u9_psd_Sd += S_u9[kk] * (psid_r[kk] + Sd_r[kk]);",
        "",
        "if (j != i) {", True,
        "T dq_jki = u1_psid_r + u2_psidd_r;",
        "d2tau_dq2[(vel_j*NUM_VEL + vel_k)*NUM_VEL + vel_i] = dq_jki;",
        "d2tau_dq2[(vel_j*NUM_VEL + vel_i)*NUM_VEL + vel_k] = dq_jki;",
        "d2tau_dvdq[(vel_j*NUM_VEL + vel_k)*NUM_VEL + vel_i] = p1;",
        "d2tau_dvdq[(vel_j*NUM_VEL + vel_i)*NUM_VEL + vel_k] = u1_S_r + u2_psd_Sd;",
        "d2tau_dqd2[(vel_j*NUM_VEL + vel_k)*NUM_VEL + vel_i] = u11_S_r;",
        "d2tau_dqd2[(vel_j*NUM_VEL + vel_i)*NUM_VEL + vel_k] = u11_S_r;",
        "dM_dq[(vel_k*NUM_VEL + vel_j)*NUM_VEL + vel_i] = S_r_u12;",
        "dM_dq[(vel_j*NUM_VEL + vel_k)*NUM_VEL + vel_i] = S_r_u12;",
    ])
    self.gen_add_end_control_flow()
    self.gen_add_code_lines([
        "if (k != j) {", True,
        "d2tau_dq2[(vel_i*NUM_VEL + vel_k)*NUM_VEL + vel_j] = p2;",
        "d2tau_dq2[(vel_k*NUM_VEL + vel_i)*NUM_VEL + vel_j] = S_r_u3;",
        "d2tau_dqd2[(vel_i*NUM_VEL + vel_j)*NUM_VEL + vel_k] = -u11_S_r;",
        "d2tau_dqd2[(vel_i*NUM_VEL + vel_k)*NUM_VEL + vel_j] = -u11_S_r;",
        "d2tau_dvdq[(vel_i*NUM_VEL + vel_j)*NUM_VEL + vel_k] = S_r_u5 + u9_psd_Sd;",
        "d2tau_dvdq[(vel_k*NUM_VEL + vel_j)*NUM_VEL + vel_i] = S_r_u6;",
        "dM_dq[(vel_k*NUM_VEL + vel_i)*NUM_VEL + vel_j] = S_r_u9;",
        "dM_dq[(vel_i*NUM_VEL + vel_k)*NUM_VEL + vel_j] = S_r_u9;",
        "if (j != i) {", True,
        "d2tau_dq2[(vel_k*NUM_VEL + vel_j)*NUM_VEL + vel_i] = S_r_u3;",
        "d2tau_dqd2[(vel_k*NUM_VEL + vel_i)*NUM_VEL + vel_j] = S_r_u10;",
        "d2tau_dqd2[(vel_k*NUM_VEL + vel_j)*NUM_VEL + vel_i] = S_r_u10;",
        "d2tau_dvdq[(vel_k*NUM_VEL + vel_i)*NUM_VEL + vel_j] = S_r_u7;",
    ])
    self.gen_add_end_control_flow()
    self.gen_add_code_line("else {", True)
    self.gen_add_code_line("d2tau_dqd2[(vel_k*NUM_VEL + vel_j)*NUM_VEL + vel_i] = S_r_u4;")
    self.gen_add_end_control_flow()
    self.gen_add_end_control_flow()
    self.gen_add_code_line("else {", True)
    self.gen_add_code_line("d2tau_dqd2[(vel_i*NUM_VEL + vel_j)*NUM_VEL + vel_k] = -dot_prod<T, 6, 1, 1>(S_u2, S_r);")
    self.gen_add_end_control_flow()

    self.gen_add_end_control_flow()  # end parallel_loop over kr_idx
    self.gen_add_sync(use_thread_group)

    self.gen_add_end_control_flow()  # end tt loop
    self.gen_add_code_line("j = wf_parent[j];")
    self.gen_add_end_control_flow()  # end while j

    self.gen_add_end_control_flow()  # end pp loop

    # Aggregate IC, BC, f into parent — parallel over 78 elements (36 IC + 36 BC + 6 f).
    self.gen_add_code_line("// Bubble subtree-aggregated IC/BC/f up — parallel over 78 elements per body.")
    self.gen_add_code_line("int parent = wf_parent[i];")
    self.gen_add_code_line("if (parent >= 0) {", True)
    self.gen_add_parallel_loop("agg_idx", "78", use_thread_group)
    self.gen_add_code_lines([
        "if (agg_idx < 36) {",
        "    IC[parent*36 + agg_idx] += IC[i*36 + agg_idx];",
        "} else if (agg_idx < 72) {",
        "    int b = agg_idx - 36;",
        "    BC[parent*36 + b] += BC[i*36 + b];",
        "} else {",
        "    int r = agg_idx - 72;",
        "    f_w[parent*6 + r] += f_w[i*6 + r];",
        "}",
    ])
    self.gen_add_end_control_flow()  # end parallel_loop agg_idx
    self.gen_add_end_control_flow()  # end if (parent >= 0)
    self.gen_add_sync(use_thread_group)
    self.gen_add_end_control_flow()  # end i loop

    # Final transpose of d2tau_dvdq trailing axes: [τ, q, qd] -> [τ, qd, q].
    # Parallelize over (a_i, b_i) with each thread handling its (b_i, c_i) upper-triangle pairs.
    self.gen_add_code_line("// d2tau_dvdq was stored [τ, q, qd]; transpose trailing axes to match RBDReference convention.")
    self.gen_add_parallel_loop("ab_i", "NUM_VEL*NUM_VEL", use_thread_group)
    self.gen_add_code_line("int a_i = ab_i / NUM_VEL; int b_i = ab_i % NUM_VEL;")
    self.gen_add_code_line("for (int c_i = b_i + 1; c_i < NUM_VEL; ++c_i) {", True)
    self.gen_add_code_line("T x = d2tau_dvdq[(a_i*NUM_VEL + b_i)*NUM_VEL + c_i];")
    self.gen_add_code_line("T y = d2tau_dvdq[(a_i*NUM_VEL + c_i)*NUM_VEL + b_i];")
    self.gen_add_code_line("d2tau_dvdq[(a_i*NUM_VEL + b_i)*NUM_VEL + c_i] = y;")
    self.gen_add_code_line("d2tau_dvdq[(a_i*NUM_VEL + c_i)*NUM_VEL + b_i] = x;")
    self.gen_add_end_control_flow()
    self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)
    self.gen_add_end_function()


def gen_idsva_so_world_frame_inner_function_call(self, use_thread_group = False):
    """Emit the call to `idsva_so_world_frame_inner` mirroring the existing call helper."""
    id_so_code_start = "idsva_so_world_frame_inner<T>(s_idsva_so, s_q, s_qd, s_qdd, "
    id_so_code_middle = self.gen_insert_helpers_function_call()
    id_so_code_end = "s_temp, gravity);"
    if use_thread_group:
        id_so_code_start = id_so_code_start.replace("(", "(tgrp, ")
    self.gen_add_code_line(id_so_code_start + id_so_code_middle + id_so_code_end)


def gen_idsva_so_world_frame_kernel(self, use_thread_group = False, single_call_timing = False):
    NUM_POS = self.robot.get_num_pos()
    n = self.robot.get_num_vel()
    NJ = self.robot.get_num_joints()
    # sv2 has its own smaller shared-mem budget (`IDSVA_SO_WORLD_FRAME_DYNAMIC_SHARED_MEM_BYTES`).
    # Use its own use_global_output flag — typically False even for g1 since sv2 is much smaller.
    use_global_output = getattr(self, "idsva_so_world_frame_use_global_output",
                                getattr(self, "idsva_so_body_frame_use_global_output", NJ > SHARED_MEMORY_JOINT_THRESHOLD))
    func_params = [
        "d_idsva_so is a pointer to memory for the final result of size 4*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS*SECOND_ORDER_COORDS = " + str(4*n**3),
        "d_q_dq_u is the vector of joint positions, velocities, and accelerations",
        "stride_q_qd_u is the stride between each q, qd, u",
        "d_robotModel is the pointer to the initialized model specific helpers on the GPU (XImats, topology_helpers, etc.)",
        "gravity is the gravity constant",
        "num_timesteps is the length of the trajectory points",
    ]
    func_notes = ["world-frame IDSVA-SO kernel: clean single-pass reference path."]
    func_def_start = "void idsva_so_world_frame_kernel(T *d_idsva_so, const T *d_q_qd_u, const int stride_q_qd_u, "
    func_def_end = "const robotModel<T> *d_robotModel, const T gravity, const int NUM_TIMESTEPS) {"
    func_def = func_def_start + func_def_end
    if single_call_timing:
        func_def = func_def.replace("kernel(", "kernel_single_timing(")
    self.gen_add_func_doc("Computes IDSVA-SO via the world-frame single-pass formulation", func_notes, func_params, None)
    self.gen_add_code_line("template <typename T, int RESOURCE_TIER = TIER_PERF>")
    self.gen_add_code_line("__global__")
    self.gen_add_code_line("__launch_bounds__(tier_max_threads<RESOURCE_TIER>())")
    self.gen_add_code_line(func_def, True)
    extra_t_buffers = [("s_q_qd_u", n*2 + NUM_POS)]
    if not use_global_output:
        extra_t_buffers.append(("s_idsva_so", 4*n**3))
    # sv2 uses its own (smaller) inner temp size.
    shared_mem_size = gen_idsva_so_world_frame_temp_mem_size(self) if self.robot.floating_base else self.gen_idsva_so_body_frame_inner_temp_mem_size()
    self.gen_XImats_helpers_temp_shared_memory_code(shared_mem_size, extra_t_buffers=extra_t_buffers)
    self.gen_add_code_line(f"T *s_q = s_q_qd_u; T *s_qd = &s_q_qd_u[{NUM_POS}]; T *s_qdd = &s_q_qd_u[{NUM_POS + n}];")
    if not single_call_timing:
        self.gen_add_parallel_loop("k", "NUM_TIMESTEPS", use_thread_group, block_level=True)
        self.gen_kernel_load_inputs("q_qd_u", "stride_q_qd_u", str(2*n + NUM_POS), use_thread_group)
        self.gen_load_update_XImats_helpers_function_call(use_thread_group)
        if use_global_output:
            self.gen_add_code_line(f"T *s_idsva_so = &d_idsva_so[k*{4*n**3}];")
        self.gen_idsva_so_world_frame_inner_function_call(use_thread_group)
        if not use_global_output:
            self.gen_kernel_save_result("idsva_so", str(4*n**3), str(4*n**3), use_thread_group)
        self.gen_add_end_control_flow()
    else:
        self.gen_kernel_load_inputs_single_timing("q_qd_u", str(NUM_POS + 2*n), use_thread_group)
        self.gen_add_code_line("// compute with NUM_TIMESTEPS as NUM_REPS for timing")
        self.gen_add_code_line("for (int rep = 0; rep < NUM_TIMESTEPS; rep++){", True)
        self.gen_anti_licm_input_reload("q_qd_u", str(NUM_POS + 2*n), use_thread_group)
        self.gen_load_update_XImats_helpers_function_call(use_thread_group)
        if use_global_output:
            self.gen_add_code_line("T *s_idsva_so = d_idsva_so;")
        self.gen_idsva_so_world_frame_inner_function_call(use_thread_group)
        self.gen_add_end_control_flow()
        if not use_global_output:
            self.gen_kernel_save_result_single_timing("idsva_so", str(4*n**3), use_thread_group)
    self.gen_add_end_function()


def gen_idsva_so_world_frame_host(self, mode = 0):
    single_call_timing = True if mode == 1 else False
    compute_only = True if mode == 2 else False
    func_params = [
        "hd_data is the packaged input and output pointers",
        "d_robotModel is the pointer to the initialized model specific helpers on the GPU",
        "gravity is the gravity constant",
        "num_timesteps is the length of the trajectory points",
        "streams are pointers to CUDA streams",
    ]
    func_def_start = "void idsva_so_world_frame_host(gridData<T, KIND> *hd_data, const robotModel<T> *d_robotModel, const T gravity, const int num_timesteps,"
    func_def_end =   "                      const dim3 block_dimms, const dim3 thread_dimms, cudaStream_t *streams) {"
    if single_call_timing:
        func_def_start = func_def_start.replace("(", "_single_timing(")
        func_def_end = "              " + func_def_end
    if compute_only:
        func_def_start = func_def_start.replace("(", "_compute_only(")
        func_def_end = "             " + func_def_end.replace(", cudaStream_t *streams", "")
    self.gen_add_func_doc("Compute IDSVA-SO via the world-frame single-pass formulation", [], func_params, None)
    self.gen_add_code_line("template <typename T, gridDataKind KIND = GRID_DATA_ALL>")
    self.gen_add_code_line("__host__")
    self.gen_add_code_line(func_def_start)
    self.gen_add_code_line(func_def_end, True)
    self.gen_add_code_line("static_assert(KIND == GRID_DATA_ALL || KIND == GRID_DATA_DYNAMICS, \"idsva_so_world_frame_host requires all-data or dynamics gridData\");")
    func_call_start = "idsva_so_world_frame_kernel<T><<<block_dimms,thread_dimms,IDSVA_SO_WORLD_FRAME_DYNAMIC_SHARED_MEM_BYTES<T>()>>>(hd_data->d_idsva_so," + \
        "hd_data->d_q_qd_u,stride_q_qd,"
    func_call_end = "d_robotModel,gravity,num_timesteps);"
    if single_call_timing:
        func_call_start = func_call_start.replace("kernel<T>", "kernel_single_timing<T>")
    self.gen_add_code_line("int stride_q_qd = Q_QD_U_STRIDE;")
    if not compute_only:
        self.gen_add_code_lines([
            "// start code with memory transfer",
            "gpuErrchk(cudaMemcpyAsync(hd_data->d_q_qd_u,hd_data->h_q_qd_u,stride_q_qd*" + ("num_timesteps*" if not single_call_timing else "") + "sizeof(T),cudaMemcpyHostToDevice,streams[0]));",
            "gpuErrchkKernel();",
        ])
    self.gen_add_code_line("// then call the kernel")
    func_call_code = [f'{func_call_start}{func_call_end}']
    # See gen_idsva_so_body_frame_host for the same fix: sync between launch and
    # clock_gettime(end) is required for real single-call timing, and
    # compute_only needs a sync so callers' batch timers see actual
    # kernel-completion time (not just async-launch overhead).
    if single_call_timing:
        func_call_code.insert(0, "struct timespec start, end; clock_gettime(CLOCK_MONOTONIC,&start);")
        func_call_code.append("gpuErrchkKernel();")
        func_call_code.append("clock_gettime(CLOCK_MONOTONIC,&end);")
    self.gen_add_code_line("gpuErrchk(grid_check_dynamic_shared_memory_bytes(\"idsva_so_world_frame\", IDSVA_SO_WORLD_FRAME_DYNAMIC_SHARED_MEM_BYTES<T>()));")
    self.gen_add_code_lines(func_call_code)
    if not compute_only:
        self.gen_add_code_lines([
            "// finally transfer the result back",
            "gpuErrchk(cudaMemcpy(hd_data->h_idsva_so,hd_data->d_idsva_so,SECOND_ORDER_TENSOR_SIZE*" + ("num_timesteps*" if not single_call_timing else "") + "sizeof(T),cudaMemcpyDeviceToHost));",
            "gpuErrchkKernel();",
        ])
    else:
        self.gen_add_code_line("gpuErrchkKernel();")
    # Label kept distinct from the original IDSVA_SO so a future bench
    # that times both can keep them separate. The parser doesn't know
    # this label today; if/when the bench wires it up, add a
    # _GRID_SINGLE_LABELS entry.
    if single_call_timing:
        from ..algo_registry import single_call_printf_line
        self.gen_add_code_line(single_call_printf_line("idsva_so_world_frame"))
    self.gen_add_end_function()


def gen_idsva_so_world_frame(self, use_thread_group = False):
    """Emit the complete world-frame IDSVA-SO path: inner, kernel, host wrappers.

    Co-exists with the existing `gen_idsva_so_body_frame` emission. Gated by the
    `enable_idsva_so_world_frame` flag in `gen_all_code`.
    """
    self.gen_idsva_so_world_frame_inner(use_thread_group)
    self.gen_idsva_so_world_frame_kernel(use_thread_group, single_call_timing=False)
    self.gen_idsva_so_world_frame_kernel(use_thread_group, single_call_timing=True)
    self.gen_idsva_so_world_frame_host(0)
    self.gen_idsva_so_world_frame_host(1)
    self.gen_idsva_so_world_frame_host(2)


def gen_idsva_so_dispatcher_host(self, mode = 0):
    """Emit `grid::idsva_so` — a host wrapper that calls the perf-winning
    variant for this robot's base type. Picked at codegen time: body_frame
    for fixed-base, world_frame for floating-base. Both inners produce
    numerically equivalent output; this is purely a perf optimization.
    Body/world host wrappers remain individually callable for direct
    comparison.

    Measured on sm_120 / RTX 5090 (2026-05-18 sweep):
      - iiwa14 (NV=7,  fixed):    body 7-9x   faster than world
      - go2    (NV=12, fixed):    body 7-10x  faster than world
      - g1     (NV=29, fixed):    world 6-15% faster than body
      - iiwa14 (NV=6,  floating): world 7x   faster than body
      - go2    (NV=18, floating): world 7x   faster than body
      - g1     (NV=35, floating): world 20x  faster than body

    The "body for fixed, world for floating" rule is the safe choice — it
    preserves the large wins at the common low-DOF fixed-base case
    (iiwa14, go2) and the large wins on every floating-base case. The
    g1_fixed regression is small (~15%) and isolated to a single
    high-DOF data point; refining the dispatcher with a NV threshold is a
    worthwhile follow-up once more high-DOF fixed-base robots exist in
    the manifest.

    Regular and compute_only modes forward to the underlying host wrapper.
    Single-timing mode inlines its own clock_gettime + kernel launch so the
    printf label says "IDSVA_SO" (not "IDSVA_SO_BODY_FRAME"/"IDSVA_SO_WORLD_FRAME"),
    which is what the bench's timing parser keys on for this row.
    """
    single_call_timing = (mode == 1)
    compute_only = (mode == 2)

    frame_suffix = "world_frame" if self.robot.floating_base else "body_frame"
    smem_macro = f"IDSVA_SO_{frame_suffix.upper()}_DYNAMIC_SHARED_MEM_BYTES<T>()"

    func_def_start = "void idsva_so(gridData<T, KIND> *hd_data, const robotModel<T> *d_robotModel, const T gravity, const int num_timesteps,"
    func_def_end =   "                      const dim3 block_dimms, const dim3 thread_dimms, cudaStream_t *streams) {"
    if single_call_timing:
        func_def_start = func_def_start.replace("(", "_single_timing(")
        func_def_end = "              " + func_def_end
    if compute_only:
        func_def_start = func_def_start.replace("(", "_compute_only(")
        func_def_end = "             " + func_def_end.replace(", cudaStream_t *streams", "")

    self.gen_add_func_doc(
        f"Dispatching IDSVA-SO wrapper (forwards to {frame_suffix} at codegen time)",
        [f"Body wins ~30x on fixed-base; world wins 2-4x on floating-base. Same numbers either way."],
        ["hd_data is the packaged input and output pointers",
         "d_robotModel is the pointer to the initialized model specific helpers on the GPU (XImats, topology_helpers, etc.)",
         "gravity is the gravity constant",
         "num_timesteps is the length of the trajectory points we need to compute over (or overloaded as test_iters for timing)",
         "streams are pointers to CUDA streams for async memory transfers (if needed)"],
        None,
    )
    self.gen_add_code_line("template <typename T, gridDataKind KIND = GRID_DATA_ALL>")
    self.gen_add_code_line("__host__")
    self.gen_add_code_line(func_def_start)
    self.gen_add_code_line(func_def_end, True)

    if single_call_timing:
        # Inline our own single-timing block so the printf label is "IDSVA_SO"
        # (rather than the underlying body_frame/world_frame label that the
        # delegate's _single_timing wrapper would print). Mirrors the structure
        # of `gen_idsva_so_body_frame_host(mode=1)`.
        # body_frame_kernel takes d_workspace (grav-shim spill); world_frame_kernel
        # does not (no spill needed). Match the underlying signature.
        kernel_name = f"idsva_so_{frame_suffix}_kernel_single_timing"
        kernel_workspace_arg = "" if self.robot.floating_base else "hd_data->d_workspace,"
        self.gen_add_code_line("static_assert(KIND == GRID_DATA_ALL || KIND == GRID_DATA_DYNAMICS, \"idsva_so requires all-data or dynamics gridData\");")
        self.gen_add_code_line("int stride_q_qd = Q_QD_U_STRIDE;")
        self.gen_add_code_lines([
            "// start code with memory transfer",
            "gpuErrchk(cudaMemcpyAsync(hd_data->d_q_qd_u,hd_data->h_q_qd_u,stride_q_qd*sizeof(T),cudaMemcpyHostToDevice,streams[0]));",
            "gpuErrchkKernel();",
        ])
        self.gen_add_code_line("// then call the kernel")
        self.gen_add_code_line(f"gpuErrchk(grid_check_dynamic_shared_memory_bytes(\"idsva_so\", {smem_macro}));")
        self.gen_add_code_line("struct timespec start, end; clock_gettime(CLOCK_MONOTONIC,&start);")
        self.gen_add_code_line(
            f"{kernel_name}<T><<<block_dimms,thread_dimms,{smem_macro}>>>(hd_data->d_idsva_so,{kernel_workspace_arg}hd_data->d_q_qd_u,stride_q_qd,d_robotModel,gravity,num_timesteps);"
        )
        self.gen_add_code_line("gpuErrchkKernel();")
        self.gen_add_code_line("clock_gettime(CLOCK_MONOTONIC,&end);")
        self.gen_add_code_lines([
            "// finally transfer the result back",
            "gpuErrchk(cudaMemcpy(hd_data->h_idsva_so,hd_data->d_idsva_so,SECOND_ORDER_TENSOR_SIZE*sizeof(T),cudaMemcpyDeviceToHost));",
            "gpuErrchkKernel();",
        ])
        from ..algo_registry import single_call_printf_line
        self.gen_add_code_line(single_call_printf_line("idsva_so"))
    else:
        # Regular and compute_only modes forward to the underlying host wrapper.
        # The underlying wrapper handles all memcpy + launch + sync correctly,
        # and neither mode emits a single-call printf label — so the label
        # collision doesn't apply here.
        target = f"idsva_so_{frame_suffix}_host"
        if compute_only:
            target += "_compute_only"
        forward_args = "hd_data, d_robotModel, gravity, num_timesteps, block_dimms, thread_dimms"
        if not compute_only:
            forward_args += ", streams"
        self.gen_add_code_line(f"{target}<T, KIND>({forward_args});")
    self.gen_add_end_function()


def gen_idsva_so_dispatcher(self):
    """Emit all three dispatcher host modes."""
    self.gen_idsva_so_dispatcher_host(0)
    self.gen_idsva_so_dispatcher_host(1)
    self.gen_idsva_so_dispatcher_host(2)
