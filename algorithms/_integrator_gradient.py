"""Time-integrator gradient codegen.

Computes `dAB = d x_{k+1} / d(x, u)` as a 2n × 3n matrix (column-major)
where x = [q; qd] (size 2n) and u (size n). For Euler:

    top n rows of dAB:
        d(q_kp1) / dq    = I_n
        d(q_kp1) / dqd   = dt * I_n
        d(q_kp1) / du    = 0
    bottom n rows of dAB:
        d(qd_kp1) / dq   = dt * dqdd/dq
        d(qd_kp1) / dqd  = I_n + dt * dqdd/dqd
        d(qd_kp1) / du   = dt * dqdd/du  (and dqdd/du = Minv)

The `dqdd/dq`, `dqdd/dqd` blocks come straight from the FD gradient's
`s_df_du` output; `Minv` is the SYMMETRIC_UPPER triangular n×n already
populated by `forward_dynamics_gradient_inner_python`.

A boolean compile-time flag `COMPUTE_X_KP1` requests the value too —
when true the kernel also writes `d_x_kp1` (no extra RBD work, just an
extra parallel loop using `s_qdd` and `s_q`/`s_qd`).
"""

from ._integrator import _integrator_type_token, _max_stages_in_use


# Per-integrator Butcher coefficients used by the multi-stage gradient.
# For each multi-stage IT, we record:
#   c_i: stage offsets used to build the i+1-th point (p_{i+1} = x + c_i*dt*xdot_{i+1})
#        — TrajoptPlant convention so c is one entry per stage transition, length N-1.
#   b_i: final combination weights (length N).
# Stage 1 uses no offset (we set "c_0 = 0" in the unified formula so the chain
# rule reduces correctly).
_INTEGRATOR_BUTCHER = {
    # IT: (stage_count, [c_1..c_{N-1}], [b_1..b_N])
    "MIDPOINT": (2, [0.5],            [0.0, 1.0]),
    "RK3":      (3, [0.5, 0.75],      [2.0/9.0, 3.0/9.0, 4.0/9.0]),
    "RK4":      (4, [0.5, 0.5, 1.0],  [1.0/6.0, 2.0/6.0, 2.0/6.0, 1.0/6.0]),
}


def gen_integrator_gradient_inner_temp_mem_size(self):
    # Identical to FD-gradient's inner mem requirement; the dAB assembly is
    # a single parallel loop over shared inputs that already exist.
    return self.gen_forward_dynamics_gradient_inner_temp_mem_size()


def gen_integrator_gradient_dAB_assembly(self, integrator_type="IT", use_thread_group=False,
                                          s_dAB_name="s_dAB",
                                          s_df_du_name="s_df_du",
                                          s_Minv_name="s_Minv"):
    """Emit the parallel loop that fills s_dAB from s_df_du and s_Minv.

    Layout: s_dAB is 2n × 3n in COLUMN-major (matches s_df_du's convention).
    Address: s_dAB[col * 2n + row].
    """
    n = self.robot.get_num_vel()
    twoN = 2 * n
    nn = n * n
    self.gen_add_parallel_loop("ind", str(twoN * 3 * n), use_thread_group)
    self.gen_add_code_line("int row = ind % " + str(twoN) + ";")
    self.gen_add_code_line("int col = ind / " + str(twoN) + ";")
    tok = _integrator_type_token(integrator_type)
    # ----- EULER -----
    self.gen_add_code_line("if constexpr (" + tok + " == IntegratorType::EULER) {", True)
    self.gen_add_code_line("T val = static_cast<T>(0);")
    self.gen_add_code_line("if (col < " + str(n) + ") {")
    self.gen_add_code_line("    // d/dq column")
    self.gen_add_code_line("    if (row < " + str(n) + ") {")
    self.gen_add_code_line("        val = (row == col) ? static_cast<T>(1) : static_cast<T>(0);")
    self.gen_add_code_line("    } else {")
    self.gen_add_code_line("        int i_local = row - " + str(n) + ";")
    self.gen_add_code_line("        val = dt * " + s_df_du_name + "[col * " + str(n) + " + i_local];")
    self.gen_add_code_line("    }")
    self.gen_add_code_line("} else if (col < " + str(2 * n) + ") {")
    self.gen_add_code_line("    // d/dqd column")
    self.gen_add_code_line("    int j_local = col - " + str(n) + ";")
    self.gen_add_code_line("    if (row < " + str(n) + ") {")
    self.gen_add_code_line("        val = (row == j_local) ? dt : static_cast<T>(0);")
    self.gen_add_code_line("    } else {")
    self.gen_add_code_line("        int i_local = row - " + str(n) + ";")
    self.gen_add_code_line("        T diag = (i_local == j_local) ? static_cast<T>(1) : static_cast<T>(0);")
    self.gen_add_code_line("        val = diag + dt * " + s_df_du_name + "[" + str(nn) + " + j_local * " + str(n) + " + i_local];")
    self.gen_add_code_line("    }")
    self.gen_add_code_line("} else {")
    self.gen_add_code_line("    // d/du column   (dqdd/du = Minv)")
    self.gen_add_code_line("    int j_local = col - " + str(2 * n) + ";")
    self.gen_add_code_line("    if (row < " + str(n) + ") {")
    self.gen_add_code_line("        val = static_cast<T>(0);")
    self.gen_add_code_line("    } else {")
    self.gen_add_code_line("        int i_local = row - " + str(n) + ";")
    self.gen_add_code_line("        int midx = (i_local <= j_local) * (j_local * " + str(n) + " + i_local) + (i_local > j_local) * (i_local * " + str(n) + " + j_local);")
    self.gen_add_code_line("        val = dt * " + s_Minv_name + "[midx];")
    self.gen_add_code_line("    }")
    self.gen_add_code_line("}")
    self.gen_add_code_line(s_dAB_name + "[ind] = val;")
    self.gen_add_end_control_flow()  # end if constexpr EULER
    # ----- SEMI-IMPLICIT EULER -----
    # v_{k+1} = v + dt * qdd
    # q_{k+1} = q + dt * v_{k+1} = q + dt*v + dt^2 * qdd
    # Gradient:
    #   d(q_kp1)/dq  = I + dt^2 * dqdd/dq
    #   d(q_kp1)/dqd = dt * I + dt^2 * dqdd/dqd
    #   d(q_kp1)/du  = dt^2 * dqdd/du = dt^2 * Minv
    #   d(qd_kp1)/dq  = dt * dqdd/dq
    #   d(qd_kp1)/dqd = I + dt * dqdd/dqd
    #   d(qd_kp1)/du  = dt * dqdd/du = dt * Minv
    self.gen_add_code_line("else if constexpr (" + tok + " == IntegratorType::SEMI_IMPLICIT_EULER) {", True)
    self.gen_add_code_line("T val = static_cast<T>(0);")
    self.gen_add_code_line("T dt2 = dt * dt;")
    self.gen_add_code_line("if (col < " + str(n) + ") {")
    self.gen_add_code_line("    // d/dq column")
    self.gen_add_code_line("    if (row < " + str(n) + ") {")
    self.gen_add_code_line("        T diag = (row == col) ? static_cast<T>(1) : static_cast<T>(0);")
    self.gen_add_code_line("        val = diag + dt2 * " + s_df_du_name + "[col * " + str(n) + " + row];")
    self.gen_add_code_line("    } else {")
    self.gen_add_code_line("        int i_local = row - " + str(n) + ";")
    self.gen_add_code_line("        val = dt * " + s_df_du_name + "[col * " + str(n) + " + i_local];")
    self.gen_add_code_line("    }")
    self.gen_add_code_line("} else if (col < " + str(2 * n) + ") {")
    self.gen_add_code_line("    // d/dqd column")
    self.gen_add_code_line("    int j_local = col - " + str(n) + ";")
    self.gen_add_code_line("    if (row < " + str(n) + ") {")
    self.gen_add_code_line("        T diag = (row == j_local) ? dt : static_cast<T>(0);")
    self.gen_add_code_line("        val = diag + dt2 * " + s_df_du_name + "[" + str(nn) + " + j_local * " + str(n) + " + row];")
    self.gen_add_code_line("    } else {")
    self.gen_add_code_line("        int i_local = row - " + str(n) + ";")
    self.gen_add_code_line("        T diag = (i_local == j_local) ? static_cast<T>(1) : static_cast<T>(0);")
    self.gen_add_code_line("        val = diag + dt * " + s_df_du_name + "[" + str(nn) + " + j_local * " + str(n) + " + i_local];")
    self.gen_add_code_line("    }")
    self.gen_add_code_line("} else {")
    self.gen_add_code_line("    // d/du column   (dqdd/du = Minv)")
    self.gen_add_code_line("    int j_local = col - " + str(2 * n) + ";")
    self.gen_add_code_line("    if (row < " + str(n) + ") {")
    self.gen_add_code_line("        int midx = (row <= j_local) * (j_local * " + str(n) + " + row) + (row > j_local) * (row * " + str(n) + " + j_local);")
    self.gen_add_code_line("        val = dt2 * " + s_Minv_name + "[midx];")
    self.gen_add_code_line("    } else {")
    self.gen_add_code_line("        int i_local = row - " + str(n) + ";")
    self.gen_add_code_line("        int midx = (i_local <= j_local) * (j_local * " + str(n) + " + i_local) + (i_local > j_local) * (i_local * " + str(n) + " + j_local);")
    self.gen_add_code_line("        val = dt * " + s_Minv_name + "[midx];")
    self.gen_add_code_line("    }")
    self.gen_add_code_line("}")
    self.gen_add_code_line(s_dAB_name + "[ind] = val;")
    self.gen_add_end_control_flow()  # end if constexpr SI_EULER
    self.gen_add_code_line("else {", True)
    self.gen_add_code_line("static_assert(" + tok + " == IntegratorType::EULER || " + tok + " == IntegratorType::SEMI_IMPLICIT_EULER,")
    self.gen_add_code_line("              \"Integrator gradient type not yet implemented (Midpoint/RK3/RK4 need multi-stage chain rule).\");")
    self.gen_add_end_control_flow()
    self.gen_add_end_control_flow()  # end parallel loop


def gen_integrator_gradient_multistage(self, use_thread_group=False, compute_x_kp1=False):
    """Emit the multi-stage gradient body inline.

    Drives N stages of forward-dynamics-gradient at intermediate states with
    the TrajoptPlant point construction (p_{i+1} = x + c_i*dt*xdot_{i+1}; xdot
    has the original v in its first n slots, so q part of p only depends on
    original (q, qd)). Computes per-stage `D_qdd_i = ∂qdd_i / ∂(q, qd, u)` of
    shape (n × 3n, column-major) via the unified recurrence:

        D_qdd_1     = [J_qq_1 | J_qv_1 | Minv_1]                  (no chain)
        D_qdd_{i+1} = base_{i+1}(c_i, dt) + c_i*dt * J_qv_{i+1} @ D_qdd_i

    where base_{i+1} has block structure (per column c, block index = c // n):
        block_q  : J_qq_{i+1}[:, c]
        block_qd : c_i*dt * J_qq_{i+1}[:, c-n] + J_qv_{i+1}[:, c-n]
        block_u  : Minv_{i+1}[:, c-2n]
    (For stage 1, c_0 = 0 unifies the recurrence: chain term vanishes and
    base reduces to [J_qq_1 | J_qv_1 | Minv_1].)

    Mutates s_q, s_qd in shared memory across stages — caller must use the
    kernel-level non-const pointers (the per-stage XImats helper is also
    re-derived from the freshly-mutated s_q).
    """
    n = self.robot.get_num_vel()
    max_stages = _max_stages_in_use()
    three_n = 3 * n
    nn = n * n

    # Save original q, qd so we can rebuild p_{i+1} on later stages.
    self.gen_add_code_line("// --- multi-stage gradient: save original q, qd ---")
    self.gen_add_parallel_loop("ind", str(n), use_thread_group)
    self.gen_add_code_line("s_q_orig[ind] = s_q[ind];")
    self.gen_add_code_line("s_qd_orig[ind] = s_qd[ind];")
    self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)

    for stage_idx in range(max_stages):
        stage_num = stage_idx + 1  # 1-indexed for readability
        # Gate stages beyond what each integrator type uses.
        gating = " || ".join(
            "IT == IntegratorType::" + name
            for name, (cnt, _, _) in _INTEGRATOR_BUTCHER.items() if cnt >= stage_num
        )
        self.gen_add_code_line(f"// --- multi-stage gradient: stage {stage_num} ---")
        self.gen_add_code_line(f"if constexpr ({gating}) {{", True)

        if stage_idx > 0:
            # Need to overwrite s_q, s_qd with p_i values and update XImats.
            # Each multi-stage IT may use a different c_{stage_idx-1} value;
            # the one selected here depends on IT.
            offset_branches = []
            for name, (cnt, c_list, _) in _INTEGRATOR_BUTCHER.items():
                if cnt >= stage_num:
                    offset_branches.append(
                        f"(IT == IntegratorType::{name}) ? static_cast<T>({c_list[stage_idx - 1]})"
                    )
            self.gen_add_code_line(
                "constexpr T c_offset = " + " : ".join(offset_branches) + " : static_cast<T>(0);"
            )
            # Build p.q = q_orig + c*dt*qd_orig, p.qd = qd_orig + c*dt*qdd_{prev}.
            # Prior stage's qdd lives in s_stage_grad_qdd at offset (stage_idx - 1) * n.
            prev_offset = (stage_idx - 1) * n
            self.gen_add_parallel_loop("ind", str(n), use_thread_group)
            self.gen_add_code_line(f"s_q[ind] = s_q_orig[ind] + c_offset * dt * s_qd_orig[ind];")
            self.gen_add_code_line(f"s_qd[ind] = s_qd_orig[ind] + c_offset * dt * s_stage_grad_qdd[{prev_offset} + ind];")
            self.gen_add_end_control_flow()
            self.gen_add_sync(use_thread_group)
            # Update XImats for the new s_q.
            self.gen_load_update_XImats_helpers_function_call(use_thread_group)
            self.gen_add_sync(use_thread_group)

        # Run FD-gradient at this stage's (s_q, s_qd, s_u).
        # After this call: s_qdd, s_Minv, s_df_du = stage `stage_num` values.
        self.gen_forward_dynamics_gradient_inner_python(
            use_thread_group=use_thread_group,
            use_qdd_Minv_input=False,
            s_df_du_name="s_df_du",
        )
        self.gen_add_sync(use_thread_group)

        # Always save this stage's qdd into s_stage_grad_qdd[stage_idx * n].
        # The next stage's FD-grad will overwrite s_qdd, so we need this
        # snapshot to (a) build p_{stage+1} on the next iteration and
        # (b) assemble x_{k+1} at the end when compute_x_kp1.
        self.gen_add_parallel_loop("ind", str(n), use_thread_group)
        self.gen_add_code_line(
            f"s_stage_grad_qdd[{stage_idx * n} + ind] = s_qdd[ind];"
        )
        self.gen_add_end_control_flow()
        self.gen_add_sync(use_thread_group)

        # Compute D_qdd_{stage_num} into s_D_qdd_stage[stage_idx * (n*3n)].
        # For stage 1 we use the unified formula with c_0 = 0 (i.e. no chain).
        if stage_idx == 0:
            c_prev_expr = "static_cast<T>(0)"
        else:
            offset_branches = []
            for name, (cnt, c_list, _) in _INTEGRATOR_BUTCHER.items():
                if cnt >= stage_num:
                    offset_branches.append(
                        f"(IT == IntegratorType::{name}) ? static_cast<T>({c_list[stage_idx - 1]})"
                    )
            c_prev_expr = " : ".join(offset_branches) + " : static_cast<T>(0)"
        self.gen_add_code_line(f"constexpr T c_prev_s{stage_num} = {c_prev_expr};")
        self.gen_add_code_line(
            f"T *s_D_qdd_cur = &s_D_qdd_stage[{stage_idx} * {n * three_n}];"
        )
        # The chain-rule input is D_qdd_{i-1}; for stage 1 we don't need it.
        if stage_idx > 0:
            self.gen_add_code_line(
                f"T *s_D_qdd_prev = &s_D_qdd_stage[{(stage_idx - 1) * n * three_n}];"
            )
        self.gen_add_parallel_loop("ind", str(n * three_n), use_thread_group)
        self.gen_add_code_line(f"int r = ind % {n};")
        self.gen_add_code_line(f"int c = ind / {n};")
        # base term per block.
        self.gen_add_code_line("T base = static_cast<T>(0);")
        self.gen_add_code_line(f"if (c < {n}) {{")
        self.gen_add_code_line(f"    base = s_df_du[c * {n} + r];                       // J_qq[r, c]")
        self.gen_add_code_line(f"}} else if (c < {2 * n}) {{")
        self.gen_add_code_line(f"    int cc = c - {n};")
        self.gen_add_code_line(f"    base = c_prev_s{stage_num} * dt * s_df_du[cc * {n} + r] + s_df_du[{nn} + cc * {n} + r];  // c*dt*J_qq + J_qv")
        self.gen_add_code_line("} else {")
        self.gen_add_code_line(f"    int cc = c - {2 * n};")
        # s_Minv is SYMMETRIC_UPPER triangular n × n.
        self.gen_add_code_line(f"    int midx = (r <= cc) * (cc * {n} + r) + (r > cc) * (r * {n} + cc);")
        self.gen_add_code_line("    base = s_Minv[midx];                                  // Minv[r, cc]")
        self.gen_add_code_line("}")
        # Chain-rule contribution.
        if stage_idx > 0:
            self.gen_add_code_line("T chain = static_cast<T>(0);")
            self.gen_add_code_line(f"for (int k = 0; k < {n}; ++k) {{")
            # J_qv[r, k] = s_df_du[n*n + k*n + r]; D_qdd_prev[k, c] = s_D_qdd_prev[c*n + k]
            self.gen_add_code_line(
                f"    chain += s_df_du[{nn} + k * {n} + r] * s_D_qdd_prev[c * {n} + k];"
            )
            self.gen_add_code_line("}")
            self.gen_add_code_line(
                f"s_D_qdd_cur[c * {n} + r] = base + c_prev_s{stage_num} * dt * chain;"
            )
        else:
            self.gen_add_code_line(f"s_D_qdd_cur[c * {n} + r] = base;")
        self.gen_add_end_control_flow()
        self.gen_add_sync(use_thread_group)

        self.gen_add_end_control_flow()  # close `if constexpr (gating)`

    # ----- Final assembly of dAB and optional x_kp1 -----
    self.gen_add_code_line("// --- multi-stage gradient: assemble final dAB ---")
    self.gen_add_parallel_loop("ind", str(2 * n * three_n), use_thread_group)
    self.gen_add_code_line(f"int row = ind % {2 * n};")
    self.gen_add_code_line(f"int col = ind / {2 * n};")
    self.gen_add_code_line("T val = static_cast<T>(0);")
    self.gen_add_code_line(f"if (row < {n}) {{")
    # Top half: q_{k+1} = q + dt*qd (same as Euler).
    self.gen_add_code_line(f"    if (col < {n}) {{")
    self.gen_add_code_line("        val = (row == col) ? static_cast<T>(1) : static_cast<T>(0);")
    self.gen_add_code_line(f"    }} else if (col < {2 * n}) {{")
    self.gen_add_code_line(f"        val = (row == col - {n}) ? dt : static_cast<T>(0);")
    self.gen_add_code_line("    } else {")
    self.gen_add_code_line("        val = static_cast<T>(0);")
    self.gen_add_code_line("    }")
    self.gen_add_code_line("} else {")
    self.gen_add_code_line(f"    int r = row - {n};")
    self.gen_add_code_line("    // bottom half: qd_{k+1} = qd + dt * sum(b_i * qdd_i)")
    self.gen_add_code_line(f"    T identity_term = ((col >= {n} && col < {2 * n}) && (r == col - {n})) ? static_cast<T>(1) : static_cast<T>(0);")
    self.gen_add_code_line("    T accel_term = static_cast<T>(0);")
    for name, (cnt, _, b_list) in _INTEGRATOR_BUTCHER.items():
        self.gen_add_code_line(f"    if constexpr (IT == IntegratorType::{name}) {{")
        for i, b in enumerate(b_list):
            if b == 0:
                continue
            self.gen_add_code_line(
                f"        accel_term += static_cast<T>({b}) * s_D_qdd_stage[{i * n * three_n} + col * {n} + r];"
            )
        self.gen_add_code_line("    }")
    self.gen_add_code_line("    val = identity_term + dt * accel_term;")
    self.gen_add_code_line("}")
    self.gen_add_code_line("s_dAB[ind] = val;")
    self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)

    # Optionally also assemble x_{k+1}.
    if compute_x_kp1:
        self.gen_add_code_line("// --- multi-stage gradient: assemble x_{k+1} ---")
        self.gen_add_parallel_loop("ind", str(2 * n), use_thread_group)
        self.gen_add_code_line(f"if (ind < {n}) {{")
        self.gen_add_code_line(f"    s_x_kp1[ind] = s_q_orig[ind] + dt * s_qd_orig[ind];")
        self.gen_add_code_line("} else {")
        self.gen_add_code_line(f"    int j = ind - {n};")
        self.gen_add_code_line("    T accel = static_cast<T>(0);")
        for name, (cnt, _, b_list) in _INTEGRATOR_BUTCHER.items():
            self.gen_add_code_line(f"    if constexpr (IT == IntegratorType::{name}) {{")
            for i, b in enumerate(b_list):
                if b == 0:
                    continue
                # All stage qdd's now live in s_stage_grad_qdd[i*n..(i+1)*n].
                src = f"s_stage_grad_qdd[{i * n} + j]"
                self.gen_add_code_line(f"        accel += static_cast<T>({b}) * {src};")
            self.gen_add_code_line("    }")
        self.gen_add_code_line(f"    s_x_kp1[ind] = s_qd_orig[j] + dt * accel;")
        self.gen_add_code_line("}")
        self.gen_add_end_control_flow()
        self.gen_add_sync(use_thread_group)


def gen_integrator_gradient_inner_python(self, use_thread_group=False, compute_x_kp1=False,
                                          integrator_type="IT", s_dAB_name="s_dAB",
                                          s_x_kp1_name="s_x_kp1"):
    """Compose: FD gradient (sets s_Minv, s_qdd, s_dc_du, s_df_du) → dAB assembly.

    Mirrors the FD-gradient inner_python with default (non-spill) arguments —
    spill machinery for ID_DU is not exercised by the integrator path on
    typical fixed-base robots (iiwa14 / kuka). If a future large robot
    triggers the spill, copy the spill-aware variant here and pass through.
    """
    n = self.robot.get_num_vel()
    # Run FD gradient, which writes s_df_du (n*2n) and leaves s_Minv, s_qdd in shared.
    # use_qdd_Minv_input=False → we compute Minv/qdd inline (we always want both fresh).
    self.gen_forward_dynamics_gradient_inner_python(
        use_thread_group=use_thread_group,
        use_qdd_Minv_input=False,
        s_df_du_name="s_df_du",
    )
    self.gen_add_sync(use_thread_group)
    # Assemble dAB.
    self.gen_integrator_gradient_dAB_assembly(
        integrator_type=integrator_type,
        use_thread_group=use_thread_group,
        s_dAB_name=s_dAB_name,
    )
    # Optionally also build x_kp1 (the value) from the s_qdd that FD-gradient
    # populated. Saves a redundant FD pass for the both-at-once API.
    if compute_x_kp1:
        self.gen_add_sync(use_thread_group)
        self.gen_integrator_finish_function_call(
            integrator_type=integrator_type,
            use_thread_group=use_thread_group,
            updated_var_names=dict(s_x_kp1_name=s_x_kp1_name),
        )


def gen_integrator_gradient_device(self, use_thread_group=False, compute_x_kp1=False):
    n = self.robot.get_num_vel()
    suffix = "_with_x_kp1" if compute_x_kp1 else ""
    extra_out_params = (["s_x_kp1 is a pointer to memory for the next state (size 2*NUM_VEL)"]
                        if compute_x_kp1 else [])
    func_params = ["s_dAB is a pointer to memory for [A | B] of size 2*NUM_VEL*3*NUM_VEL (column-major)"] + extra_out_params + [
        "s_q is the vector of joint positions",
        "s_qd is the vector of joint velocities",
        "s_u is the vector of joint input torques",
        "d_robotModel is the pointer to the initialized model specific helpers on the GPU",
        "gravity is the gravity constant",
        "dt is the integration timestep",
    ]
    sig_x_kp1 = "T *s_x_kp1, " if compute_x_kp1 else ""
    func_def_start = "void integrator_gradient" + suffix + "_device(T *s_dAB, " + sig_x_kp1 + "const T *s_q, const T *s_qd, const T *s_u, "
    func_def_end = "const robotModel<T> *d_robotModel, const T gravity, const T dt) {"
    if use_thread_group:
        func_def_start = func_def_start.replace("(", "(cgrps::thread_group tgrp, ", 1)
        func_params.insert(0, "tgrp is the handle to the thread_group running this function")
    self.gen_add_func_doc("Computes the gradient of the integrator step" +
                          (" and the next state x_{k+1}" if compute_x_kp1 else ""),
                          [], func_params, None)
    self.gen_add_code_line("template <typename T, IntegratorType IT = IntegratorType::EULER>")
    self.gen_add_code_line("__device__")
    self.gen_add_code_line(func_def_start + func_def_end, True)
    # Allocate the per-call shared scratch (mirrors FD-gradient device).
    inner_temp_size = self.gen_integrator_gradient_inner_temp_mem_size()
    extra_t_buffers = [("s_vaf", 18 * n), ("s_dc_du", n * 2 * n), ("s_df_du", n * 2 * n),
                       ("s_Minv", n * n), ("s_qdd", n)]
    self.gen_XImats_helpers_temp_shared_memory_code(
        inner_temp_size, extra_t_buffers=extra_t_buffers, include_linalg_scratch=True,
    )
    self.gen_load_update_XImats_helpers_function_call(use_thread_group)
    self.gen_integrator_gradient_inner_python(
        use_thread_group=use_thread_group,
        compute_x_kp1=compute_x_kp1,
        integrator_type="IT",
        s_dAB_name="s_dAB",
        s_x_kp1_name="s_x_kp1",
    )
    self.gen_add_end_function()


def gen_integrator_gradient_kernel(self, use_thread_group=False, compute_x_kp1=False, single_call_timing=False):
    n = self.robot.get_num_vel()
    suffix = "_with_x_kp1" if compute_x_kp1 else ""
    func_params = ["d_dAB is a pointer to memory for [A | B] of size 2*NUM_VEL*3*NUM_VEL per timestep"]
    if compute_x_kp1:
        func_params.append("d_x_kp1 is a pointer to memory for the next state (size 2*NUM_VEL per timestep)")
    func_params += [
        "d_q_qd_u is the packed joint positions, velocities, and input torques",
        "stride_q_qd_u is the stride between each (q, qd, u) tuple in d_q_qd_u",
        "d_robotModel is the pointer to the initialized model specific helpers on the GPU",
        "gravity is the gravity constant",
        "dt is the integration timestep",
        "num_timesteps is the length of the trajectory (or overloaded as test_iters for timing)",
    ]
    sig_x_kp1 = "T *d_x_kp1, " if compute_x_kp1 else ""
    func_def_start = "void integrator_gradient" + suffix + "_kernel(T *d_dAB, " + sig_x_kp1 + "const T *d_q_qd_u, const int stride_q_qd_u, "
    func_def_end = "const robotModel<T> *d_robotModel, const T gravity, const T dt, const int NUM_TIMESTEPS) {"
    func_def = func_def_start + func_def_end
    if single_call_timing:
        func_def = func_def.replace("kernel(", "kernel_single_timing(")
    self.gen_add_func_doc("Computes the gradient of the integrator step per timestep" +
                          (" and the next state x_{k+1}" if compute_x_kp1 else ""),
                          [], func_params, None)
    self.gen_add_code_line("template <typename T, IntegratorType IT = IntegratorType::EULER, int RESOURCE_TIER = TIER_PERF>")
    self.gen_add_code_line("__global__")
    self.gen_add_code_line("__launch_bounds__(tier_max_threads<RESOURCE_TIER>())")
    self.gen_add_code_line(func_def, True)
    inner_temp_size = self.gen_integrator_gradient_inner_temp_mem_size()
    fb = self.robot.floating_base
    max_stages = _max_stages_in_use()
    extra_t_buffers = [
        ("s_q_qd_u", 3 * n + fb),
        ("s_dAB", 2 * n * 3 * n),
        ("s_df_du", n * 2 * n),
        ("s_dc_du", n * 2 * n),
        ("s_vaf", 18 * n),
        ("s_Minv", n * n),
        ("s_qdd", n),
        # Multi-stage scratch — allocated for every IT (single-stage just doesn't use it).
        ("s_q_orig", n),
        ("s_qd_orig", n),
        ("s_stage_grad_qdd", max_stages * n),
        ("s_D_qdd_stage", max_stages * n * 3 * n),
    ]
    if compute_x_kp1:
        extra_t_buffers.append(("s_x_kp1", 2 * n + fb))  # = nq + nv
    self.gen_XImats_helpers_temp_shared_memory_code(
        inner_temp_size, extra_t_buffers=extra_t_buffers, include_linalg_scratch=True,
    )
    self.gen_add_code_line(
        "T *s_q = s_q_qd_u; T *s_qd = &s_q_qd_u[" + str(n + fb) + "]; T *s_u = &s_q_qd_u[" + str(2 * n + fb) + "];"
    )
    if use_thread_group:
        self.gen_add_code_line("cgrps::thread_group tgrp = TBD;")
    def _emit_body():
        # Dispatch single-stage vs multi-stage at compile time on IT.
        self.gen_add_code_line(
            "if constexpr (IT == IntegratorType::EULER || IT == IntegratorType::SEMI_IMPLICIT_EULER) {", True
        )
        self.gen_integrator_gradient_inner_python(
            use_thread_group=use_thread_group,
            compute_x_kp1=compute_x_kp1,
            integrator_type="IT",
            s_dAB_name="s_dAB",
            s_x_kp1_name="s_x_kp1",
        )
        self.gen_add_end_control_flow()
        self.gen_add_code_line("else {", True)
        self.gen_integrator_gradient_multistage(
            use_thread_group=use_thread_group,
            compute_x_kp1=compute_x_kp1,
        )
        self.gen_add_end_control_flow()

    if not single_call_timing:
        self.gen_add_parallel_loop("k", "NUM_TIMESTEPS", use_thread_group, block_level=True)
        self.gen_kernel_load_inputs("q_qd_u", "stride_q_qd_u", str(3 * n + fb), use_thread_group)
        self.gen_add_code_line("// compute")
        self.gen_load_update_XImats_helpers_function_call(use_thread_group)
        _emit_body()
        self.gen_add_sync(use_thread_group)
        self.gen_kernel_save_result("dAB", str(2 * n * 3 * n), str(2 * n * 3 * n), use_thread_group)
        if compute_x_kp1:
            self.gen_kernel_save_result("x_kp1", str(2 * n + fb), str(2 * n + fb), use_thread_group)
        self.gen_add_end_control_flow()
    else:
        input_count = 3 * n + fb
        self.gen_kernel_load_inputs_single_timing("q_qd_u", str(input_count))
        self.gen_add_code_line("// compute with NUM_TIMESTEPS as NUM_REPS for timing")
        self.gen_add_code_line("for (int rep = 0; rep < NUM_TIMESTEPS; rep++){", True)
        self.gen_anti_licm_input_reload("q_qd_u", str(input_count), use_thread_group, feedback_from="dAB")
        self.gen_load_update_XImats_helpers_function_call(use_thread_group)
        _emit_body()
        self.gen_anti_licm_output_write("dAB")
        self.gen_add_end_control_flow()
        self.gen_kernel_save_result_single_timing("dAB", str(2 * n * 3 * n), use_thread_group)
        if compute_x_kp1:
            self.gen_kernel_save_result_single_timing("x_kp1", str(2 * n + fb), use_thread_group)
    self.gen_add_end_function()


def gen_integrator_gradient_host(self, mode=0, compute_x_kp1=False):
    single_call_timing = mode == 1
    compute_only = mode == 2
    suffix = "_with_x_kp1" if compute_x_kp1 else ""
    base_name = "integrator_gradient" + suffix
    func_params = ["hd_data is the packaged input and output pointers",
                   "d_robotModel is the pointer to the initialized model specific helpers on the GPU",
                   "gravity is the gravity constant",
                   "dt is the integration timestep",
                   "num_timesteps is the length of the trajectory (or overloaded as test_iters for timing)",
                   "streams are pointers to CUDA streams for async memory transfers (if needed)"]
    func_def_start = "void " + base_name + "(gridData<T, KIND> *hd_data, const robotModel<T> *d_robotModel, const T gravity, const T dt, const int num_timesteps,"
    func_def_end = "                  const dim3 block_dimms, const dim3 thread_dimms, cudaStream_t *streams) {"
    if single_call_timing:
        func_def_start = func_def_start.replace("(", "_single_timing(", 1)
        func_def_end = "              " + func_def_end
    if compute_only:
        func_def_start = func_def_start.replace("(", "_compute_only(", 1)
        func_def_end = "             " + func_def_end.replace(", cudaStream_t *streams", "")
    self.gen_add_func_doc(
        "Run the integrator gradient (default Euler) per timestep" +
        (" and also write x_{k+1}" if compute_x_kp1 else ""),
        [], func_params, None,
    )
    self.gen_add_code_line("template <typename T, IntegratorType IT = IntegratorType::EULER, gridDataKind KIND = GRID_DATA_ALL>")
    self.gen_add_code_line("__host__")
    self.gen_add_code_line(func_def_start)
    self.gen_add_code_line(func_def_end, True)
    self.gen_add_code_line("static_assert(KIND == GRID_DATA_ALL || KIND == GRID_DATA_DYNAMICS, \"" + base_name + " requires all-data or dynamics gridData\");")
    kernel_args_x_kp1 = "hd_data->d_x_kp1," if compute_x_kp1 else ""
    func_call_start = (base_name + "_kernel<T, IT><<<block_dimms,thread_dimms,INTEGRATOR_DU_DYNAMIC_SHARED_MEM_BYTES<T>()>>>(" +
                       "hd_data->d_dAB," + kernel_args_x_kp1 + "hd_data->d_q_qd_u,stride_q_qd_u,")
    func_call_end = "d_robotModel,gravity,dt,num_timesteps);"
    if single_call_timing:
        func_call_start = func_call_start.replace("kernel<T, IT>", "kernel_single_timing<T, IT>")
    self.gen_add_code_line("int stride_q_qd_u = 3*NUM_JOINTS;")
    if not compute_only:
        self.gen_add_code_lines([
            "// start code with memory transfer",
            "gpuErrchk(cudaMemcpyAsync(hd_data->d_q_qd_u,hd_data->h_q_qd_u,stride_q_qd_u*" +
                ("num_timesteps*" if not single_call_timing else "") + "sizeof(T),cudaMemcpyHostToDevice,streams[0]));",
            "gpuErrchkKernel();",
        ])
    self.gen_add_code_line("// then call the kernel")
    func_call_code = [func_call_start + func_call_end, "gpuErrchkKernel();"]
    if single_call_timing:
        func_call_code.insert(0, "struct timespec start, end; clock_gettime(CLOCK_MONOTONIC,&start);")
        func_call_code.append("clock_gettime(CLOCK_MONOTONIC,&end);")
    self.gen_add_code_line("gpuErrchk(grid_check_dynamic_shared_memory_bytes(\"" + base_name + "\", INTEGRATOR_DU_DYNAMIC_SHARED_MEM_BYTES<T>()));")
    self.gen_add_code_lines(func_call_code)
    if not compute_only:
        self.gen_add_code_lines([
            "// finally transfer the result back",
            "gpuErrchk(cudaMemcpy(hd_data->h_dAB,hd_data->d_dAB,2*NUM_JOINTS*3*NUM_JOINTS*" +
                ("num_timesteps*" if not single_call_timing else "") + "sizeof(T),cudaMemcpyDeviceToHost));",
            "gpuErrchkKernel();",
        ])
        if compute_x_kp1:
            self.gen_add_code_lines([
                "gpuErrchk(cudaMemcpy(hd_data->h_x_kp1,hd_data->d_x_kp1,(NUM_POS + NUM_VEL)*" +
                    ("num_timesteps*" if not single_call_timing else "") + "sizeof(T),cudaMemcpyDeviceToHost));",
                "gpuErrchkKernel();",
            ])
    if single_call_timing:
        from ..algo_registry import single_call_printf_line
        # both variants reuse the same printf key for now
        self.gen_add_code_line(single_call_printf_line("integrator_gradient" if not compute_x_kp1 else "integrator_with_gradient"))
    self.gen_add_end_function()


def gen_integrator_gradient(self, use_thread_group=False):
    # Gradient-only device + kernel + host
    self.gen_integrator_gradient_device(use_thread_group, compute_x_kp1=False)
    self.gen_integrator_gradient_kernel(use_thread_group, compute_x_kp1=False, single_call_timing=True)
    self.gen_integrator_gradient_kernel(use_thread_group, compute_x_kp1=False, single_call_timing=False)
    self.gen_integrator_gradient_host(0, compute_x_kp1=False)
    self.gen_integrator_gradient_host(1, compute_x_kp1=False)
    self.gen_integrator_gradient_host(2, compute_x_kp1=False)
    # Gradient + x_kp1 (both-at-once) device + kernel + host
    self.gen_integrator_gradient_device(use_thread_group, compute_x_kp1=True)
    self.gen_integrator_gradient_kernel(use_thread_group, compute_x_kp1=True, single_call_timing=True)
    self.gen_integrator_gradient_kernel(use_thread_group, compute_x_kp1=True, single_call_timing=False)
    self.gen_integrator_gradient_host(0, compute_x_kp1=True)
    self.gen_integrator_gradient_host(1, compute_x_kp1=True)
    self.gen_integrator_gradient_host(2, compute_x_kp1=True)
