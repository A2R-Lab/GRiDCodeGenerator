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
    fb = self.robot.floating_base
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
    if fb:
        # Floating-base top-nv rows: d(q_kp1)/dq = dIntegrate_q(q, dt*qd).
        # SE(3) Adjoint is block-diagonal: top-left 6x6 from s_dInt_q_6x6,
        # identity for the revolute joint block.
        self.gen_add_code_line("        if (row < 6 && col < 6) {")
        self.gen_add_code_line("            val = s_dInt_q_6x6[row * 6 + col];")
        self.gen_add_code_line("        } else if (row >= 6 && col >= 6) {")
        self.gen_add_code_line("            val = (row == col) ? static_cast<T>(1) : static_cast<T>(0);")
        self.gen_add_code_line("        } else { val = static_cast<T>(0); }")
    else:
        self.gen_add_code_line("        val = (row == col) ? static_cast<T>(1) : static_cast<T>(0);")
    self.gen_add_code_line("    } else {")
    self.gen_add_code_line("        int i_local = row - " + str(n) + ";")
    self.gen_add_code_line("        val = dt * " + s_df_du_name + "[col * " + str(n) + " + i_local];")
    self.gen_add_code_line("    }")
    self.gen_add_code_line("} else if (col < " + str(2 * n) + ") {")
    self.gen_add_code_line("    // d/dqd column")
    self.gen_add_code_line("    int j_local = col - " + str(n) + ";")
    self.gen_add_code_line("    if (row < " + str(n) + ") {")
    if fb:
        # Top-nv rows of d/dqd: dt * dIntegrate_v(q, dt*qd).
        self.gen_add_code_line("        if (row < 6 && j_local < 6) {")
        self.gen_add_code_line("            val = dt * s_dInt_v_6x6[row * 6 + j_local];")
        self.gen_add_code_line("        } else if (row >= 6 && j_local >= 6) {")
        self.gen_add_code_line("            val = (row == j_local) ? dt : static_cast<T>(0);")
        self.gen_add_code_line("        } else { val = static_cast<T>(0); }")
    else:
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
    if fb:
        # Floating-base SI-Euler: v_new = qd + dt*qdd; q_new = integrate(q, dt*v_new).
        #   bottom rows  = [dvdq | dvdv | dvdu] = [dt*J_qq | I + dt*J_qv | dt*Minv]
        #   top rows     = [dInt_q + dt*dInt_v@dvdq | dt*dInt_v@dvdv | dt*dInt_v@dvdu]
        # dInt_q / dInt_v are evaluated at v_dt = dt*v_new (precomputed into the
        # 6x6 blocks by gen_integrator_gradient_inner_python). dInt is block-diag:
        # the 6x6 free-flyer block plus identity on the revolute joints, so the
        # dInt_v @ dvdX matmul only mixes rows < 6.
        self.gen_add_code_line("if (col < " + str(n) + ") {")
        self.gen_add_code_line("    // d/dq column. dvdq[k,c] = dt*J_qq[k,c] = dt*s_df_du[c*n + k].")
        self.gen_add_code_line("    int c = col;")
        self.gen_add_code_line("    if (row < " + str(n) + ") {")
        self.gen_add_code_line("        T dInt_q_term = (row < 6 && c < 6) ? s_dInt_q_6x6[row * 6 + c]")
        self.gen_add_code_line("                       : ((row >= 6 && c >= 6) ? ((row == c) ? static_cast<T>(1) : static_cast<T>(0)) : static_cast<T>(0));")
        self.gen_add_code_line("        T mm;")
        self.gen_add_code_line("        if (row < 6) { mm = static_cast<T>(0); for (int k = 0; k < 6; ++k) mm += s_dInt_v_6x6[row * 6 + k] * (dt * " + s_df_du_name + "[c * " + str(n) + " + k]); }")
        self.gen_add_code_line("        else { mm = dt * " + s_df_du_name + "[c * " + str(n) + " + row]; }")
        self.gen_add_code_line("        val = dInt_q_term + dt * mm;")
        self.gen_add_code_line("    } else {")
        self.gen_add_code_line("        int i_local = row - " + str(n) + ";")
        self.gen_add_code_line("        val = dt * " + s_df_du_name + "[c * " + str(n) + " + i_local];")
        self.gen_add_code_line("    }")
        self.gen_add_code_line("} else if (col < " + str(2 * n) + ") {")
        self.gen_add_code_line("    // d/dqd column. dvdv[k,c] = (k==c) + dt*J_qv[k,c].")
        self.gen_add_code_line("    int c = col - " + str(n) + ";")
        self.gen_add_code_line("    if (row < " + str(n) + ") {")
        self.gen_add_code_line("        T mm;")
        self.gen_add_code_line("        if (row < 6) { mm = static_cast<T>(0); for (int k = 0; k < 6; ++k) { T dvdv_k = ((k == c) ? static_cast<T>(1) : static_cast<T>(0)) + dt * " + s_df_du_name + "[" + str(nn) + " + c * " + str(n) + " + k]; mm += s_dInt_v_6x6[row * 6 + k] * dvdv_k; } }")
        self.gen_add_code_line("        else { mm = ((row == c) ? static_cast<T>(1) : static_cast<T>(0)) + dt * " + s_df_du_name + "[" + str(nn) + " + c * " + str(n) + " + row]; }")
        self.gen_add_code_line("        val = dt * mm;")
        self.gen_add_code_line("    } else {")
        self.gen_add_code_line("        int i_local = row - " + str(n) + ";")
        self.gen_add_code_line("        T diag = (i_local == c) ? static_cast<T>(1) : static_cast<T>(0);")
        self.gen_add_code_line("        val = diag + dt * " + s_df_du_name + "[" + str(nn) + " + c * " + str(n) + " + i_local];")
        self.gen_add_code_line("    }")
        self.gen_add_code_line("} else {")
        self.gen_add_code_line("    // d/du column. dvdu[k,c] = dt*Minv[k,c] (SYMMETRIC_UPPER).")
        self.gen_add_code_line("    int c = col - " + str(2 * n) + ";")
        self.gen_add_code_line("    if (row < " + str(n) + ") {")
        self.gen_add_code_line("        T mm;")
        self.gen_add_code_line("        if (row < 6) { mm = static_cast<T>(0); for (int k = 0; k < 6; ++k) { int midx = (k <= c) * (c * " + str(n) + " + k) + (k > c) * (k * " + str(n) + " + c); mm += s_dInt_v_6x6[row * 6 + k] * (dt * " + s_Minv_name + "[midx]); } }")
        self.gen_add_code_line("        else { int midx = (row <= c) * (c * " + str(n) + " + row) + (row > c) * (row * " + str(n) + " + c); mm = dt * " + s_Minv_name + "[midx]; }")
        self.gen_add_code_line("        val = dt * mm;")
        self.gen_add_code_line("    } else {")
        self.gen_add_code_line("        int i_local = row - " + str(n) + ";")
        self.gen_add_code_line("        int midx = (i_local <= c) * (c * " + str(n) + " + i_local) + (i_local > c) * (i_local * " + str(n) + " + c);")
        self.gen_add_code_line("        val = dt * " + s_Minv_name + "[midx];")
        self.gen_add_code_line("    }")
        self.gen_add_code_line("}")
    else:
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
    self.gen_add_code_line("              \"dAB assembly handles single-stage IT only; Midpoint/RK3/RK4 are routed through gen_integrator_gradient_multistage.\");")
    self.gen_add_end_control_flow()
    self.gen_add_end_control_flow()  # end parallel loop


def gen_integrator_gradient_multistage(self, use_thread_group=False, compute_x_kp1=False,
                                       d_temp_spill_name="nullptr", temp_spill_flag_name="false"):
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
    nq = self.robot.get_num_pos()
    fb = self.robot.floating_base
    max_stages = _max_stages_in_use()
    three_n = 3 * n
    nn = n * n

    # Save original q (nq entries — floating-base carries the 7-element pose
    # prefix), qd (n) so we can rebuild p_{i+1} on later stages.
    self.gen_add_code_line("// --- multi-stage gradient: save original q, qd ---")
    self.gen_add_parallel_loop("ind", str(nq), use_thread_group)
    self.gen_add_code_line("s_q_orig[ind] = s_q[ind];")
    self.gen_add_end_control_flow()
    self.gen_add_parallel_loop("ind", str(n), use_thread_group)
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
            # Build p.q = integrate(q_orig, c*dt*qd_orig), p.qd = qd_orig + c*dt*qdd_{prev}.
            # Prior stage's qdd lives in s_stage_grad_qdd at offset (stage_idx - 1) * n.
            prev_offset = (stage_idx - 1) * n
            if fb:
                # Floating-base q-update is the SE(3) Lie retract over the full
                # nq layout; the velocity update is the plain Euler add.
                self.gen_add_serial_ops(use_thread_group)
                self.gen_add_code_line(f"T v_scaled[{n}];")
                self.gen_add_code_line(f"for (int i = 0; i < {n}; ++i) v_scaled[i] = c_offset * dt * s_qd_orig[i];")
                self.gen_add_code_line(f"grid_integrate_floating_q<T, {nq}>(s_q_orig, v_scaled, s_q);")
                self.gen_add_end_control_flow()
                self.gen_add_parallel_loop("ind", str(n), use_thread_group)
                self.gen_add_code_line(f"s_qd[ind] = s_qd_orig[ind] + c_offset * dt * s_stage_grad_qdd[{prev_offset} + ind];")
                self.gen_add_end_control_flow()
            else:
                self.gen_add_parallel_loop("ind", str(n), use_thread_group)
                self.gen_add_code_line(f"s_q[ind] = s_q_orig[ind] + c_offset * dt * s_qd_orig[ind];")
                self.gen_add_code_line(f"s_qd[ind] = s_qd_orig[ind] + c_offset * dt * s_stage_grad_qdd[{prev_offset} + ind];")
                self.gen_add_end_control_flow()
            self.gen_add_sync(use_thread_group)
            # Update XImats for the new s_q.
            self.gen_load_update_XImats_helpers_function_call(use_thread_group)
            self.gen_add_sync(use_thread_group)
            if fb:
                # Per-stage SE(3) dIntegrate blocks at v_dt = c*dt*qd_orig (the
                # q-perturbation increment for p.q = integrate(q_orig, c*dt*qd_orig)).
                # Reused buffers s_dInt_*_6x6 — consumed in this stage's D_qdd loop
                # below before the next stage overwrites them.
                self.gen_add_serial_ops(use_thread_group)
                self.gen_add_code_line(f"T v_dt_stage[{n}];")
                self.gen_add_code_line(f"for (int i = 0; i < {n}; ++i) v_dt_stage[i] = c_offset * dt * s_qd_orig[i];")
                self.gen_add_code_line("grid_dIntegrate_q_block<T>(v_dt_stage, s_dInt_q_6x6);")
                self.gen_add_code_line("grid_dIntegrate_v_block<T>(v_dt_stage, s_dInt_v_6x6);")
                self.gen_add_end_control_flow()
                self.gen_add_sync(use_thread_group)

        # Run FD-gradient at this stage's (s_q, s_qd, s_u).
        # After this call: s_qdd, s_Minv, s_df_du = stage `stage_num` values.
        self.gen_forward_dynamics_gradient_inner_python(
            use_thread_group=use_thread_group,
            use_qdd_Minv_input=False,
            s_df_du_name="s_df_du",
            d_temp_spill_name=d_temp_spill_name,
            temp_spill_flag_name=temp_spill_flag_name,
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
        # base term per block. For floating-base on stage > 1 the J_qq columns
        # are projected through the SE(3) dIntegrate blocks (block-diagonal: 6x6
        # free-flyer block + identity joints), mirroring J_qq_i @ dInt in Python.
        project = fb and stage_idx > 0
        self.gen_add_code_line("T base = static_cast<T>(0);")
        self.gen_add_code_line(f"if (c < {n}) {{")
        if project:
            self.gen_add_code_line("    // (J_qq @ dInt_q)[r, c]")
            self.gen_add_code_line("    if (c < 6) {")
            self.gen_add_code_line(f"        for (int k = 0; k < 6; ++k) base += s_df_du[k * {n} + r] * s_dInt_q_6x6[k * 6 + c];")
            self.gen_add_code_line("    } else {")
            self.gen_add_code_line(f"        base = s_df_du[c * {n} + r];")
            self.gen_add_code_line("    }")
        else:
            self.gen_add_code_line(f"    base = s_df_du[c * {n} + r];                       // J_qq[r, c]")
        self.gen_add_code_line(f"}} else if (c < {2 * n}) {{")
        self.gen_add_code_line(f"    int cc = c - {n};")
        if project:
            self.gen_add_code_line("    // c*dt*(J_qq @ dInt_v)[r, cc] + J_qv[r, cc]")
            self.gen_add_code_line("    T proj = static_cast<T>(0);")
            self.gen_add_code_line("    if (cc < 6) {")
            self.gen_add_code_line(f"        for (int k = 0; k < 6; ++k) proj += s_df_du[k * {n} + r] * s_dInt_v_6x6[k * 6 + cc];")
            self.gen_add_code_line("    } else {")
            self.gen_add_code_line(f"        proj = s_df_du[cc * {n} + r];")
            self.gen_add_code_line("    }")
            self.gen_add_code_line(f"    base = c_prev_s{stage_num} * dt * proj + s_df_du[{nn} + cc * {n} + r];")
        else:
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
    # q_{k+1} = integrate(q, dt*qd) (Euler-style for every RK variant), so the
    # top rows are [dInt_q | dt*dInt_v | 0] at v_dt = dt*qd. For fixed-base these
    # reduce to [I | dt*I | 0].
    if fb:
        self.gen_add_serial_ops(use_thread_group)
        self.gen_add_code_line(f"T v_dt_final[{n}];")
        self.gen_add_code_line(f"for (int i = 0; i < {n}; ++i) v_dt_final[i] = dt * s_qd_orig[i];")
        self.gen_add_code_line("grid_dIntegrate_q_block<T>(v_dt_final, s_dInt_q_6x6);")
        self.gen_add_code_line("grid_dIntegrate_v_block<T>(v_dt_final, s_dInt_v_6x6);")
        self.gen_add_end_control_flow()
        self.gen_add_sync(use_thread_group)
    self.gen_add_code_line("// --- multi-stage gradient: assemble final dAB ---")
    self.gen_add_parallel_loop("ind", str(2 * n * three_n), use_thread_group)
    self.gen_add_code_line(f"int row = ind % {2 * n};")
    self.gen_add_code_line(f"int col = ind / {2 * n};")
    self.gen_add_code_line("T val = static_cast<T>(0);")
    self.gen_add_code_line(f"if (row < {n}) {{")
    if fb:
        self.gen_add_code_line(f"    if (col < {n}) {{")
        self.gen_add_code_line("        if (row < 6 && col < 6) { val = s_dInt_q_6x6[row * 6 + col]; }")
        self.gen_add_code_line("        else if (row >= 6 && col >= 6) { val = (row == col) ? static_cast<T>(1) : static_cast<T>(0); }")
        self.gen_add_code_line("        else { val = static_cast<T>(0); }")
        self.gen_add_code_line(f"    }} else if (col < {2 * n}) {{")
        self.gen_add_code_line(f"        int j = col - {n};")
        self.gen_add_code_line("        if (row < 6 && j < 6) { val = dt * s_dInt_v_6x6[row * 6 + j]; }")
        self.gen_add_code_line("        else if (row >= 6 && j >= 6) { val = (row == j) ? dt : static_cast<T>(0); }")
        self.gen_add_code_line("        else { val = static_cast<T>(0); }")
        self.gen_add_code_line("    } else {")
        self.gen_add_code_line("        val = static_cast<T>(0);")
        self.gen_add_code_line("    }")
    else:
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
        # v_{k+1} = qd + dt * sum(b_i * qdd_i); stage qdds live in s_stage_grad_qdd.
        self.gen_add_parallel_loop("ind", str(n), use_thread_group)
        self.gen_add_code_line("T accel = static_cast<T>(0);")
        for name, (cnt, _, b_list) in _INTEGRATOR_BUTCHER.items():
            self.gen_add_code_line(f"if constexpr (IT == IntegratorType::{name}) {{")
            for i, b in enumerate(b_list):
                if b == 0:
                    continue
                self.gen_add_code_line(f"    accel += static_cast<T>({b}) * s_stage_grad_qdd[{i * n} + ind];")
            self.gen_add_code_line("}")
        v_out_index = f"{nq} + ind" if fb else f"{n} + ind"
        self.gen_add_code_line(f"s_x_kp1[{v_out_index}] = s_qd_orig[ind] + dt * accel;")
        self.gen_add_end_control_flow()
        self.gen_add_sync(use_thread_group)
        # q_{k+1} = integrate(q, dt*qd) (Euler-style q-update for every RK variant).
        if fb:
            self.gen_add_serial_ops(use_thread_group)
            self.gen_add_code_line(f"T v_scaled_x[{n}];")
            self.gen_add_code_line(f"for (int i = 0; i < {n}; ++i) v_scaled_x[i] = dt * s_qd_orig[i];")
            self.gen_add_code_line(f"grid_integrate_floating_q<T, {nq}>(s_q_orig, v_scaled_x, s_x_kp1);")
            self.gen_add_end_control_flow()
        else:
            self.gen_add_parallel_loop("ind", str(n), use_thread_group)
            self.gen_add_code_line("s_x_kp1[ind] = s_q_orig[ind] + dt * s_qd_orig[ind];")
            self.gen_add_end_control_flow()
        self.gen_add_sync(use_thread_group)


def gen_integrator_gradient_inner_python(self, use_thread_group=False, compute_x_kp1=False,
                                          integrator_type="IT", s_dAB_name="s_dAB",
                                          s_x_kp1_name="s_x_kp1",
                                          d_temp_spill_name="nullptr", temp_spill_flag_name="false"):
    """Compose: FD gradient (sets s_Minv, s_qdd, s_dc_du, s_df_du) → dAB assembly.

    This is the single-stage path (Euler / SI-Euler). For floating-base it also
    computes the 6x6 SE(3) dIntegrate blocks (s_dInt_q_6x6, s_dInt_v_6x6) at the
    q-update increment — dt*qd for Euler, dt*v_new for SI-Euler — and reads them
    in the dAB top-nv rows. The SI-Euler floating top rows additionally fold in
    the dInt_v @ dv/dX matmul (see the SEMI_IMPLICIT_EULER branch in
    gen_integrator_gradient_dAB_assembly). Multi-stage (Midpoint/RK3/RK4) goes
    through gen_integrator_gradient_multistage instead.
    """
    fb = self.robot.floating_base
    n = self.robot.get_num_vel()
    self.gen_forward_dynamics_gradient_inner_python(
        use_thread_group=use_thread_group,
        use_qdd_Minv_input=False,
        s_df_du_name="s_df_du",
        d_temp_spill_name=d_temp_spill_name,
        temp_spill_flag_name=temp_spill_flag_name,
    )
    self.gen_add_sync(use_thread_group)
    if fb:
        # Precompute the SE(3) dIntegrate blocks at the q-update increment v_dt.
        # Euler:    q_new = integrate(q, dt*qd)            -> v_dt = dt*qd
        # SI-Euler: q_new = integrate(q, dt*v_new), where  -> v_dt = dt*(qd + dt*qdd)
        #           v_new = qd + dt*qdd  (s_qdd holds qdd after the FD gradient).
        self.gen_add_serial_ops(use_thread_group)
        self.gen_add_code_line(f"T v_dt_for_dInt[{n}];")
        self.gen_add_code_line("if constexpr (" + _integrator_type_token(integrator_type) + " == IntegratorType::SEMI_IMPLICIT_EULER) {")
        self.gen_add_code_line(f"    for (int i = 0; i < {n}; ++i) v_dt_for_dInt[i] = dt * (s_qd[i] + dt * s_qdd[i]);")
        self.gen_add_code_line("} else {")
        self.gen_add_code_line(f"    for (int i = 0; i < {n}; ++i) v_dt_for_dInt[i] = dt * s_qd[i];")
        self.gen_add_code_line("}")
        self.gen_add_code_line("grid_dIntegrate_q_block<T>(v_dt_for_dInt, s_dInt_q_6x6);")
        self.gen_add_code_line("grid_dIntegrate_v_block<T>(v_dt_for_dInt, s_dInt_v_6x6);")
        self.gen_add_end_control_flow()
        self.gen_add_sync(use_thread_group)
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


def gen_integrator_gradient_full_inner_function_call(self, use_thread_group=False, compute_x_kp1=False,
                                                     scratch_in_smem_expr="true",
                                                     use_da_df_spill_expr="false",
                                                     d_workspace_pool_name="nullptr",
                                                     d_temp_spill_name="nullptr"):
    """Emit the call to `integrator_gradient[_with_x_kp1]_full_inner`. Arg order MUST
    match the def in gen_integrator_gradient_full_inner. The FD-grad inner POOL
    placement region (d_workspace) and the id_du da_df band spill region
    (d_temp_spill) default to nullptr (unused under the matching if-constexpr); the
    kernel passes real pointers per tier. s_D_qdd_stage / s_dAB remain SEPARATE
    caller-placed pointers — they are threaded through unchanged."""
    suffix = "_with_x_kp1" if compute_x_kp1 else ""
    fname = "integrator_gradient" + suffix + "_full_inner"
    tmpl = "<T, IT, " + scratch_in_smem_expr + ", " + use_da_df_spill_expr + ">"
    start = fname + tmpl + "(s_dAB, "
    if compute_x_kp1:
        start += "s_x_kp1, "
    start += ("s_q, s_qd, s_u, s_df_du, s_dc_du, s_vaf, s_Minv, s_qdd, "
              "s_q_orig, s_qd_orig, s_stage_grad_qdd, s_D_qdd_stage, "
              "s_dInt_q_6x6, s_dInt_v_6x6, ")
    middle = self.gen_insert_helpers_function_call()
    end = ("s_temp, " + d_workspace_pool_name + ", " + d_temp_spill_name + ", "
           + "d_robotModel, gravity, dt);")
    if use_thread_group:
        start = start.replace("(", "(tgrp, ", 1)
    self.gen_add_code_line(start + middle + end)


def gen_integrator_gradient_full_inner(self, use_thread_group=False, compute_x_kp1=False):
    """Emit `integrator_gradient[_with_x_kp1]_full_inner` — the whole integrator
    gradient orchestration as ONE inner that OWNS its FD-grad scratch (s_temp) pool
    placement (inner-owns-placement; mirrors gen_inverse_dynamics_gradient_full_inner /
    gen_fdsva_so_full_inner). It wraps, in order:
      [repoint s_temp] -> load_update_XImats -> (compile-time IT dispatch)
        single-stage: gen_integrator_gradient_inner_python (Euler / SI-Euler)
        multi-stage : gen_integrator_gradient_multistage    (Midpoint / RK3 / RK4)
    Because the s_temp repoint happens at the very top, EVERY consumer below —
    including the XImats helper's sincos scratch and the per-stage XImats refresh in
    the multi-stage path — follows the placement, so the kernel never repoints
    s_temp from the outside.

    TWO independent template flags:
      SCRATCH_IN_SMEM : the shared FD-grad inner s_temp pool lives in smem (true) or
                        routes the WHOLE pool to d_workspace (false; the rung-2
                        whole-pool global-temp path, formerly the kernel's line-744
                        repoint). Dominant lever on big floating humanoids.
      USE_DA_DF_SPILL : the FD-grad inner's id_du band selectively spills its
                        da_dq..fxvi band to d_temp_spill (rung 1). Threaded through
                        the stable gen_forward_dynamics_gradient_inner_python
                        composition surface to the id_du band sub-inner.

    Pointer params are caller-supplied (the kernel decides where the OUTPUT s_dAB
    and the multi-band scratch s_D_qdd_stage live — smem or de-aliased workspace
    sub-offsets — and hands in the spill regions): only the FD-grad inner s_temp
    POOL placement is the inner's call. s_q / s_qd are NON-const: the multi-stage
    path MUTATES them in shared memory across RK stages (and re-derives the per-stage
    XImats from the freshly-mutated s_q)."""
    n = self.robot.get_num_vel()
    fb = self.robot.floating_base
    suffix = "_with_x_kp1" if compute_x_kp1 else ""
    fname = "integrator_gradient" + suffix + "_full_inner"
    func_params = [
        "s_dAB is the output [A | B] buffer (caller places); size 2*NUM_VEL*3*NUM_VEL = " + str(2 * n * 3 * n),
    ]
    if compute_x_kp1:
        func_params.append("s_x_kp1 is the next-state output (caller places); size NUM_POS + NUM_VEL")
    func_params += [
        "s_q is the vector of joint positions (NON-const: mutated across RK stages)",
        "s_qd is the vector of joint velocities (NON-const: mutated across RK stages)",
        "s_u is the vector of joint input torques",
        "s_df_du / s_dc_du / s_vaf / s_Minv / s_qdd are FD-grad in/out scratch (caller places)",
        "s_q_orig / s_qd_orig / s_stage_grad_qdd / s_D_qdd_stage are multi-stage scratch (caller places; s_D_qdd_stage is a de-aliased workspace band when spilled)",
        "s_dInt_q_6x6 / s_dInt_v_6x6 are the floating-base SE(3) dIntegrate blocks (unused fixed-base)",
        "s_temp is the FD-grad inner scratch pool (used when SCRATCH_IN_SMEM)",
        "d_workspace is the global scratch pool the FD-grad inner s_temp routes to (used when !SCRATCH_IN_SMEM)",
        "d_temp_spill is the id_du da_df band spill region (used when USE_DA_DF_SPILL)",
        "d_robotModel holds XImats/topology; gravity is the gravity constant; dt is the timestep",
    ]
    func_def_start = "void " + fname + "(T *s_dAB, "
    if compute_x_kp1:
        func_def_start += "T *s_x_kp1, "
    func_def_middle = ("T *s_q, T *s_qd, const T *s_u, T *s_df_du, T *s_dc_du, T *s_vaf, "
                       "T *s_Minv, T *s_qdd, T *s_q_orig, T *s_qd_orig, T *s_stage_grad_qdd, "
                       "T *s_D_qdd_stage, T *s_dInt_q_6x6, T *s_dInt_v_6x6, ")
    func_def_end = ("T *s_temp, T *d_workspace, T *d_temp_spill, "
                    "const robotModel<T> *d_robotModel, const T gravity, const T dt) {")
    if use_thread_group:
        func_def_start += "cgrps::thread_group tgrp, "
        func_params.insert(0, "tgrp is the handle to the thread_group running this function")
    func_def_middle, func_params = self.gen_insert_helpers_func_def_params(func_def_middle, func_params, -2)
    func_def = func_def_start + func_def_middle + func_def_end
    self.gen_add_func_doc("integrator gradient orchestration as a single inner-owns-placement device function",
                          ["Owns the FD-grad inner s_temp pool placement; the repoint covers every consumer below (incl. the XImats helper's sincos scratch and the per-stage XImats refresh)"],
                          func_params, None)
    self.gen_add_code_line("template <typename T, IntegratorType IT = IntegratorType::EULER, bool SCRATCH_IN_SMEM = true, bool USE_DA_DF_SPILL = false>")
    # __forceinline__ so the whole orchestration inlines into the calling kernel.
    # Under -rdc a separate __device__ wrapper keeps its RBD callees as distinct
    # functions whose regcount must fit the kernel's launch_bounds budget -> ptxas
    # regcount error. Inlining folds them into the kernel. See _fdsva_so.py:295-300.
    self.gen_add_code_line("__device__ __forceinline__")
    self.gen_add_code_line(func_def, True)
    # Inner owns the FD-grad pool placement; the repoint covers every consumer below
    # (incl. the XImats helper's sincos scratch and the multi-stage per-stage XImats
    # refresh), so no caller-side repoint. This is the migrated kernel line-744 case.
    self.gen_add_code_line("if constexpr(!SCRATCH_IN_SMEM){ s_temp = d_workspace; } else { (void)d_workspace; }")
    # XImats helper INSIDE the inner AFTER the repoint (no shared-helper edit) so its
    # sincos scratch follows the placed s_temp pool.
    self.gen_load_update_XImats_helpers_function_call(use_thread_group)
    # Compile-time IT dispatch: single-stage (Euler / SI-Euler) vs multi-stage
    # (Midpoint / RK3 / RK4). Per-rung band flags are passed as 'true'/'false'
    # literals through the stable FD-grad _inner_python composition surface.
    spill_flag = "USE_DA_DF_SPILL"
    self.gen_add_code_line(
        "if constexpr (IT == IntegratorType::EULER || IT == IntegratorType::SEMI_IMPLICIT_EULER) {", True
    )
    self.gen_integrator_gradient_inner_python(
        use_thread_group=use_thread_group,
        compute_x_kp1=compute_x_kp1,
        integrator_type="IT",
        s_dAB_name="s_dAB",
        s_x_kp1_name="s_x_kp1",
        d_temp_spill_name="d_temp_spill",
        temp_spill_flag_name=spill_flag,
    )
    self.gen_add_end_control_flow()
    self.gen_add_code_line("else {", True)
    self.gen_integrator_gradient_multistage(
        use_thread_group=use_thread_group,
        compute_x_kp1=compute_x_kp1,
        d_temp_spill_name="d_temp_spill",
        temp_spill_flag_name=spill_flag,
    )
    self.gen_add_end_control_flow()
    self.gen_add_end_function()


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
                       ("s_Minv", n * n), ("s_qdd", n),
                       # Floating-base SE(3) dIntegrate blocks (unused for fixed-base).
                       ("s_dInt_q_6x6", 36), ("s_dInt_v_6x6", 36)]
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
    # d_workspace holds the L2-pinned scratch the s_D_qdd_stage buffer spills to
    # at LITE/MINIMAL (mirrors forward_dynamics_gradient_kernel's d_workspace).
    func_def_start = ("void integrator_gradient" + suffix + "_kernel(T *d_dAB, " + sig_x_kp1 +
                      "unsigned char *d_workspace, const T *d_q_qd_u, const int stride_q_qd_u, ")
    func_def_end = "const robotModel<T> *d_robotModel, const T gravity, const T dt, const int NUM_TIMESTEPS) {"
    func_def = func_def_start + func_def_end
    if single_call_timing:
        func_def = func_def.replace("kernel(", "kernel_single_timing(")
    func_params.insert(1 if not compute_x_kp1 else 2,
                       "d_workspace is the L2-pinned global scratch for the spilled s_D_qdd_stage (LITE/MINIMAL tiers)")
    self.gen_add_func_doc("Computes the gradient of the integrator step per timestep" +
                          (" and the next state x_{k+1}" if compute_x_kp1 else ""),
                          [], func_params, None)
    self.gen_add_code_line("template <typename T, IntegratorType IT = IntegratorType::EULER, int RESOURCE_TIER = GRID_DEFAULT_RESOURCE_TIER>")
    self.gen_add_code_line("__global__")
    # Pin launch_bounds to MAX_PERF_LEVEL_THREADS (PERF cap), NOT tier_max_threads: the
    # integrator gradient is register-bound by its RBD callees, so the LITE/MINIMAL
    # thread bump starves them and ptxas errors under -rdc=true. Tier behavior here
    # is the s_D_qdd_stage smem spill, which is independent of launch_bounds.
    self.gen_add_code_line("__launch_bounds__(MAX_PERF_LEVEL_THREADS)")
    self.gen_add_code_line(func_def, True)
    inner_temp_full = self.gen_integrator_gradient_inner_temp_mem_size()
    inner_temp_selective = max(self.gen_direct_minv_inner_temp_mem_size(),
                               self.gen_inverse_dynamics_gradient_temp_layout()["selective_shared_count"])
    fb = self.robot.floating_base
    max_stages = _max_stages_in_use()
    d_qdd_count = max_stages * n * 3 * n

    def _emit_body(dqdd_in_smem, dab_in_smem, inner_level):
        # Surgical per-tier body. The 3 distinct buffers spill independently to
        # SEPARATE non-aliasing d_workspace sub-offsets (the integrator gradient
        # never runs concurrently with id_du/fd_du/fdsva_so, so it reuses those
        # sections) — the de-aliased multi-band layout is preserved, NOT collapsed:
        #   - s_D_qdd_stage (max_stages*nv*3nv) -> Dqdd region (offset 0) when !dqdd_in_smem  [caller-placed]
        #   - s_dAB output (2nv*3nv)            -> dAB region              when !dab_in_smem  [caller-placed]
        #   - the FD-grad inner s_temp POOL (the integrator_gradient_full_inner OWNS
        #     this placement via its SCRATCH_IN_SMEM template flag):
        #       inner_level 0: full smem (SCRATCH_IN_SMEM=true);
        #       1: da_df-band SELECTIVE spill (s_temp shrinks, only the id_du band
        #          leaves smem to d_temp_spill; SCRATCH_IN_SMEM=true, USE_DA_DF_SPILL=true);
        #       2: whole inner POOL -> inner region (SCRATCH_IN_SMEM=false; the inner
        #          repoints s_temp=d_workspace at its top — this is the migrated
        #          former kernel line-744 repoint).
        # The hot scaffold (s_dc_du / s_vaf / s_Minv) always stays in smem. The kernel
        # only slices the band base pointers + passes the per-rung flags as literals.
        inner_temp_size = (inner_temp_full if inner_level == 0
                           else (inner_temp_selective if inner_level == 1 else 0))
        extra_t_buffers = [("s_q_qd_u", 3 * n + fb)]
        if dab_in_smem:
            extra_t_buffers.append(("s_dAB", 2 * n * 3 * n))
        extra_t_buffers += [
            ("s_df_du", n * 2 * n),
            ("s_dc_du", n * 2 * n),
            ("s_vaf", 18 * n),
            ("s_Minv", n * n),
            ("s_qdd", n),
            # Multi-stage scratch — allocated for every IT (single-stage just doesn't use it).
            # s_q_orig holds the full nq pose (floating-base adds the quaternion slot).
            ("s_q_orig", n + fb),
            ("s_qd_orig", n),
            ("s_stage_grad_qdd", max_stages * n),
        ]
        if dqdd_in_smem:
            extra_t_buffers.append(("s_D_qdd_stage", d_qdd_count))
        extra_t_buffers += [
            # Floating-base 6x6 SE(3) dIntegrate blocks (Euler single-stage path).
            # For fixed-base these stay unused.
            ("s_dInt_q_6x6", 36),
            ("s_dInt_v_6x6", 36),
        ]
        if compute_x_kp1:
            extra_t_buffers.append(("s_x_kp1", 2 * n + fb))  # = nq + nv
        self.gen_XImats_helpers_temp_shared_memory_code(
            inner_temp_size, extra_t_buffers=extra_t_buffers, include_linalg_scratch=True,
        )
        self.gen_add_code_line(
            "T *s_q = s_q_qd_u; T *s_qd = &s_q_qd_u[" + str(n + fb) + "]; T *s_u = &s_q_qd_u[" + str(2 * n + fb) + "];"
        )
        # The kernel only SLICES the de-aliased workspace band base pointers and
        # passes per-rung flags; the full_inner OWNS the FD-grad inner s_temp pool
        # placement (the rung-2 whole-pool repoint is its SCRATCH_IN_SMEM=false path).
        # Per-rung flags are passed as 'true'/'false' literals.
        scratch_in_smem_expr = "false" if inner_level == 2 else "true"
        spill_flag = "true" if inner_level == 1 else "false"
        # d_temp_spill is the id_du da_df band region (rung 1); always declared so the
        # full_inner call can reference it (nullptr unless inner_level==1).
        self.gen_add_code_line("T *d_temp_spill = nullptr;")
        if use_thread_group:
            self.gen_add_code_line("cgrps::thread_group tgrp = TBD;")

        def _emit_spill_pointers(slot_expr):
            # Slice the de-aliased multi-band workspace base pointers. The 3 distinct
            # buffers (s_D_qdd_stage, s_dAB, the FD-grad inner pool) keep SEPARATE
            # non-aliasing workspace sub-offsets; only the FD-grad inner POOL repoint
            # moved into the full_inner (its SCRATCH_IN_SMEM=false path). The kernel
            # passes that pool base via d_workspace_pool_name below.
            if not dqdd_in_smem:
                self.gen_add_code_line(
                    "T *s_D_qdd_stage = reinterpret_cast<T *>(&d_workspace[" + slot_expr + "]);"
                )
            if not dab_in_smem:
                self.gen_add_code_line(
                    "T *s_dAB = reinterpret_cast<T *>(&d_workspace[" + slot_expr + " + GRID_INTEGRATOR_DU_DAB_OFFSET_BYTES<T>()]);"
                )
            if inner_level == 1:
                self.gen_add_code_line(
                    "d_temp_spill = reinterpret_cast<T *>(&d_workspace[" + slot_expr + " + GRID_INTEGRATOR_DU_INNER_OFFSET_BYTES<T>()]);"
                )

        def _emit_full_inner_call(slot_expr):
            # The FD-grad inner pool base (only consumed by the inner when
            # SCRATCH_IN_SMEM=false; nullptr otherwise).
            pool_name = ("reinterpret_cast<T *>(&d_workspace[" + slot_expr + " + GRID_INTEGRATOR_DU_INNER_OFFSET_BYTES<T>()])"
                         if inner_level == 2 else "nullptr")
            spill_name = "d_temp_spill" if inner_level == 1 else "nullptr"
            self.gen_integrator_gradient_full_inner_function_call(
                use_thread_group=use_thread_group,
                compute_x_kp1=compute_x_kp1,
                scratch_in_smem_expr=scratch_in_smem_expr,
                use_da_df_spill_expr=spill_flag,
                d_workspace_pool_name=pool_name,
                d_temp_spill_name=spill_name,
            )

        if not single_call_timing:
            self.gen_add_parallel_loop("k", "NUM_TIMESTEPS", use_thread_group, block_level=True)
            self.gen_kernel_load_inputs("q_qd_u", "stride_q_qd_u", str(3 * n + fb), use_thread_group)
            self.gen_add_code_line("// compute — the orchestration inner owns its FD-grad s_temp pool placement")
            _emit_spill_pointers("k * GRID_WORKSPACE_BYTES_PER_TIMESTEP<T>()")
            _emit_full_inner_call("k * GRID_WORKSPACE_BYTES_PER_TIMESTEP<T>()")
            self.gen_add_sync(use_thread_group)
            self.gen_kernel_save_result("dAB", str(2 * n * 3 * n), str(2 * n * 3 * n), use_thread_group)
            if compute_x_kp1:
                self.gen_kernel_save_result("x_kp1", str(2 * n + fb), str(2 * n + fb), use_thread_group)
            self.gen_add_end_control_flow()
        else:
            input_count = 3 * n + fb
            self.gen_kernel_load_inputs_single_timing("q_qd_u", str(input_count))
            _emit_spill_pointers("0")
            self.gen_add_code_line("// compute with NUM_TIMESTEPS as NUM_REPS for timing")
            self.gen_add_code_line("for (int rep = 0; rep < NUM_TIMESTEPS; rep++){", True)
            self.gen_anti_licm_input_reload("q_qd_u", str(input_count), use_thread_group, feedback_from="dAB")
            _emit_full_inner_call("0")
            self.gen_anti_licm_output_write("dAB")
            self.gen_add_end_control_flow()
            self.gen_kernel_save_result_single_timing("dAB", str(2 * n * 3 * n), use_thread_group)
            if compute_x_kp1:
                self.gen_kernel_save_result_single_timing("x_kp1", str(2 * n + fb), use_thread_group)

    # Per-tier surgical placement (perf, lite, minimal). When all three rungs
    # agree (small robots that fit at PERF), emit a single body; otherwise gate
    # per tier on RESOURCE_TIER.
    picks = getattr(self, "integrator_du_spill_tier_3way", (0, 0, 0))
    dqdd_smem = getattr(self, "integrator_du_dqdd_in_smem_per_tier", (True, True, True))
    dab_smem = getattr(self, "integrator_du_dab_in_smem_per_tier", (True, True, True))
    inner_lvl = getattr(self, "integrator_du_inner_level_per_tier", (0, 0, 0))
    if picks[0] == picks[1] == picks[2]:
        _emit_body(dqdd_smem[0], dab_smem[0], inner_lvl[0])
    else:
        for tier_idx, tier_name in enumerate(("TIER_PERF", "TIER_LITE", "TIER_MINIMAL")):
            head = ("if constexpr (RESOURCE_TIER == " + tier_name + ") {") if tier_idx == 0 else \
                   ("else if constexpr (RESOURCE_TIER == " + tier_name + ") {")
            self.gen_add_code_line(head, True)
            _emit_body(dqdd_smem[tier_idx], dab_smem[tier_idx], inner_lvl[tier_idx])
            self.gen_add_end_control_flow()
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
                       "hd_data->d_dAB," + kernel_args_x_kp1 + "hd_data->d_workspace,hd_data->d_q_qd_u,stride_q_qd_u,")
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
    workspace_bytes = ("GRID_WORKSPACE_BYTES_PER_TIMESTEP<T>()" if single_call_timing
                       else "GRID_WORKSPACE_BYTES_PER_TIMESTEP<T>()*static_cast<size_t>(num_timesteps)")
    self.gen_add_code_line("if (GRID_INTEGRATOR_DU_USES_WORKSPACE) {gpuErrchk(grid_begin_l2_persisting(0, hd_data->d_workspace, " + workspace_bytes + "));}")
    self.gen_add_code_lines(func_call_code)
    self.gen_add_code_line("if (GRID_INTEGRATOR_DU_USES_WORKSPACE) {gpuErrchk(grid_end_l2_persisting(0));}")
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
    # Inner-owns-placement orchestration inners (one per output kind). Emitted
    # before the kernels that call them. Mirrors id_du / fdsva_so full_inner.
    self.gen_integrator_gradient_full_inner(use_thread_group, compute_x_kp1=False)
    self.gen_integrator_gradient_full_inner(use_thread_group, compute_x_kp1=True)
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
