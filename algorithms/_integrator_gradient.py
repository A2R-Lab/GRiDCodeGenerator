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
    fd_grad = self.gen_forward_dynamics_gradient_inner_temp_mem_size()
    # MUJOCO_OUTPUT epilogue (floating only) reuses the (dead) FD-grad s_temp pool
    # as the 2n*3n mjx-output scratch band; guarantee the pool can hold it. For
    # every robot tested fd_grad >> 2n*3n so this max is a no-op; it only bumps the
    # pool for a hypothetical very-large-nv floating robot. Fixed-base is unchanged
    # (no epilogue), so its codegen stays byte-identical.
    if self.robot.floating_base:
        n = self.robot.get_num_vel()
        return max(fd_grad, 2 * n * 3 * n)
    return fd_grad


def gen_integrator_gradient_dAB_assembly(self, integrator_type="IT",
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
    self.gen_add_parallel_loop("ind", str(twoN * 3 * n))
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


def _emit_integrator_gradient_mjx_output(self, integrator_type, s_mjx_scratch="s_temp"):
    """Emit the MuJoCo (mjx) output-convention epilogue for the integrator gradient,
    transforming the pin dAB = [A | B] held in ``s_dAB`` to the mjx convention IN
    PLACE. Floating-base, single-stage (Euler / SI-Euler) only. Runs at the END of
    the single-stage path where the FD-grad ``s_temp`` pool is DEAD (reused as the
    2n*3n mjx-output scratch band) and ``s_qdd`` / ``s_qd`` / ``s_u`` / ``s_q`` are
    all live (no recompute needed, like fdsva_so's epilogue).

    ``s_dAB`` is 2n*3n COLUMN-major (``s_dAB[col*2n + row]``); output rows are
    ``[q_{k+1} tangent(n); qd_{k+1}(n)]``, input cols ``[dq(n) | dqd(n) | du(n)]``.
    Transcribed verbatim from docs/open-tasks/mjx_proto/proto_integ_grad_mjx.py
    (validated <1e-15 vs integrator_gradient_pin_to_mjx). The assembly:
      * BOTTOM rows reframe as a VELOCITY output (G applied to base-linear rows; the
        3 base-rotation q-cols pick up g_dot @ qd_{k+1,pin}); columns reframe by the
        input-conversion Jacobians (Ginv on base-linear cols; the _cross_cols
        velocity/force couplings on the 3 base-rotation cols).
      * TOP rows: angular + joint output rows are the pin top block reframed; the 3
        base-LINEAR output rows are overwritten with the mjx GLOBAL-add tangent
        (identity on the base-linear q col + dt*W, W=qd for euler, W=qd_{k+1} for si).
    """
    n = self.robot.get_num_vel()
    twoN = 2 * n
    si = "(" + _integrator_type_token(integrator_type) + " == IntegratorType::SEMI_IMPLICIT_EULER)"
    self.gen_add_code_line("// === mjx output convention (floating-base integrator gradient) ===")
    self.gen_add_code_line("T *s_mjx = " + s_mjx_scratch + ";   // 2n*3n mjx output band (dead FD-grad pool)")
    self.gen_add_code_line("const bool si_mjx = " + si + ";")
    # ---- single-thread assembly (correctness-first; nv small) ----
    # The inner C for/if blocks below carry their own balanced braces as literal
    # text; only the outer threadIdx guard uses the codegen indent tracker
    # (gen_add_end_control_flow). s_dAB is col-major P(o,c)=s_dAB[c*2n+o].
    self.gen_add_code_line("if (threadIdx.x == 0 && threadIdx.y == 0) {", True)
    for line in [
        # R (row-major R[3*i+j]) from the xyzw base quaternion s_q[3..6].
        "T qx = s_q[3], qy = s_q[4], qz = s_q[5], qw = s_q[6];",
        "T xx = qx*qx, yy = qy*qy, zz = qz*qz;",
        "T xy = qx*qy, xz = qx*qz, yz = qy*qz, wx = qw*qx, wy = qw*qy, wz = qw*qz;",
        "T R[9];",
        "R[0] = static_cast<T>(1) - static_cast<T>(2)*(yy+zz); R[1] = static_cast<T>(2)*(xy-wz);                    R[2] = static_cast<T>(2)*(xz+wy);",
        "R[3] = static_cast<T>(2)*(xy+wz);                    R[4] = static_cast<T>(1) - static_cast<T>(2)*(xx+zz); R[5] = static_cast<T>(2)*(yz-wx);",
        "R[6] = static_cast<T>(2)*(xz-wy);                    R[7] = static_cast<T>(2)*(yz+wx);                    R[8] = static_cast<T>(1) - static_cast<T>(2)*(xx+yy);",
        "T vlin[3] = {s_qd[0], s_qd[1], s_qd[2]};",
        "T ulin[3] = {s_u[0],  s_u[1],  s_u[2]};",
        "// qd_{k+1,pin}[lin] = (qd + dt*qdd)[lin] (for the g_dot velocity-output term)",
        "T qk1lin[3] = {s_qd[0] + dt*s_qdd[0], s_qd[1] + dt*s_qdd[1], s_qd[2] + dt*s_qdd[2]};",
        "// two output halves (top=0, bot=n); each column block reframed then G-row-rotated.",
        "for (int half = 0; half < 2; ++half) {",
        "  int base = half * " + str(n) + ";",
        "  // --- q-block (cols 0..n-1) ---",
        "  for (int c = 0; c < " + str(n) + "; ++c) {",
        "    T inner[" + str(n) + "];",
        "    for (int r = 0; r < " + str(n) + "; ++r) {",
        "      int o = base + r;",
        # (M @ Ginv)[r,c]: base-linear cols mix via Ginv[j,c]=R[c,j]; else identity.
        "      if (c < 3) { inner[r] = s_dAB[0*" + str(twoN) + "+o]*R[c*3+0] + s_dAB[1*" + str(twoN) + "+o]*R[c*3+1] + s_dAB[2*" + str(twoN) + "+o]*R[c*3+2]; }",
        "      else { inner[r] = s_dAB[c*" + str(twoN) + "+o]; }",
        "    }",
        # _cross_cols couplings on ANG cols (3..5): Jv_q on qd-block, Ju_q on u-block.
        "    if (c >= 3 && c < 6) {",
        "      int a = c - 3;",
        "      T ea[3] = {static_cast<T>(0),static_cast<T>(0),static_cast<T>(0)}; ea[a] = static_cast<T>(1);",
        "      T jc[3] = {ea[1]*vlin[2]-ea[2]*vlin[1], ea[2]*vlin[0]-ea[0]*vlin[2], ea[0]*vlin[1]-ea[1]*vlin[0]};",
        "      T uc[3] = {ea[1]*ulin[2]-ea[2]*ulin[1], ea[2]*ulin[0]-ea[0]*ulin[2], ea[0]*ulin[1]-ea[1]*ulin[0]};",
        "      for (int r = 0; r < " + str(n) + "; ++r) {",
        "        int o = base + r;",
        "        inner[r] += s_dAB[(" + str(n) + "+0)*" + str(twoN) + "+o]*(-jc[0]) + s_dAB[(" + str(n) + "+1)*" + str(twoN) + "+o]*(-jc[1]) + s_dAB[(" + str(n) + "+2)*" + str(twoN) + "+o]*(-jc[2]);",
        "        inner[r] += s_dAB[(" + str(2*n) + "+0)*" + str(twoN) + "+o]*(-uc[0]) + s_dAB[(" + str(2*n) + "+1)*" + str(twoN) + "+o]*(-uc[1]) + s_dAB[(" + str(2*n) + "+2)*" + str(twoN) + "+o]*(-uc[2]);",
        "      }",
        "    }",
        # G applied to base-linear rows: R @ inner[0:3].
        "    T l0 = inner[0], l1 = inner[1], l2 = inner[2];",
        "    inner[0] = R[0]*l0 + R[1]*l1 + R[2]*l2;",
        "    inner[1] = R[3]*l0 + R[4]*l1 + R[5]*l2;",
        "    inner[2] = R[6]*l0 + R[7]*l1 + R[8]*l2;",
        # g_dot @ qd_{k+1,pin} on the bottom-half ang q-cols (velocity output).
        "    if (half == 1 && c >= 3 && c < 6) {",
        "      int a = c - 3;",
        "      T ea[3] = {static_cast<T>(0),static_cast<T>(0),static_cast<T>(0)}; ea[a] = static_cast<T>(1);",
        "      T sk[3] = {ea[1]*qk1lin[2]-ea[2]*qk1lin[1], ea[2]*qk1lin[0]-ea[0]*qk1lin[2], ea[0]*qk1lin[1]-ea[1]*qk1lin[0]};",
        "      inner[0] += R[0]*sk[0] + R[1]*sk[1] + R[2]*sk[2];",
        "      inner[1] += R[3]*sk[0] + R[4]*sk[1] + R[5]*sk[2];",
        "      inner[2] += R[6]*sk[0] + R[7]*sk[1] + R[8]*sk[2];",
        "    }",
        "    for (int r = 0; r < " + str(n) + "; ++r) s_mjx[c*" + str(twoN) + "+(base+r)] = inner[r];",
        "  }",
        # --- qd-block (cols n..2n-1) and u-block (cols 2n..3n-1): plain M @ Ginv then G rows ---
    ]:
        self.gen_add_code_line(line)
    for off in (n, 2 * n):
        for line in [
            "  for (int cl = 0; cl < " + str(n) + "; ++cl) {",
            "    int c = " + str(off) + " + cl;",
            "    T inner[" + str(n) + "];",
            "    for (int r = 0; r < " + str(n) + "; ++r) {",
            "      int o = base + r;",
            "      if (cl < 3) { inner[r] = s_dAB[(" + str(off) + "+0)*" + str(twoN) + "+o]*R[cl*3+0] + s_dAB[(" + str(off) + "+1)*" + str(twoN) + "+o]*R[cl*3+1] + s_dAB[(" + str(off) + "+2)*" + str(twoN) + "+o]*R[cl*3+2]; }",
            "      else { inner[r] = s_dAB[c*" + str(twoN) + "+o]; }",
            "    }",
            "    T l0 = inner[0], l1 = inner[1], l2 = inner[2];",
            "    inner[0] = R[0]*l0 + R[1]*l1 + R[2]*l2;",
            "    inner[1] = R[3]*l0 + R[4]*l1 + R[5]*l2;",
            "    inner[2] = R[6]*l0 + R[7]*l1 + R[8]*l2;",
            "    for (int r = 0; r < " + str(n) + "; ++r) s_mjx[c*" + str(twoN) + "+(base+r)] = inner[r];",
            "  }",
        ]:
            self.gen_add_code_line(line)
    self.gen_add_code_line("}")  # end half loop
    # ---- base-LINEAR output rows of the TOP half: mjx GLOBAL-add tangent ----
    for line in [
        "for (int a = 0; a < 3; ++a) {",
        "  int o = a;",
        "  for (int c = 0; c < " + str(3 * n) + "; ++c) s_mjx[c*" + str(twoN) + "+o] = static_cast<T>(0);",
        "  s_mjx[a*" + str(twoN) + "+o] = static_cast<T>(1);   // identity on base-linear q col a",
        "  if (si_mjx) {",
        "    for (int c = 0; c < " + str(3 * n) + "; ++c) s_mjx[c*" + str(twoN) + "+o] += dt * s_mjx[c*" + str(twoN) + "+(" + str(n) + "+a)];",
        "  } else {",
        "    s_mjx[(" + str(n) + "+a)*" + str(twoN) + "+o] += dt;   // dt on qd base-linear col a",
        "  }",
        "}",
    ]:
        self.gen_add_code_line(line)
    self.gen_add_end_control_flow()  # if threadIdx == 0
    self.gen_add_sync()
    # ---- copy the mjx band back over s_dAB (block-parallel) ----
    self.gen_add_parallel_loop("ind", str(twoN * 3 * n))
    self.gen_add_code_line("s_dAB[ind] = s_mjx[ind];")
    self.gen_add_end_control_flow()
    self.gen_add_sync()


def gen_integrator_gradient_multistage(self, compute_x_kp1=False,
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
    self.gen_add_parallel_loop("ind", str(nq))
    self.gen_add_code_line("s_q_orig[ind] = s_q[ind];")
    self.gen_add_end_control_flow()
    self.gen_add_parallel_loop("ind", str(n))
    self.gen_add_code_line("s_qd_orig[ind] = s_qd[ind];")
    self.gen_add_end_control_flow()
    self.gen_add_sync()

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
                self.gen_add_serial_ops()
                self.gen_add_code_line(f"T v_scaled[{n}];")
                self.gen_add_code_line(f"for (int i = 0; i < {n}; ++i) v_scaled[i] = c_offset * dt * s_qd_orig[i];")
                self.gen_add_code_line(f"grid_integrate_floating_q<T, {nq}>(s_q_orig, v_scaled, s_q);")
                self.gen_add_end_control_flow()
                self.gen_add_parallel_loop("ind", str(n))
                self.gen_add_code_line(f"s_qd[ind] = s_qd_orig[ind] + c_offset * dt * s_stage_grad_qdd[{prev_offset} + ind];")
                self.gen_add_end_control_flow()
            else:
                self.gen_add_parallel_loop("ind", str(n))
                self.gen_add_code_line(f"s_q[ind] = s_q_orig[ind] + c_offset * dt * s_qd_orig[ind];")
                self.gen_add_code_line(f"s_qd[ind] = s_qd_orig[ind] + c_offset * dt * s_stage_grad_qdd[{prev_offset} + ind];")
                self.gen_add_end_control_flow()
            self.gen_add_sync()
            # Update XImats for the new s_q.
            self.gen_load_update_XImats_helpers_function_call()
            self.gen_add_sync()
            if fb:
                # Per-stage SE(3) dIntegrate blocks at v_dt = c*dt*qd_orig (the
                # q-perturbation increment for p.q = integrate(q_orig, c*dt*qd_orig)).
                # Reused buffers s_dInt_*_6x6 — consumed in this stage's D_qdd loop
                # below before the next stage overwrites them.
                self.gen_add_serial_ops()
                self.gen_add_code_line(f"T v_dt_stage[{n}];")
                self.gen_add_code_line(f"for (int i = 0; i < {n}; ++i) v_dt_stage[i] = c_offset * dt * s_qd_orig[i];")
                self.gen_add_code_line("grid_dIntegrate_q_block<T>(v_dt_stage, s_dInt_q_6x6);")
                self.gen_add_code_line("grid_dIntegrate_v_block<T>(v_dt_stage, s_dInt_v_6x6);")
                self.gen_add_end_control_flow()
                self.gen_add_sync()

        # Run FD-gradient at this stage's (s_q, s_qd, s_u).
        # After this call: s_qdd, s_Minv, s_df_du = stage `stage_num` values.
        self.gen_forward_dynamics_gradient_inner_python(
            use_qdd_Minv_input=False,
            s_df_du_name="s_df_du",
            d_temp_spill_name=d_temp_spill_name,
            temp_spill_flag_name=temp_spill_flag_name,
            d_f_ext_name="nullptr",  # integrator gradient does not thread external forces
        )
        self.gen_add_sync()

        # Always save this stage's qdd into s_stage_grad_qdd[stage_idx * n].
        # The next stage's FD-grad will overwrite s_qdd, so we need this
        # snapshot to (a) build p_{stage+1} on the next iteration and
        # (b) assemble x_{k+1} at the end when compute_x_kp1.
        self.gen_add_parallel_loop("ind", str(n))
        self.gen_add_code_line(
            f"s_stage_grad_qdd[{stage_idx * n} + ind] = s_qdd[ind];"
        )
        self.gen_add_end_control_flow()
        self.gen_add_sync()

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
        self.gen_add_parallel_loop("ind", str(n * three_n))
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
        self.gen_add_sync()

        self.gen_add_end_control_flow()  # close `if constexpr (gating)`

    # ----- Final assembly of dAB and optional x_kp1 -----
    # q_{k+1} = integrate(q, dt*qd) (Euler-style for every RK variant), so the
    # top rows are [dInt_q | dt*dInt_v | 0] at v_dt = dt*qd. For fixed-base these
    # reduce to [I | dt*I | 0].
    if fb:
        self.gen_add_serial_ops()
        self.gen_add_code_line(f"T v_dt_final[{n}];")
        self.gen_add_code_line(f"for (int i = 0; i < {n}; ++i) v_dt_final[i] = dt * s_qd_orig[i];")
        self.gen_add_code_line("grid_dIntegrate_q_block<T>(v_dt_final, s_dInt_q_6x6);")
        self.gen_add_code_line("grid_dIntegrate_v_block<T>(v_dt_final, s_dInt_v_6x6);")
        self.gen_add_end_control_flow()
        self.gen_add_sync()
    self.gen_add_code_line("// --- multi-stage gradient: assemble final dAB ---")
    self.gen_add_parallel_loop("ind", str(2 * n * three_n))
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
    self.gen_add_sync()

    # Optionally also assemble x_{k+1}.
    if compute_x_kp1:
        self.gen_add_code_line("// --- multi-stage gradient: assemble x_{k+1} ---")
        # v_{k+1} = qd + dt * sum(b_i * qdd_i); stage qdds live in s_stage_grad_qdd.
        self.gen_add_parallel_loop("ind", str(n))
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
        self.gen_add_sync()
        # q_{k+1} = integrate(q, dt*qd) (Euler-style q-update for every RK variant).
        if fb:
            self.gen_add_serial_ops()
            self.gen_add_code_line(f"T v_scaled_x[{n}];")
            self.gen_add_code_line(f"for (int i = 0; i < {n}; ++i) v_scaled_x[i] = dt * s_qd_orig[i];")
            self.gen_add_code_line(f"grid_integrate_floating_q<T, {nq}>(s_q_orig, v_scaled_x, s_x_kp1);")
            self.gen_add_end_control_flow()
        else:
            self.gen_add_parallel_loop("ind", str(n))
            self.gen_add_code_line("s_x_kp1[ind] = s_q_orig[ind] + dt * s_qd_orig[ind];")
            self.gen_add_end_control_flow()
        self.gen_add_sync()


def gen_integrator_gradient_inner_python(self, compute_x_kp1=False,
                                          integrator_type="IT", s_dAB_name="s_dAB",
                                          s_x_kp1_name="s_x_kp1",
                                          d_temp_spill_name="nullptr", temp_spill_flag_name="false",
                                          mujoco_output_expr=None):
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
        use_qdd_Minv_input=False,
        s_df_du_name="s_df_du",
        d_temp_spill_name=d_temp_spill_name,
        temp_spill_flag_name=temp_spill_flag_name,
        d_f_ext_name="nullptr",  # integrator gradient does not thread external forces
    )
    self.gen_add_sync()
    if fb:
        # Precompute the SE(3) dIntegrate blocks at the q-update increment v_dt.
        # Euler:    q_new = integrate(q, dt*qd)            -> v_dt = dt*qd
        # SI-Euler: q_new = integrate(q, dt*v_new), where  -> v_dt = dt*(qd + dt*qdd)
        #           v_new = qd + dt*qdd  (s_qdd holds qdd after the FD gradient).
        self.gen_add_serial_ops()
        self.gen_add_code_line(f"T v_dt_for_dInt[{n}];")
        self.gen_add_code_line("if constexpr (" + _integrator_type_token(integrator_type) + " == IntegratorType::SEMI_IMPLICIT_EULER) {")
        self.gen_add_code_line(f"    for (int i = 0; i < {n}; ++i) v_dt_for_dInt[i] = dt * (s_qd[i] + dt * s_qdd[i]);")
        self.gen_add_code_line("} else {")
        self.gen_add_code_line(f"    for (int i = 0; i < {n}; ++i) v_dt_for_dInt[i] = dt * s_qd[i];")
        self.gen_add_code_line("}")
        self.gen_add_code_line("grid_dIntegrate_q_block<T>(v_dt_for_dInt, s_dInt_q_6x6);")
        self.gen_add_code_line("grid_dIntegrate_v_block<T>(v_dt_for_dInt, s_dInt_v_6x6);")
        self.gen_add_end_control_flow()
        self.gen_add_sync()
    self.gen_integrator_gradient_dAB_assembly(
        integrator_type=integrator_type,
        s_dAB_name=s_dAB_name,
    )
    # Optionally also build x_kp1 (the value) from the s_qdd that FD-gradient
    # populated. Saves a redundant FD pass for the both-at-once API.
    if compute_x_kp1:
        self.gen_add_sync()
        self.gen_integrator_finish_function_call(
            integrator_type=integrator_type,
            updated_var_names=dict(s_x_kp1_name=s_x_kp1_name),
        )
    # ---- mjx output-convention epilogue (floating single-stage only) ----
    # Runs at the END where s_dAB holds the finalized pin gradient and the FD-grad
    # s_temp pool is DEAD (reused as the 2n*3n mjx scratch). s_qdd / s_qd / s_u /
    # s_q are all live (no recompute, like fdsva_so). The x_kp1 value (if built) is
    # converted to mjx by the GLOBAL-add retract + qd reframe + quat reorder.
    if fb and mujoco_output_expr is not None:
        self.gen_add_code_line("if constexpr (" + mujoco_output_expr + ") {", True)
        _emit_integrator_gradient_mjx_output(self, integrator_type, s_mjx_scratch="s_temp")
        if compute_x_kp1:
            # x_kp1 = [q (nq); qd (nv)]: base-position GLOBAL add (mjx retract), qd
            # output reframe by G, base quaternion xyzw->wxyz. s_q still holds the
            # ORIGINAL base position (integrate wrote OUT-of-place into s_x_kp1);
            # s_qd[0:3] is the raw mjx global base-linear velocity used for the step.
            self.gen_mjx_retract(s_x_kp1_name, "s_q", "s_qd", "dt")
            # qd_{k+1,mjx} = G qd_{k+1,pin}: rotate the base-linear velocity rows by R.
            # Parenthesize the qd-block base so the helper's [i] indexing binds to the
            # offset pointer (s_x_kp1 + NQ), not to the literal (operator precedence).
            self.gen_mjx_base_rotate("(" + s_x_kp1_name + " + " + str(self.robot.get_num_pos()) + ")", q_name="s_q")
            self.gen_add_code_lines([
                "// mjx x_kp1: base quaternion xyzw->wxyz (inverse of input reorder)",
                "if (threadIdx.x == 0 && threadIdx.y == 0) {", True,
                "T qw_out = " + s_x_kp1_name + "[6];",
                s_x_kp1_name + "[6] = " + s_x_kp1_name + "[5]; " + s_x_kp1_name + "[5] = " + s_x_kp1_name + "[4]; "
                    + s_x_kp1_name + "[4] = " + s_x_kp1_name + "[3]; " + s_x_kp1_name + "[3] = qw_out;",
            ])
            self.gen_add_end_control_flow()
            self.gen_add_sync()
        self.gen_add_end_control_flow()


def gen_integrator_gradient_device_function_call(self, compute_x_kp1=False,
                                                     scratch_in_smem_expr="true",
                                                     use_da_df_spill_expr="false",
                                                     d_workspace_pool_name="nullptr",
                                                     d_temp_spill_name="nullptr",
                                                     mujoco_output_expr=None):
    """Emit the call to `integrator_gradient_device` / `integrator_with_gradient_device`. Arg order MUST
    match the def in gen_integrator_gradient_device. The FD-grad inner POOL
    placement region (d_workspace) and the inverse_dynamics_gradient da_df band spill region
    (d_temp_spill) default to nullptr (unused under the matching if-constexpr); the
    kernel passes real pointers per tier. s_D_qdd_stage / s_dAB remain SEPARATE
    caller-placed pointers — they are threaded through unchanged. mujoco_output_expr
    (floating only) appends the trailing MUJOCO_OUTPUT template flag."""
    fname = ("integrator_with_gradient" if compute_x_kp1 else "integrator_gradient") + "_device"
    tmpl_flags = "<T, IT, " + scratch_in_smem_expr + ", " + use_da_df_spill_expr
    tmpl = (tmpl_flags + ", " + mujoco_output_expr + ">") if mujoco_output_expr is not None else (tmpl_flags + ">")
    start = fname + tmpl + "(s_dAB, "
    if compute_x_kp1:
        start += "s_x_kp1, "
    start += ("s_q, s_qd, s_u, s_df_du, s_dc_du, s_vaf, s_Minv, s_qdd, "
              "s_q_orig, s_qd_orig, s_stage_grad_qdd, s_D_qdd_stage, "
              "s_dInt_q_6x6, s_dInt_v_6x6, ")
    middle = self.gen_insert_helpers_function_call()
    end = ("s_temp, " + d_workspace_pool_name + ", " + d_temp_spill_name + ", "
           + "d_robotModel, gravity, dt);")
    self.gen_add_code_line(start + middle + end)


def gen_integrator_gradient_device(self, compute_x_kp1=False):
    """Emit `integrator_gradient_device` (compute_x_kp1=False) / `integrator_with_gradient_device` (compute_x_kp1=True) — the whole integrator
    gradient orchestration as ONE inner that OWNS its FD-grad scratch (s_temp) pool
    placement (inner-owns-placement; mirrors gen_inverse_dynamics_gradient_device /
    gen_fdsva_so_device). It wraps, in order:
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
      USE_DA_DF_SPILL : the FD-grad inner's inverse_dynamics_gradient band selectively spills its
                        da_dq..fxvi band to d_temp_spill (rung 1). Threaded through
                        the stable gen_forward_dynamics_gradient_inner_python
                        composition surface to the inverse_dynamics_gradient band sub-inner.

    Pointer params are caller-supplied (the kernel decides where the OUTPUT s_dAB
    and the multi-band scratch s_D_qdd_stage live — smem or de-aliased workspace
    sub-offsets — and hands in the spill regions): only the FD-grad inner s_temp
    POOL placement is the inner's call. s_q / s_qd are NON-const: the multi-stage
    path MUTATES them in shared memory across RK stages (and re-derives the per-stage
    XImats from the freshly-mutated s_q)."""
    n = self.robot.get_num_vel()
    fb = self.robot.floating_base
    fname = ("integrator_with_gradient" if compute_x_kp1 else "integrator_gradient") + "_device"
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
        "d_temp_spill is the inverse_dynamics_gradient da_df band spill region (used when USE_DA_DF_SPILL)",
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
    func_def_middle, func_params = self.gen_insert_helpers_func_def_params(func_def_middle, func_params, -2)
    func_def = func_def_start + func_def_middle + func_def_end
    self.gen_add_func_doc("integrator gradient orchestration as a single inner-owns-placement device function",
                          ["Owns the FD-grad inner s_temp pool placement; the repoint covers every consumer below (incl. the XImats helper's sincos scratch and the per-stage XImats refresh)"],
                          func_params, None)
    # MUJOCO_OUTPUT (floating only): compile-time mjx output-convention flag, LAST so
    # existing positional <T,IT,SCRATCH,SPILL> call sites are unaffected; default
    # false if-constexpr-elides the epilogue -> byte-identical PTX on the pin path.
    mjx_device = self.robot.floating_base
    if mjx_device:
        self.gen_add_code_line("template <typename T, IntegratorType IT = IntegratorType::EULER, bool SCRATCH_IN_SMEM = true, bool USE_DA_DF_SPILL = false, bool MUJOCO_OUTPUT = false>")
    else:
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
    self.gen_load_update_XImats_helpers_function_call()
    # Compile-time IT dispatch: single-stage (Euler / SI-Euler) vs multi-stage
    # (Midpoint / RK3 / RK4). Per-rung band flags are passed as 'true'/'false'
    # literals through the stable FD-grad _inner_python composition surface.
    spill_flag = "USE_DA_DF_SPILL"
    self.gen_add_code_line(
        "if constexpr (IT == IntegratorType::EULER || IT == IntegratorType::SEMI_IMPLICIT_EULER) {", True
    )
    self.gen_integrator_gradient_inner_python(
        compute_x_kp1=compute_x_kp1,
        integrator_type="IT",
        s_dAB_name="s_dAB",
        s_x_kp1_name="s_x_kp1",
        d_temp_spill_name="d_temp_spill",
        temp_spill_flag_name=spill_flag,
        mujoco_output_expr=("MUJOCO_OUTPUT" if mjx_device else None),
    )
    self.gen_add_end_control_flow()
    self.gen_add_code_line("else {", True)
    if mjx_device:
        # Multi-stage RK + mjx is a clean-break deferral (the 2nd-order-free chain
        # rule still needs the mjx input/output reparam derived per stage). Refuse
        # rather than emit a silently-wrong tensor.
        self.gen_add_code_line("static_assert(!MUJOCO_OUTPUT, "
                               "\"integrator_gradient MUJOCO_OUTPUT is single-stage (EULER/SEMI_IMPLICIT_EULER) only; \"")
        self.gen_add_code_line("              \"multi-stage RK mjx is deferred.\");")
    self.gen_integrator_gradient_multistage(
        compute_x_kp1=compute_x_kp1,
        d_temp_spill_name="d_temp_spill",
        temp_spill_flag_name=spill_flag,
    )
    self.gen_add_end_control_flow()
    self.gen_add_end_function()


def gen_integrator_gradient_kernel(self, compute_x_kp1=False, single_call_timing=False):
    n = self.robot.get_num_vel()
    nq = self.robot.get_num_pos()
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
    func_def_start = ("void " + ("integrator_with_gradient" if compute_x_kp1 else "integrator_gradient") + "_kernel(T *d_dAB, " + sig_x_kp1 +
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
    # MUJOCO_OUTPUT (floating only): compile-time mjx flag appended LAST after
    # RESOURCE_TIER so existing <T,IT,TIER> call sites are unaffected; default false
    # if-constexpr-elides the input-convert + dAB epilogue -> byte-identical pin PTX.
    mjx_kernel = self.robot.floating_base
    if mjx_kernel:
        self.gen_add_code_line("template <typename T, IntegratorType IT = IntegratorType::EULER, int RESOURCE_TIER = GRID_DEFAULT_RESOURCE_TIER, bool MUJOCO_OUTPUT = false>")
    else:
        self.gen_add_code_line("template <typename T, IntegratorType IT = IntegratorType::EULER, int RESOURCE_TIER = GRID_DEFAULT_RESOURCE_TIER>")
    self.gen_add_code_line("__global__")
    # Pin launch_bounds to MAX_PERF_LEVEL_THREADS (PERF cap), NOT tier_max_threads: the
    # integrator gradient is register-bound by its RBD callees, so the LITE/MINIMAL
    # thread bump starves them and ptxas errors under -rdc=true. Tier behavior here
    # is the s_D_qdd_stage smem spill, which is independent of launch_bounds.
    self.gen_add_code_line("__launch_bounds__(MAX_PERF_LEVEL_THREADS)")
    self.gen_add_code_line(func_def, True)
    inner_temp_full = self.gen_integrator_gradient_inner_temp_mem_size()
    # The da_df-band SELECTIVE inner level only exists for the SPARSE (non-mimic)
    # inverse_dynamics_gradient inner. The MIMIC inner is a dense serial fold that
    # ignores USE_DA_DF_SPILL and always writes its full pool, so it cannot shrink:
    # size the "selective" pool to the full inner there (matches GCG.py's
    # _integrator_gradient_inner_selective). Mimic never emits inner_level 1 (see GCG's
    # mimic-aware integrator_gradient_inner_level_per_tier), so this only guards against a
    # mis-sized pool if that invariant ever changes.
    inner_temp_selective = (inner_temp_full if self.robot_has_mimic_joints()
                            else max(self.gen_minv_inner_temp_mem_size(),
                                     self.gen_inverse_dynamics_gradient_temp_layout()["selective_shared_count"]))
    fb = self.robot.floating_base
    max_stages = _max_stages_in_use()
    d_qdd_count = max_stages * n * 3 * n

    def _emit_body(dqdd_in_smem, dab_in_smem, inner_level):
        # Surgical per-tier body. The 3 distinct buffers spill independently to
        # SEPARATE non-aliasing d_workspace sub-offsets (the integrator gradient
        # never runs concurrently with inverse_dynamics_gradient/forward_dynamics_gradient/fdsva_so, so it reuses those
        # sections) — the de-aliased multi-band layout is preserved, NOT collapsed:
        #   - s_D_qdd_stage (max_stages*nv*3nv) -> Dqdd region (offset 0) when !dqdd_in_smem  [caller-placed]
        #   - s_dAB output (2nv*3nv)            -> dAB region              when !dab_in_smem  [caller-placed]
        #   - the FD-grad inner s_temp POOL (the integrator_gradient_device OWNS
        #     this placement via its SCRATCH_IN_SMEM template flag):
        #       inner_level 0: full smem (SCRATCH_IN_SMEM=true);
        #       1: da_df-band SELECTIVE spill (s_temp shrinks, only the inverse_dynamics_gradient band
        #          leaves smem to d_temp_spill; SCRATCH_IN_SMEM=true, USE_DA_DF_SPILL=true);
        #       2: whole inner POOL -> inner region (SCRATCH_IN_SMEM=false; the inner
        #          repoints s_temp=d_workspace at its top).
        # The hot scaffold (s_dc_du / s_vaf / s_Minv) always stays in smem. The kernel
        # only slices the band base pointers + passes the per-rung flags as literals.
        inner_temp_size = (inner_temp_full if inner_level == 0
                           else (inner_temp_selective if inner_level == 1 else 0))
        # s_vaf is body-indexed (NB bodies, stride 6). For a MIMIC robot (fixed
        # base) NB > nv, so the composed FD-grad inner's ID sub-inner writes
        # 18*NB entries — size it 18*NB to keep those writes from overflowing
        # into the adjacent s_Minv/s_qdd buffers (mirrors forward_dynamics_gradient's kernel sizing
        # in _forward_dynamics_gradient.py). Non-mimic keeps 18*n (byte-identical;
        # floating nv > NB).
        vaf_cnt = 18 * (self.robot.get_num_joints() if self.robot_has_mimic_joints() else n)
        # Canonical per-timestep INPUT packing (mirrors id/crba/aba/forward_dynamics):
        # q, qd, u each in a NUM_JOINTS(=nq)-wide slot at stride 3*nq; slice qd at nq,
        # u at 2*nq. For a FIXED base nq==nv -> 3*nq == 3*nv+fb byte-identical; for a
        # FLOATING base nq=nv+1 the old nv-strided u offset (2*nv+fb) under-read by
        # nq-nv and mis-sliced u -- the floating B=1 + batch input bug. The dAB OUTPUT
        # is genuinely tangent-space: 2*nv x 3*nv real values, per-timestep stride
        # 2*nv*3*nv (the d_dAB buffer + binding transfer are 2*nv*3*nv-sized), so the
        # save stride stays 2*nv*3*nv. The x_kp1 output is [q (nq); qd (nv)] = nq+nv.
        input_count = 3 * nq
        extra_t_buffers = [("s_q_qd_u", input_count)]
        if dab_in_smem:
            extra_t_buffers.append(("s_dAB", 2 * n * 3 * n))
        extra_t_buffers += [
            ("s_df_du", n * 2 * n),
            ("s_dc_du", n * 2 * n),
            ("s_vaf", vaf_cnt),
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
            "T *s_q = s_q_qd_u; T *s_qd = &s_q_qd_u[" + str(nq) + "]; T *s_u = &s_q_qd_u[" + str(2 * nq) + "];"
        )
        # The kernel only SLICES the de-aliased workspace band base pointers and
        # passes per-rung flags; the device OWNS the FD-grad inner s_temp pool
        # placement (the rung-2 whole-pool repoint is its SCRATCH_IN_SMEM=false path).
        # Per-rung flags are passed as 'true'/'false' literals.
        scratch_in_smem_expr = "false" if inner_level == 2 else "true"
        spill_flag = "true" if inner_level == 1 else "false"
        # d_temp_spill is the inverse_dynamics_gradient da_df band region (rung 1); always declared so the
        # device call can reference it (nullptr unless inner_level==1).
        self.gen_add_code_line("T *d_temp_spill = nullptr; (void)d_temp_spill;")

        def _emit_spill_pointers(slot_expr):
            # Slice the de-aliased multi-band workspace base pointers. The 3 distinct
            # buffers (s_D_qdd_stage, s_dAB, the FD-grad inner pool) keep SEPARATE
            # non-aliasing workspace sub-offsets; only the FD-grad inner POOL repoint
            # moved into the device (its SCRATCH_IN_SMEM=false path). The kernel
            # passes that pool base via d_workspace_pool_name below.
            if not dqdd_in_smem:
                self.gen_add_code_line(
                    "T *s_D_qdd_stage = reinterpret_cast<T *>(&d_workspace[" + slot_expr + "]);"
                )
            if not dab_in_smem:
                self.gen_add_code_line(
                    "T *s_dAB = reinterpret_cast<T *>(&d_workspace[" + slot_expr + " + GRID_INTEGRATOR_GRADIENT_DAB_OFFSET_BYTES<T>()]);"
                )
            if inner_level == 1:
                self.gen_add_code_line(
                    "d_temp_spill = reinterpret_cast<T *>(&d_workspace[" + slot_expr + " + GRID_INTEGRATOR_GRADIENT_INNER_OFFSET_BYTES<T>()]);"
                )

        def _emit_device_call(slot_expr):
            # The FD-grad inner pool base (only consumed by the inner when
            # SCRATCH_IN_SMEM=false; nullptr otherwise).
            pool_name = ("reinterpret_cast<T *>(&d_workspace[" + slot_expr + " + GRID_INTEGRATOR_GRADIENT_INNER_OFFSET_BYTES<T>()])"
                         if inner_level == 2 else "nullptr")
            spill_name = "d_temp_spill" if inner_level == 1 else "nullptr"
            self.gen_integrator_gradient_device_function_call(
                compute_x_kp1=compute_x_kp1,
                scratch_in_smem_expr=scratch_in_smem_expr,
                use_da_df_spill_expr=spill_flag,
                d_workspace_pool_name=pool_name,
                d_temp_spill_name=spill_name,
                mujoco_output_expr=("MUJOCO_OUTPUT" if mjx_kernel else None),
            )

        # MUJOCO_OUTPUT: convert the mjx-frame inputs (quat wxyz->xyzw, base-linear
        # velocity R^T, force R^T) to the pin frame so the RBD callees + the dAB
        # epilogue see pin quantities. The qdd is the kernel's own output (not an
        # input), so only q/qd/u are converted. Must follow the load + precede the
        # XImats build inside the device call. No-op on the pin path.
        def _emit_mjx_input_convert():
            if mjx_kernel:
                self.gen_add_code_line("if constexpr (MUJOCO_OUTPUT) {", True)
                self.gen_mjx_input_convert(q_name="s_q", qd_name="s_qd", u_name="s_u")
                self.gen_add_end_control_flow()

        if not single_call_timing:
            self.gen_add_parallel_loop("k", "NUM_TIMESTEPS", block_level=True)
            self.gen_kernel_load_inputs("q_qd_u",str(input_count),stride="stride_q_qd_u")
            _emit_mjx_input_convert()
            self.gen_add_code_line("// compute — the orchestration inner owns its FD-grad s_temp pool placement")
            _emit_spill_pointers("k * GRID_WORKSPACE_BYTES_PER_TIMESTEP<T>()")
            _emit_device_call("k * GRID_WORKSPACE_BYTES_PER_TIMESTEP<T>()")
            self.gen_add_sync()
            self.gen_kernel_save_result("dAB",str(2 * n * 3 * n),stride=str(2 * n * 3 * n))
            if compute_x_kp1:
                self.gen_kernel_save_result("x_kp1",str(nq + n),stride=str(nq + n))
            self.gen_add_end_control_flow()
        else:
            self.gen_kernel_load_inputs("q_qd_u",str(input_count))
            _emit_mjx_input_convert()
            _emit_spill_pointers("0")
            self.gen_add_code_line("// compute with NUM_TIMESTEPS as NUM_REPS for timing")
            self.gen_add_code_line("for (int rep = 0; rep < NUM_TIMESTEPS; rep++){", True)
            self.gen_anti_licm_input_reload("q_qd_u", str(input_count), feedback_from="dAB")
            _emit_device_call("0")
            self.gen_anti_licm_output_write("dAB")
            self.gen_add_end_control_flow()
            self.gen_kernel_save_result("dAB",str(2 * n * 3 * n))
            if compute_x_kp1:
                self.gen_kernel_save_result("x_kp1",str(nq + n))

    # Per-tier surgical placement (perf, lite, minimal). When all three rungs
    # agree (small robots that fit at PERF), emit a single body; otherwise gate
    # per tier on RESOURCE_TIER.
    # NOTE: this is the one tier-dispatch site that does NOT route through the
    # shared gen_tier_dispatch helper (B+C §1.2): the collapse predicate keys on
    # `picks` alone, but the body indexes THREE parallel per-tier tuples
    # (dqdd_smem / dab_smem / inner_lvl) by tier position — not by the pick value
    # — so the helper's value-based emit_body_fn(pick) contract doesn't fit.
    # Kept bespoke (the plan explicitly allows this for the irregular sites).
    picks = getattr(self, "integrator_gradient_spill_tier_3way", (0, 0, 0))
    dqdd_smem = getattr(self, "integrator_gradient_dqdd_in_smem_per_tier", (True, True, True))
    dab_smem = getattr(self, "integrator_gradient_dab_in_smem_per_tier", (True, True, True))
    inner_lvl = getattr(self, "integrator_gradient_inner_level_per_tier", (0, 0, 0))
    if picks[0] == picks[1] == picks[2]:
        _emit_body(dqdd_smem[0], dab_smem[0], inner_lvl[0])
    else:
        for tier_idx, tier_name in enumerate(("TIER_SHARED", "TIER_LITE", "TIER_MINIMAL")):
            head = ("if constexpr (RESOURCE_TIER == " + tier_name + ") {") if tier_idx == 0 else \
                   ("else if constexpr (RESOURCE_TIER == " + tier_name + ") {")
            self.gen_add_code_line(head, True)
            _emit_body(dqdd_smem[tier_idx], dab_smem[tier_idx], inner_lvl[tier_idx])
            self.gen_add_end_control_flow()
    self.gen_add_end_function()


def gen_integrator_gradient_host(self, mode=0, compute_x_kp1=False):
    single_call_timing = mode == 1
    compute_only = mode == 2
    base_name = ("integrator_with_gradient" if compute_x_kp1 else "integrator_gradient")
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
    # MUJOCO_OUTPUT (floating only) host flag, LAST: forwarded to the kernel launch
    # (naming IT + the tier positionally to reach the trailing flag). Default false
    # -> byte-identical pin codegen. Binding calls grid::integrator_gradient<T, IT,
    # KIND, /*MUJOCO_OUTPUT=*/true>.
    mjx_host = self.robot.floating_base
    if mjx_host:
        self.gen_add_code_line("template <typename T, IntegratorType IT = IntegratorType::EULER, gridDataKind KIND = GRID_DATA_ALL, bool MUJOCO_OUTPUT = false>")
    else:
        self.gen_add_code_line("template <typename T, IntegratorType IT = IntegratorType::EULER, gridDataKind KIND = GRID_DATA_ALL>")
    self.gen_add_code_line("__host__")
    self.gen_add_code_line(func_def_start)
    self.gen_add_code_line(func_def_end, True)
    self.gen_add_code_line("static_assert(KIND == GRID_DATA_ALL || KIND == GRID_DATA_DYNAMICS, \"" + base_name + " requires all-data or dynamics gridData\");")
    kernel_args_x_kp1 = "hd_data->d_x_kp1," if compute_x_kp1 else ""
    kernel_tmpl = ("<T, IT, GRID_DEFAULT_RESOURCE_TIER, MUJOCO_OUTPUT>" if mjx_host else "<T, IT>")
    func_call_start = (base_name + "_kernel" + kernel_tmpl + "<<<block_dimms,thread_dimms,INTEGRATOR_DU_DYNAMIC_SHARED_MEM_BYTES<T>()>>>(" +
                       "hd_data->d_dAB," + kernel_args_x_kp1 + "hd_data->d_workspace,hd_data->d_q_qd_u,stride_q_qd_u,")
    func_call_end = "d_robotModel,gravity,dt,num_timesteps);"
    if single_call_timing:
        func_call_start = func_call_start.replace(base_name + "_kernel" + kernel_tmpl,
                                                  base_name + "_kernel_single_timing" + kernel_tmpl)
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
    self.gen_add_code_line("if (GRID_INTEGRATOR_GRADIENT_USES_WORKSPACE) {gpuErrchk(grid_begin_l2_persisting(0, hd_data->d_workspace, " + workspace_bytes + "));}")
    self.gen_add_code_lines(func_call_code)
    self.gen_add_code_line("if (GRID_INTEGRATOR_GRADIENT_USES_WORKSPACE) {gpuErrchk(grid_end_l2_persisting(0));}")
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


# Floating-base SE(3) scratch the hessian inner carves from the (post-fdsva_so)
# s_temp pool: w[6] + dInt_v[36] + d2Int_qv[216] + d2Int_vv[216].
FLOATING_HESSIAN_SE3_SCRATCH = 6 + 36 + 216 + 216


def gen_integrator_hessian_device_floating(self):
    """Emit the FLOATING-BASE body of integrator_hessian_device (after the
    EULER/SI-EULER static_assert). Composes fdsva_so_device (the four 2nd-order
    forward-dynamics blocks D2qdd in s_df2, plus the first-order s_df_du / s_Minv
    / s_qdd) with the SE(3) retract second derivatives, exactly mirroring
    RBDReference._PlantMixin.plant_step_hessian (floating branch):

      velocity rows : H[nv+i, a, b] = dt * D2qdd[i, b, a]   (a/b TRANSPOSED vs
                      the fixed sweep: perturb in a, gradient-column in b).
      position rows, EULER (nonzero only when a is in the qd block, a=nv+a'):
            H[i, nv+a', b in q ] = dt   * d2Int_qv[i, b , a']
            H[i, nv+a', b in qd] = dt^2 * d2Int_vv[i, b', a']  (b=nv+b')
        everything else 0; NOT symmetrized (the FD-of-pinocchio ground truth is
        asymmetric -- only the qd perturbation axis carries the cross term).
      position rows, SI-EULER (Vgrad = dv_{k+1}/dz = [dt*fd_dq | I+dt*fd_dqd | dt*Minv]):
            t1 (b in q): dt  * sum_c d2Int_qv[i,b,c] * Vgrad[c,a]
            t2:          dt^2* sum_{m,c} d2Int_vv[i,m,c] * Vgrad[c,a] * Vgrad[m,b]
            t3:          dt^2* sum_m dInt_v[i,m] * D2qdd[m,a,b]
        where the SE(3) blocks (dInt_v, d2Int_qv, d2Int_vv) are nonzero ONLY in
        the free-flyer 6x6(x6) corner; revolute rows/cols of dInt_v are identity
        (so t3's m-sum picks up the i-th D2qdd block directly for i>=6).

    All d2Int / dInt SE(3) blocks are finite-differenced in DOUBLE (see
    grid_d2Integrate_block) then stored as T, so a float32 kernel matches the
    float64 oracle. The blocks are computed ONCE block-cooperatively into the
    post-fdsva_so s_temp pool, then a single fully-parallel sweep over the
    2*nv*3*nv*3*nv output cells reads them (each thread owns one output cell)."""
    n = self.robot.get_num_vel()
    nz = 3 * n
    nn = n * n
    nnn = n * n * n
    # fdsva_so (D2qdd blocks in s_df2; first-order s_df_du=[fd_dq|fd_dqd], s_Minv,
    # s_qdd). Same composition + flags as the fixed-base path.
    self.gen_fdsva_so_device_function_call(
        scratch_in_smem_expr="SCRATCH_IN_SMEM",
        fd_grad_use_spill_expr="FD_GRAD_USE_SPILL",
        contract_in_smem_expr="CONTRACT_IN_SMEM",
        d_workspace_pool_name="d_workspace",
        d_fd_grad_spill_name="d_fd_grad_spill",
        s_fdsva_temp_name="s_fdsva_temp")
    self.gen_add_sync()
    self.gen_add_code_line("T *d2a_dqdq = s_df2;")
    self.gen_add_code_line("T *d2a_dvdq = &s_df2[" + str(nnn) + "];")
    self.gen_add_code_line("T *d2a_dvdv = &s_df2[" + str(2 * nnn) + "];")
    self.gen_add_code_line("T *d2a_dtdq = &s_df2[" + str(3 * nnn) + "];")
    self.gen_add_code_line("const bool si = (IT == IntegratorType::SEMI_IMPLICIT_EULER);")
    self.gen_add_code_line("const T dt2 = dt * dt;")
    # The fdsva_so pool (s_temp; routed to d_workspace when !SCRATCH_IN_SMEM) is
    # free after the call returns -- carve the tiny SE(3) scratch from its front.
    self.gen_add_code_line("// SE(3) retract scratch (block-shared, carved from the freed fdsva_so pool).")
    self.gen_add_code_line("T *s_se3 = (SCRATCH_IN_SMEM) ? s_temp : d_workspace;")
    self.gen_add_code_line("T *s_w       = &s_se3[0];")
    self.gen_add_code_line("T *s_dInt_v  = &s_se3[6];          // 6x6 first-order dIntegrate_v")
    self.gen_add_code_line("T *s_d2Int_qv = &s_se3[6 + 36];    // 6x6x6 [o*36 + j*6 + k]")
    self.gen_add_code_line("T *s_d2Int_vv = &s_se3[6 + 36 + 216];")
    # w = the q-update increment (free-flyer 6): Euler dt*qd; SI-Euler dt*(qd+dt*qdd).
    self.gen_add_serial_ops()
    self.gen_add_code_line("for (int m = 0; m < 6; ++m) s_w[m] = si ? (dt * (s_qd[m] + dt * s_qdd[m])) : (dt * s_qd[m]);")
    self.gen_add_end_control_flow()
    self.gen_add_sync()
    # Compute the SE(3) blocks once, block-cooperatively. dInt_v (6x6) is one
    # thread; the two 6x6x6 d2Int tensors are FD'd in double (each thread owns
    # the full direction-k stencil over the 6 free-flyer directions).
    self.gen_add_serial_ops()
    self.gen_add_code_line("double w_d[6]; for (int m = 0; m < 6; ++m) w_d[m] = static_cast<double>(s_w[m]);")
    self.gen_add_code_line("double dIv_d[36]; grid_dIntegrate_v_block<double>(w_d, dIv_d);")
    self.gen_add_code_line("for (int m = 0; m < 36; ++m) s_dInt_v[m] = static_cast<T>(dIv_d[m]);")
    self.gen_add_code_line("grid_d2Integrate_block<T, true >(s_w, s_d2Int_qv);")
    self.gen_add_code_line("grid_d2Integrate_block<T, false>(s_w, s_d2Int_vv);")
    self.gen_add_end_control_flow()
    self.gen_add_sync()
    # ---- one fully-parallel sweep over the 2*nv*nz*nz output cells ----
    # D2qdd[m,p,q] (z-block lookup, same structure as the fixed-base assembly).
    # Defined as a device lambda so the velocity-row transpose and the SI-Euler
    # t3 contraction both read it without recompute.
    self.gen_add_code_line("auto D2qdd = [&](int m, int p, int q) -> T {")
    self.gen_add_code_line("    int pblk = p / " + str(n) + ", qblk = q / " + str(n) + ";")
    self.gen_add_code_line("    int pj = p % " + str(n) + ", qk = q % " + str(n) + ";")
    self.gen_add_code_line("    if (pblk == 0 && qblk == 0) return d2a_dqdq[m*" + str(nn) + " + pj*" + str(n) + " + qk];")
    self.gen_add_code_line("    if (pblk == 1 && qblk == 1) return d2a_dvdv[m*" + str(nn) + " + pj*" + str(n) + " + qk];")
    self.gen_add_code_line("    if (pblk == 1 && qblk == 0) return d2a_dvdq[m*" + str(nn) + " + pj*" + str(n) + " + qk];")
    self.gen_add_code_line("    if (pblk == 0 && qblk == 1) return d2a_dvdq[m*" + str(nn) + " + qk*" + str(n) + " + pj];")
    self.gen_add_code_line("    if (pblk == 2 && qblk == 0) return d2a_dtdq[m*" + str(nn) + " + pj*" + str(n) + " + qk];")
    self.gen_add_code_line("    if (pblk == 0 && qblk == 2) return d2a_dtdq[m*" + str(nn) + " + qk*" + str(n) + " + pj];")
    self.gen_add_code_line("    return static_cast<T>(0);")
    self.gen_add_code_line("};")
    # Vgrad[c, axis] = dv_{k+1}/dz: a in q -> dt*fd_dq[c,a]; a in qd -> (c==a')+dt*fd_dqd[c,a'];
    #                 a in u -> dt*Minv[c,a']. fd_dq/fd_dqd are column-major (s_df_du[col*n+row]);
    #                 Minv is SYMMETRIC_UPPER. Only used by SI-Euler.
    self.gen_add_code_line("auto Vgrad = [&](int c, int axis) -> T {")
    self.gen_add_code_line("    int blk = axis / " + str(n) + ", a2 = axis % " + str(n) + ";")
    self.gen_add_code_line("    if (blk == 0) return dt * s_df_du[a2*" + str(n) + " + c];")
    self.gen_add_code_line("    if (blk == 1) return ((c == a2) ? static_cast<T>(1) : static_cast<T>(0)) + dt * s_df_du[" + str(nn) + " + a2*" + str(n) + " + c];")
    self.gen_add_code_line("    int midx = (c <= a2) * (a2*" + str(n) + " + c) + (c > a2) * (c*" + str(n) + " + a2);")
    self.gen_add_code_line("    return dt * s_Minv[midx];")
    self.gen_add_code_line("};")
    # dInt_v[i,m]: free-flyer 6x6 corner from s_dInt_v; identity on the revolute block.
    self.gen_add_code_line("auto dIntv = [&](int i, int m) -> T {")
    self.gen_add_code_line("    if (i < 6 && m < 6) return s_dInt_v[i*6 + m];")
    self.gen_add_code_line("    return (i == m) ? static_cast<T>(1) : static_cast<T>(0);")
    self.gen_add_code_line("};")
    self.gen_add_parallel_loop("ind", str(2 * n * nz * nz))
    self.gen_add_code_line("int o = ind / " + str(nz * nz) + ";")
    self.gen_add_code_line("int a = (ind / " + str(nz) + ") % " + str(nz) + ";")
    self.gen_add_code_line("int b = ind % " + str(nz) + ";")
    self.gen_add_code_line("T val = static_cast<T>(0);")
    self.gen_add_code_line("if (o >= " + str(n) + ") {", True)
    self.gen_add_code_line("// velocity rows: dt * D2qdd[i, b, a]  (a/b transposed vs the fixed sweep).")
    self.gen_add_code_line("int i = o - " + str(n) + ";")
    self.gen_add_code_line("val = dt * D2qdd(i, b, a);")
    self.gen_add_end_control_flow()
    self.gen_add_code_line("else if (!si) {", True)
    self.gen_add_code_line("// position rows, EULER: nonzero only when a is in the qd block.")
    self.gen_add_code_line("int i = o;")
    self.gen_add_code_line("if (a >= " + str(n) + " && a < " + str(2 * n) + " && i < 6) {", True)
    self.gen_add_code_line("int aL = a - " + str(n) + ";")
    self.gen_add_code_line("if (b < " + str(n) + ") { if (b < 6 && aL < 6) val = dt * s_d2Int_qv[i*36 + b*6 + aL]; }")
    self.gen_add_code_line("else if (b < " + str(2 * n) + ") { int bL = b - " + str(n) + "; if (bL < 6 && aL < 6) val = dt2 * s_d2Int_vv[i*36 + bL*6 + aL]; }")
    self.gen_add_end_control_flow()
    self.gen_add_end_control_flow()
    self.gen_add_code_line("else {", True)
    self.gen_add_code_line("// position rows, SI-EULER: t1 + t2 + t3.")
    self.gen_add_code_line("int i = o;")
    self.gen_add_code_line("T acc = static_cast<T>(0);")
    self.gen_add_code_line("if (i < 6) {", True)
    self.gen_add_code_line("// t1 (b in q only): dt * sum_c d2Int_qv[i,b,c] * Vgrad[c,a].")
    self.gen_add_code_line("if (b < " + str(n) + " && b < 6) { for (int c = 0; c < 6; ++c) acc += dt * s_d2Int_qv[i*36 + b*6 + c] * Vgrad(c, a); }")
    self.gen_add_code_line("// t2: dt^2 * sum_{m,c} d2Int_vv[i,m,c] * Vgrad[c,a] * Vgrad[m,b].")
    self.gen_add_code_line("for (int m = 0; m < 6; ++m) { T vmb = Vgrad(m, b); if (vmb != static_cast<T>(0)) for (int c = 0; c < 6; ++c) acc += dt2 * s_d2Int_vv[i*36 + m*6 + c] * Vgrad(c, a) * vmb; }")
    self.gen_add_end_control_flow()
    self.gen_add_code_line("// t3: dt^2 * sum_m dInt_v[i,m] * D2qdd[m,a,b]. dInt_v is block-diagonal")
    self.gen_add_code_line("//     (6x6 free-flyer corner + identity revolute), so the m-sum is the 6")
    self.gen_add_code_line("//     free-flyer rows plus the single identity term m==i for i>=6.")
    self.gen_add_code_line("if (i < 6) { for (int m = 0; m < 6; ++m) acc += dt2 * dIntv(i, m) * D2qdd(m, a, b); }")
    self.gen_add_code_line("else { acc += dt2 * D2qdd(i, a, b); }")
    self.gen_add_code_line("val = acc;")
    self.gen_add_end_control_flow()
    self.gen_add_code_line("s_d2AB[ind] = val;")
    self.gen_add_end_control_flow()
    self.gen_add_sync()


def gen_integrator_hessian_device_function_call(self,
                                                scratch_in_smem_expr="true",
                                                fd_grad_use_spill_expr="false",
                                                contract_in_smem_expr="true",
                                                d_workspace_pool_name="nullptr",
                                                d_fd_grad_spill_name="nullptr",
                                                s_fdsva_temp_name="nullptr"):
    """Emit the call to `integrator_hessian_device`. Arg order MUST match the def
    in gen_integrator_hessian_device. The fdsva_so spill/pool regions default to
    nullptr (unused under the SHARED tier's all-smem placement)."""
    tmpl = ("<T, IT, " + scratch_in_smem_expr + ", " + fd_grad_use_spill_expr
            + ", " + contract_in_smem_expr + ">")
    start = ("integrator_hessian_device" + tmpl
             + "(s_d2AB, s_df2, s_idsva_so, s_Minv, s_df_du, s_qdd, s_q, s_qd, s_u, ")
    middle = self.gen_insert_helpers_function_call()
    end = ("s_temp, " + d_workspace_pool_name + ", " + d_fd_grad_spill_name + ", "
           + s_fdsva_temp_name + ", d_robotModel, gravity, dt);")
    self.gen_add_code_line(start + middle + end)


def gen_integrator_hessian_device(self):
    """Emit `integrator_hessian_device` — the second-order sensitivity of the
    integrator step x_{k+1} = [q; v] as ONE inner that composes
    `fdsva_so_device` (producing the four 2nd-order forward-dynamics blocks in
    s_df2) and assembles the s_d2AB surface (the plant_step_hessian output).

    Output s_d2AB has shape (2*NUM_VEL, 3*NUM_VEL, 3*NUM_VEL) in row-major
    (C-order) flat layout: H[o*nz*nz + a*nz + b] = d^2 x_{k+1}[o] / dz[a] dz[b],
    z = [dq(nv); dqd(nv); du(nv)], output rows = [position-tangent(nv); velocity(nv)].

    Assembly (mirrors RBDReference._PlantMixin._d2qdd_tangent + plant_step_hessian):
    the full forward-dynamics Hessian D2[i,a,b] = d^2 qdd[i]/dz[a]dz[b] is built
    from the fdsva_so blocks in s_df2 = [d2a_dqdq | d2a_dvdq | d2a_dvdv | d2a_dtdq]
    (each nv^3, laid out [i*nv*nv + j*nv + k]) with the z-block structure
        [q,q]=d2a_dqdq; [qd,qd]=d2a_dvdv; [qd,q]=d2a_dvdq, [q,qd]=d2a_dvdq^T(jk);
        [u,q]=d2a_dtdq, [q,u]=d2a_dtdq^T(jk); all u-u / u-qd / qd-u blocks = 0.
    Velocity rows (bottom nv) = dt*D2 (both Euler and SI-Euler). Position rows
    (top nv) = 0 (Euler) or dt*dt*D2 (SI-Euler, q_{k+1}=q+dt*v_{k+1}).

    Scope: EULER + SEMI_IMPLICIT_EULER, fixed-base. Multi-stage RK (the 2nd-order
    chain rule) and floating-base (the SE(3) retract Hessian) static_assert out
    (clean-break: no silently-wrong tensor). SCRATCH_IN_SMEM=true is the SHARED
    (full-smem PERF) tier; the fdsva_so spill flags are threaded through for
    later tier work but default to the all-smem placement."""
    n = self.robot.get_num_vel()
    nz = 3 * n
    func_params = [
        "s_d2AB is the output Hessian (2*NUM_VEL x 3*NUM_VEL x 3*NUM_VEL, row-major); size " + str(2 * n * nz * nz),
        "s_df2/s_idsva_so/s_Minv/s_df_du/s_qdd are fdsva_so in/out scratch (caller places)",
        "s_q/s_qd/s_u are the joint positions, velocities, and input torques",
        "s_temp is the fdsva_so shared scratch pool (used when SCRATCH_IN_SMEM)",
        "d_workspace/d_fd_grad_spill/s_fdsva_temp are fdsva_so spill regions (SHARED tier: nullptr)",
        "d_robotModel holds XImats/topology; gravity is the gravity constant; dt is the timestep",
    ]
    func_def_start = "void integrator_hessian_device(T *s_d2AB, "
    func_def_middle = ("T *s_df2, T *s_idsva_so, T *s_Minv, T *s_df_du, T *s_qdd, "
                       "const T *s_q, const T *s_qd, const T *s_u, ")
    func_def_end = ("T *s_temp, T *d_workspace, T *d_fd_grad_spill, T *s_fdsva_temp, "
                    "const robotModel<T> *d_robotModel, const T gravity, const T dt) {")
    func_def_middle, func_params = self.gen_insert_helpers_func_def_params(func_def_middle, func_params, -2)
    func_def = func_def_start + func_def_middle + func_def_end
    self.gen_add_func_doc("integrator hessian (plant_step_hessian s_d2AB surface): composes fdsva_so + dt-scaled assembly",
                          [], func_params, None)
    self.gen_add_code_line("template <typename T, IntegratorType IT = IntegratorType::EULER, "
                           "bool SCRATCH_IN_SMEM = true, bool FD_GRAD_USE_SPILL = false, bool CONTRACT_IN_SMEM = true>")
    self.gen_add_code_line("__device__ __forceinline__")
    self.gen_add_code_line(func_def, True)
    # Clean-break deferrals: only single-stage Euler / SI-Euler on a fixed base.
    self.gen_add_code_line("static_assert(IT == IntegratorType::EULER || IT == IntegratorType::SEMI_IMPLICIT_EULER,")
    self.gen_add_code_line("    \"integrator_hessian_device: only EULER / SEMI_IMPLICIT_EULER are supported \"")
    self.gen_add_code_line("    \"(multi-stage RK 2nd-order chain rule is deferred; see f1_plant_step_hessian_plan.md).\");")
    if self.robot.floating_base:
        gen_integrator_hessian_device_floating(self)
        self.gen_add_end_function()
        return
    # The four 2nd-order forward-dynamics blocks (fdsva_so output), each nv^3.
    self.gen_fdsva_so_device_function_call(
        scratch_in_smem_expr="SCRATCH_IN_SMEM",
        fd_grad_use_spill_expr="FD_GRAD_USE_SPILL",
        contract_in_smem_expr="CONTRACT_IN_SMEM",
        d_workspace_pool_name="d_workspace",
        d_fd_grad_spill_name="d_fd_grad_spill",
        s_fdsva_temp_name="s_fdsva_temp")
    self.gen_add_sync()
    self.gen_add_code_line("T *d2a_dqdq = s_df2;")
    self.gen_add_code_line("T *d2a_dvdq = &s_df2[" + str(n * n * n) + "];")
    self.gen_add_code_line("T *d2a_dvdv = &s_df2[" + str(2 * n * n * n) + "];")
    self.gen_add_code_line("T *d2a_dtdq = &s_df2[" + str(3 * n * n * n) + "];")
    # SI-Euler also carries the dt^2*D2 position rows (top nv); Euler leaves them 0.
    self.gen_add_code_line("const bool si = (IT == IntegratorType::SEMI_IMPLICIT_EULER);")
    self.gen_add_code_line("const T dt2 = dt * dt;")
    # Assemble + dt-scale in one fully-parallel sweep over the 2*nv*nz*nz output
    # cells (max in-block parallelism; each cell is an independent scatter). For
    # output cell (o, a, b): o<nv selects a position row (SI-Euler dt^2, Euler 0),
    # o>=nv a velocity row (dt). The (a,b) z-block picks which fdsva_so block (and
    # its i,j,k -> [i*nv*nv + j*nv + k] index) contributes, else 0.
    self.gen_add_parallel_loop("ind", str(2 * n * nz * nz))
    self.gen_add_code_line("int o = ind / " + str(nz * nz) + ";")
    self.gen_add_code_line("int a = (ind / " + str(nz) + ") % " + str(nz) + ";")
    self.gen_add_code_line("int b = ind % " + str(nz) + ";")
    self.gen_add_code_line("int i = o % " + str(n) + ";  // qdd component index for this output row")
    self.gen_add_code_line("T d2 = static_cast<T>(0);")
    # z-block lookup: a,b in {q:[0,nv), qd:[nv,2nv), u:[2nv,3nv)}. The fdsva_so
    # blocks store [d/(velocity|torque), d/q] (j-index first), so the transposed
    # off-diagonal blocks swap which of (a,b) supplies j vs k.
    self.gen_add_code_line("int ablk = a / " + str(n) + "; int bblk = b / " + str(n) + ";")
    self.gen_add_code_line("int aj = a % " + str(n) + "; int bk = b % " + str(n) + ";")
    self.gen_add_code_line("if (ablk == 0 && bblk == 0)      d2 = d2a_dqdq[i*" + str(n * n) + " + aj*" + str(n) + " + bk];  // [q,q]")
    self.gen_add_code_line("else if (ablk == 1 && bblk == 1) d2 = d2a_dvdv[i*" + str(n * n) + " + aj*" + str(n) + " + bk];  // [qd,qd]")
    self.gen_add_code_line("else if (ablk == 1 && bblk == 0) d2 = d2a_dvdq[i*" + str(n * n) + " + aj*" + str(n) + " + bk];  // [qd,q]")
    self.gen_add_code_line("else if (ablk == 0 && bblk == 1) d2 = d2a_dvdq[i*" + str(n * n) + " + bk*" + str(n) + " + aj];  // [q,qd] = [qd,q]^T(jk)")
    self.gen_add_code_line("else if (ablk == 2 && bblk == 0) d2 = d2a_dtdq[i*" + str(n * n) + " + aj*" + str(n) + " + bk];  // [u,q]")
    self.gen_add_code_line("else if (ablk == 0 && bblk == 2) d2 = d2a_dtdq[i*" + str(n * n) + " + bk*" + str(n) + " + aj];  // [q,u] = [u,q]^T(jk)")
    self.gen_add_code_line("// all u-u / u-qd / qd-u blocks are identically zero (qdd linear in u).")
    self.gen_add_code_line("T scale = (o < " + str(n) + ") ? (si ? dt2 : static_cast<T>(0)) : dt;")
    self.gen_add_code_line("s_d2AB[ind] = scale * d2;")
    self.gen_add_end_control_flow()
    self.gen_add_sync()
    self.gen_add_end_function()


def gen_integrator_gradient(self):
    # Canonical _device (orchestrator: owns s_temp placement; called from kernel).
    # One per output kind (gradient-only vs gradient + x_kp1).
    self.gen_integrator_gradient_device(compute_x_kp1=False)
    self.gen_integrator_gradient_device(compute_x_kp1=True)
    # Gradient-only kernel + host
    self.gen_integrator_gradient_kernel(compute_x_kp1=False, single_call_timing=True)
    self.gen_integrator_gradient_kernel(compute_x_kp1=False, single_call_timing=False)
    self.gen_integrator_gradient_host(0, compute_x_kp1=False)
    self.gen_integrator_gradient_host(1, compute_x_kp1=False)
    self.gen_integrator_gradient_host(2, compute_x_kp1=False)
    # Gradient + x_kp1 (both-at-once) kernel + host
    self.gen_integrator_gradient_kernel(compute_x_kp1=True, single_call_timing=True)
    self.gen_integrator_gradient_kernel(compute_x_kp1=True, single_call_timing=False)
    self.gen_integrator_gradient_host(0, compute_x_kp1=True)
    self.gen_integrator_gradient_host(1, compute_x_kp1=True)
    self.gen_integrator_gradient_host(2, compute_x_kp1=True)
