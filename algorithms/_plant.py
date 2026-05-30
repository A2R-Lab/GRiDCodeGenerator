"""Plant namespace codegen (T6).

Emits a SIBLING `namespace grid_plant { ... }` block AFTER the `grid`
namespace closes. Everything here is ADDITIVE: it composes the already-emitted
`grid::` device functions (integrator value/gradient, end-effector pose +
Jacobian) and the vendored GLASS `grid_linalg_*` primitives. ZERO edits are made
to any existing `grid::` emit path.

State convention (matches the integrator + GATO/PDDP references):
    x = [q (NUM_POS); qd (NUM_VEL)]   (fixed-base: NUM_POS == NUM_VEL == n)
    u = control torques (NUM_VEL)
The plant-gradient block order is [A | B] = [dx_{k+1}/dx | dx_{k+1}/du], the
same `s_dAB` surface `grid::integrator_gradient` already produces (2n x 3n,
column-major).

Cost convention (matches PDDP TrajoptCost / GATO trackingcost):
    quadratic cost  = 1/2 * r^T diag(W) r           (the 1/2 scaling)
    gradient        = diag(W) r
    Gauss-Newton hessian (the RATIFIED choice):
        quadratic   -> diag(W)
        ee-position -> J_p^T W J_p
The true analytic 2nd-order hessian (cost curvature folded with the integrator
Hessian) is intentionally NOT emitted — grid has no analytic 2nd-order
integrator and adding one would violate "additive". It is left as a labeled
TODO at each relevant site.

Barriers (log-barrier, per GATO jointBarrier):
    b(x)  = -mu * ( log(x - lower) + log(upper - x) )
    b'(x) = -mu * ( 1/(x - lower) - 1/(upper - x) )
    b''(x)=  mu * ( 1/(x - lower)^2 + 1/(upper - x)^2 )
Bounds are passed in as explicit `s_lower` / `s_upper` pointers (the URDF parser
only carries position limits today, so velocity/torque bounds must be supplied
by the caller). An `isfinite` guard skips any side whose bound is +/-inf, so an
unbounded joint contributes EXACTLY zero to value/gradient/hessian.
"""


# ---------------------------------------------------------------------------
# Plant step (value) + plant step gradient — thin wrappers over the integrator.
# ---------------------------------------------------------------------------

def gen_plant_step(self):
    """`plant_step` — thin wrapper over `grid::integrator_device` (value).

    x_{k+1} = integrator(x_k, u_k, dt). s_x is [q; qd]; we slice q/qd and call
    the integrator's auto-allocating device wrapper (which owns its own scratch).
    """
    nq = self.robot.get_num_pos()
    func_params = [
        "s_x_kp1 is the next state output (size NUM_POS + NUM_VEL)",
        "s_x is the current state [q (NUM_POS); qd (NUM_VEL)]",
        "s_u is the control torque vector (size NUM_VEL)",
        "d_robotModel is the GPU model helpers (XImats, topology, ...)",
        "gravity is the gravity constant",
        "dt is the integration timestep",
    ]
    self.gen_add_func_doc("Plant step: x_{k+1} = integrator(x_k, u_k, dt) (thin wrapper over grid::integrator_device)",
                          [], func_params, None)
    self.gen_add_code_line("template <typename T, grid::IntegratorType IT = grid::IntegratorType::EULER>")
    self.gen_add_code_line("__device__")
    self.gen_add_code_line("void plant_step(T *s_x_kp1, const T *s_x, const T *s_u, "
                           "const grid::robotModel<T> *d_robotModel, const T gravity, const T dt) {", True)
    self.gen_add_code_line("const T *s_q  = s_x;")
    self.gen_add_code_line("const T *s_qd = &s_x[" + str(nq) + "];")
    self.gen_add_code_line("grid::integrator_device<T, IT>(s_x_kp1, s_q, s_qd, s_u, d_robotModel, gravity, dt);")
    self.gen_add_end_function()


def gen_plant_step_gradient(self, with_value=False):
    """`plant_step_gradient[_and_value]` — thin wrapper over
    `grid::integrator_gradient_device` (the [A|B] = s_dAB surface).

    Pass-through: the emitted s_dAB IS grid's integrator gradient (this is what
    the equivalence test asserts). When `with_value`, also returns x_{k+1} via
    the both-at-once integrator gradient device (no extra RBD work).

    The caller supplies the FD-grad scratch buffers the integrator-gradient
    device needs (it is an inner-owns-placement orchestrator); we forward them
    straight through. s_q / s_qd are NON-const because the multi-stage RK path
    mutates them in place across stages (see grid::integrator_gradient_device).

    Plant HESSIAN: per the ratified decision the plant-step hessian used by the
    cost layer is the Gauss-Newton outer product of the COST gradient, assembled
    in the cost hessian functions below. The true second-order integrator
    Hessian (d^2 x_{k+1} / d(x,u)^2, a 2n x 3n x 3n tensor) is NOT emitted here.
    // TODO(plant-2nd-order): grid has no analytic 2nd-order integrator; a true
    // plant Hessian would require one. Left out deliberately (additive-only).
    """
    suffix = "_and_value" if with_value else ""
    fname = "plant_step_gradient" + suffix
    func_params = [
        "s_dAB is the [A | B] output (2*NUM_VEL x 3*NUM_VEL, column-major)",
    ]
    if with_value:
        func_params.append("s_x_kp1 is the next-state output (size NUM_POS + NUM_VEL)")
    func_params += [
        "s_x is the current state [q; qd] (NON-const: mutated across RK stages)",
        "s_u is the control torque vector (size NUM_VEL)",
        "s_df_du / s_dc_du / s_vaf / s_Minv / s_qdd are FD-grad in/out scratch (caller-placed)",
        "s_q_orig / s_qd_orig / s_stage_grad_qdd / s_D_qdd_stage are multi-stage scratch (caller-placed)",
        "s_dInt_q_6x6 / s_dInt_v_6x6 are floating-base SE(3) dIntegrate blocks (unused fixed-base)",
        "s_temp / d_workspace / d_temp_spill are the integrator-gradient scratch arenas (caller-placed)",
        "d_robotModel / gravity / dt as for plant_step",
    ]
    nq = self.robot.get_num_pos()
    self.gen_add_func_doc("Plant step gradient [A|B]" + (" + value" if with_value else "") +
                          " (thin wrapper over grid::integrator_gradient" + ("_with_x_kp1" if with_value else "") +
                          "_device — pass-through)",
                          [], func_params, None)
    self.gen_add_code_line("template <typename T, grid::IntegratorType IT = grid::IntegratorType::EULER, "
                           "bool SCRATCH_IN_SMEM = true, bool USE_DA_DF_SPILL = false>")
    self.gen_add_code_line("__device__")
    sig = "void " + fname + "(T *s_dAB, "
    if with_value:
        sig += "T *s_x_kp1, "
    # The middle params end with the SE(3) dInt blocks; the shared XImats /
    # topology helpers are injected here (same mechanism the integrator-gradient
    # device uses) so this stays a faithful pass-through. The s_temp / workspace
    # / spill arenas + (model, gravity, dt) close the signature.
    sig_middle = ("T *s_x, const T *s_u, T *s_df_du, T *s_dc_du, T *s_vaf, T *s_Minv, T *s_qdd, "
                  "T *s_q_orig, T *s_qd_orig, T *s_stage_grad_qdd, T *s_D_qdd_stage, "
                  "T *s_dInt_q_6x6, T *s_dInt_v_6x6, ")
    sig_middle, func_params = self.gen_insert_helpers_func_def_params(sig_middle, func_params, -1)
    sig_end = ("T *s_temp, T *d_workspace, T *d_temp_spill, "
               "const grid::robotModel<T> *d_robotModel, const T gravity, const T dt) {")
    self.gen_add_code_line(sig + sig_middle + sig_end, True)
    self.gen_add_code_line("T *s_q  = s_x;")
    self.gen_add_code_line("T *s_qd = &s_x[" + str(nq) + "];")
    inner = "grid::integrator_gradient" + ("_with_x_kp1" if with_value else "") + "_device" \
            "<T, IT, SCRATCH_IN_SMEM, USE_DA_DF_SPILL>(s_dAB, "
    if with_value:
        inner += "s_x_kp1, "
    inner_middle = ("s_q, s_qd, s_u, s_df_du, s_dc_du, s_vaf, s_Minv, s_qdd, "
                    "s_q_orig, s_qd_orig, s_stage_grad_qdd, s_D_qdd_stage, "
                    "s_dInt_q_6x6, s_dInt_v_6x6, ")
    inner_helpers = self.gen_insert_helpers_function_call()
    inner_end = ("s_temp, d_workspace, d_temp_spill, d_robotModel, gravity, dt);")
    self.gen_add_code_line(inner + inner_middle + inner_helpers + inner_end)
    self.gen_add_end_function()


# ---------------------------------------------------------------------------
# Quadratic state / input cost (value, gradient, GN-diag hessian, fused).
# ---------------------------------------------------------------------------

def _gen_quadratic_cost_family(self, which):
    """Emit the quadratic cost family for `which` in {"state", "input"}.

    state cost: r = x - x_des (size NX = NUM_POS + NUM_VEL), weights s_Q (NX).
    input cost: r = u - u_des (size NU = NUM_VEL),           weights s_R (NU).
    Cost = 1/2 * sum_i W_i r_i^2. Gradient_i = W_i r_i. GN hessian = diag(W).
    """
    nq = self.robot.get_num_pos()
    nv = self.robot.get_num_vel()
    if which == "state":
        size = nq + nv
        var, des, w = "s_x", "s_x_des", "s_Q"
        base = "quadratic_state_cost"
        size_doc = "NUM_POS + NUM_VEL = " + str(size)
    else:
        size = nv
        var, des, w = "s_u", "s_u_des", "s_R"
        base = "quadratic_input_cost"
        size_doc = "NUM_VEL = " + str(size)
    N = str(size)

    # ---- value: cost = 1/2 sum W_i (var_i - des_i)^2, accumulated into s_out[0] ----
    self.gen_add_func_doc(
        base + ": value = 1/2 * sum_i " + w + "[i] * (" + var + "[i] - " + des + "[i])^2",
        ["Block-cooperative: each thread accumulates its strided terms into s_scratch, then a serial reduction writes s_out[0].",
         "s_scratch must hold at least " + N + " elements."],
        ["s_out is the scalar cost output (s_out[0])",
         var + " is the current value (size " + size_doc + ")",
         des + " is the desired/target value (size " + size_doc + ")",
         w + " is the diagonal weight vector (size " + size_doc + ")",
         "s_scratch is shared scratch of size >= " + N],
        None)
    self.gen_add_code_line("template <typename T>")
    self.gen_add_code_line("__device__")
    self.gen_add_code_line("void " + base + "(T *s_out, const T *" + var + ", const T *" + des +
                           ", const T *" + w + ", T *s_scratch) {", True)
    self.gen_add_parallel_loop("i", N)
    self.gen_add_code_line("T r = " + var + "[i] - " + des + "[i];")
    self.gen_add_code_line("s_scratch[i] = static_cast<T>(0.5) * " + w + "[i] * r * r;")
    self.gen_add_end_control_flow()
    self.gen_add_sync()
    self.gen_add_serial_ops()
    self.gen_add_code_line("T acc = static_cast<T>(0);")
    self.gen_add_code_line("for (int i = 0; i < " + N + "; ++i) acc += s_scratch[i];")
    self.gen_add_code_line("s_out[0] = acc;")
    self.gen_add_end_control_flow()
    self.gen_add_end_function()

    # ---- gradient: g_i = W_i (var_i - des_i), written with `mode` (set or add) ----
    self.gen_add_func_doc(
        base + "_gradient: g[i] = " + w + "[i] * (" + var + "[i] - " + des + "[i])",
        ["ACCUMULATE=false overwrites s_grad; ACCUMULATE=true adds into it (for fusing into a packed [x;u] gradient)."],
        ["s_grad is the gradient output (size " + size_doc + ")",
         var + " / " + des + " / " + w + " as in the value function"],
        None)
    self.gen_add_code_line("template <typename T, bool ACCUMULATE = false>")
    self.gen_add_code_line("__device__")
    self.gen_add_code_line("void " + base + "_gradient(T *s_grad, const T *" + var + ", const T *" + des +
                           ", const T *" + w + ") {", True)
    self.gen_add_parallel_loop("i", N)
    self.gen_add_code_line("T g = " + w + "[i] * (" + var + "[i] - " + des + "[i]);")
    self.gen_add_code_line("if (ACCUMULATE) { s_grad[i] += g; } else { s_grad[i] = g; }")
    self.gen_add_end_control_flow()
    self.gen_add_end_function()

    # ---- GN-diag hessian: H = diag(W), column-major n x n ----
    self.gen_add_func_doc(
        base + "_hessian: Gauss-Newton hessian = diag(" + w + ") (RATIFIED: GN outer product; for a quadratic cost this is exactly diag(W))",
        ["Writes a dense column-major " + N + " x " + N + " matrix; off-diagonal entries are zero.",
         "ACCUMULATE=false overwrites; ACCUMULATE=true adds into the diagonal of an existing block."],
        ["s_hess is the dense hessian output (size " + N + "*" + N + ", column-major)",
         w + " is the diagonal weight vector (size " + size_doc + ")"],
        None)
    self.gen_add_code_line("template <typename T, bool ACCUMULATE = false>")
    self.gen_add_code_line("__device__")
    self.gen_add_code_line("void " + base + "_hessian(T *s_hess, const T *" + w + ") {", True)
    self.gen_add_parallel_loop("ind", str(size * size))
    self.gen_add_code_line("int row = ind % " + N + ";")
    self.gen_add_code_line("int col = ind / " + N + ";")
    self.gen_add_code_line("T h = (row == col) ? " + w + "[row] : static_cast<T>(0);")
    self.gen_add_code_line("if (ACCUMULATE) { s_hess[ind] += h; } else { s_hess[ind] = h; }")
    self.gen_add_end_control_flow()
    self.gen_add_end_function()

    # ---- fused value + grad + hess ----
    self.gen_add_func_doc(
        base + "_value_grad_hess: fused value + gradient + GN-diag hessian in one pass",
        ["Convenience fusion of the three functions above; same conventions and ACCUMULATE semantics for grad/hess."],
        ["s_out / s_grad / s_hess are the three outputs",
         var + " / " + des + " / " + w + " / s_scratch as above"],
        None)
    self.gen_add_code_line("template <typename T, bool ACCUMULATE = false>")
    self.gen_add_code_line("__device__")
    self.gen_add_code_line("void " + base + "_value_grad_hess(T *s_out, T *s_grad, T *s_hess, "
                           "const T *" + var + ", const T *" + des + ", const T *" + w + ", T *s_scratch) {", True)
    self.gen_add_code_line(base + "<T>(s_out, " + var + ", " + des + ", " + w + ", s_scratch);")
    self.gen_add_code_line(base + "_gradient<T, ACCUMULATE>(s_grad, " + var + ", " + des + ", " + w + ");")
    self.gen_add_code_line(base + "_hessian<T, ACCUMULATE>(s_hess, " + w + ");")
    self.gen_add_end_function()


def gen_quadratic_state_cost(self):
    _gen_quadratic_cost_family(self, "state")


def gen_quadratic_input_cost(self):
    _gen_quadratic_cost_family(self, "input")


# ---------------------------------------------------------------------------
# End-effector position cost (value, gradient wrt x=[q;qd], GN hessian J_p^T W J_p).
# ---------------------------------------------------------------------------

def gen_ee_pos_cost(self):
    """ee_pos_cost family. p(q) = grid::end_effector_pose (rows 0..2 of the 6-pose);
    J_p = rows 0..2 of grid::end_effector_pose_gradient (layout
    s_deePos[6*NV*ee + 6*vi + row]). Templated on `int EE = 0`.

        r       = p(q) - p_des                              (3-vector)
        value   = 1/2 * sum_{r} W[r] * r[r]^2
        grad_q  = J_p^T W r        (size NV)
        grad_x  = [grad_q ; 0]     (the qd-block is exactly zero)
        GN hess = J_p^T W J_p      (NV x NV block; the qd rows/cols are zero)

    W is a 3-vector of per-axis position weights. The true second-order term
    (sum_r W[r] r[r] * d^2 p_r/dq^2, i.e. folding the EE Hessian) is dropped per
    the ratified Gauss-Newton decision.
    // TODO(ee-2nd-order): add the W*r weighted EE-Hessian term for a true Newton
    // hessian once/if an analytic d2ee path is wired in here.
    """
    nq = self.robot.get_num_pos()
    nv = self.robot.get_num_vel()
    num_ees = self.robot.get_total_leaf_nodes()
    nx = nq + nv

    # ---- value ----
    self.gen_add_func_doc(
        "ee_pos_cost: value = 1/2 * sum_r W[r] * (p_r(q) - p_des_r)^2 over the 3 position axes",
        ["Calls grid::end_effector_pose_device for p(q) (auto-allocating; owns its scratch).",
         "EE selects which end-effector (0.." + str(num_ees - 1) + ").",
         "s_eePos must hold 6*NUM_EE; s_scratch unused here but kept for signature uniformity."],
        ["s_out is the scalar cost output (s_out[0])",
         "s_q is the joint position vector (size NUM_POS)",
         "s_p_des is the desired EE position (3-vector)",
         "s_W is the per-axis position weight (3-vector)",
         "s_eePos is scratch for the 6*NUM_EE pose (the position is rows 0..2 of EE block)",
         "d_robotModel is the GPU model helpers"],
        None)
    self.gen_add_code_line("template <typename T, int EE = 0>")
    self.gen_add_code_line("__device__")
    self.gen_add_code_line("void ee_pos_cost(T *s_out, const T *s_q, const T *s_p_des, const T *s_W, "
                           "T *s_eePos, const grid::robotModel<T> *d_robotModel) {", True)
    self.gen_add_code_line("grid::end_effector_pose_device<T>(s_eePos, s_q, d_robotModel);")
    self.gen_add_sync()
    self.gen_add_serial_ops()
    self.gen_add_code_line("T acc = static_cast<T>(0);")
    self.gen_add_code_line("#pragma unroll")
    self.gen_add_code_line("for (int r = 0; r < 3; ++r) { T e = s_eePos[6*EE + r] - s_p_des[r]; acc += static_cast<T>(0.5) * s_W[r] * e * e; }")
    self.gen_add_code_line("s_out[0] = acc;")
    self.gen_add_end_control_flow()
    self.gen_add_end_function()

    # ---- gradient wrt x = [q; qd] (qd block is zero) ----
    self.gen_add_func_doc(
        "ee_pos_cost_gradient: grad_x = [J_p^T W (p - p_des) ; 0], over x = [q; qd]",
        ["Calls grid::end_effector_pose_device (for p) and grid::end_effector_pose_gradient_device (for J_p).",
         "J_p = rows 0..2 of s_deePos, layout s_deePos[6*NUM_VEL*ee + 6*vi + row].",
         "The qd-block of the gradient (entries NUM_VEL.." + str(nx - 1) + ") is set to exactly zero.",
         "ACCUMULATE=false overwrites s_grad; true adds (for fusing with a state-cost gradient)."],
        ["s_grad is the gradient over x (size NUM_POS + NUM_VEL = " + str(nx) + ")",
         "s_q / s_p_des / s_W / d_robotModel as above",
         "s_eePos is 6*NUM_EE pose scratch; s_deePos is 6*NUM_VEL*NUM_EE Jacobian scratch"],
        None)
    self.gen_add_code_line("template <typename T, int EE = 0, bool ACCUMULATE = false>")
    self.gen_add_code_line("__device__")
    self.gen_add_code_line("void ee_pos_cost_gradient(T *s_grad, const T *s_q, const T *s_p_des, const T *s_W, "
                           "T *s_eePos, T *s_deePos, const grid::robotModel<T> *d_robotModel) {", True)
    self.gen_add_code_line("grid::end_effector_pose_device<T>(s_eePos, s_q, d_robotModel);")
    self.gen_add_code_line("grid::end_effector_pose_gradient_device<T>(s_deePos, s_q, d_robotModel);")
    self.gen_add_sync()
    # grad_q[i] = sum_r J_p[r,i] * W[r] * (p_r - p_des_r)
    self.gen_add_parallel_loop("i", str(nv))
    self.gen_add_code_line("T g = static_cast<T>(0);")
    self.gen_add_code_line("#pragma unroll")
    self.gen_add_code_line("for (int r = 0; r < 3; ++r) {")
    self.gen_add_code_line("    T Jri = s_deePos[6*" + str(nv) + "*EE + 6*i + r];")
    self.gen_add_code_line("    T e   = s_eePos[6*EE + r] - s_p_des[r];")
    self.gen_add_code_line("    g += Jri * s_W[r] * e;")
    self.gen_add_code_line("}")
    self.gen_add_code_line("if (ACCUMULATE) { s_grad[i] += g; } else { s_grad[i] = g; }")
    self.gen_add_end_control_flow()
    # qd-block is exactly zero (only meaningful for the non-accumulate path).
    self.gen_add_code_line("if (!ACCUMULATE) {", True)
    self.gen_add_parallel_loop("i", str(nv))
    self.gen_add_code_line("s_grad[" + str(nq) + " + i] = static_cast<T>(0);")
    self.gen_add_end_control_flow()
    self.gen_add_end_control_flow()
    self.gen_add_end_function()

    # ---- GN hessian J_p^T W J_p over the q-block of x (qd rows/cols zero) ----
    self.gen_add_func_doc(
        "ee_pos_cost_hessian: Gauss-Newton hessian = J_p^T W J_p in the q-block of the x-hessian",
        ["RATIFIED GN choice: H = J_p^T diag(W) J_p (the W*r weighted EE-Hessian term is dropped).",
         "Dense column-major NX x NX (NX = NUM_POS + NUM_VEL = " + str(nx) + "); only the top-left NUM_VEL x NUM_VEL q-block is non-zero.",
         "ACCUMULATE=false overwrites the whole NX x NX block; true adds the q-block into an existing hessian."],
        ["s_hess is the dense x-hessian output (size " + str(nx) + "*" + str(nx) + ", column-major)",
         "s_q / s_W / d_robotModel as above; s_deePos is 6*NUM_VEL*NUM_EE Jacobian scratch"],
        None)
    self.gen_add_code_line("template <typename T, int EE = 0, bool ACCUMULATE = false>")
    self.gen_add_code_line("__device__")
    self.gen_add_code_line("void ee_pos_cost_hessian(T *s_hess, const T *s_q, const T *s_W, "
                           "T *s_deePos, const grid::robotModel<T> *d_robotModel) {", True)
    self.gen_add_code_line("grid::end_effector_pose_gradient_device<T>(s_deePos, s_q, d_robotModel);")
    self.gen_add_sync()
    # H[i,j] = sum_r J_p[r,i] * W[r] * J_p[r,j], column-major over the full NX x NX
    # block (zero outside the NUM_VEL x NUM_VEL q-block).
    self.gen_add_parallel_loop("ind", str(nx * nx))
    self.gen_add_code_line("int row = ind % " + str(nx) + ";")
    self.gen_add_code_line("int col = ind / " + str(nx) + ";")
    self.gen_add_code_line("T h = static_cast<T>(0);")
    self.gen_add_code_line("if (row < " + str(nv) + " && col < " + str(nv) + ") {")
    self.gen_add_code_line("    #pragma unroll")
    self.gen_add_code_line("    for (int r = 0; r < 3; ++r) {")
    self.gen_add_code_line("        T Jri = s_deePos[6*" + str(nv) + "*EE + 6*row + r];")
    self.gen_add_code_line("        T Jrj = s_deePos[6*" + str(nv) + "*EE + 6*col + r];")
    self.gen_add_code_line("        h += Jri * s_W[r] * Jrj;")
    self.gen_add_code_line("    }")
    self.gen_add_code_line("}")
    self.gen_add_code_line("if (ACCUMULATE) { s_hess[ind] += h; } else { s_hess[ind] = h; }")
    self.gen_add_end_control_flow()
    self.gen_add_end_function()


# ---------------------------------------------------------------------------
# Log-barriers (joint position / velocity / torque). Explicit bound pointers.
# ---------------------------------------------------------------------------

def _gen_barrier_helpers(self):
    """Emit the scalar one-sided-safe log-barrier helpers (value/grad/hess).

    Each `isfinite`-guards both sides so a +/-inf bound contributes zero. A tiny
    margin floor keeps the log/reciprocal finite at (and just outside) the
    boundary without changing the interior value materially.
    """
    self.gen_add_func_doc(
        "Scalar log-barrier helpers: b = -mu*(log(x-lo)+log(hi-x)); isfinite-guarded so an inf bound contributes zero",
        ["Margin floored at 1e-10 (value) / 1e-6 (grad,hess) to stay finite at the boundary."])
    self.gen_add_code_lines([
        "template <typename T> __device__ inline T grid_plant_log_barrier(T x, T lo, T hi, T mu) {",
        "    T b = static_cast<T>(0);",
        "    if (isfinite(lo)) { T d = x - lo; d = (d <= static_cast<T>(1e-10)) ? static_cast<T>(1e-10) : d; b -= log(d); }",
        "    if (isfinite(hi)) { T d = hi - x; d = (d <= static_cast<T>(1e-10)) ? static_cast<T>(1e-10) : d; b -= log(d); }",
        "    return mu * b;",
        "}",
        "",
        "template <typename T> __device__ inline T grid_plant_log_barrier_grad(T x, T lo, T hi, T mu) {",
        "    T g = static_cast<T>(0);",
        "    const T eps = static_cast<T>(1e-6);",
        "    if (isfinite(lo)) { T d = x - lo; T a = (d < static_cast<T>(0)) ? -d : d; if (a < eps) a = eps; d = (d < static_cast<T>(0)) ? -a : a; g -= static_cast<T>(1) / d; }",
        "    if (isfinite(hi)) { T d = hi - x; T a = (d < static_cast<T>(0)) ? -d : d; if (a < eps) a = eps; d = (d < static_cast<T>(0)) ? -a : a; g += static_cast<T>(1) / d; }",
        "    return mu * g;",
        "}",
        "",
        "template <typename T> __device__ inline T grid_plant_log_barrier_hess(T x, T lo, T hi, T mu) {",
        "    T h = static_cast<T>(0);",
        "    const T eps = static_cast<T>(1e-6);",
        "    if (isfinite(lo)) { T d = x - lo; T a = (d < static_cast<T>(0)) ? -d : d; if (a < eps) a = eps; h += static_cast<T>(1) / (a * a); }",
        "    if (isfinite(hi)) { T d = hi - x; T a = (d < static_cast<T>(0)) ? -d : d; if (a < eps) a = eps; h += static_cast<T>(1) / (a * a); }",
        "    return mu * h;",
        "}",
        "",
    ])


def _gen_one_barrier(self, base, doc_target, count, slice_offset):
    """Emit `<base>` (value), `<base>_gradient`, `<base>_hessian` for `count`
    bounded DOFs reading `s_var[slice_offset + i]`, writing grad/hess into the
    matching slice of a packed buffer (so torque/velocity barriers add into the
    right rows of a [x;u] gradient/hessian). All barriers ADD (+=) into outputs.

    `slice_offset` is the row/col offset of this barrier's DOFs inside the
    packed gradient/hessian (0 for a standalone buffer; NUM_POS for the qd block
    of an x-gradient; etc.). The hessian is written into a dense column-major
    `stride x stride` block — caller passes the block leading dimension.
    """
    N = str(count)
    self.gen_add_func_doc(
        base + ": value += -mu * sum_i ( log(x_i - lo_i) + log(hi_i - x_i) ) over " + doc_target,
        ["Block-cooperative reduction into s_scratch then a serial add into s_out[0].",
         "s_lower / s_upper are explicit bound vectors (size " + N + "); an inf entry skips that side (isfinite guard).",
         "s_scratch must hold at least " + N + " elements."],
        ["s_out is the scalar barrier cost (added into s_out[0])",
         "s_var is the variable vector (this barrier reads s_var[" + str(slice_offset) + " + i])",
         "s_lower / s_upper are the per-DOF bounds (size " + N + ")",
         "mu is the barrier weight",
         "s_scratch is shared scratch of size >= " + N],
        None)
    self.gen_add_code_line("template <typename T>")
    self.gen_add_code_line("__device__")
    self.gen_add_code_line("void " + base + "(T *s_out, const T *s_var, const T *s_lower, const T *s_upper, "
                           "const T mu, T *s_scratch) {", True)
    self.gen_add_parallel_loop("i", N)
    self.gen_add_code_line("s_scratch[i] = grid_plant_log_barrier<T>(s_var[" + str(slice_offset) + " + i], s_lower[i], s_upper[i], mu);")
    self.gen_add_end_control_flow()
    self.gen_add_sync()
    self.gen_add_serial_ops()
    self.gen_add_code_line("T acc = static_cast<T>(0);")
    self.gen_add_code_line("for (int i = 0; i < " + N + "; ++i) acc += s_scratch[i];")
    self.gen_add_code_line("s_out[0] += acc;")
    self.gen_add_end_control_flow()
    self.gen_add_end_function()

    # gradient (always adds into the right slice)
    self.gen_add_func_doc(
        base + "_gradient: s_grad[GRAD_OFFSET + i] += -mu*(1/(x_i-lo_i) - 1/(hi_i-x_i))",
        ["Adds into the packed gradient at GRAD_OFFSET (a template arg so velocity/torque barriers land in the qd / u rows).",
         "isfinite-guarded per side; an unbounded DOF adds exactly zero."],
        ["s_grad is the packed gradient output (added into)",
         "s_var / s_lower / s_upper / mu as in the value function"],
        None)
    self.gen_add_code_line("template <typename T, int VAR_OFFSET = " + str(slice_offset) + ", int GRAD_OFFSET = 0>")
    self.gen_add_code_line("__device__")
    self.gen_add_code_line("void " + base + "_gradient(T *s_grad, const T *s_var, const T *s_lower, const T *s_upper, const T mu) {", True)
    self.gen_add_parallel_loop("i", N)
    self.gen_add_code_line("s_grad[GRAD_OFFSET + i] += grid_plant_log_barrier_grad<T>(s_var[VAR_OFFSET + i], s_lower[i], s_upper[i], mu);")
    self.gen_add_end_control_flow()
    self.gen_add_end_function()

    # hessian (adds onto the diagonal of a dense column-major block)
    self.gen_add_func_doc(
        base + "_hessian: s_hess[(HESS_OFFSET+i)*HESS_STRIDE + (HESS_OFFSET+i)] += mu*(1/(x_i-lo_i)^2 + 1/(hi_i-x_i)^2)",
        ["Adds the barrier curvature onto the DIAGONAL of a dense column-major HESS_STRIDE x HESS_STRIDE block.",
         "HESS_OFFSET places it in the qd / u block; HESS_STRIDE is the block leading dimension.",
         "isfinite-guarded per side; an unbounded DOF adds exactly zero."],
        ["s_hess is the dense column-major hessian (added into)",
         "s_var / s_lower / s_upper / mu as above"],
        None)
    self.gen_add_code_line("template <typename T, int HESS_STRIDE, int VAR_OFFSET = " + str(slice_offset) + ", int HESS_OFFSET = 0>")
    self.gen_add_code_line("__device__")
    self.gen_add_code_line("void " + base + "_hessian(T *s_hess, const T *s_var, const T *s_lower, const T *s_upper, const T mu) {", True)
    self.gen_add_parallel_loop("i", N)
    self.gen_add_code_line("int d = HESS_OFFSET + i;")
    self.gen_add_code_line("s_hess[d * HESS_STRIDE + d] += grid_plant_log_barrier_hess<T>(s_var[VAR_OFFSET + i], s_lower[i], s_upper[i], mu);")
    self.gen_add_end_control_flow()
    self.gen_add_end_function()


def gen_plant_barriers(self):
    """Emit the three log-barriers. URDF parses only position limits today, so
    velocity/torque barriers read caller-supplied explicit bound pointers — no
    URDFParser change. Default VAR_OFFSETs place each on the right block of the
    packed state x = [q; qd]:
      joint_position_barrier reads x[0..NUM_POS)           (q block)
      joint_velocity_barrier reads x[NUM_POS..NUM_POS+NUM_VEL)  (qd block)
      joint_torque_barrier   reads u[0..NUM_VEL)           (standalone u buffer)
    """
    nq = self.robot.get_num_pos()
    nv = self.robot.get_num_vel()
    _gen_barrier_helpers(self)
    _gen_one_barrier(self, "joint_position_barrier", "the NUM_POS q-DOFs", nq, 0)
    _gen_one_barrier(self, "joint_velocity_barrier", "the NUM_VEL qd-DOFs", nv, nq)
    _gen_one_barrier(self, "joint_torque_barrier",   "the NUM_VEL u-DOFs",  nv, 0)


# ---------------------------------------------------------------------------
# Top-level emit: open the sibling namespace and gate each sub-emit on deps.
# ---------------------------------------------------------------------------

def gen_grid_plant(self, algorithms):
    """Emit the sibling `namespace grid_plant { ... }` block. Called AFTER the
    `grid` namespace closes. Each sub-emit is gated on the `grid::` deps it
    composes being present in `algorithms`; when a dep is missing we emit a
    `#warning`-style comment instead of an undefined call.
    """
    self.gen_add_code_line("")
    self.gen_add_func_doc("Plant namespace: cost / constraint / plant-step primitives composed over grid::")
    self.gen_add_code_line("namespace " + self.file_namespace + "_plant {", True)

    # Quadratic costs have no grid:: dep — always emit.
    self.gen_quadratic_state_cost()
    self.gen_quadratic_input_cost()

    # Barriers have no grid:: dep — always emit.
    self.gen_plant_barriers()

    # Plant step needs the integrator value.
    if "integrator" in algorithms:
        self.gen_plant_step()
    else:
        self.gen_add_code_line("// [grid_plant] plant_step skipped: requires the 'integrator' algorithm (grid::integrator_device) — not generated.")

    # Plant step gradient needs the integrator gradient.
    if ("integrator_gradient" in algorithms) or ("integrator_with_gradient" in algorithms):
        self.gen_plant_step_gradient(with_value=False)
        self.gen_plant_step_gradient(with_value=True)
    else:
        self.gen_add_code_line("// [grid_plant] plant_step_gradient[_and_value] skipped: requires 'integrator_gradient' (grid::integrator_gradient_device) — not generated.")

    # EE position cost needs both ee_pose and ee_pose_gradient.
    if ("ee_pose" in algorithms) and ("ee_pose_gradient" in algorithms):
        self.gen_ee_pos_cost()
    else:
        self.gen_add_code_line("// [grid_plant] ee_pos_cost skipped: requires both 'ee_pose' and 'ee_pose_gradient' (grid::end_effector_pose[_gradient]_device) — not generated.")

    self.gen_add_end_control_flow()  # close namespace grid_plant
