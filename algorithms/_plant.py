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
                          " (thin wrapper over grid::" + ("integrator_with_gradient" if with_value else "integrator_gradient") +
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
    inner = "grid::" + ("integrator_with_gradient" if with_value else "integrator_gradient") + "_device" \
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
    s_end_effector_pose_gradient[6*NV*ee + 6*vi + row]). Templated on `int EE = 0`.

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
         "s_end_effector_pose must hold 6*NUM_EE; s_scratch unused here but kept for signature uniformity."],
        ["s_out is the scalar cost output (s_out[0])",
         "s_q is the joint position vector (size NUM_POS)",
         "s_p_des is the desired EE position (3-vector)",
         "s_W is the per-axis position weight (3-vector)",
         "s_end_effector_pose is scratch for the 6*NUM_EE pose (the position is rows 0..2 of EE block)",
         "d_robotModel is the GPU model helpers"],
        None)
    self.gen_add_code_line("template <typename T, int EE = 0>")
    self.gen_add_code_line("__device__")
    self.gen_add_code_line("void ee_pos_cost(T *s_out, const T *s_q, const T *s_p_des, const T *s_W, "
                           "T *s_end_effector_pose, const grid::robotModel<T> *d_robotModel) {", True)
    self.gen_add_code_line("grid::end_effector_pose_device<T>(s_end_effector_pose, s_q, d_robotModel);")
    self.gen_add_sync()
    self.gen_add_serial_ops()
    self.gen_add_code_line("T acc = static_cast<T>(0);")
    self.gen_add_code_line("#pragma unroll")
    self.gen_add_code_line("for (int r = 0; r < 3; ++r) { T e = s_end_effector_pose[6*EE + r] - s_p_des[r]; acc += static_cast<T>(0.5) * s_W[r] * e * e; }")
    self.gen_add_code_line("s_out[0] = acc;")
    self.gen_add_end_control_flow()
    self.gen_add_end_function()

    # ---- gradient wrt x = [q; qd] (qd block is zero) ----
    self.gen_add_func_doc(
        "ee_pos_cost_gradient: grad_x = [J_p^T W (p - p_des) ; 0], over x = [q; qd]",
        ["Calls grid::end_effector_pose_device (for p) and grid::end_effector_pose_gradient_device (for J_p).",
         "J_p = rows 0..2 of s_end_effector_pose_gradient, layout s_end_effector_pose_gradient[6*NUM_VEL*ee + 6*vi + row].",
         "The qd-block of the gradient (entries NUM_VEL.." + str(nx - 1) + ") is set to exactly zero.",
         "ACCUMULATE=false overwrites s_grad; true adds (for fusing with a state-cost gradient)."],
        ["s_grad is the gradient over x (size NUM_POS + NUM_VEL = " + str(nx) + ")",
         "s_q / s_p_des / s_W / d_robotModel as above",
         "s_end_effector_pose is 6*NUM_EE pose scratch; s_end_effector_pose_gradient is 6*NUM_VEL*NUM_EE Jacobian scratch"],
        None)
    self.gen_add_code_line("template <typename T, int EE = 0, bool ACCUMULATE = false>")
    self.gen_add_code_line("__device__")
    self.gen_add_code_line("void ee_pos_cost_gradient(T *s_grad, const T *s_q, const T *s_p_des, const T *s_W, "
                           "T *s_end_effector_pose, T *s_end_effector_pose_gradient, const grid::robotModel<T> *d_robotModel) {", True)
    self.gen_add_code_line("grid::end_effector_pose_device<T>(s_end_effector_pose, s_q, d_robotModel);")
    self.gen_add_code_line("grid::end_effector_pose_gradient_device<T>(s_end_effector_pose_gradient, s_q, d_robotModel);")
    self.gen_add_sync()
    # grad_q[i] = sum_r J_p[r,i] * W[r] * (p_r - p_des_r)
    self.gen_add_parallel_loop("i", str(nv))
    self.gen_add_code_line("T g = static_cast<T>(0);")
    self.gen_add_code_line("#pragma unroll")
    self.gen_add_code_line("for (int r = 0; r < 3; ++r) {")
    self.gen_add_code_line("    T Jri = s_end_effector_pose_gradient[6*" + str(nv) + "*EE + 6*i + r];")
    self.gen_add_code_line("    T e   = s_end_effector_pose[6*EE + r] - s_p_des[r];")
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
         "s_q / s_W / d_robotModel as above; s_end_effector_pose_gradient is 6*NUM_VEL*NUM_EE Jacobian scratch"],
        None)
    self.gen_add_code_line("template <typename T, int EE = 0, bool ACCUMULATE = false>")
    self.gen_add_code_line("__device__")
    self.gen_add_code_line("void ee_pos_cost_hessian(T *s_hess, const T *s_q, const T *s_W, "
                           "T *s_end_effector_pose_gradient, const grid::robotModel<T> *d_robotModel) {", True)
    self.gen_add_code_line("grid::end_effector_pose_gradient_device<T>(s_end_effector_pose_gradient, s_q, d_robotModel);")
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
    self.gen_add_code_line("        T Jri = s_end_effector_pose_gradient[6*" + str(nv) + "*EE + 6*row + r];")
    self.gen_add_code_line("        T Jrj = s_end_effector_pose_gradient[6*" + str(nv) + "*EE + 6*col + r];")
    self.gen_add_code_line("        h += Jri * s_W[r] * Jrj;")
    self.gen_add_code_line("    }")
    self.gen_add_code_line("}")
    self.gen_add_code_line("if (ACCUMULATE) { s_hess[ind] += h; } else { s_hess[ind] = h; }")
    self.gen_add_end_control_flow()
    self.gen_add_end_function()


# ---------------------------------------------------------------------------
# CoM-tracking cost (R3 plant hook). Direct clone of ee_pos_cost with the EE
# position/Jacobian replaced by the CoM position / CoM Jacobian from
# grid::com_device (which writes [p_com(3); J_com(3 x NUM_VEL, column-major)]).
#   r       = p_com(q) - p_des                          (3-vector)
#   value   = 1/2 sum_r W[r] r[r]^2
#   grad_x  = [J_com^T W r ; 0]
#   GN hess = J_com^T W J_com  (q-block of the NX x NX hessian; qd rows/cols 0)
# ---------------------------------------------------------------------------

def gen_com_cost(self):
    nq = self.robot.get_num_pos()
    nv = self.robot.get_num_vel()
    nx = nq + nv
    com_out = 3 + 3 * nv  # com (3) + Jcom (3 x nv) layout from grid::com_device
    # ---- value ----
    self.gen_add_func_doc(
        "com_cost: value = 1/2 sum_r W[r] (p_com_r(q) - p_des_r)^2 over the 3 CoM axes",
        ["Calls grid::com_device for [p_com; J_com] (auto-allocating; owns its scratch).",
         "s_com scratch must hold 3 + 3*NUM_VEL (the com device output)."],
        ["s_out scalar cost", "s_q joint positions", "s_p_des desired CoM (3)",
         "s_W per-axis weight (3)", "s_com scratch (3 + 3*NUM_VEL)", "d_robotModel GPU model helpers"],
        None)
    self.gen_add_code_line("template <typename T>")
    self.gen_add_code_line("__device__")
    self.gen_add_code_line("void com_cost(T *s_out, const T *s_q, const T *s_p_des, const T *s_W, "
                           "T *s_com, const grid::robotModel<T> *d_robotModel) {", True)
    self.gen_add_code_line("grid::com_device<T>(s_com, s_q, d_robotModel);")
    self.gen_add_sync()
    self.gen_add_serial_ops()
    self.gen_add_code_line("T acc = static_cast<T>(0);")
    self.gen_add_code_line("#pragma unroll")
    self.gen_add_code_line("for (int r = 0; r < 3; ++r) { T e = s_com[r] - s_p_des[r]; acc += static_cast<T>(0.5) * s_W[r] * e * e; }")
    self.gen_add_code_line("s_out[0] = acc;")
    self.gen_add_end_control_flow()
    self.gen_add_end_function()

    # ---- gradient wrt x = [q; qd] (qd block zero) ----
    self.gen_add_func_doc(
        "com_cost_gradient: grad_x = [J_com^T W (p_com - p_des) ; 0]",
        ["J_com = rows of s_com starting at offset 3, layout s_com[3 + 3*vi + r] (3 x NUM_VEL column-major)."],
        ["s_grad gradient over x (" + str(nx) + ")", "s_q / s_p_des / s_W / d_robotModel as above",
         "s_com scratch (3 + 3*NUM_VEL)"], None)
    self.gen_add_code_line("template <typename T, bool ACCUMULATE = false>")
    self.gen_add_code_line("__device__")
    self.gen_add_code_line("void com_cost_gradient(T *s_grad, const T *s_q, const T *s_p_des, const T *s_W, "
                           "T *s_com, const grid::robotModel<T> *d_robotModel) {", True)
    self.gen_add_code_line("grid::com_device<T>(s_com, s_q, d_robotModel);")
    self.gen_add_sync()
    self.gen_add_parallel_loop("i", str(nv))
    self.gen_add_code_line("T g = static_cast<T>(0);")
    self.gen_add_code_line("#pragma unroll")
    self.gen_add_code_line("for (int r = 0; r < 3; ++r) { T Jri = s_com[3 + 3*i + r]; T e = s_com[r] - s_p_des[r]; g += Jri * s_W[r] * e; }")
    self.gen_add_code_line("if (ACCUMULATE) { s_grad[i] += g; } else { s_grad[i] = g; }")
    self.gen_add_end_control_flow()
    self.gen_add_code_line("if (!ACCUMULATE) {", True)
    # Zero the entire non-q-gradient tail [nv, nx): the meaningful gradient occupies the
    # first nv slots [0, nv), everything after must be zero (matches the GN hessian
    # convention below, nonzero only on [0,nv)x[0,nv)). Zeroing [nq, nq+nv) left [nv, nq)
    # UNINITIALIZED for floating-base robots (nq>nv) -> stale shared mem (go2 nq=19,nv=18
    # left s_grad[18] stale). For fixed-base (nq==nv) this is byte-identical to the old loop.
    self.gen_add_parallel_loop("i", str(nq))
    self.gen_add_code_line("s_grad[" + str(nv) + " + i] = static_cast<T>(0);")
    self.gen_add_end_control_flow()
    self.gen_add_end_control_flow()
    self.gen_add_end_function()

    # ---- GN hessian J_com^T W J_com over the q-block of x ----
    self.gen_add_func_doc(
        "com_cost_hessian: Gauss-Newton hessian = J_com^T diag(W) J_com in the q-block of the x-hessian",
        ["Dense column-major NX x NX; only the top-left NUM_VEL x NUM_VEL q-block is non-zero."],
        ["s_hess dense x-hessian (" + str(nx*nx) + ")", "s_q / s_W / d_robotModel as above",
         "s_com scratch (3 + 3*NUM_VEL)"], None)
    self.gen_add_code_line("template <typename T, bool ACCUMULATE = false>")
    self.gen_add_code_line("__device__")
    self.gen_add_code_line("void com_cost_hessian(T *s_hess, const T *s_q, const T *s_W, "
                           "T *s_com, const grid::robotModel<T> *d_robotModel) {", True)
    self.gen_add_code_line("grid::com_device<T>(s_com, s_q, d_robotModel);")
    self.gen_add_sync()
    self.gen_add_parallel_loop("ind", str(nx * nx))
    self.gen_add_code_line("int row = ind % " + str(nx) + "; int col = ind / " + str(nx) + ";")
    self.gen_add_code_line("T h = static_cast<T>(0);")
    self.gen_add_code_line("if (row < " + str(nv) + " && col < " + str(nv) + ") {")
    self.gen_add_code_line("    #pragma unroll")
    self.gen_add_code_line("    for (int r = 0; r < 3; ++r) { T Jri = s_com[3 + 3*row + r]; T Jrj = s_com[3 + 3*col + r]; h += Jri * s_W[r] * Jrj; }")
    self.gen_add_code_line("}")
    self.gen_add_code_line("if (ACCUMULATE) { s_hess[ind] += h; } else { s_hess[ind] = h; }")
    self.gen_add_end_control_flow()
    self.gen_add_end_function()


# ---------------------------------------------------------------------------
# Centroidal-momentum-tracking cost (R2 plant hook). The momentum h = A(q) qd
# has velocity-Jacobian J_h = A (the CMM), so this is the same value/grad/GN-hess
# pattern with p->h, J->A but the variable is x=[q;qd] and h depends on qd
# linearly: grad_qd = A^T W r, GN-hess qd-block = A^T W A (the q-dependence of A
# is dropped, Gauss-Newton style, matching the ee_pos_cost ratified choice).
# grid::ccrba_device writes [A (6 x NUM_VEL, column-major); h (6)].
# ---------------------------------------------------------------------------

def gen_momentum_cost(self):
    nq = self.robot.get_num_pos()
    nv = self.robot.get_num_vel()
    nx = nq + nv
    self.gen_add_func_doc(
        "momentum_cost: value = 1/2 sum_r W[r] (h_r(q,qd) - h_des_r)^2 over the 6 centroidal components",
        ["Calls grid::ccrba_device for [A; h] (auto-allocating; owns its scratch).",
         "s_ccrba scratch must hold 6*NUM_VEL + 6 (the ccrba device output)."],
        ["s_out scalar cost", "s_q / s_qd joint position/velocity", "s_h_des desired momentum (6)",
         "s_W per-component weight (6)", "s_ccrba scratch (6*NUM_VEL + 6)", "d_robotModel GPU model helpers"],
        None)
    self.gen_add_code_line("template <typename T>")
    self.gen_add_code_line("__device__")
    self.gen_add_code_line("void momentum_cost(T *s_out, const T *s_q, const T *s_qd, const T *s_h_des, const T *s_W, "
                           "T *s_ccrba, const grid::robotModel<T> *d_robotModel) {", True)
    self.gen_add_code_line("grid::ccrba_device<T>(s_ccrba, s_q, s_qd, d_robotModel);")
    self.gen_add_sync()
    self.gen_add_serial_ops()
    self.gen_add_code_line("T acc = static_cast<T>(0);")
    self.gen_add_code_line("#pragma unroll")
    self.gen_add_code_line("for (int r = 0; r < 6; ++r) { T e = s_ccrba[" + str(6*nv) + " + r] - s_h_des[r]; acc += static_cast<T>(0.5) * s_W[r] * e * e; }")
    self.gen_add_code_line("s_out[0] = acc;")
    self.gen_add_end_control_flow()
    self.gen_add_end_function()

    # gradient wrt x = [q; qd]; the q-block is dropped (GN on A), qd-block = A^T W r
    self.gen_add_func_doc(
        "momentum_cost_gradient: grad_x = [0 ; A^T W (h - h_des)] (q-block dropped, GN on A)",
        ["A = s_ccrba[r + 6*vi] (6 x NUM_VEL column-major); h = s_ccrba[6*NUM_VEL + r]."],
        ["s_grad gradient over x (" + str(nx) + ")", "s_q / s_qd / s_h_des / s_W / d_robotModel as above",
         "s_ccrba scratch (6*NUM_VEL + 6)"], None)
    self.gen_add_code_line("template <typename T, bool ACCUMULATE = false>")
    self.gen_add_code_line("__device__")
    self.gen_add_code_line("void momentum_cost_gradient(T *s_grad, const T *s_q, const T *s_qd, const T *s_h_des, const T *s_W, "
                           "T *s_ccrba, const grid::robotModel<T> *d_robotModel) {", True)
    self.gen_add_code_line("grid::ccrba_device<T>(s_ccrba, s_q, s_qd, d_robotModel);")
    self.gen_add_sync()
    self.gen_add_code_line("if (!ACCUMULATE) {", True)
    self.gen_add_parallel_loop("i", str(nq))
    self.gen_add_code_line("s_grad[i] = static_cast<T>(0);")
    self.gen_add_end_control_flow()
    self.gen_add_end_control_flow()
    self.gen_add_parallel_loop("i", str(nv))
    self.gen_add_code_line("T g = static_cast<T>(0);")
    self.gen_add_code_line("#pragma unroll")
    self.gen_add_code_line("for (int r = 0; r < 6; ++r) { T Ari = s_ccrba[r + 6*i]; T e = s_ccrba[" + str(6*nv) + " + r] - s_h_des[r]; g += Ari * s_W[r] * e; }")
    self.gen_add_code_line("if (ACCUMULATE) { s_grad[" + str(nq) + " + i] += g; } else { s_grad[" + str(nq) + " + i] = g; }")
    self.gen_add_end_control_flow()
    self.gen_add_end_function()

    # GN hessian A^T W A in the qd-block of the NX x NX hessian
    self.gen_add_func_doc(
        "momentum_cost_hessian: Gauss-Newton hessian = A^T diag(W) A in the qd-block of the x-hessian",
        ["Dense column-major NX x NX; only the bottom-right NUM_VEL x NUM_VEL qd-block is non-zero."],
        ["s_hess dense x-hessian (" + str(nx*nx) + ")", "s_q / s_qd / s_W / d_robotModel as above",
         "s_ccrba scratch (6*NUM_VEL + 6)"], None)
    self.gen_add_code_line("template <typename T, bool ACCUMULATE = false>")
    self.gen_add_code_line("__device__")
    self.gen_add_code_line("void momentum_cost_hessian(T *s_hess, const T *s_q, const T *s_qd, const T *s_W, "
                           "T *s_ccrba, const grid::robotModel<T> *d_robotModel) {", True)
    self.gen_add_code_line("grid::ccrba_device<T>(s_ccrba, s_q, s_qd, d_robotModel);")
    self.gen_add_sync()
    self.gen_add_parallel_loop("ind", str(nx * nx))
    self.gen_add_code_line("int row = ind % " + str(nx) + "; int col = ind / " + str(nx) + ";")
    self.gen_add_code_line("T h = static_cast<T>(0);")
    self.gen_add_code_line("if (row >= " + str(nq) + " && col >= " + str(nq) + ") {")
    self.gen_add_code_line("    int vi = row - " + str(nq) + "; int vj = col - " + str(nq) + ";")
    self.gen_add_code_line("    #pragma unroll")
    self.gen_add_code_line("    for (int r = 0; r < 6; ++r) { T Ari = s_ccrba[r + 6*vi]; T Arj = s_ccrba[r + 6*vj]; h += Ari * s_W[r] * Arj; }")
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
# Kernels + host wrappers (the binding layer: G1 grid_plant Python surface).
#
# Each plant kernel runs ONE BLOCK PER TIMESTEP (block-level grid-stride loop)
# and calls the auto-allocating `grid_plant::`/`grid::` device functions, which
# own the WHOLE `extern __shared__` arena. The kernel therefore passes GLOBAL
# device pointers for state / input / output straight through to the device
# function — the device fn reads/writes through whatever pointer it is given,
# so the heavy RBD scratch stays internal and we never collide with it.
#
# Reduction-style cost / barrier device functions additionally take an
# `s_scratch` shared buffer; the kernel declares a tiny `__shared__` array for
# it (these kernels do no RBD arena work, so a static shared array is fine).
#
# The host wrappers take raw device pointers for the plant-specific in/out
# buffers (desired states, weights, bounds, scalar outputs). The grid_rbd C ABI
# (wrapper_template.cu) allocates those device buffers and stages H<->D copies.
# ---------------------------------------------------------------------------

def gen_plant_step_kernel(self):
    """`plant_step_kernel` — one block/timestep, calls grid_plant::plant_step.

    Inputs/outputs are global; the device fn (-> grid::integrator_device) owns
    the shared arena. Reuses grid::INTEGRATOR_DYNAMIC_SHARED_MEM_BYTES.
    """
    nx = self.robot.get_num_pos() + self.robot.get_num_vel()
    nv = self.robot.get_num_vel()
    self.gen_add_func_doc("plant_step kernel: x_{k+1} = integrator(x_k, u_k, dt) per timestep",
                          [],
                          ["d_x_kp1 is the next-state output (NUM_POS+NUM_VEL per timestep)",
                           "d_x is the packed current state [q; qd] (NUM_POS+NUM_VEL per timestep)",
                           "d_u is the packed control torque (NUM_VEL per timestep)",
                           "stride_x / stride_u are the per-timestep strides",
                           "d_robotModel / gravity / dt as for plant_step",
                           "NUM_TIMESTEPS is the batch size"], None)
    self.gen_add_code_line("template <typename T, grid::IntegratorType IT = grid::IntegratorType::EULER>")
    self.gen_add_code_line("__global__")
    self.gen_add_code_line("void plant_step_kernel(T *d_x_kp1, const T *d_x, const T *d_u, "
                           "const int stride_x, const int stride_u, "
                           "const grid::robotModel<T> *d_robotModel, const T gravity, const T dt, const int NUM_TIMESTEPS) {", True)
    self.gen_add_parallel_loop("k", "NUM_TIMESTEPS", block_level=True)
    self.gen_add_code_line("plant_step<T, IT>(&d_x_kp1[k*" + str(nx) + "], &d_x[k*stride_x], &d_u[k*stride_u], d_robotModel, gravity, dt);")
    self.gen_add_sync()
    self.gen_add_end_control_flow()
    self.gen_add_end_function()


def gen_plant_step_gradient_kernel(self):
    """`plant_step_gradient_kernel` — one block/timestep [A|B] = s_dAB.

    plant_step_gradient is an inner-owns-placement orchestrator (it forwards to
    grid::integrator_gradient_device, which OWNS its FD-grad s_temp pool
    placement). The kernel therefore sets up the WHOLE scratch arena in shared
    memory itself — mirroring grid::integrator_gradient_kernel's PERF/full-smem
    body (every band in smem: SCRATCH_IN_SMEM=true, no workspace/spill) — then
    calls grid_plant::plant_step_gradient (the thin pass-through). The emitted
    s_dAB is byte-identical to grid::integrator_gradient's. Inputs/outputs are
    global; the launch reserves grid::INTEGRATOR_DU_DYNAMIC_SHARED_MEM_BYTES.

    x = [q (NUM_POS); qd (NUM_VEL)] (NON-const: the multi-stage RK path mutates
    q/qd across stages in smem); u = control torque (NUM_VEL)."""
    n = self.robot.get_num_vel()
    fb = 1 if self.robot.floating_base else 0
    nx = self.robot.get_num_pos() + self.robot.get_num_vel()
    from ._integrator import _max_stages_in_use
    max_stages = _max_stages_in_use()
    d_qdd_count = max_stages * n * 3 * n
    inner_temp_full = self.gen_integrator_gradient_inner_temp_mem_size()
    # s_vaf is body-indexed (NB bodies, stride 6); size by NB for mimic robots
    # (NB>nv) so the composed FD-grad ID sub-inner's 18*NB writes don't overflow
    # the adjacent buffers (mirrors integrator_gradient_kernel's vaf sizing).
    vaf_cnt = 18 * (self.robot.get_num_joints() if self.robot_has_mimic_joints() else n)
    self.gen_add_func_doc("plant_step_gradient kernel: [A|B] = integrator_gradient([q;qd], u, dt) per timestep "
                          "(full-smem scratch arena; pass-through to grid::integrator_gradient_device)",
                          [],
                          ["d_dAB is the [A|B] output (2*NUM_VEL*3*NUM_VEL per timestep, column-major)",
                           "d_x is the packed current state [q; qd] (NUM_POS+NUM_VEL per timestep)",
                           "d_u is the packed control torque (NUM_VEL per timestep)",
                           "stride_x / stride_u are the per-timestep strides",
                           "d_robotModel / gravity / dt as for plant_step",
                           "NUM_TIMESTEPS is the batch size"], None)
    self.gen_add_code_line("template <typename T, grid::IntegratorType IT = grid::IntegratorType::EULER>")
    self.gen_add_code_line("__global__")
    self.gen_add_code_line("__launch_bounds__(grid::MAX_PERF_LEVEL_THREADS)")
    self.gen_add_code_line("void plant_step_gradient_kernel(T *d_dAB, const T *d_x, const T *d_u, "
                           "const int stride_x, const int stride_u, "
                           "const grid::robotModel<T> *d_robotModel, const T gravity, const T dt, const int NUM_TIMESTEPS) {", True)
    # The shared-arena machinery (grid_align_up / grid_arena_ptr / the *_BYTES
    # helpers) + IntegratorType / robotModel live in `namespace grid`; this
    # kernel is in the sibling `grid_plant`, so pull them into scope. The only
    # grid_plant symbol the body references (plant_step_gradient) is found by
    # enclosing-namespace lookup.
    self.gen_add_code_line("using namespace grid;")
    # Whole scratch arena in shared memory (PERF/full-smem; mirrors the
    # integrator_gradient_kernel _emit_body(dqdd_in_smem=True, dab_in_smem=True,
    # inner_level=0) layout). s_x (= s_q;s_qd) + s_u are staged from global; the
    # FD-grad bands, the dAB output, and the multi-stage scratch all live here.
    extra_t_buffers = [
        ("s_x", nx),
        ("s_u", n),
        ("s_dAB", 2 * n * 3 * n),
        ("s_df_du", n * 2 * n),
        ("s_dc_du", n * 2 * n),
        ("s_vaf", vaf_cnt),
        ("s_Minv", n * n),
        ("s_qdd", n),
        ("s_q_orig", n + fb),
        ("s_qd_orig", n),
        ("s_stage_grad_qdd", max_stages * n),
        ("s_D_qdd_stage", d_qdd_count),
        ("s_dInt_q_6x6", 36),
        ("s_dInt_v_6x6", 36),
    ]
    self.gen_XImats_helpers_temp_shared_memory_code(
        inner_temp_full, extra_t_buffers=extra_t_buffers, include_linalg_scratch=True)
    self.gen_add_parallel_loop("k", "NUM_TIMESTEPS", block_level=True)
    # stage x (q;qd) + u into smem (plant_step_gradient mutates q/qd across stages).
    self.gen_add_parallel_loop("ind", str(nx))
    self.gen_add_code_line("s_x[ind] = d_x[k*stride_x + ind];")
    self.gen_add_end_control_flow()
    self.gen_add_parallel_loop("ind", str(n))
    self.gen_add_code_line("s_u[ind] = d_u[k*stride_u + ind];")
    self.gen_add_end_control_flow()
    self.gen_add_sync()
    # All scratch is smem (SCRATCH_IN_SMEM=true, no spill): pass nullptr for the
    # global workspace + spill regions. The orchestrator owns its s_temp pool.
    self.gen_add_code_line("plant_step_gradient<T, IT, true, false>("
                           "s_dAB, s_x, s_u, s_df_du, s_dc_du, s_vaf, s_Minv, s_qdd, "
                           "s_q_orig, s_qd_orig, s_stage_grad_qdd, s_D_qdd_stage, "
                           "s_dInt_q_6x6, s_dInt_v_6x6, "
                           + self.gen_insert_helpers_function_call()
                           + "s_temp, nullptr, nullptr, d_robotModel, gravity, dt);")
    self.gen_add_sync()
    self.gen_add_serial_ops()
    self.gen_add_code_line("for (int ind = 0; ind < " + str(2 * n * 3 * n) + "; ++ind) d_dAB[k*" + str(2 * n * 3 * n) + " + ind] = s_dAB[ind];")
    self.gen_add_end_control_flow()
    self.gen_add_sync()
    self.gen_add_end_control_flow()
    self.gen_add_end_function()


def gen_com_cost_kernel(self):
    """`com_cost_kernel` — one block/timestep value+grad_x+GN-hess_x.

    Calls grid_plant::com_cost[_gradient/_hessian], which call the
    auto-allocating grid::com_device. Global in/out; reuses
    grid::COM_DYNAMIC_SHARED_MEM_BYTES for the launch smem. Mirrors
    gen_ee_pos_cost_kernel (CoM position/Jacobian in place of EE)."""
    nq = self.robot.get_num_pos()
    nv = self.robot.get_num_vel()
    nx = nq + nv
    com_out = 3 + 3 * nv  # grid::com_device output: [p_com(3); J_com(3 x NV)]
    self.gen_add_func_doc("com_cost_kernel: value + grad_x + GN hess_x per timestep",
                          [], ["d_out scalar cost (1 per timestep)",
                               "d_grad grad over x (" + str(nx) + " per timestep)",
                               "d_hess dense col-major x-hessian (" + str(nx*nx) + " per timestep)",
                               "d_q joint positions (NUM_POS per timestep)",
                               "d_p_des desired CoM position (3 per timestep)",
                               "d_W per-axis weight (3 per timestep)",
                               "d_com global scratch (3 + 3*NUM_VEL per timestep)",
                               "NUM_TIMESTEPS is the batch size"], None)
    self.gen_add_code_line("template <typename T>")
    self.gen_add_code_line("__global__")
    self.gen_add_code_line("void com_cost_kernel(T *d_out, T *d_grad, T *d_hess, "
                           "const T *d_q, const T *d_p_des, const T *d_W, T *d_com, "
                           "const grid::robotModel<T> *d_robotModel, const int NUM_TIMESTEPS) {", True)
    self.gen_add_parallel_loop("k", "NUM_TIMESTEPS", block_level=True)
    self.gen_add_code_line("const T *s_q = &d_q[k*" + str(nq) + "]; const T *s_p_des = &d_p_des[k*3]; const T *s_W = &d_W[k*3];")
    self.gen_add_code_line("T *s_com = &d_com[k*" + str(com_out) + "];")
    self.gen_add_code_line("com_cost<T>(&d_out[k], s_q, s_p_des, s_W, s_com, d_robotModel);")
    self.gen_add_sync()
    self.gen_add_code_line("com_cost_gradient<T>(&d_grad[k*" + str(nx) + "], s_q, s_p_des, s_W, s_com, d_robotModel);")
    self.gen_add_sync()
    self.gen_add_code_line("com_cost_hessian<T>(&d_hess[k*" + str(nx*nx) + "], s_q, s_W, s_com, d_robotModel);")
    self.gen_add_sync()
    self.gen_add_end_control_flow()
    self.gen_add_end_function()


def gen_momentum_cost_kernel(self):
    """`momentum_cost_kernel` — one block/timestep value+grad_x+GN-hess_x.

    Calls grid_plant::momentum_cost[_gradient/_hessian], which call the
    auto-allocating grid::ccrba_device. Global in/out; reuses
    grid::CCRBA_DYNAMIC_SHARED_MEM_BYTES for the launch smem. Reads q + qd (the
    momentum h = A qd depends on qd)."""
    nq = self.robot.get_num_pos()
    nv = self.robot.get_num_vel()
    nx = nq + nv
    ccrba_out = 6 * nv + 6  # grid::ccrba_device output: [A(6 x NV); h(6)]
    self.gen_add_func_doc("momentum_cost_kernel: value + grad_x + GN hess_x per timestep",
                          [], ["d_out scalar cost (1 per timestep)",
                               "d_grad grad over x (" + str(nx) + " per timestep)",
                               "d_hess dense col-major x-hessian (" + str(nx*nx) + " per timestep)",
                               "d_q / d_qd joint positions/velocities (NUM_POS / NUM_VEL per timestep)",
                               "d_h_des desired centroidal momentum (6 per timestep)",
                               "d_W per-component weight (6 per timestep)",
                               "d_ccrba global scratch (6*NUM_VEL + 6 per timestep)",
                               "NUM_TIMESTEPS is the batch size"], None)
    self.gen_add_code_line("template <typename T>")
    self.gen_add_code_line("__global__")
    self.gen_add_code_line("void momentum_cost_kernel(T *d_out, T *d_grad, T *d_hess, "
                           "const T *d_q, const T *d_qd, const T *d_h_des, const T *d_W, T *d_ccrba, "
                           "const grid::robotModel<T> *d_robotModel, const int NUM_TIMESTEPS) {", True)
    self.gen_add_parallel_loop("k", "NUM_TIMESTEPS", block_level=True)
    self.gen_add_code_line("const T *s_q = &d_q[k*" + str(nq) + "]; const T *s_qd = &d_qd[k*" + str(nv) + "];")
    self.gen_add_code_line("const T *s_h_des = &d_h_des[k*6]; const T *s_W = &d_W[k*6];")
    self.gen_add_code_line("T *s_ccrba = &d_ccrba[k*" + str(ccrba_out) + "];")
    self.gen_add_code_line("momentum_cost<T>(&d_out[k], s_q, s_qd, s_h_des, s_W, s_ccrba, d_robotModel);")
    self.gen_add_sync()
    self.gen_add_code_line("momentum_cost_gradient<T>(&d_grad[k*" + str(nx) + "], s_q, s_qd, s_h_des, s_W, s_ccrba, d_robotModel);")
    self.gen_add_sync()
    self.gen_add_code_line("momentum_cost_hessian<T>(&d_hess[k*" + str(nx*nx) + "], s_q, s_qd, s_W, s_ccrba, d_robotModel);")
    self.gen_add_sync()
    self.gen_add_end_control_flow()
    self.gen_add_end_function()


def gen_quadratic_cost_kernel(self, which):
    """`<base>_kernel` — one block/timestep value+grad+GN-diag-hess.

    No RBD arena needed; a small static `__shared__` reduction scratch suffices.
    """
    nq = self.robot.get_num_pos()
    nv = self.robot.get_num_vel()
    if which == "state":
        size = nq + nv
        var, des, w = "x", "x_des", "Q"
        base = "quadratic_state_cost"
    else:
        size = nv
        var, des, w = "u", "u_des", "R"
        base = "quadratic_input_cost"
    N = str(size)
    self.gen_add_func_doc(base + "_kernel: value + gradient + GN-diag hessian per timestep",
                          [], ["d_out scalar cost (1 per timestep)",
                               "d_grad gradient (" + N + " per timestep)",
                               "d_hess dense col-major hessian (" + N + "*" + N + " per timestep)",
                               "d_" + var + " / d_" + des + " / d_" + w + " inputs (" + N + " per timestep)",
                               "NUM_TIMESTEPS is the batch size"], None)
    self.gen_add_code_line("template <typename T>")
    self.gen_add_code_line("__global__")
    self.gen_add_code_line("void " + base + "_kernel(T *d_out, T *d_grad, T *d_hess, "
                           "const T *d_" + var + ", const T *d_" + des + ", const T *d_" + w + ", const int NUM_TIMESTEPS) {", True)
    self.gen_add_code_line("__shared__ T s_scratch[" + N + "];")
    self.gen_add_parallel_loop("k", "NUM_TIMESTEPS", block_level=True)
    self.gen_add_code_line(base + "_value_grad_hess<T>(&d_out[k], &d_grad[k*" + N + "], &d_hess[k*" + str(size*size) + "], "
                           "&d_" + var + "[k*" + N + "], &d_" + des + "[k*" + N + "], &d_" + w + "[k*" + N + "], s_scratch);")
    self.gen_add_sync()
    self.gen_add_end_control_flow()
    self.gen_add_end_function()


def gen_ee_pos_cost_kernel(self):
    """`ee_pos_cost_kernel` — one block/timestep value+grad_x+GN-hess_x.

    Calls grid_plant::ee_pos_cost[_gradient/_hessian], which call the
    auto-allocating grid::end_effector_pose[_gradient]_device. Global in/out;
    reuses grid::END_EFFECTOR_POSE_DYNAMIC_SHARED_MEM_BYTES for the launch smem.
    """
    nq = self.robot.get_num_pos()
    nv = self.robot.get_num_vel()
    num_ees = self.robot.get_total_leaf_nodes()
    nx = nq + nv
    self.gen_add_func_doc("ee_pos_cost_kernel: value + grad_x + GN hess_x per timestep (EE=0)",
                          [], ["d_out scalar cost (1 per timestep)",
                               "d_grad grad over x (" + str(nx) + " per timestep)",
                               "d_hess dense col-major x-hessian (" + str(nx*nx) + " per timestep)",
                               "d_q joint positions (NUM_POS per timestep)",
                               "d_p_des desired EE position (3 per timestep)",
                               "d_W per-axis weight (3 per timestep)",
                               "d_end_effector_pose / d_end_effector_pose_gradient global scratch (6*NUM_EES / 6*NUM_VEL*NUM_EES per timestep)",
                               "NUM_TIMESTEPS is the batch size"], None)
    self.gen_add_code_line("template <typename T, int EE = 0>")
    self.gen_add_code_line("__global__")
    self.gen_add_code_line("void ee_pos_cost_kernel(T *d_out, T *d_grad, T *d_hess, "
                           "const T *d_q, const T *d_p_des, const T *d_W, T *d_end_effector_pose, T *d_end_effector_pose_gradient, "
                           "const grid::robotModel<T> *d_robotModel, const int NUM_TIMESTEPS) {", True)
    self.gen_add_parallel_loop("k", "NUM_TIMESTEPS", block_level=True)
    self.gen_add_code_line("const T *s_q = &d_q[k*" + str(nq) + "]; const T *s_p_des = &d_p_des[k*3]; const T *s_W = &d_W[k*3];")
    self.gen_add_code_line("T *s_end_effector_pose = &d_end_effector_pose[k*" + str(6*num_ees) + "]; T *s_end_effector_pose_gradient = &d_end_effector_pose_gradient[k*" + str(6*nv*num_ees) + "];")
    self.gen_add_code_line("ee_pos_cost<T, EE>(&d_out[k], s_q, s_p_des, s_W, s_end_effector_pose, d_robotModel);")
    self.gen_add_sync()
    self.gen_add_code_line("ee_pos_cost_gradient<T, EE>(&d_grad[k*" + str(nx) + "], s_q, s_p_des, s_W, s_end_effector_pose, s_end_effector_pose_gradient, d_robotModel);")
    self.gen_add_sync()
    self.gen_add_code_line("ee_pos_cost_hessian<T, EE>(&d_hess[k*" + str(nx*nx) + "], s_q, s_W, s_end_effector_pose_gradient, d_robotModel);")
    self.gen_add_sync()
    self.gen_add_end_control_flow()
    self.gen_add_end_function()


def gen_barrier_kernel(self, base, count, var_offset):
    """`<base>_kernel` — one block/timestep value+grad+hess-diag for a barrier.

    `count` bounded DOFs, reading the var slice at `var_offset`. Grad/hess are
    written into standalone packed buffers (offset 0) since the Python surface
    returns the per-DOF gradient/hessian-diagonal directly.
    """
    N = str(count)
    self.gen_add_func_doc(base + "_kernel: value + grad + hess-diagonal per timestep",
                          [], ["d_out scalar barrier cost (1 per timestep)",
                               "d_grad per-DOF gradient (" + N + " per timestep)",
                               "d_hess_diag per-DOF hessian diagonal (" + N + " per timestep)",
                               "d_var variable vector (" + N + " per timestep)",
                               "d_lower / d_upper per-DOF bounds (" + N + " per timestep)",
                               "mu barrier weight; NUM_TIMESTEPS is the batch size"], None)
    self.gen_add_code_line("template <typename T>")
    self.gen_add_code_line("__global__")
    self.gen_add_code_line("void " + base + "_kernel(T *d_out, T *d_grad, T *d_hess_diag, "
                           "const T *d_var, const T *d_lower, const T *d_upper, const T mu, const int NUM_TIMESTEPS) {", True)
    self.gen_add_code_line("__shared__ T s_scratch[" + N + "];")
    self.gen_add_parallel_loop("k", "NUM_TIMESTEPS", block_level=True)
    # zero the scalar out (the value fn ADDS into it), then run the three fns.
    self.gen_add_serial_ops()
    self.gen_add_code_line("d_out[k] = static_cast<T>(0);")
    self.gen_add_end_control_flow()
    self.gen_add_sync()
    self.gen_add_code_line("const T *s_var = &d_var[k*" + N + "]; const T *s_lo = &d_lower[k*" + N + "]; const T *s_hi = &d_upper[k*" + N + "];")
    # value reads s_var[VAR_OFFSET + i]; we pass an offset-0 view, so call with VAR_OFFSET=0 semantics.
    self.gen_add_code_line(base + "<T>(&d_out[k], s_var, s_lo, s_hi, mu, s_scratch);")
    self.gen_add_sync()
    # grad/hess templates default VAR_OFFSET to the packed slice offset; we pass an
    # offset-0 view of s_var and want offset-0 writes, so force both offsets to 0.
    self.gen_add_parallel_loop("i", N)
    self.gen_add_code_line("d_grad[k*" + N + " + i] = grid_plant_log_barrier_grad<T>(s_var[i], s_lo[i], s_hi[i], mu);")
    self.gen_add_code_line("d_hess_diag[k*" + N + " + i] = grid_plant_log_barrier_hess<T>(s_var[i], s_lo[i], s_hi[i], mu);")
    self.gen_add_end_control_flow()
    self.gen_add_sync()
    self.gen_add_end_control_flow()
    self.gen_add_end_function()


def gen_plant_kernels(self, algorithms):
    """Emit the binding-layer kernels for the plant value surface."""
    self.gen_quadratic_cost_kernel("state")
    self.gen_quadratic_cost_kernel("input")
    nq = self.robot.get_num_pos()
    nv = self.robot.get_num_vel()
    gen_barrier_kernel(self, "joint_position_barrier", nq, 0)
    gen_barrier_kernel(self, "joint_velocity_barrier", nv, nq)
    gen_barrier_kernel(self, "joint_torque_barrier",   nv, 0)
    if "integrator" in algorithms:
        gen_plant_step_kernel(self)
        # Signal to the binding layer (wrapper_template.cu) that plant_step exists.
        self.gen_add_code_line("#define GRID_PLANT_HAS_STEP 1")
    if ("integrator_gradient" in algorithms) or ("integrator_with_gradient" in algorithms):
        gen_plant_step_gradient_kernel(self)
        self.gen_add_code_line("#define GRID_PLANT_HAS_STEP_GRADIENT 1")
    if ("end_effector_pose" in algorithms) and ("end_effector_pose_gradient" in algorithms):
        gen_ee_pos_cost_kernel(self)
        self.gen_add_code_line("#define GRID_PLANT_HAS_EE_COST 1")
    # CoM / centroidal-momentum cost kernels emit whenever their device fns do
    # (gated identically to gen_com_cost/gen_momentum_cost in gen_grid_plant:
    # require grid::com_device + grid::ccrba_device, non-mimic).
    centroidal_ok = ("com" in algorithms and "ccrba" in algorithms
                     and not self.robot_has_mimic_joints())
    if centroidal_ok:
        gen_com_cost_kernel(self)
        self.gen_add_code_line("#define GRID_PLANT_HAS_COM_COST 1")
        gen_momentum_cost_kernel(self)
        self.gen_add_code_line("#define GRID_PLANT_HAS_MOMENTUM_COST 1")


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
    if ("end_effector_pose" in algorithms) and ("end_effector_pose_gradient" in algorithms):
        self.gen_ee_pos_cost()
    else:
        self.gen_add_code_line("// [grid_plant] ee_pos_cost skipped: requires both 'end_effector_pose' and 'end_effector_pose_gradient' (grid::end_effector_pose[_gradient]_device) — not generated.")

    # CoM-tracking / centroidal-momentum-tracking costs need the centroidal
    # kinematics-domain device fns (grid::com_device / grid::ccrba_device), which
    # are emitted ONLY when their `com` / `ccrba` keys are selected (and the robot
    # is non-mimic). Gating on `end_effector_pose` alone was wrong: a profile that
    # pulls in ee_pose for some OTHER reason (e.g. the frame_jacobian family, whose
    # normalization adds end_effector_pose) but does not request com/ccrba would
    # emit com_cost/momentum_cost referencing undefined grid::com_device/ccrba_device.
    centroidal_ok = ("com" in algorithms and "ccrba" in algorithms
                     and not self.robot_has_mimic_joints())
    if centroidal_ok:
        gen_com_cost(self)
        gen_momentum_cost(self)
    else:
        self.gen_add_code_line("// [grid_plant] com_cost/momentum_cost skipped: require grid::com_device/ccrba_device (need 'com'+'ccrba', non-mimic).")

    # Binding layer (G1): emit the per-timestep kernels that wrap the device
    # functions above, so the grid_rbd Python/C-ABI surface can launch them.
    self.gen_plant_kernels(algorithms)

    self.gen_add_end_control_flow()  # close namespace grid_plant
