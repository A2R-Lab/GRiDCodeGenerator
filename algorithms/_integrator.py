"""Time-integrator value codegen.

Mirrors `_forward_dynamics.py` (inner / device / kernel / host layers) but
emits a single time step `x_{k+1} = integrator(x_k, u_k, dt; f_dyn)` where
`f_dyn` is the existing forward dynamics. State is `x = [q (nq); qd (nv)]`
(fixed-base, so `nq == nv == n`), control `u` is size `n`, output `x_{k+1}`
is size `2n`. dt is a per-call scalar threaded through device/kernel/host.

The integrator type is selected at compile time by an `IntegratorType IT`
template parameter. Only `EULER` is wired up here; the `_dispatch` helper
below is structured so adding semi-implicit Euler / Midpoint / RK3 / RK4
later is purely additive (one extra `if constexpr (IT == ...)` branch).
"""


# integrator name <-> codegen-side string constant
_INTEGRATOR_TYPES = ("EULER", "SEMI_IMPLICIT_EULER", "MIDPOINT", "RK3", "RK4")

# Number of forward-dynamics evaluations each integrator type requires.
# Used at codegen time to size shared-memory buffers (per-stage qdd) and to
# guide which stage-computation branches are emitted.
_STAGE_COUNT = {
    "EULER": 1,
    "SEMI_IMPLICIT_EULER": 1,
    "MIDPOINT": 2,
    "RK3": 3,
    "RK4": 4,
}


def _max_stages_in_use():
    """Maximum stage count among all currently-emitted integrator types.

    For now, the kernel statically allocates per-stage scratch sized for the
    most-expensive integrator (RK4). This keeps shared-memory layout simple
    and the cost is small (a handful of extra n-sized buffers).
    """
    return max(_STAGE_COUNT.values())


def _integrator_type_token(integrator_type):
    """Either a known enum value (compile-time enum) or a raw template
    parameter passthrough (e.g. "IT" inside a `template <..., IntegratorType IT>`
    scope)."""
    if integrator_type in _INTEGRATOR_TYPES:
        return "IntegratorType::" + integrator_type
    return integrator_type


def gen_integrator_inner_temp_mem_size(self):
    # Integrator's only extra scratch is the FD itself; the assembly step is
    # in-place over a parallel loop with no additional storage.
    return self.gen_forward_dynamics_inner_temp_mem_size()


def gen_integrator_finish_function_call(self, integrator_type="IT", use_thread_group=False, updated_var_names=None):
    var_names = dict(
        s_x_kp1_name="s_x_kp1",
        s_q_name="s_q",
        s_qd_name="s_qd",
        s_qdd_name="s_qdd",
        dt_name="dt",
    )
    if updated_var_names is not None:
        for key, value in updated_var_names.items():
            var_names[key] = value
    code = ("integrator_finish<T, " + _integrator_type_token(integrator_type) + ">(" +
            var_names["s_x_kp1_name"] + ", " +
            var_names["s_q_name"] + ", " +
            var_names["s_qd_name"] + ", " +
            var_names["s_qdd_name"] + ", " +
            var_names["dt_name"] + ");")
    if use_thread_group:
        code = code.replace("(", "(tgrp, ", 1)
    self.gen_add_code_line(code)


def gen_integrator_finish(self, use_thread_group=False):
    """Emit a templated `integrator_finish<T, IntegratorType IT>` device function.

    For EULER:
        x_{k+1}[i]    = q[i]  + dt * qd[i]    for i in [0, n)     // q + dt*qd
        x_{k+1}[n+i]  = qd[i] + dt * qdd[i]   for i in [0, n)     // qd + dt*qdd
    Assumes the underlying forward dynamics has already populated s_qdd.
    """
    n = self.robot.get_num_vel()
    func_params = ["s_x_kp1 is a pointer to memory for the next state (size 2*NUM_VEL)",
                   "s_q is the vector of joint positions",
                   "s_qd is the vector of joint velocities",
                   "s_qdd is the vector of joint accelerations (output of forward_dynamics)",
                   "dt is the integration timestep"]
    func_def = "void integrator_finish(T *s_x_kp1, const T *s_q, const T *s_qd, const T *s_qdd, const T dt) {"
    func_notes = ["Assumes s_qdd is already computed for the current (s_q, s_qd, s_u)",
                  "Does not internally sync the thread group, so it should be called after all threads have finished computing their values"]
    if use_thread_group:
        func_def = func_def.replace("(", "(cgrps::thread_group tgrp, ", 1)
        func_params.insert(0, "tgrp is the handle to the thread_group running this function")
    self.gen_add_func_doc("Finish the integrator step: write x_{k+1} from (q, qd, qdd) per the integrator type",
                          func_notes, func_params, None)
    self.gen_add_code_line("template <typename T, IntegratorType IT>")
    self.gen_add_code_line("__device__")
    self.gen_add_code_line(func_def, True)

    # Single parallel loop over [0, 2n) that handles both halves.
    self.gen_add_parallel_loop("ind", str(2 * n), use_thread_group)
    self.gen_add_code_line("if constexpr (IT == IntegratorType::EULER) {", True)
    self.gen_add_code_line("if (ind < " + str(n) + ") {")
    self.gen_add_code_line("    s_x_kp1[ind] = s_q[ind] + dt * s_qd[ind];")
    self.gen_add_code_line("} else {")
    self.gen_add_code_line("    int j = ind - " + str(n) + ";")
    self.gen_add_code_line("    s_x_kp1[ind] = s_qd[j] + dt * s_qdd[j];")
    self.gen_add_code_line("}")
    self.gen_add_end_control_flow()  # end if constexpr EULER
    # Semi-implicit Euler: v_{k+1} = v + dt*qdd; q_{k+1} = q + dt*v_{k+1}.
    self.gen_add_code_line("else if constexpr (IT == IntegratorType::SEMI_IMPLICIT_EULER) {", True)
    self.gen_add_code_line("if (ind < " + str(n) + ") {")
    self.gen_add_code_line("    int j = ind;")
    self.gen_add_code_line("    T vkp1 = s_qd[j] + dt * s_qdd[j];")
    self.gen_add_code_line("    s_x_kp1[ind] = s_q[j] + dt * vkp1;")
    self.gen_add_code_line("} else {")
    self.gen_add_code_line("    int j = ind - " + str(n) + ";")
    self.gen_add_code_line("    s_x_kp1[ind] = s_qd[j] + dt * s_qdd[j];")
    self.gen_add_code_line("}")
    self.gen_add_end_control_flow()  # end if constexpr SI_EULER
    # Midpoint/RK3/RK4 use the multi-stage finish path (driven from
    # integrator_inner directly), which calls `integrator_finish_multistage`.
    # The single-stage finish here is only used by Euler / Semi-Implicit Euler.
    self.gen_add_code_line("else {", True)
    self.gen_add_code_line("static_assert(IT == IntegratorType::EULER || IT == IntegratorType::SEMI_IMPLICIT_EULER,")
    self.gen_add_code_line("              \"integrator_finish is only for single-stage integrators; multi-stage uses inner directly.\");")
    self.gen_add_end_control_flow()
    self.gen_add_end_control_flow()  # end parallel loop
    self.gen_add_end_function()


def gen_integrator_inner_function_call(self, integrator_type="IT", use_thread_group=False, updated_var_names=None):
    var_names = dict(
        s_x_kp1_name="s_x_kp1",
        s_q_name="s_q",
        s_qd_name="s_qd",
        s_u_name="s_u",
        s_qdd_name="s_qdd",
        s_stage_qdd_name="s_stage_qdd",
        s_stage_point_name="s_stage_point",
        d_robotModel_name="d_robotModel",
        s_temp_name="s_temp",
        dt_name="dt",
        gravity_name="gravity",
    )
    if updated_var_names is not None:
        for key, value in updated_var_names.items():
            var_names[key] = value
    code_start = ("integrator_inner<T, " + _integrator_type_token(integrator_type) + ">(" +
                  var_names["s_x_kp1_name"] + ", " +
                  var_names["s_q_name"] + ", " +
                  var_names["s_qd_name"] + ", " +
                  var_names["s_u_name"] + ", " +
                  var_names["s_qdd_name"] + ", " +
                  var_names["s_stage_qdd_name"] + ", " +
                  var_names["s_stage_point_name"] + ", ")
    code_end = (var_names["d_robotModel_name"] + ", " +
                var_names["s_temp_name"] + ", " +
                var_names["gravity_name"] + ", " +
                var_names["dt_name"] + ");")
    code_middle = self.gen_insert_helpers_function_call()
    if use_thread_group:
        code_start = code_start.replace("(", "(tgrp, ", 1)
    self.gen_add_code_line(code_start + code_middle + code_end)


def gen_integrator_inner(self, use_thread_group=False):
    """Templated inner: invokes forward_dynamics_inner(es) then either the
    single-stage integrator_finish or a multi-stage weighted assembly.

    Caller owns:
      - `s_qdd`: stage-1 qdd output (size n) — always used.
      - `s_stage_qdd`: stages 2..N qdd outputs (size (max_stages-1)*n) — only
        used for multi-stage integrators (Midpoint/RK3/RK4).
      - `s_stage_point`: intermediate state scratch (size (max_stages-1)*2n) —
        only used for multi-stage integrators.
    For Euler/SI-Euler, `s_stage_qdd` / `s_stage_point` are allocated but
    never touched.
    """
    n = self.robot.get_num_vel()
    max_stages = _max_stages_in_use()
    extra_qdd_count = (max_stages - 1) * n
    extra_point_count = (max_stages - 1) * 2 * n
    func_params = ["s_x_kp1 is a pointer to memory for the next state (size 2*NUM_VEL)",
                   "s_q is the vector of joint positions",
                   "s_qd is the vector of joint velocities",
                   "s_u is the vector of joint input torques",
                   "s_qdd is shared memory for the stage-1 joint accelerations (size NUM_VEL)",
                   "s_stage_qdd is shared memory for stages 2..N qdd outputs (size " + str(extra_qdd_count) + ")",
                   "s_stage_point is shared memory for stages 2..N intermediate (q,qd) states (size " + str(extra_point_count) + ")",
                   "s_temp is the pointer to the shared memory needed of size: " +
                       str(self.gen_integrator_inner_temp_mem_size()),
                   "gravity is the gravity constant",
                   "dt is the integration timestep"]
    func_def_start = ("void integrator_inner(T *s_x_kp1, const T *s_q, const T *s_qd, const T *s_u, "
                      "T *s_qdd, T *s_stage_qdd, T *s_stage_point, ")
    # d_robotModel is needed by multi-stage integrators to recompute s_XImats
    # at intermediate states. For single-stage (Euler / SI Euler) it's unused.
    func_def_end = "const robotModel<T> *d_robotModel, T *s_temp, const T gravity, const T dt) {"
    func_def_start, func_params = self.gen_insert_helpers_func_def_params(func_def_start, func_params, -3)
    func_notes = ["Assumes s_XImats is updated already for the current s_q",
                  "For Midpoint/RK3/RK4, re-runs forward_dynamics at intermediate states and weights stage qdd outputs."]
    if use_thread_group:
        func_def_start = func_def_start.replace("(", "(cgrps::thread_group tgrp, ", 1)
        func_params.insert(0, "tgrp is the handle to the thread_group running this function")
    self.gen_add_func_doc("Computes a single integrator step (x_{k+1} = integrator(x_k, u_k, dt))",
                          func_notes, func_params, None)
    self.gen_add_code_line("template <typename T, IntegratorType IT>")
    self.gen_add_code_line("__device__")
    self.gen_add_code_line(func_def_start + func_def_end, True)

    # Stage 1: always run forward dynamics on (q, qd, u).
    self.gen_forward_dynamics_inner_function_call(use_thread_group)
    self.gen_add_sync(use_thread_group)

    # Single-stage branch — Euler / Semi-Implicit Euler.
    self.gen_add_code_line("if constexpr (IT == IntegratorType::EULER || IT == IntegratorType::SEMI_IMPLICIT_EULER) {", True)
    self.gen_integrator_finish_function_call(integrator_type="IT", use_thread_group=use_thread_group)
    self.gen_add_end_control_flow()

    # Multi-stage branch — emits each subsequent stage in turn, with an
    # if-constexpr to gate which stages actually run for which IT.
    # All multi-stage IT values share the same stage-driver structure with
    # different Butcher coefficients selected at compile time.
    self.gen_add_code_line("else {", True)
    # Aliases for stage-scratch slices.
    self.gen_add_code_line("T *s_qdd_2 = &s_stage_qdd[0];")
    self.gen_add_code_line("T *s_p1_q  = &s_stage_point[0];")
    self.gen_add_code_line("T *s_p1_qd = &s_stage_point[" + str(n) + "];")
    if max_stages >= 3:
        self.gen_add_code_line("T *s_qdd_3 = &s_stage_qdd[" + str(n) + "];")
        self.gen_add_code_line("T *s_p2_q  = &s_stage_point[" + str(2 * n) + "];")
        self.gen_add_code_line("T *s_p2_qd = &s_stage_point[" + str(3 * n) + "];")
    if max_stages >= 4:
        self.gen_add_code_line("T *s_qdd_4 = &s_stage_qdd[" + str(2 * n) + "];")
        self.gen_add_code_line("T *s_p3_q  = &s_stage_point[" + str(4 * n) + "];")
        self.gen_add_code_line("T *s_p3_qd = &s_stage_point[" + str(5 * n) + "];")

    # ----- Stage 2 (Midpoint / RK3 / RK4): p1 = x + c1*dt*[qd; qdd_1] -----
    # Midpoint: c1 = 0.5. RK3: c1 = 0.5. RK4: c1 = 0.5. (All three use 0.5
    # for stage 2's offset.)
    self.gen_add_code_line("constexpr T c1 = static_cast<T>(0.5);")
    self.gen_add_parallel_loop("ind", str(n), use_thread_group)
    self.gen_add_code_line("s_p1_q[ind]  = s_q[ind]  + c1 * dt * s_qd[ind];")
    self.gen_add_code_line("s_p1_qd[ind] = s_qd[ind] + c1 * dt * s_qdd[ind];")
    self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)
    # IMPORTANT: re-derive s_XImats for the stage-2 configuration before
    # invoking FD — the helper was last populated for s_q (stage 1).
    self.gen_load_update_XImats_helpers_function_call(use_thread_group, updated_var_names=dict(s_q_name="s_p1_q"))
    self.gen_add_sync(use_thread_group)
    # FD at p1.
    self.gen_forward_dynamics_inner_function_call(use_thread_group, updated_var_names=dict(
        s_q_name="s_p1_q", s_qd_name="s_p1_qd", s_qdd_name="s_qdd_2",
    ))
    self.gen_add_sync(use_thread_group)

    # ----- Stage 3 (RK3 / RK4) -----
    if max_stages >= 3:
        self.gen_add_code_line("if constexpr (IT == IntegratorType::RK3 || IT == IntegratorType::RK4) {", True)
        # TrajoptPlant convention: xdot_i = [qd; qdd_i] (note: uses original qd,
        # NOT the stage-i velocity). So p_2 = xk + c2*dt*xdot_2 means
        # p_2.q = q + c2*dt*qd, p_2.qd = qd + c2*dt*qdd_2.
        # RK3: c2 = 0.75 (point2 = xk + 0.75*dt*xdot_2)
        # RK4: c2 = 0.5  (point2 = xk + 0.5*dt*xdot_2)
        self.gen_add_code_line("constexpr T c2 = (IT == IntegratorType::RK3) ? static_cast<T>(0.75) : static_cast<T>(0.5);")
        self.gen_add_parallel_loop("ind", str(n), use_thread_group)
        self.gen_add_code_line("s_p2_q[ind]  = s_q[ind]  + c2 * dt * s_qd[ind];")
        self.gen_add_code_line("s_p2_qd[ind] = s_qd[ind] + c2 * dt * s_qdd_2[ind];")
        self.gen_add_end_control_flow()
        self.gen_add_sync(use_thread_group)
        self.gen_load_update_XImats_helpers_function_call(use_thread_group, updated_var_names=dict(s_q_name="s_p2_q"))
        self.gen_add_sync(use_thread_group)
        self.gen_forward_dynamics_inner_function_call(use_thread_group, updated_var_names=dict(
            s_q_name="s_p2_q", s_qd_name="s_p2_qd", s_qdd_name="s_qdd_3",
        ))
        self.gen_add_sync(use_thread_group)
        self.gen_add_end_control_flow()

    # ----- Stage 4 (RK4 only) -----
    if max_stages >= 4:
        self.gen_add_code_line("if constexpr (IT == IntegratorType::RK4) {", True)
        self.gen_add_code_line("constexpr T c3 = static_cast<T>(1.0);")
        self.gen_add_parallel_loop("ind", str(n), use_thread_group)
        self.gen_add_code_line("s_p3_q[ind]  = s_q[ind]  + c3 * dt * s_qd[ind];")
        self.gen_add_code_line("s_p3_qd[ind] = s_qd[ind] + c3 * dt * s_qdd_3[ind];")
        self.gen_add_end_control_flow()
        self.gen_add_sync(use_thread_group)
        self.gen_load_update_XImats_helpers_function_call(use_thread_group, updated_var_names=dict(s_q_name="s_p3_q"))
        self.gen_add_sync(use_thread_group)
        self.gen_forward_dynamics_inner_function_call(use_thread_group, updated_var_names=dict(
            s_q_name="s_p3_q", s_qd_name="s_p3_qd", s_qdd_name="s_qdd_4",
        ))
        self.gen_add_sync(use_thread_group)
        self.gen_add_end_control_flow()

    # ----- Final assembly: x_{k+1} = xk + dt * sum(b_i * xdot_i) -----
    # In TrajoptPlant's convention, xdot_i = [qd; qdd_i] with the SAME qd for
    # every stage. So:
    #   x_{k+1}.q  = q + dt * (sum b_i) * qd  = q + dt * qd     (b's always sum to 1)
    #   x_{k+1}.qd = qd + dt * sum(b_i * qdd_i)
    # That keeps the q update Euler-style across all multi-stage variants.
    self.gen_add_code_line("// final assembly: q_{k+1} = q + dt*qd; qd_{k+1} = qd + dt*sum(b_i * qdd_i)")
    self.gen_add_parallel_loop("ind", str(2 * n), use_thread_group)
    self.gen_add_code_line("if (ind < " + str(n) + ") {")
    self.gen_add_code_line("    s_x_kp1[ind] = s_q[ind] + dt * s_qd[ind];")
    self.gen_add_code_line("} else {")
    self.gen_add_code_line("    int j = ind - " + str(n) + ";")
    self.gen_add_code_line("    T accel = static_cast<T>(0);")
    self.gen_add_code_line("    if constexpr (IT == IntegratorType::MIDPOINT) {")
    self.gen_add_code_line("        accel = s_qdd_2[j];")
    self.gen_add_code_line("    } else if constexpr (IT == IntegratorType::RK3) {")
    self.gen_add_code_line("        constexpr T b1 = static_cast<T>(2.0/9.0);")
    self.gen_add_code_line("        constexpr T b2 = static_cast<T>(3.0/9.0);")
    self.gen_add_code_line("        constexpr T b3 = static_cast<T>(4.0/9.0);")
    self.gen_add_code_line("        accel = b1 * s_qdd[j] + b2 * s_qdd_2[j] + b3 * s_qdd_3[j];")
    self.gen_add_code_line("    } else if constexpr (IT == IntegratorType::RK4) {")
    self.gen_add_code_line("        constexpr T b1 = static_cast<T>(1.0/6.0);")
    self.gen_add_code_line("        constexpr T b2 = static_cast<T>(2.0/6.0);")
    self.gen_add_code_line("        constexpr T b3 = static_cast<T>(2.0/6.0);")
    self.gen_add_code_line("        constexpr T b4 = static_cast<T>(1.0/6.0);")
    self.gen_add_code_line("        accel = b1 * s_qdd[j] + b2 * s_qdd_2[j] + b3 * s_qdd_3[j] + b4 * s_qdd_4[j];")
    self.gen_add_code_line("    }")
    self.gen_add_code_line("    s_x_kp1[ind] = s_qd[j] + dt * accel;")
    self.gen_add_code_line("}")
    self.gen_add_end_control_flow()
    self.gen_add_end_control_flow()  # end else (multi-stage)
    self.gen_add_end_function()


def gen_integrator_device(self, use_thread_group=False):
    n = self.robot.get_num_vel()
    func_params = ["s_x_kp1 is a pointer to memory for the next state (size 2*NUM_VEL)",
                   "s_q is the vector of joint positions",
                   "s_qd is the vector of joint velocities",
                   "s_u is the vector of joint input torques",
                   "d_robotModel is the pointer to the initialized model specific helpers on the GPU (XImats, topology_helpers, etc.)",
                   "gravity is the gravity constant",
                   "dt is the integration timestep"]
    func_def_start = "void integrator_device(T *s_x_kp1, const T *s_q, const T *s_qd, const T *s_u, "
    func_def_end = "const robotModel<T> *d_robotModel, const T gravity, const T dt) {"
    if use_thread_group:
        func_def_start = func_def_start.replace("(", "(cgrps::thread_group tgrp, ", 1)
        func_params.insert(0, "tgrp is the handle to the thread_group running this function")
    self.gen_add_func_doc("Computes a single integrator step using the precomputed robotModel",
                          [], func_params, None)
    self.gen_add_code_line("template <typename T, IntegratorType IT = IntegratorType::EULER>")
    self.gen_add_code_line("__device__")
    self.gen_add_code_line(func_def_start + func_def_end, True)
    shared_mem_size = self.gen_integrator_inner_temp_mem_size()
    max_stages = _max_stages_in_use()
    extra_t_buffers = [
        ("s_qdd", n),
        ("s_stage_qdd", (max_stages - 1) * n),
        ("s_stage_point", (max_stages - 1) * 2 * n),
    ]
    self.gen_XImats_helpers_temp_shared_memory_code(shared_mem_size, extra_t_buffers=extra_t_buffers, include_linalg_scratch=True)
    self.gen_load_update_XImats_helpers_function_call(use_thread_group)
    self.gen_integrator_inner_function_call(integrator_type="IT", use_thread_group=use_thread_group)
    self.gen_add_end_function()


def gen_integrator_kernel(self, use_thread_group=False, single_call_timing=False):
    n = self.robot.get_num_vel()
    func_params = ["d_x_kp1 is a pointer to memory for the next state (size 2*NUM_VEL per timestep)",
                   "d_q_qd_u is the packed joint positions, velocities, and input torques",
                   "stride_q_qd_u is the stride between each (q, qd, u) tuple in d_q_qd_u",
                   "d_robotModel is the pointer to the initialized model specific helpers on the GPU",
                   "gravity is the gravity constant",
                   "dt is the integration timestep",
                   "num_timesteps is the length of the trajectory (or overloaded as test_iters for timing)"]
    func_def_start = "void integrator_kernel(T *d_x_kp1, const T *d_q_qd_u, const int stride_q_qd_u, "
    func_def_end = "const robotModel<T> *d_robotModel, const T gravity, const T dt, const int NUM_TIMESTEPS) {"
    func_def = func_def_start + func_def_end
    if single_call_timing:
        func_def = func_def.replace("kernel(", "kernel_single_timing(")
    self.gen_add_func_doc("Computes a single integrator step per timestep (Euler by default)",
                          [], func_params, None)
    self.gen_add_code_line("template <typename T, IntegratorType IT = IntegratorType::EULER, int RESOURCE_TIER = TIER_PERF>")
    self.gen_add_code_line("__global__")
    self.gen_add_code_line("__launch_bounds__(tier_max_threads<RESOURCE_TIER>())")
    self.gen_add_code_line(func_def, True)
    shared_mem_size = self.gen_integrator_inner_temp_mem_size()
    fb = self.robot.floating_base  # 0 for fixed-base
    max_stages = _max_stages_in_use()
    extra_t_buffers = [
        ("s_q_qd_u", 3 * n + fb),
        ("s_qdd", n),
        ("s_stage_qdd", (max_stages - 1) * n),
        ("s_stage_point", (max_stages - 1) * 2 * n),
        ("s_x_kp1", 2 * n),
    ]
    self.gen_XImats_helpers_temp_shared_memory_code(shared_mem_size, extra_t_buffers=extra_t_buffers, include_linalg_scratch=True)
    self.gen_add_code_line(
        "T *s_q = s_q_qd_u; T *s_qd = &s_q_qd_u[" + str(n + fb) + "]; T *s_u = &s_q_qd_u[" + str(2 * n + fb) + "];"
    )
    if use_thread_group:
        self.gen_add_code_line("cgrps::thread_group tgrp = TBD;")
    if not single_call_timing:
        self.gen_add_parallel_loop("k", "NUM_TIMESTEPS", use_thread_group, block_level=True)
        self.gen_kernel_load_inputs("q_qd_u", "stride_q_qd_u", str(3 * n + fb), use_thread_group)
        self.gen_add_code_line("// compute")
        self.gen_load_update_XImats_helpers_function_call(use_thread_group)
        self.gen_integrator_inner_function_call(integrator_type="IT", use_thread_group=use_thread_group)
        self.gen_add_sync(use_thread_group)
        self.gen_kernel_save_result("x_kp1", str(2 * n), str(2 * n), use_thread_group)
        self.gen_add_end_control_flow()
    else:
        input_count = 3 * n + fb
        self.gen_kernel_load_inputs_single_timing("q_qd_u", str(input_count))
        self.gen_add_code_line("// compute with NUM_TIMESTEPS as NUM_REPS for timing")
        self.gen_add_code_line("for (int rep = 0; rep < NUM_TIMESTEPS; rep++){", True)
        self.gen_anti_licm_input_reload("q_qd_u", str(input_count), use_thread_group, feedback_from="x_kp1")
        self.gen_load_update_XImats_helpers_function_call(use_thread_group)
        self.gen_integrator_inner_function_call(integrator_type="IT", use_thread_group=use_thread_group)
        self.gen_anti_licm_output_write("x_kp1")
        self.gen_add_end_control_flow()
        self.gen_kernel_save_result_single_timing("x_kp1", str(2 * n), use_thread_group)
    self.gen_add_end_function()


def gen_integrator_host(self, mode=0):
    single_call_timing = mode == 1
    compute_only = mode == 2
    func_params = ["hd_data is the packaged input and output pointers",
                   "d_robotModel is the pointer to the initialized model specific helpers on the GPU",
                   "gravity is the gravity constant",
                   "dt is the integration timestep",
                   "num_timesteps is the length of the trajectory (or overloaded as test_iters for timing)",
                   "streams are pointers to CUDA streams for async memory transfers (if needed)"]
    func_def_start = "void integrator(gridData<T, KIND> *hd_data, const robotModel<T> *d_robotModel, const T gravity, const T dt, const int num_timesteps,"
    func_def_end = "                  const dim3 block_dimms, const dim3 thread_dimms, cudaStream_t *streams) {"
    if single_call_timing:
        func_def_start = func_def_start.replace("(", "_single_timing(", 1)
        func_def_end = "              " + func_def_end
    if compute_only:
        func_def_start = func_def_start.replace("(", "_compute_only(", 1)
        func_def_end = "             " + func_def_end.replace(", cudaStream_t *streams", "")
    self.gen_add_func_doc("Run a single integrator step (default Euler) per timestep",
                          [], func_params, None)
    self.gen_add_code_line("template <typename T, IntegratorType IT = IntegratorType::EULER, gridDataKind KIND = GRID_DATA_ALL>")
    self.gen_add_code_line("__host__")
    self.gen_add_code_line(func_def_start)
    self.gen_add_code_line(func_def_end, True)
    self.gen_add_code_line("static_assert(KIND == GRID_DATA_ALL || KIND == GRID_DATA_DYNAMICS, \"integrator requires all-data or dynamics gridData\");")
    func_call_start = "integrator_kernel<T, IT><<<block_dimms,thread_dimms,INTEGRATOR_DYNAMIC_SHARED_MEM_BYTES<T>()>>>(hd_data->d_x_kp1,hd_data->d_q_qd_u,stride_q_qd_u,"
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
    self.gen_add_code_line("gpuErrchk(grid_check_dynamic_shared_memory_bytes(\"integrator\", INTEGRATOR_DYNAMIC_SHARED_MEM_BYTES<T>()));")
    self.gen_add_code_lines(func_call_code)
    if not compute_only:
        self.gen_add_code_lines([
            "// finally transfer the result back",
            "gpuErrchk(cudaMemcpy(hd_data->h_x_kp1,hd_data->d_x_kp1,2*NUM_JOINTS*" +
                ("num_timesteps*" if not single_call_timing else "") + "sizeof(T),cudaMemcpyDeviceToHost));",
            "gpuErrchkKernel();",
        ])
    if single_call_timing:
        from ..algo_registry import single_call_printf_line
        self.gen_add_code_line(single_call_printf_line("integrator"))
    self.gen_add_end_function()


def gen_integrator(self, use_thread_group=False):
    # Emit finish + inner (templated on IT), then EULER-typed device/kernel/host.
    self.gen_integrator_finish(use_thread_group)
    self.gen_integrator_inner(use_thread_group)
    self.gen_integrator_device(use_thread_group)
    self.gen_integrator_kernel(use_thread_group, single_call_timing=True)
    self.gen_integrator_kernel(use_thread_group, single_call_timing=False)
    self.gen_integrator_host(0)
    self.gen_integrator_host(1)
    self.gen_integrator_host(2)
