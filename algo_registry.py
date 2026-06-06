"""Single source of truth for GRiD algorithm metadata.

ONE verbose canonical name per algorithm, used IDENTICALLY as: registry `key` ==
emitted `grid::` symbol == bench key == algorithm_list token == printf label. There
are NO short aliases or legacy labels here — this is a deliberate clean break.

Each entry binds:
  - `key`: the JSON output key (lowercase) — also the parser's lowercased label match
  - `display`: the verbose label rendered in benchmark reports
  - `section`: the section heading in the report (Core Dynamics / Gradients / ...)

The codegen's printf label is `key.upper()` unless `printf_label` is overridden.
The parser then lowercases the matched line and looks it up against `key`.

Consumers:
  - `test/benchmarks/timing_parser.py` builds its label→key maps from this list.
  - `test/benchmarks/generate_report.py` builds ALGO_DISPLAY and ALGO_SECTIONS.
  - (Future) per-algo `gen_*_host` codegen functions read `printf_label` here.

Adding a new algo touches exactly this file (plus the corresponding codegen
emission of the kernel + the bench's `#if GRID_HAS_X` measure wrappers).
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class AlgoEntry:
    key: str
    display: str
    section: str
    printf_label_override: str | None = None

    @property
    def printf_label(self) -> str:
        return self.printf_label_override or self.key.upper()


ALGO_REGISTRY: tuple[AlgoEntry, ...] = (
    # Core Dynamics
    AlgoEntry("inverse_dynamics",     "Inverse Dynamics (RNEA / Recursive Newton-Euler Algorithm)",
              "Core Dynamics"),
    AlgoEntry("minv",                 "Minv (M⁻¹, computed directly)",
              "Core Dynamics"),
    AlgoEntry("forward_dynamics",     "Forward Dynamics (Minv+RNEA)",
              "Core Dynamics"),
    AlgoEntry("aba",                  "ABA (Articulated Body Algorithm)",
              "Core Dynamics"),
    AlgoEntry("crba",                 "CRBA",
              "Core Dynamics"),

    # Gradients
    AlgoEntry("inverse_dynamics_gradient", "Inverse Dynamics Gradient (∂ID/∂q,v)", "Gradients"),
    AlgoEntry("forward_dynamics_gradient", "Forward Dynamics Gradient (∂FD/∂q,v)", "Gradients"),
    AlgoEntry("f_ext_gradient",       "F_EXT_GRAD (∂tau/∂fext=-Jᵀ, ∂q̈/∂fext=M⁻¹Jᵀ)",
              "Gradients"),
    AlgoEntry("f_ext_gradient_dq",    "F_EXT_GRADIENT_DQ (∂(inverse_dynamics_gradient)/∂fext=-∂Jᵀ/∂q, fixed base)",
              "Gradients"),
    AlgoEntry("inverse_dynamics_regressor", "Inverse Dynamics Regressor (Joint-torque Y; tau=Y·π, ∂tau/∂π)",
              "Gradients"),
    AlgoEntry("forward_dynamics_parameter_gradient", "Forward Dynamics Parameter Gradient (∂q̈/∂π = -M⁻¹·Y)",
              "Gradients"),
    AlgoEntry("kinetic_energy_regressor", "Kinetic Energy Regressor (KE = y_KE·π, length 10·NB)",
              "Gradients"),
    AlgoEntry("potential_energy_regressor", "Potential Energy Regressor (PE = y_PE·π, length 10·NB)",
              "Gradients"),

    # Integrators
    AlgoEntry("integrator",           "Integrator (x_{k+1})",              "Integrators"),
    AlgoEntry("integrator_gradient",  "Integrator_Gradient (∂x_{k+1}/∂x,u)",
              "Integrators"),
    AlgoEntry("integrator_with_gradient",
              "Integrator_With_Gradient (x_{k+1} + ∂x_{k+1}/∂x,u)",
              "Integrators"),
    AlgoEntry("integrator_hessian", "Integrator_Hessian (∂²x_{k+1}/∂z², z=[q,qd,u]; plant_step_hessian s_d2AB)",
              "Integrators"),

    # Kinematics
    AlgoEntry("end_effector_pose",              "END_EFFECTOR_POSE",                 "Kinematics"),
    AlgoEntry("end_effector_pose_gradient",     "END_EFFECTOR_POSE_GRADIENT (Jacobian)", "Kinematics"),
    AlgoEntry("end_effector_pose_hessian",      "END_EFFECTOR_POSE_HESSIAN (2nd-order EE Jacobian)", "Kinematics"),
    AlgoEntry("frame_jacobian",       "FRAME_JACOBIAN (general-frame J: LOCAL/WORLD/LWA)", "Kinematics"),
    AlgoEntry("frame_jacobian_dot",   "FRAME_JACOBIAN_DOT (time derivative Jdot of the general-frame J)", "Kinematics"),
    AlgoEntry("osc_inertia",          "OSC_INERTIA (operational-space inertia Lambda = (J Minv J^T)^-1)", "Kinematics"),
    AlgoEntry("end_effector_pose_runtime",          "END_EFFECTOR_POSE_RUNTIME (runtime target/offset pose [xyz;rpy])", "Kinematics"),
    AlgoEntry("end_effector_pose_gradient_runtime", "END_EFFECTOR_POSE_GRADIENT_RUNTIME (runtime target/offset pose Jacobian)", "Kinematics"),

    # Second-Order
    AlgoEntry("idsva_so",             "IDSVA_SO (dispatched: body for fixed, world for floating)",
              "Second-Order"),
    AlgoEntry("idsva_so_body_frame",  "IDSVA_SO_BODY_FRAME (2nd-order ID, body-frame)",
              "Second-Order"),
    AlgoEntry("idsva_so_world_frame", "IDSVA_SO_WORLD_FRAME (2nd-order ID, world-frame)",
              "Second-Order"),
    AlgoEntry("fdsva_so",             "FDSVA_SO (2nd-order FD)",
              "Second-Order"),

    # Centroidal / Energy / CoM (G2 quick-wins, R1-R3)
    AlgoEntry("generalized_gravity", "Generalized Gravity g(q)=RNEA(q,0,0)",
              "Centroidal"),
    AlgoEntry("nonlinear_effects",   "Nonlinear Effects c(q,qd)=RNEA(q,qd,0)",
              "Centroidal"),
    AlgoEntry("energy",              "Energy (KE/PE/mechanical)",
              "Centroidal"),
    AlgoEntry("com",                 "CoM + CoM Jacobian",
              "Centroidal"),
    AlgoEntry("ccrba",               "CCRBA (A, h)",
              "Centroidal"),
    AlgoEntry("coriolis_matrix",     "Coriolis Matrix C(q,q̇)",
              "Centroidal"),
    AlgoEntry("dccrba",              "dCCRBA (∂A/∂q tensor, 6×NV×NV)",
              "Centroidal"),
    AlgoEntry("cmm_time_variation",  "CMM Time Variation (Ȧ, 6×NV)",
              "Centroidal"),

    # Plant (T6): cost / constraint / plant-step primitives emitted in the
    # sibling `grid_plant` namespace. No standalone benchmarked kernel — this
    # entry exists so the family has a registry key (display / sectioning) and
    # so future per-primitive bench wrappers can reference it.
    AlgoEntry("plant",                "Plant (cost/constraint/step primitives)",
              "Plant"),
)


def build_single_label_map() -> dict[str, str]:
    """e.g. {'single call inverse_dynamics': 'inverse_dynamics', 'single call idsva_so_body_frame': 'idsva_so_body_frame', ...}."""
    return {f"single call {e.printf_label.lower()}": e.key for e in ALGO_REGISTRY}


def build_batch_with_mem_label_map() -> dict[str, str]:
    """e.g. {'inverse_dynamics with memory': 'inverse_dynamics', 'fdsva_so with memory': 'fdsva_so', ...}."""
    return {f"{e.printf_label.lower()} with memory": e.key for e in ALGO_REGISTRY}


def build_batch_compute_only_label_map() -> dict[str, str]:
    """e.g. {'inverse_dynamics compute only': 'inverse_dynamics', 'fdsva_so compute only': 'fdsva_so', ...}."""
    return {f"{e.printf_label.lower()} compute only": e.key for e in ALGO_REGISTRY}


def build_display_map() -> dict[str, str]:
    """key → display string, for generate_report.ALGO_DISPLAY."""
    return {e.key: e.display for e in ALGO_REGISTRY}


def build_sections_map() -> dict[str, list[str]]:
    """section → [keys], preserving registry order, for generate_report.ALGO_SECTIONS."""
    out: dict[str, list[str]] = {}
    for entry in ALGO_REGISTRY:
        out.setdefault(entry.section, []).append(entry.key)
    return out


_BY_KEY = {e.key: e for e in ALGO_REGISTRY}


def printf_label_for(key: str) -> str:
    """Return the codegen printf label for an algorithm key.

    Used by per-algo `gen_*_host` codegen functions so the emitted
    `printf("Single Call X ...")` string is derived from this registry
    instead of hardcoded in each algorithm file.
    """
    try:
        return _BY_KEY[key].printf_label
    except KeyError:
        raise KeyError(f"{key!r} not in ALGO_REGISTRY — add it to algo_registry.py") from None


def single_call_printf_line(key: str) -> str:
    """Return the full C++ printf statement codegen should emit for single-call timing."""
    return (
        f'printf("Single Call {printf_label_for(key)} %fus\\n",'
        'time_delta_us_timespec(start,end)/static_cast<double>(num_timesteps));'
    )
