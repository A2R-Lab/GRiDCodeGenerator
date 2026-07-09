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
    # Collision (W3): baked sphere data + config_free emitted in the sibling
    # `grid_collision` namespace (composed over grid::multi_target_position + the
    # static SDF geometry header). No standalone benchmarked kernel — registry key
    # only, like `plant`.
    AlgoEntry("collision",            "Collision (config_free / spherized-URDF)",
              "Collision"),
)


# ─────────────────────────────────────────────────────────────────────────────
# Per-algo DESCRIPTOR table (item M, Step 0). One row per algorithm capturing the
# launch-config + kernel-attribute METADATA that is otherwise scattered across
# GRiDCodeGenerator.py. This step generates NOTHING — it is the parity safety net
# (test/test_algo_descriptor_parity.py asserts the table reproduces the live
# LAUNCH_CONFIG_ALGO_TO_SYMBOL dict and the KERNEL_ATTR_MANIFEST metadata exactly).
# Later steps (per docs/open-tasks/design_descriptor_table_spec.md) extend the
# schema with the arena/spill closures and DRIVE those sites from these rows.
#
# Only the IRREGULAR fields are stored per row; everything regular is derived:
#   - bytes_macro defaults to "<KEY.upper()>_DYNAMIC_SHARED_MEM_BYTES" — overridden
#     for the 4 algos that share a sibling's macro (integrator_gradient /
#     integrator_with_gradient -> INTEGRATOR_DU; generalized_gravity /
#     nonlinear_effects -> INVERSE_DYNAMICS_BIAS).
#   - gate_attr defaults to None (registration falls back to membership in
#     generated_algorithms) — overridden for the 7 opt-in / conditionally-emitted
#     kernels that carry an explicit `generate_*` / `_*_emitted` flag.
#   - has_kernel_attr defaults True — False for the 3 registry keys with no own
#     cudaFuncSetAttribute entry: `idsva_so` (a dispatch alias for the body/world
#     kernels), `integrator_hessian` (no standalone benchmarked kernel yet), and
#     `plant` (cost/constraint primitives, no __global__).
#   - autotune_keys defaults to () — set for the 17 algos that carry a baked
#     launch_cfg<> (the bench-abbreviated JSON key(s) that map to this grid symbol).


@dataclass(frozen=True)
class AlgoDescriptor:
    key: str                              # == AlgoEntry.key (the join)
    autotune_keys: tuple[str, ...] = ()   # bench JSON keys -> LAUNCH_CONFIG_ALGO_TO_SYMBOL
    has_kernel_attr: bool = True          # has its own KERNEL_ATTR_MANIFEST entry
    gate_attr: str | None = None          # explicit generate_*/_*_emitted gate, else None
    bytes_macro_stem: str | None = None   # override the default *_DYNAMIC_SHARED_MEM_BYTES stem

    @property
    def carries_launch_cfg(self) -> bool:
        return bool(self.autotune_keys)

    @property
    def enum_name(self) -> str:
        return "GRID_ALGO_" + self.key.upper()

    @property
    def bytes_macro(self) -> str:
        """Full `NAME<T>()` string as it appears in KERNEL_ATTR_MANIFEST."""
        stem = self.bytes_macro_stem or (self.key.upper() + "_DYNAMIC_SHARED_MEM_BYTES")
        return stem + "<T>()"


ALGO_DESCRIPTORS: tuple[AlgoDescriptor, ...] = (
    # Core Dynamics
    AlgoDescriptor("inverse_dynamics", autotune_keys=("id",)),
    AlgoDescriptor("minv", autotune_keys=("minv",)),
    AlgoDescriptor("forward_dynamics", autotune_keys=("fd",)),
    AlgoDescriptor("aba", autotune_keys=("aba",)),
    AlgoDescriptor("crba", autotune_keys=("crba",)),

    # Gradients
    AlgoDescriptor("inverse_dynamics_gradient", autotune_keys=("id_du",),
                   gate_attr="generate_inverse_dynamics_gradient"),
    AlgoDescriptor("forward_dynamics_gradient", autotune_keys=("fd_du",),
                   gate_attr="generate_forward_dynamics_gradient"),
    AlgoDescriptor("f_ext_gradient"),
    AlgoDescriptor("f_ext_gradient_dq", gate_attr="_f_ext_gradient_dq_emitted"),
    AlgoDescriptor("inverse_dynamics_regressor"),
    AlgoDescriptor("forward_dynamics_parameter_gradient"),
    AlgoDescriptor("kinetic_energy_regressor"),
    AlgoDescriptor("potential_energy_regressor"),

    # Kinematics
    AlgoDescriptor("end_effector_pose", autotune_keys=("ee_pose",)),
    AlgoDescriptor("end_effector_pose_gradient", autotune_keys=("ee_pose_gradient",)),
    AlgoDescriptor("end_effector_pose_hessian", autotune_keys=("ee_pose_hessian",),
                   gate_attr="generate_end_effector_pose_hessian"),
    AlgoDescriptor("frame_jacobian"),
    AlgoDescriptor("frame_jacobian_dot"),
    AlgoDescriptor("osc_inertia"),
    AlgoDescriptor("end_effector_pose_runtime"),
    AlgoDescriptor("end_effector_pose_gradient_runtime"),

    # Second-Order
    AlgoDescriptor("idsva_so", autotune_keys=("idsva_so",), has_kernel_attr=False),
    AlgoDescriptor("idsva_so_body_frame", autotune_keys=("idsva_so_body_frame",),
                   gate_attr="generate_idsva_so_body_frame"),
    AlgoDescriptor("idsva_so_world_frame", autotune_keys=("idsva_so_world_frame",),
                   gate_attr="generate_idsva_so_world_frame"),
    AlgoDescriptor("fdsva_so", autotune_keys=("fdsva_so",), gate_attr="generate_fdsva_so"),

    # Integrators — kept AFTER Second-Order so `[d for d in ALGO_DESCRIPTORS if
    # d.carries_launch_cfg]` reproduces the launch-config enum order (integrator
    # group LAST), which Step 1 relies on for a byte-identical grid.cuh. (The
    # ALGO_REGISTRY order above is the report/section order and is unaffected;
    # the descriptor test asserts key-equality as a set/dict, order-independent.)
    AlgoDescriptor("integrator", autotune_keys=("integrator",)),
    AlgoDescriptor("integrator_gradient", autotune_keys=("integrator_gradient",),
                   bytes_macro_stem="INTEGRATOR_DU_DYNAMIC_SHARED_MEM_BYTES"),
    AlgoDescriptor("integrator_with_gradient", autotune_keys=("integrator_with_gradient",),
                   bytes_macro_stem="INTEGRATOR_DU_DYNAMIC_SHARED_MEM_BYTES"),
    AlgoDescriptor("integrator_hessian", has_kernel_attr=False),

    # Centroidal / Energy / CoM
    AlgoDescriptor("generalized_gravity",
                   bytes_macro_stem="INVERSE_DYNAMICS_BIAS_DYNAMIC_SHARED_MEM_BYTES"),
    AlgoDescriptor("nonlinear_effects",
                   bytes_macro_stem="INVERSE_DYNAMICS_BIAS_DYNAMIC_SHARED_MEM_BYTES"),
    AlgoDescriptor("energy"),
    AlgoDescriptor("com"),
    AlgoDescriptor("ccrba"),
    AlgoDescriptor("coriolis_matrix"),
    AlgoDescriptor("dccrba"),
    AlgoDescriptor("cmm_time_variation"),

    # Plant (no standalone kernel)
    AlgoDescriptor("plant", has_kernel_attr=False),
    # Collision (no standalone benchmarked kernel; config_free lives in grid_collision)
    AlgoDescriptor("collision", has_kernel_attr=False),
)


_BY_KEY_DESCRIPTOR = {d.key: d for d in ALGO_DESCRIPTORS}


def descriptor_for(key: str) -> AlgoDescriptor:
    try:
        return _BY_KEY_DESCRIPTOR[key]
    except KeyError:
        raise KeyError(f"{key!r} not in ALGO_DESCRIPTORS — add it to algo_registry.py") from None


def build_launch_config_algo_to_symbol() -> dict[str, str]:
    """Reconstruct LAUNCH_CONFIG_ALGO_TO_SYMBOL (bench JSON key -> grid symbol)
    from the descriptor rows. Parity-checked against the live dict in
    test/test_algo_descriptor_parity.py; the eventual Step-1 generation will
    consume this directly."""
    return {k: d.key for d in ALGO_DESCRIPTORS for k in d.autotune_keys}


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
