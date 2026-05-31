"""Single source of truth for GRiD algorithm metadata.

Each entry binds:
  - `key`: the JSON output key (lowercase) — also the parser's lowercased label match
  - `display`: the verbose label rendered in benchmark reports
  - `section`: the section heading in the report (Core Dynamics / Gradients / ...)
  - `legacy_labels`: alternate (case-insensitive) labels accepted by the parser,
    for backwards compatibility with older generated binaries

The codegen's printf label is `key.upper()` unless `printf_label` is overridden.
The parser then lowercases the matched line and looks it up against `key` (or one of
the legacy labels).

Consumers:
  - `test/benchmarks/timing_parser.py` builds its label→key maps from this list.
  - `test/benchmarks/generate_report.py` builds ALGO_DISPLAY and ALGO_SECTIONS.
  - (Future) per-algo `gen_*_host` codegen functions read `printf_label` here.

Adding a new algo touches exactly this file (plus the corresponding codegen
emission of the kernel + the bench's `#if GRID_HAS_X` measure wrappers).
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class AlgoEntry:
    key: str
    display: str
    section: str
    legacy_labels: tuple[str, ...] = ()
    printf_label_override: str | None = None

    @property
    def printf_label(self) -> str:
        return self.printf_label_override or self.key.upper()


ALGO_REGISTRY: tuple[AlgoEntry, ...] = (
    # Core Dynamics
    AlgoEntry("id",                   "ID (Inverse Dynamics)",
              "Core Dynamics", legacy_labels=("inverse dynamics",)),
    AlgoEntry("minv",                 "Minv (M⁻¹)",
              "Core Dynamics", legacy_labels=("minv (direct)",)),
    AlgoEntry("fd",                   "FD (Minv+RNEA)",
              "Core Dynamics", legacy_labels=("forward dynamics",)),
    AlgoEntry("aba",                  "ABA (Articulated Body)",
              "Core Dynamics", legacy_labels=("aba (articulated body)",)),
    AlgoEntry("crba",                 "CRBA",
              "Core Dynamics"),

    # Gradients
    AlgoEntry("id_du",                "ID_DU (∂ID/∂q,v)",                  "Gradients"),
    AlgoEntry("fd_du",                "FD_DU (∂FD/∂q,v)",                  "Gradients"),
    AlgoEntry("f_ext_gradient",       "F_EXT_GRAD (∂tau/∂fext=-Jᵀ, ∂q̈/∂fext=M⁻¹Jᵀ)",
              "Gradients"),
    AlgoEntry("f_ext_gradient_dq",    "F_EXT_GRAD_DQ (∂(id_du)/∂fext=-∂Jᵀ/∂q, fixed base)",
              "Gradients"),
    AlgoEntry("regressor",            "REGRESSOR (Joint-torque Y; tau=Y·π, ∂tau/∂π)",
              "Gradients", legacy_labels=("inverse_dynamics_regressor", "joint torque regressor")),
    AlgoEntry("fd_parameter_gradient", "FD_PARAM_GRAD (∂q̈/∂π = -M⁻¹·Y)",
              "Gradients", legacy_labels=("fd_parameter_gradient",)),

    # Integrators
    AlgoEntry("integrator",           "Integrator (x_{k+1})",              "Integrators"),
    AlgoEntry("integrator_gradient",  "Integrator_Gradient (∂x_{k+1}/∂x,u)",
              "Integrators"),
    AlgoEntry("integrator_with_gradient",
              "Integrator_With_Gradient (x_{k+1} + ∂x_{k+1}/∂x,u)",
              "Integrators"),

    # Kinematics
    AlgoEntry("ee_pose",              "EE_POSE",                           "Kinematics",
              legacy_labels=("eepos",)),
    AlgoEntry("ee_pose_gradient",     "EE_POSE_GRADIENT (Jacobian)",       "Kinematics",
              legacy_labels=("deepos",)),
    AlgoEntry("ee_pose_hessian",      "EE_POSE_HESSIAN (2nd-order EE Jacobian)", "Kinematics"),
    AlgoEntry("frame_jacobian",       "FRAME_JACOBIAN (general-frame J: LOCAL/WORLD/LWA)", "Kinematics"),

    # Second-Order
    AlgoEntry("idsva_so",             "IDSVA_SO (dispatched: body for fixed, world for floating)",
              "Second-Order"),
    AlgoEntry("idsva_so_body_frame",  "IDSVA_SO_BODY_FRAME (2nd-order ID, body-frame)",
              "Second-Order", legacy_labels=("id_so",)),
    AlgoEntry("idsva_so_world_frame", "IDSVA_SO_WORLD_FRAME (2nd-order ID, world-frame)",
              "Second-Order"),
    AlgoEntry("fdsva_so",             "FDSVA_SO (2nd-order FD)",
              "Second-Order", legacy_labels=("fd_so",)),

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

    # Plant (T6): cost / constraint / plant-step primitives emitted in the
    # sibling `grid_plant` namespace. No standalone benchmarked kernel — this
    # entry exists so the family has a registry key (display / sectioning) and
    # so future per-primitive bench wrappers can reference it.
    AlgoEntry("plant",                "Plant (cost/constraint/step primitives)",
              "Plant"),
)


def _all_label_forms(entry: AlgoEntry) -> tuple[str, ...]:
    """Every lowercased label string that parses to this entry's key."""
    forms = [entry.printf_label.lower()]
    forms.extend(label.lower() for label in entry.legacy_labels)
    return tuple(forms)


def build_single_label_map() -> dict[str, str]:
    """e.g. {'single call id': 'id', 'single call idsva_so_body_frame': 'idsva_so_body_frame', ...}."""
    out: dict[str, str] = {}
    for entry in ALGO_REGISTRY:
        for form in _all_label_forms(entry):
            out[f"single call {form}"] = entry.key
    return out


def build_batch_with_mem_label_map() -> dict[str, str]:
    """e.g. {'id with memory': 'id', 'fdsva_so with memory': 'fdsva_so', ...}."""
    out: dict[str, str] = {}
    for entry in ALGO_REGISTRY:
        for form in _all_label_forms(entry):
            out[f"{form} with memory"] = entry.key
    return out


def build_batch_compute_only_label_map() -> dict[str, str]:
    """e.g. {'id compute only': 'id', 'fdsva_so compute only': 'fdsva_so', ...}."""
    out: dict[str, str] = {}
    for entry in ALGO_REGISTRY:
        for form in _all_label_forms(entry):
            out[f"{form} compute only"] = entry.key
    return out


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
