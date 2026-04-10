"""Compile a PBPK topology + compound into an ODE right-hand side function.

The ``build_compound_pbpk_params`` function assembles every numeric
parameter the ODE needs (Kp per tissue, fu_b, CLint_liver_L_h, CL_renal)
via :class:`charon.core.parameter_bridge.ParameterBridge` — so the full
IVIVE audit trail is still produced by the bridge on every Pipeline run.

``build_rhs`` closes over the topology + params and returns a
``rhs(t, y) -> dy`` function suitable for :func:`scipy.integrate.solve_ivp`.

Units (all interfaces):
  * amounts A : mg
  * volumes V : L
  * concentrations C : mg/L (equivalent to ng/uL)
  * flows Q : L/h
  * clearances : L/h
  * time : h
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np

from charon.core.parameter_bridge import ParameterBridge
from charon.core.schema import CompoundConfig
from charon.pbpk.kp_calculator import compute_all_kp
from charon.pbpk.topology import PBPKTopology, PORTAL_TISSUES

_COMPOUND_TYPES = ("neutral", "acid", "base", "zwitterion")


@dataclass(frozen=True)
class CompoundPBPKParams:
    """Resolved, PBPK-ready compound parameters for a specific topology."""

    name: str
    molecular_weight: float
    logp: float
    pka_acid: float | None
    pka_base: float | None
    compound_type: str
    fu_p: float
    bp_ratio: float
    fu_b: float
    clint_liver_L_h: float
    cl_renal_L_h: float
    kp_by_tissue: dict[str, float]


def infer_compound_type(pka_acid: float | None, pka_base: float | None) -> str:
    """Classify a compound from its pKa values (matches Sisyphus thresholds).

    Rules:
      - acidic pKa < 7.0 AND basic pKa > 8.0  -> zwitterion
      - acidic pKa < 7.0 only                  -> acid
      - basic pKa > 8.0 only                   -> base
      - otherwise                               -> neutral
    """
    is_acidic = pka_acid is not None and pka_acid < 7.0
    is_basic = pka_base is not None and pka_base > 8.0
    if is_acidic and is_basic:
        return "zwitterion"
    if is_acidic:
        return "acid"
    if is_basic:
        return "base"
    return "neutral"


def _get_property_value(prop) -> float | None:
    if prop is None:
        return None
    return float(prop.value)


def _require(value: float | None, label: str) -> float:
    if value is None:
        raise ValueError(
            f"{label} is required to build PBPK parameters but was not provided"
        )
    return value


def build_compound_pbpk_params(
    compound: CompoundConfig,
    topology: PBPKTopology,
    bridge: ParameterBridge,
    *,
    compound_type: str | None = None,
    clint_system: Literal["HLM", "hepatocytes"] = "HLM",
    liver_model: str = "well_stirred",
    override_cl_renal_L_h: float | None = None,
) -> CompoundPBPKParams:
    """Assemble all PBPK ODE parameters for a given compound and species.

    Parameters
    ----------
    compound : CompoundConfig
        Validated compound with at minimum: logP, fu_p, fu_inc, BP ratio,
        CLint and molecular weight.
    topology : PBPKTopology
        Loaded species topology (tissue volumes, flows, compositions).
    bridge : ParameterBridge
        Active IVIVE bridge instance; used to run CLint->CLh and derive
        CLint_liver_L_h with a full ConversionLog audit.
    compound_type : str, optional
        Override the pKa-based compound_type classification.  Useful for
        weak bases (e.g. midazolam, pKa_base=6.2) where the Sisyphus
        threshold of >8.0 would misclassify the drug as neutral.
    clint_system : "HLM" or "hepatocytes"
        Which in-vitro system the CLint value comes from.  HLM applies
        fu_inc correction; hepatocytes do not.
    liver_model : str
        Liver model name passed to the bridge (for the audit trail only —
        the PBPK ODE consumes CLint_liver_L_h regardless).
    override_cl_renal_L_h : float, optional
        If supplied, use this value instead of running the bridge's renal
        clearance estimator.  Enables compounds with experimental renal CL.
    """
    props = compound.properties
    logp = _require(_get_property_value(props.physicochemical.logp), "logp")
    pka_acid = _get_property_value(props.physicochemical.pka_acid)
    pka_base = _get_property_value(props.physicochemical.pka_base)
    fu_p = _require(_get_property_value(props.binding.fu_p), "fu_p")
    fu_inc = _require(_get_property_value(props.binding.fu_inc), "fu_inc")
    bp_ratio = _require(_get_property_value(props.binding.bp_ratio), "bp_ratio")
    clint = _require(
        _get_property_value(props.metabolism.clint_uL_min_mg), "clint_uL_min_mg"
    )
    mw = compound.molecular_weight
    if mw is None:
        raise ValueError(
            f"compound {compound.name!r} has no molecular_weight; required for PBPK"
        )

    resolved_type = compound_type or infer_compound_type(pka_acid, pka_base)
    if resolved_type not in _COMPOUND_TYPES:
        raise ValueError(
            f"compound_type must be one of {_COMPOUND_TYPES}, got {resolved_type!r}"
        )

    pka_for_kp: float | None
    if resolved_type == "acid":
        pka_for_kp = pka_acid
    elif resolved_type in ("base", "zwitterion"):
        pka_for_kp = pka_base
    else:
        pka_for_kp = None

    tissue_comps = {
        name: node.composition for name, node in topology.tissues.items()
    }

    kp_by_tissue = compute_all_kp(
        logp=logp,
        pka=pka_for_kp,
        compound_type=resolved_type,
        tissue_compositions=tissue_comps,
        plasma_composition=topology.plasma_composition,
        method="rodgers_rowland",
    )

    # Hepatic IVIVE with full audit via ParameterBridge.
    hep = bridge.clint_to_clh(
        clint=clint,
        fu_inc=fu_inc,
        fu_p=fu_p,
        bp_ratio=bp_ratio,
        system=clint_system,
        mppgl=40.0,
        hepatocellularity=120.0,
        liver_weight_g=topology.liver_weight_g,
        qh_L_h=topology.tissues["liver"].blood_flow_L_h,
        model=liver_model,
    )
    clint_liver_L_h = hep.clint_liver_L_h
    assert clint_liver_L_h is not None

    if override_cl_renal_L_h is not None:
        cl_renal_L_h = float(override_cl_renal_L_h)
    else:
        cl_renal_prop = _get_property_value(props.renal.clrenal_L_h)
        if cl_renal_prop is not None:
            cl_renal_L_h = cl_renal_prop
        else:
            cl_renal_L_h = bridge.assign_renal_clearance(
                fu_p=fu_p,
                gfr_mL_min=topology.gfr_mL_min,
            )

    fu_b = fu_p / bp_ratio

    return CompoundPBPKParams(
        name=compound.name,
        molecular_weight=float(mw),
        logp=logp,
        pka_acid=pka_acid,
        pka_base=pka_base,
        compound_type=resolved_type,
        fu_p=fu_p,
        bp_ratio=bp_ratio,
        fu_b=fu_b,
        clint_liver_L_h=clint_liver_L_h,
        cl_renal_L_h=cl_renal_L_h,
        kp_by_tissue=dict(kp_by_tissue),
    )


def build_rhs(
    topology: PBPKTopology,
    params: CompoundPBPKParams,
    *,
    infusion_rate_mg_per_h: float = 0.0,
    infusion_duration_h: float = 0.0,
):
    """Return a closure ``rhs(t, y) -> dy`` for the PBPK ODE system.

    State vector layout (length 2 + N_tissues):
      y[0] = A_venous      (mg)
      y[1] = A_arterial    (mg)
      y[2:] = A_tissue_i   (mg) in topology.tissue_names() order

    The returned function is a plain Python callable with no shared state,
    safe to pass to :func:`scipy.integrate.solve_ivp`.  The Jacobian is
    left for scipy to estimate via finite differences; for 17 states the
    BDF setup cost is negligible.

    Parameters
    ----------
    infusion_rate_mg_per_h : float
        Source term added to the venous pool while ``t <
        infusion_duration_h``.  Zero for IV bolus.
    infusion_duration_h : float
        End of the infusion in hours.  Only used when
        ``infusion_rate_mg_per_h > 0``.
    """
    tissue_names = topology.tissue_names()
    n_tissues = len(tissue_names)

    volumes = np.array(
        [topology.tissues[name].volume_L for name in tissue_names],
        dtype=np.float64,
    )
    flows = np.array(
        [topology.tissues[name].blood_flow_L_h for name in tissue_names],
        dtype=np.float64,
    )
    kp = np.array(
        [params.kp_by_tissue[name] for name in tissue_names],
        dtype=np.float64,
    )

    drains_to = [topology.tissues[name].drains_to for name in tissue_names]

    lung_idx = tissue_names.index("lung")
    liver_idx = tissue_names.index("liver")
    kidney_idx = tissue_names.index("kidney")
    portal_indices = np.array(
        [tissue_names.index(p) for p in PORTAL_TISSUES], dtype=np.int64
    )

    bp = params.bp_ratio
    fu_b = params.fu_b
    clint_liver = params.clint_liver_L_h
    cl_renal = params.cl_renal_L_h

    v_ven = topology.venous_volume_L
    v_art = topology.arterial_volume_L
    q_co = topology.cardiac_output_L_h
    q_ha = topology.hepatic_artery_L_h
    q_liver_total = flows[liver_idx]

    def rhs(t: float, y: np.ndarray) -> np.ndarray:
        dy = np.empty_like(y)

        a_ven = y[0]
        a_art = y[1]
        a_tissues = y[2:]

        c_ven = a_ven / v_ven
        c_art = a_art / v_art

        c_tissue = a_tissues / volumes
        c_blood_out = c_tissue * bp / kp

        dy_tissues = np.zeros_like(a_tissues)
        for i in range(n_tissues):
            if i == lung_idx or i == liver_idx:
                continue
            dy_tissues[i] = flows[i] * (c_art - c_blood_out[i])

        dy_tissues[lung_idx] = q_co * (c_ven - c_blood_out[lung_idx])

        dy_tissues[kidney_idx] -= cl_renal * c_blood_out[kidney_idx] / bp

        portal_inflow = 0.0
        for pi in portal_indices:
            portal_inflow += flows[pi] * c_blood_out[pi]
        hepatic_elim = clint_liver * fu_b * c_blood_out[liver_idx]
        dy_tissues[liver_idx] = (
            q_ha * c_art
            + portal_inflow
            - q_liver_total * c_blood_out[liver_idx]
            - hepatic_elim
        )

        venous_inflow = 0.0
        for i in range(n_tissues):
            if i == lung_idx:
                continue
            if drains_to[i] == "venous":
                if i == liver_idx:
                    venous_inflow += q_liver_total * c_blood_out[liver_idx]
                else:
                    venous_inflow += flows[i] * c_blood_out[i]
        dy_ven = venous_inflow - q_co * c_ven

        if infusion_rate_mg_per_h > 0.0 and t < infusion_duration_h:
            dy_ven += infusion_rate_mg_per_h

        dy_art = q_co * c_blood_out[lung_idx] - q_co * c_art

        dy[0] = dy_ven
        dy[1] = dy_art
        dy[2:] = dy_tissues
        return dy

    return rhs
