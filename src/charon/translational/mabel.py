"""MABEL (Minimum Anticipated Biological Effect Level) dose projection.

Two sub-estimates using apparent PK (CL/F for oral, CL for IV):
  1. Single-dose (Cmax): dose = Css_total × Vd_apparent
  2. Steady-state (Css): dose = Css_total × CL_apparent × tau
MRSD = min(single, steady) / safety_factor.

Unit conventions:
    target_kd_nM      : nM
    molecular_weight  : g/mol
    css_u_mg_L        : nM × g/mol × 1e-6  →  mg/L
    css_total_mg_L    : css_u_mg_L / fu_p   →  mg/L
    dose_cmax_mg      : css_total_mg_L × Vd_apparent_L  →  mg
    dose_ss_mg        : css_total_mg_L × CL_apparent_L_h × tau_h  →  mg
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class MABELResult:
    target_kd_nM: float
    target_conc_mg_L: float        # total plasma Css (mg/L)
    fu_p: float
    molecular_weight: float
    dose_cmax_mg: float
    dose_ss_mg: float
    safety_factor: float
    tau_h: float
    mrsd_mg: float
    limiting_approach: str         # "cmax" or "steady_state"


def compute_mabel(
    *,
    target_kd_nM: float,
    molecular_weight: float,
    fu_p: float,
    cl_apparent_L_h: float,
    vd_apparent_L: float,
    safety_factor: float = 10.0,
    tau_h: float = 24.0,
) -> MABELResult:
    """Project MABEL-based FIH dose from target binding affinity and PK.

    Args:
        target_kd_nM: Target binding affinity Kd (nM). Must be > 0.
        molecular_weight: Molecular weight of compound (g/mol).
        fu_p: Plasma unbound fraction (0–1]. Must be > 0.
        cl_apparent_L_h: Apparent clearance CL/F (oral) or CL (IV) in L/h.
        vd_apparent_L: Apparent volume of distribution Vd/F (or Vss) in L.
        safety_factor: Divisor applied to the more conservative dose
            estimate. Default = 10.
        tau_h: Dosing interval in hours. Default = 24 (once daily).

    Returns:
        MABELResult with both dose estimates, the limiting approach, and
        mrsd_mg.

    Hand-calc example (Cmax approach):
        css_u  = 10 nM × 325.77 g/mol × 1e-6  = 3.2577e-3 mg/L
        css_t  = 3.2577e-3 / 0.03             = 0.10859 mg/L
        dose   = 0.10859 × Vd_L

    The lower of dose_cmax and dose_ss is selected as the conservative
    starting point before applying safety_factor.
    """
    if target_kd_nM <= 0:
        raise ValueError(f"target_kd_nM must be > 0, got {target_kd_nM}")
    if fu_p <= 0:
        raise ValueError(f"fu_p must be > 0, got {fu_p}")

    # nM × g/mol × 1e-6 = mg/L
    css_u_mg_L = target_kd_nM * molecular_weight * 1e-6
    css_total_mg_L = css_u_mg_L / fu_p

    dose_cmax = css_total_mg_L * vd_apparent_L
    dose_ss = css_total_mg_L * cl_apparent_L_h * tau_h

    if dose_cmax <= dose_ss:
        limiting = "cmax"
        mrsd = dose_cmax / safety_factor
    else:
        limiting = "steady_state"
        mrsd = dose_ss / safety_factor

    return MABELResult(
        target_kd_nM=target_kd_nM,
        target_conc_mg_L=css_total_mg_L,
        fu_p=fu_p,
        molecular_weight=molecular_weight,
        dose_cmax_mg=dose_cmax,
        dose_ss_mg=dose_ss,
        safety_factor=safety_factor,
        tau_h=tau_h,
        mrsd_mg=mrsd,
        limiting_approach=limiting,
    )
