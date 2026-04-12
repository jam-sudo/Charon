"""PAD (Pharmacologically Active Dose) projection.

    dose = Ceff × tau × CL_apparent
    MRSD = dose / safety_factor

Unit conventions:
    target_ceff_nM    : nM
    molecular_weight  : g/mol
    ceff_mg_L         : nM × g/mol × 1e-6  →  mg/L
    auc_target        : ceff_mg_L × tau_h  →  mg·h/L
    dose_mg           : auc_target × cl_apparent_L_h  →  mg
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class PADResult:
    target_ceff_nM: float
    target_conc_mg_L: float        # ceff in mg/L
    auc_target_mg_h_L: float       # target AUC over tau
    cl_apparent_L_h: float
    dose_mg: float
    safety_factor: float
    tau_h: float
    mrsd_mg: float


def compute_pad(
    *,
    target_ceff_nM: float,
    molecular_weight: float,
    cl_apparent_L_h: float,
    safety_factor: float = 10.0,
    tau_h: float = 24.0,
) -> PADResult:
    """Project PAD-based FIH dose from efficacious concentration and PK.

    Args:
        target_ceff_nM: Efficacious plasma concentration Ceff (nM). Must
            be > 0. Typically EC50 or IC50 from pharmacology.
        molecular_weight: Molecular weight (g/mol).
        cl_apparent_L_h: Apparent clearance CL/F (oral) or CL (IV) in
            L/h. Must be > 0.
        safety_factor: Divisor applied to the projected dose. Default = 10.
        tau_h: Dosing interval in hours. Default = 24 (once daily).

    Returns:
        PADResult with dose_mg, mrsd_mg, and intermediate quantities.

    Hand-calc (tau=24):
        ceff   = 100 nM × 325.77 g/mol × 1e-6  = 0.032577 mg/L
        AUC    = 0.032577 × 24                 = 0.78185 mg·h/L
        dose   = 0.78185 × 18.6               = 14.542 mg
        MRSD   = 14.542 / 10                  = 1.4542 mg
    """
    if target_ceff_nM <= 0:
        raise ValueError(f"target_ceff_nM must be > 0, got {target_ceff_nM}")
    if cl_apparent_L_h <= 0:
        raise ValueError(f"cl_apparent_L_h must be > 0, got {cl_apparent_L_h}")

    # nM × g/mol × 1e-6 = mg/L
    ceff_mg_L = target_ceff_nM * molecular_weight * 1e-6
    auc_target = ceff_mg_L * tau_h
    dose = auc_target * cl_apparent_L_h
    mrsd = dose / safety_factor

    return PADResult(
        target_ceff_nM=target_ceff_nM,
        target_conc_mg_L=ceff_mg_L,
        auc_target_mg_h_L=auc_target,
        cl_apparent_L_h=cl_apparent_L_h,
        dose_mg=dose,
        safety_factor=safety_factor,
        tau_h=tau_h,
        mrsd_mg=mrsd,
    )
