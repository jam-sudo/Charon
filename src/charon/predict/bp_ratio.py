"""Blood:plasma ratio (B:P) estimation.

Design note (Sprint 2): Charon intentionally does NOT train an ML RBP
model. Available RBP data (``adme_reference.csv`` 153 rows) is all in
the conformal calibration set, and using it for training would leak
the holdout. Omega's XGBoost RBP baseline also reported R² ≈ -0.08 on
scaffold CV, which is worse than a constant-mean predictor. Therefore
Charon uses an empirical formula grounded in Hct physiology.

Formula (Poulin & Theil 2002 simplification):
    B:P = 1 + Hct × (fu_p × Kp_rbc − 1)

where
    Hct   = 0.45 (standard human hematocrit)
    Kp_rbc = erythrocyte:plasma water partition, defaulted to 1.0 for
             neutral drugs. Basic drugs (pKa > 8) are slightly favoured
             by erythrocytes (Kp_rbc > 1) and acidic drugs the opposite.

The result is clipped to the physiological range [0.5, 3.0].
"""

from __future__ import annotations

import math

HUMAN_HCT: float = 0.45
BP_MIN: float = 0.5
BP_MAX: float = 3.0

# Compound-type defaults for Kp_rbc (rough physiological priors).
_KP_RBC_DEFAULTS: dict[str, float] = {
    "neutral": 1.0,
    "acid": 0.7,
    "base": 1.3,
    "zwitterion": 1.0,
}


def predict_bp_ratio(
    fu_p: float,
    compound_type: str = "neutral",
    hct: float = HUMAN_HCT,
    kp_rbc: float | None = None,
) -> float:
    """Estimate the blood:plasma concentration ratio from fu_p.

    Args:
        fu_p: Plasma unbound fraction, in [0, 1].
        compound_type: One of ``"neutral"``, ``"acid"``, ``"base"``,
            ``"zwitterion"``. Selects the default ``kp_rbc`` prior.
        hct: Hematocrit. Defaults to the standard adult value 0.45.
        kp_rbc: Override for erythrocyte:plasma water partition. If
            ``None``, the compound-type default is used.

    Returns:
        B:P ratio, clipped to [0.5, 3.0].

    Raises:
        ValueError: If ``fu_p`` is outside [0, 1], ``hct`` is outside
            (0, 1), or inputs are not finite.
    """
    if not math.isfinite(fu_p) or fu_p < 0.0 or fu_p > 1.0:
        raise ValueError(f"fu_p must be in [0, 1], got {fu_p}")
    if not math.isfinite(hct) or hct <= 0.0 or hct >= 1.0:
        raise ValueError(f"hct must be in (0, 1), got {hct}")

    if kp_rbc is None:
        kp_rbc = _KP_RBC_DEFAULTS.get(compound_type, 1.0)
    if not math.isfinite(kp_rbc) or kp_rbc <= 0.0:
        raise ValueError(f"kp_rbc must be positive and finite, got {kp_rbc}")

    bp = 1.0 + hct * (fu_p * kp_rbc - 1.0)
    return float(max(BP_MIN, min(BP_MAX, bp)))
