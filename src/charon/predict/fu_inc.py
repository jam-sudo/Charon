"""Microsomal unbound fraction (fu_inc) estimation.

Uses the Austin et al. (2002) correlation, which is equivalent to the
Hallifax & Houston (2006) parameterization at 1 mg/mL microsomal protein:

    fu_inc = 1 / (1 + 10^(0.072 * logP^2 + 0.067 * logP - 1.126))

References:
    Austin RP et al. (2002). "The influence of nonspecific microsomal
    binding on apparent intrinsic clearance, and its prediction from
    physicochemical properties." Drug Metab Dispos 30(12):1497-1503.

    Hallifax D, Houston JB (2006). "Binding of drugs to hepatic
    microsomes: comment and assessment of current prediction methodology
    with recommendation for improvement." Drug Metab Dispos 34(4):724-726.

CRITICAL DOMAIN RULE:
    fu_inc correction applies to HLM (human liver microsomes) only. It
    MUST NOT be applied to hepatocyte CLint values, which already reflect
    cellular (intact) binding. Applying fu_inc to hepatocyte data will
    overpredict clearance.
"""

from __future__ import annotations

import math


def predict_fu_inc(logp: float) -> float:
    """Predict the microsomal unbound fraction from logP.

    Args:
        logp: Octanol-water partition coefficient (dimensionless).

    Returns:
        fu_inc in (0, 1]. Lower values indicate more non-specific
        binding, typical of lipophilic compounds.

    Raises:
        ValueError: If ``logp`` is not finite.

    Notes:
        Derived for 1 mg/mL HLM incubations. Do not apply to hepatocyte
        intrinsic clearance data. The correlation has an upper bound of
        1.0 and asymptotically approaches 0.0 for logP → +∞.
    """
    if not math.isfinite(logp):
        raise ValueError(f"logp must be finite, got {logp}")
    log_ratio = 0.072 * (logp ** 2) + 0.067 * logp - 1.126
    return 1.0 / (1.0 + 10.0 ** log_ratio)
