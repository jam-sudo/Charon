"""Renal clearance estimation (Layer 1 wrapper).

Delegates the actual calculation to
``charon.core.parameter_bridge.ParameterBridge.assign_renal_clearance``.
This module exists solely for API symmetry with the other ``predict/``
modules so that ``predict.__init__`` can expose a flat set of
``predict_*`` functions.
"""

from __future__ import annotations

from charon.core.parameter_bridge import ParameterBridge
from charon.core.units import HUMAN_GFR_ML_MIN

_bridge = ParameterBridge()


def estimate_renal_clearance(
    fu_p: float,
    gfr_mL_min: float = HUMAN_GFR_ML_MIN,
    is_active_secretion: bool = False,
    net_secretion_factor: float = 0.0,
) -> float:
    """Estimate renal clearance in L/h.

    CL_renal = fu_p × GFR_L_h × (1 + net_secretion_factor)

    Args:
        fu_p: Plasma unbound fraction, in [0, 1].
        gfr_mL_min: Glomerular filtration rate in mL/min.
        is_active_secretion: Set ``True`` to annotate that the compound
            is thought to undergo active tubular secretion. The flag is
            advisory — quantitative secretion contribution must come
            from ``net_secretion_factor``.
        net_secretion_factor: Multiplicative secretion factor; 0.0 means
            filtration only. A value of 0.5 means 50% additional
            clearance from secretion.

    Returns:
        Renal clearance in L/h.

    Raises:
        ValueError: If ``fu_p`` is outside [0, 1].
    """
    return _bridge.assign_renal_clearance(
        fu_p=fu_p,
        gfr_mL_min=gfr_mL_min,
        is_active_secretion=is_active_secretion,
        net_secretion_factor=net_secretion_factor,
    )
