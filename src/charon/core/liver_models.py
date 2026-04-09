"""Three hepatic clearance models for IVIVE predictions."""

from __future__ import annotations

import math


def _validate_inputs(qh: float, fu_b: float, clint_liver: float) -> None:
    """Validate common inputs shared by all liver models."""
    if qh <= 0:
        raise ValueError(f"qh must be > 0, got {qh}")
    if fu_b < 0 or fu_b > 1:
        raise ValueError(f"fu_b must be in [0, 1], got {fu_b}")
    if clint_liver < 0:
        raise ValueError(f"clint_liver must be >= 0, got {clint_liver}")


def well_stirred(qh: float, fu_b: float, clint_liver: float) -> float:
    """
    Well-stirred liver model (industry standard for IVIVE).

    CLh = (Qh * fu_b * CLint) / (Qh + fu_b * CLint)

    Most common model.  Assumes instantaneous equilibration, homogeneous
    mixing.  For high CLint: CLh approaches Qh (flow-limited).  For low
    CLint: CLh approaches fu_b * CLint.

    Returns CLh in L/h.
    """
    _validate_inputs(qh, fu_b, clint_liver)

    if clint_liver == 0.0:
        return 0.0

    clh = (qh * fu_b * clint_liver) / (qh + fu_b * clint_liver)
    # Clamp to Qh in case of floating-point overshoot.
    return min(clh, qh)


def parallel_tube(qh: float, fu_b: float, clint_liver: float) -> float:
    """
    Parallel-tube liver model.

    CLh = Qh * (1 - exp(-fu_b * CLint / Qh))

    Assumes plug flow through liver.  More conservative than well-stirred
    for high extraction.  For high CLint: CLh approaches Qh.  For low
    CLint: converges with well-stirred model.

    Returns CLh in L/h.
    """
    _validate_inputs(qh, fu_b, clint_liver)

    if clint_liver == 0.0:
        return 0.0

    clh = qh * (1.0 - math.exp(-fu_b * clint_liver / qh))
    return min(clh, qh)


def dispersion(
    qh: float,
    fu_b: float,
    clint_liver: float,
    dn: float = 0.17,
) -> float:
    """
    Dispersion liver model (Roberts & Rowland, 1986).

    Intermediate between well-stirred and parallel-tube.  DN is the axial
    dispersion number (default 0.17).  For DN->0: approaches parallel-tube.
    For DN->inf: approaches well-stirred.

    Formula:
      RN = fu_b * CLint / Qh
      a  = sqrt(1 + 4 * RN * DN)
      E  = 1 - 4*a / ((1+a)^2 * exp((a-1)/(2*DN))
                       - (1-a)^2 * exp(-(a+1)/(2*DN)))
      CLh = Qh * E

    Returns CLh in L/h.
    """
    _validate_inputs(qh, fu_b, clint_liver)

    if clint_liver == 0.0:
        return 0.0

    rn = fu_b * clint_liver / qh
    a = math.sqrt(1.0 + 4.0 * rn * dn)

    exp_pos_arg = (a - 1.0) / (2.0 * dn)
    exp_neg_arg = -(a + 1.0) / (2.0 * dn)

    # (1+a)^2 * exp((a-1)/(2*DN))  --  this term grows with large CLint
    # (1-a)^2 * exp(-(a+1)/(2*DN)) --  this term vanishes for large CLint
    term_pos = (1.0 + a) ** 2
    term_neg = (1.0 - a) ** 2

    # Handle potential overflow in exp_pos_arg.  Python's math.exp raises
    # OverflowError for very large arguments.  In that limit the extraction
    # ratio E -> 1, so CLh -> Qh.
    try:
        exp_pos = math.exp(exp_pos_arg)
    except OverflowError:
        return qh

    # exp_neg_arg is always <= 0, so exp() is in (0, 1] -- no overflow risk.
    exp_neg = math.exp(exp_neg_arg)

    denominator = term_pos * exp_pos - term_neg * exp_neg

    # Guard against zero denominator (should not happen for valid inputs,
    # but be defensive).
    if denominator == 0.0:
        return qh

    extraction = 1.0 - (4.0 * a) / denominator
    # Clamp extraction to [0, 1].
    extraction = max(0.0, min(1.0, extraction))

    clh = qh * extraction
    return min(clh, qh)


_MODEL_REGISTRY: dict[str, object] = {
    "well_stirred": well_stirred,
    "parallel_tube": parallel_tube,
    "dispersion": dispersion,
}


def get_liver_model(name: str):
    """
    Factory function.  Returns the liver model function by name.

    Valid names: "well_stirred", "parallel_tube", "dispersion"
    Raises ValueError for unknown names.
    """
    try:
        return _MODEL_REGISTRY[name]
    except KeyError:
        valid = ", ".join(sorted(_MODEL_REGISTRY))
        raise ValueError(
            f"Unknown liver model {name!r}. Valid models: {valid}"
        ) from None
