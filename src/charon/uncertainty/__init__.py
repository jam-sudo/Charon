"""Charon Layer 4 — uncertainty quantification."""

from charon.uncertainty.dose_range import UncertaintyResult, compute_dose_range
from charon.uncertainty.propagation import PropagationResult, propagate, override_compound
from charon.uncertainty.sampling import (
    SamplingResult,
    build_param_specs,
    generate_lhs_samples,
)
from charon.uncertainty.sobol import compute_sensitivity, compute_sensitivity_with_r2

__all__ = [
    "PropagationResult",
    "SamplingResult",
    "UncertaintyResult",
    "build_param_specs",
    "compute_dose_range",
    "compute_sensitivity",
    "compute_sensitivity_with_r2",
    "generate_lhs_samples",
    "override_compound",
    "propagate",
]
