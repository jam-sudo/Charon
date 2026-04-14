"""Unit tests for dose CI aggregation (dose_range.py).

All numerical claims are verified with hand calculations in the docstrings.
"""

from __future__ import annotations

import numpy as np
import pytest

from charon.uncertainty.dose_range import UncertaintyResult, compute_dose_range


# ---------------------------------------------------------------------------
# Fixtures / helpers
# ---------------------------------------------------------------------------

def _simple_sensitivity() -> dict[str, float]:
    return {"fu_p": 0.70, "clint_uL_min_mg": 0.20, "mppgl": 0.10}


def _make_result(doses, *, sensitivity=None, param_names=None, parameter_matrix=None):
    if sensitivity is None:
        sensitivity = _simple_sensitivity()
    if param_names is None:
        param_names = tuple(sensitivity.keys())
    return compute_dose_range(
        np.array(doses, dtype=float),
        sensitivity=sensitivity,
        param_names=param_names,
        parameter_matrix=parameter_matrix,
    )


# ---------------------------------------------------------------------------
# Test 1: Geometric mean point estimate
# ---------------------------------------------------------------------------

def test_geometric_mean():
    """doses = [10, 100] → geomean ≈ 31.62.

    Hand calculation:
        log_doses = [ln(10), ln(100)] = [2.3026, 4.6052]
        mean_log   = (2.3026 + 4.6052) / 2 = 3.4539
        geomean    = exp(3.4539) = 31.623  ✓
    """
    result = _make_result([10.0, 100.0])
    assert abs(result.point_estimate_mg - 31.623) < 0.01, (
        f"Geometric mean should be ≈31.62, got {result.point_estimate_mg:.4f}"
    )


# ---------------------------------------------------------------------------
# Test 2: CI ordering
# ---------------------------------------------------------------------------

def test_ci_ordering():
    """lower < point_estimate < upper for a typical multi-sample distribution."""
    rng = np.random.default_rng(7)
    doses = rng.lognormal(mean=np.log(50), sigma=0.8, size=200)
    result = _make_result(doses)

    assert result.ci_90_lower_mg < result.point_estimate_mg, (
        "5th pct should be below the geometric mean"
    )
    assert result.point_estimate_mg < result.ci_90_upper_mg, (
        "Geometric mean should be below 95th pct"
    )


# ---------------------------------------------------------------------------
# Test 3: Confidence = HIGH for narrow CI
# ---------------------------------------------------------------------------

def test_confidence_high():
    """A tight lognormal distribution (sigma=0.1) gives ci_ratio < 3 → HIGH.

    Hand estimation:
        sigma=0.1 → 90% CI ≈ exp(±1.645*0.1) = exp(±0.165)
        ratio ≈ exp(0.33) ≈ 1.39 << 3  → HIGH  ✓
    """
    rng = np.random.default_rng(3)
    doses = rng.lognormal(mean=np.log(30), sigma=0.10, size=500)
    result = _make_result(doses)

    assert result.confidence == "HIGH", (
        f"Expected HIGH confidence, got {result.confidence} "
        f"(ci_ratio={result.ci_ratio:.2f})"
    )


# ---------------------------------------------------------------------------
# Test 4: Confidence = LOW for wide CI
# ---------------------------------------------------------------------------

def test_confidence_low():
    """A very wide lognormal distribution (sigma=2.0) gives ci_ratio ≥ 10 → LOW.

    Hand estimation:
        sigma=2.0 → 90% CI ≈ exp(±1.645*2.0) = exp(±3.29)
        ratio ≈ exp(6.58) ≈ 718 >> 10  → LOW  ✓
    """
    rng = np.random.default_rng(5)
    doses = rng.lognormal(mean=np.log(50), sigma=2.0, size=500)
    result = _make_result(doses)

    assert result.confidence == "LOW", (
        f"Expected LOW confidence, got {result.confidence} "
        f"(ci_ratio={result.ci_ratio:.2f})"
    )


# ---------------------------------------------------------------------------
# Test 5: Recommendation mentions the top parameter
# ---------------------------------------------------------------------------

def test_recommendation_string():
    """When fu_p has the highest sensitivity, it appears in the recommendation."""
    sensitivity = {"fu_p": 0.70, "clint_uL_min_mg": 0.20, "mppgl": 0.10}
    rng = np.random.default_rng(9)
    doses = rng.lognormal(mean=np.log(20), sigma=0.5, size=100)
    result = _make_result(doses, sensitivity=sensitivity)

    assert "fu_p" in result.recommendation, (
        f"'fu_p' should appear in recommendation: {result.recommendation!r}"
    )
    assert result.limiting_parameter == "fu_p"


# ---------------------------------------------------------------------------
# Test 6: n_successful counts only positive doses
# ---------------------------------------------------------------------------

def test_n_successful():
    """3 positive doses (negatives filtered) → n_successful=3, n_samples=5."""
    doses = np.array([10.0, 20.0, -5.0, 0.0, 50.0])
    result = _make_result(doses)

    assert result.n_successful == 3, (
        f"Only 3 positive values → n_successful=3, got {result.n_successful}"
    )
    assert result.n_samples == 5, (
        f"Total samples (including invalid) should be 5, got {result.n_samples}"
    )


# ---------------------------------------------------------------------------
# Test 7: Raises on all-invalid doses
# ---------------------------------------------------------------------------

def test_raises_on_no_valid_doses():
    """ValueError when all doses are non-positive."""
    doses = np.array([-1.0, 0.0, -2.0])
    with pytest.raises(ValueError, match="No valid"):
        _make_result(doses)
