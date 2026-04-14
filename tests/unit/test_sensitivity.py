"""Unit tests for SRC-based sensitivity analysis (sobol.py).

All tests use synthetic data where the ground truth is known analytically.
"""

from __future__ import annotations

import numpy as np
import pytest

from charon.uncertainty.sobol import compute_sensitivity, compute_sensitivity_with_r2

RNG = np.random.default_rng(0)
N = 1000  # number of samples used in generative tests


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------

def _make_matrix(n: int, p: int, seed: int = 42) -> np.ndarray:
    """Return an (n, p) matrix of i.i.d. standard-normal columns."""
    rng = np.random.default_rng(seed)
    return rng.standard_normal((n, p))


# ---------------------------------------------------------------------------
# Test 1: Single dominant parameter
# ---------------------------------------------------------------------------

def test_single_dominant_param():
    """X0 contributes >>80% when dose = exp(2*X0 + 0.1*X1 + 0.01*X2).

    Hand calculation:
        β_std ≈ [2, 0.1, 0.01]  (all params have σ≈1 before standardisation)
        SRC²₀ = 4 / (4 + 0.01 + 0.0001) ≈ 0.997
        So importance["A"] >> 0.80  ✓
    """
    X = _make_matrix(N, 3)
    doses = np.exp(2.0 * X[:, 0] + 0.1 * X[:, 1] + 0.01 * X[:, 2])
    names = ("A", "B", "C")
    importance = compute_sensitivity(X, doses, names)

    assert importance["A"] > 0.80, (
        f"A should dominate: got {importance['A']:.3f}"
    )


# ---------------------------------------------------------------------------
# Test 2: Equal contribution
# ---------------------------------------------------------------------------

def test_equal_contribution():
    """dose = exp(X0 + X1) → A ≈ B (within 0.15).

    Hand calculation:
        β_std ≈ [1, 1]
        SRC²₀ = SRC²₁ = 0.5  → importance ≈ 0.5 each  ✓
    """
    X = _make_matrix(N, 2)
    doses = np.exp(X[:, 0] + X[:, 1])
    names = ("A", "B")
    importance = compute_sensitivity(X, doses, names)

    assert abs(importance["A"] - importance["B"]) < 0.15, (
        f"A and B should be equal: A={importance['A']:.3f}, B={importance['B']:.3f}"
    )


# ---------------------------------------------------------------------------
# Test 3: All parameter names returned
# ---------------------------------------------------------------------------

def test_returns_all_params():
    """4 input parameters → 4 keys in output dict."""
    X = _make_matrix(N, 4)
    doses = np.exp(X[:, 0] + X[:, 1] + X[:, 2] + X[:, 3])
    names = ("fu_p", "clint_uL_min_mg", "mppgl", "bp_ratio")
    importance = compute_sensitivity(X, doses, names)

    assert set(importance.keys()) == set(names)
    assert len(importance) == 4


# ---------------------------------------------------------------------------
# Test 4: Values sum to approximately 1
# ---------------------------------------------------------------------------

def test_values_sum_near_one():
    """Normalised SRC² indices must sum to ≈1 (within 0.05).

    For any linear model in log-dose space the normalisation is exact;
    numerical noise from lstsq is negligible.
    """
    X = _make_matrix(N, 3)
    doses = np.exp(0.5 * X[:, 0] + 1.5 * X[:, 1] + 0.3 * X[:, 2])
    names = ("alpha", "beta", "gamma")
    importance = compute_sensitivity(X, doses, names)

    total = sum(importance.values())
    assert abs(total - 1.0) < 0.05, f"Sum of importances should be 1.0, got {total:.4f}"


# ---------------------------------------------------------------------------
# Test 5: R² returned and meaningful
# ---------------------------------------------------------------------------

def test_r_squared_returned():
    """For a linear log-dose model R² should be close to 1 (>0.5).

    With N=1000 i.i.d. normal X and dose = exp(β·X) the in-sample R²
    of the standardised regression should be very high (close to 1).
    """
    X = _make_matrix(N, 2)
    doses = np.exp(2.0 * X[:, 0] + 1.0 * X[:, 1])
    names = ("p1", "p2")
    importance, r2 = compute_sensitivity_with_r2(X, doses, names)

    assert r2 > 0.5, f"R² should be > 0.5 for a linear model, got {r2:.4f}"
    # Sanity: importance keys match
    assert set(importance.keys()) == {"p1", "p2"}
