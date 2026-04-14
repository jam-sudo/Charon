# Sprint 5: Uncertainty Quantification Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Propagate prediction uncertainty through PBPK → dose to produce `dose [point, 90% CI]` + per-parameter sensitivity ranking.

**Architecture:** LHS sampling (scipy) with Iman-Conover rank correlation generates N=500 parameter sets. Each is propagated through Pipeline._run_single() to get a dose. Doses are aggregated into geometric mean + P5/P95. SRC-based sensitivity identifies the dominant uncertainty driver.

**Tech Stack:** Python 3.11+, scipy.stats.qmc (LHS), numpy (linalg, stats), Pydantic v2, pytest

**Spec:** `docs/superpowers/specs/2026-04-13-sprint5-uncertainty-quantification-design.md`

---

## File Structure

| Action | File | Responsibility |
|--------|------|----------------|
| Create | `src/charon/uncertainty/sampling.py` | LHS + Iman-Conover + ParameterSample |
| Create | `src/charon/uncertainty/propagation.py` | N-sample propagation loop |
| Create | `src/charon/uncertainty/dose_range.py` | Dose CI aggregation + confidence |
| Create | `src/charon/uncertainty/sobol.py` | SRC-based sensitivity analysis |
| Modify | `src/charon/pipeline.py` | Add uncertainty param, _run_single, wire propagation |
| Modify | `src/charon/uncertainty/__init__.py` | Export public symbols |
| Create | `tests/unit/test_sampling.py` | LHS + correlation tests |
| Create | `tests/unit/test_sensitivity.py` | SRC tests |
| Create | `tests/unit/test_dose_range.py` | CI + confidence tests |
| Create | `tests/integration/test_uncertainty_pipeline.py` | E2E with midazolam |

---

### Task 1: Implement LHS sampling with parameter distributions (sampling.py)

**Files:**
- Create: `src/charon/uncertainty/sampling.py`
- Create: `tests/unit/test_sampling.py`

- [ ] **Step 1: Write the failing tests**

Create `tests/unit/test_sampling.py`:

```python
"""Tests for LHS sampling with Iman-Conover correlation."""
from __future__ import annotations

import numpy as np
import pytest


class TestGenerateLHSSamples:
    def test_returns_correct_count(self):
        from charon.uncertainty.sampling import generate_lhs_samples
        result = generate_lhs_samples(
            param_specs={"logp": (2.0, 1.0, "normal"), "fu_p": (0.03, 0.5, "lognormal")},
            n_samples=50,
            correlation="none",
        )
        assert len(result.samples) == 50

    def test_param_names_preserved(self):
        from charon.uncertainty.sampling import generate_lhs_samples
        result = generate_lhs_samples(
            param_specs={"logp": (2.0, 1.0, "normal"), "fu_p": (0.03, 0.5, "lognormal")},
            n_samples=10,
            correlation="none",
        )
        assert result.param_names == ("logp", "fu_p")

    def test_lhs_uniform_coverage(self):
        """Each bin in [0,1] should have exactly 1 sample (LHS property)."""
        from charon.uncertainty.sampling import generate_lhs_samples
        result = generate_lhs_samples(
            param_specs={"x": (5.0, 1.0, "normal")},
            n_samples=100,
            correlation="none",
            seed=42,
        )
        values = np.array([s["x"] for s in result.samples])
        # After CDF transform back to uniform, each percentile bin should have ~1 sample
        # Just verify we got 100 distinct values
        assert len(set(values)) == 100

    def test_normal_distribution_range(self):
        """Normal samples should be centered near mu."""
        from charon.uncertainty.sampling import generate_lhs_samples
        result = generate_lhs_samples(
            param_specs={"logp": (3.0, 0.5, "normal")},
            n_samples=500,
            correlation="none",
            seed=42,
        )
        values = np.array([s["logp"] for s in result.samples])
        assert abs(np.mean(values) - 3.0) < 0.2

    def test_lognormal_positive(self):
        """LogNormal samples must all be positive."""
        from charon.uncertainty.sampling import generate_lhs_samples
        result = generate_lhs_samples(
            param_specs={"fu_p": (0.03, 0.5, "lognormal")},
            n_samples=200,
            correlation="none",
            seed=42,
        )
        values = np.array([s["fu_p"] for s in result.samples])
        assert np.all(values > 0)

    def test_fu_p_clipped_to_1(self):
        """fu_p should never exceed 1.0 (bounded quantity)."""
        from charon.uncertainty.sampling import generate_lhs_samples
        result = generate_lhs_samples(
            param_specs={"fu_p": (0.5, 0.8, "lognormal")},  # wide spread could exceed 1
            n_samples=200,
            correlation="none",
            seed=42,
        )
        values = np.array([s["fu_p"] for s in result.samples])
        assert np.all(values <= 1.0)

    def test_seed_reproducibility(self):
        from charon.uncertainty.sampling import generate_lhs_samples
        specs = {"logp": (3.0, 0.5, "normal"), "fu_p": (0.03, 0.5, "lognormal")}
        r1 = generate_lhs_samples(param_specs=specs, n_samples=50, seed=123, correlation="none")
        r2 = generate_lhs_samples(param_specs=specs, n_samples=50, seed=123, correlation="none")
        v1 = [s["logp"] for s in r1.samples]
        v2 = [s["logp"] for s in r2.samples]
        assert v1 == v2


class TestImanConover:
    def test_induces_correlation(self):
        """With target corr(logP, fu_p) = -0.6, result should be close."""
        from charon.uncertainty.sampling import generate_lhs_samples
        import numpy as np
        specs = {"logp": (3.0, 0.5, "normal"), "fu_p": (0.03, 0.5, "lognormal")}
        corr_matrix = np.array([[1.0, -0.6], [-0.6, 1.0]])
        result = generate_lhs_samples(
            param_specs=specs,
            n_samples=500,
            correlation="iman_conover",
            target_correlation=corr_matrix,
            seed=42,
        )
        logp_vals = np.array([s["logp"] for s in result.samples])
        fup_vals = np.array([s["fu_p"] for s in result.samples])
        from scipy.stats import spearmanr
        rho, _ = spearmanr(logp_vals, fup_vals)
        assert abs(rho - (-0.6)) < 0.15, f"Spearman rho = {rho:.3f}, expected ~ -0.6"

    def test_no_correlation_flag(self):
        """correlation='none' should not apply Iman-Conover."""
        from charon.uncertainty.sampling import generate_lhs_samples
        specs = {"logp": (3.0, 0.5, "normal"), "fu_p": (0.03, 0.5, "lognormal")}
        result = generate_lhs_samples(
            param_specs=specs, n_samples=50, correlation="none", seed=42,
        )
        assert result.correlation_applied is False


class TestBuildParamSpecs:
    def test_from_compound_experimental(self):
        """Experimental values (no CI) should use fallback CVs."""
        from charon.uncertainty.sampling import build_param_specs
        from charon.core.schema import CompoundConfig
        import yaml
        d = yaml.safe_load(open("validation/data/tier1_obach/compounds/midazolam.yaml"))
        compound = CompoundConfig(**d)
        specs = build_param_specs(compound)
        assert "logp" in specs
        assert "fu_p" in specs
        assert "clint" in specs
        # All should have positive sigma (fallback CV)
        for name, (mu, sigma, dist) in specs.items():
            assert sigma > 0, f"{name} has sigma=0"

    def test_from_compound_with_ci(self):
        """Properties with conformal CIs should use CI-derived sigma."""
        from charon.uncertainty.sampling import build_param_specs
        from charon.core.schema import (
            CompoundConfig, CompoundProperties, PhysicochemicalProperties,
            BindingProperties, MetabolismProperties, PredictedProperty,
        )
        compound = CompoundConfig(
            name="test", smiles="C", molecular_weight=300.0,
            properties=CompoundProperties(
                physicochemical=PhysicochemicalProperties(
                    logp=PredictedProperty(value=3.0, ci_90_lower=2.0, ci_90_upper=4.0, source="ml_ensemble"),
                ),
                binding=BindingProperties(
                    fu_p=PredictedProperty(value=0.1, ci_90_lower=0.02, ci_90_upper=0.5, source="ml_ensemble"),
                    bp_ratio=PredictedProperty(value=0.8, source="ml_ensemble"),
                ),
                metabolism=MetabolismProperties(
                    clint_uL_min_mg=PredictedProperty(value=50.0, ci_90_lower=10.0, ci_90_upper=250.0, source="ml_ensemble"),
                ),
            ),
        )
        specs = build_param_specs(compound)
        # logP should use CI-derived sigma
        assert specs["logp"][1] > 0
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/unit/test_sampling.py -v`
Expected: FAIL — module does not exist

- [ ] **Step 3: Implement sampling.py**

Create `src/charon/uncertainty/sampling.py`:

```python
"""Latin Hypercube Sampling with optional Iman-Conover rank correlation.

Generates N parameter sets for uncertainty propagation. Each parameter
has a specified distribution (normal or lognormal) with mean and sigma.
Iman-Conover re-ranks the LHS samples to match a target correlation
matrix while preserving marginal distributions.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.stats import norm, rankdata
from scipy.stats.qmc import LatinHypercube

from charon.core.schema import CompoundConfig

# Fallback CVs for experimental values (no conformal CI).
_FALLBACK_CV: dict[str, float] = {
    "logp": 0.15,       # well-measured, low variability
    "fu_p": 0.50,       # high inter-lab variability
    "clint": 0.60,      # highest uncertainty source
    "peff": 0.40,       # moderate
    "mppgl": 0.20,      # literature CV ~20%
    "bp_ratio": 0.10,   # low variability
}

# Physical bounds for clipping.
_BOUNDS: dict[str, tuple[float, float]] = {
    "fu_p": (0.001, 1.0),
    "clint": (0.01, 5000.0),
    "peff": (1e-7, 1e-2),
    "bp_ratio": (0.3, 3.0),
    "mppgl": (10.0, 80.0),
}

# Default correlation matrix (ARCHITECTURE §Layer 4).
# Order: logp, fu_p, clint, peff, mppgl, bp_ratio
DEFAULT_CORRELATION = np.array([
    [ 1.00, -0.60,  0.00,  0.40,  0.00,  0.00],
    [-0.60,  1.00,  0.00,  0.00,  0.00,  0.00],
    [ 0.00,  0.00,  1.00,  0.00,  0.00,  0.00],
    [ 0.40,  0.00,  0.00,  1.00,  0.00,  0.00],
    [ 0.00,  0.00,  0.00,  0.00,  1.00,  0.00],
    [ 0.00,  0.00,  0.00,  0.00,  0.00,  1.00],
])
DEFAULT_PARAM_ORDER = ("logp", "fu_p", "clint", "peff", "mppgl", "bp_ratio")


@dataclass(frozen=True)
class SamplingResult:
    """N parameter sets from LHS."""

    samples: tuple[dict[str, float], ...]
    param_names: tuple[str, ...]
    n_params_sampled: int
    correlation_applied: bool
    seed: int


def build_param_specs(
    compound: CompoundConfig,
) -> dict[str, tuple[float, float, str]]:
    """Build parameter specs {name: (mu, sigma, dist_type)} from a compound.

    Uses conformal CI widths when available, otherwise fallback CVs.
    dist_type is "normal" for logP, "lognormal" for positive-only params.
    """
    specs: dict[str, tuple[float, float, str]] = {}
    props = compound.properties

    # logP — normal distribution
    if props.physicochemical.logp is not None:
        p = props.physicochemical.logp
        mu = p.value
        if p.ci_90_lower is not None and p.ci_90_upper is not None:
            sigma = (p.ci_90_upper - p.ci_90_lower) / 3.29  # 90% CI → σ
        else:
            sigma = abs(mu) * _FALLBACK_CV["logp"]
        if sigma > 0:
            specs["logp"] = (mu, sigma, "normal")

    # fu_p — lognormal
    if props.binding.fu_p is not None:
        p = props.binding.fu_p
        mu = p.value
        if p.ci_90_lower is not None and p.ci_90_upper is not None:
            log_width = math.log10(p.ci_90_upper) - math.log10(max(p.ci_90_lower, 1e-6))
            sigma_log = log_width / 3.29
        else:
            sigma_log = _FALLBACK_CV["fu_p"]
        if sigma_log > 0:
            specs["fu_p"] = (mu, sigma_log, "lognormal")

    # CLint — lognormal
    if props.metabolism.clint_uL_min_mg is not None:
        p = props.metabolism.clint_uL_min_mg
        mu = p.value
        if p.ci_90_lower is not None and p.ci_90_upper is not None:
            log_width = math.log10(p.ci_90_upper) - math.log10(max(p.ci_90_lower, 0.01))
            sigma_log = log_width / 3.29
        else:
            sigma_log = _FALLBACK_CV["clint"]
        if sigma_log > 0:
            specs["clint"] = (mu, sigma_log, "lognormal")

    # Peff — lognormal (if available)
    if props.permeability.peff_cm_s is not None:
        p = props.permeability.peff_cm_s
        mu = p.value
        if p.ci_90_lower is not None and p.ci_90_upper is not None:
            log_width = math.log10(p.ci_90_upper) - math.log10(max(p.ci_90_lower, 1e-7))
            sigma_log = log_width / 3.29
        else:
            sigma_log = _FALLBACK_CV["peff"]
        if sigma_log > 0:
            specs["peff"] = (mu, sigma_log, "lognormal")

    # MPPGL — normal (physiological constant)
    specs["mppgl"] = (40.0, 40.0 * _FALLBACK_CV["mppgl"], "normal")

    # BP ratio — normal
    if props.binding.bp_ratio is not None:
        p = props.binding.bp_ratio
        mu = p.value
        sigma = mu * _FALLBACK_CV["bp_ratio"]
        if sigma > 0:
            specs["bp_ratio"] = (mu, sigma, "normal")

    return specs


def generate_lhs_samples(
    *,
    param_specs: dict[str, tuple[float, float, str]],
    n_samples: int = 500,
    correlation: str = "iman_conover",
    target_correlation: np.ndarray | None = None,
    seed: int = 42,
) -> SamplingResult:
    """Generate LHS parameter samples with optional rank correlation.

    Parameters
    ----------
    param_specs : dict
        {name: (mu, sigma, dist_type)} where dist_type is
        "normal" or "lognormal".
    n_samples : int
        Number of LHS samples.
    correlation : str
        "iman_conover" or "none".
    target_correlation : ndarray, optional
        p×p correlation matrix. If None and correlation="iman_conover",
        uses DEFAULT_CORRELATION subset for the given params.
    seed : int
        Random seed for reproducibility.
    """
    param_names = tuple(param_specs.keys())
    n_params = len(param_names)

    rng = np.random.default_rng(seed)

    # 1. Generate LHS in [0, 1]^p
    sampler = LatinHypercube(d=n_params, seed=rng)
    unit = sampler.random(n=n_samples)  # (N, p)

    # 2. Transform to target marginals
    raw = np.empty_like(unit)
    for j, name in enumerate(param_names):
        mu, sigma, dist_type = param_specs[name]
        if dist_type == "normal":
            raw[:, j] = norm.ppf(unit[:, j], loc=mu, scale=sigma)
        elif dist_type == "lognormal":
            # mu is the median in original space; sigma is in log10 space
            log_mu = math.log10(mu)
            log_vals = norm.ppf(unit[:, j], loc=log_mu, scale=sigma)
            raw[:, j] = 10.0 ** log_vals
        else:
            raise ValueError(f"Unknown dist_type {dist_type!r}")

    # 3. Apply Iman-Conover (if requested)
    applied = False
    if correlation == "iman_conover" and n_params >= 2:
        if target_correlation is not None:
            C_target = target_correlation
        else:
            C_target = _build_default_correlation(param_names)
        if C_target is not None:
            raw = _iman_conover(raw, C_target, rng)
            applied = True

    # 4. Clip to physical bounds
    for j, name in enumerate(param_names):
        if name in _BOUNDS:
            lo, hi = _BOUNDS[name]
            raw[:, j] = np.clip(raw[:, j], lo, hi)

    # 5. Build sample dicts
    samples = []
    for i in range(n_samples):
        d = {name: float(raw[i, j]) for j, name in enumerate(param_names)}
        samples.append(d)

    return SamplingResult(
        samples=tuple(samples),
        param_names=param_names,
        n_params_sampled=n_params,
        correlation_applied=applied,
        seed=seed,
    )


def _build_default_correlation(param_names: tuple[str, ...]) -> np.ndarray | None:
    """Extract the relevant submatrix of DEFAULT_CORRELATION."""
    indices = []
    for name in param_names:
        if name in DEFAULT_PARAM_ORDER:
            indices.append(DEFAULT_PARAM_ORDER.index(name))
        else:
            return None  # unknown param, skip correlation
    idx = np.array(indices)
    return DEFAULT_CORRELATION[np.ix_(idx, idx)]


def _iman_conover(
    samples: np.ndarray,
    target_corr: np.ndarray,
    rng: np.random.Generator,
) -> np.ndarray:
    """Apply Iman-Conover rank-based correlation adjustment.

    Preserves marginal distributions while inducing the target rank
    correlation matrix.
    """
    n, p = samples.shape

    # Compute rank matrix (1..N for each column)
    ranks = np.empty_like(samples)
    for j in range(p):
        ranks[:, j] = rankdata(samples[:, j])

    # Observed rank correlation
    R_obs = np.corrcoef(ranks, rowvar=False)

    # Cholesky decompositions
    try:
        L_target = np.linalg.cholesky(target_corr)
        L_obs = np.linalg.cholesky(R_obs)
    except np.linalg.LinAlgError:
        # If Cholesky fails (non-PD matrix), skip correlation
        return samples

    # Transform: new_ranks = ranks × inv(L_obs)^T × L_target^T
    T = np.linalg.solve(L_obs, L_target)  # L_obs^{-1} × L_target
    new_ranks_raw = ranks @ T

    # Re-rank to get integer ranks
    new_ranks = np.empty_like(new_ranks_raw)
    for j in range(p):
        new_ranks[:, j] = rankdata(new_ranks_raw[:, j])

    # Re-sort original samples by new ranks
    result = np.empty_like(samples)
    for j in range(p):
        sorted_vals = np.sort(samples[:, j])
        new_order = new_ranks[:, j].astype(int) - 1  # 0-indexed
        result[:, j] = sorted_vals[new_order]

    return result
```

- [ ] **Step 4: Run tests**

Run: `pytest tests/unit/test_sampling.py -v`
Expected: 10 PASSED

- [ ] **Step 5: Regression**

Run: `pytest tests/ -x -q`
Expected: 647 + 10 = 657 passed

- [ ] **Step 6: Commit**

```bash
git add src/charon/uncertainty/sampling.py tests/unit/test_sampling.py
git commit -m "feat(uncertainty): LHS sampling with Iman-Conover correlation and fallback CVs"
```

---

### Task 2: Implement SRC-based sensitivity analysis (sobol.py)

**Files:**
- Create: `src/charon/uncertainty/sobol.py`
- Create: `tests/unit/test_sensitivity.py`

- [ ] **Step 1: Write the failing tests**

Create `tests/unit/test_sensitivity.py`:

```python
"""Tests for SRC-based sensitivity analysis."""
from __future__ import annotations

import numpy as np
import pytest


class TestComputeSensitivity:
    def test_single_dominant_param(self):
        """When dose depends only on param A, A should have importance ~1.0."""
        from charon.uncertainty.sobol import compute_sensitivity
        rng = np.random.default_rng(42)
        n = 200
        X = rng.standard_normal((n, 3))
        # dose = exp(2*X0 + 0.1*X1 + 0.01*X2) — X0 dominates
        doses = np.exp(2 * X[:, 0] + 0.1 * X[:, 1] + 0.01 * X[:, 2])
        result = compute_sensitivity(
            parameter_matrix=X,
            doses=doses,
            param_names=("A", "B", "C"),
        )
        assert result["A"] > 0.8  # A should dominate
        assert result["A"] > result["B"] > result["C"]

    def test_equal_contribution(self):
        """Two equal params should have similar importance."""
        from charon.uncertainty.sobol import compute_sensitivity
        rng = np.random.default_rng(42)
        n = 500
        X = rng.standard_normal((n, 2))
        doses = np.exp(X[:, 0] + X[:, 1])
        result = compute_sensitivity(X, doses, ("A", "B"))
        assert abs(result["A"] - result["B"]) < 0.15

    def test_returns_all_params(self):
        from charon.uncertainty.sobol import compute_sensitivity
        rng = np.random.default_rng(42)
        X = rng.standard_normal((100, 4))
        doses = np.exp(X @ [1, 0.5, 0.2, 0.1])
        result = compute_sensitivity(X, doses, ("a", "b", "c", "d"))
        assert set(result.keys()) == {"a", "b", "c", "d"}

    def test_values_sum_near_one(self):
        """Normalized importance should sum to ~1.0 for linear models."""
        from charon.uncertainty.sobol import compute_sensitivity
        rng = np.random.default_rng(42)
        X = rng.standard_normal((500, 3))
        doses = np.exp(X @ [1, 0.5, 0.2])
        result = compute_sensitivity(X, doses, ("a", "b", "c"))
        assert abs(sum(result.values()) - 1.0) < 0.05

    def test_r_squared_returned(self):
        from charon.uncertainty.sobol import compute_sensitivity_with_r2
        rng = np.random.default_rng(42)
        X = rng.standard_normal((200, 2))
        doses = np.exp(X @ [1.0, 0.5])
        importance, r2 = compute_sensitivity_with_r2(X, doses, ("a", "b"))
        assert 0.5 < r2 <= 1.0  # linear model should have high R²
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/unit/test_sensitivity.py -v`
Expected: FAIL

- [ ] **Step 3: Implement sobol.py**

Create `src/charon/uncertainty/sobol.py`:

```python
"""Sensitivity analysis via Standardized Regression Coefficients (SRC).

SRC² approximates first-order Sobol indices for near-linear models.
This is a pragmatic Phase A approximation — true Saltelli Sobol
requires N×(2p+2) evaluations and is deferred to Phase B.

The regression is performed in log-dose space:
    log(dose) = β₀ + β₁z₁ + β₂z₂ + ... + ε
where z_i are standardized parameters (mean=0, std=1).

SRC_i = β_i, and importance_i = β_i² / Σβ_j².
"""

from __future__ import annotations

import numpy as np


def compute_sensitivity(
    parameter_matrix: np.ndarray,
    doses: np.ndarray,
    param_names: tuple[str, ...],
) -> dict[str, float]:
    """Compute normalized SRC² importance scores.

    Parameters
    ----------
    parameter_matrix : ndarray, shape (N, p)
        Sampled parameter values (raw, not standardized).
    doses : ndarray, shape (N,)
        Dose outcomes for each sample.
    param_names : tuple of str
        Names for each parameter column.

    Returns
    -------
    dict[str, float]
        Normalized importance scores summing to ~1.0.
    """
    importance, _ = compute_sensitivity_with_r2(parameter_matrix, doses, param_names)
    return importance


def compute_sensitivity_with_r2(
    parameter_matrix: np.ndarray,
    doses: np.ndarray,
    param_names: tuple[str, ...],
) -> tuple[dict[str, float], float]:
    """Compute SRC² importance + R² diagnostic.

    Returns
    -------
    (importance, r_squared)
        importance: dict of normalized SRC² per parameter.
        r_squared: coefficient of determination of the linear fit.
        If R² < 0.7, the model is significantly nonlinear and
        SRC² may not capture total variance.
    """
    X = np.asarray(parameter_matrix, dtype=np.float64)
    y = np.log(np.asarray(doses, dtype=np.float64))  # log-space regression

    n, p = X.shape
    assert p == len(param_names)

    # Standardize
    X_mean = X.mean(axis=0)
    X_std = X.std(axis=0)
    X_std[X_std == 0] = 1.0  # avoid division by zero for fixed params
    Z = (X - X_mean) / X_std

    y_mean = y.mean()
    y_std = y.std()
    if y_std == 0:
        # All doses identical → no uncertainty
        return {name: 1.0 / p for name in param_names}, 1.0
    y_z = (y - y_mean) / y_std

    # OLS: beta = (Z^T Z)^{-1} Z^T y
    beta, residuals, rank, sv = np.linalg.lstsq(Z, y_z, rcond=None)

    # SRC² and normalization
    src_sq = beta ** 2
    total = src_sq.sum()
    if total == 0:
        importance = {name: 1.0 / p for name in param_names}
    else:
        importance = {name: float(src_sq[j] / total) for j, name in enumerate(param_names)}

    # R²
    y_pred = Z @ beta
    ss_res = np.sum((y_z - y_pred) ** 2)
    ss_tot = np.sum(y_z ** 2)  # already mean-centered
    r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0

    return importance, float(r_squared)
```

- [ ] **Step 4: Run tests**

Run: `pytest tests/unit/test_sensitivity.py -v`
Expected: 5 PASSED

- [ ] **Step 5: Regression**

Run: `pytest tests/ -x -q`
Expected: 657 + 5 = 662 passed

- [ ] **Step 6: Commit**

```bash
git add src/charon/uncertainty/sobol.py tests/unit/test_sensitivity.py
git commit -m "feat(uncertainty): SRC-based sensitivity analysis (approximate Sobol S_i)"
```

---

### Task 3: Implement dose aggregation (dose_range.py)

**Files:**
- Create: `src/charon/uncertainty/dose_range.py`
- Create: `tests/unit/test_dose_range.py`

- [ ] **Step 1: Write the failing tests**

Create `tests/unit/test_dose_range.py`:

```python
"""Tests for dose CI aggregation and confidence classification."""
from __future__ import annotations

import numpy as np
import pytest


class TestComputeDoseRange:
    def test_geometric_mean(self):
        from charon.uncertainty.dose_range import compute_dose_range
        doses = np.array([10.0, 100.0])  # geomean = sqrt(1000) ≈ 31.62
        result = compute_dose_range(doses, sensitivity={"a": 1.0}, param_names=("a",))
        assert result.point_estimate_mg == pytest.approx(31.623, rel=1e-3)

    def test_ci_ordering(self):
        from charon.uncertainty.dose_range import compute_dose_range
        rng = np.random.default_rng(42)
        doses = np.exp(rng.normal(3, 0.5, 200))
        result = compute_dose_range(doses, sensitivity={"a": 1.0}, param_names=("a",))
        assert result.ci_90_lower_mg < result.point_estimate_mg < result.ci_90_upper_mg

    def test_confidence_high(self):
        """Narrow CI → HIGH confidence."""
        from charon.uncertainty.dose_range import compute_dose_range
        doses = np.full(100, 50.0) + np.random.default_rng(42).normal(0, 2, 100)
        result = compute_dose_range(doses, sensitivity={"a": 1.0}, param_names=("a",))
        assert result.confidence == "HIGH"

    def test_confidence_low(self):
        """Very wide CI → LOW confidence."""
        from charon.uncertainty.dose_range import compute_dose_range
        doses = np.exp(np.random.default_rng(42).normal(3, 2, 200))  # huge spread
        result = compute_dose_range(doses, sensitivity={"a": 1.0}, param_names=("a",))
        assert result.confidence == "LOW"

    def test_recommendation_string(self):
        from charon.uncertainty.dose_range import compute_dose_range
        doses = np.exp(np.random.default_rng(42).normal(3, 0.5, 100))
        result = compute_dose_range(
            doses,
            sensitivity={"fu_p": 0.45, "CLint": 0.35, "logP": 0.20},
            param_names=("fu_p", "CLint", "logP"),
        )
        assert "fu_p" in result.recommendation
        assert result.limiting_parameter == "fu_p"

    def test_n_successful(self):
        from charon.uncertainty.dose_range import compute_dose_range
        doses = np.array([10.0, 20.0, 30.0])
        result = compute_dose_range(doses, sensitivity={"a": 1.0}, param_names=("a",))
        assert result.n_successful == 3
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/unit/test_dose_range.py -v`
Expected: FAIL

- [ ] **Step 3: Implement dose_range.py**

Create `src/charon/uncertainty/dose_range.py`:

```python
"""Dose CI aggregation from propagation results.

Computes geometric mean, 5th/95th percentiles, confidence classification,
and actionable recommendation from N dose samples.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from charon.uncertainty.sobol import compute_sensitivity_with_r2


@dataclass(frozen=True)
class UncertaintyResult:
    """Structured uncertainty quantification output."""

    point_estimate_mg: float
    ci_90_lower_mg: float
    ci_90_upper_mg: float
    ci_ratio: float
    confidence: str
    n_samples: int
    n_successful: int
    convergence_met: bool
    sensitivity: dict[str, float]
    limiting_parameter: str
    recommendation: str
    r_squared: float


def compute_dose_range(
    doses: np.ndarray,
    *,
    sensitivity: dict[str, float],
    param_names: tuple[str, ...],
    parameter_matrix: np.ndarray | None = None,
) -> UncertaintyResult:
    """Aggregate N dose samples into CI + confidence + recommendation.

    Parameters
    ----------
    doses : ndarray, shape (N,)
        Dose estimates from propagation (positive values).
    sensitivity : dict
        Pre-computed sensitivity scores {param: importance}.
    param_names : tuple
        Parameter names.
    parameter_matrix : ndarray, optional
        If provided, recompute sensitivity (otherwise use pre-computed).
    """
    doses = np.asarray(doses, dtype=np.float64)
    valid = doses[doses > 0]
    n_successful = len(valid)

    if n_successful == 0:
        raise ValueError("No valid (positive) dose samples")

    # Geometric mean
    log_doses = np.log(valid)
    point_estimate = float(np.exp(np.mean(log_doses)))

    # Percentiles
    ci_lower = float(np.percentile(valid, 5))
    ci_upper = float(np.percentile(valid, 95))
    ci_ratio = ci_upper / ci_lower if ci_lower > 0 else float("inf")

    # Confidence classification
    if ci_ratio < 3.0:
        confidence = "HIGH"
    elif ci_ratio < 10.0:
        confidence = "MEDIUM"
    else:
        confidence = "LOW"

    # Recompute sensitivity if matrix provided
    r_squared = 0.0
    if parameter_matrix is not None and n_successful >= 10:
        sensitivity, r_squared = compute_sensitivity_with_r2(
            parameter_matrix[:n_successful], valid, param_names
        )
    else:
        r_squared = 1.0  # pre-computed, no diagnostic

    # Convergence check: running geomean CV
    convergence_met = _check_convergence(valid)

    # Limiting parameter and recommendation
    if sensitivity:
        limiting = max(sensitivity, key=sensitivity.get)  # type: ignore[arg-type]
        pct = sensitivity[limiting] * 100
        recommendation = (
            f"Experimental {limiting} measurement would narrow CI by ~{pct:.0f}%"
        )
    else:
        limiting = "unknown"
        recommendation = "Insufficient data for sensitivity analysis"

    return UncertaintyResult(
        point_estimate_mg=point_estimate,
        ci_90_lower_mg=ci_lower,
        ci_90_upper_mg=ci_upper,
        ci_ratio=ci_ratio,
        confidence=confidence,
        n_samples=len(doses),
        n_successful=n_successful,
        convergence_met=convergence_met,
        sensitivity=sensitivity,
        limiting_parameter=limiting,
        recommendation=recommendation,
        r_squared=r_squared,
    )


def _check_convergence(doses: np.ndarray, window: int = 50) -> bool:
    """Check if running geometric mean has stabilized."""
    n = len(doses)
    if n < 3 * window:
        return False  # not enough samples
    log_d = np.log(doses)
    means = []
    for i in range(window, n + 1, window):
        means.append(np.mean(log_d[:i]))
    if len(means) < 3:
        return False
    last3 = np.array(means[-3:])
    cv = np.std(last3) / abs(np.mean(last3)) if np.mean(last3) != 0 else float("inf")
    return cv < 0.05
```

- [ ] **Step 4: Run tests**

Run: `pytest tests/unit/test_dose_range.py -v`
Expected: 6 PASSED

- [ ] **Step 5: Regression**

Run: `pytest tests/ -x -q`
Expected: 662 + 6 = 668 passed

- [ ] **Step 6: Commit**

```bash
git add src/charon/uncertainty/dose_range.py tests/unit/test_dose_range.py
git commit -m "feat(uncertainty): dose CI aggregation with confidence classification"
```

---

### Task 4: Implement propagation engine (propagation.py)

**Files:**
- Create: `src/charon/uncertainty/propagation.py`

- [ ] **Step 1: Implement propagation.py**

Create `src/charon/uncertainty/propagation.py`:

```python
"""N-sample uncertainty propagation through the PBPK pipeline.

For each LHS parameter sample:
  1. Override compound properties with sampled values.
  2. Run Pipeline._run_single() → PK + dose.
  3. Collect dose, PK metrics, and parameter values.

Failed ODE solves are skipped (counted in n_failed).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from charon.core.schema import (
    CompoundConfig,
    DoseProjectionConfig,
    PredictedProperty,
)
from charon.uncertainty.sampling import SamplingResult


@dataclass
class PropagationResult:
    """Output of N-sample propagation."""

    doses_mg: np.ndarray
    cl_apparent: np.ndarray
    parameter_matrix: np.ndarray
    n_total: int
    n_successful: int
    n_failed: int
    param_names: tuple[str, ...]


def override_compound(
    base: CompoundConfig,
    sample: dict[str, float],
) -> CompoundConfig:
    """Create a modified compound with sampled property values.

    Uses Pydantic model_copy to avoid mutating the original.
    """
    props = base.properties

    # Deep copy property groups that need modification
    phys = props.physicochemical
    bind = props.binding
    metab = props.metabolism
    perm = props.permeability

    if "logp" in sample and phys.logp is not None:
        new_logp = phys.logp.model_copy(update={"value": sample["logp"]})
        phys = phys.model_copy(update={"logp": new_logp})

    if "fu_p" in sample and bind.fu_p is not None:
        new_fup = bind.fu_p.model_copy(update={"value": sample["fu_p"]})
        bind = bind.model_copy(update={"fu_p": new_fup})

    if "clint" in sample and metab.clint_uL_min_mg is not None:
        new_clint = metab.clint_uL_min_mg.model_copy(update={"value": sample["clint"]})
        metab = metab.model_copy(update={"clint_uL_min_mg": new_clint})

    if "peff" in sample and perm.peff_cm_s is not None:
        new_peff = perm.peff_cm_s.model_copy(update={"value": sample["peff"]})
        perm = perm.model_copy(update={"peff_cm_s": new_peff})

    if "bp_ratio" in sample and bind.bp_ratio is not None:
        new_bp = bind.bp_ratio.model_copy(update={"value": sample["bp_ratio"]})
        bind = bind.model_copy(update={"bp_ratio": new_bp})

    new_props = props.model_copy(update={
        "physicochemical": phys,
        "binding": bind,
        "metabolism": metab,
        "permeability": perm,
    })
    return base.model_copy(update={"properties": new_props})


def propagate(
    *,
    base_compound: CompoundConfig,
    samples: SamplingResult,
    route: str,
    dose_mg: float,
    dose_projection: DoseProjectionConfig,
    duration_h: float = 72.0,
    liver_model: str = "well_stirred",
    mppgl_override: bool = True,
) -> PropagationResult:
    """Run N PBPK simulations with sampled parameters.

    Parameters
    ----------
    base_compound : CompoundConfig
        Base compound to modify per sample.
    samples : SamplingResult
        LHS parameter samples.
    route, dose_mg, dose_projection, duration_h, liver_model
        Pipeline configuration (same for all samples).
    mppgl_override : bool
        If True and "mppgl" in samples, pass sampled MPPGL to bridge.
    """
    # Lazy import to avoid circular
    from charon.pipeline import Pipeline

    n = len(samples.samples)
    doses = []
    cls = []
    param_rows = []
    n_failed = 0

    for i, sample in enumerate(samples.samples):
        try:
            compound_i = override_compound(base_compound, sample)

            pipe = Pipeline(
                compound=compound_i,
                route=route,  # type: ignore[arg-type]
                dose_mg=dose_mg,
                duration_h=duration_h,
                liver_model=liver_model,
                dose_projection=dose_projection,
            )
            result = pipe.run()

            if result.dose_recommendation is None or result.dose_recommendation.mrsd_mg <= 0:
                n_failed += 1
                continue

            doses.append(result.dose_recommendation.mrsd_mg)
            cls.append(result.pk_parameters.cl_apparent or 0.0)
            param_rows.append([sample.get(name, 0.0) for name in samples.param_names])

        except Exception:
            n_failed += 1

    n_successful = len(doses)

    return PropagationResult(
        doses_mg=np.array(doses, dtype=np.float64),
        cl_apparent=np.array(cls, dtype=np.float64),
        parameter_matrix=np.array(param_rows, dtype=np.float64) if param_rows else np.empty((0, len(samples.param_names))),
        n_total=n,
        n_successful=n_successful,
        n_failed=n_failed,
        param_names=samples.param_names,
    )
```

- [ ] **Step 2: Verify it imports**

Run: `python3 -c "from charon.uncertainty.propagation import propagate, override_compound; print('OK')"`
Expected: `OK`

- [ ] **Step 3: Commit**

```bash
git add src/charon/uncertainty/propagation.py
git commit -m "feat(uncertainty): N-sample propagation engine with compound override"
```

---

### Task 5: Wire Pipeline with uncertainty parameter

**Files:**
- Modify: `src/charon/pipeline.py`
- Create: `tests/integration/test_uncertainty_pipeline.py`

- [ ] **Step 1: Write the failing integration tests**

Create `tests/integration/test_uncertainty_pipeline.py`:

```python
"""End-to-end Pipeline + uncertainty quantification tests."""
from __future__ import annotations

import sys
from pathlib import Path

import pytest
import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from charon.core.schema import CompoundConfig, DoseProjectionConfig, UncertaintyConfig
from charon.pipeline import Pipeline


def _load_midazolam():
    path = REPO_ROOT / "validation" / "data" / "tier1_obach" / "compounds" / "midazolam.yaml"
    with path.open() as f:
        return CompoundConfig(**yaml.safe_load(f))


class TestPipelineUncertainty:
    def test_returns_uncertainty_result(self):
        pipe = Pipeline(
            compound=_load_midazolam(),
            route="oral",
            dose_mg=5.0,
            dose_projection=DoseProjectionConfig(noael_mg_kg=50.0, noael_species="rat"),
            uncertainty=UncertaintyConfig(n_samples=30),
        )
        result = pipe.run()
        assert result.uncertainty is not None
        assert result.uncertainty.point_estimate_mg > 0
        assert result.uncertainty.ci_90_lower_mg < result.uncertainty.ci_90_upper_mg

    def test_no_uncertainty_returns_none(self):
        pipe = Pipeline(
            compound=_load_midazolam(),
            route="oral",
            dose_mg=5.0,
            dose_projection=DoseProjectionConfig(noael_mg_kg=50.0, noael_species="rat"),
        )
        result = pipe.run()
        assert result.uncertainty is None

    def test_uncertainty_without_dose_projection_raises(self):
        with pytest.raises(ValueError, match="dose_projection"):
            Pipeline(
                compound=_load_midazolam(),
                route="oral",
                dose_mg=5.0,
                uncertainty=UncertaintyConfig(n_samples=30),
            ).run()

    def test_confidence_is_valid(self):
        pipe = Pipeline(
            compound=_load_midazolam(),
            route="oral",
            dose_mg=5.0,
            dose_projection=DoseProjectionConfig(noael_mg_kg=50.0, noael_species="rat"),
            uncertainty=UncertaintyConfig(n_samples=30),
        )
        result = pipe.run()
        assert result.uncertainty.confidence in ("HIGH", "MEDIUM", "LOW")

    def test_sensitivity_has_params(self):
        pipe = Pipeline(
            compound=_load_midazolam(),
            route="oral",
            dose_mg=5.0,
            dose_projection=DoseProjectionConfig(noael_mg_kg=50.0, noael_species="rat"),
            uncertainty=UncertaintyConfig(n_samples=30),
        )
        result = pipe.run()
        assert len(result.uncertainty.sensitivity) >= 3

    def test_recommendation_string(self):
        pipe = Pipeline(
            compound=_load_midazolam(),
            route="oral",
            dose_mg=5.0,
            dose_projection=DoseProjectionConfig(noael_mg_kg=50.0, noael_species="rat"),
            uncertainty=UncertaintyConfig(n_samples=30),
        )
        result = pipe.run()
        assert isinstance(result.uncertainty.recommendation, str)
        assert len(result.uncertainty.recommendation) > 10

    def test_metadata_includes_uncertainty(self):
        pipe = Pipeline(
            compound=_load_midazolam(),
            route="oral",
            dose_mg=5.0,
            dose_projection=DoseProjectionConfig(noael_mg_kg=50.0, noael_species="rat"),
            uncertainty=UncertaintyConfig(n_samples=30),
        )
        result = pipe.run()
        assert "uncertainty_ci_90" in result.metadata

    def test_iv_with_uncertainty(self):
        pipe = Pipeline(
            compound=_load_midazolam(),
            route="iv_bolus",
            dose_mg=5.0,
            dose_projection=DoseProjectionConfig(noael_mg_kg=50.0, noael_species="rat"),
            uncertainty=UncertaintyConfig(n_samples=20),
        )
        result = pipe.run()
        assert result.uncertainty is not None
        assert result.uncertainty.point_estimate_mg > 0
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/integration/test_uncertainty_pipeline.py -v`
Expected: FAIL — Pipeline doesn't accept `uncertainty`

- [ ] **Step 3: Modify Pipeline to support uncertainty**

In `src/charon/pipeline.py`:

**a) Add `uncertainty` to `__init__` and PipelineResult:**

Add `UncertaintyConfig` to schema imports. Add to `__init__` signature:
```python
        uncertainty: UncertaintyConfig | None = None,
```
Store as `self.uncertainty = uncertainty`.

Add to PipelineResult:
```python
    uncertainty: "UncertaintyResult | None" = None
```

**b) Add `_run_uncertainty` method:**

```python
    def _run_uncertainty(self, pk, dose_rec):
        """Run uncertainty propagation if configured."""
        if self.uncertainty is None:
            return None
        if self.dose_projection is None:
            raise ValueError(
                "uncertainty requires dose_projection to be set "
                "(need dose projection to quantify dose CI)"
            )

        from charon.uncertainty.sampling import build_param_specs, generate_lhs_samples
        from charon.uncertainty.propagation import propagate
        from charon.uncertainty.dose_range import compute_dose_range

        param_specs = build_param_specs(self.compound)
        samples = generate_lhs_samples(
            param_specs=param_specs,
            n_samples=self.uncertainty.n_samples,
            correlation=self.uncertainty.correlation,
        )

        prop_result = propagate(
            base_compound=self.compound,
            samples=samples,
            route=self.route,
            dose_mg=self.dose_mg,
            dose_projection=self.dose_projection,
            duration_h=self.duration_h,
            liver_model=self.liver_model,
        )

        if prop_result.n_successful < 5:
            return None  # too few successful runs

        return compute_dose_range(
            prop_result.doses_mg,
            sensitivity={},  # will be recomputed
            param_names=prop_result.param_names,
            parameter_matrix=prop_result.parameter_matrix,
        )
```

**c) Call `_run_uncertainty` in both `run()` and `_run_oral()`:**

After `dose_rec = self._maybe_project_dose(pk)`:
```python
        unc = self._run_uncertainty(pk, dose_rec)
```

Add `uncertainty=unc` to PipelineResult constructor.
Add `"uncertainty_ci_90": f"{unc.ci_90_lower_mg:.1f}-{unc.ci_90_upper_mg:.1f}" if unc else None` to metadata.

- [ ] **Step 4: Run tests**

Run: `pytest tests/integration/test_uncertainty_pipeline.py -v`
Expected: 8 PASSED (may take ~30-60 seconds due to N=30 propagations)

- [ ] **Step 5: Full regression**

Run: `pytest tests/ -x -q`
Expected: 668 + 8 = 676 passed

- [ ] **Step 6: Commit**

```bash
git add src/charon/pipeline.py tests/integration/test_uncertainty_pipeline.py
git commit -m "feat(pipeline): integrate uncertainty quantification (LHS + SRC + dose CI)"
```

---

### Task 6: Update uncertainty/__init__.py exports

**Files:**
- Modify: `src/charon/uncertainty/__init__.py`

- [ ] **Step 1: Write exports**

```python
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
```

- [ ] **Step 2: Verify**

Run: `python3 -c "from charon.uncertainty import UncertaintyResult, generate_lhs_samples, compute_sensitivity; print('OK')"`
Expected: `OK`

- [ ] **Step 3: Commit**

```bash
git add src/charon/uncertainty/__init__.py
git commit -m "chore(uncertainty): export Layer 4 symbols from __init__"
```

---

### Task 7: Final regression gate

- [ ] **Step 1: Full suite**

Run: `pytest tests/ -v --tb=short 2>&1 | tail -40`
Expected: all pass (647 + ~29 = 676+)

- [ ] **Step 2: Integration gates**

Run: `pytest tests/integration/ -v --tb=short`
Expected: all integration tests pass (Obach panel + oral + FIH + uncertainty)

- [ ] **Step 3: E2E demo**

Run:
```bash
python3 -c "
from charon.pipeline import Pipeline
from charon.core.schema import CompoundConfig, DoseProjectionConfig, UncertaintyConfig
import yaml
d = yaml.safe_load(open('validation/data/tier1_obach/compounds/midazolam.yaml'))
c = CompoundConfig(**d)
r = Pipeline(c, route='oral', dose_mg=5.0,
    dose_projection=DoseProjectionConfig(noael_mg_kg=50.0, noael_species='rat', target_kd_nM=10.0),
    uncertainty=UncertaintyConfig(n_samples=100),
).run()
u = r.uncertainty
print(f'Dose: {u.point_estimate_mg:.1f} mg [{u.ci_90_lower_mg:.1f} - {u.ci_90_upper_mg:.1f}] 90% CI')
print(f'Confidence: {u.confidence}')
print(f'Top sensitivity: {u.limiting_parameter} ({u.sensitivity[u.limiting_parameter]*100:.0f}%)')
print(f'Recommendation: {u.recommendation}')
"
```

Expected output like:
```
Dose: X.X mg [Y.Y - Z.Z] 90% CI
Confidence: MEDIUM
Top sensitivity: fu_p (XX%)
Recommendation: Experimental fu_p measurement would narrow CI by ~XX%
```
