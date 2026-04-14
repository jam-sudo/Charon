# Sprint 5 — Uncertainty Quantification (LHS + Sensitivity + Dose CI)

**Date:** 2026-04-13
**Scope:** LHS sampling with Iman-Conover correlation, N-sample PBPK
propagation, dose 90% confidence interval, SRC-based sensitivity
analysis, Pipeline integration.
**Prior state:** Sprint 4a complete (647 tests, SMILES → oral PK → FIH
dose). Layer 1 conformal CIs exist but are not propagated to dose.
**Expected scale:** ~15 tasks.

---

## 1. Goal

Propagate prediction uncertainty through the entire pipeline to produce:

1. Dose recommendation with 90% CI: `45 mg [28–72 mg]`.
2. Per-parameter sensitivity ranking: "fu_p accounts for 42% of dose
   variance."
3. Actionable recommendation: "Experimental fu_p measurement would
   narrow CI by ~40%."

**Non-goals:**

- True Saltelli Sobol indices (Phase B — requires N×(2p+2) evaluations)
- Population variability / virtual trial (Phase B)
- Multiprocessing / parallelization (Phase B)
- Monte Carlo method (LHS only; MC stub deferred)
- Correlation matrix estimation from training data (use literature values)

---

## 2. Architecture

```
Layer 1 Conformal CIs + Literature CVs
            ↓
    [1] LHS Sampler (N=500 parameter sets)
            ↓ Iman-Conover correlation
    [2] Propagation Engine (N × Pipeline.run())
            ↓
    [3] Dose Aggregator (geometric mean, percentiles)
            ↓
    [4] Sensitivity Analyzer (SRC² ≈ Sobol S_i)
            ↓
    UncertaintyResult (CI, sensitivity, recommendation)
```

---

## 3. Sampling (sampling.py)

### 3.1 Parameters to sample

| Parameter | Distribution | Range source | Transform |
|-----------|-------------|--------------|-----------|
| logP | Normal(μ=pred, σ) | Conformal CI: σ = (log_upper - log_lower) / 3.29 | linear |
| fu_p | LogNormal | Conformal CI (log-space symmetric) | log10 (CLAUDE.md §7.10) |
| CLint | LogNormal | Conformal CI (log-space symmetric) | log10 |
| Peff | LogNormal | Conformal CI (log-space symmetric) | log10 |
| MPPGL | Normal(40, 8) | CV ~20% literature | linear |
| BP ratio | Normal(μ=pred, σ=0.1×pred) | CV ~10% assumed | linear |

If a conformal CI is not available for a property (e.g., experimental
value with no CI), that parameter is held fixed (not sampled).

### 3.2 LHS algorithm

Standard scipy LHS:
```python
from scipy.stats.qmc import LatinHypercube
sampler = LatinHypercube(d=n_params, seed=rng)
unit_samples = sampler.random(n=n_samples)  # (N, p) in [0,1]
# Transform each column to target marginal distribution
```

### 3.3 Iman-Conover rank correlation

Induces a target rank correlation matrix on the LHS samples while
preserving each parameter's marginal distribution.

Algorithm:
1. Generate N×p uncorrelated LHS samples (already done in 3.2).
2. Compute the rank matrix R of the samples.
3. Compute the Cholesky decomposition of the target correlation
   matrix: L_target = cholesky(C_target).
4. Compute the Cholesky of the observed rank correlation matrix:
   L_observed = cholesky(corr(R)).
5. Transform ranks: R_new = R × L_observed⁻¹ᵀ × L_targetᵀ.
6. Re-sort each column of the original samples by the new rank order.

### 3.4 Target correlation matrix

From ARCHITECTURE §Layer 4 (literature values):

```
         logP   fu_p  CLint  Peff  MPPGL  BP
logP     1.00  -0.60  0.00   0.40  0.00  0.00
fu_p    -0.60   1.00  0.00   0.00  0.00  0.00
CLint    0.00   0.00  1.00   0.00  0.00  0.00
Peff     0.40   0.00  0.00   1.00  0.00  0.00
MPPGL    0.00   0.00  0.00   0.00  1.00  0.00
BP       0.00   0.00  0.00   0.00  0.00  1.00
```

### 3.5 Interface

```python
@dataclass(frozen=True)
class ParameterSample:
    """One sampled parameter set."""
    logp: float
    fu_p: float
    clint_uL_min_mg: float
    peff_cm_s: float | None
    mppgl: float
    bp_ratio: float

@dataclass(frozen=True)
class SamplingResult:
    """N parameter sets from LHS + Iman-Conover."""
    samples: tuple[ParameterSample, ...]
    n_params_sampled: int
    correlation_applied: bool
    seed: int

def generate_lhs_samples(
    *,
    compound: CompoundConfig,
    n_samples: int = 500,
    correlation: str = "iman_conover",
    seed: int = 42,
) -> SamplingResult:
```

The function reads conformal CIs from the compound's
PredictedProperty.ci_90_lower / ci_90_upper to derive per-parameter
distribution widths. Properties without CIs use literature CVs as
fallback.

---

## 4. Propagation (propagation.py)

### 4.1 Core loop

```python
def propagate_uncertainty(
    *,
    base_compound: CompoundConfig,
    samples: SamplingResult,
    route: str,
    dose_mg: float,
    dose_projection: DoseProjectionConfig,
    duration_h: float = 72.0,
    liver_model: str = "well_stirred",
) -> PropagationResult:
```

For each sample:
1. Clone base_compound with sampled property overrides.
2. Run `Pipeline(compound_i, route, dose_mg, dose_projection).run()`.
3. Collect: mrsd_mg, cl_apparent, auc, fa, fg, fh, bioavailability.
4. If ODE solver fails for a sample → skip, increment n_failed.

### 4.2 Compound override mechanism

```python
def _override_compound(
    base: CompoundConfig,
    sample: ParameterSample,
) -> CompoundConfig:
```

Uses Pydantic model_copy(update=...) to create a modified compound
with sampled logP, fu_p, CLint, Peff, BP. MPPGL override is passed
to the Pipeline via a separate mechanism (kwargs or bridge override).

### 4.3 Output

```python
@dataclass
class PropagationResult:
    doses_mg: np.ndarray           # (N_successful,)
    cl_apparent: np.ndarray
    auc_0_inf: np.ndarray
    fa_values: np.ndarray
    fg_values: np.ndarray
    fh_values: np.ndarray
    parameter_matrix: np.ndarray   # (N_successful, n_params) for SRC
    n_total: int
    n_successful: int
    n_failed: int
    param_names: tuple[str, ...]
```

---

## 5. Dose Aggregation (dose_range.py)

### 5.1 Statistics

```python
point_estimate = geometric_mean(doses)
ci_90_lower = np.percentile(doses, 5)
ci_90_upper = np.percentile(doses, 95)
ci_ratio = ci_90_upper / ci_90_lower
```

### 5.2 Confidence classification

| ci_ratio | Confidence | Meaning |
|----------|-----------|---------|
| < 3 | HIGH | Narrow CI, predictions reliable |
| 3–10 | MEDIUM | Moderate uncertainty, usable with caution |
| ≥ 10 | LOW | Wide CI, experimental data recommended |

### 5.3 Convergence check

Running geometric mean of doses computed at every 50-sample increment.
If CV of the last 3 increments < 5%, convergence met.

### 5.4 Interface

```python
@dataclass(frozen=True)
class UncertaintyResult:
    point_estimate_mg: float
    ci_90_lower_mg: float
    ci_90_upper_mg: float
    ci_ratio: float
    confidence: str                  # HIGH / MEDIUM / LOW
    n_samples: int
    n_successful: int
    convergence_met: bool
    sensitivity: dict[str, float]    # param → SRC² (normalized)
    limiting_parameter: str
    recommendation: str
    pk_summary: dict[str, dict]      # per-PK-metric: mean, P5, P95

def compute_dose_uncertainty(
    propagation: PropagationResult,
) -> UncertaintyResult:
```

---

## 6. Sensitivity Analysis (sobol.py)

### 6.1 SRC method

Standardized Regression Coefficients approximate first-order Sobol
indices for near-linear models. PBPK dose response is approximately
monotonic in each parameter, so SRC² ≈ S_i.

```python
# Standardize inputs and output
Z = (X - X.mean(0)) / X.std(0)  # (N, p) standardized params
y = (log_doses - log_doses.mean()) / log_doses.std()

# OLS regression
beta = np.linalg.lstsq(Z, y, rcond=None)[0]  # (p,)
src_squared = beta ** 2
# Normalize to sum to ~1 (R² of the regression)
importance = src_squared / src_squared.sum()
```

### 6.2 Output

```python
def compute_sensitivity(
    parameter_matrix: np.ndarray,
    doses: np.ndarray,
    param_names: tuple[str, ...],
) -> dict[str, float]:
```

Returns `{"fu_p": 0.42, "CLint": 0.31, ...}` — normalized importance
scores that sum to ~1.0 (up to R² of the linear fit).

### 6.3 Actionable recommendation

The parameter with highest SRC² gets a recommendation:
```
"Experimental {param} measurement would narrow CI by ~{importance*100:.0f}%"
```

### 6.4 Naming convention

The output field is called `sensitivity` (not `sobol_indices`) to be
honest about the method. The spec and report can state "approximate
first-order sensitivity indices via standardized regression" rather
than claiming Sobol analysis.

---

## 7. Pipeline Integration

### 7.1 Pipeline changes

Add `uncertainty: UncertaintyConfig | None = None` to Pipeline.__init__.

When `uncertainty` is not None AND `dose_projection` is not None:
1. Generate LHS samples.
2. Run propagation (N × Pipeline inner loop).
3. Compute dose CI + sensitivity.
4. Attach UncertaintyResult to PipelineResult.

When `uncertainty` is set but `dose_projection` is None → ValueError
(uncertainty needs dose projection to quantify dose CI).

### 7.2 PipelineResult extension

```python
@dataclass
class PipelineResult:
    ...
    uncertainty: UncertaintyResult | None = None
```

### 7.3 Inner loop design

The propagation loop calls a lightweight internal method that skips
the uncertainty layer (to avoid recursion). This is a private
`_run_single(compound_i)` method that runs one PBPK simulation +
dose projection without uncertainty:

```python
def _run_single(self, compound: CompoundConfig) -> PipelineResult:
    """Single deterministic run (no uncertainty recursion)."""
    # Same as run() but without uncertainty propagation
```

---

## 8. Schema Changes

### 8.1 UncertaintyConfig — already exists, no changes needed

```python
class UncertaintyConfig(BaseModel):
    method: Literal["lhs", "monte_carlo"] = "lhs"
    n_samples: int = 500
    correlation: Literal["iman_conover", "none"] = "iman_conover"
```

### 8.2 PipelineResult — add uncertainty field

Already shown in §7.2.

---

## 9. File Changes

### New files

| File | Purpose |
|------|---------|
| `src/charon/uncertainty/sampling.py` | LHS + Iman-Conover + ParameterSample |
| `src/charon/uncertainty/propagation.py` | N-sample propagation loop |
| `src/charon/uncertainty/dose_range.py` | Dose CI aggregation + confidence |
| `src/charon/uncertainty/sobol.py` | SRC-based sensitivity analysis |
| `tests/unit/test_sampling.py` | LHS + correlation tests |
| `tests/unit/test_sensitivity.py` | SRC computation tests |
| `tests/unit/test_dose_range.py` | CI aggregation tests |
| `tests/integration/test_uncertainty_pipeline.py` | E2E Pipeline + uncertainty |

### Modified files

| File | Changes |
|------|---------|
| `src/charon/pipeline.py` | Add uncertainty param, _run_single, propagation call |
| `src/charon/uncertainty/__init__.py` | Export public symbols |

---

## 10. Validation

### 10.1 Acceptance criteria

| Criterion | Target |
|-----------|--------|
| LHS produces N samples with near-uniform marginals | KS test p > 0.05 |
| Iman-Conover achieves target correlation ±0.15 | rank correlation check |
| N=50 propagation completes without crash | Smoke test |
| Dose CI: lower < point < upper | Ordering invariant |
| Confidence classification matches ci_ratio | HIGH/MEDIUM/LOW correct |
| SRC: fu_p or CLint is top sensitivity for midazolam | Known from Omega |
| N=500 propagation completes in < 300 seconds | Performance gate |
| All 647 existing tests pass | Regression |
| New tests: ≥20 | Coverage |

### 10.2 Midazolam E2E test

```python
Pipeline(
    compound=midazolam,
    route="oral",
    dose_mg=5.0,
    dose_projection=DoseProjectionConfig(noael_mg_kg=50.0, noael_species="rat"),
    uncertainty=UncertaintyConfig(n_samples=100),  # reduced for test speed
).run()
```

Verify: UncertaintyResult has finite CI, sensitivity dict has ≥3
parameters, recommendation string is non-empty.

---

## 11. Known Limitations

### 11.1 SRC vs true Sobol

SRC² approximates first-order Sobol S_i under linearity. For highly
nonlinear dose-response (e.g., extreme fu_p < 0.01 where CLh becomes
very sensitive), SRC may underestimate interaction effects. The R² of
the regression serves as a diagnostic: if R² < 0.7, warn that
nonlinearity is significant and SRC² may not sum to total variance.

### 11.2 ODE failures

Some parameter combinations (extreme fu_p + extreme CLint) may cause
BDF solver divergence. Failed samples are skipped. If >10% fail,
a warning is emitted. CI is computed from successful samples only.

### 11.3 Performance

N=500 at ~200ms/solve = ~100 seconds. Acceptable for CLI but not for
interactive dashboards. Multiprocessing deferred to Phase B.

### 11.4 Conformal CI as distribution width

Using conformal CIs (90% coverage) to derive sampling distribution
widths assumes the CI reflects a normal-ish distribution. This is
approximate — conformal gives marginal coverage, not a proper
posterior. The CI width is used as a pragmatic measure of prediction
uncertainty, not a Bayesian credible interval.

---

## 12. Definition of Done

- [ ] LHS sampling produces N correlated parameter sets.
- [ ] Iman-Conover achieves target correlations ±0.15.
- [ ] Propagation runs N simulations and collects dose/PK arrays.
- [ ] Dose CI computed: geometric mean + P5/P95.
- [ ] Confidence classified as HIGH/MEDIUM/LOW from ci_ratio.
- [ ] SRC sensitivity computed and normalized.
- [ ] Recommendation string generated from top sensitivity parameter.
- [ ] Pipeline(uncertainty=...) returns UncertaintyResult.
- [ ] Midazolam E2E produces finite CI with sensible sensitivity ranking.
- [ ] All 647 existing tests pass (regression).
- [ ] New tests: ≥20.
