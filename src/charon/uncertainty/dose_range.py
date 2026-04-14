"""Dose CI aggregation from uncertainty propagation results.

Takes an array of Monte Carlo dose samples (output of the PBPK/translational
pipeline run under LHS sampling) and computes:

  - Geometric-mean point estimate
  - 90% CI (5th–95th percentile)
  - CI ratio (upper / lower) as a spread metric
  - Confidence classification: HIGH (<3×), MEDIUM (<10×), LOW (≥10×)
  - SRC-based sensitivity indices (which parameters drive uncertainty)
  - Convergence check (rolling CV of cumulative mean < 5%)
  - Limiting parameter and prioritisation recommendation

Units
-----
All dose values are assumed to be in **mg** (mg/kg or mg total — caller
specifies context).  The point estimate and CI bounds are returned in the
same unit as the input.

References
----------
FDA Guidance for Industry: Estimating the Maximum Safe Starting Dose in
Initial Clinical Trials (2005), pp. 7–9.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from charon.uncertainty.sobol import compute_sensitivity_with_r2


@dataclass(frozen=True)
class UncertaintyResult:
    """Immutable summary of dose uncertainty analysis.

    Attributes
    ----------
    point_estimate_mg : float
        Geometric mean of the valid dose distribution (mg).
    ci_90_lower_mg : float
        5th percentile of valid doses (mg).
    ci_90_upper_mg : float
        95th percentile of valid doses (mg).
    ci_ratio : float
        ci_90_upper_mg / ci_90_lower_mg.  Dimensionless spread metric.
    confidence : str
        "HIGH" (ratio < 3), "MEDIUM" (ratio < 10), or "LOW" (ratio ≥ 10).
    n_samples : int
        Total number of samples passed (including invalid/negative).
    n_successful : int
        Number of valid (positive) samples used for statistics.
    convergence_met : bool
        True if rolling cumulative-mean CV < 5% over last 3 windows.
    sensitivity : dict[str, float]
        SRC² importance per parameter (values sum to ≈1).
    limiting_parameter : str
        Parameter with highest sensitivity index.
    recommendation : str
        Human-readable prioritisation suggestion.
    r_squared : float
        R² of the SRC regression in log-dose space.  Values <0.7 indicate
        substantial non-linearity; SRC indices should be treated as
        approximate.
    """

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
    """Aggregate Monte Carlo dose samples into an UncertaintyResult.

    Parameters
    ----------
    doses : np.ndarray, shape (n_samples,)
        Dose samples in mg.  Negative and zero values are filtered out.
    sensitivity : dict[str, float]
        Pre-computed sensitivity dict (used as fallback when
        *parameter_matrix* is None or too small to recompute).
    param_names : tuple[str, ...]
        Parameter names matching columns in *parameter_matrix*.
    parameter_matrix : np.ndarray or None, shape (n_samples, n_params)
        If provided *and* n_successful >= 10, SRC indices are recomputed
        from the valid-dose sub-matrix, overriding the *sensitivity* arg.

    Returns
    -------
    UncertaintyResult

    Raises
    ------
    ValueError
        If no valid (positive) dose samples are present.

    Notes
    -----
    Hand-calculation check (test_geometric_mean):
        doses = [10, 100]
        log_doses = [log(10), log(100)] = [2.303, 4.605]
        mean_log = 3.454
        geomean = exp(3.454) ≈ 31.62  ✓

    CI check:
        5th pct of [10, 100] = 10.0 (only 2 samples)
        95th pct = 100.0
        ci_ratio = 10.0 → "LOW"
    """
    doses = np.asarray(doses, dtype=np.float64)
    valid = doses[doses > 0]
    n_successful = int(len(valid))

    if n_successful == 0:
        raise ValueError("No valid (positive) dose samples — cannot compute CI")

    # --- Point estimate: geometric mean ---
    log_doses = np.log(valid)
    point_estimate = float(np.exp(np.mean(log_doses)))

    # --- 90% CI (5th–95th percentile) ---
    ci_lower = float(np.percentile(valid, 5))
    ci_upper = float(np.percentile(valid, 95))
    ci_ratio = ci_upper / ci_lower if ci_lower > 0 else float("inf")

    # --- Confidence classification ---
    if ci_ratio < 3.0:
        confidence = "HIGH"
    elif ci_ratio < 10.0:
        confidence = "MEDIUM"
    else:
        confidence = "LOW"

    # --- SRC sensitivity (recompute if matrix provided and enough samples) ---
    r_squared = 0.0
    if parameter_matrix is not None and n_successful >= 10:
        valid_mask = doses > 0
        valid_matrix = np.asarray(parameter_matrix, dtype=np.float64)[valid_mask]
        sensitivity, r_squared = compute_sensitivity_with_r2(
            valid_matrix, valid, param_names
        )
    else:
        # Use provided sensitivity; R² not computable without matrix
        r_squared = 1.0

    # --- Convergence check ---
    convergence_met = _check_convergence(valid)

    # --- Limiting parameter and recommendation ---
    if sensitivity:
        limiting = max(sensitivity, key=sensitivity.__getitem__)
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
        n_samples=int(len(doses)),
        n_successful=n_successful,
        convergence_met=convergence_met,
        sensitivity=sensitivity,
        limiting_parameter=limiting,
        recommendation=recommendation,
        r_squared=r_squared,
    )


def _check_convergence(doses: np.ndarray, window: int = 50) -> bool:
    """Check whether the running geometric mean has converged.

    Convergence criterion: the coefficient of variation (CV) of the last
    three cumulative-mean values falls below 5%.

    Parameters
    ----------
    doses : np.ndarray
        Valid (positive) dose samples in chronological order.
    window : int
        Window size for computing cumulative means.  Defaults to 50.

    Returns
    -------
    bool
        True if converged, False if not enough samples or CV ≥ 5%.

    Notes
    -----
    Requires at least 3 * window samples for a meaningful convergence
    assessment.
    """
    n = len(doses)
    if n < 3 * window:
        return False

    log_d = np.log(doses)
    means = []
    for i in range(window, n + 1, window):
        means.append(float(np.mean(log_d[:i])))

    if len(means) < 3:
        return False

    last3 = np.array(means[-3:])
    mean_val = np.mean(last3)
    if mean_val == 0:
        return False
    cv = np.std(last3) / abs(mean_val)
    return bool(cv < 0.05)
