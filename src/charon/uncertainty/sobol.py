"""Sensitivity analysis via Standardized Regression Coefficients (SRC).

SRC² approximates first-order Sobol indices for near-linear models.
Regression is performed in log-dose space:

    log(dose) = β₀ + β₁z₁ + β₂z₂ + ... + ε

where z_i are standardized parameters (zero mean, unit variance).
Importance is normalised so all indices sum to 1:

    importance_i = β_i² / Σβ_j²

R² of the regression is returned as a linearity diagnostic.  When
R² < 0.7 the model is substantially non-linear and the SRC indices
should be interpreted with caution.

References
----------
Saltelli, A. et al. (2008). Global Sensitivity Analysis: The Primer.
Wiley. Chapter 4.
"""

from __future__ import annotations

import numpy as np


def compute_sensitivity(
    parameter_matrix: np.ndarray,
    doses: np.ndarray,
    param_names: tuple[str, ...],
) -> dict[str, float]:
    """Compute SRC-based sensitivity indices (importance scores sum to 1).

    Parameters
    ----------
    parameter_matrix : np.ndarray, shape (n_samples, n_params)
        Matrix of sampled parameter values (one row per sample).
    doses : np.ndarray, shape (n_samples,)
        Corresponding dose outputs (positive values only used).
    param_names : tuple[str, ...]
        Parameter names in the same column order as *parameter_matrix*.

    Returns
    -------
    dict[str, float]
        Mapping of parameter name -> normalised SRC² importance [0, 1].
        Values sum to approximately 1.0.
    """
    importance, _ = compute_sensitivity_with_r2(parameter_matrix, doses, param_names)
    return importance


def compute_sensitivity_with_r2(
    parameter_matrix: np.ndarray,
    doses: np.ndarray,
    param_names: tuple[str, ...],
) -> tuple[dict[str, float], float]:
    """Compute SRC-based sensitivity indices and regression R².

    Parameters
    ----------
    parameter_matrix : np.ndarray, shape (n_samples, n_params)
        Matrix of sampled parameter values.
    doses : np.ndarray, shape (n_samples,)
        Corresponding dose outputs (positive values only used).
    param_names : tuple[str, ...]
        Parameter names in the same column order as *parameter_matrix*.

    Returns
    -------
    importance : dict[str, float]
        Normalised SRC² importance per parameter.
    r_squared : float
        Coefficient of determination of the standardised regression.
        High values (>0.7) indicate that SRC indices are reliable.

    Notes
    -----
    Regression is performed on *all* provided rows (caller is responsible
    for pre-filtering invalid samples before passing to this function).
    Log-clipping at 1e-30 prevents log(0) errors.

    Hand-calculation check (single dominant param example):
        dose = exp(2*X0 + 0.1*X1 + 0.01*X2)
        β_standardised ≈ [2σ₀, 0.1σ₁, 0.01σ₂] / y_std
        SRC²₀ ≈ (2σ₀)² / [(2σ₀)² + (0.1σ₁)² + (0.01σ₂)²] >> 0.8
    """
    X = np.asarray(parameter_matrix, dtype=np.float64)
    y = np.log(np.maximum(np.asarray(doses, dtype=np.float64), 1e-30))

    n, p = X.shape
    if p != len(param_names):
        raise ValueError(
            f"parameter_matrix has {p} columns but {len(param_names)} names given"
        )

    # --- Standardise X ---
    X_mean = X.mean(axis=0)
    X_std = X.std(axis=0)
    X_std[X_std == 0] = 1.0          # avoid div-by-zero for constant params
    Z = (X - X_mean) / X_std

    # --- Standardise y ---
    y_mean = y.mean()
    y_std = y.std()
    if y_std == 0:
        # All doses identical → uniform importance
        importance = {name: 1.0 / p for name in param_names}
        return importance, 1.0
    y_z = (y - y_mean) / y_std

    # --- OLS regression (no intercept needed after standardisation) ---
    beta, _residuals, _rank, _sv = np.linalg.lstsq(Z, y_z, rcond=None)

    # --- Normalised SRC² ---
    src_sq = beta ** 2
    total = src_sq.sum()
    if total == 0:
        importance = {name: 1.0 / p for name in param_names}
    else:
        importance = {
            name: float(src_sq[j] / total) for j, name in enumerate(param_names)
        }

    # --- R² ---
    y_pred = Z @ beta
    ss_res = np.sum((y_z - y_pred) ** 2)
    ss_tot = np.sum(y_z ** 2)
    r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0

    return importance, float(r_squared)
