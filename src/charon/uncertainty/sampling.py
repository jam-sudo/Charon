"""Latin Hypercube Sampling with optional Iman-Conover rank correlation.

Generates N parameter sets for uncertainty propagation through PBPK
simulation.  Each parameter can follow a normal or lognormal marginal
distribution, with spread derived from conformal prediction CIs (when
available) or pharmacological fallback CVs.

Units
-----
All parameters are sampled in their native units.  Lognormal parameters
use log10-space internally (``log_mu = log10(mu)``, ``sigma`` is in
log10 units).  Physical bounds are enforced post-sampling:

  - fu_p:              [0.001, 1.0]
  - clint_uL_min_mg:   [0.01, 5000]
  - peff_cm_s:         [1e-7, 0.1]
  - bp_ratio:          [0.3, 5.0]
  - mppgl:             [10, 100]
  - logp:              [-5, 10]

References
----------
Iman, R.L. & Conover, W.J. (1982). A distribution-free approach to inducing
rank correlation among input variables. Communications in Statistics -
Simulation and Computation, 11(3), 311-334.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Literal

import numpy as np
from numpy.typing import NDArray
from scipy.stats import norm
from scipy.stats.qmc import LatinHypercube

from charon.core.schema import CompoundConfig

# ---------------------------------------------------------------------------
# Fallback CVs for experimental values without CIs
# ---------------------------------------------------------------------------

_FALLBACK_CV: dict[str, float] = {
    "logp": 0.15,
    "fu_p": 0.50,
    "clint_uL_min_mg": 0.60,
    "peff_cm_s": 0.40,
    "mppgl": 0.20,
    "bp_ratio": 0.10,
}

# Distribution type per parameter (normal vs lognormal)
_DIST_TYPE: dict[str, Literal["normal", "lognormal"]] = {
    "logp": "normal",
    "fu_p": "lognormal",
    "clint_uL_min_mg": "lognormal",
    "peff_cm_s": "lognormal",
    "mppgl": "normal",
    "bp_ratio": "normal",
}

# Physical bounds per parameter: (lower, upper)
_PHYSICAL_BOUNDS: dict[str, tuple[float, float]] = {
    "logp": (-5.0, 10.0),
    "fu_p": (0.001, 1.0),
    "clint_uL_min_mg": (0.01, 5000.0),
    "peff_cm_s": (1e-7, 0.1),
    "mppgl": (10.0, 100.0),
    "bp_ratio": (0.3, 5.0),
}

# Default target rank-correlation matrix (ordered as listed)
_CORR_PARAMS_ORDER = ("logp", "fu_p", "clint_uL_min_mg", "peff_cm_s", "mppgl", "bp_ratio")

_DEFAULT_TARGET_CORR = np.array([
    # logP   fu_p  CLint  Peff  MPPGL  BP
    [1.00, -0.60,  0.00,  0.40,  0.00,  0.00],  # logP
    [-0.60,  1.00,  0.00,  0.00,  0.00,  0.00],  # fu_p
    [0.00,  0.00,  1.00,  0.00,  0.00,  0.00],  # CLint
    [0.40,  0.00,  0.00,  1.00,  0.00,  0.00],  # Peff
    [0.00,  0.00,  0.00,  0.00,  1.00,  0.00],  # MPPGL
    [0.00,  0.00,  0.00,  0.00,  0.00,  1.00],  # BP
])

# CI width divisor: 90% CI = ±1.645σ → width = 3.29σ
_CI_WIDTH_DIVISOR = 3.29


# ---------------------------------------------------------------------------
# SamplingResult
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class SamplingResult:
    """Immutable container for LHS sampling output.

    Attributes
    ----------
    samples : tuple[dict[str, float], ...]
        Each element is one parameter set (name -> value).
    param_names : tuple[str, ...]
        Ordered parameter names.
    n_params_sampled : int
        Number of distinct parameters sampled.
    correlation_applied : bool
        Whether Iman-Conover correlation was applied.
    seed : int
        Random seed used for reproducibility.
    """

    samples: tuple[dict[str, float], ...]
    param_names: tuple[str, ...]
    n_params_sampled: int
    correlation_applied: bool
    seed: int


# ---------------------------------------------------------------------------
# build_param_specs
# ---------------------------------------------------------------------------

def build_param_specs(
    compound: CompoundConfig,
) -> dict[str, tuple[float, float, str]]:
    """Build parameter specifications from compound properties.

    For each ADMET parameter, determines (mu, sigma, distribution_type):

    - If a 90% CI is available (from conformal prediction), sigma is derived
      from CI width:
        - Normal:    sigma = (upper - lower) / 3.29
        - Lognormal: sigma_log = (log10(upper) - log10(lower)) / 3.29
    - If no CI (experimental values), fallback CVs are used:
        - Normal:    sigma = mu * CV
        - Lognormal: sigma_log = CV * log10(e)  [CV in original space → log10 space]

    MPPGL is always included: Normal(40, 8).
    BP ratio is included when available.

    Parameters
    ----------
    compound : CompoundConfig
        Compound with predicted or experimental properties.

    Returns
    -------
    dict[str, tuple[float, float, str]]
        Keys are parameter names.
        Values are (mu, sigma, distribution_type) where:
        - For "normal": mu and sigma are in native units
        - For "lognormal": mu is in native units, sigma is in log10 space
    """
    specs: dict[str, tuple[float, float, str]] = {}
    props = compound.properties

    # --- logP ---
    logp_prop = props.physicochemical.logp
    if logp_prop is not None:
        mu = logp_prop.value
        if logp_prop.ci_90_lower is not None and logp_prop.ci_90_upper is not None:
            sigma = (logp_prop.ci_90_upper - logp_prop.ci_90_lower) / _CI_WIDTH_DIVISOR
        else:
            sigma = abs(mu) * _FALLBACK_CV["logp"] if mu != 0 else 0.15
        specs["logp"] = (mu, sigma, "normal")

    # --- fu_p ---
    fup_prop = props.binding.fu_p
    if fup_prop is not None:
        mu = fup_prop.value
        if fup_prop.ci_90_lower is not None and fup_prop.ci_90_upper is not None:
            sigma_log = (
                math.log10(fup_prop.ci_90_upper) - math.log10(fup_prop.ci_90_lower)
            ) / _CI_WIDTH_DIVISOR
        else:
            # CV in original space -> sigma in log10 space
            sigma_log = _FALLBACK_CV["fu_p"] * math.log10(math.e)
        specs["fu_p"] = (mu, sigma_log, "lognormal")

    # --- CLint ---
    clint_prop = props.metabolism.clint_uL_min_mg
    if clint_prop is not None:
        mu = clint_prop.value
        if clint_prop.ci_90_lower is not None and clint_prop.ci_90_upper is not None:
            sigma_log = (
                math.log10(clint_prop.ci_90_upper) - math.log10(clint_prop.ci_90_lower)
            ) / _CI_WIDTH_DIVISOR
        else:
            sigma_log = _FALLBACK_CV["clint_uL_min_mg"] * math.log10(math.e)
        specs["clint_uL_min_mg"] = (mu, sigma_log, "lognormal")

    # --- Peff ---
    peff_prop = props.permeability.peff_cm_s
    if peff_prop is not None:
        mu = peff_prop.value
        if peff_prop.ci_90_lower is not None and peff_prop.ci_90_upper is not None:
            sigma_log = (
                math.log10(peff_prop.ci_90_upper) - math.log10(peff_prop.ci_90_lower)
            ) / _CI_WIDTH_DIVISOR
        else:
            sigma_log = _FALLBACK_CV["peff_cm_s"] * math.log10(math.e)
        specs["peff_cm_s"] = (mu, sigma_log, "lognormal")

    # --- BP ratio ---
    bp_prop = props.binding.bp_ratio
    if bp_prop is not None:
        mu = bp_prop.value
        if bp_prop.ci_90_lower is not None and bp_prop.ci_90_upper is not None:
            sigma = (bp_prop.ci_90_upper - bp_prop.ci_90_lower) / _CI_WIDTH_DIVISOR
        else:
            sigma = mu * _FALLBACK_CV["bp_ratio"]
        specs["bp_ratio"] = (mu, sigma, "normal")

    # --- MPPGL (always included, physiological constant) ---
    specs["mppgl"] = (40.0, 8.0, "normal")

    return specs


# ---------------------------------------------------------------------------
# Iman-Conover rank correlation
# ---------------------------------------------------------------------------

def _iman_conover(
    samples: NDArray[np.float64],
    target_corr: NDArray[np.float64],
    rng: np.random.Generator,
) -> NDArray[np.float64]:
    """Apply Iman-Conover rank reordering to induce target rank correlation.

    Algorithm (Iman & Conover 1982):
      1. Compute rank matrix R from the input samples, centre it.
      2. Compute current rank correlation C = corr(R*).
      3. Cholesky decompose C = P P^T and target T = Q Q^T.
      4. Transform: R_new = R* @ P^{-T} @ Q^T
         This maps the current correlation structure to the target.
      5. Re-sort each column of original samples by the new rank order.

    Parameters
    ----------
    samples : ndarray of shape (n_samples, n_params)
        Samples from marginal distributions (already transformed from [0,1]).
    target_corr : ndarray of shape (n_params, n_params)
        Desired Spearman rank correlation matrix (symmetric, positive-definite).
    rng : numpy.random.Generator
        Not used directly but kept for API consistency.

    Returns
    -------
    ndarray of shape (n_samples, n_params)
        Re-ordered samples with induced rank correlation structure.
    """
    n_samples, n_params = samples.shape

    # Step 1: Compute rank matrix (ties broken by first occurrence), centre it
    rank_matrix = np.empty_like(samples)
    for j in range(n_params):
        order = np.argsort(samples[:, j])
        rank_matrix[order, j] = np.arange(1, n_samples + 1, dtype=np.float64)
    rank_centered = rank_matrix - rank_matrix.mean(axis=0)

    # Step 2: Current rank correlation C
    rank_cov = (rank_centered.T @ rank_centered) / (n_samples - 1)
    rank_std = np.sqrt(np.diag(rank_cov))
    rank_corr = rank_cov / np.outer(rank_std, rank_std)
    rank_corr = (rank_corr + rank_corr.T) / 2.0
    np.fill_diagonal(rank_corr, 1.0)

    # Step 3: Cholesky decompositions — C = P P^T, T_target = Q Q^T
    P = np.linalg.cholesky(rank_corr)       # lower-triangular
    Q = np.linalg.cholesky(target_corr)      # lower-triangular

    # Step 4: Transform R_new = R* @ P^{-T} @ Q^T
    # P^{-T} = inv(P^T) = inv(P)^T
    P_inv_T = np.linalg.inv(P).T
    new_ranks = rank_centered @ P_inv_T @ Q.T

    # Step 5: Re-sort each column of original samples by new rank order
    result = np.empty_like(samples)
    for j in range(n_params):
        sorted_col = np.sort(samples[:, j])
        new_order = np.argsort(np.argsort(new_ranks[:, j]))
        result[:, j] = sorted_col[new_order]

    return result


# ---------------------------------------------------------------------------
# Main sampling function
# ---------------------------------------------------------------------------

def generate_lhs_samples(
    *,
    param_specs: dict[str, tuple[float, float, str]],
    n_samples: int = 500,
    correlation: Literal["iman_conover", "none"] = "iman_conover",
    target_correlation: NDArray[np.float64] | None = None,
    seed: int = 42,
) -> SamplingResult:
    """Generate Latin Hypercube samples with optional Iman-Conover correlation.

    Parameters
    ----------
    param_specs : dict[str, tuple[float, float, str]]
        Keys are parameter names. Values are (mu, sigma, distribution_type)
        where distribution_type is "normal" or "lognormal".
        For lognormal: mu is in native units, sigma is in log10 space.
    n_samples : int
        Number of parameter sets to generate.
    correlation : {"iman_conover", "none"}
        Whether to apply Iman-Conover rank reordering.
    target_correlation : ndarray or None
        Custom target correlation matrix. If None, uses the default
        6x6 matrix from ARCHITECTURE.md (logP-fu_p: -0.6, logP-Peff: 0.4).
    seed : int
        Random seed for reproducibility.

    Returns
    -------
    SamplingResult
        Frozen dataclass with samples and metadata.
    """
    param_names = tuple(param_specs.keys())
    n_params = len(param_names)

    if n_params == 0:
        return SamplingResult(
            samples=(),
            param_names=(),
            n_params_sampled=0,
            correlation_applied=False,
            seed=seed,
        )

    # --- Step 1: Generate LHS samples in [0, 1] ---
    rng = np.random.default_rng(seed)
    sampler = LatinHypercube(d=n_params, seed=rng)
    unit_samples = sampler.random(n=n_samples)  # shape (n_samples, n_params)

    # --- Step 2: Transform to target marginal distributions ---
    samples = np.empty_like(unit_samples)
    for j, name in enumerate(param_names):
        mu, sigma, dist = param_specs[name]
        if dist == "normal":
            samples[:, j] = norm.ppf(unit_samples[:, j], loc=mu, scale=sigma)
        elif dist == "lognormal":
            # mu is native-scale value; sigma is log10-space std dev
            log_mu = math.log10(mu)
            log_vals = norm.ppf(unit_samples[:, j], loc=log_mu, scale=sigma)
            samples[:, j] = np.power(10.0, log_vals)
        else:
            raise ValueError(f"Unknown distribution type {dist!r} for {name}")

    # --- Step 3: Apply Iman-Conover correlation (optional) ---
    correlation_applied = False
    if correlation == "iman_conover" and n_params >= 2:
        # Build target correlation sub-matrix for the active parameters
        if target_correlation is not None:
            corr_matrix = target_correlation
        else:
            corr_matrix = _build_submatrix(param_names)

        if corr_matrix is not None:
            samples = _iman_conover(samples, corr_matrix, rng)
            correlation_applied = True

    # --- Step 4: Clip to physical bounds ---
    for j, name in enumerate(param_names):
        if name in _PHYSICAL_BOUNDS:
            lo, hi = _PHYSICAL_BOUNDS[name]
            samples[:, j] = np.clip(samples[:, j], lo, hi)

    # --- Step 5: Pack into list of dicts ---
    sample_dicts = tuple(
        {name: float(samples[i, j]) for j, name in enumerate(param_names)}
        for i in range(n_samples)
    )

    return SamplingResult(
        samples=sample_dicts,
        param_names=param_names,
        n_params_sampled=n_params,
        correlation_applied=correlation_applied,
        seed=seed,
    )


def _build_submatrix(
    param_names: tuple[str, ...],
) -> NDArray[np.float64] | None:
    """Extract the sub-matrix of the default correlation for active params.

    Returns an identity matrix if no params overlap with the default set,
    which effectively means "no correlation" (Iman-Conover with I is a no-op).
    """
    n = len(param_names)
    sub = np.eye(n)

    # Map active parameter names to their indices in the default matrix
    idx_map: dict[str, int] = {}
    for name in param_names:
        if name in _CORR_PARAMS_ORDER:
            idx_map[name] = _CORR_PARAMS_ORDER.index(name)

    # Fill in off-diagonal entries from the default matrix
    for i, name_i in enumerate(param_names):
        for j, name_j in enumerate(param_names):
            if i != j and name_i in idx_map and name_j in idx_map:
                sub[i, j] = _DEFAULT_TARGET_CORR[idx_map[name_i], idx_map[name_j]]

    return sub
