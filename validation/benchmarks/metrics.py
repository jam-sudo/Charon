"""Benchmark metrics for Layer 1-3 validation suites.

Metrics
-------
fold_error       : max(pred/obs, obs/pred), always >= 1.0
aafe             : geometric mean of fold errors (absolute average fold error)
within_n_fold    : fraction of predictions with fold_error <= n
mae              : mean absolute error
rmse             : root mean squared error
pearson_r        : Pearson product-moment correlation coefficient
within_abs_diff  : fraction of predictions with |pred-obs| <= threshold
"""

from __future__ import annotations

import math
from collections.abc import Iterable


def fold_error(predicted: float, observed: float) -> float:
    """Two-sided fold error. Always >= 1.0."""
    if predicted <= 0 or observed <= 0:
        raise ValueError(
            f"fold_error requires positive values, got "
            f"predicted={predicted}, observed={observed}"
        )
    return max(predicted / observed, observed / predicted)


def aafe(predicted: Iterable[float], observed: Iterable[float]) -> float:
    """Absolute Average Fold Error (geometric mean of fold errors).

    Definition: AAFE = 10^( mean(|log10(pred/obs)|) ).
    """
    preds = list(predicted)
    obs = list(observed)
    if len(preds) != len(obs):
        raise ValueError(
            f"aafe: length mismatch ({len(preds)} vs {len(obs)})"
        )
    if not preds:
        raise ValueError("aafe: empty input")
    log_fold_errors = []
    for p, o in zip(preds, obs):
        if p <= 0 or o <= 0:
            raise ValueError(
                f"aafe requires positive values, got pred={p}, obs={o}"
            )
        log_fold_errors.append(abs(math.log10(p / o)))
    mean_log_fold = sum(log_fold_errors) / len(log_fold_errors)
    return 10 ** mean_log_fold


def within_n_fold(
    predicted: Iterable[float],
    observed: Iterable[float],
    *,
    n: float = 2.0,
) -> float:
    """Fraction of predictions with fold_error <= n (range [0, 1])."""
    preds = list(predicted)
    obs = list(observed)
    if len(preds) != len(obs):
        raise ValueError(
            f"within_n_fold: length mismatch ({len(preds)} vs {len(obs)})"
        )
    if not preds:
        raise ValueError("within_n_fold: empty input")
    hits = 0
    for p, o in zip(preds, obs):
        if fold_error(p, o) <= n:
            hits += 1
    return hits / len(preds)


def mae(predicted: Iterable[float], observed: Iterable[float]) -> float:
    """Mean Absolute Error: mean(|pred - obs|).

    Parameters
    ----------
    predicted : iterable of float
    observed  : iterable of float

    Returns
    -------
    float
        Non-negative MAE value.
    """
    preds = list(predicted)
    obs = list(observed)
    if len(preds) != len(obs):
        raise ValueError(
            f"mae: length mismatch ({len(preds)} vs {len(obs)})"
        )
    if not preds:
        raise ValueError("mae: empty input")
    return sum(abs(p - o) for p, o in zip(preds, obs)) / len(preds)


def rmse(predicted: Iterable[float], observed: Iterable[float]) -> float:
    """Root Mean Squared Error: sqrt(mean((pred - obs)^2)).

    Parameters
    ----------
    predicted : iterable of float
    observed  : iterable of float

    Returns
    -------
    float
        Non-negative RMSE value.
    """
    preds = list(predicted)
    obs = list(observed)
    if len(preds) != len(obs):
        raise ValueError(
            f"rmse: length mismatch ({len(preds)} vs {len(obs)})"
        )
    if not preds:
        raise ValueError("rmse: empty input")
    return math.sqrt(sum((p - o) ** 2 for p, o in zip(preds, obs)) / len(preds))


def pearson_r(predicted: Iterable[float], observed: Iterable[float]) -> float:
    """Pearson product-moment correlation coefficient.

    Returns 0.0 when n < 2 or when either series has zero standard deviation.

    Parameters
    ----------
    predicted : iterable of float
    observed  : iterable of float

    Returns
    -------
    float
        Correlation coefficient in [-1, 1], or 0.0 for degenerate cases.
    """
    preds = list(predicted)
    obs = list(observed)
    if len(preds) != len(obs):
        raise ValueError(
            f"pearson_r: length mismatch ({len(preds)} vs {len(obs)})"
        )
    if not preds:
        raise ValueError("pearson_r: empty input")
    n = len(preds)
    if n < 2:
        return 0.0
    mean_p = sum(preds) / n
    mean_o = sum(obs) / n
    deviations_p = [p - mean_p for p in preds]
    deviations_o = [o - mean_o for o in obs]
    ss_p = sum(d ** 2 for d in deviations_p)
    ss_o = sum(d ** 2 for d in deviations_o)
    if ss_p == 0.0 or ss_o == 0.0:
        return 0.0
    cov = sum(dp * do for dp, do in zip(deviations_p, deviations_o))
    return cov / math.sqrt(ss_p * ss_o)


def within_abs_diff(
    predicted: Iterable[float],
    observed: Iterable[float],
    *,
    threshold: float = 0.5,
) -> float:
    """Fraction of predictions where |pred - obs| <= threshold (range [0, 1]).

    Parameters
    ----------
    predicted : iterable of float
    observed  : iterable of float
    threshold : float, keyword-only
        Maximum allowed absolute difference (default 0.5).

    Returns
    -------
    float
        Fraction in [0, 1].
    """
    preds = list(predicted)
    obs = list(observed)
    if len(preds) != len(obs):
        raise ValueError(
            f"within_abs_diff: length mismatch ({len(preds)} vs {len(obs)})"
        )
    if not preds:
        raise ValueError("within_abs_diff: empty input")
    hits = sum(1 for p, o in zip(preds, obs) if abs(p - o) <= threshold)
    return hits / len(preds)
