"""Benchmark metrics for Layer 1-3 validation suites.

Metrics
-------
fold_error    : max(pred/obs, obs/pred), always >= 1.0
aafe          : geometric mean of fold errors (absolute average fold error)
within_n_fold : fraction of predictions with fold_error <= n
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
