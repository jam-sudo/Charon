"""Tests for benchmark aggregation helpers."""

from __future__ import annotations

import math
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from validation.benchmarks.layer2_human_pk import (  # noqa: E402
    PanelRow,
    PanelSummary,
    aggregate_summary,
)


def _row(key: str, fold_cl: float, fold_vss: float, fold_thalf: float,
         strict: bool, override_tissues=None):
    """Convenience: synthesize a PanelRow from target fold errors.

    `observed` is always 1.0 for all metrics; `predicted` is chosen so
    that max(pred/obs, obs/pred) == target fold.
    """
    def pred_from_fold(f: float) -> float:
        # Always return pred > obs so fold = pred/obs = f
        return f
    predicted = {
        "cl_L_h": pred_from_fold(fold_cl),
        "vss_L": pred_from_fold(fold_vss),
        "t_half_h": pred_from_fold(fold_thalf),
    }
    observed = {"cl_L_h": 1.0, "vss_L": 1.0, "t_half_h": 1.0}
    return PanelRow(
        key=key,
        predicted=predicted,
        observed=observed,
        fold={"cl_L_h": fold_cl, "vss_L": fold_vss, "t_half_h": fold_thalf},
        pass_2_fold={
            "cl_L_h": fold_cl <= 2.0,
            "vss_L": fold_vss <= 2.0,
            "t_half_h": fold_thalf <= 2.0,
        },
        override_tissues=override_tissues or [],
        strict_targets=strict,
        mode="no_override",
    )


class TestAggregateSummary:
    def test_single_row(self):
        rows = [_row("theo", 1.0, 1.0, 1.0, strict=True)]
        s = aggregate_summary(rows, mode="no_override")
        assert s.n == 1
        assert s.mode == "no_override"
        assert s.aafe["cl_L_h"] == pytest.approx(1.0)
        assert s.aafe["vss_L"] == pytest.approx(1.0)
        assert s.aafe["t_half_h"] == pytest.approx(1.0)
        assert s.within_2_fold["cl_L_h"] == pytest.approx(1.0)
        assert s.strict_failures == 0

    def test_aafe_geometric_mean(self):
        rows = [
            _row("a", fold_cl=2.0, fold_vss=1.0, fold_thalf=1.0, strict=False),
            _row("b", fold_cl=8.0, fold_vss=1.0, fold_thalf=1.0, strict=False),
        ]
        s = aggregate_summary(rows, mode="no_override")
        # geometric mean of (2, 8) = 4
        assert s.aafe["cl_L_h"] == pytest.approx(4.0)

    def test_within_2_fold_fraction(self):
        rows = [
            _row("a", 1.5, 1.0, 1.0, strict=False),
            _row("b", 3.0, 1.0, 1.0, strict=False),
        ]
        s = aggregate_summary(rows, mode="no_override")
        assert s.within_2_fold["cl_L_h"] == pytest.approx(0.5)
        assert s.within_3_fold["cl_L_h"] == pytest.approx(1.0)

    def test_strict_failure_count(self):
        rows = [
            _row("a", 1.5, 1.0, 1.0, strict=True),     # PASS all metrics
            _row("b", 3.0, 1.0, 1.0, strict=True),     # FAIL cl_L_h fold
            _row("c", 1.0, 5.0, 1.0, strict=False),    # non-strict; doesn't count
        ]
        s = aggregate_summary(rows, mode="no_override")
        assert s.strict_failures == 1

    def test_mode_label_preserved(self):
        rows = [_row("a", 1.0, 1.0, 1.0, strict=False)]
        s_no = aggregate_summary(rows, mode="no_override")
        s_yes = aggregate_summary(rows, mode="with_override")
        assert s_no.mode == "no_override"
        assert s_yes.mode == "with_override"

    def test_empty_rows_returns_zero_n(self):
        s = aggregate_summary([], mode="no_override")
        assert s.n == 0
        # AAFE on empty is undefined; implementation returns NaN
        assert math.isnan(s.aafe["cl_L_h"])
