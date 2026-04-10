"""End-to-end Obach panel smoke test.

Loads the real panel.yaml, runs the full benchmark, and asserts:
  - All 10 compounds produce finite positive PK values in both modes.
  - theophylline strict 2-fold gate passes (regression invariant).
  - Panel AAFE_Vss with-override stays under the 5.0 sanity floor.
  - ≥2 compounds have recorded kp_overrides (DoD §2 check).

The actual AAFE values are `record_property`'d into pytest results
so CI artifacts preserve the measurement.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from validation.benchmarks.layer2_human_pk import (  # noqa: E402
    DEFAULT_PANEL_PATH,
    run_benchmark,
)


@pytest.fixture(scope="module")
def benchmark_result():
    summaries, rows = run_benchmark(DEFAULT_PANEL_PATH)
    return summaries, rows


class TestObachPanelSmoke:
    def test_all_compounds_produce_finite_pk_no_override(self, benchmark_result):
        _, rows = benchmark_result
        for r in rows["no_override"]:
            for metric in ("cl_L_h", "vss_L", "t_half_h"):
                assert math.isfinite(r.predicted[metric]), (
                    f"{r.key} {metric} non-finite: {r.predicted[metric]}"
                )
                assert r.predicted[metric] > 0, (
                    f"{r.key} {metric} non-positive: {r.predicted[metric]}"
                )

    def test_all_compounds_produce_finite_pk_with_override(self, benchmark_result):
        _, rows = benchmark_result
        for r in rows["with_override"]:
            for metric in ("cl_L_h", "vss_L", "t_half_h"):
                assert math.isfinite(r.predicted[metric]), (
                    f"{r.key} {metric} non-finite: {r.predicted[metric]}"
                )
                assert r.predicted[metric] > 0, (
                    f"{r.key} {metric} non-positive: {r.predicted[metric]}"
                )

    def test_theophylline_strict_gate_passes(self, benchmark_result):
        _, rows = benchmark_result
        theo_rows = [r for r in rows["with_override"] if r.key == "theophylline"]
        assert len(theo_rows) == 1
        theo = theo_rows[0]
        assert theo.strict_targets is True
        for metric in ("cl_L_h", "vss_L", "t_half_h"):
            assert theo.pass_2_fold[metric], (
                f"theophylline {metric} fold={theo.fold[metric]:.2f} > 2.0"
            )

    def test_panel_aafe_vss_below_sanity_floor(self, benchmark_result, record_property):
        summaries, _ = benchmark_result
        aafe_vss_yes = summaries["with_override"].aafe["vss_L"]
        aafe_vss_no = summaries["no_override"].aafe["vss_L"]
        record_property("aafe_vss_with_override", aafe_vss_yes)
        record_property("aafe_vss_no_override", aafe_vss_no)
        assert aafe_vss_yes < 5.0, (
            f"Session post-condition failed: AAFE_Vss with-override = "
            f"{aafe_vss_yes:.2f} exceeds 5.0 sanity floor"
        )

    def test_panel_aafe_cl_below_sanity_floor(self, benchmark_result, record_property):
        summaries, _ = benchmark_result
        aafe_cl_yes = summaries["with_override"].aafe["cl_L_h"]
        aafe_cl_no = summaries["no_override"].aafe["cl_L_h"]
        record_property("aafe_cl_with_override", aafe_cl_yes)
        record_property("aafe_cl_no_override", aafe_cl_no)
        # CL sanity floor: 5.0 (same as Vss, generous)
        assert aafe_cl_yes < 5.0

    def test_at_least_two_compounds_have_verified_overrides(self, benchmark_result):
        """DoD §2: ≥2 compounds must carry empirical Kp overrides."""
        _, rows = benchmark_result
        compounds_with_overrides = {
            r.key for r in rows["with_override"] if r.override_tissues
        }
        assert len(compounds_with_overrides) >= 2, (
            f"DoD §2 violation: only {len(compounds_with_overrides)} "
            f"compound(s) have verified overrides: {compounds_with_overrides}"
        )

    def test_within_2_fold_reported(self, benchmark_result, record_property):
        summaries, _ = benchmark_result
        for mode in ("no_override", "with_override"):
            for metric in ("cl_L_h", "vss_L", "t_half_h"):
                frac = summaries[mode].within_2_fold[metric]
                record_property(f"within_2_fold__{mode}__{metric}", frac)
                assert 0.0 <= frac <= 1.0
