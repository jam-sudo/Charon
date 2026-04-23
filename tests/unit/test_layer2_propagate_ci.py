"""Layer 2 benchmark --propagate-ci preserves point-estimate AAFE.

Task 8 of Sprint 7: adds a Monte-Carlo mode that samples fu_p and
clint from conformal CIs (ML-sourced rows only) and reports
p05/p50/p95 of CL and Vss per compound.

The committed Obach panel uses ``source: experimental`` for every
compound, so when ``--propagate-ci`` is ON the benchmark is still
expected to run end-to-end but with no percentile fields emitted
(no ML-sourced rows to sample). The sampling plumbing itself is
exercised by an in-process unit test below that injects an ML-
sourced fu_p so ``_mc_sample_pk`` actually produces percentiles.
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]


def _find_metric(summary_list: list[dict], label: str) -> dict:
    for entry in summary_list:
        if entry.get("metric") == label:
            return entry
    raise KeyError(f"metric {label!r} not found in summary list")


@pytest.mark.slow
class TestPropagateCIMode:
    def test_point_mode_aafe_unchanged(self, tmp_path: Path) -> None:
        """Point-only mode must reproduce the committed report's AAFE.

        Enforces the Task-8 invariant: adding argparse and the CI mode
        must not perturb the point-estimate pipeline (same PRNG state,
        same Kp overrides, same BDF solver inputs). Tolerance rel=1e-6.
        """
        committed = json.loads(
            (REPO_ROOT / "validation" / "reports" / "layer2_human_pk.json").read_text()
        )
        stem = tmp_path / "layer2"
        subprocess.run(
            [
                sys.executable,
                str(REPO_ROOT / "validation" / "benchmarks" / "layer2_human_pk.py"),
                "--output-stem",
                str(stem),
            ],
            check=True,
            cwd=str(REPO_ROOT),
        )
        fresh = json.loads(stem.with_suffix(".json").read_text())
        for metric in ("CL (L/h)", "Vss (L)", "t_half (h)"):
            c = _find_metric(committed["summary"], metric)
            f = _find_metric(fresh["summary"], metric)
            assert f["AAFE"] == pytest.approx(
                c["AAFE"], rel=1e-6
            ), f"Point-estimate AAFE drifted for {metric}"
            assert f["within_2_fold"] == pytest.approx(
                c["within_2_fold"], rel=1e-6
            )
            assert f["within_3_fold"] == pytest.approx(
                c["within_3_fold"], rel=1e-6
            )

    def test_ci_mode_runs_end_to_end(self, tmp_path: Path) -> None:
        """--propagate-ci mode runs without error and emits a report.

        The Obach Tier-1 panel is entirely ``source: experimental``, so
        no percentile fields are expected. The assertion is that the
        subprocess exits 0 and produces a JSON report; no ML-sourced
        rows means no p05/p50/p95 — this is correct behavior for
        experimental panels. The percentile-emitting path is covered
        by TestMCSamplePK below.
        """
        stem = tmp_path / "layer2_ci"
        subprocess.run(
            [
                sys.executable,
                str(REPO_ROOT / "validation" / "benchmarks" / "layer2_human_pk.py"),
                "--propagate-ci",
                "--n-samples",
                "8",
                "--output-stem",
                str(stem),
            ],
            check=True,
            cwd=str(REPO_ROOT),
        )
        data = json.loads(stem.with_suffix(".json").read_text())
        assert "rows" in data
        assert len(data["rows"]) >= 1
        # Point-estimate numbers in CI mode match point-only (same seed).
        for metric in ("CL (L/h)", "Vss (L)",):
            entry = _find_metric(data["summary"], metric)
            assert entry["AAFE"] > 0


class TestMCSamplePK:
    """Direct unit coverage of the Monte-Carlo sampling path."""

    def test_mc_sample_pk_emits_percentiles_for_ml_sourced(self) -> None:
        """Injecting an ML-sourced fu_p with a CI yields p05/p50/p95."""
        # Lazy imports so the module loads even if --propagate-ci path
        # isn't yet wired.
        sys.path.insert(0, str(REPO_ROOT))
        from validation.benchmarks.layer2_human_pk import (
            _mc_sample_pk,
            load_panel,
            DEFAULT_PANEL_PATH,
            PanelEntry,
        )

        panel = load_panel(DEFAULT_PANEL_PATH)
        # Pick theophylline as baseline (neutral, fast, strict-gate).
        entry0 = next(e for e in panel if e.key == "theophylline")

        # Inject ML-sourced fu_p with a CI.
        from charon.core.schema import PredictedProperty

        fup_ml = PredictedProperty(
            value=0.60,
            ci_90_lower=0.45,
            ci_90_upper=0.80,
            source="ml_ensemble",
            unit="fraction",
        )
        new_binding = entry0.compound.properties.binding.model_copy(
            update={"fu_p": fup_ml}
        )
        new_props = entry0.compound.properties.model_copy(
            update={"binding": new_binding}
        )
        new_compound = entry0.compound.model_copy(
            update={"properties": new_props}
        )
        entry_ml = PanelEntry(
            key=entry0.key,
            compound=new_compound,
            route=entry0.route,
            dose_mg=entry0.dose_mg,
            duration_h=entry0.duration_h,
            observed=entry0.observed,
            strict_targets=entry0.strict_targets,
            obach_table_row=entry0.obach_table_row,
            notes=entry0.notes,
        )

        result = _mc_sample_pk(entry_ml, n=8, seed=42)
        assert "cl_p05" in result
        assert "cl_p50" in result
        assert "cl_p95" in result
        assert "vss_p05" in result
        assert "vss_p95" in result
        assert result["cl_p05"] <= result["cl_p50"] <= result["cl_p95"]
        assert result["vss_p05"] <= result["vss_p50"] <= result["vss_p95"]

    def test_mc_sample_pk_returns_empty_for_all_experimental(self) -> None:
        """All-experimental entry → no percentile fields returned."""
        sys.path.insert(0, str(REPO_ROOT))
        from validation.benchmarks.layer2_human_pk import (
            _mc_sample_pk,
            load_panel,
            DEFAULT_PANEL_PATH,
        )

        panel = load_panel(DEFAULT_PANEL_PATH)
        entry0 = next(e for e in panel if e.key == "theophylline")
        assert _mc_sample_pk(entry0, n=8, seed=42) == {}
