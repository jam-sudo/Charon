"""Unit tests for layer3_fih_dose benchmark helpers."""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]


class TestFoldError:
    def test_symmetric_fold_error(self):
        from validation.benchmarks.layer3_fih_dose import fold_error
        assert fold_error(10.0, 5.0) == pytest.approx(2.0)
        assert fold_error(5.0, 10.0) == pytest.approx(2.0)

    def test_identity(self):
        from validation.benchmarks.layer3_fih_dose import fold_error
        assert fold_error(7.0, 7.0) == pytest.approx(1.0)


class TestSanityFloor:
    def test_pass_when_pred_below(self):
        from validation.benchmarks.layer3_fih_dose import passes_floor
        assert passes_floor(predicted_mg=1.0, approved_starting_mg=2.0)

    def test_pass_at_equal(self):
        from validation.benchmarks.layer3_fih_dose import passes_floor
        assert passes_floor(predicted_mg=2.0, approved_starting_mg=2.0)

    def test_fail_above(self):
        from validation.benchmarks.layer3_fih_dose import passes_floor
        assert not passes_floor(predicted_mg=3.0, approved_starting_mg=2.0)


@pytest.mark.slow
class TestEndToEndRun:
    def test_script_runs_and_emits_report(self, tmp_path):
        stem = tmp_path / "l3"
        result = subprocess.run(
            [
                sys.executable,
                str(REPO_ROOT / "validation" / "benchmarks" / "layer3_fih_dose.py"),
                "--output-stem", str(stem),
            ],
            check=False,
            capture_output=True,
            text=True,
        )
        # Whether exit 0 or 2, the report must be written.
        data = json.loads(stem.with_suffix(".json").read_text())
        assert "extra_sections" in data
        sections = data["extra_sections"]
        gold_key = next(k for k in sections if k.lower().startswith("gold"))
        floor_key = next(k for k in sections if k.lower().startswith("sanity"))
        floor_rows = sections[floor_key]
        all_pass = all(row["pass_floor"] for row in floor_rows)
        assert result.returncode == 0, (
            f"Layer 3 benchmark exited {result.returncode}. Sanity floor failures: "
            f"{[r['compound'] for r in floor_rows if not r['pass_floor']]}"
        )
        assert all_pass
