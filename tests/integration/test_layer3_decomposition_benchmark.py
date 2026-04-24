"""Integration tests for the Sprint 10 decomposition orchestrator.

Runs the full 12-compound Tier A decomposition and verifies structural
properties of the output.
"""

from __future__ import annotations

import json
import math
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "validation" / "benchmarks" / "layer3_ivive_decomposition.py"


@pytest.fixture(scope="module")
def decomposition_json(tmp_path_factory) -> dict:
    """Run the orchestrator once and cache the JSON output for all tests."""
    out_stem = tmp_path_factory.mktemp("sprint10") / "decomp"
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "--output-stem", str(out_stem)],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        timeout=600,
    )
    assert result.returncode == 0, (
        f"orchestrator failed: stdout={result.stdout} stderr={result.stderr}"
    )
    json_path = out_stem.with_suffix(".json")
    assert json_path.exists()
    return json.loads(json_path.read_text())


def test_decomposition_covers_all_12_tier_a(decomposition_json):
    panel_path = REPO_ROOT / "validation" / "data" / "fih_reference" / "panel.yaml"
    panel = yaml.safe_load(panel_path.read_text())["panel"]
    tier_a = [c["name"] for c in panel["compounds"] if c["tier"] == "gold"]
    rows = decomposition_json["extra_sections"]["Per-compound decomposition"]
    names_in_rows = {r["compound"] for r in rows}
    assert names_in_rows == set(tier_a), (
        f"missing: {set(tier_a) - names_in_rows}, "
        f"extra: {names_in_rows - set(tier_a)}"
    )
    assert len(rows) == 12


def test_decomposition_log_additivity_per_compound(decomposition_json):
    """For each of 12 compounds: log10 of signed observed fold equals
    the log10 sum of (liver, route, residual) signed factors within rtol=1e-6."""
    rows = decomposition_json["extra_sections"]["Per-compound decomposition"]
    for row in rows:
        lhs = math.log10(row["fold_observed_signed"])
        rhs = (
            math.log10(row["fold_liver_model_signed"])
            + math.log10(row["fold_route_bias"])
            + math.log10(row["fold_residual_signed"])
        )
        assert lhs == pytest.approx(rhs, rel=1e-6, abs=1e-9), (
            f"{row['compound']}: log-additivity violated "
            f"(lhs={lhs}, rhs={rhs})"
        )


def test_decomposition_every_row_numeric(decomposition_json):
    """No NaN / null in required numeric fields."""
    required = {
        "fold_observed",
        "fold_liver_model",
        "fold_route_bias",
        "fold_residual",
    }
    rows = decomposition_json["extra_sections"]["Per-compound decomposition"]
    for row in rows:
        for key in required:
            v = row[key]
            assert isinstance(v, (int, float))
            assert not math.isnan(v), f"{row['compound']}: {key} is NaN"


def test_liver_model_whatif_produces_nonzero_attribution(decomposition_json):
    """Regression guard: afb2464 had aggregate_pct_liver_model = 0% because
    Pipeline(liver_model=X) is a no-op for the PBPK ODE. Analytical scaling
    (a4c5b5d+) must produce non-zero liver-model attribution for a panel
    whose compounds span clint values from near-zero to high extraction.
    """
    summary = decomposition_json["summary"]
    assert summary["aggregate_pct_liver_model"] > 0.0, (
        "aggregate_pct_liver_model is 0% — "
        "likely regression to per-model Pipeline calls (Pipeline.liver_model "
        "is a no-op for the PBPK ODE)"
    )


def test_liver_model_diverges_for_majority_of_panel(decomposition_json):
    """At least 5 of 12 Tier A compounds should show non-trivial liver-model
    divergence (|log10(fold_liver_model)| > 0.001). This threshold is a
    regression guard against the Sprint 10 bug class (aggregate_pct_liver_model
    = 0% because Pipeline.liver_model was a no-op), NOT a pharmacological
    claim. The exact count depends on simulation route:
      - Sprint 10 (iv_bolus): ~10/12 diverged (non-flow-limited IV MRSDs)
      - Sprint 11 (oral):     ~7/12 diverged (more flow-limited under oral)
    Compounds with near-zero CLint (lisinopril, diazepam) always converge.
    """
    rows = decomposition_json["extra_sections"]["Per-compound decomposition"]
    n_divergent = sum(
        1 for r in rows
        if r["fold_liver_model"] > 0 and abs(math.log10(r["fold_liver_model"])) > 0.001
    )
    assert n_divergent >= 5, (
        f"Only {n_divergent}/12 compounds show liver-model divergence; "
        f"expected >= 5 on this panel. If 0, likely regression to "
        f"Pipeline.liver_model no-op (see Sprint 10 bug class)."
    )
