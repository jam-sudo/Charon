"""Layer 1 benchmark exposes per-property conformal coverage."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

REPORT = Path(__file__).resolve().parents[2] / "validation" / "reports" / "layer1_admet.json"


@pytest.mark.skipif(not REPORT.exists(), reason="Layer 1 report not generated yet")
class TestCoverageSection:
    def test_json_has_conformal_section(self):
        data = json.loads(REPORT.read_text())
        assert "extra_sections" in data
        coverage_rows = data["extra_sections"].get("Conformal coverage", [])
        assert len(coverage_rows) >= 2
        props = {r["property"] for r in coverage_rows}
        assert {"fup", "clint_hepatocyte"} <= props
        for row in coverage_rows:
            assert row["n_samples"] >= 100
            assert 0.0 < row["empirical_coverage"] <= 1.0
            assert row["mean_fold_error"] >= 1.0

    def test_markdown_has_coverage_section(self):
        txt = REPORT.with_suffix(".md").read_text()
        assert "## Conformal coverage" in txt
