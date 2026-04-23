"""End-to-end: `charon report <smi>` Markdown/JSON contains Layer 1 CI.

Guards against silent regression where ``predict_properties`` somewhere
in the CLI chain gets called with ``CONFORMAL_OFF`` or the conformal
intervals get stripped before hitting the report.

Notes on CLI surface (adapted from the Task 5 spec to match reality):

- The dose flag is ``--dose`` (not ``--dose-mg``).
- ``charon report`` requires a dose-projection target; we supply
  ``--noael 5 --noael-species rat`` (MRSD = min of three methods).
- ``export_report`` always writes both a ``.md`` and a sibling ``.json``,
  so ``--json`` is not required.
- The JSON payload mirrors ``ReportData``: ADME entries live under the
  top-level ``properties`` dict keyed by property name (``fu_p``,
  ``clint_uL_min_mg`` ...), not under ``compound.adme``.  Each entry
  carries ``value`` / ``ci_lower`` / ``ci_upper`` / ``unit`` / ``source``.
"""
from __future__ import annotations

import json
import subprocess
import sys

CAFFEINE = "Cn1cnc2n(C)c(=O)n(C)c(=O)c12"


def _run_report(out_md) -> subprocess.CompletedProcess:
    return subprocess.run(
        [
            sys.executable, "-m", "charon.cli.main", "report",
            CAFFEINE,
            "--route", "iv_bolus",
            "--dose", "1.0",
            "--noael", "5",
            "--noael-species", "rat",
            "--output", str(out_md),
        ],
        capture_output=True,
        text=True,
        check=False,
    )


class TestReportCIPresent:
    def test_markdown_report_has_90pct_ci(self, tmp_path):
        out_md = tmp_path / "caffeine.md"
        result = _run_report(out_md)
        assert result.returncode == 0, result.stderr
        txt = out_md.read_text()
        # ADME table header renders a "90% CI" column.
        assert "90% CI" in txt
        # And the rendered rows contain at least one bracketed interval.
        assert "[" in txt and "]" in txt

    def test_json_report_has_ci_fields(self, tmp_path):
        out_md = tmp_path / "caffeine.md"
        out_json = tmp_path / "caffeine.json"
        result = _run_report(out_md)
        assert result.returncode == 0, result.stderr
        assert out_json.exists(), "JSON sibling was not written"

        data = json.loads(out_json.read_text())
        props = data["properties"]
        assert "fu_p" in props, f"fu_p missing from properties; keys={list(props)}"

        fup = props["fu_p"]
        assert fup.get("ci_lower") is not None, "fu_p ci_lower is None (CI stripped?)"
        assert fup.get("ci_upper") is not None, "fu_p ci_upper is None (CI stripped?)"
        # Sanity: interval must bracket the point prediction.
        assert fup["ci_lower"] <= fup["value"] <= fup["ci_upper"]

        # Also assert CLint carries a CI (core of the Sprint 7 wiring).
        clint = props.get("clint_uL_min_mg")
        assert clint is not None, "clint_uL_min_mg missing from properties"
        assert clint.get("ci_lower") is not None, "clint ci_lower is None"
        assert clint.get("ci_upper") is not None, "clint ci_upper is None"
