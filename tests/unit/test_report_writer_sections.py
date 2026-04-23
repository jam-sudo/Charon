"""Tests for emit_report extra_sections support."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from validation.benchmarks.report_writer import emit_report  # noqa: E402


def _base_payload():
    return {
        "title": "Test",
        "panel": "test",
        "date_utc": "2026-04-23T00:00:00+00:00",
        "summary": {"x": 1},
        "rows": [{"name": "a", "value": 1.0}],
        "notes": ["note"],
    }


def test_extra_sections_rendered_in_markdown(tmp_path):
    payload = _base_payload()
    payload["extra_sections"] = {
        "Conformal coverage": [
            {"property": "fup", "n_samples": 153, "empirical_coverage": 0.902},
            {"property": "clint_hepatocyte", "n_samples": 1441, "empirical_coverage": 0.890},
        ]
    }
    md_path, json_path = emit_report(payload, stem=tmp_path / "rep")
    md = md_path.read_text()
    assert "## Conformal coverage" in md
    assert "fup" in md
    assert "clint_hepatocyte" in md


def test_extra_sections_preserved_in_json(tmp_path):
    payload = _base_payload()
    payload["extra_sections"] = {"Foo": [{"a": 1}]}
    md_path, json_path = emit_report(payload, stem=tmp_path / "rep")
    data = json.loads(json_path.read_text())
    assert data["extra_sections"]["Foo"] == [{"a": 1}]


def test_rows_optional_when_extra_sections_present(tmp_path):
    """Benchmarks like Layer 3 have no single ``rows`` list."""
    payload = _base_payload()
    payload["rows"] = []  # allow empty
    payload["extra_sections"] = {
        "Gold results": [{"compound": "midazolam", "fold": 1.2}],
        "Sanity floor": [{"compound": "warfarin", "pass_floor": True}],
    }
    md_path, _ = emit_report(payload, stem=tmp_path / "rep")
    md = md_path.read_text()
    assert "## Gold results" in md
    assert "## Sanity floor" in md
