"""Unit tests for validation/benchmarks/report_writer.py."""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from validation.benchmarks.report_writer import emit_report  # noqa: E402

# ---------------------------------------------------------------------------
# Minimal valid payload fixture
# ---------------------------------------------------------------------------

MINIMAL_PAYLOAD: dict = {
    "title": "Test Report",
    "panel": "test_panel",
    "date_utc": "2026-04-15T00:00:00Z",
    "summary": {},
    "rows": [],
    "notes": ["All checks passed."],
}


def _full_payload(tmp_path: Path) -> dict:
    """Return a rich payload with optional keys for broader tests."""
    return {
        "title": "Layer 1 ADMET Validation",
        "panel": "obach_1999",
        "date_utc": "2026-04-15T12:00:00Z",
        "summary": [
            {"metric": "AAFE_CLint", "value": 2.31, "target": "<2.5", "met": True}
        ],
        "rows": [
            {
                "compound": "midazolam",
                "pred_CLint": 15.2,
                "obs_CLint": 14.0,
                "fold_error": 1.09,
            },
        ],
        "notes": ["Reference: Obach 1999.", "Scaffold CV used."],
        "disclaimer": "For research use only.",
        "targets": [
            {"metric": "AAFE", "target": "<2.5", "met": True}
        ],
        "excluded": [
            {"compound": "compound_X", "reason": "no reference data"}
        ],
    }


# ---------------------------------------------------------------------------
# Test 1 — both files created with correct suffixes
# ---------------------------------------------------------------------------

class TestCreatesFiles:
    def test_creates_md_and_json(self, tmp_path: Path):
        stem = tmp_path / "subdir" / "report"
        md_path, json_path = emit_report(MINIMAL_PAYLOAD, stem=stem)

        assert md_path.exists(), "Markdown file should exist"
        assert json_path.exists(), "JSON file should exist"
        assert md_path.suffix == ".md"
        assert json_path.suffix == ".json"
        assert md_path.is_absolute()
        assert json_path.is_absolute()

    def test_creates_parent_dirs(self, tmp_path: Path):
        stem = tmp_path / "a" / "b" / "c" / "report"
        md_path, json_path = emit_report(MINIMAL_PAYLOAD, stem=stem)
        assert md_path.exists()
        assert json_path.exists()


# ---------------------------------------------------------------------------
# Test 2 — missing required key raises ValueError
# ---------------------------------------------------------------------------

class TestMissingRequiredKey:
    def test_missing_required_key_raises(self, tmp_path: Path):
        bad_payload = dict(MINIMAL_PAYLOAD)
        del bad_payload["rows"]
        with pytest.raises(ValueError, match="missing required"):
            emit_report(bad_payload, stem=tmp_path / "report")

    def test_all_required_keys_checked(self, tmp_path: Path):
        """Each required key individually triggers ValueError when absent."""
        required_keys = ["title", "panel", "date_utc", "summary", "rows", "notes"]
        for key in required_keys:
            bad = dict(MINIMAL_PAYLOAD)
            del bad[key]
            with pytest.raises(ValueError, match="missing required"):
                emit_report(bad, stem=tmp_path / f"report_{key}")


# ---------------------------------------------------------------------------
# Test 3 — empty rows and empty summary allowed
# ---------------------------------------------------------------------------

class TestEmptyRowsAllowed:
    def test_empty_rows_allowed(self, tmp_path: Path):
        payload = dict(MINIMAL_PAYLOAD)
        payload["rows"] = []
        payload["summary"] = {}
        md_path, json_path = emit_report(payload, stem=tmp_path / "report")
        assert md_path.exists()
        assert json_path.exists()


# ---------------------------------------------------------------------------
# Test 4 — pipe character escaped in Markdown
# ---------------------------------------------------------------------------

class TestPipeEscape:
    def test_pipe_escape_in_md(self, tmp_path: Path):
        payload = dict(MINIMAL_PAYLOAD)
        payload["rows"] = [{"compound": "a|b", "value": 1.0}]
        md_path, _ = emit_report(payload, stem=tmp_path / "report")

        md_text = md_path.read_text(encoding="utf-8")
        assert r"a\|b" in md_text, "Pipe in cell must be escaped as \\|"
        assert "a|b" not in md_text.replace(r"a\|b", ""), (
            "Unescaped pipe should not appear as standalone cell content"
        )


# ---------------------------------------------------------------------------
# Test 5 — NaN renders as "-" in Markdown
# ---------------------------------------------------------------------------

class TestNanRendering:
    def test_nan_renders_as_dash_in_md(self, tmp_path: Path):
        payload = dict(MINIMAL_PAYLOAD)
        payload["rows"] = [{"compound": "testcpd", "pred": float("nan")}]
        md_path, _ = emit_report(payload, stem=tmp_path / "report")

        md_text = md_path.read_text(encoding="utf-8")
        # The NaN cell must render as "-"
        assert " - " in md_text or "| - " in md_text or "- |" in md_text, (
            "NaN should render as '-' in Markdown table"
        )
        assert "nan" not in md_text.lower(), "Raw 'nan' must not appear in Markdown"


# ---------------------------------------------------------------------------
# Test 6 — inf becomes null in JSON
# ---------------------------------------------------------------------------

class TestInfInJson:
    def test_inf_becomes_null_in_json(self, tmp_path: Path):
        payload = dict(MINIMAL_PAYLOAD)
        payload["rows"] = [{"compound": "testcpd", "value": float("inf")}]
        _, json_path = emit_report(payload, stem=tmp_path / "report")

        raw = json_path.read_text(encoding="utf-8")
        # Must be valid JSON (no Infinity literals)
        data = json.loads(raw)
        rows = data["rows"]
        assert rows[0]["value"] is None, "inf must be serialised as null in JSON"


# ---------------------------------------------------------------------------
# Test 7 — JSON roundtrip preserves structure
# ---------------------------------------------------------------------------

class TestJsonRoundtrip:
    def test_json_roundtrip(self, tmp_path: Path):
        payload = _full_payload(tmp_path)
        _, json_path = emit_report(payload, stem=tmp_path / "report")

        data = json.loads(json_path.read_text(encoding="utf-8"))
        assert data["title"] == "Layer 1 ADMET Validation"
        assert data["panel"] == "obach_1999"
        assert isinstance(data["rows"], list)
        assert len(data["rows"]) == 1
        assert data["rows"][0]["compound"] == "midazolam"
        # summary preserved (list form)
        assert isinstance(data["summary"], list)
        assert data["summary"][0]["metric"] == "AAFE_CLint"


# ---------------------------------------------------------------------------
# Test 8 — disclaimer rendered as blockquote when present
# ---------------------------------------------------------------------------

class TestDisclaimerRendering:
    def test_disclaimer_rendered_when_present(self, tmp_path: Path):
        payload = dict(MINIMAL_PAYLOAD)
        payload["disclaimer"] = "For research use only."
        md_path, _ = emit_report(payload, stem=tmp_path / "report")

        md_text = md_path.read_text(encoding="utf-8")
        assert "For research use only." in md_text
        # Should appear in a blockquote ("> " prefix)
        assert "> " in md_text, "Disclaimer must be rendered as a blockquote"

    def test_no_disclaimer_section_when_absent(self, tmp_path: Path):
        payload = dict(MINIMAL_PAYLOAD)
        assert "disclaimer" not in payload
        md_path, _ = emit_report(payload, stem=tmp_path / "report")

        md_text = md_path.read_text(encoding="utf-8")
        # No blockquote line should appear
        lines_with_blockquote = [ln for ln in md_text.splitlines() if ln.startswith(">")]
        assert not lines_with_blockquote, "No blockquote expected without disclaimer"
