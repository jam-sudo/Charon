"""Sprint 16 unit tests for report_writer.py history-preservation feature.

Tests the sentinel marker preservation: when an existing .md file contains
`<!-- BEGIN_PRESERVED_HISTORY -->`, content from the marker through EOF is
preserved across `emit_report` calls.
"""

from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from validation.benchmarks.report_writer import (  # noqa: E402
    _HISTORY_MARKER,
    emit_report,
)


def _minimal_payload(title: str = "Test Report") -> dict:
    return {
        "title": title,
        "panel": "test_panel",
        "date_utc": "2026-04-15T00:00:00Z",
        "summary": {},
        "rows": [],
        "notes": ["test note"],
    }


def test_emit_report_preserves_history_round_trip(tmp_path: Path) -> None:
    """User content after marker is preserved across regeneration."""
    stem = tmp_path / "report"

    # Run 1: write fresh report (no marker, no history yet)
    emit_report(_minimal_payload("Run 1"), stem=stem)
    md_path = stem.with_suffix(".md")

    # Append marker + user history
    appended = (
        "\n"
        + _HISTORY_MARKER
        + "\n\n"
        + "## §1 Sprint X narrative\n\n"
        + "Manually authored historical content goes here.\n"
    )
    md_path.write_text(md_path.read_text(encoding="utf-8") + appended, encoding="utf-8")

    # Run 2: regenerate with different payload
    emit_report(_minimal_payload("Run 2"), stem=stem)

    final_content = md_path.read_text(encoding="utf-8")

    # New title from Run 2 must be present
    assert "Run 2" in final_content, "Run 2 title missing — generated content not written"
    # Run 1 title must be absent (data zone overwritten)
    assert "Run 1" not in final_content, "Run 1 title leaked — data zone not regenerated"
    # User history must be preserved
    assert "## §1 Sprint X narrative" in final_content, "User narrative section missing"
    assert "Manually authored historical content goes here." in final_content, (
        "User narrative body missing"
    )
    # Marker must be present exactly once
    assert final_content.count(_HISTORY_MARKER) == 1, (
        f"Expected marker once, got {final_content.count(_HISTORY_MARKER)}"
    )


def test_emit_report_no_marker_means_no_preservation(tmp_path: Path) -> None:
    """Without the marker, existing content is wiped (current behavior)."""
    stem = tmp_path / "report"

    # Run 1: fresh write
    emit_report(_minimal_payload("Run 1"), stem=stem)
    md_path = stem.with_suffix(".md")

    # Append stale content WITHOUT marker
    md_path.write_text(
        md_path.read_text(encoding="utf-8") + "\n## Stale\n\nOld content.\n",
        encoding="utf-8",
    )

    # Run 2: regenerate
    emit_report(_minimal_payload("Run 2"), stem=stem)

    final_content = md_path.read_text(encoding="utf-8")

    assert "Run 2" in final_content
    assert "## Stale" not in final_content, "Stale content survived without marker"
    assert "Old content." not in final_content
    assert _HISTORY_MARKER not in final_content, (
        "Marker should not be auto-added when not previously present"
    )


def test_emit_report_marker_only_idempotent(tmp_path: Path) -> None:
    """Marker present but no following content: regen is idempotent."""
    stem = tmp_path / "report"

    # Run 1: write + append marker only
    emit_report(_minimal_payload("Run 1"), stem=stem)
    md_path = stem.with_suffix(".md")
    md_path.write_text(
        md_path.read_text(encoding="utf-8") + "\n" + _HISTORY_MARKER + "\n",
        encoding="utf-8",
    )

    content_after_run1_with_marker = md_path.read_text(encoding="utf-8")

    # Run 2 with same payload should be idempotent
    emit_report(_minimal_payload("Run 1"), stem=stem)
    content_after_run2 = md_path.read_text(encoding="utf-8")

    # Marker still present exactly once
    assert content_after_run2.count(_HISTORY_MARKER) == 1
    # File content stable across regen with same payload + marker
    assert content_after_run1_with_marker == content_after_run2, (
        "Regen with marker-only history should be byte-identical to prior state"
    )


def test_emit_report_marker_first_occurrence_wins(tmp_path: Path) -> None:
    """If marker appears multiple times, first occurrence is the boundary."""
    stem = tmp_path / "report"

    emit_report(_minimal_payload("Run 1"), stem=stem)
    md_path = stem.with_suffix(".md")

    # Append two markers with content between them
    appended = (
        "\n"
        + _HISTORY_MARKER
        + "\n\n## A\n\nFirst section.\n\n"
        + _HISTORY_MARKER
        + "\n\n## B\n\nSecond section.\n"
    )
    md_path.write_text(md_path.read_text(encoding="utf-8") + appended, encoding="utf-8")

    # Regenerate
    emit_report(_minimal_payload("Run 2"), stem=stem)

    final_content = md_path.read_text(encoding="utf-8")

    # Both narrative sections preserved (extracted from first marker to EOF)
    assert "## A" in final_content
    assert "## B" in final_content
    assert "First section." in final_content
    assert "Second section." in final_content
    # Both markers preserved
    assert final_content.count(_HISTORY_MARKER) == 2
    # Run 2 generated content present
    assert "Run 2" in final_content
    assert "Run 1" not in final_content
