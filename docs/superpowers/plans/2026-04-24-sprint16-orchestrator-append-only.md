# Sprint 16 — Layer 3 Report Orchestrator Append-Only Refactor Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add sentinel marker preservation (`<!-- BEGIN_PRESERVED_HISTORY -->`) to `validation/benchmarks/report_writer.py:emit_report()` so manually-appended sprint narrative sections survive benchmark regeneration.

**Architecture:** Single-function refactor. Marker presence in existing `.md` is the opt-in signal (no API change). When marker is found, content from marker through EOF is preserved; otherwise current overwrite behavior unchanged. Backwards compatible — Layer 1/2 reports unaffected. Layer 3 reports get marker manually added once.

**Tech Stack:** Python 3.11+, pytest, pathlib. No new dependencies.

**Spec:** `docs/superpowers/specs/2026-04-24-sprint16-orchestrator-append-only-design.md`

---

## File Structure

| File | Change |
|---|---|
| `validation/benchmarks/report_writer.py` | Modify — add `_HISTORY_MARKER` constant + `_extract_preserved_history()` helper + integrate into `emit_report()` (~25 LOC) |
| `tests/unit/test_report_writer_preservation.py` | NEW — 4 unit tests (~110 LOC) |
| `validation/reports/layer3_fih_dose.md` | Modify — insert marker before first sprint narrative section (1 line) |
| `validation/reports/layer3_ivive_decomposition.md` | Modify — insert marker before §9 Sprint 13 narrative (1 line) |

**Unchanged (reaffirming):**
- All `src/charon/` code
- `validation/benchmarks/layer3_fih_dose.py`, `layer3_ivive_decomposition.py`, `layer1_admet.py`, `layer2_human_pk.py`
- All compound YAMLs, `panel.yaml`
- All other unit + integration tests

---

## Current `report_writer.py` shape (for reference)

The relevant existing code (lines 247-291) is the `emit_report()` function:

```python
def emit_report(payload: dict, *, stem: Path) -> tuple[Path, Path]:
    """Emit benchmark report as ``{stem}.md`` and ``{stem}.json``."""
    # Validate required keys
    missing = [k for k in _REQUIRED_KEYS if k not in payload]
    if missing:
        raise ValueError(...)

    stem = Path(stem).resolve()
    stem.parent.mkdir(parents=True, exist_ok=True)

    md_path = stem.with_suffix(".md")
    json_path = stem.with_suffix(".json")

    # Write Markdown
    md_content = _render_markdown(payload)
    md_path.write_text(md_content, encoding="utf-8")  # ← LINE 284, the overwrite

    # Write JSON
    sanitised = _sanitise(payload)
    json_content = json.dumps(sanitised, indent=2, allow_nan=False)
    json_path.write_text(json_content, encoding="utf-8")

    return md_path, json_path
```

The Sprint 16 change inserts marker preservation logic between `_render_markdown()` and `md_path.write_text()`.

---

## Task 1: Implement marker preservation in `report_writer.py`

**Files:**
- Create: `tests/unit/test_report_writer_preservation.py`
- Modify: `validation/benchmarks/report_writer.py`

This task uses TDD: write failing tests first, then implement.

- [ ] **Step 1: Write the 4 failing tests**

Create `tests/unit/test_report_writer_preservation.py` with this exact content:

```python
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
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
cd /home/jam/Charon/.worktrees/sprint16-orchestrator-append-only
pytest tests/unit/test_report_writer_preservation.py -v
```

Expected: 4 tests fail with `ImportError: cannot import name '_HISTORY_MARKER' from 'validation.benchmarks.report_writer'`. This confirms the helper and constant don't yet exist.

- [ ] **Step 3: Add `_HISTORY_MARKER` constant + `_extract_preserved_history` helper**

Modify `validation/benchmarks/report_writer.py`. Find this section near the top of the file (right after the docstring closing `"""` and the `from __future__ import annotations` import block, around line 47):

```python
# ---------------------------------------------------------------------------
# Required payload keys
# ---------------------------------------------------------------------------

_REQUIRED_KEYS: tuple[str, ...] = (
    "title",
    "panel",
    "date_utc",
    "summary",
    "rows",
    "notes",
)
```

Add immediately AFTER the `_REQUIRED_KEYS` block (before the `# Cell formatting helpers` separator):

```python


# ---------------------------------------------------------------------------
# History preservation (Sprint 16)
# ---------------------------------------------------------------------------

# Sentinel marker for opt-in history preservation. When an existing .md file
# contains this marker, emit_report preserves content from the marker through
# end-of-file across regenerations. Place it at the boundary between
# orchestrator-generated content (above) and user-managed sprint narratives
# (below). Marker absence = current overwrite behavior (backwards compatible).
_HISTORY_MARKER: str = "<!-- BEGIN_PRESERVED_HISTORY -->"


def _extract_preserved_history(md_path: Path) -> str:
    """Return content from history marker to EOF (inclusive of marker), or ''.

    Returns empty string if the file does not exist OR does not contain the
    marker. The marker itself is included in the returned content so it
    persists across regenerations.

    Parameters
    ----------
    md_path:
        Path to the existing markdown report file.

    Returns
    -------
    str
        Preserved content (marker + everything after it) or ``""`` if no
        preservation should occur.
    """
    if not md_path.exists():
        return ""
    content = md_path.read_text(encoding="utf-8")
    idx = content.find(_HISTORY_MARKER)
    if idx == -1:
        return ""
    return content[idx:]
```

- [ ] **Step 4: Integrate preservation into `emit_report()`**

Find the `emit_report()` function (line 247 in current code). Locate this block (around line 282-284):

```python
    md_path = stem.with_suffix(".md")
    json_path = stem.with_suffix(".json")

    # Write Markdown
    md_content = _render_markdown(payload)
    md_path.write_text(md_content, encoding="utf-8")
```

Replace it with:

```python
    md_path = stem.with_suffix(".md")
    json_path = stem.with_suffix(".json")

    # Write Markdown — preserve history below the sentinel marker if present
    preserved = _extract_preserved_history(md_path)
    md_content = _render_markdown(payload)
    if preserved:
        md_content = md_content.rstrip() + "\n\n" + preserved
    md_path.write_text(md_content, encoding="utf-8")
```

- [ ] **Step 5: Run the 4 new tests to verify they pass**

```bash
pytest tests/unit/test_report_writer_preservation.py -v
```

Expected: 4/4 PASS. Output should show all four test functions with PASSED status.

If any fail: re-read the test assertions vs implementation. Common issues:
- Marker not exported (check `_HISTORY_MARKER` is module-level, not inside a function)
- `_extract_preserved_history` returns wrong value (should include marker itself in returned string)
- `emit_report` doesn't append `preserved` correctly (should rstrip new content first, then `\n\n` + preserved)

- [ ] **Step 6: Run existing report_writer tests to verify no regression**

```bash
pytest tests/unit/test_report_writer.py tests/unit/test_report_writer_sections.py -q
```

Expected: all pass. The pre-existing test suite uses fresh `tmp_path` fixtures, so they don't have markers in their files → preservation doesn't trigger → behavior unchanged.

If any pre-existing test fails: STOP, report BLOCKED. Likely cause: introduced a bug in the rstrip / newline handling that affects basic write behavior.

- [ ] **Step 7: Commit**

```bash
cd /home/jam/Charon/.worktrees/sprint16-orchestrator-append-only
git add validation/benchmarks/report_writer.py tests/unit/test_report_writer_preservation.py
git commit -m "feat(sprint16): add sentinel marker history preservation to emit_report"
```

---

## Task 2: Add markers to Layer 3 reports

**Files:**
- Modify: `validation/reports/layer3_fih_dose.md`
- Modify: `validation/reports/layer3_ivive_decomposition.md`

The marker is added once, manually, immediately before the first sprint narrative section in each Layer 3 report. After this task, future orchestrator regenerations will preserve all narrative content below the marker.

- [ ] **Step 1: Locate insertion point in `layer3_fih_dose.md`**

The current end-of-orchestrator-content is after the `## Notes` section. Find where the first sprint narrative begins:

```bash
cd /home/jam/Charon/.worktrees/sprint16-orchestrator-append-only
grep -n "^## " validation/reports/layer3_fih_dose.md
```

Look for the first `## Sprint` section heading (or `## §`). The marker goes on its own line, with blank lines around it, immediately BEFORE that heading.

- [ ] **Step 2: Insert marker into `layer3_fih_dose.md`**

Use the Edit tool to insert the marker. The exact context depends on the current file (which may have Sprint 15 narrative as the first appended section). Find the line that starts the FIRST sprint narrative section (e.g., `## Sprint 15 comparison...`). Insert above it:

Pattern: locate the last line of the orchestrator-managed `## Notes` section content, then insert the marker on the next line. Example:

If the file currently has:
```
## Notes

- Decomposition: ...
- ...

## Sprint 15 comparison (CYP2D6 enhancement for propranolol — partial closure)
```

Edit to insert `<!-- BEGIN_PRESERVED_HISTORY -->` between the last Notes bullet and the first sprint section:

```
## Notes

- Decomposition: ...
- ...

<!-- BEGIN_PRESERVED_HISTORY -->

## Sprint 15 comparison (CYP2D6 enhancement for propranolol — partial closure)
```

Use the Edit tool with `old_string` matching enough surrounding context to be unique (typically the last note bullet + blank line + first sprint heading) and `new_string` containing the marker line + blank line inserted.

- [ ] **Step 3: Verify marker is in place in `layer3_fih_dose.md`**

```bash
grep -n "BEGIN_PRESERVED_HISTORY" validation/reports/layer3_fih_dose.md
```

Expected: exactly one match. The line number should be after `## Notes` and before any `## Sprint` heading.

- [ ] **Step 4: Locate insertion point in `layer3_ivive_decomposition.md`**

```bash
grep -n "^## " validation/reports/layer3_ivive_decomposition.md
```

Find the first sprint narrative section (likely `## §9. Sprint 13...` or `## §10 Sprint 15 audit...`). The marker goes immediately before that.

- [ ] **Step 5: Insert marker into `layer3_ivive_decomposition.md`**

Same procedure as Step 2 but for this file. Insert `<!-- BEGIN_PRESERVED_HISTORY -->` between the orchestrator-managed end (after the last `## Notes` bullet) and the first sprint narrative section.

- [ ] **Step 6: Verify marker is in place in `layer3_ivive_decomposition.md`**

```bash
grep -n "BEGIN_PRESERVED_HISTORY" validation/reports/layer3_ivive_decomposition.md
```

Expected: exactly one match. Line number should be after `## Notes` and before any `## §N` or `## Sprint` heading.

- [ ] **Step 7: Commit**

```bash
git add validation/reports/layer3_fih_dose.md validation/reports/layer3_ivive_decomposition.md
git commit -m "chore(sprint16): add history preservation markers to Layer 3 reports"
```

---

## Task 3: Verification — regenerate Layer 3 benchmarks + full suite

**Files:**
- Re-generate (transparently — no code change in benchmark scripts): `validation/reports/layer3_fih_dose.{md,json}`, `validation/reports/layer3_ivive_decomposition.{md,json}`

This task verifies that running the orchestrators after Tasks 1+2 preserves the narrative content as designed.

- [ ] **Step 1: Capture pre-regen narrative section count**

```bash
cd /home/jam/Charon/.worktrees/sprint16-orchestrator-append-only

# Count narrative sections (## Sprint or ## §) in fih_dose.md
echo "fih_dose narrative section count BEFORE regen:"
grep -c "^## \(Sprint\|§\)" validation/reports/layer3_fih_dose.md

echo "decomposition narrative section count BEFORE regen:"
grep -c "^## \(Sprint\|§\)" validation/reports/layer3_ivive_decomposition.md
```

Record the two counts (one for each file).

- [ ] **Step 2: Run Layer 3 benchmarks**

```bash
python3 validation/benchmarks/layer3_fih_dose.py
python3 validation/benchmarks/layer3_ivive_decomposition.py
```

Expected output:
- `[OK] Sanity floor: 12/12 pass.`
- `[OK] Decomposition wrote ... (12 compounds)`

No errors. The benchmark scripts are unchanged — they go through the modified `emit_report()` which now sees the marker in the existing files and preserves content.

- [ ] **Step 3: Verify narrative sections preserved**

```bash
echo "fih_dose narrative section count AFTER regen:"
grep -c "^## \(Sprint\|§\)" validation/reports/layer3_fih_dose.md

echo "decomposition narrative section count AFTER regen:"
grep -c "^## \(Sprint\|§\)" validation/reports/layer3_ivive_decomposition.md

echo "Marker count in fih_dose:"
grep -c "BEGIN_PRESERVED_HISTORY" validation/reports/layer3_fih_dose.md

echo "Marker count in decomposition:"
grep -c "BEGIN_PRESERVED_HISTORY" validation/reports/layer3_ivive_decomposition.md
```

Expected:
- Narrative section counts MATCH the pre-regen counts from Step 1
- Marker count = 1 in each file (still exactly one marker)

If counts differ: STOP, report BLOCKED. Investigate whether marker placement was wrong (e.g., narrative ended up in the data zone) or `_extract_preserved_history` has a bug.

- [ ] **Step 4: Verify the file structure (data zone above marker, history below)**

```bash
# Check that marker comes AFTER the orchestrator-managed sections
python3 -c "
from pathlib import Path
for fname in ['validation/reports/layer3_fih_dose.md', 'validation/reports/layer3_ivive_decomposition.md']:
    content = Path(fname).read_text()
    marker_idx = content.find('<!-- BEGIN_PRESERVED_HISTORY -->')
    notes_idx = content.find('## Notes')
    print(f'{fname}:')
    print(f'  Notes at offset: {notes_idx}')
    print(f'  Marker at offset: {marker_idx}')
    assert marker_idx > notes_idx, 'marker must come AFTER Notes section'
    # Verify no narrative section between Notes and marker
    between = content[notes_idx:marker_idx]
    has_narrative_in_data_zone = any(
        h in between for h in ['## Sprint', '## §']
    )
    assert not has_narrative_in_data_zone, 'narrative section ended up in data zone'
    print(f'  Structure OK')
"
```

Expected: both files print `Structure OK`, no assertion failures.

- [ ] **Step 5: Run full pytest suite**

```bash
pytest -q 2>&1 | tail -8
```

Expected: **945 passed** (941 Sprint 15 baseline + 4 new preservation tests).

If `test_point_mode_aafe_unchanged` fails (known flaky from prior sprints): retry once.

```bash
pip install -e ".[dev]" --quiet
pytest tests/unit/test_layer2_propagate_ci.py::TestPropagateCIMode::test_point_mode_aafe_unchanged -v
```

If still failing or any other test fails: STOP, report BLOCKED.

- [ ] **Step 6: Commit regenerated reports (if any data changes occurred)**

```bash
git status validation/reports/
git diff --stat validation/reports/
```

If `layer3_fih_dose.{md,json}` or `layer3_ivive_decomposition.{md,json}` show changes (e.g., timestamp updated or other minor data tweaks):

```bash
git add validation/reports/layer3_fih_dose.md \
        validation/reports/layer3_fih_dose.json \
        validation/reports/layer3_ivive_decomposition.md \
        validation/reports/layer3_ivive_decomposition.json
git commit -m "chore(sprint16): regenerate Layer 3 reports after marker insertion (verify preservation)"
```

If no changes (regen produced byte-identical output because marker preservation worked perfectly): skip the commit. Report this is the case.

**Do NOT commit:** `validation/reports/layer2_human_pk.{md,json}` drift (pre-existing unrelated drift; same exclusion as Sprints 12-15).

- [ ] **Step 7: Final branch state check**

```bash
git log --oneline main..HEAD
```

Expected: 2 or 3 commits depending on whether Step 6 produced a regen commit:
- 2 commits if regen was idempotent: `feat(sprint16) ...` + `chore(sprint16) ... markers`
- 3 commits if regen produced data updates: above + `chore(sprint16) ... regenerate Layer 3 reports`

```bash
git status
```

Expected: working tree clean except known `validation/reports/layer2_human_pk.*` drift.

---

## Self-Review Notes (for the implementer)

- **Marker placement matters.** In Task 2, the marker MUST go after `## Notes` (the last orchestrator-managed section) and BEFORE the first sprint narrative. If you place it inside the data zone (e.g., after `## Summary`), the next regen will overwrite the marker placement (because the orchestrator regenerates the full data zone). The marker effectively MUST live in the user-managed zone.

- **No code changes in benchmark scripts.** Tasks 1+2 are sufficient. The benchmark scripts (`layer3_fih_dose.py`, `layer3_ivive_decomposition.py`) call `emit_report()` unchanged — preservation is fully transparent. If you find yourself touching these files, STOP.

- **Honest verification in Task 3.** The Step 1 → Step 3 narrative section count comparison is the key proof that preservation works. If counts differ, the test isn't subjective — something is broken.

- **Idempotent regen is OK.** If Task 3 Step 6 produces no diff, that's a GOOD sign — it means the orchestrator-generated content was already current AND the marker preservation worked perfectly. No regen commit needed.

- **Backwards compat is the design contract.** Task 1's existing-test regression (Step 6) is critical. If `test_report_writer.py` tests fail, the implementation has changed default behavior (marker absent should mean overwrite, just like before).

- **Test-first discipline.** Task 1 uses TDD: write tests, verify they FAIL with the right error, implement, verify they PASS. This proves both that the tests are real (not no-ops) and that the implementation correctly addresses what was tested.

- **Commit count estimate:** 2 (if regen idempotent) or 3 (if regen produces data updates) feature commits + merge.

- **Estimated final test count:** 941 + 4 = **945**.

- **Worktree:** This plan should be executed in a worktree at `.worktrees/sprint16-orchestrator-append-only` on branch `feature/sprint16-orchestrator-append-only`. The using-git-worktrees skill creates this before Task 1.
