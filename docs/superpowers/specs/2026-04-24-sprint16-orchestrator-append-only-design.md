# Sprint 16 — Layer 3 report orchestrator append-only refactor

**Status:** Design approved (brainstorming with standing autonomy), pending user review of written spec before plan.

## 1. Context

Sprint 15 (merged `c0e4672`, 2026-04-24) encountered an architectural fragility flagged in its narrative + ticket: the Layer 3 report orchestrators (`validation/benchmarks/layer3_fih_dose.py`, `layer3_ivive_decomposition.py`) overwrite their entire `.md` output via the shared `validation/benchmarks/report_writer.py:emit_report()` function. This wipes any manually-appended sprint narrative sections every time a benchmark is regenerated.

Concrete losses observed during Sprint 15:
- `layer3_fih_dose.md`: Sprint 13 close-but-not-quite comparison table and Sprint 14 diazepam audit subsection were dropped in Sprint 15's regeneration (Task 4). Sprint 15 narrative + Sprint 14 reconciliation were appended fresh; prior history is now only retrievable via `git log -p` of the .md file.
- `layer3_ivive_decomposition.md`: Sprint 13 §9 narrative and Sprint 15 audit §10 were both wiped during Sprint 15 Task 4's regeneration. The Sprint 15 implementer manually restored §9 + §10 from `git show HEAD:` before appending §11 (per the Sprint 15 implementation report).

This manual-restoration workflow is error-prone (a future sprint could easily forget) and inconsistent (different .md files now have different historical depth). The fragility is present in `report_writer.py:284`:

```python
md_path.write_text(md_content, encoding="utf-8")
```

`emit_report()` is shared by all 4 layer benchmarks (`layer1_admet.py`, `layer2_human_pk.py`, `layer3_fih_dose.py`, `layer3_ivive_decomposition.py`), so a single fix benefits the entire validation suite.

## 2. Goal

Add **sentinel marker preservation** to `emit_report()`:

- Define constant `_HISTORY_MARKER = "<!-- BEGIN_PRESERVED_HISTORY -->"`.
- Before writing the new `.md` content, read the existing `.md` (if present), search for the marker, and capture everything from the marker through end-of-file.
- Append captured content to the newly-generated content, then write.
- If marker is absent (or file is absent), behavior is identical to current code (backwards compatible — Layer 1/2 reports unaffected).

Add the marker manually (once) to existing Layer 3 reports immediately before their first sprint narrative section. This converts the existing narratives from "orphan content" to "preserved history".

## 3. Non-goals

- Restoring lost Sprint 13/14 narrative content from `git log` (data-stale orphans; would create internal inconsistency with current Sprint 15 data tables — see Sprint 15 spec §3).
- Adding marker to Layer 1/2 reports (no narrative content currently — defer to future need).
- Sidecar history files (Approach B from brainstorming, rejected for separation overhead).
- Structured `extra_appendix` payload field (Approach C, rejected for parsing complexity).
- New API parameter to `emit_report()` (marker presence is the opt-in signal).
- JSON file preservation (`.json` is pure data — narratives only live in `.md`).
- END marker (`<!-- END_PRESERVED_HISTORY -->`) — YAGNI; preserved content always extends to EOF.
- Modifications to `validation/benchmarks/layer3_*.py` benchmark scripts (preservation is fully transparent at the `emit_report()` level).
- CLI flag to disable preservation (manual file deletion suffices for emergency reset).

## 4. Architecture

Single-function refactor. The contract:

```
[orchestrator-generated content][\n\n][marker][\n][user-managed history]
^                              ^         ^                              ^
|---- regenerated each run ----|         |---- preserved across runs --|
                                  ^
                              boundary
```

`emit_report()` flow with preservation:

1. Validate payload keys (existing).
2. Read existing `.md` file (if present); call `_extract_preserved_history()`.
3. Render new Markdown from payload (`_render_markdown()`).
4. If preserved content non-empty: concatenate `new_md.rstrip() + "\n\n" + preserved`.
5. Write `.md`. Write `.json` (unchanged).

`_extract_preserved_history(md_path: Path) -> str`:

```python
def _extract_preserved_history(md_path: Path) -> str:
    """Return content from history marker to EOF (inclusive of marker), or ''."""
    if not md_path.exists():
        return ""
    content = md_path.read_text(encoding="utf-8")
    idx = content.find(_HISTORY_MARKER)
    if idx == -1:
        return ""
    return content[idx:]
```

Edge cases handled correctly by this logic:
| Scenario | preserved | New file content |
|---|---|---|
| File doesn't exist | `""` | new only |
| File exists, no marker | `""` | new only (existing wiped — current behavior) |
| File exists, marker at EOF, no content after | `marker` | `new + marker` (idempotent) |
| File exists, marker mid-file with content after | `marker + history` | `new + marker + history` |
| File exists, multiple markers | first marker + everything after (incl. duplicate markers) | `new + (preserved that contains both markers)` |
| User removes marker between regens | `""` | new only (orphan history wiped — documented contract) |

## 5. Component design

### 5.1 Constant
```python
_HISTORY_MARKER = "<!-- BEGIN_PRESERVED_HISTORY -->"
```
HTML comment — invisible in rendered Markdown, visible in raw `.md`. Sufficiently unique that orchestrator-generated content (data tables, summary rows) won't contain it accidentally.

### 5.2 New helper
`_extract_preserved_history(md_path: Path) -> str` (see §4). Module-private.

### 5.3 Modified function
`emit_report(payload, *, stem) -> tuple[Path, Path]` — signature unchanged. Body adds marker preservation step before `md_path.write_text(...)`.

## 6. Files touched (exhaustive)

**Modified:**
- `validation/benchmarks/report_writer.py` — add constant, helper, integrate into `emit_report()` (~25 LOC)
- `validation/reports/layer3_fih_dose.md` — insert marker line before first narrative section (1 line)
- `validation/reports/layer3_ivive_decomposition.md` — insert marker line before §9 (Sprint 13) (1 line)

**New:**
- `tests/unit/test_report_writer_preservation.py` — 4 unit tests (~80 LOC)

**Unchanged (reaffirming):**
- All `src/charon/` code
- `validation/benchmarks/layer3_fih_dose.py`, `layer3_ivive_decomposition.py`
- `validation/benchmarks/layer1_admet.py`, `layer2_human_pk.py`
- All compound YAMLs, panel.yaml
- All other unit + integration tests

## 7. Marker placement convention

For Layer 3 reports, marker goes immediately after the last orchestrator-generated section (`## Notes`) and before the first manually-authored sprint narrative:

```markdown
## Notes

- Decomposition: fold_observed = fold_liver_model * fold_route_bias * fold_residual.
- ...

<!-- BEGIN_PRESERVED_HISTORY -->

## §9. Sprint 13 — UGT/CYP2C9 correction for diclofenac

After `hepatic_clint_multiplier: 3.5` added to diclofenac.yaml...
```

Future maintainers: append new sprint narratives below the marker. Do not add anything between the marker and `## Notes` — that content lives in the orchestrator's generated zone and will be wiped.

## 8. Testing plan

### 8.1 New unit tests in `tests/unit/test_report_writer_preservation.py`

1. **`test_emit_report_preserves_history_round_trip`**: 
   - Write report once via `emit_report` with minimal payload to a tmp path
   - Manually append `_HISTORY_MARKER + "\n\n## Test History\n\nUser content"` to the .md
   - Call `emit_report` again with a DIFFERENT payload
   - Assert: new payload's content is in the file, AND `## Test History` is preserved

2. **`test_emit_report_no_marker_means_no_preservation`**: 
   - Write report once
   - Manually append `## Stale\nContent` (without marker)
   - Call `emit_report` again
   - Assert: stale content is gone (current behavior preserved)

3. **`test_emit_report_marker_only_idempotent`**: 
   - Write report
   - Append `_HISTORY_MARKER + "\n"` (marker only, no following content)
   - Call `emit_report` twice with same payload
   - Assert: file content stable; marker present once at end

4. **`test_emit_report_marker_first_occurrence_wins`**: 
   - Write report
   - Append `_HISTORY_MARKER + "\n## A\n" + _HISTORY_MARKER + "\n## B\n"`
   - Call `emit_report` again
   - Assert: both `## A` and `## B` and both markers preserved (extracted from first marker to EOF)

### 8.2 Manual verification (post-Task 2)

After adding markers to Layer 3 reports:

```bash
python3 validation/benchmarks/layer3_fih_dose.py
python3 validation/benchmarks/layer3_ivive_decomposition.py
```

Expected:
- `[OK] Sanity floor: 12/12 pass.`
- `[OK] Decomposition wrote ... (12 compounds)`
- Both `.md` files preserve all `## Sprint N` and `## §N` sections after the marker

```bash
# Diff before/after regen — only the data section above marker should change
git diff validation/reports/layer3_fih_dose.md
git diff validation/reports/layer3_ivive_decomposition.md
```

Expected: changes confined to data section above marker; preserved history section below marker is byte-identical (or only differs in `Generated:` timestamp if there were data updates).

### 8.3 Regression

- Sprint 12/13/14/15 integration tests (6 total) — pass unchanged
- `test_fih_pipeline.py::test_tier_a_oral_pipeline_runs_to_mrsd` — pass
- Layer 1, 2 benchmarks — unaffected (no marker, no preservation)
- Existing report_writer tests (if any) — pass

### 8.4 Full suite

- Sprint 15 baseline: 941 tests
- Sprint 16: 941 + 4 new unit tests = **945 expected**

## 9. Success criteria

- `emit_report()` preserves content after `_HISTORY_MARKER` in existing `.md` files.
- 4 unit tests pass (round-trip, no-marker, marker-only, multi-marker).
- Layer 3 reports have marker added; regenerating either orchestrator no longer loses sprint narrative.
- Layer 1/2 reports unchanged (no markers, current behavior).
- 945 tests pass.
- No code in `src/charon/` modified.

## 10. Risks

| Risk | Mitigation |
|---|---|
| Marker accidentally appears in orchestrator-generated content | Marker is sufficiently unique (HTML comment with explicit name); document that orchestrator-generated content must not contain the literal marker |
| User places marker incorrectly (e.g., before `## Results`) → corrupts data section | Convention documented; user error visible immediately on next regen (file structure looks wrong) |
| Forgetting to add marker to a new report needing preservation | Default behavior unchanged (no marker = current); developer must intentionally add marker, which is the opt-in signal |
| `_extract_preserved_history` with very large file | `Path.read_text()` loads entire file; benchmark `.md` files are <100KB so memory is non-issue |
| User adds content BETWEEN `## Notes` and marker (in data zone) | Wiped on regen — documented contract; data zone is orchestrator-managed |
| Marker text changes in future (versioning) | Refactor required at that time; current marker should be stable for foreseeable future |
| Test file grows large | 4 tests are independent and small; ~80 LOC total |

## 11. Deliverables (minimum)

- 1 file modified (`report_writer.py`): constant + helper + integration (~25 LOC)
- 1 new test file (`test_report_writer_preservation.py`): 4 tests
- 2 markers added to Layer 3 `.md` files (1 line each)
- Verification: orchestrator regen confirms preservation
- 945 tests pass

## 12. Estimated test count

Sprint 15 ended at 941. Sprint 16 adds 4 unit tests. **Expected: 945.**

## 13. Sprint 17+ candidates (for context, not implementation)

After Sprint 16 closes the orchestrator fragility:
- Sprint 17: lisinopril 4.13x renal CL refinement (PEPT1 transporter modeling)
- Sprint 18: diclofenac 3.10x close-but-not-quite — investigate UGT-specific kinetics or extend literature multiplier search
- Sprint 19+: per-CYP2D6 ML retraining for propranolol 4.82x (architectural)
- Sprint 3b Sessions 2b/2c: formulation + precipitation (older carryover)

These are NOT Sprint 16 scope — Sprint 16 is purely infrastructure.
