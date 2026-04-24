# Sprint 14 — diazepam parameter audit (honest investigation)

**Status:** Design approved (brainstorming, standing autonomy), pending user spec review before plan.

## 1. Context

Sprint 13 (merged `8d4c812`, 2026-04-24) brought diclofenac fold from 10.23x to 3.10x — close but honest outcome. Tier A within-3x remained at 8/12 = 66.7% (§8 still PASSED).

Remaining residuals after Sprint 13:
- propranolol 28.55x — PBPK/ACAT architectural issue (oral F≈0.80 vs literature 0.26); not a multiplier fix, deferred to a larger investigation.
- **diazepam 4.91x** — Sprint 14 target.
- lisinopril 4.13x — renal-dominant, `hepatic_clint_multiplier` not applicable.

Diazepam is NOT systematically under-predicted by HLM IVIVE (Obach 1999 ratio ~1.5-2x, reference-drug status). Applying the Sprint 12/13 multiplier pattern without literature justification would be over-tuning. Sprint 14 instead audits diazepam's stored parameter values against primary literature and corrects any that are wrong.

## 2. Goal

Verify diazepam's `target_ceff_nM` (panel.yaml) and `clint_uL_min_mg`, `fu_p`, `bp_ratio` (compound YAML) against primary literature. Correct any value that disagrees with literature consensus by > 1.5x. Rerun benchmarks. Accept whatever honest outcome emerges — improvement, null result, or a literature-justified conservative multiplier.

## 3. Non-goals

- Applying an empirical multiplier to diazepam without primary-literature justification (would be over-tuning — diazepam is a reference-drug for IVIVE accuracy).
- Extended clearance / low-fu_p special-handling infrastructure (architectural; deferred).
- Changes to schema, parameter_bridge, pipeline, or any code file.
- Other compounds (propranolol / lisinopril are separate residuals).

## 4. Audit checklist

For diazepam, verify each of these values via WebSearch + literature cross-check:

### 4.1 `panel.yaml` diazepam entry
- `target_ceff_nM` — primary references:
  - **Greenblatt 1980** J Clin Pharmacol 20(11-12):616 — sedative/anxiolytic therapeutic range ~200-500 ng/mL (= 702-1755 nM at MW 284.74) for total Cp. Free Cp (at fu_p=0.013) ~10-23 nM.
  - **Mandelli 1978** Clin Pharmacokinet 3(1):72 — therapeutic window total 200-600 ng/mL.
  - Cross-check which (free or total) `target_ceff_nM` is intended to be; if mismatched units, that's the bug.

### 4.2 `diazepam.yaml` metabolism / binding
- `clint_uL_min_mg` — Obach 1999 DMD 27:1350 Table 2 lists diazepam HLM CLint.
- `fu_p` — literature range 0.013-0.016 across studies (Obach 1999, Benet 1996). Current value should be within this band.
- `bp_ratio` — diazepam BP ratio ~0.60-0.75. Current value should be within this band.

### 4.3 Expected discovery outcomes

- **Most likely:** `target_ceff_nM` is near the low end (~10-30 nM free, or ~300 nM total) — if current value matches one of the valid conventions (free vs total), no change. If it's obviously wrong (e.g., 1.0 nM), that's the bug.
- **Secondary:** `clint_uL_min_mg` matches Obach 1999; no change.
- **Tertiary:** `fu_p`/`bp_ratio` within literature band; no change.

## 5. Action protocol

For each audited value:

1. **Verified correct** (within 1.5x of primary literature) — no change.
2. **Off by >1.5x with primary-literature support for different value** — correct the value, cite the primary source, document the change in a YAML comment.
3. **Off but no clear primary reference** — flag but do not change; document uncertainty.

If NO values need correction:
- Sprint 14 reports honest null result: "diazepam residual 4.91x is irreducible at the current framework; closure requires extended-clearance model or low-fu_p special handling (Sprint 15+ scope)."

If ≥1 value corrected:
- Rerun benchmarks.
- New diazepam fold recorded.
- Other 11 compounds should show zero delta (only diazepam's YAML + its panel entry change).

If audit supports a conservative multiplier (i.e., primary literature explicitly reports ~2x underprediction for diazepam like Obach-documented values):
- Apply `hepatic_clint_multiplier: 2.0` with citation.
- Note explicitly in YAML and report: "conservative 2.0x at the lower end of Obach 1999 literature 1.5-2.5x range for CYP3A4/2C19 substrates; diazepam is a reference-drug so factor is applied conservatively."

## 6. Files touched (exhaustive, all branches)

**Conditional (depending on audit outcome):**
- `validation/data/fih_reference/panel.yaml` — if `target_ceff_nM` corrected
- `validation/data/tier1_obach/compounds/diazepam.yaml` — if CLint/fu_p/bp_ratio corrected OR if multiplier applied
- `validation/reports/layer3_fih_dose.{md,json}` — regenerated if any data changed
- `validation/reports/layer3_ivive_decomposition.{md,json}` — regenerated if any data changed

**Always:**
- `docs/superpowers/sprint10-ivive-bias-ticket.md` — Sprint 14 audit findings appended

**Unchanged:**
- All source code (schema, parameter_bridge, pipeline, etc.)
- All tests (no new tests — audit is literature-only)
- All other compound YAMLs

## 7. Testing plan

No new tests. Existing 938 tests should still pass regardless of outcome:
- If data corrected: tests pass (schema validation accepts corrected values)
- If null result: tests pass (no changes to testable code)
- If multiplier applied: diazepam MRSD changes but no integration test exists specifically for it (and Sprint 12/13 tests for atorvastatin/diclofenac continue passing)

If the implementer sees any test failure after their edits, that indicates a broken YAML — fix the YAML syntax.

## 8. Success criteria

Whichever is honest:

- **Best:** Data correction found → diazepam fold improves → within 3x possibly → 9/12 = 75%
- **Expected:** All values verify correct → honest null result documented → fold unchanged at ~4.9x, Tier A stays 8/12
- **Acceptable fallback:** Conservative literature-backed 2.0x multiplier applied → diazepam fold drops ~2x, likely within 3x → 9/12 if successful

## 9. Risks

| Risk | Mitigation |
|---|---|
| WebSearch inconclusive for diazepam literature values | Use Obach 1999 as authoritative source; if Obach value matches YAML, no change |
| `target_ceff_nM` convention (free vs total) ambiguous in panel.yaml | Investigate which convention is used across other compounds; use same convention |
| Audit finds no issue → "waste" of sprint effort | Honest null result IS a valid scientific outcome. Sprint 10 was also research-only. |
| Applying 2.0x multiplier is borderline over-tuning | Explicit caveat + only if literature (Obach 1999) directly supports ~2x for diazepam |

## 10. Deliverables summary (minimum)

- Sprint 14 audit findings document (could be a report file or ticket section)
- 0-3 YAML corrections with citations
- Possibly 1 multiplier application (with conservative caveat)
- Regenerated benchmarks (if data changed)
- Ticket update

## 11. Estimated test count

938 → 938 (no new tests).
