# Sprint 17 — lisinopril parameter audit + framework-limit closure

**Status:** Design approved (brainstorming with standing autonomy), pending user review of written spec before plan.

## 1. Context

Sprint 16 closure flagged lisinopril 4.13x as `Sprint 17 candidate; renal CL refinement`. Subsequent context exploration during Sprint 17 brainstorming revealed that "renal CL refinement" is a misframing: lisinopril already has experimental `clrenal_L_h: 5.0 L/h` (Beermann 1988) used directly via the schema override path (`src/charon/pbpk/ode_compiler.py:252`), bypassing IVIVE entirely. Decomposition (`validation/reports/layer3_ivive_decomposition.md`) confirms:

- `fold_observed = 4.126`
- `fold_liver_model = 1` (CLint ≈ 0, hepatic CL negligible)
- `fold_route_bias = 1` (already on oral)
- `fold_residual = 4.126` (entire fold is unexplained by decomposition)
- `cl_renal_L_h = 5.0` (matches input directly)
- `F_pred = 0.2424` ≈ `F_obs = 0.25` (Beermann 1988 range 0.25-0.29)

Both PK observables (CL_total ≈ 5.32 L/h pred vs 5.1 L/h obs; F ≈ 0.24 pred vs 0.25 obs) are accurate. The 4.13x residual is therefore not a pharmacology prediction error — it must arise from PAD methodology assumptions (`target_ceff_nM = 170 nM`, `safety_factor = 10`) versus the FDA Prinivil label starting dose (10 mg PO OD). Sprint 10 §4 already audited `target_ceff_nM = 170` and found it plausible; no follow-up parameter audit has been done.

Sprint 17 is the formal audit-first close-out following the Sprint 14 (diazepam) honest-null pattern. The expected outcome is Branch C — declaring lisinopril 4.13x as a benchmark methodology gap (PAD with SF=10 vs FDA chronic dose label), not a pharmacology gap.

## 2. Goal

Verify every lisinopril parameter against primary literature; publish a per-parameter audit table; sensitivity-sweep the dominant levers (target_ceff_nM, safety_factor, Peff); reach a Branch A/B/C decision with empirical evidence; close the Sprint 16 lisinopril ticket with appropriate language.

This is a **research-only audit sprint**. No `src/charon/` code is touched (Sprint 14 precedent).

## 3. Non-goals

- PAD methodology refactor (`safety_factor`, target_ceff_nM redefinition for renal-dominant compounds) — too broad; would change benchmark scoring across all 12 Tier A compounds.
- PEPT1 / PEPT2 transporter infrastructure (intestinal / renal di-tripeptide transporters). Lisinopril `F_pred ≈ F_obs` already; explicit PEPT1 modeling cannot move the dial. Defer as future architectural sprint candidate (separate from this scope).
- OAT3 modeling for lisinopril renal handling — `cl_renal ≈ GFR × fu_p = 5.4 L/h` matches observed 5.0 L/h, indicating filtration-only behaviour with negligible net secretion/reabsorption. No transporter override would help.
- Reference-dose change in `validation/data/fih_reference/panel.yaml` (10 mg → 5 mg lower-bound) — that is benchmark methodology, not pharmacology.
- Updating `target_ceff_nM` to a Cavg-based redefinition — Sprint 10 §4 already verified 170 nM as plausible Cmax-based target; redefining cross-cuts the panel.
- §8 closure — Sprint 16 closed §8 PASSED at 8/12 = 66.7%. Sprint 17 expected outcome (Branch C) leaves §8 unchanged.

## 4. Architecture

Single-script audit + documentation update. No production code or tests modified.

```
scripts/sprint17_audit.py              [NEW]   F-decomposition + parameter audit + sensitivity sweeps
docs/superpowers/sprint10-ivive-bias-ticket.md  [MODIFIED]  Sprint 17 section + lisinopril resolution
validation/reports/layer3_ivive_decomposition.md  [MODIFIED]  §12 Sprint 17 audit narrative
validation/reports/layer3_fih_dose.md  [MODIFIED]  Sprint 17 closure narrative
```

The Layer 3 report modifications append below the `<!-- BEGIN_PRESERVED_HISTORY -->` markers introduced by Sprint 16; this is the first sprint to actually exercise that infrastructure for narrative additions. No regeneration of the Layer 3 reports is needed for Branch C (no parameter change → no benchmark drift).

For Branches A/B (parameter correction): YAML edit + integration test + regen via `validation/benchmarks/layer3_*.py`. The Sprint 16 marker mechanism preserves §9-§11 + new §12 across regen.

## 5. Audit script design

`scripts/sprint17_audit.py` (~150-200 LOC, idempotent stdout, no file writes):

### 5.1 Step 1 — Charon F-decomposition

- Load `validation/data/tier1_obach/compounds/lisinopril.yaml`
- Run `Pipeline(lisinopril, route="oral", dose_mg=1.0, dose_projection=DoseProjectionConfig(target_ceff_nM=170.0, safety_factor=10.0, tau_h=24.0))`
- Extract: `Fa, Fg, Fh, F_oral, clint_liver_L_h, cl_renal_L_h, CL_total, MRSD_PAD`
- Print "Decomposition" table with predicted vs analytical-expected vs literature-range columns

### 5.2 Step 2 — Parameter audit table

Per parameter: charon value | literature primary source | literature value/range | within range? | source PMID/DOI

| Parameter | Charon | Source (primary) | Literature range | OK? |
|---|---|---|---|---|
| `clrenal_L_h` | 5.0 | Beermann 1988 (PMID:3284731 or equivalent) | 4.5-5.4 (~95% of CL_total = 5.1) | ? |
| `fu_p` | 0.75 | Lancaster 1988 | 0.70-0.80 (low binding) | ? |
| `bp_ratio` | 0.85 | Beermann 1988 | 0.80-0.90 | ? |
| `clint_uL_min_mg` | 0.1 | Beermann 1988 (negligible hepatic) | < 1.0 | ? |
| `target_ceff_nM` | 170 | Beermann 1988 / panel notes | 150-200 nM (Cmax 70 ng/mL × 1000/MW=405) | ? (Sprint 10 §4 confirmed plausible) |
| `peff_cm_s` | 0.3e-4 | Knutter 2008 (cited; primary Peff NOT located) | LOW CONFIDENCE | flag |
| `F_oral_obs` | 0.25 (panel) | Beermann 1988 | 0.25-0.29 | ? |

Each row's "OK?" column populated empirically (✓ in range / ✗ out of range / ⚠ confidence flag).

### 5.3 Step 3 — Sensitivity sweeps

Three orthogonal sweeps to characterise residual structure:

**3a. target_ceff_nM** ∈ {85, 170, 340, 700} (0.5×, 1×, 2×, 4× current 170):
- For each: print MRSD_PAD, fold vs ref=10 mg
- Find m_close_3x and m_close_1x (smallest target_ceff_nM that brings fold within 3x / 1x)

**3b. safety_factor** ∈ {3, 5, 10}:
- For each: MRSD_PAD, fold vs ref=10 mg
- Document SF=3 vs SF=10 dose ratio

**3c. peff_cm_s** ∈ {0.1, 0.3, 0.5, 1.0, 2.0} × 1e-4:
- For each: F_oral predicted, MRSD_PAD, fold
- Identify whether passive Peff alone can explain the residual (expected: no — observed F=0.25 already matches predicted)

### 5.4 Step 4 — Branch decision rationale

Print decision tree outcome based on Steps 1-3 results:
- **Branch A** triggered if: audit finds a parameter outside literature range AND correction brings fold within 3x
- **Branch B** triggered if: parameter correction supportable but fold remains > 3x
- **Branch C** triggered if: all parameters within literature ranges AND no correction realistically closes residual; classify as "framework methodology gap" (PAD/SF=10 vs FDA label start dose)

Prior probability based on Sprint 10 §4 + decomposition: **Branch C dominant.**

## 6. Files touched (exhaustive)

### 6.1 New
- `scripts/sprint17_audit.py` — audit script (~150-200 LOC)

### 6.2 Modified — Branch C (expected)
- `docs/superpowers/sprint10-ivive-bias-ticket.md` — Sprint 17 section + lisinopril ticket resolution language
- `validation/reports/layer3_ivive_decomposition.md` — append §12 Sprint 17 audit narrative below sentinel marker
- `validation/reports/layer3_fih_dose.md` — append Sprint 17 closure narrative below sentinel marker

### 6.3 Modified — Branch A/B (contingent)
- `validation/data/tier1_obach/compounds/lisinopril.yaml` — parameter correction + citation block (Sprint 12/13/15 pattern)
- `tests/integration/test_lisinopril_<correction>.py` — NEW, mirrors Sprint 12/13/15 integration test pattern (3 tests, ratio band)
- `validation/reports/layer3_ivive_decomposition.{md,json}` — regen with corrected parameters
- `validation/reports/layer3_fih_dose.{md,json}` — regen with corrected parameters

### 6.4 Unchanged (reaffirming)
- `src/charon/` — entire production code untouched
- `validation/benchmarks/*.py` — orchestrators untouched
- `validation/benchmarks/report_writer.py` — Sprint 16 infrastructure exercised, not modified
- All other compound YAMLs

## 7. Branch decision protocol

Strict gate after Task 1 (audit run): **user reviews stdout output, confirms Branch (A/B/C), then Task 3 proceeds branch-specific.**

If a parameter is found outside literature range (Branch A trigger):
- Correction must use literature midpoint or anchored value (CLAUDE.md §6.5 honesty)
- Cited primary sources (≥2 if possible)
- Multiplier-template only if applicable (lisinopril has minimal hepatic CL; not the right pattern)
- More likely correction: direct YAML value update (e.g., fu_p, clrenal_L_h)

If close-but-not-quite (Branch B):
- Document partial improvement honestly; do NOT inflate to force §8 movement
- Sprint 13/15 precedent — accept residual outside 3x with explicit literature anchor

If honest null (Branch C — expected):
- No code change
- Three documentation locations updated (ticket, decomposition report §12, fih_dose report)
- "Framework methodology gap, not pharmacology gap" framing — analogous to Sprint 14 diazepam "framework-limited (low fu_p well-stirred sensitivity)" language

## 8. Testing plan

### 8.1 Branch C (no production code change)

- Existing 945 tests must pass unchanged (`pytest tests/` baseline)
- Sprint 17 audit script runs idempotently (manual verification: stdout stable across reruns)
- No new unit/integration tests added
- Test count: **945 (unchanged)**

### 8.2 Branch A/B (parameter correction)

- New integration test `tests/integration/test_lisinopril_<correction>.py` with 3 tests (Sprint 12/13/15 pattern):
  - `test_lisinopril_<param>_change_applied` — YAML override loaded
  - `test_lisinopril_<param>_dose_recommendation_changes` — MRSD shifts by expected ratio band
  - `test_lisinopril_<param>_within_literature_range` — corrected value sanity check
- 945 + 3 = **948** test count expected
- All Sprint 12/13/15 integration tests + Tier A regression remain pass (no cross-compound effect; YAML touch is isolated to lisinopril)

### 8.3 Regen verification (Branches A/B only)

- `python3 validation/benchmarks/layer3_fih_dose.py` → fold update reflects YAML correction
- `python3 validation/benchmarks/layer3_ivive_decomposition.py` → decomposition row update
- Both reports preserve §9-§11 narrative (Sprint 16 marker preservation verified — first usage)
- Manual diff inspection: only lisinopril row + Sprint 17 narrative changes

## 9. Success criteria

1. `scripts/sprint17_audit.py` exists, runs idempotently, prints clear decomposition + audit table + sensitivity sweep + Branch decision
2. Each parameter has a primary literature source explicitly cited (Beermann 1988, Lancaster 1988, etc.) with status (in-range / out-of-range / low-confidence-flagged)
3. Branch decision documented with empirical evidence (sensitivity sweep numerical results)
4. Sprint 16 lisinopril ticket entry updated with resolution language
5. Layer 3 reports preserve all prior narrative (Sprint 16 marker proven in production)
6. Test suite remains green (945 for Branch C, 948 for Branch A/B)
7. §8 status updated correctly (Branch C: 8/12 unchanged; Branch A: potentially 9/12 = 75.0%)
8. CLAUDE.md §6.5 honesty discipline preserved (no inflated multipliers)

## 10. Risks

| Risk | Mitigation |
|---|---|
| Audit finds parameter mid-correction zone (Branch B partial) | Sprint 13/15 precedent: accept honest close-but-not-quite, document anchor |
| Audit script fails to import Pipeline | Sprint 14 precedent: confirm imports first; structure mirrors `sprint15_audit.py` |
| Regen wipes Sprint 16 markers (regression) | Sprint 16 marker tests still in suite; if wipe occurs it's a Sprint 16 regression — fix at that level, not here |
| User overrides Branch C and pushes for Branch A under unverified parameter correction | CLAUDE.md §6.5 must be cited; Sprint 14 precedent shows null is acceptable |
| Peff "LOW CONFIDENCE" finding triggers desire to "fix" Peff with no primary citation | Branch boundaries strict — if no primary lit, no correction; document low-confidence flag instead |
| Branch A/B ratio band too narrow → integration test flaky | Use Sprint 12/13/15 ratio band convention ([N×0.6, N×1.5] for fold ratio) |

## 11. Deliverables (minimum, Branch C)

- 1 new file (`scripts/sprint17_audit.py`)
- 3 modified docs (ticket, decomposition report §12, fih_dose report)
- 0 modified production code
- 0 modified tests
- 945 tests pass

## 12. Estimated test count

**Branch C (expected):** 945 (Sprint 16 baseline unchanged).
**Branch A/B (contingent):** 948 (3 integration tests added).

## 13. Sprint 18+ candidates (for context, not implementation)

After Sprint 17 closes the lisinopril ticket:
- Sprint 18: propranolol 4.82x — per-CYP2D6 ML retraining (architectural)
- Sprint 3b 2b/2c: oral formulation + precipitation kinetics (older carryover)
- Future: PEPT1/PEPT2 transporter infrastructure as architectural sprint (multi-compound benefit; not justified for lisinopril alone)

These are NOT Sprint 17 scope.
