# Sprint 15 — propranolol CYP2D6 IVIVE correction (audit + multiplier)

**Status:** Design approved (brainstorming with standing autonomy), pending user review of written spec before plan.

## 1. Context

The Sprint 10–14 cycle (2026-04-23 to 2026-04-24) brought §8 to PASS at 8/12 = 66.7% Tier A within-3x. Remaining residuals after Sprint 14:

- **propranolol 28.55x** — largest; Sprint 15 target
- diclofenac 3.10x — Sprint 13 close-but-not-quite (literature midpoint multiplier 3.5)
- diazepam 4.91x — Sprint 14 honest null (parameters verified, framework-limited)
- lisinopril 4.13x — renal-dominant, hepatic multiplier inappropriate

The Sprint 14 ticket characterized propranolol as: *"ACAT oral F computation architectural gap (ACAT likely gives F≈0.80 vs literature 0.26). Not a multiplier fix."* This characterization predates F-decomposition analysis. A direct analytical decomposition (re-derived in Sprint 15 brainstorming) suggests a different root cause:

| Component | Charon (analytical) | Notes |
|---|---|---|
| Fa | ~0.96 | Peff = 2.9e-4 cm/s → near-complete absorption |
| Fg | ~1.0 | Non-CYP3A4 substrate; gut metabolism negligible |
| Fh | ~0.93 | Qh / (Qh + fu_b · CLint_liver) = 99.45 / 107.74 |
| F_oral | ~0.89 | vs literature F = 0.26 → 3.4x overprediction |
| CLh_ws | 7.65 L/h (matches report) | vs literature in vivo CL ≈ 50 L/h → 6.5x underprediction |

The F-overprediction traces to CLint underprediction. Charon stores CLint_HLM = 13.6 µL/min/mg (Obach 1999 Table 2), but propranolol is a known CYP2D6 + CYP1A2 substrate with documented HLM IVIVE bias for high-extraction lipophilic bases (Ito & Houston 2005, Hallifax & Houston 2010).

**This is the same Sprint 12/13 pattern, not an ACAT architectural gap.** Sprint 15 verifies the decomposition empirically (Task 1), looks up the literature multiplier (Task 2), and applies it conditionally (Tasks 3–4). The Sprint 14 ticket characterization will be reconciled in Sprint 15 narrative.

## 2. Goal

Two-phase audit-then-act:

1. **Audit (Task 1):** Verify F-decomposition by capturing Charon's actual `Fa, Fg, Fh, F_oral` from `PKParameters` output for propranolol oral simulation. Run an empirical multiplier sweep (m ∈ {1, 2, 3, 5, 8, 12, 15, 20, 25, 30}) to produce a fold-vs-m curve.
2. **Literature (Task 2):** WebSearch + primary-literature cross-check for propranolol-specific or CYP2D6-class in vivo / in vitro CLint ratio. Choose conservative midpoint.

Conditional on audit + literature outcomes, apply `hepatic_clint_multiplier` to `propranolol.yaml` (Task 3) and regenerate benchmarks (Task 4). Accept whichever honest outcome emerges.

## 3. Non-goals

- **ACAT code changes.** Fa is expected normal; investigation only if Task 1 reveals Fa < 0.85.
- **Schema / parameter_bridge / pipeline / ode_compiler code.** Sprint 12 infrastructure reused (the field, kwarg, and call site already exist).
- **Other compounds.** Tier A panel touches one entry only.
- **CYP2D6 EM/PM polymorphism modeling.** Phase B scope.
- **Extended clearance model.** Sprint 16 candidate if Sprint 15 yields close-but-not-quite or null.
- **panel.yaml changes.** propranolol target_ceff_nM unchanged (audit may flag but not edit).

## 4. Audit checklist (Task 1)

For propranolol oral simulation (current `propranolol.yaml`, no multiplier), capture from `PKParameters`:

| Parameter | Expected (analytical) | Acceptable range | Action if outside |
|---|---|---|---|
| `fa` | ~0.96 | [0.85, 1.00] | If < 0.85 → ACAT issue; switch to Branch C, flag Sprint 16 |
| `fg` | ~1.0 | [0.90, 1.00] | If < 0.90 → unexpected gut metabolism; investigate |
| `fh` | ~0.93 | [0.85, 0.97] | If < 0.85 → CL not underpredicted; multiplier inappropriate |
| `bioavailability` (F_oral) | ~0.89 | [0.78, 0.96] | If < 0.78 → reconsider mechanism |
| `cl_apparent / fa·fg·fh` (CLh check) | 7.65 L/h | [6.5, 9.0] | Sanity match to decomposition table |

### 4.1 Multiplier sweep

For each m ∈ {1, 2, 3, 5, 8, 12, 15, 20, 25, 30}, override `hepatic_clint_multiplier` in-memory, run the FIH pipeline, record predicted MRSD and fold. Produce a table:

```
m  | MRSD_mg | fold_observed | within_3x?
1  |   0.3502 |     28.55     |   no
2  |   ...    |     ...       |   ...
...
30 |   ...    |     ...       |   ...
```

This curve is the primary empirical evidence. Identify `m_close_3x` (smallest m bringing fold ≤ 3x) and `m_close_1x` (m bringing fold closest to 1.0).

### 4.2 Output

Audit subsection appended to `validation/reports/layer3_ivive_decomposition.md` as `## §10 Sprint 15 audit (propranolol F-decomposition)`. Includes F-decomposition table + multiplier sweep curve. Markdown only; no JSON required. Output location matches Sprint 14 precedent (audit appended to existing report).

## 5. Literature multiplier search (Task 2)

WebSearch + primary-citation verification for propranolol IVIVE bias:

1. **Obach 1999 DMD 27:1350** — Table 3 (predicted vs observed CL_int back-calculation). Look up propranolol-specific row.
2. **Ito & Houston 2005 Pharm Res 22:103** — basic-drug HLM IVIVE survey; CYP2D6 / high-extraction subset average bias.
3. **Hallifax & Houston 2010 Drug Metab Pharmacokinet 25:74** — well-stirred extrapolation bias correction recommendations; includes ER-stratified bias factors.
4. **Riley 2005 Drug Metab Dispos 33:1304** — average HLM IVIVE bias factor (commonly cited 3x but ER-dependent).
5. **Tucker 2010 J Clin Pharmacol 50:1029** *(or earlier propranolol PK reviews)* — alternative source for in vivo CLint estimate.

### 5.1 Action protocol by literature outcome

| Literature support | Multiplier (apply) | Branch | Sprint 15 outcome |
|---|---|---|---|
| Cited propranolol-specific or class ratio ≥ 12x (≥2 sources) | midpoint of cited range (likely ~15) | A | Likely within-3x → §8 = 9/12 = 75% |
| Cited ratio 5–10x (≥2 sources) | midpoint (likely ~7) | B | Close-but-not-quite (Sprint 13 pattern); §8 = 8/12 |
| Cited ratio 3–5x | midpoint (likely ~4) | B | Substantially improved but well outside 3x; honest reporting |
| Cited ratio < 3x or absent / contradictory | NO multiplier | C | Null result; §8 = 8/12; Sprint 16 architectural flagged |

**Honesty discipline (CLAUDE.md §6.5):** do NOT exceed cited range to force §8 closure. If only Branch B is justified, accept it.

## 6. Architecture

No code changes. Sprint 15 reuses Sprint 12 infrastructure unchanged:

- `MetabolismProperties.hepatic_clint_multiplier: PredictedProperty | None = None` (schema, validated > 0)
- `ParameterBridge.clint_to_clh(..., clint_multiplier: float | None = None)` (kwarg)
- `ode_compiler.py:233` Pipeline forwarding (covers IV + oral; single call site)

Sprint 15 = **data + test + report narrative + ticket** only.

## 7. Files touched

**Conditional (Branches A or B — multiplier applied):**
- `validation/data/tier1_obach/compounds/propranolol.yaml` — add `hepatic_clint_multiplier` block under `metabolism` with literature citation method string
- `validation/reports/layer3_fih_dose.{md,json}` — regenerated
- `validation/reports/layer3_ivive_decomposition.{md,json}` — regenerated
- `tests/integration/test_propranolol_cyp2d6_enhancement.py` — new (3 tests)

**Always:**
- `docs/superpowers/sprint10-ivive-bias-ticket.md` — Sprint 15 section appended (with Sprint 14 reconciliation)
- Audit findings (Task 1 output) — written to `validation/reports/layer3_ivive_decomposition.md` Sprint 15 subsection or sibling file

**Unchanged (reaffirming):**
- All `src/charon/` code
- All other compound YAMLs
- `validation/data/fih_reference/panel.yaml`
- All other tests

## 8. Testing plan

### 8.1 New integration tests (Branches A or B only)

`tests/integration/test_propranolol_cyp2d6_enhancement.py` — mirrors Sprint 12/13 patterns:

1. **`test_propranolol_mrsd_with_multiplier_larger_than_without`** — strip multiplier in-memory → compare MRSDs; assert `enhanced > base`.
2. **`test_propranolol_mrsd_ratio_within_literature_range`** — assert `enhanced / base` within `[4, 25]`. The wider band reflects larger literature uncertainty for CYP2D6 vs OATP/UGT (Sprint 12 used [4, 12]; Sprint 13 used [2, 5]).
3. **`test_propranolol_enhanced_clint_liver_metadata_scales`** — assert `clint_liver_L_h` enhanced/base equals exactly the multiplier value (rtol 1e-3).

### 8.2 Regression

- `test_atorvastatin_oatp_enhancement.py` (3 tests) — pass unchanged
- `test_diclofenac_ugt_enhancement.py` (3 tests) — pass unchanged
- `test_fih_pipeline.py::test_tier_a_oral_pipeline_runs_to_mrsd` — 12/12 pass
- `test_layer3_fih_dose.py` Tier B sanity floor — 12/12 pass
- Sprint 10 decomposition log-additivity invariant — holds

### 8.3 Full suite

- Branch A or B: 938 → 941 (+3 integration tests)
- Branch C: 938 → 938 (no new tests)

## 9. Success criteria

- **Branch A (best):** propranolol fold ≤ 3x → §8 = 9/12 = 75%. Likely if literature documents ≥ 12x ratio.
- **Branch B (acceptable):** fold improves substantially (e.g., 28.55 → 4–6x) but stays outside 3x. §8 stays 8/12. Sprint 13 honest close-but-not-quite reporting.
- **Branch C (null):** literature absent / contradictory OR audit failure. No YAML change. §8 stays 8/12. Sprint 16 architectural investigation flagged.
- **Regression:** zero change on the other 11 Tier A compounds (multiplier scoped to propranolol).
- **Sprint 14 reconciliation:** Sprint 15 narrative explicitly states that Sprint 14 ticket's "not a multiplier fix" claim was unverified at the time and is corrected here.

## 10. Risks

| Risk | Mitigation |
|---|---|
| WebSearch yields no propranolol-specific ratio | Fall back to CYP2D6 class average (Hallifax-Houston 2010); honest range disclosure in YAML method string |
| Cited ratio supports only Branch B (close-but-not-quite) | Accept per §6.5; do NOT inflate above cited range; reuse Sprint 13 close-but-not-quite reporting template |
| Multiplier 15–20x feels "large" relative to Sprint 12 (8x) and Sprint 13 (3.5x) | Require ≥ 2 independent literature sources before applying; cite all in YAML method string; Sprint 15 narrative explains CYP2D6 IVIVE bias is documented to be larger than UGT or OATP1B1 systems |
| Audit Task 1 reveals Fa < 0.85 or Fg < 0.90 | Switch to Branch C; flag Sprint 16 ACAT investigation |
| Multiplier sweep contradicts back-calc estimate | Use empirical sweep (more reliable); reconcile with literature in narrative |
| CYP2D6 EM/PM polymorphism — F_lit = 0.26 is EM-only | Document in YAML; use median EM literature data (consistent with Wood 1978 and standard reviews) |
| Sprint 14 ticket reversal could appear inconsistent across reports | Be explicit: Sprint 14's claim was unverified at the time; Sprint 15 audit corrects it; both reports remain in repo for transparency |
| `PKParameters` Fa/Fg/Fh not directly exposed in current Pipeline result API | `PKParameters` from `pk_extract.py:215-227` includes `fa`, `fg`, `fh`, `bioavailability` fields. Pipeline result includes simulation, params, topology — sufficient to call `extract_oral_pk_parameters` directly in audit script |
| Test count regression / flaky test | None expected from this sprint; flaky `test_point_mode_aafe_unchanged` is pre-existing and unrelated |

## 11. Deliverables (minimum)

- 1 audit script run (Task 1) producing F-decomposition table + multiplier sweep markdown subsection
- 1 literature search summary (Task 2) with primary citations + chosen branch + chosen multiplier value (or no-multiplier decision)
- Conditional: 1 YAML edit (multiplier block), 1 integration test file, 2 regenerated reports
- 1 ticket Sprint 15 section with Sprint 14 reconciliation

## 12. Estimated test count

Sprint 14 ended at 938. Branches A/B → 941 (+3 integration tests). Branch C → 938 unchanged.

## 13. Predicted outcome (numerical, analytical)

For propranolol with multiplier `m`, using PAD ∝ CL_total / F_oral with renal CL = 0.1 L/h:

| m | CLint_liver_L_h | CLh_ws_L_h | Fh | F_oral | MRSD_mg | fold | within 3x |
|---|---|---|---|---|---|---|---|
| 1 (current) | 51.0 | 7.65 | 0.926 | 0.89 | 0.350 | 28.55 | no |
| 5 | 255 | 30.5 | 0.665 | 0.638 | 1.95 | 5.13 | no |
| 8 | 408 | 41.4 | 0.560 | 0.538 | 3.13 | 3.20 | no (just) |
| 12 | 612 | 50.4 | 0.470 | 0.451 | 4.59 | 2.18 | yes |
| 15 | 765 | 55.3 | 0.444 | 0.427 | 5.27 | 1.90 | yes |
| 18 | 918 | 59.6 | 0.385 | 0.370 | 6.26 | 1.60 | yes |
| 20 | 1020 | 62.1 | 0.375 | 0.360 | 6.95 | 1.44 | yes |

(Empirical sweep in Task 1 will replace this prediction.)

**Decision boundary:** m ≈ 9 brings fold to ~3x. Branches A (lit ≥ 12x) closes well within 3x. Branch B (lit 5–10x) stays outside 3x but improves dramatically.

---

## Appendix A — Sprint 14 ticket reconciliation language (draft)

To be inserted in Sprint 15 ticket section:

> Sprint 14 (diazepam audit) ticket characterized propranolol's residual as *"ACAT oral F computation architectural gap... Not a multiplier fix."* This characterization was a hypothesis based on a quick mental model, not on F-decomposition data. Sprint 15 audit (Task 1) explicitly captures `Fa, Fg, Fh` from the pipeline output and finds:
>
> - Fa ≈ [audit value] — ACAT functioning normally for high-permeability propranolol
> - Fg ≈ [audit value] — non-CYP3A4 substrate, no gut metabolism issue
> - Fh ≈ [audit value] — analytically traceable to CLint underprediction
>
> The gap is therefore the same Sprint 12/13 multiplier-template pattern (HLM CLint underprediction for IVIVE-known substrate classes), not ACAT architecture. Sprint 15 applies the multiplier template with literature citation per CLAUDE.md §6.5 honesty discipline.
