# Sprint 13 — UGT/CYP2C9 IVIVE correction for diclofenac

**Status:** Design approved (brainstorming with standing autonomy), pending user spec review before plan.

## 1. Context

Sprint 12 (merged `50c6c28`, 2026-04-24) closed atorvastatin's OATP1B1 residual by adding the generic `MetabolismProperties.hepatic_clint_multiplier` field and setting `value: 8.0` for atorvastatin only. That pushed Tier A within-3x from 7/12 (58.3%) to 8/12 (66.7%), achieving the §8 target ≥60%.

After Sprint 12, the remaining Tier A residuals (symmetric fold-error) are:
- **diclofenac: 10.23x** — CYP2C9 + UGT2B7 mixed clearance; HLM CLint systematically under-predicts UGT-mediated path
- propranolol: 28.55x — CYP2D6 + extensive first-pass (Sprint 14 target)
- diazepam: 4.91x — very low fu_p (0.013) well-stirred sensitivity
- lisinopril: 4.13x — non-hepatic elimination + low Peff

Diclofenac is the next achievable closure: its IVIVE gap is well-characterised (Miners 2006 Br J Clin Pharmacol; Rowland 2013 Drug Metab Rev) and the Sprint 12 template applies unchanged — same schema field, new YAML value, new literature citation, regenerate benchmarks.

## 2. Goal

Apply `hepatic_clint_multiplier` to diclofenac.yaml with a literature-sourced value for the UGT2B7 + CYP2C9 IVIVE gap. Rerun Layer 3 Tier A benchmark. Target: diclofenac fold ≤ 3x → Tier A within-3x 8/12 → 9/12 = 75%. Honest reporting per CLAUDE.md §6.5 if literature-supported multiplier does not fully close the gap.

## 3. Non-goals

- **Other UGT substrates.** Tier A panel has diclofenac as the only UGT-dominated compound.
- **`fm_UGT` schema field / per-enzyme phenotyping.** The current multiplier treats the full hepatic clearance path uniformly; per-enzyme splitting is not needed to close the single-compound gap. Future UGT scope would warrant separate treatment.
- **UGT-specific ML retraining.** Requires training data beyond what's available.
- **Schema / ParameterBridge / Pipeline changes.** Infrastructure was completed in Sprint 12; Sprint 13 is data-only + tests.
- **Any compound YAML beyond diclofenac.**
- **Panel.yaml changes.**

## 4. Literature basis for the multiplier

Published analyses of UGT-substrate IVIVE (via HLM with UDPGA activation):

- **Miners 2006 (Br J Clin Pharmacol 62:16):** systematic UGT2B7 substrate under-prediction averaging ~2.5-4x across the Obach panel.
- **Rowland 2013 (Drug Metab Rev 45:381):** compilation of in vivo/in vitro CLint_u ratios; diclofenac-specific cited values in the 3-4x range.
- **Obach 1999 (DMD 27:1350) Table 3 row-level ratios:** diclofenac in vivo/in vitro ~3-4x (classic under-prediction drug in the Obach panel).

**Planned value: 3.5** (midpoint of cited 3-4x range). Per Sprint 12's pattern, the implementer verifies via WebSearch before commit; if verified literature converges on a different value inside 2-5x, use verified value.

### 4.1 Expected outcome

- Current fold: 10.23 (diclofenac MRSD predicted ~4.89 mg vs reference 50 mg; 50/4.89 ≈ 10.2)
- With 3.5x multiplier: MRSD ~4.89 × 3.5 = ~17 mg (fold 50/17 ≈ 2.9, within 3x) — expected best-case
- With 3.0x multiplier: fold ~3.4 (just outside 3x)
- With 2.5x multiplier: fold ~4.1 (outside 3x)

Given the expected-case estimate of 2.9x sits right at the 3x boundary, outcome is **not guaranteed**. The honesty clause applies — do not inflate the multiplier beyond cited range to force a pass.

## 5. Architecture

No code changes. Uses infrastructure built in Sprint 12:

- Schema field `MetabolismProperties.hepatic_clint_multiplier: PredictedProperty | None = None` (already exists)
- `ParameterBridge.clint_to_clh(..., clint_multiplier=...)` kwarg (already exists)
- `ode_compiler.py:233` forwards `metab.hepatic_clint_multiplier.value` to ParameterBridge (already wired)

Sprint 13 adds data + test + report narrative only.

## 6. Files touched (exhaustive)

**Modified:**
- `validation/data/tier1_obach/compounds/diclofenac.yaml` — add `hepatic_clint_multiplier` block under `metabolism`
- `validation/reports/layer3_fih_dose.{md,json}` — regenerated
- `validation/reports/layer3_ivive_decomposition.{md,json}` — regenerated
- `docs/superpowers/sprint10-ivive-bias-ticket.md` — append Sprint 13 status

**New:**
- `tests/integration/test_diclofenac_ugt_enhancement.py` — mirrors `test_atorvastatin_oatp_enhancement.py` pattern; before/after MRSD comparison + clint_liver scaling

**Unchanged (reaffirming):**
- `src/charon/core/schema.py`
- `src/charon/core/parameter_bridge.py`
- `src/charon/pipeline.py`
- `src/charon/pbpk/ode_compiler.py`
- All other compound YAMLs
- `validation/data/fih_reference/panel.yaml`
- `tests/unit/test_decomposition.py`, `test_parameter_bridge.py`, etc.

## 7. Testing plan

### 7.1 New integration tests

`tests/integration/test_diclofenac_ugt_enhancement.py` — mirrors the Sprint 12 atorvastatin test structure:

1. `test_diclofenac_mrsd_with_multiplier_larger_than_without` — strip multiplier in-memory → compare MRSDs; assert enhanced > base.
2. `test_diclofenac_mrsd_ratio_within_literature_range` — assert ratio in [2.0, 5.0] (broader than atorvastatin's [4, 12] because diclofenac's baseline is lower-extraction so scaling is more linear but smaller).
3. `test_diclofenac_enhanced_clint_liver_metadata_scales` — assert `clint_liver_L_h` enhanced / base is exactly the multiplier value (within rtol 1e-3).

### 7.2 Regression

- Existing `test_atorvastatin_oatp_enhancement.py` still passes (atorvastatin unchanged).
- `test_fih_pipeline.py::test_tier_a_oral_pipeline_runs_to_mrsd` still 12/12 pass.
- `test_layer3_fih_dose.py` Tier B sanity floor still 12/12 pass.
- Sprint 10 decomposition log-additivity invariant still holds.

### 7.3 Full suite

- Expected count: 935 (Sprint 12 baseline) + 3 new integration tests = 938.

## 8. Success criteria

- **Primary target (stretch):** diclofenac fold ≤ 3.0x → Tier A within-3x = 9/12 = 75%.
- **Minimum acceptable:** diclofenac fold improves substantially from 10.23x (by at least the factor of the literature multiplier, well-stirred curvature permitting). If it lands in 3.1-4.0x range (still outside 3x), the sprint is reported as "meaningful improvement, §8 still PASSED at 8/12 = 66.7%, diclofenac approaches 3x boundary".
- **Regression:** zero change on the other 11 Tier A compounds.
- **No test failures.**

## 9. Risks

| Risk | Mitigation |
|---|---|
| Literature multiplier supports only 2-3x → diclofenac stays above 3x | Honest report; Sprint 13 is still progress. Do NOT tune multiplier beyond cited range. |
| Well-stirred curvature (if diclofenac is high-extraction): multiplier 3.5x yields < 3.5x CLh scaling | Expected effect small for diclofenac (fu_p=0.005, extreme low — low extraction, near-linear scaling expected). Verified by integration test range. |
| ConversionLog audit trail for diclofenac looks different from atorvastatin | Same `clint_enhancement` step mechanism; tested by integration test. No divergence expected. |
| diclofenac literature misinterpreted as pure UGT when CYP2C9 dominates → wrong mechanism framing | Method string explicitly cites both CYP2C9 + UGT2B7 pathways; empirical factor covers both. |

## 10. Deliverables summary

- 1 YAML edit (diclofenac.yaml) with Miners 2006 / Rowland 2013 / Obach 1999 citation
- 1 new integration test file (~3 tests)
- 2 regenerated report files with Sprint 13 narrative section (`## §9` in decomposition md; "Sprint 13 comparison" in fih_dose md)
- 1 ticket status section appended

## 11. Estimated test count

Sprint 12 ended at 935. Sprint 13 adds 3 integration tests. Expected end state: 938.
