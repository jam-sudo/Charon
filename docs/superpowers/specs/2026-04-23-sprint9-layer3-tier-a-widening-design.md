# Sprint 9 — Layer 3 Tier A Panel Widening

**Date**: 2026-04-23
**Status**: Design (awaiting user approval)
**Owner**: Charon maintainers
**Supersedes**: none
**Related specs**:
- `2026-04-23-sprint7-conformal-integration-and-layer3-design.md` (introduced Tier A at n=5)

---

## 1. Context & Motivation

Sprint 7 delivered the Layer 3 FIH dose benchmark with a 5-compound Tier A gold panel. Result: 3/5 (60%) within 3-fold of published FIH / approved-start doses — **exactly on the CLAUDE.md §8 acceptance boundary** ("within 3-fold for ≥60% of compounds"). The Sprint 7 final review explicitly flagged this thin margin and recommended widening the panel to ≥10 compounds in a follow-up sprint.

Sprint 9 closes that follow-up. Expanded panel:

- **5 existing Tier A compounds** (unchanged results) — midazolam, warfarin, propranolol, verapamil, omeprazole.
- **4 Obach promotions** — theophylline, diclofenac, diazepam, metoprolol. These have FDA-label starting doses and existing `validation/data/tier1_obach/compounds/*.yaml` entries (Tier 1 experimental ADME data). Promoting them from Tier B sanity floor to Tier A gold is free data.
- **3 non-Obach additions** — acetaminophen, lisinopril, atorvastatin. Chosen to broaden elimination-pathway coverage (UGT/renal/CYP3A4-with-transporters) beyond the hepatic-CYP bias of the Obach-12 panel. Each requires a new compound YAML with primary-literature ADME values.

Total Tier A gold: **12 compounds**. Boundary resilience: a single fold-error flip now shifts the fraction by ~8 percentage points instead of 20, so the §8 gate is no longer single-point-of-failure.

### Why not ≥15?

Sprint 7 follow-up said "≥10". Curating non-Obach compound YAMLs is primary-literature work (Obach 2008 human PK predictions; FDA label / briefing docs; published therapeutic Cp). Three new compounds is a manageable batch; five or more inflates the data-validation burden without materially improving boundary resilience beyond n=12.

### What this design excludes (and why)

- **No oral route migration.** Compound YAMLs still lack Papp/Peff. Keeping `route: iv_bolus` per Sprint 7 — 1/F bias is documented and Tier A is report-only (Tier B gate is untouched).
- **No Tier B expansion.** The sanity floor is already 12 Obach compounds, 12/12 pass. Widening it is separate work.
- **No classifier changes.** Sprint 8's Tier 3 is orthogonal — Tier A compounds all have experimental CLint overrides, so Tier 3 never fires on them.
- **No benchmark-code changes.** `validation/benchmarks/layer3_fih_dose.py` is panel-driven — it reads `panel.yaml` and runs. No script logic changes.

---

## 2. Goals

1. Expand Tier A gold from 5 to 12 compounds in `validation/data/fih_reference/panel.yaml`.
2. Create 3 new compound YAMLs (`acetaminophen.yaml`, `lisinopril.yaml`, `atorvastatin.yaml`) in `validation/data/tier1_obach/compounds/` with primary-literature ADME values.
3. Tier A accuracy criterion (spec §8: "within 3-fold for ≥60% of compounds") must pass with comfortable margin (≥67% target; fail at <60%).
4. Schema test updates to reflect new minimum compound count.
5. README Layer 3 subsection updated with new n=12 numbers.
6. Regenerated `validation/reports/layer3_fih_dose.{md,json}` with the full 12-compound Tier A table.

## 3. Non-Goals

- Not changing Layer 3 benchmark script logic.
- Not switching any compound to oral route.
- Not rebalancing Tier B sanity floor.
- Not improving IVIVE accuracy for the out-of-tolerance compounds (that's research-grade, Sprint 9+).
- Not adding compounds whose FDA starting dose is unclear or whose ADME primary literature is incomplete.

---

## 4. Compound Dossier

### 4.1 Obach promotions (4 compounds)

Already have `validation/data/tier1_obach/compounds/*.yaml` with experimental ADME from Obach 1999 Table 2. `panel.yaml` currently lists them in Tier B only. Add Tier A gold entries:

| Compound | FIH reference | Route | `target_ceff_nM` | Therapeutic Cp rationale | Source |
|---|---|---|---|---|---|
| theophylline | 100 mg PO | iv_bolus | 55000 | Cp ~10 µg/mL, MW 180 | FDA label initial 100–200 mg |
| diclofenac | 50 mg PO | iv_bolus | 5000 | Cp ~1.5 µg/mL, MW 296 | FDA Voltaren label 50 mg PO TID |
| diazepam | 2 mg PO | iv_bolus | 1800 | Cp ~500 ng/mL, MW 285 | FDA Valium label 2–10 mg initial |
| metoprolol | 50 mg PO | iv_bolus | 380 | Cp ~100 ng/mL, MW 267 | FDA Lopressor label 50 mg PO BID |

### 4.2 New compounds (3 compounds)

Each new compound YAML carries Obach-style experimental ADME (source-cited) so Tier 1 overrides short-circuit any Tier 2/3 ML path.

**acetaminophen** (paracetamol)
- Route: iv_bolus; `reference_fih_mg`: 500; `target_ceff_nM`: 66000 (Cp ~10 µg/mL, MW 151).
- ADME: MW 151.16; logP 0.46; fu_p 0.80 (~20% protein bound); CLint_exp 1.8 µL/min/10^6 cells (Obach 2008 hepatocyte); bp_ratio 1.0; CLrenal 2.0 L/h.
- Source citations: FDA OTC monograph (325–650 mg q4-6h); Obach 2008 CL=21.6 L/h Vss=63 L; Prescott 1980 fu_p.
- Rationale: dominant elimination via UGT (sulfate/glucuronide) + minor CYP2E1. Broadens beyond CYP-focused panel.

**lisinopril**
- Route: iv_bolus; `reference_fih_mg`: 10; `target_ceff_nM`: 170 (Cp ~70 ng/mL, MW 405).
- ADME: MW 405.49; logP -1.22; fu_p 0.75 (~25% protein bound); CLint_exp 0 µL/min (not hepatically metabolised — `source: experimental`, explicitly zero); bp_ratio 0.85; CLrenal 5.0 L/h (~100% renal).
- Source citations: FDA Prinivil label 5–10 mg initial PO; Beermann 1988 human PK.
- Rationale: **renal-dominant elimination** (CLint=0 forces the pipeline to use CLrenal only). Stress-tests whether PAD path handles a renally-cleared drug correctly.

**atorvastatin**
- Route: iv_bolus; `reference_fih_mg`: 10; `target_ceff_nM`: 5.3 (Cp ~3 ng/mL, MW 558.6).
- ADME: MW 558.64; logP 6.36; fu_p 0.02 (~98% protein bound); CLint_exp 128 µL/min/mg HLM (note: **microsomal unit** — include `fu_inc: 0.07` for Austin correction per CLAUDE.md §6a); bp_ratio 0.61; CLrenal 0.0 L/h (<2% renal).
- Source citations: FDA Lipitor label 10 mg PO starting; Obach 2008 CL=39 L/h Vss=381 L; Lennernäs 2003 transporter review (OATP1B1 uptake — acknowledged limitation).
- Rationale: CYP3A4 + OATP1B1 uptake transporter. Charon doesn't model transporters explicitly, so this compound exercises the "Tier 1 override saves us" path — it's an honest stress test.

**Curation burden**: each new YAML is ~20 lines. Values cited inline with DOI-addressable references. No new math.

### 4.3 Expected fold-error outcomes

Literature-anchored predictions (informative, not binding):

| Compound | Prior Tier 2 CL_pred | Obach CL_obs | Fold | Within 3x? | Notes |
|---|---|---|---|---|---|
| midazolam | 5.77 | 21 | 3.6× | NO | existing; 1.73 was reported, discrepancy from PAD-path |
| warfarin | 0.78 | 0.19 | 4.1× | NO | existing; 1.04 reported — PAD gives different answer |
| propranolol | 4.72 | 50 | 10.6× | NO | existing; 36x reported (1/F dominant) |
| verapamil | 12.27 | 60 | 4.9× | NO | existing; 8x reported |
| omeprazole | 19.25 | 34.2 | 1.8× | YES | existing; 1.50 reported |
| theophylline | 3.26 | 2.9 | 1.1× | **YES** | Obach's easiest compound for IVIVE |
| diclofenac | 1.14 | 16 | 14× | NO | Obach's worst IVIVE (unusual CYP2C9/UGT) |
| diazepam | 0.026 | 1.6 | 62× | NO | very-low-fu_p disaster |
| metoprolol | 17.54 | 63 | 3.6× | borderline | may flip either way |
| acetaminophen | (new) | 21.6 | ? | likely YES | high fu_p, straightforward UGT |
| lisinopril | (new) | 5.1 | ? | likely YES | CLint=0 forces CLrenal path |
| atorvastatin | (new) | 39 | ? | ? | CYP3A4 + OATP — could overpredict |

**Boundary-hit projection**: 3-5 within 3x out of 12 = 25-42%. **This could FAIL §8.**

If the projection is correct, Sprint 9 will surface the honest truth that Charon's point-estimate IVIVE accuracy on this panel is below the §8 target. That is an acceptable outcome *as long as* Sprint 9 reports the number correctly. The alternative — cherry-picking only easy compounds — is worse.

**Contingency**: if actual within-3x ratio ≥ 60%, declare Sprint 9 success and close the Sprint 7 follow-up. If it's <60%, commit the report anyway (honest) and file a Sprint 10 ticket for IVIVE bias investigation. **Do not silently massage the panel to hit 60%**.

---

## 5. Module-by-Module Changes

### New files

| File | Responsibility |
|---|---|
| `validation/data/tier1_obach/compounds/acetaminophen.yaml` | Experimental ADME for acetaminophen. |
| `validation/data/tier1_obach/compounds/lisinopril.yaml` | Experimental ADME for lisinopril (CLint=0, renal-dominant). |
| `validation/data/tier1_obach/compounds/atorvastatin.yaml` | Experimental ADME for atorvastatin (high binding, CYP3A4, transporter caveat). |

### Modified files

| File | Change |
|---|---|
| `validation/data/fih_reference/panel.yaml` | +7 gold compounds (4 Obach promotions + 3 new). Tier B list untouched. |
| `tests/unit/test_fih_panel_schema.py` | `test_at_least_five_gold` → `test_at_least_twelve_gold`; `test_top_level_fields` threshold update from 15 → 22 compounds total. |
| `validation/reports/layer3_fih_dose.{md,json}` | Regenerated. Gold table grows to 12 rows. |
| `README.md` | Layer 3 subsection: update "(n=5)" → "(n=12)", update within-3x/within-10x counts. |

### Untouched

- `validation/benchmarks/layer3_fih_dose.py` — panel-driven, no code change.
- `validation/data/fih_reference/panel.yaml` Tier B section — unchanged.
- `models/*` — no retraining.
- Regression goldens — compound YAMLs added do not enter the regression pipeline (regression uses Pipeline.from_smiles, not panel-driven benchmarks).

---

## 6. Testing Strategy

### 6.1 Schema

`tests/unit/test_fih_panel_schema.py` updates:

- `test_top_level_fields`: total compounds ≥ 22 (was ≥ 15).
- `test_at_least_twelve_gold`: gold ≥ 12 (was ≥ 5).
- `test_gold_has_reference_fih` / `test_gold_has_target_ceff`: unchanged, will now iterate over 12 instead of 5.
- `test_each_compound_has_property_yaml`: new names (`acetaminophen`, `lisinopril`, `atorvastatin`) must resolve to existing files.

### 6.2 Compound YAML validation

Each new YAML must pass `charon.core.compound_config.load_compound_yaml` — the existing loader will reject malformed entries. Add one explicit smoke test per new compound:

```python
@pytest.mark.parametrize("name", ["acetaminophen", "lisinopril", "atorvastatin"])
def test_new_compound_loads(name):
    path = Path("validation/data/tier1_obach/compounds") / f"{name}.yaml"
    compound = load_compound_yaml(path)
    assert compound.name == name
    assert compound.molecular_weight > 0
    assert compound.properties.metabolism.clint_uL_min_mg is not None or \
           compound.properties.renal.clrenal_L_h is not None, \
           f"{name} must declare at least one clearance path"
```

### 6.3 Benchmark end-to-end

Re-run `python validation/benchmarks/layer3_fih_dose.py`. Expectations:

- Tier B: 12/12 pass (unchanged).
- Tier A: 12 rows. Report fold-error, within_3x, within_10x for each.
- Exit 0 unless Tier B breaks (it shouldn't).

### 6.4 Regression (advisory, not gating)

Run `pytest tests/regression/ -v`. No change expected — regression tests don't hit the FIH panel.

### 6.5 Acceptance criterion

Sprint 9 is successful iff:

- All 7 new / promoted compounds have legitimate FDA-sourced `reference_fih_mg` and Obach-style experimental ADME in their YAML.
- Panel schema tests updated and green.
- Tier B 12/12 pass.
- Tier A report committed with actual numbers.
- README updated to cite actual results.
- If Tier A within-3x ≥ 60%: §8 target PASS recorded; Sprint 7 follow-up closed.
- If Tier A within-3x < 60%: §8 target reported as FAIL, Sprint 10 follow-up ticket filed, commit proceeds. **Honest failure is acceptable; silent reshuffling to pass is not.**

---

## 7. Risks & Mitigations

| Risk | Likelihood | Mitigation |
|---|---|---|
| New-compound ADME values contain a transcription error (e.g., fu_p for wrong species) | Medium | Each YAML cites primary-literature DOI inline; spec reviewer checks vs Obach 2008 Table 2 directly. |
| Tier A within-3x drops below 60% after expansion | **High** | Honest reporting per §6.5. Sprint 10 follow-up for IVIVE bias correction. Not a code defect — a validation truth. |
| `fu_inc` required for atorvastatin's microsomal CLint entry makes the CYP3A4 pathway unusual | Medium | Document `fu_inc: 0.07` explicitly; test via `pytest tests/` that the compound YAML parses and runs through `parameter_bridge` without the Omega fu_p-double-apply bug (CLAUDE.md §0b-1). |
| atorvastatin OATP1B1 uptake not modelled → IVIVE overpredicts CL | High | Flag explicitly in the compound YAML `notes:` field. Acceptance criteria #5 is honest reporting — if atorvastatin fails within-10x, document; don't drop it. |
| lisinopril CLint=0 triggers divide-by-zero somewhere in the pipeline | Low | Existing Obach compounds don't test this, but the well-stirred liver model with CLint=0 returns CLh=0 (not an error) — verified in parameter_bridge tests. Still, add an assertion in the smoke test. |
| `pytest.approx` tolerance on the report test is too tight for the refreshed numbers | Low | Task-time fix: loosen regression snapshots if needed, but only for the compounds changed. |

---

## 8. Out-of-Scope Follow-ups

- **Sprint 10**: IVIVE bias investigation for highly-bound low-fu_p compounds (warfarin, diazepam) — research-grade.
- **Sprint 11**: Add Papp/Peff to compound YAMLs; switch Layer 3 panel to `oral` route for compounds whose clinical reference is oral.
- **Later**: OATP/MRP transporter plumbing (affects atorvastatin, statins class, many kinase inhibitors).
- **Later**: ULARGE-scale panel expansion (FDA briefing-doc ingestion at scale) — requires a sustained curation program.

---

## 9. Timeline Estimate

- **Day 1**: Curate 3 new compound YAMLs with primary-literature citations. Smoke-test they load via `load_compound_yaml`.
- **Day 2**: Update `panel.yaml` (7 new gold entries). Update schema tests. Run benchmark, inspect fold-errors.
- **Day 3**: Update README. Commit reports. If within-3x ≥ 60%, merge; if <60%, file Sprint 10 ticket and merge with honest numbers.

Total: ~3 focused days. No new code, no new math, no new ML. Pure data curation + panel refresh.
