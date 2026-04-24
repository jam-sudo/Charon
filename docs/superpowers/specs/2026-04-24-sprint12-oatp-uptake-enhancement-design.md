# Sprint 12 — OATP1B1 hepatic uptake enhancement (minimal)

**Status:** Design approved (brainstorming with standing autonomy), pending user spec review before plan.

## 1. Context

Sprint 11 (merged `8ef3eb6`, 2026-04-24) migrated Tier A to oral route and closed the 42% route_bias aggregate Sprint 10 diagnosed. Outcome: within-3x = 7/12 = 58.3%, just 2% shy of §8 target ≥60%. The residual breakdown after oral migration:

- atorvastatin 4.64x (OATP1B1 unmodelled; largest per-compound mechanistic gap)
- propranolol 28.55x (CYP2D6 / F-correction partial)
- diclofenac 10.23x (UGT/CYP2C9 IVIVE underprediction)
- diazepam 4.91x (very low fu_p well-stirred sensitivity)
- lisinopril 4.13x (non-hepatic elimination)

Atorvastatin is the single most achievable residual to close via a mechanistic fix. Published extended-clearance analyses (Izumi 2018 Drug Metab Pharmacokinet 33(2):57–73; Barton 2013 Pharm Res 30:1077) consistently report that in vivo CLint for OATP1B1 substrates exceeds in vitro (hepatocyte or HLM) CLint by ~5–12x — the "IVIVE gap" for active-uptake substrates. Closing this gap for atorvastatin is expected to move its fold-error within 3x and cross §8 target 60%.

## 2. Goal

Add one optional field to the compound schema (`hepatic_clint_multiplier`) that enhances CLint_liver before the liver extraction model is applied, for the narrow case of OATP1B1-substrate compounds whose passive IVIVE under-predicts clearance. Curate the multiplier for atorvastatin from primary literature. Rerun Layer 3 Tier A benchmark and confirm §8 target crossing.

## 3. Non-goals

- **Full extended clearance model (PS_uptake + PS_metabolism separation).** Major PBPK ODE refactor; deferred to Sprint 15+ only if multiple OATP substrates enter future panels.
- **Other OATP substrates.** The current Tier A panel contains atorvastatin as the only documented OATP1B1 substrate. The field is generic (any compound may set it) but no other YAML is updated in this sprint.
- **Pgp substrate handling / efflux transporters.** Orthogonal mechanism, separate sprint.
- **Biliary clearance modeling.** Not currently represented in Charon PBPK.
- **Schema migration for `oatp_substrate: PredictedProperty | None`** field that exists but is unused. Leaving it untouched.
- **Changes to liver_models.py, PBPK ODE, uncertainty propagation.**

## 4. Architecture

### 4.1 Schema extension

Add one optional field to `MetabolismProperties` in `src/charon/core/schema.py`:

```python
class MetabolismProperties(BaseModel):
    primary_cyp: str | None = None
    secondary_cyp: str | None = None
    clint_uL_min_mg: PredictedProperty | None = None
    fm_cyp3a4: float | None = None
    hepatic_clint_multiplier: PredictedProperty | None = None  # NEW
```

Where `PredictedProperty` is the existing schema type (value + source + method). Default `None` preserves all current behaviour — no existing YAML needs modification.

Validator: `hepatic_clint_multiplier.value > 0` when present. A value of 1.0 is functionally a no-op but allowed (signals "explicitly no enhancement applied, not missing data").

### 4.2 ParameterBridge signature extension

`src/charon/core/parameter_bridge.py` — `ParameterBridge.clint_to_clh` gains one optional kwarg:

```python
def clint_to_clh(
    self,
    clint_uL_min_mg: float,
    ...,
    clint_multiplier: float | None = None,  # NEW
) -> HepaticClearance:
```

Behaviour:
- If `clint_multiplier is None` or exactly 1.0: no change from current behaviour.
- If `clint_multiplier > 0` and `!= 1.0`: after computing `clint_liver_L_h` (Step 3 of the IVIVE flow, per CLAUDE.md §6a), multiply by the factor and log the step:

```python
if clint_multiplier is not None and clint_multiplier != 1.0:
    if clint_multiplier <= 0:
        raise ValueError(f"clint_multiplier must be > 0, got {clint_multiplier}")
    enhanced = clint_liver_L_h * clint_multiplier
    steps.append(ConversionStep(
        name="clint_enhancement",
        value=enhanced,
        unit="L/h",
        formula=f"CLint_liver * multiplier = {clint_liver_L_h:.4g} * {clint_multiplier} = {enhanced:.4g}",
        notes=f"Empirical enhancement for uptake-limited substrates (e.g. OATP1B1). "
              f"Corrects passive IVIVE under-prediction."
    ))
    clint_liver_L_h = enhanced
```

The enhanced value then feeds into the liver extraction model exactly as an ordinary CLint_liver would. The enhancement is pharmacologically approximate (treats uptake as if it were additional metabolism) — rigour-level commentary is in Risk 9.3 below.

### 4.3 Pipeline wiring

`src/charon/pipeline.py` — the `_build_base_params` (or equivalent internal method) reads `compound.properties.metabolism.hepatic_clint_multiplier` and passes its `.value` to `ParameterBridge.clint_to_clh(clint_multiplier=...)`. If the field is `None`, pass `None` explicitly.

### 4.4 Atorvastatin YAML edit

`validation/data/tier1_obach/compounds/atorvastatin.yaml` — add to `properties.metabolism`:

```yaml
metabolism:
  ...
  hepatic_clint_multiplier:
    value: 8.0
    source: literature
    unit: ratio
    method: "Izumi 2018 Drug Metab Pharmacokinet 33(2):57 / Barton 2013 Pharm Res 30:1077 — in vivo/in vitro CLint_u ratio for OATP1B1 substrates 5-12x; 8 = midpoint consistent with atorvastatin published analyses"
```

Implementation-time verification: the implementer performs WebSearch against `"atorvastatin OATP1B1 IVIVE enhancement factor"` or `"atorvastatin in vivo in vitro CLint ratio"`. If literature median/median is outside 5–12x, use the verified value and cite.

### 4.5 ConversionLog capture

The `clint_enhancement` step, when present, appears in the `intermediate_steps` list of the `ConversionLog` returned by `clint_to_clh`. This preserves Charon's audit-trail requirement (CLAUDE.md §6g). Reports reading the ConversionLog will automatically surface the enhancement factor.

## 5. Testing plan

### 5.1 Unit tests (new)

`tests/unit/test_parameter_bridge.py` (or equivalent existing file):

1. `test_clint_multiplier_default_none_matches_current_behaviour` — with `clint_multiplier=None`, result is bit-for-bit identical to current behaviour for a representative input (midazolam-like values). Regression guard against accidental semantic drift.

2. `test_clint_multiplier_equal_one_is_noop` — `clint_multiplier=1.0` produces identical CLh to `None` case (within `rtol=1e-12`).

3. `test_clint_multiplier_scales_clint_liver` — `clint_multiplier=8.0` produces a CLint_liver 8x larger than the non-enhanced case (measured from ConversionLog `intermediate_steps`).

4. `test_clint_multiplier_raises_on_non_positive` — `clint_multiplier=0` or negative → ValueError.

5. `test_conversion_log_contains_enhancement_step` — when multiplier > 0 and != 1.0, intermediate_steps contains a `name="clint_enhancement"` entry with expected formula string.

### 5.2 Integration tests

`tests/integration/test_atorvastatin_oatp_enhancement.py` (new) or extended `test_fih_pipeline.py`:

1. `test_atorvastatin_oatp_enhancement_improves_mrsd` — compare atorvastatin MRSD with and without the multiplier (programmatically toggle by loading the YAML, zeroing the multiplier, and running Pipeline). Assert: with-multiplier MRSD is strictly larger than without-multiplier MRSD, and the ratio is within `[5, 12]` (matching the multiplier value of 8 ± logarithmic tolerance).

2. `test_tier_a_oral_smoke_includes_atorvastatin_enhanced` — the existing Sprint 11 smoke (`test_tier_a_oral_pipeline_runs_to_mrsd`) continues to pass. No regression on any compound.

### 5.3 Schema tests

Existing `test_compound_config.py` continues to pass (additive optional field on `MetabolismProperties`). If a test asserts schema completeness of compound YAMLs, that test should NOT require `hepatic_clint_multiplier` to be present (it's optional).

## 6. Success criteria

- **Atorvastatin fold-error within 3x.** Sprint 11: 4.64x. Sprint 12 target: ≤ 3.0x.
- **Tier A within-3x ≥ 8/12 = 66.7%** (§8 target ≥60% PASSED).
- **No regression on the other 11 Tier A compounds.** Each compound's fold-error in Sprint 12 matches Sprint 11's within ±5% numerical noise tolerance.
- **Tier B sanity floor 12/12** unchanged.
- All ~930 tests pass.
- **Honesty clause:** if atorvastatin's post-enhancement fold-error is worse than 3x or better than 1/3x (over-correction), report honestly per CLAUDE.md §6.5. Do not tune the multiplier to hit the target — use the literature-curated value and accept the outcome.

## 7. Files touched

**Modified:**
- `src/charon/core/schema.py` (+1 optional field on `MetabolismProperties`)
- `src/charon/core/parameter_bridge.py` (+1 kwarg, +1 ConversionStep emission)
- `src/charon/pipeline.py` (read + pass through)
- `validation/data/tier1_obach/compounds/atorvastatin.yaml` (+hepatic_clint_multiplier)
- Regenerated: `validation/reports/layer3_fih_dose.{md,json}`
- Regenerated: `validation/reports/layer3_ivive_decomposition.{md,json}`

**New:**
- `tests/integration/test_atorvastatin_oatp_enhancement.py` (or extension of existing integration test file)

**Unchanged:**
- `src/charon/core/liver_models.py`
- `src/charon/pbpk/*.py`
- `src/charon/translational/*.py`
- `validation/data/fih_reference/panel.yaml`
- Any compound YAML other than atorvastatin

## 8. Risks

| Risk | Mitigation |
|---|---|
| Enhancement factor varies widely across literature (5-12x) → atorvastatin over/under-shoots | Use median (8.0); WebSearch-verify during curation; document range in method string; accept honest outcome whichever way |
| "Multiplying CLint" conflates uptake with metabolism | Document in docstring + method string that this is empirical — correct numerical effect for uptake-limited substrates, not a claim of mechanistic accuracy |
| Appears compound-specific tuning since only atorvastatin uses it | Schema field is generic; multiplier has clear literature provenance; applies to *any* future OATP/NTCP substrate |
| Over-correction causes regression (atorvastatin now over-predicted >3x) | Success criteria explicitly allow honest reporting without retuning |
| Sprint 11's tests break due to schema change | Schema addition is strictly optional (`= None`); no existing YAML requires modification |
| Pipeline doesn't read the new field → silent bug | Integration test checks: with-multiplier MRSD > without-multiplier MRSD |

## 9. Alternate designs considered and rejected

### 9.1 Full extended clearance model (Option B from brainstorming)

Split liver compartment into blood and tissue sub-compartments, with separate PS_uptake and CL_metabolism parameters. Would require:
- New PBPK state variables (liver_tissue)
- PS_uptake parameter curation for every compound (or default from passive Peff × surface area)
- Stiff ODE re-validation against Sprint 3b-2a's 12-compound panel
- Schema changes broader than one field

**Rejected:** Scope vastly exceeds Sprint 12's goal (close atorvastatin's 4.64x residual and hit §8). Extended clearance is a future architectural sprint warranted only when multiple substrates demand the rigour.

### 9.2 Kp override (Option C)

Use existing `empirical_kp_by_tissue` mechanism to set liver Kp high for atorvastatin, simulating accumulation.

**Rejected — pharmacologically wrong:** Kp affects Vss and t1/2 but NOT systemic CL in the PAD MRSD formula. A high liver Kp would not change MRSD in a compartmental PBPK where the elimination rate depends on CLint and fu_b at the liver-blood interface. Would not fix the 5x underprediction.

### 9.3 Treating the multiplier as mechanistically accurate

Reframing the field as `oatp_uptake_cl_multiplier` or similar would suggest the multiplier represents a physically meaningful uptake clearance. It does not — it is an empirical IVIVE correction factor.

**Rejected:** The name `hepatic_clint_multiplier` is chosen for honesty. Docstring and method-string citation explicitly state it is empirical. Future extended-clearance work (Sprint 15+) may introduce `oatp_ps_uptake_ml_min_10_6_cells` or similar mechanistically-meaningful parameters and deprecate the multiplier.

## 10. Follow-ups (out of Sprint 12 scope)

- Sprint 13 — UGT / CYP2C9 calibration refresh (diclofenac 10.23x)
- Sprint 14 — CYP2D6 / beta-blocker F-correction path (propranolol 28.55x)
- Sprint 15+ — extended clearance model (if more OATP/NTCP substrates enter panel)
- Documentation — add OATP multiplier to the compound-YAML authoring guide if one exists

## 11. Estimated test count

Sprint 11 ended at 925. Sprint 12 adds 5–7 tests. Expected end state: ~930–932.
