# Sprint 11 — Tier A oral route migration (Papp/Peff curation)

**Status:** Design approved (brainstorming with standing autonomy), pending user spec review before plan.

## 1. Context

Sprint 9 (merged `ac61b72`) delivered the honest Layer 3 Tier A widening to n=12, reporting within-3x = 5/12 = 41.7% (§8 target ≥60% FAILED). Sprint 10 (merged `cb3ada0`) partitioned the fold-error via a three-factor signed log-additive decomposition. Aggregate attribution:

- liver_model choice: 4.9%
- **1/F route bias: 42.1%**
- residual: 137.2% (atorvastatin OATP-dominant)

The route-bias component is a *route-comparison artefact*, not an IVIVE error: the production pipeline simulates IV bolus while 11/12 Tier A clinical reference doses are oral-route. Sprint 11 removes this artefact by migrating the simulation path to oral. Expected upper-bound per-compound improvement is the full 1/F factor (e.g., propranolol 36.3x → ~9.4x, verapamil 8.1x → ~1.8x, atorvastatin 70.9x → ~10x).

ACAT + Fg was validated in Sprint 3b Session 2a (commit 2026-04-12; midazolam 1.01x, felodipine 1.26x, nifedipine 1.06x vs observed Fg), so the oral pipeline is production-ready. The only missing piece for Tier A oral MRSD is Papp/Peff curation on the 11 panel compounds that do not yet carry it.

## 2. Goal

Produce a Layer 3 Tier A benchmark run with simulation route = `oral`, predicated on curated literature Peff (or Papp with `papp_to_peff()` conversion), and report the honest within-3x outcome.

## 3. Non-goals

- **Tier B oral migration.** Tier B stays `iv_bolus`; its sanity-floor test already passes 12/12. Panel will carry mixed routes (Tier A oral, Tier B iv_bolus) — this is acceptable.
- **OATP1B1 modelling** (atorvastatin residual). Deferred to Sprint 12.
- **UGT / CYP2C9 calibration refresh** (diclofenac residual). Deferred to Sprint 13.
- **Solubility / formulation / precipitation**. ACAT oral works without explicit solubility (Sprint 3b-2a validated). Formulation is Sprint 3b-2b/2c scope.
- **Schema changes to `CompoundConfig`.** `permeability.peff_cm_s` and `permeability.papp_nm_s` fields already exist (`src/charon/core/schema.py:197-198`).
- **Sprint 10 decomposition library rewrite.** The orchestrator is updated minimally to track `simulation_route` per row; the pure-function library (`decomposition.py`) is untouched.

## 4. Success criteria

- 11 Tier A compound YAMLs carry a curated Peff (preferred) or Papp value with source-attributable metadata (`source`, `method`, citation).
- `panel.yaml` Tier A entries have `route: oral`; Tier B entries unchanged.
- Smoke: all 12 Tier A compounds complete `Pipeline(route='oral')` without exception or NaN.
- `layer3_fih_dose` benchmark produces a fresh within-3x number on the oral panel.
- Target: within-3x ≥ 6/12 = 50% (strong progress). Stretch: ≥ 7/12 = 58%, approaching §8 ≥60%.
- Decomposition rerun shows `aggregate_pct_route_bias` collapses to near 0% (sanity check that the migration actually eliminated the artefact).
- Regression: Tier B sanity-floor 12/12 still passes; Sprint 10 log-additivity invariant still holds.
- **Honesty clause (per CLAUDE.md §6.5):** if within-3x does not improve or improves below 50%, report the outcome verbatim. Do not revert, do not tune Peff values to hit the target.

## 5. Architecture

### 5.1 Per-compound YAML curation

Each of 11 Tier A compounds (midazolam already has it) gains a `properties.permeability` block:

```yaml
properties:
  permeability:
    peff_cm_s:
      value: <float>
      source: literature | compilation | ml_ensemble
      method: "<citation string>"
```

Or when only Papp is available:

```yaml
properties:
  permeability:
    papp_nm_s:
      value: <float>
      source: literature | compilation
      method: "<citation string>"
```

Fallback chain per compound:

1. **Obach 2008 Table 2** (8 panel compounds: warfarin, propranolol, verapamil, omeprazole, theophylline, diclofenac, diazepam, metoprolol) — well-curated Peff or Papp.
2. **Primary literature** for 3 Sprint 9 additions (acetaminophen, lisinopril, atorvastatin).
3. **Varma 2012 J Med Chem compilation** as secondary literature fallback.
4. **ML prediction via `predict_properties`** as last resort — flagged `source: ml_ensemble` so the report can surface which compounds lean on ML-predicted permeability.

Each value's citation string is auditable; no placeholder values.

### 5.2 Panel route flip

`validation/data/fih_reference/panel.yaml`:

- Tier A 12 entries: `route: iv_bolus` → `route: oral`
- `route_note` block: replace the existing IV-vs-oral caveat with a short note explaining that Tier A is now oral-route, Tier B remains iv_bolus for compatibility with existing sanity-floor tests.
- No other fields changed in panel.yaml. `target_ceff_nM`, `reference_fih_mg`, `tier`, etc., all preserved.

### 5.3 Benchmark re-run

`validation/benchmarks/layer3_fih_dose.py` does not change. It already reads `route` per-entry and runs `Pipeline(route=entry["route"], ...)`. Re-running it after the panel flip produces the new Tier A fold-error table. Commit the regenerated report alongside the compound YAML changes for a reviewable diff.

### 5.4 Decomposition orchestrator update

`validation/benchmarks/layer3_ivive_decomposition.py` needs a small addition: track `simulation_route` per row and use it to short-circuit route-bias to 1.0 when simulation and reference routes match.

Concretely, in `run_panel`:

```python
simulation_route = entry["route"]  # From panel.yaml Tier A entry
reference_route = bioav_row["fih_reference_route"]
if simulation_route == reference_route:
    # No 1/F artefact when simulating the same route as the clinical reference.
    fold_route, flags = 1.0, ()
else:
    fold_route, flags = compute_route_bias_factor(
        route_ref=reference_route, f_lit=bioav_row["f_oral"]
    )
```

Then feed `fold_route` + `flags` into the existing `decompose_fold_error` call by overriding its computed `fold_route_bias` post-hoc, OR by introducing a new thin wrapper that takes precomputed route_bias.

**Concrete choice:** Add an optional kwarg `route_bias_override: float | None = None` to `decompose_fold_error`. When provided, it bypasses `compute_route_bias_factor`. This keeps the library pure-function and testable; the orchestrator supplies the override.

Also add `simulation_route` and `reference_route` columns to each row dict so the report is legible.

### 5.5 Schema test update

`tests/unit/test_fih_panel_schema.py` adds two guards:

1. Every Tier A compound's YAML carries `properties.permeability.peff_cm_s` OR `properties.permeability.papp_nm_s`. If neither → test fails.
2. Every Tier A compound's panel entry has `route: oral`. If not → test fails.

### 5.6 Smoke test

`tests/integration/test_fih_pipeline.py` (existing; extend or add a new test):

Run `Pipeline(compound, route='oral', dose_mg=1.0, dose_projection=...)` for all 12 Tier A compounds. Assert: `result.dose_recommendation is not None` and `mrsd_mg > 0` for each. No NaN, no exceptions.

## 6. Data flow

```
Obach 2008 / primary lit / Varma 2012 / ML fallback
    │
    ▼
11 compound YAMLs gain permeability block
    │
    ▼
panel.yaml Tier A route flipped to oral
    │
    ▼
layer3_fih_dose.py rerun → new fold-error table
    │
    ▼
layer3_ivive_decomposition.py rerun (with simulation_route)
    │
    ▼
Report narrative comparison: Sprint 9 (iv) vs Sprint 11 (oral)
```

## 7. Files touched (exhaustive)

**Modified:**
- `validation/data/tier1_obach/compounds/acetaminophen.yaml`
- `validation/data/tier1_obach/compounds/atorvastatin.yaml`
- `validation/data/tier1_obach/compounds/diazepam.yaml`
- `validation/data/tier1_obach/compounds/diclofenac.yaml`
- `validation/data/tier1_obach/compounds/lisinopril.yaml`
- `validation/data/tier1_obach/compounds/metoprolol.yaml`
- `validation/data/tier1_obach/compounds/omeprazole.yaml`
- `validation/data/tier1_obach/compounds/propranolol.yaml`
- `validation/data/tier1_obach/compounds/theophylline.yaml`
- `validation/data/tier1_obach/compounds/verapamil.yaml`
- `validation/data/tier1_obach/compounds/warfarin.yaml`
- `validation/data/fih_reference/panel.yaml` (Tier A route + route_note)
- `validation/benchmarks/layer3_ivive_decomposition.py` (simulation_route tracking; route_bias override)
- `src/charon/translational/decomposition.py` (optional `route_bias_override` kwarg)
- `tests/unit/test_fih_panel_schema.py` (+2 guards)
- `tests/unit/test_decomposition.py` (+1-2 tests for override kwarg)
- `tests/integration/test_fih_pipeline.py` (+1 smoke test for Tier A oral)
- `validation/reports/layer3_fih_dose.{md,json}` (regenerated)
- `validation/reports/layer3_ivive_decomposition.{md,json}` (regenerated)

**Not modified:**
- `src/charon/pipeline.py`
- `src/charon/pbpk/acat.py`, `ode_compiler.py`, `solver.py`
- `src/charon/translational/pad.py`, `dose_projector.py`
- `src/charon/core/liver_models.py`, `schema.py`
- `validation/data/fih_reference/bioavailability.csv` (Sprint 10 artefact; no changes needed — `fih_reference_route` already captures the clinical reference route)

## 8. Risks

| Risk | Mitigation |
|---|---|
| ACAT numerical failure on low-Peff compound (lisinopril, F=0.25) | Smoke test catches this pre-benchmark. If fails: flag with `ml_permeability` source + explicit note in report; do not block release. |
| Peff literature disagrees across sources | Use median of cited range; notes column captures range. |
| within-3x does not improve or worsens | Report honestly per §6.5. Document which compounds regressed (if any) in report §4. Sprint 12/13 scope determined by residual pattern. |
| Atorvastatin OATP dominates even after 1/F removal (expected residual ~10x) | Document in report §3. Atorvastatin is a known Sprint 12 target — not a Sprint 11 failure. |
| Some compounds have substantially different Cmax vs. simulation (high Vd distributional lag) | PAD path uses steady-state approximation; this is inherent to the framework, noted in Sprint 3b-2a. |
| Tier B sanity floor breaks after orchestrator changes | Regression-tested by rerunning existing gated test. |

## 9. Testing plan

### 9.1 Unit tests (new / extended)

- `tests/unit/test_fih_panel_schema.py`:
  - `test_tier_a_compounds_have_permeability` — assert Peff or Papp in every Tier A YAML
  - `test_tier_a_panel_route_is_oral` — assert every Tier A panel entry has `route: oral`
- `tests/unit/test_decomposition.py`:
  - `test_decompose_with_route_bias_override_ignores_f_lit` — override kwarg bypasses `compute_route_bias_factor`

### 9.2 Integration tests (new)

- `tests/integration/test_fih_pipeline.py`:
  - `test_tier_a_oral_all_12_compounds_run` — parameterised smoke test over 12 Tier A compounds; assert `mrsd_mg > 0` and non-NaN

### 9.3 Regression

- Existing `test_layer3_fih_dose.py` continues to pass (Tier B gated)
- `test_layer3_decomposition_benchmark.py` — log-additivity still holds on the new (oral) simulation data
- Full suite: `pytest -q` green

## 10. Deliverables summary

- 11 modified compound YAMLs with curated permeability
- 1 modified `panel.yaml` (Tier A route flipped)
- 2 regenerated validation reports
- 1 updated orchestrator + 1 updated decomposition library (`route_bias_override`)
- ~4 new tests
- Sprint 11 status note appended to `docs/superpowers/sprint10-ivive-bias-ticket.md` (documents which residuals remain and which Sprint 12/13 tickets are justified)

## 11. Estimated test count

Sprint 10 ended at 909 tests. Sprint 11 adds ~4 tests. Expected end state: ~913 tests. All pass.
