# Sprint 7 — Conformal Pipeline Integration + Layer 3 FIH Dose Benchmark

**Date**: 2026-04-23
**Status**: Design (awaiting user approval)
**Owner**: Charon maintainers
**Supersedes**: none
**Related specs**: `2026-04-15-phase-a-finish-design.md` (predecessor)

---

## 1. Context & Motivation

Phase A finish (Sprint previous) delivered a "honest MVP" README and runnable
Layer 1/2 benchmark reports. A deeper audit before declaring Phase A complete
surfaced two load-bearing gaps:

### Gap A — Conformal prediction is built but unwired

CLAUDE.md §6k declares conformal prediction intervals a mandatory feature and
§6j makes the same declaration for CLint uncertainty. The code layer confirms
the intent:

- `src/charon/predict/conformal.py` (273 lines): `ConformalPredictor` with
  `calibrate()`, `calibrate_from_oof()`, `get_interval()`, physical bounds,
  coverage warning.
- `src/charon/predict/__init__.py:98-162`: `predict_properties(..., conformal=None)`
  explicitly attaches `ci_90_lower` / `ci_90_upper` to `fup` and
  `clint_hepatocyte` when a calibrated predictor is supplied.
- `src/charon/uncertainty/sampling.py:162-210`: Sobol/LHS sampling reads the
  attached CI widths to set per-parameter sigmas.

Production call sites reveal the gap:

- `cli/main.py:59`: `predict_properties(args.smiles)` — no conformal passed.
- `pipeline.py:104`: `predict_properties(smiles)` — no conformal passed.
- `validation/benchmarks/layer1_admet.py:171`: same.
- `scripts/train_clint.py`: computes OOF residuals in memory (line 147) but
  **never persists them**, so `calibrate_from_oof` has no input to work from.
- Grep for `ConformalPredictor` outside tests returns zero production
  call sites.

Net effect: every predicted `fup` and `clint_hepatocyte` ships with
`ci_90_lower = ci_90_upper = None`; Layer 4 uncertainty silently falls back
to hardcoded default sigmas; the generated report shows no honest Layer 1
confidence bands. The product's headline uncertainty claim is not backed by
running code.

### Gap B — Layer 3 FIH dose benchmark does not exist

`validation/benchmarks/layer3_fih_dose.py` was committed in the Phase A
finish sprint as an empty-body stub with a deferral rationale. No dataset,
no script, no report.

Layer 3 is the product's actual output — `recommend` subcommand returns an
FIH dose, and `ARCHITECTURE.md §8` sets the acceptance bar at "within 3-fold
for ≥60% of compounds". Without a benchmark, the product's headline claim
is unverified, and Phase A cannot honestly be called complete.

### Why these two together

Both gaps are validation-story gaps: "we said we do X, do we actually?" Both
sit at the boundary between trained models and user-visible output. Both
block a credible Phase A closeout. Tackling them in the same sprint keeps
the scope small enough to execute in one pass while producing a coherent
story: *after this sprint, every user-facing prediction carries a CI, and
the end-to-end dose claim is measured*.

### What this design excludes (and why)

- **Tier 3 CLint classification fallback** (§6j). Self-contained feature,
  independently valuable, non-blocking. Bundling inflates scope without
  compounding benefit. Defer to Sprint 8.
- **IVIVE bias correction** for highly bound / low-CLint compounds (warfarin,
  diazepam). Research-grade problem, belongs in a dedicated investigation
  sprint — not in an integration sprint.
- **BP-ratio conformal calibration**. The BP predictor is empirical
  (Hct-based formula), not ML, so conformal is conceptually the wrong
  tool — an empirical uncertainty budget is the right approach and belongs
  in a separate BP-characterisation task.

---

## 2. Goals

1. **Persist CLint OOF residuals at training time** so `calibrate_from_oof`
   has a reproducible on-disk input.
2. **Ship a `ConformalPredictor.load_default()`** factory that returns a
   calibrated predictor for `fup` and `clint_hepatocyte` with a JSON cache
   to avoid re-calibrating on every run.
3. **Auto-attach `ci_90_lower` / `ci_90_upper`** to `fup` and
   `clint_hepatocyte` in every path that runs `predict_properties` — CLI,
   `Pipeline.from_smiles`, benchmarks. Existing explicit `conformal=None`
   opt-out preserved for tests.
4. **Extend Layer 1 benchmark report** with empirical coverage and mean
   CI width columns, demonstrating that conformal actually ships.
5. **Extend Layer 2 benchmark report** to propagate CLint CI into per-compound
   CL CI (Monte-Carlo over the log-space CI), showing honest Layer 2 bands
   — without changing point-estimate AAFE computation.
6. **Implement Layer 3 FIH dose benchmark**: 5 gold FIH cases with
   FDA-briefing-document-sourced starting doses + full Obach-12 sanity
   floor check (`MRSD ≤ approved starting dose`).
7. **Commit Layer 3 report** (`.md` + `.json`) alongside Layer 1/2.

## 3. Non-Goals

- Not improving point-estimate AAFE for CL/Vss/t½.
- Not implementing Tier 3 classification fallback.
- Not adding BP-ratio or Peff conformal calibration.
- Not adding a new training dataset for fup/clint.
- Not adding Layer 3 to the Phase A "strict" gate — this benchmark reports
  only; a strict target can be enabled in a later sprint once baseline is known.

---

## 4. Architecture

### 4.1 Conformal integration — lifecycle

```
Training (one-shot, manual):
  scripts/train_clint.py
    └─ fits XGBoost
    └─ computes oof_log (scaffold 5-fold)
    └─ NEW: np.save("models/xgboost_clint_oof_residuals.npy",
                     abs(oof_log - y_true_log))
  scripts/train_fup.py
    └─ NEW: same persistence for fup (enables OOF-based fup calibration
            as a cross-check against the reference-CSV calibration)

First import (runtime):
  ConformalPredictor.load_default()
    1. Try load models/conformal_cache.json
    2. If missing OR source files newer than cache:
       - Instantiate ADMETPredictor
       - calibrate(predictor)  # fup via adme_reference.csv
       - calibrate_from_oof("clint_hepatocyte",
                            np.load("xgboost_clint_oof_residuals.npy"))
       - Serialize _reports dict (quantile_log10, factor, coverage, n)
         to conformal_cache.json
    3. Return calibrated predictor

Subsequent calls (within a process):
  Module-level singleton cached in predict/__init__.py.
  CLI / Pipeline / benchmark all get the same instance.

User-supplied SMILES:
  predict_properties(smiles)
    ├─ conformal = _default_conformal()   # lazy singleton
    ├─ adme = predictor.predict(smiles)
    ├─ fup_lo, fup_hi = conformal.get_interval("fup", adme.fup)
    ├─ clint_lo, clint_hi = conformal.get_interval("clint_hepatocyte",
                                                    adme.clint_hepatocyte)
    └─ return CompoundProperties with ci_90 attached
```

### 4.2 JSON cache format

```json
{
  "schema_version": 1,
  "generated_utc": "2026-04-23T12:34:56+00:00",
  "source_hashes": {
    "adme_reference.csv": "sha256:...",
    "xgboost_clint.json": "sha256:...",
    "xgboost_clint_oof_residuals.npy": "sha256:..."
  },
  "coverage_target": 0.90,
  "reports": {
    "fup": {
      "n_samples": 153,
      "quantile_log10": 0.412,
      "factor": 2.58,
      "empirical_coverage": 0.902,
      "median_fold_error": 1.74,
      "mean_fold_error": 1.91,
      "warning": null
    },
    "clint_hepatocyte": { ... }
  }
}
```

Invalidated on any source hash mismatch. No pickle — strict JSON for
cross-process / cross-platform determinism.

### 4.3 Layer 2 CI propagation (minimal, not a full re-design)

Layer 2 benchmark currently reports a point estimate per compound. New
behaviour: for each compound, draw N=100 log-space samples of
`clint_hepatocyte` (if ML-sourced — skip if experimental override) and
`fup` (same rule) from their conformal intervals, run the existing
simulate-then-extract pipeline, record CL as `[p05, p50, p95]`.

Key constraints:
- Only propagate CI on properties whose source is `ml_ensemble`. Experimental
  overrides retain point estimates — the user asserted a known value.
- Point-estimate AAFE (reported against `p50`) must numerically equal the
  current report to prevent silent regression. Assertion in test.
- Sampling is optional (`--propagate-ci` CLI flag on the benchmark);
  default report still shows point-only to keep CI runtime bounded.

### 4.4 Layer 3 FIH dose benchmark

Two-tier validation protocol:

**Tier A — gold FIH comparison (5 compounds):**

| Drug | Reference source | Published FIH start dose |
|---|---|---|
| Midazolam | NDA 20-942 | 1 mg IV |
| Warfarin | FDA label, initial dose | 2 mg PO |
| Propranolol | NDA 16-418 | 10 mg PO |
| Verapamil | NDA 18-817 | 40 mg PO |
| Omeprazole | NDA 19-810 | 20 mg PO |

**Reference-source caveat**: Each compound's published FIH dose is captured
in `panel.yaml` with an inline source citation. Where the primary FDA
briefing document is not publicly retrievable, the approved-label starting
dose (a conservative proxy — real FIH doses are typically at or below the
approved starting dose) is used with an explicit `source_type: label_start`
tag. The benchmark report surfaces this distinction so readers can discount
the fold-error accordingly.

Charon computes MRSD for each. Metric: fold-error = `max(pred/ref, ref/pred)`.
Acceptance: ≥ 3/5 within 3-fold, ≥ 4/5 within 10-fold.

**Tier B — sanity floor (Obach-12 full panel):**

For all 12 Obach compounds with a documented approved-label starting dose,
require `Charon_MRSD ≤ approved_starting_dose`. This is a one-sided sanity
check (a safe MRSD is always at or below the real clinical starting dose).
Acceptance: 12/12 pass.

Dataset lives in `validation/data/fih_reference/panel.yaml` with per-compound
sources cited.

---

## 5. Module-by-Module Changes

### Modified

| File | Change |
|---|---|
| `scripts/train_clint.py` | `np.save("models/xgboost_clint_oof_residuals.npy", residuals)` after OOF CV; update metadata JSON to record residual file hash. |
| `scripts/train_fup.py` | Same persistence for fup (used as cross-check, not default). |
| `src/charon/predict/conformal.py` | Add `load_default(cache_path=None)` classmethod; add `save_cache()` / `_load_cache()` helpers; add source-hash invalidation; add `clint_hepatocyte` to `_PROPERTY_BOUNDS` (already present — confirm). |
| `src/charon/predict/__init__.py` | Replace `conformal=None` default with lazy module-level `_default_conformal()`; respect explicit `None` as opt-out; update docstring. |
| `validation/benchmarks/layer1_admet.py` | Add coverage + mean CI width columns; report `ConformalPredictor.load_default()` used. |
| `validation/benchmarks/layer2_human_pk.py` | Add `--propagate-ci` flag; when set, N=100 MC over log-CI for ML-sourced `fup`/`clint`; output p05/p50/p95 for CL and Vss. |
| `src/charon/report/narrative.py` | If `fup.ci_90_lower`/`fup.ci_90_upper` present, render "90% CI [lo, hi]" in ADME table rows. |
| `validation/reports/layer1_admet.{md,json}` | Regenerated. |
| `validation/reports/layer2_human_pk.{md,json}` | Regenerated (point-only mode unchanged; new ci mode produces second report variant). |

### New

| File | Responsibility |
|---|---|
| `models/xgboost_clint_oof_residuals.npy` | Persisted scaffold-CV residuals. |
| `models/xgboost_fup_oof_residuals.npy` | Same for fup. |
| `models/conformal_cache.json` | Serialized CoverageReport dict. |
| `validation/benchmarks/layer3_fih_dose.py` | Runnable Layer 3 benchmark script. |
| `validation/data/fih_reference/panel.yaml` | 5 gold + 12 sanity-floor compounds with references. |
| `validation/reports/layer3_fih_dose.md` | Generated report. |
| `validation/reports/layer3_fih_dose.json` | Generated report. |
| `tests/unit/test_conformal_default.py` | Tests for `load_default`, cache round-trip, hash invalidation. |
| `tests/unit/test_predict_properties_conformal.py` | Confirms `predict_properties(smi)` attaches CI by default. |
| `tests/unit/test_layer3_fih_dose.py` | Unit tests for Layer 3 benchmark logic. |
| `tests/integration/test_conformal_pipeline.py` | End-to-end: SMILES → report carries Layer 1 CI. |

### Explicitly untouched

- `src/charon/core/schema.py` — `PredictedProperty.ci_90_*` already exists.
- `src/charon/uncertainty/sampling.py` — already reads `ci_90_*`, no change needed.
- `src/charon/pbpk/` — Layer 2 CI propagation happens in the *benchmark*, not
  in the core pipeline.

---

## 6. Layer 3 FIH Dose Benchmark — Detail

### 6.1 Input YAML schema

```yaml
panel:
  name: charon_sprint7_fih
  version: 1
  compounds:
    - name: midazolam
      smiles: "..."
      tier: gold
      reference_fih_mg: 1.0
      route: iv_bolus
      source: "NDA 20-942 briefing document, p.14"
      source_notes: "Phase 1 healthy volunteer study, 1 mg IV bolus"

    - name: propranolol
      smiles: "..."
      tier: gold
      reference_fih_mg: 10.0
      route: oral
      source: "NDA 16-418"

    - name: warfarin
      smiles: "..."
      tier: sanity_floor
      approved_starting_dose_mg: 2.0
      route: oral
      source: "FDA label, initial dose 2-5 mg"
      # Tier B pass criterion: Charon_MRSD_mg <= 2.0
```

### 6.2 Computation

```python
for compound in panel:
    props = load_compound_yaml(compound.name)  # reuses tier1_obach YAMLs
    pipe = Pipeline(props, route=compound.route, dose_mg=1.0, ...)
    result = pipe.run()
    mrsd_mg = result.dose_recommendation.mrsd_mg

    if compound.tier == "gold":
        fold = max(mrsd_mg / compound.reference_fih_mg,
                   compound.reference_fih_mg / mrsd_mg)
        row = {..., "fold_error": fold, "within_3x": fold <= 3.0,
               "within_10x": fold <= 10.0}
    else:  # sanity_floor
        pass_floor = mrsd_mg <= compound.approved_starting_dose_mg
        row = {..., "pass_floor": pass_floor}
```

### 6.3 Report structure

- Summary: Tier A pass rates (3x, 10x), Tier B pass rate.
- Gold table: compound, MRSD_pred, FIH_ref, fold, 3x pass, 10x pass.
- Sanity-floor table: compound, MRSD_pred, approved_start, pass.
- Notes: what was measured, what wasn't, data source list.

### 6.4 Acceptance criteria (reported, not gated)

| Criterion | Target | Action if missed |
|---|---|---|
| Gold Tier A within-3x | ≥ 3/5 | Report as is; open finding issue. |
| Gold Tier A within-10x | ≥ 4/5 | Report as is; open finding issue. |
| Sanity Tier B pass | 12/12 | **Fail CI** — MRSD exceeding clinical starting dose is a safety defect. |

---

## 7. Testing Strategy

### 7.1 Unit tests (new)

- `test_conformal_default.py`:
  - Load fresh (no cache) → produces cache file with correct schema.
  - Load second time → uses cache (no re-calibration) — asserted by
    comparing fup call counts via a spy.
  - Source file hash change → cache invalidated → re-calibration happens.
  - Cache corrupted → fall back to full calibration with warning.
- `test_predict_properties_conformal.py`:
  - Default call attaches CI on fup and clint.
  - Explicit `conformal=None` opt-out produces no CI.
  - Values respect physical bounds (fup upper ≤ 1.0).
- `test_layer3_fih_dose.py`:
  - Fold-error computation symmetric.
  - Sanity-floor passes for deliberately safe MRSD, fails for unsafe.
  - Report payload conforms to `report_writer` schema.

### 7.2 Integration tests (new)

- `test_conformal_pipeline.py`:
  - `charon report CCO --route iv_bolus --dose-mg 1` produces Markdown
    containing `"90% CI"` substring for fup and clint rows.
  - JSON report has numeric `ci_90_lower` / `ci_90_upper` on those
    properties.

### 7.3 Regression tests

- Existing `tests/regression/test_known_drugs.py` golden snapshots
  **will change** because reports now include CI. Update the 5 golden files
  in the same commit as the integration change, with the diff visible in
  the PR.
- Layer 1/2 AAFE values must remain numerically unchanged (point
  estimates untouched by conformal attachment). Assertion added to the
  benchmark test file.

### 7.4 Validation benchmark re-run

- Layer 1: coverage reported ≥ 85% (warning floor); sample counts stable.
- Layer 2 point-only mode: AAFE values bit-identical to pre-sprint report.
- Layer 2 ci mode: new variant report committed.
- Layer 3: baseline report committed — whatever the numbers are, they are
  the honest starting baseline.

---

## 8. Acceptance Criteria

Phase A closeout criteria met iff all of:

1. `charon report <smi> ...` output contains `90% CI` substrings on fup
   and clint rows.
2. `pytest tests/` passes, 790 → ~820 tests.
3. `validation/benchmarks/layer1_admet.py` run produces report with
   coverage column, committed.
4. `validation/benchmarks/layer2_human_pk.py --propagate-ci` run produces
   a second report variant, committed.
5. `validation/benchmarks/layer3_fih_dose.py` runs end-to-end, report
   committed.
6. Tier B sanity-floor: 12/12 pass. **Any failure blocks sprint completion**
   and requires explicit user approval to proceed (with a filed follow-up
   ticket). An unsafe MRSD is treated as a regulatory-class defect, not a
   best-effort reporting item.
7. Tier A within-3x: reported number captured in README "validation"
   section.
8. No new `# type: ignore` or `# noqa` added.
9. No change to CLAUDE.md MUST-NOT rules violated (audit during review).

---

## 9. Risks & Mitigations

| Risk | Likelihood | Mitigation |
|---|---|---|
| OOF residuals from existing trained models are lost (not persisted) and re-training perturbs point estimates | High | Retrain clint+fup as part of the sprint; commit new model + metadata + residuals atomically; compare new vs old AAFE. Drift <5% auto-accepted; drift ≥5% stops the sprint and surfaces to user — likely caused by non-determinism and a fixed seed may be needed before proceeding. |
| Module-level singleton fights with test isolation | Medium | Expose `reset_default_conformal()` for test teardown; `conformal=None` opt-out preserved. |
| Layer 2 CI propagation adds runtime (N=100 × 12 compounds × 72h BDF) | Medium | Keep point-only the default; gate behind `--propagate-ci`; document wall time. |
| FDA briefing doc numbers not easily verifiable | Medium | Cite each reference compound inline in YAML; a later curation pass can swap for more authoritative values without touching code. |
| Cache corruption in CI environments with stale models/ | Low | Strict hash check; re-calibration is idempotent; logs clearly say which path was taken. |
| Discovering Tier B sanity failures (unsafe MRSD) on some Obach compound | Medium | Honest result; fail the check, surface as a known IVIVE bias issue, file follow-up investigation. Better now than after Phase B. |

---

## 10. Out-of-Scope Follow-ups

Documented for handoff:

- **Sprint 8**: Tier 3 CLint classification fallback (§6j Tier 3).
- **Sprint 9**: IVIVE bias investigation for highly bound compounds.
- **Later**: BP-ratio empirical uncertainty budget; Peff prediction + CI.
- **Later**: Expand Layer 3 gold panel beyond 5 compounds.

---

## 11. Timeline Estimate

- **Day 1**: OOF residual persistence + retrain clint/fup + `load_default()` + JSON cache + tests.
- **Day 2**: `predict_properties` default wiring + CLI/pipeline integration + golden snapshot updates + integration tests.
- **Day 3**: Layer 1/2 benchmark updates + reports + Layer 3 dataset curation.
- **Day 4**: `layer3_fih_dose.py` implementation + tests + reports.
- **Day 5**: Verification, README update, PR.

Total: ~5 focused days (one working week). No net-new research, no new
ML training architecture, no new scientific claim — all integration work
on existing pieces.
