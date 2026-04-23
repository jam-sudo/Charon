# Sprint 8 — CLint Tier 3 Classification Fallback

**Date**: 2026-04-23
**Status**: Design (awaiting user approval)
**Owner**: Charon maintainers
**Supersedes**: none
**Related specs**:
- `2026-04-23-sprint7-conformal-integration-and-layer3-design.md` (predecessor; Tier 2 + conformal integration)

---

## 1. Context & Motivation

CLAUDE.md §6j declares a three-tier strategy for CLint:

- **Tier 1** — experimental override (user-supplied). Implemented.
- **Tier 2** — XGBoost regression with conformal CI. Implemented in Sprint 2 / integrated in Sprint 7.
- **Tier 3** — classification fallback (Low/Med/High buckets) for compounds outside the applicability domain. **Not yet implemented.**

The spec justifies Tier 3 explicitly: "SMILES-only CLint prediction cannot reliably achieve AAFE < 2.0, so for structurally novel scaffolds a probabilistic bucket label is safer than a precise-looking point estimate with a huge interval." Without Tier 3, novel compounds silently get a Tier 2 point estimate whose conformal CI is honest but wide (factor 6.74×), and users have no clear signal to obtain experimental CLint.

Sprint 8 closes this gap.

### Empirical verification (done during brainstorming)

1. **Bucket balance** — Spec buckets [<10, 10-50, >50] µL/min/10^6 cells divide the training set into 33% / 38% / 29% (1495 samples). Quantile 0.33 = 10.0, quantile 0.67 = 40.7 → spec boundaries are a natural tercile split. No rebalancing needed.

2. **Classifier feasibility** — Scaffold 5-fold CV of XGBoost multi-class (3 classes, 200 trees, depth 5, `random_state=42`) on 1476 compounds: **53.6% accuracy, Macro F1 = 0.539**. Per-class F1 {Low 0.58, Med 0.50, High 0.54}. Critically, 2-off confusions (Low↔High) are rare: 33/480 (6.9%) and 57/428 (13.3%). Most confusions are adjacent-bucket — the safety-relevant property.

3. **Layer 0 AD is dead code** — `GuardrailChecker` in `src/charon/core/guardrails.py` is defined but never instantiated in production. `_reference_fps` defaults to `None`, which returns `"UNKNOWN"` always. Wiring Layer 0 AD site-wide is a separate infra project. Sprint 8 ships a **CLint-local AD** instead: max Tanimoto similarity against the 1476-compound CLint training set.

4. **AD=LOW prevalence** — Under the 0.3 threshold inherited from `guardrails.py`: 4/14 Obach compounds (dextromethorphan, omeprazole, theophylline, verapamil) and 50% of the 151-compound adme_reference set fire as LOW. This is a **feature, not a bug**: diversity-curated reference panels are expected to span applicability domain. Tier 3 firing on half of novel scaffolds is the intended behavior; the README must document this clearly.

---

## 2. Goals

1. **Ship a CLint 3-class classifier** (`models/xgboost_clint_classifier.json`) with scaffold-CV Macro F1 ≥ 0.50 as the hard gate.
2. **Add CLint-local AD** (`predict/clint_ad.py`) decoupled from Layer 0.
3. **Auto-trigger Tier 3** when AD_max < 0.3 or when the caller passes `force_tier3=True` to `predict_properties`.
4. **Propagate Tier 3 output** as a `PredictedProperty` with `source="classification"`, a CRITICAL flag, a bucket-range CI, and a new `classifier_probs` field carrying `{low, med, high}` probabilities.
5. **Extend Layer 4 uncertainty** to sample categorical × log-uniform-within-bucket for classification-sourced CLint.
6. **Show CRITICAL warnings** in the Markdown/JSON report when any property is classification-sourced.
7. **Document AD-LOW prevalence** in README + report narrative ("expected on ~50% of diverse compounds; experimental CLint resolves it").

## 3. Non-Goals

- **Not wiring Layer 0 `GuardrailChecker` site-wide.** Sprint 8 ships CLint-local AD only; general Layer 0 AD is a separate effort.
- **Not adding per-compound local uncertainty for Tier 2.** Conformal gives a global quantile; per-compound confidence estimation is research-grade (e.g. tree-leaf statistics, ensemble disagreement) and out of scope.
- **Not rebalancing bucket boundaries.** [<10, 10-50, >50] is balanced and matches CLAUDE.md §6j.
- **Not adding classification fallback to `fu_p` or other Layer 1 properties.** CLint is the MUST-NOT-recommend-without-data property per §6j; others are Tier 2 only.
- **Not adding new training data.** Same 1495 compounds, same features.

---

## 4. Architecture

### 4.1 Tier 3 lifecycle

```
Training (one-shot, manual):
  scripts/train_clint_classifier.py
    - load data/training/clint_merged.csv (tdc_hep + chembl)
    - exclude validation InChIKey-14 set
    - bucket clint_hep into {Low=0, Med=1, High=2}
    - scaffold 5-fold CV → compute Macro F1 + confusion matrix
    - fail if Macro F1 < 0.50
    - fit final model on full data
    - save models/xgboost_clint_classifier.json
    - append entry to models/model_metadata.json
    - persist OOF class probabilities for reliability diagrams

First predict_properties call:
  clint_ad.load_default() returns a ClintLocalAD instance
    - loads the same 1476-compound Morgan FP set cached to disk
    - cache file: models/clint_ad_fingerprints.npz

predict_properties(smiles, force_tier3=False):
  1. Run Tier 2 (existing XGBoost regression) + conformal CIs.
  2. Compute AD_max = max Tanimoto vs CLint training FPs.
  3. Decision:
     - AD_max >= 0.5 (HIGH): Tier 2 only, flag="ad_high"
     - AD_max >= 0.3 (MODERATE): Tier 2 + flag="ad_moderate_caution"
     - AD_max <  0.3 (LOW)  OR  force_tier3: Tier 3
  4. If Tier 3:
     - ClassifierPredictor.predict_proba(smiles) → {low, med, high}
     - point_estimate = bucket_center(argmax_bucket)
     - ci_90 = (bucket_lower, bucket_upper) of the argmax bucket
     - source = "classification"
     - flag = "clint_tier3_classification; CRITICAL — experimental measurement essential; AD_max={x:.3f}"
     - classifier_probs = {"low": p_low, "med": p_med, "high": p_high}
     - Attach to CompoundProperties.metabolism.clint_uL_min_mg
```

### 4.2 Bucket definitions (frozen)

```python
CLINT_BUCKETS = {
    "low":  {"range": (0.1, 10.0),   "center": 3.0,    "label": 0},
    "med":  {"range": (10.0, 50.0),  "center": 22.0,   "label": 1},
    "high": {"range": (50.0, 1000.0),"center": 200.0,  "label": 2},
}
```

Centers are log-space geometric midpoints: `exp((ln lo + ln hi) / 2)`. Rationale: CLint is log-distributed; arithmetic midpoint skews high.

### 4.3 New/modified types

```python
# charon/core/schema.py — add "classification" to SourceType literal
SourceType = Literal[
    "ml_ensemble", "ml_pka", "correlation", "derived",
    "physiological", "experimental", "literature",
    "classification",   # NEW
]

# charon/core/schema.py — add classifier_probs field to PredictedProperty
class PredictedProperty(BaseModel):
    value: float
    ci_90_lower: float | None = None
    ci_90_upper: float | None = None
    source: SourceType
    unit: str | None = None
    method: str | None = None
    flag: str | None = None
    classifier_probs: dict[str, float] | None = None   # NEW

    @model_validator(mode="after")
    def _probs_sum_to_one_if_present(self):
        if self.classifier_probs is not None:
            s = sum(self.classifier_probs.values())
            if not (0.99 <= s <= 1.01):
                raise ValueError(f"classifier_probs must sum to ~1, got {s}")
        return self
```

### 4.4 Layer 4 uncertainty sampling

`src/charon/uncertainty/sampling.py` currently samples `clint` as log-normal from conformal CI (lines 183-197). Add a branch:

```python
if clint_prop.source == "classification" and clint_prop.classifier_probs is not None:
    probs = clint_prop.classifier_probs
    samples = _sample_classification_clint(probs, n_samples, rng)
elif clint_prop.ci_90_lower is not None and clint_prop.ci_90_upper is not None:
    samples = _sample_lognormal_from_ci(clint_prop, n_samples, rng)
else:
    samples = np.full(n_samples, clint_prop.value)


def _sample_classification_clint(probs, n, rng):
    """Categorical over {low, med, high} then log-uniform within bucket."""
    bucket_names = ["low", "med", "high"]
    ranges = [(0.1, 10.0), (10.0, 50.0), (50.0, 1000.0)]
    p = np.array([probs[b] for b in bucket_names])
    labels = rng.choice(3, size=n, p=p)
    samples = np.empty(n)
    for i, lbl in enumerate(labels):
        lo, hi = ranges[lbl]
        samples[i] = 10.0 ** rng.uniform(np.log10(lo), np.log10(hi))
    return samples
```

### 4.5 Report narrative

When `classifier_probs` is present, the ADME table for CLint shows:

```markdown
| clint_uL_min_mg | 3.0  | Low bucket [0.1, 10.0] | µL/min/10^6 cells | classification |
```

Plus a new top-level **`## CRITICAL warnings`** section listing every classification-sourced property with its full flag text and actionable guidance ("CLint for this compound is classification-only because the structure is outside Charon's applicability domain. Experimental measurement is required before any dose decision.").

### 4.6 CLint-local AD — `predict/clint_ad.py`

```python
_DEFAULT_FINGERPRINT_CACHE = Path(__file__).resolve().parents[3] / "models" / "clint_ad_fingerprints.npz"
_DEFAULT_TRAINING_CSV = Path(__file__).resolve().parents[3] / "data" / "training" / "clint_merged.csv"
AD_LOW_THRESHOLD = 0.3
AD_MODERATE_THRESHOLD = 0.5


class ClintLocalAD:
    """Tanimoto-based applicability domain for CLint predictions."""

    def __init__(self, reference_fps: list | None = None):
        self._reference_fps = reference_fps or []

    @classmethod
    def load_default(cls, cache_path: Path | None = None) -> "ClintLocalAD":
        """Load or build the cached CLint training FP set."""
        ...

    def max_similarity(self, smiles: str) -> float:
        """Max Tanimoto similarity to training FPs. NaN if SMILES invalid."""
        ...

    def classify(self, max_sim: float) -> Literal["HIGH", "MODERATE", "LOW"]:
        if max_sim >= AD_MODERATE_THRESHOLD: return "HIGH"
        if max_sim >= AD_LOW_THRESHOLD: return "MODERATE"
        return "LOW"
```

Cache is a `.npz` with the bit-vector-style packed Morgan FPs and the source-CSV SHA-256 for invalidation. Same pattern as Sprint 7 `conformal_cache.json`.

---

## 5. Module-by-Module Changes

### New

| File | Responsibility |
|---|---|
| `scripts/train_clint_classifier.py` | Train + persist multi-class CLint classifier. |
| `src/charon/predict/clint_ad.py` | CLint-local AD class + load_default factory + cache. |
| `src/charon/predict/clint_classifier.py` | Wraps loaded XGBoost classifier, exposes `predict_proba`. |
| `models/xgboost_clint_classifier.json` | Trained classifier weights. |
| `models/clint_ad_fingerprints.npz` | Cached training-set Morgan FPs for AD. |
| `tests/unit/test_clint_ad.py` | Similarity, thresholds, cache. |
| `tests/unit/test_clint_classifier.py` | `predict_proba` output shape + normalization. |
| `tests/unit/test_tier3_integration.py` | End-to-end `predict_properties` branching on AD. |
| `tests/unit/test_uncertainty_classification_sampling.py` | Categorical × log-uniform sampler. |

### Modified

| File | Change |
|---|---|
| `src/charon/core/schema.py` | Add `"classification"` to `SourceType` Literal; add `classifier_probs: dict[str, float] \| None` to `PredictedProperty` with validator. |
| `src/charon/predict/__init__.py` | Add `force_tier3: bool = False` kwarg to `predict_properties`; insert AD check + Tier 3 branch after existing ADME prediction. |
| `src/charon/uncertainty/sampling.py` | Add classification-sampling branch before log-normal branch for CLint. |
| `src/charon/report/narrative.py` | Render `classifier_probs` in ADME table; emit new `## CRITICAL warnings` section. |
| `src/charon/report/collector.py` | Pass `classifier_probs` + full `flag` through to payload. |
| `README.md` | Add "CLint Tier 3 and AD" subsection documenting ~50% AD=LOW prevalence on diverse compounds. |
| `models/model_metadata.json` | Add `xgboost_clint_classifier` entry with classifier metrics. |

### Explicitly untouched

- `src/charon/core/guardrails.py` — Layer 0 AD stays as-is (separate effort).
- `src/charon/pbpk/` — PBPK consumes a scalar CLint; point estimate flows normally.
- `src/charon/translational/` — dose projection unchanged.
- Regression goldens — refreshed if narrative changes force a golden update.

---

## 6. Testing Strategy

### 6.1 Unit

- **`test_clint_ad.py`** — known in-domain compound returns HIGH; known out-of-domain (e.g., a sugar, NaCl) returns LOW; cache round-trip; hash invalidation.
- **`test_clint_classifier.py`** — `predict_proba(smi)` returns 3-element dict summing to 1.0; argmax matches a known Obach compound's true bucket.
- **`test_uncertainty_classification_sampling.py`** — sample N=10000, per-bucket fraction matches input probs within 2σ; log-uniform coverage within each bucket.
- **Schema** — `PredictedProperty` rejects `classifier_probs` that don't sum to ~1.

### 6.2 Integration

- **`test_tier3_integration.py`**:
  - A known-diverse SMILES (e.g. a sugar or heavily modified peptide) → `predict_properties` returns `clint.source == "classification"` + CRITICAL flag + `classifier_probs` present.
  - A known in-domain SMILES (caffeine) → Tier 2 (source = `ml_ensemble`) unchanged.
  - `force_tier3=True` on an in-domain compound → Tier 3 anyway.
  - End-to-end `charon report` CLI on an AD=LOW compound → Markdown contains `## CRITICAL warnings` + `Low bucket` text.

### 6.3 Regression

- Re-run `pytest tests/regression/test_known_drugs.py`. None of the 5 reference drugs should flip from Tier 2 to Tier 3 (they are all in-AD). If any do, **investigate** — do not silently update the golden.
- Obach-12 panel: dextromethorphan, omeprazole, theophylline, verapamil are AD=LOW. Their behavior in benchmarks will change when classification fires. Re-run Layer 2/3 benchmarks, review the diff manually, and commit as intentional.

### 6.4 Benchmark

- Layer 2 benchmark (Obach): expected to change for the 4 AD=LOW compounds (`clint_uL_min_mg` now classification-sourced). The Layer 2 benchmark already ships Tier 1 experimental overrides for every compound, so the **point-estimate** AAFE should be unchanged. Layer 4 Monte-Carlo may diverge; document.

---

## 7. Acceptance Criteria

Sprint 8 complete iff all of:

1. `pytest tests/` passes, ~837 → ~880 tests.
2. `scripts/train_clint_classifier.py` runs with scaffold-CV Macro F1 ≥ 0.50 (hard gate).
3. `charon report` on a known AD=LOW compound contains `## CRITICAL warnings` and `Low bucket` substrings.
4. In-domain compounds (caffeine, midazolam) keep `source="ml_ensemble"` — no Tier 2 → Tier 3 regression.
5. Layer 4 uncertainty produces CLint samples with the correct bucket distribution (within 2σ at N=10000).
6. README "CLint Tier 3 and AD" subsection present, stating the ~50% AD-LOW prevalence expectation.
7. No CLAUDE.md MUST-NOT violations (no fu_inc on hepatocytes, no conformal in linear space, etc.).
8. `PredictedProperty.classifier_probs` validator rejects malformed inputs.
9. Schema change to `SourceType` does not break existing tests (updated to accept the new literal).

---

## 8. Risks & Mitigations

| Risk | Likelihood | Mitigation |
|---|---|---|
| Scaffold-CV Macro F1 drops below 0.50 after minor training-data change | Low | Pin `random_state=42`; alert on drift > 0.03 in training-run gate. |
| Adding `classifier_probs` field to `PredictedProperty` breaks existing serialized JSON consumers | Medium | Field is optional (None default); all existing tests continue to pass because absent field serializes as `null` — not a schema break. |
| Tier 3 fires unexpectedly on a production-important in-domain compound due to fingerprint quirk | Medium | `force_tier3=False` override lets users bypass; Obach-12 regression run catches common regressions before merge. |
| Categorical × log-uniform sampler produces physically implausible CLint values at bucket extremes (0.1 or 1000) | Low | Bucket endpoints match CLint physical clip range `[0.1, 1000]`; no extrapolation. |
| AD=LOW prevalence on adme_reference (50%) surprises users and erodes trust | Medium | Document clearly in README and report narrative; explicitly frame as "feature, not bug." |
| `PredictedProperty.source` Literal expansion invalidates existing Pydantic-generated schemas (e.g., OpenAPI) | Low | No Pydantic-backed API consumer exists today; document as a schema evolution in the commit. |
| Training data bucket imbalance post-retraining shifts class probabilities | Low | Sprint 8 does not add training data. Bucket counts are frozen. |

---

## 9. Out-of-Scope Follow-ups

- Per-compound uncertainty signal for Tier 2 (local density / ensemble disagreement).
- Layer 0 `GuardrailChecker` site-wide wiring + project-level reference set.
- Classification fallback for other Layer 1 properties (`fu_p`, `bp_ratio`).
- Calibration plot / reliability diagram shipped with the report.
- Rebalancing training data (adding peptide / macrocycle / natural-product rows).
- Tier 3 for species other than human.

---

## 10. Timeline Estimate

- **Day 1**: `train_clint_classifier.py` + model training + metadata + OOF probs persistence + Macro F1 gate verification.
- **Day 2**: `predict/clint_ad.py` + `predict/clint_classifier.py` + FP cache + load_default.
- **Day 3**: Schema changes (`SourceType`, `classifier_probs`) + `predict_properties` branching + `force_tier3` kwarg + integration test.
- **Day 4**: Layer 4 sampling branch + uncertainty test + report narrative CRITICAL section + regression sweep.
- **Day 5**: README + benchmarks re-run + commit hygiene + PR.

Total: ~5 focused days. No new ML training architecture, no new scientific claim beyond "classifier Macro F1 = 0.54"; all integration on existing Charon patterns established in Sprint 7.
