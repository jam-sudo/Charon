# Phase A Finish — Validation Infrastructure + Honest README

**Date**: 2026-04-15
**Status**: Design (awaiting user approval)
**Owner**: Charon maintainers
**Supersedes**: none

---

## 1. Context & Motivation

Sprint 6 (Report + CLI) closed Phase A feature implementation. Layers 0–4,
the 5 CLI subcommands, and the regulatory Markdown/JSON report are all in
place, with ~760 passing tests. What is **not** yet in place is:

1. A **public, runnable validation benchmark suite** that third parties can
   point at a dataset and reproduce.
2. A **snapshot of current accuracy** against the Obach 1999 panel with
   committed `.md` + `.json` reports under source control.
3. **Two reference compounds (omeprazole, dextromethorphan)** that
   `CLAUDE.md §8` explicitly calls out as required validation drugs.
4. A **regression guard** that prevents silent model drift for the 5
   reference drugs.
5. A **README** that honestly states what Charon does, what it does not
   do, what the current accuracy is, and when it should not be used.

Before any Phase B work or large tech-debt sprint begins, these pieces
need to exist so that future changes can be measured against a stable
baseline.

### Why this design changed after deep review

An initial architecture proposal assumed three benchmark scripts
(Layer 1, Layer 2, Layer 3) could all be populated from existing data,
and that the README could pitch Charon as an open alternative to Simcyp.
Investigation of the validation directories showed:

- `validation/data/tier2_drugs_at_fda/`, `tier3_literature/`, and
  `pk_db/` are **empty**. Layer 3 FIH dose validation has no dataset.
- Obach 1999 Table 6 (the 29-compound in-vivo panel) does **not**
  contain omeprazole or dextromethorphan. Supplementary primary
  literature is required for those two compound YAMLs.
- The current Obach panel snapshot (`with_overrides.txt`) reports
  `AAFE_CL = 4.89`, `AAFE_Vss = 3.32`, `within_3x_CL = 30%`. Phase A
  targets (§8 of CLAUDE.md) are **not met** for most compounds; only
  theophylline, antipyrine, and caffeine are close.

The design therefore narrows to a **validation-infrastructure + honest-
reporting** sprint. Accuracy investigation is explicitly deferred to the
subsequent tech-debt step.

---

## 2. Goals

1. Provide a runnable, reproducible **Layer 1 ADMET benchmark** that
   computes per-property accuracy metrics from the conformal calibration
   set, with explicit disclaimers about what the numbers mean.
2. Extend the existing **Layer 2 human PK benchmark** to emit committed
   `.md` + `.json` reports using a shared writer, without changing its
   existing logic or strict-gate semantics.
3. Add **omeprazole** and **dextromethorphan** compound YAMLs, sourced
   from supplementary primary literature with explicit citations and
   phenotype caveats.
4. Add a **regression test suite** for the 5 CLAUDE.md reference drugs
   that guards against unintended model drift.
5. Write a **README** that is MVP-honest: documents current accuracy
   limits prominently, avoids marketing language, and points to
   reproducible reports.
6. Leave **Layer 3 FIH dose benchmark** as an explicitly deferred stub
   with a documented rationale.

## 3. Non-goals (out of scope for this sprint)

- Investigating or fixing accuracy gaps on any Obach compound. AAFE
  improvement work belongs to the subsequent tech-debt sprint.
- Building the `tier2_drugs_at_fda/` dataset. This is a separate
  literature curation project.
- Training new ML models, re-tuning hyperparameters, or changing the
  conformal calibration procedure.
- Implementing the Sprint 5 ethanol `fu_p` edge case fix or the
  CLint Tier 3 classification fallback (both tracked for tech-debt).
- Any Phase B feature (population variability, DDI, API).
- Any Phase C feature (dashboard, Pareto optimization).
- Replacing the existing `charon.report` package. This sprint's
  `report_writer.py` lives under `validation/benchmarks/` and is
  intentionally separate — it serves aggregate benchmark output,
  not per-compound FIH regulatory reports.

## 4. Success criteria

- `python3 validation/benchmarks/layer1_admet.py` exits 0 and writes
  `validation/reports/layer1_admet.md` and `.json`.
- `python3 validation/benchmarks/layer2_human_pk.py` exits with its
  existing strict-gate semantics and writes
  `validation/reports/layer2_human_pk.md` and `.json`.
- `pytest tests/regression/test_known_drugs.py -v` passes on all 5
  reference drugs, with ±25% tolerance against committed golden
  baselines, in under 30 seconds.
- `pytest tests/` full suite (existing ~760 + new) passes with no
  regressions.
- `README.md` exists, contains actual (not placeholder) AAFE numbers
  captured during this sprint, and links to the committed report files.
- `validation/reports/` and `tests/regression/golden_outputs/` are
  git-tracked and not accidentally ignored.

## 5. Architecture & file layout

### 5.1 New files

```
validation/
  data/tier1_obach/compounds/
    omeprazole.yaml                       [NEW]
    dextromethorphan.yaml                 [NEW]
  benchmarks/
    layer1_admet.py                       [NEW ~200 lines]
    report_writer.py                      [NEW ~120 lines]
  reports/                                [NEW directory, git-tracked]
    layer1_admet.md
    layer1_admet.json
    layer2_human_pk.md
    layer2_human_pk.json

tests/regression/
  test_known_drugs.py                     [NEW ~180 lines]
  golden_outputs/
    midazolam_baseline.json               [NEW]
    warfarin_baseline.json                [NEW]
    diclofenac_baseline.json              [NEW]
    omeprazole_baseline.json              [NEW]
    dextromethorphan_baseline.json        [NEW]

tests/unit/
  test_layer1_benchmark.py                [NEW ~100 lines]
  test_report_writer.py                   [NEW ~80 lines]

README.md                                  [NEW ~350 lines]
```

### 5.2 Modified files

```
validation/data/tier1_obach/panel.yaml    [+2 entries, strict_targets=false]
validation/benchmarks/layer2_human_pk.py  [+emit_report call at end]
validation/benchmarks/metrics.py          [+mae, rmse, pearson_r, within_abs_diff]
validation/benchmarks/layer3_fih_dose.py  [0 lines → docstring stub explaining deferral]
.gitignore                                 [verify validation/reports/ not ignored]
```

### 5.3 Dependency graph

```
compound YAMLs (omeprazole, DM) ──┐
                                   ├─> layer2_human_pk.py ──┐
panel.yaml (+2 entries) ──────────┤                         ├─> validation/reports/*
                                   │                         │
adme_reference.csv ────────────────┼─> layer1_admet.py ─────┘
                                   │        │
metrics.py (expanded) ─────────────┤        │
                                   │        │
report_writer.py ──────────────────┘        │
                                            ▼
                                         README.md
                                            │
                                            ▼
                         tests/regression/test_known_drugs.py
                                            │
                                            ▼
                         tests/regression/golden_outputs/*.json
```

### 5.4 Module boundaries

- **`report_writer.py`** has a single public function `emit_report`.
  It converts a well-typed payload dict into a `.md` + `.json` pair
  on disk. It does not know about PBPK, ADME, or compound semantics.
  It is not used by the production `charon.report` package.
- **`layer1_admet.py`** is a standalone script. It imports
  `charon.predict.predict_properties` and the shared
  `metrics.py` and `report_writer.py`. It does not share state with
  `layer2_human_pk.py`.
- **`test_known_drugs.py`** imports `charon.Pipeline` and the
  compound SMILES directly as module constants. It does not load
  compound YAMLs. This keeps regression tests independent of YAML
  schema changes.

---

## 6. Layer 1 benchmark (`layer1_admet.py`)

### 6.1 Properties evaluated

| Property | Charon source | `adme_reference.csv` column | Source tag in report |
|---|---|---|---|
| `logP` | RDKit Crippen cLogP (descriptor) | `logP` | `descriptor (RDKit Crippen)` |
| `fu_p` | XGBoost ADMETPredictor | `fup` | `ML model (XGBoost)` |
| `bp_ratio` | empirical formula using pKa prior | `rbp` | `empirical formula` |
| `peff_cm_s` | derived from Papp or proxy | `peff_cm_s` | `derived` |

**CLint is explicitly excluded.** Rationale: Charon predicts hepatocyte
CLint in μL/min/10^6 cells. `adme_reference.csv` stores
`clint_3a4_uL_min_pmol` which is recombinant CYP3A4 activity in μL/min
per pmol of enzyme. Converting between them requires CYP3A4 abundance
(nominal 155 pmol/mg HLM) *and* a fraction-metabolized assumption,
introducing enough modeling noise that the AAFE number would not
reflect the underlying ML model. The report's `notes` field documents
this exclusion.

### 6.2 Metric policy per property

- **`logP`**: AAFE is not defined for signed log quantities. Report
  `MAE`, `RMSE`, `R²` in natural units plus `within_0.5_log_pct`
  (fraction of rows with `|pred − obs| ≤ 0.5`). The `MAE < 0.5`
  target is a **sprint-local convention** (not from CLAUDE.md §8),
  consistent with typical medicinal-chemistry tolerance for Crippen
  cLogP accuracy.
- **`fu_p`**: Report `AAFE`, `within_2fold_pct`, `within_3fold_pct`.
  Guard `min_positive = 1e-4`.
- **`bp_ratio`**: Same as `fu_p`. Guard `min_positive = 0.1`.
- **`peff_cm_s`**: Same as `fu_p`. Guard `min_positive = 1e-8`.

### 6.3 Error handling

The script must not crash on any input row. All error modes are
funnelled to an `excluded` list in the payload:

- SMILES parse failure (`ValueError` from `predict_properties`) →
  row excluded wholesale, reason `"SMILES invalid"`.
- Any other exception from `predict_properties` → row excluded,
  reason `"prediction failed: <exc type>"`.
- Missing observed value for a property → that property skipped for
  that row (the row is still counted for other properties).
- `pred ≤ min_positive` or `obs ≤ min_positive` → that property
  skipped for that row, reason `"non-positive value"`.
- Prediction returns `None` for a field → property skipped,
  reason `"prediction missing"`.

Exit codes:

- `0` on normal completion, regardless of whether AAFE targets are met.
- `2` on environmental failure (CSV missing, permission denied,
  model weights missing).

### 6.4 Report payload schema

```python
{
    "title": "Charon Layer 1 ADMET Benchmark",
    "panel": "adme_reference.csv (n=153)",
    "date_utc": "2026-04-15T...",
    "disclaimer": "Evaluated on the conformal calibration set...",
    "summary": {
        "logP":     {"source_tag": "descriptor (RDKit Crippen)",
                     "n": 153, "n_valid": 150,
                     "mae": float, "rmse": float, "r2": float,
                     "within_0.5_log_pct": float},
        "fu_p":     {"source_tag": "ML model (XGBoost)",
                     "n": 153, "n_valid": 148,
                     "aafe": float,
                     "within_2fold_pct": float,
                     "within_3fold_pct": float},
        "bp_ratio": {...},
        "peff_cm_s": {...},
    },
    "rows": [
        {"name": "midazolam", "property": "logP",
         "pred": 3.80, "obs": 3.89, "abs_diff": 0.09},
        {"name": "midazolam", "property": "fu_p",
         "pred": 0.035, "obs": 0.035, "fold_error": 1.00},
        ...
    ],
    "excluded": [
        {"name": "compound_X", "reason": "SMILES invalid"},
        {"name": "compound_Y", "property": "fu_p",
         "reason": "observed fu_p <= 1e-4"},
    ],
    "targets": {
        "fu_p":     {"metric": "aafe", "target": 2.0, "met": bool},
        "bp_ratio": {"metric": "aafe", "target": 2.0, "met": bool},
        "peff_cm_s":{"metric": "aafe", "target": 2.0, "met": bool},
        "logP":     {"metric": "mae",  "target": 0.5, "met": bool},
    },
    "notes": [
        "3 compounds excluded due to SMILES parse errors",
        "CLint not evaluated (unit mismatch)",
    ],
}
```

### 6.5 CLI

```
python3 validation/benchmarks/layer1_admet.py
python3 validation/benchmarks/layer1_admet.py --csv <path>
python3 validation/benchmarks/layer1_admet.py --limit <N>
```

`--limit N` processes the first N rows only (development/debug). No
`--json-only` flag is added in this sprint.

### 6.6 Unit tests

`tests/unit/test_layer1_benchmark.py` (~100 lines):

- Fixture CSV with 6 rows: 2 normal drugs (midazolam, warfarin),
  1 invalid SMILES, 1 with `fu_p = 0`, 1 with missing `logP` column,
  1 with normal values.
- Test that `run_benchmark(mini_csv)` returns a payload dict matching
  the schema.
- Test that invalid SMILES row appears in `excluded`.
- Test that fu_p = 0 row is skipped from fu_p summary but counted
  for logP.
- Test that summary counts match hand-calculated values.
- `report_writer.emit_report` is mocked — this test does not touch
  the filesystem.

---

## 7. `report_writer.py` — shared report emitter

### 7.1 Public API

```python
def emit_report(payload: dict, *, stem: Path) -> tuple[Path, Path]:
    """Write payload to <stem>.md and <stem>.json.

    Returns (md_path, json_path), both absolute and resolved.

    Raises:
        ValueError: if required payload keys missing.
        OSError: if parent dir cannot be created.
    """
```

### 7.2 Payload contract

Required keys: `title`, `panel`, `date_utc`, `summary`, `rows`,
`notes`.

Optional keys: `disclaimer`, `targets`, `excluded`.

Any other keys are ignored by the Markdown renderer but preserved
verbatim in the JSON output.

### 7.3 Markdown rendering

```
# {title}

**Generated**: {date_utc}
**Panel**: {panel}

> **Disclaimer**: {disclaimer}     (only if present)

## Summary

| ... pivoted per-metric table from payload["summary"] ... |

## Per-compound results

| ... table from payload["rows"] ... |

## Targets

| Metric | Target | Value | Status |
| ...

## Excluded

| Name | Property | Reason |
| ...

## Notes

- {note}
- ...
```

Rendering rules:

- Numeric formatting: 4 significant figures by default; `fold_error
  ≥ 10` → `.0f`; values `< 0.1` → scientific notation.
- `NaN` and `±inf` render as `"-"`.
- Table cells escape `|` as `\|`, `\` as `\\`, newlines as a space
  (same convention as `charon.report.narrative._md_cell`).
- Status column uses ASCII markers `[PASS]`, `[FAIL]`, `[ - ]`
  (no emoji — consistent with Sprint 6 reviewer feedback).

### 7.4 JSON rendering

```python
json.dumps(sanitised_payload, indent=2, allow_nan=False,
           default=_json_default)
```

`_json_default` handles:

- `numpy.floating` → `float`
- `numpy.integer` → `int`
- `datetime` → ISO8601 string
- `Path` → POSIX string

Non-finite floats are replaced with `None` before serialisation
(reuses the pattern from `charon.report.export._sanitize_for_json`,
but reimplemented here to avoid cross-package dependency).

### 7.5 Unit tests

`tests/unit/test_report_writer.py` (~80 lines):

- Missing required key raises `ValueError`.
- Empty `rows` and `summary` are allowed.
- Pipe and backslash in cell values are escaped correctly.
- `NaN` and `inf` in payload render as `"-"` in Markdown and
  `null` in JSON.
- Round-trip: `json.loads(written_file) == sanitised(payload)`.
- Unicode compound names (`ω-agatoxin`) render correctly.

---

## 8. Layer 2 integration

### 8.1 Changes to `layer2_human_pk.py`

Minimal diff (~20 lines):

```python
# imports
from datetime import datetime, timezone
from validation.benchmarks.report_writer import emit_report

REPORTS_DIR = REPO_ROOT / "validation" / "reports"

def _build_payload(rows, summary, panel_name):
    return {
        "title": "Charon Layer 2 Human PBPK Benchmark",
        "panel": f"Obach 1999 Tier-1 Panel (n={len(rows)})",
        "date_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "summary": {
            "cl_L_h":   {"aafe": summary.aafe_cl, ...},
            "vss_L":    {"aafe": summary.aafe_vss, ...},
            "t_half_h": {"aafe": summary.aafe_t_half, ...},
        },
        "rows": [_panel_row_to_dict(r) for r in rows],
        "targets": {
            "cl_L_h":   {"metric": "aafe", "target": 2.5, "met": summary.aafe_cl < 2.5},
            "vss_L":    {"metric": "aafe", "target": 3.0, "met": summary.aafe_vss < 3.0},
        },
        "notes": [...],
    }

def main() -> int:
    # (existing logic unchanged)
    ...
    payload = _build_payload(rows, summary, panel_name)
    emit_report(payload, stem=REPORTS_DIR / "layer2_human_pk")
    return 0 if all_strict_pass else 1
```

### 8.2 Invariants preserved

- Existing stdout table format remains byte-identical (used by
  `with_overrides.txt` snapshots as readback).
- Strict-gate exit code semantics unchanged.
- Existing `baseline_rr_only.txt` and `with_overrides.txt` files
  are **not deleted**. They remain for git-history diffing.
- Existing panel YAML loading and Pipeline configuration code
  unchanged.

### 8.3 Effect of new compounds on summary

Adding omeprazole and dextromethorphan to `panel.yaml` (with
`strict_targets: false`) grows n from 10 to 12. Neither compound is
expected to meet the 2-fold CL gate. Summary AAFE will shift
(probably higher); this is the intended outcome — the report
reflects real current state.

---

## 9. Compound YAMLs

> **DATA VERIFICATION PROTOCOL**: All compound YAML values MUST be
> verified against their cited primary sources during implementation
> (Task 3 and Task 4). The implementer should use ChEMBL MCP tools
> (compound_search, get_admet, get_bioactivity) and web search to
> cross-check SMILES (via InChIKey match), MW, and each pharmacological
> parameter. Any value that cannot be confirmed against its cited paper
> should be flagged with `# VERIFY` in the YAML and documented in the
> commit message. Values below are best estimates; several were
> corrected after ChEMBL cross-validation (see §9.4 errata).

### 9.1 `omeprazole.yaml`

```yaml
name: omeprazole
# ChEMBL1503 canonical SMILES (sulfoxide ylide notation):
smiles: "COc1ccc2[nH]c([S+]([O-])Cc3ncc(C)c(OC)c3C)nc2c1"
# Alternative (double-bond S=O): COc1ccc2[nH]c(S(=O)Cc3ncc(C)c(OC)c3C)nc2c1
# Both produce same InChIKey SUBDBMMJDZJVOS-UHFFFAOYSA-N (racemic)
molecular_weight: 345.42   # ChEMBL verified
source: literature
properties:
  physicochemical:
    logp:
      value: 2.23
      source: literature
      method: "DrugBank DB00338 (experimental logP); ChEMBL ALogP=2.9 (calculated)"
    pka_acid:
      value: 8.8
      source: literature
      method: "Andersson T, Clin Pharmacokinet 1996;31(1):9-28"
    pka_base:
      value: 4.0
      source: literature
      method: "Andersson 1996 (pyridine nitrogen)"
    compound_type: ampholyte
  binding:
    fu_p:
      value: 0.05
      source: literature
      unit: fraction
      method: "FDA label + Andersson 1996 (~95% bound, consensus)"
    fu_inc:
      value: 0.97
      source: literature
      unit: fraction
      method: "Austin 2002 estimate from logP=2.23; VERIFY against Lopez-Rodriguez 2010 Xenobiotica 40:617 if direct measurement available"
    bp_ratio:
      value: 0.59
      source: literature
      unit: ratio
      method: "DrugBank DB00338"
  metabolism:
    clint_uL_min_mg:
      value: 110.0
      source: literature
      unit: uL/min/mg
      method: "VERIFY: check Lopez-Rodriguez 2010 Table 2 or Naritomi 2001 Drug Metab Dispos 29:1316 for HLM CLint"
      # If Lopez-Rodriguez 2010 reports Vmax/Km rather than direct CLint,
      # derive CLint = Vmax/Km and note the conversion.
    fm_cyp2c19: 0.8
    fm_cyp3a4: 0.2
  renal:
    clrenal_L_h:
      value: 0.0
      source: literature
      method: "Andersson 1996 (negligible urinary excretion of parent)"
observed_iv:
  dose_mg: 40.0
  cl_L_h: 34.2        # ~0.5 L/min, high extraction, EM population
  vss_L: 24.5          # ~0.35 L/kg × 70 kg
  t_half_h: 0.9
  source: "Andersson T, Clin Pharmacokinet 1996;31(1):9-28 (EM population mean)"
  phenotype_caveat: "EM population; CYP2C19 PMs ~3x lower CL (not modelled)"
```

### 9.2 `dextromethorphan.yaml`

```yaml
name: dextromethorphan
# ChEMBL52440 canonical SMILES (verified stereochemistry):
smiles: "COc1ccc2c(c1)[C@]13CCCC[C@@H]1[C@H](C2)N(C)CC3"
# InChIKey: MKXZASYAUGDDCJ-NJAFHUGGSA-N (dextro enantiomer)
molecular_weight: 271.40   # ChEMBL verified (C18H25NO)
source: literature
properties:
  physicochemical:
    logp:
      # CORRECTED: original spec had 4.11 which is too high.
      # ChEMBL ALogP = 3.38 (calculated, Wildman-Crippen).
      # Literature experimental range: 2.7-3.4 (method-dependent).
      # VERIFY against DrugBank DB00514 experimental logP during implementation.
      value: 3.0
      source: literature
      method: "Literature consensus (range 2.7-3.4); VERIFY DB00514"
    pka_base:
      # Some sources cite 8.3, others 9.2. VERIFY primary source.
      value: 8.3
      source: literature
      method: "Literature (tertiary amine morphinan); range 8.3-9.2 in sources"
    compound_type: base
  binding:
    fu_p:
      value: 0.68
      source: literature
      unit: fraction
      method: "Bolaji 1993 Br J Clin Pharmacol 36:131 (equilibrium dialysis)"
    fu_inc:
      # Recalculated with corrected logP ~3.0 (Austin 2002 formula).
      # Exact value depends on final logP; recalculate during implementation.
      value: 0.80
      source: literature
      unit: fraction
      method: "Austin 2002 estimate; recalculate with verified logP"
    bp_ratio:
      value: 1.04
      source: literature
      unit: ratio
      method: "DrugBank DB00514"
  metabolism:
    clint_uL_min_mg:
      # CORRECTED: Bolaji 1993 reports protein binding and in vivo PK,
      # NOT HLM CLint. Need proper HLM microsomal CLint source.
      # Candidates: Obach 1999 Table 2, Ito & Houston 1997, or
      # Hallifax & Houston 2006. VERIFY during implementation.
      value: 85.0
      source: literature
      unit: uL/min/mg
      method: "VERIFY: check Obach 1999 Table 2 or Houston group publications for HLM CLint"
    fm_cyp2d6: 0.85
    fm_cyp3a4: 0.15
  renal:
    clrenal_L_h:
      value: 0.0
      source: literature
      method: "negligible urinary excretion of parent"
observed_iv:
  # CORRECTED: 80 L/h exceeds typical EM range (40-70 L/h).
  # Bolaji 1993 reports in vivo CL from IV dosing in EMs.
  # Schadel 1995 J Clin Psychopharmacol 15:263 also reports EM PK.
  # VERIFY exact values against primary sources.
  dose_mg: 30.0
  cl_L_h: 50.0          # conservative middle of EM range (40-70 L/h)
  vss_L: 250.0          # ~3.6 L/kg, plausible for basic lipophilic drug
  t_half_h: 3.5         # adjusted for lower CL (t½ = 0.693 × Vss / CL)
  source: "Schadel 1995 + Bolaji 1993 (EM population); VERIFY exact values"
  phenotype_caveat: "EM population; CYP2D6 PMs ~10x lower CL (not modelled)"
```

### 9.4 Errata from post-design ChEMBL cross-validation

The following corrections were made after querying ChEMBL (CHEMBL1503,
CHEMBL52440) and cross-referencing citations:

1. **Omeprazole citation**: `Andersson 1990 Clin Pharmacokinet 18:93`
   does not appear to exist as cited. Corrected to `Andersson T, Clin
   Pharmacokinet 1996;31(1):9-28` — the canonical omeprazole PK review.
2. **DM SMILES**: Original had `[C@@]1([H])[C@@H](C2)` which differs
   from ChEMBL canonical `[C@@H]1[C@H](C2)` at one stereocenter.
   Corrected to ChEMBL52440 canonical (InChIKey verified).
3. **DM logP**: 4.11 → 3.0 (range 2.7–3.4). Original value is outside
   the published experimental range and above ChEMBL ALogP (3.38).
4. **DM CLint attribution**: Bolaji 1993 reports protein binding and
   in vivo PK, not HLM CLint. CLint source needs verification against
   Obach 1999 or Houston-group publications.
5. **DM observed CL**: 80 → 50 L/h. 80 approaches hepatic blood flow
   (~90 L/h) which implies >100% hepatic extraction — physiologically
   unlikely. EM range is 40–70 L/h.
6. **DM pKa**: 9.2 → 8.3. Multiple sources report 8.3; 9.2 may
   reflect a different measurement method. flagged for verification.
7. **DM fu_inc**: 0.85 → 0.80. Recalculated with corrected logP ≈ 3.0
   (Austin 2002 formula); exact value depends on final logP.
8. **DM t_half**: 3.0 → 3.5 h. Adjusted for consistency with CL=50
   and Vss=250 (t½ = 0.693 × 250 / 50 = 3.47 h).

### 9.3 `panel.yaml` entries

```yaml
  - key: omeprazole
    compound_file: compounds/omeprazole.yaml
    route: iv_bolus
    dose_mg: 40.0
    duration_h: 24.0
    observed: {cl_L_h: 34.2, vss_L: 24.5, t_half_h: 0.9}
    strict_targets: false
    notes: >
      Ampholyte, CYP2C19 dominant. Not in Obach 1999; supplementary
      literature source. EM population means; PM phenotype not modelled.

  - key: dextromethorphan
    compound_file: compounds/dextromethorphan.yaml
    route: iv_bolus
    dose_mg: 30.0
    duration_h: 48.0
    observed: {cl_L_h: 50.0, vss_L: 250.0, t_half_h: 3.5}
    strict_targets: false
    notes: >
      Basic, lipophilic, CYP2D6 dominant. Not in Obach 1999. High Vss,
      phenotype-sensitive. EM population means; PM phenotype not modelled.
      Observed PK values need primary-source verification (see §9.4).
```

---

## 10. Regression tests

### 10.1 Purpose

The regression suite guards against **unintended model drift**. It does
not validate accuracy against literature — that is the job of the
Layer 2 benchmark. It only checks: "if we change the model, does the
pipeline still produce output within ±25% of the last committed
snapshot?"

### 10.2 Structure

`tests/regression/test_known_drugs.py` (~180 lines):

- 5 reference drugs as module-level `REFERENCE_DRUGS` tuple.
- Each drug has: `key`, `smiles`, `route`, `dose_mg`,
  `noael_mg_kg`, `noael_species`, `target_kd_nM`,
  `target_ceff_nM`, `phenotype_caveat`.
- Two parameterised tests:
  - `test_pipeline_runs(drug)`: runs pipeline, asserts `cmax > 0`,
    `auc_0_inf > 0`, `mrsd_mg` is finite and positive (if a
    dose recommendation was produced).
  - `test_dose_within_snapshot(drug)`: compares current primary
    metric to committed golden baseline with fold-error tolerance
    of 1.25 (±25%).

### 10.3 Primary metric selection

```
if dose_recommendation is not None:
    primary = ("mrsd_mg", rec.mrsd_mg)
elif auc_0_inf > 0:
    primary = ("auc_0_inf", pk.auc_0_inf)
else:
    primary = ("cmax", pk.cmax)
```

This allows omeprazole and dextromethorphan to be tested without
fabricated MABEL or NOAEL parameters (neither has public IND data).

### 10.4 Golden baseline schema

`tests/regression/golden_outputs/<drug>_baseline.json`:

```json
{
  "drug": "midazolam",
  "smiles": "Cc1ncc2n1...",
  "charon_version": "0.1.0-alpha",
  "primary_metric": {"name": "mrsd_mg", "value": 1.234},
  "secondary_metrics": {
    "cmax_ug_L": 45.6,
    "auc_ug_h_L": 123.4,
    "t_half_h": 3.5
  },
  "generated_utc": "2026-04-15T12:34:56+00:00",
  "caveat": null
}
```

`caveat` is populated with a string like `"CYP polymorphism; baseline
may drift"` for omeprazole and dextromethorphan; otherwise `null`.

### 10.5 Baseline regeneration workflow

```bash
UPDATE_BASELINES=1 pytest tests/regression/test_known_drugs.py
```

When `UPDATE_BASELINES=1` is set:

- `test_dose_within_snapshot` writes the current value as the new
  baseline and calls `pytest.skip(f"Updated baseline for {drug.key}")`.
- This is the only way baselines are updated. Reviewers can inspect
  the diff in `golden_outputs/*.json` to verify the update is
  intentional.

If a baseline file is missing (not covered by `UPDATE_BASELINES`),
the test **fails** explicitly with a message pointing at the
regeneration command. This prevents silent CI passes when a new
drug is added without a baseline.

### 10.6 Timing budget

Pipeline run without uncertainty: ~3–5 s per drug.
Total for 5 drugs × 2 tests: ~20–25 s.

Target: **regression suite completes in under 30 seconds** on a
cold cache. The suite runs with `uncertainty=None` to avoid the
LHS 500-sample overhead.

---

## 11. README structure

### 11.1 Table of contents (fixed order)

1. Charon — one-paragraph description
2. Project status (Phase A MVP)
3. **Known limitations** (prominent — directly after TOC)
4. Install
5. Quickstart (CLI + Python API)
6. Reference drug quick-check
7. Validation status (Layer 1 / Layer 2 / Layer 3-deferred)
8. Architecture overview (pointer to `ARCHITECTURE.md`)
9. Development
10. Roadmap (Phase B, Phase C)
11. Citing Charon
12. License
13. Acknowledgments

### 11.2 Tone rules

- First-person plural (`we`, `our`) when describing project
  decisions.
- **Banned phrases**: "alternative to", "replaces", "industry-
  leading", "state-of-the-art", "best-in-class", "production-ready".
- Concrete numbers always. `"AAFE ≈ 4.9 on Obach n=12"` is
  acceptable; `"limited accuracy on some compounds"` is not.
- Failed gates appear in README tables as `✗` or `[FAIL]`. No
  whitewashing.
- Relative links only (e.g., `[Layer 2 report](validation/reports/layer2_human_pk.md)`).

### 11.3 Known limitations section (content outline)

The limitations section lists, at minimum, the following items, each
as a single bullet or short paragraph:

- Accuracy vs Phase A targets (actual AAFE numbers, pointed at the
  Layer 2 report).
- CLint Tier 2 ML limit (scaffold-CV AAFE ≈ 2.5, experimental
  override recommended).
- Layer 3 FIH dose validation is deferred.
- No CYP phenotype modelling.
- Conformal CI is marginal, not conditional (§6d of CLAUDE.md).
- Single-subject only (no population variability).
- PBPK ODE reproducibility is `rtol=1e-6`, not bit-identical.
- R&R `Kp` overpredicts adipose for lipophilic weak bases;
  empirical override supported per compound.
- Sprint 5 LHS edge case for `fu_p ≈ 1.0` compounds.
- No DDI, no metabolite tracking, no explicit transporter modelling
  beyond P-gp proxy.
- Not FDA-cleared. Research / educational tool only.
- Salt-form MW and salt correction is the user's responsibility.

### 11.4 Validation-status tables

Layer 1 table (actual values filled during implementation):

| Property | n | Metric | Value | Target | Met |
|---|---|---|---|---|---|
| logP (Crippen) | ~150 | MAE | TBD | <0.5 | TBD |
| fu_p (XGBoost) | ~148 | AAFE | TBD | <2.0 | TBD |
| bp_ratio (formula) | ~148 | AAFE | TBD | <2.0 | TBD |
| peff_cm_s (derived) | ~140 | AAFE | TBD | <2.0 | TBD |
| CLint | — | — | excluded | — | — |

Link: `validation/reports/layer1_admet.md`

Layer 2 table (actual values filled during implementation):

| Metric | AAFE | within-2× | within-3× | Target AAFE | Met |
|---|---|---|---|---|---|
| CL | TBD | TBD | TBD | <2.5 | TBD |
| Vss | TBD | TBD | TBD | <3.0 | TBD |
| t½ | TBD | TBD | TBD | — | — |

Link: `validation/reports/layer2_human_pk.md`

Layer 3: `DEFERRED. No public tier2_drugs_at_fda dataset has been
built. FIH dose accuracy is not yet validated.`

### 11.5 TBD handling during implementation

All `TBD` fields in the README are substituted with real values
during Task 13 (README write) after Task 9 (benchmark run) has
produced the actual numbers. The plan's Task 13 step explicitly
blocks until Task 9 reports exist.

---

## 12. Implementation workflow

### 12.1 Task list (17 tasks)

Phase 1 — infrastructure and data:

| # | Task |
|---|---|
| 1 | Implement `validation/benchmarks/report_writer.py` + unit tests |
| 2 | Extend `validation/benchmarks/metrics.py` with `mae`, `rmse`, `pearson_r`, `within_abs_diff` |
| 3 | Write `validation/data/tier1_obach/compounds/omeprazole.yaml` — verify all values via ChEMBL/web search against cited sources (§9.1 + §9.4 errata) |
| 4 | Write `validation/data/tier1_obach/compounds/dextromethorphan.yaml` — verify all values via ChEMBL/web search, especially logP, CLint source, pKa (§9.2 + §9.4 errata) |
| 5 | Extend `validation/data/tier1_obach/panel.yaml` with 2 new entries |

Phase 2 — benchmark scripts:

| # | Task |
|---|---|
| 6 | Implement `validation/benchmarks/layer1_admet.py` + unit tests |
| 7 | Modify `validation/benchmarks/layer2_human_pk.py` to call `emit_report` at end of `main` |
| 8 | Replace `validation/benchmarks/layer3_fih_dose.py` empty file with docstring-only stub documenting deferral |
| 9 | Run both benchmarks; commit `validation/reports/*.md` and `*.json` |

Phase 3 — regression suite:

| # | Task |
|---|---|
| 10 | Implement `tests/regression/test_known_drugs.py` |
| 11 | Regenerate baselines: `UPDATE_BASELINES=1 pytest …`; commit `golden_outputs/*.json` |
| 12 | Verify regression suite passes with default tolerance on all 5 drugs |

Phase 4 — documentation:

| # | Task |
|---|---|
| 13 | Write `README.md` substituting real numbers from Task 9 reports |
| 14 | Verify `.gitignore` does not accidentally exclude `validation/reports/` or `golden_outputs/` |
| 15 | Final sweep: remove all `TBD`/`FIXME`, verify all relative links |

Phase 5 — closing:

| # | Task |
|---|---|
| 16 | Run full `pytest tests/` (existing ~760 + new); confirm zero regressions |
| 17 | Push all commits to `origin/main` |

### 12.2 Dependency graph

```
T1 ──┬──> T6
      └──> T7
T2 ───> T6
T3 ──┐
T4 ──┼──> T5 ──> T7
T6, T7 ──> T9 ──> T13 ──> T15 ──> T16 ──> T17
T8 ──> T13 (deferral text references T8 stub)
T10 ──> T11 ──> T12 ──> T16
T14 ──> T16
```

No parallelisable branches within a subagent dispatch. T3 and T4
are independent and may be dispatched in one prompt.

### 12.3 Subagent allocation (subagent-driven-development skill)

- **Large new files** (T1, T6, T10, T13): implementer subagent +
  spec reviewer + code-quality reviewer.
- **Small modifications** (T2, T5, T7, T8, T14): single implementer
  subagent with self-review only.
- **YAML data** (T3, T4): implementer verifies citations and YAML
  schema; no code review step.
- **Non-code tasks** (T9, T11, T17): direct execution in the
  controller session (these are running commands, not writing
  code).

### 12.4 Commit policy

- One commit per task (atomic).
- Task 9 and Task 11 produce `chore:` commits with generated
  artefacts only (no code change).
- Task 13 (README) goes in after Task 9 so the README can quote
  real numbers.
- Task 17 is the single `git push`.

### 12.5 Failure handling

A task is complete only when:

- For code changes: `pytest tests/` is fully green.
- For YAML additions: `python -c "import yaml; yaml.safe_load(open(<path>))"` succeeds.
- For benchmark runs: exit 0 and both output files exist.
- For README: all relative links resolve and no `TBD` remains.

If a task fails at review, the controller dispatches a fix
subagent with specific instructions, then re-runs the relevant
reviewer. Rollback to the previous commit is used only if the
fix attempts are not converging.

---

## 13. Open questions and risks

### 13.1 Risks

- **Layer 2 AAFE may worsen after adding omeprazole and DM.** The
  summary `AAFE_CL` could go from 4.89 to 6+. This is an expected
  and intended outcome — the report reflects reality. README text
  is written after Task 9 so it quotes actual numbers, not
  pre-sprint expectations.
- **Golden baseline regeneration is one-shot.** If a baseline
  turns out to be wrong on first capture (e.g., a stale model
  file), Task 12 will fail and the baseline must be regenerated.
  This is fine — `UPDATE_BASELINES=1` exists for exactly this.
- **adme_reference.csv unit handling.** The benchmark excludes
  CLint; this must be clearly documented in the report payload and
  unit tests so a future contributor doesn't re-enable CLint
  comparison and get garbage numbers.
- **Citation accuracy for omeprazole / DM YAMLs.** Values are
  sourced from secondary/tertiary sources (DrugBank, textbooks) as
  well as primary literature. Each field carries a `method:`
  citation string. Reviewers should verify the citations exist but
  need not re-derive the numbers.
- **Phase A is declared MVP, not validated-to-target.** The
  README wording needs user approval before it ships publicly. If
  the user disagrees with the honest-reporting tone, the sprint
  produces everything except the README and stops for re-design.

### 13.2 Explicitly deferred

- Layer 3 FIH dose benchmark script (requires dataset).
- Layer 2 accuracy investigation (belongs to tech-debt sprint).
- Sprint 5 ethanol `fu_p` fix (tech-debt sprint).
- CLint Tier 3 classification fallback (tech-debt sprint).
- tier2_drugs_at_fda dataset curation.
- CYP phenotype modelling.

---

## 14. Approval gate

This spec is ready for user review. After approval:

1. This spec is committed to git.
2. The `writing-plans` skill is invoked to generate a detailed
   step-by-step implementation plan.
3. Implementation proceeds via `subagent-driven-development` skill
   (one task, one subagent, two-stage review).

No code changes have been made yet.
