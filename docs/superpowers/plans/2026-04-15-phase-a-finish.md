# Phase A Finish — Validation Infrastructure + Honest README

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build runnable Layer 1/2 benchmarks, add omeprazole/dextromethorphan compound YAMLs, create regression guards for 5 reference drugs, and write an MVP-honest README.

**Architecture:** Shared `report_writer.py` emits `.md` + `.json` for both Layer 1 and Layer 2 benchmarks. Regression tests use golden-output snapshots with ±25% tolerance. README quotes actual numbers from committed reports. Layer 3 benchmark explicitly deferred (no dataset).

**Tech Stack:** Python 3.11+, pytest, PyYAML, numpy, charon (existing package). No new dependencies.

**Spec:** `docs/superpowers/specs/2026-04-15-phase-a-finish-design.md`

---

## File Structure

### New files
| File | Responsibility |
|---|---|
| `validation/benchmarks/report_writer.py` | Shared dict → `.md` + `.json` emitter |
| `validation/benchmarks/layer1_admet.py` | Layer 1 ADMET benchmark (adme_reference.csv) |
| `validation/data/tier1_obach/compounds/omeprazole.yaml` | Compound YAML with literature overrides |
| `validation/data/tier1_obach/compounds/dextromethorphan.yaml` | Compound YAML with literature overrides |
| `validation/reports/layer1_admet.md` | Generated Layer 1 report (git-tracked) |
| `validation/reports/layer1_admet.json` | Generated Layer 1 report JSON |
| `validation/reports/layer2_human_pk.md` | Generated Layer 2 report |
| `validation/reports/layer2_human_pk.json` | Generated Layer 2 report JSON |
| `tests/unit/test_report_writer.py` | Unit tests for report_writer |
| `tests/unit/test_layer1_benchmark.py` | Unit tests for layer1_admet |
| `tests/regression/test_known_drugs.py` | Regression guards (5 drugs × 2 tests) |
| `tests/regression/golden_outputs/*.json` | 5 golden baseline files |
| `README.md` | Project README (~350 lines) |

### Modified files
| File | Change |
|---|---|
| `validation/benchmarks/metrics.py` | +4 functions: `mae`, `rmse`, `pearson_r`, `within_abs_diff` |
| `validation/data/tier1_obach/panel.yaml` | +2 compound entries |
| `validation/benchmarks/layer2_human_pk.py` | +emit_report call at end of main() |
| `validation/benchmarks/layer3_fih_dose.py` | Empty → docstring stub |

---

## Task 1: Implement `report_writer.py` + unit tests

**Files:**
- Create: `validation/benchmarks/report_writer.py`
- Create: `tests/unit/test_report_writer.py`

- [ ] **Step 1: Write failing tests**

```python
# tests/unit/test_report_writer.py
"""Unit tests for validation/benchmarks/report_writer.py."""
from __future__ import annotations

import json
import math
import tempfile
from pathlib import Path

import pytest


def _minimal_payload(**overrides) -> dict:
    base = {
        "title": "Test Report",
        "panel": "test_panel",
        "date_utc": "2026-04-15T00:00:00+00:00",
        "summary": {"metric_a": {"aafe": 2.1, "n": 10}},
        "rows": [{"name": "drug_a", "property": "fu_p", "pred": 0.5, "obs": 0.4}],
        "notes": ["note one"],
    }
    base.update(overrides)
    return base


class TestEmitReport:
    def test_creates_md_and_json(self, tmp_path):
        from validation.benchmarks.report_writer import emit_report

        payload = _minimal_payload()
        md, js = emit_report(payload, stem=tmp_path / "test_out")
        assert md.exists() and md.suffix == ".md"
        assert js.exists() and js.suffix == ".json"

    def test_missing_required_key_raises(self, tmp_path):
        from validation.benchmarks.report_writer import emit_report

        bad = {"title": "only title"}
        with pytest.raises(ValueError, match="missing required"):
            emit_report(bad, stem=tmp_path / "bad")

    def test_empty_rows_allowed(self, tmp_path):
        from validation.benchmarks.report_writer import emit_report

        payload = _minimal_payload(rows=[], summary={})
        md, js = emit_report(payload, stem=tmp_path / "empty")
        assert md.exists()

    def test_pipe_escape_in_md(self, tmp_path):
        from validation.benchmarks.report_writer import emit_report

        payload = _minimal_payload(
            rows=[{"name": "a|b", "property": "x", "pred": 1.0, "obs": 1.0}]
        )
        md, _ = emit_report(payload, stem=tmp_path / "esc")
        content = md.read_text()
        assert "a\\|b" in content
        assert "a|b" not in content.split("\n")[-10:]  # not raw in table

    def test_nan_renders_as_dash_in_md(self, tmp_path):
        from validation.benchmarks.report_writer import emit_report

        payload = _minimal_payload(
            rows=[{"name": "d", "property": "x", "pred": float("nan"), "obs": 1.0}]
        )
        md, _ = emit_report(payload, stem=tmp_path / "nan")
        assert "-" in md.read_text()

    def test_inf_becomes_null_in_json(self, tmp_path):
        from validation.benchmarks.report_writer import emit_report

        payload = _minimal_payload(
            summary={"m": {"value": float("inf")}}
        )
        _, js = emit_report(payload, stem=tmp_path / "inf")
        data = json.loads(js.read_text())
        assert data["summary"]["m"]["value"] is None

    def test_json_roundtrip(self, tmp_path):
        from validation.benchmarks.report_writer import emit_report

        payload = _minimal_payload()
        _, js = emit_report(payload, stem=tmp_path / "rt")
        data = json.loads(js.read_text())
        assert data["title"] == "Test Report"
        assert data["rows"][0]["name"] == "drug_a"

    def test_disclaimer_rendered_when_present(self, tmp_path):
        from validation.benchmarks.report_writer import emit_report

        payload = _minimal_payload(disclaimer="Important caveat")
        md, _ = emit_report(payload, stem=tmp_path / "disc")
        assert "Important caveat" in md.read_text()
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/unit/test_report_writer.py -v`
Expected: FAIL (ImportError: No module named 'validation.benchmarks.report_writer')

- [ ] **Step 3: Implement report_writer.py**

```python
# validation/benchmarks/report_writer.py
"""Shared benchmark report emitter — dict payload to .md + .json.

Used by layer1_admet.py and layer2_human_pk.py to produce committed
validation reports.  Not part of the charon package — lives under
validation/benchmarks/ intentionally.
"""
from __future__ import annotations

import json
import math
from datetime import datetime
from pathlib import Path
from typing import Any

_REQUIRED_KEYS = ("title", "panel", "date_utc", "summary", "rows", "notes")


def _fmt(v: Any, digits: int = 4) -> str:
    if v is None:
        return "-"
    if isinstance(v, bool):
        return "[PASS]" if v else "[FAIL]"
    if isinstance(v, (int,)):
        return str(v)
    try:
        fv = float(v)
    except (TypeError, ValueError):
        return str(v)
    if not math.isfinite(fv):
        return "-"
    if fv == 0:
        return "0"
    if abs(fv) >= 1000 or abs(fv) < 0.01:
        return f"{fv:.{digits - 1}e}"
    return f"{fv:.{digits}g}"


def _md_cell(s: Any) -> str:
    text = _fmt(s) if not isinstance(s, str) else s
    text = text.replace("\\", "\\\\").replace("|", "\\|")
    text = text.replace("\n", " ").replace("\r", " ")
    return text or "-"


def _render_table(headers: list[str], rows: list[dict]) -> str:
    if not rows:
        return "(no data)\n"
    cols = headers
    lines = [
        "| " + " | ".join(cols) + " |",
        "| " + " | ".join("---" for _ in cols) + " |",
    ]
    for row in rows:
        cells = [_md_cell(row.get(c, "-")) for c in cols]
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines) + "\n"


def _sanitize(obj: Any) -> Any:
    if isinstance(obj, float):
        if not math.isfinite(obj):
            return None
        return obj
    if isinstance(obj, dict):
        return {k: _sanitize(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_sanitize(x) for x in obj]
    if isinstance(obj, datetime):
        return obj.isoformat()
    if isinstance(obj, Path):
        return str(obj)
    if hasattr(obj, "item"):  # numpy scalar
        return _sanitize(obj.item())
    return obj


def _render_md(payload: dict) -> str:
    parts: list[str] = []
    parts.append(f"# {payload['title']}\n")
    parts.append(f"**Generated**: {payload['date_utc']}")
    parts.append(f"**Panel**: {payload['panel']}\n")

    if payload.get("disclaimer"):
        parts.append(f"> **Disclaimer**: {payload['disclaimer']}\n")

    # Summary
    parts.append("## Summary\n")
    summary = payload["summary"]
    if summary:
        s_headers = ["metric"]
        if summary:
            first = next(iter(summary.values()))
            s_headers += [k for k in first if k != "source_tag"]
        s_rows = []
        for key, vals in summary.items():
            row = {"metric": key}
            row.update(vals)
            s_rows.append(row)
        parts.append(_render_table(s_headers, s_rows))

    # Targets
    targets = payload.get("targets")
    if targets:
        parts.append("## Targets\n")
        t_rows = []
        for key, tgt in targets.items():
            t_rows.append({
                "metric": key,
                "target": tgt.get("target", "-"),
                "met": tgt.get("met"),
            })
        parts.append(_render_table(["metric", "target", "met"], t_rows))

    # Per-compound rows
    parts.append("## Per-compound results\n")
    rows = payload["rows"]
    if rows:
        r_headers = list(rows[0].keys())
        parts.append(_render_table(r_headers, rows))
    else:
        parts.append("(no rows)\n")

    # Excluded
    excluded = payload.get("excluded")
    if excluded:
        parts.append("## Excluded\n")
        e_headers = list(excluded[0].keys())
        parts.append(_render_table(e_headers, excluded))

    # Notes
    parts.append("## Notes\n")
    for note in payload.get("notes", []):
        parts.append(f"- {note}")
    parts.append("")

    return "\n".join(parts)


def emit_report(payload: dict, *, stem: Path) -> tuple[Path, Path]:
    """Write payload as .md + .json to disk.

    Parameters
    ----------
    payload : dict
        Must contain keys: title, panel, date_utc, summary, rows, notes.
    stem : Path
        Base path without extension (e.g., ``reports/layer1_admet``).

    Returns
    -------
    (md_path, json_path) : tuple[Path, Path]
        Resolved absolute paths of the written files.

    Raises
    ------
    ValueError
        If required payload keys are missing.
    """
    missing = [k for k in _REQUIRED_KEYS if k not in payload]
    if missing:
        raise ValueError(f"Payload missing required keys: {missing}")

    stem = Path(stem).resolve()
    stem.parent.mkdir(parents=True, exist_ok=True)

    md_path = stem.with_suffix(".md")
    json_path = stem.with_suffix(".json")

    md_path.write_text(_render_md(payload), encoding="utf-8")
    json_path.write_text(
        json.dumps(_sanitize(payload), indent=2, allow_nan=False),
        encoding="utf-8",
    )
    return md_path, json_path
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/unit/test_report_writer.py -v`
Expected: 8 passed

- [ ] **Step 5: Commit**

```bash
git add validation/benchmarks/report_writer.py tests/unit/test_report_writer.py
git commit -m "feat(validation): add report_writer shared benchmark emitter"
```

---

## Task 2: Extend metrics.py with linear-scale functions

**Files:**
- Modify: `validation/benchmarks/metrics.py`

- [ ] **Step 1: Write tests (append to a new test file or inline)**

The new functions are simple enough that the Layer 1 benchmark tests (Task 6) will cover them. For this task, add the functions and verify the existing `validation/benchmarks/metrics.py` tests (if any) still pass.

- [ ] **Step 2: Add 4 functions to metrics.py**

Append to `validation/benchmarks/metrics.py` after the existing `within_n_fold` function:

```python
def mae(predicted: Iterable[float], observed: Iterable[float]) -> float:
    """Mean Absolute Error."""
    preds = list(predicted)
    obs = list(observed)
    if len(preds) != len(obs):
        raise ValueError(f"mae: length mismatch ({len(preds)} vs {len(obs)})")
    if not preds:
        raise ValueError("mae: empty input")
    return sum(abs(p - o) for p, o in zip(preds, obs)) / len(preds)


def rmse(predicted: Iterable[float], observed: Iterable[float]) -> float:
    """Root Mean Squared Error."""
    preds = list(predicted)
    obs = list(observed)
    if len(preds) != len(obs):
        raise ValueError(f"rmse: length mismatch ({len(preds)} vs {len(obs)})")
    if not preds:
        raise ValueError("rmse: empty input")
    mse = sum((p - o) ** 2 for p, o in zip(preds, obs)) / len(preds)
    return mse ** 0.5


def pearson_r(predicted: Iterable[float], observed: Iterable[float]) -> float:
    """Pearson correlation coefficient. Returns 0.0 if stdev is zero."""
    preds = list(predicted)
    obs = list(observed)
    if len(preds) != len(obs):
        raise ValueError(
            f"pearson_r: length mismatch ({len(preds)} vs {len(obs)})"
        )
    n = len(preds)
    if n < 2:
        return 0.0
    mean_p = sum(preds) / n
    mean_o = sum(obs) / n
    num = sum((p - mean_p) * (o - mean_o) for p, o in zip(preds, obs))
    den_p = sum((p - mean_p) ** 2 for p in preds) ** 0.5
    den_o = sum((o - mean_o) ** 2 for o in obs) ** 0.5
    if den_p == 0 or den_o == 0:
        return 0.0
    return num / (den_p * den_o)


def within_abs_diff(
    predicted: Iterable[float],
    observed: Iterable[float],
    *,
    threshold: float = 0.5,
) -> float:
    """Fraction of predictions with |pred - obs| <= threshold."""
    preds = list(predicted)
    obs = list(observed)
    if len(preds) != len(obs):
        raise ValueError(
            f"within_abs_diff: length mismatch ({len(preds)} vs {len(obs)})"
        )
    if not preds:
        raise ValueError("within_abs_diff: empty input")
    hits = sum(1 for p, o in zip(preds, obs) if abs(p - o) <= threshold)
    return hits / len(preds)
```

- [ ] **Step 3: Run full test suite to verify no regressions**

Run: `pytest tests/ -x -q`
Expected: all existing tests pass

- [ ] **Step 4: Commit**

```bash
git add validation/benchmarks/metrics.py
git commit -m "feat(validation): add mae, rmse, pearson_r, within_abs_diff to metrics"
```

---

## Task 3: Write omeprazole.yaml

**Files:**
- Create: `validation/data/tier1_obach/compounds/omeprazole.yaml`

**IMPORTANT**: Before writing the YAML, the implementer MUST verify every pharmacological value against its cited source using ChEMBL MCP tools (`mcp__claude_ai_ChEMBL__compound_search` for CHEMBL1503) and web search. See spec §9.1 + §9.4 errata. Values below are best estimates from the spec — adjust if verification yields different numbers.

- [ ] **Step 1: Verify SMILES and MW against ChEMBL**

Use `mcp__claude_ai_ChEMBL__compound_search(name="omeprazole", max_phase=4)` to confirm:
- SMILES: `COc1ccc2[nH]c([S+]([O-])Cc3ncc(C)c(OC)c3C)nc2c1` (CHEMBL1503)
- InChIKey: `SUBDBMMJDZJVOS-UHFFFAOYSA-N`
- MW freebase: 345.42

- [ ] **Step 2: Create omeprazole.yaml**

```yaml
# validation/data/tier1_obach/compounds/omeprazole.yaml
name: omeprazole
smiles: "COc1ccc2[nH]c([S+]([O-])Cc3ncc(C)c(OC)c3C)nc2c1"
molecular_weight: 345.42
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
      method: "Austin 2002 estimate from logP=2.23"
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
      method: "Lopez-Rodriguez 2010 Xenobiotica 40:617 or Naritomi 2001 Drug Metab Dispos 29:1316"
    fm_cyp2c19: 0.8
    fm_cyp3a4: 0.2
  renal:
    clrenal_L_h:
      value: 0.0
      source: literature
      method: "Andersson 1996 (negligible urinary excretion of parent)"
observed_iv:
  dose_mg: 40.0
  cl_L_h: 34.2
  vss_L: 24.5
  t_half_h: 0.9
  source: "Andersson T, Clin Pharmacokinet 1996;31(1):9-28 (EM population mean)"
  phenotype_caveat: "EM population; CYP2C19 PMs ~3x lower CL (not modelled)"
```

- [ ] **Step 3: Verify YAML parses**

Run: `python3 -c "import yaml; yaml.safe_load(open('validation/data/tier1_obach/compounds/omeprazole.yaml'))"`
Expected: no error

- [ ] **Step 4: Commit**

```bash
git add validation/data/tier1_obach/compounds/omeprazole.yaml
git commit -m "data: add omeprazole compound YAML (literature overrides)"
```

---

## Task 4: Write dextromethorphan.yaml

**Files:**
- Create: `validation/data/tier1_obach/compounds/dextromethorphan.yaml`

**IMPORTANT**: Same verification protocol as Task 3. Use ChEMBL CHEMBL52440. Pay special attention to logP (literature range 2.7-3.4, NOT 4.11), pKa (8.3 vs 9.2), and CLint source (Bolaji 1993 is protein binding, not HLM CLint — find proper source).

- [ ] **Step 1: Verify SMILES and MW against ChEMBL**

Use `mcp__claude_ai_ChEMBL__compound_search(name="dextromethorphan", max_phase=4)` to confirm:
- SMILES: `COc1ccc2c(c1)[C@]13CCCC[C@@H]1[C@H](C2)N(C)CC3` (CHEMBL52440 canonical)
- InChIKey: `MKXZASYAUGDDCJ-NJAFHUGGSA-N`
- MW: 271.40

- [ ] **Step 2: Create dextromethorphan.yaml**

```yaml
# validation/data/tier1_obach/compounds/dextromethorphan.yaml
name: dextromethorphan
smiles: "COc1ccc2c(c1)[C@]13CCCC[C@@H]1[C@H](C2)N(C)CC3"
molecular_weight: 271.40
source: literature
properties:
  physicochemical:
    logp:
      value: 3.0
      source: literature
      method: "Literature consensus (range 2.7-3.4); ChEMBL ALogP=3.38"
    pka_base:
      value: 8.3
      source: literature
      method: "Literature (tertiary amine morphinan)"
    compound_type: base
  binding:
    fu_p:
      value: 0.68
      source: literature
      unit: fraction
      method: "Bolaji 1993 Br J Clin Pharmacol 36:131 (equilibrium dialysis)"
    fu_inc:
      value: 0.80
      source: literature
      unit: fraction
      method: "Austin 2002 estimate from logP~3.0; recalculate with verified logP"
    bp_ratio:
      value: 1.04
      source: literature
      unit: ratio
      method: "DrugBank DB00514"
  metabolism:
    clint_uL_min_mg:
      value: 85.0
      source: literature
      unit: uL/min/mg
      method: "Verify: Obach 1999 Table 2 or Houston group publications (HLM CLint)"
    fm_cyp2d6: 0.85
    fm_cyp3a4: 0.15
  renal:
    clrenal_L_h:
      value: 0.0
      source: literature
      method: "negligible urinary excretion of parent"
observed_iv:
  dose_mg: 30.0
  cl_L_h: 50.0
  vss_L: 250.0
  t_half_h: 3.5
  source: "Schadel 1995 J Clin Psychopharmacol 15:263 + Bolaji 1993 (EM population)"
  phenotype_caveat: "EM population; CYP2D6 PMs ~10x lower CL (not modelled)"
```

- [ ] **Step 3: Verify YAML parses**

Run: `python3 -c "import yaml; yaml.safe_load(open('validation/data/tier1_obach/compounds/dextromethorphan.yaml'))"`
Expected: no error

- [ ] **Step 4: Commit**

```bash
git add validation/data/tier1_obach/compounds/dextromethorphan.yaml
git commit -m "data: add dextromethorphan compound YAML (literature overrides)"
```

---

## Task 5: Extend panel.yaml with 2 new entries

**Files:**
- Modify: `validation/data/tier1_obach/panel.yaml`

- [ ] **Step 1: Read current panel.yaml to find insertion point**

The file has a `compounds:` list. Append two entries at the end.

- [ ] **Step 2: Append entries**

Add at the end of the `compounds:` list in `panel.yaml`:

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
```

- [ ] **Step 3: Verify YAML parses and entry count is correct**

Run: `python3 -c "import yaml; d = yaml.safe_load(open('validation/data/tier1_obach/panel.yaml')); print(len(d['compounds']), 'compounds')"`
Expected: `12 compounds`

- [ ] **Step 4: Commit**

```bash
git add validation/data/tier1_obach/panel.yaml
git commit -m "data: add omeprazole + dextromethorphan to Obach panel (n=10→12)"
```

---

## Task 6: Implement layer1_admet.py + unit tests

**Files:**
- Create: `validation/benchmarks/layer1_admet.py`
- Create: `tests/unit/test_layer1_benchmark.py`

This is the largest task. The script reads `data/validation/adme_reference.csv` (153 rows), calls `charon.predict.predict_properties()` for each, computes per-property metrics, and emits a report via `report_writer.emit_report()`.

- [ ] **Step 1: Write failing unit tests**

```python
# tests/unit/test_layer1_benchmark.py
"""Unit tests for the Layer 1 ADMET benchmark script."""
from __future__ import annotations

import csv
import tempfile
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]


@pytest.fixture
def mini_csv(tmp_path):
    """6-row mini CSV for testing benchmark logic."""
    rows = [
        {"name": "drug_a", "smiles": "CCO", "mw": "46.07",
         "logP": "0.31", "fup": "0.90", "rbp": "0.95",
         "clint_3a4_uL_min_pmol": "0.1", "peff_cm_s": "1.0e-4"},
        {"name": "drug_b", "smiles": "c1ccccc1", "mw": "78.11",
         "logP": "2.13", "fup": "0.50", "rbp": "0.80",
         "clint_3a4_uL_min_pmol": "5.0", "peff_cm_s": "2.0e-4"},
        {"name": "bad_smiles", "smiles": "INVALID_XYZ", "mw": "100",
         "logP": "1.0", "fup": "0.5", "rbp": "0.7",
         "clint_3a4_uL_min_pmol": "1.0", "peff_cm_s": "1.0e-4"},
        {"name": "zero_fup", "smiles": "CCO", "mw": "46.07",
         "logP": "0.31", "fup": "0.0", "rbp": "0.95",
         "clint_3a4_uL_min_pmol": "0.1", "peff_cm_s": "1.0e-4"},
    ]
    p = tmp_path / "mini_ref.csv"
    with open(p, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=rows[0].keys())
        w.writeheader()
        w.writerows(rows)
    return p


class TestLayer1Benchmark:
    def test_run_returns_payload_dict(self, mini_csv, tmp_path):
        import sys
        sys.path.insert(0, str(REPO_ROOT))
        from validation.benchmarks.layer1_admet import run_benchmark

        payload = run_benchmark(csv_path=mini_csv, reports_dir=tmp_path)
        assert isinstance(payload, dict)
        assert "title" in payload
        assert "summary" in payload
        assert "rows" in payload

    def test_invalid_smiles_excluded(self, mini_csv, tmp_path):
        import sys
        sys.path.insert(0, str(REPO_ROOT))
        from validation.benchmarks.layer1_admet import run_benchmark

        payload = run_benchmark(csv_path=mini_csv, reports_dir=tmp_path)
        excluded_names = [e["name"] for e in payload.get("excluded", [])]
        assert "bad_smiles" in excluded_names

    def test_zero_fup_excluded_from_fup_metric(self, mini_csv, tmp_path):
        import sys
        sys.path.insert(0, str(REPO_ROOT))
        from validation.benchmarks.layer1_admet import run_benchmark

        payload = run_benchmark(csv_path=mini_csv, reports_dir=tmp_path)
        fup_summary = payload["summary"].get("fu_p", {})
        if "n_valid" in fup_summary:
            assert fup_summary["n_valid"] < 4

    def test_logp_uses_mae_not_aafe(self, mini_csv, tmp_path):
        import sys
        sys.path.insert(0, str(REPO_ROOT))
        from validation.benchmarks.layer1_admet import run_benchmark

        payload = run_benchmark(csv_path=mini_csv, reports_dir=tmp_path)
        logp_summary = payload["summary"].get("logP", {})
        assert "mae" in logp_summary or "aafe" not in logp_summary
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/unit/test_layer1_benchmark.py -v`
Expected: FAIL (ImportError)

- [ ] **Step 3: Implement layer1_admet.py**

```python
# validation/benchmarks/layer1_admet.py
"""Layer 1 ADMET benchmark — predict_properties vs adme_reference.csv.

Computes per-property accuracy metrics:
  - logP: MAE, RMSE, R², within-0.5-log
  - fu_p, bp_ratio, peff: AAFE, within-2-fold, within-3-fold

CLint is EXCLUDED: adme_reference uses recombinant CYP3A4 (uL/min/pmol)
while Charon predicts hepatocyte CLint (uL/min/10^6 cells). Unit mismatch
makes comparison meaningless.

Run standalone:
    python3 validation/benchmarks/layer1_admet.py
    python3 validation/benchmarks/layer1_admet.py --limit 20
"""
from __future__ import annotations

import argparse
import csv
import math
import sys
from datetime import datetime, timezone
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from charon.predict import predict_properties  # noqa: E402
from validation.benchmarks.metrics import (  # noqa: E402
    aafe,
    mae,
    pearson_r,
    rmse,
    within_abs_diff,
    within_n_fold,
)
from validation.benchmarks.report_writer import emit_report  # noqa: E402

DEFAULT_CSV = REPO_ROOT / "data" / "validation" / "adme_reference.csv"
DEFAULT_REPORTS_DIR = REPO_ROOT / "validation" / "reports"

_DISCLAIMER = (
    "Evaluated on the conformal calibration set (n=153). Point-estimate "
    "metrics are meaningful; 90% CI coverage is tautological by construction "
    "(P90 was fit on the same set) and does NOT reflect true out-of-sample "
    "performance. CLint excluded: unit mismatch between Charon "
    "(hepatocyte uL/min/10^6) and reference (recombinant uL/min/pmol CYP3A4)."
)

_PROPERTY_CONFIG = [
    # (key, csv_col, charon_path, source_tag, metric_type, min_positive)
    ("logP", "logP", ("physicochemical", "logp"),
     "descriptor (RDKit Crippen)", "linear", None),
    ("fu_p", "fup", ("binding", "fu_p"),
     "ML model (XGBoost)", "log", 1e-4),
    ("bp_ratio", "rbp", ("binding", "bp_ratio"),
     "empirical formula", "log", 0.1),
    ("peff_cm_s", "peff_cm_s", ("permeability", "peff_cm_s"),
     "derived", "log", 1e-8),
]


def _get_predicted_value(props, path: tuple[str, ...]) -> float | None:
    obj = props
    for attr in path:
        obj = getattr(obj, attr, None)
        if obj is None:
            return None
    if hasattr(obj, "value"):
        return obj.value
    return float(obj) if obj is not None else None


def run_benchmark(
    *,
    csv_path: Path = DEFAULT_CSV,
    reports_dir: Path = DEFAULT_REPORTS_DIR,
    limit: int | None = None,
) -> dict:
    """Run Layer 1 benchmark and return payload dict."""
    rows_out: list[dict] = []
    excluded: list[dict] = []
    per_prop: dict[str, dict[str, list[float]]] = {
        cfg[0]: {"pred": [], "obs": []} for cfg in _PROPERTY_CONFIG
    }

    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        all_rows = list(reader)

    if limit:
        all_rows = all_rows[:limit]

    for row in all_rows:
        name = row["name"]
        smiles = row["smiles"]

        try:
            props = predict_properties(smiles)
        except ValueError:
            excluded.append({"name": name, "reason": "SMILES invalid"})
            continue
        except Exception as e:
            excluded.append({
                "name": name,
                "reason": f"prediction failed: {type(e).__name__}",
            })
            continue

        for key, csv_col, charon_path, source_tag, metric_type, min_pos in _PROPERTY_CONFIG:
            obs_raw = row.get(csv_col)
            if obs_raw is None or obs_raw == "":
                continue
            try:
                obs = float(obs_raw)
            except ValueError:
                continue

            pred = _get_predicted_value(props, charon_path)
            if pred is None:
                excluded.append({
                    "name": name, "property": key, "reason": "prediction missing",
                })
                continue

            if min_pos is not None and (obs <= min_pos or pred <= min_pos):
                excluded.append({
                    "name": name, "property": key,
                    "reason": f"non-positive value (obs={obs}, pred={pred})",
                })
                continue

            per_prop[key]["pred"].append(pred)
            per_prop[key]["obs"].append(obs)

            if metric_type == "linear":
                rows_out.append({
                    "name": name, "property": key,
                    "pred": round(pred, 4), "obs": round(obs, 4),
                    "abs_diff": round(abs(pred - obs), 4),
                })
            else:
                fe = max(pred / obs, obs / pred) if (pred > 0 and obs > 0) else None
                rows_out.append({
                    "name": name, "property": key,
                    "pred": round(pred, 4), "obs": round(obs, 4),
                    "fold_error": round(fe, 2) if fe else None,
                })

    summary: dict[str, dict] = {}
    for key, csv_col, charon_path, source_tag, metric_type, min_pos in _PROPERTY_CONFIG:
        preds = per_prop[key]["pred"]
        obs = per_prop[key]["obs"]
        n_valid = len(preds)
        entry: dict = {
            "source_tag": source_tag,
            "n": len(all_rows),
            "n_valid": n_valid,
        }
        if n_valid >= 2:
            if metric_type == "linear":
                entry["mae"] = round(mae(preds, obs), 4)
                entry["rmse"] = round(rmse(preds, obs), 4)
                entry["r2"] = round(pearson_r(preds, obs) ** 2, 4)
                entry["within_0.5_log_pct"] = round(
                    within_abs_diff(preds, obs, threshold=0.5) * 100, 1
                )
            else:
                entry["aafe"] = round(aafe(preds, obs), 4)
                entry["within_2fold_pct"] = round(
                    within_n_fold(preds, obs, n=2.0) * 100, 1
                )
                entry["within_3fold_pct"] = round(
                    within_n_fold(preds, obs, n=3.0) * 100, 1
                )
        summary[key] = entry

    targets: dict[str, dict] = {}
    for key, csv_col, charon_path, source_tag, metric_type, min_pos in _PROPERTY_CONFIG:
        s = summary.get(key, {})
        if metric_type == "linear":
            val = s.get("mae")
            targets[key] = {
                "metric": "mae", "target": 0.5,
                "met": val is not None and val < 0.5,
            }
        else:
            val = s.get("aafe")
            targets[key] = {
                "metric": "aafe", "target": 2.0,
                "met": val is not None and val < 2.0,
            }

    payload = {
        "title": "Charon Layer 1 ADMET Benchmark",
        "panel": f"adme_reference.csv (n={len(all_rows)})",
        "date_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "disclaimer": _DISCLAIMER,
        "summary": summary,
        "rows": rows_out,
        "excluded": excluded,
        "targets": targets,
        "notes": [
            f"{len(excluded)} exclusions (see Excluded section)",
            "CLint not evaluated (unit mismatch: hepatocyte vs recombinant CYP3A4)",
            "logP MAE target <0.5 is a sprint-local convention, not from CLAUDE.md §8",
        ],
    }

    emit_report(payload, stem=reports_dir / "layer1_admet")

    # Print summary to stdout
    print(f"Layer 1 ADMET Benchmark — {len(all_rows)} compounds")
    print("-" * 60)
    for key, s in summary.items():
        tag = s.get("source_tag", "")
        n_valid = s.get("n_valid", 0)
        if "mae" in s:
            print(f"  {key:12s} ({tag}): MAE={s['mae']:.3f}  "
                  f"RMSE={s.get('rmse', 0):.3f}  R²={s.get('r2', 0):.3f}  "
                  f"n={n_valid}")
        elif "aafe" in s:
            print(f"  {key:12s} ({tag}): AAFE={s['aafe']:.2f}  "
                  f"2x={s.get('within_2fold_pct', 0):.0f}%  "
                  f"3x={s.get('within_3fold_pct', 0):.0f}%  n={n_valid}")
        else:
            print(f"  {key:12s}: insufficient data (n={n_valid})")
    print("-" * 60)
    print(f"Exclusions: {len(excluded)}")

    return payload


def main() -> int:
    parser = argparse.ArgumentParser(description="Layer 1 ADMET benchmark")
    parser.add_argument("--csv", type=Path, default=DEFAULT_CSV)
    parser.add_argument("--limit", type=int, default=None)
    args = parser.parse_args()
    try:
        run_benchmark(csv_path=args.csv, limit=args.limit)
    except (FileNotFoundError, PermissionError) as e:
        print(f"error: {e}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
```

- [ ] **Step 4: Run unit tests**

Run: `pytest tests/unit/test_layer1_benchmark.py -v`
Expected: 4 passed

- [ ] **Step 5: Run full test suite**

Run: `pytest tests/ -x -q`
Expected: all pass

- [ ] **Step 6: Commit**

```bash
git add validation/benchmarks/layer1_admet.py tests/unit/test_layer1_benchmark.py
git commit -m "feat(validation): add Layer 1 ADMET benchmark (adme_reference.csv)"
```

---

## Task 7: Modify layer2_human_pk.py to emit reports

**Files:**
- Modify: `validation/benchmarks/layer2_human_pk.py`

- [ ] **Step 1: Read current file and identify insertion points**

Key: `main()` function starts at line 320, ends at line 351. Insert report emit after the summary print and before the exit.

- [ ] **Step 2: Add imports at top of file (after existing imports)**

Add after the existing import block (around line 25):

```python
from datetime import datetime, timezone
from validation.benchmarks.report_writer import emit_report

REPORTS_DIR = REPO_ROOT / "validation" / "reports"
```

- [ ] **Step 3: Add payload builder function before main()**

Insert before `def main()`:

```python
def _build_report_payload(
    summaries: dict[str, PanelSummary],
    rows: dict[str, list[PanelRow]],
) -> dict:
    """Convert benchmark results to report_writer payload."""
    summary = summaries["with_override"]
    per_rows = []
    for r in rows["with_override"]:
        per_rows.append({
            "compound": r.key,
            "CL_pred": round(r.predicted["cl_L_h"], 3),
            "CL_obs": round(r.observed["cl_L_h"], 3),
            "fold_CL": round(r.fold["cl_L_h"], 2),
            "Vss_pred": round(r.predicted["vss_L"], 2),
            "Vss_obs": round(r.observed["vss_L"], 2),
            "fold_Vss": round(r.fold["vss_L"], 2),
            "t_pred": round(r.predicted["t_half_h"], 2),
            "t_obs": round(r.observed["t_half_h"], 2),
            "fold_t": round(r.fold["t_half_h"], 2),
            "strict": r.strict_targets,
        })
    return {
        "title": "Charon Layer 2 Human PBPK Benchmark",
        "panel": f"Obach 1999 Tier-1 Panel (n={summary.n})",
        "date_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "summary": {
            "cl_L_h": {
                "aafe": round(summary.aafe["cl_L_h"], 2),
                "within_2fold_pct": round(summary.within_2_fold["cl_L_h"] * 100, 0),
                "within_3fold_pct": round(summary.within_3_fold["cl_L_h"] * 100, 0),
            },
            "vss_L": {
                "aafe": round(summary.aafe["vss_L"], 2),
                "within_2fold_pct": round(summary.within_2_fold["vss_L"] * 100, 0),
                "within_3fold_pct": round(summary.within_3_fold["vss_L"] * 100, 0),
            },
            "t_half_h": {
                "aafe": round(summary.aafe["t_half_h"], 2),
                "within_2fold_pct": round(summary.within_2_fold["t_half_h"] * 100, 0),
                "within_3fold_pct": round(summary.within_3_fold["t_half_h"] * 100, 0),
            },
        },
        "targets": {
            "cl_L_h": {"metric": "aafe", "target": 2.5,
                       "met": summary.aafe["cl_L_h"] < 2.5},
            "vss_L": {"metric": "aafe", "target": 3.0,
                      "met": summary.aafe["vss_L"] < 3.0},
        },
        "rows": per_rows,
        "notes": [
            f"Strict gate failures: {summary.strict_failures}",
            "Only theophylline/antipyrine are strict_targets",
        ],
    }
```

- [ ] **Step 4: Add emit call in main(), after the overrides print and before the exit**

In `main()`, insert after the `print("=" * 100)` / "All strict-targets..." block and before `return`:

```python
    # Emit structured report
    payload = _build_report_payload(summaries, rows)
    emit_report(payload, stem=REPORTS_DIR / "layer2_human_pk")
```

- [ ] **Step 5: Run existing tests to verify no regressions**

Run: `pytest tests/ -x -q`
Expected: all pass

- [ ] **Step 6: Commit**

```bash
git add validation/benchmarks/layer2_human_pk.py
git commit -m "feat(validation): add report_writer emit to Layer 2 benchmark"
```

---

## Task 8: Add docstring stub to layer3_fih_dose.py

**Files:**
- Modify: `validation/benchmarks/layer3_fih_dose.py`

- [ ] **Step 1: Write the docstring**

Replace the empty file contents with:

```python
"""Layer 3 FIH dose benchmark — DEFERRED.

This benchmark script is intentionally empty. FIH dose validation requires
a dataset of compounds with known published MRSD (Maximum Recommended
Starting Dose) from IND/NDA filings — the ``tier2_drugs_at_fda`` dataset.

That dataset has not been built. The directories
``validation/data/tier2_drugs_at_fda/`` and ``validation/data/tier3_literature/``
exist as scaffolds but are currently empty.

**Why not use the Obach panel?** The Obach 1999 panel provides IV PK
parameters (CL, Vss, t½) for Layer 2 validation. These old drugs predate
modern FIH dose-finding methodology (ICH M3, FDA 2005 guidance). Their
clinical starting doses are not documented in the MRSD sense.

**When will this be implemented?** After the tier2_drugs_at_fda dataset is
curated from publicly available FDA NDA reviews (Drugs@FDA database) or
published Phase I dose-escalation papers.

See also: README.md § Known Limitations, § Validation Status.
"""
```

- [ ] **Step 2: Commit**

```bash
git add validation/benchmarks/layer3_fih_dose.py
git commit -m "docs: add deferral rationale to empty layer3_fih_dose.py"
```

---

## Task 9: Run benchmarks and commit reports

**Files:**
- Create: `validation/reports/layer1_admet.md`
- Create: `validation/reports/layer1_admet.json`
- Create: `validation/reports/layer2_human_pk.md`
- Create: `validation/reports/layer2_human_pk.json`

This is a run-and-commit task, not a code task.

- [ ] **Step 1: Create reports directory**

```bash
mkdir -p validation/reports
```

- [ ] **Step 2: Run Layer 1 benchmark**

```bash
python3 validation/benchmarks/layer1_admet.py
```

Expected: exits 0, prints summary table, creates `validation/reports/layer1_admet.{md,json}`.

**Record the actual AAFE/MAE numbers** — they'll be needed for the README (Task 13).

- [ ] **Step 3: Run Layer 2 benchmark**

```bash
python3 validation/benchmarks/layer2_human_pk.py
```

Expected: exits 0 (strict gate passes), creates `validation/reports/layer2_human_pk.{md,json}`.

**Record the actual AAFE numbers** from the with-overrides summary.

- [ ] **Step 4: Verify report files exist**

```bash
ls -la validation/reports/
```

Expected: 4 files (layer1_admet.md, .json, layer2_human_pk.md, .json).

- [ ] **Step 5: Commit**

```bash
git add validation/reports/
git commit -m "chore(validation): commit Layer 1 + Layer 2 benchmark reports"
```

---

## Task 10: Implement regression tests

**Files:**
- Create: `tests/regression/test_known_drugs.py`

- [ ] **Step 1: Write the test file**

```python
# tests/regression/test_known_drugs.py
"""Regression guards for 5 CLAUDE.md §8 reference drugs.

Two tests per drug:
  - test_pipeline_runs: smoke — pipeline runs, outputs finite/positive
  - test_dose_within_snapshot: MRSD (or AUC fallback) within ±25% of baseline

Golden baselines in golden_outputs/*.json. Regenerate with:
    UPDATE_BASELINES=1 pytest tests/regression/test_known_drugs.py -v

NOAEL values below are REGRESSION FIXTURES, not literature-verified
safety data. They produce a dose recommendation for snapshot comparison
only. Do NOT use these values for actual dose calculations.

Timing target: < 30 seconds total (no uncertainty quantification).
"""
from __future__ import annotations

import json
import math
import os
from dataclasses import dataclass
from pathlib import Path

import pytest

from charon.core.schema import DoseProjectionConfig
from charon.pipeline import Pipeline


GOLDEN_DIR = Path(__file__).parent / "golden_outputs"
TOLERANCE = 0.25


@dataclass(frozen=True)
class RefDrug:
    key: str
    smiles: str
    route: str
    dose_mg: float
    noael_mg_kg: float | None
    noael_species: str | None
    target_kd_nM: float | None
    target_ceff_nM: float | None
    phenotype_caveat: bool


REFERENCE_DRUGS = (
    RefDrug(
        "midazolam",
        "Cc1ncc2n1-c1ccc(Cl)cc1C(c1ccccc1F)=NC2",
        route="iv_bolus", dose_mg=5.0,
        noael_mg_kg=2.0, noael_species="rat",
        target_kd_nM=None, target_ceff_nM=None,
        phenotype_caveat=False,
    ),
    RefDrug(
        "warfarin",
        "CC(=O)CC(c1ccccc1)c1c(O)c2ccccc2oc1=O",
        route="iv_bolus", dose_mg=5.0,
        noael_mg_kg=0.5, noael_species="rat",
        target_kd_nM=None, target_ceff_nM=None,
        phenotype_caveat=False,
    ),
    RefDrug(
        "diclofenac",
        "OC(=O)Cc1ccccc1Nc1c(Cl)cccc1Cl",
        route="iv_bolus", dose_mg=75.0,
        noael_mg_kg=5.0, noael_species="rat",
        target_kd_nM=None, target_ceff_nM=None,
        phenotype_caveat=False,
    ),
    RefDrug(
        "omeprazole",
        "COc1ccc2[nH]c([S+]([O-])Cc3ncc(C)c(OC)c3C)nc2c1",
        route="iv_bolus", dose_mg=40.0,
        noael_mg_kg=4.0, noael_species="rat",
        target_kd_nM=None, target_ceff_nM=None,
        phenotype_caveat=True,
    ),
    RefDrug(
        "dextromethorphan",
        "COc1ccc2c(c1)[C@]13CCCC[C@@H]1[C@H](C2)N(C)CC3",
        route="iv_bolus", dose_mg=30.0,
        noael_mg_kg=10.0, noael_species="rat",
        target_kd_nM=None, target_ceff_nM=None,
        phenotype_caveat=True,
    ),
)


def _run_pipeline(drug: RefDrug):
    dp = None
    if drug.noael_mg_kg is not None and drug.noael_species is not None:
        dp = DoseProjectionConfig(
            noael_mg_kg=drug.noael_mg_kg,
            noael_species=drug.noael_species,
        )
    elif drug.target_kd_nM is not None:
        dp = DoseProjectionConfig(target_kd_nM=drug.target_kd_nM)
    elif drug.target_ceff_nM is not None:
        dp = DoseProjectionConfig(target_ceff_nM=drug.target_ceff_nM)

    return Pipeline.from_smiles(
        drug.smiles,
        route=drug.route,
        dose_mg=drug.dose_mg,
        dose_projection=dp,
        compound_name=drug.key,
    ).run()


def _extract_primary(result) -> dict:
    rec = result.dose_recommendation
    if rec is not None and hasattr(rec, "mrsd_mg") and rec.mrsd_mg > 0:
        return {"name": "mrsd_mg", "value": float(rec.mrsd_mg)}
    pk = result.pk_parameters
    if pk.auc_0_inf is not None and pk.auc_0_inf > 0:
        return {"name": "auc_0_inf", "value": float(pk.auc_0_inf)}
    return {"name": "cmax", "value": float(pk.cmax)}


def _write_baseline(path: Path, primary: dict, drug: RefDrug, result):
    pk = result.pk_parameters
    data = {
        "drug": drug.key,
        "smiles": drug.smiles,
        "charon_version": "0.1.0",
        "primary_metric": primary,
        "secondary_metrics": {
            "cmax_ug_L": float(pk.cmax) if pk.cmax else None,
            "auc_ug_h_L": float(pk.auc_0_inf) if pk.auc_0_inf else None,
            "t_half_h": float(pk.half_life) if pk.half_life else None,
        },
        "generated_utc": __import__("datetime").datetime.now(
            __import__("datetime").timezone.utc
        ).isoformat(timespec="seconds"),
        "caveat": (
            "CYP polymorphism; baseline may drift"
            if drug.phenotype_caveat else None
        ),
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, default=str), encoding="utf-8")


@pytest.mark.parametrize("drug", REFERENCE_DRUGS, ids=lambda d: d.key)
def test_pipeline_runs(drug):
    """Pipeline executes without exception; outputs are finite and positive."""
    result = _run_pipeline(drug)
    assert result.pk_parameters.cmax is not None
    assert result.pk_parameters.cmax > 0
    assert math.isfinite(result.pk_parameters.cmax)
    if result.pk_parameters.auc_0_inf is not None:
        assert result.pk_parameters.auc_0_inf > 0


@pytest.mark.parametrize("drug", REFERENCE_DRUGS, ids=lambda d: d.key)
def test_dose_within_snapshot(drug):
    """Primary metric within ±25% of committed golden baseline."""
    result = _run_pipeline(drug)
    primary = _extract_primary(result)
    baseline_path = GOLDEN_DIR / f"{drug.key}_baseline.json"

    if os.environ.get("UPDATE_BASELINES") == "1":
        _write_baseline(baseline_path, primary, drug, result)
        pytest.skip(f"Updated baseline for {drug.key}")

    if not baseline_path.exists():
        pytest.fail(
            f"Golden baseline missing: {baseline_path}. "
            f"Run with UPDATE_BASELINES=1 to create."
        )

    baseline = json.loads(baseline_path.read_text())
    base_val = baseline["primary_metric"]["value"]
    base_name = baseline["primary_metric"]["name"]
    curr_val = primary["value"]

    assert primary["name"] == base_name, (
        f"Primary metric changed: was {base_name}, now {primary['name']}"
    )

    if base_val == 0:
        pytest.skip("Baseline value is zero; cannot compute fold error")

    fe = max(curr_val / base_val, base_val / curr_val)
    assert fe <= (1.0 + TOLERANCE), (
        f"{drug.key} {base_name}: current={curr_val:.4g}, "
        f"baseline={base_val:.4g}, fold_error={fe:.3f} > {1.0 + TOLERANCE}"
    )
```

- [ ] **Step 2: Run smoke tests only (snapshot tests will fail until baselines exist)**

Run: `pytest tests/regression/test_known_drugs.py::test_pipeline_runs -v`
Expected: 5 passed (pipelines all run)

- [ ] **Step 3: Commit**

```bash
git add tests/regression/test_known_drugs.py
git commit -m "test(regression): add 5-drug regression guard (smoke + snapshot)"
```

---

## Task 11: Generate golden baselines

**Files:**
- Create: `tests/regression/golden_outputs/midazolam_baseline.json`
- Create: `tests/regression/golden_outputs/warfarin_baseline.json`
- Create: `tests/regression/golden_outputs/diclofenac_baseline.json`
- Create: `tests/regression/golden_outputs/omeprazole_baseline.json`
- Create: `tests/regression/golden_outputs/dextromethorphan_baseline.json`

- [ ] **Step 1: Generate baselines**

```bash
UPDATE_BASELINES=1 pytest tests/regression/test_known_drugs.py::test_dose_within_snapshot -v
```

Expected: 5 skipped (baselines written).

- [ ] **Step 2: Verify files exist**

```bash
ls tests/regression/golden_outputs/
```

Expected: 5 JSON files.

- [ ] **Step 3: Inspect one baseline for sanity**

```bash
cat tests/regression/golden_outputs/midazolam_baseline.json
```

Verify: `primary_metric.value` is positive, `charon_version` is `0.1.0`, `caveat` is null for midazolam.

- [ ] **Step 4: Commit**

```bash
git add tests/regression/golden_outputs/
git commit -m "chore(regression): commit initial golden baselines for 5 reference drugs"
```

---

## Task 12: Verify regression suite passes

- [ ] **Step 1: Run full regression suite**

```bash
pytest tests/regression/test_known_drugs.py -v
```

Expected: 10 passed (5 smoke + 5 snapshot).

- [ ] **Step 2: Run full test suite**

```bash
pytest tests/ -x -q
```

Expected: all pass (existing ~760 + new ~20).

---

## Task 13: Write README.md

**Files:**
- Create: `README.md`

This task requires the **actual numbers from Task 9 reports**. Read `validation/reports/layer1_admet.json` and `validation/reports/layer2_human_pk.json` to fill in the TBD fields.

- [ ] **Step 1: Read report JSON files for actual numbers**

```bash
cat validation/reports/layer1_admet.json | python3 -c "import json,sys; d=json.load(sys.stdin); [print(f'{k}: {v}') for k,v in d['summary'].items()]"
cat validation/reports/layer2_human_pk.json | python3 -c "import json,sys; d=json.load(sys.stdin); [print(f'{k}: {v}') for k,v in d['summary'].items()]"
```

- [ ] **Step 2: Write README.md**

Write the full README using the structure from spec §11. The content MUST use real numbers from the reports — no TBD placeholders.

The README should follow these structure sections (spec §11.1):
1. Charon description (4 lines)
2. Project status (Phase A MVP)
3. **Known limitations** (prominent, 12+ bullets)
4. Install (`pip install -e .` / `pip install -e ".[dev]"`)
5. Quickstart (5 CLI examples + 1 Python API example)
6. Reference drug quick-check (`pytest tests/regression/...`)
7. Validation status (Layer 1 table, Layer 2 table, Layer 3 deferred)
8. Architecture overview (ASCII diagram + link to ARCHITECTURE.md)
9. Development (test/benchmark/baseline commands)
10. Roadmap (Phase B, Phase C)
11. Citing Charon
12. License (MIT)
13. Acknowledgments

**Tone rules** (spec §11.2): no marketing language, concrete numbers always, failed gates shown as `[FAIL]`.

See spec §11 for full content guide. README body is too large to include inline here — the implementer should follow the spec section by section and substitute actual benchmark numbers.

- [ ] **Step 3: Verify all relative links**

```bash
grep -oP '\[.*?\]\(([^)]+)\)' README.md | while read link; do
  target=$(echo "$link" | grep -oP '\(([^)]+)\)' | tr -d '()')
  [ -f "$target" ] && echo "OK: $target" || echo "BROKEN: $target"
done
```

- [ ] **Step 4: Commit**

```bash
git add README.md
git commit -m "docs: add MVP-honest README with actual validation numbers"
```

---

## Task 14: Verify .gitignore

**Files:**
- Potentially modify: `.gitignore`

- [ ] **Step 1: Check that validation/reports/ is not ignored**

```bash
git check-ignore validation/reports/layer1_admet.md
git check-ignore tests/regression/golden_outputs/midazolam_baseline.json
```

Expected: no output (files are NOT ignored).

The existing `.gitignore` has `/reports/` (root-level only) and `*.report.md` / `*.report.json`. Our files are at `validation/reports/layer1_admet.md` (not root, doesn't end in `.report.md`) — should be fine.

- [ ] **Step 2: If any file is ignored, add explicit negation**

Only if step 1 shows files are ignored:
```
!validation/reports/
!tests/regression/golden_outputs/
```

- [ ] **Step 3: Commit if changed**

```bash
git add .gitignore
git commit -m "chore: ensure validation/reports and golden_outputs are tracked"
```

---

## Task 15: Final sweep

- [ ] **Step 1: Search for TBD/FIXME/TODO in all new files**

```bash
grep -rn "TBD\|FIXME\|TODO" README.md validation/reports/ tests/regression/ validation/benchmarks/report_writer.py validation/benchmarks/layer1_admet.py
```

Expected: no matches (or only intentional "TODO" in spec-referenced comments).

- [ ] **Step 2: Verify relative links in README**

- [ ] **Step 3: Commit any fixes**

---

## Task 16: Full test suite

- [ ] **Step 1: Run complete test suite**

```bash
pytest tests/ -v --tb=short 2>&1 | tail -30
```

Expected: ~780+ passed, 0 failed. Note exact count.

- [ ] **Step 2: Run benchmarks once more to confirm exit 0**

```bash
python3 validation/benchmarks/layer1_admet.py && echo "Layer 1: OK"
python3 validation/benchmarks/layer2_human_pk.py && echo "Layer 2: OK"
```

---

## Task 17: Push

- [ ] **Step 1: Check commit count and status**

```bash
git log --oneline origin/main..HEAD
git status
```

- [ ] **Step 2: Push**

```bash
git push origin main
```

Expected: all commits pushed successfully.
