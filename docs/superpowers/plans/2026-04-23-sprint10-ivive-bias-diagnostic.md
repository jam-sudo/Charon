# Sprint 10 — IVIVE bias diagnostic Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Produce a diagnostic report decomposing each Tier A fold-error into `fold_liver_model × fold_route_bias × fold_residual`, without changing production code.

**Architecture:** Pure-function library (`translational/decomposition.py`) + curated data file (`bioavailability.csv`) + thin orchestrator (`validation/benchmarks/layer3_ivive_decomposition.py`) that drives `Pipeline(..., liver_model=...)` three times per compound and emits `{md, json}` report via existing `emit_report`.

**Tech Stack:** Python 3.11+, pytest, PyYAML, dataclasses (stdlib). No new dependencies. Existing `Pipeline` already parameterises `liver_model`, so no production code changes are required.

**Spec:** `docs/superpowers/specs/2026-04-23-sprint10-ivive-bias-diagnostic-design.md`

---

## File Structure

**New files (all under `feature/sprint10-ivive-bias-diagnostic` branch):**

| Path | Responsibility |
|---|---|
| `src/charon/translational/decomposition.py` | Pure functions: `decompose_fold_error`, `select_best_alternate_liver_model`, `compute_route_bias_factor` + `DecompositionResult` dataclass. No I/O. |
| `validation/data/fih_reference/bioavailability.csv` | Curated literature-F table: 12 Tier A compounds × {f_oral, fih_reference_route, f_source, f_doi_or_pmid, notes}. |
| `validation/benchmarks/layer3_ivive_decomposition.py` | Orchestrator. Runs `Pipeline` × 3 liver models per compound, calls `decompose_fold_error`, writes `{stem}.md` + `{stem}.json`. |
| `validation/reports/layer3_ivive_decomposition.md` | Committed output (generated, checked-in so it's diffable). |
| `validation/reports/layer3_ivive_decomposition.json` | Committed output. |
| `tests/unit/test_decomposition.py` | 7 unit tests (pure-function math). |
| `tests/integration/test_layer3_decomposition_benchmark.py` | 3 integration tests (end-to-end orchestrator). |
| `docs/superpowers/sprint10-followup-lisinopril-ceff.md` | Follow-up ticket (if lisinopril literature check recommends a change). |

**Files NOT touched (reaffirming spec §3 / §10):** `src/charon/core/liver_models.py`, `src/charon/translational/pad.py`, `src/charon/translational/dose_projector.py`, `validation/data/fih_reference/panel.yaml`, everything under `src/charon/predict/`, `src/charon/pbpk/`, `src/charon/uncertainty/`.

**File-path convention note:** Existing tests live flat under `tests/unit/` (no `translational/` subdir despite spec §7.1's path). Plan uses `tests/unit/test_decomposition.py` to follow the established pattern.

---

## Task 1: Curate `bioavailability.csv` + schema test

**Files:**
- Create: `validation/data/fih_reference/bioavailability.csv`
- Create: `tests/unit/test_bioavailability_csv.py`

- [ ] **Step 1: Write the failing schema test**

```python
# tests/unit/test_bioavailability_csv.py
"""Schema & coverage tests for the Tier A bioavailability CSV."""

from __future__ import annotations

import csv
from pathlib import Path

import pytest
import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]
CSV_PATH = REPO_ROOT / "validation" / "data" / "fih_reference" / "bioavailability.csv"
PANEL_PATH = REPO_ROOT / "validation" / "data" / "fih_reference" / "panel.yaml"

REQUIRED_COLS = {
    "compound",
    "fih_reference_route",
    "f_oral",
    "f_source",
    "f_doi_or_pmid",
    "notes",
}


def _load_rows() -> list[dict]:
    with CSV_PATH.open() as fh:
        return list(csv.DictReader(fh))


def test_csv_exists():
    assert CSV_PATH.exists(), f"missing {CSV_PATH}"


def test_csv_columns_match_schema():
    rows = _load_rows()
    assert rows, "CSV has no rows"
    assert set(rows[0].keys()) == REQUIRED_COLS


def test_csv_covers_all_tier_a_compounds():
    panel = yaml.safe_load(PANEL_PATH.read_text())["panel"]
    tier_a = {c["name"] for c in panel["compounds"] if c["tier"] == "gold"}
    csv_compounds = {r["compound"] for r in _load_rows()}
    missing = tier_a - csv_compounds
    assert not missing, f"Tier A compounds missing from bioavailability.csv: {sorted(missing)}"


def test_csv_route_values_constrained():
    allowed = {"oral", "iv"}
    for row in _load_rows():
        assert row["fih_reference_route"] in allowed, (
            f"{row['compound']}: unexpected route "
            f"{row['fih_reference_route']!r}; allowed={sorted(allowed)}"
        )


def test_csv_f_oral_numeric_or_blank_iv():
    for row in _load_rows():
        f_raw = row["f_oral"].strip()
        if row["fih_reference_route"] == "iv":
            assert f_raw == "", (
                f"{row['compound']}: IV-reference rows must have blank f_oral, "
                f"got {f_raw!r}"
            )
        else:
            value = float(f_raw)
            assert 0.0 < value <= 1.0, (
                f"{row['compound']}: f_oral must be in (0, 1], got {value}"
            )


def test_csv_every_oral_row_has_citation():
    for row in _load_rows():
        if row["fih_reference_route"] == "oral":
            assert row["f_source"].strip(), f"{row['compound']}: f_source blank"
            assert row["f_doi_or_pmid"].strip(), (
                f"{row['compound']}: f_doi_or_pmid blank"
            )
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
pytest tests/unit/test_bioavailability_csv.py -v
```

Expected: all six tests FAIL with `FileNotFoundError` on missing CSV.

- [ ] **Step 3: Curate the CSV**

Cross-check each value against `panel.yaml`'s `source` field to determine whether the clinical reference dose was IV or oral.

Panel sources (from `validation/data/fih_reference/panel.yaml`):
- midazolam: "single-dose IV Phase 1" → IV
- warfarin: "Coumadin label initial dose" → oral (Coumadin is oral tablet)
- propranolol: label start → oral
- verapamil: label start → oral
- omeprazole: label start → oral
- theophylline: label start → oral
- diclofenac: label start → oral
- diazepam: FIH or label start → need to check (assume oral; diazepam IV 5 mg is a separate indication)
- metoprolol: label start → oral
- acetaminophen: label start → oral
- lisinopril: label start → oral
- atorvastatin: label start → oral

Write `validation/data/fih_reference/bioavailability.csv`:

```csv
compound,fih_reference_route,f_oral,f_source,f_doi_or_pmid,notes
midazolam,iv,,Heizmann 1983,PMID:6352420,IV reference dose; f_oral 0.30-0.60 range reported in oral studies but not applicable to IV FIH
warfarin,oral,1.0,Breckenridge 1970,PMID:5441976,Racemate near-complete absorption; F ~100% cited
propranolol,oral,0.26,Wood 1978,PMID:565544,Extensive first-pass; F 15-40% range; 0.26 = median
verapamil,oral,0.22,Eichelbaum 1981,PMID:7296790,High first-pass; F 10-35% range; 0.22 = median single-dose
omeprazole,oral,0.40,Andersson 1996,PMID:8737132,Single-dose F 0.30-0.40; increases to ~0.70 on repeat dosing
theophylline,oral,1.0,Hendeles 1977,PMID:872579,Near-complete absorption; F ~96-100%
diclofenac,oral,0.54,Willis 1979,PMID:523773,F 0.54 after oral tablet; subject to enterohepatic recirculation
diazepam,oral,0.93,Greenblatt 1980,PMID:7389287,F 0.90-1.00 oral vs IV cross-over
metoprolol,oral,0.50,Regardh 1980,PMID:7002442,F 0.40-0.60 range; 0.50 = midpoint of extensive metabolizers
acetaminophen,oral,0.88,Rawlins 1977,PMID:860545,F 0.80-0.90 oral; some first-pass sulfation
lisinopril,oral,0.25,Beermann 1988,PMID:3048950,Low GI absorption; F 0.25-0.29 reported range
atorvastatin,oral,0.14,Lennernas 2003,PMID:12965954,F ~0.14 reflects OATP1B1 hepatic uptake + CYP3A4 first-pass
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
pytest tests/unit/test_bioavailability_csv.py -v
```

Expected: 6 PASS.

- [ ] **Step 5: Commit**

```bash
git add validation/data/fih_reference/bioavailability.csv tests/unit/test_bioavailability_csv.py
git commit -m "data(sprint10): curated literature F table for Tier A (bioavailability.csv)"
```

---

## Task 2: `decomposition.py` library + unit tests (TDD)

**Files:**
- Create: `src/charon/translational/decomposition.py`
- Create: `tests/unit/test_decomposition.py`

**Design note — SIGNED factors throughout.** The log-additive invariant

```
log10(fold_observed_signed) = log10(fold_liver_model_signed)
                            + log10(fold_route_bias)
                            + log10(fold_residual_signed)
```

only holds when factors are *signed ratios* (can be < 1 meaning "under-predicted" or > 1 meaning "over-predicted"), not symmetric `max(x, 1/x)` forms. `DecompositionResult` therefore exposes signed values. A small helper `to_symmetric(signed)` is provided for any reporting that prefers the `max(x, 1/x)` convention (the benchmark `fold_error` metric is symmetric).

- [ ] **Step 1: Write the failing unit tests**

```python
# tests/unit/test_decomposition.py
"""Unit tests for charon.translational.decomposition.

All factors are signed (can be < 1 or > 1). Log-additive invariant:

    log10(fold_observed_signed) = log10(fold_liver_model_signed)
                                + log10(fold_route_bias)
                                + log10(fold_residual_signed)
"""

from __future__ import annotations

import math

import pytest

from charon.translational.decomposition import (
    DecompositionResult,
    compute_route_bias_factor,
    decompose_fold_error,
    select_best_alternate_liver_model,
    to_symmetric,
)


# ---------------------------------------------------------------------------
# to_symmetric helper
# ---------------------------------------------------------------------------

def test_to_symmetric_above_one():
    assert to_symmetric(2.0) == 2.0


def test_to_symmetric_below_one():
    assert to_symmetric(0.5) == pytest.approx(2.0, rel=1e-12)


def test_to_symmetric_exactly_one():
    assert to_symmetric(1.0) == 1.0


# ---------------------------------------------------------------------------
# select_best_alternate_liver_model
# ---------------------------------------------------------------------------

def test_select_best_alternate_picks_closest_to_reference():
    # ws=10, pt=5, disp=15; ref=5 -> pt is closest (|log10(5/5)| = 0)
    name, mrsd = select_best_alternate_liver_model(
        ws=10.0, pt=5.0, disp=15.0, reference=5.0
    )
    assert name == "parallel_tube"
    assert mrsd == 5.0


def test_select_best_alternate_picks_dispersion_when_closer():
    # ws=2, pt=1, disp=4; ref=5 -> |log10(disp/ref)|≈0.097, |log10(pt/ref)|≈0.699
    name, mrsd = select_best_alternate_liver_model(
        ws=2.0, pt=1.0, disp=4.0, reference=5.0
    )
    assert name == "dispersion"
    assert mrsd == 4.0


def test_select_best_alternate_falls_back_to_well_stirred_when_no_improvement():
    # ws=5 (exact match), pt=100, disp=0.5; alternates both worse than ws
    name, mrsd = select_best_alternate_liver_model(
        ws=5.0, pt=100.0, disp=0.5, reference=5.0
    )
    assert name == "well_stirred"
    assert mrsd == 5.0


# ---------------------------------------------------------------------------
# compute_route_bias_factor
# ---------------------------------------------------------------------------

def test_route_bias_iv_reference_returns_unity():
    # IV reference: no 1/F bias exists, factor=1.0, no flags
    factor, flags = compute_route_bias_factor(route_ref="iv", f_lit=None)
    assert factor == 1.0
    assert flags == []


def test_route_bias_iv_reference_ignores_f_lit():
    # Even if f_lit supplied, IV route means factor=1.0
    factor, flags = compute_route_bias_factor(route_ref="iv", f_lit=0.3)
    assert factor == 1.0
    assert flags == []


def test_route_bias_oral_with_known_f():
    # propranolol F=0.26 -> factor = 1/0.26 = 3.846...
    factor, flags = compute_route_bias_factor(route_ref="oral", f_lit=0.26)
    assert factor == pytest.approx(1.0 / 0.26, rel=1e-9)
    assert flags == []


def test_route_bias_oral_with_unknown_f_flags():
    # Oral reference but no F literature -> factor=1.0 + flag
    factor, flags = compute_route_bias_factor(route_ref="oral", f_lit=None)
    assert factor == 1.0
    assert flags == ["f_unknown"]


# ---------------------------------------------------------------------------
# decompose_fold_error (end-to-end, SIGNED invariant)
# ---------------------------------------------------------------------------

def test_additivity_invariant_synthetic():
    """Signed-factor log-additivity.

    Hand-calc: ws=10, pt=5, disp=15, ref=5, F=0.5 (oral)
      fold_observed_signed = 10 / 5 = 2.0
      best alt: pt (exact match at 5.0). fold_liver_model_signed = 10 / 5 = 2.0
      fold_route_bias = 1 / 0.5 = 2.0
      fold_residual_signed = 2.0 / (2.0 * 2.0) = 0.5
      Invariant: log10(2.0) == log10(2.0) + log10(2.0) + log10(0.5)
                 0.30103     == 0.30103 + 0.30103 + (-0.30103) = 0.30103 ✓
    """
    result = decompose_fold_error(
        mrsd_ws=10.0,
        mrsd_pt=5.0,
        mrsd_disp=15.0,
        f_lit=0.5,
        route_ref="oral",
        fih_reference_mg=5.0,
    )
    log_lhs = math.log10(result.fold_observed_signed)
    log_rhs = (
        math.log10(result.fold_liver_model_signed)
        + math.log10(result.fold_route_bias)
        + math.log10(result.fold_residual_signed)
    )
    assert log_lhs == pytest.approx(log_rhs, rel=1e-9, abs=1e-9)


def test_decompose_concrete_signed_values():
    """Exact numeric check on every signed factor for the canonical case."""
    result = decompose_fold_error(
        mrsd_ws=10.0,
        mrsd_pt=5.0,
        mrsd_disp=15.0,
        f_lit=0.5,
        route_ref="oral",
        fih_reference_mg=5.0,
    )
    assert result.fold_observed_signed == pytest.approx(2.0, rel=1e-12)
    assert result.fold_liver_model_signed == pytest.approx(2.0, rel=1e-12)
    assert result.fold_route_bias == pytest.approx(2.0, rel=1e-12)
    assert result.fold_residual_signed == pytest.approx(0.5, rel=1e-12)
    assert result.best_alt_model_name == "parallel_tube"
    assert result.flags == []


def test_decompose_iv_reference_no_route_factor():
    """IV reference -> route bias factor must be exactly 1.0, zero log contribution."""
    result = decompose_fold_error(
        mrsd_ws=2.0,
        mrsd_pt=1.0,
        mrsd_disp=3.0,
        f_lit=None,
        route_ref="iv",
        fih_reference_mg=1.0,
    )
    assert result.fold_route_bias == 1.0
    assert "f_unknown" not in result.flags


def test_decompose_well_stirred_is_best_returns_unity_liver_factor():
    """When well_stirred is closest to reference, liver_model factor = 1.0 exactly."""
    # ws=5.0 (exact), pt=100, disp=0.5 -> alternates both worse
    result = decompose_fold_error(
        mrsd_ws=5.0,
        mrsd_pt=100.0,
        mrsd_disp=0.5,
        f_lit=None,
        route_ref="iv",
        fih_reference_mg=5.0,
    )
    assert result.fold_liver_model_signed == 1.0
    assert result.best_alt_model_name == "well_stirred"


def test_decompose_symmetric_report_from_signed():
    """to_symmetric applied to signed fields yields the benchmark-style fold."""
    result = decompose_fold_error(
        mrsd_ws=10.0,
        mrsd_pt=5.0,
        mrsd_disp=15.0,
        f_lit=0.5,
        route_ref="oral",
        fih_reference_mg=5.0,
    )
    # signed 2.0 -> symmetric 2.0; signed 0.5 -> symmetric 2.0
    assert to_symmetric(result.fold_observed_signed) == pytest.approx(2.0)
    assert to_symmetric(result.fold_residual_signed) == pytest.approx(2.0)
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
pytest tests/unit/test_decomposition.py -v
```

Expected: all tests FAIL with `ModuleNotFoundError: charon.translational.decomposition`.

- [ ] **Step 3: Implement the library**

```python
# src/charon/translational/decomposition.py
"""Sprint 10 research-only: fold-error decomposition for Layer 3 Tier A panel.

No production wiring — `charon.pipeline.Pipeline` does not import this module.
Everything is pure-function; I/O is the orchestrator's job
(`validation/benchmarks/layer3_ivive_decomposition.py`).

Signed log-additive decomposition:

    log10(fold_observed_signed) = log10(fold_liver_model_signed)
                                + log10(fold_route_bias)
                                + log10(fold_residual_signed)

where:

- `fold_observed_signed = mrsd_ws / fih_reference_mg` — the production
  (well_stirred) prediction relative to the clinical reference dose, in
  signed form (<1 means under-predicted, >1 means over-predicted).

- `fold_liver_model_signed = mrsd_ws / mrsd_best_alt` where best_alt is
  whichever of {parallel_tube, dispersion} minimises |log10(mrsd/reference)|.
  If neither alternate is closer to the reference than well_stirred, the
  factor is exactly 1.0 (no attributable improvement).

- `fold_route_bias = 1 / F_literature` for oral-reference compounds (a
  route-comparison artefact, not an IVIVE error — the production pipeline
  predicts IV MRSD and the reference is an oral dose). 1.0 for IV
  references or when F is unknown (flagged). Always >= 1.

- `fold_residual_signed = fold_observed_signed / (liver * route_bias)` —
  the unexplained remainder (transporter / non-hepatic / UGT / model-gap).

The helper `to_symmetric(signed)` converts any signed factor to the
benchmark-style symmetric fold `max(x, 1/x)` for reporting.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Literal


@dataclass(frozen=True)
class DecompositionResult:
    """Result of a single-compound fold-error decomposition (signed factors).

    Attributes:
        fold_observed_signed: mrsd_ws / fih_reference_mg.
        fold_liver_model_signed: mrsd_ws / mrsd_best_alt (1.0 if no
            alternate improves the prediction).
        fold_route_bias: 1/F_literature for oral references, 1.0 for IV
            or unknown-F.
        fold_residual_signed: observed / (liver * route_bias).
        best_alt_model_name: "parallel_tube", "dispersion", or "well_stirred".
        flags: list of string flags, e.g. ["f_unknown"].
    """

    fold_observed_signed: float
    fold_liver_model_signed: float
    fold_route_bias: float
    fold_residual_signed: float
    best_alt_model_name: str
    flags: list[str] = field(default_factory=list)


def to_symmetric(signed: float) -> float:
    """Convert a signed fold factor to the benchmark-style symmetric form.

    max(x, 1/x); 1.0 for x == 1.0. Returns math.inf for non-positive input.
    """
    if signed <= 0:
        return math.inf
    return max(signed, 1.0 / signed)


def select_best_alternate_liver_model(
    ws: float,
    pt: float,
    disp: float,
    reference: float,
) -> tuple[str, float]:
    """Return (name, mrsd) of the liver model closest to the reference.

    Closeness = minimum |log10(mrsd / reference)|. If neither alternate
    (parallel_tube, dispersion) is closer to reference than well_stirred,
    returns ("well_stirred", ws) — meaning no attributable improvement
    from liver-model choice.

    Tie-break when alternates equally close: prefer parallel_tube (lexical).
    """
    if reference <= 0:
        raise ValueError(f"reference must be > 0, got {reference}")
    candidates = {
        "well_stirred": ws,
        "parallel_tube": pt,
        "dispersion": disp,
    }
    distances = {
        name: abs(math.log10(v / reference)) if v > 0 else math.inf
        for name, v in candidates.items()
    }
    ws_dist = distances["well_stirred"]
    alt_sorted = sorted(
        [("parallel_tube", distances["parallel_tube"]),
         ("dispersion", distances["dispersion"])],
        key=lambda kv: kv[1],
    )
    best_alt_name, best_alt_dist = alt_sorted[0]
    if best_alt_dist < ws_dist:
        return best_alt_name, candidates[best_alt_name]
    return "well_stirred", ws


def compute_route_bias_factor(
    route_ref: Literal["iv", "oral"],
    f_lit: float | None,
) -> tuple[float, list[str]]:
    """Return (factor, flags) for the 1/F route-mismatch attribution.

    IV reference → (1.0, []). Oral + known F → (1/F, []). Oral + unknown F
    → (1.0, ["f_unknown"]).
    """
    if route_ref not in ("iv", "oral"):
        raise ValueError(
            f"route_ref must be 'iv' or 'oral', got {route_ref!r}"
        )
    if route_ref == "iv":
        return 1.0, []
    if f_lit is None:
        return 1.0, ["f_unknown"]
    if f_lit <= 0 or f_lit > 1.0:
        raise ValueError(f"f_lit must be in (0, 1], got {f_lit}")
    return 1.0 / f_lit, []


def decompose_fold_error(
    mrsd_ws: float,
    mrsd_pt: float,
    mrsd_disp: float,
    f_lit: float | None,
    route_ref: Literal["iv", "oral"],
    fih_reference_mg: float,
) -> DecompositionResult:
    """Decompose Tier A observed fold-error into three multiplicative factors (signed).

    Args:
        mrsd_ws: MRSD predicted with well_stirred liver model (production default).
        mrsd_pt: MRSD predicted with parallel_tube.
        mrsd_disp: MRSD predicted with dispersion.
        f_lit: Literature bioavailability (0 < F <= 1) for oral route_ref;
            None for IV or unknown.
        route_ref: 'iv' or 'oral' — the route of the clinical FIH reference
            dose (NOT the simulation route, which is always iv_bolus per
            Sprint 7's panel.yaml).
        fih_reference_mg: Reference FIH dose in mg.

    Returns:
        DecompositionResult with signed factors satisfying:
            log10(fold_observed_signed) ==
                log10(fold_liver_model_signed)
              + log10(fold_route_bias)
              + log10(fold_residual_signed)
    """
    if mrsd_ws <= 0 or fih_reference_mg <= 0:
        raise ValueError(
            f"mrsd_ws and fih_reference_mg must be > 0, "
            f"got {mrsd_ws}, {fih_reference_mg}"
        )

    fold_obs_signed = mrsd_ws / fih_reference_mg

    best_alt_name, best_alt_mrsd = select_best_alternate_liver_model(
        ws=mrsd_ws, pt=mrsd_pt, disp=mrsd_disp, reference=fih_reference_mg
    )
    if best_alt_name == "well_stirred":
        fold_liver_signed = 1.0
    else:
        fold_liver_signed = mrsd_ws / best_alt_mrsd

    fold_route, flags = compute_route_bias_factor(
        route_ref=route_ref, f_lit=f_lit
    )

    fold_residual_signed = fold_obs_signed / (fold_liver_signed * fold_route)

    return DecompositionResult(
        fold_observed_signed=fold_obs_signed,
        fold_liver_model_signed=fold_liver_signed,
        fold_route_bias=fold_route,
        fold_residual_signed=fold_residual_signed,
        best_alt_model_name=best_alt_name,
        flags=list(flags),
    )
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
pytest tests/unit/test_decomposition.py -v
```

Expected: 14 PASS (3 to_symmetric + 3 select_best + 4 route_bias + 4 decompose_fold_error). Zero failures.

- [ ] **Step 5: Commit**

```bash
git add src/charon/translational/decomposition.py tests/unit/test_decomposition.py
git commit -m "feat(sprint10): decomposition library for Layer 3 fold-error attribution"
```

---

## Task 3: Orchestrator script — per-compound × 3 liver models

**Files:**
- Create: `validation/benchmarks/layer3_ivive_decomposition.py`
- Create: `tests/integration/test_layer3_decomposition_benchmark.py`

- [ ] **Step 1: Write the failing integration test**

```python
# tests/integration/test_layer3_decomposition_benchmark.py
"""Integration tests for the Sprint 10 decomposition orchestrator.

Runs the full 12-compound Tier A decomposition and verifies structural
properties of the output.
"""

from __future__ import annotations

import json
import math
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "validation" / "benchmarks" / "layer3_ivive_decomposition.py"


@pytest.fixture(scope="module")
def decomposition_json(tmp_path_factory) -> dict:
    """Run the orchestrator once and cache the JSON output for all tests."""
    out_stem = tmp_path_factory.mktemp("sprint10") / "decomp"
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "--output-stem", str(out_stem)],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        timeout=600,
    )
    assert result.returncode == 0, (
        f"orchestrator failed: stdout={result.stdout} stderr={result.stderr}"
    )
    json_path = out_stem.with_suffix(".json")
    assert json_path.exists()
    return json.loads(json_path.read_text())


def test_decomposition_covers_all_12_tier_a(decomposition_json):
    panel_path = REPO_ROOT / "validation" / "data" / "fih_reference" / "panel.yaml"
    panel = yaml.safe_load(panel_path.read_text())["panel"]
    tier_a = [c["name"] for c in panel["compounds"] if c["tier"] == "gold"]
    rows = decomposition_json["extra_sections"]["Per-compound decomposition"]
    names_in_rows = {r["compound"] for r in rows}
    assert names_in_rows == set(tier_a), (
        f"missing: {set(tier_a) - names_in_rows}, "
        f"extra: {names_in_rows - set(tier_a)}"
    )
    assert len(rows) == 12


def test_decomposition_log_additivity_per_compound(decomposition_json):
    """For each of 12 compounds: log10 of signed observed fold equals
    the log10 sum of (liver, route, residual) signed factors within rtol=1e-6."""
    rows = decomposition_json["extra_sections"]["Per-compound decomposition"]
    for row in rows:
        lhs = math.log10(row["fold_observed_signed"])
        rhs = (
            math.log10(row["fold_liver_model_signed"])
            + math.log10(row["fold_route_bias"])
            + math.log10(row["fold_residual_signed"])
        )
        assert lhs == pytest.approx(rhs, rel=1e-6, abs=1e-9), (
            f"{row['compound']}: log-additivity violated "
            f"(lhs={lhs}, rhs={rhs})"
        )


def test_decomposition_every_row_numeric(decomposition_json):
    """No NaN / null in required numeric fields."""
    required = {
        "fold_observed",
        "fold_liver_model",
        "fold_route_bias",
        "fold_residual",
    }
    rows = decomposition_json["extra_sections"]["Per-compound decomposition"]
    for row in rows:
        for key in required:
            v = row[key]
            assert isinstance(v, (int, float))
            assert not math.isnan(v), f"{row['compound']}: {key} is NaN"
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
pytest tests/integration/test_layer3_decomposition_benchmark.py -v
```

Expected: all three FAIL with `FileNotFoundError` on missing script.

- [ ] **Step 3: Write the orchestrator**

```python
# validation/benchmarks/layer3_ivive_decomposition.py
"""Sprint 10 — Tier A IVIVE fold-error decomposition.

See docs/superpowers/specs/2026-04-23-sprint10-ivive-bias-diagnostic-design.md.

For each of 12 Tier A compounds:
  1. Run Pipeline with liver_model in {well_stirred, parallel_tube, dispersion}.
  2. Load literature F from bioavailability.csv.
  3. Call decompose_fold_error().
  4. Aggregate into a report.

This is research-only — no production code changes, no gating.
Exit 0 on success, 1 on data errors.
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from datetime import datetime, timezone
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from charon import Pipeline  # noqa: E402
from charon.core.schema import CompoundConfig, DoseProjectionConfig  # noqa: E402
from charon.translational.decomposition import (  # noqa: E402
    decompose_fold_error,
    to_symmetric,
)
from validation.benchmarks.report_writer import emit_report  # noqa: E402

DEFAULT_PANEL = REPO_ROOT / "validation" / "data" / "fih_reference" / "panel.yaml"
DEFAULT_BIOAVAILABILITY = (
    REPO_ROOT / "validation" / "data" / "fih_reference" / "bioavailability.csv"
)
COMPOUNDS_DIR = REPO_ROOT / "validation" / "data" / "tier1_obach" / "compounds"
DEFAULT_STEM = REPO_ROOT / "validation" / "reports" / "layer3_ivive_decomposition"

LIVER_MODELS = ("well_stirred", "parallel_tube", "dispersion")


def _load_compound(name: str) -> CompoundConfig:
    path = COMPOUNDS_DIR / f"{name}.yaml"
    if not path.exists():
        raise FileNotFoundError(f"compound YAML not found: {path}")
    return CompoundConfig.model_validate(yaml.safe_load(path.read_text()))


def _compute_mrsd(entry: dict, liver_model: str) -> float:
    """Single MRSD via PAD path with the specified liver model."""
    compound = _load_compound(entry["name"])
    pipe = Pipeline(
        compound,
        route=entry["route"],
        dose_mg=1.0,
        liver_model=liver_model,
        dose_projection=DoseProjectionConfig(
            target_ceff_nM=float(entry["target_ceff_nM"]),
            safety_factor=10.0,
            tau_h=24.0,
        ),
    )
    result = pipe.run()
    if result.dose_recommendation is None:
        raise RuntimeError(f"No dose recommendation for {entry['name']}")
    return float(result.dose_recommendation.mrsd_mg)


def _load_bioavailability(path: Path) -> dict[str, dict]:
    with path.open() as fh:
        rows = list(csv.DictReader(fh))
    out: dict[str, dict] = {}
    for r in rows:
        f_raw = r["f_oral"].strip()
        out[r["compound"]] = {
            "fih_reference_route": r["fih_reference_route"],
            "f_oral": float(f_raw) if f_raw else None,
            "f_source": r["f_source"],
            "f_doi_or_pmid": r["f_doi_or_pmid"],
            "notes": r["notes"],
        }
    return out


def run_panel(panel_path: Path, bioav_path: Path) -> dict:
    panel = yaml.safe_load(panel_path.read_text())["panel"]
    bioav = _load_bioavailability(bioav_path)

    rows: list[dict] = []
    tier_a = [c for c in panel["compounds"] if c["tier"] == "gold"]
    for entry in tier_a:
        name = entry["name"]
        if name not in bioav:
            raise RuntimeError(
                f"{name} missing from bioavailability.csv — rerun Task 1 curation"
            )
        mrsds = {m: _compute_mrsd(entry, m) for m in LIVER_MODELS}
        bioav_row = bioav[name]
        result = decompose_fold_error(
            mrsd_ws=mrsds["well_stirred"],
            mrsd_pt=mrsds["parallel_tube"],
            mrsd_disp=mrsds["dispersion"],
            f_lit=bioav_row["f_oral"],
            route_ref=bioav_row["fih_reference_route"],
            fih_reference_mg=float(entry["reference_fih_mg"]),
        )
        rows.append({
            "compound": name,
            "mrsd_ws_mg": mrsds["well_stirred"],
            "mrsd_pt_mg": mrsds["parallel_tube"],
            "mrsd_disp_mg": mrsds["dispersion"],
            "reference_fih_mg": float(entry["reference_fih_mg"]),
            "fih_reference_route": bioav_row["fih_reference_route"],
            "f_lit": bioav_row["f_oral"],
            # Signed factors (log-additivity invariant holds on these).
            "fold_observed_signed": result.fold_observed_signed,
            "fold_liver_model_signed": result.fold_liver_model_signed,
            "fold_route_bias": result.fold_route_bias,
            "fold_residual_signed": result.fold_residual_signed,
            # Symmetric forms for human-readable reporting.
            "fold_observed": to_symmetric(result.fold_observed_signed),
            "fold_liver_model": to_symmetric(result.fold_liver_model_signed),
            "fold_residual": to_symmetric(result.fold_residual_signed),
            "best_alt_model": result.best_alt_model_name,
            "flags": ",".join(result.flags) if result.flags else "-",
            "f_source": bioav_row["f_source"],
            "notes": bioav_row["notes"],
        })

    # Aggregate attribution percentages (use symmetric forms — always >= 1,
    # so abs(log10) is equivalent to abs(log10(signed))).
    def _pct(key: str) -> float:
        num = sum(abs(math.log10(r[key])) for r in rows if r[key] > 0)
        denom = sum(
            abs(math.log10(r["fold_observed"])) for r in rows if r["fold_observed"] > 0
        )
        return (num / denom * 100.0) if denom > 0 else 0.0

    summary = {
        "n_compounds": len(rows),
        "aggregate_pct_liver_model": _pct("fold_liver_model"),
        "aggregate_pct_route_bias": _pct("fold_route_bias"),
        "aggregate_pct_residual": _pct("fold_residual"),
    }

    # Sort rows by residual descending (worst-unexplained first) for report.
    rows_sorted = sorted(rows, key=lambda r: -r["fold_residual"])

    return {
        "title": "Charon Sprint 10 — Tier A IVIVE Fold-Error Decomposition",
        "panel": panel["name"],
        "date_utc": datetime.now(timezone.utc).isoformat(),
        "summary": summary,
        "rows": [],
        "extra_sections": {
            "Per-compound decomposition": rows_sorted,
        },
        "notes": [
            "Decomposition: fold_observed = fold_liver_model * fold_route_bias * fold_residual.",
            "fold_liver_model: ws/best_alt if alternate improves prediction, else 1.0.",
            "fold_route_bias: 1/F for oral-reference compounds, 1.0 for IV or unknown.",
            "fold_residual: unexplained remainder (transporter, non-hepatic, UGT, model-gap).",
            "Sorted by fold_residual descending (worst-unexplained first).",
            "Aggregate %: 100 * sum(|log10(factor)|) / sum(|log10(fold_observed)|).",
            "Research only — no production code changes.",
        ],
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel", type=Path, default=DEFAULT_PANEL)
    ap.add_argument("--bioavailability", type=Path, default=DEFAULT_BIOAVAILABILITY)
    ap.add_argument("--output-stem", type=Path, default=DEFAULT_STEM)
    args = ap.parse_args(argv)

    try:
        payload = run_panel(args.panel, args.bioavailability)
    except RuntimeError as exc:
        print(f"[FAIL] {exc}", file=sys.stderr)
        return 1

    emit_report(payload, stem=args.output_stem)
    print(
        f"[OK] Decomposition wrote {args.output_stem}.md + .json "
        f"({payload['summary']['n_compounds']} compounds)"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
```

- [ ] **Step 4: Run integration tests**

```bash
pytest tests/integration/test_layer3_decomposition_benchmark.py -v
```

Expected: 3 PASS. Takes ~2-3 minutes (12 compounds × 3 liver models = 36 Pipeline runs).

- [ ] **Step 5: Commit**

```bash
git add validation/benchmarks/layer3_ivive_decomposition.py tests/integration/test_layer3_decomposition_benchmark.py
git commit -m "feat(sprint10): decomposition orchestrator for Tier A 12-compound panel"
```

---

## Task 4: Generate report + review lisinopril ceff

**Files:**
- Generate: `validation/reports/layer3_ivive_decomposition.{md,json}` (produced by orchestrator)
- Create: `docs/superpowers/sprint10-followup-lisinopril-ceff.md` (if warranted)

- [ ] **Step 1: Run the orchestrator and commit the output**

```bash
python validation/benchmarks/layer3_ivive_decomposition.py
```

Expected stdout: `[OK] Decomposition wrote ... (12 compounds)`.

Inspect `validation/reports/layer3_ivive_decomposition.md`:

```bash
less validation/reports/layer3_ivive_decomposition.md
```

- [ ] **Step 2: Add narrative sections to the report**

The orchestrator's `emit_report` produces the table and summary but not free-text narrative. Amend the report post-hoc by appending these sections directly to the markdown file (this is a one-off research output, not a regenerated artefact).

Append to `validation/reports/layer3_ivive_decomposition.md`:

```markdown

## §3. Pattern analysis

### High-F-variance oral β-blockers
Compounds: propranolol (F=0.26), metoprolol (F=0.50), verapamil (F=0.22).
Expected dominant factor: `fold_route_bias`. Sprint 11 (Papp/oral migration)
is directly responsive to this category.

### Transporter-limited
Compound: atorvastatin (F=0.14; OATP1B1 hepatic uptake not modelled).
Expected signature: large `fold_route_bias` + non-trivial residual. The
residual quantifies the OATP unmodelled gap.

### Non-hepatic elimination
Compound: lisinopril (CLint ≈ 0, CLrenal-dominant).
Expected signature: residual dominates. See §4 for ceff review.

### CYP2C9 / UGT substrates
Compound: diclofenac (CYP2C9 + glucuronidation).
Expected: mixed — partial route bias + residual from IVIVE underprediction
(Obach-known pattern).

### Very-low fu_p (well-stirred sensitivity)
Compounds: warfarin (fu_p 0.01), diazepam (fu_p 0.013), atorvastatin (fu_p 0.02).
Expected: `fold_liver_model` carries signal if parallel_tube or dispersion
closer to reference.

## §4. Lisinopril target_ceff_nM literature review

Current panel value: `target_ceff_nM = 170 nM`.

Beermann 1988 (PMID:3048950) reports lisinopril steady-state Cmax after 10 mg
p.o. at approximately 60–80 ng/mL (MW 405.49 → ~148–197 nM). Gomez 1985
reports Cp_ss ~90 ng/mL (~222 nM) after 20 mg daily.

Recommendation: **170 nM is plausible** (within the 148–222 nM literature
band for 10–20 mg daily). No change to panel.yaml recommended from this
diagnostic. The lisinopril fold-error is driven by non-hepatic elimination
attribution, not an incorrect ceff target.

No follow-up ticket filed for this compound in Sprint 10.

## §5. Sprint 11+ priority ranking

Ranked by count of Tier A compounds whose **largest attributable factor**
is addressable by each candidate:

1. **Sprint 11 (Papp/Peff + oral route migration)** — directly removes
   `fold_route_bias` for every oral-reference compound (11/12). Upper-bound
   improvement per compound: full removal of 1/F factor.
2. **Sprint 12 (OATP1B1 plumbing)** — attributable to atorvastatin residual
   only (1/12 compound). Large per-compound effect (70x fold → potentially
   ~10x).
3. **Sprint 13 (UGT-specific calibration)** — diclofenac residual (1/12
   compound), plus acetaminophen if UGT accuracy is borderline.
4. **Liver-model selection policy** — if `fold_liver_model` explains
   meaningful aggregate %, production could auto-select among {ws, pt, disp}
   based on fu_p; if aggregate % is < 5%, de-prioritise.

Actual numeric percentages are in the summary at the top of this report —
update this prose ranking once numbers are in.
```

- [ ] **Step 3: Verify lisinopril literature conclusion**

If Beermann 1988 / Gomez 1985 disagree with the 170 nM value beyond the
148–222 nM band, instead create:

```bash
cat > docs/superpowers/sprint10-followup-lisinopril-ceff.md <<'EOF'
# Sprint 10 follow-up — lisinopril target_ceff_nM review

Current panel value: 170 nM.
Literature-recommended range (Sprint 10 review): <replace with finding>.
Proposed action: <change panel.yaml to X nM | keep current>.
EOF
git add docs/superpowers/sprint10-followup-lisinopril-ceff.md
```

Otherwise skip (the §4 conclusion suffices).

- [ ] **Step 4: Commit report + narrative**

```bash
git add validation/reports/layer3_ivive_decomposition.md validation/reports/layer3_ivive_decomposition.json
git commit -m "chore(sprint10): Tier A IVIVE decomposition report (n=12)"
```

---

## Task 5: Update Sprint 10 ticket status + run full test suite

**Files:**
- Modify: `docs/superpowers/sprint10-ivive-bias-ticket.md`

- [ ] **Step 1: Update ticket status footer**

Edit `docs/superpowers/sprint10-ivive-bias-ticket.md` and append a new section at the bottom:

```markdown

## Status — Sprint 10 diagnostic merged (2026-04-23)

Diagnostic produced at `validation/reports/layer3_ivive_decomposition.md`.
Key findings (fill from report summary after orchestrator run):

- Aggregate attribution: liver_model = N%, route_bias = N%, residual = N%.
- Sprint 11 (Papp/oral) addresses <N>/12 compounds as dominant factor.
- Sprint 12 (OATP) addresses 1/12 (atorvastatin) as dominant residual.

Remediation tickets:
- Sprint 11 — Papp/Peff curation + oral route migration (next).
- Sprint 12 — OATP1B1 plumbing (deferred).
```

- [ ] **Step 2: Run the full test suite**

```bash
pytest -q
```

Expected: ~905 tests pass (882 baseline + 14 unit decomposition + 6 unit bioavailability-csv + 3 integration). Zero failures. No new warnings beyond the pre-existing `pytest.mark.slow` warnings.

If any existing test fails, **stop** and report — do not proceed. Sprint 10 is not supposed to touch production paths, so any breakage is a signal the orchestrator accidentally wrote to shared state (it should not).

- [ ] **Step 3: Commit ticket status update**

```bash
git add docs/superpowers/sprint10-ivive-bias-ticket.md
git commit -m "docs(sprint10): update ticket with diagnostic findings summary"
```

---

## Self-Review Notes (for the implementer)

- **Type consistency check:** `decompose_fold_error` signature in Task 2 is `(mrsd_ws, mrsd_pt, mrsd_disp, f_lit, route_ref, fih_reference_mg)`. The Task 3 orchestrator must call it with these exact keyword args. `DecompositionResult` fields used by Task 3: `fold_observed_signed`, `fold_liver_model_signed`, `fold_route_bias`, `fold_residual_signed`, `best_alt_model_name`, `flags`. The orchestrator also imports `to_symmetric` to produce the symmetric reporting columns.

- **Liver-model parameter propagation:** `Pipeline(liver_model=...)` is already wired end-to-end (verified: `src/charon/pipeline.py:66`, `:189`; `src/charon/pbpk/ode_compiler.py:107`). No changes needed in production paths.

- **Panel route vs. FIH reference route:** `panel.yaml` has `route: iv_bolus` uniformly (Sprint 7 constraint). `bioavailability.csv`'s `fih_reference_route` is an independent annotation describing the *clinical reference dose*. Never confuse these.

- **Tier B is out of scope:** Sprint 10 only touches Tier A (gold, n=12). Tier B (sanity_floor, n=12) is untouched — it is a one-sided gated check that Sprint 9 kept passing.

- **Atomic write not required:** the orchestrator writes once via `emit_report`; it is not a cache, so no tmp+rename pattern. (Contrast with Sprint 7's `conformal.py` cache.)

- **No new dependencies:** `csv` is stdlib. `dataclasses` is stdlib. No changes to `pyproject.toml`.

- **Commit count estimate:** 5 feature commits (Tasks 1–5) + merge commit.
