# Sprint 6 — Report Generation + CLI Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Turn the working Charon pipeline into a shippable tool with a regulatory-oriented FIH dose rationale report and a 5-command CLI.

**Architecture:** `PipelineResult → Collector (flatten to ReportData) → Narrative (render Markdown) → Export (write .md and .json)`. CLI (`cli/main.py`, argparse) wires `predict` / `simulate` / `translate` / `recommend` / `report` subcommands through the existing `Pipeline`. Zero new runtime dependencies — argparse and f-string only.

**Tech Stack:** Python 3.11, argparse (stdlib), dataclasses (stdlib), json (stdlib), pytest. Reuses existing `charon.pipeline.Pipeline`, `charon.predict.predict_properties`, and schema/dataclasses.

**Spec:** `docs/superpowers/specs/2026-04-14-sprint6-report-cli-design.md`

---

## File Structure

| File | Purpose | Task |
|------|---------|------|
| `src/charon/report/collector.py` | `ReportData` dataclass + `collect()` (pure data flattening) | 1, 2, 3 |
| `src/charon/report/narrative.py` | `render_report()` + 9 section renderers (pure f-string) | 4, 5, 6, 7 |
| `src/charon/report/export.py` | `export_markdown`, `export_json`, `export_report` (file I/O) | 8 |
| `src/charon/report/__init__.py` | Public symbol exports | 8 |
| `src/charon/cli/main.py` | argparse dispatcher + 5 subcommand handlers | 9, 10, 11 |
| `src/charon/cli/__init__.py` | Re-export `main` | 9 |
| `pyproject.toml` | Add `[project.scripts] charon = ...` | 9 |
| `tests/unit/test_collector.py` | Collector unit tests | 1, 2, 3 |
| `tests/unit/test_narrative.py` | Narrative unit tests | 4, 5, 6, 7 |
| `tests/unit/test_export.py` | Export unit tests | 8 |
| `tests/unit/test_cli.py` | CLI unit tests (via `main([...])` capsys) | 9, 10, 11 |
| `tests/integration/test_report_e2e.py` | End-to-end SMILES → report | 12 |

Totals: 6 source files, 5 test files, 1 config file. All files < 400 lines.

---

## Preflight

- [ ] **Preflight 1: Confirm baseline tests pass**

Run: `pytest tests/ --tb=no -q 2>&1 | tail -5`
Expected: `678 passed`.

- [ ] **Preflight 2: Confirm empty scaffolds exist**

Run: `ls -la src/charon/report/ src/charon/cli/`
Expected: `collector.py`, `narrative.py`, `export.py`, `figures.py`, `__init__.py` (all empty), `templates/` dir with 3 empty YAMLs; `cli/main.py` and `cli/__init__.py` empty.

---

## Task 1: ReportData dataclass + minimal collect() skeleton

**Files:**
- Create: `src/charon/report/collector.py`
- Test: `tests/unit/test_collector.py`

Goal: `ReportData` frozen dataclass with all required fields and a `collect()` function that populates identity + route + dose + metadata/timestamp from a minimal `PipelineResult`. Property flattening, IVIVE, PK table, dose/uncertainty are stubbed to empty defaults and filled in subsequent tasks.

- [ ] **Step 1: Write the failing test**

Create `tests/unit/test_collector.py`:

```python
"""Unit tests for charon.report.collector."""

from __future__ import annotations

from dataclasses import is_dataclass

import numpy as np
import pytest

from charon.core.schema import (
    CompoundConfig,
    CompoundProperties,
    PKParameters,
)
from charon.pbpk.solver import SimulationResult
from charon.pipeline import PipelineResult
from charon.report.collector import ReportData, collect


def _make_minimal_result(
    *,
    name: str = "test-compound",
    smiles: str = "CCO",
    route: str = "iv_bolus",
    dose_mg: float = 10.0,
) -> PipelineResult:
    """Build a PipelineResult with just enough data for collector tests."""
    compound = CompoundConfig(
        name=name,
        smiles=smiles,
        molecular_weight=46.07,
        source="predicted",
        properties=CompoundProperties(),
    )
    time_h = np.array([0.0, 1.0, 2.0, 4.0, 8.0, 24.0], dtype=float)
    cp_plasma = np.array([10.0, 8.0, 6.0, 4.0, 2.0, 0.5], dtype=float)
    cp_blood = cp_plasma * 0.9
    sim = SimulationResult(
        time_h=time_h,
        cp_plasma=cp_plasma,
        cp_blood=cp_blood,
        tissue_trajectories={},
        solver_method="BDF",
        solver_nfev=42,
    )
    pk = PKParameters(
        cmax=10.0,
        tmax=0.0,
        auc_0_inf=50.0,
        half_life=4.0,
        cl_apparent=0.2,
        vss=1.5,
    )
    return PipelineResult(
        compound=compound,
        pk_parameters=pk,
        time_h=time_h,
        cp_plasma=cp_plasma,
        cp_blood=cp_blood,
        simulation=sim,
        metadata={
            "species": "human",
            "route": route,
            "dose_mg": dose_mg,
            "duration_h": 24.0,
            "liver_model": "well_stirred",
        },
    )


def test_report_data_is_frozen_dataclass():
    assert is_dataclass(ReportData)
    assert ReportData.__dataclass_params__.frozen is True  # type: ignore[attr-defined]


def test_collect_identity_fields():
    result = _make_minimal_result(name="ethanol", smiles="CCO", dose_mg=10.0)
    data = collect(result)
    assert data.compound_name == "ethanol"
    assert data.smiles == "CCO"
    assert data.molecular_weight == pytest.approx(46.07)
    assert data.source == "predicted"
    assert data.route == "iv_bolus"
    assert data.dose_mg == pytest.approx(10.0)
    assert data.duration_h == pytest.approx(24.0)


def test_collect_timestamp_is_iso8601():
    result = _make_minimal_result()
    data = collect(result)
    assert "T" in data.timestamp
    assert len(data.timestamp) >= 19


def test_collect_respects_explicit_timestamp():
    result = _make_minimal_result()
    data = collect(result, timestamp="2026-04-15T12:00:00+00:00")
    assert data.timestamp == "2026-04-15T12:00:00+00:00"


def test_collect_warnings_default_empty():
    result = _make_minimal_result()
    data = collect(result)
    assert data.warnings == []


def test_collect_warnings_passthrough():
    result = _make_minimal_result()
    data = collect(result, warnings=["high uncertainty"])
    assert data.warnings == ["high uncertainty"]
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/unit/test_collector.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'charon.report.collector'` (or similar import error — the file is empty).

- [ ] **Step 3: Write minimal implementation**

Create `src/charon/report/collector.py`:

```python
"""Sprint 6 — Report data collection.

Flattens a :class:`charon.pipeline.PipelineResult` into a
:class:`ReportData` dataclass suitable for consumption by the narrative
renderer and the JSON exporter.  The collector performs no numerical
computation beyond indexing and copying — any derived metric must be
computed upstream in the pipeline.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime, timezone

from charon.pipeline import PipelineResult


@dataclass(frozen=True)
class ReportData:
    """Structured data for the Sprint 6 regulatory report.

    All fields are plain Python values (no numpy arrays, no Pydantic
    models).  This makes the dataclass trivially JSON-serialisable.
    """

    # Identity
    compound_name: str
    smiles: str
    molecular_weight: float | None
    source: str
    compound_type: str | None

    # Layer 1: ADME properties, flattened
    properties: dict[str, dict]

    # Bridge: IVIVE audit fields (pulled verbatim from metadata)
    ivive_summary: dict

    # Layer 2: PK parameters + sampled Cp-time table
    pk_params: dict[str, float | None]
    pk_table: list[dict]
    route: str
    dose_mg: float
    duration_h: float

    # Layer 3 / Layer 4 (optional)
    dose_recommendation: dict | None
    uncertainty: dict | None

    # Layer 0 + run metadata
    warnings: list[str] = field(default_factory=list)
    metadata: dict = field(default_factory=dict)
    timestamp: str = ""
    charon_version: str = "0.1.0"


def _iso_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def collect(
    result: PipelineResult,
    *,
    warnings: list[str] | None = None,
    timestamp: str | None = None,
) -> ReportData:
    """Flatten a ``PipelineResult`` into a ``ReportData`` snapshot.

    Subsequent tasks extend this function to populate ADME properties,
    IVIVE summary, PK table, dose recommendation, and uncertainty.
    """
    md = result.metadata or {}
    compound = result.compound

    return ReportData(
        compound_name=compound.name,
        smiles=compound.smiles,
        molecular_weight=compound.molecular_weight,
        source=compound.source,
        compound_type=compound.properties.physicochemical.compound_type,
        properties={},
        ivive_summary={},
        pk_params={},
        pk_table=[],
        route=str(md.get("route", "")),
        dose_mg=float(md.get("dose_mg", 0.0)),
        duration_h=float(md.get("duration_h", 0.0)),
        dose_recommendation=None,
        uncertainty=None,
        warnings=list(warnings) if warnings else [],
        metadata=dict(md),
        timestamp=timestamp or _iso_now(),
    )
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/unit/test_collector.py -v`
Expected: 6 passed.

- [ ] **Step 5: Commit**

```bash
git add src/charon/report/collector.py tests/unit/test_collector.py
git commit -m "feat(report): ReportData dataclass + minimal collect() skeleton"
```

---

## Task 2: Flatten ADME properties + IVIVE summary

**Files:**
- Modify: `src/charon/report/collector.py`
- Test: `tests/unit/test_collector.py`

Goal: `collect()` walks `result.compound.properties` emitting a flat
`{name: {value, ci_lower, ci_upper, unit, source, flag}}` dict, and builds
an `ivive_summary` dict from `result.metadata`.

- [ ] **Step 1: Write the failing test**

Append to `tests/unit/test_collector.py`:

```python
from charon.core.schema import (
    BindingProperties,
    MetabolismProperties,
    PhysicochemicalProperties,
    PredictedProperty,
)


def _make_result_with_properties() -> PipelineResult:
    compound = CompoundConfig(
        name="midazolam",
        smiles="Cc1ncc(n1C)[C@@H](c2ccccc2)OC(=O)N",
        molecular_weight=325.77,
        source="experimental",
        properties=CompoundProperties(
            physicochemical=PhysicochemicalProperties(
                logp=PredictedProperty(
                    value=3.89,
                    ci_90_lower=3.4,
                    ci_90_upper=4.3,
                    source="ml_ensemble",
                    unit="log",
                ),
                compound_type="base",
            ),
            binding=BindingProperties(
                fu_p=PredictedProperty(
                    value=0.032,
                    ci_90_lower=0.02,
                    ci_90_upper=0.05,
                    source="ml_ensemble",
                    unit="fraction",
                ),
            ),
            metabolism=MetabolismProperties(
                clint_uL_min_mg=PredictedProperty(
                    value=93.0,
                    source="experimental",
                    unit="uL/min/mg",
                    flag="clint_tier2_ml",
                ),
            ),
        ),
    )
    time_h = np.array([0.0, 1.0, 2.0])
    cp = np.array([1.0, 0.5, 0.25])
    sim = SimulationResult(
        time_h=time_h,
        cp_plasma=cp,
        cp_blood=cp,
        tissue_trajectories={},
        solver_method="BDF",
        solver_nfev=10,
    )
    return PipelineResult(
        compound=compound,
        pk_parameters=PKParameters(),
        time_h=time_h,
        cp_plasma=cp,
        cp_blood=cp,
        simulation=sim,
        metadata={
            "route": "oral",
            "dose_mg": 5.0,
            "duration_h": 24.0,
            "liver_model": "well_stirred",
            "clint_liver_L_h": 348.75,
            "cl_renal_L_h": 0.23,
            "fu_b": 0.071,
            "compound_type": "base",
        },
    )


def test_collect_properties_flattened():
    result = _make_result_with_properties()
    data = collect(result)
    # logp
    assert "logp" in data.properties
    assert data.properties["logp"]["value"] == pytest.approx(3.89)
    assert data.properties["logp"]["ci_lower"] == pytest.approx(3.4)
    assert data.properties["logp"]["ci_upper"] == pytest.approx(4.3)
    assert data.properties["logp"]["source"] == "ml_ensemble"
    # fu_p
    assert "fu_p" in data.properties
    assert data.properties["fu_p"]["value"] == pytest.approx(0.032)
    # clint
    assert "clint_uL_min_mg" in data.properties
    assert data.properties["clint_uL_min_mg"]["flag"] == "clint_tier2_ml"


def test_collect_compound_type():
    result = _make_result_with_properties()
    data = collect(result)
    assert data.compound_type == "base"


def test_collect_ivive_summary_from_metadata():
    result = _make_result_with_properties()
    data = collect(result)
    assert data.ivive_summary["clint_liver_L_h"] == pytest.approx(348.75)
    assert data.ivive_summary["cl_renal_L_h"] == pytest.approx(0.23)
    assert data.ivive_summary["fu_b"] == pytest.approx(0.071)
    assert data.ivive_summary["liver_model"] == "well_stirred"


def test_collect_skips_none_properties():
    result = _make_minimal_result()
    data = collect(result)
    # No properties populated in the minimal fixture
    assert data.properties == {}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/unit/test_collector.py::test_collect_properties_flattened -v`
Expected: FAIL — `properties` is still `{}`.

- [ ] **Step 3: Write minimal implementation**

Edit `src/charon/report/collector.py`. Add these helpers above `collect()`:

```python
from charon.core.schema import PredictedProperty


_PROPERTY_FIELDS: list[tuple[str, str, str]] = [
    # (category, attribute_on_category, report_key)
    ("physicochemical", "logp", "logp"),
    ("physicochemical", "pka_acid", "pka_acid"),
    ("physicochemical", "pka_base", "pka_base"),
    ("physicochemical", "solubility_ug_ml", "solubility_ug_ml"),
    ("binding", "fu_p", "fu_p"),
    ("binding", "fu_inc", "fu_inc"),
    ("binding", "bp_ratio", "bp_ratio"),
    ("metabolism", "clint_uL_min_mg", "clint_uL_min_mg"),
    ("permeability", "papp_nm_s", "papp_nm_s"),
    ("permeability", "peff_cm_s", "peff_cm_s"),
    ("safety", "herg_ic50_uM", "herg_ic50_uM"),
    ("renal", "clrenal_L_h", "clrenal_L_h"),
]


def _flatten_property(prop: PredictedProperty) -> dict:
    return {
        "value": float(prop.value),
        "ci_lower": None if prop.ci_90_lower is None else float(prop.ci_90_lower),
        "ci_upper": None if prop.ci_90_upper is None else float(prop.ci_90_upper),
        "unit": prop.unit,
        "source": prop.source,
        "flag": prop.flag,
        "method": prop.method,
    }


def _flatten_properties(props) -> dict[str, dict]:
    out: dict[str, dict] = {}
    for category, attr, key in _PROPERTY_FIELDS:
        section = getattr(props, category, None)
        if section is None:
            continue
        prop = getattr(section, attr, None)
        if prop is None:
            continue
        out[key] = _flatten_property(prop)
    return out


_IVIVE_KEYS = (
    "clint_liver_L_h",
    "cl_renal_L_h",
    "fu_b",
    "liver_model",
    "compound_type",
    "clint_gut_L_h",
)


def _ivive_summary_from_metadata(md: dict) -> dict:
    return {k: md[k] for k in _IVIVE_KEYS if k in md}
```

Then in the `collect()` call, replace:

```python
        properties={},
        ivive_summary={},
```

with:

```python
        properties=_flatten_properties(compound.properties),
        ivive_summary=_ivive_summary_from_metadata(md),
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/unit/test_collector.py -v`
Expected: 10 passed.

- [ ] **Step 5: Commit**

```bash
git add src/charon/report/collector.py tests/unit/test_collector.py
git commit -m "feat(report): flatten ADME properties and IVIVE summary in collect()"
```

---

## Task 3: Flatten PK parameters, PK table, dose recommendation, uncertainty

**Files:**
- Modify: `src/charon/report/collector.py`
- Test: `tests/unit/test_collector.py`

Goal: `collect()` populates `pk_params`, `pk_table` (canonical timepoint sampling), `dose_recommendation`, and `uncertainty`.

- [ ] **Step 1: Write the failing test**

Append to `tests/unit/test_collector.py`:

```python
from dataclasses import asdict as dc_asdict

from charon.translational.dose_projector import FIHDoseRecommendation
from charon.translational.hed import HEDResult
from charon.uncertainty.dose_range import UncertaintyResult


def test_collect_pk_params_flattened():
    result = _make_minimal_result()
    data = collect(result)
    assert data.pk_params["cmax"] == pytest.approx(10.0)
    assert data.pk_params["auc_0_inf"] == pytest.approx(50.0)
    assert data.pk_params["cl_apparent"] == pytest.approx(0.2)
    assert data.pk_params["vss"] == pytest.approx(1.5)
    assert "bioavailability" in data.pk_params


def test_collect_pk_table_uses_canonical_timepoints():
    result = _make_minimal_result()
    data = collect(result)
    # Canonical points ≤ max(time_h=24.0): 0, 0.25, 0.5, 1, 2, 4, 6, 8, 12, 24
    times = [row["time_h"] for row in data.pk_table]
    assert 0.0 in times
    assert 1.0 in times
    assert 24.0 in times
    # 48 and 72 are beyond duration
    assert 48.0 not in times
    assert 72.0 not in times
    # Monotonic and ≤ 12 entries
    assert times == sorted(times)
    assert len(data.pk_table) <= 12


def test_collect_pk_table_cp_values():
    result = _make_minimal_result()
    data = collect(result)
    first = data.pk_table[0]
    assert "cp_plasma_ug_L" in first
    assert "cp_blood_ug_L" in first
    assert isinstance(first["cp_plasma_ug_L"], float)


def test_collect_dose_recommendation_flattened():
    result = _make_minimal_result()
    hed = HEDResult(
        hed_mg_kg=8.06,
        mrsd_mg=56.45,
        species="rat",
        km_animal=6.2,
        km_human=37.0,
        noael_mg_kg=50.0,
        safety_factor=10.0,
        body_weight_kg=70.0,
    )
    rec = FIHDoseRecommendation(
        mrsd_mg=56.45,
        limiting_method="hed",
        hed=hed,
        mabel=None,
        pad=None,
        safety_factor=10.0,
        salt_factor=1.0,
        route="oral",
        rationale="HED: 56.45 mg",
    )
    result.dose_recommendation = rec
    data = collect(result)
    assert data.dose_recommendation is not None
    assert data.dose_recommendation["mrsd_mg"] == pytest.approx(56.45)
    assert data.dose_recommendation["limiting_method"] == "hed"
    assert data.dose_recommendation["hed"]["mrsd_mg"] == pytest.approx(56.45)
    assert data.dose_recommendation["mabel"] is None


def test_collect_uncertainty_flattened():
    result = _make_minimal_result()
    unc = UncertaintyResult(
        point_estimate_mg=5.0,
        ci_90_lower_mg=2.0,
        ci_90_upper_mg=12.0,
        ci_ratio=6.0,
        confidence="MEDIUM",
        n_samples=100,
        n_successful=95,
        convergence_met=True,
        sensitivity={"clint": 0.7, "fu_p": 0.2, "logp": 0.1},
        limiting_parameter="clint",
        recommendation="Experimental clint measurement would narrow CI by ~70%",
        r_squared=0.85,
    )
    result.uncertainty = unc
    data = collect(result)
    assert data.uncertainty is not None
    assert data.uncertainty["point_estimate_mg"] == pytest.approx(5.0)
    assert data.uncertainty["confidence"] == "MEDIUM"
    assert data.uncertainty["sensitivity"]["clint"] == pytest.approx(0.7)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/unit/test_collector.py -v`
Expected: 5 new tests FAIL — `pk_params` is empty, `dose_recommendation` is None, etc.

- [ ] **Step 3: Write minimal implementation**

Edit `src/charon/report/collector.py`. Add imports and helpers:

```python
from dataclasses import asdict as _dc_asdict

import numpy as np

from charon.core.schema import PKParameters
from charon.translational.dose_projector import FIHDoseRecommendation
from charon.uncertainty.dose_range import UncertaintyResult


_CANONICAL_TIMEPOINTS_H: tuple[float, ...] = (
    0.0, 0.25, 0.5, 1.0, 2.0, 4.0, 6.0, 8.0, 12.0, 24.0, 48.0, 72.0,
)


def _flatten_pk_params(pk: PKParameters) -> dict[str, float | None]:
    return {
        "cmax": pk.cmax,
        "tmax": pk.tmax,
        "auc_0_inf": pk.auc_0_inf,
        "auc_0_24": pk.auc_0_24,
        "half_life": pk.half_life,
        "cl_apparent": pk.cl_apparent,
        "vss": pk.vss,
        "bioavailability": pk.bioavailability,
        "fa": pk.fa,
        "fg": pk.fg,
        "fh": pk.fh,
    }


def _sample_pk_table(
    time_h: np.ndarray,
    cp_plasma: np.ndarray,
    cp_blood: np.ndarray,
) -> list[dict]:
    t_arr = np.asarray(time_h, dtype=float)
    if t_arr.size == 0:
        return []
    t_max = float(t_arr[-1])
    rows: list[dict] = []
    seen_idx: set[int] = set()
    for t in _CANONICAL_TIMEPOINTS_H:
        if t > t_max:
            break
        idx = int(np.searchsorted(t_arr, t))
        if idx >= t_arr.size:
            idx = t_arr.size - 1
        # searchsorted gives left insertion; pick the nearest of idx-1/idx
        if idx > 0 and abs(t_arr[idx - 1] - t) <= abs(t_arr[idx] - t):
            idx -= 1
        if idx in seen_idx:
            continue
        seen_idx.add(idx)
        rows.append(
            {
                "time_h": float(t_arr[idx]),
                "cp_plasma_ug_L": float(cp_plasma[idx]),
                "cp_blood_ug_L": float(cp_blood[idx]),
            }
        )
    return rows


def _flatten_dose_recommendation(
    rec: FIHDoseRecommendation | None,
) -> dict | None:
    if rec is None:
        return None
    return rec.model_dump()


def _flatten_uncertainty(unc: UncertaintyResult | None) -> dict | None:
    if unc is None:
        return None
    return _dc_asdict(unc)
```

Then in `collect()`, replace:

```python
        pk_params={},
        pk_table=[],
        ...
        dose_recommendation=None,
        uncertainty=None,
```

with:

```python
        pk_params=_flatten_pk_params(result.pk_parameters),
        pk_table=_sample_pk_table(result.time_h, result.cp_plasma, result.cp_blood),
        ...
        dose_recommendation=_flatten_dose_recommendation(result.dose_recommendation),
        uncertainty=_flatten_uncertainty(result.uncertainty),
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/unit/test_collector.py -v`
Expected: 15 passed.

- [ ] **Step 5: Commit**

```bash
git add src/charon/report/collector.py tests/unit/test_collector.py
git commit -m "feat(report): flatten PK params, PK table, dose rec, uncertainty"
```

---

## Task 4: Narrative — header, executive summary, compound profile (sections 1–2)

**Files:**
- Create: `src/charon/report/narrative.py`
- Test: `tests/unit/test_narrative.py`

Goal: Formatting helpers (`format_value`) and three renderers: `_render_header`, `_render_executive_summary`, `_render_compound_profile`. No full `render_report()` yet.

- [ ] **Step 1: Write the failing test**

Create `tests/unit/test_narrative.py`:

```python
"""Unit tests for charon.report.narrative."""

from __future__ import annotations

import pytest

from charon.report.collector import ReportData
from charon.report.narrative import (
    _render_compound_profile,
    _render_executive_summary,
    _render_header,
    format_value,
)


def _make_data(**overrides) -> ReportData:
    base = dict(
        compound_name="midazolam",
        smiles="Cc1ncc(n1C)[C@@H](c2ccccc2)OC(=O)N",
        molecular_weight=325.77,
        source="experimental",
        compound_type="base",
        properties={},
        ivive_summary={},
        pk_params={},
        pk_table=[],
        route="oral",
        dose_mg=5.0,
        duration_h=24.0,
        dose_recommendation=None,
        uncertainty=None,
        warnings=[],
        metadata={"species": "human"},
        timestamp="2026-04-15T00:00:00+00:00",
        charon_version="0.1.0",
    )
    base.update(overrides)
    return ReportData(**base)


def test_format_value_none():
    assert format_value(None) == "-"


def test_format_value_small():
    assert format_value(3.89) == "3.89"


def test_format_value_large_scientific():
    s = format_value(1234567.0)
    assert "e" in s or "E" in s or "1.23" in s


def test_format_value_tiny_scientific():
    s = format_value(0.000123)
    assert "e" in s or "E" in s


def test_render_header_contains_name():
    data = _make_data()
    out = _render_header(data)
    assert "midazolam" in out
    assert out.startswith("# ")


def test_render_header_contains_timestamp():
    data = _make_data()
    out = _render_header(data)
    assert "2026-04-15" in out


def test_render_executive_summary_without_dose():
    data = _make_data()
    out = _render_executive_summary(data)
    assert "## 1. Executive Summary" in out
    # No dose rec -> pipeline-only note
    assert "No FIH dose projection" in out or "not run" in out.lower()


def test_render_executive_summary_with_dose_no_uncertainty():
    data = _make_data(
        dose_recommendation={
            "mrsd_mg": 56.45,
            "limiting_method": "hed",
            "route": "oral",
            "salt_factor": 1.0,
            "safety_factor": 10.0,
            "rationale": "HED",
            "hed": {"mrsd_mg": 56.45},
            "mabel": None,
            "pad": None,
        }
    )
    out = _render_executive_summary(data)
    assert "56.45" in out
    assert "HED" in out or "hed" in out.lower()


def test_render_executive_summary_with_uncertainty():
    data = _make_data(
        dose_recommendation={
            "mrsd_mg": 56.45,
            "limiting_method": "hed",
            "route": "oral",
            "salt_factor": 1.0,
            "safety_factor": 10.0,
            "rationale": "HED",
            "hed": {"mrsd_mg": 56.45},
            "mabel": None,
            "pad": None,
        },
        uncertainty={
            "point_estimate_mg": 50.0,
            "ci_90_lower_mg": 20.0,
            "ci_90_upper_mg": 120.0,
            "confidence": "LOW",
            "ci_ratio": 6.0,
            "n_samples": 100,
            "n_successful": 95,
            "sensitivity": {},
            "limiting_parameter": "clint",
            "recommendation": "",
            "convergence_met": True,
            "r_squared": 0.9,
        },
    )
    out = _render_executive_summary(data)
    assert "20" in out and "120" in out  # CI bounds
    assert "LOW" in out


def test_render_compound_profile_lists_identity():
    data = _make_data()
    out = _render_compound_profile(data)
    assert "## 2. Compound Profile" in out
    assert "Cc1ncc" in out  # SMILES
    assert "325.77" in out
    assert "base" in out
    assert "experimental" in out
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/unit/test_narrative.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'charon.report.narrative'`.

- [ ] **Step 3: Write minimal implementation**

Create `src/charon/report/narrative.py`:

```python
"""Sprint 6 — Markdown narrative rendering.

Pure f-string based rendering of a :class:`ReportData` into a
regulatory-oriented Markdown document.  No template engine, no new
runtime dependencies.
"""

from __future__ import annotations

from charon.report.collector import ReportData


def format_value(v: float | int | None, *, digits: int = 4) -> str:
    """Format a numeric value for a Markdown table cell.

    - ``None`` → ``"-"``
    - ``abs(v) >= 1000`` or ``abs(v) < 0.01`` (and not zero) → scientific
      with ``digits-1`` significant figures.
    - otherwise → fixed-ish with ``digits`` significant figures.
    """
    if v is None:
        return "-"
    try:
        fv = float(v)
    except (TypeError, ValueError):
        return str(v)
    if fv == 0.0:
        return "0"
    a = abs(fv)
    if a >= 1000.0 or a < 0.01:
        return f"{fv:.{digits - 1}e}"
    return f"{fv:.{digits}g}"


def _render_header(data: ReportData) -> str:
    return (
        f"# FIH Dose Rationale Report — {data.compound_name}\n\n"
        f"*Generated: {data.timestamp} · "
        f"Charon {data.charon_version}*"
    )


def _render_executive_summary(data: ReportData) -> str:
    lines: list[str] = ["## 1. Executive Summary", ""]
    rec = data.dose_recommendation
    unc = data.uncertainty

    if rec is None:
        lines.append(
            "No FIH dose projection was run for this pipeline execution "
            "(no NOAEL, target Kd, or target Ceff provided)."
        )
    else:
        mrsd = format_value(rec.get("mrsd_mg"))
        method = str(rec.get("limiting_method", "?")).upper()
        route = rec.get("route", data.route)
        headline = (
            f"**Recommended FIH starting dose: {mrsd} mg** "
            f"({route}, limiting method: {method})"
        )
        if unc is not None:
            lo = format_value(unc.get("ci_90_lower_mg"))
            hi = format_value(unc.get("ci_90_upper_mg"))
            conf = unc.get("confidence", "?")
            headline += (
                f"  \n90% CI: [{lo} – {hi}] mg  ·  confidence: **{conf}**"
            )
        lines.append(headline)
        lines.append("")
        lines.append(
            f"Compound: {data.compound_name} ({data.smiles}). "
            f"Source: {data.source}. Route: {data.route}, "
            f"simulated dose: {format_value(data.dose_mg)} mg."
        )

    return "\n".join(lines)


def _render_compound_profile(data: ReportData) -> str:
    lines = ["## 2. Compound Profile", ""]
    lines.append(f"- **Name:** {data.compound_name}")
    lines.append(f"- **SMILES:** `{data.smiles}`")
    lines.append(
        f"- **Molecular weight:** {format_value(data.molecular_weight)} g/mol"
    )
    lines.append(f"- **Compound type:** {data.compound_type or '-'}")
    lines.append(f"- **Source:** {data.source}")
    return "\n".join(lines)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/unit/test_narrative.py -v`
Expected: 10 passed.

- [ ] **Step 5: Commit**

```bash
git add src/charon/report/narrative.py tests/unit/test_narrative.py
git commit -m "feat(report): narrative header, executive summary, compound profile"
```

---

## Task 5: Narrative — ADME table and IVIVE section (sections 3–4)

**Files:**
- Modify: `src/charon/report/narrative.py`
- Test: `tests/unit/test_narrative.py`

- [ ] **Step 1: Write the failing test**

Append to `tests/unit/test_narrative.py`:

```python
from charon.report.narrative import _render_adme_table, _render_ivive_audit


def test_render_adme_table_header():
    data = _make_data(
        properties={
            "logp": {
                "value": 3.89,
                "ci_lower": 3.4,
                "ci_upper": 4.3,
                "unit": "log",
                "source": "ml_ensemble",
                "flag": None,
                "method": None,
            },
        }
    )
    out = _render_adme_table(data)
    assert "## 3. ADME Predictions" in out
    assert "| Property" in out
    assert "| --- " in out or "|---" in out
    assert "logp" in out
    assert "3.89" in out


def test_render_adme_table_skips_missing():
    data = _make_data(properties={})
    out = _render_adme_table(data)
    # Header still emitted
    assert "## 3. ADME Predictions" in out
    # But no rows
    assert "logp" not in out
    assert "fu_p" not in out


def test_render_adme_table_ci_dash_when_missing():
    data = _make_data(
        properties={
            "fu_p": {
                "value": 0.032,
                "ci_lower": None,
                "ci_upper": None,
                "unit": "fraction",
                "source": "experimental",
                "flag": None,
                "method": None,
            }
        }
    )
    out = _render_adme_table(data)
    # The row should have a dash in the CI column
    row_lines = [ln for ln in out.split("\n") if ln.startswith("| fu_p")]
    assert len(row_lines) == 1
    # CI column is the third pipe-separated field
    cells = [c.strip() for c in row_lines[0].split("|")[1:-1]]
    assert cells[2] == "-"


def test_render_ivive_audit_narrative():
    data = _make_data(
        ivive_summary={
            "clint_liver_L_h": 348.75,
            "cl_renal_L_h": 0.23,
            "fu_b": 0.071,
            "liver_model": "well_stirred",
            "compound_type": "base",
        },
        pk_params={"cl_apparent": 24.1},
    )
    out = _render_ivive_audit(data)
    assert "## 4. IVIVE" in out
    assert "well_stirred" in out
    assert "348.75" in out or "348.8" in out or "3.49e+02" in out
    assert "0.071" in out or "7.1e-02" in out
    assert "24.1" in out  # cl_apparent


def test_render_ivive_audit_handles_missing_metadata():
    data = _make_data(ivive_summary={}, pk_params={})
    out = _render_ivive_audit(data)
    assert "## 4. IVIVE" in out
    # Should not crash; contains placeholder
    assert "-" in out or "not available" in out.lower()
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/unit/test_narrative.py -v -k "adme or ivive"`
Expected: FAIL — `ImportError: cannot import name '_render_adme_table'`.

- [ ] **Step 3: Write minimal implementation**

Append to `src/charon/report/narrative.py`:

```python
_ADME_ROW_ORDER: tuple[str, ...] = (
    "logp",
    "pka_acid",
    "pka_base",
    "solubility_ug_ml",
    "fu_p",
    "fu_inc",
    "bp_ratio",
    "clint_uL_min_mg",
    "papp_nm_s",
    "peff_cm_s",
    "herg_ic50_uM",
    "clrenal_L_h",
)


def _ci_cell(entry: dict) -> str:
    lo = entry.get("ci_lower")
    hi = entry.get("ci_upper")
    if lo is None or hi is None:
        return "-"
    return f"[{format_value(lo)}, {format_value(hi)}]"


def _render_adme_table(data: ReportData) -> str:
    lines = ["## 3. ADME Predictions", ""]
    if not data.properties:
        lines.append("*No ADME properties recorded.*")
        return "\n".join(lines)

    lines.append("| Property | Value | 90% CI | Unit | Source | Flag |")
    lines.append("| --- | --- | --- | --- | --- | --- |")
    for key in _ADME_ROW_ORDER:
        entry = data.properties.get(key)
        if entry is None:
            continue
        lines.append(
            "| {name} | {val} | {ci} | {unit} | {src} | {flag} |".format(
                name=key,
                val=format_value(entry.get("value")),
                ci=_ci_cell(entry),
                unit=entry.get("unit") or "-",
                src=entry.get("source") or "-",
                flag=entry.get("flag") or "-",
            )
        )
    return "\n".join(lines)


def _render_ivive_audit(data: ReportData) -> str:
    lines = ["## 4. IVIVE & Hepatic Clearance", ""]
    ivive = data.ivive_summary or {}
    if not ivive:
        lines.append("*IVIVE summary not available in run metadata.*")
        return "\n".join(lines)

    model = ivive.get("liver_model", "-")
    clint_liver = ivive.get("clint_liver_L_h")
    cl_renal = ivive.get("cl_renal_L_h")
    fu_b = ivive.get("fu_b")
    cl_app = data.pk_params.get("cl_apparent")

    lines.append(
        f"CLint was scaled via the **{model}** model. Whole-liver "
        f"CLint = {format_value(clint_liver)} L/h. "
        f"fu_b = {format_value(fu_b)}. "
        f"Renal CL = {format_value(cl_renal)} L/h. "
        f"In-vivo apparent CL from simulation = "
        f"{format_value(cl_app)} L/h."
    )
    lines.append("")
    lines.append(f"- **Liver model:** {model}")
    lines.append(f"- **CLint (whole liver):** {format_value(clint_liver)} L/h")
    lines.append(f"- **fu_b:** {format_value(fu_b)}")
    lines.append(f"- **Renal CL:** {format_value(cl_renal)} L/h")
    lines.append(f"- **In-vivo apparent CL:** {format_value(cl_app)} L/h")
    return "\n".join(lines)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/unit/test_narrative.py -v`
Expected: 15 passed.

- [ ] **Step 5: Commit**

```bash
git add src/charon/report/narrative.py tests/unit/test_narrative.py
git commit -m "feat(report): ADME table and IVIVE audit renderers"
```

---

## Task 6: Narrative — PK results and dose projection (sections 5–6)

**Files:**
- Modify: `src/charon/report/narrative.py`
- Test: `tests/unit/test_narrative.py`

- [ ] **Step 1: Write the failing test**

Append to `tests/unit/test_narrative.py`:

```python
from charon.report.narrative import _render_pk_results, _render_dose_projection


def test_render_pk_results_parameters_table():
    data = _make_data(
        pk_params={
            "cmax": 120.0,
            "tmax": 1.0,
            "auc_0_inf": 500.0,
            "auc_0_24": 480.0,
            "half_life": 3.5,
            "cl_apparent": 24.1,
            "vss": 95.0,
            "bioavailability": 0.44,
            "fa": 0.95,
            "fg": 0.57,
            "fh": 0.82,
        },
        pk_table=[
            {"time_h": 0.0, "cp_plasma_ug_L": 0.0, "cp_blood_ug_L": 0.0},
            {"time_h": 1.0, "cp_plasma_ug_L": 120.0, "cp_blood_ug_L": 108.0},
            {"time_h": 4.0, "cp_plasma_ug_L": 60.0, "cp_blood_ug_L": 54.0},
        ],
    )
    out = _render_pk_results(data)
    assert "## 5. PK Simulation Results" in out
    assert "120" in out  # Cmax
    assert "0.44" in out or "4.4e-01" in out  # F
    assert "Fa" in out and "Fg" in out and "Fh" in out
    # Cp-time table
    assert "| Time" in out
    assert "1" in out and "4" in out


def test_render_pk_results_handles_missing_params():
    data = _make_data(pk_params={"cmax": None, "auc_0_inf": None}, pk_table=[])
    out = _render_pk_results(data)
    assert "## 5. PK Simulation Results" in out
    # Dash for missing values
    assert "-" in out


def test_render_dose_projection_no_rec():
    data = _make_data()
    out = _render_dose_projection(data)
    assert "## 6. FIH Dose Projection" in out
    assert "not run" in out.lower() or "no " in out.lower()


def test_render_dose_projection_hed_only():
    data = _make_data(
        dose_recommendation={
            "mrsd_mg": 56.45,
            "limiting_method": "hed",
            "route": "oral",
            "safety_factor": 10.0,
            "salt_factor": 1.0,
            "rationale": "HED: 56.45 mg\nLimiting: hed",
            "hed": {"mrsd_mg": 56.45},
            "mabel": None,
            "pad": None,
        }
    )
    out = _render_dose_projection(data)
    assert "## 6. FIH Dose Projection" in out
    assert "| Method" in out
    # Three rows: HED, MABEL, PAD
    assert "HED" in out
    assert "MABEL" in out
    assert "PAD" in out
    assert "insufficient inputs" in out.lower()
    # Rationale present
    assert "Limiting: hed" in out
    # Limiting method marked
    assert "**56.45**" in out or "**56.5**" in out or "**5.6e+01**" in out
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/unit/test_narrative.py -v -k "pk_results or dose_projection"`
Expected: FAIL — `ImportError: cannot import name '_render_pk_results'`.

- [ ] **Step 3: Write minimal implementation**

Append to `src/charon/report/narrative.py`:

```python
_PK_PARAM_LABELS: tuple[tuple[str, str, str], ...] = (
    ("cmax", "Cmax", "ug/L"),
    ("tmax", "Tmax", "h"),
    ("auc_0_inf", "AUC(0-inf)", "ug*h/L"),
    ("auc_0_24", "AUC(0-24)", "ug*h/L"),
    ("half_life", "t1/2", "h"),
    ("cl_apparent", "CL apparent", "L/h"),
    ("vss", "Vss", "L"),
    ("bioavailability", "F", "-"),
    ("fa", "Fa", "-"),
    ("fg", "Fg", "-"),
    ("fh", "Fh", "-"),
)


def _render_pk_results(data: ReportData) -> str:
    lines = [
        "## 5. PK Simulation Results",
        "",
        f"Route: **{data.route}**  ·  "
        f"Dose: **{format_value(data.dose_mg)} mg**  ·  "
        f"Duration: **{format_value(data.duration_h)} h**",
        "",
        "### PK Parameters",
        "",
        "| Parameter | Value | Unit |",
        "| --- | --- | --- |",
    ]
    for key, label, unit in _PK_PARAM_LABELS:
        val = data.pk_params.get(key)
        lines.append(f"| {label} | {format_value(val)} | {unit} |")

    lines.append("")
    lines.append("### Concentration-Time Profile (canonical timepoints)")
    lines.append("")
    if not data.pk_table:
        lines.append("*No time-course data available.*")
        return "\n".join(lines)

    lines.append("| Time (h) | Cp plasma (ug/L) | Cp blood (ug/L) |")
    lines.append("| --- | --- | --- |")
    for row in data.pk_table:
        lines.append(
            f"| {format_value(row['time_h'])} "
            f"| {format_value(row['cp_plasma_ug_L'])} "
            f"| {format_value(row['cp_blood_ug_L'])} |"
        )
    return "\n".join(lines)


def _dose_method_row(
    name_display: str,
    method_key: str,
    rec: dict,
    limiting: str,
) -> str:
    sub = rec.get(method_key)
    if sub is None:
        return f"| {name_display} | - | - | insufficient inputs |"
    mrsd = sub.get("mrsd_mg")
    val_cell = format_value(mrsd)
    if method_key == limiting:
        val_cell = f"**{val_cell}**"
    inputs = {
        "hed": f"NOAEL={format_value(sub.get('noael_mg_kg'))} mg/kg "
               f"({sub.get('species', '-')})",
        "mabel": f"Kd={format_value(sub.get('target_kd_nM'))} nM",
        "pad": f"Ceff={format_value(sub.get('target_ceff_nM'))} nM",
    }.get(method_key, "-")
    return f"| {name_display} | {val_cell} | {inputs} | computed |"


def _render_dose_projection(data: ReportData) -> str:
    lines = ["## 6. FIH Dose Projection", ""]
    rec = data.dose_recommendation
    if rec is None:
        lines.append(
            "*FIH dose projection was not run for this pipeline execution.*"
        )
        return "\n".join(lines)

    limiting = str(rec.get("limiting_method", ""))
    lines.append("| Method | MRSD (mg) | Inputs | Status |")
    lines.append("| --- | --- | --- | --- |")
    lines.append(_dose_method_row("HED", "hed", rec, limiting))
    lines.append(_dose_method_row("MABEL", "mabel", rec, limiting))
    lines.append(_dose_method_row("PAD", "pad", rec, limiting))
    lines.append("")
    lines.append(
        f"Safety factor: **{format_value(rec.get('safety_factor'))}**  ·  "
        f"Salt factor: **{format_value(rec.get('salt_factor'))}**  ·  "
        f"Limiting method: **{limiting.upper()}**"
    )
    rationale = rec.get("rationale") or ""
    if rationale:
        lines.append("")
        lines.append("```")
        lines.append(rationale)
        lines.append("```")
    return "\n".join(lines)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/unit/test_narrative.py -v`
Expected: 19 passed.

- [ ] **Step 5: Commit**

```bash
git add src/charon/report/narrative.py tests/unit/test_narrative.py
git commit -m "feat(report): PK results and dose projection renderers"
```

---

## Task 7: Narrative — uncertainty, limitations, appendix + render_report()

**Files:**
- Modify: `src/charon/report/narrative.py`
- Test: `tests/unit/test_narrative.py`

- [ ] **Step 1: Write the failing test**

Append to `tests/unit/test_narrative.py`:

```python
from charon.report.narrative import (
    _render_appendix,
    _render_limitations,
    _render_uncertainty,
    render_report,
)


def test_render_uncertainty_none_returns_empty():
    data = _make_data(uncertainty=None)
    out = _render_uncertainty(data)
    assert out == "" or out.strip() == ""


def test_render_uncertainty_populated():
    data = _make_data(
        uncertainty={
            "point_estimate_mg": 50.0,
            "ci_90_lower_mg": 20.0,
            "ci_90_upper_mg": 120.0,
            "ci_ratio": 6.0,
            "confidence": "MEDIUM",
            "n_samples": 100,
            "n_successful": 95,
            "convergence_met": True,
            "sensitivity": {"clint": 0.7, "fu_p": 0.2, "logp": 0.1},
            "limiting_parameter": "clint",
            "recommendation": "Experimental clint measurement would narrow CI by ~70%",
            "r_squared": 0.85,
        }
    )
    out = _render_uncertainty(data)
    assert "## 7. Uncertainty Analysis" in out
    assert "50" in out and "20" in out and "120" in out
    assert "MEDIUM" in out
    assert "clint" in out
    assert "70" in out  # 70.0%
    assert "0.85" in out  # R²


def test_render_uncertainty_low_r_squared_warning():
    data = _make_data(
        uncertainty={
            "point_estimate_mg": 50.0,
            "ci_90_lower_mg": 20.0,
            "ci_90_upper_mg": 120.0,
            "ci_ratio": 6.0,
            "confidence": "MEDIUM",
            "n_samples": 100,
            "n_successful": 95,
            "convergence_met": True,
            "sensitivity": {"clint": 0.5},
            "limiting_parameter": "clint",
            "recommendation": "",
            "r_squared": 0.5,
        }
    )
    out = _render_uncertainty(data)
    assert "non-linear" in out.lower() or "nonlinear" in out.lower() or "warn" in out.lower()


def test_render_limitations_always_has_boilerplate():
    data = _make_data()
    out = _render_limitations(data)
    assert "## 8. Limitations" in out
    assert "well-stirred" in out.lower() or "well stirred" in out.lower()
    assert "IR" in out or "immediate-release" in out.lower()


def test_render_limitations_appends_warnings():
    data = _make_data(warnings=["Applicability domain: LOW"])
    out = _render_limitations(data)
    assert "Applicability domain: LOW" in out


def test_render_limitations_appends_property_flags():
    data = _make_data(
        properties={
            "clint_uL_min_mg": {
                "value": 93.0,
                "ci_lower": None,
                "ci_upper": None,
                "unit": "uL/min/mg",
                "source": "ml_ensemble",
                "flag": "clint_tier2_ml",
                "method": None,
            }
        }
    )
    out = _render_limitations(data)
    assert "clint_tier2_ml" in out


def test_render_appendix_has_metadata_and_version():
    data = _make_data(metadata={"species": "human", "solver_method": "BDF"})
    out = _render_appendix(data)
    assert "## 9. Appendix" in out
    assert "0.1.0" in out  # charon_version
    assert "BDF" in out
    assert "2026-04-15" in out


def test_render_report_contains_all_sections():
    data = _make_data(
        properties={
            "logp": {
                "value": 3.89,
                "ci_lower": 3.4,
                "ci_upper": 4.3,
                "unit": "log",
                "source": "ml_ensemble",
                "flag": None,
                "method": None,
            }
        },
        ivive_summary={"liver_model": "well_stirred", "fu_b": 0.071},
        pk_params={"cmax": 120.0, "cl_apparent": 24.1, "auc_0_inf": 500.0},
        pk_table=[{"time_h": 1.0, "cp_plasma_ug_L": 120.0, "cp_blood_ug_L": 108.0}],
        dose_recommendation={
            "mrsd_mg": 56.45,
            "limiting_method": "hed",
            "route": "oral",
            "safety_factor": 10.0,
            "salt_factor": 1.0,
            "rationale": "ok",
            "hed": {"mrsd_mg": 56.45},
            "mabel": None,
            "pad": None,
        },
    )
    out = render_report(data)
    for section in (
        "# FIH Dose Rationale Report",
        "## 1. Executive Summary",
        "## 2. Compound Profile",
        "## 3. ADME Predictions",
        "## 4. IVIVE",
        "## 5. PK Simulation Results",
        "## 6. FIH Dose Projection",
        "## 8. Limitations",
        "## 9. Appendix",
    ):
        assert section in out, f"missing section: {section}"
    # No uncertainty section when uncertainty is None
    assert "## 7. Uncertainty Analysis" not in out


def test_render_report_includes_uncertainty_when_present():
    data = _make_data(
        uncertainty={
            "point_estimate_mg": 50.0,
            "ci_90_lower_mg": 20.0,
            "ci_90_upper_mg": 120.0,
            "ci_ratio": 6.0,
            "confidence": "MEDIUM",
            "n_samples": 100,
            "n_successful": 95,
            "convergence_met": True,
            "sensitivity": {"clint": 0.7},
            "limiting_parameter": "clint",
            "recommendation": "",
            "r_squared": 0.9,
        }
    )
    out = render_report(data)
    assert "## 7. Uncertainty Analysis" in out


def test_render_report_no_raw_none():
    data = _make_data()
    out = render_report(data)
    # We should not have literal "None" appearing as a value
    # (the format_value helper replaces None with "-").
    # Some legitimate text may contain the word "none", so only check
    # for standalone occurrences inside table cells.
    assert "| None |" not in out
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/unit/test_narrative.py -v -k "uncertainty or limitations or appendix or render_report"`
Expected: FAIL — `_render_uncertainty`, `render_report` etc. do not exist.

- [ ] **Step 3: Write minimal implementation**

Append to `src/charon/report/narrative.py`:

```python
def _render_uncertainty(data: ReportData) -> str:
    unc = data.uncertainty
    if unc is None:
        return ""
    lines = ["## 7. Uncertainty Analysis", ""]
    point = format_value(unc.get("point_estimate_mg"))
    lo = format_value(unc.get("ci_90_lower_mg"))
    hi = format_value(unc.get("ci_90_upper_mg"))
    conf = unc.get("confidence", "?")
    lines.append(
        f"**Dose:** {point} mg  ·  "
        f"**90% CI:** [{lo} – {hi}] mg  ·  "
        f"**Confidence:** {conf}"
    )
    lines.append("")

    sens = unc.get("sensitivity") or {}
    if sens:
        lines.append("### Parameter Sensitivity (SRC²)")
        lines.append("")
        lines.append("| Rank | Parameter | Importance |")
        lines.append("| --- | --- | --- |")
        sorted_sens = sorted(sens.items(), key=lambda kv: kv[1], reverse=True)
        for rank, (name, imp) in enumerate(sorted_sens, start=1):
            pct = float(imp) * 100.0
            lines.append(f"| {rank} | {name} | {pct:.1f}% |")
        lines.append("")

    r2 = unc.get("r_squared")
    if r2 is not None:
        lines.append(f"SRC regression R² = {format_value(r2)}")
        if float(r2) < 0.7:
            lines.append(
                "> **Warning:** R² < 0.7 indicates substantial nonlinear "
                "behavior; SRC² values underestimate interaction effects and "
                "should be treated as first-order approximations only."
            )
        lines.append("")

    rec_text = unc.get("recommendation") or ""
    if rec_text:
        lines.append(f"> {rec_text}")
        lines.append("")

    lines.append(
        f"Samples: {unc.get('n_successful', 0)} successful / "
        f"{unc.get('n_samples', 0)} total  ·  "
        f"converged: {unc.get('convergence_met', False)}"
    )
    return "\n".join(lines)


_LIMITATIONS_BOILERPLATE = (
    "- Hepatic clearance uses the well-stirred liver model by default; "
    "parallel-tube or dispersion models may be more appropriate for "
    "high-extraction compounds.\n"
    "- Oral absorption assumes an immediate-release (IR) formulation; "
    "modified-release kinetics are not modeled in Phase A.\n"
    "- Conformal 90% prediction intervals provide marginal coverage, "
    "not conditional coverage. Actual coverage on out-of-domain compounds "
    "may be substantially lower than 90%.\n"
    "- Rodgers & Rowland Kp prediction may overestimate tissue partitioning "
    "for highly-bound weak bases and lipophilic acids.\n"
    "- Dose projection uses apparent PK (CL/F for oral, CL for IV). No "
    "separate bioavailability correction is applied."
)


def _render_limitations(data: ReportData) -> str:
    lines = ["## 8. Limitations & Caveats", "", _LIMITATIONS_BOILERPLATE, ""]
    extras: list[str] = []
    for w in data.warnings:
        extras.append(f"- {w}")
    flag_seen: set[str] = set()
    for key, entry in data.properties.items():
        flag = entry.get("flag") if isinstance(entry, dict) else None
        if flag and flag not in flag_seen:
            flag_seen.add(flag)
            extras.append(f"- `{key}`: flag = `{flag}`")
    if extras:
        lines.append("Additional run-specific notes:")
        lines.append("")
        lines.extend(extras)
    return "\n".join(lines)


def _render_appendix(data: ReportData) -> str:
    lines = ["## 9. Appendix", ""]
    lines.append(f"- **Charon version:** {data.charon_version}")
    lines.append(f"- **Timestamp:** {data.timestamp}")
    lines.append("")
    lines.append("```yaml")
    # Render metadata as sorted key: value lines (stable ordering).
    for key in sorted(data.metadata.keys()):
        val = data.metadata[key]
        lines.append(f"{key}: {val}")
    lines.append("```")
    return "\n".join(lines)


_SECTION_RENDERERS = (
    _render_header,
    _render_executive_summary,
    _render_compound_profile,
    _render_adme_table,
    _render_ivive_audit,
    _render_pk_results,
    _render_dose_projection,
    _render_uncertainty,
    _render_limitations,
    _render_appendix,
)


def render_report(data: ReportData) -> str:
    """Render a full Markdown report from a :class:`ReportData`."""
    parts: list[str] = []
    for fn in _SECTION_RENDERERS:
        chunk = fn(data)
        if chunk:
            parts.append(chunk)
    return "\n\n".join(parts) + "\n"
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/unit/test_narrative.py -v`
Expected: 28 passed.

- [ ] **Step 5: Commit**

```bash
git add src/charon/report/narrative.py tests/unit/test_narrative.py
git commit -m "feat(report): uncertainty, limitations, appendix + render_report"
```

---

## Task 8: Export module + report/__init__.py

**Files:**
- Create: `src/charon/report/export.py`
- Modify: `src/charon/report/__init__.py`
- Test: `tests/unit/test_export.py`

- [ ] **Step 1: Write the failing test**

Create `tests/unit/test_export.py`:

```python
"""Unit tests for charon.report.export."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from charon.report.collector import ReportData
from charon.report.export import (
    export_json,
    export_markdown,
    export_report,
)


def _make_data() -> ReportData:
    return ReportData(
        compound_name="test",
        smiles="CCO",
        molecular_weight=46.07,
        source="predicted",
        compound_type="neutral",
        properties={
            "logp": {
                "value": -0.31,
                "ci_lower": None,
                "ci_upper": None,
                "unit": "log",
                "source": "ml_ensemble",
                "flag": None,
                "method": None,
            }
        },
        ivive_summary={"liver_model": "well_stirred", "fu_b": 0.95},
        pk_params={"cmax": 10.0, "auc_0_inf": 50.0},
        pk_table=[{"time_h": 0.0, "cp_plasma_ug_L": 10.0, "cp_blood_ug_L": 9.0}],
        route="iv_bolus",
        dose_mg=10.0,
        duration_h=24.0,
        dose_recommendation=None,
        uncertainty=None,
        warnings=[],
        metadata={"species": "human", "solver_nfev": 42},
        timestamp="2026-04-15T00:00:00+00:00",
        charon_version="0.1.0",
    )


def test_export_markdown_writes_file(tmp_path: Path):
    data = _make_data()
    out = tmp_path / "report.md"
    path = export_markdown(data, out)
    assert path == out.resolve()
    text = out.read_text()
    assert "# FIH Dose Rationale Report" in text
    assert "test" in text
    # At least 9 sections
    assert text.count("## ") >= 8  # Sections 1-6, 8, 9 (Section 7 skipped)


def test_export_json_writes_valid_json(tmp_path: Path):
    data = _make_data()
    out = tmp_path / "report.json"
    path = export_json(data, out)
    assert path == out.resolve()
    parsed = json.loads(out.read_text())
    assert parsed["compound_name"] == "test"
    assert parsed["smiles"] == "CCO"
    assert parsed["pk_params"]["cmax"] == 10.0
    assert parsed["ivive_summary"]["liver_model"] == "well_stirred"


def test_export_json_has_all_top_level_keys(tmp_path: Path):
    data = _make_data()
    out = tmp_path / "r.json"
    export_json(data, out)
    parsed = json.loads(out.read_text())
    expected_keys = {
        "compound_name", "smiles", "molecular_weight", "source",
        "compound_type", "properties", "ivive_summary", "pk_params",
        "pk_table", "route", "dose_mg", "duration_h",
        "dose_recommendation", "uncertainty", "warnings", "metadata",
        "timestamp", "charon_version",
    }
    assert expected_keys.issubset(parsed.keys())


def test_export_json_full_profile_optional(tmp_path: Path):
    data = _make_data()
    out = tmp_path / "r.json"
    export_json(
        data,
        out,
        include_full_profile=True,
        full_profile={
            "time_h": [0.0, 1.0, 2.0],
            "cp_plasma_ug_L": [10.0, 8.0, 6.0],
            "cp_blood_ug_L": [9.0, 7.2, 5.4],
        },
    )
    parsed = json.loads(out.read_text())
    assert "full_profile" in parsed
    assert parsed["full_profile"]["time_h"] == [0.0, 1.0, 2.0]


def test_export_report_writes_both_md_and_json(tmp_path: Path):
    data = _make_data()
    md_path, json_path = export_report(data, tmp_path / "run.md")
    assert md_path.exists()
    assert json_path.exists()
    assert md_path.suffix == ".md"
    assert json_path.suffix == ".json"
    assert md_path.stem == json_path.stem == "run"


def test_export_report_handles_no_suffix(tmp_path: Path):
    data = _make_data()
    md_path, json_path = export_report(data, tmp_path / "run")
    assert md_path.exists()
    assert json_path.exists()
    assert md_path.name == "run.md"
    assert json_path.name == "run.json"


def test_export_json_handles_numpy_values(tmp_path: Path):
    data = _make_data()
    # Simulate a numpy value sneaking through (e.g. np.float64 inside metadata)
    data_with_np = ReportData(
        **{
            **{k: getattr(data, k) for k in data.__dataclass_fields__},
            "metadata": {"solver_nfev": np.int64(42), "cl": np.float64(1.23)},
        }
    )
    out = tmp_path / "np.json"
    export_json(data_with_np, out)
    parsed = json.loads(out.read_text())
    assert parsed["metadata"]["solver_nfev"] == 42
    assert parsed["metadata"]["cl"] == pytest.approx(1.23)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/unit/test_export.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'charon.report.export'`.

- [ ] **Step 3: Write minimal implementation**

Create `src/charon/report/export.py`:

```python
"""Sprint 6 — Report file exporters.

Writes a :class:`ReportData` to a Markdown file and a sibling JSON file.
No external dependencies: uses the ``json`` stdlib module and a small
``_json_default`` helper to handle numpy scalars.
"""

from __future__ import annotations

import json
from dataclasses import asdict
from pathlib import Path
from typing import Any

import numpy as np

from charon.report.collector import ReportData
from charon.report.narrative import render_report


def _json_default(obj: Any) -> Any:
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.floating):
        return float(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, Path):
        return str(obj)
    raise TypeError(f"Unsupported type for JSON serialization: {type(obj)!r}")


def export_markdown(data: ReportData, path: Path) -> Path:
    """Render *data* as Markdown and write to *path*.

    Returns the resolved absolute path.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(render_report(data), encoding="utf-8")
    return path.resolve()


def export_json(
    data: ReportData,
    path: Path,
    *,
    include_full_profile: bool = False,
    full_profile: dict | None = None,
) -> Path:
    """Serialise *data* to JSON at *path*.  Returns the resolved path."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = asdict(data)
    if include_full_profile and full_profile is not None:
        payload["full_profile"] = full_profile
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=False, default=_json_default),
        encoding="utf-8",
    )
    return path.resolve()


def export_report(
    data: ReportData,
    output: Path,
    *,
    full_profile: dict | None = None,
) -> tuple[Path, Path]:
    """Write both the Markdown and JSON siblings of a report.

    If *output* ends in ``.md`` the JSON path replaces the suffix with
    ``.json``; otherwise both ``.md`` and ``.json`` are appended to the
    given path.
    """
    output = Path(output)
    if output.suffix.lower() == ".md":
        md_path = output
        json_path = output.with_suffix(".json")
    else:
        md_path = output.with_name(output.name + ".md")
        json_path = output.with_name(output.name + ".json")

    md_resolved = export_markdown(data, md_path)
    json_resolved = export_json(
        data,
        json_path,
        include_full_profile=full_profile is not None,
        full_profile=full_profile,
    )
    return md_resolved, json_resolved
```

Create/update `src/charon/report/__init__.py`:

```python
"""Charon Sprint 6 — Report collection, rendering, and export."""

from charon.report.collector import ReportData, collect
from charon.report.export import (
    export_json,
    export_markdown,
    export_report,
)
from charon.report.narrative import render_report

__all__ = [
    "ReportData",
    "collect",
    "export_json",
    "export_markdown",
    "export_report",
    "render_report",
]
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/unit/test_export.py tests/unit/test_collector.py tests/unit/test_narrative.py -v`
Expected: all passing.

- [ ] **Step 5: Commit**

```bash
git add src/charon/report/export.py src/charon/report/__init__.py tests/unit/test_export.py
git commit -m "feat(report): export module for Markdown + JSON + __init__ surface"
```

---

## Task 9: CLI scaffolding — main dispatcher + `predict` and `simulate` subcommands

**Files:**
- Create: `src/charon/cli/main.py`
- Modify: `src/charon/cli/__init__.py`
- Modify: `pyproject.toml`
- Test: `tests/unit/test_cli.py`

Goal: argparse dispatcher, `predict` and `simulate` subcommands end-to-end, `pyproject.toml` entry point.

- [ ] **Step 1: Write the failing test**

Create `tests/unit/test_cli.py`:

```python
"""Unit tests for charon.cli.main (invoked directly for speed)."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from charon.cli.main import main


def test_main_help_lists_subcommands(capsys):
    with pytest.raises(SystemExit) as exc:
        main(["--help"])
    assert exc.value.code == 0
    out = capsys.readouterr().out
    for sub in ("predict", "simulate", "translate", "recommend", "report"):
        assert sub in out


def test_main_no_args_exits_nonzero(capsys):
    rc = main([])
    assert rc != 0


def test_predict_subcommand_basic(capsys, monkeypatch):
    """`charon predict <smiles>` calls predict_properties and prints a table."""
    from charon.core.schema import (
        BindingProperties,
        CompoundProperties,
        PhysicochemicalProperties,
        PredictedProperty,
    )

    def fake_predict(smiles, **_):
        return CompoundProperties(
            physicochemical=PhysicochemicalProperties(
                logp=PredictedProperty(value=0.5, source="ml_ensemble", unit="log"),
            ),
            binding=BindingProperties(
                fu_p=PredictedProperty(
                    value=0.2, source="ml_ensemble", unit="fraction"
                ),
            ),
        )

    monkeypatch.setattr("charon.cli.main.predict_properties", fake_predict)
    rc = main(["predict", "CCO"])
    assert rc == 0
    out = capsys.readouterr().out
    assert "logp" in out
    assert "0.5" in out
    assert "fu_p" in out


def test_predict_json_flag(capsys, monkeypatch):
    from charon.core.schema import CompoundProperties

    monkeypatch.setattr(
        "charon.cli.main.predict_properties",
        lambda smiles, **_: CompoundProperties(),
    )
    rc = main(["predict", "CCO", "--json"])
    assert rc == 0
    out = capsys.readouterr().out
    parsed = json.loads(out)
    assert isinstance(parsed, dict)


def test_predict_invalid_smiles(capsys, monkeypatch):
    def boom(smiles, **_):
        raise ValueError("bad SMILES")

    monkeypatch.setattr("charon.cli.main.predict_properties", boom)
    rc = main(["predict", "not-a-smiles"])
    assert rc == 1
    err = capsys.readouterr().err
    assert "Invalid" in err or "bad SMILES" in err


def test_simulate_subcommand_runs_pipeline(capsys, monkeypatch):
    """`charon simulate` runs Pipeline.from_smiles(...).run() and prints PK."""
    import numpy as np

    from charon.core.schema import (
        CompoundConfig,
        CompoundProperties,
        PKParameters,
    )
    from charon.pbpk.solver import SimulationResult
    from charon.pipeline import PipelineResult

    def fake_from_smiles(smiles, **kwargs):
        class _FakePipeline:
            def run(self):
                t = np.array([0.0, 1.0, 2.0])
                cp = np.array([5.0, 2.5, 1.25])
                return PipelineResult(
                    compound=CompoundConfig(
                        name=smiles,
                        smiles=smiles,
                        molecular_weight=46.07,
                        source="predicted",
                        properties=CompoundProperties(),
                    ),
                    pk_parameters=PKParameters(cmax=5.0, auc_0_inf=12.0),
                    time_h=t,
                    cp_plasma=cp,
                    cp_blood=cp,
                    simulation=SimulationResult(
                        time_h=t,
                        cp_plasma=cp,
                        cp_blood=cp,
                        tissue_trajectories={},
                        solver_method="BDF",
                        solver_nfev=7,
                    ),
                    metadata={
                        "route": kwargs.get("route"),
                        "dose_mg": kwargs.get("dose_mg"),
                        "duration_h": kwargs.get("duration_h", 72.0),
                        "liver_model": "well_stirred",
                    },
                )
        return _FakePipeline()

    monkeypatch.setattr(
        "charon.cli.main.Pipeline.from_smiles", staticmethod(fake_from_smiles)
    )
    rc = main(["simulate", "CCO", "--route", "iv_bolus", "--dose", "10"])
    assert rc == 0
    out = capsys.readouterr().out
    assert "Cmax" in out or "cmax" in out
    assert "5" in out
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/unit/test_cli.py -v`
Expected: FAIL — `charon.cli.main` has no `main` attribute.

- [ ] **Step 3: Write minimal implementation**

Create `src/charon/cli/main.py`:

```python
"""Charon CLI — five subcommands sharing a single dispatcher.

Subcommands:

    predict    Layer 0 + Layer 1 only (ADME prediction)
    simulate   Full pipeline up to Layer 2 (PBPK)
    translate  Pipeline + Layer 3 (FIH dose projection)
    recommend  translate + optional Layer 4 (uncertainty quantification)
    report     recommend + write Markdown/JSON report files

Invocation:

    charon <subcommand> <SMILES> [options...]
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import asdict
from typing import Sequence

from charon.core.schema import DoseProjectionConfig, UncertaintyConfig
from charon.pipeline import Pipeline
from charon.predict import predict_properties
from charon.report import (
    ReportData,
    collect,
    export_report,
    render_report,
)


# ---------------------------------------------------------------------------
# Small helpers (printers and error exits)
# ---------------------------------------------------------------------------


def _fail(msg: str, code: int = 1) -> int:
    print(f"error: {msg}", file=sys.stderr)
    return code


def _print_table(rows: list[list[str]], headers: list[str]) -> None:
    print("| " + " | ".join(headers) + " |")
    print("| " + " | ".join("---" for _ in headers) + " |")
    for row in rows:
        print("| " + " | ".join(str(c) for c in row) + " |")


def _format_value(v: float | int | None) -> str:
    if v is None:
        return "-"
    a = abs(float(v))
    if a == 0.0:
        return "0"
    if a >= 1000.0 or a < 0.01:
        return f"{float(v):.3e}"
    return f"{float(v):.4g}"


# ---------------------------------------------------------------------------
# Subcommand: predict
# ---------------------------------------------------------------------------


def _cmd_predict(args: argparse.Namespace) -> int:
    try:
        props = predict_properties(args.smiles)
    except ValueError as e:
        return _fail(f"Invalid SMILES: {e}")
    except Exception as e:  # pragma: no cover - defensive
        return _fail(f"{type(e).__name__}: {e}", code=2)

    if args.json:
        print(json.dumps(props.model_dump(), indent=2, default=str))
        return 0

    rows: list[list[str]] = []
    for category_name in (
        "physicochemical",
        "binding",
        "metabolism",
        "permeability",
        "renal",
        "safety",
    ):
        cat = getattr(props, category_name, None)
        if cat is None:
            continue
        for attr_name in cat.model_fields:
            val = getattr(cat, attr_name, None)
            if val is None:
                continue
            if hasattr(val, "value"):
                rows.append(
                    [
                        attr_name,
                        _format_value(val.value),
                        val.unit or "-",
                        val.source or "-",
                    ]
                )
    if not rows:
        print("(no properties predicted)")
        return 0
    _print_table(rows, ["Property", "Value", "Unit", "Source"])
    return 0


# ---------------------------------------------------------------------------
# Subcommand: simulate
# ---------------------------------------------------------------------------


def _build_pipeline_from_smiles(args: argparse.Namespace, **extra) -> Pipeline:
    return Pipeline.from_smiles(
        args.smiles,
        route=args.route,
        dose_mg=args.dose,
        species=args.species,
        duration_h=args.duration,
        infusion_duration_h=getattr(args, "infusion_duration", 0.0),
        liver_model=args.liver_model,
        compound_name=args.compound_name or args.smiles,
        **extra,
    )


def _print_pk_summary(data: ReportData) -> None:
    rows = []
    labels = [
        ("cmax", "Cmax", "ug/L"),
        ("tmax", "Tmax", "h"),
        ("auc_0_inf", "AUC(0-inf)", "ug*h/L"),
        ("half_life", "t1/2", "h"),
        ("cl_apparent", "CL apparent", "L/h"),
        ("vss", "Vss", "L"),
        ("bioavailability", "F", "-"),
        ("fa", "Fa", "-"),
        ("fg", "Fg", "-"),
        ("fh", "Fh", "-"),
    ]
    for key, label, unit in labels:
        val = data.pk_params.get(key)
        rows.append([label, _format_value(val), unit])
    _print_table(rows, ["Parameter", "Value", "Unit"])


def _cmd_simulate(args: argparse.Namespace) -> int:
    try:
        pipe = _build_pipeline_from_smiles(args)
        result = pipe.run()
    except ValueError as e:
        return _fail(f"Invalid input: {e}")
    except Exception as e:  # pragma: no cover - defensive
        return _fail(f"{type(e).__name__}: {e}", code=2)

    data = collect(result)
    if args.json:
        print(json.dumps(asdict(data), indent=2, default=str))
        return 0
    print(f"Pipeline: {args.smiles}  |  route={args.route}  |  dose={args.dose} mg")
    print()
    _print_pk_summary(data)
    return 0


# ---------------------------------------------------------------------------
# Argument parser
# ---------------------------------------------------------------------------


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="charon",
        description="Charon: SMILES → FIH dose recommendation pipeline.",
    )
    sub = p.add_subparsers(dest="command", metavar="COMMAND")

    # Shared option helper
    def _add_common(sp: argparse.ArgumentParser) -> None:
        sp.add_argument(
            "--species", default="human", choices=["human", "rat", "dog", "monkey"]
        )
        sp.add_argument(
            "--liver-model",
            default="well_stirred",
            choices=["well_stirred", "parallel_tube", "dispersion"],
            dest="liver_model",
        )
        sp.add_argument("--compound-name", default=None, dest="compound_name")
        sp.add_argument("--json", action="store_true", default=False)
        sp.add_argument("-q", "--quiet", action="store_true", default=False)

    # predict
    sp_pred = sub.add_parser("predict", help="ADME prediction (Layer 0+1 only)")
    sp_pred.add_argument("smiles")
    _add_common(sp_pred)
    sp_pred.set_defaults(func=_cmd_predict)

    # simulate
    sp_sim = sub.add_parser(
        "simulate", help="PK simulation (Layer 0+1+2, no dose projection)"
    )
    sp_sim.add_argument("smiles")
    sp_sim.add_argument(
        "--route",
        required=True,
        choices=["iv_bolus", "iv_infusion", "oral"],
    )
    sp_sim.add_argument("--dose", type=float, required=True, help="dose in mg")
    sp_sim.add_argument("--duration", type=float, default=72.0)
    sp_sim.add_argument(
        "--infusion-duration", type=float, default=0.0, dest="infusion_duration"
    )
    _add_common(sp_sim)
    sp_sim.set_defaults(func=_cmd_simulate)

    return p


def main(argv: Sequence[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    if not getattr(args, "command", None):
        parser.print_help(sys.stderr)
        return 1
    func = args.func
    return func(args)


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
```

Update `src/charon/cli/__init__.py`:

```python
"""Charon CLI entry point."""

from charon.cli.main import main

__all__ = ["main"]
```

Update `pyproject.toml` — add after the `[project]` block (before `[project.optional-dependencies]`):

```toml
[project.scripts]
charon = "charon.cli.main:main"
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/unit/test_cli.py -v`
Expected: 6 passed.

- [ ] **Step 5: Commit**

```bash
git add src/charon/cli/main.py src/charon/cli/__init__.py pyproject.toml tests/unit/test_cli.py
git commit -m "feat(cli): argparse dispatcher with predict and simulate subcommands"
```

---

## Task 10: CLI — `translate` and `recommend` subcommands

**Files:**
- Modify: `src/charon/cli/main.py`
- Test: `tests/unit/test_cli.py`

- [ ] **Step 1: Write the failing test**

Append to `tests/unit/test_cli.py`:

```python
def _install_fake_pipeline_with_dose(monkeypatch, *, uncertainty=None):
    """Install a fake Pipeline.from_smiles that returns a PipelineResult with
    a dose_recommendation (and optional uncertainty)."""
    import numpy as np

    from charon.core.schema import CompoundConfig, CompoundProperties, PKParameters
    from charon.pbpk.solver import SimulationResult
    from charon.pipeline import PipelineResult
    from charon.translational.dose_projector import FIHDoseRecommendation
    from charon.translational.hed import HEDResult

    def fake_from_smiles(smiles, **kwargs):
        class _FakePipe:
            def __init__(self):
                self.dose_projection = kwargs.get("dose_projection")
                self.uncertainty = kwargs.get("uncertainty")

            def run(self):
                t = np.array([0.0, 1.0, 2.0])
                cp = np.array([5.0, 2.5, 1.25])
                rec = None
                if self.dose_projection is not None:
                    hed = HEDResult(
                        hed_mg_kg=8.06,
                        mrsd_mg=56.45,
                        species="rat",
                        km_animal=6.2,
                        km_human=37.0,
                        noael_mg_kg=50.0,
                        safety_factor=10.0,
                        body_weight_kg=70.0,
                    )
                    rec = FIHDoseRecommendation(
                        mrsd_mg=56.45,
                        limiting_method="hed",
                        hed=hed,
                        mabel=None,
                        pad=None,
                        safety_factor=10.0,
                        salt_factor=1.0,
                        route=kwargs.get("route", "oral"),
                        rationale="HED 56.45 mg",
                    )
                return PipelineResult(
                    compound=CompoundConfig(
                        name=smiles,
                        smiles=smiles,
                        molecular_weight=46.07,
                        source="predicted",
                        properties=CompoundProperties(),
                    ),
                    pk_parameters=PKParameters(
                        cmax=5.0, auc_0_inf=12.0, half_life=2.0, cl_apparent=0.2
                    ),
                    time_h=t,
                    cp_plasma=cp,
                    cp_blood=cp,
                    simulation=SimulationResult(
                        time_h=t, cp_plasma=cp, cp_blood=cp,
                        tissue_trajectories={}, solver_method="BDF",
                        solver_nfev=7,
                    ),
                    metadata={
                        "route": kwargs.get("route"),
                        "dose_mg": kwargs.get("dose_mg"),
                        "duration_h": kwargs.get("duration_h", 72.0),
                        "liver_model": "well_stirred",
                    },
                    dose_recommendation=rec,
                    uncertainty=uncertainty,
                )
        return _FakePipe()

    monkeypatch.setattr(
        "charon.cli.main.Pipeline.from_smiles", staticmethod(fake_from_smiles)
    )


def test_translate_requires_at_least_one_target(capsys, monkeypatch):
    _install_fake_pipeline_with_dose(monkeypatch)
    rc = main(["translate", "CCO", "--route", "oral", "--dose", "5"])
    # No NOAEL / Kd / Ceff -> should fail
    assert rc != 0
    err = capsys.readouterr().err
    assert (
        "noael" in err.lower() or "target" in err.lower() or "required" in err.lower()
    )


def test_translate_with_noael(capsys, monkeypatch):
    _install_fake_pipeline_with_dose(monkeypatch)
    rc = main(
        [
            "translate",
            "CCO",
            "--route", "oral",
            "--dose", "5",
            "--noael", "50",
            "--noael-species", "rat",
        ]
    )
    assert rc == 0
    out = capsys.readouterr().out
    assert "56.45" in out or "56.5" in out or "5.6" in out
    assert "HED" in out


def test_recommend_without_uncertainty(capsys, monkeypatch):
    _install_fake_pipeline_with_dose(monkeypatch)
    rc = main(
        [
            "recommend",
            "CCO",
            "--route", "oral",
            "--dose", "5",
            "--noael", "50",
            "--noael-species", "rat",
        ]
    )
    assert rc == 0
    out = capsys.readouterr().out
    assert "56.45" in out or "56.5" in out


def test_recommend_with_uncertainty_flag(capsys, monkeypatch):
    from charon.uncertainty.dose_range import UncertaintyResult

    unc = UncertaintyResult(
        point_estimate_mg=56.0,
        ci_90_lower_mg=25.0,
        ci_90_upper_mg=125.0,
        ci_ratio=5.0,
        confidence="MEDIUM",
        n_samples=50,
        n_successful=48,
        convergence_met=False,
        sensitivity={"clint_uL_min_mg": 0.6, "fu_p": 0.3, "logp": 0.1},
        limiting_parameter="clint_uL_min_mg",
        recommendation="Experimental CLint measurement would narrow CI by ~60%",
        r_squared=0.82,
    )
    _install_fake_pipeline_with_dose(monkeypatch, uncertainty=unc)
    rc = main(
        [
            "recommend",
            "CCO",
            "--route", "oral",
            "--dose", "5",
            "--noael", "50",
            "--noael-species", "rat",
            "--uncertainty",
            "--n-samples", "50",
        ]
    )
    assert rc == 0
    out = capsys.readouterr().out
    assert "25" in out and "125" in out
    assert "MEDIUM" in out
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/unit/test_cli.py -v -k "translate or recommend"`
Expected: FAIL — subcommands not registered.

- [ ] **Step 3: Write minimal implementation**

Edit `src/charon/cli/main.py`. Add the two subcommand handlers after `_cmd_simulate`:

```python
# ---------------------------------------------------------------------------
# Subcommand: translate
# ---------------------------------------------------------------------------


def _build_dose_projection(args: argparse.Namespace) -> DoseProjectionConfig | None:
    has_hed = args.noael is not None and args.noael_species is not None
    has_mabel = args.target_kd is not None
    has_pad = args.target_ceff is not None
    if not (has_hed or has_mabel or has_pad):
        return None
    return DoseProjectionConfig(
        noael_mg_kg=args.noael,
        noael_species=args.noael_species,
        target_kd_nM=args.target_kd,
        target_ceff_nM=args.target_ceff,
        safety_factor=args.safety_factor,
        tau_h=args.tau,
        body_weight_kg=args.body_weight,
    )


def _print_dose_recommendation(data: ReportData) -> None:
    rec = data.dose_recommendation
    if rec is None:
        print("(no dose recommendation produced)")
        return
    print(rec.get("rationale", ""))
    print()
    rows: list[list[str]] = []
    for key, label in (("hed", "HED"), ("mabel", "MABEL"), ("pad", "PAD")):
        sub = rec.get(key)
        if sub is None:
            rows.append([label, "-", "insufficient inputs"])
        else:
            mark = " <-" if rec.get("limiting_method") == key else ""
            rows.append([label, _format_value(sub.get("mrsd_mg")), f"computed{mark}"])
    _print_table(rows, ["Method", "MRSD (mg)", "Status"])


def _cmd_translate(args: argparse.Namespace) -> int:
    dp = _build_dose_projection(args)
    if dp is None:
        return _fail(
            "translate requires at least one target: "
            "--noael + --noael-species, --target-kd, or --target-ceff"
        )
    try:
        pipe = _build_pipeline_from_smiles(args, dose_projection=dp)
        result = pipe.run()
    except ValueError as e:
        return _fail(f"Invalid input: {e}")
    except Exception as e:  # pragma: no cover
        return _fail(f"{type(e).__name__}: {e}", code=2)

    data = collect(result)
    if args.json:
        print(json.dumps(asdict(data), indent=2, default=str))
        return 0
    _print_dose_recommendation(data)
    return 0


# ---------------------------------------------------------------------------
# Subcommand: recommend
# ---------------------------------------------------------------------------


def _print_uncertainty_summary(data: ReportData) -> None:
    unc = data.uncertainty
    if unc is None:
        return
    lo = _format_value(unc.get("ci_90_lower_mg"))
    hi = _format_value(unc.get("ci_90_upper_mg"))
    point = _format_value(unc.get("point_estimate_mg"))
    conf = unc.get("confidence", "?")
    print()
    print(
        f"Uncertainty: {point} mg [{lo} – {hi}] 90% CI  ·  confidence: {conf}"
    )
    sens = unc.get("sensitivity") or {}
    if sens:
        top = max(sens, key=sens.get)
        print(f"Top sensitivity: {top} ({sens[top] * 100:.1f}%)")
    rec_txt = unc.get("recommendation") or ""
    if rec_txt:
        print(rec_txt)


def _cmd_recommend(args: argparse.Namespace) -> int:
    dp = _build_dose_projection(args)
    if dp is None:
        return _fail(
            "recommend requires at least one target: "
            "--noael + --noael-species, --target-kd, or --target-ceff"
        )
    unc_cfg: UncertaintyConfig | None = None
    if args.uncertainty:
        unc_cfg = UncertaintyConfig(n_samples=args.n_samples)
    try:
        pipe = _build_pipeline_from_smiles(
            args, dose_projection=dp, uncertainty=unc_cfg
        )
        result = pipe.run()
    except ValueError as e:
        return _fail(f"Invalid input: {e}")
    except Exception as e:  # pragma: no cover
        return _fail(f"{type(e).__name__}: {e}", code=2)

    data = collect(result)
    if args.json:
        print(json.dumps(asdict(data), indent=2, default=str))
        return 0
    _print_dose_recommendation(data)
    _print_uncertainty_summary(data)
    return 0
```

Then in `_build_parser()`, add these options and subcommands before `return p`:

```python
    # Shared dose-projection options
    def _add_dose_opts(sp: argparse.ArgumentParser) -> None:
        sp.add_argument("--noael", type=float, default=None)
        sp.add_argument("--noael-species", default=None, dest="noael_species")
        sp.add_argument("--target-kd", type=float, default=None, dest="target_kd")
        sp.add_argument(
            "--target-ceff", type=float, default=None, dest="target_ceff"
        )
        sp.add_argument(
            "--safety-factor", type=float, default=10.0, dest="safety_factor"
        )
        sp.add_argument("--tau", type=float, default=24.0)
        sp.add_argument(
            "--body-weight", type=float, default=70.0, dest="body_weight"
        )

    # translate
    sp_tr = sub.add_parser(
        "translate", help="Pipeline + FIH dose projection (Layers 0-3)"
    )
    sp_tr.add_argument("smiles")
    sp_tr.add_argument(
        "--route", required=True, choices=["iv_bolus", "iv_infusion", "oral"]
    )
    sp_tr.add_argument("--dose", type=float, required=True, help="dose in mg")
    sp_tr.add_argument("--duration", type=float, default=72.0)
    sp_tr.add_argument(
        "--infusion-duration", type=float, default=0.0, dest="infusion_duration"
    )
    _add_dose_opts(sp_tr)
    _add_common(sp_tr)
    sp_tr.set_defaults(func=_cmd_translate)

    # recommend
    sp_rc = sub.add_parser(
        "recommend",
        help="Pipeline + dose projection + optional uncertainty (Layers 0-4)",
    )
    sp_rc.add_argument("smiles")
    sp_rc.add_argument(
        "--route", required=True, choices=["iv_bolus", "iv_infusion", "oral"]
    )
    sp_rc.add_argument("--dose", type=float, required=True, help="dose in mg")
    sp_rc.add_argument("--duration", type=float, default=72.0)
    sp_rc.add_argument(
        "--infusion-duration", type=float, default=0.0, dest="infusion_duration"
    )
    _add_dose_opts(sp_rc)
    sp_rc.add_argument("--uncertainty", action="store_true", default=False)
    sp_rc.add_argument("--n-samples", type=int, default=500, dest="n_samples")
    sp_rc.add_argument("--seed", type=int, default=42)
    _add_common(sp_rc)
    sp_rc.set_defaults(func=_cmd_recommend)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/unit/test_cli.py -v`
Expected: 10 passed.

- [ ] **Step 5: Commit**

```bash
git add src/charon/cli/main.py tests/unit/test_cli.py
git commit -m "feat(cli): translate and recommend subcommands with dose projection"
```

---

## Task 11: CLI — `report` subcommand

**Files:**
- Modify: `src/charon/cli/main.py`
- Test: `tests/unit/test_cli.py`

- [ ] **Step 1: Write the failing test**

Append to `tests/unit/test_cli.py`:

```python
def test_report_writes_md_and_json(tmp_path, monkeypatch):
    _install_fake_pipeline_with_dose(monkeypatch)
    out_base = tmp_path / "run"
    rc = main(
        [
            "report",
            "CCO",
            "--route", "oral",
            "--dose", "5",
            "--noael", "50",
            "--noael-species", "rat",
            "--output", str(out_base) + ".md",
        ]
    )
    assert rc == 0
    md = tmp_path / "run.md"
    js = tmp_path / "run.json"
    assert md.exists()
    assert js.exists()
    md_text = md.read_text()
    assert "# FIH Dose Rationale Report" in md_text
    assert "CCO" in md_text
    parsed = json.loads(js.read_text())
    assert parsed["route"] == "oral"
    assert parsed["dose_mg"] == pytest.approx(5.0)


def test_report_accepts_path_without_suffix(tmp_path, monkeypatch):
    _install_fake_pipeline_with_dose(monkeypatch)
    out_base = tmp_path / "no_suffix"
    rc = main(
        [
            "report",
            "CCO",
            "--route", "oral",
            "--dose", "5",
            "--noael", "50",
            "--noael-species", "rat",
            "--output", str(out_base),
        ]
    )
    assert rc == 0
    assert (tmp_path / "no_suffix.md").exists()
    assert (tmp_path / "no_suffix.json").exists()


def test_report_requires_output(capsys, monkeypatch):
    _install_fake_pipeline_with_dose(monkeypatch)
    with pytest.raises(SystemExit):
        # argparse exits for missing required --output
        main(
            [
                "report", "CCO",
                "--route", "oral", "--dose", "5",
                "--noael", "50", "--noael-species", "rat",
            ]
        )
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/unit/test_cli.py -v -k "report"`
Expected: FAIL — `report` subcommand not registered.

- [ ] **Step 3: Write minimal implementation**

Edit `src/charon/cli/main.py`. Add after `_cmd_recommend`:

```python
# ---------------------------------------------------------------------------
# Subcommand: report
# ---------------------------------------------------------------------------


def _cmd_report(args: argparse.Namespace) -> int:
    from pathlib import Path as _P

    dp = _build_dose_projection(args)
    if dp is None:
        return _fail(
            "report requires at least one target: "
            "--noael + --noael-species, --target-kd, or --target-ceff"
        )
    unc_cfg: UncertaintyConfig | None = None
    if args.uncertainty:
        unc_cfg = UncertaintyConfig(n_samples=args.n_samples)

    if not args.quiet:
        print(f"[1/4] Running pipeline for {args.smiles}...")
    try:
        pipe = _build_pipeline_from_smiles(
            args, dose_projection=dp, uncertainty=unc_cfg
        )
        result = pipe.run()
    except ValueError as e:
        return _fail(f"Invalid input: {e}")
    except Exception as e:  # pragma: no cover
        return _fail(f"{type(e).__name__}: {e}", code=2)

    if not args.quiet:
        if args.uncertainty:
            print(f"[2/4] Uncertainty propagation ({args.n_samples} samples)...")
        print("[3/4] Collecting report data...")
    data = collect(result)

    full_profile: dict | None = None
    if args.include_full_profile:
        import numpy as _np
        full_profile = {
            "time_h": _np.asarray(result.time_h, dtype=float).tolist(),
            "cp_plasma_ug_L": _np.asarray(result.cp_plasma, dtype=float).tolist(),
            "cp_blood_ug_L": _np.asarray(result.cp_blood, dtype=float).tolist(),
        }

    if not args.quiet:
        print("[4/4] Writing report...")
    md_path, json_path = export_report(data, _P(args.output), full_profile=full_profile)
    print()
    print("Wrote:")
    print(f"  {md_path}")
    print(f"  {json_path}")
    return 0
```

Then in `_build_parser()`, add the `report` subcommand:

```python
    # report
    sp_rp = sub.add_parser(
        "report",
        help="Full pipeline + write Markdown/JSON report files",
    )
    sp_rp.add_argument("smiles")
    sp_rp.add_argument(
        "--route", required=True, choices=["iv_bolus", "iv_infusion", "oral"]
    )
    sp_rp.add_argument("--dose", type=float, required=True, help="dose in mg")
    sp_rp.add_argument("--duration", type=float, default=72.0)
    sp_rp.add_argument(
        "--infusion-duration", type=float, default=0.0, dest="infusion_duration"
    )
    _add_dose_opts(sp_rp)
    sp_rp.add_argument("--uncertainty", action="store_true", default=False)
    sp_rp.add_argument("--n-samples", type=int, default=500, dest="n_samples")
    sp_rp.add_argument("--seed", type=int, default=42)
    sp_rp.add_argument(
        "--output", required=True, help="output path (.md appended if absent)"
    )
    sp_rp.add_argument(
        "--include-full-profile",
        action="store_true",
        default=False,
        dest="include_full_profile",
    )
    _add_common(sp_rp)
    sp_rp.set_defaults(func=_cmd_report)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/unit/test_cli.py -v`
Expected: 13 passed.

- [ ] **Step 5: Commit**

```bash
git add src/charon/cli/main.py tests/unit/test_cli.py
git commit -m "feat(cli): report subcommand writes Markdown + JSON siblings"
```

---

## Task 12: End-to-end integration test + full-suite regression gate

**Files:**
- Create: `tests/integration/test_report_e2e.py`

Goal: Run the real pipeline through `report` on a small smoke input and verify both files are written and parseable. No mocks. Then re-run the whole test suite to confirm no regression.

- [ ] **Step 1: Write the failing test**

Create `tests/integration/test_report_e2e.py`:

```python
"""End-to-end integration tests for Sprint 6 report + CLI.

These tests run the real Charon pipeline on a small well-behaved
compound and verify that the generated report contains the expected
structure.  Kept fast by using a tiny uncertainty sample count.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from charon.cli.main import main


# ---------------------------------------------------------------------------
# Fixture: use ethanol (CCO) — fastest known valid SMILES in the suite
# ---------------------------------------------------------------------------

_SMILES = "CCO"


def test_cli_simulate_ethanol_iv_bolus(capsys):
    rc = main(
        [
            "simulate",
            _SMILES,
            "--route", "iv_bolus",
            "--dose", "10",
            "--duration", "24",
        ]
    )
    out = capsys.readouterr().out
    assert rc == 0
    assert "Cmax" in out
    assert "CL apparent" in out


def test_cli_report_e2e_writes_md_and_json(tmp_path: Path):
    out_base = tmp_path / "ethanol_report"
    rc = main(
        [
            "report",
            _SMILES,
            "--route", "iv_bolus",
            "--dose", "10",
            "--duration", "24",
            "--noael", "500",
            "--noael-species", "rat",
            "--output", str(out_base) + ".md",
            "--quiet",
        ]
    )
    assert rc == 0

    md_path = tmp_path / "ethanol_report.md"
    json_path = tmp_path / "ethanol_report.json"
    assert md_path.exists()
    assert json_path.exists()

    md_text = md_path.read_text()
    for section in (
        "# FIH Dose Rationale Report",
        "## 1. Executive Summary",
        "## 2. Compound Profile",
        "## 3. ADME Predictions",
        "## 4. IVIVE",
        "## 5. PK Simulation Results",
        "## 6. FIH Dose Projection",
        "## 8. Limitations",
        "## 9. Appendix",
    ):
        assert section in md_text, f"missing section: {section}"

    # Uncertainty was not requested, so section 7 should not appear
    assert "## 7. Uncertainty Analysis" not in md_text

    # JSON sanity
    payload = json.loads(json_path.read_text())
    assert payload["smiles"] == _SMILES
    assert payload["route"] == "iv_bolus"
    assert payload["dose_mg"] == pytest.approx(10.0)
    assert payload["dose_recommendation"] is not None
    assert payload["dose_recommendation"]["mrsd_mg"] > 0


def test_cli_report_e2e_with_uncertainty(tmp_path: Path):
    out_base = tmp_path / "ethanol_unc"
    rc = main(
        [
            "report",
            _SMILES,
            "--route", "iv_bolus",
            "--dose", "10",
            "--duration", "24",
            "--noael", "500",
            "--noael-species", "rat",
            "--uncertainty",
            "--n-samples", "30",  # small for test speed
            "--output", str(out_base),
            "--quiet",
        ]
    )
    assert rc == 0

    md_path = tmp_path / "ethanol_unc.md"
    json_path = tmp_path / "ethanol_unc.json"
    assert md_path.exists()
    assert json_path.exists()

    md_text = md_path.read_text()
    assert "## 7. Uncertainty Analysis" in md_text

    payload = json.loads(json_path.read_text())
    assert payload["uncertainty"] is not None
    assert payload["uncertainty"]["ci_90_lower_mg"] > 0
    assert payload["uncertainty"]["ci_90_upper_mg"] > payload["uncertainty"]["ci_90_lower_mg"]
```

- [ ] **Step 2: Run the new integration test**

Run: `pytest tests/integration/test_report_e2e.py -v`
Expected: 3 passed. If the uncertainty test is slow (> 60 s), reduce `--n-samples` to 20. If `simulate` fails due to a missing ADMET model file, skip the test with `pytest.importorskip` on `charon.predict.admet_ensemble` — but in this repo the models are checked in.

- [ ] **Step 3: Run full regression**

Run: `pytest tests/ --tb=no -q 2>&1 | tail -5`
Expected: `713 passed` (678 baseline + ~35 new tests). Adjust +/- a few for incidental tests.

- [ ] **Step 4: If any pre-existing test regressed, diagnose without destructive actions**

If any test that previously passed is now failing, read the failing test + the symbol it references and trace the regression to a change in this sprint. Fix forward (new commit). Do NOT `git reset` or rewrite history.

- [ ] **Step 5: Commit**

```bash
git add tests/integration/test_report_e2e.py
git commit -m "test(report): end-to-end CLI report integration"
```

---

## Definition of Done

All checkboxes above completed, plus:

- [ ] `pytest tests/` shows all tests green (no regressions).
- [ ] `charon --help` lists exactly the five subcommands.
- [ ] `charon report CCO --route iv_bolus --dose 10 --noael 500 --noael-species rat --output /tmp/demo.md` produces both `/tmp/demo.md` and `/tmp/demo.json`.
- [ ] The Markdown report contains sections 1, 2, 3, 4, 5, 6, 8, 9 (and 7 when uncertainty is enabled).
- [ ] Empty scaffolds `src/charon/report/figures.py` and `src/charon/report/templates/*.yaml` remain untouched.
- [ ] Spec section 9 acceptance criteria all verified via the test suite.
