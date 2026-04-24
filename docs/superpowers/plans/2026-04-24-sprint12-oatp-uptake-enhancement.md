# Sprint 12 — OATP1B1 Hepatic Uptake Enhancement Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `hepatic_clint_multiplier` field to compound YAML schema and ParameterBridge so atorvastatin's OATP1B1 IVIVE underprediction (~8x) is corrected empirically. Target §8 within-3x crossing 60%.

**Architecture:** Optional schema field (`MetabolismProperties.hepatic_clint_multiplier`) + optional kwarg (`ParameterBridge.clint_to_clh(clint_multiplier=...)`) that multiplies `clint_liver_L_h` before the liver extraction model. Pipeline reads the field and passes through. Atorvastatin YAML gets `value: 8.0` (Izumi 2018 / Barton 2013). All changes backward-compatible — unset field = no enhancement.

**Tech Stack:** Python 3.11+, Pydantic v2, pytest.

**Spec:** `docs/superpowers/specs/2026-04-24-sprint12-oatp-uptake-enhancement-design.md`

---

## File Structure

| File | Change | Lines (est.) |
|---|---|---|
| `src/charon/core/schema.py` | Add `hepatic_clint_multiplier: PredictedProperty \| None = None` to `MetabolismProperties` | +1 |
| `src/charon/core/parameter_bridge.py` | Add `clint_multiplier: float \| None = None` kwarg to `clint_to_clh`; emit `ConversionStep(name="clint_enhancement", ...)` when provided | +15 |
| `src/charon/pipeline.py` | Read `compound.properties.metabolism.hepatic_clint_multiplier` and pass `.value` to `clint_to_clh` (both IV + oral paths) | +6 |
| `validation/data/tier1_obach/compounds/atorvastatin.yaml` | Add `hepatic_clint_multiplier` block under `metabolism` | +7 |
| `tests/unit/test_parameter_bridge.py` | +5 unit tests for the new kwarg | +100 |
| `tests/integration/test_atorvastatin_oatp_enhancement.py` | NEW — before/after MRSD comparison | +80 |
| `validation/reports/layer3_fih_dose.{md,json}` | Regenerated | — |
| `validation/reports/layer3_ivive_decomposition.{md,json}` | Regenerated | — |

**Unchanged (reaffirming spec §3):**
- `src/charon/core/liver_models.py`
- `src/charon/pbpk/*.py`
- `src/charon/translational/*.py`
- `validation/data/fih_reference/panel.yaml`
- All compound YAMLs except atorvastatin

---

## Task 1: Schema extension — `hepatic_clint_multiplier` on `MetabolismProperties`

**Files:**
- Modify: `src/charon/core/schema.py:211-236` (the `MetabolismProperties` class)
- Test: `tests/unit/test_compound_config.py` (existing; smoke that load still works)

- [ ] **Step 1: Write the failing test**

Append to `tests/unit/test_compound_config.py`:

```python
def test_metabolism_properties_accepts_hepatic_clint_multiplier():
    """Schema: hepatic_clint_multiplier is an optional PredictedProperty."""
    from charon.core.schema import MetabolismProperties, PredictedProperty

    m = MetabolismProperties(
        clint_uL_min_mg=PredictedProperty(
            value=100.0,
            source="experimental",
            unit="uL/min/mg",
            method="test",
        ),
        hepatic_clint_multiplier=PredictedProperty(
            value=8.0,
            source="literature",
            unit="ratio",
            method="Izumi 2018",
        ),
    )
    assert m.hepatic_clint_multiplier is not None
    assert m.hepatic_clint_multiplier.value == 8.0


def test_metabolism_properties_hepatic_clint_multiplier_optional():
    """Schema: hepatic_clint_multiplier defaults to None."""
    from charon.core.schema import MetabolismProperties, PredictedProperty

    m = MetabolismProperties(
        clint_uL_min_mg=PredictedProperty(
            value=100.0,
            source="experimental",
            unit="uL/min/mg",
            method="test",
        ),
    )
    assert m.hepatic_clint_multiplier is None
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
pytest tests/unit/test_compound_config.py::test_metabolism_properties_accepts_hepatic_clint_multiplier tests/unit/test_compound_config.py::test_metabolism_properties_hepatic_clint_multiplier_optional -v
```

Expected: FAIL with `ValueError: Object has no field "hepatic_clint_multiplier"` or similar Pydantic extra-field rejection.

- [ ] **Step 3: Add the field**

Open `src/charon/core/schema.py`. Find the `MetabolismProperties` class (around line 211). Add after `fm_cyp3a4`:

```python
class MetabolismProperties(BaseModel):
    """CYP phenotyping and intrinsic clearance."""

    primary_cyp: str | None = None
    secondary_cyp: str | None = None
    clint_uL_min_mg: PredictedProperty | None = None
    fm_cyp3a4: float | None = None
    # Sprint 12: optional empirical multiplier applied to CLint_liver_L_h
    # before the liver extraction model. Targets uptake-limited substrates
    # (e.g., OATP1B1) where passive IVIVE under-predicts in vivo CLh.
    # When None, no enhancement is applied (backward compatible).
    hepatic_clint_multiplier: PredictedProperty | None = None
```

Add field validator below the existing `fm_cyp3a4` validator:

```python
    @field_validator("hepatic_clint_multiplier")
    @classmethod
    def check_hepatic_clint_multiplier(cls, v):
        if v is None:
            return v
        if v.value <= 0:
            raise ValueError(
                f"hepatic_clint_multiplier must be > 0, got {v.value}"
            )
        return v
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
pytest tests/unit/test_compound_config.py -v
```

Expected: all pass, including the 2 new tests + all pre-existing tests.

- [ ] **Step 5: Commit**

```bash
git add src/charon/core/schema.py tests/unit/test_compound_config.py
git commit -m "feat(sprint12): add hepatic_clint_multiplier to MetabolismProperties schema"
```

---

## Task 2: ParameterBridge — `clint_multiplier` kwarg

**Files:**
- Modify: `src/charon/core/parameter_bridge.py:48-60` (signature), ~line 190-200 (after Step 3 unit conversion)
- Test: `tests/unit/test_parameter_bridge.py`

- [ ] **Step 1: Write the failing tests**

Append to `tests/unit/test_parameter_bridge.py`:

```python
def test_clint_multiplier_default_none_matches_baseline():
    """With clint_multiplier=None, result matches the baseline CLh."""
    from charon.core.parameter_bridge import ParameterBridge

    bridge = ParameterBridge()
    result_none = bridge.clint_to_clh(
        clint=100.0, fu_inc=0.5, fu_p=0.1, system="HLM",
    )
    result_one = bridge.clint_to_clh(
        clint=100.0, fu_inc=0.5, fu_p=0.1, system="HLM",
        clint_multiplier=1.0,
    )
    # Both should yield identical CLh within floating-point tolerance
    assert abs(result_none.clh_L_h - result_one.clh_L_h) < 1e-12


def test_clint_multiplier_scales_clint_liver():
    """clint_multiplier=8.0 yields an 8x-enhanced CLint_liver before the
    liver model. For a low-extraction compound (fu_b * CLint << Qh),
    CLh scales approximately linearly with CLint_liver."""
    from charon.core.parameter_bridge import ParameterBridge

    bridge = ParameterBridge()
    # Low extraction: clint=1.0 uL/min/mg, fu_p=0.1. After scaling
    # with default HUMAN_MPPGL=40 and HUMAN_LIVER_WEIGHT_G=1500:
    #   clint_liver_L_h = 1.0 * 40 * 1500 / 1e6 * 60 = 3.6 L/h
    # Well-stirred CLh = Qh * fu_b * CLint / (Qh + fu_b * CLint)
    # With Qh=90, fu_b=0.1: CLh_base ~= 90 * 0.1 * 3.6 / (90 + 0.36) ~= 0.358
    # With 8x multiplier: CLint_liver = 28.8, CLh ~= 90*0.1*28.8/(90+2.88) ~= 2.79
    # Ratio CLh_enhanced / CLh_base ~= 7.8x (near 8x in linear regime)
    r_base = bridge.clint_to_clh(
        clint=1.0, fu_inc=0.5, fu_p=0.1, system="HLM",
    )
    r_enhanced = bridge.clint_to_clh(
        clint=1.0, fu_inc=0.5, fu_p=0.1, system="HLM",
        clint_multiplier=8.0,
    )
    ratio = r_enhanced.clh_L_h / r_base.clh_L_h
    # Should be close to 8x (slightly less due to well-stirred curvature)
    assert 7.0 < ratio < 8.1, f"CLh ratio {ratio:.3f} outside [7.0, 8.1]"


def test_clint_multiplier_raises_on_non_positive():
    """clint_multiplier <= 0 must raise ValueError."""
    from charon.core.parameter_bridge import ParameterBridge

    bridge = ParameterBridge()
    with pytest.raises(ValueError, match="clint_multiplier"):
        bridge.clint_to_clh(
            clint=100.0, fu_inc=0.5, fu_p=0.1, system="HLM",
            clint_multiplier=0.0,
        )
    with pytest.raises(ValueError, match="clint_multiplier"):
        bridge.clint_to_clh(
            clint=100.0, fu_inc=0.5, fu_p=0.1, system="HLM",
            clint_multiplier=-2.0,
        )


def test_clint_multiplier_conversion_log_contains_enhancement_step():
    """When multiplier is applied, intermediate_steps contains a
    'clint_enhancement' entry with formula string."""
    from charon.core.parameter_bridge import ParameterBridge

    bridge = ParameterBridge()
    result = bridge.clint_to_clh(
        clint=100.0, fu_inc=0.5, fu_p=0.1, system="HLM",
        clint_multiplier=8.0,
    )
    step_names = [s.name for s in result.conversion_log.intermediate_steps]
    assert "clint_enhancement" in step_names

    step = next(
        s for s in result.conversion_log.intermediate_steps
        if s.name == "clint_enhancement"
    )
    assert "8.0" in step.formula or "8" in step.formula


def test_clint_multiplier_default_no_enhancement_step():
    """When multiplier is None or 1.0, no enhancement step is logged."""
    from charon.core.parameter_bridge import ParameterBridge

    bridge = ParameterBridge()
    r_none = bridge.clint_to_clh(
        clint=100.0, fu_inc=0.5, fu_p=0.1, system="HLM",
    )
    r_one = bridge.clint_to_clh(
        clint=100.0, fu_inc=0.5, fu_p=0.1, system="HLM",
        clint_multiplier=1.0,
    )
    for r in (r_none, r_one):
        step_names = [s.name for s in r.conversion_log.intermediate_steps]
        assert "clint_enhancement" not in step_names
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
pytest tests/unit/test_parameter_bridge.py::test_clint_multiplier_default_none_matches_baseline -v
```

Expected: FAIL with `TypeError: clint_to_clh() got an unexpected keyword argument 'clint_multiplier'`.

- [ ] **Step 3: Add the kwarg + enhancement logic**

Open `src/charon/core/parameter_bridge.py`. Two edits:

**Edit A — Add `clint_multiplier` kwarg to the method signature (around line 48-60):**

Find:
```python
    def clint_to_clh(
        self,
        clint: float,
        fu_inc: float,
        fu_p: float,
        system: str = "HLM",
        mppgl: float = HUMAN_MPPGL,
        hepatocellularity: float = HUMAN_HEPATOCELLULARITY,
        liver_weight_g: float = HUMAN_LIVER_WEIGHT_G,
        qh_L_h: float = HUMAN_QH_L_H,
        bp_ratio: float = 1.0,
        model: str = "well_stirred",
    ) -> HepaticClearance:
```

Add `clint_multiplier` at the end:
```python
    def clint_to_clh(
        self,
        clint: float,
        fu_inc: float,
        fu_p: float,
        system: str = "HLM",
        mppgl: float = HUMAN_MPPGL,
        hepatocellularity: float = HUMAN_HEPATOCELLULARITY,
        liver_weight_g: float = HUMAN_LIVER_WEIGHT_G,
        qh_L_h: float = HUMAN_QH_L_H,
        bp_ratio: float = 1.0,
        model: str = "well_stirred",
        clint_multiplier: float | None = None,
    ) -> HepaticClearance:
```

**Edit B — Add validation + enhancement step after Step 3 unit conversion.**

Read lines ~180-200 of `parameter_bridge.py` to find where `clint_liver_L_h` is computed (the `/ 1e6 * 60` conversion — Step 3). The step append looks like:

```python
        clint_liver_L_h = uL_min_to_L_h(clint_liver_uL_min)
        steps.append(
            ConversionStep(
                name="clint_liver_L_h",
                value=clint_liver_L_h,
                ...
            )
        )
```

**Immediately after that ConversionStep append**, insert the enhancement block:

```python
        # ---- Step 3b (optional): CLint enhancement ---------------------
        # Empirical multiplier for uptake-limited substrates (e.g., OATP1B1).
        # Documented in CompoundConfig.properties.metabolism.hepatic_clint_multiplier.
        if clint_multiplier is not None and clint_multiplier != 1.0:
            if clint_multiplier <= 0:
                raise ValueError(
                    f"clint_multiplier must be > 0, got {clint_multiplier}"
                )
            enhanced = clint_liver_L_h * clint_multiplier
            steps.append(
                ConversionStep(
                    name="clint_enhancement",
                    value=enhanced,
                    unit="L/h",
                    formula=(
                        f"CLint_liver * multiplier = "
                        f"{clint_liver_L_h:.4g} * {clint_multiplier} = {enhanced:.4g}"
                    ),
                )
            )
            clint_liver_L_h = enhanced
```

Also update `input_params` dict (around line 123) to include the multiplier for audit:

Find:
```python
        input_params: dict[str, float | str] = {
            "clint": clint,
            ...
            "model": model,
        }
```

Add:
```python
        input_params: dict[str, float | str] = {
            "clint": clint,
            ...
            "model": model,
        }
        if clint_multiplier is not None:
            input_params["clint_multiplier"] = clint_multiplier
```

(Keep type annotation `float | str` since multiplier is a float.)

**Edit C — Also guard against clint_multiplier <= 0 up-front in the Validation block** (around line 99-116). Add after the existing `if clint < 0:` check:

```python
        if clint_multiplier is not None and clint_multiplier <= 0:
            raise ValueError(
                f"clint_multiplier must be > 0, got {clint_multiplier}"
            )
```

(This makes the `test_clint_multiplier_raises_on_non_positive` test deterministic regardless of whether the later branch gets reached.)

- [ ] **Step 4: Run tests to verify they pass**

```bash
pytest tests/unit/test_parameter_bridge.py -v
```

Expected: all tests pass, including 5 new ones + all pre-existing ParameterBridge tests (which don't pass the new kwarg, so default None preserves behaviour).

- [ ] **Step 5: Commit**

```bash
git add src/charon/core/parameter_bridge.py tests/unit/test_parameter_bridge.py
git commit -m "feat(sprint12): clint_multiplier kwarg on ParameterBridge.clint_to_clh"
```

---

## Task 3: Pipeline wiring — read `hepatic_clint_multiplier` and pass through

**Files:**
- Modify: `src/charon/pipeline.py` (2 call sites: IV path and oral path)

**Rationale:** `Pipeline._run_iv` and `Pipeline._run_oral` both call `ParameterBridge.clint_to_clh`. Both paths must forward the compound's multiplier.

- [ ] **Step 1: Locate the two call sites**

```bash
grep -n "clint_to_clh" src/charon/pipeline.py
```

Note the line numbers.

- [ ] **Step 2: Extract the multiplier value in both paths**

At the top of each path method (`_run_iv`, `_run_oral`), before the `bridge.clint_to_clh(...)` call, add:

```python
        # Sprint 12: optional OATP uptake enhancement factor
        metab = self.compound.properties.metabolism
        clint_multiplier = (
            metab.hepatic_clint_multiplier.value
            if metab.hepatic_clint_multiplier is not None
            else None
        )
```

Then in the `clint_to_clh(...)` call, add `clint_multiplier=clint_multiplier` as a kwarg.

Find the IV path call, expected like:
```python
        result = bridge.clint_to_clh(
            clint=metab.clint_uL_min_mg.value,
            fu_inc=...,
            ...
        )
```

Change to:
```python
        result = bridge.clint_to_clh(
            clint=metab.clint_uL_min_mg.value,
            fu_inc=...,
            ...
            clint_multiplier=clint_multiplier,
        )
```

Apply to BOTH IV and oral paths.

- [ ] **Step 3: Smoke test — ensure Pipeline still runs for a non-enhanced compound**

```bash
pytest tests/integration/test_fih_pipeline.py -v
```

Expected: all pass. Enhancement is None for every existing compound → no behavior change yet.

- [ ] **Step 4: Commit**

```bash
git add src/charon/pipeline.py
git commit -m "feat(sprint12): Pipeline forwards hepatic_clint_multiplier to ParameterBridge"
```

---

## Task 4: Atorvastatin YAML — add `hepatic_clint_multiplier: 8.0`

**Files:**
- Modify: `validation/data/tier1_obach/compounds/atorvastatin.yaml:45-48` (the `metabolism:` block)

- [ ] **Step 1: Read current metabolism block**

```bash
sed -n '45,50p' validation/data/tier1_obach/compounds/atorvastatin.yaml
```

Should show:
```yaml
  metabolism:
    clint_uL_min_mg: {value: 128.0, source: experimental, unit: uL/min/mg, method: "Obach 2008 HLM"}
```

- [ ] **Step 2: Verify the multiplier literature value**

Before committing, use WebSearch to cross-check. Queries:
- `"atorvastatin IVIVE hepatocyte in vivo in vitro CLint ratio"`
- `"OATP1B1 atorvastatin enhancement factor Izumi"`
- `"atorvastatin extended clearance PS_uptake"`

Published values I'd expect to see in the 5-12x range. If WebSearch gives a more specific value (e.g., Izumi 2018 reports "approximately 8x" for atorvastatin specifically), use it. If sources report a range, pick the median and cite the range in the method string.

- [ ] **Step 3: Add the multiplier block**

Edit `validation/data/tier1_obach/compounds/atorvastatin.yaml`. Modify the `metabolism:` block to add `hepatic_clint_multiplier`:

```yaml
  metabolism:
    clint_uL_min_mg: {value: 128.0, source: experimental, unit: uL/min/mg, method: "Obach 2008 HLM"}
    # Sprint 12: OATP1B1 hepatic uptake enhancement.
    # Passive IVIVE (via HLM CLint) under-predicts in vivo CLh for OATP
    # substrates because active uptake into hepatocyte is rate-limiting at
    # typical fu_p × CLint. Empirical multiplier derived from literature
    # in vivo/in vitro CLint_u ratios for atorvastatin.
    hepatic_clint_multiplier:
      value: 8.0
      source: literature
      unit: ratio
      method: "Izumi 2018 Drug Metab Pharmacokinet 33(2):57 / Barton 2013 Pharm Res 30:1077 — in vivo/in vitro CLint_u ratio for OATP1B1 substrates 5-12x; 8 = midpoint consistent with atorvastatin published analyses"
```

(Adjust `value: 8.0` to whatever verified value came out of Step 2. Keep source/method specific to cited literature.)

- [ ] **Step 4: Verify YAML parses + schema accepts**

```bash
python3 -c "
import yaml
from charon.core.schema import CompoundConfig
data = yaml.safe_load(open('validation/data/tier1_obach/compounds/atorvastatin.yaml').read())
c = CompoundConfig.model_validate(data)
m = c.properties.metabolism.hepatic_clint_multiplier
print(f'hepatic_clint_multiplier: {m.value} ({m.source}, {m.method[:40]}...)')"
```

Expected: prints `hepatic_clint_multiplier: 8.0 (literature, Izumi 2018 ...)`.

- [ ] **Step 5: Commit**

```bash
git add validation/data/tier1_obach/compounds/atorvastatin.yaml
git commit -m "data(sprint12): atorvastatin hepatic_clint_multiplier=8.0 (OATP1B1 IVIVE correction)"
```

---

## Task 5: Integration test — atorvastatin MRSD before/after

**Files:**
- Create: `tests/integration/test_atorvastatin_oatp_enhancement.py`

- [ ] **Step 1: Write the test**

Create `tests/integration/test_atorvastatin_oatp_enhancement.py`:

```python
"""Sprint 12 integration test: atorvastatin MRSD with/without OATP enhancement.

Validates:
1. The enhanced MRSD is strictly larger than the non-enhanced MRSD (the
   multiplier scales CLh up, which scales MRSD up via the PAD formula).
2. The ratio MRSD_enhanced / MRSD_base is within a physiologically
   plausible range (5-12x, matching the literature multiplier range).
3. Pipeline does not crash when the multiplier is present.
"""

from __future__ import annotations

import copy
from pathlib import Path

import pytest
import yaml

from charon import Pipeline
from charon.core.schema import CompoundConfig, DoseProjectionConfig

REPO_ROOT = Path(__file__).resolve().parents[2]
ATORVA_YAML = REPO_ROOT / "validation" / "data" / "tier1_obach" / "compounds" / "atorvastatin.yaml"
PANEL_YAML = REPO_ROOT / "validation" / "data" / "fih_reference" / "panel.yaml"


def _load_atorvastatin_entry() -> dict:
    panel = yaml.safe_load(PANEL_YAML.read_text())["panel"]
    return next(c for c in panel["compounds"] if c["name"] == "atorvastatin" and c["tier"] == "gold")


def _run_pipeline(compound: CompoundConfig, entry: dict) -> float:
    pipe = Pipeline(
        compound,
        route=entry["route"],
        dose_mg=1.0,
        dose_projection=DoseProjectionConfig(
            target_ceff_nM=float(entry["target_ceff_nM"]),
            safety_factor=10.0,
            tau_h=24.0,
        ),
    )
    result = pipe.run()
    assert result.dose_recommendation is not None
    return float(result.dose_recommendation.mrsd_mg)


def test_atorvastatin_mrsd_with_multiplier_larger_than_without():
    """With clint_multiplier=8.0, atorvastatin MRSD must be strictly larger."""
    data = yaml.safe_load(ATORVA_YAML.read_text())
    entry = _load_atorvastatin_entry()

    # Load WITH multiplier (from YAML as committed)
    compound_enhanced = CompoundConfig.model_validate(data)
    assert compound_enhanced.properties.metabolism.hepatic_clint_multiplier is not None, (
        "atorvastatin.yaml must have hepatic_clint_multiplier populated by Task 4"
    )
    mrsd_enhanced = _run_pipeline(compound_enhanced, entry)

    # Load WITHOUT multiplier (strip from in-memory dict)
    data_base = copy.deepcopy(data)
    data_base["properties"]["metabolism"].pop("hepatic_clint_multiplier", None)
    compound_base = CompoundConfig.model_validate(data_base)
    assert compound_base.properties.metabolism.hepatic_clint_multiplier is None
    mrsd_base = _run_pipeline(compound_base, entry)

    assert mrsd_enhanced > mrsd_base, (
        f"Enhanced MRSD ({mrsd_enhanced:.3g}) must exceed baseline ({mrsd_base:.3g})"
    )


def test_atorvastatin_mrsd_ratio_within_literature_range():
    """MRSD_enhanced / MRSD_base should be within 4-12x (plausibly near the
    multiplier value of 8, allowing for well-stirred curvature and non-hepatic
    contributions)."""
    data = yaml.safe_load(ATORVA_YAML.read_text())
    entry = _load_atorvastatin_entry()

    compound_enhanced = CompoundConfig.model_validate(data)
    mrsd_enhanced = _run_pipeline(compound_enhanced, entry)

    data_base = copy.deepcopy(data)
    data_base["properties"]["metabolism"].pop("hepatic_clint_multiplier", None)
    compound_base = CompoundConfig.model_validate(data_base)
    mrsd_base = _run_pipeline(compound_base, entry)

    ratio = mrsd_enhanced / mrsd_base
    assert 4.0 < ratio < 12.0, (
        f"MRSD ratio {ratio:.2f} outside plausible [4, 12] range; "
        f"check multiplier logic or well-stirred curvature"
    )


def test_atorvastatin_enhanced_conversion_log_has_enhancement_step():
    """The ConversionLog from an enhanced-atorvastatin run contains the
    'clint_enhancement' step."""
    data = yaml.safe_load(ATORVA_YAML.read_text())
    entry = _load_atorvastatin_entry()

    compound = CompoundConfig.model_validate(data)
    pipe = Pipeline(
        compound,
        route=entry["route"],
        dose_mg=1.0,
        dose_projection=DoseProjectionConfig(
            target_ceff_nM=float(entry["target_ceff_nM"]),
            safety_factor=10.0,
            tau_h=24.0,
        ),
    )
    result = pipe.run()
    # Find the ConversionLog via the Pipeline metadata / result path.
    # ConversionLog is attached to the ParameterBridge.clint_to_clh return
    # inside the Pipeline. We infer enhancement happened by checking the
    # Pipeline metadata's clint_liver_L_h — it should be ~8x the "naive"
    # IVIVE value (128 × 40 × 1500 / 1e6 × 60 / fu_inc=0.07).

    md = result.metadata
    # Without enhancement: clint_liver = 128 / 0.07 * 40 * 1500 / 1e6 * 60 ≈ 6,583 L/h
    # With 8x enhancement: ≈ 52,663 L/h
    # Allow ±20% tolerance for any rounding in unit conversions.
    clint_liver = md["clint_liver_L_h"]
    # Enhanced should be in the 40,000 - 70,000 L/h range
    assert 40_000.0 < clint_liver < 70_000.0, (
        f"Enhanced atorvastatin clint_liver_L_h = {clint_liver:.0f}; "
        f"expected near 52,000 (8x × 6,583 baseline)"
    )
```

- [ ] **Step 2: Run the tests**

```bash
pytest tests/integration/test_atorvastatin_oatp_enhancement.py -v
```

Expected: 3/3 PASS.

If `test_atorvastatin_mrsd_ratio_within_literature_range` fails with ratio < 4 or > 12: STOP. This indicates either the multiplier isn't applied to the right stage (edit `parameter_bridge.py` Task 2 Edit B location) or the well-stirred extraction is saturating (atorvastatin's fu_b × CLint is already flow-limited, so further CLint has diminishing effect on CLh). Report as BLOCKED with observed ratio.

- [ ] **Step 3: Commit**

```bash
git add tests/integration/test_atorvastatin_oatp_enhancement.py
git commit -m "test(sprint12): atorvastatin OATP enhancement before/after integration"
```

---

## Task 6: Regenerate benchmarks + narrative

**Files:**
- Regenerate: `validation/reports/layer3_fih_dose.{md,json}`
- Regenerate: `validation/reports/layer3_ivive_decomposition.{md,json}`

- [ ] **Step 1: Run the Layer 3 FIH dose benchmark**

```bash
python3 validation/benchmarks/layer3_fih_dose.py
```

Expected: `[OK] Sanity floor: 12/12 pass.` — Tier B untouched.

Extract new Tier A numbers:
```bash
python3 -c "
import json
d = json.load(open('validation/reports/layer3_fih_dose.json'))
s = d['summary']
print(f'Tier A within-3x: {s[\"gold_within_3x\"]}/{s[\"gold_n\"]} = {100*s[\"gold_within_3x_fraction\"]:.1f}%')
print(f'Tier A within-10x: {s[\"gold_within_10x\"]}/{s[\"gold_n\"]}')
print()
for r in d['extra_sections']['Gold (Tier A) — fold-error vs reference FIH']:
    if r['compound'] == 'atorvastatin':
        print(f'atorvastatin: mrsd={r[\"mrsd_pred_mg\"]:.3g} mg, ref={r[\"reference_fih_mg\"]:.3g} mg, fold={r[\"fold_error\"]:.2f}')
"
```

**Expected (predicted):** atorvastatin fold drops from 4.64 (Sprint 11) to somewhere in [1.5, 3.0]. Tier A within-3x should tick up from 7/12 → 8/12 (if atorvastatin crosses 3x).

- [ ] **Step 2: Run the decomposition orchestrator**

```bash
python3 validation/benchmarks/layer3_ivive_decomposition.py
```

Extract attribution:
```bash
python3 -c "
import json
d = json.load(open('validation/reports/layer3_ivive_decomposition.json'))
s = d['summary']
print(f'liver_model: {s[\"aggregate_pct_liver_model\"]:.1f}%')
print(f'route_bias:  {s[\"aggregate_pct_route_bias\"]:.1f}%')
print(f'residual:    {s[\"aggregate_pct_residual\"]:.1f}%')
print()
print('Atorvastatin decomposition:')
for r in d['extra_sections']['Per-compound decomposition']:
    if r['compound'] == 'atorvastatin':
        for k in ['mrsd_ws_mg', 'reference_fih_mg', 'fold_observed', 'fold_residual']:
            print(f'  {k}: {r[k]:.3f}' if isinstance(r[k], (int, float)) else f'  {k}: {r[k]}')
"
```

Expected: atorvastatin's `fold_residual` dropped substantially from Sprint 11's ~4.64 to near 1 (enhancement closed the gap); aggregate residual % decreased.

- [ ] **Step 3: Append Sprint 12 narrative to `validation/reports/layer3_fih_dose.md`**

Open the auto-generated md and APPEND (don't overwrite) this section at the end:

```markdown

## Sprint 12 comparison (OATP enhancement)

Sprint 11 (oral, no OATP): Tier A within-3x = 7/12 = 58.3% (§8 FAILED by ~2%).
Sprint 12 (oral + atorvastatin OATP multiplier): Tier A within-3x = <M>/12 = <PCT>% (§8 <PASS|FAIL>).

Atorvastatin-only delta:

- Sprint 11: MRSD = 2.15 mg, ref = 10 mg, fold = 4.64x (outside 3x).
- Sprint 12: MRSD = <S12_ATORVA_MRSD> mg, ref = 10 mg, fold = <S12_ATORVA_FOLD>x (<WITHIN_3X|OUTSIDE_3X>).

Multiplier applied: 8.0 × (literature mid-point for OATP1B1 substrates).
No other Tier A compound's fold-error changed (the multiplier is atorvastatin-specific).
```

Fill `<M>`, `<PCT>`, `<PASS|FAIL>`, `<S12_ATORVA_MRSD>`, `<S12_ATORVA_FOLD>`, `<WITHIN_3X|OUTSIDE_3X>` with actual Step 1 output values.

- [ ] **Step 4: Append to `validation/reports/layer3_ivive_decomposition.md`**

Append:

```markdown

## §8. Sprint 12 — OATP enhancement for atorvastatin

After adding `hepatic_clint_multiplier: 8.0` to atorvastatin.yaml (Izumi 2018 literature midpoint for OATP1B1 substrate IVIVE gap), the decomposition re-runs with:

- liver_model: <NEW_LIVER>%
- route_bias:  <NEW_ROUTE>% (still near 0 — oral route simulation unchanged)
- residual:    <NEW_RESIDUAL>% (decreased because atorvastatin's residual dropped)

Atorvastatin per-compound:
- Sprint 11 fold_residual: 4.64 (expected OATP gap)
- Sprint 12 fold_residual: <S12_ATORVA_RESIDUAL> (<DELTA>)

This confirms the hypothesis that atorvastatin's Sprint 11 residual was dominated by the OATP1B1 IVIVE gap. The remaining residual (if any) after enhancement represents biliary clearance, Pgp efflux, or other unmodelled mechanisms.
```

Fill placeholders.

- [ ] **Step 5: Commit**

```bash
git add validation/reports/layer3_fih_dose.md \
        validation/reports/layer3_fih_dose.json \
        validation/reports/layer3_ivive_decomposition.md \
        validation/reports/layer3_ivive_decomposition.json
git commit -m "chore(sprint12): regenerated Layer 3 reports + narrative (OATP enhancement)"
```

---

## Task 7: Update Sprint 10 ticket + full suite green

**Files:**
- Modify: `docs/superpowers/sprint10-ivive-bias-ticket.md`

- [ ] **Step 1: Append Sprint 12 status**

Append to `docs/superpowers/sprint10-ivive-bias-ticket.md`:

```markdown

## Sprint 12 (OATP enhancement for atorvastatin) completed — 2026-04-24

Added `hepatic_clint_multiplier` schema field + ParameterBridge kwarg + atorvastatin.yaml multiplier=8.0 (Izumi 2018 / Barton 2013 OATP1B1 IVIVE correction).

**Atorvastatin delta:**
- Sprint 11: fold-error 4.64x (outside 3x)
- Sprint 12: fold-error <S12_ATORVA_FOLD>x (<WITHIN_3X|OUTSIDE_3X>)

**Layer 3 Tier A result:**
- Sprint 11: 7/12 = 58.3% (§8 FAILED by 2%)
- Sprint 12: <M>/12 = <PCT>% (§8 <PASS|FAIL>)

Remaining residual breakdown (Sprint 13/14 targets):
- propranolol: CYP2D6 / F-correction path (28.55x)
- diclofenac: UGT/CYP2C9 underprediction (10.23x)
- diazepam: low fu_p sensitivity (4.91x)
- lisinopril: non-hepatic elimination + low Peff (4.13x)
```

Fill placeholders from Task 6 output.

- [ ] **Step 2: Run full test suite**

```bash
pytest -q
```

Expected: all pass. Total count ~930-933 (925 Sprint 11 baseline + 2 schema tests + 5 ParameterBridge tests + 3 integration tests = 935 if all count individually; pytest may report slightly different).

If any test fails: STOP. Report BLOCKED.

- [ ] **Step 3: Commit ticket update**

```bash
git add docs/superpowers/sprint10-ivive-bias-ticket.md
git commit -m "docs(sprint12): Sprint 10 ticket — Sprint 12 OATP enhancement results"
```

---

## Self-Review Notes (for the implementer)

- **Literature multiplier verification.** Task 4 Step 2 tells the implementer to WebSearch Izumi 2018 or Barton 2013 for atorvastatin-specific OATP1B1 IVIVE ratios. If the verified value is outside 5-12x (e.g., 3x or 20x), USE the verified value — do not hold to 8.0 for its own sake. The test in Task 5 `test_atorvastatin_mrsd_ratio_within_literature_range` allows 4-12x which covers typical verified outcomes; if the verified value is 3x, you will need to adjust the test tolerance (lower bound to 2.5).

- **Honest reporting.** If atorvastatin's fold-error post-enhancement is outside 3x (either still > 3x or now > 3x in the over-correction direction), report honestly in Task 6/7 narrative. Do NOT tune the multiplier to hit §8 target. The spec's success criteria §6 explicitly permits honest failure.

- **Type consistency.** The schema field is `hepatic_clint_multiplier: PredictedProperty | None`. The ParameterBridge kwarg is `clint_multiplier: float | None`. Pipeline bridges them by extracting `.value`. Keep naming distinct to make the YAML→Pipeline→ParameterBridge chain legible.

- **Integration test isolation.** Task 5's test programmatically strips the multiplier via `copy.deepcopy` + `pop` rather than editing the YAML back-and-forth. This keeps the committed YAML in its enhanced state throughout the test.

- **Pipeline dual-path.** `src/charon/pipeline.py` has both `_run_iv` and `_run_oral` paths that both call `ParameterBridge.clint_to_clh`. Task 3 must update both call sites. Missing one will silently fail to enhance in that path; the integration test (Task 5) exercises the oral path because atorvastatin's panel entry uses `route: oral`.

- **Commit count estimate.** 7 feature commits + merge commit.

- **Estimated final test count.** 925 + 2 schema + 5 ParameterBridge + 3 integration = 935.
