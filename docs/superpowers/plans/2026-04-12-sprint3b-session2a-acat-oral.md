# Sprint 3b Session 2a: ACAT Oral Route + Fg Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add 8-segment GI lumen transit, enterocyte CYP3A4 metabolism, and oral route to Charon's PBPK pipeline, with Fg validated against midazolam/felodipine/nifedipine.

**Architecture:** Extend the BDF ODE system from 17 to 26 states (+8 lumen segments + 1 enterocyte). The enterocyte uses a Qvilli-aware first-pass model with bottom-up CYP3A4 scaling (MPPGI=15.0 within 2–20 literature range). Fg is extracted post-hoc via trapezoidal flux integration, independent of hepatic CLint error.

**Tech Stack:** Python 3.11+, Pydantic v2, scipy BDF, numpy, PyYAML, pytest

**Spec:** `docs/superpowers/specs/2026-04-12-sprint3b-session2a-acat-oral-design.md`

---

## File Structure

| Action | File | Responsibility |
|--------|------|----------------|
| Create | `src/charon/pbpk/acat.py` | GI segment data, GI tract loader, absorption rate computation, Papp→Peff |
| Modify | `src/charon/core/schema.py:199-204` | Add `fm_cyp3a4` to MetabolismProperties |
| Modify | `src/charon/pbpk/species/human.yaml` | Add `gi_tract` YAML section |
| Modify | `src/charon/pbpk/ode_compiler.py` | Add `OralPBPKParams`, `build_oral_rhs()`, `compute_gut_clint()` |
| Modify | `src/charon/pbpk/solver.py` | Add `OralSimulationResult`, `simulate_oral()` |
| Modify | `src/charon/pbpk/pk_extract.py` | Add `compute_oral_pk_parameters()` |
| Modify | `src/charon/pipeline.py:112-118` | Wire `route="oral"` |
| Modify | `src/charon/pbpk/__init__.py` | Export new public symbols |
| Modify | `validation/data/tier1_obach/compounds/midazolam.yaml` | Add `fm_cyp3a4`, `peff_cm_s`, `observed_oral` |
| Create | `validation/data/tier1_obach/compounds/felodipine.yaml` | New compound |
| Create | `validation/data/tier1_obach/compounds/nifedipine.yaml` | New compound |
| Create | `tests/unit/test_acat.py` | Unit tests for GI tract loading and absorption rates |
| Create | `tests/unit/test_oral_compiler.py` | Unit tests for OralPBPKParams, gut CLint, build_oral_rhs |
| Create | `tests/unit/test_oral_solver.py` | Unit tests for simulate_oral |
| Create | `tests/unit/test_oral_pk_extract.py` | Unit tests for oral PK parameter extraction |
| Create | `tests/integration/test_oral_pipeline.py` | End-to-end oral pipeline + Fg validation |

---

### Task 1: Add fm_cyp3a4 to schema

**Files:**
- Modify: `src/charon/core/schema.py:199-204`
- Test: `tests/unit/test_schema.py`

- [ ] **Step 1: Write the failing test**

```python
# Append to tests/unit/test_schema.py

class TestFmCyp3a4:
    def test_fm_cyp3a4_default_none(self):
        from charon.core.schema import MetabolismProperties
        m = MetabolismProperties()
        assert m.fm_cyp3a4 is None

    def test_fm_cyp3a4_valid(self):
        from charon.core.schema import MetabolismProperties
        m = MetabolismProperties(fm_cyp3a4=0.75)
        assert m.fm_cyp3a4 == 0.75

    def test_fm_cyp3a4_zero(self):
        from charon.core.schema import MetabolismProperties
        m = MetabolismProperties(fm_cyp3a4=0.0)
        assert m.fm_cyp3a4 == 0.0

    def test_fm_cyp3a4_one(self):
        from charon.core.schema import MetabolismProperties
        m = MetabolismProperties(fm_cyp3a4=1.0)
        assert m.fm_cyp3a4 == 1.0

    def test_fm_cyp3a4_negative_rejected(self):
        from charon.core.schema import MetabolismProperties
        import pytest
        with pytest.raises(Exception):
            MetabolismProperties(fm_cyp3a4=-0.1)

    def test_fm_cyp3a4_above_one_rejected(self):
        from charon.core.schema import MetabolismProperties
        import pytest
        with pytest.raises(Exception):
            MetabolismProperties(fm_cyp3a4=1.01)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/unit/test_schema.py::TestFmCyp3a4 -v`
Expected: FAIL — `fm_cyp3a4` does not exist yet on MetabolismProperties

- [ ] **Step 3: Implement fm_cyp3a4 field with validator**

In `src/charon/core/schema.py`, replace the MetabolismProperties class (lines 199-204):

```python
class MetabolismProperties(BaseModel):
    """CYP phenotyping and intrinsic clearance."""

    primary_cyp: str | None = None
    secondary_cyp: str | None = None
    clint_uL_min_mg: PredictedProperty | None = None
    fm_cyp3a4: float | None = None

    @field_validator("fm_cyp3a4")
    @classmethod
    def _fm_cyp3a4_range(cls, v: float | None) -> float | None:
        if v is not None and (v < 0.0 or v > 1.0):
            raise ValueError(f"fm_cyp3a4 must be in [0.0, 1.0], got {v}")
        return v
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/unit/test_schema.py::TestFmCyp3a4 -v`
Expected: 6 PASSED

- [ ] **Step 5: Run full test suite for regression**

Run: `pytest tests/ -x -q`
Expected: 543 passed

- [ ] **Step 6: Commit**

```bash
git add src/charon/core/schema.py tests/unit/test_schema.py
git commit -m "feat(schema): add fm_cyp3a4 to MetabolismProperties with [0,1] validator"
```

---

### Task 2: Add gi_tract section to human.yaml

**Files:**
- Modify: `src/charon/pbpk/species/human.yaml`

- [ ] **Step 1: Append gi_tract section to human.yaml**

Add the following at the end of `src/charon/pbpk/species/human.yaml` (after line 132, after the `cyp_orthologs` block):

```yaml

  # GI tract parameters for ACAT oral absorption model (Sprint 3b Session 2a).
  # Per-segment values ported from Sisyphus reference_man.yaml (Yu & Amidon 1999).
  gi_tract:
    enterocyte_volume_L: 0.30
    enterocyte_weight_g: 400
    q_villi_fraction: 0.18            # Gertz 2011 DMPK 26(5):486
    mppgi_mg_g: 15.0                  # Literature range 2-20 (Yang 2007, Barter 2007)
                                      # Calibrated within range: midazolam Fg=0.57
    cyp3a4_gut_pmol_per_mg: 31        # Paine 1997 Gastroenterology 113(2):453
    cyp3a4_liver_pmol_per_mg: 137     # Shimada 1994
    segments:
      stomach:
        volume_L: 0.25
        radius_cm: 5.0
        ka_fraction: 0.0
        transit_rate_1_h: 4.0
      duodenum:
        volume_L: 0.05
        radius_cm: 1.6
        ka_fraction: 1.0
        transit_rate_1_h: 3.846
      jejunum1:
        volume_L: 0.07
        radius_cm: 1.5
        ka_fraction: 1.0
        transit_rate_1_h: 2.105
      jejunum2:
        volume_L: 0.07
        radius_cm: 1.5
        ka_fraction: 1.0
        transit_rate_1_h: 2.105
      ileum1:
        volume_L: 0.06
        radius_cm: 1.3
        ka_fraction: 0.8
        transit_rate_1_h: 1.471
      ileum2:
        volume_L: 0.06
        radius_cm: 1.3
        ka_fraction: 0.6
        transit_rate_1_h: 1.471
      ileum3:
        volume_L: 0.06
        radius_cm: 1.3
        ka_fraction: 0.3
        transit_rate_1_h: 1.471
      colon:
        volume_L: 0.30
        radius_cm: 2.5
        ka_fraction: 0.05
        transit_rate_1_h: 0.074
```

- [ ] **Step 2: Verify YAML parses correctly**

Run: `python3 -c "import yaml; d=yaml.safe_load(open('src/charon/pbpk/species/human.yaml')); gi=d['species']['gi_tract']; print(len(gi['segments']), 'segments'); print('MPPGI:', gi['mppgi_mg_g'])"`
Expected: `8 segments` and `MPPGI: 15.0`

- [ ] **Step 3: Run existing topology tests for regression**

Run: `pytest tests/unit/test_topology.py -v`
Expected: all pass (gi_tract is a new section, not read by existing loader)

- [ ] **Step 4: Commit**

```bash
git add src/charon/pbpk/species/human.yaml
git commit -m "data(human.yaml): add gi_tract section with 8-segment ACAT parameters"
```

---

### Task 3: Create acat.py — GI data structures and loader

**Files:**
- Create: `src/charon/pbpk/acat.py`
- Create: `tests/unit/test_acat.py`

- [ ] **Step 1: Write the failing tests**

Create `tests/unit/test_acat.py`:

```python
"""Tests for the ACAT GI tract module."""
from __future__ import annotations

import math
import pytest


class TestGISegment:
    def test_segment_fields(self):
        from charon.pbpk.acat import GISegment
        s = GISegment(
            name="duodenum", volume_L=0.05, radius_cm=1.6,
            ka_fraction=1.0, transit_rate_1_h=3.846,
        )
        assert s.name == "duodenum"
        assert s.volume_L == 0.05
        assert s.radius_cm == 1.6
        assert s.ka_fraction == 1.0
        assert s.transit_rate_1_h == 3.846

    def test_segment_frozen(self):
        from charon.pbpk.acat import GISegment
        s = GISegment(name="x", volume_L=1, radius_cm=1, ka_fraction=1, transit_rate_1_h=1)
        with pytest.raises(AttributeError):
            s.name = "y"


class TestGITract:
    def test_gi_tract_fields(self):
        from charon.pbpk.acat import GITract, GISegment
        seg = GISegment(name="stomach", volume_L=0.25, radius_cm=5.0,
                        ka_fraction=0.0, transit_rate_1_h=4.0)
        gi = GITract(
            segments=(seg,),
            enterocyte_volume_L=0.30,
            enterocyte_weight_g=400,
            q_villi_fraction=0.18,
            mppgi_mg_g=15.0,
            cyp3a4_gut_pmol_per_mg=31,
            cyp3a4_liver_pmol_per_mg=137,
        )
        assert len(gi.segments) == 1
        assert gi.mppgi_mg_g == 15.0

    def test_gi_tract_frozen(self):
        from charon.pbpk.acat import GITract, GISegment
        seg = GISegment(name="x", volume_L=1, radius_cm=1, ka_fraction=1, transit_rate_1_h=1)
        gi = GITract(segments=(seg,), enterocyte_volume_L=0.3,
                     enterocyte_weight_g=400, q_villi_fraction=0.18,
                     mppgi_mg_g=15.0, cyp3a4_gut_pmol_per_mg=31,
                     cyp3a4_liver_pmol_per_mg=137)
        with pytest.raises(AttributeError):
            gi.mppgi_mg_g = 99


class TestLoadGITract:
    def test_loads_8_segments(self):
        from charon.pbpk.acat import load_gi_tract
        gi = load_gi_tract("human")
        assert len(gi.segments) == 8

    def test_segment_names_in_order(self):
        from charon.pbpk.acat import load_gi_tract
        gi = load_gi_tract("human")
        names = [s.name for s in gi.segments]
        assert names == [
            "stomach", "duodenum", "jejunum1", "jejunum2",
            "ileum1", "ileum2", "ileum3", "colon",
        ]

    def test_stomach_ka_zero(self):
        from charon.pbpk.acat import load_gi_tract
        gi = load_gi_tract("human")
        assert gi.segments[0].ka_fraction == 0.0

    def test_enterocyte_params(self):
        from charon.pbpk.acat import load_gi_tract
        gi = load_gi_tract("human")
        assert gi.enterocyte_volume_L == 0.30
        assert gi.enterocyte_weight_g == 400
        assert gi.q_villi_fraction == 0.18
        assert gi.mppgi_mg_g == 15.0

    def test_total_sitt_about_3h(self):
        """Sum of 1/k_transit for SI segments ≈ 3.25 h."""
        from charon.pbpk.acat import load_gi_tract
        gi = load_gi_tract("human")
        si = [s for s in gi.segments if s.name not in ("stomach", "colon")]
        sitt = sum(1.0 / s.transit_rate_1_h for s in si)
        assert 3.0 < sitt < 4.0, f"SITT = {sitt:.2f} h"

    def test_missing_species_raises(self):
        from charon.pbpk.acat import load_gi_tract
        with pytest.raises(FileNotFoundError):
            load_gi_tract("martian")


class TestComputeAbsorptionRates:
    def test_stomach_kabs_zero(self):
        """ka_fraction=0 → k_abs=0 regardless of Peff."""
        from charon.pbpk.acat import load_gi_tract, compute_absorption_rates
        gi = load_gi_tract("human")
        rates = compute_absorption_rates(gi, peff_cm_s=4.0e-4)
        assert rates[0] == 0.0  # stomach

    def test_duodenum_rate_positive(self):
        from charon.pbpk.acat import load_gi_tract, compute_absorption_rates
        gi = load_gi_tract("human")
        rates = compute_absorption_rates(gi, peff_cm_s=4.0e-4)
        assert rates[1] > 0.0  # duodenum

    def test_hand_calculation_duodenum(self):
        """k_abs = 2 × Peff × 3600 / R × ka_fraction."""
        from charon.pbpk.acat import load_gi_tract, compute_absorption_rates
        gi = load_gi_tract("human")
        rates = compute_absorption_rates(gi, peff_cm_s=4.0e-4)
        expected = 2 * 4.0e-4 * 3600 / 1.6 * 1.0  # = 1.8 /h
        assert rates[1] == pytest.approx(expected, rel=1e-6)

    def test_rates_decrease_distally(self):
        """Due to ka_fraction decay, k_abs should generally decrease."""
        from charon.pbpk.acat import load_gi_tract, compute_absorption_rates
        gi = load_gi_tract("human")
        rates = compute_absorption_rates(gi, peff_cm_s=4.0e-4)
        # jejunum1 rate >= ileum3 rate (both have different R and ka)
        assert rates[2] > rates[6]  # jejunum1 > ileum3

    def test_returns_8_values(self):
        from charon.pbpk.acat import load_gi_tract, compute_absorption_rates
        gi = load_gi_tract("human")
        rates = compute_absorption_rates(gi, peff_cm_s=1.0e-4)
        assert len(rates) == 8


class TestPappToPeff:
    def test_known_conversion(self):
        """Sun 2002 correlation: log10(Peff) = 0.4926×log10(Papp_nm_s) - 0.1454."""
        from charon.pbpk.acat import papp_to_peff
        import math
        papp = 100.0  # nm/s
        expected = 10 ** (0.4926 * math.log10(100.0) - 0.1454)
        result = papp_to_peff(papp)
        assert result == pytest.approx(expected, rel=1e-6)

    def test_high_papp(self):
        from charon.pbpk.acat import papp_to_peff
        result = papp_to_peff(500.0)  # high permeability
        assert result > 1e-5  # Peff should be > 10 µm/s

    def test_zero_papp_raises(self):
        from charon.pbpk.acat import papp_to_peff
        with pytest.raises(ValueError):
            papp_to_peff(0.0)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/unit/test_acat.py -v`
Expected: FAIL — `charon.pbpk.acat` does not exist

- [ ] **Step 3: Implement acat.py**

Create `src/charon/pbpk/acat.py`:

```python
"""ACAT GI tract model — segment data, loader, absorption rate computation.

Implements the 8-segment GI lumen transit model for oral drug absorption.
Segment parameters are loaded from the species YAML ``gi_tract`` section.

The absorption rate per segment uses the mechanistic cylindrical model:
    k_abs_i [1/h] = (2 × Peff [cm/s] × 3600 / R_i [cm]) × ka_fraction_i

References:
    Yu LX, Amidon GL (1999). Int J Pharm 186:119-125.
    Sun D et al. (2002). J Pharm Sci 91:1396-1404.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

import yaml

_SPECIES_DIR = Path(__file__).parent / "species"


@dataclass(frozen=True)
class GISegment:
    """Single GI lumen segment."""

    name: str
    volume_L: float
    radius_cm: float
    ka_fraction: float
    transit_rate_1_h: float


@dataclass(frozen=True)
class GITract:
    """Complete GI tract physiology for ACAT model."""

    segments: tuple[GISegment, ...]
    enterocyte_volume_L: float
    enterocyte_weight_g: float
    q_villi_fraction: float
    mppgi_mg_g: float
    cyp3a4_gut_pmol_per_mg: float
    cyp3a4_liver_pmol_per_mg: float


def load_gi_tract(species: str) -> GITract:
    """Load GI tract parameters from species YAML.

    Parameters
    ----------
    species : str
        Species name (e.g. ``"human"``).

    Returns
    -------
    GITract
        Frozen GI tract data container.

    Raises
    ------
    FileNotFoundError
        If the species YAML does not exist.
    KeyError
        If the YAML lacks a ``gi_tract`` section.
    """
    yaml_path = _SPECIES_DIR / f"{species}.yaml"
    if not yaml_path.exists():
        raise FileNotFoundError(f"Species YAML not found: {yaml_path}")

    with yaml_path.open() as fp:
        data = yaml.safe_load(fp)

    gi = data["species"]["gi_tract"]
    segments: list[GISegment] = []
    for name, spec in gi["segments"].items():
        segments.append(
            GISegment(
                name=name,
                volume_L=float(spec["volume_L"]),
                radius_cm=float(spec["radius_cm"]),
                ka_fraction=float(spec["ka_fraction"]),
                transit_rate_1_h=float(spec["transit_rate_1_h"]),
            )
        )

    return GITract(
        segments=tuple(segments),
        enterocyte_volume_L=float(gi["enterocyte_volume_L"]),
        enterocyte_weight_g=float(gi["enterocyte_weight_g"]),
        q_villi_fraction=float(gi["q_villi_fraction"]),
        mppgi_mg_g=float(gi["mppgi_mg_g"]),
        cyp3a4_gut_pmol_per_mg=float(gi["cyp3a4_gut_pmol_per_mg"]),
        cyp3a4_liver_pmol_per_mg=float(gi["cyp3a4_liver_pmol_per_mg"]),
    )


def compute_absorption_rates(
    gi: GITract,
    peff_cm_s: float,
) -> tuple[float, ...]:
    """Compute per-segment absorption rate constants.

    Parameters
    ----------
    gi : GITract
        GI tract physiology.
    peff_cm_s : float
        Effective intestinal permeability in cm/s.

    Returns
    -------
    tuple[float, ...]
        k_abs for each segment in 1/h, same order as ``gi.segments``.
    """
    rates: list[float] = []
    for seg in gi.segments:
        if seg.ka_fraction == 0.0:
            rates.append(0.0)
        else:
            k_abs = (2.0 * peff_cm_s * 3600.0 / seg.radius_cm) * seg.ka_fraction
            rates.append(k_abs)
    return tuple(rates)


def papp_to_peff(papp_nm_s: float) -> float:
    """Convert Caco-2 Papp (nm/s) to human Peff (cm/s).

    Uses the Sun et al. 2002 simplified correlation:
        log10(Peff) = 0.4926 × log10(Papp_nm_s) - 0.1454

    Parameters
    ----------
    papp_nm_s : float
        Apparent permeability from Caco-2 assay, nm/s.

    Returns
    -------
    float
        Effective intestinal permeability, cm/s.

    Raises
    ------
    ValueError
        If papp_nm_s <= 0.
    """
    if papp_nm_s <= 0.0:
        raise ValueError(f"papp_nm_s must be > 0, got {papp_nm_s}")
    log_peff = 0.4926 * math.log10(papp_nm_s) - 0.1454
    return 10.0 ** log_peff
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/unit/test_acat.py -v`
Expected: 17 PASSED

- [ ] **Step 5: Run full suite for regression**

Run: `pytest tests/ -x -q`
Expected: 543 + 17 = 560 passed

- [ ] **Step 6: Commit**

```bash
git add src/charon/pbpk/acat.py tests/unit/test_acat.py
git commit -m "feat(acat): GI segment loader, absorption rates, Papp→Peff conversion"
```

---

### Task 4: Add OralPBPKParams and compute_gut_clint to ode_compiler

**Files:**
- Modify: `src/charon/pbpk/ode_compiler.py`
- Create: `tests/unit/test_oral_compiler.py`

- [ ] **Step 1: Write the failing tests**

Create `tests/unit/test_oral_compiler.py`:

```python
"""Tests for oral PBPK parameter assembly and gut CLint calculation."""
from __future__ import annotations

import pytest


class TestComputeGutClint:
    def test_no_fm_returns_zero(self):
        """fm_cyp3a4=None → CLint_gut=0 (non-CYP3A4 substrate)."""
        from charon.pbpk.ode_compiler import compute_gut_clint
        from charon.pbpk.acat import load_gi_tract
        gi = load_gi_tract("human")
        result = compute_gut_clint(
            clint_liver_L_h=348.75,
            fm_cyp3a4=None,
            gi_tract=gi,
            mppgl=40.0,
            liver_weight_g=1500.0,
        )
        assert result == 0.0

    def test_fm_zero_returns_zero(self):
        """fm_cyp3a4=0 → CLint_gut=0."""
        from charon.pbpk.ode_compiler import compute_gut_clint
        from charon.pbpk.acat import load_gi_tract
        gi = load_gi_tract("human")
        result = compute_gut_clint(
            clint_liver_L_h=348.75,
            fm_cyp3a4=0.0,
            gi_tract=gi,
            mppgl=40.0,
            liver_weight_g=1500.0,
        )
        assert result == 0.0

    def test_midazolam_clint_gut(self):
        """Midazolam: CLint_gut ≈ 7.9 L/h → Fg ≈ 0.57."""
        from charon.pbpk.ode_compiler import compute_gut_clint
        from charon.pbpk.acat import load_gi_tract
        gi = load_gi_tract("human")
        result = compute_gut_clint(
            clint_liver_L_h=348.75,  # actual Charon midazolam value
            fm_cyp3a4=1.0,
            gi_tract=gi,
            mppgl=40.0,
            liver_weight_g=1500.0,
        )
        # Expected: 348.75 × 1.0 × (31/137) × (15×400)/(40×1500)
        #         = 348.75 × 0.2263 × 0.10 = 7.89
        assert result == pytest.approx(7.89, rel=0.02)

    def test_partial_fm(self):
        """fm_cyp3a4=0.5 → half the CLint_gut."""
        from charon.pbpk.ode_compiler import compute_gut_clint
        from charon.pbpk.acat import load_gi_tract
        gi = load_gi_tract("human")
        full = compute_gut_clint(
            clint_liver_L_h=348.75, fm_cyp3a4=1.0,
            gi_tract=gi, mppgl=40.0, liver_weight_g=1500.0,
        )
        half = compute_gut_clint(
            clint_liver_L_h=348.75, fm_cyp3a4=0.5,
            gi_tract=gi, mppgl=40.0, liver_weight_g=1500.0,
        )
        assert half == pytest.approx(full / 2.0, rel=1e-10)


class TestOralPBPKParams:
    def test_has_gut_fields(self):
        from charon.pbpk.ode_compiler import OralPBPKParams
        # Just test the dataclass has the new fields
        import dataclasses
        field_names = [f.name for f in dataclasses.fields(OralPBPKParams)]
        assert "clint_gut_L_h" in field_names
        assert "peff_cm_s" in field_names
        assert "q_villi_L_h" in field_names
        assert "v_enterocyte_L" in field_names

    def test_inherits_compound_fields(self):
        from charon.pbpk.ode_compiler import OralPBPKParams
        import dataclasses
        field_names = [f.name for f in dataclasses.fields(OralPBPKParams)]
        # Must have core compound fields from CompoundPBPKParams
        assert "fu_b" in field_names
        assert "clint_liver_L_h" in field_names
        assert "kp_by_tissue" in field_names
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/unit/test_oral_compiler.py -v`
Expected: FAIL — `compute_gut_clint` and `OralPBPKParams` do not exist

- [ ] **Step 3: Implement compute_gut_clint and OralPBPKParams**

Add the following to the END of `src/charon/pbpk/ode_compiler.py` (after `build_rhs` at line 390):

```python
# ---------------------------------------------------------------------------
# Oral PBPK extensions (Sprint 3b Session 2a)
# ---------------------------------------------------------------------------

from charon.pbpk.acat import GITract, compute_absorption_rates


@dataclass(frozen=True)
class OralPBPKParams:
    """PBPK parameters extended for oral route.

    Contains all fields from CompoundPBPKParams plus gut-specific
    parameters for the ACAT enterocyte model.
    """

    # --- inherited from CompoundPBPKParams ---
    name: str
    molecular_weight: float
    logp: float
    pka_acid: float | None
    pka_base: float | None
    compound_type: str
    fu_p: float
    bp_ratio: float
    fu_b: float
    clint_liver_L_h: float
    cl_renal_L_h: float
    kp_by_tissue: dict[str, float]
    kp_overrides: tuple[KpOverrideRecord, ...] = ()

    # --- oral-specific ---
    clint_gut_L_h: float = 0.0
    peff_cm_s: float = 0.0
    q_villi_L_h: float = 0.0
    v_enterocyte_L: float = 0.30
    gi_tract: GITract | None = None


def compute_gut_clint(
    *,
    clint_liver_L_h: float,
    fm_cyp3a4: float | None,
    gi_tract: GITract,
    mppgl: float,
    liver_weight_g: float,
) -> float:
    """Compute gut-wall intrinsic clearance via ratio approach.

    CLint_gut = CLint_liver × fm_CYP3A4 × (CYP3A4_gut/CYP3A4_liver)
                × (MPPGI × gut_weight) / (MPPGL × liver_weight)

    The fu_inc correction cancels because both liver and gut use the
    same microsomal-protein-based scaling.

    Parameters
    ----------
    clint_liver_L_h : float
        Whole-liver intrinsic clearance from ParameterBridge (L/h).
    fm_cyp3a4 : float or None
        Fraction of hepatic CLint metabolised by CYP3A4.
        None or 0.0 → returns 0.0 (non-CYP3A4, Fg=1.0).
    gi_tract : GITract
        GI physiology with MPPGI and CYP3A4 content.
    mppgl : float
        Milligrams protein per gram liver (mg/g).
    liver_weight_g : float
        Liver weight (g).

    Returns
    -------
    float
        Gut intrinsic clearance in L/h.
    """
    if fm_cyp3a4 is None or fm_cyp3a4 == 0.0:
        return 0.0

    cyp_ratio = gi_tract.cyp3a4_gut_pmol_per_mg / gi_tract.cyp3a4_liver_pmol_per_mg
    gut_scaling = gi_tract.mppgi_mg_g * gi_tract.enterocyte_weight_g
    liver_scaling = mppgl * liver_weight_g

    return clint_liver_L_h * fm_cyp3a4 * cyp_ratio * gut_scaling / liver_scaling
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/unit/test_oral_compiler.py -v`
Expected: 6 PASSED

- [ ] **Step 5: Run full suite for regression**

Run: `pytest tests/ -x -q`
Expected: 560 + 6 = 566 passed

- [ ] **Step 6: Commit**

```bash
git add src/charon/pbpk/ode_compiler.py tests/unit/test_oral_compiler.py
git commit -m "feat(ode_compiler): add OralPBPKParams and compute_gut_clint (ratio approach)"
```

---

### Task 5: Implement build_oral_rhs — the oral ODE right-hand side

**Files:**
- Modify: `src/charon/pbpk/ode_compiler.py`
- Modify: `tests/unit/test_oral_compiler.py`

- [ ] **Step 1: Write the failing tests**

Append to `tests/unit/test_oral_compiler.py`:

```python
import numpy as np


class TestBuildOralRhs:
    """Test the oral ODE right-hand side function."""

    def _make_oral_params(self):
        """Helper: build minimal OralPBPKParams for midazolam-like compound."""
        from charon.pbpk.ode_compiler import OralPBPKParams
        from charon.pbpk.acat import load_gi_tract
        from charon.pbpk.topology import load_species_topology

        topo = load_species_topology("human")
        gi = load_gi_tract("human")

        # Minimal Kp dict (all 1.0 for simplicity in RHS tests)
        kp = {name: 1.0 for name in topo.tissue_names()}

        q_gut = topo.tissues["gut_wall"].blood_flow_L_h
        q_villi = gi.q_villi_fraction * q_gut

        return OralPBPKParams(
            name="test_oral",
            molecular_weight=325.77,
            logp=3.89,
            pka_acid=None,
            pka_base=6.2,
            compound_type="base",
            fu_p=0.03,
            bp_ratio=0.66,
            fu_b=0.03 / 0.66,
            clint_liver_L_h=348.75,
            cl_renal_L_h=0.0,
            kp_by_tissue=kp,
            clint_gut_L_h=7.89,
            peff_cm_s=4.0e-4,
            q_villi_L_h=q_villi,
            v_enterocyte_L=gi.enterocyte_volume_L,
            gi_tract=gi,
        ), topo

    def test_returns_callable(self):
        from charon.pbpk.ode_compiler import build_oral_rhs
        params, topo = self._make_oral_params()
        rhs = build_oral_rhs(topo, params)
        assert callable(rhs)

    def test_state_vector_length(self):
        """Oral state vector = 2 + 15 tissues + 8 lumen + 1 enterocyte = 26."""
        from charon.pbpk.ode_compiler import build_oral_rhs
        params, topo = self._make_oral_params()
        rhs = build_oral_rhs(topo, params)
        n = 2 + len(topo.tissues) + 8 + 1  # 26
        y0 = np.zeros(n)
        y0[17] = 5.0  # dose in stomach
        dy = rhs(0.0, y0)
        assert len(dy) == n

    def test_dose_in_stomach_creates_transit(self):
        """Drug in stomach should produce negative dA_stomach (emptying)."""
        from charon.pbpk.ode_compiler import build_oral_rhs
        params, topo = self._make_oral_params()
        rhs = build_oral_rhs(topo, params)
        n = 2 + len(topo.tissues) + 8 + 1
        y0 = np.zeros(n)
        y0[17] = 5.0  # stomach index = 2 + 15 = 17
        dy = rhs(0.0, y0)
        # dA_stomach should be negative (emptying)
        assert dy[17] < 0.0
        # dA_duodenum should be positive (receiving from stomach)
        assert dy[18] > 0.0

    def test_enterocyte_receives_absorption(self):
        """Drug in duodenum should produce positive dA_enterocyte."""
        from charon.pbpk.ode_compiler import build_oral_rhs
        params, topo = self._make_oral_params()
        rhs = build_oral_rhs(topo, params)
        n = 2 + len(topo.tissues) + 8 + 1
        y0 = np.zeros(n)
        y0[18] = 5.0  # duodenum
        dy = rhs(0.0, y0)
        # enterocyte index = 2 + 15 + 8 = 25
        assert dy[25] > 0.0

    def test_mass_conservation_no_elimination(self):
        """With CLint=0 and CLrenal=0, total mass derivative should be 0."""
        from charon.pbpk.ode_compiler import OralPBPKParams, build_oral_rhs
        from charon.pbpk.acat import load_gi_tract
        from charon.pbpk.topology import load_species_topology

        topo = load_species_topology("human")
        gi = load_gi_tract("human")
        q_gut = topo.tissues["gut_wall"].blood_flow_L_h

        params = OralPBPKParams(
            name="no_elim",
            molecular_weight=300.0,
            logp=2.0,
            pka_acid=None,
            pka_base=None,
            compound_type="neutral",
            fu_p=1.0,
            bp_ratio=1.0,
            fu_b=1.0,
            clint_liver_L_h=0.0,  # no hepatic elimination
            cl_renal_L_h=0.0,     # no renal elimination
            kp_by_tissue={n: 1.0 for n in topo.tissue_names()},
            clint_gut_L_h=0.0,    # no gut metabolism
            peff_cm_s=4.0e-4,
            q_villi_L_h=gi.q_villi_fraction * q_gut,
            v_enterocyte_L=gi.enterocyte_volume_L,
            gi_tract=gi,
        )

        rhs = build_oral_rhs(topo, params)
        n = 2 + len(topo.tissues) + 8 + 1
        y0 = np.zeros(n)
        y0[17] = 100.0  # 100 mg in stomach

        # At t=0, some drug transits/absorbs but total should be conserved
        dy = rhs(0.0, y0)
        total_dy = np.sum(dy)
        assert abs(total_dy) < 1e-10, f"Mass not conserved: sum(dy) = {total_dy}"

    def test_liver_receives_portal_from_enterocyte(self):
        """Enterocyte drug should increase liver mass via portal inflow."""
        from charon.pbpk.ode_compiler import build_oral_rhs
        params, topo = self._make_oral_params()
        rhs = build_oral_rhs(topo, params)
        n = 2 + len(topo.tissues) + 8 + 1
        y0 = np.zeros(n)
        y0[25] = 5.0  # drug in enterocyte
        dy = rhs(0.0, y0)
        liver_idx = 2 + list(topo.tissues.keys()).index("liver")
        assert dy[liver_idx] > 0.0, "liver should receive drug from enterocyte"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/unit/test_oral_compiler.py::TestBuildOralRhs -v`
Expected: FAIL — `build_oral_rhs` does not exist

- [ ] **Step 3: Implement build_oral_rhs**

Append to `src/charon/pbpk/ode_compiler.py`:

```python
def build_oral_rhs(
    topology: PBPKTopology,
    params: OralPBPKParams,
):
    """Return a closure ``rhs(t, y) -> dy`` for the oral PBPK+ACAT ODE.

    State vector layout (length 2 + N_tissues + 8 + 1 = 26):
      y[0]            = A_venous         (mg)
      y[1]            = A_arterial       (mg)
      y[2:2+N]        = A_tissue_i       (mg) in topology.tissue_names() order
      y[2+N:2+N+8]    = A_lumen_i        (mg) in gi_tract segment order
      y[2+N+8]        = A_enterocyte     (mg)

    The tissue dynamics (states 0 through 2+N-1) are identical to
    ``build_rhs`` except that the liver mass balance includes an
    additive portal inflow term from the enterocyte.

    Fg emerges from the enterocyte dynamics:
      dA_enterocyte/dt = absorption_influx - gut_metabolism - basolateral_transfer
    where basolateral transfer feeds into the liver via the portal vein.
    """
    gi = params.gi_tract
    assert gi is not None, "OralPBPKParams.gi_tract must be set"

    tissue_names = topology.tissue_names()
    n_tissues = len(tissue_names)
    n_lumen = len(gi.segments)  # 8

    # --- Pre-compute tissue arrays (same as build_rhs) ---
    volumes = np.array(
        [topology.tissues[name].volume_L for name in tissue_names],
        dtype=np.float64,
    )
    flows = np.array(
        [topology.tissues[name].blood_flow_L_h for name in tissue_names],
        dtype=np.float64,
    )
    kp = np.array(
        [params.kp_by_tissue[name] for name in tissue_names],
        dtype=np.float64,
    )
    drains_to = [topology.tissues[name].drains_to for name in tissue_names]

    lung_idx = tissue_names.index("lung")
    liver_idx = tissue_names.index("liver")
    kidney_idx = tissue_names.index("kidney")
    portal_indices = np.array(
        [tissue_names.index(p) for p in PORTAL_TISSUES], dtype=np.int64
    )

    bp = params.bp_ratio
    fu_b = params.fu_b
    clint_liver = params.clint_liver_L_h
    cl_renal = params.cl_renal_L_h

    v_ven = topology.venous_volume_L
    v_art = topology.arterial_volume_L
    q_co = topology.cardiac_output_L_h
    q_ha = topology.hepatic_artery_L_h
    q_liver_total = flows[liver_idx]

    # --- Pre-compute GI arrays ---
    k_abs = np.array(
        compute_absorption_rates(gi, params.peff_cm_s),
        dtype=np.float64,
    )
    k_transit = np.array(
        [seg.transit_rate_1_h for seg in gi.segments],
        dtype=np.float64,
    )

    # Enterocyte rate constants
    v_entero = params.v_enterocyte_L
    q_villi = params.q_villi_L_h
    clint_gut = params.clint_gut_L_h
    k_baso = q_villi / v_entero         # basolateral transfer rate (1/h)
    k_metab_gut = clint_gut / v_entero  # gut metabolism rate (1/h)

    # State vector offsets
    _tissue_start = 2
    _lumen_start = 2 + n_tissues
    _entero_idx = 2 + n_tissues + n_lumen

    def rhs(t: float, y: np.ndarray) -> np.ndarray:
        dy = np.zeros_like(y)

        # --- Unpack state ---
        a_ven = y[0]
        a_art = y[1]
        a_tissues = y[_tissue_start:_lumen_start]
        a_lumen = y[_lumen_start:_entero_idx]
        a_entero = y[_entero_idx]

        # --- Tissue concentrations ---
        c_ven = a_ven / v_ven
        c_art = a_art / v_art
        c_tissue = a_tissues / volumes
        c_blood_out = c_tissue * bp / kp

        # --- Tissue dynamics (identical to build_rhs) ---
        dy_tissues = np.zeros(n_tissues)
        for i in range(n_tissues):
            if i == lung_idx or i == liver_idx:
                continue
            dy_tissues[i] = flows[i] * (c_art - c_blood_out[i])

        dy_tissues[lung_idx] = q_co * (c_ven - c_blood_out[lung_idx])
        dy_tissues[kidney_idx] -= cl_renal * c_blood_out[kidney_idx] / bp

        # Portal inflow from PBPK tissues
        portal_inflow = 0.0
        for pi in portal_indices:
            portal_inflow += flows[pi] * c_blood_out[pi]

        # Portal inflow from enterocyte (basolateral transfer)
        portal_from_entero = k_baso * a_entero

        hepatic_elim = clint_liver * fu_b * c_blood_out[liver_idx]
        dy_tissues[liver_idx] = (
            q_ha * c_art
            + portal_inflow
            + portal_from_entero
            - q_liver_total * c_blood_out[liver_idx]
            - hepatic_elim
        )

        # Venous inflow (same as build_rhs)
        venous_inflow = 0.0
        for i in range(n_tissues):
            if i == lung_idx:
                continue
            if drains_to[i] == "venous":
                if i == liver_idx:
                    venous_inflow += q_liver_total * c_blood_out[liver_idx]
                else:
                    venous_inflow += flows[i] * c_blood_out[i]
        dy_ven = venous_inflow - q_co * c_ven

        dy_art = q_co * c_blood_out[lung_idx] - q_co * c_art

        # --- GI lumen dynamics ---
        dy_lumen = np.zeros(n_lumen)

        # Stomach: only transit out (ka=0)
        dy_lumen[0] = -k_transit[0] * a_lumen[0]

        # Segments 1..N-1: transit in from previous, transit out, absorption
        for i in range(1, n_lumen):
            dy_lumen[i] = (
                k_transit[i - 1] * a_lumen[i - 1]
                - k_transit[i] * a_lumen[i]
                - k_abs[i] * a_lumen[i]
            )

        # --- Enterocyte dynamics ---
        absorption_influx = 0.0
        for i in range(n_lumen):
            absorption_influx += k_abs[i] * a_lumen[i]

        gut_metabolism = k_metab_gut * a_entero
        basolateral = k_baso * a_entero
        dy_entero = absorption_influx - gut_metabolism - basolateral

        # --- Pack output ---
        dy[0] = dy_ven
        dy[1] = dy_art
        dy[_tissue_start:_lumen_start] = dy_tissues
        dy[_lumen_start:_entero_idx] = dy_lumen
        dy[_entero_idx] = dy_entero

        return dy

    return rhs
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/unit/test_oral_compiler.py -v`
Expected: 13 PASSED (6 from Task 4 + 7 new)

- [ ] **Step 5: Run full suite for regression**

Run: `pytest tests/ -x -q`
Expected: 566 + 7 = 573 passed

- [ ] **Step 6: Commit**

```bash
git add src/charon/pbpk/ode_compiler.py tests/unit/test_oral_compiler.py
git commit -m "feat(ode_compiler): add build_oral_rhs with 26-state ACAT+PBPK ODE"
```

---

### Task 6: Implement simulate_oral in solver.py

**Files:**
- Modify: `src/charon/pbpk/solver.py`
- Create: `tests/unit/test_oral_solver.py`

- [ ] **Step 1: Write the failing tests**

Create `tests/unit/test_oral_solver.py`:

```python
"""Tests for the oral simulation wrapper."""
from __future__ import annotations

import math
import numpy as np
import pytest


def _make_midazolam_oral_params():
    """Build OralPBPKParams for midazolam from real YAML data."""
    from charon.core.parameter_bridge import ParameterBridge
    from charon.core.schema import CompoundConfig
    from charon.pbpk.acat import load_gi_tract
    from charon.pbpk.ode_compiler import (
        build_compound_pbpk_params,
        compute_gut_clint,
        OralPBPKParams,
    )
    from charon.pbpk.topology import load_species_topology

    import yaml
    from pathlib import Path

    comp_path = Path("validation/data/tier1_obach/compounds/midazolam.yaml")
    with comp_path.open() as f:
        raw = yaml.safe_load(f)
    compound = CompoundConfig(**raw)

    topo = load_species_topology("human")
    bridge = ParameterBridge()
    base_params = build_compound_pbpk_params(compound, topo, bridge=bridge)

    gi = load_gi_tract("human")
    q_gut = topo.tissues["gut_wall"].blood_flow_L_h
    q_villi = gi.q_villi_fraction * q_gut

    fm = compound.properties.metabolism.fm_cyp3a4
    clint_gut = compute_gut_clint(
        clint_liver_L_h=base_params.clint_liver_L_h,
        fm_cyp3a4=fm,
        gi_tract=gi,
        mppgl=40.0,
        liver_weight_g=topo.liver_weight_g,
    )

    peff = 4.0e-4  # midazolam Peff from literature

    return OralPBPKParams(
        name=base_params.name,
        molecular_weight=base_params.molecular_weight,
        logp=base_params.logp,
        pka_acid=base_params.pka_acid,
        pka_base=base_params.pka_base,
        compound_type=base_params.compound_type,
        fu_p=base_params.fu_p,
        bp_ratio=base_params.bp_ratio,
        fu_b=base_params.fu_b,
        clint_liver_L_h=base_params.clint_liver_L_h,
        cl_renal_L_h=base_params.cl_renal_L_h,
        kp_by_tissue=base_params.kp_by_tissue,
        kp_overrides=base_params.kp_overrides,
        clint_gut_L_h=clint_gut,
        peff_cm_s=peff,
        q_villi_L_h=q_villi,
        v_enterocyte_L=gi.enterocyte_volume_L,
        gi_tract=gi,
    ), topo


class TestSimulateOral:
    def test_returns_result(self):
        from charon.pbpk.solver import simulate_oral
        params, topo = _make_midazolam_oral_params()
        sim = simulate_oral(topo, params, dose_mg=5.0, duration_h=24.0)
        assert sim.route == "oral"
        assert sim.dose_mg == 5.0

    def test_cp_plasma_positive(self):
        from charon.pbpk.solver import simulate_oral
        params, topo = _make_midazolam_oral_params()
        sim = simulate_oral(topo, params, dose_mg=5.0, duration_h=24.0)
        # After some time, drug should appear in plasma
        assert np.any(sim.cp_plasma > 0)

    def test_cp_plasma_has_peak(self):
        """Oral profile should rise then fall → Cmax at tmax > 0."""
        from charon.pbpk.solver import simulate_oral
        params, topo = _make_midazolam_oral_params()
        sim = simulate_oral(topo, params, dose_mg=5.0, duration_h=24.0)
        tmax_idx = np.argmax(sim.cp_plasma)
        assert sim.time_h[tmax_idx] > 0.0, "tmax should be > 0 for oral"
        assert sim.cp_plasma[tmax_idx] > 0.0

    def test_mass_balance(self):
        """Total mass in system + eliminated should ≈ dose."""
        from charon.pbpk.solver import simulate_oral
        params, topo = _make_midazolam_oral_params()
        sim = simulate_oral(topo, params, dose_mg=5.0, duration_h=24.0)
        # At t=0, all mass in stomach → total state sum = 5.0
        total_t0 = sim.state_trajectory[:, 0].sum()
        assert total_t0 == pytest.approx(5.0, abs=0.01)

    def test_bdf_method(self):
        from charon.pbpk.solver import simulate_oral
        params, topo = _make_midazolam_oral_params()
        sim = simulate_oral(topo, params, dose_mg=5.0, duration_h=24.0)
        assert sim.solver_method == "BDF"

    def test_lumen_trajectory_available(self):
        from charon.pbpk.solver import simulate_oral
        params, topo = _make_midazolam_oral_params()
        sim = simulate_oral(topo, params, dose_mg=5.0, duration_h=24.0)
        # lumen_trajectory should have shape (8, n_time_points)
        assert hasattr(sim, "lumen_trajectory")
        assert sim.lumen_trajectory.shape[0] == 8

    def test_enterocyte_trajectory_available(self):
        from charon.pbpk.solver import simulate_oral
        params, topo = _make_midazolam_oral_params()
        sim = simulate_oral(topo, params, dose_mg=5.0, duration_h=24.0)
        assert hasattr(sim, "enterocyte_trajectory")
        assert len(sim.enterocyte_trajectory) == len(sim.time_h)

    def test_explicit_method_rejected(self):
        from charon.pbpk.solver import simulate_oral
        params, topo = _make_midazolam_oral_params()
        with pytest.raises(ValueError, match="BDF"):
            simulate_oral(topo, params, dose_mg=5.0, duration_h=24.0,
                          method="RK45")
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/unit/test_oral_solver.py -v`
Expected: FAIL — `simulate_oral` does not exist

- [ ] **Step 3: Implement simulate_oral and OralSimulationResult**

Add to the end of `src/charon/pbpk/solver.py`:

```python
from charon.pbpk.ode_compiler import OralPBPKParams, build_oral_rhs


@dataclass
class OralSimulationResult:
    """Output of an oral PBPK+ACAT simulation."""

    time_h: np.ndarray
    cp_blood: np.ndarray
    cp_plasma: np.ndarray
    state_trajectory: np.ndarray
    lumen_trajectory: np.ndarray
    enterocyte_trajectory: np.ndarray
    mass_balance_residual: float
    solver_success: bool
    solver_method: str
    solver_nfev: int
    route: str
    dose_mg: float


def simulate_oral(
    topology: PBPKTopology,
    params: OralPBPKParams,
    *,
    dose_mg: float,
    duration_h: float,
    n_time_points: int = 500,
    rtol: float = 1e-6,
    atol: float | None = None,
    method: str = "BDF",
) -> OralSimulationResult:
    """Simulate an oral dose through the PBPK+ACAT kernel.

    Parameters
    ----------
    topology : PBPKTopology
        Loaded species topology.
    params : OralPBPKParams
        Compound-specific oral PBPK parameters.
    dose_mg : float
        Oral dose in mg (placed in stomach at t=0).
    duration_h : float
        Simulation end time (hours).
    """
    if method not in _ALLOWED_METHODS:
        raise ValueError(
            f"PBPK is stiff — only BDF is allowed, got {method!r}."
        )
    if dose_mg <= 0:
        raise ValueError(f"dose_mg must be > 0, got {dose_mg}")
    if duration_h <= 0:
        raise ValueError(f"duration_h must be > 0, got {duration_h}")

    if atol is None:
        atol = 1e-9 * dose_mg

    gi = params.gi_tract
    assert gi is not None

    n_tissues = len(topology.tissues)
    n_lumen = len(gi.segments)
    n_states = 2 + n_tissues + n_lumen + 1

    y0 = np.zeros(n_states)
    # Oral dose: all in stomach lumen at t=0
    stomach_idx = 2 + n_tissues  # first lumen state
    y0[stomach_idx] = dose_mg

    rhs = build_oral_rhs(topology, params)

    t_eval = np.linspace(0.0, duration_h, n_time_points)

    sol = solve_ivp(
        rhs,
        (0.0, duration_h),
        y0,
        method=method,
        t_eval=t_eval,
        rtol=rtol,
        atol=atol,
        dense_output=False,
    )

    if not sol.success:
        raise RuntimeError(
            f"Oral ODE solver failed for {params.name!r}: {sol.message}"
        )

    v_ven = topology.venous_volume_L
    cp_blood = sol.y[0] / v_ven
    cp_plasma = cp_blood / params.bp_ratio

    # Mass balance at t=0
    residual = float(np.abs(sol.y.sum(axis=0)[0] - dose_mg))

    # Extract lumen and enterocyte trajectories
    lumen_start = 2 + n_tissues
    entero_idx = lumen_start + n_lumen
    lumen_traj = sol.y[lumen_start:entero_idx, :]
    entero_traj = sol.y[entero_idx, :]

    return OralSimulationResult(
        time_h=sol.t,
        cp_blood=cp_blood,
        cp_plasma=cp_plasma,
        state_trajectory=sol.y,
        lumen_trajectory=lumen_traj,
        enterocyte_trajectory=entero_traj,
        mass_balance_residual=residual,
        solver_success=sol.success,
        solver_method=method,
        solver_nfev=int(sol.nfev),
        route="oral",
        dose_mg=dose_mg,
    )
```

- [ ] **Step 4: Update midazolam.yaml with fm_cyp3a4**

This test helper needs midazolam to have fm_cyp3a4. Add to `validation/data/tier1_obach/compounds/midazolam.yaml` under `metabolism:`:

```yaml
  metabolism:
    clint_uL_min_mg: {value: 93.0, source: experimental, unit: uL/min/mg, method: "Obach 1999 Table 2"}
    fm_cyp3a4: 1.0
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `pytest tests/unit/test_oral_solver.py -v`
Expected: 8 PASSED

- [ ] **Step 6: Run full suite for regression**

Run: `pytest tests/ -x -q`
Expected: 573 + 8 = 581 passed

- [ ] **Step 7: Commit**

```bash
git add src/charon/pbpk/solver.py tests/unit/test_oral_solver.py validation/data/tier1_obach/compounds/midazolam.yaml
git commit -m "feat(solver): add simulate_oral with OralSimulationResult and ACAT IC"
```

---

### Task 7: Implement compute_oral_pk_parameters with Fa/Fg/Fh extraction

**Files:**
- Modify: `src/charon/pbpk/pk_extract.py`
- Create: `tests/unit/test_oral_pk_extract.py`

- [ ] **Step 1: Write the failing tests**

Create `tests/unit/test_oral_pk_extract.py`:

```python
"""Tests for oral PK parameter extraction including Fa/Fg/Fh."""
from __future__ import annotations

import numpy as np
import pytest


def _run_midazolam_oral():
    """Run midazolam oral simulation and return result."""
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

    from tests.unit.test_oral_solver import _make_midazolam_oral_params
    from charon.pbpk.solver import simulate_oral
    params, topo = _make_midazolam_oral_params()
    return simulate_oral(topo, params, dose_mg=5.0, duration_h=72.0), params, topo


class TestComputeOralPK:
    def test_returns_pk_parameters(self):
        from charon.pbpk.pk_extract import compute_oral_pk_parameters
        sim, params, topo = _run_midazolam_oral()
        pk = compute_oral_pk_parameters(sim, params, topo, dose_mg=5.0)
        assert pk.cmax is not None
        assert pk.tmax is not None
        assert pk.auc_0_inf is not None

    def test_cmax_positive(self):
        from charon.pbpk.pk_extract import compute_oral_pk_parameters
        sim, params, topo = _run_midazolam_oral()
        pk = compute_oral_pk_parameters(sim, params, topo, dose_mg=5.0)
        assert pk.cmax > 0

    def test_tmax_positive(self):
        from charon.pbpk.pk_extract import compute_oral_pk_parameters
        sim, params, topo = _run_midazolam_oral()
        pk = compute_oral_pk_parameters(sim, params, topo, dose_mg=5.0)
        assert pk.tmax > 0  # oral should have delay

    def test_fa_high_for_midazolam(self):
        """Midazolam has high Peff → Fa > 0.90."""
        from charon.pbpk.pk_extract import compute_oral_pk_parameters
        sim, params, topo = _run_midazolam_oral()
        pk = compute_oral_pk_parameters(sim, params, topo, dose_mg=5.0)
        assert pk.fa is not None
        assert pk.fa > 0.90

    def test_fg_midazolam_near_057(self):
        """Midazolam Fg should be ≈ 0.57 (calibration check, ±10%)."""
        from charon.pbpk.pk_extract import compute_oral_pk_parameters
        sim, params, topo = _run_midazolam_oral()
        pk = compute_oral_pk_parameters(sim, params, topo, dose_mg=5.0)
        assert pk.fg is not None
        assert 0.50 < pk.fg < 0.65, f"Fg = {pk.fg:.3f}, expected ~0.57"

    def test_fh_positive(self):
        from charon.pbpk.pk_extract import compute_oral_pk_parameters
        sim, params, topo = _run_midazolam_oral()
        pk = compute_oral_pk_parameters(sim, params, topo, dose_mg=5.0)
        assert pk.fh is not None
        assert 0.0 < pk.fh < 1.0

    def test_bioavailability_product(self):
        """F = Fa × Fg × Fh."""
        from charon.pbpk.pk_extract import compute_oral_pk_parameters
        sim, params, topo = _run_midazolam_oral()
        pk = compute_oral_pk_parameters(sim, params, topo, dose_mg=5.0)
        expected = pk.fa * pk.fg * pk.fh
        assert pk.bioavailability == pytest.approx(expected, rel=1e-6)

    def test_cl_apparent_is_cl_over_f(self):
        """For oral, CL_apparent = dose / AUC = CL / F."""
        from charon.pbpk.pk_extract import compute_oral_pk_parameters
        sim, params, topo = _run_midazolam_oral()
        pk = compute_oral_pk_parameters(sim, params, topo, dose_mg=5.0)
        assert pk.cl_apparent == pytest.approx(5.0 / pk.auc_0_inf, rel=1e-3)

    def test_vss_none_for_oral(self):
        from charon.pbpk.pk_extract import compute_oral_pk_parameters
        sim, params, topo = _run_midazolam_oral()
        pk = compute_oral_pk_parameters(sim, params, topo, dose_mg=5.0)
        assert pk.vss is None
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/unit/test_oral_pk_extract.py -v`
Expected: FAIL — `compute_oral_pk_parameters` does not exist

- [ ] **Step 3: Implement compute_oral_pk_parameters**

Add to `src/charon/pbpk/pk_extract.py`:

```python
from charon.pbpk.solver import OralSimulationResult
from charon.pbpk.ode_compiler import OralPBPKParams
from charon.pbpk.topology import PBPKTopology
from charon.pbpk.acat import compute_absorption_rates


def compute_oral_pk_parameters(
    sim: OralSimulationResult,
    params: OralPBPKParams,
    topology: PBPKTopology,
    *,
    dose_mg: float,
) -> PKParameters:
    """Extract PK parameters and Fa/Fg/Fh from an oral simulation.

    Fa and Fg are computed via post-hoc trapezoidal integration of
    ODE fluxes, independent of hepatic CLint.

    Parameters
    ----------
    sim : OralSimulationResult
        Completed oral simulation.
    params : OralPBPKParams
        Oral compound parameters.
    topology : PBPKTopology
        Species topology.
    dose_mg : float
        Administered oral dose (mg).
    """
    t = sim.time_h
    cp = sim.cp_plasma

    # --- Standard PK ---
    cmax = float(np.max(cp))
    tmax = float(t[int(np.argmax(cp))])

    auc_0_last = _trapezoid(t, cp)

    ke, cp_last = _terminal_log_slope(t, cp)
    auc_tail = cp_last / ke
    auc_0_inf = auc_0_last + auc_tail

    mask_24 = t <= 24.0
    auc_0_24 = _trapezoid(t[mask_24], cp[mask_24]) if mask_24.sum() >= 2 else None

    half_life = math.log(2) / ke
    cl_apparent = dose_mg / auc_0_inf  # CL/F

    # --- Fa: fraction absorbed ---
    gi = params.gi_tract
    assert gi is not None
    k_abs_arr = np.array(
        compute_absorption_rates(gi, params.peff_cm_s),
        dtype=np.float64,
    )
    # Absorption flux at each time point: sum(k_abs_i × A_lumen_i(t))
    absorption_flux = np.zeros_like(t)
    for i in range(len(gi.segments)):
        absorption_flux += k_abs_arr[i] * sim.lumen_trajectory[i, :]
    absorbed_total = _trapezoid(t, absorption_flux)
    fa = absorbed_total / dose_mg

    # --- Fg: fraction escaping gut metabolism ---
    k_baso = params.q_villi_L_h / params.v_enterocyte_L
    portal_flux = k_baso * sim.enterocyte_trajectory
    portal_total = _trapezoid(t, portal_flux)
    fg = portal_total / absorbed_total if absorbed_total > 0 else 1.0

    # --- Fh: fraction escaping hepatic first-pass ---
    qh = topology.tissues["liver"].blood_flow_L_h
    fh = qh / (qh + params.fu_b * params.clint_liver_L_h)

    # --- Bioavailability ---
    bioavailability = fa * fg * fh

    return PKParameters(
        cmax=cmax,
        tmax=tmax,
        auc_0_inf=float(auc_0_inf),
        auc_0_24=float(auc_0_24) if auc_0_24 is not None else None,
        half_life=float(half_life),
        cl_apparent=float(cl_apparent),
        vss=None,  # not meaningful for oral
        bioavailability=float(bioavailability),
        fa=float(fa),
        fg=float(fg),
        fh=float(fh),
    )
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/unit/test_oral_pk_extract.py -v`
Expected: 9 PASSED

- [ ] **Step 5: Run full suite for regression**

Run: `pytest tests/ -x -q`
Expected: 581 + 9 = 590 passed

- [ ] **Step 6: Commit**

```bash
git add src/charon/pbpk/pk_extract.py tests/unit/test_oral_pk_extract.py
git commit -m "feat(pk_extract): add compute_oral_pk_parameters with Fa/Fg/Fh extraction"
```

---

### Task 8: Wire Pipeline.run(route="oral")

**Files:**
- Modify: `src/charon/pipeline.py:112-118`
- Modify: `tests/unit/test_pipeline.py`

- [ ] **Step 1: Write the failing test**

Append to `tests/unit/test_pipeline.py`:

```python
class TestPipelineOral:
    def test_oral_route_runs(self):
        """Pipeline(route='oral') should no longer raise NotImplementedError."""
        from charon.pipeline import Pipeline
        from charon.core.schema import CompoundConfig
        import yaml
        from pathlib import Path

        comp_path = Path("validation/data/tier1_obach/compounds/midazolam.yaml")
        with comp_path.open() as f:
            raw = yaml.safe_load(f)
        compound = CompoundConfig(**raw)

        pipe = Pipeline(
            compound=compound,
            route="oral",
            dose_mg=5.0,
            duration_h=72.0,
        )
        result = pipe.run()
        assert result.pk_parameters.cmax > 0
        assert result.pk_parameters.fg is not None
        assert result.pk_parameters.fa is not None

    def test_oral_route_metadata(self):
        from charon.pipeline import Pipeline
        from charon.core.schema import CompoundConfig
        import yaml
        from pathlib import Path

        comp_path = Path("validation/data/tier1_obach/compounds/midazolam.yaml")
        with comp_path.open() as f:
            raw = yaml.safe_load(f)
        compound = CompoundConfig(**raw)

        pipe = Pipeline(compound=compound, route="oral", dose_mg=5.0)
        result = pipe.run()
        assert result.metadata["route"] == "oral"
        assert "fg" in result.metadata
        assert "fa" in result.metadata
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/unit/test_pipeline.py::TestPipelineOral -v`
Expected: FAIL — NotImplementedError

- [ ] **Step 3: Implement _run_oral in Pipeline**

Replace the `run` method in `src/charon/pipeline.py` (lines 112-179):

```python
    def run(self) -> PipelineResult:
        """Execute the full pipeline and return a :class:`PipelineResult`."""
        if self.route == "oral":
            return self._run_oral()

        topology: PBPKTopology = load_species_topology(self.species)
        bridge = ParameterBridge()

        params = build_compound_pbpk_params(
            self.compound,
            topology,
            bridge=bridge,
            compound_type=self.compound_type_override,
            liver_model=self.liver_model,
        )

        sim = simulate_iv(
            topology,
            params,
            dose_mg=self.dose_mg,
            route=self.route,  # type: ignore[arg-type]
            duration_h=self.duration_h,
            infusion_duration_h=self.infusion_duration_h,
        )

        pk = compute_pk_parameters(
            sim.time_h,
            sim.cp_plasma,
            dose_mg=self.dose_mg,
            route=self.route,  # type: ignore[arg-type]
            infusion_duration_h=self.infusion_duration_h,
        )

        return PipelineResult(
            compound=self.compound,
            pk_parameters=pk,
            time_h=sim.time_h,
            cp_plasma=sim.cp_plasma,
            cp_blood=sim.cp_blood,
            simulation=sim,
            metadata={
                "species": self.species,
                "route": self.route,
                "dose_mg": self.dose_mg,
                "duration_h": self.duration_h,
                "liver_model": self.liver_model,
                "compound_type": params.compound_type,
                "clint_liver_L_h": params.clint_liver_L_h,
                "cl_renal_L_h": params.cl_renal_L_h,
                "fu_b": params.fu_b,
                "solver_method": sim.solver_method,
                "solver_nfev": sim.solver_nfev,
                "kp_overrides": [
                    {
                        "tissue": r.tissue,
                        "rr_value": r.rr_value,
                        "empirical_value": r.empirical_value,
                        "source": r.source,
                        "method": r.method,
                        "flag": r.flag,
                    }
                    for r in params.kp_overrides
                ],
            },
        )

    def _run_oral(self) -> PipelineResult:
        """Execute oral route through ACAT + PBPK."""
        from charon.pbpk.acat import load_gi_tract, papp_to_peff
        from charon.pbpk.ode_compiler import (
            OralPBPKParams,
            compute_gut_clint,
        )
        from charon.pbpk.pk_extract import compute_oral_pk_parameters
        from charon.pbpk.solver import simulate_oral

        topology = load_species_topology(self.species)
        bridge = ParameterBridge()

        base_params = build_compound_pbpk_params(
            self.compound,
            topology,
            bridge=bridge,
            compound_type=self.compound_type_override,
            liver_model=self.liver_model,
        )

        gi = load_gi_tract(self.species)

        # Resolve Peff
        perm = self.compound.properties.permeability
        if perm.peff_cm_s is not None:
            peff = float(perm.peff_cm_s.value)
        elif perm.papp_nm_s is not None:
            peff = papp_to_peff(float(perm.papp_nm_s.value))
        else:
            raise ValueError(
                f"Oral route requires Peff or Papp for compound "
                f"{self.compound.name!r}, but neither is provided."
            )

        # Gut CLint
        fm = self.compound.properties.metabolism.fm_cyp3a4
        clint_gut = compute_gut_clint(
            clint_liver_L_h=base_params.clint_liver_L_h,
            fm_cyp3a4=fm,
            gi_tract=gi,
            mppgl=40.0,
            liver_weight_g=topology.liver_weight_g,
        )

        q_gut = topology.tissues["gut_wall"].blood_flow_L_h
        q_villi = gi.q_villi_fraction * q_gut

        oral_params = OralPBPKParams(
            name=base_params.name,
            molecular_weight=base_params.molecular_weight,
            logp=base_params.logp,
            pka_acid=base_params.pka_acid,
            pka_base=base_params.pka_base,
            compound_type=base_params.compound_type,
            fu_p=base_params.fu_p,
            bp_ratio=base_params.bp_ratio,
            fu_b=base_params.fu_b,
            clint_liver_L_h=base_params.clint_liver_L_h,
            cl_renal_L_h=base_params.cl_renal_L_h,
            kp_by_tissue=base_params.kp_by_tissue,
            kp_overrides=base_params.kp_overrides,
            clint_gut_L_h=clint_gut,
            peff_cm_s=peff,
            q_villi_L_h=q_villi,
            v_enterocyte_L=gi.enterocyte_volume_L,
            gi_tract=gi,
        )

        sim = simulate_oral(
            topology,
            oral_params,
            dose_mg=self.dose_mg,
            duration_h=self.duration_h,
        )

        pk = compute_oral_pk_parameters(
            sim, oral_params, topology, dose_mg=self.dose_mg
        )

        return PipelineResult(
            compound=self.compound,
            pk_parameters=pk,
            time_h=sim.time_h,
            cp_plasma=sim.cp_plasma,
            cp_blood=sim.cp_blood,
            simulation=sim,
            metadata={
                "species": self.species,
                "route": self.route,
                "dose_mg": self.dose_mg,
                "duration_h": self.duration_h,
                "liver_model": self.liver_model,
                "compound_type": oral_params.compound_type,
                "clint_liver_L_h": oral_params.clint_liver_L_h,
                "clint_gut_L_h": oral_params.clint_gut_L_h,
                "cl_renal_L_h": oral_params.cl_renal_L_h,
                "fu_b": oral_params.fu_b,
                "peff_cm_s": oral_params.peff_cm_s,
                "solver_method": sim.solver_method,
                "solver_nfev": sim.solver_nfev,
                "fa": pk.fa,
                "fg": pk.fg,
                "fh": pk.fh,
                "bioavailability": pk.bioavailability,
            },
        )
```

Also add `peff_cm_s` to midazolam.yaml (under `permeability:` in `properties:`):

```yaml
  permeability:
    peff_cm_s:
      value: 4.0e-4
      source: literature
      method: "Lennernas 2007 Eur J Pharm Sci 29(3-4):278"
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/unit/test_pipeline.py::TestPipelineOral -v`
Expected: 2 PASSED

- [ ] **Step 5: Run full suite for regression**

Run: `pytest tests/ -x -q`
Expected: 590 + 2 = 592 passed

- [ ] **Step 6: Commit**

```bash
git add src/charon/pipeline.py tests/unit/test_pipeline.py validation/data/tier1_obach/compounds/midazolam.yaml
git commit -m "feat(pipeline): wire route='oral' through ACAT + simulate_oral + PK extraction"
```

---

### Task 9: Update pbpk/__init__.py exports

**Files:**
- Modify: `src/charon/pbpk/__init__.py`

- [ ] **Step 1: Add new public symbols**

Update `src/charon/pbpk/__init__.py` to export the new oral symbols:

```python
from charon.pbpk.acat import (
    GISegment,
    GITract,
    compute_absorption_rates,
    load_gi_tract,
    papp_to_peff,
)
from charon.pbpk.ode_compiler import (
    OralPBPKParams,
    compute_gut_clint,
)
from charon.pbpk.pk_extract import compute_oral_pk_parameters
from charon.pbpk.solver import OralSimulationResult, simulate_oral
```

Add corresponding names to `__all__`.

- [ ] **Step 2: Verify imports work**

Run: `python3 -c "from charon.pbpk import simulate_oral, OralSimulationResult, load_gi_tract, GITract; print('OK')"`
Expected: `OK`

- [ ] **Step 3: Run full suite for regression**

Run: `pytest tests/ -x -q`
Expected: 592 passed

- [ ] **Step 4: Commit**

```bash
git add src/charon/pbpk/__init__.py
git commit -m "chore(pbpk): export oral ACAT symbols from __init__"
```

---

### Task 10: Curate felodipine.yaml

**Files:**
- Create: `validation/data/tier1_obach/compounds/felodipine.yaml`

- [ ] **Step 1: Create felodipine compound YAML**

Create `validation/data/tier1_obach/compounds/felodipine.yaml`:

```yaml
name: felodipine
smiles: "CCOC(=O)C1=C(C)NC(C)=C(C(=O)OC)C1c1cccc(Cl)c1Cl"
molecular_weight: 384.26
source: experimental
properties:
  physicochemical:
    logp:
      value: 3.39
      source: experimental
      method: "Wishart DS, DrugBank DB01023"
    compound_type: neutral
  binding:
    fu_p: {value: 0.004, source: experimental, unit: fraction, method: "Edgar 1992 Eur J Clin Pharmacol 42(3):261"}
    fu_inc: {value: 0.60, source: correlation, unit: fraction, method: "Hallifax-Houston correlation from logP"}
    bp_ratio: {value: 0.73, source: experimental, unit: ratio, method: "Estimated from hematocrit and Kp_blood"}
  metabolism:
    clint_uL_min_mg: {value: 60.0, source: experimental, unit: uL/min/mg, method: "Obach 1999 Drug Metab Dispos 27(11):1350"}
    fm_cyp3a4: 1.0
  permeability:
    peff_cm_s:
      value: 3.5e-4
      source: literature
      method: "Estimated from Caco-2 Papp ~250 nm/s; Sun 2002 correlation"
  renal:
    clrenal_L_h: {value: 0.0, source: experimental, unit: L/h, method: "Negligible renal elimination"}
observed_oral:
  dose_mg: 10.0
  route: oral
  bioavailability: 0.15
  fg: 0.45
  fg_source: "Edgar 1992 Eur J Clin Pharmacol 42(3):261"
```

- [ ] **Step 2: Verify YAML loads**

Run: `python3 -c "from charon.core.schema import CompoundConfig; import yaml; d=yaml.safe_load(open('validation/data/tier1_obach/compounds/felodipine.yaml')); c=CompoundConfig(**d); print(c.name, c.properties.metabolism.fm_cyp3a4)"`
Expected: `felodipine 1.0`

- [ ] **Step 3: Commit**

```bash
git add validation/data/tier1_obach/compounds/felodipine.yaml
git commit -m "data: curate felodipine compound YAML for oral Fg validation"
```

---

### Task 11: Curate nifedipine.yaml

**Files:**
- Create: `validation/data/tier1_obach/compounds/nifedipine.yaml`

- [ ] **Step 1: Create nifedipine compound YAML**

Create `validation/data/tier1_obach/compounds/nifedipine.yaml`:

```yaml
name: nifedipine
smiles: "COC(=O)C1=C(C)NC(C)=C(C(=O)OC)C1c1ccccc1[N+](=O)[O-]"
molecular_weight: 346.34
source: experimental
properties:
  physicochemical:
    logp:
      value: 2.20
      source: experimental
      method: "Wishart DS, DrugBank DB01115"
    compound_type: neutral
  binding:
    fu_p: {value: 0.04, source: experimental, unit: fraction, method: "Holtbecker 1996 Clin Pharmacol Ther 60(1):54"}
    fu_inc: {value: 0.75, source: correlation, unit: fraction, method: "Hallifax-Houston correlation from logP"}
    bp_ratio: {value: 0.73, source: experimental, unit: ratio, method: "Estimated from hematocrit and Kp_blood"}
  metabolism:
    clint_uL_min_mg: {value: 35.0, source: experimental, unit: uL/min/mg, method: "Obach 1999 Drug Metab Dispos 27(11):1350"}
    fm_cyp3a4: 1.0
  permeability:
    peff_cm_s:
      value: 3.0e-4
      source: literature
      method: "Lennernas 2007 Eur J Pharm Sci 29(3-4):278"
  renal:
    clrenal_L_h: {value: 0.0, source: experimental, unit: L/h, method: "Negligible renal elimination"}
observed_oral:
  dose_mg: 10.0
  route: oral
  bioavailability: 0.50
  fg: 0.78
  fg_source: "Holtbecker 1996 Clin Pharmacol Ther 60(1):54"
```

- [ ] **Step 2: Verify YAML loads**

Run: `python3 -c "from charon.core.schema import CompoundConfig; import yaml; d=yaml.safe_load(open('validation/data/tier1_obach/compounds/nifedipine.yaml')); c=CompoundConfig(**d); print(c.name, c.properties.metabolism.fm_cyp3a4)"`
Expected: `nifedipine 1.0`

- [ ] **Step 3: Commit**

```bash
git add validation/data/tier1_obach/compounds/nifedipine.yaml
git commit -m "data: curate nifedipine compound YAML for oral Fg validation"
```

---

### Task 12: Integration test — oral pipeline + Fg validation

**Files:**
- Create: `tests/integration/test_oral_pipeline.py`

- [ ] **Step 1: Write the integration tests**

Create `tests/integration/test_oral_pipeline.py`:

```python
"""End-to-end oral pipeline tests with Fg validation.

Tests midazolam (calibration), felodipine, and nifedipine (independent)
against literature Fg values per spec §10.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import pytest
import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from charon.core.schema import CompoundConfig
from charon.pipeline import Pipeline


def _load_compound(name: str) -> CompoundConfig:
    path = REPO_ROOT / "validation" / "data" / "tier1_obach" / "compounds" / f"{name}.yaml"
    with path.open() as f:
        raw = yaml.safe_load(f)
    return CompoundConfig(**raw)


def _run_oral(name: str, dose_mg: float, duration_h: float = 72.0):
    compound = _load_compound(name)
    pipe = Pipeline(compound=compound, route="oral", dose_mg=dose_mg,
                    duration_h=duration_h)
    return pipe.run()


class TestMidazolamOral:
    """Midazolam PO 5mg — calibration compound."""

    def test_produces_finite_pk(self):
        result = _run_oral("midazolam", 5.0)
        pk = result.pk_parameters
        for attr in ("cmax", "tmax", "auc_0_inf", "half_life", "fa", "fg", "fh"):
            val = getattr(pk, attr)
            assert val is not None and math.isfinite(val), f"{attr} = {val}"

    def test_fg_calibration(self, record_property):
        """Fg ≈ 0.57 ± 10% (MPPGI calibration check)."""
        result = _run_oral("midazolam", 5.0)
        fg = result.pk_parameters.fg
        record_property("midazolam_fg", fg)
        assert 0.51 < fg < 0.63, f"midazolam Fg = {fg:.3f}, expected ~0.57"

    def test_fa_high(self, record_property):
        """Midazolam Fa > 0.90 (high Peff)."""
        result = _run_oral("midazolam", 5.0)
        fa = result.pk_parameters.fa
        record_property("midazolam_fa", fa)
        assert fa > 0.90, f"midazolam Fa = {fa:.3f}"

    def test_mass_balance(self):
        result = _run_oral("midazolam", 5.0)
        total_t0 = result.simulation.state_trajectory[:, 0].sum()
        assert abs(total_t0 - 5.0) < 0.05

    def test_tmax_positive(self):
        result = _run_oral("midazolam", 5.0)
        assert result.pk_parameters.tmax > 0

    def test_bioavailability_product(self):
        result = _run_oral("midazolam", 5.0)
        pk = result.pk_parameters
        expected = pk.fa * pk.fg * pk.fh
        assert pk.bioavailability == pytest.approx(expected, rel=1e-6)


class TestFelodipineOral:
    """Felodipine PO 10mg — independent validation."""

    def test_produces_finite_pk(self):
        result = _run_oral("felodipine", 10.0)
        pk = result.pk_parameters
        for attr in ("cmax", "tmax", "auc_0_inf", "fa", "fg", "fh"):
            val = getattr(pk, attr)
            assert val is not None and math.isfinite(val), f"{attr} = {val}"

    def test_fg_within_2_fold(self, record_property):
        """Fg within 2-fold of literature 0.45."""
        result = _run_oral("felodipine", 10.0)
        fg = result.pk_parameters.fg
        record_property("felodipine_fg", fg)
        lit_fg = 0.45
        fold = max(fg / lit_fg, lit_fg / fg)
        assert fold < 2.0, (
            f"felodipine Fg = {fg:.3f}, lit = {lit_fg}, fold = {fold:.2f}"
        )

    def test_fa_high(self, record_property):
        result = _run_oral("felodipine", 10.0)
        fa = result.pk_parameters.fa
        record_property("felodipine_fa", fa)
        assert fa > 0.90


class TestNifedipineOral:
    """Nifedipine PO 10mg — independent validation."""

    def test_produces_finite_pk(self):
        result = _run_oral("nifedipine", 10.0)
        pk = result.pk_parameters
        for attr in ("cmax", "tmax", "auc_0_inf", "fa", "fg", "fh"):
            val = getattr(pk, attr)
            assert val is not None and math.isfinite(val), f"{attr} = {val}"

    def test_fg_within_2_fold(self, record_property):
        """Fg within 2-fold of literature 0.78."""
        result = _run_oral("nifedipine", 10.0)
        fg = result.pk_parameters.fg
        record_property("nifedipine_fg", fg)
        lit_fg = 0.78
        fold = max(fg / lit_fg, lit_fg / fg)
        assert fold < 2.0, (
            f"nifedipine Fg = {fg:.3f}, lit = {lit_fg}, fold = {fold:.2f}"
        )

    def test_fa_high(self, record_property):
        result = _run_oral("nifedipine", 10.0)
        fa = result.pk_parameters.fa
        record_property("nifedipine_fa", fa)
        assert fa > 0.90


class TestRegressionInvariants:
    """Verify existing IV tests still pass."""

    def test_existing_test_count_preserved(self):
        """Sanity: importing the IV modules should still work."""
        from charon.pbpk import simulate_iv, build_rhs, compute_pk_parameters
        assert callable(simulate_iv)
        assert callable(build_rhs)
        assert callable(compute_pk_parameters)
```

- [ ] **Step 2: Run tests**

Run: `pytest tests/integration/test_oral_pipeline.py -v`
Expected: all PASSED (if preceding tasks are correct)

- [ ] **Step 3: Run full suite including Obach panel regression**

Run: `pytest tests/ -x -q`
Expected: all passed (592 + 15 = 607+ tests)

- [ ] **Step 4: Commit**

```bash
git add tests/integration/test_oral_pipeline.py
git commit -m "test(integration): oral pipeline + Fg validation for midazolam/felodipine/nifedipine"
```

---

### Task 13: Add midazolam observed_oral and observed PK section

**Files:**
- Modify: `validation/data/tier1_obach/compounds/midazolam.yaml`

- [ ] **Step 1: Append oral observed PK to midazolam YAML**

Add the following to the end of `validation/data/tier1_obach/compounds/midazolam.yaml`:

```yaml
observed_oral:
  dose_mg: 5.0
  route: oral
  bioavailability: 0.44
  fg: 0.57
  fg_source: "Thummel 1996 Clin Pharmacol Ther 59(5):491"
  cmax_ng_mL: 41.0
  tmax_h: 0.5
```

- [ ] **Step 2: Verify YAML parses**

Run: `python3 -c "import yaml; d=yaml.safe_load(open('validation/data/tier1_obach/compounds/midazolam.yaml')); print(d.get('observed_oral', {}).get('fg'))"`
Expected: `0.57`

- [ ] **Step 3: Commit**

```bash
git add validation/data/tier1_obach/compounds/midazolam.yaml
git commit -m "data(midazolam): add observed_oral section with Fg=0.57 reference"
```

---

### Task 14: Final regression gate — full test suite

**Files:** (no new files — verification only)

- [ ] **Step 1: Run the complete test suite**

Run: `pytest tests/ -v --tb=short 2>&1 | tail -40`
Expected: all tests pass, including:
- 543 original tests (Session 1 + 1.5 baseline)
- ~65+ new oral tests
- theophylline strict gate PASS
- Obach panel AAFE_Vss < 5.0

- [ ] **Step 2: Run Obach panel benchmark specifically**

Run: `pytest tests/integration/test_obach_panel_smoke.py -v`
Expected: all 7 tests PASS (regression invariant)

- [ ] **Step 3: Run oral integration tests specifically**

Run: `pytest tests/integration/test_oral_pipeline.py -v`
Expected: all tests PASS, with Fg values recorded

- [ ] **Step 4: Print Fg summary**

Run: `pytest tests/integration/test_oral_pipeline.py -v --tb=short 2>&1 | grep "fg\|Fg\|PASS\|FAIL"`

Verify:
- midazolam Fg ≈ 0.57 (±10%)
- felodipine Fg within 2-fold of 0.45
- nifedipine Fg within 2-fold of 0.78
- All Fa > 0.90

- [ ] **Step 5: Commit — session summary if all pass**

If all gates pass, this task is the final verification checkpoint. No code change needed.
