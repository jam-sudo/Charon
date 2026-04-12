# Sprint 4a: FIH Dose Projection Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add three FIH dose projection methods (HED, MABEL, PAD) that take human PK from Layer 2 + user inputs and produce a conservative dose recommendation (MRSD = min of all methods).

**Architecture:** Each method is an independent pure function in its own file under `translational/`. A `DoseProjector` coordinator runs available methods, selects the minimum, and returns `FIHDoseRecommendation`. Pipeline gains an optional `dose_projection` parameter.

**Tech Stack:** Python 3.11+, Pydantic v2, pytest, dataclasses

**Spec:** `docs/superpowers/specs/2026-04-12-sprint4a-fih-dose-projection-design.md`

---

## File Structure

| Action | File | Responsibility |
|--------|------|----------------|
| Modify | `src/charon/core/schema.py:314-321` | Add `tau_h`, `body_weight_kg` to DoseProjectionConfig |
| Create | `src/charon/translational/hed.py` | HED computation (BSA scaling, FDA Km) |
| Create | `src/charon/translational/mabel.py` | MABEL computation (PK-driven, Cmax + Css) |
| Create | `src/charon/translational/pad.py` | PAD computation (efficacy-driven) |
| Create | `src/charon/translational/dose_projector.py` | Coordinator: runs methods, MRSD selection, rationale |
| Modify | `src/charon/pipeline.py:53-65` | Add `dose_projection` param, call DoseProjector |
| Modify | `src/charon/pipeline.py:38-47` | Add `dose_recommendation` to PipelineResult |
| Modify | `src/charon/translational/__init__.py` | Export public symbols |
| Create | `tests/unit/test_hed.py` | HED unit tests |
| Create | `tests/unit/test_mabel.py` | MABEL unit tests |
| Create | `tests/unit/test_pad.py` | PAD unit tests |
| Create | `tests/unit/test_dose_projector.py` | Coordinator unit tests |
| Create | `tests/integration/test_fih_pipeline.py` | E2E Pipeline + dose projection |

---

### Task 1: Extend DoseProjectionConfig in schema.py

**Files:**
- Modify: `src/charon/core/schema.py:314-321`
- Test: `tests/unit/test_schema.py`

- [ ] **Step 1: Write the failing test**

Append to `tests/unit/test_schema.py`:

```python
class TestDoseProjectionConfig:
    def test_default_tau_h(self):
        from charon.core.schema import DoseProjectionConfig
        d = DoseProjectionConfig()
        assert d.tau_h == 24.0

    def test_default_body_weight(self):
        from charon.core.schema import DoseProjectionConfig
        d = DoseProjectionConfig()
        assert d.body_weight_kg == 70.0

    def test_custom_values(self):
        from charon.core.schema import DoseProjectionConfig
        d = DoseProjectionConfig(
            noael_mg_kg=50.0,
            noael_species="rat",
            target_kd_nM=10.0,
            tau_h=12.0,
            body_weight_kg=60.0,
        )
        assert d.tau_h == 12.0
        assert d.body_weight_kg == 60.0
        assert d.noael_mg_kg == 50.0
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/unit/test_schema.py::TestDoseProjectionConfig -v`
Expected: FAIL — `tau_h` does not exist

- [ ] **Step 3: Update DoseProjectionConfig**

In `src/charon/core/schema.py`, replace `DoseProjectionConfig` (lines 314-321):

```python
class DoseProjectionConfig(BaseModel):
    """First-in-human dose-projection inputs."""

    noael_mg_kg: float | None = None
    noael_species: str | None = None
    safety_factor: float = 10.0
    target_kd_nM: float | None = None
    target_ceff_nM: float | None = None
    tau_h: float = 24.0
    body_weight_kg: float = 70.0
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/unit/test_schema.py::TestDoseProjectionConfig -v`
Expected: 3 PASSED

- [ ] **Step 5: Run full suite for regression**

Run: `pytest tests/ -x -q`
Expected: 612 passed

- [ ] **Step 6: Commit**

```bash
git add src/charon/core/schema.py tests/unit/test_schema.py
git commit -m "feat(schema): add tau_h and body_weight_kg to DoseProjectionConfig"
```

---

### Task 2: Implement HED (hed.py)

**Files:**
- Create: `src/charon/translational/hed.py`
- Create: `tests/unit/test_hed.py`

- [ ] **Step 1: Write the failing tests**

Create `tests/unit/test_hed.py`:

```python
"""Tests for HED (Human Equivalent Dose) computation."""
from __future__ import annotations

import pytest


class TestComputeHED:
    def test_rat_noael_50(self):
        """Hand calc: HED=50×6.2/37=8.378, MRSD=8.378×70/10=58.649."""
        from charon.translational.hed import compute_hed
        result = compute_hed(
            noael_mg_kg=50.0,
            noael_species="rat",
            safety_factor=10.0,
            body_weight_kg=70.0,
        )
        assert result.hed_mg_kg == pytest.approx(50.0 * 6.2 / 37.0, rel=1e-6)
        assert result.mrsd_mg == pytest.approx(50.0 * 6.2 / 37.0 * 70.0 / 10.0, rel=1e-6)

    def test_dog_noael_10(self):
        """Hand calc: HED=10×20/37=5.405, MRSD=5.405×70/10=37.838."""
        from charon.translational.hed import compute_hed
        result = compute_hed(noael_mg_kg=10.0, noael_species="dog")
        assert result.hed_mg_kg == pytest.approx(10.0 * 20.0 / 37.0, rel=1e-6)
        assert result.mrsd_mg == pytest.approx(10.0 * 20.0 / 37.0 * 70.0 / 10.0, rel=1e-6)

    def test_monkey_noael_25(self):
        """Monkey Km=12."""
        from charon.translational.hed import compute_hed
        result = compute_hed(noael_mg_kg=25.0, noael_species="monkey")
        assert result.hed_mg_kg == pytest.approx(25.0 * 12.0 / 37.0, rel=1e-6)

    def test_mouse(self):
        """Mouse Km=3."""
        from charon.translational.hed import compute_hed
        result = compute_hed(noael_mg_kg=100.0, noael_species="mouse")
        assert result.km_animal == 3.0

    def test_case_insensitive_species(self):
        from charon.translational.hed import compute_hed
        r1 = compute_hed(noael_mg_kg=50.0, noael_species="Rat")
        r2 = compute_hed(noael_mg_kg=50.0, noael_species="RAT")
        assert r1.mrsd_mg == r2.mrsd_mg

    def test_unknown_species_raises(self):
        from charon.translational.hed import compute_hed
        with pytest.raises(ValueError, match="species"):
            compute_hed(noael_mg_kg=50.0, noael_species="hamster")

    def test_custom_safety_factor(self):
        from charon.translational.hed import compute_hed
        r3 = compute_hed(noael_mg_kg=50.0, noael_species="rat", safety_factor=3.0)
        r10 = compute_hed(noael_mg_kg=50.0, noael_species="rat", safety_factor=10.0)
        assert r3.mrsd_mg == pytest.approx(r10.mrsd_mg * 10.0 / 3.0, rel=1e-6)

    def test_custom_body_weight(self):
        from charon.translational.hed import compute_hed
        r = compute_hed(noael_mg_kg=50.0, noael_species="rat", body_weight_kg=50.0)
        assert r.body_weight_kg == 50.0
        assert r.mrsd_mg == pytest.approx(50.0 * 6.2 / 37.0 * 50.0 / 10.0, rel=1e-6)

    def test_result_fields(self):
        from charon.translational.hed import compute_hed, HEDResult
        r = compute_hed(noael_mg_kg=50.0, noael_species="rat")
        assert isinstance(r, HEDResult)
        assert r.noael_mg_kg == 50.0
        assert r.noael_species == "rat"
        assert r.km_animal == 6.2
        assert r.km_human == 37.0
        assert r.safety_factor == 10.0

    def test_negative_noael_raises(self):
        from charon.translational.hed import compute_hed
        with pytest.raises(ValueError):
            compute_hed(noael_mg_kg=-1.0, noael_species="rat")
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/unit/test_hed.py -v`
Expected: FAIL — module does not exist

- [ ] **Step 3: Implement hed.py**

Create `src/charon/translational/hed.py`:

```python
"""HED (Human Equivalent Dose) computation via FDA BSA scaling.

Converts a preclinical NOAEL to a human equivalent dose using the
body surface area (BSA) correction factor Km, per FDA Guidance for
Industry (2005) "Estimating the Maximum Safe Starting Dose in Initial
Clinical Trials for Therapeutics in Adult Healthy Volunteers."

    HED [mg/kg] = NOAEL [mg/kg] × (Km_animal / Km_human)
    MRSD [mg] = HED × body_weight_kg / safety_factor
"""

from __future__ import annotations

from dataclasses import dataclass

# FDA 2005 Guidance Km values (kg/m²) — body weight / BSA ratio.
KM_BY_SPECIES: dict[str, float] = {
    "human": 37.0,
    "rat": 6.2,
    "mouse": 3.0,
    "dog": 20.0,
    "monkey": 12.0,
    "rabbit": 12.0,
    "guinea_pig": 8.0,
}


@dataclass(frozen=True)
class HEDResult:
    """Result of a single HED dose projection."""

    noael_mg_kg: float
    noael_species: str
    km_animal: float
    km_human: float
    hed_mg_kg: float
    body_weight_kg: float
    safety_factor: float
    mrsd_mg: float


def compute_hed(
    *,
    noael_mg_kg: float,
    noael_species: str,
    safety_factor: float = 10.0,
    body_weight_kg: float = 70.0,
) -> HEDResult:
    """Compute Human Equivalent Dose from a preclinical NOAEL.

    Parameters
    ----------
    noael_mg_kg : float
        No-Observed-Adverse-Effect Level in mg/kg (preclinical).
    noael_species : str
        Species name (case-insensitive). Must be in KM_BY_SPECIES.
    safety_factor : float
        Applied to HED to get MRSD. Default 10 (FDA standard).
    body_weight_kg : float
        Human body weight for mg/kg → mg conversion. Default 70 kg.

    Returns
    -------
    HEDResult

    Raises
    ------
    ValueError
        If species is unknown or noael_mg_kg <= 0.
    """
    if noael_mg_kg <= 0:
        raise ValueError(f"noael_mg_kg must be > 0, got {noael_mg_kg}")

    species_key = noael_species.lower()
    if species_key not in KM_BY_SPECIES:
        raise ValueError(
            f"Unknown species {noael_species!r}. "
            f"Supported: {sorted(KM_BY_SPECIES.keys())}"
        )

    km_animal = KM_BY_SPECIES[species_key]
    km_human = KM_BY_SPECIES["human"]

    hed_mg_kg = noael_mg_kg * (km_animal / km_human)
    mrsd_mg = hed_mg_kg * body_weight_kg / safety_factor

    return HEDResult(
        noael_mg_kg=noael_mg_kg,
        noael_species=species_key,
        km_animal=km_animal,
        km_human=km_human,
        hed_mg_kg=hed_mg_kg,
        body_weight_kg=body_weight_kg,
        safety_factor=safety_factor,
        mrsd_mg=mrsd_mg,
    )
```

- [ ] **Step 4: Run tests**

Run: `pytest tests/unit/test_hed.py -v`
Expected: 10 PASSED

- [ ] **Step 5: Regression**

Run: `pytest tests/ -x -q`
Expected: 612 + 10 = 622 passed

- [ ] **Step 6: Commit**

```bash
git add src/charon/translational/hed.py tests/unit/test_hed.py
git commit -m "feat(translational): implement HED computation with FDA BSA scaling"
```

---

### Task 3: Implement MABEL (mabel.py)

**Files:**
- Create: `src/charon/translational/mabel.py`
- Create: `tests/unit/test_mabel.py`

- [ ] **Step 1: Write the failing tests**

Create `tests/unit/test_mabel.py`:

```python
"""Tests for MABEL (Minimum Anticipated Biological Effect Level)."""
from __future__ import annotations

import math
import pytest


class TestComputeMABEL:
    def test_midazolam_hand_calc(self):
        """
        Kd=10 nM, MW=325.77, fu_p=0.03.
        Css_u = 10 × 325.77 / 1e6 = 3.258e-3 mg/L
        Css_total = 3.258e-3 / 0.03 = 0.1086 mg/L
        CL_apparent=18.6 L/h, half_life=1.7 h
        Vd_apparent = 1.7 × 18.6 / ln(2) = 45.61 L
        dose_cmax = 0.1086 × 45.61 = 4.953 mg
        dose_ss = 0.1086 × 18.6 × 24 = 48.48 mg
        MRSD = min(4.953, 48.48) / 10 = 0.4953 mg
        """
        from charon.translational.mabel import compute_mabel
        result = compute_mabel(
            target_kd_nM=10.0,
            molecular_weight=325.77,
            fu_p=0.03,
            cl_apparent_L_h=18.6,
            vd_apparent_L=1.7 * 18.6 / math.log(2),
            safety_factor=10.0,
            tau_h=24.0,
        )
        css_u = 10.0 * 325.77 / 1e6
        css_total = css_u / 0.03
        vd = 1.7 * 18.6 / math.log(2)
        dose_cmax = css_total * vd
        dose_ss = css_total * 18.6 * 24.0
        expected_mrsd = min(dose_cmax, dose_ss) / 10.0

        assert result.mrsd_mg == pytest.approx(expected_mrsd, rel=1e-4)
        assert result.limiting_approach == "cmax"

    def test_target_conc_conversion(self):
        """Verify nM → mg/L conversion."""
        from charon.translational.mabel import compute_mabel
        result = compute_mabel(
            target_kd_nM=100.0,
            molecular_weight=400.0,
            fu_p=0.5,
            cl_apparent_L_h=10.0,
            vd_apparent_L=100.0,
        )
        # Css_u = 100 × 400 / 1e6 = 0.04 mg/L
        # Css_total = 0.04 / 0.5 = 0.08 mg/L
        assert result.target_conc_mg_L == pytest.approx(0.08, rel=1e-6)

    def test_ss_limiting_with_short_tau(self):
        """With short tau, SS dose can be less than Cmax dose."""
        from charon.translational.mabel import compute_mabel
        result = compute_mabel(
            target_kd_nM=10.0,
            molecular_weight=300.0,
            fu_p=0.1,
            cl_apparent_L_h=5.0,
            vd_apparent_L=500.0,  # large Vd → large Cmax dose
            tau_h=4.0,            # short interval → small SS dose
        )
        # dose_cmax = Css × 500 = large
        # dose_ss = Css × 5 × 4 = 20×Css → smaller than Css×500
        assert result.limiting_approach == "steady_state"

    def test_result_fields(self):
        from charon.translational.mabel import compute_mabel, MABELResult
        r = compute_mabel(
            target_kd_nM=10.0,
            molecular_weight=300.0,
            fu_p=0.5,
            cl_apparent_L_h=10.0,
            vd_apparent_L=50.0,
        )
        assert isinstance(r, MABELResult)
        assert r.target_kd_nM == 10.0
        assert r.safety_factor == 10.0
        assert r.tau_h == 24.0
        assert r.dose_cmax_mg > 0
        assert r.dose_ss_mg > 0
        assert r.mrsd_mg > 0

    def test_negative_kd_raises(self):
        from charon.translational.mabel import compute_mabel
        with pytest.raises(ValueError):
            compute_mabel(
                target_kd_nM=-1.0,
                molecular_weight=300.0,
                fu_p=0.5,
                cl_apparent_L_h=10.0,
                vd_apparent_L=50.0,
            )

    def test_fu_p_zero_raises(self):
        from charon.translational.mabel import compute_mabel
        with pytest.raises(ValueError):
            compute_mabel(
                target_kd_nM=10.0,
                molecular_weight=300.0,
                fu_p=0.0,
                cl_apparent_L_h=10.0,
                vd_apparent_L=50.0,
            )
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/unit/test_mabel.py -v`
Expected: FAIL — module does not exist

- [ ] **Step 3: Implement mabel.py**

Create `src/charon/translational/mabel.py`:

```python
"""MABEL (Minimum Anticipated Biological Effect Level) dose projection.

PK-driven approach (no PD model). Two sub-estimates:
  1. Single-dose (Cmax-based): dose = Css_total × Vd_apparent
  2. Steady-state (Css-based): dose = Css_total × CL_apparent × tau

MRSD = min(single, steady) / safety_factor.

All PK parameters are "apparent" — CL/F for oral, CL for IV.
No separate bioavailability division needed.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class MABELResult:
    """Result of a MABEL dose projection."""

    target_kd_nM: float
    target_conc_mg_L: float
    fu_p: float
    molecular_weight: float
    dose_cmax_mg: float
    dose_ss_mg: float
    safety_factor: float
    tau_h: float
    mrsd_mg: float
    limiting_approach: str


def compute_mabel(
    *,
    target_kd_nM: float,
    molecular_weight: float,
    fu_p: float,
    cl_apparent_L_h: float,
    vd_apparent_L: float,
    safety_factor: float = 10.0,
    tau_h: float = 24.0,
) -> MABELResult:
    """Compute MABEL-based FIH dose from target Kd and human PK.

    Parameters
    ----------
    target_kd_nM : float
        Target binding affinity (Kd) in nM.
    molecular_weight : float
        Molecular weight in g/mol.
    fu_p : float
        Fraction unbound in plasma (0, 1].
    cl_apparent_L_h : float
        Apparent clearance: CL/F for oral, CL for IV.
    vd_apparent_L : float
        Apparent volume of distribution: Vd/F for oral, Vss for IV.
    safety_factor : float
        Default 10.
    tau_h : float
        Dosing interval in hours. Default 24.

    Returns
    -------
    MABELResult

    Raises
    ------
    ValueError
        If target_kd_nM <= 0 or fu_p <= 0.
    """
    if target_kd_nM <= 0:
        raise ValueError(f"target_kd_nM must be > 0, got {target_kd_nM}")
    if fu_p <= 0:
        raise ValueError(f"fu_p must be > 0 for MABEL (division by fu_p), got {fu_p}")

    # nM → mg/L: C [mg/L] = C [nM] × MW [g/mol] × 1e-6
    css_u_mg_L = target_kd_nM * molecular_weight * 1e-6
    css_total_mg_L = css_u_mg_L / fu_p

    dose_cmax = css_total_mg_L * vd_apparent_L
    dose_ss = css_total_mg_L * cl_apparent_L_h * tau_h

    if dose_cmax <= dose_ss:
        limiting = "cmax"
        mrsd = dose_cmax / safety_factor
    else:
        limiting = "steady_state"
        mrsd = dose_ss / safety_factor

    return MABELResult(
        target_kd_nM=target_kd_nM,
        target_conc_mg_L=css_total_mg_L,
        fu_p=fu_p,
        molecular_weight=molecular_weight,
        dose_cmax_mg=dose_cmax,
        dose_ss_mg=dose_ss,
        safety_factor=safety_factor,
        tau_h=tau_h,
        mrsd_mg=mrsd,
        limiting_approach=limiting,
    )
```

- [ ] **Step 4: Run tests**

Run: `pytest tests/unit/test_mabel.py -v`
Expected: 6 PASSED

- [ ] **Step 5: Regression**

Run: `pytest tests/ -x -q`
Expected: 622 + 6 = 628 passed

- [ ] **Step 6: Commit**

```bash
git add src/charon/translational/mabel.py tests/unit/test_mabel.py
git commit -m "feat(translational): implement MABEL dose projection (PK-driven, Cmax + Css)"
```

---

### Task 4: Implement PAD (pad.py)

**Files:**
- Create: `src/charon/translational/pad.py`
- Create: `tests/unit/test_pad.py`

- [ ] **Step 1: Write the failing tests**

Create `tests/unit/test_pad.py`:

```python
"""Tests for PAD (Pharmacologically Active Dose)."""
from __future__ import annotations

import pytest


class TestComputePAD:
    def test_hand_calc(self):
        """
        Ceff=100 nM, MW=325.77, CL_apparent=18.6 L/h, tau=24h.
        Ceff_mg_L = 100 × 325.77 / 1e6 = 0.032577 mg/L
        AUC_target = 0.032577 × 24 = 0.78185 mg·h/L
        dose = 0.78185 × 18.6 = 14.542 mg
        MRSD = 14.542 / 10 = 1.4542 mg
        """
        from charon.translational.pad import compute_pad
        result = compute_pad(
            target_ceff_nM=100.0,
            molecular_weight=325.77,
            cl_apparent_L_h=18.6,
            safety_factor=10.0,
            tau_h=24.0,
        )
        ceff_mg = 100.0 * 325.77 / 1e6
        auc = ceff_mg * 24.0
        dose = auc * 18.6
        expected = dose / 10.0
        assert result.mrsd_mg == pytest.approx(expected, rel=1e-4)

    def test_short_tau(self):
        from charon.translational.pad import compute_pad
        r24 = compute_pad(
            target_ceff_nM=100.0, molecular_weight=300.0,
            cl_apparent_L_h=10.0, tau_h=24.0,
        )
        r12 = compute_pad(
            target_ceff_nM=100.0, molecular_weight=300.0,
            cl_apparent_L_h=10.0, tau_h=12.0,
        )
        assert r12.mrsd_mg == pytest.approx(r24.mrsd_mg / 2.0, rel=1e-6)

    def test_result_fields(self):
        from charon.translational.pad import compute_pad, PADResult
        r = compute_pad(
            target_ceff_nM=50.0, molecular_weight=400.0,
            cl_apparent_L_h=15.0,
        )
        assert isinstance(r, PADResult)
        assert r.target_ceff_nM == 50.0
        assert r.dose_mg > 0
        assert r.mrsd_mg > 0
        assert r.safety_factor == 10.0

    def test_negative_ceff_raises(self):
        from charon.translational.pad import compute_pad
        with pytest.raises(ValueError):
            compute_pad(
                target_ceff_nM=-10.0, molecular_weight=300.0,
                cl_apparent_L_h=10.0,
            )

    def test_zero_cl_raises(self):
        from charon.translational.pad import compute_pad
        with pytest.raises(ValueError):
            compute_pad(
                target_ceff_nM=10.0, molecular_weight=300.0,
                cl_apparent_L_h=0.0,
            )
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/unit/test_pad.py -v`
Expected: FAIL — module does not exist

- [ ] **Step 3: Implement pad.py**

Create `src/charon/translational/pad.py`:

```python
"""PAD (Pharmacologically Active Dose) projection.

Efficacy-driven: compute the dose needed to maintain average plasma
concentration ≥ target Ceff over a dosing interval.

    Ceff [mg/L] = Ceff [nM] × MW / 1e6
    AUC_target = Ceff × tau
    dose = AUC_target × CL_apparent
    MRSD = dose / safety_factor

CL_apparent is CL/F for oral or CL for IV — already accounts for
bioavailability.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class PADResult:
    """Result of a PAD dose projection."""

    target_ceff_nM: float
    target_conc_mg_L: float
    auc_target_mg_h_L: float
    cl_apparent_L_h: float
    dose_mg: float
    safety_factor: float
    tau_h: float
    mrsd_mg: float


def compute_pad(
    *,
    target_ceff_nM: float,
    molecular_weight: float,
    cl_apparent_L_h: float,
    safety_factor: float = 10.0,
    tau_h: float = 24.0,
) -> PADResult:
    """Compute PAD-based FIH dose from target efficacious concentration.

    Parameters
    ----------
    target_ceff_nM : float
        Target efficacious total plasma concentration in nM.
    molecular_weight : float
        Molecular weight in g/mol.
    cl_apparent_L_h : float
        Apparent clearance: CL/F for oral, CL for IV.
    safety_factor : float
        Default 10.
    tau_h : float
        Dosing interval in hours. Default 24.

    Returns
    -------
    PADResult

    Raises
    ------
    ValueError
        If inputs are non-positive.
    """
    if target_ceff_nM <= 0:
        raise ValueError(f"target_ceff_nM must be > 0, got {target_ceff_nM}")
    if cl_apparent_L_h <= 0:
        raise ValueError(f"cl_apparent_L_h must be > 0, got {cl_apparent_L_h}")

    ceff_mg_L = target_ceff_nM * molecular_weight * 1e-6
    auc_target = ceff_mg_L * tau_h
    dose = auc_target * cl_apparent_L_h
    mrsd = dose / safety_factor

    return PADResult(
        target_ceff_nM=target_ceff_nM,
        target_conc_mg_L=ceff_mg_L,
        auc_target_mg_h_L=auc_target,
        cl_apparent_L_h=cl_apparent_L_h,
        dose_mg=dose,
        safety_factor=safety_factor,
        tau_h=tau_h,
        mrsd_mg=mrsd,
    )
```

- [ ] **Step 4: Run tests**

Run: `pytest tests/unit/test_pad.py -v`
Expected: 5 PASSED

- [ ] **Step 5: Regression**

Run: `pytest tests/ -x -q`
Expected: 628 + 5 = 633 passed

- [ ] **Step 6: Commit**

```bash
git add src/charon/translational/pad.py tests/unit/test_pad.py
git commit -m "feat(translational): implement PAD dose projection (efficacy-driven)"
```

---

### Task 5: Implement DoseProjector coordinator (dose_projector.py)

**Files:**
- Create: `src/charon/translational/dose_projector.py`
- Create: `tests/unit/test_dose_projector.py`

- [ ] **Step 1: Write the failing tests**

Create `tests/unit/test_dose_projector.py`:

```python
"""Tests for DoseProjector coordinator."""
from __future__ import annotations

import math
import pytest


def _make_pk():
    """Minimal PKParameters for testing."""
    from charon.core.schema import PKParameters
    return PKParameters(
        cmax=0.02,
        tmax=0.87,
        auc_0_inf=0.27,
        half_life=1.7,
        cl_apparent=18.6,
        vss=None,
        bioavailability=0.48,
        fa=0.97,
        fg=0.575,
        fh=0.86,
    )


def _make_compound():
    """Minimal CompoundConfig for testing."""
    from charon.core.schema import CompoundConfig
    return CompoundConfig(
        name="test_drug",
        smiles="C",
        molecular_weight=325.77,
        properties={
            "binding": {"fu_p": {"value": 0.03, "source": "experimental"}},
        },
    )


class TestProjectFIHDose:
    def test_hed_only(self):
        from charon.translational.dose_projector import project_fih_dose
        from charon.core.schema import DoseProjectionConfig
        config = DoseProjectionConfig(
            noael_mg_kg=50.0, noael_species="rat",
        )
        rec = project_fih_dose(
            pk=_make_pk(), compound=_make_compound(), config=config,
            route="oral",
        )
        assert rec.hed is not None
        assert rec.mabel is None
        assert rec.pad is None
        assert rec.mrsd_mg == rec.hed.mrsd_mg
        assert rec.limiting_method == "hed"

    def test_mabel_only(self):
        from charon.translational.dose_projector import project_fih_dose
        from charon.core.schema import DoseProjectionConfig
        config = DoseProjectionConfig(target_kd_nM=10.0)
        rec = project_fih_dose(
            pk=_make_pk(), compound=_make_compound(), config=config,
            route="oral",
        )
        assert rec.hed is None
        assert rec.mabel is not None
        assert rec.mrsd_mg == rec.mabel.mrsd_mg

    def test_pad_only(self):
        from charon.translational.dose_projector import project_fih_dose
        from charon.core.schema import DoseProjectionConfig
        config = DoseProjectionConfig(target_ceff_nM=100.0)
        rec = project_fih_dose(
            pk=_make_pk(), compound=_make_compound(), config=config,
            route="oral",
        )
        assert rec.hed is None
        assert rec.pad is not None
        assert rec.mrsd_mg == rec.pad.mrsd_mg

    def test_all_three_min_wins(self):
        from charon.translational.dose_projector import project_fih_dose
        from charon.core.schema import DoseProjectionConfig
        config = DoseProjectionConfig(
            noael_mg_kg=50.0, noael_species="rat",
            target_kd_nM=10.0,
            target_ceff_nM=100.0,
        )
        rec = project_fih_dose(
            pk=_make_pk(), compound=_make_compound(), config=config,
            route="oral",
        )
        assert rec.hed is not None
        assert rec.mabel is not None
        assert rec.pad is not None
        all_mrsd = [rec.hed.mrsd_mg, rec.mabel.mrsd_mg, rec.pad.mrsd_mg]
        assert rec.mrsd_mg == pytest.approx(min(all_mrsd), rel=1e-10)

    def test_no_inputs_raises(self):
        from charon.translational.dose_projector import project_fih_dose
        from charon.core.schema import DoseProjectionConfig
        config = DoseProjectionConfig()  # nothing provided
        with pytest.raises(ValueError, match="at least one"):
            project_fih_dose(
                pk=_make_pk(), compound=_make_compound(), config=config,
                route="oral",
            )

    def test_salt_correction(self):
        from charon.translational.dose_projector import project_fih_dose
        from charon.core.schema import DoseProjectionConfig, CompoundConfig, SaltForm
        config = DoseProjectionConfig(noael_mg_kg=50.0, noael_species="rat")
        compound = CompoundConfig(
            name="test_salt",
            smiles="C",
            molecular_weight=300.0,
            salt_form=SaltForm(name="HCl", mw_salt=336.0, salt_factor=300.0/336.0),
            properties={"binding": {"fu_p": {"value": 0.5, "source": "experimental"}}},
        )
        rec = project_fih_dose(
            pk=_make_pk(), compound=compound, config=config,
            route="oral",
        )
        assert rec.salt_factor < 1.0
        # MRSD should be salt-corrected (smaller than without salt)
        assert rec.salt_factor == pytest.approx(300.0 / 336.0, rel=1e-6)

    def test_rationale_contains_method_names(self):
        from charon.translational.dose_projector import project_fih_dose
        from charon.core.schema import DoseProjectionConfig
        config = DoseProjectionConfig(
            noael_mg_kg=50.0, noael_species="rat",
            target_kd_nM=10.0,
        )
        rec = project_fih_dose(
            pk=_make_pk(), compound=_make_compound(), config=config,
            route="oral",
        )
        assert "HED" in rec.rationale
        assert "MABEL" in rec.rationale

    def test_custom_safety_factor(self):
        from charon.translational.dose_projector import project_fih_dose
        from charon.core.schema import DoseProjectionConfig
        config3 = DoseProjectionConfig(noael_mg_kg=50.0, noael_species="rat", safety_factor=3.0)
        config10 = DoseProjectionConfig(noael_mg_kg=50.0, noael_species="rat", safety_factor=10.0)
        r3 = project_fih_dose(pk=_make_pk(), compound=_make_compound(), config=config3, route="oral")
        r10 = project_fih_dose(pk=_make_pk(), compound=_make_compound(), config=config10, route="oral")
        assert r3.mrsd_mg > r10.mrsd_mg
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/unit/test_dose_projector.py -v`
Expected: FAIL — module does not exist

- [ ] **Step 3: Implement dose_projector.py**

Create `src/charon/translational/dose_projector.py`:

```python
"""FIH Dose Projector — coordinator for HED, MABEL, and PAD methods.

Runs all applicable methods based on available inputs, selects the most
conservative (minimum) MRSD, and produces a structured recommendation.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Literal

from pydantic import BaseModel, ConfigDict

from charon.core.schema import (
    CompoundConfig,
    DoseProjectionConfig,
    PKParameters,
)
from charon.translational.hed import HEDResult, compute_hed
from charon.translational.mabel import MABELResult, compute_mabel
from charon.translational.pad import PADResult, compute_pad


class FIHDoseRecommendation(BaseModel):
    """Structured FIH dose recommendation."""

    model_config = ConfigDict(frozen=True)

    mrsd_mg: float
    limiting_method: str
    hed: HEDResult | None = None
    mabel: MABELResult | None = None
    pad: PADResult | None = None
    safety_factor: float
    salt_factor: float = 1.0
    route: str = "oral"
    rationale: str


def project_fih_dose(
    *,
    pk: PKParameters,
    compound: CompoundConfig,
    config: DoseProjectionConfig,
    route: str = "oral",
) -> FIHDoseRecommendation:
    """Run all applicable dose projection methods and select MRSD.

    Parameters
    ----------
    pk : PKParameters
        PK output from Layer 2 PBPK simulation.
    compound : CompoundConfig
        Compound data (for MW, fu_p, salt_form).
    config : DoseProjectionConfig
        User-provided dose projection inputs.
    route : str
        Route used in simulation ("oral", "iv_bolus", "iv_infusion").

    Returns
    -------
    FIHDoseRecommendation

    Raises
    ------
    ValueError
        If no method has sufficient inputs.
    """
    sf = config.safety_factor
    mw = compound.molecular_weight or 0.0
    salt_factor = 1.0
    if compound.salt_form is not None:
        salt_factor = compound.salt_form.salt_factor

    # Resolve apparent PK
    cl_apparent = pk.cl_apparent or 0.0
    if route == "oral" or pk.vss is None:
        half_life = pk.half_life or 1.0
        vd_apparent = half_life * cl_apparent / math.log(2) if cl_apparent > 0 else 0.0
    else:
        vd_apparent = pk.vss

    fu_p = None
    if compound.properties.binding.fu_p is not None:
        fu_p = compound.properties.binding.fu_p.value

    # --- Run methods ---
    hed_result: HEDResult | None = None
    mabel_result: MABELResult | None = None
    pad_result: PADResult | None = None
    candidates: list[tuple[str, float]] = []

    # HED
    if config.noael_mg_kg is not None and config.noael_species is not None:
        hed_result = compute_hed(
            noael_mg_kg=config.noael_mg_kg,
            noael_species=config.noael_species,
            safety_factor=sf,
            body_weight_kg=config.body_weight_kg,
        )
        candidates.append(("hed", hed_result.mrsd_mg * salt_factor))

    # MABEL
    if config.target_kd_nM is not None and fu_p is not None and fu_p > 0 and mw > 0:
        mabel_result = compute_mabel(
            target_kd_nM=config.target_kd_nM,
            molecular_weight=mw,
            fu_p=fu_p,
            cl_apparent_L_h=cl_apparent,
            vd_apparent_L=vd_apparent,
            safety_factor=sf,
            tau_h=config.tau_h,
        )
        candidates.append(("mabel", mabel_result.mrsd_mg * salt_factor))

    # PAD
    if config.target_ceff_nM is not None and mw > 0 and cl_apparent > 0:
        pad_result = compute_pad(
            target_ceff_nM=config.target_ceff_nM,
            molecular_weight=mw,
            cl_apparent_L_h=cl_apparent,
            safety_factor=sf,
            tau_h=config.tau_h,
        )
        candidates.append(("pad", pad_result.mrsd_mg * salt_factor))

    if not candidates:
        raise ValueError(
            "FIH dose projection requires at least one of: "
            "noael_mg_kg+noael_species, target_kd_nM, or target_ceff_nM"
        )

    # Select minimum
    limiting_method, mrsd_mg = min(candidates, key=lambda x: x[1])

    # Build rationale
    lines = [f"FIH dose recommendation: {mrsd_mg:.2f} mg"]
    for method_name, method_mrsd in candidates:
        label = method_name.upper()
        marker = " <- limiting" if method_name == limiting_method else ""
        lines.append(f"  {label}: {method_mrsd:.2f} mg (SF={sf}){marker}")
    for method_name in ("hed", "mabel", "pad"):
        if method_name not in [c[0] for c in candidates]:
            lines.append(f"  {method_name.upper()}: not computed (insufficient inputs)")
    if salt_factor != 1.0:
        lines.append(f"Salt correction applied: factor={salt_factor:.4f}")
    lines.append(f"Most conservative method: {limiting_method.upper()}")
    rationale = "\n".join(lines)

    return FIHDoseRecommendation(
        mrsd_mg=mrsd_mg,
        limiting_method=limiting_method,
        hed=hed_result,
        mabel=mabel_result,
        pad=pad_result,
        safety_factor=sf,
        salt_factor=salt_factor,
        route=route,
        rationale=rationale,
    )
```

- [ ] **Step 4: Run tests**

Run: `pytest tests/unit/test_dose_projector.py -v`
Expected: 8 PASSED

- [ ] **Step 5: Regression**

Run: `pytest tests/ -x -q`
Expected: 633 + 8 = 641 passed

- [ ] **Step 6: Commit**

```bash
git add src/charon/translational/dose_projector.py tests/unit/test_dose_projector.py
git commit -m "feat(translational): implement DoseProjector coordinator (MRSD = min of methods)"
```

---

### Task 6: Wire Pipeline with dose_projection

**Files:**
- Modify: `src/charon/pipeline.py`
- Create: `tests/integration/test_fih_pipeline.py`

- [ ] **Step 1: Write the failing tests**

Create `tests/integration/test_fih_pipeline.py`:

```python
"""End-to-end Pipeline + FIH dose projection tests."""
from __future__ import annotations

import math
import sys
from pathlib import Path

import pytest
import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from charon.core.schema import CompoundConfig, DoseProjectionConfig
from charon.pipeline import Pipeline


def _load_midazolam():
    path = REPO_ROOT / "validation" / "data" / "tier1_obach" / "compounds" / "midazolam.yaml"
    with path.open() as f:
        raw = yaml.safe_load(f)
    return CompoundConfig(**raw)


class TestPipelineWithDoseProjection:
    def test_oral_with_hed(self):
        compound = _load_midazolam()
        pipe = Pipeline(
            compound=compound,
            route="oral",
            dose_mg=5.0,
            dose_projection=DoseProjectionConfig(
                noael_mg_kg=50.0,
                noael_species="rat",
            ),
        )
        result = pipe.run()
        assert result.dose_recommendation is not None
        assert result.dose_recommendation.hed is not None
        assert result.dose_recommendation.mrsd_mg > 0
        # HED = 50 × 6.2/37 × 70 / 10 ≈ 58.65
        assert result.dose_recommendation.mrsd_mg == pytest.approx(58.65, rel=0.01)

    def test_oral_with_all_three(self):
        compound = _load_midazolam()
        pipe = Pipeline(
            compound=compound,
            route="oral",
            dose_mg=5.0,
            dose_projection=DoseProjectionConfig(
                noael_mg_kg=50.0,
                noael_species="rat",
                target_kd_nM=10.0,
                target_ceff_nM=100.0,
            ),
        )
        result = pipe.run()
        rec = result.dose_recommendation
        assert rec is not None
        assert rec.hed is not None
        assert rec.mabel is not None
        assert rec.pad is not None
        all_mrsd = [rec.hed.mrsd_mg, rec.mabel.mrsd_mg, rec.pad.mrsd_mg]
        # MRSD should be the minimum (most conservative)
        assert rec.mrsd_mg <= min(all_mrsd) + 0.01

    def test_no_dose_projection_returns_none(self):
        compound = _load_midazolam()
        pipe = Pipeline(compound=compound, route="oral", dose_mg=5.0)
        result = pipe.run()
        assert result.dose_recommendation is None

    def test_iv_with_hed(self):
        compound = _load_midazolam()
        pipe = Pipeline(
            compound=compound,
            route="iv_bolus",
            dose_mg=5.0,
            dose_projection=DoseProjectionConfig(
                noael_mg_kg=50.0,
                noael_species="rat",
            ),
        )
        result = pipe.run()
        assert result.dose_recommendation is not None
        assert result.dose_recommendation.hed is not None

    def test_rationale_is_string(self):
        compound = _load_midazolam()
        pipe = Pipeline(
            compound=compound,
            route="oral",
            dose_mg=5.0,
            dose_projection=DoseProjectionConfig(
                noael_mg_kg=50.0,
                noael_species="rat",
                target_kd_nM=10.0,
            ),
        )
        result = pipe.run()
        assert isinstance(result.dose_recommendation.rationale, str)
        assert "HED" in result.dose_recommendation.rationale
        assert "MABEL" in result.dose_recommendation.rationale

    def test_metadata_includes_dose(self):
        compound = _load_midazolam()
        pipe = Pipeline(
            compound=compound,
            route="oral",
            dose_mg=5.0,
            dose_projection=DoseProjectionConfig(
                noael_mg_kg=50.0,
                noael_species="rat",
            ),
        )
        result = pipe.run()
        assert "mrsd_mg" in result.metadata
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/integration/test_fih_pipeline.py -v`
Expected: FAIL — Pipeline doesn't accept `dose_projection`

- [ ] **Step 3: Modify Pipeline**

In `src/charon/pipeline.py`, make these changes:

**a) Add `dose_recommendation` to PipelineResult (line 38-47):**

```python
@dataclass
class PipelineResult:
    """Output of :meth:`Pipeline.run`."""

    compound: CompoundConfig
    pk_parameters: PKParameters
    time_h: np.ndarray
    cp_plasma: np.ndarray
    cp_blood: np.ndarray
    simulation: SimulationResult
    metadata: dict = field(default_factory=dict)
    dose_recommendation: "FIHDoseRecommendation | None" = None
```

Add this import at the top of pipeline.py (with other imports):
```python
from __future__ import annotations
```

**b) Add `dose_projection` to `Pipeline.__init__` (line 53-65):**

```python
    def __init__(
        self,
        compound: CompoundConfig,
        *,
        route: Literal["iv_bolus", "iv_infusion", "oral"],
        dose_mg: float,
        species: str = "human",
        duration_h: float = 72.0,
        infusion_duration_h: float = 0.0,
        liver_model: str = "well_stirred",
        compound_type_override: str | None = None,
        dose_projection: DoseProjectionConfig | None = None,
    ) -> None:
        self.compound = compound
        self.route = route
        self.dose_mg = dose_mg
        self.species = species
        self.duration_h = duration_h
        self.infusion_duration_h = infusion_duration_h
        self.liver_model = liver_model
        self.compound_type_override = compound_type_override
        self.dose_projection = dose_projection
```

Add `DoseProjectionConfig` to the imports from schema.

**c) Add `_maybe_project_dose` helper method:**

```python
    def _maybe_project_dose(self, pk: PKParameters) -> "FIHDoseRecommendation | None":
        """Run dose projection if configured."""
        if self.dose_projection is None:
            return None
        # Check if any method has inputs
        dp = self.dose_projection
        has_hed = dp.noael_mg_kg is not None and dp.noael_species is not None
        has_mabel = dp.target_kd_nM is not None
        has_pad = dp.target_ceff_nM is not None
        if not (has_hed or has_mabel or has_pad):
            return None

        from charon.translational.dose_projector import project_fih_dose
        return project_fih_dose(
            pk=pk,
            compound=self.compound,
            config=self.dose_projection,
            route=self.route,
        )
```

**d) Call `_maybe_project_dose` at the end of both `run()` and `_run_oral()`:**

In the IV `run()` method, before the final `return PipelineResult(...)`:

```python
        dose_rec = self._maybe_project_dose(pk)
```

Then add `dose_recommendation=dose_rec` to the PipelineResult constructor and `"mrsd_mg": dose_rec.mrsd_mg if dose_rec else None` to metadata.

Similarly in `_run_oral()`, after PK extraction:

```python
        dose_rec = self._maybe_project_dose(pk)
```

Add `dose_recommendation=dose_rec` and metadata field.

- [ ] **Step 4: Run tests**

Run: `pytest tests/integration/test_fih_pipeline.py -v`
Expected: 6 PASSED

- [ ] **Step 5: Run full regression**

Run: `pytest tests/ -x -q`
Expected: 641 + 6 = 647 passed

- [ ] **Step 6: Commit**

```bash
git add src/charon/pipeline.py tests/integration/test_fih_pipeline.py
git commit -m "feat(pipeline): integrate FIH dose projection (HED+MABEL+PAD) into Pipeline.run()"
```

---

### Task 7: Update translational/__init__.py exports

**Files:**
- Modify: `src/charon/translational/__init__.py`

- [ ] **Step 1: Write exports**

```python
"""Charon Layer 3 — translational scaling and FIH dose projection."""

from charon.translational.dose_projector import (
    FIHDoseRecommendation,
    project_fih_dose,
)
from charon.translational.hed import HEDResult, compute_hed
from charon.translational.mabel import MABELResult, compute_mabel
from charon.translational.pad import PADResult, compute_pad

__all__ = [
    "FIHDoseRecommendation",
    "HEDResult",
    "MABELResult",
    "PADResult",
    "compute_hed",
    "compute_mabel",
    "compute_pad",
    "project_fih_dose",
]
```

- [ ] **Step 2: Verify imports**

Run: `python3 -c "from charon.translational import project_fih_dose, compute_hed, compute_mabel, compute_pad; print('OK')"`
Expected: `OK`

- [ ] **Step 3: Regression**

Run: `pytest tests/ -x -q`
Expected: 647 passed

- [ ] **Step 4: Commit**

```bash
git add src/charon/translational/__init__.py
git commit -m "chore(translational): export FIH dose projection symbols from __init__"
```

---

### Task 8: Final regression gate

**Files:** (verification only)

- [ ] **Step 1: Run full test suite**

Run: `pytest tests/ -v --tb=short 2>&1 | tail -40`
Expected: all pass (612 baseline + ~35 new = 647+)

- [ ] **Step 2: Run Obach panel regression**

Run: `pytest tests/integration/test_obach_panel_smoke.py -v`
Expected: 7/7 PASS

- [ ] **Step 3: Run oral Fg validation**

Run: `pytest tests/integration/test_oral_pipeline.py -v`
Expected: 14/14 PASS

- [ ] **Step 4: Run FIH pipeline tests**

Run: `pytest tests/integration/test_fih_pipeline.py -v`
Expected: 6/6 PASS

- [ ] **Step 5: Print dose projection summary**

Run: `python3 -c "
from charon.pipeline import Pipeline
from charon.core.schema import CompoundConfig, DoseProjectionConfig
import yaml
d = yaml.safe_load(open('validation/data/tier1_obach/compounds/midazolam.yaml'))
c = CompoundConfig(**d)
r = Pipeline(c, route='oral', dose_mg=5.0, dose_projection=DoseProjectionConfig(
    noael_mg_kg=50.0, noael_species='rat', target_kd_nM=10.0, target_ceff_nM=100.0,
)).run()
print(r.dose_recommendation.rationale)
"`

Expected output showing all three methods with MRSD = min.
