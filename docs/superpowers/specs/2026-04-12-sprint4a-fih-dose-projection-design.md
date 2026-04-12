# Sprint 4a — FIH Dose Projection (HED + MABEL + PAD)

**Date:** 2026-04-12
**Scope:** Implement three FIH dose projection methods (HED from NOAEL,
MABEL from target Kd, PAD from target Ceff) using human PK from Layer 2.
No cross-species PBPK or allometric scaling.
**Prior state:** Session 2a complete (612 tests, oral + IV pipeline, Fg
validated). Layer 2 PK output provides CL, Vss, AUC, F, half_life.
**Expected scale:** ~12 tasks.

---

## 1. Goal

Add FIH (First-in-Human) dose recommendation to Charon's pipeline:

1. HED from user-provided NOAEL via FDA BSA scaling.
2. MABEL from user-provided target Kd + predicted human PK.
3. PAD from user-provided target Ceff + predicted human PK.
4. MRSD = min(available methods) — always most conservative.
5. Pipeline integration: `Pipeline.run()` with `dose_projection` config
   produces `FIHDoseRecommendation` alongside PK results.

**Non-goals for Sprint 4a:**

- Cross-species PBPK simulation (Sprint 4b)
- Allometric scaling / Rule of Exponents (Sprint 4b — needs animal PK data)
- Consensus scaling (Sprint 4b)
- Uncertainty quantification on dose (Sprint 5)
- Population variability (Phase B)

---

## 2. Architecture

```
Layer 2 PK → DoseProjector → FIHDoseRecommendation
                  ↑
         DoseProjectionConfig
         (NOAEL, species, Kd, Ceff, safety_factor)
```

DoseProjector runs each method independently, collects results, and
selects MRSD = min(all available). Methods with missing inputs are
silently skipped. At least one method must be runnable.

---

## 3. HED — Human Equivalent Dose (hed.py)

### 3.1 Formula

```
HED [mg/kg] = NOAEL [mg/kg] × (Km_animal / Km_human)
MRSD [mg] = HED × body_weight_kg / safety_factor
```

### 3.2 Km values (FDA Guidance for Industry, 2005)

| Species | Km (kg/m²) |
|---------|-----------|
| human | 37.0 |
| rat | 6.2 |
| mouse | 3.0 |
| dog | 20.0 |
| monkey | 12.0 |
| rabbit | 12.0 |
| guinea_pig | 8.0 |

Stored as a module-level dict in `hed.py`. Species name matching is
case-insensitive.

### 3.3 Salt form correction

If the compound has a salt form (`CompoundConfig.salt_form` with
`salt_factor != 1.0`), the dose is corrected:

```
mrsd_free_base = mrsd_salt × salt_factor
```

where `salt_factor = MW_free / MW_salt < 1.0`.

### 3.4 Interface

```python
@dataclass(frozen=True)
class HEDResult:
    noael_mg_kg: float
    noael_species: str
    km_animal: float
    km_human: float
    hed_mg_kg: float
    body_weight_kg: float
    safety_factor: float
    mrsd_mg: float        # before salt correction
```

```python
def compute_hed(
    *,
    noael_mg_kg: float,
    noael_species: str,
    safety_factor: float = 10.0,
    body_weight_kg: float = 70.0,
) -> HEDResult:
```

### 3.5 Validation

Hand calculation: rat NOAEL 50 mg/kg, SF=10, BW=70 kg.
- HED = 50 × (6.2 / 37) = 8.378 mg/kg
- MRSD = 8.378 × 70 / 10 = 58.65 mg

---

## 4. MABEL — Minimum Anticipated Biological Effect Level (mabel.py)

### 4.1 PK-driven approach (no PD model)

Two sub-methods, both reported, MRSD = min:

**Single-dose (Cmax-based, conservative for FIH SAD):**
```
Css_u_target [mg/L] = target_kd_nM × MW × 1e-6
Css_total [mg/L] = Css_u_target / fu_p
dose_cmax_mg = Css_total × Vd_apparent_L
mrsd_cmax = dose_cmax_mg / safety_factor
```

Where Vd_apparent:
- IV route: Vss (from PKParameters)
- Oral route: half_life × CL_apparent / ln(2)  (= Vd/F)

Since CL_apparent and Vd_apparent already incorporate F for oral,
no separate /F is needed. This keeps the formula route-agnostic.

**Steady-state (Css-based, for repeat dosing context):**
```
dose_ss_mg = Css_total × CL_apparent_L_h × tau_h
mrsd_ss = dose_ss_mg / safety_factor
```

CL_apparent = CL/F for oral (already accounts for bioavailability),
= CL for IV (F=1 by definition).

MABEL MRSD = min(mrsd_cmax, mrsd_ss).

### 4.2 Unit conversion

target_kd_nM → mg/L:
```
kd_mg_L = kd_nM × molecular_weight_g_mol × 1e-6
         = kd_nM × MW / 1e6
```

### 4.3 Interface

```python
@dataclass(frozen=True)
class MABELResult:
    target_kd_nM: float
    target_conc_mg_L: float   # Css_total (total, not unbound)
    fu_p: float
    molecular_weight: float
    dose_cmax_mg: float       # single-dose estimate
    dose_ss_mg: float         # steady-state estimate
    safety_factor: float
    tau_h: float
    mrsd_mg: float            # min(cmax, ss) / SF
    limiting_approach: str    # "cmax" or "steady_state"
```

```python
def compute_mabel(
    *,
    target_kd_nM: float,
    molecular_weight: float,
    fu_p: float,
    cl_apparent_L_h: float,   # CL/F (oral) or CL (IV)
    vd_apparent_L: float,     # Vd/F (oral) or Vss (IV)
    safety_factor: float = 10.0,
    tau_h: float = 24.0,
) -> MABELResult:
```

### 4.4 PK input resolution

All dose formulas use apparent PK parameters (already F-corrected for
oral, or raw values for IV where F=1). No separate F division needed.

| Route | cl_apparent | vd_apparent | Interpretation |
|-------|-------------|-------------|---------------|
| oral | pk.cl_apparent (= CL/F) | pk.half_life × pk.cl_apparent / ln(2) (= Vd/F) | Dose is for oral administration |
| IV | pk.cl_apparent (= CL) | pk.vss | Dose is for IV administration |

The DoseProjector resolves these from PKParameters based on the route
used in the simulation. The `route` is documented in the recommendation
so the user knows what the dose refers to.

### 4.5 Validation

midazolam: Kd=10 nM, MW=325.77, fu_p=0.03.
- Css_u = 10 × 325.77 / 1e6 = 0.003258 mg/L
- Css_total = 0.003258 / 0.03 = 0.1086 mg/L
- With CL/F and Vd from oral PK → compute dose → verify manually

---

## 5. PAD — Pharmacologically Active Dose (pad.py)

### 5.1 Formula

```
Ceff_mg_L = target_ceff_nM × MW / 1e6
AUC_target = Ceff_mg_L × tau_h
dose_mg = AUC_target × CL_apparent_L_h
mrsd_pad = dose_mg / safety_factor
```

This assumes maintaining average concentration ≥ Ceff over the dosing
interval. More conservative than peak-based.

### 5.2 Interface

```python
@dataclass(frozen=True)
class PADResult:
    target_ceff_nM: float
    target_conc_mg_L: float
    auc_target_mg_h_L: float
    cl_L_h: float
    bioavailability: float
    dose_mg: float
    safety_factor: float
    tau_h: float
    mrsd_mg: float
```

```python
def compute_pad(
    *,
    target_ceff_nM: float,
    molecular_weight: float,
    cl_apparent_L_h: float,   # CL/F (oral) or CL (IV)
    safety_factor: float = 10.0,
    tau_h: float = 24.0,
) -> PADResult:
```

---

## 6. DoseProjector — Coordinator

### 6.1 FIHDoseRecommendation

```python
class FIHDoseRecommendation(BaseModel):
    model_config = ConfigDict(frozen=True)

    mrsd_mg: float
    limiting_method: str          # "hed", "mabel", or "pad"
    hed: HEDResult | None = None
    mabel: MABELResult | None = None
    pad: PADResult | None = None
    safety_factor: float
    salt_factor: float = 1.0
    rationale: str                # human-readable summary
```

### 6.2 Logic

```python
def project_fih_dose(
    *,
    pk: PKParameters,
    compound: CompoundConfig,
    config: DoseProjectionConfig,
    body_weight_kg: float = 70.0,
) -> FIHDoseRecommendation:
```

1. If `config.noael_mg_kg` and `config.noael_species` → run HED.
2. If `config.target_kd_nM` → run MABEL using PK.
3. If `config.target_ceff_nM` → run PAD using PK.
4. If no method runnable → raise ValueError.
5. Collect results, apply salt_factor to each MRSD.
6. MRSD = min(available MRSD values).
7. Build rationale string explaining which methods ran and why the
   minimum was chosen.

### 6.3 Rationale generation

Example output:
```
FIH dose recommendation: 5.86 mg (free base)
  HED (rat NOAEL 50 mg/kg, SF=10): 58.65 mg
  MABEL (Kd=10 nM, SF=10): 5.86 mg ← limiting
  PAD: not computed (target_ceff_nM not provided)
Most conservative method: MABEL
Salt correction: none (salt_factor=1.0)
```

---

## 7. Schema Changes

### 7.1 DoseProjectionConfig — already exists

```python
class DoseProjectionConfig(BaseModel):
    noael_mg_kg: float | None = None
    noael_species: str | None = None
    safety_factor: float = 10.0
    target_kd_nM: float | None = None
    target_ceff_nM: float | None = None
```

Already in schema.py. No changes needed.

### 7.2 Add tau_h and body_weight_kg to DoseProjectionConfig

```python
class DoseProjectionConfig(BaseModel):
    noael_mg_kg: float | None = None
    noael_species: str | None = None
    safety_factor: float = 10.0
    target_kd_nM: float | None = None
    target_ceff_nM: float | None = None
    tau_h: float = 24.0               # NEW: dosing interval
    body_weight_kg: float = 70.0      # NEW: subject body weight
```

### 7.3 PipelineResult — add dose_recommendation

Add to `PipelineResult`:
```python
dose_recommendation: FIHDoseRecommendation | None = None
```

### 7.4 Pipeline — add dose_projection parameter

```python
class Pipeline:
    def __init__(self, ..., dose_projection: DoseProjectionConfig | None = None):
```

---

## 8. File Changes

### New files

| File | Purpose |
|------|---------|
| `src/charon/translational/hed.py` | HED computation (BSA scaling) |
| `src/charon/translational/mabel.py` | MABEL computation (PK-driven) |
| `src/charon/translational/pad.py` | PAD computation |
| `src/charon/translational/dose_projector.py` | Coordinator: runs methods, selects MRSD |
| `tests/unit/test_hed.py` | HED unit tests |
| `tests/unit/test_mabel.py` | MABEL unit tests |
| `tests/unit/test_pad.py` | PAD unit tests |
| `tests/unit/test_dose_projector.py` | Coordinator tests |
| `tests/integration/test_fih_pipeline.py` | E2E: Pipeline with dose_projection |

### Modified files

| File | Changes |
|------|---------|
| `src/charon/core/schema.py` | Add tau_h, body_weight_kg to DoseProjectionConfig |
| `src/charon/pipeline.py` | Add dose_projection param, call DoseProjector after PK |
| `src/charon/translational/__init__.py` | Export public symbols |

---

## 9. Validation

### 9.1 Hand calculations

**HED test case:**
- Input: rat NOAEL=50 mg/kg, SF=10, BW=70 kg
- HED = 50 × 6.2/37 = 8.378 mg/kg
- MRSD = 8.378 × 70 / 10 = 58.65 mg

**MABEL test case:**
- Input: Kd=10 nM, MW=325.77, fu_p=0.03, CL/F=18.6 L/h, half_life=1.7 h
- Css_u = 10 × 325.77 / 1e6 = 3.258e-3 mg/L
- Css_total = 3.258e-3 / 0.03 = 0.1086 mg/L
- Vd_apparent = 1.7 × 18.6 / ln(2) = 45.6 L
- dose_cmax = 0.1086 × 45.6 = 4.95 mg (before SF)
- dose_ss = 0.1086 × 18.6 × 24 = 48.48 mg (before SF)
- MABEL MRSD = min(4.95, 48.48) / 10 = 0.495 mg

**PAD test case (oral route):**
- Input: Ceff=100 nM, MW=325.77, CL_apparent=18.6 L/h (= CL/F)
- Ceff_mg_L = 100 × 325.77 / 1e6 = 0.03258 mg/L
- AUC_target = 0.03258 × 24 = 0.7819 mg·h/L
- dose = 0.7819 × 18.6 = 14.54 mg (CL_apparent already incorporates F)
- MRSD = 14.54 / 10 = 1.454 mg

### 9.2 Acceptance criteria

| Criterion | Target |
|-----------|--------|
| HED matches hand calculation to 6 digits | Numerical correctness |
| MABEL matches hand calculation | Correct unit conversion nM→mg/L |
| PAD matches hand calculation | Correct AUC→dose conversion |
| MRSD = min(available) | Most conservative selection |
| Missing inputs → method skipped | Graceful degradation |
| All 612 existing tests pass | Regression |
| Pipeline(dose_projection=...) returns recommendation | E2E |
| Salt correction applied when present | Safety |
| Unknown species raises ValueError | Input validation |

---

## 10. CLAUDE.md Safety Rules Compliance

| Rule | How addressed |
|------|--------------|
| §6i: MRSD = min(NOAEL, MABEL, PAD) | Explicit min() in coordinator |
| §6i: salt_factor correction | Applied to all MRSD values |
| §6i: all three methods reported | FIHDoseRecommendation has all results |
| §0b-5: most conservative value | min() enforced, never max() |
| §0b-6: salt form correction | dose_free = dose_salt × salt_factor |

---

## 11. Definition of Done

- [ ] `compute_hed()` matches hand calculation for rat, dog, monkey.
- [ ] `compute_mabel()` matches hand calculation with correct nM→mg/L.
- [ ] `compute_pad()` matches hand calculation.
- [ ] `project_fih_dose()` selects min(available) as MRSD.
- [ ] Missing inputs → method gracefully skipped.
- [ ] Salt correction applied when salt_form present.
- [ ] Unknown species → ValueError.
- [ ] Pipeline(dose_projection=...) returns FIHDoseRecommendation.
- [ ] midazolam E2E: oral PK + dose projection → realistic recommendation.
- [ ] All 612 existing tests pass (regression).
- [ ] New tests: ≥25 new tests.
