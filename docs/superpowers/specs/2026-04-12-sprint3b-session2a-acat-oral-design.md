# Sprint 3b Session 2a — ACAT Core + Oral Route + Fg

**Date:** 2026-04-12
**Scope:** Add 8-segment GI lumen transit, enterocyte compartment with
CYP3A4 gut-wall metabolism, oral route in Pipeline, Fg extraction, and
validation against midazolam + felodipine + nifedipine PO.
**Prior state:** Session 1 + 1.5 complete (543 tests, IV-only kernel,
AAFE_Vss=3.12, theophylline strict gate PASS).
**Expected scale:** ~15 tasks.

---

## 1. Goal

Turn the Sprint 3a IV-only PBPK kernel into a functional oral dosing
engine that:

1. Simulates oral drug absorption through an 8-segment GI transit model.
2. Computes Fg (gut-wall first-pass extraction) mechanistically within
   the ODE via a dedicated enterocyte compartment.
3. Reports Fa, Fg, Fh, and F = Fa x Fg x Fh for every oral simulation.
4. Validates Fg within 2-fold of literature values for at least 2
   independent CYP3A4 substrates.

**Non-goals for Session 2a** (deferred to 2b/2c):

- Noyes-Whitney dissolution (2b)
- BCS classification (2b)
- Food effect / fed-state parameters (2b)
- P-gp efflux transport (future)
- Enterohepatic recirculation (future)
- Non-CYP3A4 gut metabolism (future)
- Per-segment enterocyte sub-compartments (future)
- Multiple dosing / steady-state oral (future)
- Precipitation / supersaturation (2c)
- Hepatocyte-sourced CLint for gut IVIVE (future)

---

## 2. Architecture Overview

### 2.1 State vector extension

Current IV kernel: `[A_venous, A_arterial, A_tissue_1..15]` = 17 states.

Session 2a adds:

| States | Count | Description |
|--------|-------|-------------|
| `A_lumen_i` (8 segments) | 8 | Dissolved drug in each GI segment |
| `A_enterocyte` | 1 | Drug in pooled enterocyte compartment |
| **Total new** | **9** | |
| **Grand total** | **26** | Trivial for BDF solver |

Fa and Fg are extracted post-hoc from the state trajectory via
trapezoidal integration of fluxes (no additional tracker states needed).

### 2.2 ODE architecture

Oral route adds three new dynamic subsystems to the existing PBPK ODE:

```
[1] GI Lumen Transit (8 segments, first-order chain)
         ↓ absorption (Peff-based)
[2] Enterocyte Compartment (CYP3A4 metabolism vs basolateral escape)
         ↓ basolateral transfer
[3] Portal Vein → Liver (existing PBPK tissue dynamics, unmodified)
```

The existing tissue dynamics (states 0–16) are **unchanged** for oral.
The only modification to the existing liver mass balance is an additive
portal inflow term from the enterocyte.

### 2.3 Code architecture

```
build_rhs()       → IV-only RHS (17 states), UNCHANGED
build_oral_rhs()  → Oral RHS (26 states), new function
                    tissue dynamics reimplemented inline (not shared helper)
                    to preserve IV regression safety via zero code coupling
```

Rationale for duplication: `ode_compiler.py` is currently 390 LOC.
Extracting a shared tissue-dynamics helper introduces coupling risk for
a function that is the single most safety-critical code path in Charon.
Refactoring into shared helpers is deferred to 2b+ when the oral path
is stable and tested.

---

## 3. GI Lumen Transit Model

### 3.1 Segment parameters

Ported from Sisyphus `reference_man.yaml` (lines 144–286), which in
turn derives from Yu & Amidon 1999 (Int J Pharm 186:119) and ACAT
literature.

| Segment | Volume (L) | Radius (cm) | ka_fraction | Transit rate (1/h) |
|---------|-----------|-------------|-------------|-------------------|
| stomach | 0.25 | 5.0 | 0.0 | 4.0 |
| duodenum | 0.05 | 1.6 | 1.0 | 3.846 |
| jejunum1 | 0.07 | 1.5 | 1.0 | 2.105 |
| jejunum2 | 0.07 | 1.5 | 1.0 | 2.105 |
| ileum1 | 0.06 | 1.3 | 0.8 | 1.471 |
| ileum2 | 0.06 | 1.3 | 0.6 | 1.471 |
| ileum3 | 0.06 | 1.3 | 0.3 | 1.471 |
| colon | 0.30 | 2.5 | 0.05 | 0.074 |

These are **compound-independent physiological constants**. The
ka_fraction values are empirical absorption-window scaling factors
accounting for villous surface area, tight junction changes, and
regional absorption capacity. They are NOT tuned to any specific
compound (CLAUDE.md §7.3 compliant).

Transit rate validation:
- Gastric emptying: k=4.0 → t_half=10.4 min (literature fasted: 10–20 min) ✓
- Total SITT: 1/3.846 + 2/2.105 + 3/1.471 = 0.26 + 0.95 + 2.04 = 3.25 h (literature: 3–4 h) ✓
- Colonic transit: k=0.074 → t_half=9.4 h (literature: 12–36 h mean ~20 h) ✓

### 3.2 Source term (oral dosing)

Immediate-release assumption (Session 2a): entire dose placed in
stomach lumen at t=0 as dissolved drug.

```
IC: A_stomach(0) = dose_mg, all other states = 0
```

### 3.3 Transit ODE

```
dA_stomach/dt = -k_transit_stomach × A_stomach

dA_segment_i/dt = k_transit_prev × A_prev
                  - k_transit_i × A_segment_i
                  - k_abs_i × A_segment_i

# Fecal excretion: colon transit empties into a virtual fecal sink.
# Mass leaving colon via transit = k_transit_colon × A_colon is
# permanently removed from the system (no downstream segment).
```

Stomach has `k_abs = 0` (ka_fraction=0), so it only empties to
duodenum — no gastric absorption.

### 3.4 Absorption rate constant

Mechanistic cylindrical model:

```
k_abs_i [1/h] = (2 × Peff [cm/s] × 3600 / R_i [cm]) × ka_fraction_i
```

Unit verification: `[cm/s] / [cm] × [s/h] × [-] = [1/h]` ✓

Peff source priority:
1. `PermeabilityProperties.peff_cm_s` (experimental or ML-predicted)
2. Fallback from Papp: `log10(Peff) = 0.4926 × log10(Papp_nm_s) - 0.1454`
   (Sun et al. 2002, simplified Caco-2→human correlation)
3. If neither available: raise ValueError (oral simulation requires Peff)

### 3.5 Fecal excretion

Unabsorbed drug reaching the end of the colon exits the system:
`fecal_rate = k_transit_colon × A_colon`. This mass is lost from the
system and contributes to `(1 - Fa)`.

---

## 4. Enterocyte Compartment

### 4.1 Physiology

The enterocyte is a pooled absorption compartment representing the
mucosal epithelial cell layer of the entire small intestine + colon.
It sits between the GI lumen and the portal blood supply.

| Parameter | Value | Source |
|-----------|-------|--------|
| V_enterocyte | 0.30 L | Estimated mucosal epithelial volume |
| Q_villi | 0.18 × Q_gut = 10.53 L/h | Gertz 2011 DMPK 26(5):486 |
| MPPGI | 6.0 mg/g | Literature range 2–20 (Yang 2007, Barter 2007); calibrated within range against midazolam Fg=0.57 (Thummel 1996) |
| Gut enterocyte weight | 400 g | Metabolically active mucosal mass |
| CYP3A4_gut (pmol/mg) | 31 | Paine 1997 Gastroenterology 113(2):453 |
| CYP3A4_liver (pmol/mg) | 137 | Shimada 1994 |

### 4.2 ODE

```
dA_enterocyte/dt = absorption_influx - gut_metabolism - basolateral_transfer

where:
  absorption_influx = sum(k_abs_i × A_lumen_i)        [mg/h]
  gut_metabolism    = (CLint_gut / V_enterocyte) × A_enterocyte  [mg/h]
  basolateral_transfer = (Q_villi / V_enterocyte) × A_enterocyte [mg/h]
```

At steady state, `Fg = k_baso / (k_baso + k_metab)` where:
- `k_baso = Q_villi / V_enterocyte` [1/h]
- `k_metab = CLint_gut / V_enterocyte` [1/h]

### 4.3 Portal inflow modification

The existing liver mass balance gains one additive term:

```python
# Existing (unchanged):
portal_inflow = sum(Q_portal_tissue × C_blood_out_portal_tissue)

# New addition:
portal_from_enterocyte = (Q_villi / V_enterocyte) × A_enterocyte

# Liver mass balance:
dy_liver = (Q_ha × C_art
            + portal_inflow
            + portal_from_enterocyte      # ← NEW
            - Q_liver × C_liver_blood_out
            - hepatic_elim)
```

This adds mass to the portal blood without additional blood flow. This
is physically correct — absorbed drug enters the portal circulation
(analogous to IV bolus adding mass to venous blood).

### 4.4 Relationship to gut_wall tissue

The gut_wall tissue (`topology.py` PORTAL_TISSUES) continues to exist
as a standard PBPK compartment for systemic drug distribution:

- **Enterocyte**: handles first-pass absorption and CYP3A4 metabolism.
  Drug flows: lumen → enterocyte → portal vein.
- **gut_wall tissue**: handles systemic recirculation distribution
  (arterial → gut_wall → portal via blood flow, with Kp partitioning).
  No CYP3A4 metabolism on recirculating drug.

V_enterocyte (0.30 L) is ADDITIONAL to V_gut_wall (1.03 L). Combined
= 1.33 L (within literature range for total gut wall 0.9–1.4 L).

**Known simplification**: systemically recirculating drug does NOT
encounter enterocyte CYP3A4. This is a standard ACAT simplification —
CYP3A4 metabolism applies only to first-pass absorption. The effect on
total clearance is second-order (gut CYP3A4 ≈ 0.45% of hepatic
CYP3A4 content by mass).

---

## 5. Gut CLint IVIVE (CYP3A4)

### 5.1 Formula — ratio approach (fu_inc-cancelling)

```python
CLint_gut_L_h = (CLint_liver_L_h
                 × fm_cyp3a4
                 × (CYP3A4_gut_per_mg / CYP3A4_liver_per_mg)
                 × (MPPGI × gut_enterocyte_weight_g)
                 / (MPPGL × liver_weight_g))
```

**Why this works (HLM only)**: `CLint_liver_L_h` from ParameterBridge
already incorporates fu_inc correction. Both liver and gut CLint use
microsomal-protein-based scaling (MPPGL and MPPGI respectively), so
fu_inc cancels in the numerator/denominator ratio. This cancellation
**only holds when both sides use the same scaling basis** (microsomal
protein). For hepatocyte-sourced CLint, the liver uses
hepatocellularity (cells/g) while the gut analog (HPGI) is unmeasured,
so the ratio approach does not apply — see §5.4.

### 5.2 fm_CYP3A4

New field in `MetabolismProperties`:

```python
fm_cyp3a4: float | None = None
```

- `None` (default) → CLint_gut = 0 → Fg = 1.0 (non-CYP3A4 substrates)
- `1.0` → entire hepatic CLint attributed to CYP3A4 (midazolam, felodipine, nifedipine)
- `0.0–1.0` → partial CYP3A4 contribution

### 5.3 Midazolam verification

```
CLint_liver_L_h ≈ 880 L/h (Charon ParameterBridge, midazolam HLM)
fm_CYP3A4 = 1.0
CYP3A4 ratio = 31/137 = 0.226
MPPGI × gut_weight = 6.0 × 400 = 2400
MPPGL × liver_weight = 40 × 1500 = 60000

CLint_gut = 880 × 1.0 × 0.226 × 2400/60000 = 880 × 0.00904 = 7.95 L/h

k_baso = 10.53/0.30 = 35.1 /h
k_metab = 7.95/0.30 = 26.5 /h
Fg = 35.1/(35.1 + 26.5) = 0.57  ✓ (Thummel 1996: 0.57)
```

### 5.4 Session 2a restriction

Only HLM-sourced CLint is supported for gut IVIVE. Hepatocyte-sourced
CLint uses cell-based scaling (hepatocellularity), which does not have
a gut analog (HPGI — hepatocytes per gram intestine — is unmeasured).

If `clint_system == "hepatocytes"` and `fm_cyp3a4 is not None`,
`build_oral_rhs` raises `NotImplementedError` with a message explaining
the limitation.

---

## 6. Fa / Fg / Fh / F Extraction

### 6.1 Fa — fraction absorbed

Post-hoc from state trajectory at simulation end:

```python
mass_remaining_lumen = sum(A_lumen_i(t_end))
mass_in_enterocyte = A_enterocyte(t_end)
# For long simulations, both should be ~0 for high-Peff drugs

# Cumulative fecal = dose - (absorbed + remaining_lumen + remaining_enterocyte + mass_in_PBPK)
# Simpler: integrate absorption flux
absorbed_total = trapz(sum(k_abs_i × A_lumen_i(t)), t)
Fa = absorbed_total / dose_mg
```

### 6.2 Fg — fraction escaping gut metabolism

Post-hoc trapezoidal integration of enterocyte fluxes:

```python
absorption_flux = sum(k_abs_i × A_lumen_i(t))           # mg/h at each t
portal_flux = (Q_villi / V_enterocyte) × A_enterocyte(t) # mg/h at each t

absorbed_total = trapz(absorption_flux, t)
portal_total = trapz(portal_flux, t)
Fg = portal_total / absorbed_total
```

This is **independent of CLint_liver** — only enterocyte dynamics are
involved. Fg validation is uncontaminated by the known CLint hepatic
overprediction (Sprint 3b Session 1 AAFE_CL = 4.89).

### 6.3 Fh — fraction escaping hepatic first-pass

Analytical from well-stirred steady-state (no additional tracker needed):

```python
Fh = Qh / (Qh + fu_b × CLint_liver_L_h)
```

**Note**: this will underestimate Fh for compounds with CLint
overprediction (midazolam: predicted Fh ≈ 0.12 vs literature ~0.72).
This is a known Sprint 3b Session 1 limitation, not a Session 2a issue.

### 6.4 Bioavailability

```python
F = Fa × Fg × Fh
```

For IV route: F = 1.0, Fa = 1.0, Fg = 1.0, Fh = None (unchanged from
current implementation).

---

## 7. Schema Changes

### 7.1 MetabolismProperties — add fm_cyp3a4

```python
class MetabolismProperties(BaseModel):
    primary_cyp: str | None = None
    secondary_cyp: str | None = None
    clint_uL_min_mg: PredictedProperty | None = None
    fm_cyp3a4: float | None = None          # NEW: fraction CYP3A4-metabolized
```

Validator: `0.0 <= fm_cyp3a4 <= 1.0` if not None.

### 7.2 PKParameters — already has oral fields

`PKParameters` already has `cmax`, `tmax`, `auc_0_inf`, `auc_0_24`,
`half_life`, `cl_apparent`, `vss`, `bioavailability`, `fa`, `fg`, `fh`.

For oral route:
- `cl_apparent` = CL/F (apparent oral clearance)
- `vss` = None (not meaningful for oral)
- `cmax`, `tmax` = from the oral Cp-time profile

### 7.3 PKParameters — no new fields needed

The existing schema covers all oral PK outputs.

---

## 8. human.yaml — GI Tract Section

New `gi_tract` section added to `src/charon/pbpk/species/human.yaml`:

```yaml
gi_tract:
  # Enterocyte physiology (pooled for all segments)
  enterocyte_volume_L: 0.30
  enterocyte_weight_g: 400
  q_villi_fraction: 0.18              # fraction of Q_gut → villous blood
                                      # Gertz 2011 DMPK 26(5):486
  mppgi_mg_g: 6.0                     # microsomal protein per g intestine
                                      # Literature range 2-20 (Yang 2007, Barter 2007)
                                      # Calibrated within range: midazolam Fg=0.57

  # CYP3A4 content for gut-liver IVIVE ratio
  cyp3a4_gut_pmol_per_mg: 31          # Paine 1997 Gastroenterology 113(2):453
  cyp3a4_liver_pmol_per_mg: 137       # Shimada 1994

  # Per-segment GI parameters (Yu & Amidon 1999, Sisyphus port)
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

---

## 9. File Changes

### New files

| File | Purpose |
|------|---------|
| `src/charon/pbpk/acat.py` | GISegment dataclass, GITract dataclass, load_gi_tract() from YAML, compute_absorption_rates(peff), papp_to_peff() fallback |
| `validation/data/tier1_obach/compounds/felodipine.yaml` | Felodipine compound data (CYP3A4 substrate, Fg validation) |
| `validation/data/tier1_obach/compounds/nifedipine.yaml` | Nifedipine compound data (CYP3A4 substrate, Fg validation) |

### Modified files

| File | Changes |
|------|---------|
| `src/charon/core/schema.py` | Add `fm_cyp3a4` to `MetabolismProperties` with validator |
| `src/charon/pbpk/species/human.yaml` | Add `gi_tract` section (§8) |
| `src/charon/pbpk/ode_compiler.py` | Add `OralPBPKParams` dataclass, `build_oral_rhs()` function |
| `src/charon/pbpk/solver.py` | Add `simulate_oral()` function, `OralSimulationResult` dataclass |
| `src/charon/pbpk/pk_extract.py` | Add `compute_oral_pk_parameters()` with Fa/Fg/Fh extraction |
| `src/charon/pipeline.py` | Wire `route="oral"` through acat + oral_rhs + simulate_oral |
| `src/charon/pbpk/topology.py` | Add optional GI data loading to `load_species_topology()` or separate loader |
| `validation/data/tier1_obach/compounds/midazolam.yaml` | Add `fm_cyp3a4: 1.0`, `peff_cm_s`, oral observed PK |

---

## 10. Validation

### 10.1 Compounds

| Compound | Dose | fm_CYP3A4 | Lit Fg | Lit F | Role | Fg source |
|----------|------|-----------|--------|-------|------|-----------|
| midazolam | PO 5 mg | 1.0 | 0.57 | 0.35–0.44 | Calibration (MPPGI=6.0 derived from this) | Thummel 1996 Clin Pharmacol Ther 59(5):491 |
| felodipine | PO 10 mg | 1.0 | 0.45 | 0.15–0.20 | Independent validation | Edgar 1992 Eur J Clin Pharmacol 42(3):261 |
| nifedipine | PO 10 mg | 1.0 | 0.78 | 0.45–0.68 | Independent validation | Holtbecker 1996 Clin Pharmacol Ther 60(1):54 |

midazolam is calibration (circular for Fg), so the independent test
has **2 compounds** (felodipine, nifedipine) spanning a Fg range of
0.45–0.78.

### 10.2 Acceptance criteria

| Criterion | Target | Rationale |
|-----------|--------|-----------|
| Fg within 2-fold for ≥2 independent compounds | felodipine + nifedipine | Core deliverable |
| Fa > 0.90 for all 3 compounds (all high-Peff) | Sanity check on absorption model | High-Peff drugs should be fully absorbed |
| midazolam Fg ≈ 0.57 (±10%) | Calibration check | Confirms MPPGI implementation |
| All existing 543 tests pass | Regression invariant | Zero IV-path breakage |
| theophylline strict gate PASS | Regression invariant | Sprint 3b Session 1 invariant |
| Mass balance residual < 1% of dose | ODE sanity | sum(states) + eliminated ≈ dose |

### 10.3 Cross-check with IV data

midazolam IV 5 mg is already in the Obach panel. Compare:
- IV CL vs oral CL/F → back-calculate F_observed = CL_iv / (CL/F)_oral
- Verify F_predicted ≈ F_observed (within the known CLint error bounds)

### 10.4 What we do NOT validate in 2a

- Overall F accuracy (dominated by Fh error from CLint overprediction)
- Dissolution-limited absorption (IR assumption only)
- Food effect
- BCS classification

---

## 11. Compound YAML Structure — Oral Extension

### 11.1 midazolam.yaml additions

```yaml
# Existing IV fields unchanged
properties:
  metabolism:
    fm_cyp3a4: 1.0
  permeability:
    peff_cm_s:
      value: 4.0e-4
      source: literature
      method: "Lennernas 2007 Eur J Pharm Sci 29(3-4):278"

# Observed oral PK (new section alongside existing observed_pk)
observed_oral:
  dose_mg: 5.0
  route: oral
  bioavailability: 0.44
  fg: 0.57
  fg_source: "Thummel 1996 Clin Pharmacol Ther 59(5):491"
  cmax_ng_mL: 41.0
  tmax_h: 0.5
```

### 11.2 felodipine.yaml (new)

All values from primary literature, §7.3 compliant. Key properties:
- MW: 384.3, logP: 3.4, pKa_base: None (neutral)
- fu_p: 0.004, BP: 0.73
- CLint HLM: ~60 uL/min/mg
- fm_CYP3A4: 1.0
- Peff: 3.5e-4 cm/s
- Observed Fg: 0.45 (Edgar 1992)

### 11.3 nifedipine.yaml (new)

Key properties:
- MW: 346.3, logP: 2.2, pKa_base: None (neutral)
- fu_p: 0.04, BP: 0.73
- CLint HLM: ~35 uL/min/mg
- fm_CYP3A4: 1.0
- Peff: 3.0e-4 cm/s
- Observed Fg: 0.78 (Holtbecker 1996)

---

## 12. Pipeline Wiring

### 12.1 Pipeline.run(route="oral") flow

```python
def run(self) -> PipelineResult:
    if self.route == "oral":
        return self._run_oral()
    # ... existing IV path unchanged

def _run_oral(self) -> PipelineResult:
    topology = load_species_topology(self.species)
    bridge = ParameterBridge()

    # 1. Build standard PBPK params (Kp, CLint_liver, fu_b — same as IV)
    params = build_compound_pbpk_params(...)

    # 2. Load GI tract from species YAML
    gi_tract = load_gi_tract(self.species)

    # 3. Compute gut CLint (ratio approach)
    clint_gut = compute_gut_clint(params, gi_tract, self.compound)

    # 4. Build oral PBPK params (extends CompoundPBPKParams)
    oral_params = OralPBPKParams(
        **params.__dict__,
        clint_gut_L_h=clint_gut,
        gi_tract=gi_tract,
        peff_cm_s=...,  # from compound permeability properties
    )

    # 5. Resolve Peff (§12.2 priority chain)
    peff = resolve_peff(self.compound)  # raises ValueError if unavailable

    # 6. Simulate
    sim = simulate_oral(topology, oral_params, dose_mg=self.dose_mg,
                        duration_h=self.duration_h)

    # 7. Extract PK + Fa/Fg/Fh
    pk = compute_oral_pk_parameters(sim, dose_mg=self.dose_mg)

    return PipelineResult(
        compound=self.compound,
        pk_parameters=pk,
        time_h=sim.time_h,
        cp_plasma=sim.cp_plasma,
        cp_blood=sim.cp_blood,
        simulation=sim,
        metadata={...},  # route, dose, Fg, Fa, etc.
    )
```

### 12.2 Peff resolution

Priority chain:
1. `compound.properties.permeability.peff_cm_s.value` (if not None)
2. `papp_to_peff(compound.properties.permeability.papp_nm_s.value)`
   (Sun 2002 correlation)
3. `raise ValueError("oral route requires Peff or Papp")`

---

## 13. Known Limitations and Honest Gaps

### 13.1 MPPGI calibration

MPPGI = 6.0 mg/g is calibrated against midazolam Fg=0.57. This is
within the literature range (2–20) but is NOT an independently measured
physiological constant. Future experimental MPPGI measurements can
directly replace this value.

### 13.2 Pooled enterocyte

All absorbed drug pools into a single enterocyte compartment with
uniform CYP3A4 density. In reality, CYP3A4 is concentrated in the
jejunum and sparse in the ileum/colon. For drugs absorbed primarily in
the jejunum (most orally bioavailable drugs), this simplification
overstates gut metabolism of the small fraction absorbed distally.
Impact is small: colon ka_fraction=0.05 means <2% of dose is affected.

### 13.3 No recirculation gut metabolism

Drug that returns to the gut wall via systemic arterial blood does NOT
encounter enterocyte CYP3A4. This is standard ACAT practice and has
minimal impact (gut CYP3A4 content ≈ 0.45% of hepatic by mass).

### 13.4 CLint overprediction propagates to Fh and F

Sprint 3b Session 1 measured AAFE_CL = 4.89. This CLint overprediction
directly causes Fh underprediction, which makes overall F predictions
unreliable. **Fg validation must be performed independently of Fh**
using the post-hoc flux integration method (§6.2).

### 13.5 IR assumption

Session 2a assumes immediate-release (all drug dissolved at t=0).
BCS II/IV drugs with dissolution-limited absorption will have Fa
overpredicted. All three validation compounds are high-permeability
with adequate solubility for the validation doses used.

---

## 14. Definition of Done

- [ ] `Pipeline.run(route="oral")` produces a PipelineResult with
      finite positive Cp-time profile, Cmax, Tmax, AUC, and Fa/Fg/Fh.
- [ ] midazolam PO 5 mg: Fg within 10% of 0.57 (calibration check).
- [ ] felodipine PO 10 mg: Fg within 2-fold of 0.45 (independent).
- [ ] nifedipine PO 10 mg: Fg within 2-fold of 0.78 (independent).
- [ ] All 3 compounds: Fa > 0.90 (high-Peff sanity).
- [ ] Mass balance residual < 1% of dose for all oral simulations.
- [ ] All existing 543 tests pass (regression).
- [ ] theophylline strict gate PASS (regression).
- [ ] Obach panel AAFE_Vss sanity floor < 5.0 (regression).
- [ ] New oral tests: unit + integration, ≥30 new tests.
