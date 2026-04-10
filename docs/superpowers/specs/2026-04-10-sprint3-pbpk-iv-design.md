# Sprint 3 — Minimal IV-only Human PBPK — Design

**Date:** 2026-04-10
**Scope:** Charon Sprint 3 (Layer 2 PBPK), IV-only subset
**Status:** Design approved, autonomous execution mode
**Supersedes:** N/A (first Sprint 3 design doc)

---

## 1. Session Goal

Ship a working human PBPK simulation engine end-to-end **for the IV route**,
so that `Pipeline(compound, route="iv_bolus", dose_mg=5).run()` produces a
plasma concentration–time profile and PK parameters (Cmax, AUC, CL, Vss,
half-life). The session target is a clean, well-tested kernel that the
subsequent Sprint 3b session can extend with ACAT / oral / formulation.

**Out of scope** (deferred to Sprint 3b): ACAT GI segments, oral route,
dissolution/BCS/food-effect, rat/dog PBPK, solver benchmarking for
wall-time target. Phase B/C features are out of scope entirely.

## 2. Validation Target (Definition-of-Done)

**Primary compound:** midazolam IV bolus, 5 mg human dose.

**Inputs** (experimental overrides, mirroring `test_parameter_bridge.py`
fixture values so hepatic math is internally consistent):

| Parameter        | Value     | Source                              |
|------------------|-----------|-------------------------------------|
| SMILES           | `C[N+]1=C(CN2CC(=NCc3ccccc31)c4ccc(F)cc4)Cl` *(illustrative)* | RDKit parse; canonical form used |
| molecular_weight | 325.77    | PubChem                             |
| logP             | 3.89      | Obach 1999                          |
| pKa_base         | 6.2       | Tertiary amine (midazolam pKa)      |
| compound_type    | `"base"`  | Explicit override (midazolam-base)  |
| fu_p             | 0.03      | Obach 1999                          |
| fu_inc (HLM)     | 0.96      | Austin 2002-ish for logP=3.89       |
| BP ratio         | 0.66      | Obach 1999                          |
| CLint (HLM)      | 93 μL/min/mg | Obach 1999                       |
| CL_renal         | 0.0 L/h   | Midazolam negligible renal clearance|

**Observed reference** (literature, healthy adults):

- CL ≈ 21 L/h
- Vss ≈ 66 L
- t½ ≈ 3 h

**Expected PBPK prediction** (hand-calculated from inputs above, well-stirred):

- CLh (via parameter_bridge) ≈ 13.5 L/h
  - CLint_u = 93/0.96 = 96.9 μL/min/mg
  - CLint_liver = 96.9 × 40 × 1500 = 5.81e6 μL/min = 348.75 L/h
  - fu_b = 0.03/0.66 = 0.0455
  - CLh_ws = 90 × 0.0455 × 348.75 / (90 + 0.0455 × 348.75) = 13.48 L/h
- CL (observed/predicted fold error) = 21/13.5 = 1.56 → within 2-fold ✓
- Vss: driven by Σ V_i × Kp_i (R&R base branch). Expected ~50–80 L;
  within 2-fold of 66 L required.

**Pass criteria (this session):**

1. `Pipeline(midazolam_compound, route="iv_bolus", dose_mg=5).run()` produces a
   plasma Cp-time profile without ODE solver failure.
2. Predicted CL within 2-fold of 21 L/h (target ~13.5 L/h, well within 2-fold).
3. Predicted Vss within 2-fold of 66 L (target 30–130 L).
4. Mass balance: amount_eliminated + amount_in_body == dose (within 1e-6 × dose).
5. BDF solver method is enforced (asserted via test).
6. All new unit + integration tests pass; all 390 pre-existing tests still pass.
7. Coverage: `src/charon/pbpk/` (excluding empty formulation/ stubs) ≥ 80%.

## 3. Architecture

### 3.1 Compartmental topology

**Perfusion-limited, 15 tissues, single compartment per tissue.** Blood
pools: one venous, one arterial. Lung receives full cardiac output from
the venous pool and discharges into arterial.

```
   ┌─────────────── venous pool ◄──────────────────────────────┐
   │                   │                                       │
   │                   ▼                                       │
   │                 lung                                      │
   │                   │                                       │
   │                   ▼                                       │
   │              arterial pool ───────────┐                   │
   │                   │                   │                   │
   │                   │ Q_i  (parallel systemic tissues)      │
   │                   │                                       │
   │           ┌───────┼────────┬──────────┬───────┐           │
   │           ▼       ▼        ▼          ▼       ▼           │
   │        brain   heart   kidney ... (non-portal)            │
   │           │       │        │          │       │           │
   │           └───────┴────────┴──────────┴───────┘           │
   │                           │                               │
   │                           └───── to venous pool ──────────┘
   │
   │           ┌────────────┬──────────┐
   │           ▼            ▼          ▼
   │        spleen        gut_wall   pancreas  (portal pre-hepatic)
   │           │            │          │
   │           └────────────┴──────────┘
   │                        │ (portal vein mix)
   │                        ▼
   │                      liver  ──── hepatic elimination (CLint_liver × fu_b)
   │    Q_HA (from arterial)│
   │   ────────────────────►│
   │                        ▼
   └──────────────── hepatic vein → venous pool
```

**Portal / hepatic arterial split** — the liver receives:

- Hepatic artery: `Q_HA = Q_liver_total − Σ Q_portal_tissues`
- Portal vein: mixed outflow of spleen + gut_wall + pancreas

For human.yaml: Q_liver_total=0.255 CO, Q_portal=0.19 CO, Q_HA=0.065 CO.
These are consistency-checked at topology build time; mismatches raise.

**Kidney** — receives arterial blood, drains to venous, has a renal
elimination term inside the kidney compartment.

**State vector layout** (deterministic, preserved by OrderedDict):

```
y[0] = A_venous_mg
y[1] = A_arterial_mg
y[2 : 2 + N_tissues] = A_tissue_i_mg   (order from YAML)
```

Total states: 2 + len(tissues). For human YAML with 15 tissues: 17 states.

### 3.2 Mass balance equations

All equations are in **blood-flow / blood-concentration** form. Tissue state
is total drug amount A (mg). Concentration leaving tissue is a
**blood-equivalent** concentration so it can be multiplied by blood flow Q.

Let:

- `V_i` = tissue volume (L)
- `Q_i` = blood flow to tissue (L/h)
- `Kp_i` = tissue:plasma partition coefficient (from Rodgers & Rowland)
- `BP` = blood:plasma ratio
- `C_tissue_i = A_i / V_i` (total tissue concentration, mg/L)
- `C_blood_out_i = C_tissue_i × BP / Kp_i`
  (blood concentration leaving tissue, perfusion-limited equilibrium)

**Venous pool:**

```
V_ven × dC_ven/dt = Σ_(non-portal,i≠liver,i≠lung) Q_i × C_blood_out_i
                  + Q_liver_total × C_blood_out_liver
                  − Q_CO × C_ven
```

(Note: kidney is included in the non-portal sum; renal elimination happens
**inside** the kidney compartment before blood leaves it.)

**Arterial pool:**

```
V_art × dC_art/dt = Q_CO × C_blood_out_lung − Q_CO × C_art
```

**Lung:**

```
dA_lung/dt = Q_CO × C_ven − Q_CO × C_blood_out_lung
```

**Non-portal, non-liver, non-kidney systemic tissue:**

```
dA_i/dt = Q_i × C_art − Q_i × C_blood_out_i
```

**Kidney** (adds renal elimination; CL_renal is a plasma clearance, so
divide blood-equivalent by BP to get plasma-equivalent):

```
dA_kidney/dt = Q_k × C_art − Q_k × C_blood_out_kidney
             − CL_renal × (C_blood_out_kidney / BP)
```

**Portal pre-hepatic tissue (spleen, gut_wall, pancreas):**

```
dA_p/dt = Q_p × C_art − Q_p × C_blood_out_p
```

(Its outflow feeds the liver, not the venous pool.)

**Liver** (receives hepatic artery + mixed portal; hepatic elimination as
well-stirred intrinsic term):

```
portal_inflow   = Σ_(p in portal) Q_p × C_blood_out_p
hepatic_elim    = CLint_liver_L_h × fu_b × C_blood_out_liver
dA_liver/dt     = Q_HA × C_art + portal_inflow
                 − Q_liver_total × C_blood_out_liver
                 − hepatic_elim
```

**Well-stirred self-check** (derived by the liver rhs at steady state):

```
0 = Q_H × C_art − Q_H × C_liver_out − CLint × fu_b × C_liver_out
→ C_liver_out = C_art × Q_H / (Q_H + fu_b × CLint)
→ Rate_elim   = fu_b × CLint × C_liver_out
              = Q_H × fu_b × CLint / (Q_H + fu_b × CLint) × C_art
              = CLh_ws × C_art     ✓  (matches liver_models.well_stirred)
```

### 3.3 Compound → PBPK params

```
build_compound_pbpk_params(compound, topology, bridge, compound_type=None):
  # 1. Resolve compound_type: explicit override → fall back to _infer_from_pka
  # 2. For each tissue in topology:
  #      Kp[tissue] = compute_kp_rodgers_rowland(logP, pKa, compound_type, comp, plasma)
  # 3. Call bridge.clint_to_clh(...) to get HepaticClearance with full ConversionLog
  # 4. Extract CLint_liver_L_h from the conversion log (or new schema field)
  # 5. CL_renal = bridge.assign_renal_clearance(fu_p, species_gfr, ...)
  # 6. fu_b = fu_p / BP
  # 7. Return frozen CompoundPBPKParams dataclass.
```

### 3.4 Solver

`scipy.integrate.solve_ivp` with:

- `method='BDF'` (hard-required; unit-tested)
- `rtol=1e-6`, `atol=1e-9 * dose_mg` (scaled)
- `dense_output=True`, `t_eval = linspace(0, duration_h, 500)`
  (uniform is good enough for Cmax/AUC integration; a future session
   can move to log-spaced early points for better distribution-phase
   resolution if accuracy demands it)
- On `sol.success == False`: raise `RuntimeError` with solver message

**Initial conditions:**

- `iv_bolus`: `y[0] = dose_mg` (venous pool), all others zero
- `iv_infusion`: all zero; rhs adds `dose_rate_mg_per_h` to `dA_venous/dt`
  while `t < infusion_duration_h`

### 3.5 PK extraction

Inputs: time grid, plasma concentration (divide venous blood by BP).

- `Cmax` = max(Cp)
- `Tmax` = t[argmax(Cp)]
- `AUC_0_last` = trapezoidal integral on (t, Cp)
- `ke` terminal slope: log-linear regression on last 30% of samples
  (minimum 3 points); if slope ≥ 0, flag and omit tail
- `AUC_tail` = Cp_last / ke (0 if slope invalid)
- `AUC_0_inf` = AUC_0_last + AUC_tail
- `AUC_0_24` = trapezoidal on (t≤24, Cp)
- `half_life` = ln(2) / ke
- `CL` = dose_mg / AUC_0_inf
- `AUMC_0_last` = trapezoidal integral on (t, t·Cp)
- `AUMC_tail` = Cp_last × (t_last/ke + 1/ke²)
- `AUMC_0_inf` = AUMC_0_last + AUMC_tail
- `MRT` = AUMC_0_inf / AUC_0_inf
- `Vss` = CL × MRT (bolus) or CL × (MRT − T_inf/2) (infusion)

## 4. File layout

```
src/charon/
├── pipeline.py                    NEW  top-level Pipeline + PipelineResult
├── __init__.py                    UPDATE export Pipeline
├── core/
│   └── schema.py                  UPDATE add clint_liver_L_h to HepaticClearance
├── pbpk/
│   ├── __init__.py                UPDATE export new public API
│   ├── kp_calculator.py           (existing, unchanged)
│   ├── topology.py                NEW  PBPKTopology + load_species_topology()
│   ├── ode_compiler.py            NEW  CompoundPBPKParams + build_rhs()
│   ├── solver.py                  NEW  SimulationResult + simulate_iv()
│   ├── pk_extract.py              NEW  compute_pk_parameters()
│   └── species/                   (existing, unchanged)
└── ...

tests/
├── unit/
│   ├── test_topology.py           NEW  YAML loading, flow consistency, portal split
│   ├── test_ode_compiler.py       NEW  params assembly, rhs signs, mass conservation
│   ├── test_solver.py             NEW  BDF enforced, analytic 1-cpt match, mass bal
│   ├── test_pk_extract.py         NEW  Cmax/AUC/CL/Vss vs hand calcs
│   └── test_parameter_bridge.py   UPDATE verify clint_liver_L_h populated
├── integration/
│   └── test_smiles_to_pk.py       NEW  midazolam IV end-to-end validation
└── ...

validation/benchmarks/
├── metrics.py                     NEW  AAFE, fold_error, within_n_fold
└── layer2_human_pk.py             NEW  midazolam IV runnable benchmark script
```

## 5. Key data types

```python
@dataclass(frozen=True)
class TissueNode:
    name: str
    volume_L: float
    blood_flow_L_h: float
    composition: TissueComposition   # for Kp calc
    parent: Literal["arterial", "portal", "venous"]
    drains_to: Literal["venous", "liver"]
    eliminates: bool                  # True for liver and kidney

@dataclass(frozen=True)
class PBPKTopology:
    species: str
    body_weight_kg: float
    cardiac_output_L_h: float
    hematocrit: float
    venous_volume_L: float
    arterial_volume_L: float
    hepatic_artery_L_h: float         # derived: Q_liver_total - Q_portal_total
    portal_tissue_names: tuple[str, ...]
    tissues: "OrderedDict[str, TissueNode]"
    plasma_composition: TissueComposition

@dataclass(frozen=True)
class CompoundPBPKParams:
    name: str
    molecular_weight: float
    logp: float
    pka_acid: float | None
    pka_base: float | None
    compound_type: str               # neutral | acid | base | zwitterion
    fu_p: float
    bp_ratio: float
    fu_b: float
    clint_liver_L_h: float           # whole-liver intrinsic CL after IVIVE scaling
    cl_renal_L_h: float
    kp_by_tissue: dict[str, float]

@dataclass
class SimulationResult:
    time_h: np.ndarray               # (N,)
    cp_blood: np.ndarray             # (N,) venous blood concentration, mg/L
    cp_plasma: np.ndarray            # (N,) venous plasma concentration, mg/L
    state_trajectory: np.ndarray     # (17, N) for audit / mass balance
    mass_balance_residual: float     # max |sum(y) - expected| / dose
    solver_success: bool
    solver_nfev: int
    solver_method: str               # "BDF"
```

## 6. Risks and mitigations

| Risk | Impact | Mitigation |
|------|--------|------------|
| R&R Kp base-branch threshold (pKa_base < 8.0 → "neutral") misclassifies weak bases like midazolam | Wrong Kp → wrong Vss | Accept explicit `compound_type` override in `build_compound_pbpk_params`; validation test uses `"base"` explicitly |
| Kp clipping at [0.01, 50] caps adipose for highly lipophilic drugs | Underprediction of Vss | Document in design; accept for session 1; Sprint 3b can revisit |
| BDF with perfect perfusion limit → very stiff (instantaneous equilibration) → slow convergence | Wall time > 200 ms | Pre-scale atol by dose; use scipy default jacobian approximation; benchmark in test but not enforce wall time this session |
| fu_p near-zero (< 0.005) → fu_b vanishes → liver elim vanishes → mass balance appears broken | Validation drift | Parameter bridge already warns at fu_p < 0.01; PBPK params assembly re-checks |
| PK extraction terminal slope misbehaves for short simulations | Incorrect CL/Vss | Auto-extend duration check: if `Cp[-1] / Cmax > 0.05`, flag in result metadata and accept (don't fail) |
| Kidney elimination uses plasma-based CL but compartment tracks blood-equiv | Wrong scaling | Explicit `C_blood_out_kidney / BP` in kidney rhs — unit-tested |
| Hepatic elimination term double-applied with parameter_bridge CLh | Major bug | PBPK uses **CLint_liver_L_h**, NOT CLh; unit-tested that rhs does not use HepaticClearance.clh_L_h |

## 7. Autonomous execution protocol

The user has granted autonomous decision authority for this session ("자율적으로
항상 최선의 최고의 선택지를 고르며 진행"). Terminal gates enforced:

1. All 390 pre-existing tests continue to pass (no regressions).
2. All new tests pass on first `pytest` run after each module is written.
3. Pre-commit hooks (if any) pass.
4. Final commit describes deliverables concretely.
5. On any failure, diagnose root cause before retrying; do not --no-verify,
   do not amend published commits.

## 8. Non-goals (explicit)

- No ACAT / oral / dissolution / BCS / food effect.
- No rat / dog / monkey PBPK in this session (species YAMLs stay as-is).
- No solver wall-time target enforcement.
- No Layer 3 scaling, Layer 4 UQ, Layer 0 guardrails changes.
- No changes to Sprint 2 ADME prediction code.
- No Phase B/C scaffolding touched.
