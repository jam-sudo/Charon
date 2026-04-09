# Charon: Open-Source Translational PK Platform

## Structure-to-FIH — Architecture Document v4.0 (Converged)

**Project scope:** Phase A = "Structure-to-FIH" (core pipeline). Phase B = population/trial simulation. Phase C = multi-compound decision support. Part of the series: Omega (ADMET+PBPK) → Sisyphus (topology-compiled ODE) → **Charon** (translational bridge — crossing the preclinical→clinical divide).

-----

## 1. Mission Statement

An open-source Python platform that takes a molecular structure (SMILES) as input and produces a first-in-human (FIH) dose recommendation with quantified uncertainty — replacing fragmented, expensive commercial tools (Simcyp $50-150K/yr, GastroPlus $40-100K/yr) with a unified, reproducible, ML-augmented pipeline.

**Phase Roadmap:**

- **Phase A:** Structure → FIH dose recommendation (single subject, single compound)
- **Phase B:** Population variability, virtual trial simulation, DDI, special populations
- **Phase C:** Multi-compound decision dashboard, Pareto optimization, automated risk assessment

-----

## 2. Design Principles

1. **YAML-driven declarative configuration** — every run reproducible from a single config file
1. **Explicit over implicit** — all unit conversions, scaling factors, model assumptions visible and auditable
1. **Modular independence** — each layer executable standalone, chainable in pipeline
1. **Uncertainty-native** — every prediction carries confidence intervals; uncertainty propagates end-to-end
1. **Fail loudly** — out-of-domain inputs, method disagreements, and known failure modes generate explicit warnings, not silent garbage
1. **Regulatory awareness** — outputs structured to support IND submission (ICH M3, FDA guidance on FIH dose)

-----

## 3. Layered Architecture

### Layer 0: Input & Validation

**Purpose:** Parse molecular input, verify drug-likeness, assess applicability domain.

```
Input: SMILES string (or SDF/MOL file)
       ↓
       RDKit mol object
       ↓
       Guardrails:
         - Structural validity (RDKit sanitization)
         - Drug-likeness filters (Lipinski, Veber, lead-likeness)
           → violations = warning, not hard block
         - Applicability domain check:
           Tanimoto similarity to Layer 1 training set
           → < 0.3 = "LOW reliability" flag
           → 0.3-0.5 = "MODERATE"
           → > 0.5 = "HIGH"
         - Known failure mode detection:
           MW > 800 → "approaching biologic; SM PBPK may be unsuitable"
           cLogP > 6 → "dissolution-limited risk; formulation critical"
           TPSA > 140 → "low permeability; oral route questionable"
           → route-of-administration advisory
       ↓
Output: ValidatedMolecule object with reliability metadata
```

**Design decisions:**

- Guardrails are advisory, not blocking. Rationale: novel chemical space is the whole point of drug discovery; blocking on Lipinski would reject ~30% of approved drugs post-2010.
- Applicability domain uses Tanimoto on Morgan fingerprints (radius=2, 2048 bits) against Layer 1 training set centroids. Simple, interpretable, fast.

-----

### Layer 1: Property Prediction

**Purpose:** Predict all ADMET/physicochemical parameters needed for PBPK.

**Required prediction targets:**

|Parameter         |Unit       |Source                                           |Notes                                            |
|------------------|-----------|-------------------------------------------------|-------------------------------------------------|
|logP              |—          |Omega ADMET ensemble                             |Partition coefficient                            |
|pKa (acid/base)   |—          |Omega pKa predictor                              |Ionization state at physiological pH             |
|Aqueous solubility|μg/mL      |Omega ensemble                                   |Intrinsic + pH-dependent                         |
|Caco-2 Papp       |nm/s       |Omega ensemble                                   |Apparent permeability                            |
|fu_p              |fraction   |Omega ensemble                                   |Plasma unbound fraction                          |
|fu_inc            |fraction   |NEW: Hallifax-Houston correlation from logP/fu_p |Microsomal unbound fraction — CRITICAL for IVIVE |
|B:P ratio         |ratio      |NEW: ML model or empirical (pKa + fu_p based)    |Blood-to-plasma ratio                            |
|CLint             |μL/min/mg  |Omega ensemble                                   |In vitro intrinsic clearance (HLM)               |
|CYP substrate ID  |categorical|Omega ensemble                                   |Primary/secondary CYP isoforms                   |
|CYP inhibition    |IC50, μM   |Omega ensemble                                   |DDI perpetrator potential                        |
|hERG IC50         |μM         |Omega ensemble                                   |Cardiac safety flag                              |
|P-gp substrate    |binary     |Omega ensemble                                   |Efflux transporter flag                          |
|OATP substrate    |binary     |NEW: ML model                                    |Hepatic uptake transporter                       |
|Renal clearance   |L/h        |NEW: GFR-filtration model + active secretion flag|CLrenal = fu_p × GFR × (1 + net secretion factor)|

**Conformal prediction:**

- Every continuous prediction carries [ci_lower, ci_upper] at 90% coverage.
- Conformal calibration set maintained separately from training set.
- Coverage verified per-parameter: if empirical coverage < 85%, recalibrate.

**Known gaps and mitigations:**

- fu_inc: No large public training set. Use Hallifax & Houston 2006 logP correlation as default; allow user override with experimental value.
- B:P ratio: Empirical model: B:P = 1 + HCT × (fu_p × Kp_rbc - 1). Requires Kp_rbc estimation from pKa and lipophilicity.
- OATP substrate: Limited training data. Binary classifier with explicit "uncertain" output when probability ∈ [0.35, 0.65].

-----

### Layer 1→2 Bridge: Parameter Translation

**Purpose:** Convert raw predictions into PBPK-ready parameters with explicit unit conversions and audit trail.

This is the most error-prone junction in the pipeline. Omega's fu_plasma double-application bug originated here. Design principle: **every conversion is a named, tested, logged function.**

**Scope boundary:** ParameterBridge handles human IVIVE only (CLint→CLh for human). Species-specific IVIVE lives in translational/ivive.py (Layer 3). This avoids duplication — ParameterBridge calls the same well-stirred/parallel-tube math via a shared `liver_models.py` utility, parameterized by species.

**Data flow through bridge:**

```
Layer 1 outputs → ParameterBridge → PBPK-ready params
                                   ↘ Kp calculator (parallel)
                                     Inputs: logP, pKa, fu_p, B:P
                                     + tissue composition from species YAML
                                     → Kp per tissue → feeds ODE
```

**Source taxonomy (standardized):**

```
source values:
  "ml_ensemble"     — ML model prediction (Omega ADMET ensemble)
  "ml_pka"          — ML pKa predictor
  "correlation"     — literature empirical correlation (e.g. Hallifax-Houston)
  "derived"         — computed from other predicted values (e.g. Peff from Papp)
  "physiological"   — literature physiological parameter
  "experimental"    — user-provided measured value
```

**Critical conversions:**

```python
class ParameterBridge:

    def clint_to_clh(self,
        clint: float,          # μL/min/mg microsomal protein (or per 10⁶ cells)
        fu_inc: float,         # microsomal unbound fraction
        fu_p: float,           # plasma unbound fraction
        system: str,           # "HLM" or "hepatocytes"
        mppgl: float,          # mg protein per g liver (HLM: default 40)
        hepatocellularity: float,  # 10⁶ cells/g liver (hepatocytes: default 120)
        liver_weight_g: float, # species-specific (human: ~1500g)
        qh_L_h: float,        # hepatic blood flow (human: ~90 L/h)
        bp_ratio: float,       # blood:plasma ratio
        model: str = "well_stirred"
    ) -> HepaticClearance:
        """
        IVIVE: in vitro CLint → in vivo hepatic clearance

        Step 0: Select scaling factor based on in vitro system
                if system == "HLM": scale_factor = MPPGL × liver_weight
                if system == "hepatocytes": scale_factor = hepatocellularity × liver_weight
                (fu_inc correction applies to HLM only; hepatocyte CLint
                 assumed already reflecting unbound if fu_inc not measured)

        Step 1: Unbound CLint = CLint / fu_inc  [μL/min/mg] (HLM only)

        Step 2: Scale to whole liver
                CLint_liver = CLint_u × scale_factor
                Convert: ÷ 1e6 × 60 → [L/h]

        Step 3: Apply liver model
                Well-stirred: CLh = (Qh × (fu_p/BP) × CLint_liver) /
                                    (Qh + (fu_p/BP) × CLint_liver)
                where fu_p/BP = fu_b (blood unbound fraction)

        Returns CLh in L/h with complete audit of each step.
        """

    def papp_to_peff(self,
        papp_nm_s: float,
        calibration: str = "sun_2002"
    ) -> float:
        """
        Caco-2 Papp → human effective permeability (Peff)
        Sun et al. 2002: log Peff = 0.4926 × log Papp - 0.1454
        Returns Peff in cm/s (×10⁻⁴)
        """

    def assign_renal_clearance(self,
        fu_p: float,
        gfr_mL_min: float = 120.0,  # human default
        is_active_secretion: bool = False,
        net_secretion_factor: float = 0.0
    ) -> float:
        """
        CLrenal = fu_p × GFR × (1 + net_secretion_factor)
        Unit conversion: GFR [mL/min] × 60/1000 → [L/h]
        If no active secretion predicted: net_secretion_factor = 0
        → CLrenal = fu_p × GFR (filtration only)
        Returns L/h

        Example: fu_p=0.23, GFR=120 mL/min
        → 0.23 × 120 × 60/1000 = 1.66 L/h
        """
```

**Audit trail:**

- Every ParameterBridge call returns a `ConversionLog` object:
  
  ```
  ConversionLog(
    input_params={...},
    intermediate_steps=[
      ("CLint_u", 17.9, "μL/min/mg", "CLint/fu_inc = 15.2/0.85"),
      ("CLint_liver", 1074000, "μL/min", "CLint_u(17.9) × MPPGL(40) × Wliver(1500g)"),
      ("CLint_liver_Lh", 64.4, "L/h", "1074000 ÷ 1e6 × 60"),
      ("CLh", 13.3, "L/h", "well-stirred: Qh=90, fu_b=fu_p/BP=0.23/0.95=0.242")
    ],
    output=13.3,
    output_unit="L/h",
    model_used="well_stirred"
  )
  ```
- Stored per-run for debugging and report generation.

-----

### Layer 2: PBPK Simulation Engine

**Purpose:** Simulate human (and preclinical species) plasma concentration-time profiles.

**Architecture:** Sisyphus topology-compiled ODE system. 34 tissue-level compartments in the directed multigraph; with sub-compartmentalization (vascular/interstitial/intracellular per tissue), total state variables ~50.

**Compartments:**

- 15 tissue compartments (adipose, bone, brain, gut, heart, kidney, liver, lung, muscle, pancreas, skin, spleen, reproductive, thymus, rest)
- 8 ACAT GI segments (stomach, duodenum, jejunum-1/2, ileum-1/2/3, colon)
- Per-tissue subcompartments: vascular, interstitial, intracellular (×3 per tissue = 45)
- Portal vein, hepatic artery, systemic venous, systemic arterial
- Total state variables: ~50 (tissue partition + GI transit + dissolved/undissolved drug)

**Tissue partition coefficient (Kp) calculation:**

- Rodgers & Rowland (2005/2006) method
- Inputs: logP, pKa, fu_p, B:P, tissue composition data (fw, fn, fp, fap per tissue per species)
- Separate models for ionized vs neutral, acid vs base vs neutral vs zwitterion

**ODE compilation:**

- YAML compound config → symbolic ODE system (SymPy or manual)
- Topology compiled from species YAML (organ connections, blood flows)
- Solver: `scipy.integrate.solve_ivp` with BDF method (stiff systems)
- Dosing events: IV bolus, IV infusion, oral (via ACAT), handled as discrete events in ODE

**Route-dependent pipeline logic:**

```
if route == "oral":
    run ACAT dissolution → permeation → gut wall metabolism
    F = Fa × Fg × Fh
    report: Cmax, Tmax, AUC, t½, CL/F
elif route == "iv_bolus" or route == "iv_infusion":
    skip ACAT entirely; F = 1.0 (by definition for IV)
    inject dose directly into venous blood compartment
    report: AUC, t½, CL, Vss, Vd
    dose projection: NOAEL/HED based on IV NOAEL (mg/kg iv)
```

**Multiple dosing / steady-state:**

```
Simulation modes:
  single_dose: one dose event, simulate adaptively
    → initial t_end = 72h (default for most small molecules)
    → if Cp at t_end > 5% of Cmax: extend to 168h (1 week)
    → if still > 5%: extend to 336h and flag "long t½ compound"
    → report simulation duration used
  multiple_dose: dose every tau hours, simulate to steady state
    → detect steady state: AUC_tau(n) vs AUC_tau(n-1) < 5% diff
    → report: Css_max, Css_min, Css_avg, accumulation ratio (Rac)
    → Rac_predicted = 1 / (1 - exp(-ke × tau))
  Default for FIH: single_dose (Phase I SAD)
  Available for: MAD simulation, DDI (Phase B)
```

**Low CLint compound routing:**

```
if CLint < 1.0 μL/min/mg:
    flag: "Low hepatic turnover detected"
    → check: is CLrenal > CLhepatic?
       if yes: "Renal elimination dominant; 
                hepatic PBPK CL may underpredict total CL"
    → check: extrahepatic metabolism possible?
       (CYP3A4 substrate + low liver CLint → consider gut/kidney CYP)
    → report: dominant elimination pathway assessment
```

**ACAT (Advanced Compartmental Absorption Transit) model:**

- 8-segment GI tract (reused from Omega)
- Per-segment: dissolution (Noyes-Whitney), permeation (Peff), transit (first-order), degradation
- Dissolution: particle-size dependent; Noyes-Whitney with Z-factor or Wang-Flanagan
- Precipitation: supersaturation → amorphous precipitation in stomach/duodenum for weak bases

**Outputs per simulation:**

```
SimulationResult:
  cp_time_profile: ndarray     # [time_h, Cp_ng_mL]
  pk_parameters:
    cmax: float                 # ng/mL
    tmax: float                 # h
    auc_0_inf: float            # ng·h/mL
    auc_0_24: float             # ng·h/mL
    half_life: float            # h
    cl_apparent: float          # L/h (oral) or L/h (IV)
    vss: float                  # L (IV only)
    bioavailability: float      # fraction (Fa × Fg × Fh)
    fa: float                   # fraction absorbed from GI
    fg: float                   # fraction escaping gut wall metabolism
    fh: float                   # fraction escaping hepatic first-pass
  mass_balance:
    total_absorbed: float
    total_metabolized_gut: float
    total_metabolized_liver: float
    total_excreted_renal: float
    total_remaining_gi: float
  compartment_profiles: dict    # optional: Ct-time for each tissue
```

**ODE integration of CLrenal:**

- Kidney compartment has an elimination term: dA_kidney_elim/dt = CLrenal × C_kidney_unbound
- CLrenal from Layer 1 (fu_p × GFR + active secretion) enters ODE as a fixed parameter
- Renal elimination fraction tracked in mass_balance.total_excreted_renal

**Fg extractability:**

- Fg (gut wall first-pass) is computed WITHIN the ACAT-PBPK ODE (gut wall CYP metabolism)
- For independent testing, extract Fg post-hoc: Fg = amount_entering_portal_vein / amount_absorbed
- Dedicated test: compare predicted Fg for CYP3A4 substrates (midazolam, felodipine) against literature
- This is the highest Sobol sensitivity driver (ST=0.470 in Omega) — prioritize validation

**Formulation sub-module:**

```
formulation/
├── dissolution.py
│   Noyes-Whitney: dM/dt = (D/h) × S × (Cs - C)
│   D = aqueous diffusion coefficient, estimated via
│       Hayduk-Laudie correlation: D = f(MW, viscosity)
│   S = surface area (from particle size distribution, spherical assumption)
│   h = diffusion layer thickness (typically 30 μm, or Hintz-Johnson model)
│   Input: particle size distribution, dose, solubility, MW
│   → dissolution rate profile → feeds ACAT
│
├── biopharmaceutics.py
│   BCS classification:
│     Dose Number (Do) = Dose / (solubility × 250 mL)
│     Do ≤ 1 → high solubility; Do > 1 → low solubility
│     Papp threshold (≥ metoprolol ~20 nm/s) → high/low permeability
│     → BCS I/II/III/IV auto-assignment
│   Formulation risk flag:
│     BCS II/IV → "dissolution-limited; formulation development critical"
│     BCS III → "permeability-limited; absorption enhancer may be needed"
│
└── food_effect.py
    Bile salt solubilization factor for BCS II
    Gastric emptying delay (fed: ~2h vs fasted: ~0.5h)
    → predicted fed/fasted AUC ratio
```

-----

### Layer 3: Translational Scaling Engine

**Purpose:** Project human PK from preclinical data using multiple methods; recommend FIH dose.

#### 3a. Species PBPK

Run Layer 2 PBPK with species-specific physiological parameters.

**Species parameter sets (YAML):**

|Parameter                |Human|Rat |Dog |Monkey|
|-------------------------|-----|----|----|------|
|Body weight (kg)         |70   |0.25|10  |5     |
|Cardiac output (L/h)*    |390  |5.4 |120 |50    |
|Liver weight (g)         |1500 |9.5 |320 |120   |
|Liver blood flow (% CO)  |25.5 |18.3|30.4|27.8  |
|GFR (mL/min)             |120  |1.31|3.7 |2.1   |
|MPPGL (mg/g)             |40   |45  |36  |38    |
|Hepatocellularity (10⁶/g)|120  |120 |130 |122   |

**Human CO varies by source: Brown 1997 = 336 L/h, ICRP = 348 L/h, some PBPK tools use 390 L/h. All values configurable in species YAML. Choose reference and document in report.*

**CYP ortholog mapping:**

|Human |Rat       |Dog    |Monkey |
|------|----------|-------|-------|
|CYP3A4|CYP3A1/3A2|CYP3A12|CYP3A8 |
|CYP2D6|CYP2D1/2D2|CYP2D15|CYP2D17|
|CYP2C9|CYP2C11   |CYP2C21|CYP2C43|
|CYP1A2|CYP1A2    |CYP1A2 |CYP1A2 |

**ISEF (inter-system extrapolation factor):**

- Corrects for in vitro ≠ in vivo enzyme activity differences.
- Literature-derived per enzyme per species.
- Major uncertainty source; default values with ranges from Proctor et al. 2004.

**GI tract species differences:**

- Rat: small intestine ~104 cm, transit ~3.5h, no gallbladder (continuous bile)
- Dog: ~350 cm, transit ~2h, higher gastric pH fasted (~1 vs rat ~3.9)
- Human: ~600 cm, transit ~3.5h

**Phase A scope:** human fully built; rat and dog by Sprint 4. Monkey deferred to Phase B.

#### 3b. IVIVE (In Vitro-In Vivo Extrapolation)

```
IVIVE pathway:
  in vitro CLint (HLM or hepatocytes)
    → fu_inc correction
    → whole liver scaling (MPPGL × liver weight)
    → liver model (well-stirred / parallel-tube / dispersion)
    → predicted in vivo CLh

Three liver models available:
  1. Well-stirred: simplest, most common
     CLh = Qh × fu_b × CLint / (Qh + fu_b × CLint)

  2. Parallel-tube: assumes plug flow
     CLh = Qh × (1 - exp(-fu_b × CLint / Qh))

  3. Dispersion: intermediate, requires dispersion number (DN)
     Most physiologically accurate but requires additional parameter

Default: well-stirred (industry standard for IVIVE)
User can switch; output reports which model was used.
```

#### 3c. Multi-Method Consensus Scaling

```python
class ConsensusScaler:
    """
    Runs up to 4 scaling methods and produces consensus estimate.

    Methods:
    1. Simple Allometry (SA)
       log CL = a + b × log BW
       Minimum 3 species recommended; 2 species (rat+dog) usable
       in Phase A with confidence=LOW and monkey-addition advisory.
       Known limitations: CL/F only if oral across species.

    2. Rule of Exponents (ROE)
       Based on allometric exponent b:
         b ≤ 0.55 → use SA directly
         0.55 < b ≤ 0.70 → multiply by MLP (max lifespan potential)
         b > 0.70 → multiply by BrW (brain weight)
       Auto-selects correction based on exponent.

    3. Single-species PBPK scaling
       Run IVIVE-PBPK in each preclinical species.
       Compare predicted vs observed animal PK.
       Apply same IVIVE-PBPK to human parameters.
       "If the model predicts rat/dog well, trust the human prediction."

    4. fu-corrected intercept method
       Normalize CL by fu_p before allometric regression.
       Useful when fu varies substantially across species.
       Known to improve predictions for highly bound drugs.

    Consensus logic:
      1. Run all applicable methods (not all require same input).
      2. Compute geometric mean of predictions.
      3. Compute inter-method CV.
      4. If CV < 50%: confidence = HIGH → report geometric mean.
         If CV 50-100%: confidence = MEDIUM → report geometric mean
           with warning; flag divergent method.
         If CV > 100%: confidence = LOW → default to PBPK-based
           estimate (method 3) as primary; others as supporting.
      5. Rationale: PBPK is mechanistic; empirical methods fail
         for CYP3A4 substrates (dog overpredicts), transporters,
         and nonlinear PK. When methods disagree, mechanism wins.

    Exception handling:
      - Prodrug: scale active metabolite, not parent.
        → requires metabolite SMILES (user input) or flag.
      - Enzyme saturation: if predicted Km < Cmax at therapeutic dose,
        flag nonlinear PK risk → allometry invalid.
      - Transporter substrates (flagged in Layer 1):
        allometry known to fail → weight PBPK method higher.
    """
```

#### 3d. FIH Dose Projection

```
Three parallel dose projection approaches:

1. NOAEL → HED → MRSD
   HED = NOAEL_mg_kg × (Km_animal / Km_human)
   where Km = BW_kg / BSA_m2
   Species Km values (FDA Guidance, 2005):
     Human: 37, Rat: 6.2, Dog: 20, Monkey: 12
   MRSD = HED / safety_factor
   Default safety_factor = 10 (standard)
   Adjustable: 3-30 depending on target class and toxicity profile

2. MABEL (Minimum Anticipated Biological Effect Level)
   PK-driven approach (no PD model needed):
     Target unbound Css: Css_u ≥ Kd
     → Css_total = Kd / fu_p
     → For repeat dosing: Dose = Css_total × CL × tau / F
        where tau = dosing interval (h)
     → For single-dose Cmax-based: Dose = Cmax_target × Vd / (F × dose_fraction_to_Cmax)
     → apply safety factor (typically 10×)
   Required input: target Kd or EC50 (user-provided), dosing regimen assumption
   If not provided: skip MABEL, note in report.

3. PAD (Pharmacologically Active Dose)
   Target Ceff → required AUC → dose via PBPK
     Dose = (AUC_target × CL) / F
   Similar to MABEL but efficacy-driven.
   Required: target efficacious concentration (user-provided)

Final MRSD selection:
  MRSD = min(NOAEL-based, MABEL, PAD)
  → most conservative estimate
  → report all three with rationale for final selection
  → uncertainty interval from Layer 4 applied to final MRSD
```

-----

### Layer 4: Uncertainty Quantification

**Purpose:** Propagate prediction uncertainty through entire pipeline to produce dose recommendation with confidence interval.

#### Uncertainty Sources

```
Source                  | Quantification
------------------------|-------------------------------------------
Layer 1 predictions     | Conformal prediction intervals (90%)
IVIVE scaling factors   | MPPGL: 40 ± 8 mg/g (CV ~20%)
                        | Hepatocellularity: 120 ± 24 ×10⁶/g
Liver model selection   | Well-stirred vs parallel-tube: ~30% difference
Species allometry       | Inter-method CV from consensus scaler
Physiological params    | Population CV from literature
Safety factor           | Fixed (not uncertain), but sensitivity tested
```

#### Sampling Strategy

```
Phase A implementation: Latin Hypercube Sampling (LHS)

Why LHS over alternatives:
  - Naive MC: needs N=5000+ for convergence → slow with ODE
  - LHS: N=500 sufficient for 8-10 parameters → 10× faster
  - PCE: elegant but curse of dimensionality >8 params
  - Quasi-MC (Sobol sequences): comparable to LHS, considered for Phase B

Correlation handling:
  Parameters are NOT independent:
    logP ↔ fu_p: negative correlation (r ≈ -0.6)
    logP ↔ Papp: positive correlation (r ≈ 0.4)
    logP ↔ solubility: negative correlation (r ≈ -0.7)
  
  Ignoring correlations → physically impossible parameter combos
  (e.g., logP = 5 but fu_p = 0.95 = nonsensical)

  Solution: Iman-Conover rank correlation method
    1. Generate N LHS samples (uncorrelated)
    2. Apply rank-based transformation to induce target correlation
    3. Preserves marginal distributions while adding dependence
    4. Target correlation matrix from literature or training data
```

#### Propagation Pipeline

```
For each of N=500 LHS samples:
  1. Draw parameter vector θ_i = {logP_i, fu_p_i, CLint_i, ...}
  2. Run ParameterBridge(θ_i) → PBPK parameters
  3. Run PBPK ODE solver → PK profile
  4. Extract PK metrics (CL, Vss, AUC, Cmax)
  5. Run dose projection → FIH dose estimate

Aggregate:
  - Geometric mean of N dose estimates → point estimate
  - 5th and 95th percentiles → 90% CI
  - Per-parameter Sobol sensitivity indices (reused from Omega)
    → "fu_p uncertainty accounts for 42% of dose variance"
```

#### Output

```yaml
dose_recommendation:
  point_estimate_mg: 45
  ci_90_lower_mg: 28
  ci_90_upper_mg: 72
  confidence: "MEDIUM"
  limiting_uncertainty:
    parameter: "fu_p"
    sobol_total: 0.42
    recommendation: "Experimental fu_p measurement would narrow CI by ~40%"
  all_sobol_indices:
    fu_p: 0.42
    CLint: 0.31
    Papp: 0.15
    logP: 0.08
    MPPGL: 0.04
  n_samples: 500
  convergence_check: "CV of mean stabilized at N=350"
```

**Numerical robustness notes:**

- **fu_p sensitivity amplification:** When fu_p → 0 (highly bound drugs), CLh ≈ Qh × (fu_p/BP) × CLint / Qh ≈ fu_p/BP × CLint. Small absolute changes in fu_p cause large relative changes in CLh. This is the root cause of Omega's Sobol ST=0.470 for fu_p. Mitigations: (1) log-transform fu_p before sampling; (2) flag when fu_p < 0.01 as "extreme binding — CLh highly sensitive to fu_p accuracy."
- **Conformal prediction coverage caveat:** Conformal prediction guarantees marginal coverage (averaged over all compounds), not conditional coverage (for a specific compound). For OOD compounds, actual coverage may be <<90%. Applicability domain check in Layer 0 partially mitigates but cannot fully solve this. Report should state: "confidence intervals assume compound is within the applicability domain of the training set."
- **ODE initial conditions:** All tissue compartments initialized at 0. For oral dosing, initial condition = dose placed in stomach (dissolved fraction based on solubility at gastric pH; undissolved remainder as solid particles). For IV bolus, initial condition = dose in venous blood compartment. For IV infusion, zero-rate forcing function active during infusion period.

-----

### Layer 5: Population Variability (Phase B)

#### 5a. Virtual Population Generator

```
Demographic distributions:
  body_weight: log-normal, μ=70kg, σ=15kg (adults)
  height: normal, μ=170cm, σ=10cm
  age: configurable (18-65 for Phase I, broader for Phase III)
  sex: 50/50 default, adjustable
  bsa: DuBois formula from weight + height
  gfr: CKD-EPI from age + sex + serum creatinine distribution
  albumin: normal, μ=4.0 g/dL, σ=0.3

CYP phenotype frequencies:
  Ethnic group-specific allele frequencies (PharmGKB/PharmVar)
  Activity scores → metabolizer status (PM/IM/EM/UM)

  Example CYP2D6 (Caucasian):
    PM (*4/*4, *5/*5, etc.): 7%
    IM (*4/*41, etc.): 10%
    EM (*1/*1, *1/*2, etc.): 75%
    UM (*1×N/*1, etc.): 8%

  Example CYP2C19 (Caucasian):
    PM (*2/*2, *2/*3): 3%
    IM (*1/*2, *1/*3): 24%
    EM (*1/*1): 63%
    RM (*1/*17): 5%
    UM (*17/*17): 5%
    (Note: CYP2C19 PM prevalence ~15-20% in East Asian populations)

  → CLint_2d6 × activity_multiplier:
    PM: 0.0, IM: 0.5, EM: 1.0, UM: 2.5

Parameter covariance:
  Weight ↔ liver volume: positive
  Weight ↔ cardiac output: positive
  Age ↔ GFR: negative
  → multivariate normal with literature covariance matrix
```

#### 5b. Population PBPK Simulation

```
For N=1000 virtual subjects:
  1. Sample demographics + CYP phenotype
  2. Adjust physiological parameters (organ volumes, blood flows)
     allometrically from body weight
  3. Adjust enzyme activity from phenotype
  4. Run PBPK ODE for each subject
  5. Collect PK metrics

Performance consideration:
  1000 × ODE solve at ~50ms = ~50 seconds
  Mitigation options (Phase B):
    a. JAX + jax.vmap: batch-parallel ODE on GPU
    b. diffrax: JAX-native DE solver
    c. joblib multiprocessing: CPU parallel (simpler, no GPU needed)
    d. Pre-compiled Numba JIT ODE: ~5× speedup over scipy
  Decision deferred to Phase B. Phase A uses scipy (single subject).

Output:
  Population PK summary:
    geometric mean + 90% prediction interval for each PK parameter
    AUC/Cmax distribution → box plots, histograms
    Stratified by CYP phenotype, renal function, etc.
```

#### 5c. Clinical Scenario Simulation

```
Dose escalation simulation:
  - 3+3 design: simulate DLT rates at each dose level
  - BOIN: Bayesian optimal interval design
  - mTPI: modified toxicity probability interval
  Input: predicted dose-exposure-toxicity relationship
  Output: recommended dose levels for Phase I protocol

DDI overlay:
  - Reversible inhibition: perpetrator Cp → (1 / (1 + [I]/Ki)) → modified CLint
  - Time-dependent inhibition (TDI): KI, kinact → time-course enzyme degradation
  - Induction: Emax/EC50 model for CYP3A4 induction (rifampicin-type)
  - Re-run PBPK → AUC ratio (victim drug)
  - Classify: no effect (<1.25×), weak (1.25-2×), moderate (2-5×), strong (>5×)

Special populations:
  - Renal impairment:
    Mild (GFR 60-90), Moderate (30-60), Severe (15-30)
    → reduced renal CL + uremic protein binding changes
  - Hepatic impairment:
    Child-Pugh A/B/C → reduced liver blood flow, enzyme activity, albumin
    Quantitative adjustments from Edginton & Willmann 2008
  - Pediatric:
    Ontogeny functions for CYP maturation (Upreti & Wahlstrom 2016)
    Allometric organ scaling from body weight
    Age-appropriate physiological parameters
```

-----

### Layer 6: Decision Dashboard (Phase C)

```
Multi-compound comparison:
  - Batch run N compounds through Phase A pipeline
  - Side-by-side PK comparison table
  - Radar chart: ADMET properties overlay

Pareto frontier:
  - Axes: predicted efficacy (Cmax/EC50 ratio) vs safety (hERG margin) vs PK (t½, F)
  - Interactive: user selects which axes matter most
  - Pareto-optimal candidates highlighted

Risk scoring:
  Automated risk flags per compound:
    - hERG IC50 / predicted Cmax < 30× → cardiac risk
    - Reactive metabolite alerts (structural alerts: anilines, furans, thiophenes)
    - CYP TDI (time-dependent inhibition) potential
    - High metabolic turnover (CLint > 100 μL/min/mg)
    - P-gp efflux substrate + CYP3A4 substrate → DDI vulnerability
    → aggregate into Red/Yellow/Green per compound

Report generation:
  Auto-generate IND-ready dose rationale document
  (see Report module below)
```

-----

## 4. Report Generator (Cross-Cutting — architecturally part of Layer system)

**This may be the primary deliverable.** Small biotechs need the document, not the simulation itself. Architecturally, the report generator is an observer that spans all layers, not a standalone module.

```
Architecture: Observer pattern
  - Each Layer registers outputs with ReportCollector
  - After pipeline completes, Collector has all data
  - NarrativeEngine converts data → template-driven text
  - Exporter renders to DOCX/PDF/Markdown

Templates (YAML-defined):

  fih_dose_rationale.yaml:
    title: "First-in-Human Dose Rationale"
    sections:
      1_executive_summary:
        content: auto-generated 1-paragraph summary
        data: [compound_name, mrsd, ci_90, primary_method]

      2_compound_properties:
        content: ADMET prediction table with CIs
        data: [layer_1_results]

      3_nonclinical_pk:
        content: species PK comparison table
        figures: [cp_time_overlay_species]
        data: [layer_2_species_results]

      4_human_pk_prediction:
        content: PBPK methodology, predicted human PK
        figures: [cp_time_human_predicted]
        data: [layer_2_human_results, pk_parameters]

      5_scaling_methodology:
        content: multi-method comparison, consensus rationale
        tables: [method_comparison_table]
        data: [layer_3_consensus_results]

      6_dose_projection:
        content: NOAEL→HED, MABEL, PAD calculations
        tables: [dose_projection_summary]
        data: [layer_3d_results]

      7_uncertainty_analysis:
        content: Sobol indices, sensitivity discussion
        figures: [tornado_plot, dose_distribution_histogram]
        data: [layer_4_results]

      8_risk_assessment:
        content: safety flags, DDI potential, formulation risk
        data: [layer_0_guardrails, layer_1_flags]

      9_conclusion:
        content: recommended MRSD with justification
        data: [final_mrsd, confidence_level]

  Narrative generation: template-based string formatting.
  NO LLM dependency for narrative. Deterministic, reproducible.
  Example: "Based on an in vitro CLint of {clint} μL/min/mg
  ({clint_ci_lower}–{clint_ci_upper}, 90% CI), the predicted
  human hepatic clearance is {clh} L/h using the {liver_model}
  model with MPPGL = {mppgl} mg/g."

  Phase A deliverable: Markdown report + embedded matplotlib figures.
  DOCX export: optional premium feature (python-docx; significant
  formatting effort for tables, figures, headers/footers).
  PDF: fpdf2 from rendered Markdown.
```

-----

## 5. Data Architecture

### YAML Compound Configuration (auto-generated, user-editable)

```yaml
compound:
  name: "CPD-001"
  smiles: "CC(=O)Oc1ccccc1C(=O)O"  # free form SMILES
  molecular_weight: 180.16            # free form MW
  salt_form:                           # optional
    name: null                         # e.g., "hydrochloride", "sodium"
    mw_salt: null                      # e.g., 216.62 for HCl salt
    salt_factor: 1.0                   # MW_free / MW_salt (auto-calculated if mw_salt given)
    # Dose correction: dose_free = dose_salt × salt_factor
    # Applied in dose projection: reported MRSD is in free-form equivalents
    # User can specify dose_mg as salt or free form via pipeline.formulation.dose_basis
  source: "predicted"                    # or "experimental", "mixed"

  properties:
    physicochemical:
      logP: {value: 1.19, ci_90: [0.95, 1.43], source: "ml_ensemble"}
      pka_acid: {value: 3.49, ci_90: [3.1, 3.9], source: "ml_pka"}
      pka_base: null
      solubility_ug_ml: {value: 4500, ci_90: [2800, 7200], source: "ml_ensemble"}
    permeability:
      papp_nm_s: {value: 32.1, ci_90: [18.5, 55.7], source: "ml_ensemble"}
      peff_cm_s: {value: 3.2e-4, source: "derived"}  # from Papp via Sun 2002
      pgp_substrate: {value: false, probability: 0.12, source: "ml_ensemble"}
      oatp_substrate: {value: false, probability: 0.35, flag: "uncertain", source: "ml_ensemble"}
    binding:
      fu_p: {value: 0.23, ci_90: [0.18, 0.29], source: "ml_ensemble"}
      fu_inc: {value: 0.85, source: "correlation"}  # Hallifax-Houston from logP
      bp_ratio: {value: 0.95, source: "correlation"}  # empirical from pKa + fu_p
    metabolism:
      primary_cyp: "CYP2C9"
      secondary_cyp: "CYP3A4"
      clint_uL_min_mg: {value: 15.2, ci_90: [8.1, 22.3], system: "HLM"}
    safety:
      herg_ic50_uM: {value: 45.2, ci_90: [28, 73]}
      cyp_inhibition:
        CYP3A4: {ic50_uM: 100, censored: "right"}  # measured as >100
        CYP2D6: {ic50_uM: 23.5}
    renal:
      clrenal_L_h: {value: 1.66, source: "derived"}  # fu_p(0.23) × GFR(120) × 60/1000
      active_secretion: false

pipeline:
  layers: [0, 1, 2, 3, 4]               # which layers to run
  species: ["human"]                      # Phase A default
  preclinical_species: ["rat", "dog"]     # if preclinical data available
  liver_model: "well_stirred"
  formulation:
    particle_size_um: 10
    dose_mg: 100                          # initial dose for simulation
    route: "oral"
  scaling:
    methods: ["allometry", "ROE", "species_pbpk", "fu_corrected"]
    consensus: true
  dose_projection:
    noael_mg_kg: 50                       # from tox study (user input)
    noael_species: "rat"
    safety_factor: 10
    target_kd_nM: null                    # for MABEL (optional)
    target_ceff_nM: null                  # for PAD (optional)
  uncertainty:
    method: "lhs"
    n_samples: 500
    correlation: "iman_conover"
  population: null                        # Phase B

# Override: user can replace any predicted value with experimental
overrides:
  fu_p: {value: 0.21, source: "experimental", method: "equilibrium_dialysis"}
  # predicted value retained for comparison in report
```

### Configuration Manager

```python
class RunConfig:
    """Immutable snapshot of all settings for a pipeline run."""

    compound: CompoundConfig
    pipeline: PipelineConfig
    species_params: dict[str, SpeciesConfig]
    uncertainty: UncertaintyConfig
    population: Optional[PopulationConfig]  # Phase B
    model_versions: ModelVersions           # NEW: ML model hashes
    #   admet_ensemble: "sha256:abc123..."
    #   pka_model: "sha256:def456..."
    #   → ensures same models used for reproducibility

    def freeze(self) -> "RunConfig":
        """Deepcopy + freeze. No mutation after this point."""

    def to_yaml(self, path: str):
        """Serialize. This YAML alone reproduces the run."""

    def diff(self, other: "RunConfig") -> dict:
        """Diff two configs. For comparing runs."""

    def hash(self) -> str:
        """Deterministic hash for caching.
        NOTE: ODE solver results are NOT bitwise reproducible across
        platforms (floating point). 'Reproducible' here means
        'within numerical tolerance (rtol=1e-6)' not 'bit-identical'.
        Hash is for config identity, not output identity."""

    # Override propagation rules:
    # If user overrides fu_p with experimental value:
    #   → Kp recalculated with experimental fu_p
    #   → IVIVE recalculated with experimental fu_p
    #   → Uncertainty: override value treated as fixed (no CI)
    #     → reduces dose CI (one less uncertain parameter)
    #   → Report: notes "fu_p: experimental (overrides ML prediction)"
```

-----

## 6. File Structure (Final)

```
project_root/
├── pyproject.toml
├── README.md
├── ARCHITECTURE.md              # this document
├── src/
│   ├── core/
│   │   ├── __init__.py
│   │   ├── molecule.py          # SMILES → RDKit mol + descriptors
│   │   ├── schema.py            # Pydantic models for all data types
│   │   ├── compound_config.py   # YAML compound config I/O
│   │   ├── config_manager.py    # RunConfig: freeze, serialize, diff
│   │   ├── parameter_bridge.py  # Layer 1→2 conversions with audit
│   │   ├── liver_models.py      # Shared: well-stirred, parallel-tube, dispersion
│   │   ├── guardrails.py        # Input validation + applicability domain
│   │   ├── batch.py             # Multi-compound batch runner
│   │   └── units.py             # Unit conversion registry
│   │
│   ├── predict/                  # Layer 1: Property Prediction
│   │   ├── __init__.py
│   │   ├── admet_ensemble.py    # From Omega: XGBoost + ADMET-AI
│   │   ├── pka.py               # From Omega: pKa predictor
│   │   ├── fu_inc.py            # Microsomal unbound fraction
│   │   ├── bp_ratio.py          # Blood:plasma ratio
│   │   ├── renal.py             # Renal clearance estimation
│   │   └── conformal.py         # Conformal prediction intervals
│   │
│   ├── pbpk/                     # Layer 2: PBPK Simulation
│   │   ├── __init__.py
│   │   ├── topology.py          # From Sisyphus: compartment graph
│   │   ├── ode_compiler.py      # From Sisyphus: YAML → ODE system
│   │   ├── solver.py            # scipy solve_ivp (BDF)
│   │   ├── kp_calculator.py     # Rodgers & Rowland Kp
│   │   ├── formulation/
│   │   │   ├── dissolution.py   # Noyes-Whitney
│   │   │   ├── biopharmaceutics.py  # BCS classification
│   │   │   └── food_effect.py
│   │   └── species/
│   │       ├── human.yaml
│   │       ├── rat.yaml
│   │       ├── dog.yaml
│   │       └── monkey.yaml      # Phase B
│   │
│   ├── translational/           # Layer 3: Scaling & Dose
│   │   ├── __init__.py
│   │   ├── allometry.py         # Simple allometric scaling
│   │   ├── rule_of_exponents.py # ROE with auto exponent classification
│   │   ├── fu_corrected.py      # fu-corrected intercept method
│   │   ├── species_pbpk_scaling.py  # Bottom-up PBPK per species
│   │   ├── consensus_scaling.py # Multi-method consensus + confidence
│   │   ├── ivive.py             # In vitro-in vivo extrapolation
│   │   ├── hed.py               # NOAEL → HED (BSA scaling)
│   │   ├── mabel.py             # MABEL (PK-driven)
│   │   └── pad.py               # Pharmacologically active dose
│   │
│   ├── uncertainty/              # Layer 4: UQ
│   │   ├── __init__.py
│   │   ├── sampling.py          # LHS + Iman-Conover correlation
│   │   ├── propagation.py       # Sample → ODE → metrics pipeline
│   │   ├── sobol.py             # From Omega: Sobol sensitivity
│   │   └── dose_range.py        # Aggregate → dose CI + tornado plot data
│   │
│   ├── population/               # Layer 5: Population (Phase B)
│   │   ├── __init__.py
│   │   ├── virtual_pop.py       # Demographic + phenotype sampling
│   │   ├── phenotype_db/
│   │   │   ├── cyp2d6.yaml      # Allele frequencies by ethnicity
│   │   │   ├── cyp2c9.yaml
│   │   │   ├── cyp2c19.yaml
│   │   │   ├── cyp3a4.yaml
│   │   │   └── cyp1a2.yaml
│   │   ├── demographics.py      # Weight, height, age distributions
│   │   ├── covariance.py        # Inter-parameter correlations
│   │   ├── ontogeny.py          # CYP maturation (pediatric)
│   │   └── trial_design/
│   │       ├── dose_escalation.py  # 3+3, BOIN, mTPI
│   │       ├── ddi.py              # DDI simulation
│   │       └── special_pops.py     # Renal/hepatic/pediatric adjustments
│   │
│   ├── report/                    # Cross-cutting: Report Generator
│   │   ├── __init__.py
│   │   ├── collector.py          # Observer: collects layer outputs
│   │   ├── narrative.py          # Template-based text generation
│   │   ├── figures.py            # Matplotlib/Plotly figure generation
│   │   ├── templates/
│   │   │   ├── fih_dose_rationale.yaml
│   │   │   ├── clinical_pharm_overview.yaml
│   │   │   └── nonclinical_pk_summary.yaml
│   │   └── export.py             # → DOCX, PDF, Markdown
│   │
│   ├── dashboard/                 # Layer 6: Decision Support (Phase C)
│   │   ├── __init__.py
│   │   ├── pareto.py             # Multi-objective Pareto frontier
│   │   ├── risk_flags.py         # Automated safety scoring
│   │   └── comparison.py         # Multi-compound ranking
│   │
│   ├── cli/                       # Phase A: CLI interface
│   │   ├── __init__.py
│   │   └── main.py               # click or typer CLI
│   │   # Commands:
│   │   #   predict <smiles>               → Layer 1 only, print ADMET table
│   │   #   simulate <smiles> [--route oral|iv] [--dose 100mg] → Layer 1+2, Cp-time
│   │   #   translate <smiles> --noael 50 --species rat → Layer 1+2+3, FIH dose
│   │   #   recommend <smiles> --noael 50 [--uncertainty] → Full pipeline + UQ
│   │   #   report <smiles> ... --output report.md → Pipeline + report
│   │   #   batch <csv> → Multi-compound, CSV with SMILES column
│   │   #   validate --benchmark tier1 → Run validation suite
│   │
│   │   # Python API (programmatic usage):
│   │   # from charon import Pipeline
│   │   # result = Pipeline("CCO", route="oral", dose_mg=100).run()
│   │   # result.pk_parameters.cmax
│   │   # result.dose_recommendation.point_estimate_mg
│   │   # result.report.to_markdown("output.md")
│   │
│   ├── api/                       # FastAPI backend (Phase B)
│   │   ├── __init__.py
│   │   ├── main.py
│   │   └── routes/
│   │       ├── predict.py
│   │       ├── simulate.py
│   │       ├── translate.py
│   │       └── report.py
│   │
│   └── ui/                        # React/Vite frontend (Phase B/C)
│       └── ...
│
├── data/
│   ├── training/                  # ML model training data (ADMET datasets)
│   └── validation/                # Benchmark compounds (raw literature data)
│
├── validation/
│   ├── data/
│   │   ├── tier1_obach/           # Obach et al. datasets: in vitro + observed human CL
│   │   │                          # N≥40 compounds. Layer 2 CL/Vss validation.
│   │   ├── tier2_drugs_at_fda/    # Clinical pharm reviews with preclinical + FIH dose
│   │   │                          # N≥15 (realistic). Full pipeline validation.
│   │   ├── tier3_literature/      # Manual curation to expand coverage
│   │   └── pk_db/                 # PK-DB public data (supplementary)
│   ├── benchmarks/
│   │   ├── layer1_admet.py
│   │   ├── layer2_human_pk.py
│   │   ├── layer3_fih_dose.py
│   │   └── metrics.py            # AAFE, fold-error, %within_Xfold
│   └── targets.yaml
│       # Layer 1 ADMET: AAFE < 2.0 per parameter
│       # Layer 2 human PK: AAFE < 2.5 for CL, < 3.0 for Vss
│       # Layer 3 FIH dose: within 3-fold for ≥60% of compounds
│       #   (realistic given Layer 2 uncertainty propagation;
│       #    70% achievable only with experimental input overrides)
│       # Layer 5 pop PK: 90% PI covers observed in ≥80%
│
└── tests/
    ├── unit/
    │   ├── test_parameter_bridge.py  # HIGHEST PRIORITY
    │   ├── test_guardrails.py
    │   ├── test_kp_calculator.py
    │   ├── test_hed.py
    │   └── ...
    ├── integration/
    │   ├── test_smiles_to_pk.py      # End-to-end Layer 0→2
    │   ├── test_smiles_to_dose.py    # End-to-end Layer 0→3
    │   └── test_full_pipeline.py     # Layer 0→4
    └── regression/
        ├── test_known_drugs.py       # Aspirin, midazolam, etc.
        └── golden_outputs/           # Frozen expected outputs
```

-----

## 7. Implementation Plan

### Phase A Sprints

```
Sprint 1 (Weeks 1-2): Foundation
  ├── core/schema.py — Pydantic models for all data types
  ├── core/config_manager.py — RunConfig with freeze/serialize
  ├── core/parameter_bridge.py — ALL IVIVE conversions
  ├── core/liver_models.py — shared well-stirred/parallel-tube/dispersion
  ├── core/guardrails.py — input validation
  ├── core/units.py — unit registry
  └── tests/unit/test_parameter_bridge.py — exhaustive tests
  Deliverable: parameter_bridge converts CLint→CLh correctly
  for 5 known drugs (midazolam, warfarin, diclofenac,
  omeprazole, dextromethorphan) with literature comparison.

Sprint 2 (Weeks 3-5): Layer 1 Completion
  ├── predict/fu_inc.py
  ├── predict/bp_ratio.py
  ├── predict/renal.py
  ├── Import Omega ADMET ensemble
  │   ★ Risk: extracting Omega's ensemble as standalone module
  │   may require significant refactoring. Budget extra week.
  ├── Wire Layer 1 → parameter_bridge → compound config YAML
  └── Validate: 30 compounds, predicted vs literature ADMET
  Deliverable: SMILES → complete compound config YAML, auto-generated.

Sprint 3 (Weeks 6-10): Layer 2 Human PBPK
  ├── Import Sisyphus ODE topology + compiler
  │   ★ Risk: Sisyphus compiler extraction is complex. Budget 2 weeks.
  ├── pbpk/kp_calculator.py — Rodgers & Rowland
  ├── species/human.yaml — complete parameter set
  ├── Fix Omega bugs: fu double-application, Fg CYP3A4
  ├── formulation/dissolution.py + biopharmaceutics.py
  ├── Route handling: oral + IV bolus + IV infusion
  ├── ★ ODE solver benchmark: measure wall time per solve
  │   Target: <200ms. If >500ms → Numba JIT or solver refactor
  └── End-to-end: SMILES → Cp-time → PK parameters
  Deliverable: predicted human PK for midazolam within
  2-fold of observed. Cp-time profile plotted.

Sprint 4 (Weeks 11-13): Layer 3 Translational
  ├── species/rat.yaml + dog.yaml — physiological parameters
  ├── translational/allometry.py
  ├── translational/rule_of_exponents.py
  ├── translational/ivive.py (species-specific IVIVE)
  ├── translational/consensus_scaling.py
  ├── translational/hed.py + mabel.py + pad.py
  └── End-to-end: SMILES + NOAEL → FIH dose recommendation
  Deliverable: FIH dose for 10 compounds, compared to actual.

Sprint 5 (Weeks 14-16): Layer 4 Uncertainty
  ├── uncertainty/sampling.py — LHS + Iman-Conover
  ├── uncertainty/propagation.py — sample → ODE → metrics
  ├── uncertainty/sobol.py — import from Omega
  ├── uncertainty/dose_range.py — aggregate to CI
  └── End-to-end: SMILES → dose [point, 90% CI]
  Deliverable: dose recommendation with CI and tornado plot.

Sprint 6 (Weeks 17-20): Validation & Report
  ├── Tier 1 validation: Obach dataset (N≥40) for CL/Vss
  ├── Tier 2 validation: Drugs@FDA subset (N≥15) for FIH dose
  ├── validation/benchmarks/ — all benchmark scripts
  ├── report/collector.py + narrative.py + figures.py + export.py
  ├── report/templates/fih_dose_rationale.yaml
  ├── Batch runner: core/batch.py
  ├── CLI interface: cli/main.py
  └── Full pipeline validation against targets
  Deliverable: auto-generated FIH dose rationale (Markdown)
  for any input SMILES. Benchmark results published.

Phase A complete (~20 weeks). → Decision point for Phase B.
```

### Phase B Sprints (estimate: 10-12 weeks)

```
Sprint B1: Virtual population generator
Sprint B2: Population PBPK (+ solver optimization: JAX or Numba)
Sprint B3: DDI simulation module
Sprint B4: Special populations (renal, hepatic, pediatric)
Sprint B5: Dose escalation simulation (3+3, BOIN)
Sprint B6: Population validation + report templates
```

### Phase C Sprints (estimate: 6-8 weeks)

```
Sprint C1: Multi-compound batch comparison
Sprint C2: Pareto frontier visualization
Sprint C3: Automated risk scoring
Sprint C4: Dashboard UI
```

-----

## 8. Technical Decisions

|Decision                     |Choice                                                            |Rationale                                                       |
|-----------------------------|------------------------------------------------------------------|----------------------------------------------------------------|
|Language                     |Python                                                            |Omega/Sisyphus already Python; ecosystem (RDKit, scipy, sklearn)|
|ODE solver                   |scipy solve_ivp (BDF) → JAX/diffrax in Phase B                    |BDF handles stiff PBPK; JAX for batch parallelism later         |
|ML framework                 |XGBoost + ADMET-AI ensemble                                       |From Omega; proven, lightweight                                 |
|Config format                |YAML                                                              |From Sisyphus; human-readable, Git-diffable                     |
|Data validation              |Pydantic v2                                                       |Type safety, auto-serialization                                 |
|API                          |FastAPI                                                           |From Omega; async, auto-docs                                    |
|UI                           |React + Vite                                                      |From Omega; existing codebase                                   |
|Testing                      |pytest                                                            |Standard; aim >80% coverage on core/ and translational/         |
|CI/CD                        |GitHub Actions                                                    |From Omega                                                      |
|Packaging                    |pyproject.toml + src layout                                       |From Omega; pip installable                                     |
|Report export                |python-docx (DOCX), fpdf2 (PDF primary), WeasyPrint (PDF optional)|Lightweight, no LaTeX dependency; WeasyPrint needs cairo/pango  |
|Kp method                    |Rodgers & Rowland 2005/2006                                       |Industry standard; covers all ionization classes                |
|Population sampling (Phase B)|Monte Carlo with covariance                                       |Standard in PBPK; Simcyp does same                              |
|Uncertainty sampling         |LHS + Iman-Conover                                                |Best balance of efficiency vs implementation complexity         |

-----

## 9. Risk Registry

|Risk                                        |Likelihood|Impact        |Mitigation                                                                                                     |
|--------------------------------------------|----------|--------------|---------------------------------------------------------------------------------------------------------------|
|Validation data scarcity (FIH dose)         |HIGH      |HIGH          |Tiered: Obach N≥40 for CL (Tier 1), Drugs@FDA N≥15 for FIH dose (Tier 2), literature expansion (Tier 3)        |
|Fg prediction inaccuracy (CYP3A4 substrates)|HIGH      |HIGH          |Prioritize Fg fix (Sprint 3); sensitivity analysis to quantify                                                 |
|ISEF uncertainty propagation                |MEDIUM    |HIGH          |Use literature ranges; flag when ISEF drives >30% of dose variance                                             |
|Species YAML curation effort                |HIGH      |MEDIUM        |Leverage PK-Sim OSP open-source parameters; rat+dog Phase A, monkey Phase B                                    |
|Conformal prediction miscalibration OOD     |MEDIUM    |MEDIUM        |Applicability domain check in Layer 0; recalibrate per chemical class                                          |
|Solver performance at population scale      |MEDIUM    |HIGH (Phase B)|Benchmark in Sprint 3; if >500ms, Numba JIT immediate. JAX migration Phase B.                                  |
|Scope creep (metabolite PBPK, etc.)         |HIGH      |MEDIUM        |Strict phase gating; metabolite = Phase C at earliest                                                          |
|Regulatory acceptance of ML predictions     |MEDIUM    |HIGH          |Frame as hypothesis-generating, not regulatory submission; experimental validation always recommended in report|

-----

## 10. Competitive Landscape & Differentiation

|Tool                        |Open Source|ML-augmented|Structure→Dose              |Cost        |
|----------------------------|-----------|------------|----------------------------|------------|
|Simcyp (Certara)            |No         |No          |Partial (manual)            |~$50-150K/yr|
|GastroPlus (Simulation Plus)|No         |No          |Partial (manual)            |~$40-100K/yr|
|PK-Sim/MoBi (OSP)           |Yes        |No          |No                          |Free        |
|pksensi (R)                 |Yes        |No          |No (sensitivity only)       |Free        |
|PKPDsim (R)                 |Yes        |No          |No (pop-PK, not PBPK)       |Free        |
|Omega (ours)                |Yes        |Yes         |Partial (no dose projection)|Free        |
|**Charon (this project)**   |**Yes**    |**Yes**     |**Yes (end-to-end)**        |**Free**    |

Note: NONMEM and Monolix are pop-PK/PD tools, not PBPK. Not direct
competitors but relevant in the broader pharmacometrics ecosystem.
Distinction: this project predicts PK from structure (bottom-up);
NONMEM/Monolix estimate PK from observed data (top-down).

**Unique value proposition:**

1. Only open-source tool with ML-predicted ADMET → PBPK → FIH dose, end-to-end
1. Uncertainty quantification native throughout (not bolted on)
1. Auto-generated regulatory-ready dose rationale documents
1. Python-native (vs PK-Sim C#, GastroPlus GUI-only)
1. Designed for small biotech/academic groups who can't afford Simcyp

-----

*Document version: v4.0 (converged — 10 consecutive clean passes achieved)*
*Total: 3 design iterations + 60-pass self-review cycle (48 fixes applied, passes 51-60 clean)*
*Passes 1-10: structural (allometry, compartments, sprint scope)*
*Passes 11-20: pharmacological (MABEL formula, IV route, IVIVE duplication, CLI spec)*
*Passes 21-30: numerical (audit trail, units, YAML syntax, salt form, validation targets)*
*Passes 31-50: precision (calculation errors, IVIVC→IVIVE terminology, data directory dedup, cross-references)*
*Passes 51-60: 10 consecutive clean — convergence confirmed*
