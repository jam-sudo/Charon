# Charon

Charon is an open-source, auditable Python platform that takes a molecular structure (SMILES) and produces a First-in-Human (FIH) dose recommendation with uncertainty quantification. It builds on the Omega ADMET ensemble and Sisyphus PBPK engine. Phase A covers single-subject, single-compound predictions from structure alone.

---

## Project status

- **Phase A**: MVP feature-complete (Layers 0--4, CLI, reports)
- **948 tests** across unit, integration, and regression suites
- **Layer 3 §8 target met** at 8/12 = 66.7% within-3x (Sprint 12 onward); Layer 1/2 accuracy targets NOT met -- see [Known limitations](#known-limitations) and [Validation status](#validation-status)
- **Phase B** (population variability, DDI, API) and **Phase C** (multi-compound dashboard) are scaffolded only -- no implementation

---

## Known limitations

This section is the most important part of this README. Read it before using any output from Charon.

- **Not FDA-cleared. Research and educational tool only.** Any FIH dose selection for real clinical trials must go through a qualified pharmacometrician with validated commercial software.
- **Accuracy vs targets**: AAFE_CL = 3.99, AAFE_Vss = 3.03 on the Obach n=12 panel (targets: 2.5 and 3.0 respectively). See [Layer 2 report](validation/reports/layer2_human_pk.md).
- **CLint Tier 2 ML limit**: scaffold-CV AAFE ~2.5. Experimental CLint override is strongly recommended for any serious use.
- **Layer 3 FIH dose validation (Sprint 17)**: 8/12 gold within 3-fold (§8 target ≥60% MET), 12/12 within 10-fold, 12/12 sanity floor pass. Computed via PAD path (target_ceff_nM, safety_factor=10) on oral routes (Sprint 11 migrated all Tier A from iv_bolus to oral). Four close-but-not-quite residuals remain (lisinopril 3.085x, diclofenac 3.096x, propranolol 4.819x, diazepam 4.910x) — each requires architectural lift (PEPT1, UGT/CYP-specific ML, extended-clearance) rather than parameter tuning. See [Validation status § Layer 3](#layer-3----fih-dose-sprint-17-n12).
- **No CYP phenotype modeling**: CYP2C19 polymorphisms (omeprazole) and CYP2D6 polymorphisms (dextromethorphan) are not modeled. Compounds primarily cleared by polymorphic CYPs will have higher prediction error.
- **Conformal CI is marginal, not conditional**: out-of-domain (OOD) compounds may have actual coverage well below the nominal 90%.
- **Single-subject only**: no population variability, no virtual trials, no special populations (Phase B scope).
- **PBPK ODE reproducibility**: numerical tolerance rtol=1e-6. Results are reproducible within that tolerance, not bit-identical across platforms.
- **Rodgers and Rowland Kp overpredicts adipose** for lipophilic weak bases. An empirical override path exists but the fundamental R&R limitation remains.
- **fu_p near 1.0 LHS edge case**: Latin Hypercube Sampling can generate fu_p > 1.0 for near-unbound compounds; these are clipped but may distort uncertainty propagation.
- **No DDI, metabolite tracking, or active transport** beyond a P-gp permeability proxy.
- **Salt form correction** is the user's responsibility. Charon does not infer salt form from SMILES.

---

## Install

Python 3.11+ required. Core dependencies: rdkit, pydantic>=2, scipy, numpy, xgboost, scikit-learn, pyyaml.

```bash
pip install -e .          # minimal install
pip install -e ".[dev]"   # with dev/test dependencies
```

---

## Quickstart

### CLI

Charon provides five subcommands. Each builds on the previous layer.

```bash
# Layer 0 + Layer 1: ADMET property prediction
charon predict "CCO"

# Through Layer 2: PBPK simulation
charon simulate "CCO" --route oral --dose 100

# Through Layer 3: FIH dose projection
charon translate "Cc1ncc2n1-c1ccc(Cl)cc1C(c1ccccc1F)=NC2" \
    --route iv_bolus --dose 5 --noael 2 --noael-species rat

# Through Layer 4: dose projection + uncertainty quantification
charon recommend "Cc1ncc2n1-c1ccc(Cl)cc1C(c1ccccc1F)=NC2" \
    --route iv_bolus --dose 5 --noael 2 --noael-species rat --uncertainty

# Full pipeline + Markdown/JSON report
charon report "Cc1ncc2n1-c1ccc(Cl)cc1C(c1ccccc1F)=NC2" \
    --route iv_bolus --dose 5 --noael 2 --noael-species rat --output midazolam
```

### Python API

```python
from charon import Pipeline

result = Pipeline.from_smiles(
    "Cc1ncc2n1-c1ccc(Cl)cc1C(c1ccccc1F)=NC2",
    route="iv_bolus",
    dose_mg=5.0,
).run()

print(f"Cmax = {result.pk_parameters.cmax:.2f} ug/L")
```

---

## Reference drug quick-check

We maintain a regression test panel of known drugs (Obach 1999 tier-1 compounds). Run it with:

```bash
pytest tests/regression/test_known_drugs.py -v
```

---

## Validation status

All numbers below are from automated benchmarks committed to this repository. We report failures honestly -- Phase A does not meet its accuracy targets.

### Layer 1 -- ADMET (adme_reference.csv, n=153)

| Property | Source | n | Metric | Value | Target | Gate |
|---|---|---|---|---|---|---|
| logP | RDKit Crippen | 151 | MAE | 0.970 | < 1.0 | [PASS] |
| fu_p | XGBoost ensemble | 151 | AAFE | 2.71 | < 2.0 | [FAIL] |
| bp_ratio | Empirical formula | 151 | AAFE | 1.35 | < 2.0 | [PASS] |
| peff_cm_s | Derived | 0 | -- | not evaluated | -- | -- |
| CLint | -- | -- | -- | excluded (unit mismatch) | -- | -- |

Full results: [Layer 1 report](validation/reports/layer1_admet.md)

CLint was excluded from the Layer 1 benchmark because Charon predicts in hepatocyte units (uL/min/10^6 cells) while the reference set uses recombinant units (uL/min/pmol CYP3A4). This is not a workaround -- the units are fundamentally different assay systems.

### Layer 2 -- Human PK (Obach 1999 tier-1, n=12, with R&R Kp overrides)

| Metric | AAFE | within 2-fold | within 3-fold | Target AAFE | Gate |
|---|---|---|---|---|---|
| CL (L/h) | 3.99 | 33% | 42% | < 2.5 | [FAIL] |
| Vss (L) | 3.03 | 42% | 50% | < 3.0 | [FAIL] |
| t_half (h) | 9.39 | 17% | 17% | -- | -- |

Full results: [Layer 2 report](validation/reports/layer2_human_pk.md)

### Layer 3 -- FIH dose (Sprint 17, n=12)

Two-tier benchmark against `validation/data/fih_reference/panel.yaml`. All Tier A compounds run as `oral` (Sprint 11 migrated from iv_bolus); reference doses are FDA-label starting doses for the same route, eliminating the Sprint 9 `1/F` bias.

| Tier | Metric | Result | §8 target |
| --- | --- | --- | --- |
| Gold (n=12) | Within-3-fold of reference FIH dose | **8/12 (66.7%)** | >= 60% **[PASS]** |
| Gold (n=12) | Within-10-fold of reference FIH dose | 12/12 (100%) | -- |
| Sanity (n=12) | MRSD <= approved starting dose | 12/12 (100%) | gated [PASS] |

Panel composition: 5 core (Sprint 7) + 4 Obach promotions (theophylline, diclofenac, diazepam, metoprolol) + 3 elimination-diversity additions (acetaminophen UGT, lisinopril renal, atorvastatin CYP3A4+OATP).

**§8 target MET** at 8/12 = 66.7%, achieved through Sprints 10-17:
- Sprint 10: IVIVE-bias diagnostic decomposition (`docs/superpowers/sprint10-ivive-bias-ticket.md`).
- Sprint 11: oral route migration → +2 within-3x (5/12 → 7/12); route_bias collapsed 42% → 0%.
- Sprint 12: atorvastatin OATP1B1 enhancement (`hepatic_clint_multiplier: 8.0`) → +1 within-3x (7/12 → 8/12, **§8 first met**).
- Sprints 13-17: close-but-not-quite refinements for the 4 remaining compounds. All are now within 5x; full closure requires architectural sprints (per-CYP ML retraining, PEPT1 transporter, extended-clearance).

**Tier A passes (within-3x, n=8):** warfarin (1.02x), theophylline (1.02x), acetaminophen (1.21x), midazolam (1.46x), omeprazole (1.60x), atorvastatin (1.70x), metoprolol (2.13x), verapamil (2.49x).

**Tier A close-but-not-quite (n=4, all within 5x):**
- lisinopril 3.085x — Sprint 17 Branch B (Peff back-calibration; PEPT1 transporter not modeled in ACAT)
- diclofenac 3.096x — Sprint 13 Branch B (UGT/CYP2C9 multiplier 3.5, literature midpoint)
- propranolol 4.819x — Sprint 15 Branch B (CYP2D6 IVIVE multiplier 6.0, literature midpoint)
- diazepam 4.910x — Sprint 14 honest null (low fu_p well-stirred sensitivity, framework-limited)

MRSD computed via PAD path with `safety_factor=10`. CLAUDE.md §6.5 honesty discipline: no multiplier inflated above cited literature range to force closure.

Full results: [Layer 3 report](validation/reports/layer3_fih_dose.md). Per-compound IVIVE decomposition: [Layer 3 IVIVE decomposition](validation/reports/layer3_ivive_decomposition.md).

### Uncertainty

Every `fu_p` and `clint_hepatocyte` prediction ships with a 90% log-space conformal confidence interval by default (calibrated against `adme_reference.csv` n=153 for `fu_p` and scaffold-CV OOF residuals n=1441 for `clint_hepatocyte`). Empirical in-sample coverage: both ~0.900. See [Layer 1 report § Conformal coverage](validation/reports/layer1_admet.md). CIs flow through `predict_properties` -> `Pipeline.from_smiles` -> `charon report` automatically; pass `CONFORMAL_OFF` to opt out.

### CLint Tier 3 and applicability domain

Charon's CLint prediction uses a three-tier strategy (CLAUDE.md §6j):

- **Tier 1** -- experimental override from compound YAML. Always preferred when available.
- **Tier 2** -- XGBoost regression (scaffold-CV AAFE 2.49). Ships with a 90% log-space conformal interval.
- **Tier 3** -- 3-class classifier (Low `[0.1, 10)` / Med `[10, 50)` / High `[50, 1000]` uL/min/10^6 cells, scaffold-CV Macro F1 = 0.54) for compounds whose CLint-local applicability domain is LOW (max Tanimoto < 0.3 vs 1476 training compounds).

A **CRITICAL warnings** section appears in the report whenever Tier 3 fires. The classifier emits a categorical distribution `{low, med, high}` that Layer 4 uncertainty propagates as categorical x log-uniform-within-bucket via an inverse-CDF that preserves Latin-Hypercube stratification. The point estimate is the log-geometric-mean center of the most-likely bucket (Low=1.0, Med=22.36, High=223.6 uL/min/10^6 cells); the CI is the bucket range.

**AD=LOW prevalence is expected to be high on structurally diverse compound sets.** On the 153-compound `adme_reference.csv`, 50% fire as LOW; on the Obach-12 panel, 4 (dextromethorphan, omeprazole, theophylline, verapamil) do. This is by design -- the CLint training distribution is narrower than general drug-like space, and Tier 3 exists precisely to surface this honestly rather than pretend the Tier 2 regression is trustworthy on novel scaffolds. For any Tier 3 output, obtain an experimental CLint measurement before any dose decision.

Pass `force_tier3=True` to `predict_properties` to route a compound through Tier 3 regardless of AD (useful for testing).

---

## Architecture overview

Charon is a 6-layer pipeline. Each layer is independently testable and chains into the next.

```
SMILES
  |
  v
[Layer 0]  Input Validation + Guardrails
  |          -> ValidatedMolecule + reliability flags
  v
[Layer 1]  Property Prediction (ADMET ensemble + conformal CI)
  |          -> logP, pKa, fu_p, CLint, Papp, hERG, ...
  v
[ParameterBridge]  Layer 1 -> 2 Translation   ** most error-prone junction **
  |          -> CLh (IVIVE), Peff, CLrenal, Kp per tissue
  |          -> every conversion logged in ConversionLog
  v
[Layer 2]  PBPK Simulation (ODE, ~50 state variables, BDF solver)
  |          -> Cp-time profile, PK parameters (Cmax, AUC, t1/2, CL)
  v
[Layer 3]  Translational Scaling (allometry, consensus, HED/MABEL/PAD)
  |          -> FIH dose recommendation
  v
[Layer 4]  Uncertainty Quantification (LHS, Sobol, dose CI)
  |          -> dose [point estimate, 90% CI] + sensitivity analysis
  v
[Report]   Auto-generated FIH dose rationale document
```

The ParameterBridge is the single most critical junction. Every IVIVE conversion is logged with full audit trail (ConversionLog) to prevent silent unit errors. See [ARCHITECTURE.md](ARCHITECTURE.md) for the full design document.

### Package layout

```
src/charon/
  core/            Foundation (schema, units, parameter_bridge, liver_models)
  predict/         Layer 1 (ADMET prediction)
  pbpk/            Layer 2 (PBPK simulation)
  translational/   Layer 3 (allometry, HED, consensus scaling)
  uncertainty/     Layer 4 (LHS, Sobol, dose CI)
  report/          Report generation
  cli/             CLI interface
  population/      Phase B scaffold (empty)
  dashboard/       Phase C scaffold (empty)
  api/             Phase B scaffold (empty)
```

---

## Development

```bash
# Full test suite (~948 tests)
pytest tests/

# Unit tests with coverage
pytest tests/unit/ -v --cov=src/charon/core --cov-report=term-missing

# Layer 1 ADMET benchmark (generates validation/reports/layer1_admet.md)
python3 validation/benchmarks/layer1_admet.py

# Layer 2 human PK benchmark (generates validation/reports/layer2_human_pk.md)
python3 validation/benchmarks/layer2_human_pk.py

# Regenerate regression baselines
UPDATE_BASELINES=1 pytest tests/regression/test_known_drugs.py
```

---

## Roadmap

- **Phase B**: population variability, virtual trial, drug-drug interactions, special populations, REST API
- **Phase C**: multi-compound dashboard, Pareto optimization

---

## Citing Charon

```
Charon v0.1.0, https://github.com/jam-sudo/Charon
```

---

## License

MIT. See [LICENSE](LICENSE).

---

## Acknowledgments

Omega ADMET ensemble, Sisyphus topology engine, Obach 1999 reference dataset, RDKit, scipy.
