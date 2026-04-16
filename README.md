# Charon

Charon is an open-source, auditable Python platform that takes a molecular structure (SMILES) and produces a First-in-Human (FIH) dose recommendation with uncertainty quantification. It builds on the Omega ADMET ensemble and Sisyphus PBPK engine. Phase A covers single-subject, single-compound predictions from structure alone.

---

## Project status

- **Phase A**: MVP feature-complete (Layers 0--4, CLI, reports)
- **790+ tests** across unit, integration, and regression suites
- **Phase A accuracy targets are NOT met** -- see [Known limitations](#known-limitations) and [Validation status](#validation-status)
- **Phase B** (population variability, DDI, API) and **Phase C** (multi-compound dashboard) are scaffolded only -- no implementation

---

## Known limitations

This section is the most important part of this README. Read it before using any output from Charon.

- **Not FDA-cleared. Research and educational tool only.** Any FIH dose selection for real clinical trials must go through a qualified pharmacometrician with validated commercial software.
- **Accuracy vs targets**: AAFE_CL = 3.99, AAFE_Vss = 3.03 on the Obach n=12 panel (targets: 2.5 and 3.0 respectively). See [Layer 2 report](validation/reports/layer2_human_pk.md).
- **CLint Tier 2 ML limit**: scaffold-CV AAFE ~2.5. Experimental CLint override is strongly recommended for any serious use.
- **Layer 3 FIH dose validation deferred**: no public tier2_drugs_at_fda dataset has been built. FIH dose accuracy is unvalidated.
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

### Layer 3 -- FIH dose

**DEFERRED.** No public tier2_drugs_at_fda dataset has been built. FIH dose accuracy has not been validated against observed clinical starting doses.

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
# Full test suite (~790 tests)
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
