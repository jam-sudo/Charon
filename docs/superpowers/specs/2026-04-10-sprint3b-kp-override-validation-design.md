# Sprint 3b Session 1 — Honest IV Kernel (Kp Override + Validation Panel) — Design

**Date:** 2026-04-10
**Scope:** Charon Sprint 3b, Session 1 of N — R&R Kp override path + Obach 1999 Layer 2 validation panel
**Status:** Design — awaiting approval
**Supersedes:** N/A
**Builds on:** `2026-04-10-sprint3-pbpk-iv-design.md` (Sprint 3a, IV kernel)

---

## 1. Session Goal

Sprint 3a shipped a working IV PBPK kernel, but only one compound
(theophylline) passes a strict gate. The R&R Kp model is known to over-predict
adipose Kp for lipophilic weak bases (midazolam Vss ≈ 6× observed), so
"the IV kernel works" is presently a claim, not a measurement.

This session turns the Sprint 3a kernel into an **honest Layer 2 PBPK
engine** by:

1. Adding a schema-level **empirical Kp override** path so tissue-level
   R&R predictions can be replaced with cited literature/experimental
   values without touching the Kp model.
2. Curating a **10-compound validation panel** from Obach 1999 covering
   neutral / acid / weak-base / strong-base chemistry, loaded from YAML
   rather than hand-constructed Python fixtures.
3. Extending the `layer2_human_pk.py` benchmark to compute and report
   **panel-level AAFE and fold-error metrics** for CL, Vss, t½ with a
   side-by-side "override vs no-override" comparison, so the structural
   limitation of R&R is quantified and left in the run log.

No Kp model is changed. No new route (oral/ACAT) is added. No new species
is validated. Every downstream sprint gains a trustworthy measurement
surface before new complexity is layered on.

---

## 2. Definition-of-Done

This session is **measurement + reporting mode**, not a hard AAFE gate
(per user decision 2026-04-10). The benchmark reports AAFE across the
panel; the CI/test gates check only the regression invariants below.

**Must hold on final commit:**

1. All 477 Sprint 3a tests still pass. No pre-existing test is deleted or
   relaxed.
2. theophylline's existing 2-fold strict gate in `test_pipeline.py::
   TestPipelineTheophyllineValidation` still passes unchanged.
3. midazolam's `TestPipelineMidazolamLimitation` still passes unchanged
   (non-gated, documents R&R limitation). A **new** testcase verifies
   the override path end-to-end using a **mechanism check**, NOT a
   specific numerical target:

   - `Pipeline(midazolam, ...).run()` with a hand-picked synthetic
     `empirical_kp_by_tissue={'adipose': {value: 10.0, ...}}` must
     return a `PipelineResult.metadata['kp_overrides']` list of length
     1 whose single entry has `tissue == "adipose"` and
     `empirical_value == 10.0`.
   - The resulting `Vss` value must be finite and strictly less than
     the no-override `Vss` from the same compound (override reduces
     overprediction; direction-of-effect check).

   This test uses synthetic override values so it is independent of
   the citation verification outcome. It proves the plumbing works.
   The actual literature-sourced override that lands in the midazolam
   YAML is tested separately at the benchmark level (DoD §4).
4. `validation/benchmarks/layer2_human_pk.py` loads
   `validation/data/tier1_obach/panel.yaml`, runs all 10 compounds to
   completion without solver failure, prints a two-pass summary table
   (no-override, with-override), and returns exit 0 unless a
   strict-targets compound fails.
5. Panel-level AAFE values for CL, Vss, t½ are printed on stdout and
   included in pytest `record_property` metadata for the integration
   smoke test, so CI artifacts preserve them.
6. Coverage on `src/charon/pbpk/ode_compiler.py` remains ≥ 85% (Sprint 3a
   baseline). New Kp override paths are covered.
7. **Session post-condition** (sanity floor, not a hard gate): if
   panel-level `AAFE_Vss` in *with-override* mode is > 5.0, the session
   is not declared complete. A > 5.0 AAFE indicates either a curation
   error or a deeper engine bug; root-cause before declaring done.

**Explicitly NOT in DoD:**

- No target like "AAFE_CL < 2.5". We are *measuring*, not *hitting* the
  ARCHITECTURE target this session. Stating that number as the measured
  panel value in the session summary is the deliverable.
- No requirement that R&R + override beat R&R alone on every metric. If
  the override injection is small (one tissue on one compound), the
  panel metric movement will be small.

**Acceptance of residual scientific gap:**
This session may conclude with a measured `AAFE_Vss` that is still
above the ARCHITECTURE target of 3.0 (the session only enforces a
5.0 sanity floor to catch catastrophic breakage). If that happens, the
session is still complete — but the *next* session must explicitly
decide whether to pursue Kp model work (Berezhkovskiy defaults revisit,
Poulin-Theil adipose branch, new Kp methods) or to accept the
limitation as documented and move on to ACAT. That decision is made
with a concrete measurement in hand, which is exactly why this session
exists.

**Override demonstration target (minimum ≥2 compounds):**
The single-compound override demo (midazolam only) is too thin — if
its citation verification fails and falls to Fallback 3
(`source: synthetic_test`), the with-override pass would be
indistinguishable from no-override and the session loses its
demonstration value. Therefore the DoD requires **at least two
compounds with verified-citation empirical Kp overrides** before
declaring the session done. Primary: midazolam (adipose, from
Björkman PBPK literature). Secondary: propranolol or diazepam
(tissue to be selected during curation based on which compound and
tissue combination has clean literature coverage — primary candidate
sources: Rodgers 2005 Table 4, Poulin & Theil 2002). If fewer than 2 verified overrides can be
committed, the session is incomplete and the implementation plan
must explicitly escalate (do not silently degrade to Fallback 3).

---

## 3. Architecture

### 3.1 Schema additions (`src/charon/core/schema.py`)

Three additive changes. All are backward-compatible: existing YAML
fixtures and Python test factories continue to work because every new
field has a default.

```python
# NEW: Distribution property container
class DistributionProperties(BaseModel):
    """Distribution-related compound properties.

    Currently holds only the empirical Kp override. Forward-compatible
    for Vss_pred, tissue fu_tissue, binding capacity etc.
    """

    empirical_kp_by_tissue: dict[str, PredictedProperty] | None = None

    @field_validator("empirical_kp_by_tissue")
    @classmethod
    def _validate_kp_values(cls, v):
        if v is None:
            return v
        if not v:
            raise ValueError(
                "empirical_kp_by_tissue must be None or a non-empty dict"
            )
        for tissue, p in v.items():
            if p.value <= 0 or p.value > 200:
                raise ValueError(
                    f"empirical_kp_by_tissue[{tissue!r}] = {p.value} "
                    f"outside physiological range (0, 200]"
                )
        return v


# ADD field to existing PhysicochemicalProperties
class PhysicochemicalProperties(BaseModel):
    logp: PredictedProperty | None = None
    pka_acid: PredictedProperty | None = None
    pka_base: PredictedProperty | None = None
    solubility_ug_ml: PredictedProperty | None = None
    compound_type: Literal[
        "neutral", "acid", "base", "zwitterion"
    ] | None = None  # NEW — enables YAML-level classification override


# ADD field to existing CompoundProperties
class CompoundProperties(BaseModel):
    physicochemical: PhysicochemicalProperties = PhysicochemicalProperties()
    permeability: PermeabilityProperties = PermeabilityProperties()
    binding: BindingProperties = BindingProperties()
    metabolism: MetabolismProperties = MetabolismProperties()
    safety: SafetyProperties = SafetyProperties()
    renal: RenalProperties = RenalProperties()
    distribution: DistributionProperties = DistributionProperties()  # NEW
```

**Design notes:**

- Tissue-name validity is NOT checked at schema time. The full tissue
  universe depends on the loaded `PBPKTopology` (which can evolve per
  species YAML). A hardcoded tissue whitelist in the schema validator
  would drift. Tissue-name checking moves to
  `build_compound_pbpk_params`, which has access to `topology.tissues`
  and raises a clear error listing valid tissues.
- Upper bound 200 is generous (R&R caps at 50; Poulin-Theil rarely
  exceeds 100 for adipose of very lipophilic drugs). It catches
  percent-vs-ratio unit errors (e.g. someone typing `value=85.0` meaning
  85% which should be 0.85) only indirectly; the `source` field on
  `PredictedProperty` is the primary documentation mechanism and carries
  the responsibility.
- `compound_type` is a `Literal` on `PhysicochemicalProperties` (not on
  `CompoundConfig`) because it is a pKa-derived classification and
  logically co-located with the pKa values.
- `ConversionStep` is NOT modified. Override audit does not flow through
  the existing `HepaticClearance.conversion_log`; it is recorded
  separately on `CompoundPBPKParams` (see 3.3).

### 3.2 `build_compound_pbpk_params` override flow

Four changes to `src/charon/pbpk/ode_compiler.py::build_compound_pbpk_params`.

**3.2.1 compound_type resolution precedence (new)**

```python
resolved_type = (
    compound_type                                          # Pipeline kwarg (highest)
    or compound.properties.physicochemical.compound_type   # YAML field
    or infer_compound_type(pka_acid, pka_base)             # fallback
)
```

This lets the 10 validation YAML files specify the classification
directly (`physicochemical.compound_type: "base"` for midazolam,
propranolol, etc.) without requiring every caller to pass a keyword
argument. The Pipeline-level kwarg still wins, preserving Sprint 3a
behavior.

**3.2.2 Empirical Kp override (new)**

After `kp_by_tissue = compute_all_kp(...)`:

```python
override_log: list[KpOverrideRecord] = []
empirical = compound.properties.distribution.empirical_kp_by_tissue
if empirical:
    valid_tissues = set(kp_by_tissue.keys())
    for tissue, p in empirical.items():
        if tissue not in valid_tissues:
            raise ValueError(
                f"empirical_kp_by_tissue[{tissue!r}] refers to a tissue "
                f"not present in topology {topology.species!r}. "
                f"Valid tissues: {sorted(valid_tissues)}"
            )
        override_log.append(
            KpOverrideRecord(
                tissue=tissue,
                rr_value=kp_by_tissue[tissue],
                empirical_value=float(p.value),
                source=p.source,
                method=p.method,
                flag=p.flag,
            )
        )
        kp_by_tissue[tissue] = float(p.value)
```

**Design notes:**

- Unknown tissue → `ValueError` at build time (not silently ignored, not
  deferred to ODE runtime). Error message lists the valid tissues from
  the current topology, which is what the user needs.
- Partial overrides are supported: override only `adipose` for midazolam,
  leave the other 14 tissues to R&R.
- `kp_override_log` captures both the R&R-original and the empirical
  value per tissue so audit trails can show "what changed, and from
  what".

**3.2.3 Return `CompoundPBPKParams` with override audit**

```python
return CompoundPBPKParams(
    ...,  # existing fields unchanged
    kp_by_tissue=dict(kp_by_tissue),
    kp_overrides=tuple(override_log),  # NEW, default ()
)
```

**3.2.4 Hardcoded mppgl/hepatocellularity — flagged AND guarded**

`build_compound_pbpk_params` currently passes `mppgl=40.0` and
`hepatocellularity=120.0` to `bridge.clint_to_clh`. These are human
values. This is correct for Sprint 3a/3b (human-only) but will be a
latent bug in Sprint 4 (rat/dog).

This session does NOT fix the hardcoding (adding species-aware
`PBPKTopology.mppgl_mg_g` / `hepatocellularity_1e6_per_g` fields would
require preclinical validation to land safely, and that validation is
explicitly out of scope). But this session DOES add a **species guard**
at the top of `build_compound_pbpk_params` to ensure non-human topology
cannot silently consume human constants:

```python
if topology.species != "human":
    raise NotImplementedError(
        f"build_compound_pbpk_params currently hardcodes human "
        f"mppgl=40.0, hepatocellularity=120.0. "
        f"species={topology.species!r} requires species-aware values; "
        f"fix scheduled for Sprint 4 (translational layer)."
    )
```

Rationale: 2 lines of code block a class of silent-failure bugs for the
small price of a clear error message. If anyone (including a future
session) tries to run the current builder on rat/dog topology, they
hit a loud wall and know exactly what needs to change.

The `KNOWN_FUTURE_WORK` section in the final commit message still
documents the fix path: "add `mppgl_mg_g` and `hepatocellularity_1e6_per_g`
to `PBPKTopology`, read from species YAML, remove guard". That is
Sprint 4 work.

### 3.3 `CompoundPBPKParams` audit additions

```python
@dataclass(frozen=True)
class KpOverrideRecord:
    """Audit record for a single tissue-level empirical Kp override."""

    tissue: str
    rr_value: float          # what R&R computed
    empirical_value: float   # what was injected
    source: str              # PredictedProperty.source — "experimental" | "literature" | ...
    method: str | None       # PredictedProperty.method — e.g. "Björkman 2001 Table 3"
    flag: str | None         # PredictedProperty.flag — optional limitation flag


@dataclass(frozen=True)
class CompoundPBPKParams:
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
    kp_overrides: tuple[KpOverrideRecord, ...] = ()   # NEW
```

`Pipeline.run()` exposes the override log in
`PipelineResult.metadata["kp_overrides"]` as a list of dicts so downstream
reporting can render "Applied 1 empirical Kp override (adipose, source:
literature, Björkman 2001)".

### 3.4 No changes to ODE, solver, or PK extraction

`build_rhs`, `simulate_iv`, and `compute_pk_parameters` see Kp via
`params.kp_by_tissue` exactly as they do in Sprint 3a. Override happens
upstream of the ODE compiler; the ODE is oblivious to whether Kp came
from R&R or from literature. This keeps the ODE kernel fully isolated
and prevents regression on Sprint 3a's mass balance + BDF behavior.

---

## 4. Validation Panel Layout

### 4.1 Directory structure

The `validation/data/tier1_obach/` directory already exists as an empty
scaffold (created during initial project setup). This session populates
it:

```
validation/data/tier1_obach/
├── README.md                        NEW — one-paragraph: "Obach 1999
│                                         Tier-1 validation panel; loaded
│                                         by layer2_human_pk.py"
├── panel.yaml                       NEW — panel metadata + observed PK
└── compounds/                       NEW
    ├── theophylline.yaml            NEW — mirrors existing Python fixture
    ├── antipyrine.yaml              NEW
    ├── caffeine.yaml                NEW
    ├── warfarin.yaml                NEW — stereoisomer resolved at curation (§4.6 #1)
    ├── diclofenac.yaml              NEW
    ├── midazolam.yaml               NEW — includes empirical_kp override
    ├── propranolol.yaml             NEW
    ├── diazepam.yaml                NEW
    ├── metoprolol.yaml              NEW
    └── verapamil.yaml               NEW
```

### 4.2 `panel.yaml` schema

Observed PK values do NOT live inside `CompoundConfig` because
`CompoundConfig` is a *prediction/measurement input* model, not a
*literature reference* model. Polluting it with `observed_pk` would bleed
validation concerns into production schema. A separate panel descriptor
keeps the two concerns cleanly separated.

```yaml
# validation/data/tier1_obach/panel.yaml
name: "Obach 1999 Tier-1 Human IV PK Validation Panel"
source: "Obach RS, Drug Metab Dispos 27(11):1350-1359, 1999"
default_duration_h: 168.0

compounds:
  - key: theophylline
    compound_file: compounds/theophylline.yaml
    route: iv_bolus
    dose_mg: 100.0
    duration_h: 168.0
    observed:
      cl_L_h: 2.9
      vss_L: 35.0
      t_half_h: 8.0
    obach_table_row: 42           # row reference in Obach Table 6
    strict_targets: true          # regression gate retained from Sprint 3a
    notes: "Neutral low-logP primary validation compound."

  - key: midazolam
    compound_file: compounds/midazolam.yaml
    route: iv_bolus
    dose_mg: 5.0
    duration_h: 168.0
    observed:
      cl_L_h: 21.0
      vss_L: 66.0
      t_half_h: 3.0
    obach_table_row: 24
    strict_targets: false         # non-gated; documents R&R limitation
    notes: "Weak base, CYP3A4. empirical_kp override in compound YAML."

  # ... 8 more entries ...
```

The `key` field is the short programmatic name used throughout the
benchmark output and test code. `compound_file` is a path relative to
`panel.yaml`'s directory. `strict_targets: true` promotes a compound to a
2-fold hard gate identical to the Sprint 3a theophylline test (the new
benchmark respects it and fails with exit 1 if any strict-targets
compound misses 2-fold on any metric).

### 4.3 Compound YAML example (theophylline)

Each file is a `CompoundConfig` literal, valid against the Pydantic
schema. Every `PredictedProperty` declares `source` and `method`
(literature citation), so the audit trail survives load → pipeline →
report.

```yaml
# validation/data/tier1_obach/compounds/theophylline.yaml
name: theophylline
smiles: "Cn1c(=O)c2[nH]cnc2n(C)c1=O"
molecular_weight: 180.17
source: experimental
properties:
  physicochemical:
    logp:
      value: -0.02
      source: experimental
      method: "Obach 1999 Table 2"
    compound_type: neutral
  binding:
    fu_p:
      value: 0.60
      source: experimental
      unit: fraction
      method: "Obach 1999 Table 2"
    fu_inc:
      value: 1.0
      source: experimental
      unit: fraction
      method: "Austin 2002 (fu_inc=1 for low logP)"
    bp_ratio:
      value: 0.85
      source: experimental
      unit: ratio
      method: "Obach 1999 Table 2"
  metabolism:
    clint_uL_min_mg:
      value: 1.8
      source: experimental
      unit: uL/min/mg
      method: "Obach 1999 Table 2"
  renal:
    clrenal_L_h:
      value: 0.1
      source: experimental
      unit: L/h
      method: "Obach 1999 Table 2"
```

### 4.4 midazolam YAML with override

The midazolam YAML is authored in two phases.

**Phase 1 (initial commit, no override yet)**: identical structure to
the theophylline example but with `pka_base`, `compound_type: base`, and
the midazolam-specific numeric values. This file parses, runs through
the pipeline, and reproduces the Sprint 3a limitation case (Vss 6×
observed).

```yaml
# validation/data/tier1_obach/compounds/midazolam.yaml  (Phase 1)
name: midazolam
smiles: "Cc1ncc2n1-c1ccc(Cl)cc1C(c1ccccc1F)=NC2"
molecular_weight: 325.77
source: experimental
properties:
  physicochemical:
    logp: {value: 3.89, source: experimental, method: "Obach 1999 Table 2"}
    pka_base: {value: 6.2, source: experimental, method: "Literature; tertiary amine"}
    compound_type: base                 # overrides threshold-based inference
  binding:
    fu_p: {value: 0.03, source: experimental, unit: fraction, method: "Obach 1999 Table 2"}
    fu_inc: {value: 0.96, source: experimental, unit: fraction, method: "Austin 2002"}
    bp_ratio: {value: 0.66, source: experimental, unit: ratio, method: "Obach 1999 Table 2"}
  metabolism:
    clint_uL_min_mg: {value: 93.0, source: experimental, unit: uL/min/mg, method: "Obach 1999 Table 2"}
  renal:
    clrenal_L_h: {value: 0.0, source: experimental, unit: L/h, method: "Negligible renal clearance"}
```

**Phase 2 (after citation verification)**: the implementation plan adds
a `distribution.empirical_kp_by_tissue.adipose` entry with a numerical
value read directly from a verified primary source (see Section 7.2).
The numerical value is **not pinned in this design spec** because
pinning a number without an open paper in front of the author would be
fabrication. The citation-verification task is a hard prerequisite
before Phase 2 is committed.

**Interaction with DoD §2 ≥2 verified-override rule**: midazolam is the
primary override compound but not the only one. The session also adds a
second verified override on one of propranolol/diazepam (see DoD §2).
The rule is: **Fallback 3 (synthetic_test) never counts toward the ≥2
verified quota**. With only two override targets in this session, that
means either compound falling to Fallback 3 brings the verified count
below 2, which triggers escalation (pause, report, ask for direction)
per Section 12 §6. The implementation plan must NOT silently commit a
synthetic override on midazolam even if the secondary compound has a
verified citation — Fallback 3 is only usable in future sessions that
carry more override targets.

### 4.5 Compound selection rationale

10 compounds across 3 R&R classification branches, all from Obach 1999.
The "R&R class" column is the classification used in the compound YAML
for `PhysicochemicalProperties.compound_type` — this is a *mechanistic*
classification, not the compound's chemical identity.

| Compound | R&R class | Chem. identity | Primary CYP | Extraction | Role |
|---|---|---|---|---|---|
| theophylline | neutral | neutral | CYP1A2 | low | primary gate, Sprint 3a carryover |
| antipyrine | neutral | neutral | CYP multiple | low | phenotyping probe, clean PK |
| caffeine | neutral | neutral | CYP1A2 | low | low-clearance case |
| diazepam | neutral¹ | weak base | CYP2C19/3A4 | low | long t½; pKa_base 3.4 → neutral at pH 7.4 |
| warfarin | acid² | acid | CYP2C9 | low | highly protein-bound acid |
| diclofenac | acid | acid | CYP2C9 | mid | moderate extraction acid |
| midazolam | base³ | weak base | CYP3A4 | high | R&R limitation case, override demo |
| propranolol | base | base | CYP2D6 | high | high extraction base |
| metoprolol | base | base | CYP2D6 | mid | β-blocker base |
| verapamil | base | base | CYP3A4 | very high | engine stress test (ER > 0.9) |

**Footnotes:**

1. **diazepam** is chemically a benzodiazepine (weak base, pKa_base ≈
   3.4). At physiological pH 7.4, the ionized fraction is ~0.01%, so
   diazepam behaves as effectively neutral. The R&R Kp calculation uses
   the neutral branch, and `compound_type: neutral` is recorded in the
   YAML. See Section 4.6 #2.
2. **warfarin** (racemate vs S-enantiomer) is resolved at curation time
   — see Section 4.6 #1. The YAML `name` field records which form the
   values correspond to.
3. **midazolam** pKa_base 6.2 is below R&R's strict "moderate base"
   threshold of 7. Sprint 3a explicitly classified it as "base" to
   exercise the base-branch equations, and the same convention is
   preserved here. This classification is a deliberate choice carried
   forward from Sprint 3a, not a silent default.

**Panel composition by R&R class:** 4 neutral, 2 acid, 4 base. The base
group is where R&R's known failure mode lives and where the override
path demonstrates value. The neutral group anchors the panel to "R&R
works here" baseline. The acid group exercises the R&R acid branch,
which Sprint 3a never tested.

**No zwitterion**: cetirizine is the obvious candidate but its IV data
in healthy adults is less standard than the others. Documented as a gap;
zwitterion coverage deferred until a reliable IV reference is identified.

### 4.6 Known curation gotchas (for the implementation plan)

Hand-authoring 10 YAMLs with per-field citations is a 2-3 hour curation
effort, not a trivial pattern-fill. The following per-compound
pharmacological judgments must be decided at curation time and recorded
as YAML comments + commit message notes:

1. **Warfarin stereoisomer.** Obach 1999 Table 2 reports a single
   "warfarin" row whose values may be racemic or S-specific. The
   curation task MUST verify which form the reported values correspond
   to, and the YAML must declare it explicitly
   (`name: "warfarin (S-enantiomer)"` or equivalent). Mis-stating the
   form is a scientific error that would poison the AAFE metric.

2. **Diazepam classification.** Diazepam has `pKa_base ≈ 3.4`, well
   below R&R's "moderate-to-strong base" threshold of 7. At
   physiological pH, diazepam is predominantly neutral, so
   `compound_type: "neutral"` (NOT "base") is pharmacologically correct
   for the R&R Kp calculation. Document the decision in the YAML and
   confirm against the `infer_compound_type` default (which would also
   say neutral since `pKa_base > 8.0` is required). This is opposite to
   midazolam (pKa_base 6.2, too low for default base threshold but
   correctly set to "base" at the YAML level because of its experimental
   behavior).

3. **Diclofenac fu_p edge case.** Diclofenac `fu_p ≈ 0.005` triggers
   the ParameterBridge "very low fu_p" warning (threshold 0.01 per
   CLAUDE.md §6j Pitfall #10). The warning is expected; the YAML must
   include a comment acknowledging it, and the benchmark must not treat
   the warning as a failure.

4. **Verapamil ER ≈ 0.95.** Verapamil runs the well-stirred liver
   model at the extreme saturation limit where `CLh ≈ Q_H`. Small
   errors in `fu_b` or `CLint_liver` cause large fold errors in CL
   prediction. Expect a wider fold error for verapamil than for the
   other compounds; this is the engine stress test, not a bug.

5. **CYP2D6 poor-metabolizer variability.** Metoprolol is a CYP2D6
   substrate; Obach reports population-mean CL. Individual variability
   is large but out of scope for point-estimate PBPK. The YAML cites
   the population-mean value explicitly.

6. **Obach table row numbering.** Obach 1999 Table 2 (in vitro) and
   Table 6 (in vivo) use different row orderings. Each compound YAML
   and panel entry must cite BOTH tables' row references
   (`obach_table_2_row` for in vitro values, `obach_table_6_row` /
   `obach_table_row` for observed CL/Vss/t½). The design's earlier
   `obach_table_row` field refers to Table 6 (observed PK).

7. **Contingency — compound unavailable.** If any one of the 10
   compounds cannot be curated from Obach 1999 with acceptable data
   quality (e.g., observed CL/Vss measured in patients, not healthy
   volunteers, or IV data missing), the implementation plan is
   permitted to **swap** the compound for a close substitute within
   the same class (e.g., replace warfarin with tolbutamide,
   metoprolol with atenolol). The swap is logged in the session
   commit message; the panel size of 10 is preserved.

These gotchas live here in the spec (not the plan) because they are
*decisions* requiring pharmacological judgment, not *tasks*. The plan
references this section when authoring the curation steps.

---

## 5. Benchmark Harness Refactor

### 5.1 Data types (new, inside `layer2_human_pk.py`)

```python
@dataclass
class PanelEntry:
    key: str
    compound: CompoundConfig         # loaded from compound_file (classification
                                     # comes from compound.properties.physicochemical.compound_type)
    route: str                       # "iv_bolus" | "iv_infusion"
    dose_mg: float
    duration_h: float
    observed: dict[str, float]       # cl_L_h, vss_L, t_half_h
    strict_targets: bool
    obach_table_row: int | None
    notes: str

@dataclass
class PanelRow:
    key: str
    predicted: dict[str, float]      # cl, vss, t_half
    observed: dict[str, float]
    fold: dict[str, float]           # per-metric fold error
    pass_2_fold: dict[str, bool]
    override_tissues: list[str]      # from PipelineResult.metadata['kp_overrides']
    mode: str                        # "no_override" | "with_override"

@dataclass
class PanelSummary:
    n: int
    mode: str
    aafe: dict[str, float]           # {cl, vss, t_half}
    within_2_fold: dict[str, float]  # fraction [0, 1] per metric
    within_3_fold: dict[str, float]
    strict_failures: int             # count of strict-target rows failing any metric
```

### 5.2 Execution flow

```
main(panel_path):
    panel = load_panel(panel_path)         # parses panel.yaml + all compound_files
    rows_no_override = []
    rows_with_override = []

    for entry in panel:
        # Pass 1: strip empirical_kp from the compound, run
        stripped = entry.compound.model_copy(
            update={'properties': _without_kp_overrides(entry.compound.properties)}
        )
        row_no = _run_one(stripped, entry, mode="no_override")
        rows_no_override.append(row_no)

        # Pass 2: run as-is (possibly with override)
        row_yes = _run_one(entry.compound, entry, mode="with_override")
        rows_with_override.append(row_yes)

    print_table(rows_no_override, title="R&R only")
    print_table(rows_with_override, title="R&R + empirical overrides")

    summary_no = aggregate(rows_no_override, mode="no_override")
    summary_yes = aggregate(rows_with_override, mode="with_override")

    print_summary(summary_no)
    print_summary(summary_yes)

    strict_failures = summary_yes.strict_failures   # gating uses the override pass
    return 0 if strict_failures == 0 else 1
```

**Two-pass rationale**: running each compound twice (with and without
empirical override) produces the side-by-side "what did the override
actually change" report. For the 8 compounds without overrides, both
passes are identical; the compounds with overrides (midazolam + at
least one more per DoD §2) are the only rows that actually move
between the two passes.

**Presentation discretion (implementation-time)**: the spec example
below prints both full tables for visual symmetry, but the
implementation is free to collapse the "with overrides" table to show
only the rows that differ from the "R&R only" baseline (plus a
"unchanged: N compounds" summary line). Either rendering is acceptable
so long as the panel-level AAFE summaries are printed for BOTH modes
and per-compound fold errors are retrievable from the returned summary
object for downstream use.

The strip function uses **nested** `model_copy` so that any future
additions to `DistributionProperties` (e.g. `vss_pred`, `tissue_fu`) are
preserved instead of being silently reset to defaults:

```python
def _without_kp_overrides(props: CompoundProperties) -> CompoundProperties:
    """Return a copy of props with empirical_kp_by_tissue cleared.

    Uses nested model_copy so other fields of DistributionProperties
    that may be added in future (e.g. Vss_pred, tissue-level fu) are
    preserved rather than silently reset to their defaults.
    """
    new_distribution = props.distribution.model_copy(
        update={'empirical_kp_by_tissue': None}
    )
    return props.model_copy(update={'distribution': new_distribution})
```

A unit test in `test_panel_loader.py` must verify that if
`DistributionProperties` ever gains a new field, the strip function
preserves it (a forward-compat guardrail test using a local subclass
or monkey-patched field suffices).

### 5.3 Stdout format

```
================================================================================
Charon Layer 2 Human PBPK Benchmark — Obach 1999 Tier-1 Panel
================================================================================

R&R only
--------------------------------------------------------------------------------
Compound         CL_pred CL_obs  fold  | Vss_pred Vss_obs fold  | t½_pred t½_obs fold
theophylline      3.26    2.90   1.12  |   37.1    35.0   1.06  |   9.5    8.0   1.19  PASS
antipyrine        ...
...
verapamil         ...
--------------------------------------------------------------------------------
Panel AAFE (R&R only)     CL: 1.82    Vss: 4.12    t½: 2.05
Within 2-fold             CL: 7/10    Vss: 4/10    t½: 6/10
Within 3-fold             CL: 9/10    Vss: 6/10    t½: 8/10

R&R + empirical overrides
--------------------------------------------------------------------------------
... same rows, midazolam Vss changed ...
--------------------------------------------------------------------------------
Panel AAFE (with override) CL: 1.82   Vss: 2.87    t½: 1.84
Within 2-fold              CL: 7/10   Vss: 6/10    t½: 7/10
Overrides applied:
  midazolam : adipose (literature, <citation>)
================================================================================
All strict-targets compounds PASS — exit 0
```

---

## 6. Testing Strategy

New test files (counts are estimates):

| File | Tests | Focus |
|---|---|---|
| `tests/unit/test_schema_distribution.py` | ~8 | `DistributionProperties` validators, `PhysicochemicalProperties.compound_type`, backward-compat round-trip of existing fixtures |
| `tests/unit/test_kp_override.py` | ~10 | `build_compound_pbpk_params` override path: no-op when empty, single-tissue replace, multi-tissue replace, unknown tissue raises with valid list, `kp_overrides` audit record shape, mass conservation post-override |
| `tests/unit/test_compound_type_resolution.py` | ~6 | Precedence: Pipeline kwarg > YAML field > inferred. All three paths exercised. |
| `tests/unit/test_panel_loader.py` | ~6 | `panel.yaml` → list[PanelEntry], missing compound_file raises, dose_mg/duration_h resolution from default + override, observed fields required |
| `tests/unit/test_panel_metrics.py` | ~6 | `aggregate()` on synthetic rows: known fold errors → expected AAFE / within_n_fold, handles mode labeling |
| `tests/unit/test_pipeline_midazolam_override.py` | ~3 | Regression: existing `TestPipelineMidazolamLimitation` unchanged; new tests use synthetic `empirical_kp_by_tissue={'adipose': 10.0}` to verify (a) `PipelineResult.metadata['kp_overrides']` contains exactly one entry with correct tissue/value, (b) `Vss` with override is strictly less than `Vss` without override (direction-of-effect), (c) no-override baseline still matches Sprint 3a `TestPipelineMidazolamLimitation` expected values. Does NOT assert a specific numerical fold-error target. |
| `tests/integration/test_obach_panel_smoke.py` | ~3 | End-to-end: all 10 compounds run without error, all PK values finite positive, panel AAFE < sanity floor (5.0 for Vss), strict-target compounds pass their 2-fold gate. Uses `record_property` so CI artifacts preserve the measured AAFE. |

Expected new tests: ~42.

**Pre-existing test handling**: Sprint 3a's `test_pipeline.py::
TestPipelineTheophyllineValidation` and `TestPipelineMidazolamLimitation`
continue to use their in-Python `CompoundConfig` factories unchanged.
The new compound_type resolution precedence (kwarg > YAML > inferred)
preserves their behavior exactly: they pass `compound_type_override`
as a Pipeline kwarg, which wins over any YAML field. The Sprint 3a
`validation/benchmarks/layer2_human_pk.py` Python factory functions
(`theophylline()`, `midazolam()`) ARE rewritten out — they are replaced
by YAML loading. This is a rewrite of a script, not modification of a
test, and has no downstream consumers beyond the benchmark itself.

**Test-writing convention (carryover from Sprint 3a)**: every
pharmacological math test includes the hand-calculation in its docstring.
For Sprint 3b the hand-calc surface is mostly in metric aggregation
(fold, AAFE) rather than new ODE math, so the docstring calcs are shorter.

**No new pharmacology math is introduced.** All new tests verify
*plumbing* (schema, override routing, YAML loading, metric aggregation)
against the existing Sprint 3a engine. This keeps the cognitive load
small and the failure modes easy to diagnose.

---

## 7. Known-Limitation Policy (Honest Output Mode)

### 7.1 Why dual-pass reporting

Reporting only the "with-override" panel would hide the R&R limitation
from future readers of the benchmark log. Reporting only the "no-override"
panel would hide the functional result the engine can achieve when given
curated inputs. Printing both leaves a permanent record of:

1. How far R&R alone gets us (scientific baseline).
2. How much an empirical Kp override helps (mechanism value).
3. Which compounds and tissues required intervention (targeted
   diagnostics for future Kp model work).

This is the honest substitute for "fix the Kp model", which is out of
scope this session.

### 7.2 Citation protocol for override values

For the midazolam adipose override, the primary target citation is
Björkman S. (2001/2002) "Prediction of drug disposition in infants and
children by means of physiologically based pharmacokinetic (PBPK)
modeling", or an equivalent peer-reviewed midazolam PBPK paper that
reports tissue Kp or V_adipose contribution explicitly.

The implementation plan includes a **citation-verification task** with
explicit completion criteria:

1. Read the actual paper (or equivalent peer-reviewed secondary source
   that reproduces the table: Rodgers & Rowland 2005/2006 review,
   Poulin & Theil 2002, or similar).
2. Record the exact numerical Kp_adipose value as printed in the
   source.
3. Record the page number, table number, and row identifier in the
   YAML `method` field as a free-text citation string (e.g.
   `"Björkman 2001 J Pharm Sci 90(7):928-941, Table 2 row 'adipose'"`).
4. If step 1 fails (paper inaccessible), move to Fallback 1.

The task is **hard-gated**: the numerical override may not be committed
until steps 1-3 are complete, or until Fallback 1/2/3 is explicitly
invoked with a recorded rationale in the commit message.

**Fallback 1**: if the Björkman paper cannot be obtained, use Rodgers et
al. 2005 Table 4 / 5 rat experimental Kp values for the closest
structural analog from their panel (diazepam is the most likely
candidate; Rodgers' panel explicitly excludes midazolam).

**Fallback 2**: if no primary citation can be verified at all, the
override demo switches to **diazepam** instead of midazolam (diazepam's
PBPK literature is older and denser). Diazepam is already in the panel,
so the switch is a re-pointing of which compound's YAML carries the
override, not a panel restructuring.

**Fallback 3**: if all citation paths fail for BOTH the primary
(midazolam) and secondary (propranolol / diazepam) override compounds,
the session is **incomplete** per DoD §2 (≥2 verified-citation
override rule). The implementation plan must escalate rather than
ship: pause execution, report the blockage, and ask for direction
(e.g. swap the override target to a different compound, extend the
citation search, accept a reduced panel of 1 verified override with
explicit DoD waiver, etc.). Fallback 3 is therefore **forbidden as a
silent fallback** in this session; it may only be invoked per-compound
when ≥2 other verified overrides exist (which is only possible in
future sessions that expand the override target list beyond two).

In practice for Sprint 3b Session 1: Fallback 3 is dead code. It is
documented so that future sessions (which may carry more override
compounds) have a mechanism for a single degraded entry, and so that
the test plumbing can exercise the `source: "synthetic_test"` code
path without requiring a fallback trigger in production.

### 7.3 "Not tuning" principle

The override values are literature-sourced or explicitly-marked
synthetic. They are **never** chosen to make Vss match observed. If a
reader wants to do that experiment they can, but the schema's `source`
field forces them to declare it. No `source: "tuned"` appears in the
panel; if it ever does, the benchmark script will print a warning row
("WARNING: tuned value used — scientific validity compromised") to make
the compromise visible.

---

## 8. Risks and Mitigations

| Risk | Impact | Mitigation |
|---|---|---|
| Obach 1999 Table value curation error (typo in CL/Vss/t½) | AAFE metrics wrong; baseline poisoned | Every numeric value cites `obach_table_row`; implementation plan has a "cross-check against primary source" task per compound; tests assert specific observed values so typos surface in the smoke test diff |
| Kp override citation (Björkman) cannot be verified | Midazolam adipose override blocked | 3-tier fallback cascade (Section 7.2): Björkman → Rodgers 2005 → swap to diazepam. DoD §2 ≥2 verified rule prevents silent all-synthetic ship — if citation cascade fails for both primary and secondary override targets, session escalates rather than degrades |
| `distribution` schema addition breaks existing YAML round-trip | Sprint 3a fixtures fail to load | `distribution: DistributionProperties = DistributionProperties()` default ensures field-absent YAML parses unchanged; explicit test round-trips the existing theophylline and midazolam Python fixtures through schema serialize → deserialize before touching `build_compound_pbpk_params` |
| Tissue-name validation drift (benchmark uses `"fat"` but topology has `"adipose"`) | Runtime error, confusing for users | `build_compound_pbpk_params` raises with the valid-tissue list included in the error message; a unit test covers the exact error text |
| Curated `compound_type: "base"` for borderline bases (pKa_base 6.2) changes Kp classification silently for anyone relying on the old "neutral" classification | Unexpected change in predicted Vss for existing users | Only the panel YAMLs set `compound_type`; no existing fixture is modified; `infer_compound_type` thresholds stay untouched; backward-compat tests verify no-change on existing fixtures |
| Benchmark wall time grows 10× with panel expansion | Slow CI | Expected to remain well under CI per-test budget (Sprint 3a theophylline integration test runs in a couple of seconds; 10 compounds in sequence should stay in the same order of magnitude). Wall time not enforced this session; if it becomes a problem, the benchmark can be marked `@pytest.mark.slow` and excluded from default CI. |
| Panel-level AAFE_Vss > 5.0 (sanity floor trip) | Session not complete; root cause unknown | DoD §7 gates this explicitly. Likely root causes: Obach curation error, compound_type mis-classification, or an actually-broken Kp model path in Sprint 3a that theophylline's neutral chemistry didn't exercise. Triage protocol in implementation plan. |
| `mppgl` / `hepatocellularity` remain hardcoded to human values | Session is clean; Sprint 4 rat/dog will need a fix | Document as known future work in a `KNOWN_FUTURE_WORK` section of the final commit message; do NOT attempt to fix here (adds risk without validation reward). |
| Two-pass benchmark runs each compound twice → doubled solver load | Slight CI wall time increase | Acceptable; alternative (caching first-pass results, or skipping pass 2 for compounds without overrides) adds complexity for negligible gain. If wall time actually becomes problematic, the "skip pass 2 when entry has no overrides" optimization is a 5-line change. |

---

## 9. Non-Goals (Explicit)

The following are deferred to later sessions. Do not touch them.

- ❌ **ACAT / oral / dissolution / BCS / food effect** — Sprint 3b Session 2
- ❌ **rat / dog / monkey PBPK validation** — Sprint 3b Session 3 or Sprint 4
- ❌ **Berezhkovskiy default transition** — design choice in this session
  explicitly rejects it as insufficient for the weak-base Vss problem
- ❌ **KP_MAX change** — stays at 50 per Sprint 3a convention
- ❌ **Poulin-Theil / hybrid Kp models** — no new Kp methods
- ❌ **Species-aware `mppgl` / `hepatocellularity`** — Sprint 4 work
- ❌ **Layer 1 ML retraining** — not touched
- ❌ **CLI changes** — Sprint 6
- ❌ **Report generator** — Sprint 6
- ❌ **Solver wall-time enforcement** — deferred
- ❌ **Phase B (population, DDI, API) directories** — not touched
- ❌ **Phase C (dashboard) directories** — not touched
- ❌ **`ConversionStep` schema modification** — audit additions flow
  through `CompoundPBPKParams.kp_overrides`, not the existing log
- ❌ **Experimental zwitterion coverage** — gap documented, no action

---

## 10. Known Future Work (Captured, Not Addressed)

Items discovered during design that must be fixed before later sprints
but are explicitly NOT in scope this session:

1. **Species-aware MPPGL / hepatocellularity.**
   `build_compound_pbpk_params` passes hardcoded `40.0` / `120.0`. Add
   `mppgl_mg_g` and `hepatocellularity_1e6_per_g` to `PBPKTopology` and
   populate from species YAML. Blocks Sprint 4.

2. **Zwitterion coverage in the Obach panel.**
   Cetirizine or ceftriaxone have IV data but curation is less clean.
   Add in Sprint 3b Session 3 or later.

3. **compound_type threshold revisit.**
   Current infer_compound_type thresholds (`pka_acid < 7.0`,
   `pka_base > 8.0`) leave midazolam (pKa_base=6.2) classified as
   neutral. A physicochemically sound threshold is `pka_base > 7.4`
   (ionized at physiological pH). Deferred because changing the default
   behavior of an ML-populated Layer 1 compound is higher-risk than
   adding YAML-level override.

4. **Automated Obach curation from a machine-readable source.**
   Hand-authoring 10 YAMLs this session is acceptable; scaling to the
   full Obach dataset (N≥40) requires a loader that consumes a CSV/JSON
   of the table. Sprint 6 validation effort.

5. **"Tuned" source audit.**
   If the `source` field on a `PredictedProperty` is ever set to
   `"tuned"`, the benchmark should print a loud warning. This session
   does not add that path; adding it in future requires a small
   extension to the `SourceType` literal and a warning pipeline.

---

## 11. File Layout Summary

```
src/charon/
├── core/
│   └── schema.py                       UPDATE  DistributionProperties (new),
│                                               PhysicochemicalProperties.compound_type,
│                                               CompoundProperties.distribution
├── pbpk/
│   └── ode_compiler.py                 UPDATE  KpOverrideRecord (new),
│                                               CompoundPBPKParams.kp_overrides,
│                                               build_compound_pbpk_params:
│                                                 - species!='human' guard
│                                                 - compound_type precedence
│                                                 - empirical override loop
│                                                 - unknown-tissue error
└── pipeline.py                         UPDATE  metadata['kp_overrides'] exposed

tests/
├── unit/
│   ├── test_schema_distribution.py           NEW
│   ├── test_kp_override.py                   NEW
│   ├── test_compound_type_resolution.py      NEW
│   ├── test_panel_loader.py                  NEW
│   ├── test_panel_metrics.py                 NEW
│   └── test_pipeline_midazolam_override.py   NEW
└── integration/
    └── test_obach_panel_smoke.py             NEW

validation/
├── data/
│   └── tier1_obach/
│       ├── README.md                         NEW
│       ├── panel.yaml                        NEW
│       └── compounds/
│           ├── theophylline.yaml             NEW
│           ├── antipyrine.yaml               NEW
│           ├── caffeine.yaml                 NEW
│           ├── warfarin.yaml                 NEW
│           ├── diclofenac.yaml               NEW
│           ├── midazolam.yaml                NEW  (with override)
│           ├── propranolol.yaml              NEW
│           ├── diazepam.yaml                 NEW
│           ├── metoprolol.yaml               NEW
│           └── verapamil.yaml                NEW
└── benchmarks/
    └── layer2_human_pk.py                    REWRITE — load panel.yaml,
                                                        two-pass metrics

docs/superpowers/
├── specs/
│   └── 2026-04-10-sprint3b-kp-override-validation-design.md   (this doc)
└── plans/
    └── 2026-04-10-sprint3b-kp-override-validation.md          (next: writing-plans)
```

---

## 12. Autonomous Execution Protocol (carryover)

The user granted autonomous execution authority for this session
("자율모드 진행"). Terminal gates:

1. All 477 pre-existing tests continue to pass. No regressions.
2. All new tests pass on first `pytest` run after each module is written.
3. Pre-commit hooks (if any) pass; do not `--no-verify`.
4. Commit messages describe deliverables concretely.
5. On any failure: diagnose root cause before retrying.
6. Do not commit override numerical values without verified citations.
   If all citation paths fail for both primary and secondary override
   compounds, **escalate rather than degrade**: pause, report the
   blockage, and ask for direction. Do NOT silently ship all-synthetic
   overrides (see DoD §2 ≥2 verified rule).
7. Session ends only when (a) all gates above pass, (b) panel AAFE_Vss
   with-override < 5.0 sanity floor, (c) session summary explicitly
   states measured panel-level AAFE values for all 10 compounds on all
   3 metrics, (d) ≥2 compounds have verified-citation empirical Kp
   overrides recorded in their YAML files.

---

## 13. Approval Checklist

- [ ] Schema additions are additive and backward-compatible
- [ ] Tissue-name validation at build time, not schema time
- [ ] Kp override audit lives on `CompoundPBPKParams`, not on
      `ConversionStep`
- [ ] `compound_type` precedence is kwarg > YAML > inferred
- [ ] Panel uses existing `validation/data/tier1_obach/` scaffold
- [ ] Observed PK lives in `panel.yaml`, not in `CompoundConfig`
- [ ] Two-pass benchmark reports both no-override and with-override
      metrics
- [ ] `_without_kp_overrides` uses nested `model_copy` to preserve
      future `DistributionProperties` fields
- [ ] `build_compound_pbpk_params` adds `species != "human"` guard
      to surface Sprint 4 TODO loudly
- [ ] Override demo requires **≥2** compounds with verified citations
      before session is complete
- [ ] Override numerical values require verified citations (3-tier
      fallback plan documented)
- [ ] Per-compound curation gotchas (diazepam classification,
      warfarin stereoisomer, diclofenac fu_p, verapamil ER) documented
- [ ] Compound-swap contingency permitted within same class
- [ ] No Kp model change, no KP_MAX change, no new Kp method
- [ ] No ACAT, no oral, no rat/dog
- [ ] `mppgl`/`hepatocellularity` hardcoded issue is **guarded** (not
      fixed) and documented as known future work
- [ ] Theophylline strict 2-fold gate retained
- [ ] Session post-condition: AAFE_Vss (with-override) < 5.0 sanity floor
- [ ] Residual scientific gap (measured AAFE may exceed ARCHITECTURE
      target) is explicitly accepted as session outcome
