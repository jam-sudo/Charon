# Sprint 3b Session 1 — Honest IV Kernel (Kp Override + Obach Panel) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a schema-level empirical Kp override path to the Sprint 3a IV PBPK kernel, curate a 10-compound Obach 1999 validation panel loaded from YAML, and extend the benchmark harness with dual-pass AAFE reporting — so Layer 2 Human PK becomes measurable rather than merely claimable.

**Architecture:** Additive schema changes (`DistributionProperties`, `PhysicochemicalProperties.compound_type`) + additive `build_compound_pbpk_params` override loop + separate `KpOverrideRecord` audit channel (no modification of `ConversionStep`). Benchmark script rewritten to load `panel.yaml` + per-compound `CompoundConfig` YAMLs, run each compound twice (no-override / with-override), and aggregate per-metric AAFE. R&R Kp model itself is NOT changed; `KP_MAX=50` stays; no new Kp methods.

**Tech Stack:** Python 3.11+, Pydantic v2, scipy (BDF solver, existing), RDKit (existing), pytest, YAML (for compound configs and panel descriptor).

**Spec:** `docs/superpowers/specs/2026-04-10-sprint3b-kp-override-validation-design.md` (commits 1a84679, e35235b, 0051429)

---

## Pre-flight Checklist

Run before Task 1 to confirm the Sprint 3a baseline is green:

- [ ] **Pre-flight step: Confirm Sprint 3a baseline passes**

```bash
cd /home/jam/Charon
pytest tests/ -q 2>&1 | tail -20
```
Expected: `477 passed` (or similar — record the number; this is the regression baseline the session must preserve).

- [ ] **Pre-flight step: Confirm current benchmark still runs**

```bash
cd /home/jam/Charon
python3 validation/benchmarks/layer2_human_pk.py 2>&1 | tail -20
```
Expected: exit 0, theophylline passes 2-fold, midazolam prints fold-errors (non-gated).

---

## Task 1: Add `DistributionProperties` schema class

**Files:**
- Modify: `src/charon/core/schema.py`
- Create: `tests/unit/test_schema_distribution.py`

- [ ] **Step 1.1: Write failing tests for `DistributionProperties`**

Create `tests/unit/test_schema_distribution.py`:

```python
"""Tests for DistributionProperties and PhysicochemicalProperties.compound_type."""

import pytest
from pydantic import ValidationError

from charon.core.schema import (
    CompoundProperties,
    DistributionProperties,
    PhysicochemicalProperties,
    PredictedProperty,
)


class TestDistributionPropertiesDefault:
    def test_default_is_none(self):
        dp = DistributionProperties()
        assert dp.empirical_kp_by_tissue is None

    def test_explicit_none(self):
        dp = DistributionProperties(empirical_kp_by_tissue=None)
        assert dp.empirical_kp_by_tissue is None


class TestDistributionPropertiesValidation:
    def _mk_kp(self, value: float) -> PredictedProperty:
        return PredictedProperty(
            value=value, source="experimental", unit="ratio"
        )

    def test_single_tissue_ok(self):
        dp = DistributionProperties(
            empirical_kp_by_tissue={"adipose": self._mk_kp(10.0)}
        )
        assert dp.empirical_kp_by_tissue is not None
        assert dp.empirical_kp_by_tissue["adipose"].value == 10.0

    def test_multi_tissue_ok(self):
        dp = DistributionProperties(
            empirical_kp_by_tissue={
                "adipose": self._mk_kp(10.0),
                "muscle": self._mk_kp(3.0),
            }
        )
        assert len(dp.empirical_kp_by_tissue) == 2

    def test_empty_dict_rejected(self):
        with pytest.raises(ValidationError, match="non-empty dict"):
            DistributionProperties(empirical_kp_by_tissue={})

    def test_zero_value_rejected(self):
        with pytest.raises(ValidationError, match="physiological range"):
            DistributionProperties(
                empirical_kp_by_tissue={"adipose": self._mk_kp(0.0)}
            )

    def test_negative_value_rejected_at_predicted_property(self):
        # PredictedProperty itself does not reject negatives, but the
        # DistributionProperties validator's (0, 200] bound does.
        with pytest.raises(ValidationError, match="physiological range"):
            DistributionProperties(
                empirical_kp_by_tissue={"adipose": self._mk_kp(-1.0)}
            )

    def test_upper_bound_200_rejected(self):
        with pytest.raises(ValidationError, match="physiological range"):
            DistributionProperties(
                empirical_kp_by_tissue={"adipose": self._mk_kp(201.0)}
            )

    def test_upper_bound_exactly_200_ok(self):
        dp = DistributionProperties(
            empirical_kp_by_tissue={"adipose": self._mk_kp(200.0)}
        )
        assert dp.empirical_kp_by_tissue["adipose"].value == 200.0
```

- [ ] **Step 1.2: Run tests to verify they fail**

```bash
pytest tests/unit/test_schema_distribution.py -v
```
Expected: FAIL with `ImportError: cannot import name 'DistributionProperties' from 'charon.core.schema'`

- [ ] **Step 1.3: Add `DistributionProperties` class to schema.py**

In `src/charon/core/schema.py`, after `class RenalProperties` and before `# Compound-level aggregates` separator, add:

```python
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
```

- [ ] **Step 1.4: Run tests to verify they pass**

```bash
pytest tests/unit/test_schema_distribution.py -v
```
Expected: 8 tests PASS.

- [ ] **Step 1.5: Run full test suite to verify no regression**

```bash
pytest tests/ -q 2>&1 | tail -10
```
Expected: 477 existing + 8 new = 485 passed (or baseline +8).

- [ ] **Step 1.6: Commit**

```bash
git add src/charon/core/schema.py tests/unit/test_schema_distribution.py
git commit -m "$(cat <<'EOF'
Sprint 3b: Add DistributionProperties schema class

Additive schema change. Introduces DistributionProperties with a single
empirical_kp_by_tissue: dict[str, PredictedProperty] | None field for
per-tissue empirical Kp overrides. Validator enforces positive values
within (0, 200] physiological range; empty dict rejected (must be None
or populated). Not yet wired into CompoundProperties — that lands in
the next task.

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 2: Add `PhysicochemicalProperties.compound_type` field

**Files:**
- Modify: `src/charon/core/schema.py`
- Modify: `tests/unit/test_schema_distribution.py`

- [ ] **Step 2.1: Append failing tests for `compound_type`**

Append to `tests/unit/test_schema_distribution.py`:

```python
class TestPhysicochemicalCompoundType:
    def test_default_is_none(self):
        pc = PhysicochemicalProperties()
        assert pc.compound_type is None

    def test_all_literal_values_accepted(self):
        for ct in ("neutral", "acid", "base", "zwitterion"):
            pc = PhysicochemicalProperties(compound_type=ct)
            assert pc.compound_type == ct

    def test_invalid_value_rejected(self):
        with pytest.raises(ValidationError):
            PhysicochemicalProperties(compound_type="polymer")  # type: ignore[arg-type]

    def test_coexists_with_logp_pka(self):
        pc = PhysicochemicalProperties(
            logp=PredictedProperty(value=3.89, source="experimental"),
            pka_base=PredictedProperty(value=6.2, source="experimental"),
            compound_type="base",
        )
        assert pc.logp.value == 3.89
        assert pc.pka_base.value == 6.2
        assert pc.compound_type == "base"
```

- [ ] **Step 2.2: Run tests to verify they fail**

```bash
pytest tests/unit/test_schema_distribution.py::TestPhysicochemicalCompoundType -v
```
Expected: FAIL with `ValidationError: Extra inputs are not permitted` for `compound_type`.

- [ ] **Step 2.3: Add `compound_type` field to `PhysicochemicalProperties`**

Locate `PhysicochemicalProperties` in `src/charon/core/schema.py` and modify:

```python
class PhysicochemicalProperties(BaseModel):
    """Lipophilicity, pKa, solubility."""

    logp: PredictedProperty | None = None
    pka_acid: PredictedProperty | None = None
    pka_base: PredictedProperty | None = None
    solubility_ug_ml: PredictedProperty | None = None
    compound_type: Literal[
        "neutral", "acid", "base", "zwitterion"
    ] | None = None
```

- [ ] **Step 2.4: Run tests to verify they pass**

```bash
pytest tests/unit/test_schema_distribution.py::TestPhysicochemicalCompoundType -v
```
Expected: 4 tests PASS.

- [ ] **Step 2.5: Run full test suite to verify no regression**

```bash
pytest tests/ -q 2>&1 | tail -10
```
Expected: baseline +12 (8 from Task 1 + 4 new).

- [ ] **Step 2.6: Commit**

```bash
git add src/charon/core/schema.py tests/unit/test_schema_distribution.py
git commit -m "$(cat <<'EOF'
Sprint 3b: Add compound_type field to PhysicochemicalProperties

Additive Literal["neutral", "acid", "base", "zwitterion"] | None field
on PhysicochemicalProperties. Enables YAML-level compound classification
override so validation compound YAMLs can pin midazolam as "base"
(pKa_base=6.2 is below infer_compound_type's >8.0 threshold) without
requiring every caller to pass a Pipeline kwarg. Not yet consumed by
build_compound_pbpk_params — wired in Task 6.

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 3: Wire `distribution` into `CompoundProperties` + backward-compat round-trip

**Files:**
- Modify: `src/charon/core/schema.py`
- Modify: `tests/unit/test_schema_distribution.py`

- [ ] **Step 3.1: Append failing backward-compat + wire-up tests**

Append to `tests/unit/test_schema_distribution.py`:

```python
class TestCompoundPropertiesDistributionField:
    def test_default_distribution_present(self):
        cp = CompoundProperties()
        assert cp.distribution is not None
        assert isinstance(cp.distribution, DistributionProperties)
        assert cp.distribution.empirical_kp_by_tissue is None

    def test_backward_compat_no_distribution_in_dict(self):
        """Existing Sprint 3a fixtures (no 'distribution' key) must parse."""
        cp = CompoundProperties.model_validate({
            "physicochemical": {
                "logp": {"value": -0.02, "source": "experimental"},
            },
            "binding": {
                "fu_p": {"value": 0.6, "source": "experimental"},
            },
        })
        assert cp.distribution.empirical_kp_by_tissue is None

    def test_explicit_distribution_round_trip(self):
        cp = CompoundProperties.model_validate({
            "distribution": {
                "empirical_kp_by_tissue": {
                    "adipose": {"value": 12.0, "source": "literature",
                                "method": "Björkman 2001"},
                },
            },
        })
        assert cp.distribution.empirical_kp_by_tissue is not None
        adipose = cp.distribution.empirical_kp_by_tissue["adipose"]
        assert adipose.value == 12.0
        assert adipose.source == "literature"
        assert adipose.method == "Björkman 2001"

    def test_yaml_like_dump_reload(self):
        """Serialize → deserialize → identical (audit trail preserved)."""
        original = CompoundProperties(
            distribution=DistributionProperties(
                empirical_kp_by_tissue={
                    "adipose": PredictedProperty(
                        value=15.0, source="literature", method="test"
                    ),
                }
            )
        )
        dumped = original.model_dump()
        reloaded = CompoundProperties.model_validate(dumped)
        assert reloaded.distribution.empirical_kp_by_tissue is not None
        assert reloaded.distribution.empirical_kp_by_tissue["adipose"].value == 15.0
```

- [ ] **Step 3.2: Run tests to verify they fail**

```bash
pytest tests/unit/test_schema_distribution.py::TestCompoundPropertiesDistributionField -v
```
Expected: FAIL — `distribution` is not a field on `CompoundProperties`.

- [ ] **Step 3.3: Add `distribution` to `CompoundProperties`**

In `src/charon/core/schema.py`, modify:

```python
class CompoundProperties(BaseModel):
    """Container for all predicted/measured compound properties."""

    physicochemical: PhysicochemicalProperties = PhysicochemicalProperties()
    permeability: PermeabilityProperties = PermeabilityProperties()
    binding: BindingProperties = BindingProperties()
    metabolism: MetabolismProperties = MetabolismProperties()
    safety: SafetyProperties = SafetyProperties()
    renal: RenalProperties = RenalProperties()
    distribution: DistributionProperties = DistributionProperties()
```

- [ ] **Step 3.4: Run tests to verify they pass**

```bash
pytest tests/unit/test_schema_distribution.py -v
```
Expected: 16 tests PASS (8 + 4 + 4).

- [ ] **Step 3.5: Run full test suite to verify Sprint 3a fixtures still load**

```bash
pytest tests/ -q 2>&1 | tail -10
```
Expected: baseline +16. Critically, `tests/unit/test_pipeline.py` (which constructs in-Python CompoundConfig fixtures for theophylline/midazolam) must still pass — the new `distribution` field defaults correctly.

- [ ] **Step 3.6: Commit**

```bash
git add src/charon/core/schema.py tests/unit/test_schema_distribution.py
git commit -m "$(cat <<'EOF'
Sprint 3b: Wire distribution field into CompoundProperties

CompoundProperties.distribution: DistributionProperties with a
default_factory so existing Sprint 3a fixtures (which have no
'distribution' key) continue to parse unchanged. Round-trip test
verifies audit fields (source, method) survive serialize/deserialize.
Schema side of the override path is now complete; ode_compiler wiring
lands in Task 7.

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 4: Add `KpOverrideRecord` dataclass + `CompoundPBPKParams.kp_overrides` field

**Files:**
- Modify: `src/charon/pbpk/ode_compiler.py`
- Create: `tests/unit/test_kp_override.py`

- [ ] **Step 4.1: Write failing tests for the dataclass**

Create `tests/unit/test_kp_override.py`:

```python
"""Tests for the Kp empirical override path in ode_compiler.

These tests cover: the KpOverrideRecord dataclass, the species guard,
compound_type precedence, the empirical override loop, and the
end-to-end interaction with CompoundPBPKParams.
"""

from __future__ import annotations

import pytest

from charon.core.parameter_bridge import ParameterBridge
from charon.core.schema import (
    BindingProperties,
    CompoundConfig,
    CompoundProperties,
    DistributionProperties,
    MetabolismProperties,
    PhysicochemicalProperties,
    PredictedProperty,
    RenalProperties,
)
from charon.pbpk.ode_compiler import (
    CompoundPBPKParams,
    KpOverrideRecord,
    build_compound_pbpk_params,
)
from charon.pbpk.topology import load_species_topology


class TestKpOverrideRecord:
    def test_construct_minimal(self):
        rec = KpOverrideRecord(
            tissue="adipose",
            rr_value=50.0,
            empirical_value=10.0,
            source="literature",
            method="Björkman 2001",
            flag=None,
        )
        assert rec.tissue == "adipose"
        assert rec.rr_value == 50.0
        assert rec.empirical_value == 10.0
        assert rec.source == "literature"
        assert rec.method == "Björkman 2001"
        assert rec.flag is None

    def test_frozen(self):
        rec = KpOverrideRecord(
            tissue="adipose",
            rr_value=50.0,
            empirical_value=10.0,
            source="literature",
            method=None,
            flag=None,
        )
        with pytest.raises(Exception):  # FrozenInstanceError
            rec.tissue = "muscle"  # type: ignore[misc]


class TestCompoundPBPKParamsKpOverridesField:
    def test_default_empty_tuple(self):
        """Existing call sites that don't set kp_overrides keep working."""
        params = CompoundPBPKParams(
            name="theophylline",
            molecular_weight=180.17,
            logp=-0.02,
            pka_acid=None,
            pka_base=None,
            compound_type="neutral",
            fu_p=0.6,
            bp_ratio=0.85,
            fu_b=0.6 / 0.85,
            clint_liver_L_h=1.0,
            cl_renal_L_h=0.1,
            kp_by_tissue={"adipose": 0.5},
        )
        assert params.kp_overrides == ()
```

- [ ] **Step 4.2: Run tests to verify they fail**

```bash
pytest tests/unit/test_kp_override.py -v
```
Expected: FAIL with `ImportError: cannot import name 'KpOverrideRecord' from 'charon.pbpk.ode_compiler'`.

- [ ] **Step 4.3: Add `KpOverrideRecord` dataclass and extend `CompoundPBPKParams`**

In `src/charon/pbpk/ode_compiler.py`, after the existing imports and before `_COMPOUND_TYPES`:

```python
@dataclass(frozen=True)
class KpOverrideRecord:
    """Audit record for a single tissue-level empirical Kp override."""

    tissue: str
    rr_value: float          # what R&R computed
    empirical_value: float   # what was injected
    source: str              # PredictedProperty.source
    method: str | None       # PredictedProperty.method (citation string)
    flag: str | None         # PredictedProperty.flag (optional)
```

Then modify `CompoundPBPKParams` to add the new field (must be last because dataclass default values follow non-default ones):

```python
@dataclass(frozen=True)
class CompoundPBPKParams:
    """Resolved, PBPK-ready compound parameters for a specific topology."""

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
```

- [ ] **Step 4.4: Run tests to verify they pass**

```bash
pytest tests/unit/test_kp_override.py::TestKpOverrideRecord tests/unit/test_kp_override.py::TestCompoundPBPKParamsKpOverridesField -v
```
Expected: 3 tests PASS.

- [ ] **Step 4.5: Run full test suite — check for regression**

```bash
pytest tests/ -q 2>&1 | tail -10
```
Expected: baseline +19. Critically, `tests/unit/test_ode_compiler.py` and `tests/unit/test_pipeline.py` must still pass — the new `kp_overrides` field has a default so existing callers are unaffected.

- [ ] **Step 4.6: Commit**

```bash
git add src/charon/pbpk/ode_compiler.py tests/unit/test_kp_override.py
git commit -m "$(cat <<'EOF'
Sprint 3b: Add KpOverrideRecord + CompoundPBPKParams.kp_overrides field

Additive dataclass + frozen field for auditing per-tissue empirical Kp
overrides. Default empty tuple preserves existing call sites. Override
loop and Pipeline metadata exposure follow in Tasks 7 and 8.

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 5: Add species guard to `build_compound_pbpk_params`

**Files:**
- Modify: `src/charon/pbpk/ode_compiler.py`
- Modify: `tests/unit/test_kp_override.py`

- [ ] **Step 5.1: Append failing species-guard test**

Append to `tests/unit/test_kp_override.py`:

```python
from collections import OrderedDict
from charon.pbpk.kp_calculator import TissueComposition
from charon.pbpk.topology import PBPKTopology, TissueNode


def _make_fake_rat_topology() -> PBPKTopology:
    """Construct a minimal valid rat topology for the species-guard test.

    Real rat.yaml is empty in Sprint 3a; we construct the dataclass
    directly to exercise the species!='human' guard without touching
    disk.
    """
    comp = TissueComposition(fn=0.01, fp=0.005, fw=0.8, pH=7.0)
    tissues = OrderedDict([
        ("lung", TissueNode(name="lung", volume_L=0.002,
                            blood_flow_L_h=5.4, composition=comp,
                            drains_to="arterial")),
        ("liver", TissueNode(name="liver", volume_L=0.0095,
                             blood_flow_L_h=1.0, composition=comp,
                             drains_to="venous")),
        ("kidney", TissueNode(name="kidney", volume_L=0.002,
                              blood_flow_L_h=0.8, composition=comp,
                              drains_to="venous")),
        ("gut_wall", TissueNode(name="gut_wall", volume_L=0.01,
                                blood_flow_L_h=0.4, composition=comp,
                                drains_to="liver")),
        ("spleen", TissueNode(name="spleen", volume_L=0.0007,
                              blood_flow_L_h=0.2, composition=comp,
                              drains_to="liver")),
        ("pancreas", TissueNode(name="pancreas", volume_L=0.0008,
                                blood_flow_L_h=0.05, composition=comp,
                                drains_to="liver")),
    ])
    plasma = TissueComposition(fn=0.0023, fp=0.00163, fw=0.945, pH=7.4)
    return PBPKTopology(
        species="rat",
        body_weight_kg=0.25,
        cardiac_output_L_h=5.4,
        hematocrit=0.46,
        venous_volume_L=0.013,
        arterial_volume_L=0.0054,
        hepatic_artery_L_h=0.35,
        tissues=tissues,
        plasma_composition=plasma,
        gfr_mL_min=1.31,
        liver_weight_g=9.5,
    )


def _make_minimal_compound() -> CompoundConfig:
    return CompoundConfig(
        name="test_compound",
        smiles="CCO",
        molecular_weight=46.07,
        source="experimental",
        properties=CompoundProperties(
            physicochemical=PhysicochemicalProperties(
                logp=PredictedProperty(value=-0.31, source="experimental"),
            ),
            binding=BindingProperties(
                fu_p=PredictedProperty(value=1.0, source="experimental", unit="fraction"),
                fu_inc=PredictedProperty(value=1.0, source="experimental", unit="fraction"),
                bp_ratio=PredictedProperty(value=1.0, source="experimental", unit="ratio"),
            ),
            metabolism=MetabolismProperties(
                clint_uL_min_mg=PredictedProperty(value=1.0, source="experimental", unit="uL/min/mg"),
            ),
            renal=RenalProperties(
                clrenal_L_h=PredictedProperty(value=0.0, source="experimental", unit="L/h"),
            ),
        ),
    )


class TestSpeciesGuard:
    def test_human_topology_accepted(self):
        topo = load_species_topology("human")
        compound = _make_minimal_compound()
        bridge = ParameterBridge()
        # Should NOT raise
        params = build_compound_pbpk_params(compound, topo, bridge=bridge)
        assert params.compound_type == "neutral"

    def test_rat_topology_rejected(self):
        topo = _make_fake_rat_topology()
        compound = _make_minimal_compound()
        bridge = ParameterBridge()
        with pytest.raises(NotImplementedError, match="species='rat'"):
            build_compound_pbpk_params(compound, topo, bridge=bridge)

    def test_error_message_mentions_sprint_4(self):
        topo = _make_fake_rat_topology()
        compound = _make_minimal_compound()
        bridge = ParameterBridge()
        with pytest.raises(NotImplementedError, match="Sprint 4"):
            build_compound_pbpk_params(compound, topo, bridge=bridge)
```

- [ ] **Step 5.2: Run tests to verify they fail**

```bash
pytest tests/unit/test_kp_override.py::TestSpeciesGuard -v
```
Expected: FAIL. `test_human_topology_accepted` passes (no guard yet), `test_rat_topology_rejected` and `test_error_message_mentions_sprint_4` fail because the rat topology is accepted and returns normally (or fails with a different error later).

- [ ] **Step 5.3: Add species guard at the top of `build_compound_pbpk_params`**

In `src/charon/pbpk/ode_compiler.py`, modify `build_compound_pbpk_params` — add the guard immediately after the signature and before `props = compound.properties`:

```python
def build_compound_pbpk_params(
    compound: CompoundConfig,
    topology: PBPKTopology,
    bridge: ParameterBridge,
    *,
    compound_type: str | None = None,
    clint_system: Literal["HLM", "hepatocytes"] = "HLM",
    liver_model: str = "well_stirred",
    override_cl_renal_L_h: float | None = None,
) -> CompoundPBPKParams:
    """Assemble all PBPK ODE parameters for a given compound and species.

    [existing docstring]
    """
    if topology.species != "human":
        raise NotImplementedError(
            f"build_compound_pbpk_params currently hardcodes human "
            f"mppgl=40.0, hepatocellularity=120.0. "
            f"species={topology.species!r} requires species-aware values; "
            f"fix scheduled for Sprint 4 (translational layer)."
        )

    props = compound.properties
    # ... rest unchanged
```

- [ ] **Step 5.4: Run tests to verify they pass**

```bash
pytest tests/unit/test_kp_override.py::TestSpeciesGuard -v
```
Expected: 3 tests PASS.

- [ ] **Step 5.5: Run full test suite to verify no regression**

```bash
pytest tests/ -q 2>&1 | tail -10
```
Expected: baseline +22. Human-species callers (Sprint 3a tests) unaffected.

- [ ] **Step 5.6: Commit**

```bash
git add src/charon/pbpk/ode_compiler.py tests/unit/test_kp_override.py
git commit -m "$(cat <<'EOF'
Sprint 3b: Add species guard to build_compound_pbpk_params

Two-line defensive guard that raises NotImplementedError when called
with a non-human PBPKTopology. The hardcoded human mppgl=40 /
hepatocellularity=120 constants remain (fixing them requires
preclinical validation, out of scope) but non-human callers now hit a
loud wall with a clear Sprint 4 reference instead of silently
consuming wrong values.

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 6: Add `compound_type` resolution precedence

**Files:**
- Modify: `src/charon/pbpk/ode_compiler.py`
- Create: `tests/unit/test_compound_type_resolution.py`

- [ ] **Step 6.1: Write failing tests for precedence**

Create `tests/unit/test_compound_type_resolution.py`:

```python
"""Tests for compound_type resolution precedence.

Precedence: Pipeline kwarg > YAML field (PhysicochemicalProperties.compound_type)
            > pKa-inferred.
"""

from __future__ import annotations

import pytest

from charon.core.parameter_bridge import ParameterBridge
from charon.core.schema import (
    BindingProperties,
    CompoundConfig,
    CompoundProperties,
    MetabolismProperties,
    PhysicochemicalProperties,
    PredictedProperty,
    RenalProperties,
)
from charon.pbpk.ode_compiler import build_compound_pbpk_params
from charon.pbpk.topology import load_species_topology


def _midazolam_like(compound_type_in_yaml: str | None = None) -> CompoundConfig:
    """Midazolam-like compound: pKa_base=6.2 (neutral by infer_compound_type)."""
    pc_kwargs = dict(
        logp=PredictedProperty(value=3.89, source="experimental"),
        pka_base=PredictedProperty(value=6.2, source="experimental"),
    )
    if compound_type_in_yaml is not None:
        pc_kwargs["compound_type"] = compound_type_in_yaml  # type: ignore[assignment]
    return CompoundConfig(
        name="midazolam_like",
        smiles="Cc1ncc2n1-c1ccc(Cl)cc1C(c1ccccc1F)=NC2",
        molecular_weight=325.77,
        source="experimental",
        properties=CompoundProperties(
            physicochemical=PhysicochemicalProperties(**pc_kwargs),
            binding=BindingProperties(
                fu_p=PredictedProperty(value=0.03, source="experimental", unit="fraction"),
                fu_inc=PredictedProperty(value=0.96, source="experimental", unit="fraction"),
                bp_ratio=PredictedProperty(value=0.66, source="experimental", unit="ratio"),
            ),
            metabolism=MetabolismProperties(
                clint_uL_min_mg=PredictedProperty(value=93.0, source="experimental", unit="uL/min/mg"),
            ),
            renal=RenalProperties(
                clrenal_L_h=PredictedProperty(value=0.0, source="experimental", unit="L/h"),
            ),
        ),
    )


class TestCompoundTypePrecedence:
    def setup_method(self):
        self.topo = load_species_topology("human")
        self.bridge = ParameterBridge()

    def test_inferred_default_is_neutral_for_midazolam(self):
        """pKa_base=6.2 < 8.0 threshold → infer_compound_type returns 'neutral'."""
        compound = _midazolam_like(compound_type_in_yaml=None)
        params = build_compound_pbpk_params(compound, self.topo, bridge=self.bridge)
        assert params.compound_type == "neutral"

    def test_yaml_field_overrides_inference(self):
        """YAML sets compound_type='base' → overrides inferred 'neutral'."""
        compound = _midazolam_like(compound_type_in_yaml="base")
        params = build_compound_pbpk_params(compound, self.topo, bridge=self.bridge)
        assert params.compound_type == "base"

    def test_kwarg_overrides_yaml(self):
        """Pipeline kwarg 'acid' beats YAML 'base' (highest precedence)."""
        compound = _midazolam_like(compound_type_in_yaml="base")
        params = build_compound_pbpk_params(
            compound, self.topo, bridge=self.bridge,
            compound_type="acid",
        )
        assert params.compound_type == "acid"

    def test_kwarg_overrides_inference(self):
        """Pipeline kwarg beats inferred when no YAML field set."""
        compound = _midazolam_like(compound_type_in_yaml=None)
        params = build_compound_pbpk_params(
            compound, self.topo, bridge=self.bridge,
            compound_type="base",
        )
        assert params.compound_type == "base"

    def test_yaml_invalid_value_rejected_at_schema_layer(self):
        """Invalid compound_type is rejected by Pydantic before reaching builder."""
        from pydantic import ValidationError
        with pytest.raises(ValidationError):
            _midazolam_like(compound_type_in_yaml="polymer")  # type: ignore[arg-type]

    def test_kwarg_invalid_value_rejected_at_builder(self):
        compound = _midazolam_like(compound_type_in_yaml=None)
        with pytest.raises(ValueError, match="compound_type must be one of"):
            build_compound_pbpk_params(
                compound, self.topo, bridge=self.bridge,
                compound_type="polymer",
            )
```

- [ ] **Step 6.2: Run tests to verify they fail**

```bash
pytest tests/unit/test_compound_type_resolution.py -v
```
Expected: FAIL. `test_inferred_default_is_neutral_for_midazolam` passes (existing behavior), but `test_yaml_field_overrides_inference` fails — the YAML field is ignored because the builder doesn't read it yet.

- [ ] **Step 6.3: Implement precedence in `build_compound_pbpk_params`**

In `src/charon/pbpk/ode_compiler.py`, locate the line:

```python
    resolved_type = compound_type or infer_compound_type(pka_acid, pka_base)
```

Replace with:

```python
    resolved_type = (
        compound_type
        or compound.properties.physicochemical.compound_type
        or infer_compound_type(pka_acid, pka_base)
    )
```

- [ ] **Step 6.4: Run tests to verify they pass**

```bash
pytest tests/unit/test_compound_type_resolution.py -v
```
Expected: 6 tests PASS.

- [ ] **Step 6.5: Run full test suite to verify no regression**

```bash
pytest tests/ -q 2>&1 | tail -10
```
Expected: baseline +28. Sprint 3a `TestPipelineMidazolamLimitation` (which passes `compound_type_override="base"` via Pipeline kwarg) still works because kwarg has highest precedence.

- [ ] **Step 6.6: Commit**

```bash
git add src/charon/pbpk/ode_compiler.py tests/unit/test_compound_type_resolution.py
git commit -m "$(cat <<'EOF'
Sprint 3b: Add compound_type resolution precedence

Precedence: Pipeline kwarg > PhysicochemicalProperties.compound_type
(YAML) > infer_compound_type(pka). The kwarg path is preserved
unchanged, so Sprint 3a call sites that pass compound_type_override
continue to work. The YAML path lets validation compound YAMLs pin
their classification without caller changes.

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 7: Implement empirical Kp override loop

**Files:**
- Modify: `src/charon/pbpk/ode_compiler.py`
- Modify: `tests/unit/test_kp_override.py`

- [ ] **Step 7.1: Append failing override-loop tests**

Append to `tests/unit/test_kp_override.py`:

```python
class TestEmpiricalKpOverride:
    def setup_method(self):
        self.topo = load_species_topology("human")
        self.bridge = ParameterBridge()

    def _compound_with_override(
        self,
        overrides: dict[str, float] | None,
    ) -> CompoundConfig:
        dist = DistributionProperties(
            empirical_kp_by_tissue=None if overrides is None else {
                t: PredictedProperty(
                    value=v, source="literature", method="test-citation"
                )
                for t, v in overrides.items()
            }
        )
        return CompoundConfig(
            name="mz",
            smiles="Cc1ncc2n1-c1ccc(Cl)cc1C(c1ccccc1F)=NC2",
            molecular_weight=325.77,
            source="experimental",
            properties=CompoundProperties(
                physicochemical=PhysicochemicalProperties(
                    logp=PredictedProperty(value=3.89, source="experimental"),
                    pka_base=PredictedProperty(value=6.2, source="experimental"),
                    compound_type="base",
                ),
                binding=BindingProperties(
                    fu_p=PredictedProperty(value=0.03, source="experimental", unit="fraction"),
                    fu_inc=PredictedProperty(value=0.96, source="experimental", unit="fraction"),
                    bp_ratio=PredictedProperty(value=0.66, source="experimental", unit="ratio"),
                ),
                metabolism=MetabolismProperties(
                    clint_uL_min_mg=PredictedProperty(value=93.0, source="experimental", unit="uL/min/mg"),
                ),
                renal=RenalProperties(
                    clrenal_L_h=PredictedProperty(value=0.0, source="experimental", unit="L/h"),
                ),
                distribution=dist,
            ),
        )

    def test_no_override_path(self):
        compound = self._compound_with_override(None)
        params = build_compound_pbpk_params(compound, self.topo, bridge=self.bridge)
        assert params.kp_overrides == ()

    def test_single_tissue_override_replaces_value(self):
        compound_no = self._compound_with_override(None)
        compound_yes = self._compound_with_override({"adipose": 10.0})
        params_no = build_compound_pbpk_params(compound_no, self.topo, bridge=self.bridge)
        params_yes = build_compound_pbpk_params(compound_yes, self.topo, bridge=self.bridge)

        # All non-adipose tissues unchanged
        for tissue in params_no.kp_by_tissue:
            if tissue == "adipose":
                continue
            assert params_yes.kp_by_tissue[tissue] == pytest.approx(
                params_no.kp_by_tissue[tissue]
            )

        # Adipose is now 10.0
        assert params_yes.kp_by_tissue["adipose"] == 10.0

        # Override recorded
        assert len(params_yes.kp_overrides) == 1
        rec = params_yes.kp_overrides[0]
        assert rec.tissue == "adipose"
        assert rec.rr_value == pytest.approx(params_no.kp_by_tissue["adipose"])
        assert rec.empirical_value == 10.0
        assert rec.source == "literature"
        assert rec.method == "test-citation"

    def test_multi_tissue_override(self):
        compound = self._compound_with_override({
            "adipose": 8.0,
            "muscle": 2.5,
        })
        params = build_compound_pbpk_params(compound, self.topo, bridge=self.bridge)
        assert params.kp_by_tissue["adipose"] == 8.0
        assert params.kp_by_tissue["muscle"] == 2.5
        assert len(params.kp_overrides) == 2
        tissues_logged = {r.tissue for r in params.kp_overrides}
        assert tissues_logged == {"adipose", "muscle"}

    def test_unknown_tissue_raises_with_valid_list(self):
        compound = self._compound_with_override({"nonsense_tissue": 5.0})
        with pytest.raises(ValueError) as exc_info:
            build_compound_pbpk_params(compound, self.topo, bridge=self.bridge)
        msg = str(exc_info.value)
        assert "'nonsense_tissue'" in msg
        assert "'human'" in msg
        assert "Valid tissues:" in msg
        # Error message lists the topology's actual tissues
        assert "'adipose'" in msg

    def test_override_propagates_to_ode_mass_balance(self):
        """The override must actually change ODE output (Vss drops)."""
        from charon.pbpk.solver import simulate_iv

        topo = self.topo
        bridge = self.bridge

        compound_no = self._compound_with_override(None)
        compound_yes = self._compound_with_override({"adipose": 10.0})

        params_no = build_compound_pbpk_params(compound_no, topo, bridge=bridge)
        params_yes = build_compound_pbpk_params(compound_yes, topo, bridge=bridge)

        sim_no = simulate_iv(topo, params_no, dose_mg=5.0,
                             route="iv_bolus", duration_h=168.0)
        sim_yes = simulate_iv(topo, params_yes, dose_mg=5.0,
                              route="iv_bolus", duration_h=168.0)

        # Both solvers must succeed
        assert sim_no.solver_success
        assert sim_yes.solver_success

        # Override should reduce the adipose contribution to Vss →
        # venous Cp late-time should be HIGHER with override (less drug
        # distributed into adipose → more in blood).
        assert sim_yes.cp_plasma[-1] > sim_no.cp_plasma[-1]
```

- [ ] **Step 7.2: Run tests to verify they fail**

```bash
pytest tests/unit/test_kp_override.py::TestEmpiricalKpOverride -v
```
Expected: FAIL. Most tests fail because the override path is not implemented.

- [ ] **Step 7.3: Implement the override loop in `build_compound_pbpk_params`**

In `src/charon/pbpk/ode_compiler.py`, locate the block after `kp_by_tissue = compute_all_kp(...)` and before `# Hepatic IVIVE with full audit via ParameterBridge.`:

```python
    kp_by_tissue = compute_all_kp(
        logp=logp,
        pka=pka_for_kp,
        compound_type=resolved_type,
        tissue_compositions=tissue_comps,
        plasma_composition=topology.plasma_composition,
        method="rodgers_rowland",
    )

    # Apply empirical Kp overrides (if any).  This replaces R&R values
    # tissue-by-tissue with cited literature/experimental Kp; the ODE
    # is unaware of whether Kp came from R&R or from override.
    kp_by_tissue_mut = dict(kp_by_tissue)
    override_log: list[KpOverrideRecord] = []
    empirical = compound.properties.distribution.empirical_kp_by_tissue
    if empirical:
        valid_tissues = set(kp_by_tissue_mut.keys())
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
                    rr_value=float(kp_by_tissue_mut[tissue]),
                    empirical_value=float(p.value),
                    source=p.source,
                    method=p.method,
                    flag=p.flag,
                )
            )
            kp_by_tissue_mut[tissue] = float(p.value)
    kp_by_tissue = kp_by_tissue_mut

    # Hepatic IVIVE with full audit via ParameterBridge.
    hep = bridge.clint_to_clh(
        ...
```

Then modify the return statement to include the override log:

```python
    return CompoundPBPKParams(
        name=compound.name,
        molecular_weight=float(mw),
        logp=logp,
        pka_acid=pka_acid,
        pka_base=pka_base,
        compound_type=resolved_type,
        fu_p=fu_p,
        bp_ratio=bp_ratio,
        fu_b=fu_b,
        clint_liver_L_h=clint_liver_L_h,
        cl_renal_L_h=cl_renal_L_h,
        kp_by_tissue=dict(kp_by_tissue),
        kp_overrides=tuple(override_log),
    )
```

- [ ] **Step 7.4: Run tests to verify they pass**

```bash
pytest tests/unit/test_kp_override.py -v
```
Expected: All Task 4 + Task 5 + Task 7 tests PASS (12 total in this file).

- [ ] **Step 7.5: Run full test suite to verify no regression**

```bash
pytest tests/ -q 2>&1 | tail -10
```
Expected: baseline +37.

- [ ] **Step 7.6: Commit**

```bash
git add src/charon/pbpk/ode_compiler.py tests/unit/test_kp_override.py
git commit -m "$(cat <<'EOF'
Sprint 3b: Implement empirical Kp override loop in build_compound_pbpk_params

After compute_all_kp, iterate compound.properties.distribution.
empirical_kp_by_tissue and replace matching entries in kp_by_tissue.
Unknown tissue names raise ValueError with a sorted list of the
topology's actual tissue keys for clear diagnostics. Each override
produces a KpOverrideRecord capturing rr_value and empirical_value
alongside the source/method citation string. End-to-end test verifies
the override propagates through the ODE (adipose-override lowers the
late-time venous Cp differential as expected).

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 8: Expose `kp_overrides` in `Pipeline.run().metadata`

**Files:**
- Modify: `src/charon/pipeline.py`
- Modify: `tests/unit/test_pipeline.py`

- [ ] **Step 8.1: Append failing test to existing test_pipeline.py**

Append to `tests/unit/test_pipeline.py`:

```python
class TestPipelineMetadataKpOverrides:
    """Verify that kp_overrides flow from build_compound_pbpk_params
    into PipelineResult.metadata['kp_overrides']."""

    def _midazolam_with_override(self):
        from charon.core.schema import (
            BindingProperties, CompoundConfig, CompoundProperties,
            DistributionProperties, MetabolismProperties,
            PhysicochemicalProperties, PredictedProperty, RenalProperties,
        )
        return CompoundConfig(
            name="midazolam-test",
            smiles="Cc1ncc2n1-c1ccc(Cl)cc1C(c1ccccc1F)=NC2",
            molecular_weight=325.77,
            source="experimental",
            properties=CompoundProperties(
                physicochemical=PhysicochemicalProperties(
                    logp=PredictedProperty(value=3.89, source="experimental"),
                    pka_base=PredictedProperty(value=6.2, source="experimental"),
                    compound_type="base",
                ),
                binding=BindingProperties(
                    fu_p=PredictedProperty(value=0.03, source="experimental", unit="fraction"),
                    fu_inc=PredictedProperty(value=0.96, source="experimental", unit="fraction"),
                    bp_ratio=PredictedProperty(value=0.66, source="experimental", unit="ratio"),
                ),
                metabolism=MetabolismProperties(
                    clint_uL_min_mg=PredictedProperty(value=93.0, source="experimental", unit="uL/min/mg"),
                ),
                renal=RenalProperties(
                    clrenal_L_h=PredictedProperty(value=0.0, source="experimental", unit="L/h"),
                ),
                distribution=DistributionProperties(
                    empirical_kp_by_tissue={
                        "adipose": PredictedProperty(
                            value=10.0, source="literature",
                            method="test-citation"
                        ),
                    }
                ),
            ),
        )

    def test_metadata_contains_kp_overrides_list(self):
        from charon import Pipeline

        pipe = Pipeline(
            compound=self._midazolam_with_override(),
            route="iv_bolus",
            dose_mg=5.0,
            duration_h=168.0,
        )
        result = pipe.run()
        overrides = result.metadata["kp_overrides"]
        assert isinstance(overrides, list)
        assert len(overrides) == 1
        entry = overrides[0]
        assert entry["tissue"] == "adipose"
        assert entry["empirical_value"] == 10.0
        assert entry["source"] == "literature"
        assert entry["method"] == "test-citation"
        assert "rr_value" in entry  # R&R original preserved for audit

    def test_metadata_overrides_empty_for_no_override_compound(self):
        """Sprint 3a theophylline fixture (no distribution field) → empty list."""
        from charon import Pipeline
        from tests.unit.test_pipeline import _build_theophylline  # if exists

        # Fallback: construct locally if helper not available
        from charon.core.schema import (
            BindingProperties, CompoundConfig, CompoundProperties,
            MetabolismProperties, PhysicochemicalProperties,
            PredictedProperty, RenalProperties,
        )
        theo = CompoundConfig(
            name="theophylline",
            smiles="Cn1c(=O)c2[nH]cnc2n(C)c1=O",
            molecular_weight=180.17,
            source="experimental",
            properties=CompoundProperties(
                physicochemical=PhysicochemicalProperties(
                    logp=PredictedProperty(value=-0.02, source="experimental"),
                ),
                binding=BindingProperties(
                    fu_p=PredictedProperty(value=0.60, source="experimental", unit="fraction"),
                    fu_inc=PredictedProperty(value=1.0, source="experimental", unit="fraction"),
                    bp_ratio=PredictedProperty(value=0.85, source="experimental", unit="ratio"),
                ),
                metabolism=MetabolismProperties(
                    clint_uL_min_mg=PredictedProperty(value=1.8, source="experimental", unit="uL/min/mg"),
                ),
                renal=RenalProperties(
                    clrenal_L_h=PredictedProperty(value=0.1, source="experimental", unit="L/h"),
                ),
            ),
        )
        pipe = Pipeline(compound=theo, route="iv_bolus", dose_mg=100.0, duration_h=168.0)
        result = pipe.run()
        assert result.metadata["kp_overrides"] == []
```

Note: remove the `from tests.unit.test_pipeline import _build_theophylline` line if no such helper exists; the fallback inline construction suffices.

- [ ] **Step 8.2: Run tests to verify they fail**

```bash
pytest tests/unit/test_pipeline.py::TestPipelineMetadataKpOverrides -v
```
Expected: FAIL with `KeyError: 'kp_overrides'` — the Pipeline metadata dict does not expose this yet.

- [ ] **Step 8.3: Add `kp_overrides` to `Pipeline.run()` metadata**

In `src/charon/pipeline.py`, locate the `return PipelineResult(...)` at the end of `run()` and modify the `metadata` dict:

```python
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
```

- [ ] **Step 8.4: Run tests to verify they pass**

```bash
pytest tests/unit/test_pipeline.py::TestPipelineMetadataKpOverrides -v
```
Expected: 2 tests PASS.

- [ ] **Step 8.5: Run full test suite to verify no regression**

```bash
pytest tests/ -q 2>&1 | tail -10
```
Expected: baseline +39. Sprint 3a theophylline/midazolam limitation tests still pass — they don't inspect `metadata['kp_overrides']` so the new field is invisible to them.

- [ ] **Step 8.6: Commit**

```bash
git add src/charon/pipeline.py tests/unit/test_pipeline.py
git commit -m "$(cat <<'EOF'
Sprint 3b: Expose kp_overrides in Pipeline.run() metadata

PipelineResult.metadata['kp_overrides'] is a list of dicts (one per
applied override) serializing the KpOverrideRecord fields for
downstream report rendering. Empty list when no overrides apply,
preserving Sprint 3a behavior for compounds without
distribution.empirical_kp_by_tissue.

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 9: midazolam direction-of-effect mechanism test

**Files:**
- Modify: `tests/unit/test_pipeline.py`

This implements DoD §3 mechanism check: using synthetic override values (not literature), verify the override path end-to-end improves midazolam Vss prediction direction.

- [ ] **Step 9.1: Append mechanism test**

Append to `tests/unit/test_pipeline.py`:

```python
class TestPipelineMidazolamOverrideMechanism:
    """DoD §3 mechanism check: override path is wired end-to-end.

    Uses SYNTHETIC override values (source='experimental', method=
    'mechanism-test'). Does NOT assert a specific numerical fold-error
    target — only direction of effect and audit record shape.
    """

    def _base_midazolam(self, with_override: bool):
        from charon.core.schema import (
            BindingProperties, CompoundConfig, CompoundProperties,
            DistributionProperties, MetabolismProperties,
            PhysicochemicalProperties, PredictedProperty, RenalProperties,
        )
        dist = None
        if with_override:
            dist = DistributionProperties(
                empirical_kp_by_tissue={
                    "adipose": PredictedProperty(
                        value=10.0, source="experimental",
                        method="mechanism-test"
                    ),
                }
            )
        return CompoundConfig(
            name="midazolam",
            smiles="Cc1ncc2n1-c1ccc(Cl)cc1C(c1ccccc1F)=NC2",
            molecular_weight=325.77,
            source="experimental",
            properties=CompoundProperties(
                physicochemical=PhysicochemicalProperties(
                    logp=PredictedProperty(value=3.89, source="experimental"),
                    pka_base=PredictedProperty(value=6.2, source="experimental"),
                    compound_type="base",
                ),
                binding=BindingProperties(
                    fu_p=PredictedProperty(value=0.03, source="experimental", unit="fraction"),
                    fu_inc=PredictedProperty(value=0.96, source="experimental", unit="fraction"),
                    bp_ratio=PredictedProperty(value=0.66, source="experimental", unit="ratio"),
                ),
                metabolism=MetabolismProperties(
                    clint_uL_min_mg=PredictedProperty(value=93.0, source="experimental", unit="uL/min/mg"),
                ),
                renal=RenalProperties(
                    clrenal_L_h=PredictedProperty(value=0.0, source="experimental", unit="L/h"),
                ),
                distribution=dist or DistributionProperties(),
            ),
        )

    def test_override_reduces_vss(self):
        """Adipose Kp override reduces tissue uptake → Vss drops."""
        from charon import Pipeline

        OBSERVED_VSS = 66.0  # L (Obach 1999)

        pipe_no = Pipeline(
            compound=self._base_midazolam(with_override=False),
            route="iv_bolus", dose_mg=5.0, duration_h=168.0,
        )
        pipe_yes = Pipeline(
            compound=self._base_midazolam(with_override=True),
            route="iv_bolus", dose_mg=5.0, duration_h=168.0,
        )
        result_no = pipe_no.run()
        result_yes = pipe_yes.run()

        vss_no = result_no.pk_parameters.vss
        vss_yes = result_yes.pk_parameters.vss

        assert vss_no is not None and vss_no > 0
        assert vss_yes is not None and vss_yes > 0

        # Mechanism check: override must reduce the predicted Vss
        assert vss_yes < vss_no, (
            f"Override should reduce Vss (adipose Kp 50→10), "
            f"got vss_no={vss_no:.1f}, vss_yes={vss_yes:.1f}"
        )

        # Direction-of-effect check: override must move Vss TOWARD
        # the observed value (not past it in the wrong direction)
        fold_no = max(vss_no / OBSERVED_VSS, OBSERVED_VSS / vss_no)
        fold_yes = max(vss_yes / OBSERVED_VSS, OBSERVED_VSS / vss_yes)
        assert fold_yes <= fold_no, (
            f"Override must bring Vss closer to observed, "
            f"fold_no={fold_no:.2f}, fold_yes={fold_yes:.2f}"
        )

    def test_override_audit_shape(self):
        from charon import Pipeline

        pipe = Pipeline(
            compound=self._base_midazolam(with_override=True),
            route="iv_bolus", dose_mg=5.0, duration_h=168.0,
        )
        result = pipe.run()
        overrides = result.metadata["kp_overrides"]
        assert len(overrides) == 1
        assert overrides[0]["tissue"] == "adipose"
        assert overrides[0]["empirical_value"] == 10.0
        assert overrides[0]["source"] == "experimental"
        assert overrides[0]["method"] == "mechanism-test"
        # rr_value is whatever R&R computed; just check it's positive
        assert overrides[0]["rr_value"] > 0

    def test_no_override_matches_existing_limitation_case(self):
        """Without override, the existing Sprint 3a limitation still reproduces."""
        from charon import Pipeline

        pipe = Pipeline(
            compound=self._base_midazolam(with_override=False),
            route="iv_bolus", dose_mg=5.0, duration_h=168.0,
        )
        result = pipe.run()
        # Sprint 3a documented that no-override midazolam is within 4x
        # of observed CL (relaxed gate). This test asserts the same
        # loose bound to confirm no unrelated regression.
        OBSERVED_CL = 21.0
        cl = result.pk_parameters.cl_apparent
        assert cl is not None and cl > 0
        fold_cl = max(cl / OBSERVED_CL, OBSERVED_CL / cl)
        assert fold_cl < 4.0, f"No-override midazolam CL fold={fold_cl:.2f}"
```

- [ ] **Step 9.2: Run tests to verify they pass**

```bash
pytest tests/unit/test_pipeline.py::TestPipelineMidazolamOverrideMechanism -v
```
Expected: 3 tests PASS. All rely on machinery landed in Tasks 1-8.

- [ ] **Step 9.3: Run full test suite to verify no regression**

```bash
pytest tests/ -q 2>&1 | tail -10
```
Expected: baseline +42.

- [ ] **Step 9.4: Commit**

```bash
git add tests/unit/test_pipeline.py
git commit -m "$(cat <<'EOF'
Sprint 3b: Add midazolam override direction-of-effect mechanism test

DoD §3 mechanism check. Uses synthetic override (Kp_adipose=10.0,
source='experimental', method='mechanism-test') to verify (1) the
override path reduces predicted Vss, (2) it moves Vss TOWARD the
observed 66 L rather than past it, and (3) the audit metadata is
populated correctly. Does NOT assert a literature-specific numerical
target — that responsibility belongs to the benchmark harness, which
uses the verified-citation override values from the midazolam YAML.

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 10: Benchmark — panel dataclasses and loader

**Files:**
- Modify: `validation/benchmarks/layer2_human_pk.py` (incremental; final rewrite in Task 14)
- Create: `tests/unit/test_panel_loader.py`
- Create: `validation/data/tier1_obach/README.md` (stub)
- Create: `validation/data/tier1_obach/compounds/` directory

- [ ] **Step 10.1: Create `validation/data/tier1_obach/README.md` stub**

Create `validation/data/tier1_obach/README.md`:

```markdown
# Obach 1999 Tier-1 Human IV PK Validation Panel

Loaded by `validation/benchmarks/layer2_human_pk.py`.

**Source:** Obach RS, *Drug Metab Dispos* 27(11):1350-1359, 1999.
"Prediction of human clearance of twenty-nine drugs from hepatic
microsomal intrinsic clearance data".

This directory contains:

- `panel.yaml` — panel metadata and observed PK values (Table 6)
- `compounds/*.yaml` — per-compound `CompoundConfig` files with
  literature-sourced physicochemical / binding / metabolism / renal
  values

Each numeric value in every compound YAML cites its source (primarily
Obach 1999 Table 2 for in vitro values; Obach 1999 Table 6 for in vivo
reference PK). Override values (in compound `distribution.empirical_kp_by_tissue`)
require verified literature citations recorded in the `method` field —
see `docs/superpowers/specs/2026-04-10-sprint3b-kp-override-validation-design.md`
Section 7.2 for the citation protocol.
```

Then `mkdir -p validation/data/tier1_obach/compounds`.

```bash
mkdir -p /home/jam/Charon/validation/data/tier1_obach/compounds
```

- [ ] **Step 10.2: Write failing loader tests**

Create `tests/unit/test_panel_loader.py`:

```python
"""Tests for the Obach 1999 panel loader."""

from __future__ import annotations

import textwrap
from pathlib import Path

import pytest
import yaml

# These imports target the loader that will live in
# validation/benchmarks/layer2_human_pk.py after Task 10.
import sys
REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from validation.benchmarks.layer2_human_pk import (  # noqa: E402
    PanelEntry,
    load_panel,
)


@pytest.fixture
def tmp_panel(tmp_path: Path) -> Path:
    """Build a minimal self-contained panel dir on disk."""
    root = tmp_path / "tier1_obach"
    (root / "compounds").mkdir(parents=True)

    theo_yaml = textwrap.dedent("""
        name: theophylline
        smiles: "Cn1c(=O)c2[nH]cnc2n(C)c1=O"
        molecular_weight: 180.17
        source: experimental
        properties:
          physicochemical:
            logp: {value: -0.02, source: experimental, method: "Obach 1999 Table 2"}
            compound_type: neutral
          binding:
            fu_p: {value: 0.60, source: experimental, unit: fraction}
            fu_inc: {value: 1.0, source: experimental, unit: fraction}
            bp_ratio: {value: 0.85, source: experimental, unit: ratio}
          metabolism:
            clint_uL_min_mg: {value: 1.8, source: experimental, unit: uL/min/mg}
          renal:
            clrenal_L_h: {value: 0.1, source: experimental, unit: L/h}
    """).lstrip()
    (root / "compounds" / "theophylline.yaml").write_text(theo_yaml)

    panel_yaml = textwrap.dedent("""
        name: "Test Panel"
        source: "Obach 1999 (trimmed)"
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
            obach_table_row: 42
            strict_targets: true
            notes: "Primary gate"
    """).lstrip()
    (root / "panel.yaml").write_text(panel_yaml)

    return root / "panel.yaml"


class TestLoadPanel:
    def test_single_entry_loads(self, tmp_panel: Path):
        entries = load_panel(tmp_panel)
        assert len(entries) == 1
        e = entries[0]
        assert isinstance(e, PanelEntry)
        assert e.key == "theophylline"
        assert e.compound.name == "theophylline"
        assert e.route == "iv_bolus"
        assert e.dose_mg == 100.0
        assert e.duration_h == 168.0
        assert e.strict_targets is True
        assert e.observed == {"cl_L_h": 2.9, "vss_L": 35.0, "t_half_h": 8.0}

    def test_compound_file_relative_to_panel_dir(self, tmp_panel: Path):
        entries = load_panel(tmp_panel)
        # Compound fully loaded through Pydantic
        assert entries[0].compound.properties.physicochemical.logp.value == -0.02

    def test_missing_compound_file_raises(self, tmp_path: Path):
        root = tmp_path / "bad_panel"
        root.mkdir()
        panel_yaml = textwrap.dedent("""
            name: "Broken"
            source: "test"
            default_duration_h: 168.0
            compounds:
              - key: ghost
                compound_file: compounds/ghost.yaml
                route: iv_bolus
                dose_mg: 1.0
                duration_h: 168.0
                observed: {cl_L_h: 1.0, vss_L: 1.0, t_half_h: 1.0}
                strict_targets: false
        """).lstrip()
        (root / "panel.yaml").write_text(panel_yaml)
        with pytest.raises(FileNotFoundError, match="ghost.yaml"):
            load_panel(root / "panel.yaml")

    def test_default_duration_fallback(self, tmp_path: Path):
        """Entry without duration_h uses panel default_duration_h."""
        root = tmp_path / "panel"
        (root / "compounds").mkdir(parents=True)
        (root / "compounds" / "theo.yaml").write_text(textwrap.dedent("""
            name: theo
            smiles: "Cn1c(=O)c2[nH]cnc2n(C)c1=O"
            molecular_weight: 180.17
            source: experimental
            properties:
              physicochemical:
                logp: {value: -0.02, source: experimental}
              binding:
                fu_p: {value: 0.6, source: experimental, unit: fraction}
                fu_inc: {value: 1.0, source: experimental, unit: fraction}
                bp_ratio: {value: 0.85, source: experimental, unit: ratio}
              metabolism:
                clint_uL_min_mg: {value: 1.8, source: experimental, unit: uL/min/mg}
              renal:
                clrenal_L_h: {value: 0.1, source: experimental, unit: L/h}
        """).lstrip())
        (root / "panel.yaml").write_text(textwrap.dedent("""
            name: "Fallback"
            source: "test"
            default_duration_h: 72.0
            compounds:
              - key: theo
                compound_file: compounds/theo.yaml
                route: iv_bolus
                dose_mg: 100.0
                observed: {cl_L_h: 2.9, vss_L: 35.0, t_half_h: 8.0}
                strict_targets: false
        """).lstrip())
        entries = load_panel(root / "panel.yaml")
        assert entries[0].duration_h == 72.0

    def test_observed_keys_required(self, tmp_path: Path):
        """Missing observed.cl_L_h → clear error."""
        root = tmp_path / "panel"
        (root / "compounds").mkdir(parents=True)
        (root / "compounds" / "c.yaml").write_text(textwrap.dedent("""
            name: c
            smiles: "C"
            molecular_weight: 16.0
            source: experimental
            properties:
              physicochemical:
                logp: {value: 0.5, source: experimental}
              binding:
                fu_p: {value: 0.9, source: experimental, unit: fraction}
                fu_inc: {value: 1.0, source: experimental, unit: fraction}
                bp_ratio: {value: 1.0, source: experimental, unit: ratio}
              metabolism:
                clint_uL_min_mg: {value: 1.0, source: experimental, unit: uL/min/mg}
              renal:
                clrenal_L_h: {value: 0.0, source: experimental, unit: L/h}
        """).lstrip())
        (root / "panel.yaml").write_text(textwrap.dedent("""
            name: x
            source: x
            default_duration_h: 72.0
            compounds:
              - key: c
                compound_file: compounds/c.yaml
                route: iv_bolus
                dose_mg: 1.0
                observed: {vss_L: 1.0, t_half_h: 1.0}
                strict_targets: false
        """).lstrip())
        with pytest.raises(KeyError, match="cl_L_h"):
            load_panel(root / "panel.yaml")
```

- [ ] **Step 10.3: Run tests to verify they fail**

```bash
pytest tests/unit/test_panel_loader.py -v
```
Expected: FAIL with `ImportError: cannot import name 'PanelEntry'` from `validation.benchmarks.layer2_human_pk`.

- [ ] **Step 10.4: Add dataclasses and loader to `layer2_human_pk.py`**

We will incrementally transform `validation/benchmarks/layer2_human_pk.py`. For this task, we ADD new types and functions without removing the existing Sprint 3a `theophylline()` / `midazolam()` Python factories yet (they get removed in Task 14).

Insert at the top of `validation/benchmarks/layer2_human_pk.py`, after the existing imports:

```python
from dataclasses import dataclass
from pathlib import Path

import yaml

from charon.core.schema import CompoundConfig


@dataclass
class PanelEntry:
    key: str
    compound: CompoundConfig
    route: str
    dose_mg: float
    duration_h: float
    observed: dict[str, float]
    strict_targets: bool
    obach_table_row: int | None
    notes: str


def load_panel(panel_path: Path) -> list[PanelEntry]:
    """Load a panel.yaml file and all compound files it references.

    `panel.yaml` format::

        name: "..."
        source: "..."
        default_duration_h: 168.0
        compounds:
          - key: theophylline
            compound_file: compounds/theophylline.yaml
            route: iv_bolus
            dose_mg: 100.0
            duration_h: 168.0          # optional; falls back to default_duration_h
            observed:
              cl_L_h: 2.9
              vss_L: 35.0
              t_half_h: 8.0
            obach_table_row: 42        # optional
            strict_targets: true
            notes: "..."               # optional

    Returns a list of `PanelEntry` objects with fully-loaded
    `CompoundConfig` instances.
    """
    panel_path = Path(panel_path)
    with panel_path.open() as f:
        raw = yaml.safe_load(f)

    default_duration = float(raw.get("default_duration_h", 168.0))
    entries: list[PanelEntry] = []

    for idx, item in enumerate(raw["compounds"]):
        # Resolve compound file relative to the panel's directory
        compound_file = panel_path.parent / item["compound_file"]
        if not compound_file.exists():
            raise FileNotFoundError(
                f"panel.yaml entry #{idx} ({item.get('key', '?')}) "
                f"references missing compound file: {compound_file}"
            )
        with compound_file.open() as cf:
            compound_data = yaml.safe_load(cf)
        compound = CompoundConfig.model_validate(compound_data)

        # Observed PK keys are required (KeyError surfaces typos early)
        observed = {
            "cl_L_h": float(item["observed"]["cl_L_h"]),
            "vss_L": float(item["observed"]["vss_L"]),
            "t_half_h": float(item["observed"]["t_half_h"]),
        }

        entries.append(
            PanelEntry(
                key=str(item["key"]),
                compound=compound,
                route=str(item["route"]),
                dose_mg=float(item["dose_mg"]),
                duration_h=float(item.get("duration_h", default_duration)),
                observed=observed,
                strict_targets=bool(item["strict_targets"]),
                obach_table_row=item.get("obach_table_row"),
                notes=str(item.get("notes", "")),
            )
        )

    return entries
```

- [ ] **Step 10.5: Run tests to verify they pass**

```bash
pytest tests/unit/test_panel_loader.py -v
```
Expected: 5 tests PASS.

- [ ] **Step 10.6: Run full test suite to verify existing benchmark still works**

```bash
pytest tests/ -q 2>&1 | tail -10
python3 validation/benchmarks/layer2_human_pk.py 2>&1 | tail -10
```
Expected: test suite baseline +47. Benchmark script still exits 0 (existing main() is untouched; loader code is additive).

- [ ] **Step 10.7: Commit**

```bash
git add validation/benchmarks/layer2_human_pk.py tests/unit/test_panel_loader.py validation/data/tier1_obach/README.md
git commit -m "$(cat <<'EOF'
Sprint 3b: Add PanelEntry dataclass and load_panel() to benchmark harness

Additive: new types and loader function coexist with the existing
Sprint 3a Python factory functions (theophylline() / midazolam()),
which will be removed in Task 14 when the new main() lands. The
loader resolves compound_file paths relative to the panel directory,
supports default_duration_h fallback, and raises clear errors for
missing files or observed PK keys. Tests exercise all paths via
temporary panel fixtures.

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 11: Benchmark — `_without_kp_overrides` strip function

**Files:**
- Modify: `validation/benchmarks/layer2_human_pk.py`
- Modify: `tests/unit/test_panel_loader.py`

- [ ] **Step 11.1: Append failing strip function test**

Append to `tests/unit/test_panel_loader.py`:

```python
class TestWithoutKpOverrides:
    def test_strips_empirical_kp(self):
        from charon.core.schema import (
            CompoundProperties, DistributionProperties,
            PhysicochemicalProperties, PredictedProperty,
        )
        from validation.benchmarks.layer2_human_pk import _without_kp_overrides

        original = CompoundProperties(
            physicochemical=PhysicochemicalProperties(
                logp=PredictedProperty(value=3.89, source="experimental"),
            ),
            distribution=DistributionProperties(
                empirical_kp_by_tissue={
                    "adipose": PredictedProperty(
                        value=10.0, source="literature", method="test"
                    ),
                }
            ),
        )
        stripped = _without_kp_overrides(original)
        assert stripped.distribution.empirical_kp_by_tissue is None
        # Other fields preserved
        assert stripped.physicochemical.logp.value == 3.89

    def test_no_override_case_is_noop(self):
        from charon.core.schema import CompoundProperties, PhysicochemicalProperties, PredictedProperty
        from validation.benchmarks.layer2_human_pk import _without_kp_overrides

        original = CompoundProperties(
            physicochemical=PhysicochemicalProperties(
                logp=PredictedProperty(value=-0.02, source="experimental"),
            ),
        )
        stripped = _without_kp_overrides(original)
        assert stripped.distribution.empirical_kp_by_tissue is None
        assert stripped.physicochemical.logp.value == -0.02

    def test_stripped_distribution_is_same_type(self):
        """The strip must not degrade to a parent type."""
        from charon.core.schema import (
            CompoundProperties, DistributionProperties,
            PredictedProperty,
        )
        from validation.benchmarks.layer2_human_pk import _without_kp_overrides

        original = CompoundProperties(
            distribution=DistributionProperties(
                empirical_kp_by_tissue={
                    "adipose": PredictedProperty(
                        value=10.0, source="literature", method="test"
                    ),
                }
            ),
        )
        stripped = _without_kp_overrides(original)
        assert isinstance(stripped.distribution, DistributionProperties)
        assert stripped.distribution.empirical_kp_by_tissue is None
        # Note: `DistributionProperties` currently has only one field, so a
        # genuine forward-compat test (verifying that future sibling fields
        # survive the strip) is not authorable here. The implementation
        # uses nested model_copy to get this behavior for free when new
        # fields are added; reviewers should preserve that pattern.
```

- [ ] **Step 11.2: Run tests to verify they fail**

```bash
pytest tests/unit/test_panel_loader.py::TestWithoutKpOverrides -v
```
Expected: FAIL with `ImportError: cannot import name '_without_kp_overrides'`.

- [ ] **Step 11.3: Add `_without_kp_overrides` to `layer2_human_pk.py`**

Append to `validation/benchmarks/layer2_human_pk.py` (after the `load_panel` function):

```python
from charon.core.schema import CompoundProperties


def _without_kp_overrides(props: CompoundProperties) -> CompoundProperties:
    """Return a copy of props with empirical_kp_by_tissue cleared.

    Uses nested model_copy so other fields of DistributionProperties
    that may be added in future (e.g. Vss_pred, tissue-level fu) are
    preserved rather than silently reset to their defaults.
    """
    new_distribution = props.distribution.model_copy(
        update={"empirical_kp_by_tissue": None}
    )
    return props.model_copy(update={"distribution": new_distribution})
```

- [ ] **Step 11.4: Run tests to verify they pass**

```bash
pytest tests/unit/test_panel_loader.py::TestWithoutKpOverrides -v
```
Expected: 3 tests PASS.

- [ ] **Step 11.5: Run full test suite**

```bash
pytest tests/ -q 2>&1 | tail -10
```
Expected: baseline +50.

- [ ] **Step 11.6: Commit**

```bash
git add validation/benchmarks/layer2_human_pk.py tests/unit/test_panel_loader.py
git commit -m "$(cat <<'EOF'
Sprint 3b: Add _without_kp_overrides strip function to benchmark harness

Uses nested model_copy on props.distribution so any future fields
added to DistributionProperties survive the strip pass. This is the
helper the two-pass benchmark uses to produce the "no override"
baseline from a compound that carries overrides.

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 12: Benchmark — metrics aggregation

**Files:**
- Modify: `validation/benchmarks/layer2_human_pk.py`
- Create: `tests/unit/test_panel_metrics.py`

- [ ] **Step 12.1: Write failing aggregation tests**

Create `tests/unit/test_panel_metrics.py`:

```python
"""Tests for benchmark aggregation helpers."""

from __future__ import annotations

import math
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from validation.benchmarks.layer2_human_pk import (  # noqa: E402
    PanelRow,
    PanelSummary,
    aggregate_summary,
)


def _row(key: str, fold_cl: float, fold_vss: float, fold_thalf: float,
         strict: bool, override_tissues=None):
    """Convenience: synthesize a PanelRow from target fold errors.

    `observed` is always 1.0 for all metrics; `predicted` is chosen so
    that max(pred/obs, obs/pred) == target fold.
    """
    def pred_from_fold(f: float) -> float:
        # Always return pred > obs so fold = pred/obs = f
        return f
    predicted = {
        "cl_L_h": pred_from_fold(fold_cl),
        "vss_L": pred_from_fold(fold_vss),
        "t_half_h": pred_from_fold(fold_thalf),
    }
    observed = {"cl_L_h": 1.0, "vss_L": 1.0, "t_half_h": 1.0}
    return PanelRow(
        key=key,
        predicted=predicted,
        observed=observed,
        fold={"cl_L_h": fold_cl, "vss_L": fold_vss, "t_half_h": fold_thalf},
        pass_2_fold={
            "cl_L_h": fold_cl <= 2.0,
            "vss_L": fold_vss <= 2.0,
            "t_half_h": fold_thalf <= 2.0,
        },
        override_tissues=override_tissues or [],
        strict_targets=strict,
        mode="no_override",
    )


class TestAggregateSummary:
    def test_single_row(self):
        rows = [_row("theo", 1.0, 1.0, 1.0, strict=True)]
        s = aggregate_summary(rows, mode="no_override")
        assert s.n == 1
        assert s.mode == "no_override"
        assert s.aafe["cl_L_h"] == pytest.approx(1.0)
        assert s.aafe["vss_L"] == pytest.approx(1.0)
        assert s.aafe["t_half_h"] == pytest.approx(1.0)
        assert s.within_2_fold["cl_L_h"] == pytest.approx(1.0)
        assert s.strict_failures == 0

    def test_aafe_geometric_mean(self):
        rows = [
            _row("a", fold_cl=2.0, fold_vss=1.0, fold_thalf=1.0, strict=False),
            _row("b", fold_cl=8.0, fold_vss=1.0, fold_thalf=1.0, strict=False),
        ]
        s = aggregate_summary(rows, mode="no_override")
        # geometric mean of (2, 8) = 4
        assert s.aafe["cl_L_h"] == pytest.approx(4.0)

    def test_within_2_fold_fraction(self):
        rows = [
            _row("a", 1.5, 1.0, 1.0, strict=False),
            _row("b", 3.0, 1.0, 1.0, strict=False),
        ]
        s = aggregate_summary(rows, mode="no_override")
        assert s.within_2_fold["cl_L_h"] == pytest.approx(0.5)
        assert s.within_3_fold["cl_L_h"] == pytest.approx(1.0)

    def test_strict_failure_count(self):
        rows = [
            _row("a", 1.5, 1.0, 1.0, strict=True),     # PASS all metrics
            _row("b", 3.0, 1.0, 1.0, strict=True),     # FAIL cl_L_h fold
            _row("c", 1.0, 5.0, 1.0, strict=False),    # non-strict; doesn't count
        ]
        s = aggregate_summary(rows, mode="no_override")
        assert s.strict_failures == 1

    def test_mode_label_preserved(self):
        rows = [_row("a", 1.0, 1.0, 1.0, strict=False)]
        s_no = aggregate_summary(rows, mode="no_override")
        s_yes = aggregate_summary(rows, mode="with_override")
        assert s_no.mode == "no_override"
        assert s_yes.mode == "with_override"

    def test_empty_rows_returns_zero_n(self):
        s = aggregate_summary([], mode="no_override")
        assert s.n == 0
        # AAFE on empty is undefined; implementation returns NaN
        assert math.isnan(s.aafe["cl_L_h"])
```

- [ ] **Step 12.2: Run tests to verify they fail**

```bash
pytest tests/unit/test_panel_metrics.py -v
```
Expected: FAIL with `ImportError: cannot import name 'PanelRow'`.

- [ ] **Step 12.3: Add `PanelRow`, `PanelSummary`, and `aggregate_summary` to `layer2_human_pk.py`**

Append to `validation/benchmarks/layer2_human_pk.py`:

```python
import math

from validation.benchmarks.metrics import aafe, within_n_fold


@dataclass
class PanelRow:
    key: str
    predicted: dict[str, float]
    observed: dict[str, float]
    fold: dict[str, float]
    pass_2_fold: dict[str, bool]
    override_tissues: list[str]
    strict_targets: bool
    mode: str  # "no_override" | "with_override"


@dataclass
class PanelSummary:
    n: int
    mode: str
    aafe: dict[str, float]
    within_2_fold: dict[str, float]
    within_3_fold: dict[str, float]
    strict_failures: int


_METRICS = ("cl_L_h", "vss_L", "t_half_h")


def aggregate_summary(rows: list[PanelRow], mode: str) -> PanelSummary:
    """Compute panel-level AAFE and within_n_fold fractions from PanelRows."""
    n = len(rows)
    if n == 0:
        return PanelSummary(
            n=0,
            mode=mode,
            aafe={m: float("nan") for m in _METRICS},
            within_2_fold={m: float("nan") for m in _METRICS},
            within_3_fold={m: float("nan") for m in _METRICS},
            strict_failures=0,
        )

    aafe_by_metric: dict[str, float] = {}
    w2_by_metric: dict[str, float] = {}
    w3_by_metric: dict[str, float] = {}
    for metric in _METRICS:
        preds = [r.predicted[metric] for r in rows]
        obs = [r.observed[metric] for r in rows]
        aafe_by_metric[metric] = aafe(preds, obs)
        w2_by_metric[metric] = within_n_fold(preds, obs, n=2.0)
        w3_by_metric[metric] = within_n_fold(preds, obs, n=3.0)

    strict_failures = 0
    for r in rows:
        if not r.strict_targets:
            continue
        if not all(r.pass_2_fold[m] for m in _METRICS):
            strict_failures += 1

    return PanelSummary(
        n=n,
        mode=mode,
        aafe=aafe_by_metric,
        within_2_fold=w2_by_metric,
        within_3_fold=w3_by_metric,
        strict_failures=strict_failures,
    )
```

Note: the `import math` and `from validation.benchmarks.metrics import ...` lines go at the top of the file with the other imports if not already present.

- [ ] **Step 12.4: Run tests to verify they pass**

```bash
pytest tests/unit/test_panel_metrics.py -v
```
Expected: 6 tests PASS.

- [ ] **Step 12.5: Run full test suite**

```bash
pytest tests/ -q 2>&1 | tail -10
```
Expected: baseline +56.

- [ ] **Step 12.6: Commit**

```bash
git add validation/benchmarks/layer2_human_pk.py tests/unit/test_panel_metrics.py
git commit -m "$(cat <<'EOF'
Sprint 3b: Add PanelRow/PanelSummary dataclasses and aggregate_summary helper

aggregate_summary computes per-metric panel-level AAFE,
within-2-fold, within-3-fold fractions, and counts strict-target
compound failures for the benchmark harness. Empty-panel case
returns n=0 with NaN metrics. Wraps the existing
validation/benchmarks/metrics.py helpers (aafe, within_n_fold)
instead of duplicating them.

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 13: Curate 10 compound YAMLs from Obach 1999

This is the hand-curation task. It does not follow strict TDD because the deliverable is data, not code. Each YAML is authored from `Obach 1999 Table 2 (in vitro)` + `Table 6 (in vivo)` with the per-compound gotchas from spec §4.6 applied.

**Files:**
- Create: `validation/data/tier1_obach/compounds/theophylline.yaml`
- Create: `validation/data/tier1_obach/compounds/antipyrine.yaml`
- Create: `validation/data/tier1_obach/compounds/caffeine.yaml`
- Create: `validation/data/tier1_obach/compounds/warfarin.yaml`
- Create: `validation/data/tier1_obach/compounds/diclofenac.yaml`
- Create: `validation/data/tier1_obach/compounds/midazolam.yaml`
- Create: `validation/data/tier1_obach/compounds/propranolol.yaml`
- Create: `validation/data/tier1_obach/compounds/diazepam.yaml`
- Create: `validation/data/tier1_obach/compounds/metoprolol.yaml`
- Create: `validation/data/tier1_obach/compounds/verapamil.yaml`

**Pharmacological judgments to apply (from spec §4.6):**

1. **warfarin.yaml `name` field:** verify during curation whether Obach 1999 reports racemic, S-, or R-warfarin. Record the answer in the YAML `name` (e.g. `"warfarin (S-enantiomer)"`) and set `compound_type: acid`.

2. **diazepam.yaml `compound_type`:** pKa_base 3.4 is well below R&R's base threshold; set `compound_type: neutral` (diazepam is effectively neutral at pH 7.4). Include a YAML comment explaining the decision.

3. **midazolam.yaml `compound_type`:** pKa_base 6.2; set `compound_type: base` (Sprint 3a convention). This is a deliberate classification override carried forward from Sprint 3a. Phase 1 has NO distribution override (Phase 2 lands in Task 17).

4. **diclofenac.yaml:** expect a ParameterBridge "very low fu_p" warning (fu_p ≈ 0.005 < 0.01); add a YAML comment acknowledging it.

5. **All 10 YAMLs:** every `PredictedProperty` declares `source: experimental` and `method: "Obach 1999 Table 2"` (or equivalent primary source for values not in Obach Table 2).

- [ ] **Step 13.1: Author `theophylline.yaml` (mirror the Sprint 3a Python factory)**

Create `validation/data/tier1_obach/compounds/theophylline.yaml`:

```yaml
# Obach 1999 Tier-1 panel — theophylline
# Sprint 3a primary validation compound, carried into the Sprint 3b panel
# as the strict-gate regression anchor.
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
      method: "Austin 2002 (fu_inc≈1 for low logP)"
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

- [ ] **Step 13.2: Author the remaining 9 compound YAMLs using the template and data table below**

For each compound, create the file at `validation/data/tier1_obach/compounds/<name>.yaml` following the theophylline template structure. Use the values below. **Before committing, cross-check each row against a printed copy of Obach 1999 Tables 2 and 6 (or an accepted secondary source) and correct any discrepancies.**

> Important: the values in the table below are from common secondary sources and Obach 1999 canonical references but MUST be verified against the actual paper during implementation. Record discrepancies as per-compound commit notes.

| Compound | logP | pKa | compound_type | fu_p | fu_inc | bp_ratio | CLint (μL/min/mg) | CL_renal (L/h) | Chem note |
|---|---|---|---|---|---|---|---|---|---|
| antipyrine | 0.38 | — | neutral | 0.90 | 1.0 | 1.0 | 0.9 | 0.1 | neutral, classic probe |
| caffeine | -0.07 | — | neutral | 0.70 | 1.0 | 1.04 | 0.9 | 0.1 | neutral, CYP1A2 |
| warfarin | 3.54 | pka_acid=5.0 | acid | 0.012 | 0.17 | 0.56 | 3.1 | 0.05 | verify racemate vs S- |
| diclofenac | 4.02 | pka_acid=4.0 | acid | 0.005 | 0.17 | 0.55 | 11.0 | 0.1 | low fu_p edge case |
| midazolam | 3.89 | pka_base=6.2 | base | 0.03 | 0.96 | 0.66 | 93.0 | 0.0 | base override from Sprint 3a |
| propranolol | 3.48 | pka_base=9.5 | base | 0.13 | 0.96 | 0.80 | 13.6 | 0.1 | high ER base |
| diazepam | 2.82 | pka_base=3.4 | **neutral** | 0.013 | 0.67 | 0.58 | 0.37 | 0.0 | pKa too low for base branch |
| metoprolol | 1.88 | pka_base=9.6 | base | 0.88 | 1.0 | 1.14 | 24.0 | 0.7 | CYP2D6 population mean |
| verapamil | 3.79 | pka_base=8.9 | base | 0.09 | 0.68 | 0.77 | 140.0 | 0.0 | ER≈0.95 stress test |

**Template A — neutral compound (use for antipyrine, caffeine, diazepam).** No pKa key at all; `compound_type: neutral`:

```yaml
name: <compound_name>
smiles: "<canonical SMILES>"
molecular_weight: <MW>
source: experimental
properties:
  physicochemical:
    logp:
      value: <logP>
      source: experimental
      method: "Obach 1999 Table 2"
    compound_type: neutral
  binding:
    fu_p: {value: <fu_p>, source: experimental, unit: fraction, method: "Obach 1999 Table 2"}
    fu_inc: {value: <fu_inc>, source: experimental, unit: fraction, method: "Obach 1999 Table 2"}
    bp_ratio: {value: <bp>, source: experimental, unit: ratio, method: "Obach 1999 Table 2"}
  metabolism:
    clint_uL_min_mg: {value: <CLint>, source: experimental, unit: uL/min/mg, method: "Obach 1999 Table 2"}
  renal:
    clrenal_L_h: {value: <CL_renal>, source: experimental, unit: L/h, method: "Obach 1999 Table 2"}
```

**Template B — acid compound (use for warfarin, diclofenac).** Adds `pka_acid`; `compound_type: acid`:

```yaml
name: <compound_name>
smiles: "<canonical SMILES>"
molecular_weight: <MW>
source: experimental
properties:
  physicochemical:
    logp:
      value: <logP>
      source: experimental
      method: "Obach 1999 Table 2"
    pka_acid:
      value: <pKa>
      source: experimental
      method: "Literature (see commit)"
    compound_type: acid
  binding:
    fu_p: {value: <fu_p>, source: experimental, unit: fraction, method: "Obach 1999 Table 2"}
    fu_inc: {value: <fu_inc>, source: experimental, unit: fraction, method: "Obach 1999 Table 2"}
    bp_ratio: {value: <bp>, source: experimental, unit: ratio, method: "Obach 1999 Table 2"}
  metabolism:
    clint_uL_min_mg: {value: <CLint>, source: experimental, unit: uL/min/mg, method: "Obach 1999 Table 2"}
  renal:
    clrenal_L_h: {value: <CL_renal>, source: experimental, unit: L/h, method: "Obach 1999 Table 2"}
```

**Template C — base compound (use for midazolam, propranolol, metoprolol, verapamil).** Adds `pka_base`; `compound_type: base`:

```yaml
name: <compound_name>
smiles: "<canonical SMILES>"
molecular_weight: <MW>
source: experimental
properties:
  physicochemical:
    logp:
      value: <logP>
      source: experimental
      method: "Obach 1999 Table 2"
    pka_base:
      value: <pKa>
      source: experimental
      method: "Literature (see commit)"
    compound_type: base
  binding:
    fu_p: {value: <fu_p>, source: experimental, unit: fraction, method: "Obach 1999 Table 2"}
    fu_inc: {value: <fu_inc>, source: experimental, unit: fraction, method: "Obach 1999 Table 2"}
    bp_ratio: {value: <bp>, source: experimental, unit: ratio, method: "Obach 1999 Table 2"}
  metabolism:
    clint_uL_min_mg: {value: <CLint>, source: experimental, unit: uL/min/mg, method: "Obach 1999 Table 2"}
  renal:
    clrenal_L_h: {value: <CL_renal>, source: experimental, unit: L/h, method: "Obach 1999 Table 2"}
```

**Compound → template mapping:**
- antipyrine, caffeine → **Template A** (neutral, no pKa)
- warfarin, diclofenac → **Template B** (acid with `pka_acid`)
- midazolam, propranolol, metoprolol, verapamil → **Template C** (base with `pka_base`)
- **diazepam → Template D** (special case: pKa_base recorded for scientific
  completeness, but `compound_type: neutral` because the pKa is below
  R&R's base threshold — see spec §4.6 #2):

```yaml
name: diazepam
smiles: "CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc12"
molecular_weight: 284.74
source: experimental
properties:
  physicochemical:
    logp:
      value: 2.82
      source: experimental
      method: "Obach 1999 Table 2"
    pka_base:
      value: 3.4
      source: experimental
      method: "Literature (chemically a weak base)"
    compound_type: neutral  # pKa_base 3.4 is below R&R base threshold
                            # (≥7); diazepam is effectively neutral at
                            # pH 7.4, so the R&R neutral branch is used.
  binding:
    fu_p: {value: 0.013, source: experimental, unit: fraction, method: "Obach 1999 Table 2"}
    fu_inc: {value: 0.67, source: experimental, unit: fraction, method: "Obach 1999 Table 2"}
    bp_ratio: {value: 0.58, source: experimental, unit: ratio, method: "Obach 1999 Table 2"}
  metabolism:
    clint_uL_min_mg: {value: 0.37, source: experimental, unit: uL/min/mg, method: "Obach 1999 Table 2"}
  renal:
    clrenal_L_h: {value: 0.0, source: experimental, unit: L/h, method: "Obach 1999 Table 2"}
```

SMILES (canonical) and molecular weights:

- antipyrine: `CN1N(C(=O)C=C1C)c1ccccc1`, MW 188.23
- caffeine: `Cn1cnc2c1c(=O)n(C)c(=O)n2C`, MW 194.19
- warfarin: `CC(=O)CC(c1ccccc1)C1=C(O)c2ccccc2OC1=O`, MW 308.33
- diclofenac: `OC(=O)Cc1ccccc1Nc1c(Cl)cccc1Cl`, MW 296.15
- midazolam: `Cc1ncc2n1-c1ccc(Cl)cc1C(c1ccccc1F)=NC2`, MW 325.77
- propranolol: `CC(C)NCC(O)COc1cccc2ccccc12`, MW 259.34
- diazepam: `CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc12`, MW 284.74
- metoprolol: `COCCc1ccc(OCC(O)CNC(C)C)cc1`, MW 267.37
- verapamil: `COc1ccc(CC(C)(C#N)CCCN(C)CCc2ccc(OC)c(OC)c2)cc1OC`, MW 454.60

Special YAML annotations:

- **warfarin.yaml**: after `name:` add a YAML comment: `# TODO at curation: verify if Obach 1999 reports racemate or S-enantiomer; update name accordingly.`

  On verification, remove the TODO comment and update `name` to the specific form.

- **diazepam.yaml**: after `compound_type: neutral` add a YAML comment:

  ```yaml
    compound_type: neutral  # pKa_base 3.4 is too low for R&R base branch
                            # (threshold ≥7); diazepam behaves as neutral at
                            # pH 7.4, so R&R uses the neutral equations.
  ```

- **diclofenac.yaml**: after `fu_p:` add a YAML comment:

  ```yaml
    fu_p:
      value: 0.005  # WARNING: fu_p < 0.01 triggers the ParameterBridge
                    # low-fu_p warning per CLAUDE.md §6j Pitfall #10. Expected.
      source: experimental
      unit: fraction
      method: "Obach 1999 Table 2"
  ```

- **midazolam.yaml**: after `compound_type: base` add a YAML comment:

  ```yaml
    compound_type: base  # Sprint 3a convention: explicit override because
                         # infer_compound_type default threshold pKa_base>8.0
                         # would classify midazolam (pKa_base=6.2) as neutral.
                         # Keep in sync with tests/unit/test_pipeline.py
                         # TestPipelineMidazolamLimitation.
  ```

- [ ] **Step 13.3: Run a sanity check — load every YAML through Pydantic**

```bash
cd /home/jam/Charon
python3 - <<'EOF'
from pathlib import Path
import yaml
from charon.core.schema import CompoundConfig

root = Path("validation/data/tier1_obach/compounds")
for yml in sorted(root.glob("*.yaml")):
    with yml.open() as f:
        data = yaml.safe_load(f)
    cc = CompoundConfig.model_validate(data)
    print(f"OK  {cc.name:15}  MW={cc.molecular_weight}  type={cc.properties.physicochemical.compound_type}")
EOF
```
Expected: 10 "OK" lines, one per compound. Any Pydantic validation error → fix the offending YAML.

- [ ] **Step 13.4: Commit**

```bash
git add validation/data/tier1_obach/compounds/*.yaml
git commit -m "$(cat <<'EOF'
Sprint 3b: Add 10 Obach 1999 Tier-1 compound YAMLs

Hand-curated CompoundConfig files for the Sprint 3b validation panel:
neutral (theophylline, antipyrine, caffeine, diazepam), acid
(warfarin, diclofenac), base (midazolam, propranolol, metoprolol,
verapamil). Every PredictedProperty cites its source (Obach 1999
Table 2 for in vitro values). Per-compound pharmacological gotchas
(warfarin stereoisomer, diazepam neutral-for-R&R, diclofenac low fu_p,
midazolam base-override convention) are documented inline as YAML
comments. Phase 1 midazolam has no distribution.empirical_kp_by_tissue
override yet; Phase 2 lands after citation verification in Task 17.

Per spec §4.6: warfarin stereoisomer resolution happens at curation
time; verify the current YAML against the actual Obach 1999 paper
copy and adjust if needed.

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 14: Author `panel.yaml` and rewrite benchmark `main()`

**Files:**
- Create: `validation/data/tier1_obach/panel.yaml`
- Modify: `validation/benchmarks/layer2_human_pk.py` (remove old Python factories, install new main)
- Create: `tests/unit/test_panel_main.py`

- [ ] **Step 14.1: Author `panel.yaml`**

Create `validation/data/tier1_obach/panel.yaml`:

```yaml
# Obach 1999 Tier-1 Human IV PK Validation Panel
# Observed CL/Vss/t_half values from Obach RS, Drug Metab Dispos
# 27(11):1350-1359, 1999, Table 6 (in vivo human PK).  obach_table_row
# fields refer to Table 6 row numbers (verify against paper copy).
name: "Obach 1999 Tier-1 Human IV PK Validation Panel"
source: "Obach RS, Drug Metab Dispos 27(11):1350-1359, 1999"
default_duration_h: 168.0

compounds:
  - key: theophylline
    compound_file: compounds/theophylline.yaml
    route: iv_bolus
    dose_mg: 100.0
    duration_h: 168.0
    observed: {cl_L_h: 2.9, vss_L: 35.0, t_half_h: 8.0}
    obach_table_row: 42
    strict_targets: true
    notes: "Neutral low-logP primary validation compound."

  - key: antipyrine
    compound_file: compounds/antipyrine.yaml
    route: iv_bolus
    dose_mg: 500.0
    observed: {cl_L_h: 2.8, vss_L: 43.0, t_half_h: 11.0}
    obach_table_row: 3
    strict_targets: false
    notes: "Neutral phenotyping probe, low extraction."

  - key: caffeine
    compound_file: compounds/caffeine.yaml
    route: iv_bolus
    dose_mg: 250.0
    observed: {cl_L_h: 6.0, vss_L: 42.0, t_half_h: 4.9}
    obach_table_row: 6
    strict_targets: false
    notes: "Neutral, CYP1A2 substrate, low clearance."

  - key: warfarin
    compound_file: compounds/warfarin.yaml
    route: iv_bolus
    dose_mg: 15.0
    observed: {cl_L_h: 0.19, vss_L: 11.0, t_half_h: 37.0}
    obach_table_row: 43
    strict_targets: false
    notes: "Highly bound acid (fu_p 0.012), CYP2C9."

  - key: diclofenac
    compound_file: compounds/diclofenac.yaml
    route: iv_bolus
    dose_mg: 50.0
    observed: {cl_L_h: 16.0, vss_L: 13.0, t_half_h: 1.2}
    obach_table_row: 14
    strict_targets: false
    notes: "Acid, fu_p edge case (0.005), mid extraction."

  - key: midazolam
    compound_file: compounds/midazolam.yaml
    route: iv_bolus
    dose_mg: 5.0
    observed: {cl_L_h: 21.0, vss_L: 66.0, t_half_h: 3.0}
    obach_table_row: 24
    strict_targets: false
    notes: "Weak base, CYP3A4. R&R limitation; override demo target."

  - key: propranolol
    compound_file: compounds/propranolol.yaml
    route: iv_bolus
    dose_mg: 10.0
    observed: {cl_L_h: 50.0, vss_L: 270.0, t_half_h: 3.9}
    obach_table_row: 32
    strict_targets: false
    notes: "Base, CYP2D6, high extraction."

  - key: diazepam
    compound_file: compounds/diazepam.yaml
    route: iv_bolus
    dose_mg: 5.0
    observed: {cl_L_h: 1.6, vss_L: 77.0, t_half_h: 43.0}
    obach_table_row: 13
    strict_targets: false
    notes: "Weak base, neutral-for-R&R (pKa_base 3.4); long t_half."

  - key: metoprolol
    compound_file: compounds/metoprolol.yaml
    route: iv_bolus
    dose_mg: 10.0
    observed: {cl_L_h: 63.0, vss_L: 290.0, t_half_h: 4.0}
    obach_table_row: 23
    strict_targets: false
    notes: "Base, CYP2D6 population mean."

  - key: verapamil
    compound_file: compounds/verapamil.yaml
    route: iv_bolus
    dose_mg: 10.0
    observed: {cl_L_h: 60.0, vss_L: 350.0, t_half_h: 4.0}
    obach_table_row: 47
    strict_targets: false
    notes: "Base, CYP3A4, ER≈0.95 engine stress test."
```

> Observed values are canonical Obach 1999 references. Cross-check against the printed paper before accepting the panel as authoritative. If any discrepancy is found, correct the YAML and note it in the commit message.

- [ ] **Step 14.2: Write failing benchmark integration test**

Create `tests/unit/test_panel_main.py`:

```python
"""Smoke test: layer2_human_pk.main() runs over the real panel."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from validation.benchmarks.layer2_human_pk import (  # noqa: E402
    main,
    run_benchmark,
)


class TestBenchmarkMain:
    def test_default_panel_runs(self, capsys):
        """Run against the real panel.yaml — exit 0 if all strict pass."""
        exit_code = main()
        captured = capsys.readouterr()
        # Should print both sections
        assert "R&R only" in captured.out or "R&R only" in captured.err
        assert "with_override" in captured.out.lower() or "override" in captured.out.lower()
        # theophylline is the sole strict-gate compound — must pass
        assert exit_code == 0, f"benchmark failed, output:\n{captured.out}"

    def test_run_benchmark_returns_two_summaries(self):
        """Programmatic API returns (no_override, with_override) summaries."""
        panel_path = (
            REPO_ROOT / "validation" / "data" / "tier1_obach" / "panel.yaml"
        )
        summaries = run_benchmark(panel_path)
        assert "no_override" in summaries
        assert "with_override" in summaries
        assert summaries["no_override"].n == 10
        assert summaries["with_override"].n == 10
```

- [ ] **Step 14.3: Run tests to verify they fail**

```bash
pytest tests/unit/test_panel_main.py -v
```
Expected: FAIL with `ImportError: cannot import name 'run_benchmark'`.

- [ ] **Step 14.4: Rewrite `layer2_human_pk.py` main()**

REPLACE the contents of `validation/benchmarks/layer2_human_pk.py` with the following. Keep the PanelEntry / PanelRow / PanelSummary / load_panel / _without_kp_overrides / aggregate_summary code that was added in Tasks 10-12 (copy them into the new file). Delete the Sprint 3a Python factory functions `theophylline()`, `midazolam()`, and the old `_run_one()` function.

```python
"""Layer 2 human PBPK benchmark — Obach 1999 Tier-1 panel (10 compounds).

Loads ``validation/data/tier1_obach/panel.yaml``, runs each compound
twice (no-override and with-override), and prints per-metric panel
AAFE / within-2-fold / within-3-fold summaries. Exit 0 iff every
``strict_targets: true`` compound passes 2-fold on CL, Vss, and t½.

Run as a standalone script::

    python3 validation/benchmarks/layer2_human_pk.py
"""

from __future__ import annotations

import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from charon import Pipeline  # noqa: E402
from charon.core.schema import (  # noqa: E402
    CompoundConfig,
    CompoundProperties,
)
from validation.benchmarks.metrics import aafe, fold_error, within_n_fold  # noqa: E402


DEFAULT_PANEL_PATH = (
    REPO_ROOT / "validation" / "data" / "tier1_obach" / "panel.yaml"
)

_METRICS = ("cl_L_h", "vss_L", "t_half_h")


# ---------------------------------------------------------------------------
# Data types
# ---------------------------------------------------------------------------


@dataclass
class PanelEntry:
    key: str
    compound: CompoundConfig
    route: str
    dose_mg: float
    duration_h: float
    observed: dict[str, float]
    strict_targets: bool
    obach_table_row: int | None
    notes: str


@dataclass
class PanelRow:
    key: str
    predicted: dict[str, float]
    observed: dict[str, float]
    fold: dict[str, float]
    pass_2_fold: dict[str, bool]
    override_tissues: list[str]
    strict_targets: bool
    mode: str  # "no_override" | "with_override"


@dataclass
class PanelSummary:
    n: int
    mode: str
    aafe: dict[str, float]
    within_2_fold: dict[str, float]
    within_3_fold: dict[str, float]
    strict_failures: int


# ---------------------------------------------------------------------------
# Loader + strip helpers
# ---------------------------------------------------------------------------


def load_panel(panel_path: Path) -> list[PanelEntry]:
    """Load a panel.yaml file and all compound files it references."""
    panel_path = Path(panel_path)
    with panel_path.open() as f:
        raw = yaml.safe_load(f)

    default_duration = float(raw.get("default_duration_h", 168.0))
    entries: list[PanelEntry] = []

    for idx, item in enumerate(raw["compounds"]):
        compound_file = panel_path.parent / item["compound_file"]
        if not compound_file.exists():
            raise FileNotFoundError(
                f"panel.yaml entry #{idx} ({item.get('key', '?')}) "
                f"references missing compound file: {compound_file}"
            )
        with compound_file.open() as cf:
            compound_data = yaml.safe_load(cf)
        compound = CompoundConfig.model_validate(compound_data)

        observed = {
            "cl_L_h": float(item["observed"]["cl_L_h"]),
            "vss_L": float(item["observed"]["vss_L"]),
            "t_half_h": float(item["observed"]["t_half_h"]),
        }

        entries.append(
            PanelEntry(
                key=str(item["key"]),
                compound=compound,
                route=str(item["route"]),
                dose_mg=float(item["dose_mg"]),
                duration_h=float(item.get("duration_h", default_duration)),
                observed=observed,
                strict_targets=bool(item["strict_targets"]),
                obach_table_row=item.get("obach_table_row"),
                notes=str(item.get("notes", "")),
            )
        )

    return entries


def _without_kp_overrides(props: CompoundProperties) -> CompoundProperties:
    """Return a copy of props with empirical_kp_by_tissue cleared.

    Uses nested model_copy so other fields of DistributionProperties
    that may be added in future (e.g. Vss_pred, tissue-level fu) are
    preserved rather than silently reset to their defaults.
    """
    new_distribution = props.distribution.model_copy(
        update={"empirical_kp_by_tissue": None}
    )
    return props.model_copy(update={"distribution": new_distribution})


# ---------------------------------------------------------------------------
# Single-compound run
# ---------------------------------------------------------------------------


def _run_one(
    compound: CompoundConfig,
    entry: PanelEntry,
    mode: Literal["no_override", "with_override"],
) -> PanelRow:
    pipe = Pipeline(
        compound=compound,
        route=entry.route,  # type: ignore[arg-type]
        dose_mg=entry.dose_mg,
        duration_h=entry.duration_h,
    )
    result = pipe.run()
    pk = result.pk_parameters

    predicted = {
        "cl_L_h": float(pk.cl_apparent) if pk.cl_apparent is not None else float("nan"),
        "vss_L": float(pk.vss) if pk.vss is not None else float("nan"),
        "t_half_h": float(pk.half_life) if pk.half_life is not None else float("nan"),
    }
    fold: dict[str, float] = {}
    pass_2_fold: dict[str, bool] = {}
    for m in _METRICS:
        p = predicted[m]
        o = entry.observed[m]
        if p <= 0 or math.isnan(p):
            fold[m] = float("inf")
            pass_2_fold[m] = False
        else:
            fold[m] = fold_error(p, o)
            pass_2_fold[m] = fold[m] <= 2.0

    override_tissues = [
        ovr["tissue"] for ovr in result.metadata.get("kp_overrides", [])
    ]

    return PanelRow(
        key=entry.key,
        predicted=predicted,
        observed=entry.observed,
        fold=fold,
        pass_2_fold=pass_2_fold,
        override_tissues=override_tissues,
        strict_targets=entry.strict_targets,
        mode=mode,
    )


# ---------------------------------------------------------------------------
# Aggregation
# ---------------------------------------------------------------------------


def aggregate_summary(rows: list[PanelRow], mode: str) -> PanelSummary:
    n = len(rows)
    if n == 0:
        return PanelSummary(
            n=0,
            mode=mode,
            aafe={m: float("nan") for m in _METRICS},
            within_2_fold={m: float("nan") for m in _METRICS},
            within_3_fold={m: float("nan") for m in _METRICS},
            strict_failures=0,
        )

    aafe_by_metric: dict[str, float] = {}
    w2_by_metric: dict[str, float] = {}
    w3_by_metric: dict[str, float] = {}
    for metric in _METRICS:
        preds = [r.predicted[metric] for r in rows]
        obs = [r.observed[metric] for r in rows]
        aafe_by_metric[metric] = aafe(preds, obs)
        w2_by_metric[metric] = within_n_fold(preds, obs, n=2.0)
        w3_by_metric[metric] = within_n_fold(preds, obs, n=3.0)

    strict_failures = 0
    for r in rows:
        if not r.strict_targets:
            continue
        if not all(r.pass_2_fold[m] for m in _METRICS):
            strict_failures += 1

    return PanelSummary(
        n=n,
        mode=mode,
        aafe=aafe_by_metric,
        within_2_fold=w2_by_metric,
        within_3_fold=w3_by_metric,
        strict_failures=strict_failures,
    )


# ---------------------------------------------------------------------------
# Main execution
# ---------------------------------------------------------------------------


def run_benchmark(panel_path: Path) -> dict[str, PanelSummary]:
    """Execute the two-pass benchmark and return summaries (no I/O to disk)."""
    panel = load_panel(panel_path)

    rows_no_override: list[PanelRow] = []
    rows_with_override: list[PanelRow] = []

    for entry in panel:
        stripped_props = _without_kp_overrides(entry.compound.properties)
        stripped = entry.compound.model_copy(update={"properties": stripped_props})

        rows_no_override.append(_run_one(stripped, entry, mode="no_override"))
        rows_with_override.append(_run_one(entry.compound, entry, mode="with_override"))

    return {
        "no_override": aggregate_summary(rows_no_override, mode="no_override"),
        "with_override": aggregate_summary(rows_with_override, mode="with_override"),
    }, {
        "no_override": rows_no_override,
        "with_override": rows_with_override,
    }  # type: ignore[return-value]


def _print_panel_table(rows: list[PanelRow], title: str) -> None:
    print()
    print(title)
    print("-" * 100)
    header = (
        f"{'compound':<14}"
        f"{'CL_pred':>10}{'CL_obs':>10}{'f_CL':>8}  |"
        f"{'Vss_pred':>10}{'Vss_obs':>10}{'f_Vss':>8}  |"
        f"{'t½_pred':>10}{'t½_obs':>10}{'f_t½':>8}  "
        f"{'verdict':>8}"
    )
    print(header)
    print("-" * 100)
    for r in rows:
        verdict = "PASS" if all(r.pass_2_fold.values()) else "FAIL"
        if r.strict_targets:
            verdict = verdict + "*"
        print(
            f"{r.key:<14}"
            f"{r.predicted['cl_L_h']:>10.3f}{r.observed['cl_L_h']:>10.3f}"
            f"{r.fold['cl_L_h']:>8.2f}  |"
            f"{r.predicted['vss_L']:>10.2f}{r.observed['vss_L']:>10.2f}"
            f"{r.fold['vss_L']:>8.2f}  |"
            f"{r.predicted['t_half_h']:>10.2f}{r.observed['t_half_h']:>10.2f}"
            f"{r.fold['t_half_h']:>8.2f}  "
            f"{verdict:>8}"
        )


def _print_summary(summary: PanelSummary, title: str) -> None:
    print()
    print(f"{title}  (n={summary.n})")
    print(
        f"  AAFE         CL: {summary.aafe['cl_L_h']:.2f}   "
        f"Vss: {summary.aafe['vss_L']:.2f}   "
        f"t½: {summary.aafe['t_half_h']:.2f}"
    )
    print(
        f"  within 2x    CL: {summary.within_2_fold['cl_L_h']*100:.0f}%   "
        f"Vss: {summary.within_2_fold['vss_L']*100:.0f}%   "
        f"t½: {summary.within_2_fold['t_half_h']*100:.0f}%"
    )
    print(
        f"  within 3x    CL: {summary.within_3_fold['cl_L_h']*100:.0f}%   "
        f"Vss: {summary.within_3_fold['vss_L']*100:.0f}%   "
        f"t½: {summary.within_3_fold['t_half_h']*100:.0f}%"
    )
    print(f"  strict gate failures: {summary.strict_failures}")


def main(panel_path: Path | None = None) -> int:
    panel_path = panel_path or DEFAULT_PANEL_PATH
    summaries_and_rows = run_benchmark(panel_path)
    summaries, rows = summaries_and_rows[0], summaries_and_rows[1]  # type: ignore[index]

    print("=" * 100)
    print("Charon Layer 2 Human PBPK Benchmark — Obach 1999 Tier-1 Panel")
    print("=" * 100)

    _print_panel_table(rows["no_override"], "R&R only (no empirical overrides)")
    _print_summary(summaries["no_override"], "Panel summary — R&R only")

    _print_panel_table(rows["with_override"], "R&R + empirical overrides")
    _print_summary(summaries["with_override"], "Panel summary — with overrides")

    # Report applied overrides
    applied = [r for r in rows["with_override"] if r.override_tissues]
    if applied:
        print()
        print("Overrides applied:")
        for r in applied:
            for tissue in r.override_tissues:
                print(f"  {r.key}: {tissue}")

    print("=" * 100)
    strict_failures = summaries["with_override"].strict_failures
    if strict_failures == 0:
        print("All strict-targets compounds PASS — exit 0")
    else:
        print(f"{strict_failures} strict-targets compound(s) FAILED — exit 1")
    print("=" * 100)

    return 0 if strict_failures == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
```

Note: `run_benchmark()` returns a tuple `(summaries, rows)`. The test at Step 14.2 inspects `summaries` via subscription; adjust the test if needed. Actually, fix the return signature to be cleaner:

Replace the `run_benchmark` function with:

```python
def run_benchmark(panel_path: Path) -> tuple[dict[str, PanelSummary], dict[str, list[PanelRow]]]:
    """Execute the two-pass benchmark and return (summaries, rows)."""
    panel = load_panel(panel_path)

    rows_no_override: list[PanelRow] = []
    rows_with_override: list[PanelRow] = []

    for entry in panel:
        stripped_props = _without_kp_overrides(entry.compound.properties)
        stripped = entry.compound.model_copy(update={"properties": stripped_props})

        rows_no_override.append(_run_one(stripped, entry, mode="no_override"))
        rows_with_override.append(_run_one(entry.compound, entry, mode="with_override"))

    summaries = {
        "no_override": aggregate_summary(rows_no_override, mode="no_override"),
        "with_override": aggregate_summary(rows_with_override, mode="with_override"),
    }
    rows = {
        "no_override": rows_no_override,
        "with_override": rows_with_override,
    }
    return summaries, rows
```

Then update `main()` to unpack cleanly:

```python
    summaries, rows = run_benchmark(panel_path)
```

Also update `test_panel_main.py::test_run_benchmark_returns_two_summaries` to:

```python
    def test_run_benchmark_returns_two_summaries(self):
        panel_path = (
            REPO_ROOT / "validation" / "data" / "tier1_obach" / "panel.yaml"
        )
        summaries, rows = run_benchmark(panel_path)
        assert "no_override" in summaries
        assert "with_override" in summaries
        assert summaries["no_override"].n == 10
        assert summaries["with_override"].n == 10
        assert len(rows["no_override"]) == 10
        assert len(rows["with_override"]) == 10
```

- [ ] **Step 14.5: Run benchmark tests**

```bash
pytest tests/unit/test_panel_main.py -v
```
Expected: 2 tests PASS. If `test_default_panel_runs` fails because a non-strict compound returns non-finite PK or the strict theophylline gate trips, diagnose (check the per-compound output in `capsys.readouterr()`).

Also run the benchmark as a script:

```bash
python3 validation/benchmarks/layer2_human_pk.py 2>&1 | tail -50
```
Expected: Full dual-pass table, panel summaries printed, exit 0 (theophylline PASS). Record the measured AAFE values from the output — they become the session summary deliverable.

- [ ] **Step 14.6: Run full test suite for regression**

```bash
pytest tests/ -q 2>&1 | tail -10
```
Expected: baseline +58 (56 from T1-12 + 2 from T14).

- [ ] **Step 14.7: Commit**

```bash
git add validation/data/tier1_obach/panel.yaml validation/benchmarks/layer2_human_pk.py tests/unit/test_panel_main.py
git commit -m "$(cat <<'EOF'
Sprint 3b: Rewrite layer2_human_pk benchmark for Obach panel + dual-pass

Replaces the Sprint 3a Python factory approach (theophylline() and
midazolam() hand-constructed in-file) with panel.yaml loading + a
two-pass execution (no-override baseline + with-override pass). Main
prints per-compound fold errors, panel-level AAFE / within-2-fold /
within-3-fold summaries for both passes, and lists which tissues were
overridden. Exit 0 iff every strict_targets compound passes 2-fold
across all three metrics.

The benchmark output is the session's primary measurement deliverable:
panel AAFE for CL / Vss / t½ is reported on stdout in both modes.

Observed PK values in panel.yaml are canonical Obach 1999 Table 6
references; verify against the actual paper if discrepancies appear.

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 15: Record the no-override baseline AAFE (measurement)

This task does not write code — it runs the benchmark and records the measured panel AAFE as the **pre-override baseline**. This baseline is what Sprint 3b Session 1 exists to measure.

**Files:**
- Create: `validation/data/tier1_obach/baseline_rr_only.txt` (or similar)

- [ ] **Step 15.1: Run the benchmark and capture output**

```bash
cd /home/jam/Charon
python3 validation/benchmarks/layer2_human_pk.py 2>&1 | tee validation/data/tier1_obach/baseline_rr_only.txt
```

- [ ] **Step 15.2: Verify sanity floor**

Inspect `baseline_rr_only.txt`. Extract the "Panel summary — R&R only" AAFE_Vss value. The **sanity floor** is AAFE_Vss < 5.0 (session post-condition per DoD §7).

- If AAFE_Vss ≥ 5.0 in no-override mode: **STOP** and investigate. Likely causes (per spec §8 risks):
  - Obach curation error (typo in observed CL/Vss/t_half)
  - compound_type mis-classification (e.g. propranolol/verapamil getting wrong branch)
  - Actual engine bug surfaced by non-neutral chemistry that Sprint 3a theophylline never exercised
  - Re-read each compound YAML against Obach 1999 Table 2/6
  - Re-run `python3 -c 'from charon import Pipeline; ...'` on the worst offender and debug
- If AAFE_Vss < 5.0: proceed. The Sprint 3b session has a valid measurement baseline.

- [ ] **Step 15.3: Commit the baseline artifact**

```bash
git add validation/data/tier1_obach/baseline_rr_only.txt
git commit -m "$(cat <<'EOF'
Sprint 3b: Record no-override baseline Obach panel measurement

Captured stdout from a full benchmark run before any empirical Kp
overrides are committed to the compound YAMLs. This is the R&R-only
measurement surface — the "before" state that the with-override
comparison (Tasks 16-18) builds on. AAFE values for CL, Vss, t½ across
the 10-compound panel are preserved as a permanent session artifact.

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 16: Citation-verification task for the PRIMARY override (midazolam adipose)

This is a citation-verification task, NOT a code task. Per spec §7.2, the numerical override may not be committed until the citation is verified in an actual paper.

**Files:**
- Modify: `validation/data/tier1_obach/compounds/midazolam.yaml` (add Phase 2 override)

- [ ] **Step 16.1: Attempt to verify the primary citation**

Primary target: Björkman S. (2001/2002) "Prediction of drug disposition in infants and children by means of physiologically based pharmacokinetic (PBPK) modeling", J Pharm Sci or equivalent. The paper should report rat or human tissue Kp values for midazolam, specifically an adipose:plasma ratio.

**Steps:**

1. Search for the paper (or equivalent midazolam PBPK paper) in accessible sources:
   - Google Scholar, PubMed, the institution's library, or a published review that reproduces the table (e.g. a subsequent PBPK review citing Björkman).
2. If accessed: record (a) exact adipose Kp value, (b) the page/table/row reference, (c) whether the value is rat or human.
3. If not accessible: invoke **Fallback 1** — Rodgers et al. 2005 Table 4/5 rat experimental Kp for the closest structural analog (expect it to NOT be midazolam directly since Rodgers' panel excludes it). The closest match from Rodgers' panel is usually diazepam or another benzodiazepine — but Rodgers' panel members are explicitly listed; read the paper to identify a valid analog.
4. If Fallback 1 also fails: invoke **Fallback 2** — move the override target from midazolam to a different compound that DOES have clean literature Kp values (propranolol and diazepam are both well-characterized in Poulin & Theil 2002). Note this in the commit message.
5. If all paths fail: **STOP**. Report the blockage to the user and ask for direction. Do NOT commit a synthetic value (Fallback 3 is forbidden in Sprint 3b Session 1 per spec §7.2 and DoD §2).

- [ ] **Step 16.2: Edit `midazolam.yaml` to add the verified override**

Assuming the verification succeeds with a value `X.Y` from a paper `<citation>`, add the following block to `midazolam.yaml` at the bottom of the `properties:` section (note: YAML indentation must match the surrounding structure):

```yaml
  distribution:
    empirical_kp_by_tissue:
      adipose:
        value: <X.Y>
        source: literature
        method: "<exact citation with page/table/row>"
```

- [ ] **Step 16.3: Verify the updated YAML loads**

```bash
python3 - <<'EOF'
import yaml
from pathlib import Path
from charon.core.schema import CompoundConfig
p = Path("validation/data/tier1_obach/compounds/midazolam.yaml")
cc = CompoundConfig.model_validate(yaml.safe_load(p.read_text()))
dist = cc.properties.distribution
print("OK empirical_kp_by_tissue:", dist.empirical_kp_by_tissue)
EOF
```
Expected: the empirical_kp_by_tissue dict prints with the verified value and citation.

- [ ] **Step 16.4: Re-run the benchmark to see the override effect**

```bash
python3 validation/benchmarks/layer2_human_pk.py 2>&1 | tee /tmp/post_midazolam_override.txt
grep -A3 "with overrides" /tmp/post_midazolam_override.txt | tail -10
```
Expected: "with overrides" summary shows a different Vss_pred for midazolam; the "Overrides applied" footer lists `midazolam: adipose`.

- [ ] **Step 16.5: Commit the Phase 2 midazolam override**

```bash
git add validation/data/tier1_obach/compounds/midazolam.yaml
git commit -m "$(cat <<'EOF'
Sprint 3b: Add verified midazolam adipose Kp override (Phase 2)

empirical_kp_by_tissue.adipose value sourced from [CITATION]. First of
≥2 verified override compounds required by DoD §2.

Citation: <full citation with page/table/row>

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 17: Citation-verification task for the SECONDARY override

Target: propranolol or diazepam. Candidate primary sources: Rodgers et al. 2005 Table 4 (propranolol IS in Rodgers' panel as a base; diazepam is NOT in the 2005 base-specific paper). Poulin & Theil 2002 also has propranolol and diazepam Kp.

**Files:**
- Modify: `validation/data/tier1_obach/compounds/propranolol.yaml` OR `diazepam.yaml`

- [ ] **Step 17.1: Pick the candidate with the cleanest literature coverage**

- Default: **propranolol** (Rodgers 2005 reports rat tissue Kp for propranolol).
- Fallback: **diazepam** if propranolol's Rodgers value is ambiguous or missing.

- [ ] **Step 17.2: Verify the secondary citation**

Same protocol as Task 16: access the source, record value + citation string, fall back to an alternative source if unavailable. Target tissue is whichever Rodgers/Poulin-Theil reports reliably for the chosen compound — likely adipose or muscle.

If both propranolol AND diazepam citation paths fail: **STOP** and escalate to the user. Do not proceed with a synthetic fallback.

- [ ] **Step 17.3: Edit the chosen compound YAML to add the verified override**

Example for propranolol (actual values and citations filled from the paper):

```yaml
  distribution:
    empirical_kp_by_tissue:
      adipose:
        value: <X.Y>
        source: literature
        method: "Rodgers et al. 2005 J Pharm Sci 94(6):1259, Table 4 propranolol adipose"
```

- [ ] **Step 17.4: Verify the YAML loads**

```bash
python3 - <<'EOF'
import yaml
from pathlib import Path
from charon.core.schema import CompoundConfig
p = Path("validation/data/tier1_obach/compounds/propranolol.yaml")  # or diazepam
cc = CompoundConfig.model_validate(yaml.safe_load(p.read_text()))
print("OK:", cc.properties.distribution.empirical_kp_by_tissue)
EOF
```

- [ ] **Step 17.5: Re-run benchmark with both overrides in place**

```bash
python3 validation/benchmarks/layer2_human_pk.py 2>&1 | tee validation/data/tier1_obach/with_overrides.txt
```
Expected: "Overrides applied" section lists **both** midazolam and the secondary compound. Panel AAFE for Vss should be lower (or at least not worse) than the no-override baseline from Task 15.

- [ ] **Step 17.6: Commit**

```bash
git add validation/data/tier1_obach/compounds/<compound>.yaml validation/data/tier1_obach/with_overrides.txt
git commit -m "$(cat <<'EOF'
Sprint 3b: Add verified secondary Kp override (<compound>) + with-overrides benchmark

Second of ≥2 verified override compounds required by DoD §2. With
both midazolam and <compound> overrides in place, the benchmark's
"with_override" pass now has meaningful deltas from the R&R-only
baseline.

Citation: <full citation>

Captured the with-overrides benchmark stdout as a permanent artifact
next to the baseline file for session audit.

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 18: Integration smoke test with the full panel

**Files:**
- Create: `tests/integration/test_obach_panel_smoke.py`

- [ ] **Step 18.1: Write the integration smoke test**

Create `tests/integration/test_obach_panel_smoke.py`:

```python
"""End-to-end Obach panel smoke test.

Loads the real panel.yaml, runs the full benchmark, and asserts:
  - All 10 compounds produce finite positive PK values in both modes.
  - theophylline strict 2-fold gate passes (regression invariant).
  - Panel AAFE_Vss with-override stays under the 5.0 sanity floor.
  - ≥2 compounds have recorded kp_overrides (DoD §2 check).

The actual AAFE values are `record_property`'d into pytest results
so CI artifacts preserve the measurement.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from validation.benchmarks.layer2_human_pk import (  # noqa: E402
    DEFAULT_PANEL_PATH,
    run_benchmark,
)


@pytest.fixture(scope="module")
def benchmark_result():
    summaries, rows = run_benchmark(DEFAULT_PANEL_PATH)
    return summaries, rows


class TestObachPanelSmoke:
    def test_all_compounds_produce_finite_pk_no_override(self, benchmark_result):
        _, rows = benchmark_result
        for r in rows["no_override"]:
            for metric in ("cl_L_h", "vss_L", "t_half_h"):
                assert math.isfinite(r.predicted[metric]), (
                    f"{r.key} {metric} non-finite: {r.predicted[metric]}"
                )
                assert r.predicted[metric] > 0, (
                    f"{r.key} {metric} non-positive: {r.predicted[metric]}"
                )

    def test_all_compounds_produce_finite_pk_with_override(self, benchmark_result):
        _, rows = benchmark_result
        for r in rows["with_override"]:
            for metric in ("cl_L_h", "vss_L", "t_half_h"):
                assert math.isfinite(r.predicted[metric]), (
                    f"{r.key} {metric} non-finite: {r.predicted[metric]}"
                )
                assert r.predicted[metric] > 0, (
                    f"{r.key} {metric} non-positive: {r.predicted[metric]}"
                )

    def test_theophylline_strict_gate_passes(self, benchmark_result):
        _, rows = benchmark_result
        theo_rows = [r for r in rows["with_override"] if r.key == "theophylline"]
        assert len(theo_rows) == 1
        theo = theo_rows[0]
        assert theo.strict_targets is True
        for metric in ("cl_L_h", "vss_L", "t_half_h"):
            assert theo.pass_2_fold[metric], (
                f"theophylline {metric} fold={theo.fold[metric]:.2f} > 2.0"
            )

    def test_panel_aafe_vss_below_sanity_floor(self, benchmark_result, record_property):
        summaries, _ = benchmark_result
        aafe_vss_yes = summaries["with_override"].aafe["vss_L"]
        aafe_vss_no = summaries["no_override"].aafe["vss_L"]
        record_property("aafe_vss_with_override", aafe_vss_yes)
        record_property("aafe_vss_no_override", aafe_vss_no)
        assert aafe_vss_yes < 5.0, (
            f"Session post-condition failed: AAFE_Vss with-override = "
            f"{aafe_vss_yes:.2f} exceeds 5.0 sanity floor"
        )

    def test_panel_aafe_cl_below_sanity_floor(self, benchmark_result, record_property):
        summaries, _ = benchmark_result
        aafe_cl_yes = summaries["with_override"].aafe["cl_L_h"]
        aafe_cl_no = summaries["no_override"].aafe["cl_L_h"]
        record_property("aafe_cl_with_override", aafe_cl_yes)
        record_property("aafe_cl_no_override", aafe_cl_no)
        # CL sanity floor: 5.0 (same as Vss, generous)
        assert aafe_cl_yes < 5.0

    def test_at_least_two_compounds_have_verified_overrides(self, benchmark_result):
        """DoD §2: ≥2 compounds must carry empirical Kp overrides."""
        _, rows = benchmark_result
        compounds_with_overrides = {
            r.key for r in rows["with_override"] if r.override_tissues
        }
        assert len(compounds_with_overrides) >= 2, (
            f"DoD §2 violation: only {len(compounds_with_overrides)} "
            f"compound(s) have verified overrides: {compounds_with_overrides}"
        )

    def test_within_2_fold_reported(self, benchmark_result, record_property):
        summaries, _ = benchmark_result
        for mode in ("no_override", "with_override"):
            for metric in ("cl_L_h", "vss_L", "t_half_h"):
                frac = summaries[mode].within_2_fold[metric]
                record_property(f"within_2_fold__{mode}__{metric}", frac)
                assert 0.0 <= frac <= 1.0
```

- [ ] **Step 18.2: Run the smoke test**

```bash
pytest tests/integration/test_obach_panel_smoke.py -v
```
Expected: 7 tests PASS. If the `test_panel_aafe_vss_below_sanity_floor` fails with AAFE_Vss ≥ 5.0, the session's sanity floor tripped — investigate per Task 15 §15.2 before proceeding.

- [ ] **Step 18.3: Run full test suite for final regression**

```bash
pytest tests/ -q 2>&1 | tail -20
```
Expected: baseline +65 (58 previous new + 7 integration smoke). Record the exact count for the session summary.

- [ ] **Step 18.4: Commit**

```bash
git add tests/integration/test_obach_panel_smoke.py
git commit -m "$(cat <<'EOF'
Sprint 3b: Add Obach panel integration smoke test

End-to-end test that loads the real panel.yaml, runs the two-pass
benchmark, and asserts: (1) all 10 compounds produce finite positive
PK in both modes, (2) theophylline strict 2-fold gate holds, (3)
panel AAFE_Vss with-override stays under 5.0 sanity floor, (4)
panel AAFE_CL with-override stays under 5.0, (5) ≥2 compounds have
verified overrides (DoD §2 check), (6) within-2-fold fractions are
reported as pytest record_property metadata so CI preserves the
measured values.

Numerical AAFE values are tracked via record_property, NOT asserted
against ARCHITECTURE targets — the session is measurement mode, not
target-hitting mode.

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 19: Final regression sweep + session summary

**Files:**
- Create: `docs/superpowers/sessions/2026-04-10-sprint3b-session1-summary.md` (or similar)

- [ ] **Step 19.1: Run the complete test suite one final time**

```bash
cd /home/jam/Charon
pytest tests/ -v 2>&1 | tail -40
```
Expected: all tests pass. Record the exact count (should be Sprint 3a baseline + ~65 new).

- [ ] **Step 19.2: Run the benchmark as a standalone script one final time**

```bash
python3 validation/benchmarks/layer2_human_pk.py 2>&1 | tee /tmp/final_benchmark.txt
```
Expected: exit 0, both panel summaries printed, overrides applied list shows ≥2 compounds. Extract:

- Panel AAFE (no_override, with_override) for CL, Vss, t½
- within-2-fold fractions for both modes
- Which compounds received overrides and which tissues

- [ ] **Step 19.3: Author session summary document**

Create `docs/superpowers/sessions/2026-04-10-sprint3b-session1-summary.md`. If the `sessions/` directory doesn't exist, create it:

```bash
mkdir -p /home/jam/Charon/docs/superpowers/sessions
```

Content template (fill in with actual measured values):

```markdown
# Sprint 3b Session 1 Summary — Honest IV Kernel

**Date:** 2026-04-10
**Spec:** `docs/superpowers/specs/2026-04-10-sprint3b-kp-override-validation-design.md`
**Plan:** `docs/superpowers/plans/2026-04-10-sprint3b-kp-override-validation.md`
**Sprint 3a baseline:** 477 tests, 1 gated compound (theophylline)
**Sprint 3b final:** <count> tests, 1 gated compound (theophylline), 10-compound panel AAFE measured

## Deliverables

- Schema: `DistributionProperties.empirical_kp_by_tissue`,
  `PhysicochemicalProperties.compound_type`, `CompoundProperties.distribution`
- ode_compiler: `KpOverrideRecord` dataclass, `CompoundPBPKParams.kp_overrides`
  field, species guard, compound_type precedence, empirical override loop
- Pipeline: `PipelineResult.metadata['kp_overrides']` exposed
- Validation data: `validation/data/tier1_obach/panel.yaml` + 10 compound YAMLs
- Benchmark: rewritten `layer2_human_pk.py` with dual-pass execution
- Tests: ~65 new tests across schema, override, compound_type, panel loader,
  panel metrics, midazolam mechanism, Obach panel smoke

## Measured panel AAFE (10 compounds, Obach 1999)

| Metric | Mode: R&R only | Mode: with override |
|---|---|---|
| AAFE CL | <X.XX> | <Y.YY> |
| AAFE Vss | <X.XX> | <Y.YY> |
| AAFE t½ | <X.XX> | <Y.YY> |
| within-2-fold CL | <N>/10 | <N>/10 |
| within-2-fold Vss | <N>/10 | <N>/10 |
| within-2-fold t½ | <N>/10 | <N>/10 |

## Override compounds (≥2 verified — DoD §2 satisfied)

| Compound | Tissue | R&R value | Empirical value | Citation |
|---|---|---|---|---|
| midazolam | adipose | <X.X> | <Y.Y> | <full citation> |
| <secondary> | <tissue> | <X.X> | <Y.Y> | <full citation> |

## Sanity floor result

AAFE_Vss (with override) = <Y.YY> — **PASS** (< 5.0 sanity floor per DoD §7)

## Residual scientific gap

<If the measured AAFE_Vss > 3.0 ARCHITECTURE target: document the gap honestly.
If the next session will pursue Kp model work (Berezhkovskiy, Poulin-Theil) vs.
accept and move to ACAT: note the decision point here.>

## Known Future Work (carry forward to later sessions)

1. **Species-aware MPPGL / hepatocellularity.**
   `build_compound_pbpk_params` still hardcodes human values; guard
   raises NotImplementedError on non-human topology. Sprint 4 must add
   `mppgl_mg_g` and `hepatocellularity_1e6_per_g` to `PBPKTopology`
   (read from species YAML) and remove the guard.
2. **Zwitterion coverage.** Obach panel has no zwitterion; cetirizine
   candidate deferred.
3. **compound_type threshold revisit.** Default inference uses
   `pKa_base > 8.0`; physiologically sound threshold is 7.4. YAML-level
   override covers validation needs; defaults untouched.
4. **Automated Obach curation.** Hand-curation of 10 YAMLs was
   tractable; scaling to N≥40 needs a CSV/JSON loader.
5. **Session 3b-2: ACAT/oral pipeline** — the next logical session, now
   building on a measurement-verified IV kernel.

## Commits (this session)

<`git log --oneline 0051429..HEAD` output>
```

Fill in all `<...>` placeholders with the actual measured values from Task 19 Step 2 and the commit history.

- [ ] **Step 19.4: Commit session summary**

```bash
git add docs/superpowers/sessions/2026-04-10-sprint3b-session1-summary.md
git commit -m "$(cat <<'EOF'
Sprint 3b Session 1: Final session summary

Measured Obach 1999 Tier-1 panel AAFE for CL/Vss/t½ across 10
compounds in both R&R-only and with-empirical-overrides modes.
≥2 compounds verified (midazolam, <secondary>). Sanity floor
AAFE_Vss with-override < 5.0 passed.

Session achieves the 'honest Layer 2 PBPK engine' goal: the IV
kernel's performance is now a measured number rather than a claim.
Next session decides Kp model work vs. ACAT based on the measured
residual gap.

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Self-Review Checklist (execute at end)

After all 19 tasks complete, run this checklist:

- [ ] `pytest tests/ -q` → all tests pass, count recorded
- [ ] `python3 validation/benchmarks/layer2_human_pk.py` → exit 0, dual-pass output
- [ ] `git log --oneline 0051429..HEAD` → every task has exactly one commit
- [ ] `grep -r "TODO\|FIXME\|XXX" src/charon/core/schema.py src/charon/pbpk/ode_compiler.py src/charon/pipeline.py validation/benchmarks/layer2_human_pk.py` → no placeholders in Sprint 3b code
- [ ] Spec §13 approval checklist items all satisfied
- [ ] DoD §2 items 1-7 all satisfied
- [ ] ≥2 verified override compounds (not synthetic)
- [ ] Panel AAFE measured values recorded in session summary
- [ ] Known future work (species-aware mppgl) documented in summary

If any item fails: diagnose root cause, fix, re-run the full regression sweep, then re-check.
