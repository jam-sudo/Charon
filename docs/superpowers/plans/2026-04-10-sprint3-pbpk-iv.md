# Sprint 3 — Minimal IV-only Human PBPK — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Ship a working human PBPK simulation engine for IV bolus/infusion so that `Pipeline(compound, route="iv_bolus", dose_mg=5).run()` produces a plasma Cp-time profile and PK parameters, validated against midazolam IV within 2-fold.

**Architecture:** 15-tissue perfusion-limited PBPK with portal/hepatic split; stiff BDF ODE solver; pure-Python NumPy rhs built from species YAML topology; top-level `Pipeline` class wires Layer 0→1→ParameterBridge→Layer 2.

**Tech Stack:** Python 3.11, pydantic v2, scipy.integrate.solve_ivp (BDF), numpy, pyyaml, pytest.

**Spec:** `docs/superpowers/specs/2026-04-10-sprint3-pbpk-iv-design.md`

**Reference drugs (from test_parameter_bridge.py fixtures):**
- Midazolam: CLint=93 uL/min/mg, fu_inc=0.96, fu_p=0.03, BP=0.66, MW=325.77, logP=3.89, pKa_base=6.2

---

## File Map

**Create:**
- `src/charon/pbpk/topology.py` — PBPKTopology dataclass + YAML loader
- `src/charon/pbpk/ode_compiler.py` — CompoundPBPKParams + build_rhs
- `src/charon/pbpk/solver.py` — SimulationResult + simulate_iv
- `src/charon/pbpk/pk_extract.py` — compute_pk_parameters (Cmax/AUC/CL/Vss)
- `src/charon/pipeline.py` — Pipeline + PipelineResult
- `tests/unit/test_topology.py`
- `tests/unit/test_ode_compiler.py`
- `tests/unit/test_solver.py`
- `tests/unit/test_pk_extract.py`
- `tests/unit/test_pipeline.py`
- `tests/integration/test_smiles_to_pk.py`
- `validation/benchmarks/metrics.py`
- `validation/benchmarks/layer2_human_pk.py`
- `tests/unit/test_benchmark_metrics.py`

**Modify:**
- `src/charon/core/schema.py` — add `clint_liver_L_h` field to HepaticClearance
- `src/charon/core/parameter_bridge.py` — populate the new field
- `src/charon/pbpk/__init__.py` — exports
- `src/charon/__init__.py` — export Pipeline and PipelineResult

**Unchanged:** kp_calculator.py, species/*.yaml, core/liver_models.py, predict/*, all other Sprint 1/2 code.

---

## Task 1: Extend HepaticClearance schema with clint_liver_L_h

Add `clint_liver_L_h` to `HepaticClearance` so PBPK ODEs can consume the whole-liver intrinsic clearance without re-doing IVIVE scaling. Backward compatible: optional field, existing tests ignore it.

**Files:**
- Modify: `src/charon/core/schema.py` (HepaticClearance class)
- Modify: `src/charon/core/parameter_bridge.py` (populate new field)
- Test: `tests/unit/test_parameter_bridge.py` (new test for the field)

- [ ] **Step 1: Write the failing test**

Append to `tests/unit/test_parameter_bridge.py`:

```python
class TestHepaticClearanceClintLiverField:
    """The clint_liver_L_h field is needed by PBPK ODEs to avoid re-scaling."""

    def test_clint_liver_L_h_populated_hlm(self, bridge):
        """For HLM with the ARCHITECTURE example, CLint_liver = 64.377 L/h."""
        result = bridge.clint_to_clh(
            clint=15.2, fu_inc=0.85, fu_p=0.23, bp_ratio=0.95, qh_L_h=90.0
        )
        # Hand calc: 15.2/0.85 * 40 * 1500 / 1e6 * 60 = 64.377 L/h
        assert result.clint_liver_L_h == pytest.approx(64.377, rel=1e-3)

    def test_clint_liver_L_h_populated_hepatocytes(self, bridge):
        """For hepatocytes: CLint_u = CLint (no fu_inc correction)."""
        result = bridge.clint_to_clh(
            clint=10.0, fu_inc=0.85, fu_p=0.5, bp_ratio=1.0,
            system="hepatocytes"
        )
        # Hand calc: 10.0 * 120 * 1500 / 1e6 * 60 = 108.0 L/h
        assert result.clint_liver_L_h == pytest.approx(108.0, rel=1e-6)
```

- [ ] **Step 2: Run test to verify it fails**

```
python3 -m pytest tests/unit/test_parameter_bridge.py::TestHepaticClearanceClintLiverField -v
```

Expected: FAIL with `AttributeError` or `pydantic.ValidationError` (field missing).

- [ ] **Step 3: Add field to HepaticClearance schema**

Edit `src/charon/core/schema.py`, replace the HepaticClearance class with:

```python
class HepaticClearance(BaseModel):
    """Result of an in-vitro → in-vivo hepatic-clearance extrapolation."""

    clh_L_h: float
    extraction_ratio: float
    model_used: str
    conversion_log: ConversionLog
    clint_liver_L_h: float | None = None
    """Whole-liver intrinsic clearance after IVIVE scaling, in L/h.

    This is the value *before* the liver extraction model is applied.  It
    is the correct input for PBPK ODEs that embed well-stirred elimination
    directly (rate_elim = CLint_liver * fu_b * C_liver_blood_out).  Using
    ``clh_L_h`` instead in a PBPK rhs would double-apply the liver model.
    """
```

- [ ] **Step 4: Populate the field in parameter_bridge**

Edit `src/charon/core/parameter_bridge.py`, in `clint_to_clh`, replace the `return HepaticClearance(...)` block (lines 246-251) with:

```python
        return HepaticClearance(
            clh_L_h=clh,
            extraction_ratio=extraction_ratio,
            model_used=model,
            conversion_log=conversion_log,
            clint_liver_L_h=clint_liver_L_h,
        )
```

- [ ] **Step 5: Run test to verify it passes**

```
python3 -m pytest tests/unit/test_parameter_bridge.py::TestHepaticClearanceClintLiverField -v
```

Expected: 2 passed.

- [ ] **Step 6: Run full parameter_bridge tests**

```
python3 -m pytest tests/unit/test_parameter_bridge.py -v
```

Expected: all existing tests still pass + 2 new.

- [ ] **Step 7: Commit**

```
git add src/charon/core/schema.py src/charon/core/parameter_bridge.py tests/unit/test_parameter_bridge.py
git commit -m "Sprint 3: Add clint_liver_L_h field to HepaticClearance for PBPK ODE consumers"
```

---

## Task 2: pbpk/topology.py — PBPKTopology + YAML loader

Produces a validated, immutable PBPK topology graph from a species YAML file. Includes tissue volumes, blood flows, portal/arterial identification, and Rodgers & Rowland tissue compositions.

**Files:**
- Create: `src/charon/pbpk/topology.py`
- Create: `tests/unit/test_topology.py`

- [ ] **Step 1: Write failing tests**

Create `tests/unit/test_topology.py`:

```python
"""Unit tests for PBPK topology loading."""

import math
from collections import OrderedDict
from pathlib import Path

import pytest

from charon.pbpk.topology import (
    PBPKTopology,
    TissueNode,
    load_species_topology,
    PORTAL_TISSUES,
)


@pytest.fixture
def human_topology() -> PBPKTopology:
    return load_species_topology("human")


class TestLoadHumanTopology:
    def test_returns_pbpk_topology(self, human_topology):
        assert isinstance(human_topology, PBPKTopology)
        assert human_topology.species == "human"

    def test_body_weight(self, human_topology):
        assert human_topology.body_weight_kg == pytest.approx(70.0)

    def test_cardiac_output(self, human_topology):
        assert human_topology.cardiac_output_L_h == pytest.approx(390.0)

    def test_hematocrit(self, human_topology):
        assert human_topology.hematocrit == pytest.approx(0.45)

    def test_blood_pool_volumes(self, human_topology):
        assert human_topology.venous_volume_L == pytest.approx(3.7)
        assert human_topology.arterial_volume_L == pytest.approx(1.5)

    def test_tissues_are_ordered(self, human_topology):
        assert isinstance(human_topology.tissues, OrderedDict)
        # Deterministic ordering: YAML insertion order preserved.
        names = list(human_topology.tissues.keys())
        assert names[0] == "lung"  # lung is first in human.yaml
        assert "liver" in names
        assert "kidney" in names

    def test_has_15_tissues(self, human_topology):
        assert len(human_topology.tissues) == 15

    def test_lung_has_full_cardiac_output(self, human_topology):
        lung = human_topology.tissues["lung"]
        # Lung receives full CO (pulmonary circulation)
        assert lung.blood_flow_L_h == pytest.approx(
            human_topology.cardiac_output_L_h
        )

    def test_liver_flow_equals_arterial_plus_portal(self, human_topology):
        liver_q = human_topology.tissues["liver"].blood_flow_L_h
        portal_q = sum(
            human_topology.tissues[p].blood_flow_L_h for p in PORTAL_TISSUES
        )
        ha_q = human_topology.hepatic_artery_L_h
        assert liver_q == pytest.approx(ha_q + portal_q, rel=1e-6)

    def test_hepatic_artery_is_nonnegative(self, human_topology):
        # Q_HA = Q_liver_total - Q_portal must be >= 0 for consistency.
        assert human_topology.hepatic_artery_L_h >= 0

    def test_hepatic_artery_matches_yaml_numbers(self, human_topology):
        # human.yaml comment: liver total = 0.065 arterial + 0.19 portal CO
        expected_ha = 0.065 * human_topology.cardiac_output_L_h
        assert human_topology.hepatic_artery_L_h == pytest.approx(expected_ha, rel=1e-6)

    def test_systemic_flows_sum_to_cardiac_output(self, human_topology):
        """All arterial outflow paths must sum to cardiac output.

        Outflow categories:
          - hepatic artery to liver (Q_HA)
          - portal tissues (spleen, gut_wall, pancreas) — drain into liver
          - non-portal non-liver non-lung systemic tissues — drain into venous
        """
        portal_sum = sum(
            human_topology.tissues[p].blood_flow_L_h for p in PORTAL_TISSUES
        )
        non_portal_systemic = sum(
            node.blood_flow_L_h
            for name, node in human_topology.tissues.items()
            if name not in {"lung", "liver"} and name not in PORTAL_TISSUES
        )
        total = human_topology.hepatic_artery_L_h + portal_sum + non_portal_systemic
        assert total == pytest.approx(
            human_topology.cardiac_output_L_h, rel=1e-3
        )

    def test_tissue_nodes_have_composition(self, human_topology):
        for node in human_topology.tissues.values():
            comp = node.composition
            assert math.isfinite(comp.fn)
            assert math.isfinite(comp.fp)
            assert math.isfinite(comp.fw)
            assert math.isfinite(comp.pH)

    def test_portal_tissues_drain_to_liver(self, human_topology):
        for name in PORTAL_TISSUES:
            node = human_topology.tissues[name]
            assert node.drains_to == "liver"

    def test_non_portal_tissues_drain_to_venous(self, human_topology):
        for name, node in human_topology.tissues.items():
            if name in PORTAL_TISSUES or name == "lung" or name == "liver":
                continue
            assert node.drains_to == "venous"

    def test_lung_drains_to_arterial(self, human_topology):
        assert human_topology.tissues["lung"].drains_to == "arterial"

    def test_liver_drains_to_venous(self, human_topology):
        assert human_topology.tissues["liver"].drains_to == "venous"

    def test_plasma_composition_present(self, human_topology):
        assert human_topology.plasma_composition.pH == pytest.approx(7.40)


class TestLoadMissingSpecies:
    def test_raises_file_not_found(self):
        with pytest.raises(FileNotFoundError):
            load_species_topology("unicorn")


class TestLoadFromPath:
    def test_load_from_explicit_path(self, tmp_path):
        # Build a minimal YAML and ensure the loader can ingest it.
        yaml_content = """\
species:
  name: mini
  body_weight_kg: 10.0
  cardiac_output_L_h: 50.0
  hematocrit: 0.4
  liver:
    weight_g: 200.0
    volume_L: 0.2
    blood_flow_fraction_CO: 0.2
  kidney:
    weight_g: 40.0
    volume_L: 0.04
    blood_flow_fraction_CO: 0.1
    gfr_mL_min: 10.0
  bsa_scaling:
    km: 10.0
  blood:
    venous_volume_L: 0.5
    arterial_volume_L: 0.2
    portal_vein_volume_L: 0.01
  plasma:
    fn: 0.0023
    fp: 0.0199
    fw: 0.945
    pH: 7.40
  tissues:
    lung:
      volume_L: 0.05
      blood_flow_fraction: 1.0
      composition: { fn: 0.003, fp: 0.013, fw: 0.81, pH: 7.0 }
    liver:
      volume_L: 0.2
      blood_flow_fraction: 0.2
      composition: { fn: 0.035, fp: 0.025, fw: 0.75, pH: 7.0 }
    kidney:
      volume_L: 0.04
      blood_flow_fraction: 0.1
      composition: { fn: 0.012, fp: 0.024, fw: 0.78, pH: 7.0 }
    spleen:
      volume_L: 0.02
      blood_flow_fraction: 0.05
      composition: { fn: 0.008, fp: 0.011, fw: 0.79, pH: 7.0 }
    gut_wall:
      volume_L: 0.05
      blood_flow_fraction: 0.1
      composition: { fn: 0.016, fp: 0.018, fw: 0.72, pH: 7.0 }
    pancreas:
      volume_L: 0.01
      blood_flow_fraction: 0.02
      composition: { fn: 0.035, fp: 0.025, fw: 0.75, pH: 7.0 }
    muscle:
      volume_L: 5.0
      blood_flow_fraction: 0.5
      composition: { fn: 0.024, fp: 0.007, fw: 0.76, pH: 7.0 }
"""
        yaml_path = tmp_path / "mini.yaml"
        yaml_path.write_text(yaml_content)
        topo = load_species_topology("mini", path=yaml_path)
        assert topo.species == "mini"
        # Liver = 0.2 * 50 = 10 L/h
        assert topo.tissues["liver"].blood_flow_L_h == pytest.approx(10.0)
        # Portal = (0.05 + 0.1 + 0.02) * 50 = 8.5 L/h
        # HA = 10 - 8.5 = 1.5 L/h
        assert topo.hepatic_artery_L_h == pytest.approx(1.5, rel=1e-6)


class TestTopologyConsistencyChecks:
    def test_rejects_liver_flow_less_than_portal(self, tmp_path):
        """Portal tissues total more than liver flow → inconsistent YAML."""
        yaml_content = """\
species:
  name: broken
  body_weight_kg: 10.0
  cardiac_output_L_h: 50.0
  hematocrit: 0.4
  liver:
    weight_g: 200.0
    volume_L: 0.2
    blood_flow_fraction_CO: 0.05
  kidney:
    weight_g: 40.0
    volume_L: 0.04
    blood_flow_fraction_CO: 0.1
    gfr_mL_min: 10.0
  bsa_scaling:
    km: 10.0
  blood:
    venous_volume_L: 0.5
    arterial_volume_L: 0.2
    portal_vein_volume_L: 0.01
  plasma:
    fn: 0.0023
    fp: 0.0199
    fw: 0.945
    pH: 7.40
  tissues:
    lung:
      volume_L: 0.05
      blood_flow_fraction: 1.0
      composition: { fn: 0.003, fp: 0.013, fw: 0.81, pH: 7.0 }
    liver:
      volume_L: 0.2
      blood_flow_fraction: 0.05
      composition: { fn: 0.035, fp: 0.025, fw: 0.75, pH: 7.0 }
    kidney:
      volume_L: 0.04
      blood_flow_fraction: 0.1
      composition: { fn: 0.012, fp: 0.024, fw: 0.78, pH: 7.0 }
    spleen:
      volume_L: 0.02
      blood_flow_fraction: 0.2
      composition: { fn: 0.008, fp: 0.011, fw: 0.79, pH: 7.0 }
    gut_wall:
      volume_L: 0.05
      blood_flow_fraction: 0.3
      composition: { fn: 0.016, fp: 0.018, fw: 0.72, pH: 7.0 }
    pancreas:
      volume_L: 0.01
      blood_flow_fraction: 0.05
      composition: { fn: 0.035, fp: 0.025, fw: 0.75, pH: 7.0 }
"""
        yaml_path = tmp_path / "broken.yaml"
        yaml_path.write_text(yaml_content)
        with pytest.raises(ValueError, match="portal"):
            load_species_topology("broken", path=yaml_path)
```

- [ ] **Step 2: Run to verify failures**

```
python3 -m pytest tests/unit/test_topology.py -v
```

Expected: ImportError / ModuleNotFoundError (topology module doesn't exist yet).

- [ ] **Step 3: Implement topology.py**

Create `src/charon/pbpk/topology.py`:

```python
"""PBPK topology: load a species YAML into an immutable, validated graph.

The PBPK topology describes which tissues exist, their volumes, blood flows,
tissue compositions (for Rodgers & Rowland Kp), and how they are wired into
the vascular tree (arterial input, venous or portal drainage, hepatic
arterial split).

The topology is a pure data container — no ODE logic lives here.  The
``ode_compiler`` module consumes it to build the right-hand side.
"""

from __future__ import annotations

from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import yaml

from charon.pbpk.kp_calculator import TissueComposition

# Tissues that drain into the portal vein rather than the systemic venous
# pool.  These provide the "portal inflow" term for the liver compartment.
PORTAL_TISSUES: tuple[str, ...] = ("spleen", "gut_wall", "pancreas")

_DRAINAGE_OVERRIDES: dict[str, str] = {
    "lung": "arterial",
    "liver": "venous",
}

_SPECIES_DIR = Path(__file__).parent / "species"


@dataclass(frozen=True)
class TissueNode:
    """Single tissue compartment in the PBPK graph."""

    name: str
    volume_L: float
    blood_flow_L_h: float
    composition: TissueComposition
    drains_to: Literal["venous", "liver", "arterial"]


@dataclass(frozen=True)
class PBPKTopology:
    """Immutable species-level PBPK graph."""

    species: str
    body_weight_kg: float
    cardiac_output_L_h: float
    hematocrit: float
    venous_volume_L: float
    arterial_volume_L: float
    hepatic_artery_L_h: float
    tissues: "OrderedDict[str, TissueNode]"
    plasma_composition: TissueComposition
    gfr_mL_min: float
    liver_weight_g: float

    def tissue_names(self) -> tuple[str, ...]:
        """Deterministic tissue ordering (YAML insertion order)."""
        return tuple(self.tissues.keys())


def _tissue_comp_from_dict(comp: dict) -> TissueComposition:
    return TissueComposition(
        fn=float(comp["fn"]),
        fp=float(comp["fp"]),
        fw=float(comp["fw"]),
        pH=float(comp["pH"]),
    )


def _drainage_for(name: str) -> Literal["venous", "liver", "arterial"]:
    if name in _DRAINAGE_OVERRIDES:
        return _DRAINAGE_OVERRIDES[name]  # type: ignore[return-value]
    if name in PORTAL_TISSUES:
        return "liver"
    return "venous"


def load_species_topology(
    species: str,
    *,
    path: Path | str | None = None,
) -> PBPKTopology:
    """Load a PBPK topology from ``src/charon/pbpk/species/<name>.yaml``.

    Parameters
    ----------
    species : str
        Species name used as the YAML filename stem.
    path : Path or str, optional
        Explicit path to the YAML file.  If supplied, ``species`` is used
        only as the returned ``PBPKTopology.species`` identifier.

    Returns
    -------
    PBPKTopology
        Frozen topology container.

    Raises
    ------
    FileNotFoundError
        If the YAML file does not exist.
    ValueError
        If the topology fails internal consistency checks (portal sum
        exceeds total liver flow, missing mandatory tissues, etc.).
    """
    if path is None:
        yaml_path = _SPECIES_DIR / f"{species}.yaml"
    else:
        yaml_path = Path(path)

    if not yaml_path.exists():
        raise FileNotFoundError(
            f"Species YAML not found: {yaml_path}"
        )

    with yaml_path.open() as fp:
        data = yaml.safe_load(fp)

    spec = data["species"]
    co_L_h = float(spec["cardiac_output_L_h"])

    tissues: OrderedDict[str, TissueNode] = OrderedDict()
    for tissue_name, tissue_spec in spec["tissues"].items():
        flow_fraction = float(tissue_spec["blood_flow_fraction"])
        if tissue_name == "lung":
            q = co_L_h  # pulmonary: full CO
        else:
            q = flow_fraction * co_L_h
        node = TissueNode(
            name=tissue_name,
            volume_L=float(tissue_spec["volume_L"]),
            blood_flow_L_h=q,
            composition=_tissue_comp_from_dict(tissue_spec["composition"]),
            drains_to=_drainage_for(tissue_name),
        )
        tissues[tissue_name] = node

    # Mandatory tissues for the minimal IV-only kernel.
    for required in ("lung", "liver", "kidney", *PORTAL_TISSUES):
        if required not in tissues:
            raise ValueError(
                f"species {species!r} topology is missing required tissue {required!r}"
            )

    # Hepatic artery = Q_liver_total - Q_portal_total
    q_liver_total = tissues["liver"].blood_flow_L_h
    q_portal_total = sum(tissues[p].blood_flow_L_h for p in PORTAL_TISSUES)
    q_hepatic_artery = q_liver_total - q_portal_total
    if q_hepatic_artery < 0:
        raise ValueError(
            f"species {species!r}: sum of portal tissue flows "
            f"({q_portal_total:.3g} L/h) exceeds total liver flow "
            f"({q_liver_total:.3g} L/h); topology is inconsistent"
        )

    plasma_comp = _tissue_comp_from_dict(spec["plasma"])

    return PBPKTopology(
        species=species,
        body_weight_kg=float(spec["body_weight_kg"]),
        cardiac_output_L_h=co_L_h,
        hematocrit=float(spec["hematocrit"]),
        venous_volume_L=float(spec["blood"]["venous_volume_L"]),
        arterial_volume_L=float(spec["blood"]["arterial_volume_L"]),
        hepatic_artery_L_h=q_hepatic_artery,
        tissues=tissues,
        plasma_composition=plasma_comp,
        gfr_mL_min=float(spec["kidney"]["gfr_mL_min"]),
        liver_weight_g=float(spec["liver"]["weight_g"]),
    )
```

- [ ] **Step 4: Run the tests to verify pass**

```
python3 -m pytest tests/unit/test_topology.py -v
```

Expected: all tests pass.

- [ ] **Step 5: Commit**

```
git add src/charon/pbpk/topology.py tests/unit/test_topology.py
git commit -m "Sprint 3: Add pbpk.topology — species YAML → PBPKTopology graph"
```

---

## Task 3: pbpk/ode_compiler.py — CompoundPBPKParams + build_rhs

Takes a topology + compound + parameter bridge → computes Kp per tissue and assembles an ODE right-hand side function closed over those parameters. This is where the perfusion-limited mass balance equations live.

**Files:**
- Create: `src/charon/pbpk/ode_compiler.py`
- Create: `tests/unit/test_ode_compiler.py`

- [ ] **Step 1: Write failing tests**

Create `tests/unit/test_ode_compiler.py`:

```python
"""Unit tests for the PBPK ODE compiler.

The tests follow a fixture ladder:
  1. Build a midazolam-like CompoundPBPKParams via the bridge + topology.
  2. Inspect the Kp dictionary, fu_b, CLint_liver_L_h values.
  3. Verify rhs returns finite derivatives.
  4. Verify mass conservation when elimination is disabled.
"""

import math
from collections import OrderedDict

import numpy as np
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
from charon.pbpk.kp_calculator import TissueComposition
from charon.pbpk.ode_compiler import (
    CompoundPBPKParams,
    build_compound_pbpk_params,
    build_rhs,
    infer_compound_type,
)
from charon.pbpk.topology import (
    PBPKTopology,
    PORTAL_TISSUES,
    load_species_topology,
)


def _predicted(value: float, unit: str = "") -> PredictedProperty:
    return PredictedProperty(value=float(value), source="experimental", unit=unit)


@pytest.fixture
def human_topology() -> PBPKTopology:
    return load_species_topology("human")


@pytest.fixture
def midazolam_compound() -> CompoundConfig:
    """Experimental-override midazolam fixture matching test_parameter_bridge."""
    return CompoundConfig(
        name="midazolam",
        smiles="CC1=NCC2=NC=CN2C1(c1ccc(F)cc1)",  # approximate
        molecular_weight=325.77,
        source="experimental",
        properties=CompoundProperties(
            physicochemical=PhysicochemicalProperties(
                logp=_predicted(3.89),
                pka_base=_predicted(6.2),
            ),
            binding=BindingProperties(
                fu_p=_predicted(0.03, "fraction"),
                fu_inc=_predicted(0.96, "fraction"),
                bp_ratio=_predicted(0.66, "ratio"),
            ),
            metabolism=MetabolismProperties(
                clint_uL_min_mg=_predicted(93.0, "uL/min/mg"),
            ),
            renal=RenalProperties(
                clrenal_L_h=_predicted(0.0, "L/h"),
            ),
        ),
    )


@pytest.fixture
def midazolam_params(human_topology, midazolam_compound) -> CompoundPBPKParams:
    return build_compound_pbpk_params(
        midazolam_compound,
        human_topology,
        bridge=ParameterBridge(),
        compound_type="base",
        clint_system="HLM",
    )


class TestInferCompoundType:
    @pytest.mark.parametrize(
        "pka_acid, pka_base, expected",
        [
            (None, None, "neutral"),
            (4.5, None, "acid"),
            (None, 9.5, "base"),
            (4.5, 9.5, "zwitterion"),
            (9.0, None, "neutral"),  # pKa_acid >= 7.0 → not ionized at pH 7.4
            (None, 7.0, "neutral"),  # pKa_base <= 8.0 → not basic enough
        ],
    )
    def test_classification(self, pka_acid, pka_base, expected):
        assert infer_compound_type(pka_acid, pka_base) == expected


class TestBuildCompoundPBPKParams:
    def test_compound_type_override(self, midazolam_params):
        # We explicitly passed "base" to override the pKa threshold rule.
        assert midazolam_params.compound_type == "base"

    def test_fu_b_calculation(self, midazolam_params):
        # fu_b = fu_p / BP = 0.03 / 0.66
        assert midazolam_params.fu_b == pytest.approx(0.03 / 0.66, rel=1e-6)

    def test_clint_liver_L_h_matches_handcalc(self, midazolam_params):
        # CLint_u = 93/0.96 = 96.875 uL/min/mg
        # scale = 40 * 1500 = 60000
        # total_uL_min = 96.875 * 60000 = 5_812_500
        # L/h = 5_812_500 / 1e6 * 60 = 348.75
        assert midazolam_params.clint_liver_L_h == pytest.approx(348.75, rel=1e-3)

    def test_kp_covers_all_topology_tissues(self, midazolam_params, human_topology):
        assert set(midazolam_params.kp_by_tissue.keys()) == set(human_topology.tissues.keys())

    def test_kp_values_finite_and_positive(self, midazolam_params):
        for tissue, kp in midazolam_params.kp_by_tissue.items():
            assert math.isfinite(kp) and kp > 0, f"{tissue}: kp={kp}"

    def test_cl_renal_zero_for_midazolam(self, midazolam_params):
        assert midazolam_params.cl_renal_L_h == pytest.approx(0.0)


class TestBuildRhs:
    def test_rhs_at_initial_state(self, human_topology, midazolam_params):
        rhs = build_rhs(human_topology, midazolam_params)
        n_tissues = len(human_topology.tissues)
        y0 = np.zeros(2 + n_tissues)
        y0[0] = 5.0  # 5 mg IV bolus in venous pool
        dy = rhs(0.0, y0)
        assert dy.shape == y0.shape
        assert np.all(np.isfinite(dy))

    def test_mass_balance_no_elimination(self, human_topology, midazolam_compound):
        """With CLint=0 and CLrenal=0 the total body burden is conserved."""
        # Construct a zero-clearance variant
        compound = midazolam_compound.model_copy(deep=True)
        compound.properties.metabolism.clint_uL_min_mg = _predicted(0.0, "uL/min/mg")
        compound.properties.renal.clrenal_L_h = _predicted(0.0, "L/h")
        params = build_compound_pbpk_params(
            compound,
            human_topology,
            bridge=ParameterBridge(),
            compound_type="base",
            clint_system="HLM",
        )
        rhs = build_rhs(human_topology, params)
        # Inject an initial distribution (non-trivial so rhs has work)
        rng = np.random.default_rng(42)
        y = rng.uniform(0.1, 1.0, size=2 + len(human_topology.tissues))
        dy = rhs(0.0, y)
        # Sum of dy should be zero (no net source/sink) to machine precision
        assert abs(np.sum(dy)) < 1e-9

    def test_mass_balance_with_liver_elim(self, human_topology, midazolam_params):
        """With CLint>0 the total time derivative equals -hepatic_elim_rate."""
        rhs = build_rhs(human_topology, midazolam_params)
        n_tissues = len(human_topology.tissues)
        y = np.zeros(2 + n_tissues)
        # Put some drug in every compartment so the rhs is non-trivial
        y[0] = 2.0   # venous
        y[1] = 1.5   # arterial
        y[2:] = 0.3
        dy = rhs(0.0, y)
        # Expected total rate = -(CLint_liver * fu_b * C_liver_blood_out)
        liver_idx = list(human_topology.tissues.keys()).index("liver") + 2
        liver_node = human_topology.tissues["liver"]
        kp_liver = midazolam_params.kp_by_tissue["liver"]
        c_tissue_liver = y[liver_idx] / liver_node.volume_L
        c_blood_out_liver = c_tissue_liver * midazolam_params.bp_ratio / kp_liver
        expected_elim_rate = (
            -midazolam_params.clint_liver_L_h
            * midazolam_params.fu_b
            * c_blood_out_liver
        )
        # CLrenal = 0 so only hepatic contributes
        assert np.sum(dy) == pytest.approx(expected_elim_rate, rel=1e-6, abs=1e-9)

    def test_rhs_zero_state_returns_zero(self, human_topology, midazolam_params):
        rhs = build_rhs(human_topology, midazolam_params)
        n_tissues = len(human_topology.tissues)
        y = np.zeros(2 + n_tissues)
        dy = rhs(0.0, y)
        assert np.allclose(dy, 0.0, atol=1e-15)

    def test_renal_elimination_active_for_nonzero_cl_renal(
        self, human_topology, midazolam_compound
    ):
        """Kidney compartment consumes drug when CL_renal > 0."""
        compound = midazolam_compound.model_copy(deep=True)
        compound.properties.metabolism.clint_uL_min_mg = _predicted(0.0, "uL/min/mg")
        # Override renal clearance to 5 L/h
        compound.properties.renal.clrenal_L_h = _predicted(5.0, "L/h")
        params = build_compound_pbpk_params(
            compound,
            human_topology,
            bridge=ParameterBridge(),
            compound_type="base",
            clint_system="HLM",
            override_cl_renal_L_h=5.0,
        )
        rhs = build_rhs(human_topology, params)
        n_tissues = len(human_topology.tissues)
        y = np.zeros(2 + n_tissues)
        # Load kidney so it has blood to eliminate
        kidney_idx = list(human_topology.tissues.keys()).index("kidney") + 2
        y[kidney_idx] = 1.0
        dy = rhs(0.0, y)
        # Total should be strictly negative: only elimination is in kidney
        assert np.sum(dy) < 0.0
```

- [ ] **Step 2: Run to verify failure**

```
python3 -m pytest tests/unit/test_ode_compiler.py -v
```

Expected: ImportError (ode_compiler does not exist).

- [ ] **Step 3: Implement ode_compiler.py**

Create `src/charon/pbpk/ode_compiler.py`:

```python
"""Compile a PBPK topology + compound into an ODE right-hand side function.

The ``build_compound_pbpk_params`` function assembles every numeric
parameter the ODE needs (Kp per tissue, fu_b, CLint_liver_L_h, CL_renal)
via :class:`charon.core.parameter_bridge.ParameterBridge` — so the full
IVIVE audit trail is still produced by the bridge on every Pipeline run.

``build_rhs`` closes over the topology + params and returns a
``rhs(t, y) -> dy`` function suitable for :func:`scipy.integrate.solve_ivp`.

Units (all interfaces):
  * amounts A : mg
  * volumes V : L
  * concentrations C : mg/L (equivalent to ng/µL)
  * flows Q : L/h
  * clearances : L/h
  * time : h
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np

from charon.core.parameter_bridge import ParameterBridge
from charon.core.schema import CompoundConfig
from charon.core.units import mL_min_to_L_h
from charon.pbpk.kp_calculator import compute_all_kp
from charon.pbpk.topology import PBPKTopology, PORTAL_TISSUES

_COMPOUND_TYPES = ("neutral", "acid", "base", "zwitterion")


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


def infer_compound_type(pka_acid: float | None, pka_base: float | None) -> str:
    """Classify a compound from its pKa values (matches Sisyphus thresholds).

    Rules:
      - acidic pKa < 7.0 AND basic pKa > 8.0  → zwitterion
      - acidic pKa < 7.0 only                  → acid
      - basic pKa > 8.0 only                   → base
      - otherwise                               → neutral
    """
    is_acidic = pka_acid is not None and pka_acid < 7.0
    is_basic = pka_base is not None and pka_base > 8.0
    if is_acidic and is_basic:
        return "zwitterion"
    if is_acidic:
        return "acid"
    if is_basic:
        return "base"
    return "neutral"


def _get_property_value(prop) -> float | None:
    if prop is None:
        return None
    return float(prop.value)


def _require(value: float | None, label: str) -> float:
    if value is None:
        raise ValueError(
            f"{label} is required to build PBPK parameters but was not provided"
        )
    return value


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

    Parameters
    ----------
    compound : CompoundConfig
        Validated compound with at minimum: logP, fu_p, fu_inc, BP ratio,
        CLint and molecular weight.
    topology : PBPKTopology
        Loaded species topology (tissue volumes, flows, compositions).
    bridge : ParameterBridge
        Active IVIVE bridge instance; used to run CLint→CLh and derive
        CLint_liver_L_h with a full ConversionLog audit.
    compound_type : str, optional
        Override the pKa-based compound_type classification.  Useful for
        weak bases (e.g. midazolam, pKa_base=6.2) where the Sisyphus
        threshold of >8.0 would misclassify the drug as neutral.
    clint_system : "HLM" or "hepatocytes"
        Which in-vitro system the CLint value comes from.  HLM applies
        fu_inc correction; hepatocytes do not.
    liver_model : str
        Liver model name passed to the bridge (for the audit trail only —
        the PBPK ODE consumes CLint_liver_L_h regardless).
    override_cl_renal_L_h : float, optional
        If supplied, use this value instead of running the bridge's renal
        clearance estimator.  Enables compounds with experimental renal CL.
    """
    props = compound.properties
    logp = _require(_get_property_value(props.physicochemical.logp), "logp")
    pka_acid = _get_property_value(props.physicochemical.pka_acid)
    pka_base = _get_property_value(props.physicochemical.pka_base)
    fu_p = _require(_get_property_value(props.binding.fu_p), "fu_p")
    fu_inc = _require(_get_property_value(props.binding.fu_inc), "fu_inc")
    bp_ratio = _require(_get_property_value(props.binding.bp_ratio), "bp_ratio")
    clint = _require(
        _get_property_value(props.metabolism.clint_uL_min_mg), "clint_uL_min_mg"
    )
    mw = compound.molecular_weight
    if mw is None:
        raise ValueError(
            f"compound {compound.name!r} has no molecular_weight; required for PBPK"
        )

    resolved_type = compound_type or infer_compound_type(pka_acid, pka_base)
    if resolved_type not in _COMPOUND_TYPES:
        raise ValueError(
            f"compound_type must be one of {_COMPOUND_TYPES}, got {resolved_type!r}"
        )

    # For Kp the pKa passed depends on which branch fires.
    pka_for_kp: float | None
    if resolved_type == "acid":
        pka_for_kp = pka_acid
    elif resolved_type in ("base", "zwitterion"):
        pka_for_kp = pka_base
    else:
        pka_for_kp = None

    # Assemble tissue compositions for Kp calc.
    tissue_comps = {
        name: node.composition for name, node in topology.tissues.items()
    }

    kp_by_tissue = compute_all_kp(
        logp=logp,
        pka=pka_for_kp,
        compound_type=resolved_type,
        tissue_compositions=tissue_comps,
        plasma_composition=topology.plasma_composition,
        method="rodgers_rowland",
    )

    # Hepatic IVIVE with full audit via ParameterBridge
    hep = bridge.clint_to_clh(
        clint=clint,
        fu_inc=fu_inc,
        fu_p=fu_p,
        bp_ratio=bp_ratio,
        system=clint_system,
        mppgl=40.0 if clint_system == "HLM" else 40.0,  # HLM only
        hepatocellularity=120.0,
        liver_weight_g=topology.liver_weight_g,
        qh_L_h=topology.tissues["liver"].blood_flow_L_h,
        model=liver_model,
    )
    clint_liver_L_h = hep.clint_liver_L_h
    assert clint_liver_L_h is not None  # populated by Task 1

    # Renal clearance
    if override_cl_renal_L_h is not None:
        cl_renal_L_h = float(override_cl_renal_L_h)
    else:
        # Use the compound's own value if supplied, otherwise use bridge default.
        cl_renal_prop = _get_property_value(props.renal.clrenal_L_h)
        if cl_renal_prop is not None:
            cl_renal_L_h = cl_renal_prop
        else:
            cl_renal_L_h = bridge.assign_renal_clearance(
                fu_p=fu_p,
                gfr_mL_min=topology.gfr_mL_min,
            )

    fu_b = fu_p / bp_ratio

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
    )


def build_rhs(
    topology: PBPKTopology,
    params: CompoundPBPKParams,
    *,
    infusion_rate_mg_per_h: float = 0.0,
    infusion_duration_h: float = 0.0,
):
    """Return a closure ``rhs(t, y) -> dy`` for the PBPK ODE system.

    State vector layout (length 2 + N_tissues):
      y[0] = A_venous      (mg)
      y[1] = A_arterial    (mg)
      y[2:] = A_tissue_i   (mg) in topology.tissue_names() order

    The returned function is a plain Python callable with no shared state,
    safe to pass to :func:`scipy.integrate.solve_ivp`.  Jacobian is left
    for scipy to estimate via finite differences (PBPK with 17 states is
    small enough).

    Parameters
    ----------
    infusion_rate_mg_per_h : float
        Source term added to the venous pool while ``t <
        infusion_duration_h``.  Zero for IV bolus.
    infusion_duration_h : float
        End of the infusion in hours.  Only used when
        ``infusion_rate_mg_per_h > 0``.
    """
    tissue_names = topology.tissue_names()
    n_tissues = len(tissue_names)

    # Cache index and scalar arrays
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

    # Precompute masks (boolean) for rhs assembly
    # systemic_into_venous: non-portal, non-lung, non-liver tissues that drain to venous
    systemic_into_venous = np.zeros(n_tissues, dtype=bool)
    for i, name in enumerate(tissue_names):
        if drains_to[i] == "venous" and name not in {"liver"}:
            # lung drains to arterial; liver drains to venous but is handled separately
            systemic_into_venous[i] = True

    bp = params.bp_ratio
    fu_b = params.fu_b
    clint_liver = params.clint_liver_L_h
    cl_renal = params.cl_renal_L_h

    v_ven = topology.venous_volume_L
    v_art = topology.arterial_volume_L
    q_co = topology.cardiac_output_L_h
    q_ha = topology.hepatic_artery_L_h
    q_liver_total = flows[liver_idx]

    def rhs(t: float, y: np.ndarray) -> np.ndarray:
        dy = np.empty_like(y)

        a_ven = y[0]
        a_art = y[1]
        a_tissues = y[2:]

        c_ven = a_ven / v_ven
        c_art = a_art / v_art

        # Blood-equivalent concentration leaving each tissue (perfusion-limited)
        c_tissue = a_tissues / volumes
        c_blood_out = c_tissue * bp / kp

        # Systemic parallel tissues (draining to venous), excluding liver
        # dA_i/dt = Q_i * (C_art - C_out_i)
        dy_tissues = np.zeros_like(a_tissues)
        for i in range(n_tissues):
            if i == lung_idx or i == liver_idx:
                continue
            dy_tissues[i] = flows[i] * (c_art - c_blood_out[i])

        # Lung: receives full CO from venous, outputs to arterial
        dy_tissues[lung_idx] = q_co * (c_ven - c_blood_out[lung_idx])

        # Kidney: add renal elimination (CL_renal is plasma clearance → /BP)
        dy_tissues[kidney_idx] -= cl_renal * c_blood_out[kidney_idx] / bp

        # Liver: hepatic artery + portal inflow - total outflow - hepatic elim
        portal_inflow = 0.0
        for pi in portal_indices:
            portal_inflow += flows[pi] * c_blood_out[pi]
        hepatic_elim = clint_liver * fu_b * c_blood_out[liver_idx]
        dy_tissues[liver_idx] = (
            q_ha * c_art
            + portal_inflow
            - q_liver_total * c_blood_out[liver_idx]
            - hepatic_elim
        )

        # Venous pool: sum of all tissue outflows that drain to venous, minus CO out
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

        # Optional infusion source
        if infusion_rate_mg_per_h > 0.0 and t < infusion_duration_h:
            dy_ven += infusion_rate_mg_per_h

        # Arterial pool: fed by lung only, drained by CO
        dy_art = q_co * c_blood_out[lung_idx] - q_co * c_art

        dy[0] = dy_ven
        dy[1] = dy_art
        dy[2:] = dy_tissues
        return dy

    return rhs
```

- [ ] **Step 4: Run tests to verify pass**

```
python3 -m pytest tests/unit/test_ode_compiler.py -v
```

Expected: all pass.

- [ ] **Step 5: Commit**

```
git add src/charon/pbpk/ode_compiler.py tests/unit/test_ode_compiler.py
git commit -m "Sprint 3: Add pbpk.ode_compiler — PBPK params + perfusion-limited rhs"
```

---

## Task 4: pbpk/solver.py — SimulationResult + simulate_iv

Wraps `scipy.solve_ivp` with BDF for stiff PBPK systems. Enforces the solver method, injects the initial condition, and packages the result.

**Files:**
- Create: `src/charon/pbpk/solver.py`
- Create: `tests/unit/test_solver.py`

- [ ] **Step 1: Write failing tests**

Create `tests/unit/test_solver.py`:

```python
"""Unit tests for the PBPK solver wrapper.

Strategy:
  1. Smoke test: midazolam IV bolus runs to completion with finite output.
  2. Mass conservation: no elimination → total mass constant (numerical).
  3. Analytic 1-cpt match: collapse the model to a 1-compartment PK
     equivalent and verify exponential decay matches within rtol.
  4. BDF enforcement: simulate_iv rejects non-BDF methods.
"""

import numpy as np
import pytest
from scipy.integrate import solve_ivp

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
from charon.pbpk.solver import SimulationResult, simulate_iv
from charon.pbpk.topology import load_species_topology


def _predicted(value: float, unit: str = "") -> PredictedProperty:
    return PredictedProperty(value=float(value), source="experimental", unit=unit)


@pytest.fixture
def human():
    return load_species_topology("human")


@pytest.fixture
def midazolam_cfg():
    return CompoundConfig(
        name="midazolam",
        smiles="C",
        molecular_weight=325.77,
        source="experimental",
        properties=CompoundProperties(
            physicochemical=PhysicochemicalProperties(
                logp=_predicted(3.89),
                pka_base=_predicted(6.2),
            ),
            binding=BindingProperties(
                fu_p=_predicted(0.03, "fraction"),
                fu_inc=_predicted(0.96, "fraction"),
                bp_ratio=_predicted(0.66, "ratio"),
            ),
            metabolism=MetabolismProperties(
                clint_uL_min_mg=_predicted(93.0, "uL/min/mg"),
            ),
            renal=RenalProperties(
                clrenal_L_h=_predicted(0.0, "L/h"),
            ),
        ),
    )


@pytest.fixture
def midazolam_params(human, midazolam_cfg):
    return build_compound_pbpk_params(
        midazolam_cfg,
        human,
        bridge=ParameterBridge(),
        compound_type="base",
        override_cl_renal_L_h=0.0,
    )


class TestSimulateIvBolus:
    def test_smoke(self, human, midazolam_params):
        result = simulate_iv(
            human,
            midazolam_params,
            dose_mg=5.0,
            route="iv_bolus",
            duration_h=24.0,
        )
        assert isinstance(result, SimulationResult)
        assert result.time_h[0] == 0.0
        assert result.time_h[-1] == pytest.approx(24.0)
        assert np.all(np.isfinite(result.cp_plasma))
        assert result.cp_plasma[0] > 0.0, "Venous Cp should jump at t=0 from bolus"
        assert result.solver_success is True
        assert result.solver_method == "BDF"

    def test_bolus_mass_conservation_no_elim(self, human, midazolam_cfg):
        """Zero CLint + zero CLrenal → total drug mass conserved over time."""
        cfg = midazolam_cfg.model_copy(deep=True)
        cfg.properties.metabolism.clint_uL_min_mg = _predicted(0.0, "uL/min/mg")
        params = build_compound_pbpk_params(
            cfg, human, bridge=ParameterBridge(), compound_type="base",
            override_cl_renal_L_h=0.0,
        )
        result = simulate_iv(
            human, params, dose_mg=5.0, route="iv_bolus", duration_h=24.0
        )
        totals = result.state_trajectory.sum(axis=0)
        max_err = np.max(np.abs(totals - 5.0))
        assert max_err < 1e-3, f"Mass drift {max_err:.3e} mg exceeds tolerance"

    def test_mass_balance_residual_reported(self, human, midazolam_params):
        result = simulate_iv(
            human, midazolam_params, dose_mg=5.0,
            route="iv_bolus", duration_h=24.0,
        )
        # With elimination mass decreases — residual vs dose is not zero,
        # but the attribute should still be a finite float.
        assert np.isfinite(result.mass_balance_residual)

    def test_bdf_enforced(self, human, midazolam_params):
        with pytest.raises(ValueError, match="BDF"):
            simulate_iv(
                human, midazolam_params, dose_mg=5.0,
                route="iv_bolus", duration_h=24.0,
                method="RK45",
            )

    def test_bad_route_raises(self, human, midazolam_params):
        with pytest.raises(ValueError, match="route"):
            simulate_iv(
                human, midazolam_params, dose_mg=5.0,
                route="oral", duration_h=24.0,
            )


class TestSimulateIvInfusion:
    def test_infusion_matches_bolus_at_long_time(self, human, midazolam_params):
        """A 1-minute infusion of 5 mg should give nearly the same 24 h
        profile as a 5 mg bolus (aside from the first few minutes)."""
        bolus = simulate_iv(
            human, midazolam_params, dose_mg=5.0,
            route="iv_bolus", duration_h=24.0,
        )
        infusion = simulate_iv(
            human, midazolam_params, dose_mg=5.0,
            route="iv_infusion", duration_h=24.0,
            infusion_duration_h=1.0 / 60.0,
        )
        # Compare Cp at t = 6 h (well past the tiny infusion window)
        idx = np.argmin(np.abs(bolus.time_h - 6.0))
        cp_b = bolus.cp_plasma[idx]
        cp_i = infusion.cp_plasma[idx]
        assert cp_i == pytest.approx(cp_b, rel=0.05)

    def test_infusion_requires_positive_duration(self, human, midazolam_params):
        with pytest.raises(ValueError, match="infusion_duration"):
            simulate_iv(
                human, midazolam_params, dose_mg=5.0,
                route="iv_infusion", duration_h=24.0,
                infusion_duration_h=0.0,
            )


class TestOneCompartmentAnalyticMatch:
    """Sanity: collapse PBPK into a 1-cpt equivalent and match exponential decay.

    We do this by building an artificial "compound" with:
      - zero CL_renal
      - zero CLint (so hepatic elim is zero)
      - Kp = 1 everywhere (so tissues are in trivial equilibrium with blood)
      - BP = 1
    Then the venous pool concentration should evolve as a pure mixing system
    — i.e., after distribution, C_ven stays at dose / V_total_apparent.

    For this check we verify:
      1. The total mass stays at dose (already covered above).
      2. After a long time, C_ven converges to a steady plateau equal to
         dose / (V_ven + V_art + sum(V_tissue)).
    """

    def test_uniform_kp_gives_mass_over_total_volume(self, human, midazolam_cfg):
        # Rebuild params with zero elimination; override Kp to 1.0 everywhere.
        cfg = midazolam_cfg.model_copy(deep=True)
        cfg.properties.metabolism.clint_uL_min_mg = _predicted(0.0, "uL/min/mg")
        params = build_compound_pbpk_params(
            cfg, human, bridge=ParameterBridge(), compound_type="base",
            override_cl_renal_L_h=0.0,
        )
        # Rewrite Kp dict in-place — params is frozen, so build a replacement
        # via dataclasses.replace.
        from dataclasses import replace
        kp_uniform = {name: 1.0 for name in params.kp_by_tissue}
        params_uniform = replace(
            params,
            bp_ratio=1.0,
            fu_b=params.fu_p / 1.0,
            kp_by_tissue=kp_uniform,
        )
        result = simulate_iv(
            human, params_uniform, dose_mg=5.0,
            route="iv_bolus", duration_h=48.0,
        )
        v_total = (
            human.venous_volume_L
            + human.arterial_volume_L
            + sum(n.volume_L for n in human.tissues.values())
        )
        expected_plateau = 5.0 / v_total  # mg/L blood; plasma = same since BP=1
        cp_late = result.cp_plasma[-1]
        assert cp_late == pytest.approx(expected_plateau, rel=0.01)
```

- [ ] **Step 2: Run to verify failure**

```
python3 -m pytest tests/unit/test_solver.py -v
```

Expected: ImportError — solver module missing.

- [ ] **Step 3: Implement solver.py**

Create `src/charon/pbpk/solver.py`:

```python
"""Stiff ODE solver wrapper for PBPK simulation.

The PBPK system is stiff (fast tissue equilibration vs slow elimination),
so this module mandates ``scipy.integrate.solve_ivp(method='BDF')``.
Explicit methods (RK45, RK23, DOP853) are rejected with ValueError.

Usage
-----

>>> topology = load_species_topology("human")
>>> params = build_compound_pbpk_params(compound, topology, bridge)
>>> result = simulate_iv(topology, params, dose_mg=5.0,
...                      route="iv_bolus", duration_h=72.0)
>>> cp_plasma = result.cp_plasma  # venous plasma concentration, mg/L
>>> time_h = result.time_h
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np
from scipy.integrate import solve_ivp

from charon.pbpk.ode_compiler import CompoundPBPKParams, build_rhs
from charon.pbpk.topology import PBPKTopology

_ALLOWED_METHODS = {"BDF"}  # explicitly rejects RK45, LSODA, etc.


@dataclass
class SimulationResult:
    """Output of an IV PBPK simulation."""

    time_h: np.ndarray          # (N,)
    cp_blood: np.ndarray        # (N,) venous blood concentration, mg/L
    cp_plasma: np.ndarray       # (N,) venous plasma concentration, mg/L
    state_trajectory: np.ndarray  # (n_states, N)
    mass_balance_residual: float
    solver_success: bool
    solver_method: str
    solver_nfev: int
    route: str
    dose_mg: float
    infusion_duration_h: float


def simulate_iv(
    topology: PBPKTopology,
    params: CompoundPBPKParams,
    *,
    dose_mg: float,
    route: Literal["iv_bolus", "iv_infusion"],
    duration_h: float,
    infusion_duration_h: float = 0.0,
    n_time_points: int = 500,
    rtol: float = 1e-6,
    atol: float | None = None,
    method: str = "BDF",
) -> SimulationResult:
    """Simulate an IV dose through the human PBPK kernel.

    Parameters
    ----------
    topology : PBPKTopology
        Loaded species topology.
    params : CompoundPBPKParams
        Compound-specific PBPK parameters (Kp, fu_b, CLint_liver_L_h, etc.).
    dose_mg : float
        Dose in mg.
    route : "iv_bolus" or "iv_infusion"
        "iv_bolus" places all dose into the venous pool at t=0.
        "iv_infusion" injects at a constant rate over infusion_duration_h.
    duration_h : float
        Simulation end time (hours).
    infusion_duration_h : float, optional
        Required for "iv_infusion" (must be > 0).
    n_time_points : int
        Number of evaluation points (uniform grid on [0, duration_h]).
    rtol, atol : float
        Solver tolerances.  ``atol`` defaults to ``1e-9 * dose_mg``.
    method : str
        Must be "BDF" — explicit methods are rejected (PBPK is stiff).

    Returns
    -------
    SimulationResult
    """
    if method not in _ALLOWED_METHODS:
        raise ValueError(
            f"PBPK is stiff — only BDF is allowed, got {method!r}. "
            f"Explicit methods (RK45, DOP853, LSODA) silently produce wrong "
            f"results for PBPK systems."
        )

    if route not in ("iv_bolus", "iv_infusion"):
        raise ValueError(
            f"route must be 'iv_bolus' or 'iv_infusion', got {route!r}"
        )

    if dose_mg <= 0:
        raise ValueError(f"dose_mg must be > 0, got {dose_mg}")

    if duration_h <= 0:
        raise ValueError(f"duration_h must be > 0, got {duration_h}")

    if route == "iv_infusion" and infusion_duration_h <= 0:
        raise ValueError(
            f"iv_infusion requires infusion_duration_h > 0, got {infusion_duration_h}"
        )

    if atol is None:
        atol = 1e-9 * dose_mg

    n_tissues = len(topology.tissues)
    y0 = np.zeros(2 + n_tissues)

    if route == "iv_bolus":
        y0[0] = dose_mg
        infusion_rate = 0.0
        inf_duration = 0.0
    else:  # iv_infusion
        infusion_rate = dose_mg / infusion_duration_h
        inf_duration = infusion_duration_h

    rhs = build_rhs(
        topology,
        params,
        infusion_rate_mg_per_h=infusion_rate,
        infusion_duration_h=inf_duration,
    )

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
            f"PBPK ODE solver failed for compound {params.name!r}: {sol.message}"
        )

    v_ven = topology.venous_volume_L
    cp_blood = sol.y[0] / v_ven
    cp_plasma = cp_blood / params.bp_ratio

    # Mass balance residual: for bolus, dose should equal sum(y) + cumulative elim.
    # Here we just expose |sum(y) at t=0 - dose_mg| as a sanity check; full
    # cumulative elim tracking is deferred to a later session (needs extra state).
    residual = float(np.abs(sol.y.sum(axis=0)[0] - dose_mg))

    return SimulationResult(
        time_h=sol.t,
        cp_blood=cp_blood,
        cp_plasma=cp_plasma,
        state_trajectory=sol.y,
        mass_balance_residual=residual,
        solver_success=sol.success,
        solver_method=method,
        solver_nfev=int(sol.nfev),
        route=route,
        dose_mg=dose_mg,
        infusion_duration_h=inf_duration,
    )
```

- [ ] **Step 4: Run solver tests**

```
python3 -m pytest tests/unit/test_solver.py -v
```

Expected: all pass.

- [ ] **Step 5: Commit**

```
git add src/charon/pbpk/solver.py tests/unit/test_solver.py
git commit -m "Sprint 3: Add pbpk.solver — BDF ODE wrapper with IV bolus/infusion"
```

---

## Task 5: pbpk/pk_extract.py — PK parameter extraction

From (time, Cp_plasma) → Cmax, Tmax, AUC, half-life, CL, Vss. Used downstream by Pipeline.

**Files:**
- Create: `src/charon/pbpk/pk_extract.py`
- Create: `tests/unit/test_pk_extract.py`

- [ ] **Step 1: Write failing tests**

Create `tests/unit/test_pk_extract.py`:

```python
"""Unit tests for PK parameter extraction from Cp-time profiles."""

import math

import numpy as np
import pytest

from charon.core.schema import PKParameters
from charon.pbpk.pk_extract import compute_pk_parameters


class TestAnalyticOneCompartmentBolus:
    """Mono-exponential IV bolus: C(t) = (D/V) * exp(-ke*t).

    With V = 10 L and CL = 5 L/h → ke = 0.5 h⁻¹, t½ = ln(2)/0.5 ≈ 1.386 h.
    A 100 mg bolus gives Cp(0) = 10 mg/L.

    Integrals:
      AUC_0_inf = D / CL = 100 / 5 = 20.0 mg·h/L
      AUMC_0_inf = D / (CL * ke) = 100 / (5 * 0.5) = 40.0 mg·h²/L
      MRT = AUMC / AUC = 2.0 h  (also = 1/ke)
      Vss = CL * MRT = 5 * 2 = 10.0 L (matches V)
    """

    @pytest.fixture
    def mono_exp(self):
        ke = 0.5
        dose = 100.0
        v = 10.0
        cp0 = dose / v
        t = np.linspace(0, 48.0, 2001)
        cp = cp0 * np.exp(-ke * t)
        return t, cp, dose, ke, v

    def test_cmax_equals_cp0(self, mono_exp):
        t, cp, dose, _, v = mono_exp
        pk = compute_pk_parameters(t, cp, dose_mg=dose, route="iv_bolus")
        assert pk.cmax == pytest.approx(dose / v, rel=1e-4)
        assert pk.tmax == pytest.approx(0.0, abs=1e-6)

    def test_auc_inf(self, mono_exp):
        t, cp, dose, ke, _ = mono_exp
        pk = compute_pk_parameters(t, cp, dose_mg=dose, route="iv_bolus")
        expected_auc = dose / 5.0  # CL = 5 L/h
        assert pk.auc_0_inf == pytest.approx(expected_auc, rel=5e-3)

    def test_half_life(self, mono_exp):
        t, cp, dose, ke, _ = mono_exp
        pk = compute_pk_parameters(t, cp, dose_mg=dose, route="iv_bolus")
        expected_t12 = math.log(2) / ke
        assert pk.half_life == pytest.approx(expected_t12, rel=1e-3)

    def test_cl(self, mono_exp):
        t, cp, dose, _, _ = mono_exp
        pk = compute_pk_parameters(t, cp, dose_mg=dose, route="iv_bolus")
        assert pk.cl_apparent == pytest.approx(5.0, rel=5e-3)

    def test_vss_equals_v(self, mono_exp):
        t, cp, dose, _, v = mono_exp
        pk = compute_pk_parameters(t, cp, dose_mg=dose, route="iv_bolus")
        # Vss = CL * MRT = V (for 1-compartment)
        assert pk.vss == pytest.approx(v, rel=5e-3)

    def test_auc_0_24(self, mono_exp):
        t, cp, dose, ke, _ = mono_exp
        pk = compute_pk_parameters(t, cp, dose_mg=dose, route="iv_bolus")
        # AUC(0,24) = (D/V/ke)*(1 - exp(-ke*24))
        expected = (dose / 10.0 / ke) * (1 - math.exp(-ke * 24.0))
        assert pk.auc_0_24 == pytest.approx(expected, rel=5e-3)


class TestInfusionBenetCorrection:
    """For IV infusion the Vss formula is Vss = CL * (MRT - T_inf/2)."""

    def test_zero_order_infusion_1_compartment(self):
        ke = 0.3
        v = 20.0
        cl = ke * v  # 6 L/h
        tinf = 2.0
        dose = 60.0
        rate = dose / tinf

        # Analytic solution for 1-cpt zero-order infusion:
        #   During infusion (t <= tinf): C = rate/CL * (1 - exp(-ke t))
        #   After infusion (t > tinf): C = (rate/CL) * (1 - exp(-ke*tinf)) * exp(-ke*(t-tinf))
        t = np.linspace(0.0, 48.0, 4801)
        c_inf = (rate / cl) * (1 - np.exp(-ke * np.minimum(t, tinf)))
        c_post = (
            (rate / cl)
            * (1 - np.exp(-ke * tinf))
            * np.exp(-ke * np.maximum(t - tinf, 0.0))
        )
        # Piecewise stitched:
        cp = np.where(t <= tinf, c_inf, c_post)

        pk = compute_pk_parameters(
            t, cp, dose_mg=dose, route="iv_infusion", infusion_duration_h=tinf
        )
        assert pk.cl_apparent == pytest.approx(cl, rel=5e-3)
        # Vss for 1-cpt infusion should still equal V (ignoring the
        # tiny numerical residual).
        assert pk.vss == pytest.approx(v, rel=5e-3)


class TestEdgeCases:
    def test_flat_profile_raises(self):
        t = np.linspace(0, 10, 101)
        cp = np.full_like(t, 1.0)
        with pytest.raises(ValueError, match="terminal slope"):
            compute_pk_parameters(t, cp, dose_mg=10.0, route="iv_bolus")

    def test_mismatched_array_lengths(self):
        with pytest.raises(ValueError, match="length"):
            compute_pk_parameters(
                np.array([0, 1, 2]),
                np.array([1, 2]),
                dose_mg=5.0,
                route="iv_bolus",
            )

    def test_short_tail_falls_back_cleanly(self):
        # Only 4 samples: log-linear fit is degenerate but should not crash.
        t = np.array([0.0, 1.0, 2.0, 3.0])
        cp = np.array([10.0, 5.0, 2.5, 1.25])  # clean exponential
        pk = compute_pk_parameters(t, cp, dose_mg=100.0, route="iv_bolus")
        assert pk.half_life is not None and pk.half_life > 0
```

- [ ] **Step 2: Run to verify failure**

```
python3 -m pytest tests/unit/test_pk_extract.py -v
```

Expected: ImportError.

- [ ] **Step 3: Implement pk_extract.py**

Create `src/charon/pbpk/pk_extract.py`:

```python
"""Extract PK parameters (Cmax, AUC, CL, Vss, t½) from a Cp-time profile."""

from __future__ import annotations

import math
from typing import Literal

import numpy as np

from charon.core.schema import PKParameters


def _trapezoid(x: np.ndarray, y: np.ndarray) -> float:
    """Trapezoidal rule integral — numpy.trapz replacement for newer NumPy."""
    return float(np.trapezoid(y, x))


def _terminal_log_slope(
    t: np.ndarray, cp: np.ndarray, frac: float = 0.3, min_points: int = 3
) -> tuple[float, float]:
    """Log-linear regression on the last ``frac`` of the profile.

    Returns ``(ke, cp_last)``.  Raises ValueError if the slope is
    non-negative (drug not declining → terminal phase undetectable).
    """
    n = len(t)
    n_tail = max(min_points, int(np.ceil(frac * n)))
    n_tail = min(n_tail, n)
    t_tail = t[-n_tail:]
    cp_tail = cp[-n_tail:]
    # Guard against zero/negative concentrations
    mask = cp_tail > 0
    if mask.sum() < 2:
        raise ValueError("terminal slope cannot be fit: fewer than 2 positive samples")
    t_fit = t_tail[mask]
    log_cp = np.log(cp_tail[mask])
    slope, _intercept = np.polyfit(t_fit, log_cp, 1)
    if slope >= 0:
        raise ValueError(
            f"terminal slope is non-negative ({slope:.3e}); simulation "
            f"too short or profile not declining"
        )
    ke = -float(slope)
    return ke, float(cp_tail[mask][-1])


def compute_pk_parameters(
    time_h: np.ndarray,
    cp_plasma: np.ndarray,
    *,
    dose_mg: float,
    route: Literal["iv_bolus", "iv_infusion"],
    infusion_duration_h: float = 0.0,
) -> PKParameters:
    """Compute standard PK parameters from a plasma concentration profile.

    Parameters
    ----------
    time_h : array
        Time grid in hours (monotonically increasing).
    cp_plasma : array
        Plasma concentration at each time point (mg/L).
    dose_mg : float
        Administered dose (mg).
    route : "iv_bolus" or "iv_infusion"
        For "iv_infusion" Vss is corrected via Benet (Vss = CL*(MRT - T_inf/2)).
    infusion_duration_h : float
        Infusion duration (only used when route="iv_infusion").
    """
    t = np.asarray(time_h, dtype=np.float64)
    cp = np.asarray(cp_plasma, dtype=np.float64)
    if t.shape != cp.shape:
        raise ValueError(
            f"time and cp must have the same length, got {t.shape} vs {cp.shape}"
        )
    if t.ndim != 1:
        raise ValueError("time and cp must be 1-D arrays")
    if dose_mg <= 0:
        raise ValueError(f"dose_mg must be > 0, got {dose_mg}")
    if route not in ("iv_bolus", "iv_infusion"):
        raise ValueError(f"route must be iv_bolus or iv_infusion, got {route!r}")

    cmax = float(np.max(cp))
    tmax = float(t[int(np.argmax(cp))])

    auc_0_last = _trapezoid(t, cp)
    # AUC_0_24
    mask_24 = t <= 24.0
    if mask_24.sum() >= 2:
        auc_0_24 = _trapezoid(t[mask_24], cp[mask_24])
    else:
        auc_0_24 = None

    ke, cp_last = _terminal_log_slope(t, cp)
    auc_tail = cp_last / ke
    auc_0_inf = auc_0_last + auc_tail

    aumc_0_last = _trapezoid(t, t * cp)
    t_last = float(t[-1])
    aumc_tail = cp_last * (t_last / ke + 1.0 / (ke * ke))
    aumc_0_inf = aumc_0_last + aumc_tail

    half_life = math.log(2) / ke
    cl = dose_mg / auc_0_inf
    mrt = aumc_0_inf / auc_0_inf

    if route == "iv_bolus":
        vss = cl * mrt
    else:
        vss = cl * (mrt - infusion_duration_h / 2.0)

    return PKParameters(
        cmax=cmax,
        tmax=tmax,
        auc_0_inf=float(auc_0_inf),
        auc_0_24=None if auc_0_24 is None else float(auc_0_24),
        half_life=float(half_life),
        cl_apparent=float(cl),
        vss=float(vss),
        bioavailability=1.0,  # IV by definition
        fa=1.0,
        fg=1.0,
        fh=None,  # not separately computed for IV
    )
```

- [ ] **Step 4: Run PK-extract tests**

```
python3 -m pytest tests/unit/test_pk_extract.py -v
```

Expected: all pass.

- [ ] **Step 5: Commit**

```
git add src/charon/pbpk/pk_extract.py tests/unit/test_pk_extract.py
git commit -m "Sprint 3: Add pbpk.pk_extract — trapezoidal AUC + Benet Vss"
```

---

## Task 6: pbpk/__init__.py — public exports

Expose the Sprint 3 PBPK API cleanly.

**Files:**
- Modify: `src/charon/pbpk/__init__.py`

- [ ] **Step 1: Write the __init__ module**

Replace `src/charon/pbpk/__init__.py` contents with:

```python
"""Charon Layer 2 — physiologically-based pharmacokinetic simulation.

Public API:

    PBPKTopology, TissueNode, load_species_topology, PORTAL_TISSUES
        Species-level PBPK graph loaded from a YAML file.

    CompoundPBPKParams, build_compound_pbpk_params, build_rhs,
    infer_compound_type
        Translate a :class:`CompoundConfig` into PBPK-ready parameters
        and build the ODE right-hand side closure.

    SimulationResult, simulate_iv
        Run an IV bolus or infusion through the human PBPK kernel using
        scipy's BDF stiff ODE solver.

    compute_pk_parameters
        Extract Cmax, AUC, CL, Vss, t½ from a concentration-time profile.

    compute_kp_rodgers_rowland, compute_all_kp, TissueComposition,
    apply_berezhkovskiy_correction
        Mechanistic Kp calculators (Rodgers & Rowland 2005/2006 plus
        Berezhkovskiy 2004 correction).
"""

from __future__ import annotations

from charon.pbpk.kp_calculator import (
    TissueComposition,
    apply_berezhkovskiy_correction,
    compute_all_kp,
    compute_kp_poulin_theil,
    compute_kp_rodgers_rowland,
)
from charon.pbpk.ode_compiler import (
    CompoundPBPKParams,
    build_compound_pbpk_params,
    build_rhs,
    infer_compound_type,
)
from charon.pbpk.pk_extract import compute_pk_parameters
from charon.pbpk.solver import SimulationResult, simulate_iv
from charon.pbpk.topology import (
    PORTAL_TISSUES,
    PBPKTopology,
    TissueNode,
    load_species_topology,
)

__all__ = [
    "CompoundPBPKParams",
    "PBPKTopology",
    "PORTAL_TISSUES",
    "SimulationResult",
    "TissueComposition",
    "TissueNode",
    "apply_berezhkovskiy_correction",
    "build_compound_pbpk_params",
    "build_rhs",
    "compute_all_kp",
    "compute_kp_poulin_theil",
    "compute_kp_rodgers_rowland",
    "compute_pk_parameters",
    "infer_compound_type",
    "load_species_topology",
    "simulate_iv",
]
```

- [ ] **Step 2: Smoke test imports**

```
python3 -c "from charon.pbpk import simulate_iv, load_species_topology, build_compound_pbpk_params, compute_pk_parameters; print('ok')"
```

Expected: `ok`.

- [ ] **Step 3: Commit**

```
git add src/charon/pbpk/__init__.py
git commit -m "Sprint 3: Publish pbpk public API"
```

---

## Task 7: charon/pipeline.py — top-level Pipeline class

Wires Layer 0 (guardrails — best-effort) + Layer 1 (predict_properties) + ParameterBridge + Layer 2 (PBPK). Provides two constructors: from SMILES (runs ML prediction) and from an explicit CompoundConfig (experimental-override path).

**Files:**
- Create: `src/charon/pipeline.py`
- Create: `tests/unit/test_pipeline.py`

- [ ] **Step 1: Write failing tests**

Create `tests/unit/test_pipeline.py`:

```python
"""Unit tests for the top-level Pipeline wiring."""

import numpy as np
import pytest

from charon.core.schema import (
    BindingProperties,
    CompoundConfig,
    CompoundProperties,
    MetabolismProperties,
    PhysicochemicalProperties,
    PKParameters,
    PredictedProperty,
    RenalProperties,
)
from charon.pipeline import Pipeline, PipelineResult


def _p(value: float, unit: str = "") -> PredictedProperty:
    return PredictedProperty(value=float(value), source="experimental", unit=unit)


@pytest.fixture
def midazolam_compound() -> CompoundConfig:
    return CompoundConfig(
        name="midazolam",
        smiles="Cn1c(Cl)nc2C(=NCc3ccccc13)c4ccc(F)cc4",  # illustrative
        molecular_weight=325.77,
        source="experimental",
        properties=CompoundProperties(
            physicochemical=PhysicochemicalProperties(
                logp=_p(3.89),
                pka_base=_p(6.2),
            ),
            binding=BindingProperties(
                fu_p=_p(0.03, "fraction"),
                fu_inc=_p(0.96, "fraction"),
                bp_ratio=_p(0.66, "ratio"),
            ),
            metabolism=MetabolismProperties(
                clint_uL_min_mg=_p(93.0, "uL/min/mg"),
            ),
            renal=RenalProperties(
                clrenal_L_h=_p(0.0, "L/h"),
            ),
        ),
    )


class TestPipelineFromCompound:
    def test_run_produces_result(self, midazolam_compound):
        pipe = Pipeline(
            compound=midazolam_compound,
            route="iv_bolus",
            dose_mg=5.0,
            duration_h=24.0,
            compound_type_override="base",
        )
        result = pipe.run()
        assert isinstance(result, PipelineResult)
        assert isinstance(result.pk_parameters, PKParameters)
        assert result.pk_parameters.cl_apparent is not None
        assert result.pk_parameters.cl_apparent > 0
        assert result.pk_parameters.vss is not None
        assert result.pk_parameters.vss > 0
        assert len(result.cp_plasma) == len(result.time_h)

    def test_iv_infusion(self, midazolam_compound):
        pipe = Pipeline(
            compound=midazolam_compound,
            route="iv_infusion",
            dose_mg=5.0,
            duration_h=24.0,
            infusion_duration_h=1.0,
            compound_type_override="base",
        )
        result = pipe.run()
        assert result.pk_parameters.cl_apparent is not None

    def test_invalid_route(self, midazolam_compound):
        with pytest.raises(NotImplementedError):
            Pipeline(
                compound=midazolam_compound,
                route="oral",
                dose_mg=5.0,
            ).run()


class TestPipelineMidazolamValidation:
    """End-to-end midazolam IV: CL and Vss within 2-fold of observed.

    Observed (literature, healthy adult IV bolus):
        CL ≈ 21 L/h
        Vss ≈ 66 L
        t½ ≈ 3 h

    Expected PBPK prediction using the input fixture values:
        CLh_predicted ≈ 13.5 L/h  (well-stirred)
    """

    def test_cl_within_2_fold(self, midazolam_compound):
        pipe = Pipeline(
            compound=midazolam_compound,
            route="iv_bolus",
            dose_mg=5.0,
            duration_h=48.0,
            compound_type_override="base",
        )
        result = pipe.run()
        cl = result.pk_parameters.cl_apparent
        assert cl is not None
        observed_cl = 21.0
        fold = max(cl / observed_cl, observed_cl / cl)
        assert fold < 2.0, f"CL fold error {fold:.2f} exceeds 2.0 (pred={cl:.2f}, obs={observed_cl})"

    def test_vss_within_2_fold(self, midazolam_compound):
        pipe = Pipeline(
            compound=midazolam_compound,
            route="iv_bolus",
            dose_mg=5.0,
            duration_h=48.0,
            compound_type_override="base",
        )
        result = pipe.run()
        vss = result.pk_parameters.vss
        assert vss is not None
        observed_vss = 66.0
        fold = max(vss / observed_vss, observed_vss / vss)
        assert fold < 2.0, f"Vss fold error {fold:.2f} exceeds 2.0 (pred={vss:.2f}, obs={observed_vss})"
```

- [ ] **Step 2: Run to verify failure**

```
python3 -m pytest tests/unit/test_pipeline.py -v
```

Expected: ImportError — pipeline module missing.

- [ ] **Step 3: Implement pipeline.py**

Create `src/charon/pipeline.py`:

```python
"""Top-level Charon pipeline: SMILES or CompoundConfig → PK prediction.

Sprint 3 scope: IV bolus / IV infusion only.  Oral (ACAT) routing is
reserved for Sprint 3b.

Two constructors are supported:

1. ``Pipeline(compound=cfg, ...)``
   — use a pre-populated CompoundConfig (experimental overrides).

2. ``Pipeline.from_smiles(smiles, ...)``
   — run Layer 1 ML prediction (ADMET ensemble + pKa + fu_inc + BP) to
     build the compound, then run PBPK.

Both paths converge on the same :meth:`run` method.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal

import numpy as np

from charon.core.parameter_bridge import ParameterBridge
from charon.core.schema import (
    CompoundConfig,
    CompoundProperties,
    PKParameters,
)
from charon.pbpk.ode_compiler import build_compound_pbpk_params
from charon.pbpk.pk_extract import compute_pk_parameters
from charon.pbpk.solver import SimulationResult, simulate_iv
from charon.pbpk.topology import PBPKTopology, load_species_topology


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


class Pipeline:
    """End-to-end Charon Phase A pipeline (Sprint 3 = IV only).

    Parameters
    ----------
    compound : CompoundConfig
        Pre-populated compound (use ``from_smiles`` to build from a SMILES).
    route : "iv_bolus" | "iv_infusion"
        Dose route.  Oral is not yet implemented.
    dose_mg : float
        Dose in mg.
    species : str
        Species YAML name (default: "human").
    duration_h : float
        Simulation duration in hours (default: 72).
    infusion_duration_h : float
        Only used for iv_infusion.
    liver_model : str
        Liver model passed to ParameterBridge (for audit trail).  The PBPK
        rhs consumes CLint_liver_L_h regardless of this choice.
    compound_type_override : str, optional
        Force compound_type in Kp calculation.  Useful for weak bases
        (midazolam pKa_base=6.2) which the Sisyphus threshold misclassifies
        as neutral.
    """

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
    ) -> None:
        self.compound = compound
        self.route = route
        self.dose_mg = dose_mg
        self.species = species
        self.duration_h = duration_h
        self.infusion_duration_h = infusion_duration_h
        self.liver_model = liver_model
        self.compound_type_override = compound_type_override

    @classmethod
    def from_smiles(
        cls,
        smiles: str,
        *,
        route: Literal["iv_bolus", "iv_infusion", "oral"],
        dose_mg: float,
        species: str = "human",
        duration_h: float = 72.0,
        infusion_duration_h: float = 0.0,
        liver_model: str = "well_stirred",
        compound_name: str | None = None,
    ) -> "Pipeline":
        """Build a Pipeline by running Layer 1 ML prediction on a SMILES."""
        from charon.predict import predict_properties
        from charon.core.molecule import Molecule

        mol = Molecule(smiles)
        properties: CompoundProperties = predict_properties(smiles)
        mw = mol.molecular_weight()
        compound = CompoundConfig(
            name=compound_name or smiles,
            smiles=smiles,
            molecular_weight=mw,
            source="predicted",
            properties=properties,
        )
        return cls(
            compound=compound,
            route=route,
            dose_mg=dose_mg,
            species=species,
            duration_h=duration_h,
            infusion_duration_h=infusion_duration_h,
            liver_model=liver_model,
        )

    def run(self) -> PipelineResult:
        """Execute the full pipeline and return a :class:`PipelineResult`."""
        if self.route == "oral":
            raise NotImplementedError(
                "Oral route (ACAT) is deferred to Sprint 3b. "
                "Use 'iv_bolus' or 'iv_infusion' for now."
            )

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
            route=self.route,  # narrowed by early return above
            duration_h=self.duration_h,
            infusion_duration_h=self.infusion_duration_h,
        )

        pk = compute_pk_parameters(
            sim.time_h,
            sim.cp_plasma,
            dose_mg=self.dose_mg,
            route=self.route,
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
            },
        )
```

- [ ] **Step 4: Verify `Molecule.molecular_weight()` exists or add a fallback**

Check the existing `Molecule` class to confirm it exposes `molecular_weight()`; if it does not, use `Chem.Descriptors.MolWt(mol.rdkit_mol)` inside `from_smiles` instead. The authoritative check is:

```
python3 -c "from charon.core.molecule import Molecule; m = Molecule('CCO'); print(m.molecular_weight())"
```

If the call errors, edit `from_smiles` to compute MW directly from RDKit:

```python
from rdkit.Chem import Descriptors
mw = float(Descriptors.MolWt(mol.rdkit_mol))
```

- [ ] **Step 5: Run pipeline tests**

```
python3 -m pytest tests/unit/test_pipeline.py -v
```

Expected: all pass, including the midazolam 2-fold validation.

- [ ] **Step 6: Commit**

```
git add src/charon/pipeline.py tests/unit/test_pipeline.py
git commit -m "Sprint 3: Add charon.Pipeline — top-level IV PBPK orchestration"
```

---

## Task 8: Export Pipeline from charon/__init__.py

**Files:**
- Modify: `src/charon/__init__.py`

- [ ] **Step 1: Write the new __init__**

Replace `src/charon/__init__.py` contents with:

```python
"""Charon — open-source translational PK platform.

Phase A public API entry points:

    Pipeline, PipelineResult
        End-to-end SMILES → PK prediction (Sprint 3: IV only).
"""

from __future__ import annotations

from charon.pipeline import Pipeline, PipelineResult

__all__ = ["Pipeline", "PipelineResult"]
```

- [ ] **Step 2: Smoke test**

```
python3 -c "from charon import Pipeline, PipelineResult; print('ok')"
```

Expected: `ok`.

- [ ] **Step 3: Commit**

```
git add src/charon/__init__.py
git commit -m "Sprint 3: Export Pipeline from charon top-level namespace"
```

---

## Task 9: Integration test — SMILES → PK (tests/integration/test_smiles_to_pk.py)

Fills the empty integration test file with a midazolam end-to-end pass using experimental overrides, matching the pipeline unit test but living in the proper integration folder for regression gating.

**Files:**
- Modify: `tests/integration/test_smiles_to_pk.py`

- [ ] **Step 1: Write the integration test**

Replace `tests/integration/test_smiles_to_pk.py` contents with:

```python
"""Layer 0-1-2 integration: SMILES/compound → Cp-time → PK parameters.

Sprint 3 IV-only subset.  Oral integration tests will land in Sprint 3b.
"""

import numpy as np
import pytest

from charon import Pipeline, PipelineResult
from charon.core.schema import (
    BindingProperties,
    CompoundConfig,
    CompoundProperties,
    MetabolismProperties,
    PhysicochemicalProperties,
    PredictedProperty,
    RenalProperties,
)


def _p(value: float, unit: str = "") -> PredictedProperty:
    return PredictedProperty(value=float(value), source="experimental", unit=unit)


@pytest.fixture
def midazolam() -> CompoundConfig:
    return CompoundConfig(
        name="midazolam",
        smiles="Cc1ncc2n1-c1ccc(Cl)cc1C(c1ccccc1F)=NC2",
        molecular_weight=325.77,
        source="experimental",
        properties=CompoundProperties(
            physicochemical=PhysicochemicalProperties(
                logp=_p(3.89),
                pka_base=_p(6.2),
            ),
            binding=BindingProperties(
                fu_p=_p(0.03, "fraction"),
                fu_inc=_p(0.96, "fraction"),
                bp_ratio=_p(0.66, "ratio"),
            ),
            metabolism=MetabolismProperties(
                clint_uL_min_mg=_p(93.0, "uL/min/mg"),
            ),
            renal=RenalProperties(
                clrenal_L_h=_p(0.0, "L/h"),
            ),
        ),
    )


class TestMidazolamIvBolusE2E:
    def test_runs_and_produces_pk(self, midazolam):
        pipe = Pipeline(
            compound=midazolam,
            route="iv_bolus",
            dose_mg=5.0,
            duration_h=48.0,
            compound_type_override="base",
        )
        result: PipelineResult = pipe.run()
        assert result.pk_parameters.cmax is not None and result.pk_parameters.cmax > 0
        assert result.pk_parameters.cl_apparent is not None and result.pk_parameters.cl_apparent > 0
        assert result.pk_parameters.vss is not None and result.pk_parameters.vss > 0
        assert result.pk_parameters.half_life is not None and result.pk_parameters.half_life > 0

    def test_cl_and_vss_within_2_fold_observed(self, midazolam):
        pipe = Pipeline(
            compound=midazolam,
            route="iv_bolus",
            dose_mg=5.0,
            duration_h=48.0,
            compound_type_override="base",
        )
        result = pipe.run()
        cl, vss = result.pk_parameters.cl_apparent, result.pk_parameters.vss
        observed = {"cl": 21.0, "vss": 66.0}
        cl_fold = max(cl / observed["cl"], observed["cl"] / cl)
        vss_fold = max(vss / observed["vss"], observed["vss"] / vss)
        assert cl_fold < 2.0, f"CL fold error {cl_fold:.2f} (pred={cl:.2f}, obs={observed['cl']})"
        assert vss_fold < 2.0, f"Vss fold error {vss_fold:.2f} (pred={vss:.2f}, obs={observed['vss']})"

    def test_metadata_includes_pbpk_specifics(self, midazolam):
        pipe = Pipeline(
            compound=midazolam,
            route="iv_bolus",
            dose_mg=5.0,
            duration_h=24.0,
            compound_type_override="base",
        )
        result = pipe.run()
        meta = result.metadata
        assert meta["solver_method"] == "BDF"
        assert meta["compound_type"] == "base"
        assert meta["fu_b"] == pytest.approx(0.03 / 0.66, rel=1e-6)
        assert meta["clint_liver_L_h"] == pytest.approx(348.75, rel=1e-3)
```

- [ ] **Step 2: Run the integration test**

```
python3 -m pytest tests/integration/test_smiles_to_pk.py -v
```

Expected: all pass.

- [ ] **Step 3: Commit**

```
git add tests/integration/test_smiles_to_pk.py
git commit -m "Sprint 3: Add Layer 0-2 integration test for midazolam IV"
```

---

## Task 10: validation/benchmarks/metrics.py — AAFE / fold error / within-N-fold

Reusable benchmark utility functions. Tested in unit tests; consumed by layer2_human_pk.py.

**Files:**
- Create: `validation/benchmarks/metrics.py`
- Create: `tests/unit/test_benchmark_metrics.py`

- [ ] **Step 1: Write failing tests**

Create `tests/unit/test_benchmark_metrics.py`:

```python
"""Unit tests for validation/benchmarks/metrics.py."""

import math
import sys
from pathlib import Path

# validation/ is outside the package; import via path
REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from validation.benchmarks.metrics import (  # noqa: E402
    aafe,
    fold_error,
    within_n_fold,
)


class TestFoldError:
    def test_equal_values(self):
        assert fold_error(10.0, 10.0) == pytest.approx(1.0)

    def test_2_fold_over(self):
        assert fold_error(20.0, 10.0) == pytest.approx(2.0)

    def test_2_fold_under(self):
        assert fold_error(5.0, 10.0) == pytest.approx(2.0)

    def test_zero_predicted_raises(self):
        with pytest.raises(ValueError):
            fold_error(0.0, 10.0)

    def test_zero_observed_raises(self):
        with pytest.raises(ValueError):
            fold_error(10.0, 0.0)


class TestAAFE:
    def test_identity(self):
        predicted = [1.0, 2.0, 3.0]
        observed = [1.0, 2.0, 3.0]
        assert aafe(predicted, observed) == pytest.approx(1.0)

    def test_uniform_2_fold(self):
        predicted = [2.0, 4.0, 8.0]
        observed = [1.0, 2.0, 4.0]
        assert aafe(predicted, observed) == pytest.approx(2.0)

    def test_mixed(self):
        predicted = [2.0, 1.0]  # 2-fold over, 2-fold under
        observed = [1.0, 2.0]
        # geometric mean of (2, 2) = 2
        assert aafe(predicted, observed) == pytest.approx(2.0)


class TestWithinNFold:
    def test_all_within(self):
        predicted = [1.1, 1.5, 0.8]
        observed = [1.0, 1.0, 1.0]
        assert within_n_fold(predicted, observed, n=2.0) == pytest.approx(1.0)

    def test_half_within(self):
        predicted = [1.1, 3.5]   # 1.1-fold, 3.5-fold
        observed = [1.0, 1.0]
        assert within_n_fold(predicted, observed, n=2.0) == pytest.approx(0.5)


# Importing pytest late so the sys.path hack can run without pytest sensing.
import pytest  # noqa: E402
```

- [ ] **Step 2: Run to verify failure**

```
python3 -m pytest tests/unit/test_benchmark_metrics.py -v
```

Expected: ImportError.

- [ ] **Step 3: Implement metrics.py**

Create `validation/benchmarks/metrics.py`:

```python
"""Benchmark metrics for Layer 1-3 validation suites.

Metrics
-------
fold_error : max(pred/obs, obs/pred)
aafe       : geometric mean of fold errors (absolute average fold error)
within_n_fold : fraction of predictions with fold_error <= n
"""

from __future__ import annotations

import math
from collections.abc import Iterable


def fold_error(predicted: float, observed: float) -> float:
    """Two-sided fold error. Always >= 1.0."""
    if predicted <= 0 or observed <= 0:
        raise ValueError(
            f"fold_error requires positive values, got "
            f"predicted={predicted}, observed={observed}"
        )
    return max(predicted / observed, observed / predicted)


def aafe(predicted: Iterable[float], observed: Iterable[float]) -> float:
    """Absolute Average Fold Error (geometric mean of fold errors).

    Definition: AAFE = 10^( mean(|log10(pred/obs)|) ).
    """
    preds = list(predicted)
    obs = list(observed)
    if len(preds) != len(obs):
        raise ValueError(
            f"aafe: length mismatch ({len(preds)} vs {len(obs)})"
        )
    if not preds:
        raise ValueError("aafe: empty input")
    log_fold_errors = []
    for p, o in zip(preds, obs):
        if p <= 0 or o <= 0:
            raise ValueError(
                f"aafe requires positive values, got pred={p}, obs={o}"
            )
        log_fold_errors.append(abs(math.log10(p / o)))
    mean_log_fold = sum(log_fold_errors) / len(log_fold_errors)
    return 10 ** mean_log_fold


def within_n_fold(
    predicted: Iterable[float],
    observed: Iterable[float],
    *,
    n: float = 2.0,
) -> float:
    """Fraction of predictions with fold_error <= n (range [0, 1])."""
    preds = list(predicted)
    obs = list(observed)
    if len(preds) != len(obs):
        raise ValueError(
            f"within_n_fold: length mismatch ({len(preds)} vs {len(obs)})"
        )
    if not preds:
        raise ValueError("within_n_fold: empty input")
    hits = 0
    for p, o in zip(preds, obs):
        if fold_error(p, o) <= n:
            hits += 1
    return hits / len(preds)
```

- [ ] **Step 4: Run metric tests**

```
python3 -m pytest tests/unit/test_benchmark_metrics.py -v
```

Expected: all pass.

- [ ] **Step 5: Commit**

```
git add validation/benchmarks/metrics.py tests/unit/test_benchmark_metrics.py
git commit -m "Sprint 3: Add validation benchmark metrics (AAFE, fold_error, within_n_fold)"
```

---

## Task 11: validation/benchmarks/layer2_human_pk.py — midazolam benchmark script

Standalone runnable script that builds the midazolam fixture, runs the Pipeline, and prints a PASS/FAIL table against the CL and Vss targets.

**Files:**
- Create: `validation/benchmarks/layer2_human_pk.py`

- [ ] **Step 1: Write the benchmark script**

Create `validation/benchmarks/layer2_human_pk.py`:

```python
"""Layer 2 human PBPK benchmark — midazolam IV.

Run as a standalone script::

    python3 validation/benchmarks/layer2_human_pk.py

Prints a table comparing predicted CL / Vss to literature observed values
and flags PASS/FAIL against the Sprint 3 2-fold target.

Sprint 3 scope: 1 compound (midazolam).  Sprint 3b will expand to the full
Obach 1999 IV dataset.
"""

from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from charon import Pipeline  # noqa: E402
from charon.core.schema import (  # noqa: E402
    BindingProperties,
    CompoundConfig,
    CompoundProperties,
    MetabolismProperties,
    PhysicochemicalProperties,
    PredictedProperty,
    RenalProperties,
)
from validation.benchmarks.metrics import fold_error  # noqa: E402


def _p(value: float, unit: str = "") -> PredictedProperty:
    return PredictedProperty(value=float(value), source="experimental", unit=unit)


def midazolam_compound() -> CompoundConfig:
    """Experimental-override midazolam for Layer 2 benchmarking."""
    return CompoundConfig(
        name="midazolam",
        smiles="Cc1ncc2n1-c1ccc(Cl)cc1C(c1ccccc1F)=NC2",
        molecular_weight=325.77,
        source="experimental",
        properties=CompoundProperties(
            physicochemical=PhysicochemicalProperties(
                logp=_p(3.89),
                pka_base=_p(6.2),
            ),
            binding=BindingProperties(
                fu_p=_p(0.03, "fraction"),
                fu_inc=_p(0.96, "fraction"),
                bp_ratio=_p(0.66, "ratio"),
            ),
            metabolism=MetabolismProperties(
                clint_uL_min_mg=_p(93.0, "uL/min/mg"),
            ),
            renal=RenalProperties(clrenal_L_h=_p(0.0, "L/h")),
        ),
    )


def main() -> int:
    observed = {"cl": 21.0, "vss": 66.0, "t_half": 3.0}
    compound = midazolam_compound()
    result = Pipeline(
        compound=compound,
        route="iv_bolus",
        dose_mg=5.0,
        duration_h=48.0,
        compound_type_override="base",
    ).run()

    pk = result.pk_parameters
    pred = {
        "cl": pk.cl_apparent,
        "vss": pk.vss,
        "t_half": pk.half_life,
    }

    print("=" * 72)
    print("Layer 2 Human PBPK Benchmark — Midazolam IV bolus 5 mg")
    print("=" * 72)
    print(f"{'Metric':<12} {'Predicted':>12} {'Observed':>12} {'Fold Err':>10} {'Verdict':>10}")
    print("-" * 72)

    target_fold = 2.0
    all_pass = True
    for metric in ("cl", "vss", "t_half"):
        p = pred[metric]
        o = observed[metric]
        fe = fold_error(p, o)
        verdict = "PASS" if fe <= target_fold else "FAIL"
        if verdict == "FAIL":
            all_pass = False
        print(f"{metric:<12} {p:>12.3f} {o:>12.3f} {fe:>10.3f} {verdict:>10}")

    print("-" * 72)
    print(f"PBPK params used:")
    print(f"  compound_type         = {result.metadata['compound_type']}")
    print(f"  fu_b                  = {result.metadata['fu_b']:.6f}")
    print(f"  CLint_liver_L_h       = {result.metadata['clint_liver_L_h']:.3f}")
    print(f"  Solver method         = {result.metadata['solver_method']}")
    print(f"  Solver nfev           = {result.metadata['solver_nfev']}")
    print("=" * 72)

    return 0 if all_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
```

- [ ] **Step 2: Run the benchmark script**

```
python3 validation/benchmarks/layer2_human_pk.py
```

Expected: table printed, exit code 0, all metrics PASS.

- [ ] **Step 3: Commit**

```
git add validation/benchmarks/layer2_human_pk.py
git commit -m "Sprint 3: Add runnable layer2_human_pk benchmark (midazolam IV)"
```

---

## Task 12: Full test suite + coverage gate + final commit

- [ ] **Step 1: Run the complete test suite**

```
python3 -m pytest tests/ -v --tb=short
```

Expected: All 390 pre-existing tests + all new tests pass. Zero failures, zero errors.

- [ ] **Step 2: Coverage check on pbpk and pipeline**

```
python3 -m pytest tests/unit/test_topology.py tests/unit/test_ode_compiler.py tests/unit/test_solver.py tests/unit/test_pk_extract.py tests/unit/test_pipeline.py tests/unit/test_benchmark_metrics.py tests/integration/test_smiles_to_pk.py --cov=src/charon/pbpk --cov=src/charon/pipeline --cov-report=term-missing
```

Expected: `src/charon/pbpk` ≥ 80% coverage (excluding empty formulation/ stubs), `src/charon/pipeline.py` ≥ 80% coverage.

- [ ] **Step 3: Run the benchmark script one more time as final smoke**

```
python3 validation/benchmarks/layer2_human_pk.py
```

Expected: PASS on all rows, exit 0.

- [ ] **Step 4: Final summary commit (if anything uncommitted)**

```
git status
git log --oneline -10
```

Expected: clean tree, Sprint 3 commits visible.

---

## Self-review checklist (applied during plan writing)

1. **Spec coverage** — Every spec section mapped to a task:
   - §2 Validation target → Task 9 + Task 11
   - §3.1 Topology → Task 2
   - §3.2 Mass balance equations → Task 3
   - §3.3 Compound → PBPK params → Task 3
   - §3.4 Solver → Task 4
   - §3.5 PK extraction → Task 5
   - §4 File layout → Tasks 1-11
   - §5 Data types → Tasks 2-7
   - §6 Risks → pKa override (Task 3), Kp cap (documented), BDF enforcement (Task 4), kidney BP scaling (Task 3)

2. **Placeholder scan** — No "TBD", "TODO", "similar to", or code-less steps. All code shown in full.

3. **Type consistency** — `CompoundPBPKParams`, `PBPKTopology`, `TissueNode`, `SimulationResult`, `PipelineResult` field names are referenced identically across Tasks 2-9. `build_rhs` signature matches between ode_compiler.py, solver.py, and the call site. `load_species_topology` signature matches between topology.py and solver test fixture. `compute_pk_parameters` signature matches between pk_extract.py and pipeline.py.

4. **Gap check** — All empty files listed in design §4 are created or explicitly deferred (formulation/ stays empty per non-goals).
