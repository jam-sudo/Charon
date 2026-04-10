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
        smiles="Cc1ncc2n1-c1ccc(Cl)cc1C(c1ccccc1F)=NC2",
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
