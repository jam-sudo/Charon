"""Tests for the Kp empirical override path in ode_compiler.

These tests cover: the KpOverrideRecord dataclass, the species guard,
compound_type precedence, the empirical override loop, and the
end-to-end interaction with CompoundPBPKParams.
"""

from __future__ import annotations

from collections import OrderedDict

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
from charon.pbpk.kp_calculator import TissueComposition
from charon.pbpk.ode_compiler import (
    CompoundPBPKParams,
    KpOverrideRecord,
    build_compound_pbpk_params,
)
from charon.pbpk.topology import PBPKTopology, TissueNode, load_species_topology


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
