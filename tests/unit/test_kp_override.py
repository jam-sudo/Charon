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

        # Override reduces adipose Kp (50 -> 10), shrinking Vd.
        # Less drug sequestered in adipose -> higher plasma exposure.
        # AUC_override > AUC_baseline demonstrates the override reached the ODE.
        import numpy as np
        auc_no = np.trapz(sim_no.cp_plasma, sim_no.time_h)
        auc_yes = np.trapz(sim_yes.cp_plasma, sim_yes.time_h)
        assert auc_yes > auc_no
