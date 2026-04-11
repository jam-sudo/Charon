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


class TestKpMethodSelection:
    def setup_method(self):
        from charon.core.parameter_bridge import ParameterBridge
        from charon.pbpk.topology import load_species_topology
        self.topo = load_species_topology("human")
        self.bridge = ParameterBridge()

    def _warfarin_like(self, kp_method: str | None = None):
        from charon.core.schema import (
            BindingProperties, CompoundConfig, CompoundProperties,
            MetabolismProperties, PhysicochemicalProperties,
            PredictedProperty, RenalProperties,
        )
        pc_kwargs = dict(
            logp=PredictedProperty(value=3.54, source="experimental"),
            pka_acid=PredictedProperty(value=5.0, source="experimental"),
            compound_type="acid",
        )
        if kp_method is not None:
            pc_kwargs["kp_method"] = kp_method
        return CompoundConfig(
            name="warfarin_like",
            smiles="CC(=O)CC(c1ccccc1)C1=C(O)c2ccccc2OC1=O",
            molecular_weight=308.33,
            source="experimental",
            properties=CompoundProperties(
                physicochemical=PhysicochemicalProperties(**pc_kwargs),
                binding=BindingProperties(
                    fu_p=PredictedProperty(value=0.012, source="experimental", unit="fraction"),
                    fu_inc=PredictedProperty(value=0.17, source="experimental", unit="fraction"),
                    bp_ratio=PredictedProperty(value=0.56, source="experimental", unit="ratio"),
                ),
                metabolism=MetabolismProperties(
                    clint_uL_min_mg=PredictedProperty(value=3.1, source="experimental", unit="uL/min/mg"),
                ),
                renal=RenalProperties(
                    clrenal_L_h=PredictedProperty(value=0.05, source="experimental", unit="L/h"),
                ),
            ),
        )

    def test_default_is_rodgers_rowland(self):
        """No kp_method field → uses R&R.

        Hand-calculation:
          logP=3.54, pKa_acid=5.0, compound_type=acid
          R&R adipose Kp ≈ 22.55 (direct compute_all_kp verification)
        """
        from charon.pbpk.ode_compiler import build_compound_pbpk_params
        compound = self._warfarin_like(kp_method=None)
        params = build_compound_pbpk_params(compound, self.topo, bridge=self.bridge)
        # R&R warfarin adipose Kp ≈ 22.55 (experimentally verified via direct
        # compute_all_kp call)
        assert params.kp_by_tissue["adipose"] == pytest.approx(22.55, rel=0.01)

    def test_explicit_rodgers_rowland(self):
        """Explicit kp_method='rodgers_rowland' → same result as default."""
        from charon.pbpk.ode_compiler import build_compound_pbpk_params
        compound = self._warfarin_like(kp_method="rodgers_rowland")
        params = build_compound_pbpk_params(compound, self.topo, bridge=self.bridge)
        assert params.kp_by_tissue["adipose"] == pytest.approx(22.55, rel=0.01)

    def test_berezhkovskiy_reduces_kp(self):
        """BZ method with low fu_p should give lower Kp than R&R.

        BZ formula: Kp_bz = Kp_rr / (1 + (Kp_rr - 1) * fu_p)
        With fu_p=0.012, Kp_rr=22.55:
          denom = 1 + (22.55 - 1) * 0.012 = 1 + 0.2586 = 1.2586
          Kp_bz = 22.55 / 1.2586 ≈ 17.92
        """
        from charon.pbpk.ode_compiler import build_compound_pbpk_params
        compound_rr = self._warfarin_like(kp_method="rodgers_rowland")
        compound_bz = self._warfarin_like(kp_method="berezhkovskiy")
        params_rr = build_compound_pbpk_params(compound_rr, self.topo, bridge=self.bridge)
        params_bz = build_compound_pbpk_params(compound_bz, self.topo, bridge=self.bridge)
        # BZ should REDUCE adipose Kp (correction for highly bound drugs)
        assert params_bz.kp_by_tissue["adipose"] < params_rr.kp_by_tissue["adipose"]

    def test_berezhkovskiy_needs_fu_p(self):
        """BZ method requires fu_p to be passed to compute_all_kp.
        build_compound_pbpk_params should handle this automatically.
        """
        from charon.pbpk.ode_compiler import build_compound_pbpk_params
        compound = self._warfarin_like(kp_method="berezhkovskiy")
        # Should NOT raise — fu_p is available via compound.properties.binding.fu_p
        params = build_compound_pbpk_params(compound, self.topo, bridge=self.bridge)
        assert params.compound_type == "acid"
        assert params.kp_by_tissue["adipose"] > 0
