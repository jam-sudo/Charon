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
