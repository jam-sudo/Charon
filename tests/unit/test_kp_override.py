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
