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
