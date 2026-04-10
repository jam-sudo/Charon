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
