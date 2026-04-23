"""predict_properties branches correctly between Tier 2 and Tier 3."""
from __future__ import annotations

import pytest

from charon.predict import predict_properties
from charon.predict.conformal import reset_default_conformal


CAFFEINE = "Cn1cnc2n(C)c(=O)n(C)c(=O)c12"
THEOPHYLLINE = "Cn1c(=O)c2[nH]cnc2n(C)c1=O"
PROPRANOLOL = "CC(C)NCC(O)COc1cccc2ccccc12"


@pytest.fixture(autouse=True)
def _reset():
    reset_default_conformal()
    yield
    reset_default_conformal()


class TestADLowTriggersTier3:
    def test_theophylline_is_tier3(self):
        props = predict_properties(THEOPHYLLINE)
        clint = props.metabolism.clint_uL_min_mg
        assert clint is not None
        assert clint.source == "classification"
        assert clint.classifier_probs is not None
        assert set(clint.classifier_probs) == {"low", "med", "high"}
        assert clint.flag is not None and "CRITICAL" in clint.flag
        # Point estimate within current bucket range (theophylline -> Low most likely)
        bucket_lo = clint.ci_90_lower
        bucket_hi = clint.ci_90_upper
        assert bucket_lo is not None and bucket_hi is not None
        assert bucket_lo <= clint.value <= bucket_hi


class TestADHighStaysTier2:
    def test_propranolol_is_tier2(self):
        props = predict_properties(PROPRANOLOL)
        clint = props.metabolism.clint_uL_min_mg
        assert clint.source == "ml_ensemble"
        assert clint.classifier_probs is None


class TestADModerateIsTier2WithCaution:
    def test_caffeine_is_tier2_with_caution_flag(self):
        props = predict_properties(CAFFEINE)
        clint = props.metabolism.clint_uL_min_mg
        assert clint.source == "ml_ensemble"
        assert clint.classifier_probs is None
        assert clint.flag is not None
        assert "ad_moderate_caution" in clint.flag


class TestForceTier3:
    def test_force_tier3_on_in_domain(self):
        props_default = predict_properties(PROPRANOLOL)
        props_forced = predict_properties(PROPRANOLOL, force_tier3=True)
        assert props_default.metabolism.clint_uL_min_mg.source == "ml_ensemble"
        assert props_forced.metabolism.clint_uL_min_mg.source == "classification"
