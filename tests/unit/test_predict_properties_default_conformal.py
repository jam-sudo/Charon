"""Default conformal-attach behaviour in predict_properties."""
from __future__ import annotations

import pytest

from charon.predict import predict_properties
from charon.predict.conformal import reset_default_conformal


CAFFEINE = "Cn1cnc2n(C)c(=O)n(C)c(=O)c12"


@pytest.fixture(autouse=True)
def _reset():
    reset_default_conformal()
    yield
    reset_default_conformal()


class TestDefaultAttaches:
    def test_fup_has_ci_by_default(self):
        props = predict_properties(CAFFEINE)
        fup = props.binding.fu_p
        assert fup is not None
        assert fup.ci_90_lower is not None
        assert fup.ci_90_upper is not None
        assert fup.ci_90_lower < fup.value <= fup.ci_90_upper + 1e-9
        assert fup.ci_90_upper <= 1.0

    def test_clint_has_ci_by_default(self):
        props = predict_properties(CAFFEINE)
        clint = props.metabolism.clint_uL_min_mg
        assert clint is not None
        assert clint.ci_90_lower is not None
        assert clint.ci_90_upper is not None
        assert clint.ci_90_lower < clint.value <= clint.ci_90_upper + 1e-9


class TestExplicitOptOut:
    def test_explicit_none_sentinel_disables_ci(self):
        from charon.predict import CONFORMAL_OFF
        props = predict_properties(CAFFEINE, conformal=CONFORMAL_OFF)
        assert props.binding.fu_p.ci_90_lower is None
        assert props.metabolism.clint_uL_min_mg.ci_90_lower is None


class TestCustomPredictor:
    def test_user_supplied_predictor_used(self):
        from charon.predict.conformal import ConformalPredictor
        custom = ConformalPredictor.load_default()
        props = predict_properties(CAFFEINE, conformal=custom)
        assert props.binding.fu_p.ci_90_lower is not None
