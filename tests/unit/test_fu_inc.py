"""Unit tests for the Austin 2002 microsomal fu_inc correlation."""

from __future__ import annotations

import math

import pytest

from charon.predict.fu_inc import predict_fu_inc


class TestPredictFuInc:
    """Austin 2002 fu_inc correlation at 1 mg/mL HLM."""

    def test_polar_compound_near_unity(self):
        """A very polar compound (logP = -2) should have fu_inc > 0.9."""
        fu_inc = predict_fu_inc(logp=-2.0)
        assert fu_inc > 0.9
        assert fu_inc <= 1.0

    def test_moderate_lipophilicity(self):
        """logP = 2 should give fu_inc in (0.7, 0.95)."""
        fu_inc = predict_fu_inc(logp=2.0)
        assert 0.7 < fu_inc < 0.95

    def test_highly_lipophilic(self):
        """logP = 6 should give fu_inc < 0.05 (heavy microsomal binding)."""
        fu_inc = predict_fu_inc(logp=6.0)
        assert fu_inc < 0.05
        assert fu_inc > 0.0

    def test_monotonic_decreasing(self):
        """fu_inc should decrease monotonically with increasing logP."""
        f1 = predict_fu_inc(logp=1.0)
        f3 = predict_fu_inc(logp=3.0)
        f5 = predict_fu_inc(logp=5.0)
        assert f1 > f3 > f5

    @pytest.mark.parametrize(
        "logp",
        [-3.0, -1.0, 0.0, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 7.0],
    )
    def test_always_in_unit_interval(self, logp):
        """fu_inc output must always lie in (0, 1]."""
        fu_inc = predict_fu_inc(logp=logp)
        assert 0.0 < fu_inc <= 1.0

    def test_exact_formula_logp_3(self):
        """Hand calculation for logP=3."""
        logp = 3.0
        log_ratio = 0.072 * (logp ** 2) + 0.067 * logp - 1.126
        expected = 1.0 / (1.0 + 10.0 ** log_ratio)
        assert predict_fu_inc(logp=logp) == pytest.approx(expected, rel=1e-12)

    def test_nan_raises(self):
        with pytest.raises(ValueError, match="finite"):
            predict_fu_inc(logp=float("nan"))

    def test_posinf_raises(self):
        with pytest.raises(ValueError, match="finite"):
            predict_fu_inc(logp=float("inf"))

    def test_neginf_raises(self):
        with pytest.raises(ValueError, match="finite"):
            predict_fu_inc(logp=float("-inf"))

    def test_logp_zero_matches_formula(self):
        """At logP=0 the correlation reduces to 1 / (1 + 10^-1.126)."""
        expected = 1.0 / (1.0 + 10.0 ** (-1.126))
        assert predict_fu_inc(logp=0.0) == pytest.approx(expected, rel=1e-12)
