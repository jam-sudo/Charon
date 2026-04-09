"""Unit tests for charon.core.liver_models hepatic clearance models."""

from __future__ import annotations

import math

import pytest

from charon.core.liver_models import dispersion, get_liver_model, parallel_tube, well_stirred


# ---------------------------------------------------------------------------
# Well-stirred model
# ---------------------------------------------------------------------------


class TestWellStirred:
    """Tests for the well-stirred liver model."""

    def test_basic_calculation(self):
        """
        Hand calculation: Qh=90, fu_b=0.5, CLint=100
        CLh = (90 * 0.5 * 100) / (90 + 0.5 * 100) = 4500 / 140 = 32.142857...
        """
        result = well_stirred(qh=90.0, fu_b=0.5, clint_liver=100.0)
        assert result == pytest.approx(4500.0 / 140.0, rel=1e-6)

    def test_zero_clint(self):
        """CLint=0 -> CLh=0."""
        result = well_stirred(qh=90.0, fu_b=0.5, clint_liver=0.0)
        assert result == 0.0

    def test_flow_limited(self):
        """Very high CLint -> CLh approaches Qh."""
        result = well_stirred(qh=90.0, fu_b=1.0, clint_liver=1e8)
        assert result == pytest.approx(90.0, rel=1e-3)

    def test_low_clint_convergence(self):
        """Low CLint: CLh approximately equals fu_b * CLint."""
        qh = 90.0
        fu_b = 0.5
        clint = 0.1  # very low
        result = well_stirred(qh=qh, fu_b=fu_b, clint_liver=clint)
        approx_result = fu_b * clint
        assert result == pytest.approx(approx_result, rel=0.01)

    def test_known_values(self):
        """
        Qh=90, fu_b=0.03, CLint=50
        CLh = (90 * 0.03 * 50) / (90 + 0.03 * 50) = 135 / 91.5 = 1.4754...
        """
        result = well_stirred(qh=90.0, fu_b=0.03, clint_liver=50.0)
        assert result == pytest.approx(135.0 / 91.5, rel=1e-6)

    def test_fu_b_one(self):
        """fu_b=1.0 (fully unbound): CLh = Qh*CLint/(Qh+CLint)."""
        result = well_stirred(qh=90.0, fu_b=1.0, clint_liver=90.0)
        assert result == pytest.approx(90.0 * 90.0 / (90.0 + 90.0), rel=1e-6)

    def test_fu_b_zero(self):
        """fu_b=0.0 (fully bound): CLh = 0."""
        result = well_stirred(qh=90.0, fu_b=0.0, clint_liver=100.0)
        assert result == pytest.approx(0.0)


# ---------------------------------------------------------------------------
# Parallel-tube model
# ---------------------------------------------------------------------------


class TestParallelTube:
    """Tests for the parallel-tube liver model."""

    def test_basic_calculation(self):
        """
        Hand calc: Qh=90, fu_b=0.5, CLint=100
        CLh = 90 * (1 - exp(-0.5*100/90)) = 90 * (1 - exp(-0.5556))
        exp(-0.5556) = 0.57375...
        CLh = 90 * (1 - 0.57375) = 90 * 0.42625 = 38.362...
        """
        result = parallel_tube(qh=90.0, fu_b=0.5, clint_liver=100.0)
        expected = 90.0 * (1.0 - math.exp(-0.5 * 100.0 / 90.0))
        assert result == pytest.approx(expected, rel=1e-6)

    def test_zero_clint(self):
        """CLint=0 -> CLh=0."""
        result = parallel_tube(qh=90.0, fu_b=0.5, clint_liver=0.0)
        assert result == 0.0

    def test_flow_limited(self):
        """Very high CLint -> CLh approaches Qh."""
        result = parallel_tube(qh=90.0, fu_b=1.0, clint_liver=1e8)
        assert result == pytest.approx(90.0, rel=1e-3)

    def test_pt_ge_ws(self):
        """Parallel-tube always >= well-stirred for same inputs."""
        test_cases = [
            (90.0, 0.5, 100.0),
            (90.0, 0.03, 50.0),
            (90.0, 1.0, 200.0),
            (90.0, 0.1, 10.0),
        ]
        for qh, fu_b, clint in test_cases:
            pt_val = parallel_tube(qh=qh, fu_b=fu_b, clint_liver=clint)
            ws_val = well_stirred(qh=qh, fu_b=fu_b, clint_liver=clint)
            assert pt_val >= ws_val - 1e-12, (
                f"PT ({pt_val}) should be >= WS ({ws_val}) for "
                f"qh={qh}, fu_b={fu_b}, clint={clint}"
            )

    def test_convergence_low_clint(self):
        """At low CLint, parallel-tube approximately equals well-stirred."""
        qh = 90.0
        fu_b = 0.5
        clint = 0.01  # very low
        pt_val = parallel_tube(qh=qh, fu_b=fu_b, clint_liver=clint)
        ws_val = well_stirred(qh=qh, fu_b=fu_b, clint_liver=clint)
        assert pt_val == pytest.approx(ws_val, rel=0.01)

    def test_fu_b_zero(self):
        """fu_b=0: exp(-0) = 1, so CLh = 0."""
        result = parallel_tube(qh=90.0, fu_b=0.0, clint_liver=100.0)
        assert result == pytest.approx(0.0)


# ---------------------------------------------------------------------------
# Dispersion model
# ---------------------------------------------------------------------------


class TestDispersion:
    """Tests for the dispersion liver model."""

    def test_dn_zero_approaches_pt(self):
        """DN very small -> approaches parallel-tube result."""
        qh, fu_b, clint = 90.0, 0.5, 100.0
        disp_val = dispersion(qh=qh, fu_b=fu_b, clint_liver=clint, dn=0.001)
        pt_val = parallel_tube(qh=qh, fu_b=fu_b, clint_liver=clint)
        assert disp_val == pytest.approx(pt_val, rel=0.05)

    def test_dn_large_approaches_ws(self):
        """DN very large -> approaches well-stirred result."""
        qh, fu_b, clint = 90.0, 0.5, 100.0
        disp_val = dispersion(qh=qh, fu_b=fu_b, clint_liver=clint, dn=100.0)
        ws_val = well_stirred(qh=qh, fu_b=fu_b, clint_liver=clint)
        assert disp_val == pytest.approx(ws_val, rel=0.05)

    def test_default_dn(self):
        """Default DN=0.17 gives result between WS and PT."""
        qh, fu_b, clint = 90.0, 0.5, 100.0
        disp_val = dispersion(qh=qh, fu_b=fu_b, clint_liver=clint)
        ws_val = well_stirred(qh=qh, fu_b=fu_b, clint_liver=clint)
        pt_val = parallel_tube(qh=qh, fu_b=fu_b, clint_liver=clint)
        assert ws_val <= disp_val <= pt_val + 1e-9

    def test_zero_clint(self):
        """CLint=0 -> CLh=0."""
        result = dispersion(qh=90.0, fu_b=0.5, clint_liver=0.0)
        assert result == 0.0

    def test_flow_limited(self):
        """Very high CLint -> CLh approaches Qh."""
        result = dispersion(qh=90.0, fu_b=1.0, clint_liver=1e6, dn=0.17)
        assert result == pytest.approx(90.0, rel=0.01)

    def test_bounded_by_qh(self):
        """CLh should never exceed Qh."""
        result = dispersion(qh=90.0, fu_b=1.0, clint_liver=1e8, dn=0.17)
        assert result <= 90.0 + 1e-9


# ---------------------------------------------------------------------------
# get_liver_model factory
# ---------------------------------------------------------------------------


class TestGetLiverModel:
    """Tests for the get_liver_model factory function."""

    def test_valid_names(self):
        for name in ["well_stirred", "parallel_tube", "dispersion"]:
            fn = get_liver_model(name)
            assert callable(fn)

    def test_well_stirred_returns_correct_function(self):
        fn = get_liver_model("well_stirred")
        assert fn is well_stirred

    def test_parallel_tube_returns_correct_function(self):
        fn = get_liver_model("parallel_tube")
        assert fn is parallel_tube

    def test_dispersion_returns_correct_function(self):
        fn = get_liver_model("dispersion")
        assert fn is dispersion

    def test_invalid_name_raises(self):
        with pytest.raises(ValueError, match="Unknown liver model"):
            get_liver_model("invalid_model")

    def test_empty_string_raises(self):
        with pytest.raises(ValueError, match="Unknown liver model"):
            get_liver_model("")

    def test_case_sensitive(self):
        with pytest.raises(ValueError, match="Unknown liver model"):
            get_liver_model("Well_Stirred")


# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------


class TestInputValidation:
    """Tests for input validation across all liver models."""

    @pytest.mark.parametrize("model_fn", [well_stirred, parallel_tube, dispersion])
    def test_negative_qh_raises(self, model_fn):
        with pytest.raises(ValueError, match="qh must be > 0"):
            model_fn(qh=-10.0, fu_b=0.5, clint_liver=100.0)

    @pytest.mark.parametrize("model_fn", [well_stirred, parallel_tube, dispersion])
    def test_zero_qh_raises(self, model_fn):
        with pytest.raises(ValueError, match="qh must be > 0"):
            model_fn(qh=0.0, fu_b=0.5, clint_liver=100.0)

    @pytest.mark.parametrize("model_fn", [well_stirred, parallel_tube, dispersion])
    def test_negative_clint_raises(self, model_fn):
        with pytest.raises(ValueError, match="clint_liver must be >= 0"):
            model_fn(qh=90.0, fu_b=0.5, clint_liver=-10.0)

    @pytest.mark.parametrize("model_fn", [well_stirred, parallel_tube, dispersion])
    def test_fu_b_below_zero_raises(self, model_fn):
        with pytest.raises(ValueError, match="fu_b must be in"):
            model_fn(qh=90.0, fu_b=-0.1, clint_liver=100.0)

    @pytest.mark.parametrize("model_fn", [well_stirred, parallel_tube, dispersion])
    def test_fu_b_above_one_raises(self, model_fn):
        with pytest.raises(ValueError, match="fu_b must be in"):
            model_fn(qh=90.0, fu_b=1.1, clint_liver=100.0)

    @pytest.mark.parametrize("model_fn", [well_stirred, parallel_tube, dispersion])
    def test_fu_b_boundary_zero_allowed(self, model_fn):
        """fu_b=0 is allowed (fully bound drug)."""
        result = model_fn(qh=90.0, fu_b=0.0, clint_liver=100.0)
        assert result >= 0.0

    @pytest.mark.parametrize("model_fn", [well_stirred, parallel_tube, dispersion])
    def test_fu_b_boundary_one_allowed(self, model_fn):
        """fu_b=1 is allowed (fully unbound drug)."""
        result = model_fn(qh=90.0, fu_b=1.0, clint_liver=100.0)
        assert result > 0.0
