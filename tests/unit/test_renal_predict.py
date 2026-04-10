"""Unit tests for the Layer-1 renal clearance wrapper."""

from __future__ import annotations

import pytest

from charon.predict.renal import estimate_renal_clearance


class TestEstimateRenalClearance:
    """Thin wrapper around ``ParameterBridge.assign_renal_clearance``."""

    def test_default_gfr_fu_p_023(self):
        """fu_p=0.23, GFR=120 mL/min → 0.23 * 7.2 = 1.656 L/h."""
        cl = estimate_renal_clearance(fu_p=0.23, gfr_mL_min=120.0)
        assert cl == pytest.approx(1.656, rel=1e-6)

    def test_with_net_secretion_factor(self):
        """fu_p=0.5, default GFR, secretion_factor=0.5 → 0.5*7.2*1.5 = 5.4 L/h."""
        cl = estimate_renal_clearance(
            fu_p=0.5, net_secretion_factor=0.5
        )
        assert cl == pytest.approx(5.4, rel=1e-6)

    def test_fu_p_zero(self):
        """Fully bound drug → no filtration clearance."""
        cl = estimate_renal_clearance(fu_p=0.0)
        assert cl == pytest.approx(0.0)

    def test_fu_p_one_equals_gfr(self):
        """Unbound drug filters at GFR exactly."""
        cl = estimate_renal_clearance(fu_p=1.0, gfr_mL_min=120.0)
        assert cl == pytest.approx(7.2, rel=1e-6)

    @pytest.mark.parametrize("fu_p", [-0.01, -1.0, 1.01, 2.0])
    def test_fu_p_out_of_range_raises(self, fu_p):
        with pytest.raises(ValueError, match="fu_p"):
            estimate_renal_clearance(fu_p=fu_p)

    def test_active_secretion_flag_passthrough(self):
        """is_active_secretion=True is a flag; mathematics unchanged unless
        net_secretion_factor is set."""
        cl_flag_only = estimate_renal_clearance(
            fu_p=0.5, is_active_secretion=True
        )
        cl_no_flag = estimate_renal_clearance(fu_p=0.5)
        assert cl_flag_only == pytest.approx(cl_no_flag)

    def test_custom_gfr(self):
        """Halving GFR halves the clearance for a given fu_p."""
        cl_full = estimate_renal_clearance(fu_p=0.5, gfr_mL_min=120.0)
        cl_half = estimate_renal_clearance(fu_p=0.5, gfr_mL_min=60.0)
        assert cl_half == pytest.approx(cl_full / 2.0, rel=1e-6)
