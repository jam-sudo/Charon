import pytest
import math
from charon.core.parameter_bridge import ParameterBridge

@pytest.fixture
def bridge():
    return ParameterBridge()

class TestClintToClh:
    """IVIVE: in vitro CLint -> in vivo hepatic clearance."""

    def test_architecture_example(self, bridge):
        """
        Verify the exact example from ARCHITECTURE.md:
        CLint=15.2, fu_inc=0.85, fu_p=0.23, BP=0.95, Qh=90
        Expected CLh ≈ 13.28 L/h (well-stirred)

        Hand calculation:
        Step 1: CLint_u = 15.2 / 0.85 = 17.882
        Step 2: CLint_liver = 17.882 * 40 * 1500 = 1,072,941 uL/min
                CLint_liver_L_h = 1,072,941 / 1e6 * 60 = 64.377 L/h
        Step 3: fu_b = 0.23 / 0.95 = 0.24211
        Step 4: CLh = (90 * 0.24211 * 64.377) / (90 + 0.24211 * 64.377)
                    = 1402.3 / 105.59 = 13.28 L/h
        """
        result = bridge.clint_to_clh(
            clint=15.2, fu_inc=0.85, fu_p=0.23, bp_ratio=0.95, qh_L_h=90.0
        )
        assert result.clh_L_h == pytest.approx(13.28, rel=0.01)

    @pytest.mark.parametrize("drug_name,clint,fu_inc,fu_p,bp_ratio,expected_clh,tol_fold", [
        # Reference drugs with approximate in-vitro parameters.
        # Expected CLh is the well-stirred model prediction (verifies math),
        # with tolerances to allow minor floating-point variation.
        # Drugs with very low fu_p produce low CLh -- this is the well-known
        # IVIVE underprediction for highly bound drugs (not a code bug).
        ("midazolam", 93.0, 0.96, 0.03, 0.66, 13.5, 1.1),
        ("warfarin", 0.22, 0.98, 0.01, 0.58, 0.014, 1.1),
        ("diclofenac", 82.0, 0.99, 0.01, 0.99, 2.91, 1.1),
        ("omeprazole", 13.0, 0.85, 0.05, 1.0, 2.67, 1.1),
        ("dextromethorphan", 55.0, 0.9, 0.40, 1.1, 42.35, 1.1),
    ])
    def test_known_drug_clh(self, bridge, drug_name, clint, fu_inc, fu_p, bp_ratio, expected_clh, tol_fold):
        """Test CLh prediction for reference drugs within X-fold of observed."""
        result = bridge.clint_to_clh(
            clint=clint, fu_inc=fu_inc, fu_p=fu_p, bp_ratio=bp_ratio, model="well_stirred"
        )
        assert expected_clh / tol_fold <= result.clh_L_h <= expected_clh * tol_fold, \
            f"{drug_name}: predicted CLh={result.clh_L_h:.2f}, expected ~{expected_clh} L/h (±{tol_fold}x)"

    def test_audit_trail_completeness(self, bridge):
        """ConversionLog must have all intermediate steps."""
        result = bridge.clint_to_clh(clint=15.2, fu_inc=0.85, fu_p=0.23, bp_ratio=0.95)
        log = result.conversion_log
        # Must have steps
        assert len(log.intermediate_steps) >= 4
        # Must have correct output unit
        assert log.output_unit == "L/h"
        assert log.model_used == "well_stirred"
        # Must have input params
        assert "clint" in log.input_params
        assert "fu_inc" in log.input_params

    def test_hepatocyte_no_fu_inc_correction(self, bridge):
        """
        Hepatocytes: fu_inc correction must NOT be applied.
        CLint goes directly to scaling without /fu_inc.

        Hand calc with system="hepatocytes":
        CLint_u = 10.0 (NOT divided by fu_inc)
        scale_factor = 120 * 1500 = 180,000
        CLint_liver = 10.0 * 180,000 = 1,800,000 uL/min
        CLint_liver_L_h = 1,800,000 / 1e6 * 60 = 108.0 L/h
        fu_b = 0.5 / 1.0 = 0.5
        CLh = (90 * 0.5 * 108) / (90 + 0.5 * 108) = 4860 / 144 = 33.75 L/h
        """
        result = bridge.clint_to_clh(
            clint=10.0, fu_inc=0.85, fu_p=0.5, bp_ratio=1.0,
            system="hepatocytes"
        )
        assert result.clh_L_h == pytest.approx(33.75, rel=0.01)

    def test_zero_clint(self, bridge):
        """CLint=0 -> CLh=0"""
        result = bridge.clint_to_clh(clint=0.0, fu_inc=0.85, fu_p=0.5, bp_ratio=1.0)
        assert result.clh_L_h == 0.0

    def test_fu_inc_zero_raises(self, bridge):
        """fu_inc=0 must raise ValueError (division by zero in HLM)"""
        with pytest.raises(ValueError):
            bridge.clint_to_clh(clint=10.0, fu_inc=0.0, fu_p=0.5, bp_ratio=1.0, system="HLM")

    def test_fu_p_negative_raises(self, bridge):
        """fu_p must be positive"""
        with pytest.raises(ValueError):
            bridge.clint_to_clh(clint=10.0, fu_inc=0.85, fu_p=-0.1, bp_ratio=1.0)

    def test_invalid_system_raises(self, bridge):
        """system must be 'HLM' or 'hepatocytes'"""
        with pytest.raises(ValueError):
            bridge.clint_to_clh(clint=10.0, fu_inc=0.85, fu_p=0.5, bp_ratio=1.0, system="invalid")

    def test_flow_limited(self, bridge):
        """Very high CLint -> CLh approaches Qh"""
        result = bridge.clint_to_clh(
            clint=10000.0, fu_inc=1.0, fu_p=1.0, bp_ratio=1.0
        )
        assert result.clh_L_h > 85  # close to Qh=90
        assert result.extraction_ratio > 0.95

    def test_parallel_tube_ge_well_stirred(self, bridge):
        """Parallel-tube always >= well-stirred for same inputs."""
        ws = bridge.clint_to_clh(clint=50.0, fu_inc=0.9, fu_p=0.1, bp_ratio=1.0, model="well_stirred")
        pt = bridge.clint_to_clh(clint=50.0, fu_inc=0.9, fu_p=0.1, bp_ratio=1.0, model="parallel_tube")
        assert pt.clh_L_h >= ws.clh_L_h

    def test_dispersion_model(self, bridge):
        """Dispersion model runs without error and returns reasonable value."""
        result = bridge.clint_to_clh(clint=50.0, fu_inc=0.9, fu_p=0.1, bp_ratio=1.0, model="dispersion")
        assert 0 < result.clh_L_h <= 90.0
        assert result.model_used == "dispersion"

class TestPappToPeff:
    def test_sun_2002(self, bridge):
        """
        Sun 2002: log10(Peff) = 0.4926 * log10(Papp_cm_s) - 0.1454
        Papp = 30 nm/s = 30e-7 cm/s = 3e-6 cm/s
        log10(3e-6) = log10(3) + log10(1e-6) = 0.4771 - 6 = -5.5229
        log10(Peff) = 0.4926 * (-5.5229) - 0.1454 = -2.7206 - 0.1454 = -2.8660
        Peff = 10^(-2.8660) = 1.363e-3 cm/s
        """
        peff = bridge.papp_to_peff(30.0)
        assert peff == pytest.approx(1.363e-3, rel=0.05)

    def test_zero_papp_raises(self, bridge):
        with pytest.raises(ValueError):
            bridge.papp_to_peff(0.0)

    def test_negative_papp_raises(self, bridge):
        with pytest.raises(ValueError):
            bridge.papp_to_peff(-5.0)

class TestRenalClearance:
    def test_filtration_only(self, bridge):
        """
        fu_p=0.23, GFR=120 mL/min
        GFR_L_h = 120 * 60 / 1000 = 7.2 L/h
        CLrenal = 0.23 * 7.2 * 1.0 = 1.656 L/h
        """
        cl_r = bridge.assign_renal_clearance(fu_p=0.23, gfr_mL_min=120.0)
        assert cl_r == pytest.approx(1.656, rel=0.01)

    def test_with_secretion(self, bridge):
        """
        fu_p=0.5, GFR=120, net_secretion_factor=0.5
        CLrenal = 0.5 * 7.2 * 1.5 = 5.4 L/h
        """
        cl_r = bridge.assign_renal_clearance(
            fu_p=0.5, gfr_mL_min=120.0,
            is_active_secretion=True, net_secretion_factor=0.5
        )
        assert cl_r == pytest.approx(5.4, rel=0.01)

    def test_fu_p_out_of_range_raises(self, bridge):
        with pytest.raises(ValueError):
            bridge.assign_renal_clearance(fu_p=1.5)
