"""Unit tests for charon.core.units conversion functions and constants."""

from __future__ import annotations

import pytest

from charon.core.units import (
    HUMAN_BODY_WEIGHT_KG,
    HUMAN_BSA_KM,
    HUMAN_GFR_ML_MIN,
    HUMAN_HEPATOCELLULARITY,
    HUMAN_LIVER_WEIGHT_G,
    HUMAN_MPPGL,
    HUMAN_QH_L_H,
    L_h_to_mL_min,
    cm_s_to_nm_s,
    mg_kg_to_mg,
    mL_min_to_L_h,
    ng_mL_to_nmol_L,
    nm_s_to_cm_s,
    nmol_L_to_ng_mL,
    uL_min_to_L_h,
)


# ---------------------------------------------------------------------------
# Physiological constants
# ---------------------------------------------------------------------------


class TestConstants:
    """Verify physiological reference constants have correct values."""

    def test_gfr(self):
        assert HUMAN_GFR_ML_MIN == 120.0

    def test_liver_weight(self):
        assert HUMAN_LIVER_WEIGHT_G == 1500.0

    def test_hepatic_blood_flow(self):
        assert HUMAN_QH_L_H == 90.0

    def test_mppgl(self):
        assert HUMAN_MPPGL == 40.0

    def test_hepatocellularity(self):
        assert HUMAN_HEPATOCELLULARITY == 120.0

    def test_body_weight(self):
        assert HUMAN_BODY_WEIGHT_KG == 70.0

    def test_bsa_km(self):
        assert HUMAN_BSA_KM == 37.0


# ---------------------------------------------------------------------------
# Volume-flow conversions
# ---------------------------------------------------------------------------


class TestULMinToLH:
    """uL/min -> L/h."""

    def test_basic(self):
        # 1e6 uL/min = 1 L/min = 60 L/h
        assert uL_min_to_L_h(1e6) == pytest.approx(60.0)

    def test_zero(self):
        assert uL_min_to_L_h(0.0) == 0.0

    def test_small_value(self):
        # 1000 uL/min = 0.001 L/min = 0.06 L/h
        assert uL_min_to_L_h(1000.0) == pytest.approx(0.06)

    def test_negative_raises(self):
        with pytest.raises(ValueError, match="non-negative"):
            uL_min_to_L_h(-1.0)


class TestMLMinToLH:
    """mL/min -> L/h."""

    def test_basic(self):
        # 1000 mL/min = 1 L/min = 60 L/h
        assert mL_min_to_L_h(1000.0) == pytest.approx(60.0)

    def test_gfr_conversion(self):
        # 120 mL/min = 7.2 L/h
        assert mL_min_to_L_h(120.0) == pytest.approx(7.2)

    def test_zero(self):
        assert mL_min_to_L_h(0.0) == 0.0

    def test_negative_raises(self):
        with pytest.raises(ValueError, match="non-negative"):
            mL_min_to_L_h(-5.0)


class TestLHToMLMin:
    """L/h -> mL/min."""

    def test_basic(self):
        # 60 L/h = 1000 mL/min
        assert L_h_to_mL_min(60.0) == pytest.approx(1000.0)

    def test_hepatic_flow(self):
        # 90 L/h = 1500 mL/min
        assert L_h_to_mL_min(90.0) == pytest.approx(1500.0)

    def test_zero(self):
        assert L_h_to_mL_min(0.0) == 0.0

    def test_negative_raises(self):
        with pytest.raises(ValueError, match="non-negative"):
            L_h_to_mL_min(-1.0)


class TestFlowRoundTrip:
    """Round-trip conversions should return to the original value."""

    @pytest.mark.parametrize("value", [0.0, 1.5, 7.2, 90.0, 1500.0])
    def test_mL_min_round_trip(self, value):
        assert L_h_to_mL_min(mL_min_to_L_h(value)) == pytest.approx(value)

    @pytest.mark.parametrize("value", [0.0, 1.5, 7.2, 90.0, 1500.0])
    def test_L_h_round_trip(self, value):
        assert mL_min_to_L_h(L_h_to_mL_min(value)) == pytest.approx(value)


# ---------------------------------------------------------------------------
# Permeability conversions
# ---------------------------------------------------------------------------


class TestNmSToeCmS:
    """nm/s -> cm/s."""

    def test_basic(self):
        # 100 nm/s = 100 * 1e-7 cm/s = 1e-5 cm/s
        assert nm_s_to_cm_s(100.0) == pytest.approx(1e-5)

    def test_one(self):
        assert nm_s_to_cm_s(1.0) == pytest.approx(1e-7)

    def test_zero(self):
        assert nm_s_to_cm_s(0.0) == 0.0

    def test_negative_raises(self):
        with pytest.raises(ValueError, match="non-negative"):
            nm_s_to_cm_s(-10.0)


class TestCmSToNmS:
    """cm/s -> nm/s."""

    def test_basic(self):
        # 1e-5 cm/s = 100 nm/s
        assert cm_s_to_nm_s(1e-5) == pytest.approx(100.0)

    def test_one(self):
        assert cm_s_to_nm_s(1e-7) == pytest.approx(1.0)

    def test_zero(self):
        assert cm_s_to_nm_s(0.0) == 0.0

    def test_negative_raises(self):
        with pytest.raises(ValueError, match="non-negative"):
            cm_s_to_nm_s(-1e-5)


class TestPermeabilityRoundTrip:
    """Round-trip permeability conversions."""

    @pytest.mark.parametrize("value", [0.0, 1.0, 50.0, 500.0, 10000.0])
    def test_nm_s_round_trip(self, value):
        assert cm_s_to_nm_s(nm_s_to_cm_s(value)) == pytest.approx(value)

    @pytest.mark.parametrize("value", [0.0, 1e-7, 1e-5, 1e-3])
    def test_cm_s_round_trip(self, value):
        assert nm_s_to_cm_s(cm_s_to_nm_s(value)) == pytest.approx(value)


# ---------------------------------------------------------------------------
# Dose / mass conversions
# ---------------------------------------------------------------------------


class TestMgKgToMg:
    """mg/kg -> mg."""

    def test_standard_body_weight(self):
        # 10 mg/kg * 70 kg = 700 mg
        assert mg_kg_to_mg(10.0, 70.0) == pytest.approx(700.0)

    def test_zero_dose(self):
        assert mg_kg_to_mg(0.0, 70.0) == 0.0

    def test_zero_weight(self):
        assert mg_kg_to_mg(10.0, 0.0) == 0.0

    def test_negative_dose_raises(self):
        with pytest.raises(ValueError, match="non-negative"):
            mg_kg_to_mg(-5.0, 70.0)

    def test_negative_weight_raises(self):
        with pytest.raises(ValueError, match="non-negative"):
            mg_kg_to_mg(5.0, -70.0)


# ---------------------------------------------------------------------------
# Concentration conversions
# ---------------------------------------------------------------------------


class TestNmolLToNgML:
    """nmol/L -> ng/mL."""

    def test_basic(self):
        # 1000 nmol/L * 300 g/mol / 1000 = 300 ng/mL
        assert nmol_L_to_ng_mL(1000.0, 300.0) == pytest.approx(300.0)

    def test_small_mw(self):
        # 500 nmol/L * 100 g/mol / 1000 = 50 ng/mL
        assert nmol_L_to_ng_mL(500.0, 100.0) == pytest.approx(50.0)

    def test_zero_concentration(self):
        assert nmol_L_to_ng_mL(0.0, 300.0) == 0.0

    def test_negative_concentration_raises(self):
        with pytest.raises(ValueError, match="non-negative"):
            nmol_L_to_ng_mL(-1.0, 300.0)

    def test_negative_mw_raises(self):
        with pytest.raises(ValueError, match="non-negative"):
            nmol_L_to_ng_mL(100.0, -300.0)


class TestNgMLToNmolL:
    """ng/mL -> nmol/L."""

    def test_basic(self):
        # 300 ng/mL * 1000 / 300 g/mol = 1000 nmol/L
        assert ng_mL_to_nmol_L(300.0, 300.0) == pytest.approx(1000.0)

    def test_small_mw(self):
        # 50 ng/mL * 1000 / 100 g/mol = 500 nmol/L
        assert ng_mL_to_nmol_L(50.0, 100.0) == pytest.approx(500.0)

    def test_zero_concentration(self):
        assert ng_mL_to_nmol_L(0.0, 300.0) == 0.0

    def test_negative_concentration_raises(self):
        with pytest.raises(ValueError, match="non-negative"):
            ng_mL_to_nmol_L(-1.0, 300.0)

    def test_negative_mw_raises(self):
        with pytest.raises(ValueError, match="non-negative"):
            ng_mL_to_nmol_L(100.0, -300.0)


class TestConcentrationRoundTrip:
    """Round-trip concentration conversions."""

    @pytest.mark.parametrize(
        "nmol, mw",
        [(100.0, 200.0), (500.0, 325.0), (1000.0, 180.16), (50.0, 500.0)],
    )
    def test_nmol_L_round_trip(self, nmol, mw):
        ng = nmol_L_to_ng_mL(nmol, mw)
        assert ng_mL_to_nmol_L(ng, mw) == pytest.approx(nmol)

    @pytest.mark.parametrize(
        "ng, mw",
        [(100.0, 200.0), (50.0, 325.0), (300.0, 180.16)],
    )
    def test_ng_mL_round_trip(self, ng, mw):
        nmol = ng_mL_to_nmol_L(ng, mw)
        assert nmol_L_to_ng_mL(nmol, mw) == pytest.approx(ng)
