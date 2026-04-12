"""Tests for MABEL (Minimum Anticipated Biological Effect Level) dose projection.

All hand-calculations use:
    Kd=10 nM, MW=325.77 g/mol (midazolam-like), fu_p=0.03, CL/F=18.6 L/h
"""

import math

import pytest

from charon.translational.mabel import MABELResult, compute_mabel


# ---------------------------------------------------------------------------
# Helper constants (midazolam-inspired)
# ---------------------------------------------------------------------------
KD_NM = 10.0
MW = 325.77
FU_P = 0.03
CL_F = 18.6   # L/h  (CL/F for oral)
# Vd = t1/2 × CL / ln2; use t1/2 ~ 1.7 h for midazolam-like
T_HALF_H = 1.7
VD_L = T_HALF_H * CL_F / math.log(2)   # ≈ 45.6 L


# ---------------------------------------------------------------------------
# 1. Midazolam hand-calc (Cmax approach likely limiting for large Vd)
# ---------------------------------------------------------------------------
def test_midazolam_hand_calc() -> None:
    """Hand-calc:
        css_u   = 10e-9 mol/L × 325.77 g/mol × 1e3 mg/g  = 3.2577e-3 mg/L
                   [equivalently: 10 nM × 325.77 × 1e-6]
        css_t   = 3.2577e-3 / 0.03  = 0.10859 mg/L
        Vd      = 1.7 × 18.6 / ln2  ≈ 45.60 L
        dose_cx = 0.10859 × 45.60   ≈ 4.9517 mg
        dose_ss = 0.10859 × 18.6 × 24 ≈ 48.495 mg
        limiting = "cmax"  (dose_cx < dose_ss)
        MRSD    = 4.9517 / 10 ≈ 0.4952 mg
    """
    result = compute_mabel(
        target_kd_nM=KD_NM,
        molecular_weight=MW,
        fu_p=FU_P,
        cl_apparent_L_h=CL_F,
        vd_apparent_L=VD_L,
    )
    css_u = KD_NM * MW * 1e-6
    css_t = css_u / FU_P
    expected_dose_cx = css_t * VD_L
    expected_dose_ss = css_t * CL_F * 24.0
    assert result.dose_cmax_mg == pytest.approx(expected_dose_cx, rel=1e-5)
    assert result.dose_ss_mg == pytest.approx(expected_dose_ss, rel=1e-5)
    assert result.limiting_approach == "cmax"
    assert result.mrsd_mg == pytest.approx(expected_dose_cx / 10.0, rel=1e-5)


# ---------------------------------------------------------------------------
# 2. Target concentration unit conversion: nM → mg/L
# ---------------------------------------------------------------------------
def test_target_conc_conversion() -> None:
    """Verify nM → mg/L conversion: c_total = kd_nM × MW × 1e-6 / fu_p."""
    result = compute_mabel(
        target_kd_nM=100.0,
        molecular_weight=400.0,
        fu_p=0.1,
        cl_apparent_L_h=10.0,
        vd_apparent_L=100.0,
    )
    expected_css_u = 100.0 * 400.0 * 1e-6   # 0.04 mg/L
    expected_css_t = expected_css_u / 0.1    # 0.4 mg/L
    assert result.target_conc_mg_L == pytest.approx(expected_css_t, rel=1e-6)


# ---------------------------------------------------------------------------
# 3. Steady-state limiting when Vd is small and tau is short
# ---------------------------------------------------------------------------
def test_ss_limiting_with_small_vd_short_tau() -> None:
    """With very small Vd and short tau, SS dose < Cmax dose → SS limits.

    Hand-calc:
        Vd_small = 1.0 L
        tau      = 1.0 h
        css_t    = kd_nM × MW × 1e-6 / fu_p
        dose_cx  = css_t × 1.0
        dose_ss  = css_t × CL × 1.0
        If CL < 1 L/h, dose_ss < dose_cx → SS limiting.
    """
    result = compute_mabel(
        target_kd_nM=50.0,
        molecular_weight=300.0,
        fu_p=0.5,
        cl_apparent_L_h=0.5,   # small CL
        vd_apparent_L=50.0,    # larger Vd
        tau_h=1.0,             # short tau
    )
    # dose_ss = css_t × 0.5 × 1 < css_t × 50 = dose_cx  → SS limits
    assert result.limiting_approach == "steady_state"
    assert result.mrsd_mg == pytest.approx(result.dose_ss_mg / result.safety_factor,
                                           rel=1e-6)


# ---------------------------------------------------------------------------
# 4. All result fields present and correctly typed
# ---------------------------------------------------------------------------
def test_result_fields() -> None:
    """MABELResult must expose every declared field."""
    result = compute_mabel(
        target_kd_nM=5.0,
        molecular_weight=350.0,
        fu_p=0.05,
        cl_apparent_L_h=5.0,
        vd_apparent_L=20.0,
        safety_factor=3.0,
        tau_h=12.0,
    )
    assert isinstance(result, MABELResult)
    assert result.target_kd_nM == pytest.approx(5.0)
    assert result.fu_p == pytest.approx(0.05)
    assert result.molecular_weight == pytest.approx(350.0)
    assert result.safety_factor == pytest.approx(3.0)
    assert result.tau_h == pytest.approx(12.0)
    assert isinstance(result.dose_cmax_mg, float)
    assert isinstance(result.dose_ss_mg, float)
    assert isinstance(result.mrsd_mg, float)
    assert result.limiting_approach in ("cmax", "steady_state")


# ---------------------------------------------------------------------------
# 5. Negative Kd raises ValueError
# ---------------------------------------------------------------------------
def test_negative_kd_raises() -> None:
    """target_kd_nM <= 0 must raise ValueError."""
    with pytest.raises(ValueError, match="target_kd_nM must be > 0"):
        compute_mabel(
            target_kd_nM=-1.0,
            molecular_weight=300.0,
            fu_p=0.1,
            cl_apparent_L_h=10.0,
            vd_apparent_L=50.0,
        )
    with pytest.raises(ValueError, match="target_kd_nM must be > 0"):
        compute_mabel(
            target_kd_nM=0.0,
            molecular_weight=300.0,
            fu_p=0.1,
            cl_apparent_L_h=10.0,
            vd_apparent_L=50.0,
        )


# ---------------------------------------------------------------------------
# 6. fu_p = 0 raises ValueError
# ---------------------------------------------------------------------------
def test_fu_p_zero_raises() -> None:
    """fu_p <= 0 must raise ValueError (would cause division by zero)."""
    with pytest.raises(ValueError, match="fu_p must be > 0"):
        compute_mabel(
            target_kd_nM=10.0,
            molecular_weight=300.0,
            fu_p=0.0,
            cl_apparent_L_h=10.0,
            vd_apparent_L=50.0,
        )
