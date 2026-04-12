"""Tests for PAD (Pharmacologically Active Dose) projection.

Hand-calculations use midazolam-like reference compound:
    Ceff=100 nM, MW=325.77 g/mol, CL/F=18.6 L/h, tau=24 h
"""

import pytest

from charon.translational.pad import PADResult, compute_pad


# ---------------------------------------------------------------------------
# 1. Canonical hand-calc
# ---------------------------------------------------------------------------
def test_hand_calc() -> None:
    """Hand-calc (tau=24 h):
        ceff_mg_L  = 100 nM × 325.77 g/mol × 1e-6  = 0.032577 mg/L
        auc_target = 0.032577 × 24                  = 0.781848 mg·h/L
        dose       = 0.781848 × 18.6                = 14.5424 mg
        MRSD       = 14.5424 / 10                   = 1.45424 mg
    """
    result = compute_pad(
        target_ceff_nM=100.0,
        molecular_weight=325.77,
        cl_apparent_L_h=18.6,
    )
    expected_ceff = 100.0 * 325.77 * 1e-6        # 0.032577 mg/L
    expected_auc = expected_ceff * 24.0           # 0.781848 mg·h/L
    expected_dose = expected_auc * 18.6           # 14.5424 mg
    expected_mrsd = expected_dose / 10.0          # 1.45424 mg

    assert result.target_conc_mg_L == pytest.approx(expected_ceff, rel=1e-6)
    assert result.auc_target_mg_h_L == pytest.approx(expected_auc, rel=1e-6)
    assert result.dose_mg == pytest.approx(expected_dose, rel=1e-6)
    assert result.mrsd_mg == pytest.approx(expected_mrsd, rel=1e-6)


# ---------------------------------------------------------------------------
# 2. Half tau gives half dose
# ---------------------------------------------------------------------------
def test_short_tau_halves_dose() -> None:
    """tau=12 h should give exactly half the dose of tau=24 h (linear in tau)."""
    r24 = compute_pad(
        target_ceff_nM=50.0,
        molecular_weight=300.0,
        cl_apparent_L_h=10.0,
        tau_h=24.0,
    )
    r12 = compute_pad(
        target_ceff_nM=50.0,
        molecular_weight=300.0,
        cl_apparent_L_h=10.0,
        tau_h=12.0,
    )
    assert r12.dose_mg == pytest.approx(r24.dose_mg / 2.0, rel=1e-6)
    assert r12.mrsd_mg == pytest.approx(r24.mrsd_mg / 2.0, rel=1e-6)


# ---------------------------------------------------------------------------
# 3. All result fields present and correctly typed
# ---------------------------------------------------------------------------
def test_result_fields() -> None:
    """PADResult must expose every declared field with correct types."""
    result = compute_pad(
        target_ceff_nM=200.0,
        molecular_weight=400.0,
        cl_apparent_L_h=5.0,
        safety_factor=3.0,
        tau_h=12.0,
    )
    assert isinstance(result, PADResult)
    assert result.target_ceff_nM == pytest.approx(200.0)
    assert result.cl_apparent_L_h == pytest.approx(5.0)
    assert result.safety_factor == pytest.approx(3.0)
    assert result.tau_h == pytest.approx(12.0)
    assert isinstance(result.target_conc_mg_L, float)
    assert isinstance(result.auc_target_mg_h_L, float)
    assert isinstance(result.dose_mg, float)
    assert isinstance(result.mrsd_mg, float)
    # MRSD = dose / safety_factor
    assert result.mrsd_mg == pytest.approx(result.dose_mg / result.safety_factor,
                                           rel=1e-6)


# ---------------------------------------------------------------------------
# 4. Negative Ceff raises ValueError
# ---------------------------------------------------------------------------
def test_negative_ceff_raises() -> None:
    """target_ceff_nM <= 0 must raise ValueError."""
    with pytest.raises(ValueError, match="target_ceff_nM must be > 0"):
        compute_pad(
            target_ceff_nM=-10.0,
            molecular_weight=300.0,
            cl_apparent_L_h=10.0,
        )
    with pytest.raises(ValueError, match="target_ceff_nM must be > 0"):
        compute_pad(
            target_ceff_nM=0.0,
            molecular_weight=300.0,
            cl_apparent_L_h=10.0,
        )


# ---------------------------------------------------------------------------
# 5. Zero CL raises ValueError
# ---------------------------------------------------------------------------
def test_zero_cl_raises() -> None:
    """cl_apparent_L_h <= 0 must raise ValueError (physically meaningless)."""
    with pytest.raises(ValueError, match="cl_apparent_L_h must be > 0"):
        compute_pad(
            target_ceff_nM=50.0,
            molecular_weight=300.0,
            cl_apparent_L_h=0.0,
        )
    with pytest.raises(ValueError, match="cl_apparent_L_h must be > 0"):
        compute_pad(
            target_ceff_nM=50.0,
            molecular_weight=300.0,
            cl_apparent_L_h=-1.0,
        )
