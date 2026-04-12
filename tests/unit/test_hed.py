"""Tests for HED (Human Equivalent Dose) computation.

All hand-calculations use FDA Guidance (2005) Table 1 Km values:
    rat=6.2, mouse=3.0, dog=20.0, monkey=12.0, rabbit=12.0,
    guinea_pig=8.0, human=37.0
"""

import math

import pytest

from charon.translational.hed import HEDResult, KM_BY_SPECIES, compute_hed


# ---------------------------------------------------------------------------
# 1. Canonical rat NOAEL = 50 mg/kg
# ---------------------------------------------------------------------------
def test_rat_noael_50() -> None:
    """Hand-calc:
        HED  = 50 × (6.2 / 37) = 8.37838...
        MRSD = 8.37838 × 70 / 10 = 58.648...
    """
    result = compute_hed(noael_mg_kg=50.0, noael_species="rat")
    expected_hed = 50.0 * (6.2 / 37.0)
    expected_mrsd = expected_hed * 70.0 / 10.0
    assert result.hed_mg_kg == pytest.approx(expected_hed, rel=1e-6)
    assert result.mrsd_mg == pytest.approx(expected_mrsd, rel=1e-6)


# ---------------------------------------------------------------------------
# 2. Dog NOAEL = 10 mg/kg
# ---------------------------------------------------------------------------
def test_dog_noael_10() -> None:
    """Hand-calc:
        HED  = 10 × (20 / 37) = 5.40540...
        MRSD = 5.40540 × 70 / 10 = 37.8378...
    """
    result = compute_hed(noael_mg_kg=10.0, noael_species="dog")
    expected_hed = 10.0 * (20.0 / 37.0)
    assert result.hed_mg_kg == pytest.approx(expected_hed, rel=1e-6)
    assert result.km_animal == pytest.approx(20.0)
    assert result.km_human == pytest.approx(37.0)


# ---------------------------------------------------------------------------
# 3. Monkey NOAEL = 25 mg/kg, Km = 12
# ---------------------------------------------------------------------------
def test_monkey_noael_25() -> None:
    """Hand-calc:
        HED  = 25 × (12 / 37) = 8.10810...
    """
    result = compute_hed(noael_mg_kg=25.0, noael_species="monkey")
    expected_hed = 25.0 * (12.0 / 37.0)
    assert result.hed_mg_kg == pytest.approx(expected_hed, rel=1e-6)
    assert result.km_animal == pytest.approx(12.0)


# ---------------------------------------------------------------------------
# 4. Mouse NOAEL = 100 mg/kg, Km = 3
# ---------------------------------------------------------------------------
def test_mouse_noael_100() -> None:
    """Hand-calc:
        HED  = 100 × (3 / 37) = 8.10810...
        MRSD = 8.10810 × 70 / 10 = 56.7567...
    """
    result = compute_hed(noael_mg_kg=100.0, noael_species="mouse")
    expected_hed = 100.0 * (3.0 / 37.0)
    assert result.hed_mg_kg == pytest.approx(expected_hed, rel=1e-6)


# ---------------------------------------------------------------------------
# 5. Case-insensitive species lookup
# ---------------------------------------------------------------------------
def test_case_insensitive_species() -> None:
    """'Rat', 'RAT', 'rat' should all resolve to Km=6.2."""
    r1 = compute_hed(noael_mg_kg=10.0, noael_species="Rat")
    r2 = compute_hed(noael_mg_kg=10.0, noael_species="RAT")
    r3 = compute_hed(noael_mg_kg=10.0, noael_species="rat")
    assert r1.hed_mg_kg == pytest.approx(r2.hed_mg_kg)
    assert r1.hed_mg_kg == pytest.approx(r3.hed_mg_kg)
    assert r1.noael_species == "rat"


# ---------------------------------------------------------------------------
# 6. Unknown species raises ValueError
# ---------------------------------------------------------------------------
def test_unknown_species_raises() -> None:
    """'hamster' is not in KM_BY_SPECIES and must raise ValueError."""
    with pytest.raises(ValueError, match="Unknown species"):
        compute_hed(noael_mg_kg=10.0, noael_species="hamster")


# ---------------------------------------------------------------------------
# 7. Custom safety factor
# ---------------------------------------------------------------------------
def test_custom_safety_factor() -> None:
    """MRSD with SF=3 should be 10/3 × MRSD with SF=10 for same inputs."""
    r10 = compute_hed(noael_mg_kg=50.0, noael_species="rat", safety_factor=10.0)
    r3 = compute_hed(noael_mg_kg=50.0, noael_species="rat", safety_factor=3.0)
    assert r3.mrsd_mg == pytest.approx(r10.mrsd_mg * 10.0 / 3.0, rel=1e-6)


# ---------------------------------------------------------------------------
# 8. Custom body weight
# ---------------------------------------------------------------------------
def test_custom_body_weight() -> None:
    """BW=50 kg should give (50/70) × MRSD compared to default BW=70."""
    r70 = compute_hed(noael_mg_kg=50.0, noael_species="rat")
    r50 = compute_hed(noael_mg_kg=50.0, noael_species="rat", body_weight_kg=50.0)
    assert r50.mrsd_mg == pytest.approx(r70.mrsd_mg * 50.0 / 70.0, rel=1e-6)


# ---------------------------------------------------------------------------
# 9. All result fields are present and correctly typed
# ---------------------------------------------------------------------------
def test_result_fields() -> None:
    """HEDResult must expose every field with correct types and values."""
    result = compute_hed(noael_mg_kg=20.0, noael_species="rat",
                         safety_factor=5.0, body_weight_kg=60.0)
    assert isinstance(result, HEDResult)
    assert result.noael_mg_kg == pytest.approx(20.0)
    assert result.noael_species == "rat"
    assert result.km_animal == pytest.approx(6.2)
    assert result.km_human == pytest.approx(37.0)
    assert isinstance(result.hed_mg_kg, float)
    assert result.body_weight_kg == pytest.approx(60.0)
    assert result.safety_factor == pytest.approx(5.0)
    assert isinstance(result.mrsd_mg, float)


# ---------------------------------------------------------------------------
# 10. Negative NOAEL raises ValueError
# ---------------------------------------------------------------------------
def test_negative_noael_raises() -> None:
    """noael_mg_kg <= 0 must raise ValueError."""
    with pytest.raises(ValueError, match="noael_mg_kg must be > 0"):
        compute_hed(noael_mg_kg=-5.0, noael_species="rat")
    with pytest.raises(ValueError, match="noael_mg_kg must be > 0"):
        compute_hed(noael_mg_kg=0.0, noael_species="rat")
