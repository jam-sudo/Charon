"""Unit tests for the Kp (tissue:plasma partition) calculator.

Verifies EXACT numerical parity with Sisyphus ±1e-10 across 5 reference
drugs × 15 tissues × 2 methods. Also exercises the Berezhkovskiy
correction and input validation.
"""

from __future__ import annotations

import sys

import pytest

# Import Sisyphus source of truth for side-by-side comparison.
sys.path.insert(0, "/home/jam/Sisyphus/src")
from sisyphus.predict.ivive import (  # noqa: E402
    _compute_kp_poulin_theil as sis_pt,
    _compute_kp_rodgers_rowland as sis_rr,
    _PLASMA_COMP as sis_plasma,
    _TISSUE_COMPOSITIONS as sis_tissues,
)

from charon.pbpk.kp_calculator import (  # noqa: E402
    KP_MAX,
    KP_MIN,
    TissueComposition,
    apply_berezhkovskiy_correction,
    compute_all_kp,
    compute_kp_poulin_theil,
    compute_kp_rodgers_rowland,
)

TOL = 1e-10

# Charon-native compositions built from Sisyphus numerical values, so
# Charon's functions get Charon types while Sisyphus is fed its own.
CHARON_PLASMA = TissueComposition(
    fn=sis_plasma.fn, fp=sis_plasma.fp, fw=sis_plasma.fw, pH=sis_plasma.pH
)
CHARON_TISSUES: dict[str, TissueComposition] = {
    name: TissueComposition(fn=t.fn, fp=t.fp, fw=t.fw, pH=t.pH)
    for name, t in sis_tissues.items()
}

# Reference drugs: (name, logP, pka, compound_type)
REFERENCE_DRUGS: list[tuple[str, float, float | None, str]] = [
    ("midazolam", 3.1, 6.15, "base"),
    ("warfarin", 2.7, 5.05, "acid"),
    ("caffeine", -0.07, None, "neutral"),
    ("propranolol", 3.48, 9.45, "base"),
    ("aspirin", 1.19, 3.5, "acid"),
]


class TestRodgersRowlandParityWithSisyphus:
    """5 drugs × 15 tissues: Charon R&R ≡ Sisyphus R&R."""

    @pytest.mark.parametrize("drug", REFERENCE_DRUGS, ids=lambda d: d[0])
    def test_rr_matches_sisyphus_all_tissues(self, drug):
        name, logp, pka, compound_type = drug
        for tissue_name, sis_comp in sis_tissues.items():
            char_comp = CHARON_TISSUES[tissue_name]
            sis_kp = sis_rr(logp, pka, compound_type, sis_comp, sis_plasma)
            char_kp = compute_kp_rodgers_rowland(
                logp, pka, compound_type, char_comp, CHARON_PLASMA
            )
            assert abs(sis_kp - char_kp) < TOL, (
                f"{name}/{tissue_name}: Charon={char_kp} Sisyphus={sis_kp}"
            )


class TestPoulinTheilParityWithSisyphus:
    """5 drugs × 15 tissues: Charon P&T ≡ Sisyphus P&T."""

    @pytest.mark.parametrize("drug", REFERENCE_DRUGS, ids=lambda d: d[0])
    def test_pt_matches_sisyphus_all_tissues(self, drug):
        name, logp, pka, compound_type = drug
        for tissue_name, sis_comp in sis_tissues.items():
            char_comp = CHARON_TISSUES[tissue_name]
            sis_kp = sis_pt(logp, pka, compound_type, sis_comp, sis_plasma)
            char_kp = compute_kp_poulin_theil(
                logp, pka, compound_type, char_comp, CHARON_PLASMA
            )
            assert abs(sis_kp - char_kp) < TOL, (
                f"{name}/{tissue_name}: Charon={char_kp} Sisyphus={sis_kp}"
            )


class TestNeutralBranchIgnoresPka:
    """Neutral branch should ignore pka entirely."""

    @pytest.mark.parametrize("pka_override", [None, 3.0, 7.0, 10.5])
    def test_rr_neutral_ignores_pka(self, pka_override):
        kp_none = compute_kp_rodgers_rowland(
            logp=2.0,
            pka=None,
            compound_type="neutral",
            tissue_comp=CHARON_TISSUES["liver"],
            plasma_comp=CHARON_PLASMA,
        )
        kp_override = compute_kp_rodgers_rowland(
            logp=2.0,
            pka=pka_override,
            compound_type="neutral",
            tissue_comp=CHARON_TISSUES["liver"],
            plasma_comp=CHARON_PLASMA,
        )
        assert kp_none == pytest.approx(kp_override, abs=TOL)

    @pytest.mark.parametrize("pka_override", [None, 3.0, 7.0, 10.5])
    def test_pt_neutral_ignores_pka(self, pka_override):
        kp_none = compute_kp_poulin_theil(
            logp=2.0,
            pka=None,
            compound_type="neutral",
            tissue_comp=CHARON_TISSUES["muscle_tissue"],
            plasma_comp=CHARON_PLASMA,
        )
        kp_override = compute_kp_poulin_theil(
            logp=2.0,
            pka=pka_override,
            compound_type="neutral",
            tissue_comp=CHARON_TISSUES["muscle_tissue"],
            plasma_comp=CHARON_PLASMA,
        )
        assert kp_none == pytest.approx(kp_override, abs=TOL)


class TestKpClipping:
    """Kp always clipped to [KP_MIN, KP_MAX]."""

    @pytest.mark.parametrize("drug", REFERENCE_DRUGS, ids=lambda d: d[0])
    def test_rr_results_in_range(self, drug):
        _, logp, pka, compound_type = drug
        for tissue_comp in CHARON_TISSUES.values():
            kp = compute_kp_rodgers_rowland(
                logp, pka, compound_type, tissue_comp, CHARON_PLASMA
            )
            assert KP_MIN <= kp <= KP_MAX

    def test_extreme_high_logp_clipped(self):
        """Very high logP saturates at KP_MAX."""
        kp = compute_kp_rodgers_rowland(
            logp=15.0,
            pka=None,
            compound_type="neutral",
            tissue_comp=CHARON_TISSUES["adipose_tissue"],
            plasma_comp=CHARON_PLASMA,
        )
        assert kp == pytest.approx(KP_MAX)

    def test_extreme_low_logp_non_negative(self):
        """Very low logP drives Kp toward its lower bound (≥ KP_MIN)."""
        kp = compute_kp_rodgers_rowland(
            logp=-15.0,
            pka=None,
            compound_type="neutral",
            tissue_comp=CHARON_TISSUES["adipose_tissue"],
            plasma_comp=CHARON_PLASMA,
        )
        assert KP_MIN <= kp <= KP_MAX


class TestBerezhkovskiyCorrection:
    def test_highly_bound_reduces_kp(self):
        """fu_p=0.01 with Kp=100 should drop meaningfully below 100."""
        corrected = apply_berezhkovskiy_correction(kp=100.0, fu_p=0.01)
        assert corrected < 100.0
        assert corrected > 0.0

    def test_fu_p_one_hand_calc(self):
        """Explicit closed-form verification at fu_p=1.

        Kp_bz = Kp / (1 + (Kp - 1) * 1) = Kp / Kp = 1.
        The task brief refers to this as "identity"; mathematically it
        is a fixed point only when Kp=1. We check the closed-form.
        """
        kp = 10.0
        expected = kp / (1.0 + (kp - 1.0) * 1.0)
        assert apply_berezhkovskiy_correction(
            kp=kp, fu_p=1.0
        ) == pytest.approx(expected, abs=TOL)

    def test_kp_equals_one_is_identity(self):
        """Kp=1 is the true fixed point of the BZ formula."""
        for fu_p in (0.0, 0.1, 0.5, 1.0):
            assert apply_berezhkovskiy_correction(
                kp=1.0, fu_p=fu_p
            ) == pytest.approx(1.0, abs=TOL)

    def test_invalid_kp_raises(self):
        with pytest.raises(ValueError, match="kp"):
            apply_berezhkovskiy_correction(kp=-1.0, fu_p=0.5)

    def test_invalid_fu_p_raises(self):
        with pytest.raises(ValueError, match="fu_p"):
            apply_berezhkovskiy_correction(kp=10.0, fu_p=1.5)

    def test_nonfinite_kp_raises(self):
        with pytest.raises(ValueError, match="kp"):
            apply_berezhkovskiy_correction(kp=float("nan"), fu_p=0.5)


class TestInputValidation:
    def test_invalid_compound_type_raises(self):
        with pytest.raises(ValueError, match="compound_type"):
            compute_kp_rodgers_rowland(
                logp=2.0,
                pka=5.0,
                compound_type="nonsense",
                tissue_comp=CHARON_TISSUES["liver"],
                plasma_comp=CHARON_PLASMA,
            )

    def test_nonfinite_logp_raises(self):
        with pytest.raises(ValueError, match="logp"):
            compute_kp_rodgers_rowland(
                logp=float("inf"),
                pka=None,
                compound_type="neutral",
                tissue_comp=CHARON_TISSUES["liver"],
                plasma_comp=CHARON_PLASMA,
            )


class TestComputeAllKp:
    def test_returns_all_tissues(self):
        results = compute_all_kp(
            logp=2.7,
            pka=5.05,
            compound_type="acid",
            tissue_compositions=CHARON_TISSUES,
            plasma_composition=CHARON_PLASMA,
            method="rodgers_rowland",
        )
        assert set(results.keys()) == set(CHARON_TISSUES.keys())
        for name, kp in results.items():
            assert KP_MIN <= kp <= KP_MAX, f"{name}: {kp} out of range"

    def test_poulin_theil_method(self):
        results = compute_all_kp(
            logp=3.1,
            pka=6.15,
            compound_type="base",
            tissue_compositions=CHARON_TISSUES,
            plasma_composition=CHARON_PLASMA,
            method="poulin_theil",
        )
        assert len(results) == len(CHARON_TISSUES)

    def test_berezhkovskiy_requires_fu_p(self):
        with pytest.raises(ValueError, match="fu_p"):
            compute_all_kp(
                logp=2.0,
                pka=None,
                compound_type="neutral",
                tissue_compositions=CHARON_TISSUES,
                plasma_composition=CHARON_PLASMA,
                method="berezhkovskiy",
                fu_p=None,
            )

    def test_berezhkovskiy_reduces_basic_drug(self):
        """For a basic drug with low fu_p the BZ correction lowers Kp."""
        rr = compute_all_kp(
            logp=3.48,
            pka=9.45,
            compound_type="base",
            tissue_compositions=CHARON_TISSUES,
            plasma_composition=CHARON_PLASMA,
            method="rodgers_rowland",
        )
        bz = compute_all_kp(
            logp=3.48,
            pka=9.45,
            compound_type="base",
            tissue_compositions=CHARON_TISSUES,
            plasma_composition=CHARON_PLASMA,
            method="berezhkovskiy",
            fu_p=0.1,
        )
        # At least one tissue is strictly reduced.
        reduced = sum(1 for t in rr if bz[t] < rr[t] - 1e-12)
        assert reduced >= 1

    def test_invalid_method_raises(self):
        with pytest.raises(ValueError, match="method"):
            compute_all_kp(
                logp=2.0,
                pka=None,
                compound_type="neutral",
                tissue_compositions=CHARON_TISSUES,
                plasma_composition=CHARON_PLASMA,
                method="not_a_method",
            )
