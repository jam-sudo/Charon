"""Unit tests for the empirical blood:plasma ratio estimator."""

from __future__ import annotations

import pytest

from charon.predict.bp_ratio import BP_MAX, BP_MIN, HUMAN_HCT, predict_bp_ratio


class TestPredictBpRatio:
    """Poulin-Theil simplification: B:P = 1 + Hct × (fu_p × Kp_rbc − 1)."""

    def test_fu_p_unity_neutral_returns_one(self):
        """With fu_p=1 and Kp_rbc=1 (neutral) the formula collapses to 1."""
        bp = predict_bp_ratio(fu_p=1.0, compound_type="neutral")
        assert bp == pytest.approx(1.0)

    def test_fu_p_half_neutral_below_one(self):
        """With plasma binding present B:P dips below 1."""
        bp = predict_bp_ratio(fu_p=0.5, compound_type="neutral")
        assert bp < 1.0
        # Hand calc: 1 + 0.45 * (0.5 * 1.0 - 1) = 1 - 0.225 = 0.775
        assert bp == pytest.approx(1.0 + HUMAN_HCT * (0.5 * 1.0 - 1.0))

    def test_very_low_fu_p_clamped_above_bp_min(self):
        """Very low fu_p must not drive B:P below BP_MIN."""
        bp = predict_bp_ratio(fu_p=0.01, compound_type="neutral")
        assert bp >= BP_MIN

    def test_basic_compound_higher_than_neutral(self):
        """Basic drugs (default Kp_rbc=1.3) have slightly higher B:P."""
        bp_neutral = predict_bp_ratio(fu_p=0.5, compound_type="neutral")
        bp_base = predict_bp_ratio(fu_p=0.5, compound_type="base")
        assert bp_base > bp_neutral

    def test_acidic_compound_lower_than_neutral(self):
        """Acids (Kp_rbc=0.7) are disfavoured by erythrocytes."""
        bp_neutral = predict_bp_ratio(fu_p=0.5, compound_type="neutral")
        bp_acid = predict_bp_ratio(fu_p=0.5, compound_type="acid")
        assert bp_acid < bp_neutral

    @pytest.mark.parametrize("fu_p", [-0.01, -1.0, 1.01, 2.0])
    def test_fu_p_out_of_range_raises(self, fu_p):
        with pytest.raises(ValueError, match="fu_p"):
            predict_bp_ratio(fu_p=fu_p)

    def test_fu_p_nan_raises(self):
        with pytest.raises(ValueError, match="fu_p"):
            predict_bp_ratio(fu_p=float("nan"))

    @pytest.mark.parametrize("hct", [0.0, -0.1, 1.0, 1.5])
    def test_hct_out_of_range_raises(self, hct):
        with pytest.raises(ValueError, match="hct"):
            predict_bp_ratio(fu_p=0.5, hct=hct)

    @pytest.mark.parametrize(
        "fu_p,compound_type",
        [
            (0.0, "neutral"),
            (0.01, "acid"),
            (0.05, "base"),
            (0.3, "zwitterion"),
            (0.75, "neutral"),
            (1.0, "base"),
        ],
    )
    def test_output_in_physiological_range(self, fu_p, compound_type):
        """Output must always be in [BP_MIN, BP_MAX]."""
        bp = predict_bp_ratio(fu_p=fu_p, compound_type=compound_type)
        assert BP_MIN <= bp <= BP_MAX

    def test_kp_rbc_override(self):
        """Explicit kp_rbc overrides the compound-type default."""
        bp = predict_bp_ratio(fu_p=0.5, compound_type="base", kp_rbc=2.0)
        # 1 + 0.45 * (0.5 * 2.0 - 1.0) = 1.0
        assert bp == pytest.approx(1.0)

    def test_invalid_kp_rbc_raises(self):
        with pytest.raises(ValueError, match="kp_rbc"):
            predict_bp_ratio(fu_p=0.5, kp_rbc=0.0)
        with pytest.raises(ValueError, match="kp_rbc"):
            predict_bp_ratio(fu_p=0.5, kp_rbc=-1.0)

    def test_upper_clamp(self):
        """With a very large kp_rbc override the result clamps at BP_MAX."""
        bp = predict_bp_ratio(fu_p=1.0, compound_type="neutral", kp_rbc=100.0)
        assert bp == pytest.approx(BP_MAX)
