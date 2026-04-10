"""Unit tests for the rule-based SMARTS pKa predictor."""

from __future__ import annotations

import pytest

from charon.predict.pka import PKaResult, predict_pka


class TestPredictPka:
    """End-to-end behaviour of ``predict_pka`` on reference compounds."""

    def test_aspirin_is_acid(self):
        """Aspirin has a carboxylic acid → acid, pka_acid = 4.5."""
        result = predict_pka("CC(=O)Oc1ccccc1C(=O)O")
        assert isinstance(result, PKaResult)
        assert result.compound_type == "acid"
        assert result.pka_acid == pytest.approx(4.5)
        assert "carboxylic_acid" in result.detected_acid_groups

    def test_propranolol_is_base(self):
        """Propranolol has a secondary amine → base, pka_base = 9.5."""
        result = predict_pka("CC(C)NCC(O)COc1cccc2ccccc12")
        assert result.compound_type == "base"
        assert result.pka_base == pytest.approx(9.5)

    def test_caffeine_is_neutral(self):
        """Caffeine has no ionizable groups detected by our patterns."""
        result = predict_pka("Cn1cnc2c1c(=O)n(C)c(=O)n2C")
        assert result.compound_type == "neutral"

    def test_metformin_base_or_neutral(self):
        """Metformin: guanidine not in our SMARTS list.

        Either ``base`` (if a secondary/primary amine is detected) or
        ``neutral`` (if nothing matches) are acceptable outcomes — we do
        not claim a guanidine pKa from this rule-based predictor.
        """
        result = predict_pka("CN(C)C(=N)NC(=N)N")
        assert result.compound_type in ("base", "neutral")

    def test_invalid_smiles_raises(self):
        with pytest.raises(ValueError, match="Invalid SMILES"):
            predict_pka("not_a_valid_smiles_xyz_123")

    def test_empty_string_raises(self):
        with pytest.raises(ValueError, match="Invalid SMILES"):
            predict_pka("")

    def test_whitespace_only_raises(self):
        with pytest.raises(ValueError, match="Invalid SMILES"):
            predict_pka("   ")

    def test_pkaresult_has_all_required_fields(self):
        """PKaResult exposes the documented attribute surface."""
        result = predict_pka("CC(=O)Oc1ccccc1C(=O)O")
        assert hasattr(result, "pka_acid")
        assert hasattr(result, "pka_base")
        assert hasattr(result, "compound_type")
        assert hasattr(result, "detected_acid_groups")
        assert hasattr(result, "detected_base_groups")
        assert isinstance(result.detected_acid_groups, tuple)
        assert isinstance(result.detected_base_groups, tuple)

    def test_detected_groups_names_are_strings(self):
        result = predict_pka("CC(=O)Oc1ccccc1C(=O)O")
        for name in result.detected_acid_groups:
            assert isinstance(name, str)
        for name in result.detected_base_groups:
            assert isinstance(name, str)

    def test_compound_type_in_valid_set(self):
        """compound_type must always be one of the four valid classes."""
        result = predict_pka("CC(=O)Oc1ccccc1C(=O)O")
        assert result.compound_type in PKaResult.VALID_TYPES

    def test_frozen_dataclass(self):
        """PKaResult is frozen (immutable)."""
        result = predict_pka("CC(=O)Oc1ccccc1C(=O)O")
        with pytest.raises((AttributeError, Exception)):
            result.pka_acid = 99.0  # type: ignore[misc]
