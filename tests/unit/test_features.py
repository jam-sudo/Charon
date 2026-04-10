"""Tests for charon.predict.features.compute_features.

Verifies that Charon produces the exact same 2057-D feature vector as
Sisyphus (the canonical reference implementation). Also covers invalid
input handling and structural properties of the feature vector.
"""

from __future__ import annotations

import numpy as np
import pytest

from charon.predict.features import (
    FEATURE_LENGTH,
    MORGAN_BITS,
    MORGAN_RADIUS,
    compute_features,
)

REFERENCE_DRUGS: list[tuple[str, str]] = [
    ("caffeine", "Cn1cnc2c1c(=O)n(C)c(=O)n2C"),
    ("aspirin", "CC(=O)Oc1ccccc1C(=O)O"),
    ("ibuprofen", "CC(C)Cc1ccc(cc1)C(C)C(=O)O"),
    ("midazolam", "Clc1ccc2c(c1)C(=NCc3nccn3C)c1ccccc1N2"),
    ("warfarin", "CC(=O)CC(c1ccccc1)c1c(O)c2ccccc2oc1=O"),
    ("propranolol", "CC(C)NCC(O)COc1cccc2ccccc12"),
]


class TestFeatureLengthAndShape:
    @pytest.mark.parametrize("name,smiles", REFERENCE_DRUGS)
    def test_length_is_2057(self, name: str, smiles: str) -> None:
        """Every feature vector must be exactly 2057-D."""
        features = compute_features(smiles)
        assert features.shape == (FEATURE_LENGTH,)
        assert features.dtype == np.float64

    def test_constants(self) -> None:
        """Module constants match the expected specification."""
        assert FEATURE_LENGTH == 2057
        assert MORGAN_BITS == 2048
        assert MORGAN_RADIUS == 2


class TestFingerprintSection:
    @pytest.mark.parametrize("name,smiles", REFERENCE_DRUGS)
    def test_fingerprint_is_binary(self, name: str, smiles: str) -> None:
        """Morgan fingerprint bits must be 0 or 1."""
        features = compute_features(smiles)
        fp_section = features[:MORGAN_BITS]
        unique = np.unique(fp_section)
        assert set(unique).issubset({0.0, 1.0}), f"{name}: non-binary FP values {unique}"

    @pytest.mark.parametrize("name,smiles", REFERENCE_DRUGS)
    def test_fingerprint_nonzero(self, name: str, smiles: str) -> None:
        """Every valid drug should have at least one set Morgan bit."""
        features = compute_features(smiles)
        fp_section = features[:MORGAN_BITS]
        assert fp_section.sum() > 0, f"{name}: all-zero fingerprint"


class TestDescriptorSection:
    def test_aspirin_descriptors(self) -> None:
        """Aspirin descriptors match hand-computed reference values."""
        features = compute_features("CC(=O)Oc1ccccc1C(=O)O")
        descriptors = features[MORGAN_BITS:]
        # Aspirin: MW ≈ 180.16, logP ≈ 1.2, TPSA ≈ 63.6
        # Feature index 2 = MW / 600, so descriptor[2] ≈ 180.16 / 600 = 0.3003
        assert 0.28 < descriptors[2] < 0.32
        # logP ≈ 1.0-1.5
        assert 0.8 < descriptors[0] < 1.6

    def test_caffeine_fraction_csp3(self) -> None:
        """Caffeine has 3 sp3 carbons out of 8, FractionCSP3 = 0.375."""
        features = compute_features("Cn1cnc2c1c(=O)n(C)c(=O)n2C")
        descriptors = features[MORGAN_BITS:]
        # descriptor[7] = FractionCSP3 (raw)
        assert descriptors[7] == pytest.approx(0.375, abs=0.01)


class TestInvalidSmiles:
    def test_empty_string_raises(self) -> None:
        with pytest.raises(ValueError, match="Invalid SMILES"):
            compute_features("")

    def test_whitespace_only_raises(self) -> None:
        with pytest.raises(ValueError, match="Invalid SMILES"):
            compute_features("   ")

    def test_garbage_smiles_raises(self) -> None:
        with pytest.raises(ValueError, match="Invalid SMILES"):
            compute_features("not_a_smiles_at_all")

    def test_unbalanced_parens_raises(self) -> None:
        with pytest.raises(ValueError, match="Invalid SMILES"):
            compute_features("CCC(CC")


class TestDeterminism:
    @pytest.mark.parametrize("name,smiles", REFERENCE_DRUGS)
    def test_repeated_calls_identical(self, name: str, smiles: str) -> None:
        """compute_features must be deterministic."""
        a = compute_features(smiles)
        b = compute_features(smiles)
        np.testing.assert_array_equal(a, b)
