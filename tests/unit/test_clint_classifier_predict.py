"""Tests for ClintClassifier.predict_proba."""
from __future__ import annotations

import math

import pytest

from charon.predict.clint_classifier import (
    BUCKET_CENTERS,
    BUCKET_RANGES,
    ClintClassifier,
)


CAFFEINE = "Cn1cnc2n(C)c(=O)n(C)c(=O)c12"
THEOPHYLLINE = "Cn1c(=O)c2[nH]cnc2n(C)c1=O"


class TestBucketConstants:
    def test_bucket_ranges_cover_clint_domain(self):
        assert BUCKET_RANGES["low"] == (0.1, 10.0)
        assert BUCKET_RANGES["med"] == (10.0, 50.0)
        assert BUCKET_RANGES["high"] == (50.0, 1000.0)

    def test_centers_are_log_geometric_means(self):
        for name, (lo, hi) in BUCKET_RANGES.items():
            expected = 10.0 ** ((math.log10(lo) + math.log10(hi)) / 2)
            assert BUCKET_CENTERS[name] == pytest.approx(expected, rel=0.35), (
                f"{name}: center {BUCKET_CENTERS[name]} vs log-mean {expected:.2f}"
            )


class TestPredictProba:
    @pytest.fixture(scope="class")
    def clf(self):
        return ClintClassifier.load_default()

    def test_returns_three_key_dict_summing_to_one(self, clf):
        probs = clf.predict_proba(CAFFEINE)
        assert set(probs) == {"low", "med", "high"}
        assert all(0.0 <= v <= 1.0 for v in probs.values())
        assert abs(sum(probs.values()) - 1.0) < 1e-5

    def test_theophylline_prefers_low_bucket(self, clf):
        """Theophylline has low observed CLint (~5 uL/min/10^6) -> Low bucket."""
        probs = clf.predict_proba(THEOPHYLLINE)
        assert probs["low"] > probs["med"]
        assert probs["low"] > probs["high"]

    def test_invalid_smiles_raises(self, clf):
        with pytest.raises(ValueError):
            clf.predict_proba("not-a-smiles-!!!")
