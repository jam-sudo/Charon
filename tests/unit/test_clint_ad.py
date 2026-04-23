"""Tests for ClintLocalAD (per-compound Tanimoto applicability domain)."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from charon.predict.clint_ad import (
    AD_LOW_THRESHOLD,
    AD_MODERATE_THRESHOLD,
    ClintLocalAD,
)


CAFFEINE = "Cn1cnc2n(C)c(=O)n(C)c(=O)c12"
GLUCOSE = "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O"


class TestThresholds:
    def test_thresholds_are_strict(self):
        assert AD_LOW_THRESHOLD < AD_MODERATE_THRESHOLD

    def test_classify_boundaries(self):
        ad = ClintLocalAD(reference_fps=[])  # empty ref for pure classify test
        assert ad.classify(0.6) == "HIGH"
        assert ad.classify(0.5) == "HIGH"         # inclusive at MODERATE threshold
        assert ad.classify(0.4) == "MODERATE"
        assert ad.classify(0.3) == "MODERATE"     # inclusive at LOW threshold
        assert ad.classify(0.2999) == "LOW"
        assert ad.classify(0.0) == "LOW"
        assert ad.classify(float("nan")) == "INVALID"


class TestLoadDefault:
    def test_first_call_creates_cache(self, tmp_path):
        cache = tmp_path / "clint_ad_fingerprints.npz"
        assert not cache.exists()
        ad = ClintLocalAD.load_default(cache_path=cache)
        assert cache.exists()
        assert len(ad._reference_fps) > 1000, (
            f"Expect >1000 training compounds; got {len(ad._reference_fps)}"
        )

    def test_second_call_uses_cache_without_rebuilding(self, tmp_path):
        cache = tmp_path / "clint_ad_fingerprints.npz"
        ad1 = ClintLocalAD.load_default(cache_path=cache)
        mtime1 = cache.stat().st_mtime
        ad2 = ClintLocalAD.load_default(cache_path=cache)
        mtime2 = cache.stat().st_mtime
        assert mtime1 == mtime2, "Cache must not be rewritten on hit"
        assert len(ad1._reference_fps) == len(ad2._reference_fps)

    def test_source_hash_change_rebuilds(self, tmp_path):
        cache = tmp_path / "clint_ad_fingerprints.npz"
        ClintLocalAD.load_default(cache_path=cache)
        # Corrupt the stored source hash
        data = dict(np.load(cache, allow_pickle=True).items())
        data["source_hash"] = np.array("sha256:deadbeef")
        np.savez(cache, **data)
        ad = ClintLocalAD.load_default(cache_path=cache)
        fresh = dict(np.load(cache, allow_pickle=True).items())
        assert str(fresh["source_hash"]) != "sha256:deadbeef"


class TestSimilarity:
    @pytest.fixture(scope="class")
    def ad(self):
        return ClintLocalAD.load_default()

    def test_caffeine_has_high_similarity(self, ad):
        sim = ad.max_similarity(CAFFEINE)
        assert 0.3 <= sim <= 1.0

    def test_glucose_outside_or_boundary(self, ad):
        sim = ad.max_similarity(GLUCOSE)
        assert 0.0 <= sim <= 0.6, "glucose should not be deep in the drug AD"

    def test_invalid_smiles_returns_nan(self, ad):
        sim = ad.max_similarity("not_a_smiles_123#$%")
        assert np.isnan(sim)
