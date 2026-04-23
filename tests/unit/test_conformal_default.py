"""Tests for ConformalPredictor.load_default (cache + hash invalidation)."""
from __future__ import annotations

import json
from pathlib import Path
from unittest.mock import patch

import pytest

from charon.predict.conformal import ConformalPredictor


REPO_ROOT = Path(__file__).resolve().parents[2]


class TestLoadDefaultFreshCache:
    def test_first_call_creates_cache(self, tmp_path):
        cache = tmp_path / "conformal_cache.json"
        assert not cache.exists()
        cp = ConformalPredictor.load_default(cache_path=cache)
        assert cache.exists()
        assert cp.is_calibrated("fup")
        assert cp.is_calibrated("clint_hepatocyte")

    def test_cache_schema_has_required_fields(self, tmp_path):
        cache = tmp_path / "conformal_cache.json"
        ConformalPredictor.load_default(cache_path=cache)
        data = json.loads(cache.read_text())
        assert data["schema_version"] == 1
        assert "generated_utc" in data
        assert "source_hashes" in data
        assert set(data["reports"]) >= {"fup", "clint_hepatocyte"}
        for rpt in data["reports"].values():
            for key in ("n_samples", "quantile_log10", "factor",
                        "empirical_coverage", "median_fold_error",
                        "mean_fold_error"):
                assert key in rpt


class TestCacheReuse:
    def test_second_call_loads_cache_without_recalibration(self, tmp_path):
        cache = tmp_path / "conformal_cache.json"
        ConformalPredictor.load_default(cache_path=cache)
        with patch.object(
            ConformalPredictor, "calibrate",
            side_effect=AssertionError("calibrate() must not run on cached load"),
        ):
            with patch.object(
                ConformalPredictor, "calibrate_from_oof",
                side_effect=AssertionError("calibrate_from_oof() must not run on cached load"),
            ):
                cp = ConformalPredictor.load_default(cache_path=cache)
        assert cp.is_calibrated("fup")
        assert cp.is_calibrated("clint_hepatocyte")


class TestCacheInvalidation:
    def test_source_hash_change_recalibrates(self, tmp_path):
        cache = tmp_path / "conformal_cache.json"
        ConformalPredictor.load_default(cache_path=cache)

        data = json.loads(cache.read_text())
        data["source_hashes"]["adme_reference.csv"] = "sha256:deadbeef"
        cache.write_text(json.dumps(data))

        cp = ConformalPredictor.load_default(cache_path=cache)
        fresh = json.loads(cache.read_text())
        assert fresh["source_hashes"]["adme_reference.csv"] != "sha256:deadbeef", (
            "Expected cache to be refreshed with a real hash after mismatch"
        )

    def test_corrupted_cache_falls_back(self, tmp_path, caplog):
        cache = tmp_path / "conformal_cache.json"
        cache.write_text("{not json")
        cp = ConformalPredictor.load_default(cache_path=cache)
        assert cp.is_calibrated("fup")
        assert "falling back" in caplog.text.lower()


class TestResetSingleton:
    def test_reset_clears_module_singleton(self, tmp_path):
        from charon.predict.conformal import reset_default_conformal
        cache = tmp_path / "conformal_cache.json"
        first = ConformalPredictor.load_default(cache_path=cache)
        reset_default_conformal()
        second = ConformalPredictor.load_default(cache_path=cache)
        assert first is not second
