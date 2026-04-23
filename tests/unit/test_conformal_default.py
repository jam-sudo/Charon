"""Tests for ConformalPredictor.load_default (cache + hash invalidation)."""
from __future__ import annotations

import json
from pathlib import Path
from unittest.mock import patch

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
        assert "coverage_target" in data
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

    def test_force_recalibrate_ignores_valid_cache(self, tmp_path):
        from charon.predict.conformal import ConformalPredictor
        cache = tmp_path / "conformal_cache.json"
        ConformalPredictor.load_default(cache_path=cache)
        # Cache is valid. force_recalibrate=True should rebuild it.
        mtime_before = cache.stat().st_mtime
        import time; time.sleep(0.05)  # ensure mtime tick
        cp = ConformalPredictor.load_default(cache_path=cache, force_recalibrate=True)
        mtime_after = cache.stat().st_mtime
        assert mtime_after > mtime_before, "force_recalibrate must rewrite cache"
        assert cp.is_calibrated("fup")
        assert cp.is_calibrated("clint_hepatocyte")

    def test_coverage_mismatch_recalibrates(self, tmp_path, caplog):
        from charon.predict.conformal import ConformalPredictor
        cache = tmp_path / "conformal_cache.json"
        ConformalPredictor.load_default(cache_path=cache)

        import json as _json
        data = _json.loads(cache.read_text())
        data["coverage_target"] = 0.80  # was 0.90
        cache.write_text(_json.dumps(data))

        # Requesting default 0.90 coverage should reject 0.80 cache.
        import logging
        with caplog.at_level(logging.WARNING):
            cp = ConformalPredictor.load_default(cache_path=cache)
        assert any("coverage" in rec.message.lower() for rec in caplog.records)
        fresh = _json.loads(cache.read_text())
        assert fresh["coverage_target"] == 0.90

    def test_missing_clint_oof_logs_warning_fup_only(self, tmp_path, caplog):
        from charon.predict.conformal import ConformalPredictor
        import logging
        cache = tmp_path / "conformal_cache.json"
        nonexistent = tmp_path / "nope.npy"
        with caplog.at_level(logging.WARNING):
            cp = ConformalPredictor.load_default(
                cache_path=cache, clint_oof_path=nonexistent
            )
        assert cp.is_calibrated("fup")
        assert not cp.is_calibrated("clint_hepatocyte")
        assert any("clint" in rec.message.lower() and "missing" in rec.message.lower()
                   for rec in caplog.records)


class TestResetSingleton:
    def test_reset_clears_module_singleton(self, tmp_path):
        from charon.predict.conformal import reset_default_conformal
        cache = tmp_path / "conformal_cache.json"
        first = ConformalPredictor.load_default(cache_path=cache)
        reset_default_conformal()
        second = ConformalPredictor.load_default(cache_path=cache)
        assert first is not second
