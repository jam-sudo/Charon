"""Unit tests for the log-space conformal predictor."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from charon.predict.admet_ensemble import ADMETPredictor
from charon.predict.conformal import (
    CoverageReport,
    ConformalPredictor,
)

CAL_PATH = Path("/home/jam/Charon/data/validation/adme_reference.csv")


@pytest.fixture(scope="module")
def predictor() -> ADMETPredictor:
    return ADMETPredictor()


@pytest.fixture(scope="module")
def calibrated_conformal(predictor: ADMETPredictor) -> ConformalPredictor:
    """Module-scoped: calibrate once and share across tests."""
    cp = ConformalPredictor(CAL_PATH)
    cp.calibrate(predictor)
    return cp


class TestConstruction:
    def test_load_with_valid_path(self):
        cp = ConformalPredictor(CAL_PATH)
        assert cp is not None

    def test_accepts_string_path(self):
        cp = ConformalPredictor(str(CAL_PATH))
        assert cp is not None

    def test_invalid_coverage_low_raises(self):
        with pytest.raises(ValueError, match="coverage"):
            ConformalPredictor(CAL_PATH, coverage=0.0)

    def test_invalid_coverage_high_raises(self):
        with pytest.raises(ValueError, match="coverage"):
            ConformalPredictor(CAL_PATH, coverage=1.0)

    def test_invalid_coverage_negative_raises(self):
        with pytest.raises(ValueError, match="coverage"):
            ConformalPredictor(CAL_PATH, coverage=-0.1)


class TestCalibration:
    def test_calibrate_returns_report_for_fup(
        self, predictor: ADMETPredictor
    ):
        cp = ConformalPredictor(CAL_PATH)
        reports = cp.calibrate(predictor)
        assert "fup" in reports
        assert isinstance(reports["fup"], CoverageReport)

    def test_fup_coverage_above_threshold(
        self, calibrated_conformal: ConformalPredictor
    ):
        """Empirical fup coverage should be ≥ 85% (target 90%)."""
        report = calibrated_conformal.coverage_report()["fup"]
        assert report.empirical_coverage >= 0.85

    def test_calibration_populates_report_fields(
        self, calibrated_conformal: ConformalPredictor
    ):
        report = calibrated_conformal.coverage_report()["fup"]
        assert report.n_samples > 0
        assert report.quantile_log10 > 0.0
        assert report.factor > 1.0
        assert report.median_fold_error >= 1.0
        assert report.mean_fold_error >= 1.0

    def test_is_calibrated_true_after_calibrate(
        self, calibrated_conformal: ConformalPredictor
    ):
        assert calibrated_conformal.is_calibrated("fup") is True

    def test_is_calibrated_false_for_unknown(
        self, calibrated_conformal: ConformalPredictor
    ):
        assert calibrated_conformal.is_calibrated("nonexistent") is False

    def test_missing_calibration_file_raises(self, tmp_path: Path):
        bad_path = tmp_path / "not_a_file.csv"
        cp = ConformalPredictor(bad_path)
        with pytest.raises(FileNotFoundError):
            cp.calibrate(ADMETPredictor())


class TestGetInterval:
    def test_fup_interval_shape(
        self, calibrated_conformal: ConformalPredictor
    ):
        lo, hi = calibrated_conformal.get_interval("fup", 0.25)
        assert lo > 0.0
        assert hi <= 1.0
        assert lo <= 0.25 <= hi

    def test_fup_upper_clipped_at_one(
        self, calibrated_conformal: ConformalPredictor
    ):
        """Upper bound for fup must never exceed 1.0."""
        lo, hi = calibrated_conformal.get_interval("fup", 1.0)
        assert hi == pytest.approx(1.0)

    def test_unknown_property_raises_runtime(
        self, calibrated_conformal: ConformalPredictor
    ):
        with pytest.raises(RuntimeError, match="not calibrated"):
            calibrated_conformal.get_interval("unknown_prop", 0.5)

    def test_negative_prediction_raises(
        self, calibrated_conformal: ConformalPredictor
    ):
        with pytest.raises(ValueError, match="positive"):
            calibrated_conformal.get_interval("fup", -0.1)

    def test_zero_prediction_raises(
        self, calibrated_conformal: ConformalPredictor
    ):
        with pytest.raises(ValueError, match="positive"):
            calibrated_conformal.get_interval("fup", 0.0)

    def test_interval_brackets_point_estimate(
        self, calibrated_conformal: ConformalPredictor
    ):
        """Even with physical clipping, lo ≤ pred ≤ hi must hold."""
        for pred in (0.01, 0.05, 0.2, 0.5, 0.95):
            lo, hi = calibrated_conformal.get_interval("fup", pred)
            assert lo <= pred <= hi


class TestCalibrateFromOof:
    def test_calibrate_clint_from_oof(self):
        cp = ConformalPredictor(CAL_PATH)
        residuals = np.array([0.1, 0.2, 0.3, 0.15, 0.25, 0.18, 0.22, 0.4])
        report = cp.calibrate_from_oof("clint", residuals)
        assert isinstance(report, CoverageReport)
        assert report.property_name == "clint"
        assert report.n_samples == 8
        assert cp.is_calibrated("clint")

    def test_oof_empty_raises(self):
        cp = ConformalPredictor(CAL_PATH)
        with pytest.raises(ValueError, match="residuals"):
            cp.calibrate_from_oof("clint", np.array([]))

    def test_oof_all_nonfinite_raises(self):
        cp = ConformalPredictor(CAL_PATH)
        with pytest.raises(ValueError, match="residuals"):
            cp.calibrate_from_oof(
                "clint", np.array([float("nan"), float("inf")])
            )

    def test_oof_accepts_signed_residuals(self):
        """The method takes abs() internally, so signed input is fine."""
        cp = ConformalPredictor(CAL_PATH)
        residuals = np.array([-0.1, 0.2, -0.3, 0.15])
        report = cp.calibrate_from_oof("clint", residuals)
        assert report.quantile_log10 > 0.0


class TestCoverageReport:
    def test_coverage_report_returns_dict(
        self, calibrated_conformal: ConformalPredictor
    ):
        rep = calibrated_conformal.coverage_report()
        assert isinstance(rep, dict)
        assert "fup" in rep

    def test_coverage_report_is_snapshot(
        self, calibrated_conformal: ConformalPredictor
    ):
        """Mutating the returned dict should not affect the predictor."""
        rep = dict(calibrated_conformal.coverage_report())
        rep.clear()
        assert calibrated_conformal.is_calibrated("fup") is True
