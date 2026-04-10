"""Unit tests for the XGBoost ADME ensemble predictor.

Uses the trained models shipped at ``<repo>/models/``.
"""

from __future__ import annotations

import math
import time
from pathlib import Path

import pytest

from charon.predict.admet_ensemble import (
    ADMEPrediction,
    ADMETPredictor,
    CLINT_MAX,
    CLINT_MIN,
    FUP_MAX,
    FUP_MIN,
)
from charon.predict.features import FEATURE_LENGTH

# Reference SMILES used throughout the tests.
ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"
PROPRANOLOL = "CC(C)NCC(O)COc1cccc2ccccc12"
CAFFEINE = "Cn1cnc2c1c(=O)n(C)c(=O)n2C"


@pytest.fixture(scope="module")
def predictor() -> ADMETPredictor:
    """Module-scoped predictor — models are loaded once per session."""
    return ADMETPredictor()


class TestADMETPredictor:
    """Loads real trained models from <repo>/models/."""

    def test_predict_aspirin_returns_prediction(self, predictor: ADMETPredictor):
        result = predictor.predict(ASPIRIN)
        assert isinstance(result, ADMEPrediction)
        assert result.smiles == ASPIRIN

    def test_fup_in_physical_bounds(self, predictor: ADMETPredictor):
        result = predictor.predict(ASPIRIN)
        assert FUP_MIN <= result.fup <= FUP_MAX

    def test_clint_in_physical_bounds(self, predictor: ADMETPredictor):
        result = predictor.predict(ASPIRIN)
        assert CLINT_MIN <= result.clint_hepatocyte <= CLINT_MAX

    def test_log10_outputs_are_finite(self, predictor: ADMETPredictor):
        result = predictor.predict(ASPIRIN)
        assert isinstance(result.fup_log10, float)
        assert isinstance(result.clint_log10, float)
        assert math.isfinite(result.fup_log10)
        assert math.isfinite(result.clint_log10)

    def test_invalid_smiles_raises(self, predictor: ADMETPredictor):
        with pytest.raises(ValueError):
            predictor.predict("not_a_valid_smiles_xyz_987")

    def test_empty_smiles_raises(self, predictor: ADMETPredictor):
        with pytest.raises(ValueError):
            predictor.predict("")

    def test_missing_model_dir_raises(self, tmp_path: Path):
        """A models_dir with no model files raises FileNotFoundError."""
        missing = tmp_path / "does_not_exist"
        bad_predictor = ADMETPredictor(models_dir=missing)
        with pytest.raises(FileNotFoundError):
            bad_predictor.predict(ASPIRIN)

    def test_cache_second_call_fast(self, predictor: ADMETPredictor):
        """Second prediction should be quick (models already cached)."""
        # Warm-up: load models.
        predictor.predict(ASPIRIN)
        t0 = time.perf_counter()
        predictor.predict(ASPIRIN)
        elapsed = time.perf_counter() - t0
        # No hard wall-time guarantee but a cached prediction should be
        # well under half a second on any reasonable machine.
        assert elapsed < 1.0

    def test_models_feature_length_matches(self, predictor: ADMETPredictor):
        """Loaded models must expect exactly FEATURE_LENGTH (=2057) features."""
        fup_model = predictor.fup_model
        clint_model = predictor.clint_model
        assert getattr(fup_model, "n_features_in_", FEATURE_LENGTH) == FEATURE_LENGTH
        assert getattr(clint_model, "n_features_in_", FEATURE_LENGTH) == FEATURE_LENGTH

    def test_deterministic_same_smiles(self, predictor: ADMETPredictor):
        """Same SMILES → identical prediction (XGBoost is deterministic)."""
        r1 = predictor.predict(PROPRANOLOL)
        r2 = predictor.predict(PROPRANOLOL)
        assert r1.fup == r2.fup
        assert r1.clint_hepatocyte == r2.clint_hepatocyte
        assert r1.fup_log10 == r2.fup_log10
        assert r1.clint_log10 == r2.clint_log10

    def test_metadata_reports_paths(self, predictor: ADMETPredictor):
        meta = predictor.metadata()
        assert "models_dir" in meta
        assert "fup_path" in meta
        assert "clint_path" in meta
        assert meta["fup_exists"] is True
        assert meta["clint_exists"] is True

    def test_different_smiles_differ(self, predictor: ADMETPredictor):
        """Different compounds should give different predictions."""
        a = predictor.predict(ASPIRIN)
        c = predictor.predict(CAFFEINE)
        # Not identical across two unrelated molecules.
        assert (a.fup, a.clint_hepatocyte) != (c.fup, c.clint_hepatocyte)
