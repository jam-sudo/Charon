"""XGBoost ADME property predictor.

Loads Charon-trained XGBoost models (``models/xgboost_fup.json`` and
``models/xgboost_clint.json``) and predicts ADME properties from a
SMILES string using the shared 2057-D feature vector.

Only fup and CLint are predicted here. B:P ratio is handled separately
by :mod:`charon.predict.bp_ratio` (which uses an empirical physiological
formula because the available RBP training data is too small and
coincides with the conformal calibration set).
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from threading import Lock
from typing import Any

import numpy as np
import xgboost as xgb

from charon.predict.features import FEATURE_LENGTH, compute_features

logger = logging.getLogger(__name__)

# Physical clipping bounds for back-transformed predictions.
FUP_MIN, FUP_MAX = 0.001, 1.0
CLINT_MIN, CLINT_MAX = 0.1, 1000.0  # uL/min/10^6 cells


def _default_models_dir() -> Path:
    """Return ``<repo_root>/models`` regardless of installation layout.

    ``__file__`` is ``<repo>/src/charon/predict/admet_ensemble.py``, so
    ``parents[3]`` lands on the repository root.
    """
    return Path(__file__).resolve().parents[3] / "models"


@dataclass(frozen=True)
class ADMEPrediction:
    """Raw ML predictions with both linear and log-space views.

    The log-space values are preserved so that the conformal predictor
    can compute symmetric log-space intervals without losing numerical
    precision from the back-transform.
    """

    smiles: str
    fup: float             # fraction, clipped to [FUP_MIN, FUP_MAX]
    fup_log10: float       # raw XGBoost output: log10(fup)
    clint_hepatocyte: float   # uL/min/10^6 cells, clipped to [CLINT_MIN, CLINT_MAX]
    clint_log10: float     # raw XGBoost output: log10(clint)


class ADMETPredictor:
    """Thread-safe lazy XGBoost ADME predictor.

    Models are loaded on first use and cached. The predictor is cheap to
    instantiate (no disk I/O in ``__init__``) so it can be created per
    pipeline run and still share its cache through
    :attr:`_model_cache` if the same ``models_dir`` is re-used.
    """

    _cache_lock: Lock = Lock()
    _model_cache: dict[tuple[Path, str], xgb.XGBRegressor] = {}

    def __init__(self, models_dir: Path | None = None) -> None:
        self._models_dir = Path(models_dir) if models_dir is not None else _default_models_dir()
        self._fup_path = self._models_dir / "xgboost_fup.json"
        self._clint_path = self._models_dir / "xgboost_clint.json"

    # ------------------------------------------------------------------
    # Model loading (lazy, cached)
    # ------------------------------------------------------------------
    def _load_model(self, name: str, path: Path) -> xgb.XGBRegressor:
        key = (self._models_dir.resolve(), name)
        with ADMETPredictor._cache_lock:
            cached = ADMETPredictor._model_cache.get(key)
            if cached is not None:
                return cached
            if not path.exists():
                raise FileNotFoundError(
                    f"XGBoost model not found: {path}. Train it with "
                    f"scripts/train_{name.split('_')[-1]}.py."
                )
            model = xgb.XGBRegressor()
            model.load_model(str(path))
            # Defensive check: models must expect the canonical feature
            # length. Mismatches indicate silent data pipeline drift.
            n_features = getattr(model, "n_features_in_", None)
            if n_features is not None and n_features != FEATURE_LENGTH:
                raise ValueError(
                    f"Model {path.name} expects {n_features} features, "
                    f"but charon.predict.features produces {FEATURE_LENGTH}."
                )
            ADMETPredictor._model_cache[key] = model
            logger.debug("Loaded ADME model %s from %s", name, path)
            return model

    @property
    def fup_model(self) -> xgb.XGBRegressor:
        return self._load_model("xgboost_fup", self._fup_path)

    @property
    def clint_model(self) -> xgb.XGBRegressor:
        return self._load_model("xgboost_clint", self._clint_path)

    # ------------------------------------------------------------------
    # Prediction
    # ------------------------------------------------------------------
    def predict(self, smiles: str) -> ADMEPrediction:
        """Predict fup and CLint for a single SMILES.

        Args:
            smiles: Input SMILES (canonical or raw).

        Returns:
            :class:`ADMEPrediction` with both clipped linear values and
            raw log10 predictions.

        Raises:
            ValueError: If ``smiles`` is invalid.
            FileNotFoundError: If the required model files are missing.
        """
        features = compute_features(smiles)
        X = features.reshape(1, -1)

        fup_log10_raw = float(self.fup_model.predict(X)[0])
        fup_linear = float(np.clip(10.0 ** fup_log10_raw, FUP_MIN, FUP_MAX))

        clint_log10_raw = float(self.clint_model.predict(X)[0])
        clint_linear = float(np.clip(10.0 ** clint_log10_raw, CLINT_MIN, CLINT_MAX))

        return ADMEPrediction(
            smiles=smiles,
            fup=fup_linear,
            fup_log10=fup_log10_raw,
            clint_hepatocyte=clint_linear,
            clint_log10=clint_log10_raw,
        )

    # ------------------------------------------------------------------
    # Testability helpers
    # ------------------------------------------------------------------
    @classmethod
    def _clear_cache(cls) -> None:
        """Drop all cached models (test fixtures only)."""
        with cls._cache_lock:
            cls._model_cache.clear()

    def metadata(self) -> dict[str, Any]:
        """Return paths and existence state for introspection."""
        return {
            "models_dir": str(self._models_dir),
            "fup_path": str(self._fup_path),
            "fup_exists": self._fup_path.exists(),
            "clint_path": str(self._clint_path),
            "clint_exists": self._clint_path.exists(),
        }
