"""Loaded XGBoost classifier wrapper for CLint Tier 3."""
from __future__ import annotations

import logging
import math
from pathlib import Path

import xgboost as xgb

from charon.predict.features import compute_features

logger = logging.getLogger(__name__)

BUCKET_NAMES = ("low", "med", "high")
BUCKET_RANGES: dict[str, tuple[float, float]] = {
    "low":  (0.1, 10.0),
    "med":  (10.0, 50.0),
    "high": (50.0, 1000.0),
}


def _log_geometric_mean(lo: float, hi: float) -> float:
    return 10.0 ** ((math.log10(lo) + math.log10(hi)) / 2)


BUCKET_CENTERS: dict[str, float] = {
    name: round(_log_geometric_mean(lo, hi), 3)
    for name, (lo, hi) in BUCKET_RANGES.items()
}

_DEFAULT_MODEL_PATH = (
    Path(__file__).resolve().parents[3] / "models" / "xgboost_clint_classifier.json"
)


class ClintClassifier:
    """Wraps the persisted XGBoost 3-class classifier.

    Use ``ClintClassifier.load_default()`` to get a singleton-style instance
    pointing at ``models/xgboost_clint_classifier.json``.
    """

    def __init__(self, model: xgb.XGBClassifier) -> None:
        self._model = model

    @classmethod
    def load_default(cls, model_path: Path | None = None) -> "ClintClassifier":
        path = Path(model_path) if model_path else _DEFAULT_MODEL_PATH
        if not path.exists():
            raise FileNotFoundError(
                f"Classifier model missing at {path} — run "
                "`python scripts/train_clint_classifier.py` first."
            )
        model = xgb.XGBClassifier()
        model.load_model(str(path))
        return cls(model)

    def predict_proba(self, smiles: str) -> dict[str, float]:
        """Return probabilities over {low, med, high} CLint buckets.

        Args:
            smiles: Canonical or raw SMILES.

        Returns:
            Dict with keys ``"low"``, ``"med"``, ``"high"``; values sum to 1.

        Raises:
            ValueError: If ``smiles`` cannot be featurised (invalid/empty).
        """
        try:
            feat = compute_features(smiles)
        except ValueError as exc:
            raise ValueError(f"Cannot featurise SMILES: {smiles!r}") from exc
        # Defensive: compute_features is documented to raise, but guard for any
        # future code path that could return None rather than raising.
        if feat is None:  # pragma: no cover - defensive only
            raise ValueError(f"Cannot featurise SMILES: {smiles!r}")
        probs = self._model.predict_proba(feat.reshape(1, -1))[0]
        return {BUCKET_NAMES[i]: float(probs[i]) for i in range(3)}
