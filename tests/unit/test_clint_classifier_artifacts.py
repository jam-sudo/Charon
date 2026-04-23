"""Sanity checks for the CLint classifier artifacts produced by training."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
MODELS_DIR = REPO_ROOT / "models"


def test_classifier_model_file_exists():
    path = MODELS_DIR / "xgboost_clint_classifier.json"
    assert path.exists(), (
        f"{path.name} missing — rerun `python scripts/train_clint_classifier.py`"
    )


def test_oof_probs_shape_and_normalisation():
    arr = np.load(MODELS_DIR / "xgboost_clint_classifier_oof_probs.npy")
    assert arr.ndim == 2
    assert arr.shape[1] == 3, "Expect 3 class probabilities per row"
    assert arr.shape[0] > 100, "Expect >100 OOF samples"
    assert np.all(np.isfinite(arr))
    row_sums = arr.sum(axis=1)
    assert np.allclose(row_sums, 1.0, atol=1e-5), (
        f"OOF rows must sum to 1, got max |sum - 1| = {np.max(np.abs(row_sums - 1)):.3e}"
    )


def test_metadata_records_classifier():
    meta = json.loads((MODELS_DIR / "model_metadata.json").read_text())
    entry = meta["xgboost_clint_classifier"]
    assert entry["n_samples"] > 100
    assert entry["macro_f1"] >= 0.50, (
        f"Macro F1 gate: expected >= 0.50, got {entry['macro_f1']:.3f}"
    )
    assert entry["bucket_boundaries"] == [10.0, 50.0]
    assert set(entry["bucket_names"]) == {"low", "med", "high"}
