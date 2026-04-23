"""Sanity tests for OOF residual files produced by training scripts."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
MODELS_DIR = REPO_ROOT / "models"


@pytest.mark.parametrize("stem", ["xgboost_clint", "xgboost_fup"])
def test_oof_residual_file_exists(stem: str):
    path = MODELS_DIR / f"{stem}_oof_residuals.npy"
    assert path.exists(), (
        f"{path.name} missing — rerun `python scripts/train_{stem.split('_')[1]}.py`"
    )


@pytest.mark.parametrize("stem", ["xgboost_clint", "xgboost_fup"])
def test_oof_residual_shape_and_values(stem: str):
    arr = np.load(MODELS_DIR / f"{stem}_oof_residuals.npy")
    assert arr.ndim == 1
    assert arr.size > 100, "Expect >100 OOF samples"
    assert np.all(np.isfinite(arr))
    assert np.all(arr >= 0.0), "Residuals are |log10(pred/obs)|, always >= 0"


@pytest.mark.parametrize("stem", ["xgboost_clint", "xgboost_fup"])
def test_metadata_records_residual_path(stem: str):
    meta = json.loads((MODELS_DIR / "model_metadata.json").read_text())
    entry = meta[stem]
    assert entry["oof_residuals_path"].endswith(f"{stem}_oof_residuals.npy")
    assert "oof_residuals_sha256" in entry
