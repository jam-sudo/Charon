"""Shared helpers for Charon ADME training scripts.

These are intentionally small: canonical SMILES / InChIKey-14 helpers,
Murcko scaffold computation, scaffold-based CV metrics, and utilities for
loading the validation holdout so we never train on calibration compounds.
"""

from __future__ import annotations

import hashlib
import json
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem  # noqa: F401  (needed for fingerprint imports used elsewhere)
from rdkit.Chem.inchi import InchiToInchiKey, MolToInchi
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.model_selection import GroupKFold

# Suppress RDKit warnings during data loading (noisy for large datasets).
RDLogger.DisableLog("rdApp.warning")

REPO_ROOT = Path(__file__).resolve().parent.parent
MODELS_DIR = REPO_ROOT / "models"
DATA_TRAIN_DIR = REPO_ROOT / "data" / "training"
DATA_VAL_DIR = REPO_ROOT / "data" / "validation"
VALIDATION_CSV = DATA_VAL_DIR / "adme_reference.csv"
MODEL_METADATA_JSON = MODELS_DIR / "model_metadata.json"


log = logging.getLogger("charon.train")


# ---------------------------------------------------------------------------
# Identifier helpers
# ---------------------------------------------------------------------------


def canonical_smiles(smiles: str) -> str | None:
    """Return RDKit canonical SMILES or None on failure."""
    if not isinstance(smiles, str) or not smiles.strip():
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return Chem.MolToSmiles(mol, isomericSmiles=True)


def inchikey14(smiles: str) -> str | None:
    """Return the connectivity-only InChIKey (first 14 chars) or None."""
    if not isinstance(smiles, str) or not smiles.strip():
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    inchi = MolToInchi(mol)
    if inchi is None:
        return None
    key = InchiToInchiKey(inchi)
    if key is None:
        return None
    return key[:14]


def murcko_scaffold(smiles: str) -> str:
    """Return canonical SMILES of the Bemis-Murcko scaffold.

    Used as the grouping key for ``GroupKFold``. Compounds with invalid
    SMILES are all placed in a single ``"__invalid__"`` bucket so they
    never accidentally cross train/test boundaries.
    """
    mol = Chem.MolFromSmiles(smiles) if isinstance(smiles, str) else None
    if mol is None:
        return "__invalid__"
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    if scaffold is None:
        return "__invalid__"
    return Chem.MolToSmiles(scaffold)


# ---------------------------------------------------------------------------
# Validation holdout loading
# ---------------------------------------------------------------------------


def load_validation_keys() -> set[str]:
    """Return the set of InChIKey-14s present in ``adme_reference.csv``.

    Training data matching any of these keys must be excluded to prevent
    leaking calibration compounds into the training set.
    """
    if not VALIDATION_CSV.exists():
        raise FileNotFoundError(f"Validation CSV missing: {VALIDATION_CSV}")
    df = pd.read_csv(VALIDATION_CSV)
    keys: set[str] = set()
    for smi in df["smiles"].dropna():
        key = inchikey14(smi)
        if key:
            keys.add(key)
    return keys


# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------


@dataclass
class CVReport:
    """Aggregated metrics for a scaffold-CV training run."""

    model_name: str
    target: str                       # e.g. "log10(fup)"
    n_samples: int
    n_folds: int
    r2: float                         # 1 - SS_res / SS_tot on the target
    mae: float                        # mean absolute error on the target
    aafe: float                       # 10^(mean |log10(pred/obs)|) in original space
    pct_within_2fold: float
    pct_within_3fold: float
    training_sources: list[str]
    feature_length: int
    notes: str = ""
    extras: dict = field(default_factory=dict)

    def as_dict(self) -> dict:
        return {
            "model_name": self.model_name,
            "target": self.target,
            "n_samples": self.n_samples,
            "n_folds": self.n_folds,
            "scaffold_cv_r2": self.r2,
            "scaffold_cv_mae": self.mae,
            "scaffold_cv_aafe": self.aafe,
            "pct_within_2fold": self.pct_within_2fold,
            "pct_within_3fold": self.pct_within_3fold,
            "training_sources": self.training_sources,
            "feature_length": self.feature_length,
            "notes": self.notes,
            **self.extras,
        }


def compute_metrics(
    y_true_target: np.ndarray,
    y_pred_target: np.ndarray,
    y_true_linear: np.ndarray,
    y_pred_linear: np.ndarray,
) -> tuple[float, float, float, float, float]:
    """Compute R² / MAE (on the model target) and AAFE / %within-fold (linear)."""
    ss_res = float(np.sum((y_true_target - y_pred_target) ** 2))
    ss_tot = float(np.sum((y_true_target - y_true_target.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
    mae = float(np.mean(np.abs(y_true_target - y_pred_target)))

    eps = 1e-12
    y_true_safe = np.maximum(y_true_linear, eps)
    y_pred_safe = np.maximum(y_pred_linear, eps)
    log_ratio = np.abs(np.log10(y_pred_safe / y_true_safe))
    aafe = float(10 ** np.mean(log_ratio))
    pct2 = float(np.mean(log_ratio <= np.log10(2.0)) * 100.0)
    pct3 = float(np.mean(log_ratio <= np.log10(3.0)) * 100.0)
    return r2, mae, aafe, pct2, pct3


# ---------------------------------------------------------------------------
# Scaffold cross-validation predictions
# ---------------------------------------------------------------------------


def scaffold_oof_predictions(
    make_model,
    X: np.ndarray,
    y: np.ndarray,
    groups: Iterable[str],
    n_splits: int = 5,
) -> np.ndarray:
    """Run GroupKFold (by scaffold) and return out-of-fold predictions.

    ``make_model`` must be a zero-arg callable that constructs a fresh
    untrained XGBoost estimator; this avoids leaking state between folds.
    """
    groups_arr = np.asarray(list(groups))
    unique = len(set(groups_arr))
    if unique < n_splits:
        raise ValueError(
            f"Not enough unique scaffolds ({unique}) for {n_splits}-fold CV. "
            "Check data filtering."
        )
    kf = GroupKFold(n_splits=n_splits)
    oof = np.full(len(y), np.nan, dtype=np.float64)
    for fold_idx, (train_idx, test_idx) in enumerate(kf.split(X, y, groups=groups_arr), start=1):
        model = make_model()
        model.fit(X[train_idx], y[train_idx])
        oof[test_idx] = model.predict(X[test_idx])
        log.info("  fold %d: train=%d, test=%d", fold_idx, len(train_idx), len(test_idx))
    assert not np.any(np.isnan(oof)), "Some samples left without a fold assignment"
    return oof


# ---------------------------------------------------------------------------
# Metadata persistence
# ---------------------------------------------------------------------------


def persist_oof_residuals(
    model_stem: str,
    oof_log: np.ndarray,
    y_log: np.ndarray,
) -> tuple[str, str]:
    """Save scaffold-CV |log10(pred/obs)| residuals and return (filename, sha256 hex).

    Filename is ``{model_stem}_oof_residuals.npy`` relative to MODELS_DIR.
    Consumed by ``ConformalPredictor.calibrate_from_oof`` for CLint/fup CIs.
    """
    residuals = np.abs(oof_log - y_log)
    path = MODELS_DIR / f"{model_stem}_oof_residuals.npy"
    np.save(path, residuals)
    sha = hashlib.sha256(path.read_bytes()).hexdigest()
    return path.name, sha


def update_model_metadata(report: CVReport) -> None:
    """Merge ``report`` into ``models/model_metadata.json``."""
    MODELS_DIR.mkdir(parents=True, exist_ok=True)
    metadata: dict = {}
    if MODEL_METADATA_JSON.exists():
        try:
            metadata = json.loads(MODEL_METADATA_JSON.read_text())
        except json.JSONDecodeError:
            log.warning("model_metadata.json is corrupt, starting fresh")
            metadata = {}
    metadata[report.model_name] = report.as_dict()
    MODEL_METADATA_JSON.write_text(json.dumps(metadata, indent=2, sort_keys=True))
    log.info("Wrote metadata for %s to %s", report.model_name, MODEL_METADATA_JSON)
