"""Train Charon's 3-class CLint classifier (Sprint 8 Tier 3 fallback).

Buckets (frozen per CLAUDE.md §6j):
  Low  : [0.1, 10.0)   µL/min/10^6 cells  -> label 0
  Med  : [10.0, 50.0)  µL/min/10^6 cells  -> label 1
  High : [50.0, 1000]  µL/min/10^6 cells  -> label 2

Uses the same scaffold-CV + validation-exclusion pipeline as train_clint.py.
Gate: scaffold-CV Macro F1 >= 0.50 — aborts otherwise.

Usage:
    python3 scripts/train_clint_classifier.py
"""
from __future__ import annotations

import hashlib
import logging
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import xgboost as xgb

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "src"))
sys.path.insert(0, str(REPO_ROOT / "scripts"))

from _train_common import (  # noqa: E402
    DATA_TRAIN_DIR,
    MODELS_DIR,
    canonical_smiles,
    inchikey14,
    load_validation_keys,
    murcko_scaffold,
    update_model_metadata_dict,
)
from charon.predict.features import FEATURE_LENGTH, compute_features  # noqa: E402
from sklearn.metrics import accuracy_score, confusion_matrix, f1_score  # noqa: E402
from sklearn.model_selection import GroupKFold  # noqa: E402


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger("charon.train.clint_classifier")

CLINT_SOURCE = DATA_TRAIN_DIR / "clint_merged.csv"
MODEL_OUT = MODELS_DIR / "xgboost_clint_classifier.json"
OOF_PROBS_OUT = MODELS_DIR / "xgboost_clint_classifier_oof_probs.npy"

CLEAN_SOURCES = {"tdc_hep", "chembl"}
CLINT_MIN, CLINT_MAX = 0.1, 1000.0
BUCKET_BOUNDARIES = [10.0, 50.0]  # <10 Low, <50 Med, else High
BUCKET_NAMES = ["low", "med", "high"]
MACRO_F1_GATE = 0.50


def bucket_label(v: float) -> int:
    if v < BUCKET_BOUNDARIES[0]:
        return 0
    if v < BUCKET_BOUNDARIES[1]:
        return 1
    return 2


def make_classifier() -> xgb.XGBClassifier:
    return xgb.XGBClassifier(
        n_estimators=200,
        max_depth=5,
        learning_rate=0.08,
        subsample=0.8,
        colsample_bytree=0.6,
        random_state=42,
        n_jobs=-1,
        tree_method="hist",
        objective="multi:softmax",
        num_class=3,
    )


def load_clint() -> pd.DataFrame:
    df = pd.read_csv(CLINT_SOURCE)
    df = df[df["source"].isin(CLEAN_SOURCES)].copy()
    df["clint_hep"] = pd.to_numeric(df["clint_hep"], errors="coerce")
    df = df.dropna(subset=["clint_hep"])
    df["canonical_smiles"] = df["smiles"].apply(canonical_smiles)
    df = df.dropna(subset=["canonical_smiles"]).drop_duplicates("canonical_smiles")
    return df


def exclude_validation(df: pd.DataFrame) -> pd.DataFrame:
    val_keys = load_validation_keys()
    df = df.copy()
    df["ik14"] = df["canonical_smiles"].apply(inchikey14)
    before = len(df)
    df = df[~df["ik14"].isin(val_keys)]
    log.info("Excluded %d validation-overlapping compounds; %d remaining",
             before - len(df), len(df))
    return df


def build_features(smiles_list: list[str]) -> tuple[np.ndarray, list[int]]:
    X_rows, valid_idx = [], []
    for i, smi in enumerate(smiles_list):
        try:
            feat = compute_features(smi)
        except ValueError:
            continue
        X_rows.append(feat)
        valid_idx.append(i)
    if not X_rows:
        return np.zeros((0, FEATURE_LENGTH)), []
    return np.vstack(X_rows), valid_idx


def scaffold_oof_probs(X: np.ndarray, y: np.ndarray, groups: list[str],
                      n_splits: int = 5) -> np.ndarray:
    gkf = GroupKFold(n_splits=n_splits)
    oof = np.zeros((len(y), 3))
    for fold, (tr, te) in enumerate(gkf.split(X, y, groups=groups)):
        clf = make_classifier()
        clf.fit(X[tr], y[tr])
        oof[te] = clf.predict_proba(X[te])
        log.info("  Fold %d: %d train / %d test", fold + 1, len(tr), len(te))
    return oof


def main() -> None:
    MODELS_DIR.mkdir(parents=True, exist_ok=True)

    df = exclude_validation(load_clint())
    X, valid_idx = build_features(df["canonical_smiles"].tolist())
    smiles = [df["canonical_smiles"].iloc[i] for i in valid_idx]
    y_linear = df["clint_hep"].iloc[valid_idx].to_numpy()
    y_clipped = np.clip(y_linear, CLINT_MIN, CLINT_MAX)
    y = np.array([bucket_label(v) for v in y_clipped])

    log.info("Features: %s; class counts %s", X.shape,
             dict(zip(*np.unique(y, return_counts=True))))

    groups = [murcko_scaffold(s) for s in smiles]
    log.info("Unique scaffolds: %d", len(set(groups)))

    log.info("Scaffold 5-fold CV ...")
    oof_probs = scaffold_oof_probs(X, y, groups, n_splits=5)
    y_pred = oof_probs.argmax(axis=1)

    acc = float(accuracy_score(y, y_pred))
    macro_f1 = float(f1_score(y, y_pred, average="macro"))
    per_class_f1 = f1_score(y, y_pred, average=None).tolist()
    cm = confusion_matrix(y, y_pred).tolist()

    log.info("Accuracy=%.3f  Macro F1=%.3f", acc, macro_f1)
    log.info("Per-class F1 (low/med/high)=%s", [round(x, 3) for x in per_class_f1])
    log.info("Confusion matrix (rows=true, cols=pred): %s", cm)

    if macro_f1 < MACRO_F1_GATE:
        raise SystemExit(
            f"Macro F1 {macro_f1:.3f} below gate {MACRO_F1_GATE}. "
            "Classifier is not fit for Tier 3 release."
        )

    log.info("Training final model on full data ...")
    final = make_classifier()
    final.fit(X, y)
    final.save_model(str(MODEL_OUT))
    np.save(OOF_PROBS_OUT, oof_probs)
    oof_sha = hashlib.sha256(OOF_PROBS_OUT.read_bytes()).hexdigest()

    update_model_metadata_dict("xgboost_clint_classifier", {
        "model_name": "xgboost_clint_classifier",
        "target": "clint_hep bucket {low <10, med 10-50, high >=50} uL/min/10^6 cells",
        "n_samples": int(len(y)),
        "n_folds": 5,
        "accuracy": acc,
        "macro_f1": macro_f1,
        "per_class_f1": per_class_f1,
        "confusion_matrix": cm,
        "class_counts": {
            BUCKET_NAMES[i]: int((y == i).sum()) for i in range(3)
        },
        "training_sources": sorted(CLEAN_SOURCES),
        "feature_length": FEATURE_LENGTH,
        "bucket_boundaries": BUCKET_BOUNDARIES,
        "bucket_names": BUCKET_NAMES,
        "oof_probs_path": "xgboost_clint_classifier_oof_probs.npy",
        "oof_probs_sha256": oof_sha,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "notes": (
            "3-class classification fallback for CLint when the compound "
            "is outside the CLint-local AD (max Tanimoto < 0.3 vs "
            "training set). Gate: scaffold-CV Macro F1 >= 0.50. "
            "Buckets are frozen per CLAUDE.md §6j."
        ),
    })
    log.info("Saved model -> %s, OOF probs -> %s, metadata updated.",
             MODEL_OUT, OOF_PROBS_OUT)


if __name__ == "__main__":
    main()
