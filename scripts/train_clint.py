"""Train Charon's CLint (intrinsic clearance) XGBoost model.

Data source: ``data/training/clint_merged.csv`` (copied from Sisyphus).
Only the clean hepatocyte sources are used:

  - ``tdc_hep``   (978 rows, μL/min/10^6 cells, TDC Clearance_Hepatocyte_AZ)
  - ``chembl``    (517 rows, same units, Sisyphus ChEMBL curation)

The ``biogen_fang`` partition is excluded because its numerical scale
(mean 67,586) is inconsistent with hepatocyte μL/min/10^6 cells — it was
never successfully re-normalized by Sisyphus. The ``tdc_mic`` partition is
excluded because it is microsomal (μL/min/mg protein), a different assay.

Compounds matching the validation InChIKey-14 set are excluded.

Target: ``log10(clint_hep)`` clipped to [log10(0.1), log10(1000)].

Usage:
    python3 scripts/train_clint.py
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xgboost as xgb

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "src"))
sys.path.insert(0, str(REPO_ROOT / "scripts"))

from _train_common import (  # noqa: E402
    CVReport,
    DATA_TRAIN_DIR,
    MODELS_DIR,
    canonical_smiles,
    compute_metrics,
    inchikey14,
    load_validation_keys,
    murcko_scaffold,
    persist_oof_residuals,
    scaffold_oof_predictions,
    update_model_metadata,
)
from charon.predict.features import FEATURE_LENGTH, compute_features  # noqa: E402


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger("charon.train.clint")

CLINT_SOURCE = DATA_TRAIN_DIR / "clint_merged.csv"
MODEL_OUT = MODELS_DIR / "xgboost_clint.json"

CLEAN_SOURCES = {"tdc_hep", "chembl"}
CLINT_MIN = 0.1       # μL/min/10^6 cells
CLINT_MAX = 1000.0    # μL/min/10^6 cells


def load_clint() -> pd.DataFrame:
    """Load merged CLint file and restrict to clean hepatocyte sources."""
    if not CLINT_SOURCE.exists():
        raise FileNotFoundError(f"CLint training data not found: {CLINT_SOURCE}")
    df = pd.read_csv(CLINT_SOURCE)
    log.info("Loaded %d rows from %s", len(df), CLINT_SOURCE.name)
    log.info("Source breakdown: %s", df["source"].value_counts().to_dict())

    df = df[df["source"].isin(CLEAN_SOURCES)].copy()
    log.info("After source filter (tdc_hep + chembl): %d rows", len(df))
    assert len(df) == 1495, f"Expected 1495 rows, got {len(df)}"

    # Physical range filter: drop zeros and unreasonable extremes.
    before = len(df)
    df = df[(df["clint_hep"] > CLINT_MIN) & (df["clint_hep"] < CLINT_MAX)].copy()
    log.info("After range filter (%.1f-%.1f): %d rows (removed %d)",
             CLINT_MIN, CLINT_MAX, len(df), before - len(df))

    # Dedup by InChIKey-14 (already present) and re-canonicalize.
    df["canonical_smiles"] = df["smiles"].map(canonical_smiles)
    df = df.dropna(subset=["canonical_smiles"])
    df = df.drop_duplicates(subset=["ik14"]).reset_index(drop=True)
    log.info("After dedup: %d unique compounds", len(df))
    return df


def exclude_validation(df: pd.DataFrame) -> pd.DataFrame:
    """Drop training rows that match the validation InChIKey-14 set."""
    val_keys = load_validation_keys()
    log.info("Loaded %d validation InChIKey-14 keys", len(val_keys))
    before = len(df)
    df = df[~df["ik14"].isin(val_keys)].reset_index(drop=True)
    log.info("Excluded %d validation compounds; %d remain", before - len(df), len(df))
    return df


def build_feature_matrix(smiles: list[str]) -> tuple[np.ndarray, list[int]]:
    """Compute features, returning (X, valid_indices)."""
    X_rows: list[np.ndarray] = []
    valid_idx: list[int] = []
    for i, smi in enumerate(smiles):
        try:
            X_rows.append(compute_features(smi))
            valid_idx.append(i)
        except ValueError:
            continue
    X = np.vstack(X_rows) if X_rows else np.zeros((0, FEATURE_LENGTH))
    return X, valid_idx


def make_clint_model() -> xgb.XGBRegressor:
    """Fresh XGBoost regressor (matches Sisyphus v1 clint hyperparameters)."""
    return xgb.XGBRegressor(
        n_estimators=200,
        max_depth=5,
        learning_rate=0.08,
        subsample=0.8,
        colsample_bytree=0.6,
        random_state=42,
        n_jobs=-1,
        tree_method="hist",
        objective="reg:squarederror",
    )


def main() -> None:
    MODELS_DIR.mkdir(parents=True, exist_ok=True)

    df = load_clint()
    df = exclude_validation(df)

    X, valid_idx = build_feature_matrix(df["canonical_smiles"].tolist())
    y_linear = df["clint_hep"].to_numpy()[valid_idx]
    smiles = [df["canonical_smiles"].iloc[i] for i in valid_idx]
    log.info("Feature matrix: %s", X.shape)

    y_log = np.log10(np.clip(y_linear, CLINT_MIN, CLINT_MAX))
    groups = [murcko_scaffold(s) for s in smiles]
    log.info("Unique Murcko scaffolds: %d", len(set(groups)))

    log.info("Running scaffold 5-fold CV ...")
    oof_log = scaffold_oof_predictions(make_clint_model, X, y_log, groups, n_splits=5)
    oof_linear = np.clip(10.0 ** oof_log, CLINT_MIN, CLINT_MAX)

    # Persist |log10(pred/obs)| residuals for ConformalPredictor.calibrate_from_oof
    residuals_name, residuals_sha = persist_oof_residuals(
        "xgboost_clint", oof_log, y_log
    )
    log.info("Saved OOF residuals to models/%s", residuals_name)

    r2, mae, aafe, pct2, pct3 = compute_metrics(
        y_true_target=y_log,
        y_pred_target=oof_log,
        y_true_linear=y_linear,
        y_pred_linear=oof_linear,
    )
    log.info("Scaffold CV R²=%.3f, MAE(log10)=%.3f, AAFE=%.3f", r2, mae, aafe)
    log.info("%% within 2-fold=%.1f%%, 3-fold=%.1f%%", pct2, pct3)

    log.info("Training final model on %d compounds ...", len(y_log))
    final_model = make_clint_model()
    final_model.fit(X, y_log)
    final_model.save_model(str(MODEL_OUT))
    log.info("Saved model to %s", MODEL_OUT)

    report = CVReport(
        model_name="xgboost_clint",
        target="log10(clint_hep_uL_min_10e6_cells)",
        n_samples=len(y_log),
        n_folds=5,
        r2=r2,
        mae=mae,
        aafe=aafe,
        pct_within_2fold=pct2,
        pct_within_3fold=pct3,
        training_sources=["tdc_hep", "chembl"],
        feature_length=FEATURE_LENGTH,
        notes=(
            "Target: log10(CLint) in uL/min/10^6 cells, clipped to "
            "[0.1, 1000]. Scaffold 5-fold GroupKFold. Sources biogen_fang "
            "and tdc_mic explicitly excluded (unit mismatch)."
        ),
        extras={
            "oof_residuals_path": residuals_name,
            "oof_residuals_sha256": residuals_sha,
        },
    )
    update_model_metadata(report)


if __name__ == "__main__":
    main()
