"""Train Charon's fup (fraction unbound in plasma) XGBoost model.

Data source: TDC PPBR_AZ (``data/training/fup_ppbr_az.tab``), filtered
to Homo sapiens. Compounds matching the validation InChIKey-14 set are
excluded so adme_reference.csv remains a clean holdout.

Target: ``log10(fup)``, which symmetrises the loss across the [0.001, 1.0]
range. Metrics are reported with honest Murcko-scaffold 5-fold CV.

Usage:
    python3 scripts/train_fup.py
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
    scaffold_oof_predictions,
    update_model_metadata,
)
from charon.predict.features import FEATURE_LENGTH, compute_features  # noqa: E402


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger("charon.train.fup")

FUP_SOURCE = DATA_TRAIN_DIR / "fup_ppbr_az.tab"
MODEL_OUT = MODELS_DIR / "xgboost_fup.json"


def load_human_fup() -> pd.DataFrame:
    """Load PPBR_AZ, keep only humans, compute fup, drop unusable rows."""
    if not FUP_SOURCE.exists():
        raise FileNotFoundError(f"fup training data not found: {FUP_SOURCE}")
    df = pd.read_csv(FUP_SOURCE, sep="\t")
    log.info("Loaded %d rows from %s", len(df), FUP_SOURCE.name)

    # Species filter. Pandas strips the CSV quotes, so we match the bare name.
    df = df[df["Species"] == "Homo sapiens"].copy()
    log.info("After Homo sapiens filter: %d rows", len(df))
    assert len(df) == 1614, f"Expected 1614 human rows, got {len(df)}"

    # Strip quoted Drug_ID for cleaner logging.
    df["Drug"] = df["Drug"].astype(str).str.strip().str.strip('"')
    df["fup"] = 1.0 - df["Y"].astype(float) / 100.0

    # Physical bounds: fup must be in (0, 1]. Values at the boundaries are
    # kept (clipped) because both fully-bound and fully-unbound drugs exist.
    before = len(df)
    df = df[(df["fup"] > 0.0) & (df["fup"] <= 1.0)].copy()
    log.info("After fup range filter: %d rows (removed %d)", len(df), before - len(df))

    # Clip fup to (0.001, 0.999) for log-space training stability.
    df["fup"] = df["fup"].clip(lower=0.001, upper=0.999)

    # Canonical SMILES + dedup by canonical SMILES (duplicates can exist).
    df["canonical_smiles"] = df["Drug"].map(canonical_smiles)
    df = df.dropna(subset=["canonical_smiles"])
    df = df.drop_duplicates(subset=["canonical_smiles"]).reset_index(drop=True)
    log.info("After canonicalization + dedup: %d unique compounds", len(df))
    return df


def exclude_validation(df: pd.DataFrame) -> pd.DataFrame:
    """Drop training rows that match the validation InChIKey-14 set."""
    val_keys = load_validation_keys()
    log.info("Loaded %d validation InChIKey-14 keys", len(val_keys))
    df = df.copy()
    df["ik14"] = df["canonical_smiles"].map(inchikey14)
    df = df.dropna(subset=["ik14"])
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


def make_fup_model() -> xgb.XGBRegressor:
    """Fresh XGBoost regressor matching the Sisyphus v1 fup hyperparameters."""
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

    df = load_human_fup()
    df = exclude_validation(df)

    X, valid_idx = build_feature_matrix(df["canonical_smiles"].tolist())
    y_fup = df["fup"].to_numpy()[valid_idx]
    smiles = [df["canonical_smiles"].iloc[i] for i in valid_idx]
    log.info("Feature matrix: %s", X.shape)

    groups = [murcko_scaffold(s) for s in smiles]
    log.info("Unique Murcko scaffolds: %d", len(set(groups)))

    y_log = np.log10(y_fup)

    log.info("Running scaffold 5-fold CV ...")
    oof_log = scaffold_oof_predictions(make_fup_model, X, y_log, groups, n_splits=5)
    oof_linear = np.clip(10.0 ** oof_log, 0.001, 1.0)

    r2, mae, aafe, pct2, pct3 = compute_metrics(
        y_true_target=y_log,
        y_pred_target=oof_log,
        y_true_linear=y_fup,
        y_pred_linear=oof_linear,
    )
    log.info("Scaffold CV R²=%.3f, MAE(log10)=%.3f, AAFE=%.3f", r2, mae, aafe)
    log.info("%% within 2-fold=%.1f%%, 3-fold=%.1f%%", pct2, pct3)

    log.info("Training final model on %d compounds ...", len(y_log))
    final_model = make_fup_model()
    final_model.fit(X, y_log)
    final_model.save_model(str(MODEL_OUT))
    log.info("Saved model to %s", MODEL_OUT)

    report = CVReport(
        model_name="xgboost_fup",
        target="log10(fup)",
        n_samples=len(y_log),
        n_folds=5,
        r2=r2,
        mae=mae,
        aafe=aafe,
        pct_within_2fold=pct2,
        pct_within_3fold=pct3,
        training_sources=["tdc_ppbr_az (Homo sapiens)"],
        feature_length=FEATURE_LENGTH,
        notes=(
            "Target: log10(fup). Scaffold 5-fold GroupKFold on Murcko "
            "scaffolds. Validation set (adme_reference.csv InChIKey-14) "
            "excluded from training."
        ),
    )
    update_model_metadata(report)


if __name__ == "__main__":
    main()
