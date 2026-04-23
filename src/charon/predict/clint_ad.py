"""CLint-local applicability domain for Tier 3 gating.

This module is intentionally decoupled from ``core/guardrails.py`` (the
project-wide Layer 0 AD, which is not yet wired into production paths).
The CLint-local AD maintains its own Morgan-fingerprint reference set
built from the CLint training compounds so that the Tier 3 trigger is
precisely tuned to CLint-regression coverage.

Cache file (`models/clint_ad_fingerprints.npz`) stores:
- ``packed_fps`` : uint8 array of packed Morgan bit vectors (N, 256).
- ``source_hash`` : SHA-256 of the training CSV at build time.
- ``training_path`` : relative path string to the source CSV.
"""
from __future__ import annotations

import hashlib
import logging
import sys
from pathlib import Path
from typing import Literal

import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem

logger = logging.getLogger(__name__)

AD_LOW_THRESHOLD = 0.3
AD_MODERATE_THRESHOLD = 0.5
FP_N_BITS = 2048
FP_RADIUS = 2

_DEFAULT_CACHE_PATH = (
    Path(__file__).resolve().parents[3] / "models" / "clint_ad_fingerprints.npz"
)
_DEFAULT_TRAINING_CSV = (
    Path(__file__).resolve().parents[3] / "data" / "training" / "clint_merged.csv"
)


def _sha256(path: Path) -> str:
    return "sha256:" + hashlib.sha256(path.read_bytes()).hexdigest()


def _morgan_bit_vect(smi: str):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(mol, FP_RADIUS, nBits=FP_N_BITS)


def _pack_fps(fps: list) -> np.ndarray:
    """Pack RDKit ExplicitBitVect list into a uint8 (N, FP_N_BITS//8) array."""
    buf = np.zeros((len(fps), FP_N_BITS // 8), dtype=np.uint8)
    for i, fp in enumerate(fps):
        arr = np.zeros(FP_N_BITS, dtype=np.uint8)
        DataStructs.ConvertToNumpyArray(fp, arr)
        buf[i] = np.packbits(arr)
    return buf


def _unpack_to_bit_vect_list(packed: np.ndarray) -> list:
    out = []
    for row in packed:
        bits = np.unpackbits(row)[:FP_N_BITS].astype(np.uint8)
        fp = DataStructs.ExplicitBitVect(FP_N_BITS)
        on_idx = np.nonzero(bits)[0].tolist()
        fp.SetBitsFromList(on_idx)
        out.append(fp)
    return out


class ClintLocalAD:
    """Tanimoto-based applicability domain for CLint predictions."""

    def __init__(self, reference_fps: list | None = None) -> None:
        self._reference_fps = reference_fps or []

    @classmethod
    def load_default(
        cls,
        cache_path: Path | None = None,
        training_csv: Path | None = None,
    ) -> "ClintLocalAD":
        cache_path = Path(cache_path) if cache_path else _DEFAULT_CACHE_PATH
        training_csv = Path(training_csv) if training_csv else _DEFAULT_TRAINING_CSV

        current_hash = _sha256(training_csv)
        if cache_path.exists():
            try:
                data = dict(np.load(cache_path, allow_pickle=True).items())
                if str(data.get("source_hash")) == current_hash:
                    fps = _unpack_to_bit_vect_list(data["packed_fps"])
                    logger.info("Loaded CLint AD cache (%d FPs) from %s",
                                len(fps), cache_path)
                    return cls(reference_fps=fps)
            except (OSError, ValueError, KeyError) as exc:
                logger.warning(
                    "Corrupted CLint AD cache at %s (%s); rebuilding",
                    cache_path, exc,
                )

        fps = cls._build_reference_fps(training_csv)
        packed = _pack_fps(fps)
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        np.savez(
            cache_path,
            packed_fps=packed,
            source_hash=np.array(current_hash),
            training_path=np.array(str(training_csv.name)),
        )
        logger.info("Built + cached CLint AD (%d FPs) -> %s", len(fps), cache_path)
        return cls(reference_fps=fps)

    @staticmethod
    def _build_reference_fps(training_csv: Path) -> list:
        """Recompute the training set Morgan FPs (same exclusions as train_clint.py)."""
        scripts_dir = training_csv.resolve().parents[2] / "scripts"
        if str(scripts_dir) not in sys.path:
            sys.path.insert(0, str(scripts_dir))
        from _train_common import (  # type: ignore[import-not-found]
            canonical_smiles, inchikey14, load_validation_keys,
        )
        import pandas as pd  # local import keeps module-import light

        df = pd.read_csv(training_csv)
        df = df[df["source"].isin({"tdc_hep", "chembl"})].copy()
        df["canonical_smiles"] = df["smiles"].apply(canonical_smiles)
        df = df.dropna(subset=["canonical_smiles"]).drop_duplicates("canonical_smiles")
        val_keys = load_validation_keys()
        df["ik14"] = df["canonical_smiles"].apply(inchikey14)
        df = df[~df["ik14"].isin(val_keys)]
        fps = [_morgan_bit_vect(s) for s in df["canonical_smiles"]]
        return [fp for fp in fps if fp is not None]

    def max_similarity(self, smiles: str) -> float:
        if not self._reference_fps:
            return float("nan")
        fp = _morgan_bit_vect(smiles)
        if fp is None:
            return float("nan")
        return max(
            DataStructs.TanimotoSimilarity(fp, ref)
            for ref in self._reference_fps
        )

    def classify(
        self, max_sim: float
    ) -> Literal["HIGH", "MODERATE", "LOW", "INVALID"]:
        if not np.isfinite(max_sim):
            return "INVALID"
        if max_sim >= AD_MODERATE_THRESHOLD:
            return "HIGH"
        if max_sim >= AD_LOW_THRESHOLD:
            return "MODERATE"
        return "LOW"
