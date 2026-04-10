"""Molecular feature computation for XGBoost ADME models.

Produces a 2057-D feature vector: 2048-bit Morgan fingerprint (radius=2)
+ 9 normalized RDKit descriptors.

IMPORTANT: Feature order and normalization constants are frozen. Changing
them invalidates all trained models. This layout is verified to match
Sisyphus ``descriptors.compute_features`` exactly so trained models can be
compared across projects without silent feature-order bugs.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

FEATURE_LENGTH = 2057
MORGAN_BITS = 2048
MORGAN_RADIUS = 2


def compute_features(smiles: str) -> NDArray[np.float64]:
    """Compute 2057-D feature vector from a SMILES string.

    Feature layout (verified against Sisyphus ``src/sisyphus/descriptors.py``):

    - [0:2048]  Morgan fingerprint (radius=2, nBits=2048)
    - [2048]    ``Descriptors.MolLogP`` (Crippen cLogP, raw)
    - [2049]    ``Descriptors.TPSA`` / 200.0
    - [2050]    ``Descriptors.MolWt`` / 600.0
    - [2051]    ``Descriptors.NumHAcceptors`` / 10.0
    - [2052]    ``Descriptors.NumHDonors`` / 5.0
    - [2053]    ``Descriptors.NumRotatableBonds`` / 15.0
    - [2054]    ``Descriptors.RingCount`` / 5.0
    - [2055]    ``Descriptors.FractionCSP3`` (raw)
    - [2056]    ``Descriptors.MolMR`` / 150.0

    Args:
        smiles: Canonical or raw SMILES.

    Returns:
        ``numpy.float64`` array of length 2057.

    Raises:
        ValueError: If ``smiles`` is empty, whitespace-only, or cannot be
            parsed by RDKit.
    """
    if not smiles or not smiles.strip():
        raise ValueError(f"Invalid SMILES: {smiles!r}")

    mol = Chem.MolFromSmiles(smiles)
    if mol is None or mol.GetNumAtoms() == 0:
        raise ValueError(f"Invalid SMILES: {smiles!r}")

    # Morgan fingerprint (2048 bits, radius 2)
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=MORGAN_RADIUS, nBits=MORGAN_BITS)
    fp_array = np.zeros(MORGAN_BITS, dtype=np.float64)
    for bit in fp.GetOnBits():
        fp_array[bit] = 1.0

    # 9 RDKit descriptors (normalized)
    descriptors = np.array(
        [
            Descriptors.MolLogP(mol),                   # LogP (Crippen), raw
            Descriptors.TPSA(mol) / 200.0,              # TPSA normalized
            Descriptors.MolWt(mol) / 600.0,             # MW normalized
            Descriptors.NumHAcceptors(mol) / 10.0,      # HBA normalized
            Descriptors.NumHDonors(mol) / 5.0,          # HBD normalized
            Descriptors.NumRotatableBonds(mol) / 15.0,  # RotBonds normalized
            Descriptors.RingCount(mol) / 5.0,           # Rings normalized
            Descriptors.FractionCSP3(mol),              # sp3 fraction, raw
            Descriptors.MolMR(mol) / 150.0,             # Molar refractivity normalized
        ],
        dtype=np.float64,
    )

    return np.concatenate([fp_array, descriptors])
