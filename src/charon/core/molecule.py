"""Molecular representation using RDKit for guardrails and Layer 1."""

from __future__ import annotations

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Crippen, Descriptors, rdMolDescriptors


class Molecule:
    """Parse SMILES and compute molecular descriptors for guardrails and Layer 1."""

    def __init__(self, smiles: str) -> None:
        """
        Parse SMILES string into RDKit mol object.

        Raises ValueError if SMILES is invalid (None, empty, whitespace, or
        unparseable by RDKit).  Stores canonical SMILES on success.
        """
        if smiles is None:
            raise ValueError("SMILES must not be None")
        if not isinstance(smiles, str) or not smiles.strip():
            raise ValueError("SMILES must be a non-empty string")

        cleaned = smiles.strip()
        self._mol = Chem.MolFromSmiles(cleaned)
        if self._mol is None:
            raise ValueError(f"Invalid SMILES: {cleaned!r}")

        self._canonical_smiles = Chem.MolToSmiles(self._mol)

    @property
    def mol(self) -> Chem.Mol:
        """RDKit mol object."""
        return self._mol

    @property
    def canonical_smiles(self) -> str:
        """Canonical SMILES string."""
        return self._canonical_smiles

    def descriptors(self) -> dict[str, float]:
        """
        Compute all descriptors needed for guardrails.

        Returns dict with keys:
          MW, logP (Wildman-Crippen cLogP), TPSA, HBD, HBA,
          rotatable_bonds, aromatic_rings, heavy_atom_count

        Note: logP here is the descriptor-level cLogP from Crippen, not
        the ML prediction from Layer 1.  Used only for guardrails and
        drug-likeness checks.
        """
        m = self._mol
        return {
            "MW": Descriptors.ExactMolWt(m),
            "logP": Crippen.MolLogP(m),
            "TPSA": Descriptors.TPSA(m),
            "HBD": rdMolDescriptors.CalcNumHBD(m),
            "HBA": rdMolDescriptors.CalcNumHBA(m),
            "rotatable_bonds": rdMolDescriptors.CalcNumRotatableBonds(m),
            "aromatic_rings": rdMolDescriptors.CalcNumAromaticRings(m),
            "heavy_atom_count": m.GetNumHeavyAtoms(),
        }

    def morgan_fingerprint(self, radius: int = 2, n_bits: int = 2048):
        """
        Generate Morgan fingerprint for Tanimoto similarity.

        Uses radius=2, 2048 bits as specified in ARCHITECTURE.md.
        Returns RDKit ExplicitBitVect.
        """
        return AllChem.GetMorganFingerprintAsBitVect(
            self._mol, radius, nBits=n_bits
        )

    def tanimoto_similarity(self, other: Molecule) -> float:
        """
        Compute Tanimoto coefficient between this molecule and another.

        Uses Morgan fingerprints (radius=2, 2048 bits).
        Returns float in [0, 1].
        """
        fp_self = self.morgan_fingerprint()
        fp_other = other.morgan_fingerprint()
        return DataStructs.TanimotoSimilarity(fp_self, fp_other)
