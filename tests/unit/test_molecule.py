"""Unit tests for charon.core.molecule.Molecule class."""

from __future__ import annotations

import pytest

from charon.core.molecule import Molecule


# ---------------------------------------------------------------------------
# Construction & SMILES validation
# ---------------------------------------------------------------------------


class TestMoleculeConstruction:
    """Test Molecule creation from SMILES strings."""

    def test_aspirin_smiles(self):
        mol = Molecule("CC(=O)Oc1ccccc1C(=O)O")
        assert mol.mol is not None
        assert mol.canonical_smiles is not None

    def test_aspirin_mw(self):
        mol = Molecule("CC(=O)Oc1ccccc1C(=O)O")
        desc = mol.descriptors()
        # ExactMolWt returns monoisotopic mass, not average MW
        assert desc["MW"] == pytest.approx(180.042, abs=0.01)

    def test_ethanol(self):
        mol = Molecule("CCO")
        assert mol.mol is not None

    def test_caffeine(self):
        mol = Molecule("Cn1c(=O)c2c(ncn2C)n(C)c1=O")
        desc = mol.descriptors()
        assert desc["MW"] == pytest.approx(194.08, abs=0.1)

    def test_invalid_smiles_raises(self):
        with pytest.raises(ValueError, match="Invalid SMILES"):
            Molecule("not_a_valid_smiles_XYZ123!!!")

    def test_empty_string_raises(self):
        with pytest.raises(ValueError, match="non-empty"):
            Molecule("")

    def test_whitespace_only_raises(self):
        with pytest.raises(ValueError, match="non-empty"):
            Molecule("   ")

    def test_none_raises(self):
        with pytest.raises(ValueError, match="None"):
            Molecule(None)


# ---------------------------------------------------------------------------
# Canonical SMILES
# ---------------------------------------------------------------------------


class TestCanonicalSmiles:
    """Test canonical SMILES determinism."""

    def test_deterministic_same_input(self):
        mol1 = Molecule("CC(=O)Oc1ccccc1C(=O)O")
        mol2 = Molecule("CC(=O)Oc1ccccc1C(=O)O")
        assert mol1.canonical_smiles == mol2.canonical_smiles

    def test_deterministic_different_ordering(self):
        """Different valid SMILES for the same molecule should canonicalize identically."""
        # Ethanol written two ways
        mol1 = Molecule("CCO")
        mol2 = Molecule("OCC")
        assert mol1.canonical_smiles == mol2.canonical_smiles

    def test_whitespace_stripped(self):
        mol1 = Molecule("CCO")
        mol2 = Molecule("  CCO  ")
        assert mol1.canonical_smiles == mol2.canonical_smiles


# ---------------------------------------------------------------------------
# Descriptors
# ---------------------------------------------------------------------------


class TestDescriptors:
    """Test descriptor computation."""

    def test_expected_keys(self):
        mol = Molecule("CC(=O)Oc1ccccc1C(=O)O")
        desc = mol.descriptors()
        expected_keys = {
            "MW",
            "logP",
            "TPSA",
            "HBD",
            "HBA",
            "rotatable_bonds",
            "aromatic_rings",
            "heavy_atom_count",
        }
        assert set(desc.keys()) == expected_keys

    def test_aspirin_descriptors_reasonable(self):
        mol = Molecule("CC(=O)Oc1ccccc1C(=O)O")
        desc = mol.descriptors()
        # MW around 180
        assert 170 < desc["MW"] < 190
        # logP around 1.2 (Crippen)
        assert -1.0 < desc["logP"] < 3.0
        # TPSA: aspirin has ester and carboxylic acid
        assert desc["TPSA"] > 0
        # HBD: carboxylic OH
        assert desc["HBD"] >= 1
        # HBA: oxygens
        assert desc["HBA"] >= 2
        # One aromatic ring
        assert desc["aromatic_rings"] == 1
        # 13 heavy atoms
        assert desc["heavy_atom_count"] == 13

    def test_ethanol_no_aromatic_rings(self):
        mol = Molecule("CCO")
        desc = mol.descriptors()
        assert desc["aromatic_rings"] == 0

    def test_descriptor_values_are_numeric(self):
        mol = Molecule("CC(=O)Oc1ccccc1C(=O)O")
        desc = mol.descriptors()
        for key, val in desc.items():
            assert isinstance(val, (int, float)), f"{key} should be numeric"


# ---------------------------------------------------------------------------
# Morgan fingerprint
# ---------------------------------------------------------------------------


class TestMorganFingerprint:
    """Test Morgan fingerprint generation."""

    def test_default_parameters(self):
        mol = Molecule("CCO")
        fp = mol.morgan_fingerprint()
        assert fp.GetNumBits() == 2048

    def test_custom_parameters(self):
        mol = Molecule("CCO")
        fp = mol.morgan_fingerprint(radius=3, n_bits=1024)
        assert fp.GetNumBits() == 1024

    def test_fingerprint_not_all_zeros(self):
        mol = Molecule("CC(=O)Oc1ccccc1C(=O)O")
        fp = mol.morgan_fingerprint()
        assert fp.GetNumOnBits() > 0

    def test_different_molecules_different_fingerprints(self):
        mol1 = Molecule("CCO")
        mol2 = Molecule("CC(=O)Oc1ccccc1C(=O)O")
        fp1 = mol1.morgan_fingerprint()
        fp2 = mol2.morgan_fingerprint()
        # They should differ (not identical bit patterns)
        assert fp1 != fp2


# ---------------------------------------------------------------------------
# Tanimoto similarity
# ---------------------------------------------------------------------------


class TestTanimotoSimilarity:
    """Test Tanimoto similarity computation."""

    def test_self_similarity_is_one(self):
        mol = Molecule("CC(=O)Oc1ccccc1C(=O)O")
        assert mol.tanimoto_similarity(mol) == pytest.approx(1.0)

    def test_identical_molecules_similarity_one(self):
        mol1 = Molecule("CCO")
        mol2 = Molecule("CCO")
        assert mol1.tanimoto_similarity(mol2) == pytest.approx(1.0)

    def test_different_molecules_less_than_one(self):
        mol1 = Molecule("CCO")  # ethanol
        mol2 = Molecule("CC(=O)Oc1ccccc1C(=O)O")  # aspirin
        sim = mol1.tanimoto_similarity(mol2)
        assert 0.0 <= sim < 1.0

    def test_symmetry(self):
        mol1 = Molecule("CCO")
        mol2 = Molecule("CC(=O)Oc1ccccc1C(=O)O")
        assert mol1.tanimoto_similarity(mol2) == pytest.approx(
            mol2.tanimoto_similarity(mol1)
        )

    def test_similar_molecules_higher_score(self):
        """Structurally similar molecules should have higher Tanimoto than dissimilar ones."""
        aspirin = Molecule("CC(=O)Oc1ccccc1C(=O)O")
        salicylic = Molecule("OC(=O)c1ccccc1O")  # salicylic acid (aspirin parent)
        ethanol = Molecule("CCO")

        sim_close = aspirin.tanimoto_similarity(salicylic)
        sim_far = aspirin.tanimoto_similarity(ethanol)
        assert sim_close > sim_far
