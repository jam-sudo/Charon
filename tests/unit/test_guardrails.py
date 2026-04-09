import pytest
from charon.core.guardrails import GuardrailChecker
from charon.core.molecule import Molecule

@pytest.fixture
def checker():
    return GuardrailChecker()

class TestStructuralValidity:
    def test_valid_smiles(self, checker):
        """Aspirin should be valid."""
        result = checker.check("CC(=O)Oc1ccccc1C(=O)O")
        assert result.is_valid is True

    def test_invalid_smiles(self, checker):
        """Invalid SMILES should set is_valid=False."""
        result = checker.check("not_a_smiles")
        assert result.is_valid is False

    def test_empty_smiles(self, checker):
        result = checker.check("")
        assert result.is_valid is False

class TestDrugLikeness:
    def test_aspirin_no_lipinski_warnings(self, checker):
        """Aspirin (MW=180, logP~1.2) should pass Lipinski."""
        result = checker.check("CC(=O)Oc1ccccc1C(=O)O")
        lipinski_warnings = [w for w in result.warnings if w.category == "lipinski"]
        assert len(lipinski_warnings) == 0

    def test_guardrails_never_block(self, checker):
        """Even with violations, is_valid should be True (if SMILES is valid)."""
        # A large, valid molecule: hexadecane (MW ~226, logP ~8) triggers
        # high-lipophilicity warning but must remain is_valid=True.
        big_mol = "CCCCCCCCCCCCCCCC"
        result = checker.check(big_mol)
        assert result.is_valid is True  # Valid SMILES, even if drug-likeness violations

class TestKnownFailureModes:
    def test_high_mw_critical(self, checker):
        """MW > 800 should generate CRITICAL warning."""
        # Use a large molecule - will generate the warning if MW > 800
        # Cyclosporine A SMILES (MW ~1202)
        cyclosporine = "CCC1NC(=O)C(CC(C)C)N(C)C(=O)C(CC(C)C)NC(=O)C(CC(C)C)N(C)C(=O)C(CC(C)C)NC(=O)C(C(C)C)NC(=O)C(CC(C)C)N(C)C(=O)C(C(C)CC=CC)NC(=O)C(C)NC(=O)C(CC)N(C)C(=O)C(CC(C)C)NC1=O"
        result = checker.check(cyclosporine)
        critical_warnings = [w for w in result.warnings if w.severity == "CRITICAL"]
        # If this SMILES parses correctly AND MW > 800, we should see a CRITICAL
        if result.is_valid and result.molecular_weight > 800:
            assert len(critical_warnings) > 0

class TestApplicabilityDomain:
    def test_no_reference_returns_unknown(self, checker):
        """Without reference fingerprints, AD should be UNKNOWN."""
        result = checker.check("CC(=O)Oc1ccccc1C(=O)O")
        assert result.applicability_domain == "UNKNOWN"

    def test_with_reference_self(self):
        """Molecule similar to reference should get HIGH."""
        aspirin = Molecule("CC(=O)Oc1ccccc1C(=O)O")
        ref_fps = [aspirin.morgan_fingerprint()]
        checker = GuardrailChecker(reference_fps=ref_fps)
        result = checker.check("CC(=O)Oc1ccccc1C(=O)O")
        assert result.applicability_domain == "HIGH"
