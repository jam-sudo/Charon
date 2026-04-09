import pytest
import yaml
from pathlib import Path
from charon.core.compound_config import load_compound_config, save_compound_config, apply_overrides
from charon.core.schema import CompoundConfig, PredictedProperty

@pytest.fixture
def sample_config():
    """Create a minimal CompoundConfig for testing."""
    return CompoundConfig(
        name="aspirin",
        smiles="CC(=O)Oc1ccccc1C(=O)O",
        molecular_weight=180.16,
    )

@pytest.fixture
def full_config():
    """Create a CompoundConfig with properties filled in."""
    from charon.core.schema import (
        CompoundProperties, PhysicochemicalProperties, BindingProperties,
        MetabolismProperties, PermeabilityProperties,
    )
    return CompoundConfig(
        name="CPD-001",
        smiles="CC(=O)Oc1ccccc1C(=O)O",
        molecular_weight=180.16,
        properties=CompoundProperties(
            physicochemical=PhysicochemicalProperties(
                logp=PredictedProperty(value=1.19, ci_90_lower=0.95, ci_90_upper=1.43, source="ml_ensemble"),
            ),
            binding=BindingProperties(
                fu_p=PredictedProperty(value=0.23, ci_90_lower=0.18, ci_90_upper=0.29, source="ml_ensemble"),
            ),
            metabolism=MetabolismProperties(
                primary_cyp="CYP2C9",
                clint_uL_min_mg=PredictedProperty(value=15.2, ci_90_lower=8.1, ci_90_upper=22.3, source="ml_ensemble"),
            ),
        ),
    )

class TestYamlRoundTrip:
    def test_save_and_load(self, sample_config, tmp_path):
        """Save -> Load should produce equivalent config."""
        path = tmp_path / "test_compound.yaml"
        save_compound_config(sample_config, path)
        loaded = load_compound_config(path)
        assert loaded.name == sample_config.name
        assert loaded.smiles == sample_config.smiles
        assert loaded.molecular_weight == sample_config.molecular_weight

    def test_full_config_roundtrip(self, full_config, tmp_path):
        """Full config with properties should survive roundtrip."""
        path = tmp_path / "test_full.yaml"
        save_compound_config(full_config, path)
        loaded = load_compound_config(path)
        assert loaded.properties.binding.fu_p.value == 0.23
        assert loaded.properties.metabolism.primary_cyp == "CYP2C9"

    def test_file_not_found(self):
        with pytest.raises(FileNotFoundError):
            load_compound_config("/nonexistent/path.yaml")

class TestOverrides:
    def test_override_fu_p(self, full_config):
        """Override fu_p with experimental value."""
        overridden = apply_overrides(full_config, {
            "fu_p": {"value": 0.21, "source": "experimental", "method": "equilibrium_dialysis"}
        })
        assert overridden.properties.binding.fu_p.value == 0.21
        assert overridden.properties.binding.fu_p.source == "experimental"
        # Original should be unchanged
        assert full_config.properties.binding.fu_p.value == 0.23

    def test_unknown_override_raises(self, full_config):
        """Unknown override key should raise KeyError."""
        with pytest.raises(KeyError):
            apply_overrides(full_config, {"unknown_param": {"value": 1.0}})
