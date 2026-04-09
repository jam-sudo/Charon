import pytest
from pathlib import Path
from charon.core.config_manager import (
    create_run_config, config_to_yaml, config_from_yaml,
    diff_configs, hash_config,
)
from charon.core.schema import CompoundConfig, PipelineConfig, RunConfig

@pytest.fixture
def compound():
    return CompoundConfig(
        name="test-compound",
        smiles="CCO",
        molecular_weight=46.07,
    )

class TestCreateRunConfig:
    def test_basic_creation(self, compound):
        config = create_run_config(compound)
        assert config.compound.name == "test-compound"
        assert config.pipeline is not None

class TestYamlRoundTrip:
    def test_save_and_load(self, compound, tmp_path):
        config = create_run_config(compound)
        path = tmp_path / "run_config.yaml"
        config_to_yaml(config, path)
        loaded = config_from_yaml(path)
        assert loaded.compound.smiles == "CCO"

class TestDiff:
    def test_identical_configs(self, compound):
        c1 = create_run_config(compound)
        c2 = create_run_config(compound)
        assert diff_configs(c1, c2) == {}

    def test_different_compounds(self):
        c1 = create_run_config(CompoundConfig(name="A", smiles="CCO", molecular_weight=46.07))
        c2 = create_run_config(CompoundConfig(name="B", smiles="CCO", molecular_weight=46.07))
        d = diff_configs(c1, c2)
        assert len(d) > 0

class TestHash:
    def test_deterministic(self, compound):
        c1 = create_run_config(compound)
        c2 = create_run_config(compound)
        assert hash_config(c1) == hash_config(c2)

    def test_different_configs_different_hash(self):
        c1 = create_run_config(CompoundConfig(name="A", smiles="CCO", molecular_weight=46.07))
        c2 = create_run_config(CompoundConfig(name="B", smiles="CCO", molecular_weight=46.07))
        assert hash_config(c1) != hash_config(c2)

class TestFrozenConfig:
    def test_cannot_mutate(self, compound):
        """RunConfig(frozen=True) should prevent mutation."""
        config = create_run_config(compound)
        with pytest.raises(Exception):  # pydantic ValidationError on mutation
            config.compound = CompoundConfig(name="X", smiles="C", molecular_weight=16.0)
