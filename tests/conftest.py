"""Shared fixtures for Charon test suite."""

from __future__ import annotations

import pytest

from charon.core.parameter_bridge import ParameterBridge
from charon.core.schema import (
    CompoundConfig,
    CompoundProperties,
    PhysicochemicalProperties,
    PredictedProperty,
    SaltForm,
)


@pytest.fixture
def bridge() -> ParameterBridge:
    """A fresh ParameterBridge instance."""
    return ParameterBridge()


@pytest.fixture
def reference_drugs() -> dict:
    """Reference drug data for integration/regression tests.

    Values are approximate consensus values drawn from literature.
    Keys: drug name -> dict of PK-relevant properties.
    """
    return {
        "midazolam": {
            "smiles": "Clc1ccc2c(c1)C(=NCC(=O)n2c3ccccc3F)c4ccccc4",
            "mw": 325.77,
            "logp": 3.89,
            "fu_p": 0.03,
            "clint_hlm": 140.0,  # uL/min/mg protein
            "bp_ratio": 0.54,
        },
        "warfarin": {
            "smiles": "CC(=O)CC(c1ccccc1)c2c(O)c3ccccc3oc2=O",
            "mw": 308.33,
            "logp": 2.7,
            "fu_p": 0.005,
            "clint_hlm": 3.0,
            "bp_ratio": 0.58,
        },
        "caffeine": {
            "smiles": "Cn1c(=O)c2c(ncn2C)n(C)c1=O",
            "mw": 194.19,
            "logp": -0.07,
            "fu_p": 0.65,
            "clint_hlm": 15.0,
            "bp_ratio": 1.0,
        },
        "aspirin": {
            "smiles": "CC(=O)Oc1ccccc1C(=O)O",
            "mw": 180.16,
            "logp": 1.24,
            "fu_p": 0.50,
            "clint_hlm": 680.0,
            "bp_ratio": 0.7,
        },
    }


@pytest.fixture
def sample_compound_config() -> CompoundConfig:
    """A minimal but complete CompoundConfig for testing."""
    return CompoundConfig(
        name="test-compound",
        smiles="CC(=O)Oc1ccccc1C(=O)O",
        molecular_weight=180.16,
        source="predicted",
        properties=CompoundProperties(
            physicochemical=PhysicochemicalProperties(
                logp=PredictedProperty(value=1.24, source="ml_ensemble"),
            ),
        ),
    )


@pytest.fixture
def compound_config_with_salt() -> CompoundConfig:
    """CompoundConfig that triggers auto salt_factor calculation."""
    return CompoundConfig(
        name="salt-test",
        smiles="CC(=O)Oc1ccccc1C(=O)O",
        molecular_weight=180.16,
        salt_form=SaltForm(name="sodium salt", mw_salt=202.14),
    )
