"""Sprint 9 smoke tests: the three newly-curated compound YAMLs must load
cleanly via charon.core.compound_config.load_compound_config, and each
must declare at least one non-zero clearance path (hepatic or renal).
"""
from __future__ import annotations

from pathlib import Path

import pytest

from charon.core.compound_config import load_compound_config

REPO_ROOT = Path(__file__).resolve().parents[2]
COMPOUNDS_DIR = REPO_ROOT / "validation" / "data" / "tier1_obach" / "compounds"
NEW_NAMES = ["acetaminophen", "lisinopril", "atorvastatin"]


@pytest.mark.parametrize("name", NEW_NAMES)
def test_new_compound_yaml_loads(name: str):
    path = COMPOUNDS_DIR / f"{name}.yaml"
    assert path.exists(), f"Missing {path}"
    compound = load_compound_config(path)
    assert compound.name == name
    assert compound.molecular_weight > 0


@pytest.mark.parametrize("name", NEW_NAMES)
def test_new_compound_has_clearance_path(name: str):
    compound = load_compound_config(COMPOUNDS_DIR / f"{name}.yaml")
    clint = compound.properties.metabolism.clint_uL_min_mg
    clren = compound.properties.renal.clrenal_L_h
    has_hepatic = clint is not None and clint.value > 0
    has_renal = clren is not None and clren.value > 0
    assert has_hepatic or has_renal, (
        f"{name} must declare at least one non-zero clearance path; "
        f"clint={clint}, clrenal={clren}"
    )


@pytest.mark.parametrize("name", NEW_NAMES)
def test_new_compound_binding_present(name: str):
    compound = load_compound_config(COMPOUNDS_DIR / f"{name}.yaml")
    fup = compound.properties.binding.fu_p
    bp = compound.properties.binding.bp_ratio
    assert fup is not None and 0.0 < fup.value <= 1.0
    assert bp is not None and bp.value > 0
