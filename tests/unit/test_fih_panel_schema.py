"""Schema tests for validation/data/fih_reference/panel.yaml."""
from __future__ import annotations

from pathlib import Path

import pytest
import yaml

PANEL = Path(__file__).resolve().parents[2] / "validation" / "data" / "fih_reference" / "panel.yaml"


@pytest.fixture(scope="module")
def panel() -> dict:
    return yaml.safe_load(PANEL.read_text())["panel"]


def test_top_level_fields(panel):
    assert panel["name"] == "charon_sprint7_fih"
    assert panel["version"] == 1
    assert isinstance(panel["compounds"], list)
    assert len(panel["compounds"]) >= 22


class TestGoldTier:
    def test_at_least_twelve_gold(self, panel):
        gold = [c for c in panel["compounds"] if c["tier"] == "gold"]
        assert len(gold) >= 12

    def test_gold_has_reference_fih(self, panel):
        for c in panel["compounds"]:
            if c["tier"] == "gold":
                assert "reference_fih_mg" in c
                assert c["reference_fih_mg"] > 0
                assert "source" in c

    def test_gold_has_target_ceff(self, panel):
        """PAD path requires target_ceff_nM for MRSD computation."""
        for c in panel["compounds"]:
            if c["tier"] == "gold":
                assert "target_ceff_nM" in c
                assert c["target_ceff_nM"] > 0


class TestSanityFloorTier:
    def test_has_approved_starting_dose(self, panel):
        for c in panel["compounds"]:
            if c["tier"] == "sanity_floor":
                assert "approved_starting_dose_mg" in c
                assert c["approved_starting_dose_mg"] > 0

    def test_has_target_ceff(self, panel):
        for c in panel["compounds"]:
            if c["tier"] == "sanity_floor":
                assert "target_ceff_nM" in c
                assert c["target_ceff_nM"] > 0


class TestCompoundYamlsExist:
    def test_each_compound_has_property_yaml(self, panel):
        compounds_dir = (
            Path(__file__).resolve().parents[2]
            / "validation"
            / "data"
            / "tier1_obach"
            / "compounds"
        )
        names = {c["name"] for c in panel["compounds"]}
        for name in names:
            assert (compounds_dir / f"{name}.yaml").exists(), f"Missing compound YAML: {name}"


def test_tier_a_compounds_have_permeability():
    """Sprint 11: Tier A oral migration requires Peff or Papp in each
    compound YAML."""
    panel_path = Path(__file__).resolve().parents[2] / "validation" / "data" / "fih_reference" / "panel.yaml"
    compounds_dir = Path(__file__).resolve().parents[2] / "validation" / "data" / "tier1_obach" / "compounds"
    panel = yaml.safe_load(panel_path.read_text())["panel"]
    tier_a_names = {c["name"] for c in panel["compounds"] if c["tier"] == "gold"}
    missing = []
    for name in sorted(tier_a_names):
        yaml_path = compounds_dir / f"{name}.yaml"
        data = yaml.safe_load(yaml_path.read_text())
        perm = data.get("properties", {}).get("permeability", {})
        peff_entry = perm.get("peff_cm_s")
        papp_entry = perm.get("papp_nm_s")
        has_peff = isinstance(peff_entry, dict) and peff_entry.get("value") is not None
        has_papp = isinstance(papp_entry, dict) and papp_entry.get("value") is not None
        if not (has_peff or has_papp):
            missing.append(name)
    assert not missing, (
        f"Tier A compounds missing Peff/Papp: {missing}. "
        f"Sprint 11 requires oral-route permeability data."
    )


def test_tier_a_panel_route_is_oral():
    """Sprint 11: Tier A panel entries simulate oral route."""
    panel_path = Path(__file__).resolve().parents[2] / "validation" / "data" / "fih_reference" / "panel.yaml"
    panel = yaml.safe_load(panel_path.read_text())["panel"]
    wrong_route = [
        c["name"] for c in panel["compounds"]
        if c["tier"] == "gold" and c["route"] != "oral"
    ]
    assert not wrong_route, (
        f"Tier A compounds with non-oral route: {wrong_route}. "
        f"Sprint 11 expects route=oral for all gold-tier entries."
    )
