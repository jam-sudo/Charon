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
    assert len(panel["compounds"]) >= 15


class TestGoldTier:
    def test_at_least_five_gold(self, panel):
        gold = [c for c in panel["compounds"] if c["tier"] == "gold"]
        assert len(gold) >= 5

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
