"""Tests for the Obach 1999 panel loader."""

from __future__ import annotations

import textwrap
from pathlib import Path

import pytest
import yaml

# These imports target the loader that will live in
# validation/benchmarks/layer2_human_pk.py after Task 10.
import sys
REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from validation.benchmarks.layer2_human_pk import (  # noqa: E402
    PanelEntry,
    load_panel,
)


@pytest.fixture
def tmp_panel(tmp_path: Path) -> Path:
    """Build a minimal self-contained panel dir on disk."""
    root = tmp_path / "tier1_obach"
    (root / "compounds").mkdir(parents=True)

    theo_yaml = textwrap.dedent("""
        name: theophylline
        smiles: "Cn1c(=O)c2[nH]cnc2n(C)c1=O"
        molecular_weight: 180.17
        source: experimental
        properties:
          physicochemical:
            logp: {value: -0.02, source: experimental, method: "Obach 1999 Table 2"}
            compound_type: neutral
          binding:
            fu_p: {value: 0.60, source: experimental, unit: fraction}
            fu_inc: {value: 1.0, source: experimental, unit: fraction}
            bp_ratio: {value: 0.85, source: experimental, unit: ratio}
          metabolism:
            clint_uL_min_mg: {value: 1.8, source: experimental, unit: uL/min/mg}
          renal:
            clrenal_L_h: {value: 0.1, source: experimental, unit: L/h}
    """).lstrip()
    (root / "compounds" / "theophylline.yaml").write_text(theo_yaml)

    panel_yaml = textwrap.dedent("""
        name: "Test Panel"
        source: "Obach 1999 (trimmed)"
        default_duration_h: 168.0
        compounds:
          - key: theophylline
            compound_file: compounds/theophylline.yaml
            route: iv_bolus
            dose_mg: 100.0
            duration_h: 168.0
            observed:
              cl_L_h: 2.9
              vss_L: 35.0
              t_half_h: 8.0
            obach_table_row: 42
            strict_targets: true
            notes: "Primary gate"
    """).lstrip()
    (root / "panel.yaml").write_text(panel_yaml)

    return root / "panel.yaml"


class TestLoadPanel:
    def test_single_entry_loads(self, tmp_panel: Path):
        entries = load_panel(tmp_panel)
        assert len(entries) == 1
        e = entries[0]
        assert isinstance(e, PanelEntry)
        assert e.key == "theophylline"
        assert e.compound.name == "theophylline"
        assert e.route == "iv_bolus"
        assert e.dose_mg == 100.0
        assert e.duration_h == 168.0
        assert e.strict_targets is True
        assert e.observed == {"cl_L_h": 2.9, "vss_L": 35.0, "t_half_h": 8.0}

    def test_compound_file_relative_to_panel_dir(self, tmp_panel: Path):
        entries = load_panel(tmp_panel)
        assert entries[0].compound.properties.physicochemical.logp.value == -0.02

    def test_missing_compound_file_raises(self, tmp_path: Path):
        root = tmp_path / "bad_panel"
        root.mkdir()
        panel_yaml = textwrap.dedent("""
            name: "Broken"
            source: "test"
            default_duration_h: 168.0
            compounds:
              - key: ghost
                compound_file: compounds/ghost.yaml
                route: iv_bolus
                dose_mg: 1.0
                duration_h: 168.0
                observed: {cl_L_h: 1.0, vss_L: 1.0, t_half_h: 1.0}
                strict_targets: false
        """).lstrip()
        (root / "panel.yaml").write_text(panel_yaml)
        with pytest.raises(FileNotFoundError, match="ghost.yaml"):
            load_panel(root / "panel.yaml")

    def test_default_duration_fallback(self, tmp_path: Path):
        """Entry without duration_h uses panel default_duration_h."""
        root = tmp_path / "panel"
        (root / "compounds").mkdir(parents=True)
        (root / "compounds" / "theo.yaml").write_text(textwrap.dedent("""
            name: theo
            smiles: "Cn1c(=O)c2[nH]cnc2n(C)c1=O"
            molecular_weight: 180.17
            source: experimental
            properties:
              physicochemical:
                logp: {value: -0.02, source: experimental}
              binding:
                fu_p: {value: 0.6, source: experimental, unit: fraction}
                fu_inc: {value: 1.0, source: experimental, unit: fraction}
                bp_ratio: {value: 0.85, source: experimental, unit: ratio}
              metabolism:
                clint_uL_min_mg: {value: 1.8, source: experimental, unit: uL/min/mg}
              renal:
                clrenal_L_h: {value: 0.1, source: experimental, unit: L/h}
        """).lstrip())
        (root / "panel.yaml").write_text(textwrap.dedent("""
            name: "Fallback"
            source: "test"
            default_duration_h: 72.0
            compounds:
              - key: theo
                compound_file: compounds/theo.yaml
                route: iv_bolus
                dose_mg: 100.0
                observed: {cl_L_h: 2.9, vss_L: 35.0, t_half_h: 8.0}
                strict_targets: false
        """).lstrip())
        entries = load_panel(root / "panel.yaml")
        assert entries[0].duration_h == 72.0

    def test_observed_keys_required(self, tmp_path: Path):
        """Missing observed.cl_L_h → clear error."""
        root = tmp_path / "panel"
        (root / "compounds").mkdir(parents=True)
        (root / "compounds" / "c.yaml").write_text(textwrap.dedent("""
            name: c
            smiles: "C"
            molecular_weight: 16.0
            source: experimental
            properties:
              physicochemical:
                logp: {value: 0.5, source: experimental}
              binding:
                fu_p: {value: 0.9, source: experimental, unit: fraction}
                fu_inc: {value: 1.0, source: experimental, unit: fraction}
                bp_ratio: {value: 1.0, source: experimental, unit: ratio}
              metabolism:
                clint_uL_min_mg: {value: 1.0, source: experimental, unit: uL/min/mg}
              renal:
                clrenal_L_h: {value: 0.0, source: experimental, unit: L/h}
        """).lstrip())
        (root / "panel.yaml").write_text(textwrap.dedent("""
            name: x
            source: x
            default_duration_h: 72.0
            compounds:
              - key: c
                compound_file: compounds/c.yaml
                route: iv_bolus
                dose_mg: 1.0
                observed: {vss_L: 1.0, t_half_h: 1.0}
                strict_targets: false
        """).lstrip())
        with pytest.raises(KeyError, match="cl_L_h"):
            load_panel(root / "panel.yaml")


class TestWithoutKpOverrides:
    def test_strips_empirical_kp(self):
        from charon.core.schema import (
            CompoundProperties, DistributionProperties,
            PhysicochemicalProperties, PredictedProperty,
        )
        from validation.benchmarks.layer2_human_pk import _without_kp_overrides

        original = CompoundProperties(
            physicochemical=PhysicochemicalProperties(
                logp=PredictedProperty(value=3.89, source="experimental"),
            ),
            distribution=DistributionProperties(
                empirical_kp_by_tissue={
                    "adipose": PredictedProperty(
                        value=10.0, source="literature", method="test"
                    ),
                }
            ),
        )
        stripped = _without_kp_overrides(original)
        assert stripped.distribution.empirical_kp_by_tissue is None
        # Other fields preserved
        assert stripped.physicochemical.logp.value == 3.89

    def test_no_override_case_is_noop(self):
        from charon.core.schema import CompoundProperties, PhysicochemicalProperties, PredictedProperty
        from validation.benchmarks.layer2_human_pk import _without_kp_overrides

        original = CompoundProperties(
            physicochemical=PhysicochemicalProperties(
                logp=PredictedProperty(value=-0.02, source="experimental"),
            ),
        )
        stripped = _without_kp_overrides(original)
        assert stripped.distribution.empirical_kp_by_tissue is None
        assert stripped.physicochemical.logp.value == -0.02

    def test_stripped_distribution_is_same_type(self):
        """The strip must not degrade to a parent type."""
        from charon.core.schema import (
            CompoundProperties, DistributionProperties,
            PredictedProperty,
        )
        from validation.benchmarks.layer2_human_pk import _without_kp_overrides

        original = CompoundProperties(
            distribution=DistributionProperties(
                empirical_kp_by_tissue={
                    "adipose": PredictedProperty(
                        value=10.0, source="literature", method="test"
                    ),
                }
            ),
        )
        stripped = _without_kp_overrides(original)
        assert isinstance(stripped.distribution, DistributionProperties)
        assert stripped.distribution.empirical_kp_by_tissue is None
