"""Unit tests for PBPK topology loading."""

import math
from collections import OrderedDict
from pathlib import Path

import pytest

from charon.pbpk.topology import (
    PBPKTopology,
    TissueNode,
    load_species_topology,
    PORTAL_TISSUES,
)


@pytest.fixture
def human_topology() -> PBPKTopology:
    return load_species_topology("human")


class TestLoadHumanTopology:
    def test_returns_pbpk_topology(self, human_topology):
        assert isinstance(human_topology, PBPKTopology)
        assert human_topology.species == "human"

    def test_body_weight(self, human_topology):
        assert human_topology.body_weight_kg == pytest.approx(70.0)

    def test_cardiac_output(self, human_topology):
        assert human_topology.cardiac_output_L_h == pytest.approx(390.0)

    def test_hematocrit(self, human_topology):
        assert human_topology.hematocrit == pytest.approx(0.45)

    def test_blood_pool_volumes(self, human_topology):
        assert human_topology.venous_volume_L == pytest.approx(3.7)
        assert human_topology.arterial_volume_L == pytest.approx(1.5)

    def test_tissues_are_ordered(self, human_topology):
        assert isinstance(human_topology.tissues, OrderedDict)
        # Deterministic ordering: YAML insertion order preserved.
        names = list(human_topology.tissues.keys())
        assert names[0] == "lung"  # lung is first in human.yaml
        assert "liver" in names
        assert "kidney" in names

    def test_has_15_tissues(self, human_topology):
        assert len(human_topology.tissues) == 15

    def test_lung_has_full_cardiac_output(self, human_topology):
        lung = human_topology.tissues["lung"]
        # Lung receives full CO (pulmonary circulation)
        assert lung.blood_flow_L_h == pytest.approx(
            human_topology.cardiac_output_L_h
        )

    def test_liver_flow_equals_arterial_plus_portal(self, human_topology):
        liver_q = human_topology.tissues["liver"].blood_flow_L_h
        portal_q = sum(
            human_topology.tissues[p].blood_flow_L_h for p in PORTAL_TISSUES
        )
        ha_q = human_topology.hepatic_artery_L_h
        assert liver_q == pytest.approx(ha_q + portal_q, rel=1e-6)

    def test_hepatic_artery_is_nonnegative(self, human_topology):
        # Q_HA = Q_liver_total - Q_portal must be >= 0 for consistency.
        assert human_topology.hepatic_artery_L_h >= 0

    def test_hepatic_artery_matches_yaml_numbers(self, human_topology):
        # human.yaml comment: liver total = 0.065 arterial + 0.19 portal CO
        expected_ha = 0.065 * human_topology.cardiac_output_L_h
        assert human_topology.hepatic_artery_L_h == pytest.approx(expected_ha, rel=1e-6)

    def test_systemic_flows_sum_to_cardiac_output(self, human_topology):
        """All arterial outflow paths must sum to cardiac output.

        Outflow categories:
          - hepatic artery to liver (Q_HA)
          - portal tissues (spleen, gut_wall, pancreas) — drain into liver
          - non-portal non-liver non-lung systemic tissues — drain into venous
        """
        portal_sum = sum(
            human_topology.tissues[p].blood_flow_L_h for p in PORTAL_TISSUES
        )
        non_portal_systemic = sum(
            node.blood_flow_L_h
            for name, node in human_topology.tissues.items()
            if name not in {"lung", "liver"} and name not in PORTAL_TISSUES
        )
        total = human_topology.hepatic_artery_L_h + portal_sum + non_portal_systemic
        assert total == pytest.approx(
            human_topology.cardiac_output_L_h, rel=1e-3
        )

    def test_tissue_nodes_have_composition(self, human_topology):
        for node in human_topology.tissues.values():
            comp = node.composition
            assert math.isfinite(comp.fn)
            assert math.isfinite(comp.fp)
            assert math.isfinite(comp.fw)
            assert math.isfinite(comp.pH)

    def test_portal_tissues_drain_to_liver(self, human_topology):
        for name in PORTAL_TISSUES:
            node = human_topology.tissues[name]
            assert node.drains_to == "liver"

    def test_non_portal_tissues_drain_to_venous(self, human_topology):
        for name, node in human_topology.tissues.items():
            if name in PORTAL_TISSUES or name == "lung" or name == "liver":
                continue
            assert node.drains_to == "venous"

    def test_lung_drains_to_arterial(self, human_topology):
        assert human_topology.tissues["lung"].drains_to == "arterial"

    def test_liver_drains_to_venous(self, human_topology):
        assert human_topology.tissues["liver"].drains_to == "venous"

    def test_plasma_composition_present(self, human_topology):
        assert human_topology.plasma_composition.pH == pytest.approx(7.40)


class TestLoadMissingSpecies:
    def test_raises_file_not_found(self):
        with pytest.raises(FileNotFoundError):
            load_species_topology("unicorn")


class TestLoadFromPath:
    def test_load_from_explicit_path(self, tmp_path):
        # Build a minimal YAML and ensure the loader can ingest it.
        yaml_content = """\
species:
  name: mini
  body_weight_kg: 10.0
  cardiac_output_L_h: 50.0
  hematocrit: 0.4
  liver:
    weight_g: 200.0
    volume_L: 0.2
    blood_flow_fraction_CO: 0.2
  kidney:
    weight_g: 40.0
    volume_L: 0.04
    blood_flow_fraction_CO: 0.1
    gfr_mL_min: 10.0
  bsa_scaling:
    km: 10.0
  blood:
    venous_volume_L: 0.5
    arterial_volume_L: 0.2
    portal_vein_volume_L: 0.01
  plasma:
    fn: 0.0023
    fp: 0.0199
    fw: 0.945
    pH: 7.40
  tissues:
    lung:
      volume_L: 0.05
      blood_flow_fraction: 1.0
      composition: { fn: 0.003, fp: 0.013, fw: 0.81, pH: 7.0 }
    liver:
      volume_L: 0.2
      blood_flow_fraction: 0.2
      composition: { fn: 0.035, fp: 0.025, fw: 0.75, pH: 7.0 }
    kidney:
      volume_L: 0.04
      blood_flow_fraction: 0.1
      composition: { fn: 0.012, fp: 0.024, fw: 0.78, pH: 7.0 }
    spleen:
      volume_L: 0.02
      blood_flow_fraction: 0.05
      composition: { fn: 0.008, fp: 0.011, fw: 0.79, pH: 7.0 }
    gut_wall:
      volume_L: 0.05
      blood_flow_fraction: 0.1
      composition: { fn: 0.016, fp: 0.018, fw: 0.72, pH: 7.0 }
    pancreas:
      volume_L: 0.01
      blood_flow_fraction: 0.02
      composition: { fn: 0.035, fp: 0.025, fw: 0.75, pH: 7.0 }
    muscle:
      volume_L: 5.0
      blood_flow_fraction: 0.5
      composition: { fn: 0.024, fp: 0.007, fw: 0.76, pH: 7.0 }
"""
        yaml_path = tmp_path / "mini.yaml"
        yaml_path.write_text(yaml_content)
        topo = load_species_topology("mini", path=yaml_path)
        assert topo.species == "mini"
        # Liver = 0.2 * 50 = 10 L/h
        assert topo.tissues["liver"].blood_flow_L_h == pytest.approx(10.0)
        # Portal = (0.05 + 0.1 + 0.02) * 50 = 8.5 L/h
        # HA = 10 - 8.5 = 1.5 L/h
        assert topo.hepatic_artery_L_h == pytest.approx(1.5, rel=1e-6)


class TestTopologyConsistencyChecks:
    def test_rejects_liver_flow_less_than_portal(self, tmp_path):
        """Portal tissues total more than liver flow → inconsistent YAML."""
        yaml_content = """\
species:
  name: broken
  body_weight_kg: 10.0
  cardiac_output_L_h: 50.0
  hematocrit: 0.4
  liver:
    weight_g: 200.0
    volume_L: 0.2
    blood_flow_fraction_CO: 0.05
  kidney:
    weight_g: 40.0
    volume_L: 0.04
    blood_flow_fraction_CO: 0.1
    gfr_mL_min: 10.0
  bsa_scaling:
    km: 10.0
  blood:
    venous_volume_L: 0.5
    arterial_volume_L: 0.2
    portal_vein_volume_L: 0.01
  plasma:
    fn: 0.0023
    fp: 0.0199
    fw: 0.945
    pH: 7.40
  tissues:
    lung:
      volume_L: 0.05
      blood_flow_fraction: 1.0
      composition: { fn: 0.003, fp: 0.013, fw: 0.81, pH: 7.0 }
    liver:
      volume_L: 0.2
      blood_flow_fraction: 0.05
      composition: { fn: 0.035, fp: 0.025, fw: 0.75, pH: 7.0 }
    kidney:
      volume_L: 0.04
      blood_flow_fraction: 0.1
      composition: { fn: 0.012, fp: 0.024, fw: 0.78, pH: 7.0 }
    spleen:
      volume_L: 0.02
      blood_flow_fraction: 0.2
      composition: { fn: 0.008, fp: 0.011, fw: 0.79, pH: 7.0 }
    gut_wall:
      volume_L: 0.05
      blood_flow_fraction: 0.3
      composition: { fn: 0.016, fp: 0.018, fw: 0.72, pH: 7.0 }
    pancreas:
      volume_L: 0.01
      blood_flow_fraction: 0.05
      composition: { fn: 0.035, fp: 0.025, fw: 0.75, pH: 7.0 }
"""
        yaml_path = tmp_path / "broken.yaml"
        yaml_path.write_text(yaml_content)
        with pytest.raises(ValueError, match="portal"):
            load_species_topology("broken", path=yaml_path)
