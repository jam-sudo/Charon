"""End-to-end Pipeline + FIH dose projection tests."""
from __future__ import annotations

import sys
from pathlib import Path

import pytest
import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from charon.core.schema import CompoundConfig, DoseProjectionConfig
from charon.pipeline import Pipeline


def _load_midazolam():
    path = REPO_ROOT / "validation" / "data" / "tier1_obach" / "compounds" / "midazolam.yaml"
    with path.open() as f:
        return CompoundConfig(**yaml.safe_load(f))


class TestPipelineWithDoseProjection:
    def test_oral_with_hed(self):
        pipe = Pipeline(
            compound=_load_midazolam(), route="oral", dose_mg=5.0,
            dose_projection=DoseProjectionConfig(noael_mg_kg=50.0, noael_species="rat"),
        )
        result = pipe.run()
        assert result.dose_recommendation is not None
        assert result.dose_recommendation.hed is not None
        assert result.dose_recommendation.mrsd_mg > 0
        assert result.dose_recommendation.mrsd_mg == pytest.approx(58.65, rel=0.01)

    def test_oral_with_all_three(self):
        pipe = Pipeline(
            compound=_load_midazolam(), route="oral", dose_mg=5.0,
            dose_projection=DoseProjectionConfig(
                noael_mg_kg=50.0, noael_species="rat",
                target_kd_nM=10.0, target_ceff_nM=100.0,
            ),
        )
        result = pipe.run()
        rec = result.dose_recommendation
        assert rec.hed is not None and rec.mabel is not None and rec.pad is not None
        all_mrsd = [rec.hed.mrsd_mg, rec.mabel.mrsd_mg, rec.pad.mrsd_mg]
        assert rec.mrsd_mg <= min(all_mrsd) + 0.01

    def test_no_dose_projection_returns_none(self):
        pipe = Pipeline(compound=_load_midazolam(), route="oral", dose_mg=5.0)
        result = pipe.run()
        assert result.dose_recommendation is None

    def test_iv_with_hed(self):
        pipe = Pipeline(
            compound=_load_midazolam(), route="iv_bolus", dose_mg=5.0,
            dose_projection=DoseProjectionConfig(noael_mg_kg=50.0, noael_species="rat"),
        )
        result = pipe.run()
        assert result.dose_recommendation is not None
        assert result.dose_recommendation.hed is not None

    def test_rationale_is_string(self):
        pipe = Pipeline(
            compound=_load_midazolam(), route="oral", dose_mg=5.0,
            dose_projection=DoseProjectionConfig(noael_mg_kg=50.0, noael_species="rat", target_kd_nM=10.0),
        )
        result = pipe.run()
        assert isinstance(result.dose_recommendation.rationale, str)
        assert "HED" in result.dose_recommendation.rationale

    def test_metadata_includes_dose(self):
        pipe = Pipeline(
            compound=_load_midazolam(), route="oral", dose_mg=5.0,
            dose_projection=DoseProjectionConfig(noael_mg_kg=50.0, noael_species="rat"),
        )
        result = pipe.run()
        assert "mrsd_mg" in result.metadata
