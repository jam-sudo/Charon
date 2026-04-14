"""End-to-end Pipeline + uncertainty quantification tests."""
from __future__ import annotations
import sys
from pathlib import Path
import pytest
import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from charon.core.schema import CompoundConfig, DoseProjectionConfig, UncertaintyConfig
from charon.pipeline import Pipeline


def _load_midazolam():
    path = REPO_ROOT / "validation" / "data" / "tier1_obach" / "compounds" / "midazolam.yaml"
    with path.open() as f:
        return CompoundConfig(**yaml.safe_load(f))


class TestPipelineUncertainty:
    def test_returns_uncertainty_result(self):
        pipe = Pipeline(
            compound=_load_midazolam(), route="oral", dose_mg=5.0,
            dose_projection=DoseProjectionConfig(noael_mg_kg=50.0, noael_species="rat"),
            uncertainty=UncertaintyConfig(n_samples=30),
        )
        result = pipe.run()
        assert result.uncertainty is not None
        assert result.uncertainty.point_estimate_mg > 0
        # HED-only dose projection yields identical MRSD for all samples
        # (NOAEL/Km not sampled), so lower == upper is valid.
        assert result.uncertainty.ci_90_lower_mg <= result.uncertainty.ci_90_upper_mg

    def test_no_uncertainty_returns_none(self):
        pipe = Pipeline(
            compound=_load_midazolam(), route="oral", dose_mg=5.0,
            dose_projection=DoseProjectionConfig(noael_mg_kg=50.0, noael_species="rat"),
        )
        result = pipe.run()
        assert result.uncertainty is None

    def test_uncertainty_without_dose_projection_raises(self):
        with pytest.raises(ValueError, match="dose_projection"):
            Pipeline(
                compound=_load_midazolam(), route="oral", dose_mg=5.0,
                uncertainty=UncertaintyConfig(n_samples=30),
            ).run()

    def test_confidence_is_valid(self):
        pipe = Pipeline(
            compound=_load_midazolam(), route="oral", dose_mg=5.0,
            dose_projection=DoseProjectionConfig(noael_mg_kg=50.0, noael_species="rat"),
            uncertainty=UncertaintyConfig(n_samples=30),
        )
        result = pipe.run()
        assert result.uncertainty.confidence in ("HIGH", "MEDIUM", "LOW")

    def test_sensitivity_has_params(self):
        pipe = Pipeline(
            compound=_load_midazolam(), route="oral", dose_mg=5.0,
            dose_projection=DoseProjectionConfig(noael_mg_kg=50.0, noael_species="rat"),
            uncertainty=UncertaintyConfig(n_samples=30),
        )
        result = pipe.run()
        assert len(result.uncertainty.sensitivity) >= 3

    def test_recommendation_string(self):
        pipe = Pipeline(
            compound=_load_midazolam(), route="oral", dose_mg=5.0,
            dose_projection=DoseProjectionConfig(noael_mg_kg=50.0, noael_species="rat"),
            uncertainty=UncertaintyConfig(n_samples=30),
        )
        result = pipe.run()
        assert isinstance(result.uncertainty.recommendation, str)
        assert len(result.uncertainty.recommendation) > 10

    def test_metadata_includes_uncertainty(self):
        pipe = Pipeline(
            compound=_load_midazolam(), route="oral", dose_mg=5.0,
            dose_projection=DoseProjectionConfig(noael_mg_kg=50.0, noael_species="rat"),
            uncertainty=UncertaintyConfig(n_samples=30),
        )
        result = pipe.run()
        assert "uncertainty_ci_90" in result.metadata

    def test_iv_with_uncertainty(self):
        pipe = Pipeline(
            compound=_load_midazolam(), route="iv_bolus", dose_mg=5.0,
            dose_projection=DoseProjectionConfig(noael_mg_kg=50.0, noael_species="rat"),
            uncertainty=UncertaintyConfig(n_samples=20),
        )
        result = pipe.run()
        assert result.uncertainty is not None
        assert result.uncertainty.point_estimate_mg > 0
