"""Tests for Tier 3 classification sampling in Layer 4 uncertainty."""
from __future__ import annotations

import math

import numpy as np
import pytest

from charon.uncertainty.sampling import _sample_classification_clint


@pytest.fixture
def rng():
    return np.random.default_rng(42)


class TestBucketDistribution:
    def test_empirical_fractions_match_probs(self, rng):
        probs = {"low": 0.6, "med": 0.3, "high": 0.1}
        samples = _sample_classification_clint(probs, n=20000, rng=rng)
        assert samples.shape == (20000,)
        lo_frac = float(np.mean(samples < 10.0))
        med_frac = float(np.mean((samples >= 10.0) & (samples < 50.0)))
        hi_frac = float(np.mean(samples >= 50.0))
        assert lo_frac == pytest.approx(0.6, abs=0.02)
        assert med_frac == pytest.approx(0.3, abs=0.02)
        assert hi_frac == pytest.approx(0.1, abs=0.02)


class TestLogUniformWithinBucket:
    def test_values_lie_in_physical_range(self, rng):
        probs = {"low": 0.34, "med": 0.33, "high": 0.33}
        samples = _sample_classification_clint(probs, n=5000, rng=rng)
        assert np.all(samples >= 0.1)
        assert np.all(samples <= 1000.0)

    def test_log_uniform_within_low_bucket(self, rng):
        probs = {"low": 1.0, "med": 0.0, "high": 0.0}
        samples = _sample_classification_clint(probs, n=5000, rng=rng)
        assert samples.min() >= 0.1
        assert samples.max() < 10.0 + 1e-6
        # Median should be close to the log-geometric midpoint (~1.0)
        median = float(np.median(samples))
        assert 0.5 <= median <= 2.0


class TestEdgeCases:
    def test_zero_prob_bucket_yields_no_samples(self, rng):
        probs = {"low": 0.0, "med": 1.0, "high": 0.0}
        samples = _sample_classification_clint(probs, n=1000, rng=rng)
        assert np.all(samples >= 10.0)
        assert np.all(samples < 50.0)

    def test_missing_bucket_key_raises(self, rng):
        with pytest.raises(KeyError):
            _sample_classification_clint({"low": 0.5, "med": 0.5}, n=10, rng=rng)


class TestBuildParamSpecsClassification:
    def test_classification_source_yields_classification_tuple(self):
        from charon.core.schema import (
            BindingProperties,
            CompoundConfig,
            CompoundProperties,
            MetabolismProperties,
            PermeabilityProperties,
            PhysicochemicalProperties,
            PredictedProperty,
        )
        from charon.uncertainty.sampling import build_param_specs

        compound = CompoundConfig(
            name="test",
            smiles="CCO",
            molecular_weight=46.07,
            source="predicted",
            properties=CompoundProperties(
                physicochemical=PhysicochemicalProperties(),
                permeability=PermeabilityProperties(),
                binding=BindingProperties(),
                metabolism=MetabolismProperties(
                    clint_uL_min_mg=PredictedProperty(
                        value=3.0,
                        ci_90_lower=0.1,
                        ci_90_upper=10.0,
                        source="classification",
                        unit="uL/min/10^6 cells",
                        classifier_probs={"low": 0.7, "med": 0.2, "high": 0.1},
                    ),
                ),
            ),
        )
        specs = build_param_specs(compound)
        mu, sigma, dist = specs["clint_uL_min_mg"]
        assert dist == "classification"
        assert mu == {"low": 0.7, "med": 0.2, "high": 0.1}
        assert sigma is None
