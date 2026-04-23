"""Tests for PredictedProperty.classifier_probs + SourceType `classification`."""
from __future__ import annotations

import pytest
from pydantic import ValidationError

from charon.core.schema import PredictedProperty


class TestSourceLiteral:
    def test_classification_source_accepted(self):
        prop = PredictedProperty(
            value=3.0, source="classification", unit="uL/min/10^6 cells",
        )
        assert prop.source == "classification"

    def test_invalid_source_rejected(self):
        with pytest.raises(ValidationError):
            PredictedProperty(value=3.0, source="bogus_source")


class TestClassifierProbsField:
    def test_none_by_default(self):
        prop = PredictedProperty(value=3.0, source="ml_ensemble")
        assert prop.classifier_probs is None

    def test_accepts_summing_to_one(self):
        probs = {"low": 0.6, "med": 0.3, "high": 0.1}
        prop = PredictedProperty(
            value=3.0, source="classification", classifier_probs=probs,
        )
        assert prop.classifier_probs == probs

    def test_rejects_not_summing_to_one(self):
        with pytest.raises(ValidationError, match="classifier_probs"):
            PredictedProperty(
                value=3.0, source="classification",
                classifier_probs={"low": 0.5, "med": 0.5, "high": 0.5},
            )

    def test_allows_two_percent_tolerance(self):
        """Tolerance band accommodates float32 -> float64 rounding."""
        prop = PredictedProperty(
            value=3.0, source="classification",
            classifier_probs={"low": 0.334, "med": 0.333, "high": 0.333},
        )
        assert prop.classifier_probs is not None
