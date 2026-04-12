"""Unit tests for charon.core.schema Pydantic models."""

from __future__ import annotations

import math

import pytest
from pydantic import ValidationError

from charon.core.schema import (
    CompoundConfig,
    CompoundProperties,
    ConversionLog,
    ConversionStep,
    GuardrailWarning,
    HepaticClearance,
    PipelineConfig,
    PredictedProperty,
    RunConfig,
    SaltForm,
    ValidationResult,
)


# ---------------------------------------------------------------------------
# PredictedProperty
# ---------------------------------------------------------------------------


class TestPredictedProperty:
    """Tests for the PredictedProperty value wrapper."""

    def test_valid_creation_minimal(self):
        pp = PredictedProperty(value=3.5, source="ml_ensemble")
        assert pp.value == 3.5
        assert pp.source == "ml_ensemble"
        assert pp.ci_90_lower is None
        assert pp.ci_90_upper is None
        assert pp.unit is None
        assert pp.method is None
        assert pp.flag is None

    def test_valid_creation_full(self):
        pp = PredictedProperty(
            value=3.5,
            ci_90_lower=2.0,
            ci_90_upper=5.0,
            source="experimental",
            unit="nM",
            method="LC-MS/MS",
            flag="high_confidence",
        )
        assert pp.value == 3.5
        assert pp.ci_90_lower == 2.0
        assert pp.ci_90_upper == 5.0
        assert pp.unit == "nM"

    def test_ci_lower_exceeds_value_raises(self):
        with pytest.raises(ValidationError, match="ci_90_lower"):
            PredictedProperty(
                value=3.0,
                ci_90_lower=4.0,
                ci_90_upper=5.0,
                source="ml_ensemble",
            )

    def test_value_exceeds_ci_upper_raises(self):
        with pytest.raises(ValidationError, match="ci_90_upper"):
            PredictedProperty(
                value=6.0,
                ci_90_lower=2.0,
                ci_90_upper=5.0,
                source="ml_ensemble",
            )

    def test_ci_bounds_at_value_allowed(self):
        """lower == value == upper is valid (zero-width CI)."""
        pp = PredictedProperty(
            value=3.0,
            ci_90_lower=3.0,
            ci_90_upper=3.0,
            source="ml_ensemble",
        )
        assert pp.ci_90_lower == pp.value == pp.ci_90_upper

    def test_nan_value_raises(self):
        with pytest.raises(ValidationError, match="finite"):
            PredictedProperty(value=float("nan"), source="ml_ensemble")

    def test_inf_value_raises(self):
        with pytest.raises(ValidationError, match="finite"):
            PredictedProperty(value=float("inf"), source="ml_ensemble")

    def test_neg_inf_value_raises(self):
        with pytest.raises(ValidationError, match="finite"):
            PredictedProperty(value=float("-inf"), source="ml_ensemble")

    def test_negative_value_allowed(self):
        """Negative values are valid (e.g. logP can be negative)."""
        pp = PredictedProperty(value=-2.5, source="ml_ensemble")
        assert pp.value == -2.5

    def test_invalid_source_raises(self):
        with pytest.raises(ValidationError):
            PredictedProperty(value=1.0, source="invalid_source_type")

    def test_all_valid_source_types(self):
        for src in [
            "ml_ensemble",
            "ml_pka",
            "correlation",
            "derived",
            "physiological",
            "experimental",
        ]:
            pp = PredictedProperty(value=1.0, source=src)
            assert pp.source == src


# ---------------------------------------------------------------------------
# SaltForm
# ---------------------------------------------------------------------------


class TestSaltForm:
    """Tests for SaltForm model."""

    def test_default_salt_factor(self):
        sf = SaltForm()
        assert sf.salt_factor == 1.0
        assert sf.name is None
        assert sf.mw_salt is None

    def test_creation_with_name(self):
        sf = SaltForm(name="hydrochloride")
        assert sf.name == "hydrochloride"
        assert sf.salt_factor == 1.0

    def test_creation_with_mw_salt(self):
        sf = SaltForm(name="sodium salt", mw_salt=202.14)
        assert sf.mw_salt == 202.14
        # salt_factor remains 1.0 at SaltForm level; auto-calc deferred
        # to CompoundConfig model_validator.
        assert sf.salt_factor == 1.0

    def test_explicit_salt_factor(self):
        sf = SaltForm(name="HCl", salt_factor=0.85)
        assert sf.salt_factor == 0.85


# ---------------------------------------------------------------------------
# ConversionLog / ConversionStep
# ---------------------------------------------------------------------------


class TestConversionLog:
    """Tests for ConversionLog and ConversionStep models."""

    def test_creation_with_steps(self):
        steps = [
            ConversionStep(
                name="scale", value=60000.0, unit="mg protein", formula="40 * 1500"
            ),
            ConversionStep(
                name="CLint_liver", value=3600.0, unit="uL/min", formula="0.06 * 60000"
            ),
        ]
        log = ConversionLog(
            input_params={"clint": 10.0, "system": "HLM"},
            intermediate_steps=steps,
            output=25.0,
            output_unit="L/h",
            model_used="well_stirred",
        )
        assert len(log.intermediate_steps) == 2
        assert log.output == 25.0
        assert log.output_unit == "L/h"
        assert log.model_used == "well_stirred"
        assert log.input_params["clint"] == 10.0

    def test_empty_steps(self):
        log = ConversionLog(
            input_params={},
            intermediate_steps=[],
            output=0.0,
            output_unit="L/h",
            model_used="well_stirred",
        )
        assert log.intermediate_steps == []


# ---------------------------------------------------------------------------
# HepaticClearance
# ---------------------------------------------------------------------------


class TestHepaticClearance:
    """Tests for HepaticClearance model."""

    def test_creation(self):
        log = ConversionLog(
            input_params={"clint": 10.0},
            intermediate_steps=[],
            output=25.0,
            output_unit="L/h",
            model_used="well_stirred",
        )
        hc = HepaticClearance(
            clh_L_h=25.0,
            extraction_ratio=0.278,
            model_used="well_stirred",
            conversion_log=log,
        )
        assert hc.clh_L_h == 25.0
        assert hc.extraction_ratio == pytest.approx(0.278)
        assert hc.model_used == "well_stirred"
        assert hc.conversion_log.output == 25.0


# ---------------------------------------------------------------------------
# CompoundConfig
# ---------------------------------------------------------------------------


class TestCompoundConfig:
    """Tests for CompoundConfig model."""

    def test_minimal_creation(self):
        cc = CompoundConfig(
            name="aspirin",
            smiles="CC(=O)Oc1ccccc1C(=O)O",
        )
        assert cc.name == "aspirin"
        assert cc.smiles == "CC(=O)Oc1ccccc1C(=O)O"
        assert cc.molecular_weight is None
        assert cc.salt_form is None
        assert cc.source == "predicted"

    def test_full_creation(self):
        cc = CompoundConfig(
            name="aspirin",
            smiles="CC(=O)Oc1ccccc1C(=O)O",
            molecular_weight=180.16,
            source="experimental",
            salt_form=SaltForm(name="free acid"),
        )
        assert cc.molecular_weight == 180.16
        assert cc.source == "experimental"
        assert cc.salt_form.name == "free acid"

    def test_auto_salt_factor_calculation(self):
        """When both mw_salt and molecular_weight provided, salt_factor is auto-computed."""
        cc = CompoundConfig(
            name="test",
            smiles="CC",
            molecular_weight=180.0,
            salt_form=SaltForm(name="sodium salt", mw_salt=200.0),
        )
        assert cc.salt_form.salt_factor == pytest.approx(180.0 / 200.0)

    def test_explicit_salt_factor_not_overwritten(self):
        """When salt_factor is explicitly set to non-1.0, it should NOT be overwritten."""
        cc = CompoundConfig(
            name="test",
            smiles="CC",
            molecular_weight=180.0,
            salt_form=SaltForm(name="HCl", mw_salt=200.0, salt_factor=0.75),
        )
        # salt_factor was explicitly 0.75 (!= 1.0), so auto-calc should skip
        assert cc.salt_form.salt_factor == pytest.approx(0.75)

    def test_no_auto_salt_without_mw(self):
        """Without molecular_weight, salt_factor remains default."""
        cc = CompoundConfig(
            name="test",
            smiles="CC",
            salt_form=SaltForm(name="HCl", mw_salt=200.0),
        )
        assert cc.salt_form.salt_factor == 1.0

    def test_default_properties(self):
        cc = CompoundConfig(name="x", smiles="C")
        assert cc.properties.physicochemical.logp is None
        assert cc.properties.binding.fu_p is None


# ---------------------------------------------------------------------------
# RunConfig (frozen=True)
# ---------------------------------------------------------------------------


class TestRunConfig:
    """Tests for RunConfig immutability."""

    def test_creation(self):
        rc = RunConfig(
            compound=CompoundConfig(name="test", smiles="CC"),
        )
        assert rc.compound.name == "test"
        assert rc.pipeline.liver_model == "well_stirred"
        assert rc.model_versions == {}

    def test_frozen_prevents_mutation(self):
        rc = RunConfig(
            compound=CompoundConfig(name="test", smiles="CC"),
        )
        with pytest.raises(ValidationError):
            rc.pipeline = PipelineConfig()

    def test_frozen_prevents_field_assignment(self):
        rc = RunConfig(
            compound=CompoundConfig(name="test", smiles="CC"),
        )
        with pytest.raises(ValidationError):
            rc.model_versions = {"v1": "1.0"}


# ---------------------------------------------------------------------------
# ValidationResult
# ---------------------------------------------------------------------------


class TestValidationResult:
    """Tests for ValidationResult model."""

    def test_minimal_creation(self):
        vr = ValidationResult(
            is_valid=True,
            smiles_canonical="CC",
            molecular_weight=30.07,
        )
        assert vr.is_valid is True
        assert vr.warnings == []
        assert vr.applicability_domain == "UNKNOWN"
        assert vr.descriptors == {}

    def test_creation_with_warnings(self):
        warnings = [
            GuardrailWarning(
                category="MW",
                message="MW > 500 Da",
                severity="WARNING",
            ),
            GuardrailWarning(
                category="logP",
                message="logP > 5",
                severity="CRITICAL",
            ),
        ]
        vr = ValidationResult(
            is_valid=False,
            smiles_canonical="CC(CC)CC(CC)CC(CC)CC",
            molecular_weight=600.0,
            warnings=warnings,
            applicability_domain="LOW",
            descriptors={"MW": 600.0, "logP": 5.5},
        )
        assert vr.is_valid is False
        assert len(vr.warnings) == 2
        assert vr.warnings[0].severity == "WARNING"
        assert vr.warnings[1].severity == "CRITICAL"
        assert vr.applicability_domain == "LOW"
        assert vr.descriptors["logP"] == 5.5

    def test_invalid_applicability_domain_raises(self):
        with pytest.raises(ValidationError):
            ValidationResult(
                is_valid=True,
                smiles_canonical="CC",
                molecular_weight=30.0,
                applicability_domain="INVALID",
            )

    def test_invalid_severity_raises(self):
        with pytest.raises(ValidationError):
            GuardrailWarning(
                category="test",
                message="test",
                severity="INVALID",
            )


# ---------------------------------------------------------------------------
# MetabolismProperties.fm_cyp3a4
# ---------------------------------------------------------------------------


class TestFmCyp3a4:
    def test_fm_cyp3a4_default_none(self):
        from charon.core.schema import MetabolismProperties
        m = MetabolismProperties()
        assert m.fm_cyp3a4 is None

    def test_fm_cyp3a4_valid(self):
        from charon.core.schema import MetabolismProperties
        m = MetabolismProperties(fm_cyp3a4=0.75)
        assert m.fm_cyp3a4 == 0.75

    def test_fm_cyp3a4_zero(self):
        from charon.core.schema import MetabolismProperties
        m = MetabolismProperties(fm_cyp3a4=0.0)
        assert m.fm_cyp3a4 == 0.0

    def test_fm_cyp3a4_one(self):
        from charon.core.schema import MetabolismProperties
        m = MetabolismProperties(fm_cyp3a4=1.0)
        assert m.fm_cyp3a4 == 1.0

    def test_fm_cyp3a4_negative_rejected(self):
        from charon.core.schema import MetabolismProperties
        import pytest
        with pytest.raises(Exception):
            MetabolismProperties(fm_cyp3a4=-0.1)

    def test_fm_cyp3a4_above_one_rejected(self):
        from charon.core.schema import MetabolismProperties
        import pytest
        with pytest.raises(Exception):
            MetabolismProperties(fm_cyp3a4=1.01)
