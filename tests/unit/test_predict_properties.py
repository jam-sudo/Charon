"""Integration tests for :func:`charon.predict.predict_properties`."""

from __future__ import annotations

from pathlib import Path

import pytest

from charon.core.schema import CompoundProperties
from charon.predict import (
    ADMETPredictor,
    ConformalPredictor,
    predict_properties,
)

ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"
PROPRANOLOL = "CC(C)NCC(O)COc1cccc2ccccc12"
CAFFEINE = "Cn1cnc2c1c(=O)n(C)c(=O)n2C"

CAL_PATH = Path("/home/jam/Charon/data/validation/adme_reference.csv")


@pytest.fixture(scope="module")
def predictor() -> ADMETPredictor:
    return ADMETPredictor()


@pytest.fixture(scope="module")
def conformal(predictor: ADMETPredictor) -> ConformalPredictor:
    cp = ConformalPredictor(CAL_PATH)
    cp.calibrate(predictor)
    return cp


class TestPredictPropertiesBasics:
    def test_returns_compound_properties(self, predictor: ADMETPredictor):
        props = predict_properties(ASPIRIN, predictor=predictor)
        assert isinstance(props, CompoundProperties)

    def test_all_expected_fields_populated(
        self, predictor: ADMETPredictor
    ):
        props = predict_properties(ASPIRIN, predictor=predictor)
        assert props.physicochemical.logp is not None
        assert props.physicochemical.pka_acid is not None  # aspirin is an acid
        assert props.binding.fu_p is not None
        assert props.binding.fu_inc is not None
        assert props.binding.bp_ratio is not None
        assert props.metabolism.clint_uL_min_mg is not None
        assert props.renal.clrenal_L_h is not None

    def test_without_conformal_no_ci(self, predictor: ADMETPredictor):
        """No conformal predictor → fu_p has no CI bounds attached."""
        props = predict_properties(ASPIRIN, predictor=predictor)
        fu_p = props.binding.fu_p
        assert fu_p is not None
        assert fu_p.ci_90_lower is None
        assert fu_p.ci_90_upper is None

    def test_with_conformal_attaches_ci(
        self,
        predictor: ADMETPredictor,
        conformal: ConformalPredictor,
    ):
        """Calibrated conformal → fu_p gets CI bounds."""
        props = predict_properties(
            ASPIRIN, predictor=predictor, conformal=conformal
        )
        fu_p = props.binding.fu_p
        assert fu_p is not None
        assert fu_p.ci_90_lower is not None
        assert fu_p.ci_90_upper is not None
        assert fu_p.ci_90_lower <= fu_p.value <= fu_p.ci_90_upper


class TestSourceTags:
    """Provenance tags must match the documented surface."""

    def test_source_tags(self, predictor: ADMETPredictor):
        props = predict_properties(ASPIRIN, predictor=predictor)
        assert props.binding.fu_p.source == "ml_ensemble"
        assert props.binding.fu_inc.source == "correlation"
        assert props.binding.bp_ratio.source == "derived"
        assert props.physicochemical.pka_acid.source == "ml_pka"
        assert props.renal.clrenal_L_h.source == "derived"

    def test_clint_flag_contains_tier2(self, predictor: ADMETPredictor):
        props = predict_properties(ASPIRIN, predictor=predictor)
        clint = props.metabolism.clint_uL_min_mg
        assert clint is not None
        assert clint.flag is not None
        assert "tier2_ml" in clint.flag


class TestCompoundTypeClassification:
    def test_aspirin_is_acid(self, predictor: ADMETPredictor):
        props = predict_properties(ASPIRIN, predictor=predictor)
        bp = props.binding.bp_ratio
        assert bp is not None
        assert "compound_type=acid" in (bp.method or "")

    def test_propranolol_is_base(self, predictor: ADMETPredictor):
        props = predict_properties(PROPRANOLOL, predictor=predictor)
        bp = props.binding.bp_ratio
        assert bp is not None
        assert "compound_type=base" in (bp.method or "")

    def test_caffeine_is_neutral(self, predictor: ADMETPredictor):
        props = predict_properties(CAFFEINE, predictor=predictor)
        bp = props.binding.bp_ratio
        assert bp is not None
        assert "compound_type=neutral" in (bp.method or "")


class TestErrorHandling:
    def test_invalid_smiles_raises(self, predictor: ADMETPredictor):
        with pytest.raises(ValueError):
            predict_properties("not_a_valid_smiles_xyz_777", predictor=predictor)

    def test_empty_smiles_raises(self, predictor: ADMETPredictor):
        with pytest.raises(ValueError):
            predict_properties("", predictor=predictor)


class TestDeterminism:
    def test_same_smiles_same_properties(
        self, predictor: ADMETPredictor
    ):
        """Two identical runs on the same SMILES yield identical scalars."""
        p1 = predict_properties(ASPIRIN, predictor=predictor)
        p2 = predict_properties(ASPIRIN, predictor=predictor)
        assert p1.binding.fu_p.value == p2.binding.fu_p.value
        assert p1.binding.fu_inc.value == p2.binding.fu_inc.value
        assert p1.binding.bp_ratio.value == p2.binding.bp_ratio.value
        assert (
            p1.metabolism.clint_uL_min_mg.value
            == p2.metabolism.clint_uL_min_mg.value
        )
        assert p1.renal.clrenal_L_h.value == p2.renal.clrenal_L_h.value
        assert p1.physicochemical.logp.value == p2.physicochemical.logp.value
