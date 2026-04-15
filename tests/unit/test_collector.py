"""Unit tests for charon.report.collector."""

from __future__ import annotations

from dataclasses import is_dataclass

import numpy as np
import pytest

from charon.core.schema import (
    BindingProperties,
    CompoundConfig,
    CompoundProperties,
    MetabolismProperties,
    PhysicochemicalProperties,
    PKParameters,
    PredictedProperty,
)
from charon.pbpk.solver import SimulationResult
from charon.pipeline import PipelineResult
from charon.report.collector import ReportData, collect


def _make_minimal_result(
    *,
    name: str = "test-compound",
    smiles: str = "CCO",
    route: str = "iv_bolus",
    dose_mg: float = 10.0,
) -> PipelineResult:
    """Build a PipelineResult with just enough data for collector tests."""
    compound = CompoundConfig(
        name=name,
        smiles=smiles,
        molecular_weight=46.07,
        source="predicted",
        properties=CompoundProperties(),
    )
    time_h = np.array([0.0, 1.0, 2.0, 4.0, 8.0, 24.0], dtype=float)
    cp_plasma = np.array([10.0, 8.0, 6.0, 4.0, 2.0, 0.5], dtype=float)
    cp_blood = cp_plasma * 0.9
    # State trajectory shape = (n_states, n_time); only a sentinel is needed here.
    state_trajectory = np.zeros((3, time_h.size), dtype=float)
    sim = SimulationResult(
        time_h=time_h,
        cp_blood=cp_blood,
        cp_plasma=cp_plasma,
        state_trajectory=state_trajectory,
        mass_balance_residual=0.0,
        solver_success=True,
        solver_method="BDF",
        solver_nfev=42,
        route=route,
        dose_mg=dose_mg,
        infusion_duration_h=0.0,
    )
    pk = PKParameters(
        cmax=10.0,
        tmax=0.0,
        auc_0_inf=50.0,
        half_life=4.0,
        cl_apparent=0.2,
        vss=1.5,
    )
    return PipelineResult(
        compound=compound,
        pk_parameters=pk,
        time_h=time_h,
        cp_plasma=cp_plasma,
        cp_blood=cp_blood,
        simulation=sim,
        metadata={
            "species": "human",
            "route": route,
            "dose_mg": dose_mg,
            "duration_h": 24.0,
            "liver_model": "well_stirred",
        },
    )


def test_report_data_is_frozen_dataclass():
    assert is_dataclass(ReportData)
    assert ReportData.__dataclass_params__.frozen is True  # type: ignore[attr-defined]


def test_collect_identity_fields():
    result = _make_minimal_result(name="ethanol", smiles="CCO", dose_mg=10.0)
    data = collect(result)
    assert data.compound_name == "ethanol"
    assert data.smiles == "CCO"
    assert data.molecular_weight == pytest.approx(46.07)
    assert data.source == "predicted"
    assert data.route == "iv_bolus"
    assert data.dose_mg == pytest.approx(10.0)
    assert data.duration_h == pytest.approx(24.0)


def test_collect_timestamp_is_iso8601():
    result = _make_minimal_result()
    data = collect(result)
    assert "T" in data.timestamp
    assert len(data.timestamp) >= 19


def test_collect_respects_explicit_timestamp():
    result = _make_minimal_result()
    data = collect(result, timestamp="2026-04-15T12:00:00+00:00")
    assert data.timestamp == "2026-04-15T12:00:00+00:00"


def test_collect_warnings_default_empty():
    result = _make_minimal_result()
    data = collect(result)
    assert data.warnings == []


def test_collect_warnings_passthrough():
    result = _make_minimal_result()
    data = collect(result, warnings=["high uncertainty"])
    assert data.warnings == ["high uncertainty"]


def _make_result_with_properties() -> PipelineResult:
    compound = CompoundConfig(
        name="midazolam",
        smiles="Cc1ncc(n1C)[C@@H](c2ccccc2)OC(=O)N",
        molecular_weight=325.77,
        source="experimental",
        properties=CompoundProperties(
            physicochemical=PhysicochemicalProperties(
                logp=PredictedProperty(
                    value=3.89,
                    ci_90_lower=3.4,
                    ci_90_upper=4.3,
                    source="ml_ensemble",
                    unit="log",
                ),
                compound_type="base",
            ),
            binding=BindingProperties(
                fu_p=PredictedProperty(
                    value=0.032,
                    ci_90_lower=0.02,
                    ci_90_upper=0.05,
                    source="ml_ensemble",
                    unit="fraction",
                ),
            ),
            metabolism=MetabolismProperties(
                clint_uL_min_mg=PredictedProperty(
                    value=93.0,
                    source="experimental",
                    unit="uL/min/mg",
                    flag="clint_tier2_ml",
                ),
            ),
        ),
    )
    time_h = np.array([0.0, 1.0, 2.0], dtype=float)
    cp_plasma = np.array([1.0, 0.5, 0.25], dtype=float)
    cp_blood = cp_plasma * 0.9
    state_trajectory = np.zeros((3, time_h.size), dtype=float)
    sim = SimulationResult(
        time_h=time_h,
        cp_blood=cp_blood,
        cp_plasma=cp_plasma,
        state_trajectory=state_trajectory,
        mass_balance_residual=0.0,
        solver_success=True,
        solver_method="BDF",
        solver_nfev=42,
        route="oral",
        dose_mg=5.0,
        infusion_duration_h=0.0,
    )
    pk = PKParameters(
        cmax=1.0,
        tmax=0.0,
        auc_0_inf=2.0,
        half_life=1.5,
        cl_apparent=0.3,
        vss=1.0,
    )
    return PipelineResult(
        compound=compound,
        pk_parameters=pk,
        time_h=time_h,
        cp_plasma=cp_plasma,
        cp_blood=cp_blood,
        simulation=sim,
        metadata={
            "route": "oral",
            "dose_mg": 5.0,
            "duration_h": 24.0,
            "liver_model": "well_stirred",
            "clint_liver_L_h": 348.75,
            "cl_renal_L_h": 0.23,
            "fu_b": 0.071,
            "compound_type": "base",
        },
    )


def test_collect_properties_flattened():
    result = _make_result_with_properties()
    data = collect(result)
    # logp
    assert "logp" in data.properties
    assert data.properties["logp"]["value"] == pytest.approx(3.89)
    assert data.properties["logp"]["ci_lower"] == pytest.approx(3.4)
    assert data.properties["logp"]["ci_upper"] == pytest.approx(4.3)
    assert data.properties["logp"]["source"] == "ml_ensemble"
    # fu_p
    assert "fu_p" in data.properties
    assert data.properties["fu_p"]["value"] == pytest.approx(0.032)
    # clint
    assert "clint_uL_min_mg" in data.properties
    assert data.properties["clint_uL_min_mg"]["flag"] == "clint_tier2_ml"


def test_collect_compound_type():
    result = _make_result_with_properties()
    data = collect(result)
    assert data.compound_type == "base"


def test_collect_ivive_summary_from_metadata():
    result = _make_result_with_properties()
    data = collect(result)
    assert data.ivive_summary["clint_liver_L_h"] == pytest.approx(348.75)
    assert data.ivive_summary["cl_renal_L_h"] == pytest.approx(0.23)
    assert data.ivive_summary["fu_b"] == pytest.approx(0.071)
    assert data.ivive_summary["liver_model"] == "well_stirred"


def test_collect_skips_none_properties():
    result = _make_minimal_result()
    data = collect(result)
    # No properties populated in the minimal fixture
    assert data.properties == {}


from charon.translational.dose_projector import FIHDoseRecommendation
from charon.translational.hed import HEDResult
from charon.uncertainty.dose_range import UncertaintyResult


def test_collect_pk_params_flattened():
    result = _make_minimal_result()
    data = collect(result)
    assert data.pk_params["cmax"] == pytest.approx(10.0)
    assert data.pk_params["auc_0_inf"] == pytest.approx(50.0)
    assert data.pk_params["cl_apparent"] == pytest.approx(0.2)
    assert data.pk_params["vss"] == pytest.approx(1.5)
    assert "bioavailability" in data.pk_params


def test_collect_pk_table_uses_canonical_timepoints():
    result = _make_minimal_result()
    data = collect(result)
    # Canonical points <= max(time_h=24.0): 0, 0.25, 0.5, 1, 2, 4, 6, 8, 12, 24
    times = [row["time_h"] for row in data.pk_table]
    assert 0.0 in times
    assert 1.0 in times
    assert 24.0 in times
    # 48 and 72 are beyond duration
    assert 48.0 not in times
    assert 72.0 not in times
    # Monotonic and <= 12 entries
    assert times == sorted(times)
    assert len(data.pk_table) <= 12


def test_collect_pk_table_cp_values():
    result = _make_minimal_result()
    data = collect(result)
    first = data.pk_table[0]
    assert "cp_plasma_ug_L" in first
    assert "cp_blood_ug_L" in first
    assert isinstance(first["cp_plasma_ug_L"], float)


def test_collect_dose_recommendation_flattened():
    result = _make_minimal_result()
    hed = HEDResult(
        noael_mg_kg=50.0,
        noael_species="rat",
        km_animal=6.2,
        km_human=37.0,
        hed_mg_kg=8.06,
        body_weight_kg=70.0,
        safety_factor=10.0,
        mrsd_mg=56.45,
    )
    rec = FIHDoseRecommendation(
        mrsd_mg=56.45,
        limiting_method="hed",
        hed=hed,
        mabel=None,
        pad=None,
        safety_factor=10.0,
        salt_factor=1.0,
        route="oral",
        rationale="HED: 56.45 mg",
    )
    result.dose_recommendation = rec
    data = collect(result)
    assert data.dose_recommendation is not None
    assert data.dose_recommendation["mrsd_mg"] == pytest.approx(56.45)
    assert data.dose_recommendation["limiting_method"] == "hed"
    assert data.dose_recommendation["hed"]["mrsd_mg"] == pytest.approx(56.45)
    assert data.dose_recommendation["mabel"] is None


def test_collect_uncertainty_flattened():
    result = _make_minimal_result()
    unc = UncertaintyResult(
        point_estimate_mg=5.0,
        ci_90_lower_mg=2.0,
        ci_90_upper_mg=12.0,
        ci_ratio=6.0,
        confidence="MEDIUM",
        n_samples=100,
        n_successful=95,
        convergence_met=True,
        sensitivity={"clint": 0.7, "fu_p": 0.2, "logp": 0.1},
        limiting_parameter="clint",
        recommendation="Experimental clint measurement would narrow CI by ~70%",
        r_squared=0.85,
    )
    result.uncertainty = unc
    data = collect(result)
    assert data.uncertainty is not None
    assert data.uncertainty["point_estimate_mg"] == pytest.approx(5.0)
    assert data.uncertainty["confidence"] == "MEDIUM"
    assert data.uncertainty["sensitivity"]["clint"] == pytest.approx(0.7)
