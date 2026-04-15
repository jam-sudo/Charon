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
