"""Unit tests for charon.report.collector."""

from __future__ import annotations

from dataclasses import is_dataclass

import numpy as np
import pytest

from charon.core.schema import (
    CompoundConfig,
    CompoundProperties,
    PKParameters,
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
