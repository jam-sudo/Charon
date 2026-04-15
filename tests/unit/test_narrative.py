"""Unit tests for charon.report.narrative."""

from __future__ import annotations

import pytest

from charon.report.collector import ReportData
from charon.report.narrative import (
    _render_compound_profile,
    _render_executive_summary,
    _render_header,
    format_value,
)


def _make_data(**overrides) -> ReportData:
    base = dict(
        compound_name="midazolam",
        smiles="Cc1ncc(n1C)[C@@H](c2ccccc2)OC(=O)N",
        molecular_weight=325.77,
        source="experimental",
        compound_type="base",
        properties={},
        ivive_summary={},
        pk_params={},
        pk_table=[],
        route="oral",
        dose_mg=5.0,
        duration_h=24.0,
        dose_recommendation=None,
        uncertainty=None,
        warnings=[],
        metadata={"species": "human"},
        timestamp="2026-04-15T00:00:00+00:00",
        charon_version="0.1.0",
    )
    base.update(overrides)
    return ReportData(**base)


def test_format_value_none():
    assert format_value(None) == "-"


def test_format_value_small():
    assert format_value(3.89) == "3.89"


def test_format_value_large_scientific():
    s = format_value(1234567.0)
    assert "e" in s or "E" in s or "1.23" in s


def test_format_value_tiny_scientific():
    s = format_value(0.000123)
    assert "e" in s or "E" in s


def test_render_header_contains_name():
    data = _make_data()
    out = _render_header(data)
    assert "midazolam" in out
    assert out.startswith("# ")


def test_render_header_contains_timestamp():
    data = _make_data()
    out = _render_header(data)
    assert "2026-04-15" in out


def test_render_executive_summary_without_dose():
    data = _make_data()
    out = _render_executive_summary(data)
    assert "## 1. Executive Summary" in out
    # No dose rec -> pipeline-only note
    assert "No FIH dose projection" in out or "not run" in out.lower()


def test_render_executive_summary_with_dose_no_uncertainty():
    data = _make_data(
        dose_recommendation={
            "mrsd_mg": 56.45,
            "limiting_method": "hed",
            "route": "oral",
            "salt_factor": 1.0,
            "safety_factor": 10.0,
            "rationale": "HED",
            "hed": {"mrsd_mg": 56.45},
            "mabel": None,
            "pad": None,
        }
    )
    out = _render_executive_summary(data)
    assert "56.45" in out
    assert "HED" in out or "hed" in out.lower()


def test_render_executive_summary_with_uncertainty():
    data = _make_data(
        dose_recommendation={
            "mrsd_mg": 56.45,
            "limiting_method": "hed",
            "route": "oral",
            "salt_factor": 1.0,
            "safety_factor": 10.0,
            "rationale": "HED",
            "hed": {"mrsd_mg": 56.45},
            "mabel": None,
            "pad": None,
        },
        uncertainty={
            "point_estimate_mg": 50.0,
            "ci_90_lower_mg": 20.0,
            "ci_90_upper_mg": 120.0,
            "confidence": "LOW",
            "ci_ratio": 6.0,
            "n_samples": 100,
            "n_successful": 95,
            "sensitivity": {},
            "limiting_parameter": "clint",
            "recommendation": "",
            "convergence_met": True,
            "r_squared": 0.9,
        },
    )
    out = _render_executive_summary(data)
    assert "20" in out and "120" in out  # CI bounds
    assert "LOW" in out


def test_render_compound_profile_lists_identity():
    data = _make_data()
    out = _render_compound_profile(data)
    assert "## 2. Compound Profile" in out
    assert "Cc1ncc" in out  # SMILES
    assert "325.77" in out
    assert "base" in out
    assert "experimental" in out
