"""Unit tests for charon.report.export."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from charon.report.collector import ReportData
from charon.report.export import (
    export_json,
    export_markdown,
    export_report,
)


def _make_data() -> ReportData:
    return ReportData(
        compound_name="test",
        smiles="CCO",
        molecular_weight=46.07,
        source="predicted",
        compound_type="neutral",
        properties={
            "logp": {
                "value": -0.31,
                "ci_lower": None,
                "ci_upper": None,
                "unit": "log",
                "source": "ml_ensemble",
                "flag": None,
                "method": None,
            }
        },
        ivive_summary={"liver_model": "well_stirred", "fu_b": 0.95},
        pk_params={"cmax": 10.0, "auc_0_inf": 50.0},
        pk_table=[{"time_h": 0.0, "cp_plasma_ug_L": 10.0, "cp_blood_ug_L": 9.0}],
        route="iv_bolus",
        dose_mg=10.0,
        duration_h=24.0,
        dose_recommendation=None,
        uncertainty=None,
        warnings=[],
        metadata={"species": "human", "solver_nfev": 42},
        timestamp="2026-04-15T00:00:00+00:00",
        charon_version="0.1.0",
    )


def test_export_markdown_writes_file(tmp_path: Path):
    data = _make_data()
    out = tmp_path / "report.md"
    path = export_markdown(data, out)
    assert path == out.resolve()
    text = out.read_text()
    assert "# FIH Dose Rationale Report" in text
    assert "test" in text
    # At least 8 sections visible (1-6, 8, 9; section 7 skipped)
    assert text.count("## ") >= 8


def test_export_json_writes_valid_json(tmp_path: Path):
    data = _make_data()
    out = tmp_path / "report.json"
    path = export_json(data, out)
    assert path == out.resolve()
    parsed = json.loads(out.read_text())
    assert parsed["compound_name"] == "test"
    assert parsed["smiles"] == "CCO"
    assert parsed["pk_params"]["cmax"] == 10.0
    assert parsed["ivive_summary"]["liver_model"] == "well_stirred"


def test_export_json_has_all_top_level_keys(tmp_path: Path):
    data = _make_data()
    out = tmp_path / "r.json"
    export_json(data, out)
    parsed = json.loads(out.read_text())
    expected_keys = {
        "compound_name", "smiles", "molecular_weight", "source",
        "compound_type", "properties", "ivive_summary", "pk_params",
        "pk_table", "route", "dose_mg", "duration_h",
        "dose_recommendation", "uncertainty", "warnings", "metadata",
        "timestamp", "charon_version",
    }
    assert expected_keys.issubset(parsed.keys())


def test_export_json_full_profile_optional(tmp_path: Path):
    data = _make_data()
    out = tmp_path / "r.json"
    export_json(
        data,
        out,
        include_full_profile=True,
        full_profile={
            "time_h": [0.0, 1.0, 2.0],
            "cp_plasma_ug_L": [10.0, 8.0, 6.0],
            "cp_blood_ug_L": [9.0, 7.2, 5.4],
        },
    )
    parsed = json.loads(out.read_text())
    assert "full_profile" in parsed
    assert parsed["full_profile"]["time_h"] == [0.0, 1.0, 2.0]


def test_export_report_writes_both_md_and_json(tmp_path: Path):
    data = _make_data()
    md_path, json_path = export_report(data, tmp_path / "run.md")
    assert md_path.exists()
    assert json_path.exists()
    assert md_path.suffix == ".md"
    assert json_path.suffix == ".json"
    assert md_path.stem == json_path.stem == "run"


def test_export_report_handles_no_suffix(tmp_path: Path):
    data = _make_data()
    md_path, json_path = export_report(data, tmp_path / "run")
    assert md_path.exists()
    assert json_path.exists()
    assert md_path.name == "run.md"
    assert json_path.name == "run.json"


def test_export_json_handles_numpy_values(tmp_path: Path):
    data = _make_data()
    # Simulate a numpy value sneaking through (e.g. np.float64 inside metadata)
    data_with_np = ReportData(
        **{
            **{k: getattr(data, k) for k in data.__dataclass_fields__},
            "metadata": {"solver_nfev": np.int64(42), "cl": np.float64(1.23)},
        }
    )
    out = tmp_path / "np.json"
    export_json(data_with_np, out)
    parsed = json.loads(out.read_text())
    assert parsed["metadata"]["solver_nfev"] == 42
    assert parsed["metadata"]["cl"] == pytest.approx(1.23)


def test_export_json_non_finite_floats_become_null(tmp_path):
    data = _make_data()
    data_with_inf = ReportData(
        **{
            **{k: getattr(data, k) for k in data.__dataclass_fields__},
            "metadata": {"bad": float("inf"), "also_bad": float("nan"), "ok": 1.5},
        }
    )
    out = tmp_path / "inf.json"
    export_json(data_with_inf, out)
    text = out.read_text()
    # Strict JSON parse (raises on bare NaN/Infinity tokens)
    parsed = json.loads(text)
    assert parsed["metadata"]["bad"] is None
    assert parsed["metadata"]["also_bad"] is None
    assert parsed["metadata"]["ok"] == 1.5
    # Should not contain bare NaN / Infinity tokens
    assert "NaN" not in text
    assert "Infinity" not in text


def test_export_json_nested_numpy_array(tmp_path):
    data = _make_data()
    data_nested = ReportData(
        **{
            **{k: getattr(data, k) for k in data.__dataclass_fields__},
            "metadata": {"arr": np.array([1.0, 2.0, 3.0])},
        }
    )
    out = tmp_path / "nested.json"
    export_json(data_nested, out)
    parsed = json.loads(out.read_text())
    assert parsed["metadata"]["arr"] == [1.0, 2.0, 3.0]


def test_export_json_numpy_bool(tmp_path):
    data = _make_data()
    data_bool = ReportData(
        **{
            **{k: getattr(data, k) for k in data.__dataclass_fields__},
            "metadata": {"converged": np.bool_(True)},
        }
    )
    out = tmp_path / "bool.json"
    export_json(data_bool, out)
    parsed = json.loads(out.read_text())
    assert parsed["metadata"]["converged"] is True
