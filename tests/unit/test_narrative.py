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


import math

from charon.report.narrative import _render_adme_table, _render_ivive_audit


def test_format_value_nan_returns_dash():
    assert format_value(float("nan")) == "-"


def test_format_value_inf_returns_dash():
    assert format_value(float("inf")) == "-"
    assert format_value(float("-inf")) == "-"


def test_render_adme_table_header():
    data = _make_data(
        properties={
            "logp": {
                "value": 3.89,
                "ci_lower": 3.4,
                "ci_upper": 4.3,
                "unit": "log",
                "source": "ml_ensemble",
                "flag": None,
                "method": None,
            },
        }
    )
    out = _render_adme_table(data)
    assert "## 3. ADME Predictions" in out
    assert "| Property" in out
    assert "| --- " in out or "|---" in out
    assert "logp" in out
    assert "3.89" in out


def test_render_adme_table_skips_missing():
    data = _make_data(properties={})
    out = _render_adme_table(data)
    # Header still emitted
    assert "## 3. ADME Predictions" in out
    # But no rows
    assert "logp" not in out
    assert "fu_p" not in out


def test_render_adme_table_ci_dash_when_missing():
    data = _make_data(
        properties={
            "fu_p": {
                "value": 0.032,
                "ci_lower": None,
                "ci_upper": None,
                "unit": "fraction",
                "source": "experimental",
                "flag": None,
                "method": None,
            }
        }
    )
    out = _render_adme_table(data)
    # The row should have a dash in the CI column
    row_lines = [ln for ln in out.split("\n") if ln.startswith("| fu_p")]
    assert len(row_lines) == 1
    # CI column is the third pipe-separated field
    cells = [c.strip() for c in row_lines[0].split("|")[1:-1]]
    assert cells[2] == "-"


def test_render_ivive_audit_narrative():
    data = _make_data(
        ivive_summary={
            "clint_liver_L_h": 348.75,
            "cl_renal_L_h": 0.23,
            "fu_b": 0.071,
            "liver_model": "well_stirred",
            "compound_type": "base",
        },
        pk_params={"cl_apparent": 24.1},
    )
    out = _render_ivive_audit(data)
    assert "## 4. IVIVE" in out
    assert "well_stirred" in out
    assert "348.75" in out or "348.8" in out or "3.49e+02" in out
    assert "0.071" in out or "7.1e-02" in out
    assert "24.1" in out  # cl_apparent


def test_render_ivive_audit_handles_missing_metadata():
    data = _make_data(ivive_summary={}, pk_params={})
    out = _render_ivive_audit(data)
    assert "## 4. IVIVE" in out
    # Should not crash; contains placeholder
    assert "-" in out or "not available" in out.lower()


from charon.report.narrative import _md_cell


def test_md_cell_escapes_pipe():
    assert _md_cell("warn|bad") == "warn\\|bad"


def test_md_cell_none_returns_dash():
    assert _md_cell(None) == "-"


def test_md_cell_normalizes_newlines():
    assert _md_cell("line1\nline2") == "line1 line2"


def test_render_adme_table_escapes_pipe_in_flag():
    data = _make_data(
        properties={
            "logp": {
                "value": 3.89,
                "ci_lower": None,
                "ci_upper": None,
                "unit": "log",
                "source": "ml_ensemble",
                "flag": "warn|bad",
                "method": None,
            }
        }
    )
    out = _render_adme_table(data)
    assert "warn\\|bad" in out
    assert "| warn|bad |" not in out


from charon.report.narrative import _render_dose_projection, _render_pk_results


def test_render_pk_results_parameters_table():
    data = _make_data(
        pk_params={
            "cmax": 120.0,
            "tmax": 1.0,
            "auc_0_inf": 500.0,
            "auc_0_24": 480.0,
            "half_life": 3.5,
            "cl_apparent": 24.1,
            "vss": 95.0,
            "bioavailability": 0.44,
            "fa": 0.95,
            "fg": 0.57,
            "fh": 0.82,
        },
        pk_table=[
            {"time_h": 0.0, "cp_plasma_ug_L": 0.0, "cp_blood_ug_L": 0.0},
            {"time_h": 1.0, "cp_plasma_ug_L": 120.0, "cp_blood_ug_L": 108.0},
            {"time_h": 4.0, "cp_plasma_ug_L": 60.0, "cp_blood_ug_L": 54.0},
        ],
    )
    out = _render_pk_results(data)
    assert "## 5. PK Simulation Results" in out
    assert "120" in out  # Cmax
    assert "0.44" in out or "4.4e-01" in out  # F
    assert "Fa" in out and "Fg" in out and "Fh" in out
    # Cp-time table
    assert "| Time" in out
    assert "1" in out and "4" in out


def test_render_pk_results_handles_missing_params():
    data = _make_data(pk_params={"cmax": None, "auc_0_inf": None}, pk_table=[])
    out = _render_pk_results(data)
    assert "## 5. PK Simulation Results" in out
    # Dash for missing values
    assert "-" in out


def test_render_dose_projection_no_rec():
    data = _make_data()
    out = _render_dose_projection(data)
    assert "## 6. FIH Dose Projection" in out
    assert "not run" in out.lower() or "no " in out.lower()


def test_render_dose_projection_hed_only():
    data = _make_data(
        dose_recommendation={
            "mrsd_mg": 56.45,
            "limiting_method": "hed",
            "route": "oral",
            "safety_factor": 10.0,
            "salt_factor": 1.0,
            "rationale": "HED: 56.45 mg\nLimiting: hed",
            "hed": {"mrsd_mg": 56.45},
            "mabel": None,
            "pad": None,
        }
    )
    out = _render_dose_projection(data)
    assert "## 6. FIH Dose Projection" in out
    assert "| Method" in out
    # Three rows: HED, MABEL, PAD
    assert "HED" in out
    assert "MABEL" in out
    assert "PAD" in out
    assert "insufficient inputs" in out.lower()
    # Rationale present
    assert "Limiting: hed" in out
    # Limiting method marked
    assert "**56.45**" in out or "**56.5**" in out or "**5.6e+01**" in out


def test_md_cell_escapes_backslash_before_pipe():
    # a\\|b in the source, i.e. backslash + pipe, must not break table columns
    assert _md_cell("a\\|b") == "a\\\\\\|b"


def test_md_cell_normalizes_carriage_return():
    assert _md_cell("line1\r\nline2") == "line1  line2"


from charon.report.narrative import (
    _render_appendix,
    _render_limitations,
    _render_uncertainty,
    render_report,
)


def test_render_uncertainty_none_returns_empty():
    data = _make_data(uncertainty=None)
    out = _render_uncertainty(data)
    assert out == "" or out.strip() == ""


def test_render_uncertainty_populated():
    data = _make_data(
        uncertainty={
            "point_estimate_mg": 50.0,
            "ci_90_lower_mg": 20.0,
            "ci_90_upper_mg": 120.0,
            "ci_ratio": 6.0,
            "confidence": "MEDIUM",
            "n_samples": 100,
            "n_successful": 95,
            "convergence_met": True,
            "sensitivity": {"clint": 0.7, "fu_p": 0.2, "logp": 0.1},
            "limiting_parameter": "clint",
            "recommendation": "Experimental clint measurement would narrow CI by ~70%",
            "r_squared": 0.85,
        }
    )
    out = _render_uncertainty(data)
    assert "## 7. Uncertainty Analysis" in out
    assert "50" in out and "20" in out and "120" in out
    assert "MEDIUM" in out
    assert "clint" in out
    assert "70" in out  # 70.0%
    assert "0.85" in out  # R²


def test_render_uncertainty_low_r_squared_warning():
    data = _make_data(
        uncertainty={
            "point_estimate_mg": 50.0,
            "ci_90_lower_mg": 20.0,
            "ci_90_upper_mg": 120.0,
            "ci_ratio": 6.0,
            "confidence": "MEDIUM",
            "n_samples": 100,
            "n_successful": 95,
            "convergence_met": True,
            "sensitivity": {"clint": 0.5},
            "limiting_parameter": "clint",
            "recommendation": "",
            "r_squared": 0.5,
        }
    )
    out = _render_uncertainty(data)
    assert "non-linear" in out.lower() or "nonlinear" in out.lower() or "warn" in out.lower()


def test_render_limitations_always_has_boilerplate():
    data = _make_data()
    out = _render_limitations(data)
    assert "## 8. Limitations" in out
    assert "well-stirred" in out.lower() or "well stirred" in out.lower()
    assert "IR" in out or "immediate-release" in out.lower()


def test_render_limitations_appends_warnings():
    data = _make_data(warnings=["Applicability domain: LOW"])
    out = _render_limitations(data)
    assert "Applicability domain: LOW" in out


def test_render_limitations_appends_property_flags():
    data = _make_data(
        properties={
            "clint_uL_min_mg": {
                "value": 93.0,
                "ci_lower": None,
                "ci_upper": None,
                "unit": "uL/min/mg",
                "source": "ml_ensemble",
                "flag": "clint_tier2_ml",
                "method": None,
            }
        }
    )
    out = _render_limitations(data)
    assert "clint_tier2_ml" in out


def test_render_appendix_has_metadata_and_version():
    data = _make_data(metadata={"species": "human", "solver_method": "BDF"})
    out = _render_appendix(data)
    assert "## 9. Appendix" in out
    assert "0.1.0" in out  # charon_version
    assert "BDF" in out
    assert "2026-04-15" in out


def test_render_report_contains_all_sections():
    data = _make_data(
        properties={
            "logp": {
                "value": 3.89,
                "ci_lower": 3.4,
                "ci_upper": 4.3,
                "unit": "log",
                "source": "ml_ensemble",
                "flag": None,
                "method": None,
            }
        },
        ivive_summary={"liver_model": "well_stirred", "fu_b": 0.071},
        pk_params={"cmax": 120.0, "cl_apparent": 24.1, "auc_0_inf": 500.0},
        pk_table=[{"time_h": 1.0, "cp_plasma_ug_L": 120.0, "cp_blood_ug_L": 108.0}],
        dose_recommendation={
            "mrsd_mg": 56.45,
            "limiting_method": "hed",
            "route": "oral",
            "safety_factor": 10.0,
            "salt_factor": 1.0,
            "rationale": "ok",
            "hed": {"mrsd_mg": 56.45},
            "mabel": None,
            "pad": None,
        },
    )
    out = render_report(data)
    for section in (
        "# FIH Dose Rationale Report",
        "## 1. Executive Summary",
        "## 2. Compound Profile",
        "## 3. ADME Predictions",
        "## 4. IVIVE",
        "## 5. PK Simulation Results",
        "## 6. FIH Dose Projection",
        "## 8. Limitations",
        "## 9. Appendix",
    ):
        assert section in out, f"missing section: {section}"
    # No uncertainty section when uncertainty is None
    assert "## 7. Uncertainty Analysis" not in out


def test_render_report_includes_uncertainty_when_present():
    data = _make_data(
        uncertainty={
            "point_estimate_mg": 50.0,
            "ci_90_lower_mg": 20.0,
            "ci_90_upper_mg": 120.0,
            "ci_ratio": 6.0,
            "confidence": "MEDIUM",
            "n_samples": 100,
            "n_successful": 95,
            "convergence_met": True,
            "sensitivity": {"clint": 0.7},
            "limiting_parameter": "clint",
            "recommendation": "",
            "r_squared": 0.9,
        }
    )
    out = render_report(data)
    assert "## 7. Uncertainty Analysis" in out


def test_render_report_no_raw_none():
    data = _make_data()
    out = render_report(data)
    assert "| None |" not in out


def test_render_appendix_handles_none_metadata_value():
    data = _make_data(metadata={"species": "human", "seed": None})
    out = _render_appendix(data)
    # Should not contain raw "None" in the yaml block
    yaml_block_lines = [
        ln for ln in out.split("\n")
        if not ln.startswith("#")
        and not ln.startswith("- ")
        and not ln.startswith("```")
        and ln.strip()
    ]
    joined = "\n".join(yaml_block_lines)
    assert "None" not in joined  # yaml null, not Python None
    # species should still appear
    assert "species" in out


def test_render_appendix_handles_numpy_scalar():
    import numpy as np
    data = _make_data(metadata={"solver_nfev": np.int64(42), "residual": np.float64(1e-8)})
    out = _render_appendix(data)
    assert "42" in out
    # yaml should parse the block
    import yaml
    yaml_start = out.index("```yaml") + len("```yaml\n")
    yaml_end = out.rindex("```")
    yaml_text = out[yaml_start:yaml_end]
    parsed = yaml.safe_load(yaml_text)
    assert parsed["solver_nfev"] == 42


def test_render_ivive_audit_handles_none_liver_model():
    data = _make_data(
        ivive_summary={"liver_model": None, "fu_b": None, "clint_liver_L_h": None},
        pk_params={"cl_apparent": None},
    )
    out = _render_ivive_audit(data)
    assert "None" not in out


def test_render_dose_projection_bolds_limiting_case_insensitive():
    data = _make_data(
        dose_recommendation={
            "mrsd_mg": 56.45,
            "limiting_method": "HED",  # uppercase
            "route": "oral",
            "safety_factor": 10.0,
            "salt_factor": 1.0,
            "rationale": "HED",
            "hed": {"mrsd_mg": 56.45},
            "mabel": None,
            "pad": None,
        }
    )
    out = _render_dose_projection(data)
    # HED row MRSD should be bolded even though upstream gave "HED" not "hed"
    assert "**56.45**" in out or "**56.5**" in out
