"""Unit tests for charon.cli.main (invoked directly for speed)."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from charon.cli.main import main


def test_main_help_lists_subcommands(capsys):
    with pytest.raises(SystemExit) as exc:
        main(["--help"])
    assert exc.value.code == 0
    out = capsys.readouterr().out
    for sub in ("predict", "simulate", "translate", "recommend", "report"):
        assert sub in out


def test_main_no_args_exits_nonzero(capsys):
    rc = main([])
    assert rc != 0


def test_predict_subcommand_basic(capsys, monkeypatch):
    """`charon predict <smiles>` calls predict_properties and prints a table."""
    from charon.core.schema import (
        BindingProperties,
        CompoundProperties,
        PhysicochemicalProperties,
        PredictedProperty,
    )

    def fake_predict(smiles, **_):
        return CompoundProperties(
            physicochemical=PhysicochemicalProperties(
                logp=PredictedProperty(value=0.5, source="ml_ensemble", unit="log"),
            ),
            binding=BindingProperties(
                fu_p=PredictedProperty(
                    value=0.2, source="ml_ensemble", unit="fraction"
                ),
            ),
        )

    monkeypatch.setattr("charon.cli.main.predict_properties", fake_predict)
    rc = main(["predict", "CCO"])
    assert rc == 0
    out = capsys.readouterr().out
    assert "logp" in out
    assert "0.5" in out
    assert "fu_p" in out


def test_predict_json_flag(capsys, monkeypatch):
    from charon.core.schema import CompoundProperties

    monkeypatch.setattr(
        "charon.cli.main.predict_properties",
        lambda smiles, **_: CompoundProperties(),
    )
    rc = main(["predict", "CCO", "--json"])
    assert rc == 0
    out = capsys.readouterr().out
    parsed = json.loads(out)
    assert isinstance(parsed, dict)


def test_predict_invalid_smiles(capsys, monkeypatch):
    def boom(smiles, **_):
        raise ValueError("bad SMILES")

    monkeypatch.setattr("charon.cli.main.predict_properties", boom)
    rc = main(["predict", "not-a-smiles"])
    assert rc == 1
    err = capsys.readouterr().err
    assert "Invalid" in err or "bad SMILES" in err


def test_simulate_subcommand_runs_pipeline(capsys, monkeypatch):
    """`charon simulate` runs Pipeline.from_smiles(...).run() and prints PK."""
    import numpy as np

    from charon.core.schema import (
        CompoundConfig,
        CompoundProperties,
        PKParameters,
    )
    from charon.pbpk.solver import SimulationResult
    from charon.pipeline import PipelineResult

    def fake_from_smiles(smiles, **kwargs):
        class _FakePipeline:
            def run(self):
                t = np.array([0.0, 1.0, 2.0])
                cp = np.array([5.0, 2.5, 1.25])
                return PipelineResult(
                    compound=CompoundConfig(
                        name=smiles,
                        smiles=smiles,
                        molecular_weight=46.07,
                        source="predicted",
                        properties=CompoundProperties(),
                    ),
                    pk_parameters=PKParameters(cmax=5.0, auc_0_inf=12.0),
                    time_h=t,
                    cp_plasma=cp,
                    cp_blood=cp,
                    simulation=SimulationResult(
                        time_h=t,
                        cp_plasma=cp,
                        cp_blood=cp,
                        state_trajectory=np.zeros((3, t.size)),
                        mass_balance_residual=0.0,
                        solver_success=True,
                        solver_method="BDF",
                        solver_nfev=7,
                        route=kwargs.get("route", "iv_bolus"),
                        dose_mg=kwargs.get("dose_mg", 10.0),
                        infusion_duration_h=kwargs.get("infusion_duration_h", 0.0),
                    ),
                    metadata={
                        "route": kwargs.get("route"),
                        "dose_mg": kwargs.get("dose_mg"),
                        "duration_h": kwargs.get("duration_h", 72.0),
                        "liver_model": "well_stirred",
                    },
                )
        return _FakePipeline()

    monkeypatch.setattr(
        "charon.cli.main.Pipeline.from_smiles", staticmethod(fake_from_smiles)
    )
    rc = main(["simulate", "CCO", "--route", "iv_bolus", "--dose", "10"])
    assert rc == 0
    out = capsys.readouterr().out
    assert "Cmax" in out or "cmax" in out
    assert "5" in out


def _install_fake_pipeline_with_dose(monkeypatch, *, uncertainty=None):
    """Install a fake Pipeline.from_smiles that returns a PipelineResult with
    a dose_recommendation (and optional uncertainty)."""
    import numpy as np

    from charon.core.schema import CompoundConfig, CompoundProperties, PKParameters
    from charon.pbpk.solver import SimulationResult
    from charon.pipeline import PipelineResult
    from charon.translational.dose_projector import FIHDoseRecommendation
    from charon.translational.hed import HEDResult

    def fake_from_smiles(smiles, **kwargs):
        class _FakePipe:
            def __init__(self):
                self.dose_projection = kwargs.get("dose_projection")
                self.uncertainty = kwargs.get("uncertainty")

            def run(self):
                t = np.array([0.0, 1.0, 2.0])
                cp = np.array([5.0, 2.5, 1.25])
                rec = None
                if self.dose_projection is not None:
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
                        route=kwargs.get("route", "oral"),
                        rationale="HED 56.45 mg",
                    )
                return PipelineResult(
                    compound=CompoundConfig(
                        name=smiles,
                        smiles=smiles,
                        molecular_weight=46.07,
                        source="predicted",
                        properties=CompoundProperties(),
                    ),
                    pk_parameters=PKParameters(
                        cmax=5.0, auc_0_inf=12.0, half_life=2.0, cl_apparent=0.2
                    ),
                    time_h=t,
                    cp_plasma=cp,
                    cp_blood=cp,
                    simulation=SimulationResult(
                        time_h=t, cp_plasma=cp, cp_blood=cp,
                        state_trajectory=np.zeros((3, t.size)),
                        mass_balance_residual=0.0,
                        solver_success=True,
                        solver_method="BDF",
                        solver_nfev=7,
                        route=kwargs.get("route", "oral"),
                        dose_mg=kwargs.get("dose_mg", 5.0),
                        infusion_duration_h=kwargs.get("infusion_duration_h", 0.0),
                    ),
                    metadata={
                        "route": kwargs.get("route"),
                        "dose_mg": kwargs.get("dose_mg"),
                        "duration_h": kwargs.get("duration_h", 72.0),
                        "liver_model": "well_stirred",
                    },
                    dose_recommendation=rec,
                    uncertainty=uncertainty,
                )
        return _FakePipe()

    monkeypatch.setattr(
        "charon.cli.main.Pipeline.from_smiles", staticmethod(fake_from_smiles)
    )


def test_translate_requires_at_least_one_target(capsys, monkeypatch):
    _install_fake_pipeline_with_dose(monkeypatch)
    rc = main(["translate", "CCO", "--route", "oral", "--dose", "5"])
    # No NOAEL / Kd / Ceff -> should fail
    assert rc != 0
    err = capsys.readouterr().err
    assert (
        "noael" in err.lower() or "target" in err.lower() or "required" in err.lower()
    )


def test_translate_with_noael(capsys, monkeypatch):
    _install_fake_pipeline_with_dose(monkeypatch)
    rc = main(
        [
            "translate",
            "CCO",
            "--route", "oral",
            "--dose", "5",
            "--noael", "50",
            "--noael-species", "rat",
        ]
    )
    assert rc == 0
    out = capsys.readouterr().out
    assert "56.45" in out or "56.5" in out
    assert "HED" in out


def test_recommend_without_uncertainty(capsys, monkeypatch):
    _install_fake_pipeline_with_dose(monkeypatch)
    rc = main(
        [
            "recommend",
            "CCO",
            "--route", "oral",
            "--dose", "5",
            "--noael", "50",
            "--noael-species", "rat",
        ]
    )
    assert rc == 0
    out = capsys.readouterr().out
    assert "56.45" in out or "56.5" in out


def test_recommend_with_uncertainty_flag(capsys, monkeypatch):
    from charon.uncertainty.dose_range import UncertaintyResult

    unc = UncertaintyResult(
        point_estimate_mg=56.0,
        ci_90_lower_mg=25.0,
        ci_90_upper_mg=125.0,
        ci_ratio=5.0,
        confidence="MEDIUM",
        n_samples=50,
        n_successful=48,
        convergence_met=False,
        sensitivity={"clint_uL_min_mg": 0.6, "fu_p": 0.3, "logp": 0.1},
        limiting_parameter="clint_uL_min_mg",
        recommendation="Experimental CLint measurement would narrow CI by ~60%",
        r_squared=0.82,
    )
    _install_fake_pipeline_with_dose(monkeypatch, uncertainty=unc)
    rc = main(
        [
            "recommend",
            "CCO",
            "--route", "oral",
            "--dose", "5",
            "--noael", "50",
            "--noael-species", "rat",
            "--uncertainty",
            "--n-samples", "50",
        ]
    )
    assert rc == 0
    out = capsys.readouterr().out
    assert "25" in out and "125" in out
    assert "MEDIUM" in out


def test_report_writes_md_and_json(tmp_path, monkeypatch):
    _install_fake_pipeline_with_dose(monkeypatch)
    out_base = tmp_path / "run"
    rc = main(
        [
            "report",
            "CCO",
            "--route", "oral",
            "--dose", "5",
            "--noael", "50",
            "--noael-species", "rat",
            "--output", str(out_base) + ".md",
        ]
    )
    assert rc == 0
    md = tmp_path / "run.md"
    js = tmp_path / "run.json"
    assert md.exists()
    assert js.exists()
    md_text = md.read_text()
    assert "# FIH Dose Rationale Report" in md_text
    assert "CCO" in md_text
    parsed = json.loads(js.read_text())
    assert parsed["route"] == "oral"
    assert parsed["dose_mg"] == pytest.approx(5.0)


def test_report_accepts_path_without_suffix(tmp_path, monkeypatch):
    _install_fake_pipeline_with_dose(monkeypatch)
    out_base = tmp_path / "no_suffix"
    rc = main(
        [
            "report",
            "CCO",
            "--route", "oral",
            "--dose", "5",
            "--noael", "50",
            "--noael-species", "rat",
            "--output", str(out_base),
        ]
    )
    assert rc == 0
    assert (tmp_path / "no_suffix.md").exists()
    assert (tmp_path / "no_suffix.json").exists()


def test_report_requires_output(capsys, monkeypatch):
    _install_fake_pipeline_with_dose(monkeypatch)
    with pytest.raises(SystemExit):
        # argparse exits for missing required --output
        main(
            [
                "report", "CCO",
                "--route", "oral", "--dose", "5",
                "--noael", "50", "--noael-species", "rat",
            ]
        )
