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
