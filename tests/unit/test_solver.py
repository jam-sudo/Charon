"""Unit tests for the PBPK solver wrapper.

Strategy:
  1. Smoke test: midazolam IV bolus runs to completion with finite output.
  2. Mass conservation: no elimination → total mass constant (numerical).
  3. Analytic 1-cpt equivalent: uniform Kp → late-time plateau matches
     dose / total_volume.
  4. BDF enforcement: simulate_iv rejects non-BDF methods.
"""

import numpy as np
import pytest

from charon.core.parameter_bridge import ParameterBridge
from charon.core.schema import (
    BindingProperties,
    CompoundConfig,
    CompoundProperties,
    MetabolismProperties,
    PhysicochemicalProperties,
    PredictedProperty,
    RenalProperties,
)
from charon.pbpk.ode_compiler import build_compound_pbpk_params
from charon.pbpk.solver import SimulationResult, simulate_iv
from charon.pbpk.topology import load_species_topology


def _predicted(value: float, unit: str = "") -> PredictedProperty:
    return PredictedProperty(value=float(value), source="experimental", unit=unit)


@pytest.fixture
def human():
    return load_species_topology("human")


@pytest.fixture
def midazolam_cfg():
    return CompoundConfig(
        name="midazolam",
        smiles="Cc1ncc2n1-c1ccc(Cl)cc1C(c1ccccc1F)=NC2",
        molecular_weight=325.77,
        source="experimental",
        properties=CompoundProperties(
            physicochemical=PhysicochemicalProperties(
                logp=_predicted(3.89),
                pka_base=_predicted(6.2),
            ),
            binding=BindingProperties(
                fu_p=_predicted(0.03, "fraction"),
                fu_inc=_predicted(0.96, "fraction"),
                bp_ratio=_predicted(0.66, "ratio"),
            ),
            metabolism=MetabolismProperties(
                clint_uL_min_mg=_predicted(93.0, "uL/min/mg"),
            ),
            renal=RenalProperties(
                clrenal_L_h=_predicted(0.0, "L/h"),
            ),
        ),
    )


@pytest.fixture
def midazolam_params(human, midazolam_cfg):
    return build_compound_pbpk_params(
        midazolam_cfg,
        human,
        bridge=ParameterBridge(),
        compound_type="base",
        override_cl_renal_L_h=0.0,
    )


class TestSimulateIvBolus:
    def test_smoke(self, human, midazolam_params):
        result = simulate_iv(
            human,
            midazolam_params,
            dose_mg=5.0,
            route="iv_bolus",
            duration_h=24.0,
        )
        assert isinstance(result, SimulationResult)
        assert result.time_h[0] == 0.0
        assert result.time_h[-1] == pytest.approx(24.0)
        assert np.all(np.isfinite(result.cp_plasma))
        assert result.cp_plasma[0] > 0.0, "Venous Cp should jump at t=0 from bolus"
        assert result.solver_success is True
        assert result.solver_method == "BDF"

    def test_bolus_mass_conservation_no_elim(self, human, midazolam_cfg):
        """Zero CLint + zero CLrenal → total drug mass conserved over time."""
        cfg = midazolam_cfg.model_copy(deep=True)
        cfg.properties.metabolism.clint_uL_min_mg = _predicted(0.0, "uL/min/mg")
        params = build_compound_pbpk_params(
            cfg, human, bridge=ParameterBridge(), compound_type="base",
            override_cl_renal_L_h=0.0,
        )
        result = simulate_iv(
            human, params, dose_mg=5.0, route="iv_bolus", duration_h=24.0
        )
        totals = result.state_trajectory.sum(axis=0)
        max_err = np.max(np.abs(totals - 5.0))
        assert max_err < 1e-3, f"Mass drift {max_err:.3e} mg exceeds tolerance"

    def test_mass_balance_residual_reported(self, human, midazolam_params):
        result = simulate_iv(
            human, midazolam_params, dose_mg=5.0,
            route="iv_bolus", duration_h=24.0,
        )
        assert np.isfinite(result.mass_balance_residual)

    def test_bdf_enforced(self, human, midazolam_params):
        with pytest.raises(ValueError, match="BDF"):
            simulate_iv(
                human, midazolam_params, dose_mg=5.0,
                route="iv_bolus", duration_h=24.0,
                method="RK45",
            )

    def test_bad_route_raises(self, human, midazolam_params):
        with pytest.raises(ValueError, match="route"):
            simulate_iv(
                human, midazolam_params, dose_mg=5.0,
                route="oral", duration_h=24.0,
            )

    def test_bad_dose_raises(self, human, midazolam_params):
        with pytest.raises(ValueError, match="dose_mg"):
            simulate_iv(
                human, midazolam_params, dose_mg=0.0,
                route="iv_bolus", duration_h=24.0,
            )

    def test_bad_duration_raises(self, human, midazolam_params):
        with pytest.raises(ValueError, match="duration_h"):
            simulate_iv(
                human, midazolam_params, dose_mg=5.0,
                route="iv_bolus", duration_h=0.0,
            )


class TestSimulateIvInfusion:
    def test_infusion_matches_bolus_at_long_time(self, human, midazolam_params):
        """A 1-minute infusion of 5 mg should give nearly the same 24 h
        profile as a 5 mg bolus (aside from the first few minutes)."""
        bolus = simulate_iv(
            human, midazolam_params, dose_mg=5.0,
            route="iv_bolus", duration_h=24.0,
        )
        infusion = simulate_iv(
            human, midazolam_params, dose_mg=5.0,
            route="iv_infusion", duration_h=24.0,
            infusion_duration_h=1.0 / 60.0,
        )
        idx = np.argmin(np.abs(bolus.time_h - 6.0))
        cp_b = bolus.cp_plasma[idx]
        cp_i = infusion.cp_plasma[idx]
        assert cp_i == pytest.approx(cp_b, rel=0.05)

    def test_infusion_requires_positive_duration(self, human, midazolam_params):
        with pytest.raises(ValueError, match="infusion_duration"):
            simulate_iv(
                human, midazolam_params, dose_mg=5.0,
                route="iv_infusion", duration_h=24.0,
                infusion_duration_h=0.0,
            )


class TestUniformKpPlateau:
    """Sanity: with Kp=1 and BP=1 everywhere and zero elimination, the
    long-time venous plasma concentration converges to dose / V_total.
    """

    def test_uniform_kp_gives_mass_over_total_volume(self, human, midazolam_cfg):
        cfg = midazolam_cfg.model_copy(deep=True)
        cfg.properties.metabolism.clint_uL_min_mg = _predicted(0.0, "uL/min/mg")
        params = build_compound_pbpk_params(
            cfg, human, bridge=ParameterBridge(), compound_type="base",
            override_cl_renal_L_h=0.0,
        )
        from dataclasses import replace
        kp_uniform = {name: 1.0 for name in params.kp_by_tissue}
        params_uniform = replace(
            params,
            bp_ratio=1.0,
            fu_b=params.fu_p / 1.0,
            kp_by_tissue=kp_uniform,
        )
        result = simulate_iv(
            human, params_uniform, dose_mg=5.0,
            route="iv_bolus", duration_h=48.0,
        )
        v_total = (
            human.venous_volume_L
            + human.arterial_volume_L
            + sum(n.volume_L for n in human.tissues.values())
        )
        expected_plateau = 5.0 / v_total
        cp_late = result.cp_plasma[-1]
        assert cp_late == pytest.approx(expected_plateau, rel=0.01)
