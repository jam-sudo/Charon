"""Tests for the oral simulation wrapper."""
from __future__ import annotations

import numpy as np
import pytest
import yaml
from pathlib import Path


def _make_midazolam_oral_params():
    """Build OralPBPKParams for midazolam from real YAML data."""
    from charon.core.parameter_bridge import ParameterBridge
    from charon.core.schema import CompoundConfig
    from charon.pbpk.acat import load_gi_tract
    from charon.pbpk.ode_compiler import (
        build_compound_pbpk_params,
        compute_gut_clint,
        OralPBPKParams,
    )
    from charon.pbpk.topology import load_species_topology

    comp_path = Path("validation/data/tier1_obach/compounds/midazolam.yaml")
    with comp_path.open() as f:
        raw = yaml.safe_load(f)
    compound = CompoundConfig(**raw)

    topo = load_species_topology("human")
    bridge = ParameterBridge()
    base_params = build_compound_pbpk_params(compound, topo, bridge=bridge)

    gi = load_gi_tract("human")
    q_gut = topo.tissues["gut_wall"].blood_flow_L_h
    q_villi = gi.q_villi_fraction * q_gut

    fm = compound.properties.metabolism.fm_cyp3a4
    clint_gut = compute_gut_clint(
        clint_liver_L_h=base_params.clint_liver_L_h,
        fm_cyp3a4=fm,
        gi_tract=gi,
        mppgl=40.0,
        liver_weight_g=topo.liver_weight_g,
    )

    peff = 4.0e-4

    return OralPBPKParams(
        name=base_params.name,
        molecular_weight=base_params.molecular_weight,
        logp=base_params.logp,
        pka_acid=base_params.pka_acid,
        pka_base=base_params.pka_base,
        compound_type=base_params.compound_type,
        fu_p=base_params.fu_p,
        bp_ratio=base_params.bp_ratio,
        fu_b=base_params.fu_b,
        clint_liver_L_h=base_params.clint_liver_L_h,
        cl_renal_L_h=base_params.cl_renal_L_h,
        kp_by_tissue=base_params.kp_by_tissue,
        kp_overrides=base_params.kp_overrides,
        clint_gut_L_h=clint_gut,
        peff_cm_s=peff,
        q_villi_L_h=q_villi,
        v_enterocyte_L=gi.enterocyte_volume_L,
        gi_tract=gi,
    ), topo


class TestSimulateOral:
    def test_returns_result(self):
        from charon.pbpk.solver import simulate_oral
        params, topo = _make_midazolam_oral_params()
        sim = simulate_oral(topo, params, dose_mg=5.0, duration_h=24.0)
        assert sim.route == "oral"
        assert sim.dose_mg == 5.0

    def test_cp_plasma_positive(self):
        from charon.pbpk.solver import simulate_oral
        params, topo = _make_midazolam_oral_params()
        sim = simulate_oral(topo, params, dose_mg=5.0, duration_h=24.0)
        assert np.any(sim.cp_plasma > 0)

    def test_cp_plasma_has_peak(self):
        """Oral profile should rise then fall → Cmax at tmax > 0."""
        from charon.pbpk.solver import simulate_oral
        params, topo = _make_midazolam_oral_params()
        sim = simulate_oral(topo, params, dose_mg=5.0, duration_h=24.0)
        tmax_idx = np.argmax(sim.cp_plasma)
        assert sim.time_h[tmax_idx] > 0.0
        assert sim.cp_plasma[tmax_idx] > 0.0

    def test_mass_balance(self):
        """Total mass at t=0 should equal dose."""
        from charon.pbpk.solver import simulate_oral
        params, topo = _make_midazolam_oral_params()
        sim = simulate_oral(topo, params, dose_mg=5.0, duration_h=24.0)
        total_t0 = sim.state_trajectory[:, 0].sum()
        assert total_t0 == pytest.approx(5.0, abs=0.01)

    def test_bdf_method(self):
        from charon.pbpk.solver import simulate_oral
        params, topo = _make_midazolam_oral_params()
        sim = simulate_oral(topo, params, dose_mg=5.0, duration_h=24.0)
        assert sim.solver_method == "BDF"

    def test_lumen_trajectory_available(self):
        from charon.pbpk.solver import simulate_oral
        params, topo = _make_midazolam_oral_params()
        sim = simulate_oral(topo, params, dose_mg=5.0, duration_h=24.0)
        assert hasattr(sim, "lumen_trajectory")
        assert sim.lumen_trajectory.shape[0] == 8

    def test_enterocyte_trajectory_available(self):
        from charon.pbpk.solver import simulate_oral
        params, topo = _make_midazolam_oral_params()
        sim = simulate_oral(topo, params, dose_mg=5.0, duration_h=24.0)
        assert hasattr(sim, "enterocyte_trajectory")
        assert len(sim.enterocyte_trajectory) == len(sim.time_h)

    def test_explicit_method_rejected(self):
        from charon.pbpk.solver import simulate_oral
        params, topo = _make_midazolam_oral_params()
        with pytest.raises(ValueError, match="BDF"):
            simulate_oral(topo, params, dose_mg=5.0, duration_h=24.0,
                          method="RK45")
