"""Tests for oral PK parameter extraction including Fa/Fg/Fh."""
from __future__ import annotations

import math
import numpy as np
import pytest
import yaml
from pathlib import Path


def _run_midazolam_oral():
    """Run midazolam oral simulation and return (sim, params, topo)."""
    from charon.core.parameter_bridge import ParameterBridge
    from charon.core.schema import CompoundConfig
    from charon.pbpk.acat import load_gi_tract
    from charon.pbpk.ode_compiler import (
        build_compound_pbpk_params,
        compute_gut_clint,
        OralPBPKParams,
    )
    from charon.pbpk.solver import simulate_oral
    from charon.pbpk.topology import load_species_topology

    comp_path = Path("validation/data/tier1_obach/compounds/midazolam.yaml")
    with comp_path.open() as f:
        raw = yaml.safe_load(f)
    compound = CompoundConfig(**raw)

    topo = load_species_topology("human")
    bridge = ParameterBridge()
    base = build_compound_pbpk_params(compound, topo, bridge=bridge)

    gi = load_gi_tract("human")
    q_gut = topo.tissues["gut_wall"].blood_flow_L_h
    q_villi = gi.q_villi_fraction * q_gut

    fm = compound.properties.metabolism.fm_cyp3a4
    clint_gut = compute_gut_clint(
        clint_liver_L_h=base.clint_liver_L_h,
        fm_cyp3a4=fm,
        gi_tract=gi,
        mppgl=40.0,
        liver_weight_g=topo.liver_weight_g,
    )

    params = OralPBPKParams(
        name=base.name, molecular_weight=base.molecular_weight,
        logp=base.logp, pka_acid=base.pka_acid, pka_base=base.pka_base,
        compound_type=base.compound_type, fu_p=base.fu_p,
        bp_ratio=base.bp_ratio, fu_b=base.fu_b,
        clint_liver_L_h=base.clint_liver_L_h, cl_renal_L_h=base.cl_renal_L_h,
        kp_by_tissue=base.kp_by_tissue, kp_overrides=base.kp_overrides,
        clint_gut_L_h=clint_gut, peff_cm_s=4.0e-4,
        q_villi_L_h=q_villi, v_enterocyte_L=gi.enterocyte_volume_L,
        gi_tract=gi,
    )

    sim = simulate_oral(topo, params, dose_mg=5.0, duration_h=72.0)
    return sim, params, topo


class TestComputeOralPK:
    def test_returns_pk_parameters(self):
        from charon.pbpk.pk_extract import compute_oral_pk_parameters
        sim, params, topo = _run_midazolam_oral()
        pk = compute_oral_pk_parameters(sim, params, topo, dose_mg=5.0)
        assert pk.cmax is not None
        assert pk.tmax is not None
        assert pk.auc_0_inf is not None

    def test_cmax_positive(self):
        from charon.pbpk.pk_extract import compute_oral_pk_parameters
        sim, params, topo = _run_midazolam_oral()
        pk = compute_oral_pk_parameters(sim, params, topo, dose_mg=5.0)
        assert pk.cmax > 0

    def test_tmax_positive(self):
        from charon.pbpk.pk_extract import compute_oral_pk_parameters
        sim, params, topo = _run_midazolam_oral()
        pk = compute_oral_pk_parameters(sim, params, topo, dose_mg=5.0)
        assert pk.tmax > 0

    def test_fa_high_for_midazolam(self):
        """Midazolam has high Peff -> Fa > 0.90."""
        from charon.pbpk.pk_extract import compute_oral_pk_parameters
        sim, params, topo = _run_midazolam_oral()
        pk = compute_oral_pk_parameters(sim, params, topo, dose_mg=5.0)
        assert pk.fa is not None
        assert pk.fa > 0.90

    def test_fg_midazolam_near_057(self):
        """Midazolam Fg should be approx 0.57 (calibration check, +/-10%)."""
        from charon.pbpk.pk_extract import compute_oral_pk_parameters
        sim, params, topo = _run_midazolam_oral()
        pk = compute_oral_pk_parameters(sim, params, topo, dose_mg=5.0)
        assert pk.fg is not None
        assert 0.50 < pk.fg < 0.65, f"Fg = {pk.fg:.3f}, expected ~0.57"

    def test_fh_positive(self):
        from charon.pbpk.pk_extract import compute_oral_pk_parameters
        sim, params, topo = _run_midazolam_oral()
        pk = compute_oral_pk_parameters(sim, params, topo, dose_mg=5.0)
        assert pk.fh is not None
        assert 0.0 < pk.fh < 1.0

    def test_bioavailability_product(self):
        """F = Fa x Fg x Fh."""
        from charon.pbpk.pk_extract import compute_oral_pk_parameters
        sim, params, topo = _run_midazolam_oral()
        pk = compute_oral_pk_parameters(sim, params, topo, dose_mg=5.0)
        expected = pk.fa * pk.fg * pk.fh
        assert pk.bioavailability == pytest.approx(expected, rel=1e-6)

    def test_cl_apparent_is_cl_over_f(self):
        """For oral, CL_apparent = dose / AUC = CL / F."""
        from charon.pbpk.pk_extract import compute_oral_pk_parameters
        sim, params, topo = _run_midazolam_oral()
        pk = compute_oral_pk_parameters(sim, params, topo, dose_mg=5.0)
        assert pk.cl_apparent == pytest.approx(5.0 / pk.auc_0_inf, rel=1e-3)

    def test_vss_none_for_oral(self):
        from charon.pbpk.pk_extract import compute_oral_pk_parameters
        sim, params, topo = _run_midazolam_oral()
        pk = compute_oral_pk_parameters(sim, params, topo, dose_mg=5.0)
        assert pk.vss is None
