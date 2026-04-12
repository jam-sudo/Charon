"""Tests for oral PBPK parameter assembly and gut CLint calculation."""
from __future__ import annotations

import dataclasses
import pytest


class TestComputeGutClint:
    def test_no_fm_returns_zero(self):
        """fm_cyp3a4=None → CLint_gut=0 (non-CYP3A4 substrate)."""
        from charon.pbpk.ode_compiler import compute_gut_clint
        from charon.pbpk.acat import load_gi_tract
        gi = load_gi_tract("human")
        result = compute_gut_clint(
            clint_liver_L_h=348.75,
            fm_cyp3a4=None,
            gi_tract=gi,
            mppgl=40.0,
            liver_weight_g=1500.0,
        )
        assert result == 0.0

    def test_fm_zero_returns_zero(self):
        """fm_cyp3a4=0 → CLint_gut=0."""
        from charon.pbpk.ode_compiler import compute_gut_clint
        from charon.pbpk.acat import load_gi_tract
        gi = load_gi_tract("human")
        result = compute_gut_clint(
            clint_liver_L_h=348.75,
            fm_cyp3a4=0.0,
            gi_tract=gi,
            mppgl=40.0,
            liver_weight_g=1500.0,
        )
        assert result == 0.0

    def test_midazolam_clint_gut(self):
        """Midazolam: CLint_gut ≈ 7.9 L/h → Fg ≈ 0.57.

        Hand calculation:
            cyp_ratio = 31/137 = 0.22628
            gut_scaling = 15 * 400 = 6000
            liver_scaling = 40 * 1500 = 60000
            CLint_gut = 348.75 * 1.0 * 0.22628 * (6000/60000)
                      = 348.75 * 0.22628 * 0.1
                      = 7.89 L/h
        """
        from charon.pbpk.ode_compiler import compute_gut_clint
        from charon.pbpk.acat import load_gi_tract
        gi = load_gi_tract("human")
        result = compute_gut_clint(
            clint_liver_L_h=348.75,
            fm_cyp3a4=1.0,
            gi_tract=gi,
            mppgl=40.0,
            liver_weight_g=1500.0,
        )
        # Expected: 348.75 × 1.0 × (31/137) × (15×400)/(40×1500) = 7.89
        assert result == pytest.approx(7.89, rel=0.02)

    def test_partial_fm(self):
        """fm_cyp3a4=0.5 → half the CLint_gut."""
        from charon.pbpk.ode_compiler import compute_gut_clint
        from charon.pbpk.acat import load_gi_tract
        gi = load_gi_tract("human")
        full = compute_gut_clint(
            clint_liver_L_h=348.75, fm_cyp3a4=1.0,
            gi_tract=gi, mppgl=40.0, liver_weight_g=1500.0,
        )
        half = compute_gut_clint(
            clint_liver_L_h=348.75, fm_cyp3a4=0.5,
            gi_tract=gi, mppgl=40.0, liver_weight_g=1500.0,
        )
        assert half == pytest.approx(full / 2.0, rel=1e-10)


class TestOralPBPKParams:
    def test_has_gut_fields(self):
        from charon.pbpk.ode_compiler import OralPBPKParams
        field_names = [f.name for f in dataclasses.fields(OralPBPKParams)]
        assert "clint_gut_L_h" in field_names
        assert "peff_cm_s" in field_names
        assert "q_villi_L_h" in field_names
        assert "v_enterocyte_L" in field_names

    def test_inherits_compound_fields(self):
        from charon.pbpk.ode_compiler import OralPBPKParams
        field_names = [f.name for f in dataclasses.fields(OralPBPKParams)]
        assert "fu_b" in field_names
        assert "clint_liver_L_h" in field_names
        assert "kp_by_tissue" in field_names


import numpy as np


class TestBuildOralRhs:
    """Test the oral ODE right-hand side function."""

    def _make_oral_params(self):
        """Helper: build minimal OralPBPKParams for midazolam-like compound."""
        from charon.pbpk.ode_compiler import OralPBPKParams
        from charon.pbpk.acat import load_gi_tract
        from charon.pbpk.topology import load_species_topology

        topo = load_species_topology("human")
        gi = load_gi_tract("human")

        kp = {name: 1.0 for name in topo.tissue_names()}
        q_gut = topo.tissues["gut_wall"].blood_flow_L_h
        q_villi = gi.q_villi_fraction * q_gut

        return OralPBPKParams(
            name="test_oral",
            molecular_weight=325.77,
            logp=3.89,
            pka_acid=None,
            pka_base=6.2,
            compound_type="base",
            fu_p=0.03,
            bp_ratio=0.66,
            fu_b=0.03 / 0.66,
            clint_liver_L_h=348.75,
            cl_renal_L_h=0.0,
            kp_by_tissue=kp,
            clint_gut_L_h=7.89,
            peff_cm_s=4.0e-4,
            q_villi_L_h=q_villi,
            v_enterocyte_L=gi.enterocyte_volume_L,
            gi_tract=gi,
        ), topo

    def test_returns_callable(self):
        from charon.pbpk.ode_compiler import build_oral_rhs
        params, topo = self._make_oral_params()
        rhs = build_oral_rhs(topo, params)
        assert callable(rhs)

    def test_state_vector_length(self):
        """Oral state vector = 2 + 15 tissues + 8 lumen + 1 enterocyte = 26."""
        from charon.pbpk.ode_compiler import build_oral_rhs
        params, topo = self._make_oral_params()
        rhs = build_oral_rhs(topo, params)
        n = 2 + len(topo.tissues) + 8 + 1
        y0 = np.zeros(n)
        y0[17] = 5.0  # dose in stomach (index 2+15=17)
        dy = rhs(0.0, y0)
        assert len(dy) == n

    def test_dose_in_stomach_creates_transit(self):
        """Drug in stomach should produce negative dA_stomach (emptying)."""
        from charon.pbpk.ode_compiler import build_oral_rhs
        params, topo = self._make_oral_params()
        rhs = build_oral_rhs(topo, params)
        n = 2 + len(topo.tissues) + 8 + 1
        y0 = np.zeros(n)
        y0[17] = 5.0
        dy = rhs(0.0, y0)
        assert dy[17] < 0.0  # stomach emptying
        assert dy[18] > 0.0  # duodenum receiving

    def test_enterocyte_receives_absorption(self):
        """Drug in duodenum should produce positive dA_enterocyte."""
        from charon.pbpk.ode_compiler import build_oral_rhs
        params, topo = self._make_oral_params()
        rhs = build_oral_rhs(topo, params)
        n = 2 + len(topo.tissues) + 8 + 1
        y0 = np.zeros(n)
        y0[18] = 5.0  # duodenum
        dy = rhs(0.0, y0)
        assert dy[25] > 0.0  # enterocyte (index 2+15+8=25)

    def test_mass_conservation_no_elimination(self):
        """With CLint=0 and CLrenal=0, total mass derivative should be 0 at t=0 with dose in stomach."""
        from charon.pbpk.ode_compiler import OralPBPKParams, build_oral_rhs
        from charon.pbpk.acat import load_gi_tract
        from charon.pbpk.topology import load_species_topology

        topo = load_species_topology("human")
        gi = load_gi_tract("human")
        q_gut = topo.tissues["gut_wall"].blood_flow_L_h

        params = OralPBPKParams(
            name="no_elim",
            molecular_weight=300.0,
            logp=2.0,
            pka_acid=None,
            pka_base=None,
            compound_type="neutral",
            fu_p=1.0,
            bp_ratio=1.0,
            fu_b=1.0,
            clint_liver_L_h=0.0,
            cl_renal_L_h=0.0,
            kp_by_tissue={n: 1.0 for n in topo.tissue_names()},
            clint_gut_L_h=0.0,
            peff_cm_s=4.0e-4,
            q_villi_L_h=gi.q_villi_fraction * q_gut,
            v_enterocyte_L=gi.enterocyte_volume_L,
            gi_tract=gi,
        )

        rhs = build_oral_rhs(topo, params)
        n = 2 + len(topo.tissues) + 8 + 1
        y0 = np.zeros(n)
        y0[17] = 100.0  # 100 mg in stomach

        dy = rhs(0.0, y0)
        total_dy = np.sum(dy)
        assert abs(total_dy) < 1e-10, f"Mass not conserved: sum(dy) = {total_dy}"

    def test_liver_receives_portal_from_enterocyte(self):
        """Enterocyte drug should increase liver mass via portal inflow."""
        from charon.pbpk.ode_compiler import build_oral_rhs
        params, topo = self._make_oral_params()
        rhs = build_oral_rhs(topo, params)
        n = 2 + len(topo.tissues) + 8 + 1
        y0 = np.zeros(n)
        y0[25] = 5.0  # drug in enterocyte
        dy = rhs(0.0, y0)
        liver_idx = 2 + list(topo.tissues.keys()).index("liver")
        assert dy[liver_idx] > 0.0, "liver should receive drug from enterocyte"
