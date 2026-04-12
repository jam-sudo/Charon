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
