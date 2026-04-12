"""Tests for the ACAT GI tract module."""
from __future__ import annotations

import math
import pytest


class TestGISegment:
    def test_segment_fields(self):
        from charon.pbpk.acat import GISegment
        s = GISegment(
            name="duodenum", volume_L=0.05, radius_cm=1.6,
            ka_fraction=1.0, transit_rate_1_h=3.846,
        )
        assert s.name == "duodenum"
        assert s.volume_L == 0.05
        assert s.radius_cm == 1.6
        assert s.ka_fraction == 1.0
        assert s.transit_rate_1_h == 3.846

    def test_segment_frozen(self):
        from charon.pbpk.acat import GISegment
        s = GISegment(name="x", volume_L=1, radius_cm=1, ka_fraction=1, transit_rate_1_h=1)
        with pytest.raises(AttributeError):
            s.name = "y"


class TestGITract:
    def test_gi_tract_fields(self):
        from charon.pbpk.acat import GITract, GISegment
        seg = GISegment(name="stomach", volume_L=0.25, radius_cm=5.0,
                        ka_fraction=0.0, transit_rate_1_h=4.0)
        gi = GITract(
            segments=(seg,),
            enterocyte_volume_L=0.30,
            enterocyte_weight_g=400,
            q_villi_fraction=0.18,
            mppgi_mg_g=15.0,
            cyp3a4_gut_pmol_per_mg=31,
            cyp3a4_liver_pmol_per_mg=137,
        )
        assert len(gi.segments) == 1
        assert gi.mppgi_mg_g == 15.0

    def test_gi_tract_frozen(self):
        from charon.pbpk.acat import GITract, GISegment
        seg = GISegment(name="x", volume_L=1, radius_cm=1, ka_fraction=1, transit_rate_1_h=1)
        gi = GITract(segments=(seg,), enterocyte_volume_L=0.3,
                     enterocyte_weight_g=400, q_villi_fraction=0.18,
                     mppgi_mg_g=15.0, cyp3a4_gut_pmol_per_mg=31,
                     cyp3a4_liver_pmol_per_mg=137)
        with pytest.raises(AttributeError):
            gi.mppgi_mg_g = 99


class TestLoadGITract:
    def test_loads_8_segments(self):
        from charon.pbpk.acat import load_gi_tract
        gi = load_gi_tract("human")
        assert len(gi.segments) == 8

    def test_segment_names_in_order(self):
        from charon.pbpk.acat import load_gi_tract
        gi = load_gi_tract("human")
        names = [s.name for s in gi.segments]
        assert names == [
            "stomach", "duodenum", "jejunum1", "jejunum2",
            "ileum1", "ileum2", "ileum3", "colon",
        ]

    def test_stomach_ka_zero(self):
        from charon.pbpk.acat import load_gi_tract
        gi = load_gi_tract("human")
        assert gi.segments[0].ka_fraction == 0.0

    def test_enterocyte_params(self):
        from charon.pbpk.acat import load_gi_tract
        gi = load_gi_tract("human")
        assert gi.enterocyte_volume_L == 0.30
        assert gi.enterocyte_weight_g == 400
        assert gi.q_villi_fraction == 0.18
        assert gi.mppgi_mg_g == 15.0

    def test_total_sitt_about_3h(self):
        """Sum of 1/k_transit for SI segments ≈ 3.25 h."""
        from charon.pbpk.acat import load_gi_tract
        gi = load_gi_tract("human")
        si = [s for s in gi.segments if s.name not in ("stomach", "colon")]
        sitt = sum(1.0 / s.transit_rate_1_h for s in si)
        assert 3.0 < sitt < 4.0, f"SITT = {sitt:.2f} h"

    def test_missing_species_raises(self):
        from charon.pbpk.acat import load_gi_tract
        with pytest.raises(FileNotFoundError):
            load_gi_tract("martian")


class TestComputeAbsorptionRates:
    def test_stomach_kabs_zero(self):
        """ka_fraction=0 → k_abs=0 regardless of Peff."""
        from charon.pbpk.acat import load_gi_tract, compute_absorption_rates
        gi = load_gi_tract("human")
        rates = compute_absorption_rates(gi, peff_cm_s=4.0e-4)
        assert rates[0] == 0.0  # stomach

    def test_duodenum_rate_positive(self):
        from charon.pbpk.acat import load_gi_tract, compute_absorption_rates
        gi = load_gi_tract("human")
        rates = compute_absorption_rates(gi, peff_cm_s=4.0e-4)
        assert rates[1] > 0.0  # duodenum

    def test_hand_calculation_duodenum(self):
        """k_abs = 2 × Peff × 3600 / R × ka_fraction."""
        from charon.pbpk.acat import load_gi_tract, compute_absorption_rates
        gi = load_gi_tract("human")
        rates = compute_absorption_rates(gi, peff_cm_s=4.0e-4)
        expected = 2 * 4.0e-4 * 3600 / 1.6 * 1.0  # = 1.8 /h
        assert rates[1] == pytest.approx(expected, rel=1e-6)

    def test_rates_decrease_distally(self):
        """Due to ka_fraction decay, k_abs should generally decrease."""
        from charon.pbpk.acat import load_gi_tract, compute_absorption_rates
        gi = load_gi_tract("human")
        rates = compute_absorption_rates(gi, peff_cm_s=4.0e-4)
        assert rates[2] > rates[6]  # jejunum1 > ileum3

    def test_returns_8_values(self):
        from charon.pbpk.acat import load_gi_tract, compute_absorption_rates
        gi = load_gi_tract("human")
        rates = compute_absorption_rates(gi, peff_cm_s=1.0e-4)
        assert len(rates) == 8


class TestPappToPeff:
    def test_known_conversion(self):
        """Sun 2002 correlation: log10(Peff) = 0.4926×log10(Papp_nm_s) - 0.1454."""
        from charon.pbpk.acat import papp_to_peff
        import math
        papp = 100.0  # nm/s
        expected = 10 ** (0.4926 * math.log10(100.0) - 0.1454)
        result = papp_to_peff(papp)
        assert result == pytest.approx(expected, rel=1e-6)

    def test_high_papp(self):
        from charon.pbpk.acat import papp_to_peff
        result = papp_to_peff(500.0)  # high permeability
        assert result > 1e-5  # Peff should be > 10 µm/s

    def test_zero_papp_raises(self):
        from charon.pbpk.acat import papp_to_peff
        with pytest.raises(ValueError):
            papp_to_peff(0.0)
