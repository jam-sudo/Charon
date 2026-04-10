"""Unit tests for PK parameter extraction from Cp-time profiles."""

import math

import numpy as np
import pytest

from charon.core.schema import PKParameters
from charon.pbpk.pk_extract import compute_pk_parameters


class TestAnalyticOneCompartmentBolus:
    """Mono-exponential IV bolus: C(t) = (D/V) * exp(-ke*t).

    With V = 10 L and CL = 5 L/h → ke = 0.5 h⁻¹, t½ = ln(2)/0.5 ≈ 1.386 h.
    A 100 mg bolus gives Cp(0) = 10 mg/L.

    Integrals:
      AUC_0_inf = D / CL = 100 / 5 = 20.0 mg·h/L
      AUMC_0_inf = D / (CL * ke) = 100 / (5 * 0.5) = 40.0 mg·h²/L
      MRT = AUMC / AUC = 2.0 h  (also = 1/ke)
      Vss = CL * MRT = 5 * 2 = 10.0 L (matches V)
    """

    @pytest.fixture
    def mono_exp(self):
        ke = 0.5
        dose = 100.0
        v = 10.0
        cp0 = dose / v
        t = np.linspace(0, 48.0, 2001)
        cp = cp0 * np.exp(-ke * t)
        return t, cp, dose, ke, v

    def test_cmax_equals_cp0(self, mono_exp):
        t, cp, dose, _, v = mono_exp
        pk = compute_pk_parameters(t, cp, dose_mg=dose, route="iv_bolus")
        assert pk.cmax == pytest.approx(dose / v, rel=1e-4)
        assert pk.tmax == pytest.approx(0.0, abs=1e-6)

    def test_auc_inf(self, mono_exp):
        t, cp, dose, ke, _ = mono_exp
        pk = compute_pk_parameters(t, cp, dose_mg=dose, route="iv_bolus")
        expected_auc = dose / 5.0  # CL = 5 L/h
        assert pk.auc_0_inf == pytest.approx(expected_auc, rel=5e-3)

    def test_half_life(self, mono_exp):
        t, cp, dose, ke, _ = mono_exp
        pk = compute_pk_parameters(t, cp, dose_mg=dose, route="iv_bolus")
        expected_t12 = math.log(2) / ke
        assert pk.half_life == pytest.approx(expected_t12, rel=1e-3)

    def test_cl(self, mono_exp):
        t, cp, dose, _, _ = mono_exp
        pk = compute_pk_parameters(t, cp, dose_mg=dose, route="iv_bolus")
        assert pk.cl_apparent == pytest.approx(5.0, rel=5e-3)

    def test_vss_equals_v(self, mono_exp):
        t, cp, dose, _, v = mono_exp
        pk = compute_pk_parameters(t, cp, dose_mg=dose, route="iv_bolus")
        # Vss = CL * MRT = V (for 1-compartment)
        assert pk.vss == pytest.approx(v, rel=5e-3)

    def test_auc_0_24(self, mono_exp):
        t, cp, dose, ke, _ = mono_exp
        pk = compute_pk_parameters(t, cp, dose_mg=dose, route="iv_bolus")
        # AUC(0,24) = (D/V/ke)*(1 - exp(-ke*24))
        expected = (dose / 10.0 / ke) * (1 - math.exp(-ke * 24.0))
        assert pk.auc_0_24 == pytest.approx(expected, rel=5e-3)


class TestInfusionBenetCorrection:
    """For IV infusion the Vss formula is Vss = CL * (MRT - T_inf/2)."""

    def test_zero_order_infusion_1_compartment(self):
        ke = 0.3
        v = 20.0
        cl = ke * v  # 6 L/h
        tinf = 2.0
        dose = 60.0
        rate = dose / tinf

        t = np.linspace(0.0, 48.0, 4801)
        c_inf = (rate / cl) * (1 - np.exp(-ke * np.minimum(t, tinf)))
        c_post = (
            (rate / cl)
            * (1 - np.exp(-ke * tinf))
            * np.exp(-ke * np.maximum(t - tinf, 0.0))
        )
        cp = np.where(t <= tinf, c_inf, c_post)

        pk = compute_pk_parameters(
            t, cp, dose_mg=dose, route="iv_infusion", infusion_duration_h=tinf
        )
        assert pk.cl_apparent == pytest.approx(cl, rel=5e-3)
        # Vss for 1-cpt infusion should still equal V (ignoring the
        # tiny numerical residual).
        assert pk.vss == pytest.approx(v, rel=5e-3)


class TestEdgeCases:
    def test_flat_profile_raises(self):
        t = np.linspace(0, 10, 101)
        cp = np.full_like(t, 1.0)
        with pytest.raises(ValueError, match="terminal slope"):
            compute_pk_parameters(t, cp, dose_mg=10.0, route="iv_bolus")

    def test_mismatched_array_lengths(self):
        with pytest.raises(ValueError, match="length"):
            compute_pk_parameters(
                np.array([0.0, 1.0, 2.0]),
                np.array([1.0, 2.0]),
                dose_mg=5.0,
                route="iv_bolus",
            )

    def test_short_tail_falls_back_cleanly(self):
        # 4 samples: log-linear fit is minimum-viable and must not crash.
        t = np.array([0.0, 1.0, 2.0, 3.0])
        cp = np.array([10.0, 5.0, 2.5, 1.25])  # clean exponential
        pk = compute_pk_parameters(t, cp, dose_mg=100.0, route="iv_bolus")
        assert pk.half_life is not None and pk.half_life > 0
