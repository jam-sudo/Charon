"""End-to-end oral pipeline tests with Fg validation.

Tests midazolam (calibration), felodipine, and nifedipine (independent)
against literature Fg values per spec section 10.

Literature Fg references:
- midazolam  Fg = 0.57  (Thummel 1996 Clin Pharmacol Ther 59(5):491;
                          Gorski 1998 Clin Pharmacol Ther 64(2):133)
- felodipine Fg = 0.45  (Edgar 1992 Eur J Clin Pharmacol 42(3):261)
- nifedipine Fg = 0.78  (Holtbecker 1996 Clin Pharmacol Ther 60(1):54)
"""

from __future__ import annotations

import math
from pathlib import Path

import pytest
import yaml

from charon import Pipeline, PipelineResult
from charon.core.schema import CompoundConfig

REPO_ROOT = Path(__file__).resolve().parents[2]


def _load_compound(name: str) -> CompoundConfig:
    """Load a compound from the tier-1 Obach validation YAML."""
    path = (
        REPO_ROOT
        / "validation"
        / "data"
        / "tier1_obach"
        / "compounds"
        / f"{name}.yaml"
    )
    with path.open() as f:
        raw = yaml.safe_load(f)
    return CompoundConfig(**raw)


def _run_oral(name: str, dose_mg: float, duration_h: float = 72.0) -> PipelineResult:
    """Run oral pipeline for a named compound."""
    compound = _load_compound(name)
    pipe = Pipeline(
        compound=compound,
        route="oral",
        dose_mg=dose_mg,
        duration_h=duration_h,
    )
    return pipe.run()


class TestMidazolamOral:
    """Midazolam PO 5 mg -- calibration compound."""

    def test_produces_finite_pk(self):
        result = _run_oral("midazolam", 5.0)
        pk = result.pk_parameters
        for attr in ("cmax", "tmax", "auc_0_inf", "half_life", "fa", "fg", "fh"):
            val = getattr(pk, attr)
            assert val is not None and math.isfinite(val), f"{attr} = {val}"

    def test_fg_calibration(self, record_property):
        """Fg within 10% of 0.57 (MPPGI calibration check).

        Hand calculation:
        Literature Fg for midazolam = 0.57 (CYP3A4 substrate, fm_cyp3a4=1.0).
        Acceptable range: 0.57 * 0.90 = 0.513, 0.57 * 1.10 = 0.627.
        """
        result = _run_oral("midazolam", 5.0)
        fg = result.pk_parameters.fg
        record_property("midazolam_fg", fg)
        assert 0.51 < fg < 0.63, f"midazolam Fg = {fg:.3f}, expected ~0.57"

    def test_fa_high(self, record_property):
        """Midazolam Fa > 0.90 (high Peff, Lennernas 2007)."""
        result = _run_oral("midazolam", 5.0)
        fa = result.pk_parameters.fa
        record_property("midazolam_fa", fa)
        assert fa > 0.90, f"midazolam Fa = {fa:.3f}"

    def test_mass_balance_t0(self):
        """At t=0 all drug mass is in stomach lumen, summing to dose_mg."""
        result = _run_oral("midazolam", 5.0)
        # state_trajectory shape: (n_states, n_time_points)
        total_t0 = result.simulation.state_trajectory[:, 0].sum()
        assert abs(total_t0 - 5.0) < 0.05

    def test_mass_balance_residual(self):
        """Solver-computed mass balance residual < 1% of dose."""
        result = _run_oral("midazolam", 5.0)
        residual_pct = result.simulation.mass_balance_residual / 5.0 * 100
        assert residual_pct < 1.0, f"mass balance residual = {residual_pct:.4f}%"

    def test_tmax_positive(self):
        result = _run_oral("midazolam", 5.0)
        assert result.pk_parameters.tmax > 0

    def test_bioavailability_product(self):
        """bioavailability == Fa * Fg * Fh (exact by construction)."""
        result = _run_oral("midazolam", 5.0)
        pk = result.pk_parameters
        expected = pk.fa * pk.fg * pk.fh
        assert pk.bioavailability == pytest.approx(expected, rel=1e-6)


class TestFelodipineOral:
    """Felodipine PO 10 mg -- independent validation."""

    def test_produces_finite_pk(self):
        result = _run_oral("felodipine", 10.0)
        pk = result.pk_parameters
        for attr in ("cmax", "tmax", "auc_0_inf", "fa", "fg", "fh"):
            val = getattr(pk, attr)
            assert val is not None and math.isfinite(val), f"{attr} = {val}"

    def test_fg_within_2_fold(self, record_property):
        """Fg within 2-fold of literature 0.45 (Edgar 1992).

        Acceptable range: 0.45/2 = 0.225 to 0.45*2 = 0.90.
        """
        result = _run_oral("felodipine", 10.0)
        fg = result.pk_parameters.fg
        record_property("felodipine_fg", fg)
        lit_fg = 0.45
        fold = max(fg / lit_fg, lit_fg / fg)
        assert fold < 2.0, (
            f"felodipine Fg = {fg:.3f}, lit = {lit_fg}, fold = {fold:.2f}"
        )

    def test_fa_high(self, record_property):
        """Felodipine Fa > 0.90 (high Peff)."""
        result = _run_oral("felodipine", 10.0)
        fa = result.pk_parameters.fa
        record_property("felodipine_fa", fa)
        assert fa > 0.90


class TestNifedipineOral:
    """Nifedipine PO 10 mg -- independent validation."""

    def test_produces_finite_pk(self):
        result = _run_oral("nifedipine", 10.0)
        pk = result.pk_parameters
        for attr in ("cmax", "tmax", "auc_0_inf", "fa", "fg", "fh"):
            val = getattr(pk, attr)
            assert val is not None and math.isfinite(val), f"{attr} = {val}"

    def test_fg_within_2_fold(self, record_property):
        """Fg within 2-fold of literature 0.78 (Holtbecker 1996).

        Acceptable range: 0.78/2 = 0.39 to 0.78*2 = 1.56 (capped at 1.0).
        """
        result = _run_oral("nifedipine", 10.0)
        fg = result.pk_parameters.fg
        record_property("nifedipine_fg", fg)
        lit_fg = 0.78
        fold = max(fg / lit_fg, lit_fg / fg)
        assert fold < 2.0, (
            f"nifedipine Fg = {fg:.3f}, lit = {lit_fg}, fold = {fold:.2f}"
        )

    def test_fa_high(self, record_property):
        """Nifedipine Fa > 0.90 (high Peff)."""
        result = _run_oral("nifedipine", 10.0)
        fa = result.pk_parameters.fa
        record_property("nifedipine_fa", fa)
        assert fa > 0.90


class TestRegressionInvariants:
    """Verify existing IV modules still work after oral additions."""

    def test_iv_imports_intact(self):
        from charon.pbpk import simulate_iv, build_rhs, compute_pk_parameters

        assert callable(simulate_iv)
        assert callable(build_rhs)
        assert callable(compute_pk_parameters)
