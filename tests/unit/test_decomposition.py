"""Unit tests for charon.translational.decomposition.

All factors are signed (can be < 1 or > 1). Log-additive invariant:

    log10(fold_observed_signed) = log10(fold_liver_model_signed)
                                + log10(fold_route_bias)
                                + log10(fold_residual_signed)
"""

from __future__ import annotations

import math

import pytest

from charon.translational.decomposition import (
    DecompositionResult,
    compute_route_bias_factor,
    decompose_fold_error,
    select_best_alternate_liver_model,
    to_symmetric,
)


# ---------------------------------------------------------------------------
# to_symmetric helper
# ---------------------------------------------------------------------------

def test_to_symmetric_above_one():
    assert to_symmetric(2.0) == 2.0


def test_to_symmetric_below_one():
    assert to_symmetric(0.5) == pytest.approx(2.0, rel=1e-12)


def test_to_symmetric_exactly_one():
    assert to_symmetric(1.0) == 1.0


# ---------------------------------------------------------------------------
# select_best_alternate_liver_model
# ---------------------------------------------------------------------------

def test_select_best_alternate_picks_closest_to_reference():
    # ws=10, pt=5, disp=15; ref=5 -> pt is closest (|log10(5/5)| = 0)
    name, mrsd = select_best_alternate_liver_model(
        ws=10.0, pt=5.0, disp=15.0, reference=5.0
    )
    assert name == "parallel_tube"
    assert mrsd == 5.0


def test_select_best_alternate_picks_dispersion_when_closer():
    # ws=2, pt=1, disp=4; ref=5 -> |log10(disp/ref)|≈0.097, |log10(pt/ref)|≈0.699
    name, mrsd = select_best_alternate_liver_model(
        ws=2.0, pt=1.0, disp=4.0, reference=5.0
    )
    assert name == "dispersion"
    assert mrsd == 4.0


def test_select_best_alternate_falls_back_to_well_stirred_when_no_improvement():
    # ws=5 (exact match), pt=100, disp=0.5; alternates both worse than ws
    name, mrsd = select_best_alternate_liver_model(
        ws=5.0, pt=100.0, disp=0.5, reference=5.0
    )
    assert name == "well_stirred"
    assert mrsd == 5.0


def test_select_best_alternate_tie_break_prefers_parallel_tube():
    # pt=3, disp=3, ref=5 — both equidistant at |log10(3/5)|=0.222
    # Tie-break must prefer parallel_tube (stable sort, listed first).
    # ws=4 is still closer (|log10(4/5)|=0.097) so that wins overall;
    # construct case where pt==disp and both beat ws.
    name, mrsd = select_best_alternate_liver_model(
        ws=100.0, pt=3.0, disp=3.0, reference=5.0
    )
    assert name == "parallel_tube"
    assert mrsd == 3.0


# ---------------------------------------------------------------------------
# compute_route_bias_factor
# ---------------------------------------------------------------------------

def test_route_bias_iv_reference_returns_unity():
    # IV reference: no 1/F bias exists, factor=1.0, no flags
    factor, flags = compute_route_bias_factor(route_ref="iv", f_lit=None)
    assert factor == 1.0
    assert flags == ()


def test_route_bias_iv_reference_ignores_f_lit():
    # Even if f_lit supplied, IV route means factor=1.0
    factor, flags = compute_route_bias_factor(route_ref="iv", f_lit=0.3)
    assert factor == 1.0
    assert flags == ()


def test_route_bias_oral_with_known_f():
    # propranolol F=0.26 -> factor = 1/0.26 = 3.846...
    factor, flags = compute_route_bias_factor(route_ref="oral", f_lit=0.26)
    assert factor == pytest.approx(1.0 / 0.26, rel=1e-9)
    assert flags == ()


def test_route_bias_oral_with_unknown_f_flags():
    # Oral reference but no F literature -> factor=1.0 + flag
    factor, flags = compute_route_bias_factor(route_ref="oral", f_lit=None)
    assert factor == 1.0
    assert flags == ("f_unknown",)


# ---------------------------------------------------------------------------
# decompose_fold_error (end-to-end, SIGNED invariant)
# ---------------------------------------------------------------------------

def test_additivity_invariant_synthetic():
    """Signed-factor log-additivity.

    Hand-calc: ws=10, pt=5, disp=15, ref=5, F=0.5 (oral)
      fold_observed_signed = 10 / 5 = 2.0
      best alt: pt (exact match at 5.0). fold_liver_model_signed = 10 / 5 = 2.0
      fold_route_bias = 1 / 0.5 = 2.0
      fold_residual_signed = 2.0 / (2.0 * 2.0) = 0.5
      Invariant: log10(2.0) == log10(2.0) + log10(2.0) + log10(0.5)
                 0.30103     == 0.30103 + 0.30103 + (-0.30103) = 0.30103 ✓
    """
    result = decompose_fold_error(
        mrsd_ws=10.0,
        mrsd_pt=5.0,
        mrsd_disp=15.0,
        f_lit=0.5,
        route_ref="oral",
        fih_reference_mg=5.0,
    )
    log_lhs = math.log10(result.fold_observed_signed)
    log_rhs = (
        math.log10(result.fold_liver_model_signed)
        + math.log10(result.fold_route_bias)
        + math.log10(result.fold_residual_signed)
    )
    assert log_lhs == pytest.approx(log_rhs, rel=1e-9, abs=1e-9)


def test_decompose_concrete_signed_values():
    """Exact numeric check on every signed factor for the canonical case."""
    result = decompose_fold_error(
        mrsd_ws=10.0,
        mrsd_pt=5.0,
        mrsd_disp=15.0,
        f_lit=0.5,
        route_ref="oral",
        fih_reference_mg=5.0,
    )
    assert result.fold_observed_signed == pytest.approx(2.0, rel=1e-12)
    assert result.fold_liver_model_signed == pytest.approx(2.0, rel=1e-12)
    assert result.fold_route_bias == pytest.approx(2.0, rel=1e-12)
    assert result.fold_residual_signed == pytest.approx(0.5, rel=1e-12)
    assert result.best_alt_model_name == "parallel_tube"
    assert result.flags == ()


def test_decompose_iv_reference_no_route_factor():
    """IV reference -> route bias factor must be exactly 1.0, zero log contribution."""
    result = decompose_fold_error(
        mrsd_ws=2.0,
        mrsd_pt=1.0,
        mrsd_disp=3.0,
        f_lit=None,
        route_ref="iv",
        fih_reference_mg=1.0,
    )
    assert result.fold_route_bias == 1.0
    assert "f_unknown" not in result.flags


def test_decompose_well_stirred_is_best_returns_unity_liver_factor():
    """When well_stirred is closest to reference, liver_model factor = 1.0 exactly."""
    # ws=5.0 (exact), pt=100, disp=0.5 -> alternates both worse
    result = decompose_fold_error(
        mrsd_ws=5.0,
        mrsd_pt=100.0,
        mrsd_disp=0.5,
        f_lit=None,
        route_ref="iv",
        fih_reference_mg=5.0,
    )
    assert result.fold_liver_model_signed == 1.0
    assert result.best_alt_model_name == "well_stirred"


def test_decompose_symmetric_report_from_signed():
    """to_symmetric applied to signed fields yields the benchmark-style fold."""
    result = decompose_fold_error(
        mrsd_ws=10.0,
        mrsd_pt=5.0,
        mrsd_disp=15.0,
        f_lit=0.5,
        route_ref="oral",
        fih_reference_mg=5.0,
    )
    # signed 2.0 -> symmetric 2.0; signed 0.5 -> symmetric 2.0
    assert to_symmetric(result.fold_observed_signed) == pytest.approx(2.0)
    assert to_symmetric(result.fold_residual_signed) == pytest.approx(2.0)
