"""Unit tests for LHS sampling with Iman-Conover correlation.

Tests cover:
  1. Correct sample count
  2. Parameter name preservation
  3. LHS uniform coverage (all samples distinct)
  4. Normal distribution centering
  5. Lognormal positivity
  6. fu_p physical clipping to [0.001, 1.0]
  7. Seed reproducibility
  8. Iman-Conover correlation induction
  9. No-correlation flag
  10. build_param_specs with experimental values (fallback CVs)
  11. build_param_specs with CI-derived sigma
"""

import numpy as np
import pytest

from charon.core.schema import (
    BindingProperties,
    CompoundConfig,
    CompoundProperties,
    MetabolismProperties,
    PermeabilityProperties,
    PhysicochemicalProperties,
    PredictedProperty,
)
from charon.uncertainty.sampling import (
    SamplingResult,
    build_param_specs,
    generate_lhs_samples,
)


def _predicted(value: float, unit: str = "", source: str = "experimental") -> PredictedProperty:
    return PredictedProperty(value=float(value), source=source, unit=unit)


def _predicted_with_ci(
    value: float, lower: float, upper: float, unit: str = "", source: str = "ml_ensemble"
) -> PredictedProperty:
    return PredictedProperty(
        value=float(value),
        ci_90_lower=float(lower),
        ci_90_upper=float(upper),
        source=source,
        unit=unit,
    )


@pytest.fixture
def midazolam_experimental():
    """Midazolam with experimental values (no CIs) -> fallback CVs."""
    return CompoundConfig(
        name="midazolam",
        smiles="Cc1ncc2n1-c1ccc(Cl)cc1C(c1ccccc1F)=NC2",
        molecular_weight=325.77,
        source="experimental",
        properties=CompoundProperties(
            physicochemical=PhysicochemicalProperties(
                logp=_predicted(3.89),
            ),
            binding=BindingProperties(
                fu_p=_predicted(0.032),
                bp_ratio=_predicted(0.67),
            ),
            metabolism=MetabolismProperties(
                clint_uL_min_mg=_predicted(93.0),
            ),
            permeability=PermeabilityProperties(
                peff_cm_s=_predicted(3.5e-4),
            ),
        ),
    )


@pytest.fixture
def compound_with_ci():
    """Compound with ML-predicted properties including 90% CIs."""
    return CompoundConfig(
        name="test_compound",
        smiles="c1ccccc1",
        molecular_weight=78.11,
        source="predicted",
        properties=CompoundProperties(
            physicochemical=PhysicochemicalProperties(
                logp=_predicted_with_ci(2.0, 1.5, 2.5),
            ),
            binding=BindingProperties(
                fu_p=_predicted_with_ci(0.10, 0.03, 0.33),
                bp_ratio=_predicted_with_ci(1.0, 0.85, 1.15),
            ),
            metabolism=MetabolismProperties(
                clint_uL_min_mg=_predicted_with_ci(30.0, 5.0, 180.0),
            ),
            permeability=PermeabilityProperties(
                peff_cm_s=_predicted_with_ci(1.0e-4, 2.0e-5, 5.0e-4),
            ),
        ),
    )


# ---- Test 1: Correct sample count ------------------------------------------

def test_returns_correct_count(midazolam_experimental):
    specs = build_param_specs(midazolam_experimental)
    result = generate_lhs_samples(
        param_specs=specs, n_samples=50, correlation="none", seed=42,
    )
    assert len(result.samples) == 50
    assert isinstance(result, SamplingResult)


# ---- Test 2: Parameter names preserved --------------------------------------

def test_param_names_preserved(midazolam_experimental):
    specs = build_param_specs(midazolam_experimental)
    result = generate_lhs_samples(
        param_specs=specs, n_samples=10, correlation="none", seed=42,
    )
    assert set(result.param_names) == set(specs.keys())
    for sample in result.samples:
        assert set(sample.keys()) == set(specs.keys())


# ---- Test 3: LHS uniform coverage ------------------------------------------

def test_lhs_uniform_coverage(midazolam_experimental):
    """100 LHS samples should all be distinct (no duplicates)."""
    specs = build_param_specs(midazolam_experimental)
    result = generate_lhs_samples(
        param_specs=specs, n_samples=100, correlation="none", seed=42,
    )
    # Check that all samples are distinct for the first parameter
    first_param = result.param_names[0]
    values = [s[first_param] for s in result.samples]
    assert len(set(values)) == 100


# ---- Test 4: Normal distribution mean check ---------------------------------

def test_normal_distribution_range(midazolam_experimental):
    """Mean of 500 normal samples should be near mu.

    logP (3.89) is sampled as normal. With CV=0.15,
    sigma = 3.89 * 0.15 = 0.58. Mean of 500 should be within ~0.15 of 3.89.
    """
    specs = build_param_specs(midazolam_experimental)
    result = generate_lhs_samples(
        param_specs=specs, n_samples=500, correlation="none", seed=42,
    )
    logp_values = [s["logp"] for s in result.samples]
    mean_logp = np.mean(logp_values)
    assert abs(mean_logp - 3.89) < 0.15, f"logP mean {mean_logp} too far from 3.89"


# ---- Test 5: Lognormal positivity -------------------------------------------

def test_lognormal_positive(midazolam_experimental):
    """All lognormal-sampled parameters must be > 0."""
    specs = build_param_specs(midazolam_experimental)
    result = generate_lhs_samples(
        param_specs=specs, n_samples=200, correlation="none", seed=42,
    )
    for param in ["fu_p", "clint_uL_min_mg", "peff_cm_s"]:
        if param in result.param_names:
            values = [s[param] for s in result.samples]
            assert all(v > 0 for v in values), f"{param} has non-positive values"


# ---- Test 6: fu_p clipped to 1.0 -------------------------------------------

def test_fu_p_clipped_to_1():
    """fu_p with very wide uncertainty should never exceed 1.0.

    Hand calculation: fu_p=0.5, lognormal CV=0.50 -> sigma_log ~= 0.213
    In 1000 samples some would exceed 1.0 without clipping.
    We use a deliberately large sigma to provoke the clip.
    """
    specs = {
        "fu_p": (0.5, 0.50, "lognormal"),  # mu=0.5, sigma_log=0.50
    }
    result = generate_lhs_samples(
        param_specs=specs, n_samples=1000, correlation="none", seed=42,
    )
    fu_p_values = [s["fu_p"] for s in result.samples]
    assert all(v <= 1.0 for v in fu_p_values), "fu_p exceeds 1.0"
    assert all(v >= 0.001 for v in fu_p_values), "fu_p below 0.001"


# ---- Test 7: Seed reproducibility ------------------------------------------

def test_seed_reproducibility(midazolam_experimental):
    specs = build_param_specs(midazolam_experimental)
    r1 = generate_lhs_samples(param_specs=specs, n_samples=50, correlation="none", seed=123)
    r2 = generate_lhs_samples(param_specs=specs, n_samples=50, correlation="none", seed=123)
    for s1, s2 in zip(r1.samples, r2.samples):
        for key in s1:
            assert s1[key] == s2[key], f"Mismatch for {key}"


# ---- Test 8: Iman-Conover induces correlation -------------------------------

def test_iman_conover_induces_correlation(midazolam_experimental):
    """Spearman rho(logP, fu_p) should be approximately -0.6 +/- 0.15.

    With N=2000 and target correlation -0.6, the Iman-Conover algorithm
    should produce Spearman rank correlation close to the target.
    """
    from scipy.stats import spearmanr

    specs = build_param_specs(midazolam_experimental)
    result = generate_lhs_samples(
        param_specs=specs, n_samples=2000, correlation="iman_conover", seed=42,
    )
    logp_vals = np.array([s["logp"] for s in result.samples])
    fup_vals = np.array([s["fu_p"] for s in result.samples])
    rho, _ = spearmanr(logp_vals, fup_vals)
    assert abs(rho - (-0.6)) < 0.15, f"Spearman rho = {rho}, expected ~-0.6"
    assert result.correlation_applied is True


# ---- Test 9: No-correlation flag --------------------------------------------

def test_no_correlation_flag(midazolam_experimental):
    specs = build_param_specs(midazolam_experimental)
    result = generate_lhs_samples(
        param_specs=specs, n_samples=50, correlation="none", seed=42,
    )
    assert result.correlation_applied is False


# ---- Test 10: build_param_specs with experimental (fallback CVs) ------------

def test_build_param_specs_experimental(midazolam_experimental):
    """Experimental values without CIs should use fallback CVs.

    Hand-checked fallback CVs:
      logP: CV=0.15 -> sigma = 3.89 * 0.15 = 0.5835, distribution=normal
      fu_p: CV=0.50 -> sigma_log = 0.50 * log10(e) ≈ 0.217, distribution=lognormal
      CLint: CV=0.60, distribution=lognormal
      Peff: CV=0.40, distribution=lognormal
      BP: CV=0.10, distribution=normal
      MPPGL: always Normal(40, 8)
    """
    specs = build_param_specs(midazolam_experimental)

    # logP: normal, sigma = mu * CV = 3.89 * 0.15
    assert specs["logp"][2] == "normal"
    assert abs(specs["logp"][0] - 3.89) < 1e-6
    expected_sigma_logp = 3.89 * 0.15
    assert abs(specs["logp"][1] - expected_sigma_logp) < 0.01

    # fu_p: lognormal
    assert specs["fu_p"][2] == "lognormal"
    assert abs(specs["fu_p"][0] - 0.032) < 1e-6

    # CLint: lognormal
    assert specs["clint_uL_min_mg"][2] == "lognormal"

    # MPPGL: always included, normal
    assert "mppgl" in specs
    assert specs["mppgl"][2] == "normal"
    assert abs(specs["mppgl"][0] - 40.0) < 1e-6
    assert abs(specs["mppgl"][1] - 8.0) < 1e-6

    # BP: normal
    assert specs["bp_ratio"][2] == "normal"


# ---- Test 11: build_param_specs with CI-derived sigma -----------------------

def test_build_param_specs_with_ci(compound_with_ci):
    """Compounds with 90% CIs should derive sigma from CI width.

    Hand calculation for logP:
      CI = [1.5, 2.5], width = 1.0
      sigma = width / 3.29 = 1.0 / 3.29 ≈ 0.3040
      distribution = normal

    Hand calculation for fu_p (lognormal):
      CI = [0.03, 0.33]
      sigma_log = (log10(0.33) - log10(0.03)) / 3.29
              = (-0.4815 - (-1.5229)) / 3.29
              = 1.0414 / 3.29
              ≈ 0.3165
    """
    specs = build_param_specs(compound_with_ci)

    # logP: normal, CI-derived sigma
    assert specs["logp"][2] == "normal"
    expected_sigma_logp = (2.5 - 1.5) / 3.29
    assert abs(specs["logp"][1] - expected_sigma_logp) < 0.01

    # fu_p: lognormal, CI-derived sigma_log
    assert specs["fu_p"][2] == "lognormal"
    expected_sigma_fup = (np.log10(0.33) - np.log10(0.03)) / 3.29
    assert abs(specs["fu_p"][1] - expected_sigma_fup) < 0.01

    # CLint: lognormal, CI-derived sigma_log
    assert specs["clint_uL_min_mg"][2] == "lognormal"
    expected_sigma_clint = (np.log10(180.0) - np.log10(5.0)) / 3.29
    assert abs(specs["clint_uL_min_mg"][1] - expected_sigma_clint) < 0.01
