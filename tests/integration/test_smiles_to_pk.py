"""Layer 0-1-2 integration: SMILES/compound → Cp-time → PK parameters.

Sprint 3 IV-only subset.  Oral integration tests will land in Sprint 3b.

Primary compound: theophylline (well-behaved R&R Kp, validates the
entire IV PBPK kernel with CL / Vss / t_half all within 2-fold of
observed human values).
"""

import numpy as np
import pytest

from charon import Pipeline, PipelineResult
from charon.core.schema import (
    BindingProperties,
    CompoundConfig,
    CompoundProperties,
    MetabolismProperties,
    PhysicochemicalProperties,
    PredictedProperty,
    RenalProperties,
)


def _p(value: float, unit: str = "") -> PredictedProperty:
    return PredictedProperty(value=float(value), source="experimental", unit=unit)


@pytest.fixture
def theophylline() -> CompoundConfig:
    return CompoundConfig(
        name="theophylline",
        smiles="Cn1c(=O)c2[nH]cnc2n(C)c1=O",
        molecular_weight=180.17,
        source="experimental",
        properties=CompoundProperties(
            physicochemical=PhysicochemicalProperties(
                logp=_p(-0.02),
            ),
            binding=BindingProperties(
                fu_p=_p(0.60, "fraction"),
                fu_inc=_p(1.0, "fraction"),
                bp_ratio=_p(0.85, "ratio"),
            ),
            metabolism=MetabolismProperties(
                clint_uL_min_mg=_p(1.8, "uL/min/mg"),
            ),
            renal=RenalProperties(
                clrenal_L_h=_p(0.1, "L/h"),
            ),
        ),
    )


class TestTheophyllineIvBolusE2E:
    def test_runs_and_produces_pk(self, theophylline):
        pipe = Pipeline(
            compound=theophylline,
            route="iv_bolus",
            dose_mg=100.0,
            duration_h=168.0,
        )
        result: PipelineResult = pipe.run()
        assert result.pk_parameters.cmax is not None and result.pk_parameters.cmax > 0
        assert result.pk_parameters.cl_apparent is not None and result.pk_parameters.cl_apparent > 0
        assert result.pk_parameters.vss is not None and result.pk_parameters.vss > 0
        assert result.pk_parameters.half_life is not None and result.pk_parameters.half_life > 0

    def test_cl_vss_t_half_within_2_fold(self, theophylline):
        """All three core PK metrics within the Sprint 3 2-fold target."""
        pipe = Pipeline(
            compound=theophylline,
            route="iv_bolus",
            dose_mg=100.0,
            duration_h=168.0,
        )
        result = pipe.run()
        obs = {"cl": 2.9, "vss": 35.0, "t_half": 8.0}
        pred = {
            "cl": result.pk_parameters.cl_apparent,
            "vss": result.pk_parameters.vss,
            "t_half": result.pk_parameters.half_life,
        }
        folds = {k: max(pred[k] / obs[k], obs[k] / pred[k]) for k in obs}
        for metric, fe in folds.items():
            assert fe < 2.0, (
                f"{metric}: fold error {fe:.2f} "
                f"(pred={pred[metric]:.3f}, obs={obs[metric]:.3f})"
            )

    def test_metadata_includes_pbpk_specifics(self, theophylline):
        pipe = Pipeline(
            compound=theophylline,
            route="iv_bolus",
            dose_mg=100.0,
            duration_h=168.0,
        )
        result = pipe.run()
        meta = result.metadata
        assert meta["solver_method"] == "BDF"
        assert meta["compound_type"] == "neutral"
        assert meta["fu_b"] == pytest.approx(0.60 / 0.85, rel=1e-6)
        # CLint_liver = 1.8/1.0 * 40 * 1500 / 1e6 * 60 = 6.48 L/h
        assert meta["clint_liver_L_h"] == pytest.approx(6.48, rel=1e-3)

    def test_iv_infusion_runs(self, theophylline):
        """IV infusion over 30 minutes still produces sane PK metrics."""
        pipe = Pipeline(
            compound=theophylline,
            route="iv_infusion",
            dose_mg=100.0,
            duration_h=168.0,
            infusion_duration_h=0.5,
        )
        result = pipe.run()
        pk = result.pk_parameters
        assert pk.cl_apparent is not None
        # Benet correction: Vss should still come out close to 35 L
        assert 15.0 < pk.vss < 80.0
