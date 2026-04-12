"""Tests for DoseProjector coordinator."""
from __future__ import annotations

import pytest


def _make_pk():
    from charon.core.schema import PKParameters
    return PKParameters(
        cmax=0.02, tmax=0.87, auc_0_inf=0.27, half_life=1.7,
        cl_apparent=18.6, vss=None, bioavailability=0.48,
        fa=0.97, fg=0.575, fh=0.86,
    )


def _make_compound():
    from charon.core.schema import CompoundConfig
    return CompoundConfig(
        name="test_drug", smiles="C", molecular_weight=325.77,
        properties={"binding": {"fu_p": {"value": 0.03, "source": "experimental"}}},
    )


class TestProjectFIHDose:
    def test_hed_only(self):
        from charon.translational.dose_projector import project_fih_dose
        from charon.core.schema import DoseProjectionConfig
        config = DoseProjectionConfig(noael_mg_kg=50.0, noael_species="rat")
        rec = project_fih_dose(pk=_make_pk(), compound=_make_compound(), config=config, route="oral")
        assert rec.hed is not None
        assert rec.mabel is None
        assert rec.pad is None
        assert rec.mrsd_mg == pytest.approx(rec.hed.mrsd_mg, rel=1e-6)
        assert rec.limiting_method == "hed"

    def test_mabel_only(self):
        from charon.translational.dose_projector import project_fih_dose
        from charon.core.schema import DoseProjectionConfig
        config = DoseProjectionConfig(target_kd_nM=10.0)
        rec = project_fih_dose(pk=_make_pk(), compound=_make_compound(), config=config, route="oral")
        assert rec.hed is None
        assert rec.mabel is not None
        assert rec.mrsd_mg == pytest.approx(rec.mabel.mrsd_mg, rel=1e-6)

    def test_pad_only(self):
        from charon.translational.dose_projector import project_fih_dose
        from charon.core.schema import DoseProjectionConfig
        config = DoseProjectionConfig(target_ceff_nM=100.0)
        rec = project_fih_dose(pk=_make_pk(), compound=_make_compound(), config=config, route="oral")
        assert rec.hed is None
        assert rec.pad is not None
        assert rec.mrsd_mg == pytest.approx(rec.pad.mrsd_mg, rel=1e-6)

    def test_all_three_min_wins(self):
        from charon.translational.dose_projector import project_fih_dose
        from charon.core.schema import DoseProjectionConfig
        config = DoseProjectionConfig(
            noael_mg_kg=50.0, noael_species="rat",
            target_kd_nM=10.0, target_ceff_nM=100.0,
        )
        rec = project_fih_dose(pk=_make_pk(), compound=_make_compound(), config=config, route="oral")
        assert rec.hed is not None and rec.mabel is not None and rec.pad is not None
        all_mrsd = [rec.hed.mrsd_mg, rec.mabel.mrsd_mg, rec.pad.mrsd_mg]
        assert rec.mrsd_mg == pytest.approx(min(all_mrsd), rel=1e-6)

    def test_no_inputs_raises(self):
        from charon.translational.dose_projector import project_fih_dose
        from charon.core.schema import DoseProjectionConfig
        with pytest.raises(ValueError, match="at least one"):
            project_fih_dose(
                pk=_make_pk(), compound=_make_compound(),
                config=DoseProjectionConfig(), route="oral",
            )

    def test_salt_correction(self):
        from charon.translational.dose_projector import project_fih_dose
        from charon.core.schema import DoseProjectionConfig, CompoundConfig, SaltForm
        config = DoseProjectionConfig(noael_mg_kg=50.0, noael_species="rat")
        compound = CompoundConfig(
            name="test_salt", smiles="C", molecular_weight=300.0,
            salt_form=SaltForm(name="HCl", mw_salt=336.0, salt_factor=300.0 / 336.0),
            properties={"binding": {"fu_p": {"value": 0.5, "source": "experimental"}}},
        )
        rec = project_fih_dose(pk=_make_pk(), compound=compound, config=config, route="oral")
        assert rec.salt_factor == pytest.approx(300.0 / 336.0, rel=1e-6)
        assert rec.salt_factor < 1.0

    def test_rationale_contains_methods(self):
        from charon.translational.dose_projector import project_fih_dose
        from charon.core.schema import DoseProjectionConfig
        config = DoseProjectionConfig(noael_mg_kg=50.0, noael_species="rat", target_kd_nM=10.0)
        rec = project_fih_dose(pk=_make_pk(), compound=_make_compound(), config=config, route="oral")
        assert "HED" in rec.rationale
        assert "MABEL" in rec.rationale

    def test_custom_safety_factor(self):
        from charon.translational.dose_projector import project_fih_dose
        from charon.core.schema import DoseProjectionConfig
        r3 = project_fih_dose(
            pk=_make_pk(), compound=_make_compound(),
            config=DoseProjectionConfig(
                noael_mg_kg=50.0, noael_species="rat", safety_factor=3.0
            ),
            route="oral",
        )
        r10 = project_fih_dose(
            pk=_make_pk(), compound=_make_compound(),
            config=DoseProjectionConfig(
                noael_mg_kg=50.0, noael_species="rat", safety_factor=10.0
            ),
            route="oral",
        )
        assert r3.mrsd_mg > r10.mrsd_mg
