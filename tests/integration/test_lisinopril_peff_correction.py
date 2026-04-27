"""Sprint 17 Branch B integration tests — lisinopril Peff back-calibration.

Verifies: (1) corrected Peff loaded from YAML, (2) F_oral now matches
literature F_obs=0.25 from Beermann 1988, (3) MRSD shifts toward
literature-consistent value (fold improves from 4.13x to ~3.09x —
close-but-not-quite, expected per Branch B honesty).
"""

from __future__ import annotations

import copy
from pathlib import Path

import yaml

from charon import Pipeline
from charon.core.schema import CompoundConfig, DoseProjectionConfig

REPO_ROOT = Path(__file__).resolve().parents[2]
LISIN_YAML = REPO_ROOT / "validation" / "data" / "tier1_obach" / "compounds" / "lisinopril.yaml"
PANEL_YAML = REPO_ROOT / "validation" / "data" / "fih_reference" / "panel.yaml"

OLD_PEFF = 0.3e-4
NEW_PEFF = 2.10e-5
F_OBS_BEERMANN = 0.25  # Beermann 1988 typical adult value
F_OBS_TOLERANCE = 0.02  # +/- 0.02 around target after calibration


def _load_entry() -> dict:
    panel = yaml.safe_load(PANEL_YAML.read_text())["panel"]
    return next(c for c in panel["compounds"] if c["name"] == "lisinopril")


def _run(compound: CompoundConfig):
    entry = _load_entry()
    pipe = Pipeline(
        compound,
        route=entry["route"],
        dose_mg=1.0,
        dose_projection=DoseProjectionConfig(
            target_ceff_nM=float(entry["target_ceff_nM"]),
            safety_factor=10.0,
            tau_h=24.0,
        ),
    )
    return pipe.run()


def test_lisinopril_peff_corrected_in_yaml():
    """The corrected Peff value is loaded from YAML."""
    data = yaml.safe_load(LISIN_YAML.read_text())
    actual = data["properties"]["permeability"]["peff_cm_s"]["value"]
    assert abs(actual - NEW_PEFF) < 1e-9, f"expected {NEW_PEFF}, got {actual}"


def test_lisinopril_F_oral_matches_beermann_1988():
    """With corrected Peff, predicted F_oral matches Beermann 1988 F_obs=0.25."""
    data = yaml.safe_load(LISIN_YAML.read_text())
    compound = CompoundConfig.model_validate(data)
    result = _run(compound)
    F_pred = result.pk_parameters.bioavailability
    assert abs(F_pred - F_OBS_BEERMANN) < F_OBS_TOLERANCE, (
        f"F_oral_pred {F_pred:.4f} not within {F_OBS_TOLERANCE} of Beermann 1988 "
        f"F_obs={F_OBS_BEERMANN}"
    )


def test_lisinopril_correction_improves_fold():
    """MRSD shifts toward reference; fold improves from ~4.13x toward 3.09x.

    Branch B close-but-not-quite: improvement is real but still outside 3x.
    """
    data = yaml.safe_load(LISIN_YAML.read_text())

    old_data = copy.deepcopy(data)
    old_data["properties"]["permeability"]["peff_cm_s"]["value"] = OLD_PEFF
    old_compound = CompoundConfig.model_validate(old_data)
    new_compound = CompoundConfig.model_validate(data)

    old_result = _run(old_compound)
    new_result = _run(new_compound)

    old_mrsd = old_result.dose_recommendation.mrsd_mg
    new_mrsd = new_result.dose_recommendation.mrsd_mg
    ratio = new_mrsd / old_mrsd

    # New MRSD should be ~30-40% larger (corrects over-prediction of F).
    # Accept band [1.20, 1.55] for ratio = new_mrsd / old_mrsd.
    assert 1.20 <= ratio <= 1.55, (
        f"MRSD ratio {ratio:.3f} (new={new_mrsd:.4g}, old={old_mrsd:.4g}) "
        f"outside expected band [1.20, 1.55]"
    )

    # Fold improves from old (~4.13x) to new (~3.09x), but stays > 3x.
    ref_mg = float(_load_entry()["reference_fih_mg"])
    new_fold = ref_mg / new_mrsd
    assert 2.9 <= new_fold <= 3.3, f"new fold {new_fold:.3f}x outside [2.9, 3.3] band"
