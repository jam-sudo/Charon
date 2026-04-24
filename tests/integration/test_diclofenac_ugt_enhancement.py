"""Sprint 13 integration test: diclofenac MRSD with/without UGT enhancement.

Mirrors the Sprint 12 atorvastatin test structure. Validates:
1. Enhanced MRSD is strictly larger than baseline (multiplier > 1 scales CLh
   upward, which scales MRSD upward via the PAD formula).
2. Ratio MRSD_enhanced / MRSD_base is within [2, 5] — broader than
   atorvastatin's range because diclofenac baseline is low-extraction
   (fu_p=0.005, CLint=11 uL/min/mg) so scaling is more linear but smaller.
3. clint_liver_L_h ratio exactly matches the multiplier value.
"""

from __future__ import annotations

import copy
from pathlib import Path

import pytest
import yaml

from charon import Pipeline
from charon.core.schema import CompoundConfig, DoseProjectionConfig

REPO_ROOT = Path(__file__).resolve().parents[2]
DICLO_YAML = REPO_ROOT / "validation" / "data" / "tier1_obach" / "compounds" / "diclofenac.yaml"
PANEL_YAML = REPO_ROOT / "validation" / "data" / "fih_reference" / "panel.yaml"


def _load_diclofenac_entry() -> dict:
    panel = yaml.safe_load(PANEL_YAML.read_text())["panel"]
    return next(
        c for c in panel["compounds"]
        if c["name"] == "diclofenac" and c["tier"] == "gold"
    )


def _run_pipeline(compound: CompoundConfig, entry: dict) -> tuple[float, dict]:
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
    result = pipe.run()
    assert result.dose_recommendation is not None
    return float(result.dose_recommendation.mrsd_mg), result.metadata


def test_diclofenac_mrsd_with_multiplier_larger_than_without():
    """With clint_multiplier > 1, diclofenac MRSD must be strictly larger."""
    data = yaml.safe_load(DICLO_YAML.read_text())
    entry = _load_diclofenac_entry()

    compound_enhanced = CompoundConfig.model_validate(data)
    assert compound_enhanced.properties.metabolism.hepatic_clint_multiplier is not None, (
        "diclofenac.yaml must have hepatic_clint_multiplier populated (Task 1 prerequisite)"
    )
    mrsd_enhanced, _ = _run_pipeline(compound_enhanced, entry)

    data_base = copy.deepcopy(data)
    data_base["properties"]["metabolism"].pop("hepatic_clint_multiplier", None)
    compound_base = CompoundConfig.model_validate(data_base)
    assert compound_base.properties.metabolism.hepatic_clint_multiplier is None
    mrsd_base, _ = _run_pipeline(compound_base, entry)

    assert mrsd_enhanced > mrsd_base, (
        f"Enhanced MRSD ({mrsd_enhanced:.3g}) must exceed baseline ({mrsd_base:.3g})"
    )


def test_diclofenac_mrsd_ratio_within_literature_range():
    """MRSD_enhanced / MRSD_base should be within [2, 5] — near the multiplier
    value (3.5) in the near-linear regime for a low-extraction compound like
    diclofenac (fu_p=0.005, so low fu_b*CLint even after enhancement)."""
    data = yaml.safe_load(DICLO_YAML.read_text())
    entry = _load_diclofenac_entry()

    compound_enhanced = CompoundConfig.model_validate(data)
    mrsd_enhanced, _ = _run_pipeline(compound_enhanced, entry)

    data_base = copy.deepcopy(data)
    data_base["properties"]["metabolism"].pop("hepatic_clint_multiplier", None)
    compound_base = CompoundConfig.model_validate(data_base)
    mrsd_base, _ = _run_pipeline(compound_base, entry)

    ratio = mrsd_enhanced / mrsd_base
    assert 2.0 < ratio < 5.0, (
        f"MRSD ratio {ratio:.2f} outside plausible [2, 5] range; "
        f"check multiplier logic or extraction saturation"
    )


def test_diclofenac_enhanced_clint_liver_metadata_scales_by_multiplier():
    """clint_liver_L_h ratio (enhanced/base) should exactly match the multiplier.
    This is a pure linear scaling — ConversionStep multiplies CLint_liver before
    liver model, so the reported clint_liver_L_h in metadata is the enhanced
    value."""
    data = yaml.safe_load(DICLO_YAML.read_text())
    entry = _load_diclofenac_entry()

    compound_enhanced = CompoundConfig.model_validate(data)
    multiplier = compound_enhanced.properties.metabolism.hepatic_clint_multiplier.value
    _, md_enhanced = _run_pipeline(compound_enhanced, entry)

    data_base = copy.deepcopy(data)
    data_base["properties"]["metabolism"].pop("hepatic_clint_multiplier", None)
    compound_base = CompoundConfig.model_validate(data_base)
    _, md_base = _run_pipeline(compound_base, entry)

    clint_enhanced = float(md_enhanced["clint_liver_L_h"])
    clint_base = float(md_base["clint_liver_L_h"])
    ratio = clint_enhanced / clint_base
    assert abs(ratio - multiplier) < 1e-3, (
        f"clint_liver_L_h ratio {ratio:.6f} != multiplier {multiplier}"
    )
