"""Sprint 12 integration test: atorvastatin MRSD with/without OATP enhancement.

Validates:
1. The enhanced MRSD is strictly larger than the non-enhanced MRSD (the
   multiplier scales CLh up, which scales MRSD up via the PAD formula).
2. The ratio MRSD_enhanced / MRSD_base is within a physiologically
   plausible range (4-12x, matching the literature multiplier range).
3. Pipeline metadata shows the enhanced clint_liver_L_h.
"""

from __future__ import annotations

import copy
from pathlib import Path

import pytest
import yaml

from charon import Pipeline
from charon.core.schema import CompoundConfig, DoseProjectionConfig

REPO_ROOT = Path(__file__).resolve().parents[2]
ATORVA_YAML = REPO_ROOT / "validation" / "data" / "tier1_obach" / "compounds" / "atorvastatin.yaml"
PANEL_YAML = REPO_ROOT / "validation" / "data" / "fih_reference" / "panel.yaml"


def _load_atorvastatin_entry() -> dict:
    panel = yaml.safe_load(PANEL_YAML.read_text())["panel"]
    return next(
        c for c in panel["compounds"]
        if c["name"] == "atorvastatin" and c["tier"] == "gold"
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


def test_atorvastatin_mrsd_with_multiplier_larger_than_without():
    """With clint_multiplier=8.0, atorvastatin MRSD must be strictly larger."""
    data = yaml.safe_load(ATORVA_YAML.read_text())
    entry = _load_atorvastatin_entry()

    # Load WITH multiplier (from committed YAML)
    compound_enhanced = CompoundConfig.model_validate(data)
    assert compound_enhanced.properties.metabolism.hepatic_clint_multiplier is not None, (
        "atorvastatin.yaml must have hepatic_clint_multiplier populated (Task 4 prerequisite)"
    )
    mrsd_enhanced, _ = _run_pipeline(compound_enhanced, entry)

    # Load WITHOUT multiplier (strip from in-memory dict)
    data_base = copy.deepcopy(data)
    data_base["properties"]["metabolism"].pop("hepatic_clint_multiplier", None)
    compound_base = CompoundConfig.model_validate(data_base)
    assert compound_base.properties.metabolism.hepatic_clint_multiplier is None
    mrsd_base, _ = _run_pipeline(compound_base, entry)

    assert mrsd_enhanced > mrsd_base, (
        f"Enhanced MRSD ({mrsd_enhanced:.3g}) must exceed baseline ({mrsd_base:.3g})"
    )


def test_atorvastatin_mrsd_ratio_within_literature_range():
    """MRSD_enhanced / MRSD_base should be within 4-12x (plausibly near the
    multiplier value of 8, allowing for well-stirred curvature)."""
    data = yaml.safe_load(ATORVA_YAML.read_text())
    entry = _load_atorvastatin_entry()

    compound_enhanced = CompoundConfig.model_validate(data)
    mrsd_enhanced, _ = _run_pipeline(compound_enhanced, entry)

    data_base = copy.deepcopy(data)
    data_base["properties"]["metabolism"].pop("hepatic_clint_multiplier", None)
    compound_base = CompoundConfig.model_validate(data_base)
    mrsd_base, _ = _run_pipeline(compound_base, entry)

    ratio = mrsd_enhanced / mrsd_base
    assert 4.0 < ratio < 12.0, (
        f"MRSD ratio {ratio:.2f} outside plausible [4, 12] range; "
        f"check multiplier logic or well-stirred curvature"
    )


def test_atorvastatin_enhanced_clint_liver_metadata_scales():
    """The Pipeline metadata exposes clint_liver_L_h. Enhanced should be
    approximately 8x the non-enhanced value."""
    data = yaml.safe_load(ATORVA_YAML.read_text())
    entry = _load_atorvastatin_entry()

    compound_enhanced = CompoundConfig.model_validate(data)
    _, md_enhanced = _run_pipeline(compound_enhanced, entry)

    data_base = copy.deepcopy(data)
    data_base["properties"]["metabolism"].pop("hepatic_clint_multiplier", None)
    compound_base = CompoundConfig.model_validate(data_base)
    _, md_base = _run_pipeline(compound_base, entry)

    clint_enhanced = float(md_enhanced["clint_liver_L_h"])
    clint_base = float(md_base["clint_liver_L_h"])
    ratio = clint_enhanced / clint_base
    # Should be 8.0x within numerical tolerance
    assert 7.8 < ratio < 8.2, (
        f"clint_liver_L_h ratio {ratio:.3f} outside [7.8, 8.2]; expected ~8x"
    )
