"""Sprint 15 integration test: propranolol MRSD with/without CYP2D6 enhancement.

Mirrors Sprint 12 (atorvastatin) and Sprint 13 (diclofenac) test structures.
Validates:
1. Enhanced MRSD is strictly larger than baseline.
2. Ratio MRSD_enhanced / MRSD_base is within [4, 25] — broader than Sprint 12
   atorvastatin's [4, 12] or Sprint 13 diclofenac's [2, 5] because:
   - propranolol's literature-supported multiplier range (CYP2D6 high-extraction
     base IVIVE bias) is larger and more uncertain than OATP1B1 or UGT2B7
   - both CL and F change with multiplier (well-stirred saturation), so the
     MRSD ratio is sub-linear vs the multiplier (e.g., m=15 → MRSD ratio ~15;
     m=8 → MRSD ratio ~9)
3. clint_liver_L_h ratio (enhanced/base) exactly matches the multiplier value.
"""

from __future__ import annotations

import copy
from pathlib import Path

import pytest
import yaml

from charon import Pipeline
from charon.core.schema import CompoundConfig, DoseProjectionConfig

REPO_ROOT = Path(__file__).resolve().parents[2]
PROPRA_YAML = REPO_ROOT / "validation" / "data" / "tier1_obach" / "compounds" / "propranolol.yaml"
PANEL_YAML = REPO_ROOT / "validation" / "data" / "fih_reference" / "panel.yaml"


def _load_propranolol_entry() -> dict:
    panel = yaml.safe_load(PANEL_YAML.read_text())["panel"]
    return next(
        c for c in panel["compounds"]
        if c["name"] == "propranolol" and c["tier"] == "gold"
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


def test_propranolol_mrsd_with_multiplier_larger_than_without():
    """With clint_multiplier > 1, propranolol MRSD must be strictly larger."""
    data = yaml.safe_load(PROPRA_YAML.read_text())
    entry = _load_propranolol_entry()

    compound_enhanced = CompoundConfig.model_validate(data)
    assert compound_enhanced.properties.metabolism.hepatic_clint_multiplier is not None, (
        "propranolol.yaml must have hepatic_clint_multiplier populated (Task 2 prerequisite)"
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


def test_propranolol_mrsd_ratio_within_literature_range():
    """MRSD_enhanced / MRSD_base within [4, 25] — band reflects CYP2D6 literature
    uncertainty (5-20x typical). For high-extraction propranolol (post-multiplier),
    well-stirred saturation makes the MRSD ratio approximately equal to the
    multiplier in the moderate-extraction regime, decreasing slightly at very
    high multipliers."""
    data = yaml.safe_load(PROPRA_YAML.read_text())
    entry = _load_propranolol_entry()

    compound_enhanced = CompoundConfig.model_validate(data)
    mrsd_enhanced, _ = _run_pipeline(compound_enhanced, entry)

    data_base = copy.deepcopy(data)
    data_base["properties"]["metabolism"].pop("hepatic_clint_multiplier", None)
    compound_base = CompoundConfig.model_validate(data_base)
    mrsd_base, _ = _run_pipeline(compound_base, entry)

    ratio = mrsd_enhanced / mrsd_base
    assert 4.0 <= ratio <= 25.0, (
        f"MRSD ratio {ratio:.2f} outside expected [4, 25] range; "
        f"check multiplier value or extraction saturation"
    )


def test_propranolol_enhanced_clint_liver_metadata_scales_by_multiplier():
    """clint_liver_L_h ratio (enhanced/base) should exactly match the multiplier.
    Pure linear scaling — ConversionStep multiplies CLint_liver before the
    liver model, so the reported clint_liver_L_h in metadata is the enhanced
    value."""
    data = yaml.safe_load(PROPRA_YAML.read_text())
    entry = _load_propranolol_entry()

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
