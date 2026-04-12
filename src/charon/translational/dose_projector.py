"""FIH Dose Projector — coordinator for HED, MABEL, and PAD methods.

Runs all applicable methods based on available inputs, selects the most
conservative (minimum) MRSD, and produces a structured recommendation.
"""

from __future__ import annotations

import math

from pydantic import BaseModel, ConfigDict

from charon.core.schema import (
    CompoundConfig,
    DoseProjectionConfig,
    PKParameters,
)
from charon.translational.hed import HEDResult, compute_hed
from charon.translational.mabel import MABELResult, compute_mabel
from charon.translational.pad import PADResult, compute_pad


class FIHDoseRecommendation(BaseModel):
    model_config = ConfigDict(frozen=True)

    mrsd_mg: float
    limiting_method: str
    hed: HEDResult | None = None
    mabel: MABELResult | None = None
    pad: PADResult | None = None
    safety_factor: float
    salt_factor: float = 1.0
    route: str = "oral"
    rationale: str


def project_fih_dose(
    *,
    pk: PKParameters,
    compound: CompoundConfig,
    config: DoseProjectionConfig,
    route: str = "oral",
) -> FIHDoseRecommendation:
    """Coordinate FIH dose projection across HED, MABEL, and PAD methods.

    Runs whichever methods have sufficient inputs, then selects the most
    conservative (minimum) MRSD as the recommended starting dose per FDA
    guidance (MRSD = min of all applicable estimates).

    Args:
        pk: PK parameters from PBPK simulation or experimental data.
            cl_apparent (L/h) and half_life (h) are used for MABEL/PAD.
        compound: Compound configuration including molecular_weight,
            salt_form, and binding.fu_p.
        config: Dose-projection inputs: noael_mg_kg/noael_species for HED,
            target_kd_nM for MABEL, target_ceff_nM for PAD.
        route: Administration route; "oral" uses CL/F, "iv" uses CL.

    Returns:
        FIHDoseRecommendation with mrsd_mg (salt-corrected), limiting_method,
        individual method results, and rationale text.

    Raises:
        ValueError: If no method has sufficient inputs to run.

    Safety rule:
        MRSD is always the minimum across all computed methods.  Salt
        correction (dose_free = dose_salt × salt_factor) is applied after
        method selection so the reported mrsd_mg is in free-base equivalents.
    """
    sf = config.safety_factor
    mw = compound.molecular_weight or 0.0
    salt_factor = 1.0
    if compound.salt_form is not None:
        salt_factor = compound.salt_form.salt_factor

    # Resolve apparent PK parameters
    cl_apparent = pk.cl_apparent or 0.0
    if route == "oral" or pk.vss is None:
        half_life = pk.half_life or 1.0
        vd_apparent = (
            half_life * cl_apparent / math.log(2) if cl_apparent > 0 else 0.0
        )
    else:
        vd_apparent = pk.vss

    fu_p: float | None = None
    if compound.properties.binding.fu_p is not None:
        fu_p = compound.properties.binding.fu_p.value

    # Run available methods
    hed_result: HEDResult | None = None
    mabel_result: MABELResult | None = None
    pad_result: PADResult | None = None
    candidates: list[tuple[str, float]] = []

    if config.noael_mg_kg is not None and config.noael_species is not None:
        hed_result = compute_hed(
            noael_mg_kg=config.noael_mg_kg,
            noael_species=config.noael_species,
            safety_factor=sf,
            body_weight_kg=config.body_weight_kg,
        )
        candidates.append(("hed", hed_result.mrsd_mg * salt_factor))

    if (
        config.target_kd_nM is not None
        and fu_p is not None
        and fu_p > 0
        and mw > 0
    ):
        mabel_result = compute_mabel(
            target_kd_nM=config.target_kd_nM,
            molecular_weight=mw,
            fu_p=fu_p,
            cl_apparent_L_h=cl_apparent,
            vd_apparent_L=vd_apparent,
            safety_factor=sf,
            tau_h=config.tau_h,
        )
        candidates.append(("mabel", mabel_result.mrsd_mg * salt_factor))

    if config.target_ceff_nM is not None and mw > 0 and cl_apparent > 0:
        pad_result = compute_pad(
            target_ceff_nM=config.target_ceff_nM,
            molecular_weight=mw,
            cl_apparent_L_h=cl_apparent,
            safety_factor=sf,
            tau_h=config.tau_h,
        )
        candidates.append(("pad", pad_result.mrsd_mg * salt_factor))

    if not candidates:
        raise ValueError(
            "FIH dose projection requires at least one of: "
            "noael_mg_kg+noael_species, target_kd_nM, or target_ceff_nM"
        )

    limiting_method, mrsd_mg = min(candidates, key=lambda x: x[1])

    # Build rationale
    lines = [f"FIH dose recommendation: {mrsd_mg:.2f} mg"]
    for name, val in candidates:
        marker = " <- limiting" if name == limiting_method else ""
        lines.append(f"  {name.upper()}: {val:.2f} mg (SF={sf}){marker}")
    ran_names = {c[0] for c in candidates}
    for name in ("hed", "mabel", "pad"):
        if name not in ran_names:
            lines.append(
                f"  {name.upper()}: not computed (insufficient inputs)"
            )
    if salt_factor != 1.0:
        lines.append(f"Salt correction applied: factor={salt_factor:.4f}")
    lines.append(f"Most conservative method: {limiting_method.upper()}")
    rationale = "\n".join(lines)

    return FIHDoseRecommendation(
        mrsd_mg=mrsd_mg,
        limiting_method=limiting_method,
        hed=hed_result,
        mabel=mabel_result,
        pad=pad_result,
        safety_factor=sf,
        salt_factor=salt_factor,
        route=route,
        rationale=rationale,
    )
