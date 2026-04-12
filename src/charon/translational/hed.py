"""HED (Human Equivalent Dose) computation via FDA BSA scaling.

    HED [mg/kg] = NOAEL [mg/kg] × (Km_animal / Km_human)
    MRSD [mg] = HED × body_weight_kg / safety_factor

References:
    FDA Guidance for Industry (2005): Estimating the Maximum Safe Starting
    Dose in Initial Clinical Trials for Therapeutics in Adult Healthy
    Volunteers. Table 1 (Km values by species).
"""

from __future__ import annotations

from dataclasses import dataclass

KM_BY_SPECIES: dict[str, float] = {
    "human": 37.0,
    "rat": 6.2,
    "mouse": 3.0,
    "dog": 20.0,
    "monkey": 12.0,
    "rabbit": 12.0,
    "guinea_pig": 8.0,
}


@dataclass(frozen=True)
class HEDResult:
    noael_mg_kg: float
    noael_species: str
    km_animal: float
    km_human: float
    hed_mg_kg: float
    body_weight_kg: float
    safety_factor: float
    mrsd_mg: float


def compute_hed(
    *,
    noael_mg_kg: float,
    noael_species: str,
    safety_factor: float = 10.0,
    body_weight_kg: float = 70.0,
) -> HEDResult:
    """Compute HED and MRSD from animal NOAEL using FDA BSA scaling.

    Args:
        noael_mg_kg: Animal NOAEL in mg/kg. Must be > 0.
        noael_species: Species name (case-insensitive). Must be in
            KM_BY_SPECIES.
        safety_factor: Safety factor divisor applied to HED×BW.
            Default = 10 (FDA recommended for healthy volunteer FIH).
        body_weight_kg: Human body weight in kg. Default = 70.

    Returns:
        HEDResult with hed_mg_kg and mrsd_mg.

    Formula:
        km_animal / km_human converts animal mg/kg to human mg/kg via
        body-surface-area (BSA) normalisation.
        MRSD [mg] = (NOAEL × km_animal/km_human × BW_kg) / safety_factor
    """
    if noael_mg_kg <= 0:
        raise ValueError(f"noael_mg_kg must be > 0, got {noael_mg_kg}")
    species_key = noael_species.lower()
    if species_key not in KM_BY_SPECIES:
        raise ValueError(
            f"Unknown species {noael_species!r}. "
            f"Supported: {sorted(KM_BY_SPECIES.keys())}"
        )
    km_animal = KM_BY_SPECIES[species_key]
    km_human = KM_BY_SPECIES["human"]
    hed_mg_kg = noael_mg_kg * (km_animal / km_human)
    mrsd_mg = hed_mg_kg * body_weight_kg / safety_factor
    return HEDResult(
        noael_mg_kg=noael_mg_kg,
        noael_species=species_key,
        km_animal=km_animal,
        km_human=km_human,
        hed_mg_kg=hed_mg_kg,
        body_weight_kg=body_weight_kg,
        safety_factor=safety_factor,
        mrsd_mg=mrsd_mg,
    )
