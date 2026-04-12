"""ACAT GI tract model — segment data, loader, absorption rate computation.

Implements the 8-segment GI lumen transit model for oral drug absorption.
Segment parameters are loaded from the species YAML ``gi_tract`` section.

The absorption rate per segment uses the mechanistic cylindrical model:
    k_abs_i [1/h] = (2 × Peff [cm/s] × 3600 / R_i [cm]) × ka_fraction_i

References:
    Yu LX, Amidon GL (1999). Int J Pharm 186:119-125.
    Sun D et al. (2002). J Pharm Sci 91:1396-1404.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

import yaml

_SPECIES_DIR = Path(__file__).parent / "species"


@dataclass(frozen=True)
class GISegment:
    """Single GI lumen segment."""

    name: str
    volume_L: float
    radius_cm: float
    ka_fraction: float
    transit_rate_1_h: float


@dataclass(frozen=True)
class GITract:
    """Complete GI tract physiology for ACAT model."""

    segments: tuple[GISegment, ...]
    enterocyte_volume_L: float
    enterocyte_weight_g: float
    q_villi_fraction: float
    mppgi_mg_g: float
    cyp3a4_gut_pmol_per_mg: float
    cyp3a4_liver_pmol_per_mg: float


def load_gi_tract(species: str) -> GITract:
    """Load GI tract parameters from species YAML.

    Parameters
    ----------
    species : str
        Species name (e.g. ``"human"``).

    Returns
    -------
    GITract
        Frozen GI tract data container.

    Raises
    ------
    FileNotFoundError
        If the species YAML does not exist.
    KeyError
        If the YAML lacks a ``gi_tract`` section.
    """
    yaml_path = _SPECIES_DIR / f"{species}.yaml"
    if not yaml_path.exists():
        raise FileNotFoundError(f"Species YAML not found: {yaml_path}")

    with yaml_path.open() as fp:
        data = yaml.safe_load(fp)

    gi = data["species"]["gi_tract"]
    segments: list[GISegment] = []
    for name, spec in gi["segments"].items():
        segments.append(
            GISegment(
                name=name,
                volume_L=float(spec["volume_L"]),
                radius_cm=float(spec["radius_cm"]),
                ka_fraction=float(spec["ka_fraction"]),
                transit_rate_1_h=float(spec["transit_rate_1_h"]),
            )
        )

    return GITract(
        segments=tuple(segments),
        enterocyte_volume_L=float(gi["enterocyte_volume_L"]),
        enterocyte_weight_g=float(gi["enterocyte_weight_g"]),
        q_villi_fraction=float(gi["q_villi_fraction"]),
        mppgi_mg_g=float(gi["mppgi_mg_g"]),
        cyp3a4_gut_pmol_per_mg=float(gi["cyp3a4_gut_pmol_per_mg"]),
        cyp3a4_liver_pmol_per_mg=float(gi["cyp3a4_liver_pmol_per_mg"]),
    )


def compute_absorption_rates(
    gi: GITract,
    peff_cm_s: float,
) -> tuple[float, ...]:
    """Compute per-segment absorption rate constants.

    k_abs_i [1/h] = (2 × Peff [cm/s] × 3600 / R_i [cm]) × ka_fraction_i

    Parameters
    ----------
    gi : GITract
        GI tract physiology.
    peff_cm_s : float
        Effective intestinal permeability in cm/s.

    Returns
    -------
    tuple[float, ...]
        k_abs for each segment in 1/h, same order as ``gi.segments``.
    """
    rates: list[float] = []
    for seg in gi.segments:
        if seg.ka_fraction == 0.0:
            rates.append(0.0)
        else:
            k_abs = (2.0 * peff_cm_s * 3600.0 / seg.radius_cm) * seg.ka_fraction
            rates.append(k_abs)
    return tuple(rates)


def papp_to_peff(papp_nm_s: float) -> float:
    """Convert Caco-2 Papp (nm/s) to human Peff (cm/s).

    Uses the Sun et al. 2002 simplified correlation:
        log10(Peff) = 0.4926 × log10(Papp_nm_s) - 0.1454

    Parameters
    ----------
    papp_nm_s : float
        Apparent permeability from Caco-2 assay, nm/s.

    Returns
    -------
    float
        Effective intestinal permeability, cm/s.

    Raises
    ------
    ValueError
        If papp_nm_s <= 0.
    """
    if papp_nm_s <= 0.0:
        raise ValueError(f"papp_nm_s must be > 0, got {papp_nm_s}")
    log_peff = 0.4926 * math.log10(papp_nm_s) - 0.1454
    return 10.0 ** log_peff
