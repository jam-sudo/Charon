"""Unit conversion registry for translational PK calculations.

Provides physiological reference constants and pure, stateless conversion
functions between common pharmacokinetic unit systems.  Every converter
enforces non-negativity where the physical quantity cannot be negative
(flows, clearances, concentrations, masses) but allows negative values for
properties such as logP.
"""

from __future__ import annotations

# ---------------------------------------------------------------------------
# Physiological reference constants (human, default)
# ---------------------------------------------------------------------------

HUMAN_GFR_ML_MIN: float = 120.0
"""Glomerular filtration rate (mL/min)."""

HUMAN_LIVER_WEIGHT_G: float = 1500.0
"""Reference liver weight (g)."""

HUMAN_QH_L_H: float = 90.0
"""Hepatic blood flow (L/h)."""

HUMAN_MPPGL: float = 40.0
"""Microsomal protein per gram of liver (mg protein / g liver)."""

HUMAN_HEPATOCELLULARITY: float = 120.0
"""Hepatocellularity (10^6 cells / g liver).

The 10^6 factor is implicit in the unit convention; the stored numeric value
is 120.0 representing 120 × 10^6 cells/g liver.
"""

HUMAN_BODY_WEIGHT_KG: float = 70.0
"""Standard human body weight (kg)."""

HUMAN_BSA_KM: float = 37.0
"""Body-weight / body-surface-area ratio used for HED conversion (kg/m^2)."""


# ---------------------------------------------------------------------------
# Helper – input guard
# ---------------------------------------------------------------------------

def _check_non_negative(value: float, name: str) -> None:
    """Raise ``ValueError`` if *value* is negative.

    Parameters
    ----------
    value : float
        The numeric value to validate.
    name : str
        Human-readable label used in the error message.
    """
    if value < 0:
        raise ValueError(f"{name} must be non-negative, got {value}")


# ---------------------------------------------------------------------------
# Volume-flow conversions
# ---------------------------------------------------------------------------

def uL_min_to_L_h(value: float) -> float:
    """Convert a flow rate from µL/min to L/h.

    Formula: ``value / 1e6 * 60``

    Parameters
    ----------
    value : float
        Flow rate in µL/min.

    Returns
    -------
    float
        Flow rate in L/h.

    Raises
    ------
    ValueError
        If *value* is negative.
    """
    _check_non_negative(value, "Flow rate (µL/min)")
    return value / 1e6 * 60


def mL_min_to_L_h(value: float) -> float:
    """Convert a flow rate from mL/min to L/h.

    Formula: ``value * 60 / 1000``

    Parameters
    ----------
    value : float
        Flow rate in mL/min.

    Returns
    -------
    float
        Flow rate in L/h.

    Raises
    ------
    ValueError
        If *value* is negative.
    """
    _check_non_negative(value, "Flow rate (mL/min)")
    return value * 60 / 1000


def L_h_to_mL_min(value: float) -> float:
    """Convert a flow rate from L/h to mL/min.

    Formula: ``value * 1000 / 60``

    Parameters
    ----------
    value : float
        Flow rate in L/h.

    Returns
    -------
    float
        Flow rate in mL/min.

    Raises
    ------
    ValueError
        If *value* is negative.
    """
    _check_non_negative(value, "Flow rate (L/h)")
    return value * 1000 / 60


# ---------------------------------------------------------------------------
# Permeability conversions
# ---------------------------------------------------------------------------

def nm_s_to_cm_s(value: float) -> float:
    """Convert permeability from nm/s to cm/s.

    Formula: ``value * 1e-7``

    Parameters
    ----------
    value : float
        Permeability in nm/s.

    Returns
    -------
    float
        Permeability in cm/s.

    Raises
    ------
    ValueError
        If *value* is negative.
    """
    _check_non_negative(value, "Permeability (nm/s)")
    return value * 1e-7


def cm_s_to_nm_s(value: float) -> float:
    """Convert permeability from cm/s to nm/s.

    Formula: ``value / 1e-7``

    Parameters
    ----------
    value : float
        Permeability in cm/s.

    Returns
    -------
    float
        Permeability in nm/s.

    Raises
    ------
    ValueError
        If *value* is negative.
    """
    _check_non_negative(value, "Permeability (cm/s)")
    return value / 1e-7


# ---------------------------------------------------------------------------
# Dose / mass conversions
# ---------------------------------------------------------------------------

def mg_kg_to_mg(value: float, body_weight_kg: float) -> float:
    """Convert a weight-normalised dose (mg/kg) to absolute mass (mg).

    Formula: ``value * body_weight_kg``

    Parameters
    ----------
    value : float
        Dose in mg/kg.
    body_weight_kg : float
        Body weight in kg.

    Returns
    -------
    float
        Absolute dose in mg.

    Raises
    ------
    ValueError
        If *value* or *body_weight_kg* is negative.
    """
    _check_non_negative(value, "Dose (mg/kg)")
    _check_non_negative(body_weight_kg, "Body weight (kg)")
    return value * body_weight_kg


# ---------------------------------------------------------------------------
# Concentration conversions
# ---------------------------------------------------------------------------

def nmol_L_to_ng_mL(value: float, molecular_weight: float) -> float:
    """Convert concentration from nmol/L to ng/mL.

    Formula: ``value * molecular_weight / 1000``

    Parameters
    ----------
    value : float
        Concentration in nmol/L (nM).
    molecular_weight : float
        Molecular weight in g/mol (Da).

    Returns
    -------
    float
        Concentration in ng/mL.

    Raises
    ------
    ValueError
        If *value* or *molecular_weight* is negative.
    """
    _check_non_negative(value, "Concentration (nmol/L)")
    _check_non_negative(molecular_weight, "Molecular weight (g/mol)")
    return value * molecular_weight / 1000


def ng_mL_to_nmol_L(value: float, molecular_weight: float) -> float:
    """Convert concentration from ng/mL to nmol/L.

    Formula: ``value * 1000 / molecular_weight``

    Parameters
    ----------
    value : float
        Concentration in ng/mL.
    molecular_weight : float
        Molecular weight in g/mol (Da).

    Returns
    -------
    float
        Concentration in nmol/L (nM).

    Raises
    ------
    ValueError
        If *value* or *molecular_weight* is negative.
    """
    _check_non_negative(value, "Concentration (ng/mL)")
    _check_non_negative(molecular_weight, "Molecular weight (g/mol)")
    return value * 1000 / molecular_weight
