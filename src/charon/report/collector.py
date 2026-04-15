"""Sprint 6 — Report data collection.

Flattens a :class:`charon.pipeline.PipelineResult` into a
:class:`ReportData` dataclass suitable for consumption by the narrative
renderer and the JSON exporter.  The collector performs no numerical
computation beyond indexing and copying — any derived metric must be
computed upstream in the pipeline.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime, timezone

from charon.core.schema import PredictedProperty
from charon.pipeline import PipelineResult


_PROPERTY_FIELDS: list[tuple[str, str, str]] = [
    # (category, attribute_on_category, report_key)
    ("physicochemical", "logp", "logp"),
    ("physicochemical", "pka_acid", "pka_acid"),
    ("physicochemical", "pka_base", "pka_base"),
    ("physicochemical", "solubility_ug_ml", "solubility_ug_ml"),
    ("binding", "fu_p", "fu_p"),
    ("binding", "fu_inc", "fu_inc"),
    ("binding", "bp_ratio", "bp_ratio"),
    ("metabolism", "clint_uL_min_mg", "clint_uL_min_mg"),
    ("permeability", "papp_nm_s", "papp_nm_s"),
    ("permeability", "peff_cm_s", "peff_cm_s"),
    ("safety", "herg_ic50_uM", "herg_ic50_uM"),
    ("renal", "clrenal_L_h", "clrenal_L_h"),
]


def _flatten_property(prop: PredictedProperty) -> dict:
    return {
        "value": float(prop.value),
        "ci_lower": None if prop.ci_90_lower is None else float(prop.ci_90_lower),
        "ci_upper": None if prop.ci_90_upper is None else float(prop.ci_90_upper),
        "unit": prop.unit,
        "source": prop.source,
        "flag": prop.flag,
        "method": prop.method,
    }


def _flatten_properties(props) -> dict[str, dict]:
    out: dict[str, dict] = {}
    for category, attr, key in _PROPERTY_FIELDS:
        section = getattr(props, category, None)
        if section is None:
            continue
        prop = getattr(section, attr, None)
        if prop is None:
            continue
        out[key] = _flatten_property(prop)
    return out


_IVIVE_KEYS = (
    "clint_liver_L_h",
    "cl_renal_L_h",
    "fu_b",
    "liver_model",
    "compound_type",
    "clint_gut_L_h",
)


def _ivive_summary_from_metadata(md: dict) -> dict:
    return {k: md[k] for k in _IVIVE_KEYS if k in md}


@dataclass(frozen=True)
class ReportData:
    """Structured data for the Sprint 6 regulatory report.

    All fields are plain Python values (no numpy arrays, no Pydantic
    models).  This makes the dataclass trivially JSON-serialisable.
    """

    # Identity
    compound_name: str
    smiles: str
    molecular_weight: float | None
    source: str
    compound_type: str | None

    # Layer 1: ADME properties, flattened
    properties: dict[str, dict]

    # Bridge: IVIVE audit fields (pulled verbatim from metadata)
    ivive_summary: dict

    # Layer 2: PK parameters + sampled Cp-time table
    pk_params: dict[str, float | None]
    pk_table: list[dict]
    route: str
    dose_mg: float
    duration_h: float

    # Layer 3 / Layer 4 (optional)
    dose_recommendation: dict | None
    uncertainty: dict | None

    # Layer 0 + run metadata
    warnings: list[str] = field(default_factory=list)
    metadata: dict = field(default_factory=dict)
    timestamp: str = ""
    charon_version: str = "0.1.0"


def _iso_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def collect(
    result: PipelineResult,
    *,
    warnings: list[str] | None = None,
    timestamp: str | None = None,
) -> ReportData:
    """Flatten a ``PipelineResult`` into a ``ReportData`` snapshot.

    Subsequent tasks extend this function to populate ADME properties,
    IVIVE summary, PK table, dose recommendation, and uncertainty.
    """
    md = result.metadata or {}
    compound = result.compound

    return ReportData(
        compound_name=compound.name,
        smiles=compound.smiles,
        molecular_weight=compound.molecular_weight,
        source=compound.source,
        compound_type=compound.properties.physicochemical.compound_type,
        properties=_flatten_properties(compound.properties),
        ivive_summary=_ivive_summary_from_metadata(md),
        pk_params={},
        pk_table=[],
        route=str(md.get("route", "")),
        dose_mg=float(md.get("dose_mg", 0.0)),
        duration_h=float(md.get("duration_h", 0.0)),
        dose_recommendation=None,
        uncertainty=None,
        warnings=list(warnings) if warnings else [],
        metadata=dict(md),
        timestamp=timestamp or _iso_now(),
    )
