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

from charon.pipeline import PipelineResult


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
        properties={},
        ivive_summary={},
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
