"""Sprint 6 — Report file exporters.

Writes a :class:`ReportData` to a Markdown file and a sibling JSON file.
No external dependencies: uses the ``json`` stdlib module and a small
``_json_default`` helper to handle numpy scalars.
"""

from __future__ import annotations

import json
from dataclasses import asdict
from pathlib import Path
from typing import Any

import numpy as np

from charon.report.collector import ReportData
from charon.report.narrative import render_report


def _json_default(obj: Any) -> Any:
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.floating):
        return float(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, Path):
        return str(obj)
    raise TypeError(f"Unsupported type for JSON serialization: {type(obj)!r}")


def export_markdown(data: ReportData, path: Path) -> Path:
    """Render *data* as Markdown and write to *path*.

    Returns the resolved absolute path.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(render_report(data), encoding="utf-8")
    return path.resolve()


def export_json(
    data: ReportData,
    path: Path,
    *,
    include_full_profile: bool = False,
    full_profile: dict | None = None,
) -> Path:
    """Serialise *data* to JSON at *path*.  Returns the resolved path."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = asdict(data)
    if include_full_profile and full_profile is not None:
        payload["full_profile"] = full_profile
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=False, default=_json_default),
        encoding="utf-8",
    )
    return path.resolve()


def export_report(
    data: ReportData,
    output: Path,
    *,
    full_profile: dict | None = None,
) -> tuple[Path, Path]:
    """Write both the Markdown and JSON siblings of a report.

    If *output* ends in ``.md`` the JSON path replaces the suffix with
    ``.json``; otherwise both ``.md`` and ``.json`` are appended to the
    given path.
    """
    output = Path(output)
    if output.suffix.lower() == ".md":
        md_path = output
        json_path = output.with_suffix(".json")
    else:
        md_path = output.with_name(output.name + ".md")
        json_path = output.with_name(output.name + ".json")

    md_resolved = export_markdown(data, md_path)
    json_resolved = export_json(
        data,
        json_path,
        include_full_profile=full_profile is not None,
        full_profile=full_profile,
    )
    return md_resolved, json_resolved
