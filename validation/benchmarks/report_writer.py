"""Shared benchmark report emitter for Charon validation suites.

Writes ``.md`` and ``.json`` report files from a structured dict payload.

Public API
----------
emit_report(payload, *, stem) -> tuple[Path, Path]
    Write ``{stem}.md`` and ``{stem}.json``, creating parent directories
    as needed. Returns resolved absolute paths.

Payload schema
--------------
Required keys:
    title       str          Report heading
    panel       str          Panel identifier (e.g. "obach_1999")
    date_utc    str          ISO-8601 UTC timestamp string
    summary     dict | list  Summary rows (dict or list of dicts)
    rows        list[dict]   Per-compound result rows
    notes       list[str]    Bullet-point notes

Optional keys:
    disclaimer  str          Shown as blockquote in Markdown
    targets     list[dict]   metric/target/met table
    excluded    list[dict]   Excluded-compound table

Cell formatting
---------------
- Numeric (float/int): 4 significant figures; sci notation when >=1000 or
  <0.01; NaN/inf -> "-"; bool -> "[PASS]" / "[FAIL]"
- String: escape ``\\`` -> ``\\\\``, ``|`` -> ``\\|``, newlines -> space;
  None -> "-"

JSON sanitisation
-----------------
Non-finite floats -> None; numpy scalars -> Python native;
datetime -> isoformat; Path -> str.
"""

from __future__ import annotations

import json
import math
from datetime import datetime
from pathlib import Path
from typing import Any

# ---------------------------------------------------------------------------
# Required payload keys
# ---------------------------------------------------------------------------

_REQUIRED_KEYS: tuple[str, ...] = (
    "title",
    "panel",
    "date_utc",
    "summary",
    "rows",
    "notes",
)


# ---------------------------------------------------------------------------
# History preservation (Sprint 16)
# ---------------------------------------------------------------------------

# Sentinel marker for opt-in history preservation. When an existing .md file
# contains this marker, emit_report preserves content from the marker through
# end-of-file across regenerations. Place it at the boundary between
# orchestrator-generated content (above) and user-managed sprint narratives
# (below). Marker absence = current overwrite behavior (backwards compatible).
_HISTORY_MARKER: str = "<!-- BEGIN_PRESERVED_HISTORY -->"


def _extract_preserved_history(md_path: Path) -> str:
    """Return content from history marker to EOF (inclusive of marker), or ''.

    Returns empty string if the file does not exist OR does not contain the
    marker. The marker itself is included in the returned content so it
    persists across regenerations.

    Parameters
    ----------
    md_path:
        Path to the existing markdown report file.

    Returns
    -------
    str
        Preserved content (marker + everything after it) or ``""`` if no
        preservation should occur.
    """
    if not md_path.exists():
        return ""
    content = md_path.read_text(encoding="utf-8")
    idx = content.find(_HISTORY_MARKER)
    if idx == -1:
        return ""
    return content[idx:]


# ---------------------------------------------------------------------------
# Cell formatting helpers
# ---------------------------------------------------------------------------

def _fmt_number(v: float | int) -> str:
    """Format a numeric value to 4 significant figures for Markdown."""
    if isinstance(v, bool):
        return "[PASS]" if v else "[FAIL]"
    try:
        f = float(v)
    except (TypeError, ValueError):
        return str(v)
    if not math.isfinite(f):
        return "-"
    if f == 0.0:
        return "0.000"
    abs_f = abs(f)
    if abs_f >= 1000 or abs_f < 0.01:
        return f"{f:.3e}"
    return f"{f:.4g}"


def _fmt_cell(v: Any) -> str:
    """Format a single table cell value for Markdown."""
    if v is None:
        return "-"
    if isinstance(v, bool):
        return "[PASS]" if v else "[FAIL]"
    if isinstance(v, (int, float)):
        return _fmt_number(v)
    s = str(v)
    # Escape backslashes first, then pipes, then flatten newlines
    s = s.replace("\\", "\\\\")
    s = s.replace("|", "\\|")
    s = s.replace("\n", " ").replace("\r", " ")
    return s


# ---------------------------------------------------------------------------
# Markdown rendering helpers
# ---------------------------------------------------------------------------

def _render_table(rows: list[dict]) -> list[str]:
    """Render a list of dicts as a Markdown table. Returns lines."""
    if not rows:
        return ["*(no data)*", ""]
    headers = list(rows[0].keys())
    lines: list[str] = []
    header_row = "| " + " | ".join(headers) + " |"
    sep_row = "| " + " | ".join(["---"] * len(headers)) + " |"
    lines.append(header_row)
    lines.append(sep_row)
    for row in rows:
        cells = [_fmt_cell(row.get(h)) for h in headers]
        lines.append("| " + " | ".join(cells) + " |")
    lines.append("")
    return lines


def _render_markdown(payload: dict) -> str:
    """Render the full Markdown document from payload."""
    lines: list[str] = []

    # Title and metadata
    lines.append(f"# {payload['title']}")
    lines.append("")
    lines.append(f"**Generated:** {payload['date_utc']}")
    lines.append(f"**Panel:** {payload['panel']}")
    lines.append("")

    # Disclaimer (optional)
    if "disclaimer" in payload and payload["disclaimer"]:
        lines.append(f"> {payload['disclaimer']}")
        lines.append("")

    # Summary section
    lines.append("## Summary")
    lines.append("")
    summary = payload["summary"]
    if isinstance(summary, list) and summary:
        lines.extend(_render_table(summary))
    elif isinstance(summary, dict) and summary:
        # Convert single-level dict to two-column table
        dict_rows = [{"key": k, "value": v} for k, v in summary.items()]
        lines.extend(_render_table(dict_rows))
    else:
        lines.append("*(no summary data)*")
        lines.append("")

    # Targets section (optional)
    if "targets" in payload and payload["targets"]:
        lines.append("## Targets")
        lines.append("")
        lines.extend(_render_table(payload["targets"]))

    # Per-compound results
    lines.append("## Results")
    lines.append("")
    rows = payload.get("rows", [])
    if rows:
        lines.extend(_render_table(rows))
    else:
        lines.append("*(no results)*")
        lines.append("")

    # Extra sections (generic dict[str, list[dict] | dict])
    for section_title, section_rows in payload.get("extra_sections", {}).items():
        lines.append(f"## {section_title}")
        lines.append("")
        if isinstance(section_rows, list) and section_rows:
            lines.extend(_render_table(section_rows))
        elif isinstance(section_rows, dict) and section_rows:
            flat = []
            for k, v in section_rows.items():
                if isinstance(v, dict):
                    flat.append({"key": k, **v})
                else:
                    flat.append({"key": k, "value": v})
            lines.extend(_render_table(flat))
        else:
            lines.append("*(no data)*")
            lines.append("")

    # Excluded section (optional)
    if "excluded" in payload and payload["excluded"]:
        lines.append("## Excluded")
        lines.append("")
        lines.extend(_render_table(payload["excluded"]))

    # Notes
    lines.append("## Notes")
    lines.append("")
    for note in payload.get("notes", []):
        lines.append(f"- {note}")
    lines.append("")

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# JSON sanitisation
# ---------------------------------------------------------------------------

def _sanitise(obj: Any) -> Any:
    """Recursively sanitise payload for JSON serialisation.

    - Non-finite floats (NaN, inf, -inf) -> None
    - numpy scalar types -> Python native (int/float/bool)
    - datetime -> isoformat string
    - Path -> str
    - dict/list -> recursed
    """
    # numpy detection without importing numpy (avoids hard dependency)
    type_name = type(obj).__module__
    if type_name == "numpy" or (
        hasattr(obj, "item") and hasattr(obj, "dtype")
    ):
        # numpy scalar — convert to Python native
        return _sanitise(obj.item())

    if isinstance(obj, bool):
        return obj
    if isinstance(obj, float):
        if not math.isfinite(obj):
            return None
        return obj
    if isinstance(obj, int):
        return obj
    if isinstance(obj, str):
        return obj
    if isinstance(obj, datetime):
        return obj.isoformat()
    if isinstance(obj, Path):
        return str(obj)
    if isinstance(obj, dict):
        return {k: _sanitise(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_sanitise(v) for v in obj]
    # Fallback: attempt str conversion
    return obj


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def emit_report(payload: dict, *, stem: Path) -> tuple[Path, Path]:
    """Emit benchmark report as ``{stem}.md`` and ``{stem}.json``.

    Parameters
    ----------
    payload:
        Report data dict. Must contain all required keys:
        ``title``, ``panel``, ``date_utc``, ``summary``, ``rows``, ``notes``.
    stem:
        File path without extension. Parent directories are created
        automatically.

    Returns
    -------
    tuple[Path, Path]
        Resolved absolute paths to the Markdown and JSON files.

    Raises
    ------
    ValueError
        If any required key is absent from *payload*.
    """
    # Validate required keys
    missing = [k for k in _REQUIRED_KEYS if k not in payload]
    if missing:
        raise ValueError(
            f"missing required key(s) in payload: {', '.join(missing)}"
        )

    stem = Path(stem).resolve()
    stem.parent.mkdir(parents=True, exist_ok=True)

    md_path = stem.with_suffix(".md")
    json_path = stem.with_suffix(".json")

    # Write Markdown — preserve history below the sentinel marker if present
    preserved = _extract_preserved_history(md_path)
    md_content = _render_markdown(payload)
    if preserved:
        md_content = md_content.rstrip() + "\n\n" + preserved
    md_path.write_text(md_content, encoding="utf-8")

    # Write JSON
    sanitised = _sanitise(payload)
    json_content = json.dumps(sanitised, indent=2, allow_nan=False)
    json_path.write_text(json_content, encoding="utf-8")

    return md_path, json_path
