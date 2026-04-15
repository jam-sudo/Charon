"""Sprint 6 — Markdown narrative rendering.

Pure f-string based rendering of a :class:`ReportData` into a
regulatory-oriented Markdown document.  No template engine, no new
runtime dependencies.
"""

from __future__ import annotations

from charon.report.collector import ReportData


def format_value(v: float | int | None, *, digits: int = 4) -> str:
    """Format a numeric value for a Markdown table cell.

    - ``None`` → ``"-"``
    - ``abs(v) >= 1000`` or ``abs(v) < 0.01`` (and not zero) → scientific
      with ``digits-1`` significant figures.
    - otherwise → fixed-ish with ``digits`` significant figures.
    """
    if v is None:
        return "-"
    try:
        fv = float(v)
    except (TypeError, ValueError):
        return str(v)
    if fv == 0.0:
        return "0"
    a = abs(fv)
    if a >= 1000.0 or a < 0.01:
        return f"{fv:.{digits - 1}e}"
    return f"{fv:.{digits}g}"


def _render_header(data: ReportData) -> str:
    return (
        f"# FIH Dose Rationale Report — {data.compound_name}\n\n"
        f"*Generated: {data.timestamp} · "
        f"Charon {data.charon_version}*"
    )


def _render_executive_summary(data: ReportData) -> str:
    lines: list[str] = ["## 1. Executive Summary", ""]
    rec = data.dose_recommendation
    unc = data.uncertainty

    if rec is None:
        lines.append(
            "No FIH dose projection was run for this pipeline execution "
            "(no NOAEL, target Kd, or target Ceff provided)."
        )
    else:
        mrsd = format_value(rec.get("mrsd_mg"))
        method = str(rec.get("limiting_method", "?")).upper()
        route = rec.get("route", data.route)
        headline = (
            f"**Recommended FIH starting dose: {mrsd} mg** "
            f"({route}, limiting method: {method})"
        )
        if unc is not None:
            lo = format_value(unc.get("ci_90_lower_mg"))
            hi = format_value(unc.get("ci_90_upper_mg"))
            conf = unc.get("confidence", "?")
            headline += (
                f"  \n90% CI: [{lo} – {hi}] mg  ·  confidence: **{conf}**"
            )
        lines.append(headline)
        lines.append("")
        lines.append(
            f"Compound: {data.compound_name} ({data.smiles}). "
            f"Source: {data.source}. Route: {data.route}, "
            f"simulated dose: {format_value(data.dose_mg)} mg."
        )

    return "\n".join(lines)


def _render_compound_profile(data: ReportData) -> str:
    lines = ["## 2. Compound Profile", ""]
    lines.append(f"- **Name:** {data.compound_name}")
    lines.append(f"- **SMILES:** `{data.smiles}`")
    lines.append(
        f"- **Molecular weight:** "
        f"{format_value(data.molecular_weight, digits=5)} g/mol"
    )
    lines.append(f"- **Compound type:** {data.compound_type or '-'}")
    lines.append(f"- **Source:** {data.source}")
    return "\n".join(lines)
