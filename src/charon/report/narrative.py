"""Sprint 6 — Markdown narrative rendering.

Pure f-string based rendering of a :class:`ReportData` into a
regulatory-oriented Markdown document.  No template engine, no new
runtime dependencies.
"""

from __future__ import annotations

import math

from charon.report.collector import ReportData


def format_value(v: float | int | None, *, digits: int = 4) -> str:
    """Format a numeric value for a Markdown table cell.

    - ``None``, ``NaN``, ``±inf`` → ``"-"``
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
    if not math.isfinite(fv):
        return "-"
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


_ADME_ROW_ORDER: tuple[str, ...] = (
    "logp",
    "pka_acid",
    "pka_base",
    "solubility_ug_ml",
    "fu_p",
    "fu_inc",
    "bp_ratio",
    "clint_uL_min_mg",
    "papp_nm_s",
    "peff_cm_s",
    "herg_ic50_uM",
    "clrenal_L_h",
)


def _md_cell(s: object | None) -> str:
    """Escape a free-form string for safe use inside a Markdown table cell.

    Replaces ``|`` with ``\\|`` and newlines with spaces.  Returns ``"-"``
    for ``None``.
    """
    if s is None:
        return "-"
    return str(s).replace("|", "\\|").replace("\n", " ")


def _ci_cell(entry: dict) -> str:
    lo = entry.get("ci_lower")
    hi = entry.get("ci_upper")
    if lo is None or hi is None:
        return "-"
    return f"[{format_value(lo)}, {format_value(hi)}]"


def _render_adme_table(data: ReportData) -> str:
    lines = ["## 3. ADME Predictions", ""]
    if not data.properties:
        lines.append("*No ADME properties recorded.*")
        return "\n".join(lines)

    lines.append("| Property | Value | 90% CI | Unit | Source | Flag |")
    lines.append("| --- | --- | --- | --- | --- | --- |")
    for key in _ADME_ROW_ORDER:
        entry = data.properties.get(key)
        if entry is None:
            continue
        lines.append(
            "| {name} | {val} | {ci} | {unit} | {src} | {flag} |".format(
                name=_md_cell(key),
                val=format_value(entry.get("value")),
                ci=_ci_cell(entry),
                unit=_md_cell(entry.get("unit")),
                src=_md_cell(entry.get("source")),
                flag=_md_cell(entry.get("flag")),
            )
        )
    return "\n".join(lines)


def _render_ivive_audit(data: ReportData) -> str:
    lines = ["## 4. IVIVE & Hepatic Clearance", ""]
    ivive = data.ivive_summary or {}
    if not ivive:
        lines.append("*IVIVE summary not available in run metadata.*")
        return "\n".join(lines)

    model = ivive.get("liver_model", "-")
    clint_liver = ivive.get("clint_liver_L_h")
    cl_renal = ivive.get("cl_renal_L_h")
    fu_b = ivive.get("fu_b")
    cl_app = data.pk_params.get("cl_apparent")

    lines.append(
        f"CLint was scaled via the **{model}** model. Whole-liver "
        f"CLint = {format_value(clint_liver)} L/h. "
        f"fu_b = {format_value(fu_b)}. "
        f"Renal CL = {format_value(cl_renal)} L/h. "
        f"In-vivo apparent CL from simulation = "
        f"{format_value(cl_app)} L/h."
    )
    lines.append("")
    lines.append(f"- **Liver model:** {model}")
    lines.append(f"- **CLint (whole liver):** {format_value(clint_liver)} L/h")
    lines.append(f"- **fu_b:** {format_value(fu_b)}")
    lines.append(f"- **Renal CL:** {format_value(cl_renal)} L/h")
    lines.append(f"- **In-vivo apparent CL:** {format_value(cl_app)} L/h")
    return "\n".join(lines)


_PK_PARAM_LABELS: tuple[tuple[str, str, str], ...] = (
    ("cmax", "Cmax", "ug/L"),
    ("tmax", "Tmax", "h"),
    ("auc_0_inf", "AUC(0-inf)", "ug*h/L"),
    ("auc_0_24", "AUC(0-24)", "ug*h/L"),
    ("half_life", "t1/2", "h"),
    ("cl_apparent", "CL apparent", "L/h"),
    ("vss", "Vss", "L"),
    ("bioavailability", "F", "-"),
    ("fa", "Fa", "-"),
    ("fg", "Fg", "-"),
    ("fh", "Fh", "-"),
)


def _render_pk_results(data: ReportData) -> str:
    lines = [
        "## 5. PK Simulation Results",
        "",
        f"Route: **{_md_cell(data.route)}**  ·  "
        f"Dose: **{format_value(data.dose_mg)} mg**  ·  "
        f"Duration: **{format_value(data.duration_h)} h**",
        "",
        "### PK Parameters",
        "",
        "| Parameter | Value | Unit |",
        "| --- | --- | --- |",
    ]
    for key, label, unit in _PK_PARAM_LABELS:
        val = data.pk_params.get(key)
        lines.append(
            f"| {_md_cell(label)} | {format_value(val)} | {_md_cell(unit)} |"
        )

    lines.append("")
    lines.append("### Concentration-Time Profile (canonical timepoints)")
    lines.append("")
    if not data.pk_table:
        lines.append("*No time-course data available.*")
        return "\n".join(lines)

    lines.append("| Time (h) | Cp plasma (ug/L) | Cp blood (ug/L) |")
    lines.append("| --- | --- | --- |")
    for row in data.pk_table:
        lines.append(
            f"| {format_value(row['time_h'])} "
            f"| {format_value(row['cp_plasma_ug_L'])} "
            f"| {format_value(row['cp_blood_ug_L'])} |"
        )
    return "\n".join(lines)


def _dose_method_row(
    name_display: str,
    method_key: str,
    rec: dict,
    limiting: str,
) -> str:
    sub = rec.get(method_key)
    if sub is None:
        return f"| {_md_cell(name_display)} | - | - | insufficient inputs |"
    mrsd = sub.get("mrsd_mg")
    val_cell = format_value(mrsd)
    if method_key == limiting:
        val_cell = f"**{val_cell}**"
    inputs = {
        "hed": f"NOAEL={format_value(sub.get('noael_mg_kg'))} mg/kg "
               f"({_md_cell(sub.get('noael_species') or sub.get('species'))})",
        "mabel": f"Kd={format_value(sub.get('target_kd_nM'))} nM",
        "pad": f"Ceff={format_value(sub.get('target_ceff_nM'))} nM",
    }.get(method_key, "-")
    return f"| {_md_cell(name_display)} | {val_cell} | {inputs} | computed |"


def _render_dose_projection(data: ReportData) -> str:
    lines = ["## 6. FIH Dose Projection", ""]
    rec = data.dose_recommendation
    if rec is None:
        lines.append(
            "*FIH dose projection was not run for this pipeline execution.*"
        )
        return "\n".join(lines)

    limiting = str(rec.get("limiting_method", ""))
    lines.append("| Method | MRSD (mg) | Inputs | Status |")
    lines.append("| --- | --- | --- | --- |")
    lines.append(_dose_method_row("HED", "hed", rec, limiting))
    lines.append(_dose_method_row("MABEL", "mabel", rec, limiting))
    lines.append(_dose_method_row("PAD", "pad", rec, limiting))
    lines.append("")
    lines.append(
        f"Safety factor: **{format_value(rec.get('safety_factor'))}**  ·  "
        f"Salt factor: **{format_value(rec.get('salt_factor'))}**  ·  "
        f"Limiting method: **{limiting.upper()}**"
    )
    rationale = rec.get("rationale") or ""
    if rationale:
        lines.append("")
        lines.append("```")
        lines.append(rationale)
        lines.append("```")
    return "\n".join(lines)
