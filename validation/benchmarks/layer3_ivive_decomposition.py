"""Sprint 10 — Tier A IVIVE fold-error decomposition.

See docs/superpowers/specs/2026-04-23-sprint10-ivive-bias-diagnostic-design.md.

For each of 12 Tier A compounds:
  1. Run Pipeline once (well_stirred baseline).
  2. Analytically compute MRSD for parallel_tube and dispersion what-ifs
     via CL-proportional scaling:
       mrsd_model = mrsd_ws * (CLh_model + cl_renal) / (CLh_ws + cl_renal)
  3. Load literature F from bioavailability.csv.
  4. Call decompose_fold_error().
  5. Aggregate into a report.

Note: Pipeline.liver_model is a no-op for the PBPK ODE, which embeds
well-stirred extraction directly (CLint_liver * fu_b * C_liver_blood_out).
Running Pipeline 3× per compound would yield identical MRSDs; the
analytical what-if approach correctly reflects the model divergence.

This is research-only — no production code changes, no gating.
Exit 0 on success, 1 on data errors.
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from datetime import datetime, timezone
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from charon import Pipeline  # noqa: E402
from charon.core.liver_models import dispersion, parallel_tube, well_stirred  # noqa: E402
from charon.core.schema import CompoundConfig, DoseProjectionConfig  # noqa: E402
from charon.pbpk.topology import load_species_topology  # noqa: E402
from charon.translational.decomposition import (  # noqa: E402
    decompose_fold_error,
    to_symmetric,
)
from validation.benchmarks.report_writer import emit_report  # noqa: E402

DEFAULT_PANEL = REPO_ROOT / "validation" / "data" / "fih_reference" / "panel.yaml"
DEFAULT_BIOAVAILABILITY = (
    REPO_ROOT / "validation" / "data" / "fih_reference" / "bioavailability.csv"
)
COMPOUNDS_DIR = REPO_ROOT / "validation" / "data" / "tier1_obach" / "compounds"
DEFAULT_STEM = REPO_ROOT / "validation" / "reports" / "layer3_ivive_decomposition"


def _load_compound(name: str) -> CompoundConfig:
    path = COMPOUNDS_DIR / f"{name}.yaml"
    if not path.exists():
        raise FileNotFoundError(f"compound YAML not found: {path}")
    return CompoundConfig.model_validate(yaml.safe_load(path.read_text()))


def _compute_baseline_and_whatif(entry: dict) -> dict[str, float]:
    """Run Pipeline once (well_stirred baseline) then analytically scale MRSD
    for parallel_tube and dispersion what-ifs.

    Rationale: Charon's PBPK ODE embeds well-stirred extraction directly
    (CLint_liver * fu_b * C_liver_blood_out), so Pipeline(liver_model=X) is
    a no-op for the simulation. MRSD is APPROXIMATELY proportional to CL_total
    via the PAD path; the ODE-derived cl_apparent is not algebraically identical
    to CLh_ws + CLrenal for high-Vss compounds (multi-compartment distributional
    effects). Alternate-model MRSDs are scaled as:

        mrsd_model = mrsd_ws * (CLh_model + cl_renal) / (CLh_ws + cl_renal)

    The approximation is sub-percent for most Tier A compounds and acceptable
    for Sprint 10's research diagnostic. Note the clint_liver_L_h passed to
    liver-model functions is the whole-liver IVIVE-scaled value (L/h), NOT
    the in-vitro microsomal CLint (uL/min/mg).

    Returns a dict with mrsd_ws_mg, mrsd_pt_mg, mrsd_disp_mg,
    clh_ws_L_h, clh_pt_L_h, clh_disp_L_h, cl_renal_L_h.
    """
    compound = _load_compound(entry["name"])
    pipe = Pipeline(
        compound,
        route=entry["route"],
        dose_mg=1.0,
        liver_model="well_stirred",
        dose_projection=DoseProjectionConfig(
            target_ceff_nM=float(entry["target_ceff_nM"]),
            safety_factor=10.0,
            tau_h=24.0,
        ),
    )
    result = pipe.run()
    if result.dose_recommendation is None:
        raise RuntimeError(f"No dose recommendation for {entry['name']}")
    mrsd_ws = float(result.dose_recommendation.mrsd_mg)

    md = result.metadata
    clint_liver = float(md["clint_liver_L_h"])
    cl_renal = float(md["cl_renal_L_h"])
    fu_b = float(md["fu_b"])
    topology = load_species_topology("human")
    qh = topology.tissues["liver"].blood_flow_L_h

    clh_ws = well_stirred(qh=qh, fu_b=fu_b, clint_liver=clint_liver)
    clh_pt = parallel_tube(qh=qh, fu_b=fu_b, clint_liver=clint_liver)
    clh_disp = dispersion(qh=qh, fu_b=fu_b, clint_liver=clint_liver)

    cl_total_ws = clh_ws + cl_renal
    if cl_total_ws <= 0:
        raise RuntimeError(
            f"{entry['name']}: cl_total_ws <= 0 "
            f"(clh_ws={clh_ws}, cl_renal={cl_renal})"
        )

    mrsd_pt = mrsd_ws * (clh_pt + cl_renal) / cl_total_ws
    mrsd_disp = mrsd_ws * (clh_disp + cl_renal) / cl_total_ws

    return {
        "mrsd_ws_mg": mrsd_ws,
        "mrsd_pt_mg": mrsd_pt,
        "mrsd_disp_mg": mrsd_disp,
        "clh_ws_L_h": clh_ws,
        "clh_pt_L_h": clh_pt,
        "clh_disp_L_h": clh_disp,
        "cl_renal_L_h": cl_renal,
    }


def _load_bioavailability(path: Path) -> dict[str, dict]:
    with path.open() as fh:
        rows = list(csv.DictReader(fh))
    out: dict[str, dict] = {}
    for r in rows:
        f_raw = r["f_oral"].strip()
        out[r["compound"]] = {
            "fih_reference_route": r["fih_reference_route"],
            "f_oral": float(f_raw) if f_raw else None,
            "f_source": r["f_source"],
            "f_doi_or_pmid": r["f_doi_or_pmid"],
            "notes": r["notes"],
        }
    return out


def run_panel(panel_path: Path, bioav_path: Path) -> dict:
    panel = yaml.safe_load(panel_path.read_text())["panel"]
    bioav = _load_bioavailability(bioav_path)

    rows: list[dict] = []
    tier_a = [c for c in panel["compounds"] if c["tier"] == "gold"]
    for entry in tier_a:
        name = entry["name"]
        if name not in bioav:
            raise RuntimeError(
                f"{name} missing from bioavailability.csv — rerun Task 1 curation"
            )
        bundle = _compute_baseline_and_whatif(entry)
        mrsds = {
            "well_stirred": bundle["mrsd_ws_mg"],
            "parallel_tube": bundle["mrsd_pt_mg"],
            "dispersion": bundle["mrsd_disp_mg"],
        }
        bioav_row = bioav[name]
        result = decompose_fold_error(
            mrsd_ws=mrsds["well_stirred"],
            mrsd_pt=mrsds["parallel_tube"],
            mrsd_disp=mrsds["dispersion"],
            f_lit=bioav_row["f_oral"],
            route_ref=bioav_row["fih_reference_route"],
            fih_reference_mg=float(entry["reference_fih_mg"]),
        )
        rows.append({
            "compound": name,
            "mrsd_ws_mg": mrsds["well_stirred"],
            "mrsd_pt_mg": mrsds["parallel_tube"],
            "mrsd_disp_mg": mrsds["dispersion"],
            "clh_ws_L_h": bundle["clh_ws_L_h"],
            "clh_pt_L_h": bundle["clh_pt_L_h"],
            "clh_disp_L_h": bundle["clh_disp_L_h"],
            "cl_renal_L_h": bundle["cl_renal_L_h"],
            "reference_fih_mg": float(entry["reference_fih_mg"]),
            "fih_reference_route": bioav_row["fih_reference_route"],
            "f_lit": bioav_row["f_oral"],
            # Signed factors (log-additivity invariant holds on these).
            "fold_observed_signed": result.fold_observed_signed,
            "fold_liver_model_signed": result.fold_liver_model_signed,
            "fold_route_bias": result.fold_route_bias,
            "fold_residual_signed": result.fold_residual_signed,
            # Symmetric forms for human-readable reporting.
            "fold_observed": to_symmetric(result.fold_observed_signed),
            "fold_liver_model": to_symmetric(result.fold_liver_model_signed),
            "fold_residual": to_symmetric(result.fold_residual_signed),
            "best_alt_model": result.best_alt_model_name,
            "flags": ",".join(result.flags) if result.flags else "-",
            "f_source": bioav_row["f_source"],
            "notes": bioav_row["notes"],
        })

    # Aggregate attribution percentages (use symmetric forms — always >= 1,
    # so abs(log10) is equivalent to abs(log10(signed))).
    def _pct(key: str) -> float:
        num = sum(abs(math.log10(r[key])) for r in rows if r[key] > 0)
        denom = sum(
            abs(math.log10(r["fold_observed"])) for r in rows if r["fold_observed"] > 0
        )
        return (num / denom * 100.0) if denom > 0 else 0.0

    summary = {
        "n_compounds": len(rows),
        "aggregate_pct_liver_model": _pct("fold_liver_model"),
        "aggregate_pct_route_bias": _pct("fold_route_bias"),
        "aggregate_pct_residual": _pct("fold_residual"),
    }

    # Sort rows by residual descending (worst-unexplained first) for report.
    rows_sorted = sorted(rows, key=lambda r: -r["fold_residual"])

    return {
        "title": "Charon Sprint 10 — Tier A IVIVE Fold-Error Decomposition",
        "panel": panel["name"],
        "date_utc": datetime.now(timezone.utc).isoformat(),
        "summary": summary,
        "rows": [],
        "extra_sections": {
            "Per-compound decomposition": rows_sorted,
        },
        "notes": [
            "Decomposition: fold_observed = fold_liver_model * fold_route_bias * fold_residual.",
            "fold_liver_model: ws/best_alt if alternate improves prediction, else 1.0.",
            "fold_route_bias: 1/F for oral-reference compounds, 1.0 for IV or unknown.",
            "fold_residual: unexplained remainder (transporter, non-hepatic, UGT, model-gap).",
            "Sorted by fold_residual descending (worst-unexplained first).",
            "Aggregate %: 100 * sum(|log10(factor)|) / sum(|log10(fold_observed)|).",
            "Liver-model what-ifs are analytical (mrsd_model = mrsd_ws * CL_model/CL_ws).",
            "Pipeline.liver_model is a no-op for the PBPK ODE; only one Pipeline run per compound.",
            "Research only — no production code changes.",
        ],
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel", type=Path, default=DEFAULT_PANEL)
    ap.add_argument("--bioavailability", type=Path, default=DEFAULT_BIOAVAILABILITY)
    ap.add_argument("--output-stem", type=Path, default=DEFAULT_STEM)
    args = ap.parse_args(argv)

    try:
        payload = run_panel(args.panel, args.bioavailability)
    except RuntimeError as exc:
        print(f"[FAIL] {exc}", file=sys.stderr)
        return 1

    emit_report(payload, stem=args.output_stem)
    print(
        f"[OK] Decomposition wrote {args.output_stem}.md + .json "
        f"({payload['summary']['n_compounds']} compounds)"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
