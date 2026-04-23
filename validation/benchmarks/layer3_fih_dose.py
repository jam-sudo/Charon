"""Layer 3 FIH dose benchmark — two-tier (gold + sanity-floor).

See docs/superpowers/specs/2026-04-23-sprint7-conformal-integration-and-layer3-design.md §6.

Tier A ("gold", n=5): report-only fold-error vs FDA-sourced FIH /
approved-starting doses.

Tier B ("sanity_floor", n=12, full Obach panel): gated check —
``MRSD_pred <= approved_starting_dose_mg``. Any failure is a
regulatory-class defect and the script exits with code 2.

MRSD is computed via the PAD (Pharmacologically Active Dose) path
using ``target_ceff_nM`` from the panel entry with ``safety_factor=10``.
"""
from __future__ import annotations

import argparse
import math
import sys
from datetime import datetime, timezone
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from charon import Pipeline  # noqa: E402
from charon.core.schema import CompoundConfig, DoseProjectionConfig  # noqa: E402
from validation.benchmarks.report_writer import emit_report  # noqa: E402

DEFAULT_PANEL = REPO_ROOT / "validation" / "data" / "fih_reference" / "panel.yaml"
COMPOUNDS_DIR = REPO_ROOT / "validation" / "data" / "tier1_obach" / "compounds"
DEFAULT_STEM = REPO_ROOT / "validation" / "reports" / "layer3_fih_dose"


def fold_error(a: float, b: float) -> float:
    """Symmetric fold-error max(a/b, b/a). Non-positive inputs -> inf."""
    if a <= 0 or b <= 0:
        return math.inf
    return max(a / b, b / a)


def passes_floor(predicted_mg: float, approved_starting_mg: float) -> bool:
    """Sanity floor: predicted MRSD must not exceed approved starting dose."""
    return predicted_mg <= approved_starting_mg + 1e-9


def _load_compound(name: str) -> CompoundConfig:
    """Load a CompoundConfig from the Obach compound YAML directory."""
    path = COMPOUNDS_DIR / f"{name}.yaml"
    if not path.exists():
        raise FileNotFoundError(f"compound YAML not found: {path}")
    data = yaml.safe_load(path.read_text())
    return CompoundConfig.model_validate(data)


def _compute_mrsd(entry: dict) -> tuple[float, str]:
    """Return (MRSD in mg, limiting_method) via the PAD path.

    Uses ``target_ceff_nM`` from the panel entry — PAD translates target
    steady-state exposure into a daily dose, then applies safety_factor=10.
    """
    compound = _load_compound(entry["name"])
    pipe = Pipeline(
        compound,
        route=entry["route"],
        dose_mg=1.0,
        dose_projection=DoseProjectionConfig(
            target_ceff_nM=float(entry["target_ceff_nM"]),
            safety_factor=10.0,
            tau_h=24.0,
        ),
    )
    result = pipe.run()
    if result.dose_recommendation is None:
        raise RuntimeError(f"No dose recommendation for {entry['name']}")
    return (
        float(result.dose_recommendation.mrsd_mg),
        result.dose_recommendation.limiting_method,
    )


def run_panel(panel_path: Path) -> dict:
    panel = yaml.safe_load(panel_path.read_text())["panel"]
    gold_rows: list[dict] = []
    floor_rows: list[dict] = []
    for entry in panel["compounds"]:
        mrsd_mg, limiting = _compute_mrsd(entry)
        if entry["tier"] == "gold":
            ref = entry["reference_fih_mg"]
            fold = fold_error(mrsd_mg, ref)
            gold_rows.append({
                "compound": entry["name"],
                "route": entry["route"],
                "mrsd_pred_mg": mrsd_mg,
                "limiting_method": limiting,
                "reference_fih_mg": ref,
                "source_type": entry.get("source_type", "unknown"),
                "fold_error": fold,
                "within_3x": fold <= 3.0,
                "within_10x": fold <= 10.0,
                "source": entry["source"],
            })
        else:
            start = entry["approved_starting_dose_mg"]
            floor_rows.append({
                "compound": entry["name"],
                "route": entry["route"],
                "mrsd_pred_mg": mrsd_mg,
                "limiting_method": limiting,
                "approved_starting_dose_mg": start,
                "pass_floor": passes_floor(mrsd_mg, start),
                "source": entry["source"],
            })
    n_gold = len(gold_rows)
    n_floor = len(floor_rows)
    summary = {
        "gold_n": n_gold,
        "gold_within_3x": sum(r["within_3x"] for r in gold_rows),
        "gold_within_3x_fraction": (
            sum(r["within_3x"] for r in gold_rows) / n_gold if n_gold else 0.0
        ),
        "gold_within_10x": sum(r["within_10x"] for r in gold_rows),
        "sanity_n": n_floor,
        "sanity_pass_count": sum(r["pass_floor"] for r in floor_rows),
        "sanity_pass_fraction": (
            sum(r["pass_floor"] for r in floor_rows) / n_floor if n_floor else 0.0
        ),
        "sanity_failures": ",".join(
            r["compound"] for r in floor_rows if not r["pass_floor"]
        ) or "-",
    }
    return {
        "title": "Charon Layer 3 FIH Dose Benchmark",
        "panel": panel["name"],
        "date_utc": datetime.now(timezone.utc).isoformat(),
        "summary": summary,
        "rows": [],
        "extra_sections": {
            "Gold (Tier A) — fold-error vs reference FIH": gold_rows,
            "Sanity floor (Tier B) — MRSD <= approved starting dose": floor_rows,
        },
        "notes": [
            "Tier A (gold): fold-error vs published FIH or approved-start dose. Report-only.",
            "Tier B (sanity_floor): MRSD_pred <= approved_starting_dose_mg. GATED — any failure blocks release.",
            "MRSD computed via PAD path (target_ceff_nM) with safety_factor=10.",
            "Sources: inline per compound in panel.yaml.",
        ],
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel", type=Path, default=DEFAULT_PANEL)
    ap.add_argument(
        "--output-stem",
        type=Path,
        default=DEFAULT_STEM,
        help="Output path stem; writes {stem}.md and {stem}.json",
    )
    args = ap.parse_args(argv)

    payload = run_panel(args.panel)
    emit_report(payload, stem=args.output_stem)

    failures = payload["summary"]["sanity_failures"]
    if failures and failures != "-":
        print(
            f"[FAIL] Sanity floor: failures in {failures}",
            file=sys.stderr,
        )
        return 2
    print(
        f"[OK] Sanity floor: {payload['summary']['sanity_pass_count']}/"
        f"{payload['summary']['sanity_n']} pass."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
