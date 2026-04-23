"""Layer 1 ADMET benchmark — conformal calibration set (n=153).

Reads ``data/validation/adme_reference.csv``, calls
``charon.predict.predict_properties`` for each compound, and computes
per-property metrics (logP: linear; fu_p / bp_ratio / peff: log-space).
Results are written via ``validation.benchmarks.report_writer.emit_report``.

Properties evaluated
--------------------
logP        csv:logP      path:physicochemical.logp     metrics: MAE, RMSE, R²
fu_p        csv:fup       path:binding.fu_p             metrics: AAFE, ≤2-fold, ≤3-fold
bp_ratio    csv:rbp       path:binding.bp_ratio         metrics: AAFE, ≤2-fold, ≤3-fold
peff_cm_s   csv:peff_cm_s path:permeability.peff_cm_s  metrics: AAFE, ≤2-fold, ≤3-fold

CLint EXCLUDED: unit mismatch between Charon (hepatocyte µL/min/10^6)
and reference (recombinant µL/min/pmol CYP3A4).

Usage
-----
python3 validation/benchmarks/layer1_admet.py
python3 validation/benchmarks/layer1_admet.py --csv path/to/adme_reference.csv
python3 validation/benchmarks/layer1_admet.py --limit 20   # quick dev run
"""

from __future__ import annotations

import argparse
import csv
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from charon.predict import predict_properties  # noqa: E402
from charon.core.schema import CompoundProperties  # noqa: E402
from validation.benchmarks.metrics import (  # noqa: E402
    aafe,
    mae,
    pearson_r,
    rmse,
    within_abs_diff,
    within_n_fold,
)
from validation.benchmarks.report_writer import emit_report  # noqa: E402

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_DEFAULT_CSV = REPO_ROOT / "data" / "validation" / "adme_reference.csv"

_DISCLAIMER = (
    "Evaluated on the conformal calibration set (n=153). Point-estimate "
    "metrics are meaningful; 90% CI coverage is tautological by construction "
    "(P90 was fit on the same set) and does NOT reflect true out-of-sample "
    "performance. CLint excluded: unit mismatch between Charon (hepatocyte "
    "uL/min/10^6) and reference (recombinant uL/min/pmol CYP3A4)."
)

# Minimum positive guards for log-space metrics. Values at or below this
# threshold are excluded from that property's metric computation.
_MIN_POSITIVE: dict[str, float | None] = {
    "logp": None,       # linear metric — no positivity required
    "fup": 1e-4,
    "rbp": 0.1,
    "peff_cm_s": 1e-8,
}

# ---------------------------------------------------------------------------
# Attribute path resolution
# ---------------------------------------------------------------------------

# Each entry: (csv_column, attribute_path_tuple)
_PROPERTY_MAP: list[tuple[str, tuple[str, ...]]] = [
    ("logP",     ("physicochemical", "logp")),
    ("fup",      ("binding", "fu_p")),
    ("rbp",      ("binding", "bp_ratio")),
    ("peff_cm_s", ("permeability", "peff_cm_s")),
]


def _get_predicted_value(
    props: CompoundProperties,
    path: tuple[str, ...],
) -> float | None:
    """Walk the attribute path on ``props`` and return ``obj.value``.

    Returns ``None`` if any intermediate attribute is ``None`` or the
    ``value`` field itself is ``None``.
    """
    obj: Any = props
    for attr in path:
        try:
            obj = getattr(obj, attr)
        except AttributeError:
            return None
        if obj is None:
            return None
    # obj is now the PredictedProperty (or may be a bare float for some paths)
    if hasattr(obj, "value"):
        return obj.value if obj.value is not None else None
    # If the leaf is already a float (edge case)
    if isinstance(obj, (int, float)):
        return float(obj)
    return None


# ---------------------------------------------------------------------------
# Core benchmark runner
# ---------------------------------------------------------------------------


def run_benchmark(
    *,
    csv_path: Path | str = _DEFAULT_CSV,
    reports_dir: Path | str | None = None,
    limit: int | None = None,
) -> dict:
    """Run the Layer 1 ADMET benchmark.

    Parameters
    ----------
    csv_path:
        Path to ``adme_reference.csv``.
    reports_dir:
        Directory for output reports. If *None*, defaults to
        ``<REPO_ROOT>/validation/reports``.
    limit:
        If set, process only the first *limit* rows (useful for rapid
        development testing).

    Returns
    -------
    dict
        Payload dict passed to ``emit_report`` (and also returned).
    """
    csv_path = Path(csv_path)
    if not csv_path.exists():
        raise FileNotFoundError(f"adme_reference.csv not found: {csv_path}")

    if reports_dir is None:
        reports_dir = REPO_ROOT / "validation" / "reports"
    reports_dir = Path(reports_dir)

    # Read CSV
    with csv_path.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        all_rows = list(reader)

    if limit is not None:
        all_rows = all_rows[:limit]

    # Accumulators per property: lists of (predicted, observed) pairs
    logp_pairs: list[tuple[float, float]] = []
    fup_pairs: list[tuple[float, float]] = []
    rbp_pairs: list[tuple[float, float]] = []
    peff_pairs: list[tuple[float, float]] = []

    result_rows: list[dict] = []
    excluded: list[dict] = []

    for row in all_rows:
        name = row.get("name", "unknown")
        smiles = row.get("smiles", "").strip()

        # --- Predict ---
        try:
            props = predict_properties(smiles)
        except ValueError:
            excluded.append({"name": name, "reason": "invalid SMILES"})
            continue
        except Exception as exc:  # noqa: BLE001
            excluded.append({"name": name, "reason": f"prediction failed: {exc}"})
            continue

        # --- Per-property extraction ---
        row_entry: dict = {"name": name, "smiles": smiles}

        # logP
        obs_logp = _parse_float(row.get("logP"))
        pred_logp = _get_predicted_value(props, ("physicochemical", "logp"))
        if obs_logp is not None and pred_logp is not None:
            logp_pairs.append((pred_logp, obs_logp))
            row_entry["logP_obs"] = obs_logp
            row_entry["logP_pred"] = pred_logp
            row_entry["logP_err"] = pred_logp - obs_logp
        else:
            row_entry["logP_obs"] = None
            row_entry["logP_pred"] = None
            row_entry["logP_err"] = None

        # fu_p
        obs_fup = _parse_float(row.get("fup"))
        pred_fup = _get_predicted_value(props, ("binding", "fu_p"))
        min_fup = _MIN_POSITIVE["fup"]
        if (
            obs_fup is not None
            and pred_fup is not None
            and obs_fup > min_fup  # type: ignore[operator]
            and pred_fup > min_fup  # type: ignore[operator]
        ):
            fup_pairs.append((pred_fup, obs_fup))
            row_entry["fup_obs"] = obs_fup
            row_entry["fup_pred"] = pred_fup
        else:
            row_entry["fup_obs"] = obs_fup
            row_entry["fup_pred"] = pred_fup

        # bp_ratio
        obs_rbp = _parse_float(row.get("rbp"))
        pred_rbp = _get_predicted_value(props, ("binding", "bp_ratio"))
        min_rbp = _MIN_POSITIVE["rbp"]
        if (
            obs_rbp is not None
            and pred_rbp is not None
            and obs_rbp > min_rbp  # type: ignore[operator]
            and pred_rbp > min_rbp  # type: ignore[operator]
        ):
            rbp_pairs.append((pred_rbp, obs_rbp))
            row_entry["rbp_obs"] = obs_rbp
            row_entry["rbp_pred"] = pred_rbp
        else:
            row_entry["rbp_obs"] = obs_rbp
            row_entry["rbp_pred"] = pred_rbp

        # peff_cm_s
        obs_peff = _parse_float(row.get("peff_cm_s"))
        pred_peff = _get_predicted_value(props, ("permeability", "peff_cm_s"))
        min_peff = _MIN_POSITIVE["peff_cm_s"]
        if (
            obs_peff is not None
            and pred_peff is not None
            and obs_peff > min_peff  # type: ignore[operator]
            and pred_peff > min_peff  # type: ignore[operator]
        ):
            peff_pairs.append((pred_peff, obs_peff))
            row_entry["peff_obs"] = obs_peff
            row_entry["peff_pred"] = pred_peff
        else:
            row_entry["peff_obs"] = obs_peff
            row_entry["peff_pred"] = pred_peff

        result_rows.append(row_entry)

    # --- Compute metrics ---
    summary: list[dict] = []

    # logP — linear metrics
    if logp_pairs:
        p_logp = [x[0] for x in logp_pairs]
        o_logp = [x[1] for x in logp_pairs]
        mae_logp = mae(p_logp, o_logp)
        rmse_logp = rmse(p_logp, o_logp)
        r2_logp = pearson_r(p_logp, o_logp) ** 2
        within05_logp = within_abs_diff(p_logp, o_logp, threshold=0.5)
        summary.append({
            "property": "logP",
            "n": len(logp_pairs),
            "metric_type": "linear",
            "mae": mae_logp,
            "rmse": rmse_logp,
            "r2": r2_logp,
            "within_0.5_logP_pct": within05_logp * 100.0,
        })
        print(
            f"logP  (n={len(logp_pairs):3d}) | "
            f"MAE={mae_logp:.3f}  RMSE={rmse_logp:.3f}  "
            f"R²={r2_logp:.3f}  within±0.5: {within05_logp*100:.1f}%"
        )
    else:
        summary.append({"property": "logP", "n": 0, "note": "no data"})

    # fu_p — log metrics
    if fup_pairs:
        p_fup = [x[0] for x in fup_pairs]
        o_fup = [x[1] for x in fup_pairs]
        aafe_fup = aafe(p_fup, o_fup)
        w2_fup = within_n_fold(p_fup, o_fup, n=2.0)
        w3_fup = within_n_fold(p_fup, o_fup, n=3.0)
        summary.append({
            "property": "fu_p",
            "n": len(fup_pairs),
            "metric_type": "log",
            "aafe": aafe_fup,
            "within_2fold_pct": w2_fup * 100.0,
            "within_3fold_pct": w3_fup * 100.0,
        })
        print(
            f"fu_p  (n={len(fup_pairs):3d}) | "
            f"AAFE={aafe_fup:.2f}  "
            f"≤2-fold: {w2_fup*100:.1f}%  ≤3-fold: {w3_fup*100:.1f}%"
        )
    else:
        summary.append({"property": "fu_p", "n": 0, "note": "no data"})

    # bp_ratio — log metrics
    if rbp_pairs:
        p_rbp = [x[0] for x in rbp_pairs]
        o_rbp = [x[1] for x in rbp_pairs]
        aafe_rbp = aafe(p_rbp, o_rbp)
        w2_rbp = within_n_fold(p_rbp, o_rbp, n=2.0)
        w3_rbp = within_n_fold(p_rbp, o_rbp, n=3.0)
        summary.append({
            "property": "bp_ratio",
            "n": len(rbp_pairs),
            "metric_type": "log",
            "aafe": aafe_rbp,
            "within_2fold_pct": w2_rbp * 100.0,
            "within_3fold_pct": w3_rbp * 100.0,
        })
        print(
            f"rbp   (n={len(rbp_pairs):3d}) | "
            f"AAFE={aafe_rbp:.2f}  "
            f"≤2-fold: {w2_rbp*100:.1f}%  ≤3-fold: {w3_rbp*100:.1f}%"
        )
    else:
        summary.append({"property": "bp_ratio", "n": 0, "note": "no data"})

    # peff_cm_s — log metrics
    if peff_pairs:
        p_peff = [x[0] for x in peff_pairs]
        o_peff = [x[1] for x in peff_pairs]
        aafe_peff = aafe(p_peff, o_peff)
        w2_peff = within_n_fold(p_peff, o_peff, n=2.0)
        w3_peff = within_n_fold(p_peff, o_peff, n=3.0)
        summary.append({
            "property": "peff_cm_s",
            "n": len(peff_pairs),
            "metric_type": "log",
            "aafe": aafe_peff,
            "within_2fold_pct": w2_peff * 100.0,
            "within_3fold_pct": w3_peff * 100.0,
        })
        print(
            f"peff  (n={len(peff_pairs):3d}) | "
            f"AAFE={aafe_peff:.2f}  "
            f"≤2-fold: {w2_peff*100:.1f}%  ≤3-fold: {w3_peff*100:.1f}%"
        )
    else:
        summary.append({
            "property": "peff_cm_s",
            "n": 0,
            "note": "no predicted values available",
        })
        print("peff  (n=  0) | no predicted values (permeability module not active)")

    # --- Build payload ---
    targets = [
        {
            "property": "logP",
            "target": "MAE < 1.0",
            "met": (
                summary[0].get("mae", float("inf")) < 1.0
                if summary[0].get("n", 0) > 0
                else None
            ),
        },
        {
            "property": "fu_p",
            "target": "AAFE < 2.0",
            "met": (
                summary[1].get("aafe", float("inf")) < 2.0
                if summary[1].get("n", 0) > 0
                else None
            ),
        },
        {
            "property": "bp_ratio",
            "target": "AAFE < 2.0",
            "met": (
                summary[2].get("aafe", float("inf")) < 2.0
                if summary[2].get("n", 0) > 0
                else None
            ),
        },
        {
            "property": "peff_cm_s",
            "target": "AAFE < 2.0",
            "met": (
                summary[3].get("aafe", float("inf")) < 2.0
                if summary[3].get("n", 0) > 0
                else None
            ),
        },
    ]

    notes = [
        "CLint excluded: unit mismatch between Charon hepatocyte (µL/min/10^6) "
        "and reference recombinant CYP3A4 (µL/min/pmol).",
        "peff_cm_s: predicted only when permeability module is active (ACAT pipeline).",
        "fu_p min-positive guard: obs and pred both > 1e-4. "
        "bp_ratio guard: both > 0.1. peff guard: both > 1e-8.",
        "logP source: RDKit Crippen cLogP (derived). "
        "fu_p source: XGBoost ml_ensemble. "
        "bp_ratio source: empirical formula.",
    ]

    payload = {
        "title": "Charon Layer 1 ADMET Benchmark — adme_reference.csv (n=153)",
        "panel": "adme_reference_v1",
        "date_utc": datetime.now(timezone.utc).isoformat(),
        "disclaimer": _DISCLAIMER,
        "summary": summary,
        "rows": result_rows,
        "excluded": [{"name": e["name"], "reason": e["reason"]} for e in excluded],
        "targets": targets,
        "notes": notes,
    }

    # --- Conformal coverage section (per-property from default predictor) ---
    try:
        from charon.predict.conformal import ConformalPredictor

        cp = ConformalPredictor.load_default()
        coverage_rows = [
            {
                "property": prop,
                "n_samples": rpt.n_samples,
                "empirical_coverage": rpt.empirical_coverage,
                "quantile_log10": rpt.quantile_log10,
                "factor": rpt.factor,
                "median_fold_error": rpt.median_fold_error,
                "mean_fold_error": rpt.mean_fold_error,
                "warning": rpt.warning or "",
            }
            for prop, rpt in cp._reports.items()
        ]
        if coverage_rows:
            payload.setdefault("extra_sections", {})["Conformal coverage"] = (
                coverage_rows
            )
    except Exception as exc:  # pragma: no cover - defensive
        print(f"[warn] conformal coverage unavailable: {exc}")

    # --- Emit report ---
    stem = reports_dir / "layer1_admet"
    md_path, json_path = emit_report(payload, stem=stem)
    print(f"\nReport written: {md_path}")
    print(f"             : {json_path}")

    return payload


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------


def _parse_float(val: str | None) -> float | None:
    """Parse a CSV cell to float; return None for empty / non-numeric."""
    if val is None:
        return None
    s = val.strip()
    if not s:
        return None
    try:
        return float(s)
    except ValueError:
        return None


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Run the Layer 1 ADMET benchmark against adme_reference.csv "
            "and emit a Markdown + JSON report."
        )
    )
    parser.add_argument(
        "--csv",
        type=Path,
        default=_DEFAULT_CSV,
        help="Path to adme_reference.csv (default: data/validation/adme_reference.csv)",
    )
    parser.add_argument(
        "--reports-dir",
        type=Path,
        default=None,
        help="Output directory for reports (default: validation/reports/)",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Process only the first N rows (for fast dev testing)",
    )
    args = parser.parse_args()

    csv_path: Path = args.csv
    if not csv_path.exists():
        print(f"ERROR: CSV not found: {csv_path}", file=sys.stderr)
        return 2

    try:
        run_benchmark(
            csv_path=csv_path,
            reports_dir=args.reports_dir,
            limit=args.limit,
        )
    except FileNotFoundError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    except Exception as exc:  # noqa: BLE001
        print(f"UNEXPECTED ERROR: {exc}", file=sys.stderr)
        return 2

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
