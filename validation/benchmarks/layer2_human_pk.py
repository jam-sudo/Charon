"""Layer 2 human PBPK benchmark — Obach 1999 Tier-1 panel (10 compounds).

Loads ``validation/data/tier1_obach/panel.yaml``, runs each compound
twice (no-override and with-override), and prints per-metric panel
AAFE / within-2-fold / within-3-fold summaries. Exit 0 iff every
``strict_targets: true`` compound passes 2-fold on CL, Vss, and t½.

Run as a standalone script::

    python3 validation/benchmarks/layer2_human_pk.py
"""

from __future__ import annotations

import math
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Literal

import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

REPORTS_DIR = REPO_ROOT / "validation" / "reports"

from charon import Pipeline  # noqa: E402
from charon.core.schema import (  # noqa: E402
    CompoundConfig,
    CompoundProperties,
)
from validation.benchmarks.metrics import aafe, fold_error, within_n_fold  # noqa: E402
from validation.benchmarks.report_writer import emit_report  # noqa: E402


DEFAULT_PANEL_PATH = (
    REPO_ROOT / "validation" / "data" / "tier1_obach" / "panel.yaml"
)

_METRICS = ("cl_L_h", "vss_L", "t_half_h")


# ---------------------------------------------------------------------------
# Data types
# ---------------------------------------------------------------------------


@dataclass
class PanelEntry:
    key: str
    compound: CompoundConfig
    route: str
    dose_mg: float
    duration_h: float
    observed: dict[str, float]
    strict_targets: bool
    obach_table_row: int | None
    notes: str


@dataclass
class PanelRow:
    key: str
    predicted: dict[str, float]
    observed: dict[str, float]
    fold: dict[str, float]
    pass_2_fold: dict[str, bool]
    override_tissues: list[str]
    strict_targets: bool
    mode: str  # "no_override" | "with_override"


@dataclass
class PanelSummary:
    n: int
    mode: str
    aafe: dict[str, float]
    within_2_fold: dict[str, float]
    within_3_fold: dict[str, float]
    strict_failures: int


# ---------------------------------------------------------------------------
# Loader + strip helpers
# ---------------------------------------------------------------------------


def load_panel(panel_path: Path) -> list[PanelEntry]:
    """Load a panel.yaml file and all compound files it references."""
    panel_path = Path(panel_path)
    with panel_path.open() as f:
        raw = yaml.safe_load(f)

    default_duration = float(raw.get("default_duration_h", 168.0))
    entries: list[PanelEntry] = []

    for idx, item in enumerate(raw["compounds"]):
        compound_file = panel_path.parent / item["compound_file"]
        if not compound_file.exists():
            raise FileNotFoundError(
                f"panel.yaml entry #{idx} ({item.get('key', '?')}) "
                f"references missing compound file: {compound_file}"
            )
        with compound_file.open() as cf:
            compound_data = yaml.safe_load(cf)
        compound = CompoundConfig.model_validate(compound_data)

        observed = {
            "cl_L_h": float(item["observed"]["cl_L_h"]),
            "vss_L": float(item["observed"]["vss_L"]),
            "t_half_h": float(item["observed"]["t_half_h"]),
        }

        entries.append(
            PanelEntry(
                key=str(item["key"]),
                compound=compound,
                route=str(item["route"]),
                dose_mg=float(item["dose_mg"]),
                duration_h=float(item.get("duration_h", default_duration)),
                observed=observed,
                strict_targets=bool(item["strict_targets"]),
                obach_table_row=item.get("obach_table_row"),
                notes=str(item.get("notes", "")),
            )
        )

    return entries


def _without_kp_overrides(props: CompoundProperties) -> CompoundProperties:
    """Return a copy of props with empirical_kp_by_tissue cleared.

    Uses nested model_copy so other fields of DistributionProperties
    that may be added in future (e.g. Vss_pred, tissue-level fu) are
    preserved rather than silently reset to their defaults.
    """
    new_distribution = props.distribution.model_copy(
        update={"empirical_kp_by_tissue": None}
    )
    return props.model_copy(update={"distribution": new_distribution})


# ---------------------------------------------------------------------------
# Single-compound run
# ---------------------------------------------------------------------------


def _run_one(
    compound: CompoundConfig,
    entry: PanelEntry,
    mode: Literal["no_override", "with_override"],
) -> PanelRow:
    pipe = Pipeline(
        compound=compound,
        route=entry.route,  # type: ignore[arg-type]
        dose_mg=entry.dose_mg,
        duration_h=entry.duration_h,
    )
    result = pipe.run()
    pk = result.pk_parameters

    predicted = {
        "cl_L_h": float(pk.cl_apparent) if pk.cl_apparent is not None else float("nan"),
        "vss_L": float(pk.vss) if pk.vss is not None else float("nan"),
        "t_half_h": float(pk.half_life) if pk.half_life is not None else float("nan"),
    }
    fold: dict[str, float] = {}
    pass_2_fold: dict[str, bool] = {}
    for m in _METRICS:
        p = predicted[m]
        o = entry.observed[m]
        if p <= 0 or math.isnan(p):
            fold[m] = float("inf")
            pass_2_fold[m] = False
        else:
            fold[m] = fold_error(p, o)
            pass_2_fold[m] = fold[m] <= 2.0

    override_tissues = [
        ovr["tissue"] for ovr in result.metadata.get("kp_overrides", [])
    ]

    return PanelRow(
        key=entry.key,
        predicted=predicted,
        observed=entry.observed,
        fold=fold,
        pass_2_fold=pass_2_fold,
        override_tissues=override_tissues,
        strict_targets=entry.strict_targets,
        mode=mode,
    )


# ---------------------------------------------------------------------------
# Aggregation
# ---------------------------------------------------------------------------


def aggregate_summary(rows: list[PanelRow], mode: str) -> PanelSummary:
    n = len(rows)
    if n == 0:
        return PanelSummary(
            n=0,
            mode=mode,
            aafe={m: float("nan") for m in _METRICS},
            within_2_fold={m: float("nan") for m in _METRICS},
            within_3_fold={m: float("nan") for m in _METRICS},
            strict_failures=0,
        )

    aafe_by_metric: dict[str, float] = {}
    w2_by_metric: dict[str, float] = {}
    w3_by_metric: dict[str, float] = {}
    for metric in _METRICS:
        preds = [r.predicted[metric] for r in rows]
        obs = [r.observed[metric] for r in rows]
        aafe_by_metric[metric] = aafe(preds, obs)
        w2_by_metric[metric] = within_n_fold(preds, obs, n=2.0)
        w3_by_metric[metric] = within_n_fold(preds, obs, n=3.0)

    strict_failures = 0
    for r in rows:
        if not r.strict_targets:
            continue
        if not all(r.pass_2_fold[m] for m in _METRICS):
            strict_failures += 1

    return PanelSummary(
        n=n,
        mode=mode,
        aafe=aafe_by_metric,
        within_2_fold=w2_by_metric,
        within_3_fold=w3_by_metric,
        strict_failures=strict_failures,
    )


# ---------------------------------------------------------------------------
# Main execution
# ---------------------------------------------------------------------------


def run_benchmark(
    panel_path: Path,
) -> tuple[dict[str, PanelSummary], dict[str, list[PanelRow]]]:
    """Execute the two-pass benchmark and return (summaries, rows)."""
    panel = load_panel(panel_path)

    rows_no_override: list[PanelRow] = []
    rows_with_override: list[PanelRow] = []

    for entry in panel:
        stripped_props = _without_kp_overrides(entry.compound.properties)
        stripped = entry.compound.model_copy(update={"properties": stripped_props})

        rows_no_override.append(_run_one(stripped, entry, mode="no_override"))
        rows_with_override.append(_run_one(entry.compound, entry, mode="with_override"))

    summaries = {
        "no_override": aggregate_summary(rows_no_override, mode="no_override"),
        "with_override": aggregate_summary(rows_with_override, mode="with_override"),
    }
    rows = {
        "no_override": rows_no_override,
        "with_override": rows_with_override,
    }
    return summaries, rows


def _print_panel_table(rows: list[PanelRow], title: str) -> None:
    print()
    print(title)
    print("-" * 100)
    header = (
        f"{'compound':<14}"
        f"{'CL_pred':>10}{'CL_obs':>10}{'f_CL':>8}  |"
        f"{'Vss_pred':>10}{'Vss_obs':>10}{'f_Vss':>8}  |"
        f"{'t_pred':>10}{'t_obs':>10}{'f_t':>8}  "
        f"{'verdict':>8}"
    )
    print(header)
    print("-" * 100)
    for r in rows:
        verdict = "PASS" if all(r.pass_2_fold.values()) else "FAIL"
        if r.strict_targets:
            verdict = verdict + "*"
        print(
            f"{r.key:<14}"
            f"{r.predicted['cl_L_h']:>10.3f}{r.observed['cl_L_h']:>10.3f}"
            f"{r.fold['cl_L_h']:>8.2f}  |"
            f"{r.predicted['vss_L']:>10.2f}{r.observed['vss_L']:>10.2f}"
            f"{r.fold['vss_L']:>8.2f}  |"
            f"{r.predicted['t_half_h']:>10.2f}{r.observed['t_half_h']:>10.2f}"
            f"{r.fold['t_half_h']:>8.2f}  "
            f"{verdict:>8}"
        )


def _print_summary(summary: PanelSummary, title: str) -> None:
    print()
    print(f"{title}  (n={summary.n})")
    print(
        f"  AAFE         CL: {summary.aafe['cl_L_h']:.2f}   "
        f"Vss: {summary.aafe['vss_L']:.2f}   "
        f"t_half: {summary.aafe['t_half_h']:.2f}"
    )
    print(
        f"  within 2x    CL: {summary.within_2_fold['cl_L_h']*100:.0f}%   "
        f"Vss: {summary.within_2_fold['vss_L']*100:.0f}%   "
        f"t_half: {summary.within_2_fold['t_half_h']*100:.0f}%"
    )
    print(
        f"  within 3x    CL: {summary.within_3_fold['cl_L_h']*100:.0f}%   "
        f"Vss: {summary.within_3_fold['vss_L']*100:.0f}%   "
        f"t_half: {summary.within_3_fold['t_half_h']*100:.0f}%"
    )
    print(f"  strict gate failures: {summary.strict_failures}")


def _build_report_payload(
    summaries: dict[str, PanelSummary],
    rows: dict[str, list[PanelRow]],
) -> dict:
    """Build a report_writer-compatible payload from "with_override" results."""
    s = summaries["with_override"]
    date_utc = datetime.now(tz=timezone.utc).isoformat()

    summary_list = [
        {
            "metric": "CL (L/h)",
            "AAFE": s.aafe["cl_L_h"],
            "within_2_fold": s.within_2_fold["cl_L_h"],
            "within_3_fold": s.within_3_fold["cl_L_h"],
        },
        {
            "metric": "Vss (L)",
            "AAFE": s.aafe["vss_L"],
            "within_2_fold": s.within_2_fold["vss_L"],
            "within_3_fold": s.within_3_fold["vss_L"],
        },
        {
            "metric": "t_half (h)",
            "AAFE": s.aafe["t_half_h"],
            "within_2_fold": s.within_2_fold["t_half_h"],
            "within_3_fold": s.within_3_fold["t_half_h"],
        },
    ]

    targets = [
        {
            "metric": "CL (L/h)",
            "target": "AAFE < 2.5",
            "met": s.aafe["cl_L_h"] < 2.5,
        },
        {
            "metric": "Vss (L)",
            "target": "AAFE < 3.0",
            "met": s.aafe["vss_L"] < 3.0,
        },
    ]

    result_rows = []
    for r in rows["with_override"]:
        result_rows.append(
            {
                "compound": r.key,
                "cl_pred": r.predicted["cl_L_h"],
                "cl_obs": r.observed["cl_L_h"],
                "cl_fold": r.fold["cl_L_h"],
                "cl_pass_2x": r.pass_2_fold["cl_L_h"],
                "vss_pred": r.predicted["vss_L"],
                "vss_obs": r.observed["vss_L"],
                "vss_fold": r.fold["vss_L"],
                "vss_pass_2x": r.pass_2_fold["vss_L"],
                "t_half_pred": r.predicted["t_half_h"],
                "t_half_obs": r.observed["t_half_h"],
                "t_half_fold": r.fold["t_half_h"],
                "t_half_pass_2x": r.pass_2_fold["t_half_h"],
                "strict_target": r.strict_targets,
            }
        )

    notes = [
        f"Panel: n={s.n} compounds (Obach 1999 Tier-1)",
        "Mode: R&R + empirical Kp overrides (with_override)",
        f"Strict-gate failures: {s.strict_failures}",
        "AAFE targets: CL < 2.5, Vss < 3.0 (ARCHITECTURE.md §8)",
        "BDF solver (scipy.integrate.solve_ivp, method='BDF')",
        "fu_p applied only via fu_b = fu_p/BP in liver model (no double-application)",
    ]

    return {
        "title": "Charon Layer 2 Human PBPK Benchmark — Obach 1999 Tier-1 Panel",
        "panel": "obach_1999_tier1",
        "date_utc": date_utc,
        "summary": summary_list,
        "targets": targets,
        "rows": result_rows,
        "notes": notes,
    }


def main(panel_path: Path | None = None) -> int:
    panel_path = panel_path or DEFAULT_PANEL_PATH
    summaries, rows = run_benchmark(panel_path)

    print("=" * 100)
    print("Charon Layer 2 Human PBPK Benchmark — Obach 1999 Tier-1 Panel")
    print("=" * 100)

    _print_panel_table(rows["no_override"], "R&R only (no empirical overrides)")
    _print_summary(summaries["no_override"], "Panel summary — R&R only")

    _print_panel_table(rows["with_override"], "R&R + empirical overrides")
    _print_summary(summaries["with_override"], "Panel summary — with overrides")

    # Report applied overrides
    applied = [r for r in rows["with_override"] if r.override_tissues]
    if applied:
        print()
        print("Overrides applied:")
        for r in applied:
            for tissue in r.override_tissues:
                print(f"  {r.key}: {tissue}")

    print("=" * 100)
    strict_failures = summaries["with_override"].strict_failures
    if strict_failures == 0:
        print("All strict-targets compounds PASS — exit 0")
    else:
        print(f"{strict_failures} strict-targets compound(s) FAILED — exit 1")
    print("=" * 100)

    payload = _build_report_payload(summaries, rows)
    emit_report(payload, stem=REPORTS_DIR / "layer2_human_pk")

    return 0 if strict_failures == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
