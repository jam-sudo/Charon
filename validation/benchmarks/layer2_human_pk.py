"""Layer 2 human PBPK benchmark — Obach 1999 Tier-1 panel (10 compounds).

Loads ``validation/data/tier1_obach/panel.yaml``, runs each compound
twice (no-override and with-override), and prints per-metric panel
AAFE / within-2-fold / within-3-fold summaries. Exit 0 iff every
``strict_targets: true`` compound passes 2-fold on CL, Vss, and t½.

Run as a standalone script::

    python3 validation/benchmarks/layer2_human_pk.py
"""

from __future__ import annotations

import argparse
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
    mc_percentiles: dict[str, float] | None = None
    """Optional Monte-Carlo percentiles (cl_p05/p50/p95, vss_p05/p50/p95).

    Populated only when ``--propagate-ci`` is active AND at least one of
    fu_p / clint is ``source == "ml_ensemble"`` with a conformal CI. Left
    ``None`` in point-estimate mode so the default execution path and
    JSON shape are bit-identical to pre-Task-8 runs.
    """


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
# Monte-Carlo over conformal CIs (Sprint 7 Task 8)
# ---------------------------------------------------------------------------


def _is_ml_with_ci(prop) -> bool:
    """True if ``prop`` is an ML-sourced PredictedProperty with a 90% CI."""
    if prop is None:
        return False
    if getattr(prop, "source", None) != "ml_ensemble":
        return False
    return (
        prop.ci_90_lower is not None
        and prop.ci_90_upper is not None
        and prop.ci_90_lower > 0
        and prop.ci_90_upper > 0
    )


def _lognormal_samples(prop, n: int, rng) -> "list[float]":
    """Sample ``n`` log-normal draws whose 90% band matches ``prop``'s CI.

    Uses log10-space symmetric CI -> sigma_log = (log10(upper) -
    log10(lower)) / (2 * 1.645). Returns a numpy array in linear space.
    """
    log_lower = math.log10(prop.ci_90_lower)
    log_upper = math.log10(prop.ci_90_upper)
    mu = math.log10(prop.value)
    sigma_log = (log_upper - log_lower) / (2.0 * 1.645)
    return 10.0 ** rng.normal(loc=mu, scale=sigma_log, size=n)


def _perturb_properties(
    props: CompoundProperties,
    *,
    fup: float | None = None,
    clint: float | None = None,
) -> CompoundProperties:
    """Return a copy of ``props`` with fu_p and/or CLint overridden.

    Preserves CI fields (they're not relevant for a single MC sample).
    Sampled values are not expected to retain the original CI so those
    fields are dropped in the copy.
    """
    updates: dict = {}
    if fup is not None and props.binding.fu_p is not None:
        # Bound physically to (0, 1].
        fup_clipped = max(min(fup, 1.0), 1e-6)
        new_fup = props.binding.fu_p.model_copy(
            update={
                "value": float(fup_clipped),
                "ci_90_lower": None,
                "ci_90_upper": None,
            }
        )
        new_binding = props.binding.model_copy(update={"fu_p": new_fup})
        updates["binding"] = new_binding
    if clint is not None and props.metabolism.clint_uL_min_mg is not None:
        clint_clipped = max(float(clint), 1e-6)
        new_clint = props.metabolism.clint_uL_min_mg.model_copy(
            update={
                "value": clint_clipped,
                "ci_90_lower": None,
                "ci_90_upper": None,
            }
        )
        new_metab = props.metabolism.model_copy(
            update={"clint_uL_min_mg": new_clint}
        )
        updates["metabolism"] = new_metab
    return props.model_copy(update=updates) if updates else props


def _mc_sample_pk(
    entry: PanelEntry,
    n: int,
    seed: int = 42,
) -> dict[str, float]:
    """Monte-Carlo sample CL/Vss for one panel entry.

    Samples fu_p and clint_uL_min_mg independently from their conformal
    CIs in log10 space (only those that are ``source="ml_ensemble"``
    with a CI). Re-runs the PBPK pipeline for each sample and reports
    p05/p50/p95 of CL and Vss. Returns ``{}`` if no ML-sourced property
    has a CI (the point-estimate-only case — caller should not emit
    percentile fields for such rows).
    """
    import numpy as np

    fup_prop = entry.compound.properties.binding.fu_p
    clint_prop = entry.compound.properties.metabolism.clint_uL_min_mg

    sample_fup = _is_ml_with_ci(fup_prop)
    sample_clint = _is_ml_with_ci(clint_prop)
    if not sample_fup and not sample_clint:
        return {}

    rng = np.random.default_rng(seed)
    samples_fup = _lognormal_samples(fup_prop, n, rng) if sample_fup else None
    samples_clint = (
        _lognormal_samples(clint_prop, n, rng) if sample_clint else None
    )

    # Additional physical clip on fu_p: well-stirred requires
    # fu_b = fu_p/BP <= 1, so fu_p <= BP. Clip sampled fu_p before
    # passing through.
    bp_ratio = entry.compound.properties.binding.bp_ratio
    fup_upper = float(bp_ratio.value) if bp_ratio is not None else 1.0
    fup_upper = min(max(fup_upper, 0.01), 1.0)

    cls: list[float] = []
    vsss: list[float] = []
    for i in range(n):
        fup_i = None
        if samples_fup is not None:
            # Clip to [1e-6, BP_ratio] so fu_b stays in [0, 1].
            fup_i = max(min(float(samples_fup[i]), fup_upper), 1e-6)
        clint_i = (
            float(samples_clint[i]) if samples_clint is not None else None
        )
        props_i = _perturb_properties(
            entry.compound.properties,
            fup=fup_i,
            clint=clint_i,
        )
        compound_i = entry.compound.model_copy(update={"properties": props_i})
        pipe_i = Pipeline(
            compound=compound_i,
            route=entry.route,  # type: ignore[arg-type]
            dose_mg=entry.dose_mg,
            duration_h=entry.duration_h,
        )
        try:
            res = pipe_i.run()
        except Exception:
            # A pathological draw (e.g. numerical edge case) should
            # not abort the whole panel — skip and continue.
            continue
        pk = res.pk_parameters
        if pk.cl_apparent is None or not math.isfinite(pk.cl_apparent):
            continue
        if pk.vss is None or not math.isfinite(pk.vss):
            continue
        cls.append(float(pk.cl_apparent))
        vsss.append(float(pk.vss))

    if not cls or not vsss:
        return {}

    return {
        "cl_p05": float(np.percentile(cls, 5)),
        "cl_p50": float(np.percentile(cls, 50)),
        "cl_p95": float(np.percentile(cls, 95)),
        "vss_p05": float(np.percentile(vsss, 5)),
        "vss_p50": float(np.percentile(vsss, 50)),
        "vss_p95": float(np.percentile(vsss, 95)),
        "n_samples": len(cls),
    }


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
    *,
    propagate_ci: bool = False,
    n_samples: int = 100,
) -> tuple[dict[str, PanelSummary], dict[str, list[PanelRow]]]:
    """Execute the two-pass benchmark and return (summaries, rows).

    When ``propagate_ci`` is False (default) the execution path is
    bit-identical to pre-Task-8: no MC calls, no extra RNG state, no
    extra fields on ``PanelRow``.

    When ``propagate_ci`` is True, the with-override rows for any
    compound whose fu_p and/or clint is ``ml_ensemble``-sourced with a
    conformal CI are additionally annotated with ``mc_percentiles``
    carrying p05/p50/p95 of CL and Vss from ``n_samples`` draws.
    """
    panel = load_panel(panel_path)

    rows_no_override: list[PanelRow] = []
    rows_with_override: list[PanelRow] = []

    for entry in panel:
        stripped_props = _without_kp_overrides(entry.compound.properties)
        stripped = entry.compound.model_copy(update={"properties": stripped_props})

        rows_no_override.append(_run_one(stripped, entry, mode="no_override"))
        row_with = _run_one(entry.compound, entry, mode="with_override")
        if propagate_ci:
            mc = _mc_sample_pk(entry, n=n_samples)
            if mc:
                row_with.mc_percentiles = mc
        rows_with_override.append(row_with)

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
    *,
    propagate_ci: bool = False,
    n_samples: int = 0,
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
        row_dict = {
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
        if r.mc_percentiles:
            # Monte-Carlo CI band (Task 8). Only present for ML-sourced
            # rows when --propagate-ci was passed; absent otherwise so
            # point-mode reports remain bit-identical to pre-Task-8.
            row_dict.update(r.mc_percentiles)
        result_rows.append(row_dict)

    notes = [
        f"Panel: n={s.n} compounds (Obach 1999 Tier-1)",
        "Mode: R&R + empirical Kp overrides (with_override)",
        f"Strict-gate failures: {s.strict_failures}",
        "AAFE targets: CL < 2.5, Vss < 3.0 (ARCHITECTURE.md §8)",
        "BDF solver (scipy.integrate.solve_ivp, method='BDF')",
        "fu_p applied only via fu_b = fu_p/BP in liver model (no double-application)",
    ]
    if propagate_ci:
        n_ci_rows = sum(1 for r in rows["with_override"] if r.mc_percentiles)
        notes.append(
            f"Monte-Carlo over conformal CIs: N={n_samples} samples per "
            f"compound, {n_ci_rows}/{s.n} rows had ml_ensemble-sourced "
            f"fu_p/CLint with conformal CIs."
        )
        notes.append(
            "fu_p and CLint sampled log-normal from conformal 90% CI "
            "(sigma_log = (log10(upper) - log10(lower)) / (2 * 1.645))."
        )

    return {
        "title": "Charon Layer 2 Human PBPK Benchmark — Obach 1999 Tier-1 Panel",
        "panel": "obach_1999_tier1",
        "date_utc": date_utc,
        "summary": summary_list,
        "targets": targets,
        "rows": result_rows,
        "notes": notes,
    }


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Charon Layer 2 Human PBPK benchmark (Obach 1999 Tier-1).",
    )
    p.add_argument(
        "--propagate-ci",
        action="store_true",
        help=(
            "Monte-Carlo CL/Vss from conformal CIs (ml_ensemble-sourced "
            "properties only). Emits p05/p50/p95 per compound."
        ),
    )
    p.add_argument(
        "--n-samples",
        type=int,
        default=100,
        help="Number of MC samples per compound when --propagate-ci is set.",
    )
    p.add_argument(
        "--output-stem",
        type=Path,
        default=None,
        help=(
            "Path stem for report output; writes {stem}.md and {stem}.json. "
            "Defaults to validation/reports/layer2_human_pk."
        ),
    )
    return p.parse_args(argv)


def main(
    panel_path: Path | None = None,
    *,
    argv: list[str] | None = None,
) -> int:
    # When called programmatically with ``argv=None`` (the default),
    # default to an empty argv so argparse doesn't consume sys.argv
    # meant for the enclosing test runner / REPL. The CLI entrypoint
    # explicitly passes ``sys.argv[1:]``.
    if argv is None:
        argv = []
    args = _parse_args(argv)
    panel_path = panel_path or DEFAULT_PANEL_PATH
    summaries, rows = run_benchmark(
        panel_path,
        propagate_ci=args.propagate_ci,
        n_samples=args.n_samples,
    )

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

    if args.propagate_ci:
        mc_rows = [r for r in rows["with_override"] if r.mc_percentiles]
        print()
        print(
            f"Monte-Carlo (N={args.n_samples}) over conformal CIs: "
            f"{len(mc_rows)}/{len(rows['with_override'])} rows ML-sourced"
        )
        for r in mc_rows:
            mc = r.mc_percentiles or {}
            print(
                f"  {r.key:<14} "
                f"CL p05/p50/p95: {mc.get('cl_p05', 0):.3f}/"
                f"{mc.get('cl_p50', 0):.3f}/{mc.get('cl_p95', 0):.3f}  "
                f"Vss p05/p50/p95: {mc.get('vss_p05', 0):.2f}/"
                f"{mc.get('vss_p50', 0):.2f}/{mc.get('vss_p95', 0):.2f}"
            )

    print("=" * 100)
    strict_failures = summaries["with_override"].strict_failures
    if strict_failures == 0:
        print("All strict-targets compounds PASS — exit 0")
    else:
        print(f"{strict_failures} strict-targets compound(s) FAILED — exit 1")
    print("=" * 100)

    payload = _build_report_payload(
        summaries,
        rows,
        propagate_ci=args.propagate_ci,
        n_samples=args.n_samples,
    )
    stem = args.output_stem if args.output_stem else REPORTS_DIR / "layer2_human_pk"
    emit_report(payload, stem=stem)

    return 0 if strict_failures == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main(argv=sys.argv[1:]))
