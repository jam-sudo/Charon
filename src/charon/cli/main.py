"""Charon CLI — five subcommands sharing a single dispatcher.

Subcommands:

    predict    Layer 0 + Layer 1 only (ADME prediction)
    simulate   Full pipeline up to Layer 2 (PBPK)
    translate  Pipeline + Layer 3 (FIH dose projection)
    recommend  translate + optional Layer 4 (uncertainty quantification)
    report     recommend + write Markdown/JSON report files

Invocation:

    charon <subcommand> <SMILES> [options...]
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import asdict
from typing import Sequence

from charon.core.schema import DoseProjectionConfig, UncertaintyConfig
from charon.pipeline import Pipeline
from charon.predict import predict_properties
from charon.report import (
    ReportData,
    collect,
)
from charon.report.narrative import format_value as _format_value


# ---------------------------------------------------------------------------
# Small helpers (printers and error exits)
# ---------------------------------------------------------------------------


def _fail(msg: str, code: int = 1) -> int:
    print(f"error: {msg}", file=sys.stderr)
    return code


def _print_table(rows: list[list[str]], headers: list[str]) -> None:
    print("| " + " | ".join(headers) + " |")
    print("| " + " | ".join("---" for _ in headers) + " |")
    for row in rows:
        print("| " + " | ".join(str(c) for c in row) + " |")


# ---------------------------------------------------------------------------
# Subcommand: predict
# ---------------------------------------------------------------------------


def _cmd_predict(args: argparse.Namespace) -> int:
    try:
        props = predict_properties(args.smiles)
    except ValueError as e:
        return _fail(f"Invalid SMILES: {e}")
    except Exception as e:  # pragma: no cover - defensive
        return _fail(f"{type(e).__name__}: {e}", code=2)

    if args.json:
        print(json.dumps(props.model_dump(), indent=2, default=str))
        return 0

    rows: list[list[str]] = []
    for category_name in (
        "physicochemical",
        "binding",
        "metabolism",
        "permeability",
        "renal",
        "safety",
    ):
        cat = getattr(props, category_name, None)
        if cat is None:
            continue
        for attr_name in type(cat).model_fields:
            val = getattr(cat, attr_name, None)
            if val is None:
                continue
            if hasattr(val, "value"):
                rows.append(
                    [
                        attr_name,
                        _format_value(val.value),
                        val.unit or "-",
                        val.source or "-",
                    ]
                )
    if not rows:
        print("(no properties predicted)")
        return 0
    _print_table(rows, ["Property", "Value", "Unit", "Source"])
    return 0


# ---------------------------------------------------------------------------
# Subcommand: simulate
# ---------------------------------------------------------------------------


def _build_pipeline_from_smiles(args: argparse.Namespace, **extra) -> Pipeline:
    return Pipeline.from_smiles(
        args.smiles,
        route=args.route,
        dose_mg=args.dose,
        species=args.species,
        duration_h=args.duration,
        infusion_duration_h=getattr(args, "infusion_duration", 0.0),
        liver_model=args.liver_model,
        compound_name=args.compound_name or args.smiles,
        **extra,
    )


def _print_pk_summary(data: ReportData) -> None:
    rows = []
    labels = [
        ("cmax", "Cmax", "ug/L"),
        ("tmax", "Tmax", "h"),
        ("auc_0_inf", "AUC(0-inf)", "ug*h/L"),
        ("half_life", "t1/2", "h"),
        ("cl_apparent", "CL apparent", "L/h"),
        ("vss", "Vss", "L"),
        ("bioavailability", "F", "-"),
        ("fa", "Fa", "-"),
        ("fg", "Fg", "-"),
        ("fh", "Fh", "-"),
    ]
    for key, label, unit in labels:
        val = data.pk_params.get(key)
        rows.append([label, _format_value(val), unit])
    _print_table(rows, ["Parameter", "Value", "Unit"])


def _cmd_simulate(args: argparse.Namespace) -> int:
    try:
        pipe = _build_pipeline_from_smiles(args)
        result = pipe.run()
    except ValueError as e:
        return _fail(f"Invalid input: {e}")
    except Exception as e:  # pragma: no cover - defensive
        return _fail(f"{type(e).__name__}: {e}", code=2)

    data = collect(result)
    if args.json:
        print(json.dumps(asdict(data), indent=2, default=str))
        return 0
    print(f"Pipeline: {args.smiles}  |  route={args.route}  |  dose={args.dose} mg")
    print()
    _print_pk_summary(data)
    return 0


# ---------------------------------------------------------------------------
# Subcommand: translate
# ---------------------------------------------------------------------------


def _build_dose_projection(args: argparse.Namespace) -> DoseProjectionConfig | None:
    has_hed = args.noael is not None and args.noael_species is not None
    has_mabel = args.target_kd is not None
    has_pad = args.target_ceff is not None
    if not (has_hed or has_mabel or has_pad):
        return None
    return DoseProjectionConfig(
        noael_mg_kg=args.noael,
        noael_species=args.noael_species,
        target_kd_nM=args.target_kd,
        target_ceff_nM=args.target_ceff,
        safety_factor=args.safety_factor,
        tau_h=args.tau,
        body_weight_kg=args.body_weight,
    )


def _print_dose_recommendation(data: ReportData) -> None:
    rec = data.dose_recommendation
    if rec is None:
        print("(no dose recommendation produced)")
        return
    print(rec.get("rationale", ""))
    print()
    rows: list[list[str]] = []
    for key, label in (("hed", "HED"), ("mabel", "MABEL"), ("pad", "PAD")):
        sub = rec.get(key)
        if sub is None:
            rows.append([label, "-", "insufficient inputs"])
        else:
            mark = " <-" if rec.get("limiting_method") == key else ""
            rows.append([label, _format_value(sub.get("mrsd_mg")), f"computed{mark}"])
    _print_table(rows, ["Method", "MRSD (mg)", "Status"])


def _cmd_translate(args: argparse.Namespace) -> int:
    dp = _build_dose_projection(args)
    if dp is None:
        return _fail(
            "translate requires at least one target: "
            "--noael + --noael-species, --target-kd, or --target-ceff"
        )
    try:
        pipe = _build_pipeline_from_smiles(args, dose_projection=dp)
        result = pipe.run()
    except ValueError as e:
        return _fail(f"Invalid input: {e}")
    except Exception as e:  # pragma: no cover
        return _fail(f"{type(e).__name__}: {e}", code=2)

    data = collect(result)
    if args.json:
        print(json.dumps(asdict(data), indent=2, default=str))
        return 0
    _print_dose_recommendation(data)
    return 0


# ---------------------------------------------------------------------------
# Subcommand: recommend
# ---------------------------------------------------------------------------


def _print_uncertainty_summary(data: ReportData) -> None:
    unc = data.uncertainty
    if unc is None:
        return
    lo = _format_value(unc.get("ci_90_lower_mg"))
    hi = _format_value(unc.get("ci_90_upper_mg"))
    point = _format_value(unc.get("point_estimate_mg"))
    conf = unc.get("confidence", "?")
    print()
    print(
        f"Uncertainty: {point} mg [{lo} – {hi}] 90% CI  ·  confidence: {conf}"
    )
    sens = unc.get("sensitivity") or {}
    if sens:
        top = max(sens, key=sens.get)
        print(f"Top sensitivity: {top} ({sens[top] * 100:.1f}%)")
    rec_txt = unc.get("recommendation") or ""
    if rec_txt:
        print(rec_txt)


def _cmd_recommend(args: argparse.Namespace) -> int:
    dp = _build_dose_projection(args)
    if dp is None:
        return _fail(
            "recommend requires at least one target: "
            "--noael + --noael-species, --target-kd, or --target-ceff"
        )
    unc_cfg: UncertaintyConfig | None = None
    if args.uncertainty:
        unc_cfg = UncertaintyConfig(n_samples=args.n_samples)
    try:
        pipe = _build_pipeline_from_smiles(
            args, dose_projection=dp, uncertainty=unc_cfg
        )
        result = pipe.run()
    except ValueError as e:
        return _fail(f"Invalid input: {e}")
    except Exception as e:  # pragma: no cover
        return _fail(f"{type(e).__name__}: {e}", code=2)

    data = collect(result)
    if args.json:
        print(json.dumps(asdict(data), indent=2, default=str))
        return 0
    _print_dose_recommendation(data)
    _print_uncertainty_summary(data)
    return 0


# ---------------------------------------------------------------------------
# Argument parser
# ---------------------------------------------------------------------------


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="charon",
        description="Charon: SMILES -> FIH dose recommendation pipeline.",
    )
    sub = p.add_subparsers(dest="command", metavar="COMMAND")

    # Shared option helper
    def _add_common(sp: argparse.ArgumentParser) -> None:
        sp.add_argument(
            "--species", default="human", choices=["human", "rat", "dog", "monkey"]
        )
        sp.add_argument(
            "--liver-model",
            default="well_stirred",
            choices=["well_stirred", "parallel_tube", "dispersion"],
            dest="liver_model",
        )
        sp.add_argument("--compound-name", default=None, dest="compound_name")
        sp.add_argument("--json", action="store_true", default=False)

    # Shared dose-projection options
    def _add_dose_opts(sp: argparse.ArgumentParser) -> None:
        sp.add_argument("--noael", type=float, default=None)
        sp.add_argument("--noael-species", default=None, dest="noael_species")
        sp.add_argument("--target-kd", type=float, default=None, dest="target_kd")
        sp.add_argument(
            "--target-ceff", type=float, default=None, dest="target_ceff"
        )
        sp.add_argument(
            "--safety-factor", type=float, default=10.0, dest="safety_factor"
        )
        sp.add_argument("--tau", type=float, default=24.0)
        sp.add_argument(
            "--body-weight", type=float, default=70.0, dest="body_weight"
        )

    # predict
    sp_pred = sub.add_parser("predict", help="ADME prediction (Layer 0+1 only)")
    sp_pred.add_argument("smiles")
    _add_common(sp_pred)
    sp_pred.set_defaults(func=_cmd_predict)

    # simulate
    sp_sim = sub.add_parser(
        "simulate", help="PK simulation (Layer 0+1+2, no dose projection)"
    )
    sp_sim.add_argument("smiles")
    sp_sim.add_argument(
        "--route",
        required=True,
        choices=["iv_bolus", "iv_infusion", "oral"],
    )
    sp_sim.add_argument("--dose", type=float, required=True, help="dose in mg")
    sp_sim.add_argument("--duration", type=float, default=72.0)
    sp_sim.add_argument(
        "--infusion-duration", type=float, default=0.0, dest="infusion_duration"
    )
    _add_common(sp_sim)
    sp_sim.set_defaults(func=_cmd_simulate)

    # translate
    sp_tr = sub.add_parser(
        "translate", help="Pipeline + FIH dose projection (Layers 0-3)"
    )
    sp_tr.add_argument("smiles")
    sp_tr.add_argument(
        "--route", required=True, choices=["iv_bolus", "iv_infusion", "oral"]
    )
    sp_tr.add_argument("--dose", type=float, required=True, help="dose in mg")
    sp_tr.add_argument("--duration", type=float, default=72.0)
    sp_tr.add_argument(
        "--infusion-duration", type=float, default=0.0, dest="infusion_duration"
    )
    _add_dose_opts(sp_tr)
    _add_common(sp_tr)
    sp_tr.set_defaults(func=_cmd_translate)

    # recommend
    sp_rc = sub.add_parser(
        "recommend",
        help="Pipeline + dose projection + optional uncertainty (Layers 0-4)",
    )
    sp_rc.add_argument("smiles")
    sp_rc.add_argument(
        "--route", required=True, choices=["iv_bolus", "iv_infusion", "oral"]
    )
    sp_rc.add_argument("--dose", type=float, required=True, help="dose in mg")
    sp_rc.add_argument("--duration", type=float, default=72.0)
    sp_rc.add_argument(
        "--infusion-duration", type=float, default=0.0, dest="infusion_duration"
    )
    _add_dose_opts(sp_rc)
    sp_rc.add_argument("--uncertainty", action="store_true", default=False)
    sp_rc.add_argument("--n-samples", type=int, default=500, dest="n_samples")
    _add_common(sp_rc)
    sp_rc.set_defaults(func=_cmd_recommend)

    # report (Task 11 stub — remove this when Task 11 registers the real handler)
    sp_rp = sub.add_parser(
        "report",
        help="(Task 11) Full pipeline + write report files — not yet implemented",
    )
    sp_rp.add_argument("smiles", nargs="?")
    sp_rp.set_defaults(
        func=lambda args: _fail("report not yet implemented", code=2)
    )

    return p


def main(argv: Sequence[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    if not getattr(args, "command", None):
        parser.print_help(sys.stderr)
        return 1
    func = args.func
    return func(args)


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
