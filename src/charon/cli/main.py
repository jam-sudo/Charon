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
    export_report,
    render_report,
)


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


def _format_value(v: float | int | None) -> str:
    if v is None:
        return "-"
    a = abs(float(v))
    if a == 0.0:
        return "0"
    if a >= 1000.0 or a < 0.01:
        return f"{float(v):.3e}"
    return f"{float(v):.4g}"


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
        sp.add_argument("-q", "--quiet", action="store_true", default=False)

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

    # Placeholder stubs for Task 10-11 so that `charon --help` already lists
    # all five subcommands. Task 10 will replace these with real handlers.
    for stub in ("translate", "recommend", "report"):
        sp_stub = sub.add_parser(
            stub, help=f"(Task 10-11) {stub} subcommand - not yet implemented"
        )
        sp_stub.add_argument("smiles", nargs="?")
        sp_stub.set_defaults(
            func=lambda args: _fail(f"{args.command} not yet implemented", code=2)
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
