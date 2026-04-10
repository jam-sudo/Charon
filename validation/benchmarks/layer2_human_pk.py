"""Layer 2 human PBPK benchmark.

Runs the Charon IV PBPK kernel on a small reference panel and prints a
table comparing predicted CL / Vss / t_half to literature observed values.

Sprint 3 scope: 2 compounds
  - theophylline  : well-behaved R&R Kp (neutral, low logP) — primary validation
  - midazolam     : documented R&R overprediction case (weak base, high logP)

Sprint 3b will extend to the full Obach 1999 IV dataset and tighten the
R&R Kp calibration for weak bases via Berezhkovskiy and/or empirical
adipose overrides.

Run as a standalone script::

    python3 validation/benchmarks/layer2_human_pk.py

Exit code: 0 if all strict-target rows PASS, 1 otherwise.
"""

from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from charon import Pipeline  # noqa: E402
from charon.core.schema import (  # noqa: E402
    BindingProperties,
    CompoundConfig,
    CompoundProperties,
    MetabolismProperties,
    PhysicochemicalProperties,
    PredictedProperty,
    RenalProperties,
)
from validation.benchmarks.metrics import fold_error  # noqa: E402


def _p(value: float, unit: str = "") -> PredictedProperty:
    return PredictedProperty(value=float(value), source="experimental", unit=unit)


# ---------------------------------------------------------------------------
# Sprint 3b: Panel loader — PanelEntry dataclass and load_panel()
# The Sprint 3a BenchmarkCompound / theophylline() / midazolam() factories
# below remain unchanged and will be removed in Task 14.
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


def load_panel(panel_path: Path) -> list[PanelEntry]:
    """Load a panel.yaml file and all compound files it references.

    `panel.yaml` format::

        name: "..."
        source: "..."
        default_duration_h: 168.0
        compounds:
          - key: theophylline
            compound_file: compounds/theophylline.yaml
            route: iv_bolus
            dose_mg: 100.0
            duration_h: 168.0          # optional; falls back to default_duration_h
            observed:
              cl_L_h: 2.9
              vss_L: 35.0
              t_half_h: 8.0
            obach_table_row: 42        # optional
            strict_targets: true
            notes: "..."               # optional

    Returns a list of `PanelEntry` objects with fully-loaded
    `CompoundConfig` instances.
    """
    panel_path = Path(panel_path)
    with panel_path.open() as f:
        raw = yaml.safe_load(f)

    default_duration = float(raw.get("default_duration_h", 168.0))
    entries: list[PanelEntry] = []

    for idx, item in enumerate(raw["compounds"]):
        # Resolve compound file relative to the panel's directory
        compound_file = panel_path.parent / item["compound_file"]
        if not compound_file.exists():
            raise FileNotFoundError(
                f"panel.yaml entry #{idx} ({item.get('key', '?')}) "
                f"references missing compound file: {compound_file}"
            )
        with compound_file.open() as cf:
            compound_data = yaml.safe_load(cf)
        compound = CompoundConfig.model_validate(compound_data)

        # Observed PK keys are required (KeyError surfaces typos early)
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


# ---------------------------------------------------------------------------
# Sprint 3a: BenchmarkCompound Python factories (removed in Task 14)
# ---------------------------------------------------------------------------

@dataclass
class BenchmarkCompound:
    compound: CompoundConfig
    dose_mg: float
    duration_h: float
    compound_type_override: str | None
    observed: dict[str, float]       # cl, vss, t_half
    strict_targets: bool             # if False, all rows print but do not gate exit


def theophylline() -> BenchmarkCompound:
    cfg = CompoundConfig(
        name="theophylline",
        smiles="Cn1c(=O)c2[nH]cnc2n(C)c1=O",
        molecular_weight=180.17,
        source="experimental",
        properties=CompoundProperties(
            physicochemical=PhysicochemicalProperties(logp=_p(-0.02)),
            binding=BindingProperties(
                fu_p=_p(0.60, "fraction"),
                fu_inc=_p(1.0, "fraction"),
                bp_ratio=_p(0.85, "ratio"),
            ),
            metabolism=MetabolismProperties(clint_uL_min_mg=_p(1.8, "uL/min/mg")),
            renal=RenalProperties(clrenal_L_h=_p(0.1, "L/h")),
        ),
    )
    return BenchmarkCompound(
        compound=cfg,
        dose_mg=100.0,
        duration_h=168.0,
        compound_type_override=None,
        observed={"cl": 2.9, "vss": 35.0, "t_half": 8.0},
        strict_targets=True,
    )


def midazolam() -> BenchmarkCompound:
    cfg = CompoundConfig(
        name="midazolam",
        smiles="Cc1ncc2n1-c1ccc(Cl)cc1C(c1ccccc1F)=NC2",
        molecular_weight=325.77,
        source="experimental",
        properties=CompoundProperties(
            physicochemical=PhysicochemicalProperties(
                logp=_p(3.89),
                pka_base=_p(6.2),
            ),
            binding=BindingProperties(
                fu_p=_p(0.03, "fraction"),
                fu_inc=_p(0.96, "fraction"),
                bp_ratio=_p(0.66, "ratio"),
            ),
            metabolism=MetabolismProperties(clint_uL_min_mg=_p(93.0, "uL/min/mg")),
            renal=RenalProperties(clrenal_L_h=_p(0.0, "L/h")),
        ),
    )
    return BenchmarkCompound(
        compound=cfg,
        dose_mg=5.0,
        duration_h=168.0,
        compound_type_override="base",
        observed={"cl": 21.0, "vss": 66.0, "t_half": 3.0},
        strict_targets=False,  # known R&R over-prediction for weak bases
    )


def _run_one(bc: BenchmarkCompound, target_fold: float) -> tuple[bool, dict]:
    pipe = Pipeline(
        compound=bc.compound,
        route="iv_bolus",
        dose_mg=bc.dose_mg,
        duration_h=bc.duration_h,
        compound_type_override=bc.compound_type_override,
    )
    result = pipe.run()
    pk = result.pk_parameters
    pred = {
        "cl": pk.cl_apparent,
        "vss": pk.vss,
        "t_half": pk.half_life,
    }

    print(f"\n{'=' * 72}")
    print(f"Compound: {bc.compound.name}  (dose={bc.dose_mg} mg IV bolus)")
    print(f"Target fold: <= {target_fold}x   Strict gate: {bc.strict_targets}")
    print("-" * 72)
    print(f"{'Metric':<12} {'Predicted':>12} {'Observed':>12} {'Fold Err':>10} {'Verdict':>10}")
    print("-" * 72)
    all_pass = True
    for metric in ("cl", "vss", "t_half"):
        p = pred[metric]
        o = bc.observed[metric]
        if p is None or p <= 0:
            verdict = "ERROR"
            fe = float("inf")
            all_pass = False
        else:
            fe = fold_error(p, o)
            verdict = "PASS" if fe <= target_fold else "FAIL"
            if verdict == "FAIL":
                all_pass = False
        print(f"{metric:<12} {p if p is not None else 0:>12.3f} {o:>12.3f} {fe:>10.3f} {verdict:>10}")
    print("-" * 72)
    print(f"PBPK params:")
    print(f"  compound_type    = {result.metadata['compound_type']}")
    print(f"  fu_b             = {result.metadata['fu_b']:.6f}")
    print(f"  CLint_liver_L_h  = {result.metadata['clint_liver_L_h']:.3f}")
    print(f"  CL_renal_L_h     = {result.metadata['cl_renal_L_h']:.3f}")
    print(f"  solver_method    = {result.metadata['solver_method']}")
    print(f"  solver_nfev      = {result.metadata['solver_nfev']}")

    return all_pass, pred


def main() -> int:
    target_fold = 2.0
    print("=" * 72)
    print("Charon Layer 2 Human PBPK Benchmark")
    print("=" * 72)

    panel = [theophylline(), midazolam()]

    gated_failures = 0
    for bc in panel:
        passed, _pred = _run_one(bc, target_fold=target_fold)
        if bc.strict_targets and not passed:
            gated_failures += 1

    print("\n" + "=" * 72)
    if gated_failures == 0:
        print("ALL STRICT-TARGET COMPOUNDS PASS (2-fold on CL, Vss, t_half)")
    else:
        print(f"{gated_failures} STRICT-TARGET COMPOUND(S) FAILED")
    print("=" * 72)

    return 0 if gated_failures == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
