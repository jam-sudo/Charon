"""Sprint 15 audit — propranolol F-decomposition + multiplier sweep.

Purpose:
1. Capture Fa/Fg/Fh/F_oral from current Charon pipeline for propranolol
   (no multiplier applied) — validate against analytical estimates.
2. Run an empirical multiplier sweep (m in {1, 2, 3, 5, 8, 12, 15, 20, 25, 30})
   to produce the fold-vs-m curve. Identify m_close_3x and m_close_1x.

Output: prints two markdown tables to stdout. Pipe to file or copy into the
audit subsection of validation/reports/layer3_ivive_decomposition.md.

Usage:
    python3 scripts/sprint15_audit.py > /tmp/sprint15_audit.md
"""

from __future__ import annotations

import copy
import sys
from pathlib import Path
from typing import Any

import yaml

from charon import Pipeline
from charon.core.schema import CompoundConfig, DoseProjectionConfig

REPO_ROOT = Path(__file__).resolve().parents[1]
PROPRA_YAML = REPO_ROOT / "validation" / "data" / "tier1_obach" / "compounds" / "propranolol.yaml"
PANEL_YAML = REPO_ROOT / "validation" / "data" / "fih_reference" / "panel.yaml"


def load_propranolol_entry() -> dict[str, Any]:
    panel = yaml.safe_load(PANEL_YAML.read_text())["panel"]
    return next(
        c for c in panel["compounds"]
        if c["name"] == "propranolol" and c["tier"] == "gold"
    )


def run_pipeline_for_compound(
    compound: CompoundConfig,
    entry: dict[str, Any],
) -> tuple[float, dict[str, Any], Any]:
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
        raise RuntimeError("Pipeline did not produce a dose recommendation")
    return (
        float(result.dose_recommendation.mrsd_mg),
        result.metadata,
        result.pk_parameters,
    )


def make_compound_with_multiplier(
    base_data: dict[str, Any], m: float | None
) -> CompoundConfig:
    data = copy.deepcopy(base_data)
    metab = data["properties"]["metabolism"]
    if m is None or m == 1.0:
        metab.pop("hepatic_clint_multiplier", None)
    else:
        metab["hepatic_clint_multiplier"] = {
            "value": float(m),
            "source": "experimental",
            "unit": "ratio",
            "method": f"sprint15 audit sweep m={m}",
        }
    return CompoundConfig.model_validate(data)


def main() -> int:
    base_data = yaml.safe_load(PROPRA_YAML.read_text())
    entry = load_propranolol_entry()
    ref_mg = float(entry["reference_fih_mg"])

    # --- Section 1: F-decomposition (no multiplier applied) ---
    base_compound = make_compound_with_multiplier(base_data, None)
    assert (
        base_compound.properties.metabolism.hepatic_clint_multiplier is None
    ), "audit baseline must have NO multiplier"
    mrsd_base, meta_base, pk_base = run_pipeline_for_compound(base_compound, entry)
    fold_base = ref_mg / mrsd_base if mrsd_base > 0 else float("inf")

    print("## §10 Sprint 15 audit — propranolol F-decomposition")
    print()
    print(f"**Generated:** by `scripts/sprint15_audit.py`")
    print(f"**Reference FIH dose:** {ref_mg} mg ({entry.get('f_source', 'n/a')})")
    print()
    print("### F-decomposition (no multiplier)")
    print()
    print("| Component | Value | Analytical expected | Within range? |")
    print("|---|---:|---:|:---:|")

    def in_range(v: float, lo: float, hi: float) -> str:
        return "yes" if lo <= v <= hi else "**NO**"

    print(
        f"| Fa | {pk_base.fa:.4f} | ~0.96 | {in_range(pk_base.fa, 0.85, 1.00)} |"
    )
    print(
        f"| Fg | {pk_base.fg:.4f} | ~1.00 | {in_range(pk_base.fg, 0.90, 1.00)} |"
    )
    print(
        f"| Fh | {pk_base.fh:.4f} | ~0.93 | {in_range(pk_base.fh, 0.85, 0.97)} |"
    )
    print(
        f"| F_oral | {pk_base.bioavailability:.4f} | ~0.89 | {in_range(pk_base.bioavailability, 0.78, 0.96)} |"
    )

    clint_liver_base = float(meta_base.get("clint_liver_L_h", float("nan")))
    print(
        f"| CLint_liver_L_h | {clint_liver_base:.2f} | ~51.0 | {in_range(clint_liver_base, 45.0, 60.0)} |"
    )
    print(
        f"| MRSD_base_mg | {mrsd_base:.4g} | 0.3502 (Sprint 14) | n/a |"
    )
    print(
        f"| fold_base | {fold_base:.2f} | 28.55 (Sprint 14) | n/a |"
    )
    print()

    # --- Section 2: multiplier sweep ---
    multipliers = [1, 2, 3, 5, 8, 12, 15, 20, 25, 30]
    print("### Multiplier sweep (fold vs m)")
    print()
    print("| m | MRSD_mg | fold_observed | within_3x | F_oral | CLint_liver_L_h |")
    print("|---:|---:|---:|:---:|---:|---:|")

    sweep_rows: list[tuple[int, float, float, bool]] = []
    for m in multipliers:
        compound = make_compound_with_multiplier(base_data, m)
        mrsd_m, meta_m, pk_m = run_pipeline_for_compound(compound, entry)
        fold_m = ref_mg / mrsd_m if mrsd_m > 0 else float("inf")
        within_3x = fold_m <= 3.0
        sweep_rows.append((m, mrsd_m, fold_m, within_3x))
        clint_liv = float(meta_m.get("clint_liver_L_h", float("nan")))
        F_m = pk_m.bioavailability
        flag = "yes" if within_3x else "no"
        print(
            f"| {m} | {mrsd_m:.4g} | {fold_m:.2f} | {flag} | {F_m:.3f} | {clint_liv:.2f} |"
        )

    # Identify boundary multipliers
    m_close_3x = next((m for m, _, fold, _ in sweep_rows if fold <= 3.0), None)
    m_close_1x_pair = min(sweep_rows, key=lambda r: abs(r[2] - 1.0))
    m_close_1x = m_close_1x_pair[0]

    print()
    print("### Boundary summary")
    print()
    print(f"- **m_close_3x:** {m_close_3x if m_close_3x is not None else 'NOT REACHED in sweep'}  (smallest m bringing fold ≤ 3x)")
    print(f"- **m_close_1x:** {m_close_1x}  (m bringing fold closest to 1.0)")
    print()
    print(
        "Use this empirical curve to choose the literature multiplier (Task 2). "
        "If literature supports a value ≥ m_close_3x, Branch A applies; if literature "
        "supports a value < m_close_3x but > 3, Branch B applies (close-but-not-quite); "
        "if literature is inconclusive or < 3, Branch C applies (null)."
    )

    return 0


if __name__ == "__main__":
    sys.exit(main())
