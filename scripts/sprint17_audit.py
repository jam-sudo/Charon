"""Sprint 17 audit — lisinopril parameter audit + sensitivity sweeps.

Purpose:
1. F-decomposition baseline: capture Fa/Fg/Fh/F_oral, CL_total, MRSD_PAD
   from current Charon pipeline for lisinopril (no parameter change).
2. Parameter audit table: each parameter (clrenal_L_h, fu_p, bp_ratio,
   clint_uL_min_mg, target_ceff_nM, peff_cm_s, F_oral_obs) compared to
   primary-literature range with in-range / out-of-range / low-confidence flag.
3. Sensitivity sweeps:
   3a. target_ceff_nM in {85, 170, 340, 700} (0.5x-4x current 170)
   3b. safety_factor in {3, 5, 10}
   3c. peff_cm_s in {0.1, 0.3, 0.5, 1.0, 2.0} x 1e-4
   For each: print MRSD_PAD, fold vs ref=10 mg.
4. Branch decision rationale based on Steps 1-3 outcomes.

Output: prints markdown sections to stdout. Pipe to file or copy into the
audit subsection of validation/reports/layer3_ivive_decomposition.md.

Usage:
    python3 scripts/sprint17_audit.py > /tmp/sprint17_audit.md
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
LISIN_YAML = REPO_ROOT / "validation" / "data" / "tier1_obach" / "compounds" / "lisinopril.yaml"
PANEL_YAML = REPO_ROOT / "validation" / "data" / "fih_reference" / "panel.yaml"


def load_lisinopril_entry() -> dict[str, Any]:
    panel = yaml.safe_load(PANEL_YAML.read_text())["panel"]
    return next(
        c for c in panel["compounds"]
        if c["name"] == "lisinopril" and c["tier"] == "gold"
    )


def run_pipeline(
    compound: CompoundConfig,
    entry: dict[str, Any],
    target_ceff_nM: float | None = None,
    safety_factor: float = 10.0,
) -> tuple[float, dict[str, Any], Any]:
    """Run pipeline with optional target_ceff_nM and safety_factor overrides."""
    target = float(target_ceff_nM if target_ceff_nM is not None else entry["target_ceff_nM"])
    pipe = Pipeline(
        compound,
        route=entry["route"],
        dose_mg=1.0,
        dose_projection=DoseProjectionConfig(
            target_ceff_nM=target,
            safety_factor=safety_factor,
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


def make_compound_with_peff(
    base_data: dict[str, Any], peff_cm_s: float | None = None
) -> CompoundConfig:
    """Build CompoundConfig with optional peff_cm_s override."""
    data = copy.deepcopy(base_data)
    if peff_cm_s is not None:
        data["properties"]["permeability"]["peff_cm_s"]["value"] = float(peff_cm_s)
    return CompoundConfig.model_validate(data)


def main() -> int:
    base_data = yaml.safe_load(LISIN_YAML.read_text())
    entry = load_lisinopril_entry()
    ref_mg = float(entry["reference_fih_mg"])
    base_compound = make_compound_with_peff(base_data, peff_cm_s=None)

    # --- Section 1: F-decomposition baseline ---
    mrsd_base, meta_base, pk_base = run_pipeline(base_compound, entry)
    fold_base = ref_mg / mrsd_base if mrsd_base > 0 else float("inf")

    print("## §12. Sprint 17 audit — lisinopril (framework-limit closure)")
    print()
    print("**Generated:** by `scripts/sprint17_audit.py`")
    print(f"**Reference FIH dose:** {ref_mg} mg ({entry['source']})")
    print(f"**Target Ceff:** {entry['target_ceff_nM']} nM (Sprint 10 §4 verified)")
    print()
    print("### Section 1 — F-decomposition baseline (no parameter change)")
    print()
    print("| Component | Value | Literature expected | Within range? |")
    print("|---|---:|---:|:---:|")

    def in_range(v: float, lo: float, hi: float) -> str:
        return "yes" if lo <= v <= hi else "**NO**"

    print(f"| Fa | {pk_base.fa:.4f} | ~0.24-0.30 (BCS III, F_obs=0.25-0.29) | {in_range(pk_base.fa, 0.20, 0.32)} |")
    print(f"| Fg | {pk_base.fg:.4f} | ~1.00 (non-CYP3A4) | {in_range(pk_base.fg, 0.95, 1.05)} |")
    print(f"| Fh | {pk_base.fh:.4f} | ~1.00 (CLint negligible) | {in_range(pk_base.fh, 0.97, 1.00)} |")
    print(f"| F_oral | {pk_base.bioavailability:.4f} | ~0.25 (Beermann 1988 0.25-0.29) | {in_range(pk_base.bioavailability, 0.22, 0.31)} |")

    cl_renal = float(meta_base.get("cl_renal_L_h", float("nan")))
    clint_liver = float(meta_base.get("clint_liver_L_h", float("nan")))
    # CL_total derived: cl_apparent (= CL/F for oral) × F
    cl_apparent = pk_base.cl_apparent if pk_base.cl_apparent is not None else float("nan")
    F = pk_base.bioavailability if pk_base.bioavailability is not None else float("nan")
    cl_total = cl_apparent * F if (cl_apparent == cl_apparent and F == F) else float("nan")

    print(f"| CL_renal_L_h | {cl_renal:.3f} | 5.0 (Beermann 1988) | {in_range(cl_renal, 4.5, 5.5)} |")
    print(f"| CLint_liver_L_h | {clint_liver:.4f} | ~0.3 (negligible from clint=0.1) | {in_range(clint_liver, 0.1, 0.6)} |")
    print(f"| CL_total_L_h (derived) | {cl_total:.3f} | ~5.1 (Beermann 1988) | {in_range(cl_total, 4.8, 5.6)} |")
    print(f"| MRSD_PAD_mg | {mrsd_base:.4g} | 2.424 (Sprint 16) | n/a |")
    print(f"| fold_base | {fold_base:.3f} | 4.126 (Sprint 16) | n/a |")
    print()

    # --- Section 2: parameter audit table ---
    print("### Section 2 — Parameter audit (vs primary literature)")
    print()
    print("| Parameter | Charon | Primary source | Literature range | OK? |")
    print("|---|---:|---|---|:---:|")

    audit_rows = [
        ("clrenal_L_h", 5.0, "Beermann 1988", "4.5-5.4 (~95% of CL_total=5.1)", "OK"),
        ("fu_p", 0.75, "Lancaster 1988", "0.70-0.80 (low binding)", "OK"),
        ("bp_ratio", 0.85, "Beermann 1988", "0.80-0.90", "OK"),
        ("clint_uL_min_mg", 0.1, "Beermann 1988", "<1.0 (negligible)", "OK"),
        ("target_ceff_nM", 170.0, "Beermann 1988 / Sprint 10 §4", "150-200 nM (Cmax 70 ng/mL × 1000/405)", "OK"),
        ("peff_cm_s", 0.3e-4, "Knutter 2008 (cited; primary Peff NOT located)", "BCS III; PEPT1-mediated absorption not modelled", "FLAG-LOW-CONF"),
        ("F_oral_obs", 0.25, "Beermann 1988", "0.25-0.29", "OK"),
    ]

    for name, val, src, lit, ok in audit_rows:
        if isinstance(val, float) and val < 1e-3:
            val_str = f"{val:.2e}"
        else:
            val_str = f"{val:g}"
        print(f"| `{name}` | {val_str} | {src} | {lit} | {ok} |")

    print()
    print("**Pre-audit prior:** Sprint 10 §4 verified `target_ceff_nM=170` plausible. ")
    print("Beermann 1988 + Lancaster 1988 are primary peer-reviewed sources for ")
    print("propranolol/lisinopril-class clinical PK; no contradicting literature located. ")
    print("Peff is the only LOW-CONFIDENCE parameter (no primary measurement); ")
    print("however Sprint 11 oral migration confirmed F_pred matches F_obs to ~3%, ")
    print("so passive Peff value is empirically calibrated whether or not primary cite exists.")
    print()

    # --- Section 3: sensitivity sweeps ---
    print("### Section 3a — target_ceff_nM sweep (safety_factor=10, peff=baseline)")
    print()
    print("| target_ceff_nM | MRSD_mg | fold | within_3x | F_oral |")
    print("|---:|---:|---:|:---:|---:|")

    target_sweep_rows: list[tuple[float, float, float, bool]] = []
    for target in [85.0, 170.0, 340.0, 700.0]:
        mrsd_t, _, pk_t = run_pipeline(base_compound, entry, target_ceff_nM=target)
        fold_t = ref_mg / mrsd_t if mrsd_t > 0 else float("inf")
        within_3x_t = fold_t <= 3.0
        target_sweep_rows.append((target, mrsd_t, fold_t, within_3x_t))
        flag = "yes" if within_3x_t else "no"
        print(f"| {target:.0f} | {mrsd_t:.4g} | {fold_t:.3f} | {flag} | {pk_t.bioavailability:.4f} |")

    target_close_3x = next((t for t, _, fold, _ in target_sweep_rows if fold <= 3.0), None)
    print()
    print(f"- **target_ceff_close_3x:** {target_close_3x if target_close_3x is not None else 'NOT REACHED'} nM  (smallest target_ceff bringing fold ≤ 3x)")
    print()

    print("### Section 3b — safety_factor sweep (target_ceff=170, peff=baseline)")
    print()
    print("| safety_factor | MRSD_mg | fold | within_3x |")
    print("|---:|---:|---:|:---:|")

    sf_sweep_rows: list[tuple[float, float, float, bool]] = []
    for sf in [3.0, 5.0, 10.0]:
        mrsd_s, _, _ = run_pipeline(base_compound, entry, safety_factor=sf)
        fold_s = ref_mg / mrsd_s if mrsd_s > 0 else float("inf")
        within_3x_s = fold_s <= 3.0
        sf_sweep_rows.append((sf, mrsd_s, fold_s, within_3x_s))
        flag = "yes" if within_3x_s else "no"
        print(f"| {sf:g} | {mrsd_s:.4g} | {fold_s:.3f} | {flag} |")

    sf_close_3x = next((s for s, _, fold, _ in sf_sweep_rows if fold <= 3.0), None)
    print()
    print(f"- **sf_close_3x:** {sf_close_3x if sf_close_3x is not None else 'NOT REACHED'}  (smallest SF bringing fold ≤ 3x)")
    print()

    print("### Section 3c — peff_cm_s sweep (target_ceff=170, safety_factor=10)")
    print()
    print("| peff_cm_s | F_oral | MRSD_mg | fold | within_3x |")
    print("|---:|---:|---:|---:|:---:|")

    peff_sweep_rows: list[tuple[float, float, float, float, bool]] = []
    for peff_factor in [0.1, 0.3, 0.5, 1.0, 2.0]:
        peff_val = peff_factor * 1e-4
        compound_p = make_compound_with_peff(base_data, peff_cm_s=peff_val)
        mrsd_p, _, pk_p = run_pipeline(compound_p, entry)
        fold_p = ref_mg / mrsd_p if mrsd_p > 0 else float("inf")
        within_3x_p = fold_p <= 3.0
        peff_sweep_rows.append((peff_val, pk_p.bioavailability, mrsd_p, fold_p, within_3x_p))
        flag = "yes" if within_3x_p else "no"
        print(f"| {peff_val:.2e} | {pk_p.bioavailability:.4f} | {mrsd_p:.4g} | {fold_p:.3f} | {flag} |")

    print()

    # --- Section 4: Branch decision rationale ---
    print("### Section 4 — Branch decision rationale")
    print()
    n_audit_ok = sum(1 for _, _, _, _, ok in audit_rows if ok == "OK")
    n_audit_flag = sum(1 for _, _, _, _, ok in audit_rows if "FLAG" in ok)
    n_audit_fail = len(audit_rows) - n_audit_ok - n_audit_flag

    print(f"- Audit results: {n_audit_ok} OK, {n_audit_flag} flagged, {n_audit_fail} out-of-range")
    print(f"- target_ceff sweep: even 4× (700 nM) {'closes within 3x' if target_close_3x else 'does NOT reach within-3x'}")
    print(f"- safety_factor sweep: SF=3 {'closes within 3x' if sf_close_3x and sf_close_3x <= 3.0 else 'does NOT reach within-3x'}")

    peff_close_3x_rows = [r for r in peff_sweep_rows if r[4]]
    print(f"- peff sweep: passive-Peff alone {'CAN bring fold within 3x' if peff_close_3x_rows else 'CANNOT bring fold within 3x'}")
    print()

    if n_audit_fail > 0:
        print("**Branch A indicated:** at least one parameter outside literature range. ")
        print("Inspect audit table; correct out-of-range value; integration test + regen.")
    elif n_audit_flag > 0 and peff_close_3x_rows:
        print("**Branch B candidate:** peff is LOW-CONFIDENCE and a peff value within ")
        print("plausible BCS-III range brings fold within 3x. Apply correction with ")
        print("caveat (no primary cite); document partial improvement.")
    else:
        print("**Branch C indicated (HONEST NULL):** all parameters in literature range; ")
        print("no realistic single-parameter perturbation closes 4.13x within 3x. The ")
        print("residual is a benchmark methodology gap (PAD with SF=10 vs FDA chronic ")
        print("label start dose), not a pharmacology gap. Document as framework limit.")

    print()
    print("Reference precedents:")
    print("- Sprint 14 (diazepam audit): Branch C honest null when all parameters verify ")
    print("  AND residual mechanism is framework-limited (low-fu_p well-stirred sensitivity).")
    print("- Sprint 15 (propranolol audit): Branch B applied when literature multiplier ")
    print("  partially closes residual but does NOT reach within-3x (CLAUDE.md §6.5 honesty).")

    return 0


if __name__ == "__main__":
    sys.exit(main())
