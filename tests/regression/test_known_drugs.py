"""Regression guards for the 5 CLAUDE.md reference drugs.

WARNING: NOAEL values below are REGRESSION FIXTURES, not literature-verified
safety data.  They are chosen to exercise the dose-projection code path and
MUST NOT be used for any real clinical or preclinical decision-making.

To regenerate golden baselines after an intentional model change:

    UPDATE_BASELINES=1 pytest tests/regression/test_known_drugs.py

Target runtime: < 30 seconds total (no uncertainty quantification).

These tests verify that pipeline output does NOT drift silently after model
changes.  They do NOT validate accuracy against literature values.
"""

from __future__ import annotations

import json
import math
import os
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

import pytest

from charon.core.schema import DoseProjectionConfig
from charon.pipeline import Pipeline

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

GOLDEN_DIR = Path(__file__).parent / "golden_outputs"
GOLDEN_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# RefDrug dataclass
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class RefDrug:
    key: str
    smiles: str
    route: str
    dose_mg: float
    noael_mg_kg: float | None
    noael_species: str | None
    target_kd_nM: float | None
    target_ceff_nM: float | None
    phenotype_caveat: bool


# ---------------------------------------------------------------------------
# Reference drugs (5 CLAUDE.md reference compounds)
# ---------------------------------------------------------------------------

REFERENCE_DRUGS = (
    RefDrug(
        "midazolam",
        "Cc1ncc2n1-c1ccc(Cl)cc1C(c1ccccc1F)=NC2",
        route="iv_bolus",
        dose_mg=5.0,
        noael_mg_kg=2.0,
        noael_species="rat",
        target_kd_nM=None,
        target_ceff_nM=None,
        phenotype_caveat=False,
    ),
    RefDrug(
        "warfarin",
        "CC(=O)CC(c1ccccc1)c1c(O)c2ccccc2oc1=O",
        route="iv_bolus",
        dose_mg=5.0,
        noael_mg_kg=0.5,
        noael_species="rat",
        target_kd_nM=None,
        target_ceff_nM=None,
        phenotype_caveat=False,
    ),
    RefDrug(
        "diclofenac",
        "OC(=O)Cc1ccccc1Nc1c(Cl)cccc1Cl",
        route="iv_bolus",
        dose_mg=75.0,
        noael_mg_kg=5.0,
        noael_species="rat",
        target_kd_nM=None,
        target_ceff_nM=None,
        phenotype_caveat=False,
    ),
    RefDrug(
        "omeprazole",
        "COc1ccc2[nH]c([S+]([O-])Cc3ncc(C)c(OC)c3C)nc2c1",
        route="iv_bolus",
        dose_mg=40.0,
        noael_mg_kg=4.0,
        noael_species="rat",
        target_kd_nM=None,
        target_ceff_nM=None,
        phenotype_caveat=True,
    ),
    RefDrug(
        "dextromethorphan",
        "COc1ccc2c(c1)[C@]13CCCC[C@@H]1[C@H](C2)N(C)CC3",
        route="iv_bolus",
        dose_mg=30.0,
        noael_mg_kg=10.0,
        noael_species="rat",
        target_kd_nM=None,
        target_ceff_nM=None,
        phenotype_caveat=True,
    ),
)

# ---------------------------------------------------------------------------
# Helper: IDs for pytest parameterize
# ---------------------------------------------------------------------------


def _drug_id(drug: RefDrug) -> str:
    return drug.key


# ---------------------------------------------------------------------------
# Pipeline runner
# ---------------------------------------------------------------------------


def _run_pipeline(drug: RefDrug):
    """Build and run the pipeline for a single reference drug."""
    dp: DoseProjectionConfig | None = None
    if drug.noael_mg_kg is not None and drug.noael_species is not None:
        dp = DoseProjectionConfig(
            noael_mg_kg=drug.noael_mg_kg,
            noael_species=drug.noael_species,
            target_kd_nM=drug.target_kd_nM,
            target_ceff_nM=drug.target_ceff_nM,
        )
    return Pipeline.from_smiles(
        drug.smiles,
        route=drug.route,
        dose_mg=drug.dose_mg,
        compound_name=drug.key,
        dose_projection=dp,
    ).run()


# ---------------------------------------------------------------------------
# Primary metric extraction
# ---------------------------------------------------------------------------


def _extract_primary(result) -> dict:
    """Return the primary regression metric as {name, value}.

    Priority:
      1. mrsd_mg  (dose recommendation)
      2. auc_0_inf
      3. cmax
    """
    dr = result.dose_recommendation
    if dr is not None and math.isfinite(dr.mrsd_mg):
        return {"name": "mrsd_mg", "value": float(dr.mrsd_mg)}

    auc = result.pk_parameters.auc_0_inf
    if auc is not None and math.isfinite(auc):
        return {"name": "auc_0_inf", "value": float(auc)}

    cmax = result.pk_parameters.cmax
    if cmax is not None and math.isfinite(cmax):
        return {"name": "cmax", "value": float(cmax)}

    raise RuntimeError(
        "No finite primary metric (mrsd_mg, auc_0_inf, or cmax) in result"
    )


# ---------------------------------------------------------------------------
# Baseline I/O
# ---------------------------------------------------------------------------


def _baseline_path(drug: RefDrug) -> Path:
    return GOLDEN_DIR / f"{drug.key}.json"


def _write_baseline(path: Path, primary: dict, drug: RefDrug, result) -> None:
    """Write a golden baseline JSON file."""
    pk = result.pk_parameters

    try:
        from importlib.metadata import version as _pkg_version
        charon_version = _pkg_version("charon")
    except Exception:
        charon_version = "0.1.0"

    caveat = (
        "phenotype_caveat: ML predictions for this compound may be less reliable "
        "due to metabolic phenotype complexity (e.g. CYP2D6 polymorphism)"
        if drug.phenotype_caveat
        else None
    )

    payload = {
        "drug": drug.key,
        "smiles": drug.smiles,
        "charon_version": charon_version,
        "primary_metric": primary,
        "secondary_metrics": {
            "cmax_ug_L": round(float(pk.cmax), 4) if pk.cmax is not None else None,
            "auc_ug_h_L": round(float(pk.auc_0_inf), 4) if pk.auc_0_inf is not None else None,
            "t_half_h": round(float(pk.half_life), 4) if pk.half_life is not None else None,
        },
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "caveat": caveat,
    }
    path.write_text(json.dumps(payload, indent=2))


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("drug", REFERENCE_DRUGS, ids=_drug_id)
def test_pipeline_runs(drug: RefDrug):
    """Smoke test: pipeline completes and produces non-zero PK output.

    Checks:
    - run() returns without exception
    - cmax > 0
    - auc_0_inf > 0
    - mrsd_mg is finite if dose_recommendation is present
    """
    result = _run_pipeline(drug)
    pk = result.pk_parameters

    assert pk.cmax is not None and pk.cmax > 0, (
        f"{drug.key}: cmax must be > 0, got {pk.cmax}"
    )
    assert pk.auc_0_inf is not None and pk.auc_0_inf > 0, (
        f"{drug.key}: auc_0_inf must be > 0, got {pk.auc_0_inf}"
    )

    dr = result.dose_recommendation
    if dr is not None:
        assert math.isfinite(dr.mrsd_mg), (
            f"{drug.key}: mrsd_mg must be finite, got {dr.mrsd_mg}"
        )


@pytest.mark.parametrize("drug", REFERENCE_DRUGS, ids=_drug_id)
def test_dose_within_snapshot(drug: RefDrug):
    """Snapshot test: primary metric must match golden baseline within ±25%.

    Fold-error tolerance: 1.25 (i.e. result / baseline in [0.80, 1.25]).

    To regenerate baselines after an intentional model change:
        UPDATE_BASELINES=1 pytest tests/regression/test_known_drugs.py

    If a baseline is missing, the test fails with instructions to regenerate.
    """
    baseline_path = _baseline_path(drug)
    update_mode = os.environ.get("UPDATE_BASELINES", "").strip() == "1"

    result = _run_pipeline(drug)
    primary = _extract_primary(result)

    if update_mode:
        _write_baseline(baseline_path, primary, drug, result)
        pytest.skip(
            f"[UPDATE_BASELINES] Wrote new baseline for {drug.key}: "
            f"{primary['name']}={primary['value']:.4g}"
        )

    if not baseline_path.exists():
        pytest.fail(
            f"Golden baseline missing for '{drug.key}'.\n"
            f"  Expected: {baseline_path}\n"
            f"  To generate baselines, run:\n"
            f"    UPDATE_BASELINES=1 pytest tests/regression/test_known_drugs.py\n"
            f"  Then commit the files under tests/regression/golden_outputs/."
        )

    baseline = json.loads(baseline_path.read_text())
    expected_metric = baseline["primary_metric"]

    assert expected_metric["name"] == primary["name"], (
        f"{drug.key}: primary metric changed from "
        f"'{expected_metric['name']}' to '{primary['name']}'. "
        f"Regenerate baselines with UPDATE_BASELINES=1."
    )

    expected_val = float(expected_metric["value"])
    actual_val = float(primary["value"])

    assert expected_val > 0, f"{drug.key}: baseline value must be > 0"
    assert actual_val > 0, f"{drug.key}: computed value must be > 0"

    fold_error = max(actual_val / expected_val, expected_val / actual_val)
    assert fold_error <= 1.25, (
        f"{drug.key} [{primary['name']}]: fold error {fold_error:.3f} exceeds 1.25x.\n"
        f"  Baseline: {expected_val:.4g}\n"
        f"  Current:  {actual_val:.4g}\n"
        f"  If this change is intentional, regenerate baselines:\n"
        f"    UPDATE_BASELINES=1 pytest tests/regression/test_known_drugs.py"
    )
