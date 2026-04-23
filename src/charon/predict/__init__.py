"""Charon Layer 1 property prediction: SMILES → CompoundProperties.

Public surface:

    predict_properties(smiles, predictor=None, conformal=None)
        Orchestrator that runs the full prediction pipeline and returns
        a populated :class:`charon.core.schema.CompoundProperties`.

    ADMETPredictor, ADMEPrediction
        XGBoost ensemble for fup and CLint.

    ConformalPredictor, CoverageReport
        Log-space symmetric conformal intervals.

    predict_pka, PKaResult
        Rule-based pKa estimation.

    predict_fu_inc
        Austin 2002 microsomal unbound fraction correlation.

    predict_bp_ratio
        Empirical blood:plasma ratio.

    estimate_renal_clearance
        Thin wrapper around ParameterBridge.

    compute_features
        Shared 2057-D molecular feature vector (also used by training
        scripts).
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Literal

from charon.core.molecule import Molecule
from charon.core.schema import (
    BindingProperties,
    CompoundProperties,
    MetabolismProperties,
    PermeabilityProperties,
    PhysicochemicalProperties,
    PredictedProperty,
    RenalProperties,
    SafetyProperties,
    SourceType,
)
from charon.predict.admet_ensemble import ADMEPrediction, ADMETPredictor
from charon.predict.bp_ratio import predict_bp_ratio
from charon.predict.conformal import (
    ConformalPredictor,
    CoverageReport,
    get_default_conformal,
)
from charon.predict.features import compute_features
from charon.predict.fu_inc import predict_fu_inc
from charon.predict.pka import PKaResult, predict_pka
from charon.predict.renal import estimate_renal_clearance


class _ConformalOff:
    """Sentinel value: explicitly disable default CI attachment.

    Passed as ``conformal=CONFORMAL_OFF`` to :func:`predict_properties`
    to opt out of the module-level default predictor.
    """

    def __repr__(self) -> str:  # pragma: no cover - trivial
        return "CONFORMAL_OFF"


CONFORMAL_OFF = _ConformalOff()


__all__ = [
    "ADMEPrediction",
    "ADMETPredictor",
    "CONFORMAL_OFF",
    "CompoundProperties",
    "ConformalPredictor",
    "CoverageReport",
    "PKaResult",
    "PredictedProperty",
    "compute_features",
    "estimate_renal_clearance",
    "predict_bp_ratio",
    "predict_fu_inc",
    "predict_pka",
    "predict_properties",
]


def _predicted(
    value: float | None,
    source: SourceType,
    unit: str | None = None,
    ci_lower: float | None = None,
    ci_upper: float | None = None,
    method: str | None = None,
    flag: str | None = None,
) -> PredictedProperty | None:
    """Wrap a raw scalar into a ``PredictedProperty`` (None-safe)."""
    if value is None:
        return None
    return PredictedProperty(
        value=float(value),
        ci_90_lower=None if ci_lower is None else float(ci_lower),
        ci_90_upper=None if ci_upper is None else float(ci_upper),
        source=source,
        unit=unit,
        method=method,
        flag=flag,
    )


def predict_properties(
    smiles: str,
    predictor: ADMETPredictor | None = None,
    conformal: "ConformalPredictor | _ConformalOff | None" = None,
) -> CompoundProperties:
    """Run the full Layer 1 prediction pipeline on a SMILES.

    Steps:

    1. Parse SMILES → :class:`Molecule` (validates structure).
    2. Run rule-based pKa → pKa_acid, pKa_base, compound_type.
    3. Run XGBoost ADME → fup, CLint (hepatocyte units).
    4. Apply Austin 2002 → fu_inc from cLogP (HLM only).
    5. Apply empirical B:P formula using pKa-informed Kp_rbc prior.
    6. Estimate renal clearance = fu_p × GFR.
    7. Attach log-space conformal intervals where calibrated.
    8. Assemble :class:`CompoundProperties` with source tags.

    Args:
        smiles: Input SMILES (canonical or raw).
        predictor: Optional pre-initialised :class:`ADMETPredictor`.
            Instantiated on demand when ``None``.
        conformal: Controls CI attachment.

            * ``None`` (default): use the module-level default
              :class:`ConformalPredictor` from
              :func:`get_default_conformal`. If its calibration data are
              missing, a warning is logged and CIs are omitted.
            * :data:`CONFORMAL_OFF`: explicitly disable CI attachment
              (leaves ``ci_90_lower`` / ``ci_90_upper`` as ``None``).
            * A user-supplied :class:`ConformalPredictor`: used as-is.

    Returns:
        Fully populated :class:`CompoundProperties`.

    Raises:
        ValueError: If ``smiles`` is invalid.
        FileNotFoundError: If ADMET model files are missing.
    """
    # 1. Validate structure (raises ValueError on bad SMILES).
    mol = Molecule(smiles)
    descriptors = mol.descriptors()

    # 2. Rule-based pKa / compound_type.
    pka_result = predict_pka(smiles)

    # 3. XGBoost ADME (fup + CLint).
    predictor = predictor or ADMETPredictor()
    adme = predictor.predict(smiles)

    # 4. fu_inc from Crippen cLogP (Austin 2002).
    clogp = float(descriptors["logP"])
    fu_inc_value = predict_fu_inc(clogp)

    # 5. Empirical B:P ratio using compound_type prior.
    bp_value = predict_bp_ratio(adme.fup, compound_type=pka_result.compound_type)

    # 6. Renal clearance (filtration only by default).
    cl_renal = estimate_renal_clearance(fu_p=adme.fup)

    # 7. Conformal intervals (default = module singleton; CONFORMAL_OFF = disabled).
    fup_lo: float | None = None
    fup_hi: float | None = None
    clint_lo: float | None = None
    clint_hi: float | None = None
    if isinstance(conformal, _ConformalOff):
        active: ConformalPredictor | None = None
    elif conformal is None:
        try:
            active = get_default_conformal()
        except (FileNotFoundError, OSError) as exc:
            import logging as _logging
            _logging.getLogger(__name__).warning(
                "Default conformal unavailable (%s); continuing without CI", exc
            )
            active = None
    else:
        active = conformal

    if active is not None:
        if active.is_calibrated("fup"):
            fup_lo, fup_hi = active.get_interval("fup", adme.fup)
        if active.is_calibrated("clint_hepatocyte"):
            clint_lo, clint_hi = active.get_interval(
                "clint_hepatocyte", adme.clint_hepatocyte
            )

    # 8. Assemble Pydantic schema.
    physicochemical = PhysicochemicalProperties(
        logp=_predicted(clogp, source="derived", unit=""),
        pka_acid=_predicted(pka_result.pka_acid, source="ml_pka"),
        pka_base=_predicted(pka_result.pka_base, source="ml_pka"),
    )
    permeability = PermeabilityProperties()
    binding = BindingProperties(
        fu_p=_predicted(
            adme.fup,
            source="ml_ensemble",
            unit="fraction",
            ci_lower=fup_lo,
            ci_upper=fup_hi,
        ),
        fu_inc=_predicted(
            fu_inc_value,
            source="correlation",
            unit="fraction",
            method="Austin 2002 (HLM 1 mg/mL)",
        ),
        bp_ratio=_predicted(
            bp_value,
            source="derived",
            unit="ratio",
            method=f"empirical Hct 0.45, compound_type={pka_result.compound_type}",
        ),
    )
    metabolism = MetabolismProperties(
        clint_uL_min_mg=_predicted(
            adme.clint_hepatocyte,
            source="ml_ensemble",
            unit="uL/min/10^6 cells",
            ci_lower=clint_lo,
            ci_upper=clint_hi,
            flag=(
                "clint_tier2_ml; hepatocyte units; experimental value recommended"
            ),
        ),
    )
    safety = SafetyProperties()
    renal = RenalProperties(
        clrenal_L_h=_predicted(
            cl_renal,
            source="derived",
            unit="L/h",
            method="fu_p * GFR (filtration only)",
        ),
        active_secretion=False,
    )

    return CompoundProperties(
        physicochemical=physicochemical,
        permeability=permeability,
        binding=binding,
        metabolism=metabolism,
        safety=safety,
        renal=renal,
    )
