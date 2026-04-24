"""Layer 1-to-2 Bridge: Convert raw predictions into PBPK-ready parameters.

This is the most error-prone junction in the pipeline.  Every conversion
is a named, tested, logged function.  Every call returns a
:class:`ConversionLog` audit trail.

CRITICAL RULES:
  - ``fu_inc`` correction applies to HLM **only**, never hepatocytes.
  - ``fu_p`` enters **only** as ``fu_b = fu_p / BP`` in the liver model,
    NOT in CLint scaling.
  - Every conversion step is logged in :class:`ConversionLog`.
"""

from __future__ import annotations

import logging
import math

from charon.core.liver_models import get_liver_model
from charon.core.schema import ConversionLog, ConversionStep, HepaticClearance
from charon.core.units import (
    HUMAN_GFR_ML_MIN,
    HUMAN_HEPATOCELLULARITY,
    HUMAN_LIVER_WEIGHT_G,
    HUMAN_MPPGL,
    HUMAN_QH_L_H,
    mL_min_to_L_h,
    nm_s_to_cm_s,
    uL_min_to_L_h,
)

logger = logging.getLogger(__name__)

_VALID_SYSTEMS = {"HLM", "hepatocytes"}


class ParameterBridge:
    """Layer 1-to-2 Bridge: Convert raw predictions into PBPK-ready parameters.

    Every public method returns either a scalar result or a structured
    result containing a full :class:`ConversionLog` audit trail.
    """

    # ------------------------------------------------------------------
    # IVIVE: CLint -> CLh
    # ------------------------------------------------------------------

    def clint_to_clh(
        self,
        clint: float,
        fu_inc: float,
        fu_p: float,
        system: str = "HLM",
        mppgl: float = HUMAN_MPPGL,
        hepatocellularity: float = HUMAN_HEPATOCELLULARITY,
        liver_weight_g: float = HUMAN_LIVER_WEIGHT_G,
        qh_L_h: float = HUMAN_QH_L_H,
        bp_ratio: float = 1.0,
        model: str = "well_stirred",
        clint_multiplier: float | None = None,
    ) -> HepaticClearance:
        """IVIVE: in-vitro CLint to in-vivo hepatic clearance.

        Parameters
        ----------
        clint : float
            In-vitro intrinsic clearance.
            Units: uL/min/mg protein (HLM) or uL/min/10^6 cells (hepatocytes).
        fu_inc : float
            Microsomal unbound fraction, in (0, 1].  Applied for HLM only.
        fu_p : float
            Plasma unbound fraction, in (0, 1].
        system : str
            ``"HLM"`` or ``"hepatocytes"``.
        mppgl : float
            Microsomal protein per gram liver (mg/g).
        hepatocellularity : float
            Hepatocellularity (10^6 cells / g liver).
        liver_weight_g : float
            Liver weight in grams.
        qh_L_h : float
            Hepatic blood flow in L/h.
        bp_ratio : float
            Blood-to-plasma ratio (must be > 0).
        model : str
            Liver model name (``"well_stirred"``, ``"parallel_tube"``,
            ``"dispersion"``).

        Returns
        -------
        HepaticClearance
            Result with ``clh_L_h``, ``extraction_ratio``, ``model_used``,
            and a complete ``conversion_log`` audit trail.

        Raises
        ------
        ValueError
            If any input violates physical constraints.
        """
        # ---- Validation ------------------------------------------------
        if system not in _VALID_SYSTEMS:
            raise ValueError(
                f"system must be one of {sorted(_VALID_SYSTEMS)}, got {system!r}"
            )
        if fu_p <= 0:
            raise ValueError(f"fu_p must be > 0, got {fu_p}")
        if fu_p < 0.01:
            logger.warning(
                "fu_p=%.4f — extreme binding; CLh highly sensitive to fu_p accuracy",
                fu_p,
            )
        if system == "HLM" and fu_inc <= 0:
            raise ValueError(f"fu_inc must be > 0 for HLM system, got {fu_inc}")
        if bp_ratio <= 0:
            raise ValueError(f"bp_ratio must be > 0, got {bp_ratio}")
        if clint < 0:
            raise ValueError(f"clint must be >= 0, got {clint}")
        if clint_multiplier is not None and clint_multiplier <= 0:
            raise ValueError(
                f"clint_multiplier must be > 0, got {clint_multiplier}"
            )

        # Resolve the liver model function early so we fail fast on bad names.
        liver_model_fn = get_liver_model(model)

        # ---- Build audit trail -----------------------------------------
        steps: list[ConversionStep] = []
        input_params: dict[str, float | str] = {
            "clint": clint,
            "fu_inc": fu_inc,
            "fu_p": fu_p,
            "system": system,
            "mppgl": mppgl,
            "hepatocellularity": hepatocellularity,
            "liver_weight_g": liver_weight_g,
            "qh_L_h": qh_L_h,
            "bp_ratio": bp_ratio,
            "model": model,
        }
        if clint_multiplier is not None:
            input_params["clint_multiplier"] = clint_multiplier

        # ---- Step 0: Scaling factor ------------------------------------
        if system == "HLM":
            scale_factor = mppgl * liver_weight_g
            scale_formula = f"MPPGL({mppgl}) * liver_weight({liver_weight_g}g)"
            scale_unit = "mg protein"
        else:
            scale_factor = hepatocellularity * liver_weight_g
            scale_formula = (
                f"hepatocellularity({hepatocellularity}) "
                f"* liver_weight({liver_weight_g}g)"
            )
            scale_unit = "10^6 cells"

        steps.append(
            ConversionStep(
                name="scale_factor",
                value=scale_factor,
                unit=scale_unit,
                formula=scale_formula,
            )
        )

        # ---- Step 1: Unbound CLint -------------------------------------
        if system == "HLM":
            clint_u = clint / fu_inc
            clint_u_formula = f"CLint/fu_inc = {clint}/{fu_inc}"
            clint_u_unit = "uL/min/mg"
        else:
            clint_u = clint  # No fu_inc correction for hepatocytes
            clint_u_formula = "CLint (no fu_inc correction for hepatocytes)"
            clint_u_unit = "uL/min/10^6 cells"

        steps.append(
            ConversionStep(
                name="CLint_u",
                value=clint_u,
                unit=clint_u_unit,
                formula=clint_u_formula,
            )
        )

        # ---- Step 2: Scale to whole liver ------------------------------
        clint_liver_uL_min = clint_u * scale_factor
        steps.append(
            ConversionStep(
                name="CLint_liver",
                value=clint_liver_uL_min,
                unit="uL/min",
                formula=(
                    f"CLint_u({clint_u:.4g}) * scale_factor({scale_factor:.4g})"
                ),
            )
        )

        clint_liver_L_h = uL_min_to_L_h(clint_liver_uL_min)
        steps.append(
            ConversionStep(
                name="CLint_liver_Lh",
                value=clint_liver_L_h,
                unit="L/h",
                formula=f"{clint_liver_uL_min:.4g} / 1e6 * 60",
            )
        )

        # ---- Step 3b (optional): CLint enhancement (Sprint 12) ---------
        # Empirical multiplier for uptake-limited substrates (e.g., OATP1B1).
        # Sourced from CompoundConfig.properties.metabolism.hepatic_clint_multiplier.
        if clint_multiplier is not None and clint_multiplier != 1.0:
            enhanced = clint_liver_L_h * clint_multiplier
            steps.append(
                ConversionStep(
                    name="clint_enhancement",
                    value=enhanced,
                    unit="L/h",
                    formula=(
                        f"CLint_liver * multiplier = "
                        f"{clint_liver_L_h:.4g} * {clint_multiplier} = {enhanced:.4g}"
                    ),
                )
            )
            clint_liver_L_h = enhanced

        # ---- Step 3: fu_b (blood unbound fraction) ---------------------
        fu_b = fu_p / bp_ratio
        steps.append(
            ConversionStep(
                name="fu_b",
                value=fu_b,
                unit="fraction",
                formula=f"fu_p/bp_ratio = {fu_p}/{bp_ratio}",
            )
        )

        # ---- Step 4: Apply liver model ---------------------------------
        clh = liver_model_fn(qh=qh_L_h, fu_b=fu_b, clint_liver=clint_liver_L_h)
        steps.append(
            ConversionStep(
                name="CLh",
                value=clh,
                unit="L/h",
                formula=(
                    f"{model}: Qh={qh_L_h}, "
                    f"fu_b=fu_p/BP={fu_p}/{bp_ratio}={fu_b:.6g}, "
                    f"CLint_liver={clint_liver_L_h:.4g}"
                ),
            )
        )

        # ---- Step 5: Extraction ratio ----------------------------------
        extraction_ratio = clh / qh_L_h if qh_L_h > 0 else 0.0
        steps.append(
            ConversionStep(
                name="extraction_ratio",
                value=extraction_ratio,
                unit="fraction",
                formula=f"CLh/Qh = {clh:.4g}/{qh_L_h}",
            )
        )

        # ---- Assemble result -------------------------------------------
        conversion_log = ConversionLog(
            input_params=input_params,
            intermediate_steps=steps,
            output=clh,
            output_unit="L/h",
            model_used=model,
        )

        return HepaticClearance(
            clh_L_h=clh,
            extraction_ratio=extraction_ratio,
            model_used=model,
            conversion_log=conversion_log,
            clint_liver_L_h=clint_liver_L_h,
        )

    # ------------------------------------------------------------------
    # Permeability: Papp -> Peff
    # ------------------------------------------------------------------

    def papp_to_peff(
        self,
        papp_nm_s: float,
        calibration: str = "sun_2002",
    ) -> float:
        """Convert Caco-2 apparent permeability to human effective permeability.

        Uses the Sun et al. 2002 calibration::

            log10(Peff_cm_s) = 0.4926 * log10(Papp_cm_s) - 0.1454

        Parameters
        ----------
        papp_nm_s : float
            Apparent permeability in nm/s (as predicted by ML models).
        calibration : str
            Calibration equation to use.  Currently only ``"sun_2002"``
            is implemented.

        Returns
        -------
        float
            Human effective permeability Peff in cm/s.

        Raises
        ------
        ValueError
            If ``papp_nm_s <= 0`` or unsupported calibration.
        """
        if papp_nm_s <= 0:
            raise ValueError(f"Papp must be positive, got {papp_nm_s}")

        if calibration != "sun_2002":
            raise ValueError(
                f"Unknown calibration {calibration!r}; only 'sun_2002' is supported"
            )

        papp_cm_s = nm_s_to_cm_s(papp_nm_s)
        log_peff = 0.4926 * math.log10(papp_cm_s) - 0.1454
        return 10 ** log_peff

    # ------------------------------------------------------------------
    # Renal clearance
    # ------------------------------------------------------------------

    def assign_renal_clearance(
        self,
        fu_p: float,
        gfr_mL_min: float = HUMAN_GFR_ML_MIN,
        is_active_secretion: bool = False,
        net_secretion_factor: float = 0.0,
    ) -> float:
        """Estimate renal clearance from glomerular filtration.

        Formula::

            CLrenal = fu_p * GFR_L_h * (1 + net_secretion_factor)

        Parameters
        ----------
        fu_p : float
            Plasma unbound fraction in [0, 1].
        gfr_mL_min : float
            Glomerular filtration rate in mL/min (default 120).
        is_active_secretion : bool
            Flag indicating active tubular secretion.
        net_secretion_factor : float
            Multiplicative factor for net secretion (0 = filtration only).

        Returns
        -------
        float
            Renal clearance in L/h.

        Raises
        ------
        ValueError
            If ``fu_p`` is outside [0, 1].
        """
        if fu_p < 0 or fu_p > 1:
            raise ValueError(f"fu_p must be in [0, 1], got {fu_p}")

        gfr_L_h = mL_min_to_L_h(gfr_mL_min)
        cl_renal = fu_p * gfr_L_h * (1.0 + net_secretion_factor)

        if is_active_secretion and net_secretion_factor == 0.0:
            logger.warning(
                "Active secretion flagged but net_secretion_factor is 0.0"
            )

        return cl_renal
