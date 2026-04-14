"""N-sample propagation engine for uncertainty quantification.

Runs N PBPK simulations with sampled parameter sets, collecting dose and
PK outcomes for downstream CI aggregation (dose_range.py) and sensitivity
analysis (sobol.py).

Architecture
------------
For each sample from LHS:
  1. ``override_compound`` creates a modified CompoundConfig with sampled values.
  2. A fresh ``Pipeline`` is constructed and run.
  3. MRSD (mg) and CL_apparent (L/h) are collected from the result.
  4. Failures are logged and skipped — never propagated.

The key mapping from sampling parameter names to CompoundConfig property
paths is:

    logp             -> properties.physicochemical.logp.value
    fu_p             -> properties.binding.fu_p.value
    clint_uL_min_mg  -> properties.metabolism.clint_uL_min_mg.value
    peff_cm_s        -> properties.permeability.peff_cm_s.value
    bp_ratio         -> properties.binding.bp_ratio.value

``mppgl`` is a physiological parameter that does not live in CompoundConfig;
it is passed through the pipeline via the YAML species file default and is
not overridden here (Phase A simplification).
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Literal

import numpy as np

from charon.core.schema import (
    CompoundConfig,
    DoseProjectionConfig,
    PredictedProperty,
)

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Parameter key -> nested CompoundConfig path
# ---------------------------------------------------------------------------

# Each entry: (property_group, field_name)
# e.g. ("physicochemical", "logp") -> properties.physicochemical.logp
_PARAM_PATH: dict[str, tuple[str, str]] = {
    "logp": ("physicochemical", "logp"),
    "fu_p": ("binding", "fu_p"),
    "clint_uL_min_mg": ("binding", "fu_p"),  # placeholder, overridden below
    "peff_cm_s": ("permeability", "peff_cm_s"),
    "bp_ratio": ("binding", "bp_ratio"),
}
# Fix: clint lives under metabolism, not binding
_PARAM_PATH["clint_uL_min_mg"] = ("metabolism", "clint_uL_min_mg")


# ---------------------------------------------------------------------------
# override_compound
# ---------------------------------------------------------------------------


def override_compound(
    base: CompoundConfig,
    sample: dict[str, float],
) -> CompoundConfig:
    """Create a modified CompoundConfig with sampled parameter values.

    For each key in *sample*, updates the corresponding
    ``PredictedProperty.value`` via nested ``model_copy``.  Keys not
    present in *sample* or where the base property is ``None`` are
    skipped.

    Parameters
    ----------
    base : CompoundConfig
        The original compound configuration (not mutated).
    sample : dict[str, float]
        Mapping of parameter name -> sampled value.  Expected keys are
        a subset of: ``logp``, ``fu_p``, ``clint_uL_min_mg``,
        ``peff_cm_s``, ``bp_ratio``.

    Returns
    -------
    CompoundConfig
        A new CompoundConfig with overridden property values.

    Notes
    -----
    ``mppgl`` is a physiological constant (not a compound property) and
    is therefore not handled here.  It would need to be injected at the
    species-topology level in a future sprint.
    """
    # Work through nested Pydantic copies to preserve immutability of the
    # original.  We accumulate changes per property group and apply them
    # in a single model_copy chain.

    props = base.properties

    # Track which groups need updating
    physchem = props.physicochemical
    binding = props.binding
    metabolism = props.metabolism
    permeability = props.permeability

    physchem_changed = False
    binding_changed = False
    metabolism_changed = False
    permeability_changed = False

    for key, value in sample.items():
        if key not in _PARAM_PATH:
            # mppgl or other non-compound params — skip silently
            continue

        group_name, field_name = _PARAM_PATH[key]

        if group_name == "physicochemical":
            prop: PredictedProperty | None = getattr(physchem, field_name, None)
            if prop is not None:
                new_prop = prop.model_copy(update={"value": value})
                physchem = physchem.model_copy(update={field_name: new_prop})
                physchem_changed = True

        elif group_name == "binding":
            prop = getattr(binding, field_name, None)
            if prop is not None:
                new_prop = prop.model_copy(update={"value": value})
                binding = binding.model_copy(update={field_name: new_prop})
                binding_changed = True

        elif group_name == "metabolism":
            prop = getattr(metabolism, field_name, None)
            if prop is not None:
                new_prop = prop.model_copy(update={"value": value})
                metabolism = metabolism.model_copy(update={field_name: new_prop})
                metabolism_changed = True

        elif group_name == "permeability":
            prop = getattr(permeability, field_name, None)
            if prop is not None:
                new_prop = prop.model_copy(update={"value": value})
                permeability = permeability.model_copy(update={field_name: new_prop})
                permeability_changed = True

    # Build updated properties (only copy what changed)
    props_update: dict[str, object] = {}
    if physchem_changed:
        props_update["physicochemical"] = physchem
    if binding_changed:
        props_update["binding"] = binding
    if metabolism_changed:
        props_update["metabolism"] = metabolism
    if permeability_changed:
        props_update["permeability"] = permeability

    if props_update:
        new_props = props.model_copy(update=props_update)
        return base.model_copy(update={"properties": new_props})

    return base.model_copy()


# ---------------------------------------------------------------------------
# PropagationResult
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class PropagationResult:
    """Immutable container for N-sample propagation output.

    Attributes
    ----------
    doses_mg : np.ndarray, shape (n_successful,)
        MRSD values from successful pipeline runs.
    cl_apparent : np.ndarray, shape (n_successful,)
        Apparent clearance values from successful runs.
    parameter_matrix : np.ndarray, shape (n_successful, n_params)
        Sampled parameter values for successful runs (columns match
        *param_names* ordering).
    n_total : int
        Total number of samples attempted.
    n_successful : int
        Number of samples that completed without error.
    n_failed : int
        Number of samples that raised exceptions.
    param_names : tuple[str, ...]
        Ordered parameter names matching columns of *parameter_matrix*.
    """

    doses_mg: np.ndarray
    cl_apparent: np.ndarray
    parameter_matrix: np.ndarray
    n_total: int
    n_successful: int
    n_failed: int
    param_names: tuple[str, ...]


# ---------------------------------------------------------------------------
# propagate
# ---------------------------------------------------------------------------


def propagate(
    *,
    base_compound: CompoundConfig,
    samples: tuple[dict[str, float], ...],
    route: Literal["iv_bolus", "iv_infusion", "oral"],
    dose_mg: float,
    dose_projection: DoseProjectionConfig,
    duration_h: float = 72.0,
    liver_model: str = "well_stirred",
) -> PropagationResult:
    """Run N PBPK simulations with sampled parameters.

    For each sample in *samples*:
      1. Override compound properties with sampled values.
      2. Construct and run a full Pipeline.
      3. Collect MRSD and CL_apparent.
      4. On exception: skip, log warning, increment failure counter.

    Parameters
    ----------
    base_compound : CompoundConfig
        Base compound configuration (template for overrides).
    samples : tuple[dict[str, float], ...]
        Parameter sets from LHS sampling (one dict per sample).
    route : {"iv_bolus", "iv_infusion", "oral"}
        Administration route.
    dose_mg : float
        Dose in mg for each simulation.
    dose_projection : DoseProjectionConfig
        FIH dose projection settings (NOAEL, MABEL, PAD inputs).
    duration_h : float
        Simulation duration in hours.
    liver_model : str
        Liver extraction model name (default "well_stirred").

    Returns
    -------
    PropagationResult
        Collected dose outcomes, CL values, and parameter matrix.
    """
    # Lazy import to avoid circular dependency at module level
    from charon.pipeline import Pipeline

    n_total = len(samples)
    if n_total == 0:
        param_names = tuple(samples[0].keys()) if samples else ()
        return PropagationResult(
            doses_mg=np.array([], dtype=np.float64),
            cl_apparent=np.array([], dtype=np.float64),
            parameter_matrix=np.empty((0, 0), dtype=np.float64),
            n_total=0,
            n_successful=0,
            n_failed=0,
            param_names=param_names,
        )

    # Determine param_names from first sample
    param_names = tuple(samples[0].keys())

    dose_results: list[float] = []
    cl_results: list[float] = []
    param_rows: list[list[float]] = []
    n_failed = 0

    for i, sample in enumerate(samples):
        try:
            compound_i = override_compound(base_compound, sample)
            pipe = Pipeline(
                compound_i,
                route=route,
                dose_mg=dose_mg,
                dose_projection=dose_projection,
                duration_h=duration_h,
                liver_model=liver_model,
            )
            result = pipe.run()

            # Extract MRSD
            mrsd = None
            if result.dose_recommendation is not None:
                mrsd = result.dose_recommendation.mrsd_mg

            if mrsd is None:
                logger.warning(
                    "Sample %d/%d: no dose_recommendation.mrsd_mg — skipping",
                    i + 1,
                    n_total,
                )
                n_failed += 1
                continue

            # Extract CL_apparent (may be None for some routes)
            cl_app = 0.0
            if result.pk_parameters.cl_apparent is not None:
                cl_app = result.pk_parameters.cl_apparent

            dose_results.append(mrsd)
            cl_results.append(cl_app)
            param_rows.append([sample[k] for k in param_names])

        except Exception:
            logger.warning(
                "Sample %d/%d failed with exception", i + 1, n_total,
                exc_info=True,
            )
            n_failed += 1

    n_successful = len(dose_results)

    return PropagationResult(
        doses_mg=np.array(dose_results, dtype=np.float64),
        cl_apparent=np.array(cl_results, dtype=np.float64),
        parameter_matrix=np.array(param_rows, dtype=np.float64).reshape(
            n_successful, len(param_names)
        ),
        n_total=n_total,
        n_successful=n_successful,
        n_failed=n_failed,
        param_names=param_names,
    )
