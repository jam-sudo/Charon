"""Sprint 10 research-only: fold-error decomposition for Layer 3 Tier A panel.

No production wiring — `charon.pipeline.Pipeline` does not import this module.
Everything is pure-function; I/O is the orchestrator's job
(`validation/benchmarks/layer3_ivive_decomposition.py`).

Signed log-additive decomposition:

    log10(fold_observed_signed) = log10(fold_liver_model_signed)
                                + log10(fold_route_bias)
                                + log10(fold_residual_signed)

where:

- `fold_observed_signed = mrsd_ws / fih_reference_mg` — the production
  (well_stirred) prediction relative to the clinical reference dose, in
  signed form (<1 means under-predicted, >1 means over-predicted).

- `fold_liver_model_signed = mrsd_ws / mrsd_best_alt` where best_alt is
  whichever of {parallel_tube, dispersion} minimises |log10(mrsd/reference)|.
  If neither alternate is closer to the reference than well_stirred, the
  factor is exactly 1.0 (no attributable improvement).

- `fold_route_bias = 1 / F_literature` for oral-reference compounds (a
  route-comparison artefact, not an IVIVE error — the production pipeline
  predicts IV MRSD and the reference is an oral dose). 1.0 for IV
  references or when F is unknown (flagged). Always >= 1.

- `fold_residual_signed = fold_observed_signed / (liver * route_bias)` —
  the unexplained remainder (transporter / non-hepatic / UGT / model-gap).

The helper `to_symmetric(signed)` converts any signed factor to the
benchmark-style symmetric fold `max(x, 1/x)` for reporting.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Literal


@dataclass(frozen=True)
class DecompositionResult:
    """Result of a single-compound fold-error decomposition (signed factors).

    Attributes:
        fold_observed_signed: mrsd_ws / fih_reference_mg.
        fold_liver_model_signed: mrsd_ws / mrsd_best_alt (1.0 if no
            alternate improves the prediction).
        fold_route_bias: 1/F_literature for oral references, 1.0 for IV
            or unknown-F.
        fold_residual_signed: observed / (liver * route_bias).
        best_alt_model_name: "parallel_tube", "dispersion", or "well_stirred".
        flags: list of string flags, e.g. ["f_unknown"].
    """

    fold_observed_signed: float
    fold_liver_model_signed: float
    fold_route_bias: float
    fold_residual_signed: float
    best_alt_model_name: str
    flags: tuple[str, ...] = field(default_factory=tuple)


def to_symmetric(signed: float) -> float:
    """Convert a signed fold factor to the benchmark-style symmetric form.

    max(x, 1/x); 1.0 for x == 1.0. Returns math.inf for non-positive input.
    """
    if signed <= 0:
        return math.inf
    return max(signed, 1.0 / signed)


def select_best_alternate_liver_model(
    ws: float,
    pt: float,
    disp: float,
    reference: float,
) -> tuple[str, float]:
    """Return (name, mrsd) of the liver model closest to the reference.

    Closeness = minimum |log10(mrsd / reference)|. If neither alternate
    (parallel_tube, dispersion) is closer to reference than well_stirred,
    returns ("well_stirred", ws) — meaning no attributable improvement
    from liver-model choice.

    Tie-break when alternates equally close: `parallel_tube` wins (it is listed first in a stable-sort, not alphabetical order).
    """
    if ws <= 0:
        raise ValueError(f"ws must be > 0, got {ws}")
    if reference <= 0:
        raise ValueError(f"reference must be > 0, got {reference}")
    candidates = {
        "well_stirred": ws,
        "parallel_tube": pt,
        "dispersion": disp,
    }
    distances = {
        name: abs(math.log10(v / reference)) if v > 0 else math.inf
        for name, v in candidates.items()
    }
    ws_dist = distances["well_stirred"]
    alt_sorted = sorted(
        [("parallel_tube", distances["parallel_tube"]),
         ("dispersion", distances["dispersion"])],
        key=lambda kv: kv[1],
    )
    best_alt_name, best_alt_dist = alt_sorted[0]
    if best_alt_dist < ws_dist:
        return best_alt_name, candidates[best_alt_name]
    return "well_stirred", ws


def compute_route_bias_factor(
    route_ref: Literal["iv", "oral"],
    f_lit: float | None,
) -> tuple[float, tuple[str, ...]]:
    """Return (factor, flags) for the 1/F route-mismatch attribution.

    IV reference → (1.0, ()). Oral + known F → (1/F, ()). Oral + unknown F
    → (1.0, ("f_unknown",)).
    """
    if route_ref not in ("iv", "oral"):
        raise ValueError(
            f"route_ref must be 'iv' or 'oral', got {route_ref!r}"
        )
    if route_ref == "iv":
        return 1.0, ()
    if f_lit is None:
        return 1.0, ("f_unknown",)
    if f_lit <= 0 or f_lit > 1.0:
        raise ValueError(f"f_lit must be in (0, 1], got {f_lit}")
    return 1.0 / f_lit, ()


def decompose_fold_error(
    mrsd_ws: float,
    mrsd_pt: float,
    mrsd_disp: float,
    f_lit: float | None,
    route_ref: Literal["iv", "oral"],
    fih_reference_mg: float,
) -> DecompositionResult:
    """Decompose Tier A observed fold-error into three multiplicative factors (signed).

    Args:
        mrsd_ws: MRSD predicted with well_stirred liver model (production default).
        mrsd_pt: MRSD predicted with parallel_tube.
        mrsd_disp: MRSD predicted with dispersion.
        f_lit: Literature bioavailability (0 < F <= 1) for oral route_ref;
            None for IV or unknown.
        route_ref: 'iv' or 'oral' — the route of the clinical FIH reference
            dose (NOT the simulation route, which is always iv_bolus per
            Sprint 7's panel.yaml).
        fih_reference_mg: Reference FIH dose in mg.

    Returns:
        DecompositionResult with signed factors satisfying:
            log10(fold_observed_signed) ==
                log10(fold_liver_model_signed)
              + log10(fold_route_bias)
              + log10(fold_residual_signed)
    """
    if mrsd_ws <= 0 or fih_reference_mg <= 0:
        raise ValueError(
            f"mrsd_ws and fih_reference_mg must be > 0, "
            f"got {mrsd_ws}, {fih_reference_mg}"
        )
    if mrsd_pt <= 0 or mrsd_disp <= 0:
        raise ValueError(
            f"mrsd_pt and mrsd_disp must be > 0, got {mrsd_pt}, {mrsd_disp}"
        )

    fold_obs_signed = mrsd_ws / fih_reference_mg

    best_alt_name, best_alt_mrsd = select_best_alternate_liver_model(
        ws=mrsd_ws, pt=mrsd_pt, disp=mrsd_disp, reference=fih_reference_mg
    )
    if best_alt_name == "well_stirred":
        fold_liver_signed = 1.0
    else:
        fold_liver_signed = mrsd_ws / best_alt_mrsd

    fold_route, flags = compute_route_bias_factor(
        route_ref=route_ref, f_lit=f_lit
    )

    fold_residual_signed = fold_obs_signed / (fold_liver_signed * fold_route)

    return DecompositionResult(
        fold_observed_signed=fold_obs_signed,
        fold_liver_model_signed=fold_liver_signed,
        fold_route_bias=fold_route,
        fold_residual_signed=fold_residual_signed,
        best_alt_model_name=best_alt_name,
        flags=tuple(flags),
    )
