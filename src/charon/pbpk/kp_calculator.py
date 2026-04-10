"""Tissue:plasma partition coefficient (Kp) calculation.

Implements Rodgers & Rowland (2005/2006) and Poulin & Theil (2000/2002)
mechanistic Kp models with an optional Berezhkovskiy (2004) correction
for highly protein-bound drugs.

Ported verbatim from
``Sisyphus/src/sisyphus/predict/ivive.py:_compute_kp_rodgers_rowland``
and ``_compute_kp_poulin_theil``. The port keeps the formulas identical
(verified to ±1e-10) while adapting the types to Charon's schema:

- ``TissueComposition`` is a small frozen dataclass (fn, fp, fw, pH).
- Functions return plain ``float`` Kp values (not Sisyphus Distribution).
- CLint decomposition lives in Sprint 3, not here.

References:
    Rodgers T, Leahy D, Rowland M (2005). "Physiologically based
        pharmacokinetic modeling 1: predicting the tissue distribution
        of moderate-to-strong bases." J Pharm Sci 94:1259-76.
    Rodgers T, Rowland M (2006). "Physiologically based pharmacokinetic
        modelling 2: predicting the tissue distribution of acids, very
        weak bases, neutrals and zwitterions." J Pharm Sci 95:1238-57.
    Poulin P, Theil F-P (2002). "Prediction of pharmacokinetics prior
        to in vivo studies. 1. Mechanism-based prediction of volume of
        distribution." J Pharm Sci 91:129-56.
    Berezhkovskiy LM (2004). "Volume of distribution at steady state
        for a linear pharmacokinetic system with peripheral elimination."
        J Pharm Sci 93:1628-40.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Callable

KP_MIN = 0.01
KP_MAX = 50.0

_VALID_COMPOUND_TYPES = frozenset({"neutral", "acid", "base", "zwitterion"})
_VALID_METHODS = frozenset({"rodgers_rowland", "poulin_theil", "berezhkovskiy"})


@dataclass(frozen=True)
class TissueComposition:
    """Volume fractions and intracellular pH for a tissue.

    Attributes:
        fn: Neutral lipid fraction.
        fp: Phospholipid fraction.
        fw: Water fraction.
        pH: Intracellular pH (or plasma pH for the plasma reference).
    """

    fn: float
    fp: float
    fw: float
    pH: float


def _lipid_partition(comp: TissueComposition, partition_coeff: float) -> float:
    """Rodgers & Rowland lipid partitioning: fn*P + fp*(0.3P + 0.7)."""
    return comp.fn * partition_coeff + comp.fp * (0.3 * partition_coeff + 0.7)


def _rr_neutral(
    p: float,
    tissue_comp: TissueComposition,
    plasma_comp: TissueComposition,
) -> float:
    """R&R / P&T neutral or unknown-pKa branch."""
    kp_num = tissue_comp.fw + _lipid_partition(tissue_comp, p)
    kp_den = plasma_comp.fw + _lipid_partition(plasma_comp, p)
    return kp_num / kp_den if kp_den > 0 else 1.0


def _rr_acid(
    p: float,
    pka: float,
    tissue_comp: TissueComposition,
    plasma_comp: TissueComposition,
) -> float:
    """R&R / P&T acid branch (distribution coefficient)."""
    d_tissue = p / (1.0 + 10.0 ** (tissue_comp.pH - pka))
    d_plasma = p / (1.0 + 10.0 ** (plasma_comp.pH - pka))
    ion_ratio = (
        (1.0 + 10.0 ** (tissue_comp.pH - pka))
        / (1.0 + 10.0 ** (plasma_comp.pH - pka))
    )
    kp_num = tissue_comp.fw * ion_ratio + _lipid_partition(tissue_comp, d_tissue)
    kp_den = plasma_comp.fw + _lipid_partition(plasma_comp, d_plasma)
    return kp_num / kp_den if kp_den > 0 else 1.0


def _rr_base_or_zwitterion(
    p: float,
    pka: float,
    tissue_comp: TissueComposition,
    plasma_comp: TissueComposition,
) -> float:
    """R&R base / zwitterion branch with phospholipid binding term."""
    ion_ratio = (
        (1.0 + 10.0 ** (pka - tissue_comp.pH))
        / (1.0 + 10.0 ** (pka - plasma_comp.pH))
    )
    phospholipid_binding = tissue_comp.fp * max(ion_ratio - 1.0, 0.0)
    kp_num = (
        tissue_comp.fw * ion_ratio
        + _lipid_partition(tissue_comp, p)
        + phospholipid_binding
    )
    kp_den = plasma_comp.fw + _lipid_partition(plasma_comp, p)
    return kp_num / kp_den if kp_den > 0 else 1.0


def _pt_base_or_zwitterion(
    p: float,
    pka: float,
    tissue_comp: TissueComposition,
    plasma_comp: TissueComposition,
) -> float:
    """P&T base / zwitterion branch WITHOUT phospholipid binding.

    The only mechanistic difference between Rodgers & Rowland and
    Poulin & Theil for bases/zwitterions: P&T drops the phospholipid
    binding term, which gives more conservative (lower) Kp values.
    """
    ion_ratio = (
        (1.0 + 10.0 ** (pka - tissue_comp.pH))
        / (1.0 + 10.0 ** (pka - plasma_comp.pH))
    )
    kp_num = tissue_comp.fw * ion_ratio + _lipid_partition(tissue_comp, p)
    kp_den = plasma_comp.fw + _lipid_partition(plasma_comp, p)
    return kp_num / kp_den if kp_den > 0 else 1.0


def _validate_inputs(
    logp: float,
    compound_type: str,
    tissue_comp: TissueComposition,
    plasma_comp: TissueComposition,
) -> None:
    if not math.isfinite(logp):
        raise ValueError(f"logp must be finite, got {logp}")
    if compound_type not in _VALID_COMPOUND_TYPES:
        raise ValueError(
            f"compound_type must be one of {sorted(_VALID_COMPOUND_TYPES)}, "
            f"got {compound_type!r}"
        )
    for label, comp in (("tissue", tissue_comp), ("plasma", plasma_comp)):
        for attr in ("fn", "fp", "fw", "pH"):
            value = getattr(comp, attr)
            if not math.isfinite(value):
                raise ValueError(f"{label}.{attr} must be finite, got {value}")


def compute_kp_rodgers_rowland(
    logp: float,
    pka: float | None,
    compound_type: str,
    tissue_comp: TissueComposition,
    plasma_comp: TissueComposition,
) -> float:
    """Compute a single tissue Kp via Rodgers & Rowland.

    Neutrals and pka=None → neutral branch. Acids use the distribution
    coefficient. Bases and zwitterions add the phospholipid binding
    term. Result is clipped to ``[KP_MIN, KP_MAX]``.
    """
    _validate_inputs(logp, compound_type, tissue_comp, plasma_comp)
    p = 10.0 ** logp

    if compound_type == "neutral" or pka is None:
        kp = _rr_neutral(p, tissue_comp, plasma_comp)
    elif compound_type == "acid":
        kp = _rr_acid(p, pka, tissue_comp, plasma_comp)
    else:  # "base" or "zwitterion"
        kp = _rr_base_or_zwitterion(p, pka, tissue_comp, plasma_comp)

    return float(max(KP_MIN, min(KP_MAX, kp)))


def compute_kp_poulin_theil(
    logp: float,
    pka: float | None,
    compound_type: str,
    tissue_comp: TissueComposition,
    plasma_comp: TissueComposition,
) -> float:
    """Compute a single tissue Kp via Poulin & Theil.

    Identical to R&R for neutrals and acids. For bases/zwitterions,
    omits the phospholipid binding term — the key mechanistic difference
    — giving more conservative Kp values.
    """
    _validate_inputs(logp, compound_type, tissue_comp, plasma_comp)
    p = 10.0 ** logp

    if compound_type == "neutral" or pka is None:
        kp = _rr_neutral(p, tissue_comp, plasma_comp)
    elif compound_type == "acid":
        kp = _rr_acid(p, pka, tissue_comp, plasma_comp)
    else:  # "base" or "zwitterion"
        kp = _pt_base_or_zwitterion(p, pka, tissue_comp, plasma_comp)

    return float(max(KP_MIN, min(KP_MAX, kp)))


def apply_berezhkovskiy_correction(kp: float, fu_p: float) -> float:
    """Berezhkovskiy (2004) plasma-protein-binding correction.

    Kp_bz = Kp_rr / (1 + (Kp_rr - 1) * fu_p)

    Accounts for the fact that only the unbound drug distributes. For
    highly bound drugs (fu_p ≪ 1) this meaningfully reduces Kp and
    VDss. For fully unbound drugs (fu_p = 1) the correction is a no-op.

    The formula is numerically safe for Kp ≥ 0 and fu_p ∈ [0, 1]. If
    the denominator collapses (should not happen under the valid input
    domain) the uncorrected Kp is returned unchanged.
    """
    if not math.isfinite(kp) or kp < 0.0:
        raise ValueError(f"kp must be non-negative finite, got {kp}")
    if not math.isfinite(fu_p) or fu_p < 0.0 or fu_p > 1.0:
        raise ValueError(f"fu_p must be in [0, 1], got {fu_p}")
    denom = 1.0 + (kp - 1.0) * fu_p
    if denom < 1e-10:
        return kp
    return kp / denom


_KP_METHOD_DISPATCH: dict[str, Callable[..., float]] = {
    "rodgers_rowland": compute_kp_rodgers_rowland,
    "poulin_theil": compute_kp_poulin_theil,
}


def compute_all_kp(
    logp: float,
    pka: float | None,
    compound_type: str,
    tissue_compositions: dict[str, TissueComposition],
    plasma_composition: TissueComposition,
    method: str = "rodgers_rowland",
    fu_p: float | None = None,
) -> dict[str, float]:
    """Compute Kp for every tissue in ``tissue_compositions``.

    Args:
        logp: Octanol-water partition coefficient.
        pka: Dissociation constant, or ``None`` for neutral-only handling.
        compound_type: One of ``"neutral"``, ``"acid"``, ``"base"``,
            ``"zwitterion"``.
        tissue_compositions: Mapping of tissue name → composition.
        plasma_composition: Reference plasma composition.
        method: ``"rodgers_rowland"``, ``"poulin_theil"``, or
            ``"berezhkovskiy"`` (R&R + BZ correction).
        fu_p: Plasma unbound fraction. Required when ``method =
            "berezhkovskiy"``; otherwise ignored.

    Returns:
        Dict of ``{tissue_name: kp_value}``. Values are clipped to
        ``[KP_MIN, KP_MAX]``.
    """
    if method not in _VALID_METHODS:
        raise ValueError(
            f"method must be one of {sorted(_VALID_METHODS)}, got {method!r}"
        )
    use_bz = method == "berezhkovskiy"
    if use_bz and fu_p is None:
        raise ValueError("fu_p is required when method='berezhkovskiy'")

    base_method = "rodgers_rowland" if use_bz else method
    kp_fn = _KP_METHOD_DISPATCH[base_method]

    results: dict[str, float] = {}
    for tissue_name, tissue_comp in tissue_compositions.items():
        kp = kp_fn(logp, pka, compound_type, tissue_comp, plasma_composition)
        if use_bz:
            assert fu_p is not None  # narrow type for mypy
            kp = apply_berezhkovskiy_correction(kp, fu_p)
            kp = float(max(KP_MIN, min(KP_MAX, kp)))
        results[tissue_name] = kp
    return results
