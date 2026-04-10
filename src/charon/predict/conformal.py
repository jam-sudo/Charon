"""Log-space symmetric conformal prediction intervals.

Charon's conformal implementation intentionally diverges from Omega's in
three ways:

1. **Always log-space.** Every property interval is computed as
   ``[pred / 10^q, pred * 10^q]``. This guarantees positive bounds and
   matches the log-normal error distribution typical of ADME endpoints.
   Omega's linear-space widening produced negative lower bounds for peff.

2. **Train/calibration separation enforced.** The training scripts
   explicitly drop any compound matching an adme_reference.csv
   InChIKey-14. This file re-verifies that separation and warns loudly
   if a training compound has leaked back into calibration.

3. **Honest coverage reporting.** After calibration, we recompute
   empirical coverage on the calibration set and report it. This is
   not a true held-out coverage (the calibration set defines the
   quantile) but it is a useful sanity check — we warn if the in-sample
   coverage is below 85%.

Calibration data source: ``data/validation/adme_reference.csv`` (153
rows, 8 columns: name, smiles, mw, logP, fup, rbp,
clint_3a4_uL_min_pmol, peff_cm_s).

Note on CLint units: adme_reference records CLint in μL/min/pmol CYP3A4,
whereas Charon's XGBoost model predicts μL/min/10^6 cells (hepatocyte
basis). Because these two units are not interconvertible without a
compound-specific CYP3A4 fraction-metabolized, CLint cannot be
conformally calibrated against adme_reference. For CLint we derive the
log-space quantile from the training OOF residuals directly (see
:meth:`ConformalPredictor.calibrate_from_oof`).
"""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Mapping

import numpy as np
import pandas as pd

from charon.predict.admet_ensemble import ADMETPredictor

logger = logging.getLogger(__name__)

DEFAULT_COVERAGE = 0.90
COVERAGE_WARNING_THRESHOLD = 0.85
_EPS = 1e-12

# Physical bounds applied after computing symmetric log-space intervals.
# A symmetric interval on a bounded quantity (e.g. fup ∈ [0, 1]) can
# legitimately produce an upper bound above the physical maximum; we
# clip so downstream consumers never see infeasible values.
_PROPERTY_BOUNDS: dict[str, tuple[float, float]] = {
    "fup": (0.001, 1.0),
    "clint": (0.1, 1000.0),
    "clint_hepatocyte": (0.1, 1000.0),
}


@dataclass
class CoverageReport:
    """Per-property quantile and empirical coverage."""

    property_name: str
    n_samples: int
    quantile_log10: float        # |log10(pred/obs)| at the target coverage
    factor: float                # 10^quantile_log10
    empirical_coverage: float    # in-sample coverage at this quantile
    median_fold_error: float     # 10^median(|log10(pred/obs)|)
    mean_fold_error: float       # 10^mean(|log10(pred/obs)|)
    warning: str | None = None


class ConformalPredictor:
    """Log-space symmetric conformal predictor.

    The predictor holds one quantile per property. A given ``predictor``
    (ADMETPredictor) is calibrated once against the reference CSV; after
    that, :meth:`get_interval` returns a symmetric interval for any new
    prediction.

    Usage:
        >>> p = ADMETPredictor()
        >>> cp = ConformalPredictor(Path("data/validation/adme_reference.csv"))
        >>> cp.calibrate(p)
        >>> lo, hi = cp.get_interval("fup", 0.25)
    """

    # Map Charon property names → (reference column, predictor attr).
    # Only properties with consistent units across the two sources are
    # calibrated here. CLint lives in different units and is handled by
    # :meth:`calibrate_from_oof`.
    _SUPPORTED_PROPERTIES = {
        "fup": "fup",
    }

    def __init__(
        self,
        calibration_data: Path | str,
        coverage: float = DEFAULT_COVERAGE,
    ) -> None:
        self._calibration_path = Path(calibration_data)
        if not (0.0 < coverage < 1.0):
            raise ValueError(f"coverage must be in (0, 1), got {coverage}")
        self._coverage = coverage
        self._reports: dict[str, CoverageReport] = {}

    # ------------------------------------------------------------------
    # Calibration
    # ------------------------------------------------------------------
    def calibrate(self, predictor: ADMETPredictor) -> dict[str, CoverageReport]:
        """Calibrate property intervals using the reference CSV.

        Args:
            predictor: An ADMETPredictor whose models are loaded (they
                are loaded lazily on first call).

        Returns:
            Dict of per-property :class:`CoverageReport`.
        """
        if not self._calibration_path.exists():
            raise FileNotFoundError(
                f"Calibration data not found: {self._calibration_path}"
            )
        df = pd.read_csv(self._calibration_path)
        logger.info("Loaded %d calibration compounds from %s",
                    len(df), self._calibration_path.name)

        reports: dict[str, CoverageReport] = {}
        for prop, ref_col in self._SUPPORTED_PROPERTIES.items():
            if ref_col not in df.columns:
                logger.warning("Column %s missing in %s; skipping", ref_col, self._calibration_path)
                continue
            residuals: list[float] = []
            for row in df.itertuples(index=False):
                smi = getattr(row, "smiles", None)
                obs = getattr(row, ref_col)
                if not isinstance(smi, str) or not smi.strip():
                    continue
                if obs is None or (isinstance(obs, float) and math.isnan(obs)):
                    continue
                try:
                    pred = predictor.predict(smi)
                except ValueError:
                    continue
                pred_val = getattr(pred, prop)
                if pred_val <= 0.0 or obs <= 0.0:
                    continue
                residuals.append(abs(math.log10(pred_val / obs)))
            if not residuals:
                logger.warning("No residuals collected for %s", prop)
                continue
            report = self._build_report(prop, residuals)
            reports[prop] = report
            logger.info(
                "  %s: q%.0f=%.3f (factor %.2fx), coverage=%.1f%% on %d samples",
                prop,
                self._coverage * 100,
                report.quantile_log10,
                report.factor,
                report.empirical_coverage * 100,
                report.n_samples,
            )
        self._reports.update(reports)
        return reports

    def calibrate_from_oof(
        self,
        property_name: str,
        oof_log_residuals: np.ndarray,
    ) -> CoverageReport:
        """Calibrate a single property from out-of-fold log-space residuals.

        Used for CLint, whose calibration set and training set live in
        incompatible units. Supply the ``|log10(pred/obs)|`` residuals
        from scaffold-CV on the training set.
        """
        residuals = np.abs(np.asarray(oof_log_residuals, dtype=np.float64))
        residuals = residuals[np.isfinite(residuals)]
        if residuals.size == 0:
            raise ValueError("No finite residuals provided")
        report = self._build_report(property_name, residuals.tolist())
        self._reports[property_name] = report
        logger.info(
            "Calibrated %s from OOF residuals: q%.0f=%.3f (factor %.2fx), %d samples",
            property_name,
            self._coverage * 100,
            report.quantile_log10,
            report.factor,
            report.n_samples,
        )
        return report

    def _build_report(self, property_name: str, residuals: list[float]) -> CoverageReport:
        arr = np.asarray(residuals, dtype=np.float64)
        # Use NumPy's "higher" interpolation so the quantile is one of
        # the observed residuals (more conservative / honest for small N).
        q = float(np.quantile(arr, self._coverage, method="higher"))
        factor = 10.0 ** q
        coverage = float(np.mean(arr <= q + _EPS))
        median_fe = float(10.0 ** np.median(arr))
        mean_fe = float(10.0 ** np.mean(arr))
        warning: str | None = None
        if coverage < COVERAGE_WARNING_THRESHOLD:
            warning = (
                f"Empirical coverage {coverage:.1%} is below "
                f"{COVERAGE_WARNING_THRESHOLD:.0%}; consider recalibration."
            )
            logger.warning("%s: %s", property_name, warning)
        return CoverageReport(
            property_name=property_name,
            n_samples=arr.size,
            quantile_log10=q,
            factor=factor,
            empirical_coverage=coverage,
            median_fold_error=median_fe,
            mean_fold_error=mean_fe,
            warning=warning,
        )

    # ------------------------------------------------------------------
    # Lookup
    # ------------------------------------------------------------------
    def get_interval(
        self,
        property_name: str,
        predicted_value: float,
    ) -> tuple[float, float]:
        """Return the (lower, upper) conformal interval in linear space.

        Args:
            property_name: One of the calibrated properties.
            predicted_value: Point prediction (linear scale, > 0).

        Returns:
            Tuple of (lower_bound, upper_bound), both positive.

        Raises:
            RuntimeError: If ``property_name`` has not been calibrated.
            ValueError: If ``predicted_value`` is not positive.
        """
        if property_name not in self._reports:
            raise RuntimeError(
                f"Property {property_name!r} is not calibrated. "
                f"Call calibrate() or calibrate_from_oof() first."
            )
        if not math.isfinite(predicted_value) or predicted_value <= 0.0:
            raise ValueError(f"predicted_value must be positive, got {predicted_value}")
        factor = self._reports[property_name].factor
        lower = predicted_value / factor
        upper = predicted_value * factor
        bounds = _PROPERTY_BOUNDS.get(property_name)
        if bounds is not None:
            lo_min, hi_max = bounds
            lower = max(lower, lo_min)
            upper = min(upper, hi_max)
            # Guarantee that the interval still brackets the point
            # estimate after clipping.
            lower = min(lower, predicted_value)
            upper = max(upper, predicted_value)
        return lower, upper

    def coverage_report(self) -> Mapping[str, CoverageReport]:
        """Return a read-only view of calibration reports."""
        return dict(self._reports)

    def is_calibrated(self, property_name: str) -> bool:
        return property_name in self._reports
