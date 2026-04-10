"""Extract PK parameters (Cmax, AUC, CL, Vss, t_half) from a Cp-time profile."""

from __future__ import annotations

import math
from typing import Literal

import numpy as np

from charon.core.schema import PKParameters


def _trapezoid(x: np.ndarray, y: np.ndarray) -> float:
    """Trapezoidal rule integral (numpy 1.x / 2.x compatible)."""
    trap = getattr(np, "trapezoid", None) or getattr(np, "trapz")
    return float(trap(y, x))


def _terminal_log_slope(
    t: np.ndarray, cp: np.ndarray, frac: float = 0.3, min_points: int = 3
) -> tuple[float, float]:
    """Log-linear regression on the last ``frac`` of the profile.

    Returns ``(ke, cp_last)``.  Raises ValueError if the slope is
    non-negative (drug not declining → terminal phase undetectable).
    """
    n = len(t)
    n_tail = max(min_points, int(np.ceil(frac * n)))
    n_tail = min(n_tail, n)
    t_tail = t[-n_tail:]
    cp_tail = cp[-n_tail:]
    mask = cp_tail > 0
    if mask.sum() < 2:
        raise ValueError("terminal slope cannot be fit: fewer than 2 positive samples")
    t_fit = t_tail[mask]
    log_cp = np.log(cp_tail[mask])
    slope, _intercept = np.polyfit(t_fit, log_cp, 1)
    if slope >= 0:
        raise ValueError(
            f"terminal slope is non-negative ({slope:.3e}); simulation "
            f"too short or profile not declining"
        )
    ke = -float(slope)
    return ke, float(cp_tail[mask][-1])


def compute_pk_parameters(
    time_h: np.ndarray,
    cp_plasma: np.ndarray,
    *,
    dose_mg: float,
    route: Literal["iv_bolus", "iv_infusion"],
    infusion_duration_h: float = 0.0,
) -> PKParameters:
    """Compute standard PK parameters from a plasma concentration profile.

    Parameters
    ----------
    time_h : array
        Time grid in hours (monotonically increasing).
    cp_plasma : array
        Plasma concentration at each time point (mg/L).
    dose_mg : float
        Administered dose (mg).
    route : "iv_bolus" or "iv_infusion"
        For "iv_infusion" Vss is corrected via Benet (Vss = CL*(MRT - T_inf/2)).
    infusion_duration_h : float
        Infusion duration (only used when route="iv_infusion").
    """
    t = np.asarray(time_h, dtype=np.float64)
    cp = np.asarray(cp_plasma, dtype=np.float64)
    if t.shape != cp.shape:
        raise ValueError(
            f"time and cp must have the same length, got {t.shape} vs {cp.shape}"
        )
    if t.ndim != 1:
        raise ValueError("time and cp must be 1-D arrays")
    if dose_mg <= 0:
        raise ValueError(f"dose_mg must be > 0, got {dose_mg}")
    if route not in ("iv_bolus", "iv_infusion"):
        raise ValueError(f"route must be iv_bolus or iv_infusion, got {route!r}")

    cmax = float(np.max(cp))
    tmax = float(t[int(np.argmax(cp))])

    auc_0_last = _trapezoid(t, cp)
    mask_24 = t <= 24.0
    if mask_24.sum() >= 2:
        auc_0_24 = _trapezoid(t[mask_24], cp[mask_24])
    else:
        auc_0_24 = None

    ke, cp_last = _terminal_log_slope(t, cp)
    auc_tail = cp_last / ke
    auc_0_inf = auc_0_last + auc_tail

    aumc_0_last = _trapezoid(t, t * cp)
    t_last = float(t[-1])
    aumc_tail = cp_last * (t_last / ke + 1.0 / (ke * ke))
    aumc_0_inf = aumc_0_last + aumc_tail

    half_life = math.log(2) / ke
    cl = dose_mg / auc_0_inf
    mrt = aumc_0_inf / auc_0_inf

    if route == "iv_bolus":
        vss = cl * mrt
    else:
        vss = cl * (mrt - infusion_duration_h / 2.0)

    return PKParameters(
        cmax=cmax,
        tmax=tmax,
        auc_0_inf=float(auc_0_inf),
        auc_0_24=None if auc_0_24 is None else float(auc_0_24),
        half_life=float(half_life),
        cl_apparent=float(cl),
        vss=float(vss),
        bioavailability=1.0,
        fa=1.0,
        fg=1.0,
        fh=None,
    )
