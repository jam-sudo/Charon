"""Extract PK parameters (Cmax, AUC, CL, Vss, t_half) from a Cp-time profile."""

from __future__ import annotations

import math
from typing import Literal

import numpy as np

from charon.core.schema import PKParameters
from charon.pbpk.acat import compute_absorption_rates


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


def compute_oral_pk_parameters(
    sim,  # OralSimulationResult (not typed to avoid circular import)
    params,  # OralPBPKParams
    topology,  # PBPKTopology
    *,
    dose_mg: float,
) -> PKParameters:
    """Extract PK parameters and Fa/Fg/Fh from an oral simulation.

    Fa and Fg are computed via post-hoc trapezoidal integration of
    ODE fluxes, independent of hepatic CLint.

    Parameters
    ----------
    sim : OralSimulationResult
        Completed oral simulation.
    params : OralPBPKParams
        Oral compound parameters.
    topology : PBPKTopology
        Species topology (for Qh, used in Fh calculation).
    dose_mg : float
        Administered oral dose (mg).

    Returns
    -------
    PKParameters
        Standard PK parameters with Fa, Fg, Fh, and bioavailability.

    Notes
    -----
    Units:
    - dose_mg: mg
    - Cp: mg/L
    - AUC: mg*h/L
    - CL/F: L/h
    - Fa, Fg, Fh: dimensionless fractions (0-1)

    Fg is derived from the ratio of portal venous flux to total absorption
    flux, computed via trapezoidal integration of the ODE state trajectories.
    This makes Fg independent of hepatic CLint error (post-hoc decomposition).

    Fh uses the analytical well-stirred formula: Fh = Qh / (Qh + fu_b * CLint_liver).
    """
    t = np.asarray(sim.time_h, dtype=np.float64)
    cp = np.asarray(sim.cp_plasma, dtype=np.float64)

    # --- Standard PK ---
    cmax = float(np.max(cp))
    tmax = float(t[int(np.argmax(cp))])

    auc_0_last = _trapezoid(t, cp)

    mask_24 = t <= 24.0
    auc_0_24 = _trapezoid(t[mask_24], cp[mask_24]) if mask_24.sum() >= 2 else None

    ke, cp_last = _terminal_log_slope(t, cp)
    auc_tail = cp_last / ke
    auc_0_inf = auc_0_last + auc_tail

    half_life = math.log(2) / ke
    cl_apparent = dose_mg / auc_0_inf  # CL/F for oral

    # --- Fa: fraction absorbed (trapezoidal integration of absorption flux) ---
    gi = params.gi_tract
    k_abs_arr = np.array(
        compute_absorption_rates(gi, params.peff_cm_s),
        dtype=np.float64,
    )
    # Absorption flux at each time point: sum(k_abs_i * A_lumen_i(t))
    absorption_flux = np.zeros_like(t)
    for i in range(len(gi.segments)):
        absorption_flux += k_abs_arr[i] * sim.lumen_trajectory[i, :]
    absorbed_total = _trapezoid(t, absorption_flux)
    fa = absorbed_total / dose_mg

    # --- Fg: fraction escaping gut metabolism (independent of CLint_liver) ---
    k_baso = params.q_villi_L_h / params.v_enterocyte_L
    portal_flux = k_baso * sim.enterocyte_trajectory
    portal_total = _trapezoid(t, portal_flux)
    fg = portal_total / absorbed_total if absorbed_total > 0 else 1.0

    # --- Fh: fraction escaping hepatic first-pass (analytical well-stirred) ---
    qh = topology.tissues["liver"].blood_flow_L_h
    fh = qh / (qh + params.fu_b * params.clint_liver_L_h)

    # --- Bioavailability ---
    bioavailability = fa * fg * fh

    return PKParameters(
        cmax=cmax,
        tmax=tmax,
        auc_0_inf=float(auc_0_inf),
        auc_0_24=float(auc_0_24) if auc_0_24 is not None else None,
        half_life=float(half_life),
        cl_apparent=float(cl_apparent),
        vss=None,  # not meaningful for oral
        bioavailability=float(bioavailability),
        fa=float(fa),
        fg=float(fg),
        fh=float(fh),
    )
