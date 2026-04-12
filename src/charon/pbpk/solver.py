"""Stiff ODE solver wrapper for PBPK simulation.

The PBPK system is stiff (fast tissue equilibration vs slow elimination),
so this module mandates ``scipy.integrate.solve_ivp(method='BDF')``.
Explicit methods (RK45, RK23, DOP853) are rejected with ValueError.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np
from scipy.integrate import solve_ivp

from charon.pbpk.ode_compiler import CompoundPBPKParams, OralPBPKParams, build_rhs, build_oral_rhs
from charon.pbpk.topology import PBPKTopology

_ALLOWED_METHODS = {"BDF"}


@dataclass
class SimulationResult:
    """Output of an IV PBPK simulation."""

    time_h: np.ndarray
    cp_blood: np.ndarray
    cp_plasma: np.ndarray
    state_trajectory: np.ndarray
    mass_balance_residual: float
    solver_success: bool
    solver_method: str
    solver_nfev: int
    route: str
    dose_mg: float
    infusion_duration_h: float


def simulate_iv(
    topology: PBPKTopology,
    params: CompoundPBPKParams,
    *,
    dose_mg: float,
    route: Literal["iv_bolus", "iv_infusion"],
    duration_h: float,
    infusion_duration_h: float = 0.0,
    n_time_points: int = 500,
    rtol: float = 1e-6,
    atol: float | None = None,
    method: str = "BDF",
) -> SimulationResult:
    """Simulate an IV dose through the human PBPK kernel.

    Parameters
    ----------
    topology : PBPKTopology
        Loaded species topology.
    params : CompoundPBPKParams
        Compound-specific PBPK parameters (Kp, fu_b, CLint_liver_L_h, etc.).
    dose_mg : float
        Dose in mg.
    route : "iv_bolus" or "iv_infusion"
        "iv_bolus" places all dose into the venous pool at t=0.
        "iv_infusion" injects at a constant rate over infusion_duration_h.
    duration_h : float
        Simulation end time (hours).
    infusion_duration_h : float, optional
        Required for "iv_infusion" (must be > 0).
    n_time_points : int
        Number of evaluation points (uniform grid on [0, duration_h]).
    rtol, atol : float
        Solver tolerances.  ``atol`` defaults to ``1e-9 * dose_mg``.
    method : str
        Must be "BDF" — explicit methods are rejected (PBPK is stiff).
    """
    if method not in _ALLOWED_METHODS:
        raise ValueError(
            f"PBPK is stiff — only BDF is allowed, got {method!r}. "
            f"Explicit methods (RK45, DOP853, LSODA) silently produce wrong "
            f"results for PBPK systems."
        )

    if route not in ("iv_bolus", "iv_infusion"):
        raise ValueError(
            f"route must be 'iv_bolus' or 'iv_infusion', got {route!r}"
        )

    if dose_mg <= 0:
        raise ValueError(f"dose_mg must be > 0, got {dose_mg}")

    if duration_h <= 0:
        raise ValueError(f"duration_h must be > 0, got {duration_h}")

    if route == "iv_infusion" and infusion_duration_h <= 0:
        raise ValueError(
            f"iv_infusion requires infusion_duration_h > 0, got {infusion_duration_h}"
        )

    if atol is None:
        atol = 1e-9 * dose_mg

    n_tissues = len(topology.tissues)
    y0 = np.zeros(2 + n_tissues)

    if route == "iv_bolus":
        y0[0] = dose_mg
        infusion_rate = 0.0
        inf_duration = 0.0
    else:  # iv_infusion
        infusion_rate = dose_mg / infusion_duration_h
        inf_duration = infusion_duration_h

    rhs = build_rhs(
        topology,
        params,
        infusion_rate_mg_per_h=infusion_rate,
        infusion_duration_h=inf_duration,
    )

    t_eval = np.linspace(0.0, duration_h, n_time_points)

    sol = solve_ivp(
        rhs,
        (0.0, duration_h),
        y0,
        method=method,
        t_eval=t_eval,
        rtol=rtol,
        atol=atol,
        dense_output=False,
    )

    if not sol.success:
        raise RuntimeError(
            f"PBPK ODE solver failed for compound {params.name!r}: {sol.message}"
        )

    v_ven = topology.venous_volume_L
    cp_blood = sol.y[0] / v_ven
    cp_plasma = cp_blood / params.bp_ratio

    # mass balance residual at t=0 (sanity check on IC injection)
    residual = float(np.abs(sol.y.sum(axis=0)[0] - dose_mg))

    return SimulationResult(
        time_h=sol.t,
        cp_blood=cp_blood,
        cp_plasma=cp_plasma,
        state_trajectory=sol.y,
        mass_balance_residual=residual,
        solver_success=sol.success,
        solver_method=method,
        solver_nfev=int(sol.nfev),
        route=route,
        dose_mg=dose_mg,
        infusion_duration_h=inf_duration,
    )


@dataclass
class OralSimulationResult:
    """Output of an oral PBPK+ACAT simulation."""

    time_h: np.ndarray
    cp_blood: np.ndarray
    cp_plasma: np.ndarray
    state_trajectory: np.ndarray
    lumen_trajectory: np.ndarray
    enterocyte_trajectory: np.ndarray
    mass_balance_residual: float
    solver_success: bool
    solver_method: str
    solver_nfev: int
    route: str
    dose_mg: float


def simulate_oral(
    topology: PBPKTopology,
    params: OralPBPKParams,
    *,
    dose_mg: float,
    duration_h: float,
    n_time_points: int = 500,
    rtol: float = 1e-6,
    atol: float | None = None,
    method: str = "BDF",
) -> OralSimulationResult:
    """Simulate an oral dose through the PBPK+ACAT kernel.

    Parameters
    ----------
    topology : PBPKTopology
        Loaded species topology.
    params : OralPBPKParams
        Compound-specific oral PBPK parameters.
    dose_mg : float
        Oral dose in mg (placed in stomach at t=0).
    duration_h : float
        Simulation end time (hours).
    n_time_points : int
        Number of evaluation points (uniform grid on [0, duration_h]).
    rtol, atol : float
        Solver tolerances.  ``atol`` defaults to ``1e-9 * dose_mg``.
    method : str
        Must be "BDF" — explicit methods are rejected (PBPK is stiff).
    """
    if method not in _ALLOWED_METHODS:
        raise ValueError(
            f"PBPK is stiff — only BDF is allowed, got {method!r}. "
            f"Explicit methods (RK45, DOP853, LSODA) silently produce wrong "
            f"results for PBPK systems."
        )
    if dose_mg <= 0:
        raise ValueError(f"dose_mg must be > 0, got {dose_mg}")
    if duration_h <= 0:
        raise ValueError(f"duration_h must be > 0, got {duration_h}")

    if atol is None:
        atol = 1e-9 * dose_mg

    gi = params.gi_tract
    assert gi is not None, "OralPBPKParams.gi_tract must not be None for simulate_oral"

    n_tissues = len(topology.tissues)
    n_lumen = len(gi.segments)
    n_states = 2 + n_tissues + n_lumen + 1  # 26

    y0 = np.zeros(n_states)
    # Oral dose: all in stomach lumen at t=0
    stomach_idx = 2 + n_tissues  # first lumen state = 17 (for 15 tissues)
    y0[stomach_idx] = dose_mg

    rhs = build_oral_rhs(topology, params)

    t_eval = np.linspace(0.0, duration_h, n_time_points)

    sol = solve_ivp(
        rhs,
        (0.0, duration_h),
        y0,
        method=method,
        t_eval=t_eval,
        rtol=rtol,
        atol=atol,
        dense_output=False,
    )

    if not sol.success:
        raise RuntimeError(
            f"Oral ODE solver failed for {params.name!r}: {sol.message}"
        )

    v_ven = topology.venous_volume_L
    cp_blood = sol.y[0] / v_ven
    cp_plasma = cp_blood / params.bp_ratio

    # Mass balance at t=0
    residual = float(np.abs(sol.y.sum(axis=0)[0] - dose_mg))

    # Extract lumen and enterocyte trajectories
    lumen_start = 2 + n_tissues
    entero_idx = lumen_start + n_lumen
    lumen_traj = sol.y[lumen_start:entero_idx, :]
    entero_traj = sol.y[entero_idx, :]

    return OralSimulationResult(
        time_h=sol.t,
        cp_blood=cp_blood,
        cp_plasma=cp_plasma,
        state_trajectory=sol.y,
        lumen_trajectory=lumen_traj,
        enterocyte_trajectory=entero_traj,
        mass_balance_residual=residual,
        solver_success=sol.success,
        solver_method=method,
        solver_nfev=int(sol.nfev),
        route="oral",
        dose_mg=dose_mg,
    )
