"""Charon Layer 2 — physiologically-based pharmacokinetic simulation.

Public API:

    PBPKTopology, TissueNode, load_species_topology, PORTAL_TISSUES
        Species-level PBPK graph loaded from a YAML file.

    CompoundPBPKParams, build_compound_pbpk_params, build_rhs,
    infer_compound_type
        Translate a :class:`CompoundConfig` into PBPK-ready parameters
        and build the ODE right-hand side closure.

    SimulationResult, simulate_iv
        Run an IV bolus or infusion through the human PBPK kernel using
        scipy's BDF stiff ODE solver.

    compute_pk_parameters
        Extract Cmax, AUC, CL, Vss, t_half from a concentration-time
        profile.

    compute_kp_rodgers_rowland, compute_all_kp, TissueComposition,
    apply_berezhkovskiy_correction
        Mechanistic Kp calculators (Rodgers & Rowland 2005/2006 plus
        Berezhkovskiy 2004 correction).
"""

from __future__ import annotations

from charon.pbpk.acat import (
    GISegment,
    GITract,
    compute_absorption_rates,
    load_gi_tract,
    papp_to_peff,
)
from charon.pbpk.kp_calculator import (
    TissueComposition,
    apply_berezhkovskiy_correction,
    compute_all_kp,
    compute_kp_poulin_theil,
    compute_kp_rodgers_rowland,
)
from charon.pbpk.ode_compiler import (
    CompoundPBPKParams,
    OralPBPKParams,
    build_compound_pbpk_params,
    build_oral_rhs,
    build_rhs,
    compute_gut_clint,
    infer_compound_type,
)
from charon.pbpk.pk_extract import compute_oral_pk_parameters, compute_pk_parameters
from charon.pbpk.solver import OralSimulationResult, SimulationResult, simulate_iv, simulate_oral
from charon.pbpk.topology import (
    PORTAL_TISSUES,
    PBPKTopology,
    TissueNode,
    load_species_topology,
)

__all__ = [
    "CompoundPBPKParams",
    "GISegment",
    "GITract",
    "OralPBPKParams",
    "OralSimulationResult",
    "PBPKTopology",
    "PORTAL_TISSUES",
    "SimulationResult",
    "TissueComposition",
    "TissueNode",
    "apply_berezhkovskiy_correction",
    "build_compound_pbpk_params",
    "build_oral_rhs",
    "build_rhs",
    "compute_absorption_rates",
    "compute_all_kp",
    "compute_gut_clint",
    "compute_kp_poulin_theil",
    "compute_kp_rodgers_rowland",
    "compute_oral_pk_parameters",
    "compute_pk_parameters",
    "infer_compound_type",
    "load_gi_tract",
    "load_species_topology",
    "papp_to_peff",
    "simulate_iv",
    "simulate_oral",
]
