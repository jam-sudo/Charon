"""Top-level Charon pipeline: SMILES or CompoundConfig → PK prediction.

Sprint 3 scope: IV bolus / IV infusion only.  Oral (ACAT) routing is
reserved for Sprint 3b.

Two constructors are supported:

1. ``Pipeline(compound=cfg, ...)``
   — use a pre-populated CompoundConfig (experimental overrides).

2. ``Pipeline.from_smiles(smiles, ...)``
   — run Layer 1 ML prediction (ADMET ensemble + pKa + fu_inc + BP) to
     build the compound, then run PBPK.

Both paths converge on the same :meth:`run` method.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal

import numpy as np

from charon.core.parameter_bridge import ParameterBridge
from charon.core.schema import (
    CompoundConfig,
    CompoundProperties,
    PKParameters,
)
from charon.pbpk.ode_compiler import build_compound_pbpk_params
from charon.pbpk.pk_extract import compute_pk_parameters
from charon.pbpk.solver import SimulationResult, simulate_iv
from charon.pbpk.topology import PBPKTopology, load_species_topology


@dataclass
class PipelineResult:
    """Output of :meth:`Pipeline.run`."""

    compound: CompoundConfig
    pk_parameters: PKParameters
    time_h: np.ndarray
    cp_plasma: np.ndarray
    cp_blood: np.ndarray
    simulation: SimulationResult
    metadata: dict = field(default_factory=dict)


class Pipeline:
    """End-to-end Charon Phase A pipeline (Sprint 3 = IV only)."""

    def __init__(
        self,
        compound: CompoundConfig,
        *,
        route: Literal["iv_bolus", "iv_infusion", "oral"],
        dose_mg: float,
        species: str = "human",
        duration_h: float = 72.0,
        infusion_duration_h: float = 0.0,
        liver_model: str = "well_stirred",
        compound_type_override: str | None = None,
    ) -> None:
        self.compound = compound
        self.route = route
        self.dose_mg = dose_mg
        self.species = species
        self.duration_h = duration_h
        self.infusion_duration_h = infusion_duration_h
        self.liver_model = liver_model
        self.compound_type_override = compound_type_override

    @classmethod
    def from_smiles(
        cls,
        smiles: str,
        *,
        route: Literal["iv_bolus", "iv_infusion", "oral"],
        dose_mg: float,
        species: str = "human",
        duration_h: float = 72.0,
        infusion_duration_h: float = 0.0,
        liver_model: str = "well_stirred",
        compound_name: str | None = None,
    ) -> "Pipeline":
        """Build a Pipeline by running Layer 1 ML prediction on a SMILES."""
        from charon.core.molecule import Molecule
        from charon.predict import predict_properties

        mol = Molecule(smiles)
        descriptors = mol.descriptors()
        mw = float(descriptors["MW"])
        properties: CompoundProperties = predict_properties(smiles)
        compound = CompoundConfig(
            name=compound_name or smiles,
            smiles=smiles,
            molecular_weight=mw,
            source="predicted",
            properties=properties,
        )
        return cls(
            compound=compound,
            route=route,
            dose_mg=dose_mg,
            species=species,
            duration_h=duration_h,
            infusion_duration_h=infusion_duration_h,
            liver_model=liver_model,
        )

    def run(self) -> PipelineResult:
        """Execute the full pipeline and return a :class:`PipelineResult`."""
        if self.route == "oral":
            raise NotImplementedError(
                "Oral route (ACAT) is deferred to Sprint 3b. "
                "Use 'iv_bolus' or 'iv_infusion' for now."
            )

        topology: PBPKTopology = load_species_topology(self.species)
        bridge = ParameterBridge()

        params = build_compound_pbpk_params(
            self.compound,
            topology,
            bridge=bridge,
            compound_type=self.compound_type_override,
            liver_model=self.liver_model,
        )

        sim = simulate_iv(
            topology,
            params,
            dose_mg=self.dose_mg,
            route=self.route,  # type: ignore[arg-type]
            duration_h=self.duration_h,
            infusion_duration_h=self.infusion_duration_h,
        )

        pk = compute_pk_parameters(
            sim.time_h,
            sim.cp_plasma,
            dose_mg=self.dose_mg,
            route=self.route,  # type: ignore[arg-type]
            infusion_duration_h=self.infusion_duration_h,
        )

        return PipelineResult(
            compound=self.compound,
            pk_parameters=pk,
            time_h=sim.time_h,
            cp_plasma=sim.cp_plasma,
            cp_blood=sim.cp_blood,
            simulation=sim,
            metadata={
                "species": self.species,
                "route": self.route,
                "dose_mg": self.dose_mg,
                "duration_h": self.duration_h,
                "liver_model": self.liver_model,
                "compound_type": params.compound_type,
                "clint_liver_L_h": params.clint_liver_L_h,
                "cl_renal_L_h": params.cl_renal_L_h,
                "fu_b": params.fu_b,
                "solver_method": sim.solver_method,
                "solver_nfev": sim.solver_nfev,
            },
        )
