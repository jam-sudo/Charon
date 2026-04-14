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
    DoseProjectionConfig,
    PKParameters,
    UncertaintyConfig,
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
    dose_recommendation: "FIHDoseRecommendation | None" = None
    uncertainty: "UncertaintyResult | None" = None


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
        dose_projection: DoseProjectionConfig | None = None,
        uncertainty: UncertaintyConfig | None = None,
    ) -> None:
        self.compound = compound
        self.route = route
        self.dose_mg = dose_mg
        self.species = species
        self.duration_h = duration_h
        self.infusion_duration_h = infusion_duration_h
        self.liver_model = liver_model
        self.compound_type_override = compound_type_override
        self.dose_projection = dose_projection
        self.uncertainty = uncertainty

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

    def _run_uncertainty(self):
        """Run uncertainty propagation if configured."""
        if self.uncertainty is None:
            return None
        if self.dose_projection is None:
            raise ValueError(
                "uncertainty requires dose_projection to be set "
                "(need dose projection to quantify dose CI)"
            )
        from charon.uncertainty.sampling import build_param_specs, generate_lhs_samples
        from charon.uncertainty.propagation import propagate
        from charon.uncertainty.dose_range import compute_dose_range

        param_specs = build_param_specs(self.compound)
        sampling_result = generate_lhs_samples(
            param_specs=param_specs,
            n_samples=self.uncertainty.n_samples,
            correlation=self.uncertainty.correlation,
        )
        prop_result = propagate(
            base_compound=self.compound,
            samples=sampling_result.samples,
            route=self.route,
            dose_mg=self.dose_mg,
            dose_projection=self.dose_projection,
            duration_h=self.duration_h,
            liver_model=self.liver_model,
        )
        if prop_result.n_successful < 5:
            return None
        return compute_dose_range(
            prop_result.doses_mg,
            sensitivity={},
            param_names=prop_result.param_names,
            parameter_matrix=prop_result.parameter_matrix,
        )

    def _maybe_project_dose(self, pk: PKParameters) -> "FIHDoseRecommendation | None":
        """Run FIH dose projection if config is provided and has sufficient inputs."""
        if self.dose_projection is None:
            return None
        dp = self.dose_projection
        has_hed = dp.noael_mg_kg is not None and dp.noael_species is not None
        has_mabel = dp.target_kd_nM is not None
        has_pad = dp.target_ceff_nM is not None
        if not (has_hed or has_mabel or has_pad):
            return None
        from charon.translational.dose_projector import project_fih_dose
        return project_fih_dose(
            pk=pk, compound=self.compound, config=self.dose_projection, route=self.route,
        )

    def run(self) -> PipelineResult:
        """Execute the full pipeline and return a :class:`PipelineResult`."""
        if self.route == "oral":
            return self._run_oral()

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

        dose_rec = self._maybe_project_dose(pk)
        unc = self._run_uncertainty()

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
                "kp_overrides": [
                    {
                        "tissue": r.tissue,
                        "rr_value": r.rr_value,
                        "empirical_value": r.empirical_value,
                        "source": r.source,
                        "method": r.method,
                        "flag": r.flag,
                    }
                    for r in params.kp_overrides
                ],
                "mrsd_mg": dose_rec.mrsd_mg if dose_rec else None,
                "uncertainty_ci_90": f"{unc.ci_90_lower_mg:.1f}-{unc.ci_90_upper_mg:.1f}" if unc else None,
            },
            dose_recommendation=dose_rec,
            uncertainty=unc,
        )

    def _run_oral(self) -> PipelineResult:
        """Execute oral route through ACAT + PBPK."""
        from charon.pbpk.acat import load_gi_tract, papp_to_peff
        from charon.pbpk.ode_compiler import (
            OralPBPKParams,
            compute_gut_clint,
        )
        from charon.pbpk.pk_extract import compute_oral_pk_parameters
        from charon.pbpk.solver import simulate_oral

        topology = load_species_topology(self.species)
        bridge = ParameterBridge()

        base_params = build_compound_pbpk_params(
            self.compound,
            topology,
            bridge=bridge,
            compound_type=self.compound_type_override,
            liver_model=self.liver_model,
        )

        gi = load_gi_tract(self.species)

        # Resolve Peff
        perm = self.compound.properties.permeability
        if perm.peff_cm_s is not None:
            peff = float(perm.peff_cm_s.value)
        elif perm.papp_nm_s is not None:
            peff = papp_to_peff(float(perm.papp_nm_s.value))
        else:
            raise ValueError(
                f"Oral route requires Peff or Papp for compound "
                f"{self.compound.name!r}, but neither is provided."
            )

        # Gut CLint
        fm = self.compound.properties.metabolism.fm_cyp3a4
        clint_gut = compute_gut_clint(
            clint_liver_L_h=base_params.clint_liver_L_h,
            fm_cyp3a4=fm,
            gi_tract=gi,
            mppgl=40.0,
            liver_weight_g=topology.liver_weight_g,
        )

        q_gut = topology.tissues["gut_wall"].blood_flow_L_h
        q_villi = gi.q_villi_fraction * q_gut

        oral_params = OralPBPKParams(
            name=base_params.name,
            molecular_weight=base_params.molecular_weight,
            logp=base_params.logp,
            pka_acid=base_params.pka_acid,
            pka_base=base_params.pka_base,
            compound_type=base_params.compound_type,
            fu_p=base_params.fu_p,
            bp_ratio=base_params.bp_ratio,
            fu_b=base_params.fu_b,
            clint_liver_L_h=base_params.clint_liver_L_h,
            cl_renal_L_h=base_params.cl_renal_L_h,
            kp_by_tissue=base_params.kp_by_tissue,
            kp_overrides=base_params.kp_overrides,
            clint_gut_L_h=clint_gut,
            peff_cm_s=peff,
            q_villi_L_h=q_villi,
            v_enterocyte_L=gi.enterocyte_volume_L,
            gi_tract=gi,
        )

        sim = simulate_oral(
            topology,
            oral_params,
            dose_mg=self.dose_mg,
            duration_h=self.duration_h,
        )

        pk = compute_oral_pk_parameters(
            sim, oral_params, topology, dose_mg=self.dose_mg
        )

        dose_rec = self._maybe_project_dose(pk)
        unc = self._run_uncertainty()

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
                "compound_type": oral_params.compound_type,
                "clint_liver_L_h": oral_params.clint_liver_L_h,
                "clint_gut_L_h": oral_params.clint_gut_L_h,
                "cl_renal_L_h": oral_params.cl_renal_L_h,
                "fu_b": oral_params.fu_b,
                "peff_cm_s": oral_params.peff_cm_s,
                "solver_method": sim.solver_method,
                "solver_nfev": sim.solver_nfev,
                "fa": pk.fa,
                "fg": pk.fg,
                "fh": pk.fh,
                "bioavailability": pk.bioavailability,
                "mrsd_mg": dose_rec.mrsd_mg if dose_rec else None,
                "uncertainty_ci_90": f"{unc.ci_90_lower_mg:.1f}-{unc.ci_90_upper_mg:.1f}" if unc else None,
            },
            dose_recommendation=dose_rec,
            uncertainty=unc,
        )
