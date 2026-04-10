"""Pydantic v2 data models for the Charon translational PK pipeline.

All models use Pydantic v2 patterns exclusively (``ConfigDict``,
``field_validator``, ``model_validator``).  No v1-style ``class Config``
or ``@validator`` decorators are used anywhere.
"""

from __future__ import annotations

import math
from typing import Literal

from pydantic import BaseModel, ConfigDict, field_validator, model_validator

# ---------------------------------------------------------------------------
# Type aliases
# ---------------------------------------------------------------------------

SourceType = Literal[
    "ml_ensemble",
    "ml_pka",
    "correlation",
    "derived",
    "physiological",
    "experimental",
]
"""Taxonomy of data-source provenance."""

# ---------------------------------------------------------------------------
# Core value wrapper
# ---------------------------------------------------------------------------


class PredictedProperty(BaseModel):
    """A predicted (or measured) scalar with optional 90% confidence interval."""

    value: float
    ci_90_lower: float | None = None
    ci_90_upper: float | None = None
    source: SourceType
    unit: str | None = None
    method: str | None = None
    flag: str | None = None

    @field_validator("value")
    @classmethod
    def _value_must_be_finite(cls, v: float) -> float:
        if not math.isfinite(v):
            raise ValueError(f"value must be finite, got {v}")
        return v

    @model_validator(mode="after")
    def _ci_bounds_consistent(self) -> "PredictedProperty":
        if self.ci_90_lower is not None and self.ci_90_upper is not None:
            if self.ci_90_lower > self.value:
                raise ValueError(
                    f"ci_90_lower ({self.ci_90_lower}) must be "
                    f"<= value ({self.value})"
                )
            if self.value > self.ci_90_upper:
                raise ValueError(
                    f"value ({self.value}) must be "
                    f"<= ci_90_upper ({self.ci_90_upper})"
                )
        return self


# ---------------------------------------------------------------------------
# Salt / formulation helpers
# ---------------------------------------------------------------------------


class SaltForm(BaseModel):
    """Describes a salt form and its relationship to the free-base parent."""

    name: str | None = None
    mw_salt: float | None = None
    salt_factor: float = 1.0

    @model_validator(mode="after")
    def _auto_salt_factor(self) -> "SaltForm":
        """If ``mw_salt`` is provided the ``salt_factor`` stays as given.

        Automatic calculation of ``salt_factor`` from ``mw_salt`` and the
        parent ``molecular_weight`` is deferred to the ``CompoundConfig``
        model validator, which has access to both values.
        """
        return self


# ---------------------------------------------------------------------------
# Conversion audit trail
# ---------------------------------------------------------------------------


class ConversionStep(BaseModel):
    """Single intermediate step in a unit-conversion pipeline."""

    name: str
    value: float
    unit: str
    formula: str


class ConversionLog(BaseModel):
    """Full audit log for a multi-step unit conversion."""

    input_params: dict[str, float | str]
    intermediate_steps: list[ConversionStep]
    output: float
    output_unit: str
    model_used: str


# ---------------------------------------------------------------------------
# Hepatic clearance result
# ---------------------------------------------------------------------------


class HepaticClearance(BaseModel):
    """Result of an in-vitro → in-vivo hepatic-clearance extrapolation."""

    clh_L_h: float
    extraction_ratio: float
    model_used: str
    conversion_log: ConversionLog
    clint_liver_L_h: float | None = None
    """Whole-liver intrinsic clearance after IVIVE scaling, in L/h.

    This is the value *before* the liver extraction model is applied.  It
    is the correct input for PBPK ODEs that embed well-stirred elimination
    directly (rate_elim = CLint_liver * fu_b * C_liver_blood_out).  Using
    ``clh_L_h`` instead in a PBPK rhs would double-apply the liver model.
    """


# ---------------------------------------------------------------------------
# Guardrails & validation
# ---------------------------------------------------------------------------


class GuardrailWarning(BaseModel):
    """A single guardrail or applicability-domain warning."""

    category: str
    message: str
    severity: Literal["INFO", "WARNING", "CRITICAL"]


class ValidationResult(BaseModel):
    """Output of the compound-validation / guardrail layer."""

    is_valid: bool
    smiles_canonical: str
    molecular_weight: float
    warnings: list[GuardrailWarning] = []
    applicability_domain: Literal["HIGH", "MODERATE", "LOW", "UNKNOWN"] = "UNKNOWN"
    descriptors: dict[str, float] = {}


# ---------------------------------------------------------------------------
# Property groups
# ---------------------------------------------------------------------------


class PhysicochemicalProperties(BaseModel):
    """Lipophilicity, pKa, solubility."""

    logp: PredictedProperty | None = None
    pka_acid: PredictedProperty | None = None
    pka_base: PredictedProperty | None = None
    solubility_ug_ml: PredictedProperty | None = None
    compound_type: Literal[
        "neutral", "acid", "base", "zwitterion"
    ] | None = None


class PermeabilityProperties(BaseModel):
    """Apparent permeability and transporter substrate predictions."""

    papp_nm_s: PredictedProperty | None = None
    peff_cm_s: PredictedProperty | None = None
    pgp_substrate: PredictedProperty | None = None
    oatp_substrate: PredictedProperty | None = None


class BindingProperties(BaseModel):
    """Plasma/incubational binding and blood-plasma partitioning."""

    fu_p: PredictedProperty | None = None
    fu_inc: PredictedProperty | None = None
    bp_ratio: PredictedProperty | None = None


class MetabolismProperties(BaseModel):
    """CYP phenotyping and intrinsic clearance."""

    primary_cyp: str | None = None
    secondary_cyp: str | None = None
    clint_uL_min_mg: PredictedProperty | None = None


class SafetyProperties(BaseModel):
    """hERG liability and CYP inhibition panel."""

    herg_ic50_uM: PredictedProperty | None = None
    cyp_inhibition: dict[str, PredictedProperty] = {}


class RenalProperties(BaseModel):
    """Renal clearance parameters."""

    clrenal_L_h: PredictedProperty | None = None
    active_secretion: bool = False


class DistributionProperties(BaseModel):
    """Distribution-related compound properties.

    Currently holds only the empirical Kp override. Forward-compatible
    for Vss_pred, tissue fu_tissue, binding capacity etc.
    """

    empirical_kp_by_tissue: dict[str, PredictedProperty] | None = None

    @field_validator("empirical_kp_by_tissue")
    @classmethod
    def _validate_kp_values(cls, v):
        if v is None:
            return v
        if not v:
            raise ValueError(
                "empirical_kp_by_tissue must be None or a non-empty dict"
            )
        for tissue, p in v.items():
            if p.value <= 0 or p.value > 200:
                raise ValueError(
                    f"empirical_kp_by_tissue[{tissue!r}] = {p.value} "
                    f"outside physiological range (0, 200]"
                )
        return v


# ---------------------------------------------------------------------------
# Compound-level aggregates
# ---------------------------------------------------------------------------


class CompoundProperties(BaseModel):
    """Container for all predicted/measured compound properties."""

    physicochemical: PhysicochemicalProperties = PhysicochemicalProperties()
    permeability: PermeabilityProperties = PermeabilityProperties()
    binding: BindingProperties = BindingProperties()
    metabolism: MetabolismProperties = MetabolismProperties()
    safety: SafetyProperties = SafetyProperties()
    renal: RenalProperties = RenalProperties()


class CompoundConfig(BaseModel):
    """Top-level compound specification (identity + properties)."""

    name: str
    smiles: str
    molecular_weight: float | None = None
    salt_form: SaltForm | None = None
    source: Literal["predicted", "experimental", "mixed"] = "predicted"
    properties: CompoundProperties = CompoundProperties()

    @model_validator(mode="after")
    def _auto_salt_factor(self) -> "CompoundConfig":
        """Compute ``salt_factor`` when both ``mw_salt`` and ``molecular_weight``
        are available and ``salt_factor`` has not been explicitly set
        (i.e. still at the default of 1.0).
        """
        if (
            self.salt_form is not None
            and self.salt_form.mw_salt is not None
            and self.molecular_weight is not None
            and self.salt_form.salt_factor == 1.0
        ):
            self.salt_form.salt_factor = (
                self.molecular_weight / self.salt_form.mw_salt
            )
        return self


# ---------------------------------------------------------------------------
# Pipeline / run configuration
# ---------------------------------------------------------------------------


class FormulationConfig(BaseModel):
    """Oral-formulation parameters."""

    particle_size_um: float = 10.0
    dose_mg: float = 100.0
    route: Literal["oral", "iv_bolus", "iv_infusion"] = "oral"


class DoseProjectionConfig(BaseModel):
    """First-in-human dose-projection inputs."""

    noael_mg_kg: float | None = None
    noael_species: str | None = None
    safety_factor: float = 10.0
    target_kd_nM: float | None = None
    target_ceff_nM: float | None = None


class UncertaintyConfig(BaseModel):
    """Uncertainty/sensitivity-analysis settings."""

    method: Literal["lhs", "monte_carlo"] = "lhs"
    n_samples: int = 500
    correlation: Literal["iman_conover", "none"] = "iman_conover"


class PipelineConfig(BaseModel):
    """Controls which pipeline layers and species to run."""

    layers: list[int] = [0, 1, 2, 3, 4]
    species: list[str] = ["human"]
    preclinical_species: list[str] = []
    liver_model: Literal["well_stirred", "parallel_tube", "dispersion"] = (
        "well_stirred"
    )
    formulation: FormulationConfig = FormulationConfig()
    dose_projection: DoseProjectionConfig = DoseProjectionConfig()
    uncertainty: UncertaintyConfig = UncertaintyConfig()


# ---------------------------------------------------------------------------
# PK output parameters
# ---------------------------------------------------------------------------


class PKParameters(BaseModel):
    """Standard pharmacokinetic output parameters."""

    cmax: float | None = None
    tmax: float | None = None
    auc_0_inf: float | None = None
    auc_0_24: float | None = None
    half_life: float | None = None
    cl_apparent: float | None = None
    vss: float | None = None
    bioavailability: float | None = None
    fa: float | None = None
    fg: float | None = None
    fh: float | None = None


# ---------------------------------------------------------------------------
# Top-level run configuration (immutable)
# ---------------------------------------------------------------------------


class RunConfig(BaseModel):
    """Immutable top-level configuration for a single pipeline run."""

    model_config = ConfigDict(frozen=True)

    compound: CompoundConfig
    pipeline: PipelineConfig = PipelineConfig()
    model_versions: dict[str, str] = {}
