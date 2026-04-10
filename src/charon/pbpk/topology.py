"""PBPK topology: load a species YAML into an immutable, validated graph.

The PBPK topology describes which tissues exist, their volumes, blood flows,
tissue compositions (for Rodgers & Rowland Kp), and how they are wired into
the vascular tree (arterial input, venous or portal drainage, hepatic
arterial split).

The topology is a pure data container — no ODE logic lives here.  The
``ode_compiler`` module consumes it to build the right-hand side.
"""

from __future__ import annotations

from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import yaml

from charon.pbpk.kp_calculator import TissueComposition

# Tissues that drain into the portal vein rather than the systemic venous
# pool.  These provide the "portal inflow" term for the liver compartment.
PORTAL_TISSUES: tuple[str, ...] = ("spleen", "gut_wall", "pancreas")

_DRAINAGE_OVERRIDES: dict[str, str] = {
    "lung": "arterial",
    "liver": "venous",
}

_SPECIES_DIR = Path(__file__).parent / "species"


@dataclass(frozen=True)
class TissueNode:
    """Single tissue compartment in the PBPK graph."""

    name: str
    volume_L: float
    blood_flow_L_h: float
    composition: TissueComposition
    drains_to: Literal["venous", "liver", "arterial"]


@dataclass(frozen=True)
class PBPKTopology:
    """Immutable species-level PBPK graph."""

    species: str
    body_weight_kg: float
    cardiac_output_L_h: float
    hematocrit: float
    venous_volume_L: float
    arterial_volume_L: float
    hepatic_artery_L_h: float
    tissues: "OrderedDict[str, TissueNode]"
    plasma_composition: TissueComposition
    gfr_mL_min: float
    liver_weight_g: float

    def tissue_names(self) -> tuple[str, ...]:
        """Deterministic tissue ordering (YAML insertion order)."""
        return tuple(self.tissues.keys())


def _tissue_comp_from_dict(comp: dict) -> TissueComposition:
    return TissueComposition(
        fn=float(comp["fn"]),
        fp=float(comp["fp"]),
        fw=float(comp["fw"]),
        pH=float(comp["pH"]),
    )


def _drainage_for(name: str) -> Literal["venous", "liver", "arterial"]:
    if name in _DRAINAGE_OVERRIDES:
        return _DRAINAGE_OVERRIDES[name]  # type: ignore[return-value]
    if name in PORTAL_TISSUES:
        return "liver"
    return "venous"


def load_species_topology(
    species: str,
    *,
    path: Path | str | None = None,
) -> PBPKTopology:
    """Load a PBPK topology from ``src/charon/pbpk/species/<name>.yaml``.

    Parameters
    ----------
    species : str
        Species name used as the YAML filename stem.
    path : Path or str, optional
        Explicit path to the YAML file.  If supplied, ``species`` is used
        only as the returned ``PBPKTopology.species`` identifier.

    Returns
    -------
    PBPKTopology
        Frozen topology container.

    Raises
    ------
    FileNotFoundError
        If the YAML file does not exist.
    ValueError
        If the topology fails internal consistency checks (portal sum
        exceeds total liver flow, missing mandatory tissues, etc.).
    """
    if path is None:
        yaml_path = _SPECIES_DIR / f"{species}.yaml"
    else:
        yaml_path = Path(path)

    if not yaml_path.exists():
        raise FileNotFoundError(
            f"Species YAML not found: {yaml_path}"
        )

    with yaml_path.open() as fp:
        data = yaml.safe_load(fp)

    spec = data["species"]
    co_L_h = float(spec["cardiac_output_L_h"])

    tissues: OrderedDict[str, TissueNode] = OrderedDict()
    for tissue_name, tissue_spec in spec["tissues"].items():
        flow_fraction = float(tissue_spec["blood_flow_fraction"])
        if tissue_name == "lung":
            q = co_L_h  # pulmonary: full CO
        else:
            q = flow_fraction * co_L_h
        node = TissueNode(
            name=tissue_name,
            volume_L=float(tissue_spec["volume_L"]),
            blood_flow_L_h=q,
            composition=_tissue_comp_from_dict(tissue_spec["composition"]),
            drains_to=_drainage_for(tissue_name),
        )
        tissues[tissue_name] = node

    # Mandatory tissues for the minimal IV-only kernel.
    for required in ("lung", "liver", "kidney", *PORTAL_TISSUES):
        if required not in tissues:
            raise ValueError(
                f"species {species!r} topology is missing required tissue {required!r}"
            )

    # Hepatic artery = Q_liver_total - Q_portal_total
    q_liver_total = tissues["liver"].blood_flow_L_h
    q_portal_total = sum(tissues[p].blood_flow_L_h for p in PORTAL_TISSUES)
    q_hepatic_artery = q_liver_total - q_portal_total
    if q_hepatic_artery < 0:
        raise ValueError(
            f"species {species!r}: sum of portal tissue flows "
            f"({q_portal_total:.3g} L/h) exceeds total liver flow "
            f"({q_liver_total:.3g} L/h); topology is inconsistent"
        )

    plasma_comp = _tissue_comp_from_dict(spec["plasma"])

    return PBPKTopology(
        species=species,
        body_weight_kg=float(spec["body_weight_kg"]),
        cardiac_output_L_h=co_L_h,
        hematocrit=float(spec["hematocrit"]),
        venous_volume_L=float(spec["blood"]["venous_volume_L"]),
        arterial_volume_L=float(spec["blood"]["arterial_volume_L"]),
        hepatic_artery_L_h=q_hepatic_artery,
        tissues=tissues,
        plasma_composition=plasma_comp,
        gfr_mL_min=float(spec["kidney"]["gfr_mL_min"]),
        liver_weight_g=float(spec["liver"]["weight_g"]),
    )
