"""YAML compound-configuration I/O and experimental-override helpers.

Depends on :mod:`charon.core.schema` (Pydantic v2 models) and PyYAML.
"""

from __future__ import annotations

from pathlib import Path

import yaml

from charon.core.schema import CompoundConfig, PredictedProperty


# -- Override field -> dotted path mapping -----------------------------------

_OVERRIDE_MAP: dict[str, tuple[str, str]] = {
    "fu_p":              ("binding",        "fu_p"),
    "fu_inc":            ("binding",        "fu_inc"),
    "bp_ratio":          ("binding",        "bp_ratio"),
    "clint_uL_min_mg":   ("metabolism",     "clint_uL_min_mg"),
    "papp_nm_s":         ("permeability",   "papp_nm_s"),
    "solubility_ug_ml":  ("physicochemical", "solubility_ug_ml"),
    "logp":              ("physicochemical", "logp"),
    "herg_ic50_uM":      ("safety",         "herg_ic50_uM"),
    "clrenal_L_h":       ("renal",          "clrenal_L_h"),
}


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def load_compound_config(path: str | Path) -> CompoundConfig:
    """Load compound configuration from a YAML file.

    Parses into a Pydantic :class:`CompoundConfig` model (validation runs
    automatically on construction).

    Raises
    ------
    FileNotFoundError
        If *path* does not exist.
    pydantic.ValidationError
        On schema violations.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Compound config file not found: {path}")
    with open(path, "r") as fh:
        data = yaml.safe_load(fh)
    return CompoundConfig(**data)


def save_compound_config(config: CompoundConfig, path: str | Path) -> None:
    """Save compound configuration to a YAML file.

    Serialises via ``config.model_dump(exclude_none=True)`` for clean YAML and
    uses ``default_flow_style=False`` for human readability.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    data = config.model_dump(exclude_none=True)
    with open(path, "w") as fh:
        yaml.dump(data, fh, default_flow_style=False, sort_keys=False)


def apply_overrides(
    config: CompoundConfig,
    overrides: dict[str, dict],
) -> CompoundConfig:
    """Apply experimental overrides to a compound config.

    Parameters
    ----------
    config:
        The base :class:`CompoundConfig` to override.
    overrides:
        Mapping of property names to dicts containing at minimum ``"value"``
        and optionally ``"source"`` (defaults to ``"experimental"``),
        ``"method"``, ``"unit"``, ``"ci_90_lower"``, ``"ci_90_upper"``,
        and ``"flag"``.

    Returns
    -------
    CompoundConfig
        A **new** ``CompoundConfig`` with the overrides applied.  The original
        *config* is not modified.

    Raises
    ------
    KeyError
        If an override key is not recognised.
    """
    # Deep-copy via round-trip through model_dump so we never mutate the input.
    data = config.model_dump()

    for key, override_values in overrides.items():
        if key not in _OVERRIDE_MAP:
            raise KeyError(
                f"Unknown override key {key!r}. "
                f"Valid keys: {sorted(_OVERRIDE_MAP)}"
            )

        group, field = _OVERRIDE_MAP[key]

        # Build a PredictedProperty from the override values.
        prop_kwargs: dict = {"value": override_values["value"]}
        prop_kwargs["source"] = override_values.get("source", "experimental")
        for optional in ("method", "unit", "ci_90_lower", "ci_90_upper", "flag"):
            if optional in override_values:
                prop_kwargs[optional] = override_values[optional]

        prop = PredictedProperty(**prop_kwargs)

        # Inject into the nested dict.
        data["properties"][group][field] = prop.model_dump()

    return CompoundConfig(**data)
