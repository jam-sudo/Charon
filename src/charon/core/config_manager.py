"""RunConfig management: create, freeze, serialise, diff, and hash.

Provides helper utilities for working with :class:`RunConfig` instances
throughout the Charon pipeline lifecycle.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import yaml

from charon.core.schema import CompoundConfig, PipelineConfig, RunConfig


# ---------------------------------------------------------------------------
# Creation
# ---------------------------------------------------------------------------


def create_run_config(
    compound: CompoundConfig,
    pipeline: PipelineConfig | None = None,
    model_versions: dict[str, str] | None = None,
) -> RunConfig:
    """Create a :class:`RunConfig` from its constituent parts.

    Applies sensible defaults for any missing component:

    * *pipeline* defaults to a fresh :class:`PipelineConfig` (all defaults).
    * *model_versions* defaults to an empty dict.
    """
    return RunConfig(
        compound=compound,
        pipeline=pipeline if pipeline is not None else PipelineConfig(),
        model_versions=model_versions if model_versions is not None else {},
    )


# ---------------------------------------------------------------------------
# YAML serialisation
# ---------------------------------------------------------------------------


def config_to_yaml(config: RunConfig, path: str | Path) -> None:
    """Serialise a :class:`RunConfig` to a YAML file.

    The resulting YAML is self-contained and sufficient to reproduce the run.
    ``None`` values are excluded to keep the file concise.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    data = config.model_dump(exclude_none=True)
    with open(path, "w") as fh:
        yaml.dump(data, fh, default_flow_style=False, sort_keys=False)


def config_from_yaml(path: str | Path) -> RunConfig:
    """Deserialise a :class:`RunConfig` from a YAML file.

    Validates against the Pydantic schema on load.

    Raises
    ------
    FileNotFoundError
        If *path* does not exist.
    pydantic.ValidationError
        On schema violations.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"RunConfig file not found: {path}")
    with open(path, "r") as fh:
        data = yaml.safe_load(fh)
    return RunConfig(**data)


# ---------------------------------------------------------------------------
# Diff
# ---------------------------------------------------------------------------


def _recursive_diff(
    a: Any,
    b: Any,
    path: str = "",
) -> dict[str, Any]:
    """Walk two structures and collect leaf-level differences.

    Returns a flat dict keyed by dotted path with ``{"old": ..., "new": ...}``
    values for every field that differs between *a* and *b*.
    """
    diffs: dict[str, Any] = {}

    if isinstance(a, dict) and isinstance(b, dict):
        all_keys = set(a) | set(b)
        for key in sorted(all_keys):
            child_path = f"{path}.{key}" if path else key
            if key not in a:
                diffs[child_path] = {"old": None, "new": b[key]}
            elif key not in b:
                diffs[child_path] = {"old": a[key], "new": None}
            else:
                diffs.update(_recursive_diff(a[key], b[key], child_path))
    elif isinstance(a, list) and isinstance(b, list):
        max_len = max(len(a), len(b))
        for i in range(max_len):
            child_path = f"{path}[{i}]"
            if i >= len(a):
                diffs[child_path] = {"old": None, "new": b[i]}
            elif i >= len(b):
                diffs[child_path] = {"old": a[i], "new": None}
            else:
                diffs.update(_recursive_diff(a[i], b[i], child_path))
    else:
        if a != b:
            diffs[path] = {"old": a, "new": b}

    return diffs


def diff_configs(config_a: RunConfig, config_b: RunConfig) -> dict[str, Any]:
    """Compute a recursive diff between two :class:`RunConfig` instances.

    Returns a dict of changed fields keyed by dotted path, each containing
    ``{"old": ..., "new": ...}`` values.  An empty dict means the two configs
    are identical.
    """
    dump_a = config_a.model_dump()
    dump_b = config_b.model_dump()
    return _recursive_diff(dump_a, dump_b)


# ---------------------------------------------------------------------------
# Hashing
# ---------------------------------------------------------------------------


def hash_config(config: RunConfig) -> str:
    """Compute a deterministic SHA-256 hash of a :class:`RunConfig`.

    The config is serialised to canonical JSON (sorted keys, no extra
    whitespace) before hashing.  The same config will always produce the same
    hash.

    .. note::

       This hashes the *config*, not the simulation output.  ODE solver
       results are **not** bitwise reproducible across platforms.
       "Reproducible" means within numerical tolerance (``rtol=1e-6``), not
       bit-identical.
    """
    canonical = json.dumps(
        config.model_dump(),
        sort_keys=True,
        separators=(",", ":"),
        default=str,
    )
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()
