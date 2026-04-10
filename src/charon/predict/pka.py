"""Rule-based pKa prediction from SMILES.

Detects common ionizable functional groups via RDKit SMARTS patterns and
assigns literature-average pKa values. The intent is not to produce
precise pKa values (ChemAxon or ACD does that better) but to produce a
robust ``compound_type`` classification (acid / base / neutral /
zwitterion) and reasonable pKa_acid / pKa_base that feed the Rodgers &
Rowland Kp calculation.

SMARTS patterns are ported verbatim from
``Omega/src/omega_pbpk/prediction/pka_predictor.py`` (Omega v0.9).

compound_type classification mirrors Sisyphus ``chemistry._classify_from_pka``:
    - acidic pKa < 7.0  AND  basic pKa > 8.0  → "zwitterion"
    - acidic pKa < 7.0  only                  → "acid"
    - basic pKa  > 8.0  only                  → "base"
    - otherwise                               → "neutral"
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import ClassVar

from rdkit import Chem

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Functional group → baseline pKa (literature averages)
# ---------------------------------------------------------------------------
GROUP_PKA: dict[str, float] = {
    "carboxylic_acid": 4.5,
    "sulfonamide": 10.5,
    "imidazole": 6.8,
    "pyridine": 5.2,
    "phenol": 9.5,
    "amine_primary": 10.0,
    "amine_secondary": 9.5,
    "amine_tertiary": 9.5,
}

# Each entry: (group_name, SMARTS, "acid" | "base")
# Priority order: first match wins within each category.
# Patterns ported verbatim from Omega.
_ACID_SMARTS: list[tuple[str, str]] = [
    ("carboxylic_acid", "[CX3](=O)[OX2H1]"),
    ("sulfonamide", "[SX4](=O)(=O)[NX3H1,NX3H2]"),
    ("phenol", "[OX2H1]c"),
]

_BASE_SMARTS: list[tuple[str, str]] = [
    ("imidazole", "[nH]1ccnc1"),
    ("pyridine", "[$([n;r6;!$([nH]);!$(n~C=O);!$(n~c=O);!$(n~n)])]"),
    ("amine_secondary", "[NH1;!$(NC=O);!$(NS=O);!$([N]=[*]);!$([n])]"),
    ("amine_primary", "[NH2;!$(NC=O);!$(NS=O);!$([N]=[*]);!$([n])]"),
    ("amine_tertiary", "[NX3;H0;!$(NC=O);!$(NS=O);!$([N]=[*]);!$([n]);!$(N[OH])]([#6])([#6])[#6]"),
]


@dataclass(frozen=True)
class PKaResult:
    """Structured pKa prediction output.

    Attributes:
        pka_acid: Most acidic protonation constant detected, or ``None``.
        pka_base: Most basic protonation constant detected, or ``None``.
        compound_type: One of ``"neutral"``, ``"acid"``, ``"base"``,
            ``"zwitterion"``.
        detected_acid_groups: Names of SMARTS patterns that matched as
            acids (useful for debugging and audit trail).
        detected_base_groups: Names of SMARTS patterns that matched as
            bases.
    """

    pka_acid: float | None
    pka_base: float | None
    compound_type: str
    detected_acid_groups: tuple[str, ...] = field(default_factory=tuple)
    detected_base_groups: tuple[str, ...] = field(default_factory=tuple)

    VALID_TYPES: ClassVar[frozenset[str]] = frozenset(
        {"neutral", "acid", "base", "zwitterion"}
    )


class _SmartsCache:
    """Lazily compiled SMARTS patterns (shared across calls)."""

    _compiled_acid: list[tuple[str, Chem.Mol]] | None = None
    _compiled_base: list[tuple[str, Chem.Mol]] | None = None

    @classmethod
    def acid_patterns(cls) -> list[tuple[str, Chem.Mol]]:
        if cls._compiled_acid is None:
            cls._compiled_acid = [
                (name, Chem.MolFromSmarts(smarts))
                for name, smarts in _ACID_SMARTS
            ]
        return cls._compiled_acid

    @classmethod
    def base_patterns(cls) -> list[tuple[str, Chem.Mol]]:
        if cls._compiled_base is None:
            cls._compiled_base = [
                (name, Chem.MolFromSmarts(smarts))
                for name, smarts in _BASE_SMARTS
            ]
        return cls._compiled_base


def _detect_acid_groups(mol: Chem.Mol) -> list[str]:
    hits: list[str] = []
    for name, pattern in _SmartsCache.acid_patterns():
        if pattern is not None and mol.HasSubstructMatch(pattern):
            hits.append(name)
    return hits


def _detect_base_groups(mol: Chem.Mol) -> list[str]:
    hits: list[str] = []
    for name, pattern in _SmartsCache.base_patterns():
        if pattern is not None and mol.HasSubstructMatch(pattern):
            hits.append(name)
    return hits


def _classify(pka_acid: float | None, pka_base: float | None) -> str:
    """Sisyphus Henderson-Hasselbalch classification thresholds."""
    is_acidic = pka_acid is not None and pka_acid < 7.0
    is_basic = pka_base is not None and pka_base > 8.0
    if is_acidic and is_basic:
        return "zwitterion"
    if is_acidic:
        return "acid"
    if is_basic:
        return "base"
    return "neutral"


def predict_pka(smiles: str) -> PKaResult:
    """Predict pKa_acid, pKa_base, and compound_type from SMILES.

    The algorithm scans for ionizable functional groups via SMARTS and
    assigns the literature-average pKa for the first match in each
    category. If multiple acid or base groups are detected, the dominant
    one (lowest pKa for acids, highest pKa for bases) is reported.

    Args:
        smiles: Input SMILES.

    Returns:
        ``PKaResult`` with predicted pKa values and compound_type.

    Raises:
        ValueError: If ``smiles`` is empty, whitespace-only, or not
            parseable by RDKit.
    """
    if not smiles or not smiles.strip():
        raise ValueError(f"Invalid SMILES: {smiles!r}")
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles!r}")

    acid_hits = _detect_acid_groups(mol)
    base_hits = _detect_base_groups(mol)

    pka_acid: float | None = None
    if acid_hits:
        # Dominant acid = lowest pKa (most ionized at pH 7.4).
        pka_acid = min(GROUP_PKA[name] for name in acid_hits)

    pka_base: float | None = None
    if base_hits:
        # Dominant base = highest pKa (most protonated at pH 7.4).
        pka_base = max(GROUP_PKA[name] for name in base_hits)

    compound_type = _classify(pka_acid, pka_base)

    return PKaResult(
        pka_acid=pka_acid,
        pka_base=pka_base,
        compound_type=compound_type,
        detected_acid_groups=tuple(acid_hits),
        detected_base_groups=tuple(base_hits),
    )
