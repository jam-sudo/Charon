"""Layer 0: Input validation, drug-likeness checks, and applicability domain.

All guardrails are ADVISORY -- they generate warnings but never block
execution.  Only a genuinely unparseable SMILES sets ``is_valid=False``.
"""

from __future__ import annotations

import logging

from rdkit import DataStructs

from charon.core.molecule import Molecule
from charon.core.schema import GuardrailWarning, ValidationResult

logger = logging.getLogger(__name__)


class GuardrailChecker:
    """Layer 0: Validate molecular input, check drug-likeness, assess applicability domain."""

    def __init__(self, reference_fps: list | None = None) -> None:
        """
        Parameters
        ----------
        reference_fps : list or None
            List of Morgan fingerprint BitVects from the training set.
            If ``None``, applicability domain is reported as ``"UNKNOWN"``.
        """
        self._reference_fps = reference_fps

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def check(self, smiles: str) -> ValidationResult:
        """Run all guardrails on input SMILES.

        Returns a :class:`ValidationResult` with ``is_valid``, ``warnings``,
        ``applicability_domain``, and computed descriptors.

        Guardrails are **advisory**: drug-likeness violations produce
        warnings, never exceptions.  Only a truly invalid SMILES sets
        ``is_valid=False``.
        """
        warnings: list[GuardrailWarning] = []

        # ----------------------------------------------------------
        # 1. Structural validity
        # ----------------------------------------------------------
        try:
            mol = Molecule(smiles)
        except ValueError as exc:
            return ValidationResult(
                is_valid=False,
                smiles_canonical=smiles,
                molecular_weight=0.0,
                warnings=[
                    GuardrailWarning(
                        category="structural_validity",
                        message=f"Invalid SMILES: {exc}",
                        severity="CRITICAL",
                    )
                ],
                applicability_domain="UNKNOWN",
                descriptors={},
            )

        desc = mol.descriptors()
        canonical = mol.canonical_smiles

        # ----------------------------------------------------------
        # 2. Lipinski Rule of 5
        # ----------------------------------------------------------
        lipinski_violations = 0
        if desc["MW"] > 500:
            lipinski_violations += 1
        if desc["logP"] > 5:
            lipinski_violations += 1
        if desc["HBD"] > 5:
            lipinski_violations += 1
        if desc["HBA"] > 10:
            lipinski_violations += 1

        if lipinski_violations >= 2:
            warnings.append(
                GuardrailWarning(
                    category="lipinski",
                    message=(
                        f"Lipinski Rule of 5: {lipinski_violations} violations "
                        f"(MW={desc['MW']:.1f}, logP={desc['logP']:.2f}, "
                        f"HBD={desc['HBD']:.0f}, HBA={desc['HBA']:.0f})"
                    ),
                    severity="WARNING",
                )
            )

        # ----------------------------------------------------------
        # 3. Veber rules
        # ----------------------------------------------------------
        if desc["rotatable_bonds"] > 10:
            warnings.append(
                GuardrailWarning(
                    category="veber",
                    message=(
                        f"Veber rule violation: rotatable bonds = "
                        f"{desc['rotatable_bonds']:.0f} (limit 10)"
                    ),
                    severity="WARNING",
                )
            )
        if desc["TPSA"] > 140:
            warnings.append(
                GuardrailWarning(
                    category="veber",
                    message=(
                        f"Veber rule violation: TPSA = {desc['TPSA']:.1f} "
                        f"(limit 140)"
                    ),
                    severity="WARNING",
                )
            )

        # ----------------------------------------------------------
        # 4. Known failure modes
        # ----------------------------------------------------------
        if desc["MW"] > 800:
            warnings.append(
                GuardrailWarning(
                    category="failure_mode",
                    message=(
                        "Approaching biologic territory; small-molecule "
                        "PBPK may be unsuitable"
                    ),
                    severity="CRITICAL",
                )
            )

        if desc["logP"] > 6:
            warnings.append(
                GuardrailWarning(
                    category="failure_mode",
                    message=(
                        "High lipophilicity; dissolution-limited risk, "
                        "formulation critical"
                    ),
                    severity="WARNING",
                )
            )

        if desc["TPSA"] > 140:
            warnings.append(
                GuardrailWarning(
                    category="failure_mode",
                    message=(
                        "Low permeability predicted; oral route may be "
                        "questionable"
                    ),
                    severity="WARNING",
                )
            )

        # ----------------------------------------------------------
        # 5. Applicability domain
        # ----------------------------------------------------------
        ad = self._assess_applicability_domain(mol)

        return ValidationResult(
            is_valid=True,
            smiles_canonical=canonical,
            molecular_weight=desc["MW"],
            warnings=warnings,
            applicability_domain=ad,
            descriptors=desc,
        )

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _assess_applicability_domain(self, mol: Molecule) -> str:
        """Compute applicability domain from Tanimoto similarity to reference set."""
        if self._reference_fps is None or len(self._reference_fps) == 0:
            return "UNKNOWN"

        query_fp = mol.morgan_fingerprint()
        max_sim = max(
            DataStructs.TanimotoSimilarity(query_fp, ref_fp)
            for ref_fp in self._reference_fps
        )

        if max_sim > 0.5:
            return "HIGH"
        elif max_sim >= 0.3:
            return "MODERATE"
        else:
            return "LOW"
