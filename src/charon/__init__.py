"""Charon — open-source translational PK platform.

Phase A public API entry points:

    Pipeline, PipelineResult
        End-to-end SMILES → PK prediction (Sprint 3: IV only).
"""

from __future__ import annotations

from charon.pipeline import Pipeline, PipelineResult

__all__ = ["Pipeline", "PipelineResult"]
