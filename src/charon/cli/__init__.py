"""Charon CLI entry point package.

The public entry point is :func:`charon.cli.main.main`.  We intentionally
do NOT re-export it at the package level so that ``charon.cli.main`` keeps
resolving to the submodule (needed for ``monkeypatch.setattr`` in tests).
"""
