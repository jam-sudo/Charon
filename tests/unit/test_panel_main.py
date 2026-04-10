"""Smoke test: layer2_human_pk.main() runs over the real panel."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from validation.benchmarks.layer2_human_pk import (  # noqa: E402
    main,
    run_benchmark,
)


class TestBenchmarkMain:
    def test_default_panel_runs(self, capsys):
        """Run against the real panel.yaml — exit 0 if all strict pass."""
        exit_code = main()
        captured = capsys.readouterr()
        # Should print both sections
        assert "R&R only" in captured.out or "R&R only" in captured.err
        assert "override" in captured.out.lower()
        # theophylline is the sole strict-gate compound — must pass
        assert exit_code == 0, f"benchmark failed, output:\n{captured.out}"

    def test_run_benchmark_returns_two_summaries(self):
        """Programmatic API returns (summaries, rows) tuple."""
        panel_path = (
            REPO_ROOT / "validation" / "data" / "tier1_obach" / "panel.yaml"
        )
        summaries, rows = run_benchmark(panel_path)
        assert "no_override" in summaries
        assert "with_override" in summaries
        assert summaries["no_override"].n == 10
        assert summaries["with_override"].n == 10
        assert len(rows["no_override"]) == 10
        assert len(rows["with_override"]) == 10
