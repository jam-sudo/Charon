"""End-to-end: `charon report` on AD=LOW compound contains CRITICAL section."""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

THEOPHYLLINE = "Cn1c(=O)c2[nH]cnc2n(C)c1=O"
REPO_ROOT = Path(__file__).resolve().parents[2]


class TestTier3Report:
    def test_markdown_has_critical_section(self, tmp_path):
        out = tmp_path / "theo"
        subprocess.run(
            [
                sys.executable, "-m", "charon.cli.main", "report",
                THEOPHYLLINE,
                "--route", "iv_bolus",
                "--dose", "1.0",
                "--noael", "5", "--noael-species", "rat",
                "--output", str(out.with_suffix(".md")),
            ],
            check=True,
            capture_output=True,
        )
        md = out.with_suffix(".md").read_text()
        assert "## CRITICAL warnings" in md
        assert "clint_uL_min_mg" in md or "clint" in md.lower()
        assert "experimental measurement essential" in md

    def test_json_carries_classifier_probs(self, tmp_path):
        out = tmp_path / "theo"
        subprocess.run(
            [
                sys.executable, "-m", "charon.cli.main", "report",
                THEOPHYLLINE,
                "--route", "iv_bolus",
                "--dose", "1.0",
                "--noael", "5", "--noael-species", "rat",
                "--output", str(out.with_suffix(".md")),
            ],
            check=True,
        )
        data = json.loads(out.with_suffix(".json").read_text())

        def _walk(obj):
            if isinstance(obj, dict):
                if obj.get("source") == "classification":
                    yield obj
                for v in obj.values():
                    yield from _walk(v)
            elif isinstance(obj, list):
                for v in obj:
                    yield from _walk(v)

        classification_entries = list(_walk(data))
        assert classification_entries, "Expected at least one classification-sourced entry"
        entry = classification_entries[0]
        probs = entry["classifier_probs"]
        assert set(probs) == {"low", "med", "high"}
        assert abs(sum(probs.values()) - 1.0) < 0.01
