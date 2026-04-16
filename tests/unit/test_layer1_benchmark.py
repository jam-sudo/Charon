"""Unit tests for validation/benchmarks/layer1_admet.py.

Uses a mini in-memory CSV fixture (written to a tmp file) so that we do not
depend on the full 153-row calibration set or slow model inference for all
compounds.

Fixture rows
------------
1. ethanol (CCO)          — valid SMILES, all obs values present
2. benzene (c1ccccc1)     — valid SMILES, all obs values present
3. "INVALID_XYZ"          — unparseable SMILES → must appear in excluded
4. "zero_fup_drug"        — fup=0.0 → must be excluded from fu_p metric
                            (fails the min-positive guard 1e-4)
"""

from __future__ import annotations

import csv
import sys
import textwrap
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from validation.benchmarks.layer1_admet import run_benchmark  # noqa: E402


# ---------------------------------------------------------------------------
# Mini CSV fixture
# ---------------------------------------------------------------------------

_MINI_CSV_CONTENT = textwrap.dedent("""\
    name,smiles,mw,logP,fup,rbp,clint_3a4_uL_min_pmol,peff_cm_s
    ethanol,CCO,46.1,0.0,0.80,0.90,1.0,2e-4
    benzene,c1ccccc1,78.1,2.0,0.50,1.10,2.0,3e-4
    bad_smiles,INVALID_XYZ,100.0,1.0,0.30,0.80,1.0,1e-4
    zero_fup_drug,CC(C)C,72.1,2.0,0.0,0.80,1.0,2e-4
""")


@pytest.fixture()
def mini_csv(tmp_path: Path) -> Path:
    """Write the mini CSV fixture to a temp file and return its path."""
    csv_file = tmp_path / "mini_adme.csv"
    csv_file.write_text(_MINI_CSV_CONTENT, encoding="utf-8")
    return csv_file


@pytest.fixture()
def benchmark_payload(mini_csv: Path, tmp_path: Path) -> dict:
    """Run the benchmark once against the mini CSV and return payload."""
    payload = run_benchmark(
        csv_path=mini_csv,
        reports_dir=tmp_path / "reports",
    )
    return payload


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


class TestRunReturnShape:
    """Test 1: run_benchmark returns a dict with all required payload keys."""

    def test_run_returns_payload_dict(self, benchmark_payload: dict) -> None:
        """run_benchmark must return a dict containing the required report keys."""
        required_keys = {"title", "panel", "date_utc", "summary", "rows", "notes"}
        for key in required_keys:
            assert key in benchmark_payload, f"Missing key in payload: '{key}'"

    def test_payload_summary_is_list(self, benchmark_payload: dict) -> None:
        assert isinstance(benchmark_payload["summary"], list)

    def test_payload_rows_is_list(self, benchmark_payload: dict) -> None:
        assert isinstance(benchmark_payload["rows"], list)

    def test_payload_excluded_is_list(self, benchmark_payload: dict) -> None:
        assert isinstance(benchmark_payload.get("excluded", []), list)


class TestInvalidSmilesExcluded:
    """Test 2: a row with an invalid SMILES appears in payload['excluded']."""

    def test_invalid_smiles_excluded(self, benchmark_payload: dict) -> None:
        """'bad_smiles' row (SMILES='INVALID_XYZ') must appear in excluded list."""
        excluded_names = [e["name"] for e in benchmark_payload.get("excluded", [])]
        assert "bad_smiles" in excluded_names, (
            f"Expected 'bad_smiles' in excluded, got: {excluded_names}"
        )

    def test_invalid_smiles_not_in_results(self, benchmark_payload: dict) -> None:
        """The invalid-SMILES compound must NOT appear in the result rows."""
        result_names = [r["name"] for r in benchmark_payload.get("rows", [])]
        assert "bad_smiles" not in result_names


class TestZeroFupExcludedFromMetric:
    """Test 3: compound with fup=0.0 is excluded from the fu_p metric."""

    def test_zero_fup_excluded_from_fup_metric(self, benchmark_payload: dict) -> None:
        """The fu_p summary 'n' must not count the zero_fup_drug row.

        The mini CSV has:
          - 2 valid compounds (ethanol, benzene) with fup=0.80 and 0.50
          - 1 invalid SMILES (excluded outright)
          - 1 zero_fup_drug with fup=0.0 (excluded by min-positive guard)

        Therefore fu_p n must be <= 2 (may be lower if predictions are also
        below the guard, but certainly NOT 3 or more).
        """
        fup_summary = next(
            (s for s in benchmark_payload["summary"] if s.get("property") == "fu_p"),
            None,
        )
        assert fup_summary is not None, "fu_p entry not found in summary"
        n_fup = fup_summary.get("n", 0)
        # zero_fup_drug (fup=0.0) must not be counted; max valid pairs = 2
        assert n_fup <= 2, (
            f"fu_p n={n_fup} > 2: zero_fup_drug was not excluded by min-positive guard"
        )


class TestLogPUsesLinearMetric:
    """Test 4: logP summary uses linear metrics (MAE), not AAFE."""

    def test_logp_uses_mae_not_aafe(self, benchmark_payload: dict) -> None:
        """logP summary entry must contain 'mae' key and must NOT contain 'aafe'."""
        logp_summary = next(
            (s for s in benchmark_payload["summary"] if s.get("property") == "logP"),
            None,
        )
        assert logp_summary is not None, "logP entry not found in summary"
        if logp_summary.get("n", 0) > 0:
            assert "mae" in logp_summary, (
                "logP summary must have 'mae' key (linear metric)"
            )
            assert "aafe" not in logp_summary, (
                "logP summary must NOT have 'aafe' key (log metric)"
            )


class TestReportFilesWritten:
    """Test 5: emit_report writes .md and .json files."""

    def test_report_files_created(self, mini_csv: Path, tmp_path: Path) -> None:
        """run_benchmark must produce both .md and .json report files."""
        reports_dir = tmp_path / "reports"
        run_benchmark(csv_path=mini_csv, reports_dir=reports_dir)
        md_file = reports_dir / "layer1_admet.md"
        json_file = reports_dir / "layer1_admet.json"
        assert md_file.exists(), f"Expected markdown report at {md_file}"
        assert json_file.exists(), f"Expected JSON report at {json_file}"
