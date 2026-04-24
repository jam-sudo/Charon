"""Schema & coverage tests for the Tier A bioavailability CSV."""

from __future__ import annotations

import csv
from pathlib import Path

import pytest
import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]
CSV_PATH = REPO_ROOT / "validation" / "data" / "fih_reference" / "bioavailability.csv"
PANEL_PATH = REPO_ROOT / "validation" / "data" / "fih_reference" / "panel.yaml"

REQUIRED_COLS = {
    "compound",
    "fih_reference_route",
    "f_oral",
    "f_source",
    "f_doi_or_pmid",
    "notes",
}


def _load_rows() -> list[dict]:
    with CSV_PATH.open() as fh:
        return list(csv.DictReader(fh))


def test_csv_exists():
    assert CSV_PATH.exists(), f"missing {CSV_PATH}"


def test_csv_columns_match_schema():
    rows = _load_rows()
    assert rows, "CSV has no rows"
    assert set(rows[0].keys()) == REQUIRED_COLS


def test_csv_covers_all_tier_a_compounds():
    panel = yaml.safe_load(PANEL_PATH.read_text())["panel"]
    tier_a = {c["name"] for c in panel["compounds"] if c["tier"] == "gold"}
    csv_compounds = {r["compound"] for r in _load_rows()}
    missing = tier_a - csv_compounds
    assert not missing, f"Tier A compounds missing from bioavailability.csv: {sorted(missing)}"


def test_csv_route_values_constrained():
    allowed = {"oral", "iv"}
    for row in _load_rows():
        assert row["fih_reference_route"] in allowed, (
            f"{row['compound']}: unexpected route "
            f"{row['fih_reference_route']!r}; allowed={sorted(allowed)}"
        )


def test_csv_f_oral_numeric_or_blank_iv():
    for row in _load_rows():
        f_raw = row["f_oral"].strip()
        if row["fih_reference_route"] == "iv":
            assert f_raw == "", (
                f"{row['compound']}: IV-reference rows must have blank f_oral, "
                f"got {f_raw!r}"
            )
        else:
            value = float(f_raw)
            assert 0.0 < value <= 1.0, (
                f"{row['compound']}: f_oral must be in (0, 1], got {value}"
            )


def test_csv_every_oral_row_has_citation():
    for row in _load_rows():
        if row["fih_reference_route"] == "oral":
            assert row["f_source"].strip(), f"{row['compound']}: f_source blank"
            assert row["f_doi_or_pmid"].strip(), (
                f"{row['compound']}: f_doi_or_pmid blank"
            )
