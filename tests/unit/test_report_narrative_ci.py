"""Sprint 7 Task 4: verify the ADME table renders a 90% CI column.

The narrative ADME renderer takes a :class:`ReportData` object whose
``properties`` mapping holds entry dicts with ``ci_lower``/``ci_upper`` keys
(populated by :mod:`charon.report.collector` from conformal bounds).

Acceptance:
1. When an entry has both ``ci_lower`` and ``ci_upper`` non-``None``, the
   rendered Markdown contains a bracketed interval like ``[0.04, 0.25]`` and
   the header includes a ``90% CI`` column.
2. When either bound is ``None``, the row's CI cell shows a placeholder
   (em dash / ``-``) and no bracketed interval appears in that row.
"""

from __future__ import annotations

from charon.report.collector import ReportData
from charon.report.narrative import _render_adme_table


def _adme_entry(
    value: float,
    ci_lo: float | None,
    ci_hi: float | None,
    *,
    unit: str = "fraction",
    source: str = "ml_ensemble",
    flag: str | None = None,
) -> dict:
    """Build an entry shaped like :func:`collector._flatten_property`."""
    return {
        "value": value,
        "ci_lower": ci_lo,
        "ci_upper": ci_hi,
        "unit": unit,
        "source": source,
        "flag": flag,
        "method": None,
    }


def _make_report(properties: dict[str, dict]) -> ReportData:
    return ReportData(
        compound_name="test",
        smiles="CCO",
        molecular_weight=46.07,
        source="experimental",
        compound_type="neutral",
        properties=properties,
        ivive_summary={},
        pk_params={},
        pk_table=[],
        route="oral",
        dose_mg=1.0,
        duration_h=24.0,
        dose_recommendation=None,
        uncertainty=None,
        warnings=[],
        metadata={"species": "human"},
        timestamp="2026-04-23T00:00:00+00:00",
        charon_version="0.1.0",
    )


class TestCIRendering:
    def test_renders_ci_when_bounds_present(self):
        data = _make_report(
            {"fu_p": _adme_entry(0.10, 0.04, 0.25)}
        )
        md = _render_adme_table(data)

        # header contains 90% CI column
        assert "90% CI" in md
        # point estimate rendered
        assert "0.1" in md
        # bracketed interval appears in the fu_p row
        row = next(
            ln for ln in md.splitlines() if ln.startswith("| fu_p ")
        )
        assert "[" in row and "]" in row
        # bounds present with expected formatting
        assert "0.04" in row
        assert "0.25" in row
        assert "[0.04, 0.25]" in row

    def test_omits_ci_when_missing(self):
        data = _make_report(
            {"fu_p": _adme_entry(0.10, None, None)}
        )
        md = _render_adme_table(data)

        # header column still present (layout stable)
        assert "90% CI" in md
        # value still rendered
        assert "0.1" in md

        row = next(
            ln for ln in md.splitlines() if ln.startswith("| fu_p ")
        )
        # no bracketed interval in the fu_p row
        assert "[" not in row
        assert "]" not in row
        # CI cell is the third pipe-separated field -> placeholder dash
        cells = [c.strip() for c in row.split("|")[1:-1]]
        assert cells[2] == "-"

    def test_omits_ci_when_only_one_bound_missing(self):
        """If either bound is None we must not render a lopsided bracket."""
        data = _make_report(
            {"fu_p": _adme_entry(0.10, 0.04, None)}
        )
        md = _render_adme_table(data)

        row = next(
            ln for ln in md.splitlines() if ln.startswith("| fu_p ")
        )
        assert "[" not in row
        assert "]" not in row
        cells = [c.strip() for c in row.split("|")[1:-1]]
        assert cells[2] == "-"

    def test_header_declares_ci_column_between_value_and_unit(self):
        """The CI column must be adjacent to the Value column."""
        data = _make_report(
            {"logp": _adme_entry(3.89, 3.4, 4.3, unit="log")}
        )
        md = _render_adme_table(data)

        # Locate header row containing the column names
        header = next(
            ln for ln in md.splitlines()
            if ln.startswith("| Property") and "90% CI" in ln
        )
        cols = [c.strip() for c in header.split("|")[1:-1]]
        # Expect: Property | Value | 90% CI | Unit | Source | Flag
        assert cols[0] == "Property"
        assert cols[1] == "Value"
        assert cols[2] == "90% CI"
        assert cols[3] == "Unit"
