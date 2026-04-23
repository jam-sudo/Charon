"""Narrative: CRITICAL warnings section when a property is classification-sourced."""
from __future__ import annotations

from charon.report.narrative import render_critical_warnings


def _classification_entry():
    return {
        "name": "clint_uL_min_mg",
        "value": 3.0,
        "ci_lower": 0.1,
        "ci_upper": 10.0,
        "unit": "uL/min/10^6 cells",
        "source": "classification",
        "flag": "clint_tier3_classification; CRITICAL — experimental measurement essential; AD_max=0.286",
        "classifier_probs": {"low": 0.82, "med": 0.08, "high": 0.10},
    }


class TestCriticalWarningsSection:
    def test_emits_section_when_any_entry_is_classification(self):
        md = render_critical_warnings([_classification_entry()])
        assert "## CRITICAL warnings" in md
        assert "clint_uL_min_mg" in md
        assert "experimental measurement essential" in md
        low_label_present = "Low: 0.82" in md or "low: 0.82" in md.lower()
        assert low_label_present

    def test_empty_section_when_no_classification_entries(self):
        entry_tier2 = {
            "name": "fu_p", "value": 0.1, "ci_lower": 0.05, "ci_upper": 0.2,
            "source": "ml_ensemble", "flag": None, "classifier_probs": None,
        }
        md = render_critical_warnings([entry_tier2])
        assert md == "" or "## CRITICAL warnings" not in md
