# Sprint 15 — propranolol CYP2D6 IVIVE Correction Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Audit propranolol's F-decomposition empirically (Task 1), then conditionally apply `hepatic_clint_multiplier` to propranolol.yaml with literature citation (Tasks 2-5). Target: fold 28.55x → ≤3x (Branch A) or substantially improved (Branch B) or null result (Branch C).

**Architecture:** Audit-first sprint (mirrors Sprint 14 honesty discipline) followed by Sprint 12/13 multiplier-template application. Reuses all Sprint 12 infrastructure (`MetabolismProperties.hepatic_clint_multiplier` field + `ParameterBridge.clint_multiplier` kwarg + `ode_compiler.py` pass-through). No code changes in `src/charon/`.

**Tech Stack:** pytest, PyYAML, copy.deepcopy for in-memory YAML manipulation. WebSearch for literature verification.

**Spec:** `docs/superpowers/specs/2026-04-24-sprint15-cyp2d6-propranolol-design.md`

---

## File Structure

| File | Change |
|---|---|
| `scripts/sprint15_audit.py` | NEW — F-decomposition + multiplier sweep audit script (~120 LOC) |
| `validation/reports/layer3_ivive_decomposition.md` | Always — append `## §10 Sprint 15 audit (propranolol F-decomposition)` |
| `validation/data/tier1_obach/compounds/propranolol.yaml` | Conditional (Branches A/B) — add `hepatic_clint_multiplier` block |
| `tests/integration/test_propranolol_cyp2d6_enhancement.py` | Conditional (Branches A/B) NEW — before/after test (~115 LOC; mirrors `test_atorvastatin_oatp_enhancement.py`) |
| `validation/reports/layer3_fih_dose.{md,json}` | Conditional (Branches A/B) — regenerated with Sprint 15 narrative |
| `validation/reports/layer3_ivive_decomposition.{md,json}` | Conditional (Branches A/B) — regenerated |
| `docs/superpowers/sprint10-ivive-bias-ticket.md` | Always — append Sprint 15 status section + Sprint 14 reconciliation language |

**Unchanged (Sprint 12 infrastructure stays):**
- `src/charon/core/schema.py`, `parameter_bridge.py`, `pipeline.py`, `ode_compiler.py`
- All other compound YAMLs
- `validation/data/fih_reference/panel.yaml`
- All other tests

---

## Current propranolol values (from this worktree)

**panel.yaml (tier: gold) entry:**
- `reference_fih_mg: 10.0`
- `route: oral`
- `target_ceff_nM: <as configured>` (audit will print)
- Source: Wood 1978 / NDA 16-418 (Inderal label)

**propranolol.yaml current relevant fields:**
- `fu_p: 0.13` (Obach 1999 Table 2)
- `fu_inc: 0.96` (Obach 1999)
- `bp_ratio: 0.80` (Obach 1999)
- `clint_uL_min_mg: 13.6` (Obach 1999)
- `clrenal_L_h: 0.1` (Obach 1999)
- `peff_cm_s: 2.9e-4` (Lennernäs 1997 / Winiwarter 1998)
- `hepatic_clint_multiplier`: NOT PRESENT (Task 2 will add if Branch A/B)

**Sprint 14 outcome:** propranolol oral MRSD = 0.3502 mg, ref 10 mg, fold 28.55x.

---

## Task 1: F-decomposition audit + multiplier sweep

**Files:**
- Create: `scripts/sprint15_audit.py`
- Modify: `validation/reports/layer3_ivive_decomposition.md` (append §10)

**This task is research-only — captures empirical evidence to validate or reject the brainstorming hypothesis (Fa~0.96, Fg~1.0, Fh~0.93, F_oral~0.89) and produce the multiplier sweep curve.**

- [ ] **Step 1: Create the audit script**

Create `scripts/sprint15_audit.py` with this exact content:

```python
"""Sprint 15 audit — propranolol F-decomposition + multiplier sweep.

Purpose:
1. Capture Fa/Fg/Fh/F_oral from current Charon pipeline for propranolol
   (no multiplier applied) — validate against analytical estimates.
2. Run an empirical multiplier sweep (m in {1, 2, 3, 5, 8, 12, 15, 20, 25, 30})
   to produce the fold-vs-m curve. Identify m_close_3x and m_close_1x.

Output: prints two markdown tables to stdout. Pipe to file or copy into the
audit subsection of validation/reports/layer3_ivive_decomposition.md.

Usage:
    python3 scripts/sprint15_audit.py > /tmp/sprint15_audit.md
"""

from __future__ import annotations

import copy
import sys
from pathlib import Path
from typing import Any

import yaml

from charon import Pipeline
from charon.core.schema import CompoundConfig, DoseProjectionConfig

REPO_ROOT = Path(__file__).resolve().parents[1]
PROPRA_YAML = REPO_ROOT / "validation" / "data" / "tier1_obach" / "compounds" / "propranolol.yaml"
PANEL_YAML = REPO_ROOT / "validation" / "data" / "fih_reference" / "panel.yaml"


def load_propranolol_entry() -> dict[str, Any]:
    panel = yaml.safe_load(PANEL_YAML.read_text())["panel"]
    return next(
        c for c in panel["compounds"]
        if c["name"] == "propranolol" and c["tier"] == "gold"
    )


def run_pipeline_for_compound(
    compound: CompoundConfig,
    entry: dict[str, Any],
) -> tuple[float, dict[str, Any], Any]:
    pipe = Pipeline(
        compound,
        route=entry["route"],
        dose_mg=1.0,
        dose_projection=DoseProjectionConfig(
            target_ceff_nM=float(entry["target_ceff_nM"]),
            safety_factor=10.0,
            tau_h=24.0,
        ),
    )
    result = pipe.run()
    if result.dose_recommendation is None:
        raise RuntimeError("Pipeline did not produce a dose recommendation")
    return (
        float(result.dose_recommendation.mrsd_mg),
        result.metadata,
        result.pk_parameters,
    )


def make_compound_with_multiplier(
    base_data: dict[str, Any], m: float | None
) -> CompoundConfig:
    data = copy.deepcopy(base_data)
    metab = data["properties"]["metabolism"]
    if m is None or m == 1.0:
        metab.pop("hepatic_clint_multiplier", None)
    else:
        metab["hepatic_clint_multiplier"] = {
            "value": float(m),
            "source": "experimental",
            "unit": "ratio",
            "method": f"sprint15 audit sweep m={m}",
        }
    return CompoundConfig.model_validate(data)


def main() -> int:
    base_data = yaml.safe_load(PROPRA_YAML.read_text())
    entry = load_propranolol_entry()
    ref_mg = float(entry["reference_fih_mg"])

    # --- Section 1: F-decomposition (no multiplier applied) ---
    base_compound = make_compound_with_multiplier(base_data, None)
    assert (
        base_compound.properties.metabolism.hepatic_clint_multiplier is None
    ), "audit baseline must have NO multiplier"
    mrsd_base, meta_base, pk_base = run_pipeline_for_compound(base_compound, entry)
    fold_base = ref_mg / mrsd_base if mrsd_base > 0 else float("inf")

    print("## §10 Sprint 15 audit — propranolol F-decomposition")
    print()
    print(f"**Generated:** by `scripts/sprint15_audit.py`")
    print(f"**Reference FIH dose:** {ref_mg} mg ({entry.get('f_source', 'n/a')})")
    print()
    print("### F-decomposition (no multiplier)")
    print()
    print("| Component | Value | Analytical expected | Within range? |")
    print("|---|---:|---:|:---:|")

    def in_range(v: float, lo: float, hi: float) -> str:
        return "yes" if lo <= v <= hi else "**NO**"

    print(
        f"| Fa | {pk_base.fa:.4f} | ~0.96 | {in_range(pk_base.fa, 0.85, 1.00)} |"
    )
    print(
        f"| Fg | {pk_base.fg:.4f} | ~1.00 | {in_range(pk_base.fg, 0.90, 1.00)} |"
    )
    print(
        f"| Fh | {pk_base.fh:.4f} | ~0.93 | {in_range(pk_base.fh, 0.85, 0.97)} |"
    )
    print(
        f"| F_oral | {pk_base.bioavailability:.4f} | ~0.89 | {in_range(pk_base.bioavailability, 0.78, 0.96)} |"
    )

    clint_liver_base = float(meta_base.get("clint_liver_L_h", float("nan")))
    print(
        f"| CLint_liver_L_h | {clint_liver_base:.2f} | ~51.0 | {in_range(clint_liver_base, 45.0, 60.0)} |"
    )
    print(
        f"| MRSD_base_mg | {mrsd_base:.4g} | 0.3502 (Sprint 14) | n/a |"
    )
    print(
        f"| fold_base | {fold_base:.2f} | 28.55 (Sprint 14) | n/a |"
    )
    print()

    # --- Section 2: multiplier sweep ---
    multipliers = [1, 2, 3, 5, 8, 12, 15, 20, 25, 30]
    print("### Multiplier sweep (fold vs m)")
    print()
    print("| m | MRSD_mg | fold_observed | within_3x | F_oral | CLint_liver_L_h |")
    print("|---:|---:|---:|:---:|---:|---:|")

    sweep_rows: list[tuple[int, float, float, bool]] = []
    for m in multipliers:
        compound = make_compound_with_multiplier(base_data, m)
        mrsd_m, meta_m, pk_m = run_pipeline_for_compound(compound, entry)
        fold_m = ref_mg / mrsd_m if mrsd_m > 0 else float("inf")
        within_3x = fold_m <= 3.0
        sweep_rows.append((m, mrsd_m, fold_m, within_3x))
        clint_liv = float(meta_m.get("clint_liver_L_h", float("nan")))
        F_m = pk_m.bioavailability
        flag = "yes" if within_3x else "no"
        print(
            f"| {m} | {mrsd_m:.4g} | {fold_m:.2f} | {flag} | {F_m:.3f} | {clint_liv:.2f} |"
        )

    # Identify boundary multipliers
    m_close_3x = next((m for m, _, fold, _ in sweep_rows if fold <= 3.0), None)
    m_close_1x_pair = min(sweep_rows, key=lambda r: abs(r[2] - 1.0))
    m_close_1x = m_close_1x_pair[0]

    print()
    print("### Boundary summary")
    print()
    print(f"- **m_close_3x:** {m_close_3x if m_close_3x is not None else 'NOT REACHED in sweep'}  (smallest m bringing fold ≤ 3x)")
    print(f"- **m_close_1x:** {m_close_1x}  (m bringing fold closest to 1.0)")
    print()
    print(
        "Use this empirical curve to choose the literature multiplier (Task 2). "
        "If literature supports a value ≥ m_close_3x, Branch A applies; if literature "
        "supports a value < m_close_3x but > 3, Branch B applies (close-but-not-quite); "
        "if literature is inconclusive or < 3, Branch C applies (null)."
    )

    return 0


if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Step 2: Run the audit script**

```bash
python3 scripts/sprint15_audit.py | tee /tmp/sprint15_audit.md
```

Expected: prints the F-decomposition table + sweep table + boundary summary. No errors.

If the script raises (e.g., `result.dose_recommendation is None`): STOP. Report BLOCKED. Likely cause: pipeline path differs for propranolol; debug Pipeline manually before continuing.

- [ ] **Step 3: Verify F-decomposition is consistent with hypothesis**

Read the output. Branch decision based on F-decomposition findings:

- **All four (Fa, Fg, Fh, F_oral) within range columns are "yes":** Hypothesis confirmed. Proceed to Task 2 — literature multiplier search.
- **Fa < 0.85 OR Fg < 0.90:** Hypothesis FAILED. F-gap is in absorption/gut metabolism, not CLint. **Switch to Branch C (null) immediately.** Document in narrative; flag Sprint 16 ACAT investigation. Skip Tasks 2-3, do only Task 4 (with Branch C narrative) and Task 5.
- **Fh < 0.85:** Charon's CLh is already high relative to literature → multiplier inappropriate. Switch to Branch C.
- **Fh > 0.97:** even more extreme low CL than expected; multiplier is appropriate but value will be very large. Proceed to Task 2.

- [ ] **Step 4: Append §10 to layer3_ivive_decomposition.md**

Append the audit script output (`/tmp/sprint15_audit.md`) to `validation/reports/layer3_ivive_decomposition.md`. Use the Read tool first to confirm the report's current end (do NOT overwrite §9 Sprint 13 or earlier sections).

```bash
# Append (preserve existing content)
cat /tmp/sprint15_audit.md >> validation/reports/layer3_ivive_decomposition.md
```

Verify the append:

```bash
tail -20 validation/reports/layer3_ivive_decomposition.md
```

Expected: shows the boundary summary lines from the audit.

- [ ] **Step 5: Commit audit script + report append**

```bash
git add scripts/sprint15_audit.py validation/reports/layer3_ivive_decomposition.md
git commit -m "audit(sprint15): propranolol F-decomposition + multiplier sweep"
```

---

## Task 2: Literature multiplier search + Branch decision + propranolol.yaml edit

**Files:**
- Conditional Modify (Branches A/B): `validation/data/tier1_obach/compounds/propranolol.yaml` (the `metabolism:` block)

- [ ] **Step 1: WebSearch for propranolol IVIVE literature**

Run WebSearch queries:

- `"Obach 1999 propranolol predicted observed CLint hepatic clearance"`
- `"propranolol CYP2D6 in vivo in vitro CLint ratio HLM underprediction"`
- `"Ito Houston 2005 basic drug HLM IVIVE underprediction high extraction"`
- `"Hallifax Houston 2010 well-stirred bias correction high extraction"`
- `"propranolol intrinsic clearance in vivo Walle 1985 Routledge"`

Goals:
1. Find propranolol-specific in vivo / in vitro CLint ratio (most-cited primary source).
2. If propranolol-specific is unavailable, find CYP2D6-class average HLM bias.
3. Cross-check at least 2 independent sources before applying.

Document findings with citation strings. Note: the value Charon needs depends on the audit's `m_close_3x`. If literature supports ≥ `m_close_3x`, Branch A. If 5x ≤ literature < `m_close_3x`, Branch B. If literature < 5x, Branch C.

- [ ] **Step 2: Branch decision**

Based on Step 1 findings + Task 1's `m_close_3x`:

| Literature support | Multiplier (apply) | Branch |
|---|---|---|
| ≥ 2 sources, ratio ≥ `m_close_3x` (likely 9-12 from sweep) | midpoint of cited range (e.g., 12, 15, 18) | A |
| ≥ 2 sources, ratio in [5, m_close_3x) | midpoint (e.g., 6, 7, 8) | B |
| Cited range only 3-5x or only 1 source | midpoint (e.g., 4) | B (best-effort) |
| Only generic CYP2D6 bias mentioned, no specific number | NO multiplier | C |
| Cited range < 3x or contradictory | NO multiplier | C |

**If Branch C:** Skip Step 3-5. Skip Task 3. Proceed to Task 4 with Branch C narrative.

**If Branch A or B:** Continue with Step 3.

- [ ] **Step 3: Edit propranolol.yaml metabolism block**

Read the current file first:

```bash
sed -n '20,30p' validation/data/tier1_obach/compounds/propranolol.yaml
```

(The `metabolism:` block is around line 22-23 with just `clint_uL_min_mg`.)

Edit using the Edit tool. Target structure (substitute YOUR_M, YOUR_LIT_RANGE, YOUR_CITATIONS from Step 1-2 findings):

```yaml
  metabolism:
    clint_uL_min_mg: {value: 13.6, source: experimental, unit: uL/min/mg, method: "Obach 1999 Table 2"}
    # Sprint 15: CYP2D6 + CYP1A2 IVIVE underprediction correction.
    # Propranolol is a high-extraction lipophilic base; HLM CLint
    # systematically underpredicts in vivo CLint for CYP2D6 substrates
    # of this class (Ito & Houston 2005, Hallifax & Houston 2010).
    # Literature documents in vivo / in vitro CLint ratio in <YOUR_LIT_RANGE>x
    # range for propranolol specifically. Empirical multiplier set to the
    # cited midpoint, applied at the CLint_liver level (before well-stirred
    # extraction model) per Sprint 12/13 pattern.
    hepatic_clint_multiplier:
      value: <YOUR_M>
      source: literature
      unit: ratio
      method: "<YOUR_CITATIONS> — CYP2D6/high-extraction base IVIVE gap; <YOUR_M> = midpoint of cited <YOUR_LIT_RANGE>x range"
```

Example (if literature supports 12-18x and you choose 15):

```yaml
    hepatic_clint_multiplier:
      value: 15.0
      source: literature
      unit: ratio
      method: "Obach 1999 DMD 27:1350 Table 3 / Ito & Houston 2005 Pharm Res 22:103 / Hallifax & Houston 2010 Drug Metab Pharmacokinet 25:74 — CYP2D6/high-extraction base HLM IVIVE underprediction 12-18x; 15.0 = midpoint"
```

- [ ] **Step 4: Verify YAML parses + schema validates**

```bash
python3 -c "
import yaml
from charon.core.schema import CompoundConfig
data = yaml.safe_load(open('validation/data/tier1_obach/compounds/propranolol.yaml').read())
c = CompoundConfig.model_validate(data)
m = c.properties.metabolism.hepatic_clint_multiplier
print(f'hepatic_clint_multiplier: {m.value} ({m.source})')
print(f'method: {m.method[:120]}...')
"
```

Expected: prints the multiplier value with `(literature)` source and the method string. No errors.

- [ ] **Step 5: Run existing YAML/schema tests**

```bash
pytest tests/unit/test_fih_new_compound_yamls.py tests/unit/test_compound_config.py -q
```

Expected: all pass. The Sprint 12 schema test for `hepatic_clint_multiplier` continues to pass.

- [ ] **Step 6: Commit**

```bash
git add validation/data/tier1_obach/compounds/propranolol.yaml
git commit -m "data(sprint15): propranolol hepatic_clint_multiplier=<YOUR_M> (CYP2D6 IVIVE correction)"
```

(Replace `<YOUR_M>` with actual value, e.g., `15.0`.)

---

## Task 3: Integration test — propranolol before/after multiplier (Branches A/B only)

**Skip this task if Branch C (no multiplier applied).**

**Files:**
- Create: `tests/integration/test_propranolol_cyp2d6_enhancement.py`

- [ ] **Step 1: Create the test file**

Create `tests/integration/test_propranolol_cyp2d6_enhancement.py` with this exact content:

```python
"""Sprint 15 integration test: propranolol MRSD with/without CYP2D6 enhancement.

Mirrors Sprint 12 (atorvastatin) and Sprint 13 (diclofenac) test structures.
Validates:
1. Enhanced MRSD is strictly larger than baseline.
2. Ratio MRSD_enhanced / MRSD_base is within [4, 25] — broader than Sprint 12
   atorvastatin's [4, 12] or Sprint 13 diclofenac's [2, 5] because:
   - propranolol's literature-supported multiplier range (CYP2D6 high-extraction
     base IVIVE bias) is larger and more uncertain than OATP1B1 or UGT2B7
   - both CL and F change with multiplier (well-stirred saturation), so the
     MRSD ratio is sub-linear vs the multiplier (e.g., m=15 → MRSD ratio ~15;
     m=8 → MRSD ratio ~9)
3. clint_liver_L_h ratio (enhanced/base) exactly matches the multiplier value.
"""

from __future__ import annotations

import copy
from pathlib import Path

import pytest
import yaml

from charon import Pipeline
from charon.core.schema import CompoundConfig, DoseProjectionConfig

REPO_ROOT = Path(__file__).resolve().parents[2]
PROPRA_YAML = REPO_ROOT / "validation" / "data" / "tier1_obach" / "compounds" / "propranolol.yaml"
PANEL_YAML = REPO_ROOT / "validation" / "data" / "fih_reference" / "panel.yaml"


def _load_propranolol_entry() -> dict:
    panel = yaml.safe_load(PANEL_YAML.read_text())["panel"]
    return next(
        c for c in panel["compounds"]
        if c["name"] == "propranolol" and c["tier"] == "gold"
    )


def _run_pipeline(compound: CompoundConfig, entry: dict) -> tuple[float, dict]:
    pipe = Pipeline(
        compound,
        route=entry["route"],
        dose_mg=1.0,
        dose_projection=DoseProjectionConfig(
            target_ceff_nM=float(entry["target_ceff_nM"]),
            safety_factor=10.0,
            tau_h=24.0,
        ),
    )
    result = pipe.run()
    assert result.dose_recommendation is not None
    return float(result.dose_recommendation.mrsd_mg), result.metadata


def test_propranolol_mrsd_with_multiplier_larger_than_without():
    """With clint_multiplier > 1, propranolol MRSD must be strictly larger."""
    data = yaml.safe_load(PROPRA_YAML.read_text())
    entry = _load_propranolol_entry()

    compound_enhanced = CompoundConfig.model_validate(data)
    assert compound_enhanced.properties.metabolism.hepatic_clint_multiplier is not None, (
        "propranolol.yaml must have hepatic_clint_multiplier populated (Task 2 prerequisite)"
    )
    mrsd_enhanced, _ = _run_pipeline(compound_enhanced, entry)

    data_base = copy.deepcopy(data)
    data_base["properties"]["metabolism"].pop("hepatic_clint_multiplier", None)
    compound_base = CompoundConfig.model_validate(data_base)
    assert compound_base.properties.metabolism.hepatic_clint_multiplier is None
    mrsd_base, _ = _run_pipeline(compound_base, entry)

    assert mrsd_enhanced > mrsd_base, (
        f"Enhanced MRSD ({mrsd_enhanced:.3g}) must exceed baseline ({mrsd_base:.3g})"
    )


def test_propranolol_mrsd_ratio_within_literature_range():
    """MRSD_enhanced / MRSD_base within [4, 25] — band reflects CYP2D6 literature
    uncertainty (5-20x typical). For high-extraction propranolol (post-multiplier),
    well-stirred saturation makes the MRSD ratio approximately equal to the
    multiplier in the moderate-extraction regime, decreasing slightly at very
    high multipliers."""
    data = yaml.safe_load(PROPRA_YAML.read_text())
    entry = _load_propranolol_entry()

    compound_enhanced = CompoundConfig.model_validate(data)
    mrsd_enhanced, _ = _run_pipeline(compound_enhanced, entry)

    data_base = copy.deepcopy(data)
    data_base["properties"]["metabolism"].pop("hepatic_clint_multiplier", None)
    compound_base = CompoundConfig.model_validate(data_base)
    mrsd_base, _ = _run_pipeline(compound_base, entry)

    ratio = mrsd_enhanced / mrsd_base
    assert 4.0 <= ratio <= 25.0, (
        f"MRSD ratio {ratio:.2f} outside expected [4, 25] range; "
        f"check multiplier value or extraction saturation"
    )


def test_propranolol_enhanced_clint_liver_metadata_scales_by_multiplier():
    """clint_liver_L_h ratio (enhanced/base) should exactly match the multiplier.
    Pure linear scaling — ConversionStep multiplies CLint_liver before the
    liver model, so the reported clint_liver_L_h in metadata is the enhanced
    value."""
    data = yaml.safe_load(PROPRA_YAML.read_text())
    entry = _load_propranolol_entry()

    compound_enhanced = CompoundConfig.model_validate(data)
    multiplier = compound_enhanced.properties.metabolism.hepatic_clint_multiplier.value
    _, md_enhanced = _run_pipeline(compound_enhanced, entry)

    data_base = copy.deepcopy(data)
    data_base["properties"]["metabolism"].pop("hepatic_clint_multiplier", None)
    compound_base = CompoundConfig.model_validate(data_base)
    _, md_base = _run_pipeline(compound_base, entry)

    clint_enhanced = float(md_enhanced["clint_liver_L_h"])
    clint_base = float(md_base["clint_liver_L_h"])
    ratio = clint_enhanced / clint_base
    assert abs(ratio - multiplier) < 1e-3, (
        f"clint_liver_L_h ratio {ratio:.6f} != multiplier {multiplier}"
    )
```

- [ ] **Step 2: Run the integration tests**

```bash
pytest tests/integration/test_propranolol_cyp2d6_enhancement.py -v
```

Expected: 3/3 PASS.

Possible failures:
- Test 2 `ratio < 4 or > 25`: investigate. If `ratio < 4`, multiplier may not be applied or extraction is unexpectedly saturated. If `ratio > 25`, the multiplier value itself is too large given Charon's pipeline (PAD path may amplify nonlinearly). Report BLOCKED.
- Test 3 `ratio != multiplier`: unit conversion bug or multiplier not flowing to ode_compiler. STOP, report BLOCKED.

- [ ] **Step 3: Commit**

```bash
git add tests/integration/test_propranolol_cyp2d6_enhancement.py
git commit -m "test(sprint15): propranolol CYP2D6 enhancement before/after integration"
```

---

## Task 4: Regenerate benchmarks + Sprint 15 narrative + Sprint 14 reconciliation

**Files:**
- Conditional regenerate (Branches A/B only): `validation/reports/layer3_fih_dose.{md,json}`, `validation/reports/layer3_ivive_decomposition.{md,json}`
- Always: append Sprint 15 narrative + Sprint 14 reconciliation language to report markdown files

**Prior per-compound folds (Sprint 14, baseline for Sprint 15) for the comparison table:**

| Compound | Sprint 14 fold |
|---|---:|
| midazolam | 1.46 |
| warfarin | 1.02 |
| propranolol | 28.55 |
| verapamil | 2.49 |
| omeprazole | 1.60 |
| theophylline | 1.02 |
| diclofenac | 3.10 |
| diazepam | 4.91 |
| metoprolol | 2.13 |
| acetaminophen | 1.21 |
| lisinopril | 4.13 |
| atorvastatin | 1.70 |

Sprint 14 within-3x: 8/12 = 66.7%.

- [ ] **Step 1: If Branches A/B (multiplier applied), regenerate benchmarks**

```bash
python3 validation/benchmarks/layer3_fih_dose.py
python3 validation/benchmarks/layer3_ivive_decomposition.py
```

Expected outputs:
- `[OK] Sanity floor: 12/12 pass.` (Tier B unchanged.)
- `[OK] Decomposition wrote ... (12 compounds)`

If Branch C (no multiplier): SKIP regeneration. Benchmarks unchanged from Sprint 14.

- [ ] **Step 2: Extract Sprint 15 propranolol fold + Tier A summary** (Branches A/B only)

```bash
python3 -c "
import json
d = json.load(open('validation/reports/layer3_fih_dose.json'))
s = d['summary']
print(f'Tier A within-3x: {s[\"gold_within_3x\"]}/{s[\"gold_n\"]} = {100*s[\"gold_within_3x_fraction\"]:.1f}%')
print(f'Tier A within-10x: {s[\"gold_within_10x\"]}/{s[\"gold_n\"]}')
print()
print(f'{\"compound\":<15} {\"mrsd_mg\":>10} {\"ref_mg\":>10} {\"fold\":>8} {\"3x\":>4}')
for r in d['extra_sections']['Gold (Tier A) — fold-error vs reference FIH']:
    print(f'{r[\"compound\"]:<15} {r[\"mrsd_pred_mg\"]:>10.3g} {r[\"reference_fih_mg\"]:>10.3g} {r[\"fold_error\"]:>8.2f} {str(r[\"within_3x\"])[:4]:>4}')
"
```

Record output for narrative.

- [ ] **Step 3: Extract decomposition aggregate + propranolol residual** (Branches A/B only)

```bash
python3 -c "
import json
d = json.load(open('validation/reports/layer3_ivive_decomposition.json'))
s = d['summary']
print(f'liver_model: {s[\"aggregate_pct_liver_model\"]:.1f}%')
print(f'route_bias:  {s[\"aggregate_pct_route_bias\"]:.1f}%')
print(f'residual:    {s[\"aggregate_pct_residual\"]:.1f}%')
print()
print('Propranolol decomposition:')
for r in d['extra_sections']['Per-compound decomposition']:
    if r['compound'] == 'propranolol':
        for k in ('mrsd_ws_mg', 'reference_fih_mg', 'fold_observed', 'fold_residual'):
            v = r.get(k)
            if isinstance(v, (int, float)):
                print(f'  {k}: {v:.3f}')
            else:
                print(f'  {k}: {v}')
"
```

- [ ] **Step 4: Append Sprint 15 narrative to `validation/reports/layer3_fih_dose.md`**

Use the Read tool first to confirm existing content. Append (do NOT overwrite). Three template variants — choose based on Task 2 branch outcome:

**Branch A (within-3x closure):**

```markdown

## Sprint 15 comparison (CYP2D6 enhancement for propranolol)

Sprint 14 (8/12 = 66.7% within-3x, §8 PASSED): propranolol at 28.55x was the largest remaining residual.
Sprint 15 (propranolol `hepatic_clint_multiplier: <YOUR_M>`): Tier A within-3x = <M_COUNT>/12 = <PCT>%.

Per-compound deltas (Sprint 14 → Sprint 15):

| Compound | Sprint 14 fold | Sprint 15 fold | Δ |
|---|---:|---:|:---:|
| midazolam | 1.46 | <S15> | <~/↓/↑> |
| warfarin | 1.02 | <S15> | <...> |
| propranolol | 28.55 | <S15_PROPRA> | <PROPRA_DELTA> |
| verapamil | 2.49 | <S15> | <...> |
| omeprazole | 1.60 | <S15> | <...> |
| theophylline | 1.02 | <S15> | <...> |
| diclofenac | 3.10 | <S15> | <...> |
| diazepam | 4.91 | <S15> | <...> |
| metoprolol | 2.13 | <S15> | <...> |
| acetaminophen | 1.21 | <S15> | <...> |
| lisinopril | 4.13 | <S15> | <...> |
| atorvastatin | 1.70 | <S15> | <...> |

**Propranolol-only:**
- Sprint 14: MRSD = 0.3502 mg, fold 28.55x (outside 3x)
- Sprint 15: MRSD = <S15_PROPRA_MRSD> mg, fold <S15_PROPRA_FOLD>x (within 3x)

Multiplier applied: <YOUR_M> (<YOUR_CITATIONS> midpoint for CYP2D6/high-extraction base IVIVE gap).
Other 11 compounds expected zero delta (multiplier scoped to propranolol).

### Sprint 14 ticket reconciliation

Sprint 14 (diazepam audit) ticket characterized propranolol's residual as *"ACAT oral F computation architectural gap... Not a multiplier fix."* Sprint 15 audit (Task 1, see §10 of decomposition report) explicitly captured `Fa, Fg, Fh, F_oral` from the pipeline output and found:

- Fa ≈ <FA_VAL> — ACAT functioning normally for high-permeability propranolol
- Fg ≈ <FG_VAL> — non-CYP3A4 substrate, no gut metabolism issue
- Fh ≈ <FH_VAL> — analytically traceable to CLint underprediction in HLM
- F_oral ≈ <F_VAL> (vs literature 0.26) — overprediction by ~3.4x

The gap is therefore the **same Sprint 12/13 multiplier-template pattern** (HLM CLint underprediction for IVIVE-known substrate classes), **not an ACAT architectural issue**. Sprint 14's claim was a hypothesis based on a quick mental model and is corrected here per CLAUDE.md §6.5 honesty.
```

**Branch B (close-but-not-quite):**

```markdown

## Sprint 15 comparison (CYP2D6 enhancement for propranolol — partial closure)

Sprint 14 (8/12 = 66.7% within-3x, §8 PASSED): propranolol at 28.55x was the largest remaining residual.
Sprint 15 (propranolol `hepatic_clint_multiplier: <YOUR_M>`): Tier A within-3x = 8/12 = 66.7% (unchanged; propranolol close but stays outside 3x).

**Propranolol-only:**
- Sprint 14: MRSD = 0.3502 mg, fold 28.55x (outside 3x)
- Sprint 15: MRSD = <S15_PROPRA_MRSD> mg, fold <S15_PROPRA_FOLD>x (outside 3x by <margin>)

Multiplier applied: <YOUR_M> (<YOUR_CITATIONS> midpoint for CYP2D6/high-extraction base IVIVE gap; cited range only supports up to <UPPER_BOUND>x — applying higher multiplier would violate CLAUDE.md §6.5 honesty).

**Honest interpretation:** Substantial improvement (28.55 → <S15_PROPRA_FOLD>x) but the literature-cited multiplier range is insufficient to fully close to 3x. Further closure requires either (a) finding additional primary literature with a higher cited ratio, or (b) architectural work (per-CYP2D6 ML-recalibration of CLint, or extended-clearance modeling of in vivo CYP2D6 turnover). Both are Sprint 16+ scope.

### Sprint 14 ticket reconciliation

Sprint 14 (diazepam audit) ticket characterized propranolol's residual as *"ACAT oral F computation architectural gap... Not a multiplier fix."* Sprint 15 audit (Task 1, see §10 of decomposition report) explicitly captured `Fa, Fg, Fh, F_oral` from the pipeline output and found:

- Fa ≈ <FA_VAL> — ACAT functioning normally for high-permeability propranolol
- Fg ≈ <FG_VAL> — non-CYP3A4 substrate, no gut metabolism issue
- Fh ≈ <FH_VAL> — analytically traceable to CLint underprediction in HLM
- F_oral ≈ <F_VAL> (vs literature 0.26) — overprediction by ~3.4x

The gap is therefore the **same Sprint 12/13 multiplier-template pattern** (HLM CLint underprediction for IVIVE-known substrate classes), **not an ACAT architectural issue**. Sprint 14's claim was a hypothesis based on a quick mental model and is corrected here per CLAUDE.md §6.5 honesty.
```

**Branch C (null):**

```markdown

## Sprint 15 — propranolol audit + null result

Sprint 14 (8/12 = 66.7% within-3x): propranolol at 28.55x was the largest remaining residual.
Sprint 15 audit (Task 1, see §10 of decomposition report) examined F-decomposition empirically.

**Audit findings:**
- Fa = <FA_VAL> (expected ~0.96)
- Fg = <FG_VAL> (expected ~1.0)
- Fh = <FH_VAL> (expected ~0.93)
- F_oral = <F_VAL>

**Reason for null:** <one of:>
- "Audit revealed Fa = <FA_VAL> < 0.85, indicating ACAT/absorption is the gap, not CLint. Multiplier inappropriate. Sprint 16 ACAT investigation flagged."
- "WebSearch yielded no propranolol-specific in vivo / in vitro CLint ratio with ≥ 2 sources; generic CYP2D6 bias mentioned in reviews but no quantitative value cited. Per CLAUDE.md §6.5 honesty, no multiplier applied."
- "Cited literature range only supports < 5x for propranolol IVIVE under-prediction, well below empirical `m_close_3x = <M_3X>`. Applying small multiplier would fail to meaningfully close gap; null result accepted."

§8 target stays PASSED at 8/12 = 66.7%. Propranolol residual remains framework-limited; closure requires architectural work in Sprint 16+.

### Sprint 14 ticket reconciliation

Sprint 14 ticket characterized propranolol as *"ACAT oral F computation architectural gap. Not a multiplier fix."* Sprint 15 audit confirmed the F-gap exists but determined that <multiplier inappropriate per audit findings | literature support insufficient>. The Sprint 14 framing was approximately correct in conclusion (no multiplier to apply) but not in mechanism.
```

Fill all `<...>` placeholders with actual values from Steps 2-3 output (Branches A/B) or Task 1 audit output (Branch C).

- [ ] **Step 5: Append Sprint 15 section to `validation/reports/layer3_ivive_decomposition.md`** (Branches A/B only)

For Branch C: skip — Task 1 already wrote §10 audit; no further section needed.

For Branches A/B, append `## §11. Sprint 15 — CYP2D6 correction for propranolol`:

```markdown

## §11. Sprint 15 — CYP2D6 correction for propranolol

After `hepatic_clint_multiplier: <YOUR_M>` added to propranolol.yaml (<YOUR_CITATIONS> midpoint for CYP2D6/high-extraction base IVIVE gap):

- liver_model: <NEW_LIVER>%
- route_bias:  <NEW_ROUTE>%
- residual:    <NEW_RESIDUAL>%

Propranolol per-compound:
- Sprint 14 fold_residual: 27.49 (signed: 0.03638)
- Sprint 15 fold_residual: <S15_PROPRA_RESIDUAL> (signed: <S15_PROPRA_RESIDUAL_SIGNED>)

The multiplier treats CYP2D6 + CYP1A2 hepatic clearance as an enhanced uniform path (empirical — analogous to Sprint 12 OATP1B1 atorvastatin and Sprint 13 UGT2B7 diclofenac approaches). Sprint 15 audit (§10) confirmed the F-gap traces to CLint underprediction (Fa, Fg normal; Fh too high), validating multiplier appropriateness.

Remaining largest residuals after Sprint 15:
- diazepam 4.91x — very low fu_p well-stirred sensitivity (Sprint 14 honest null; framework-limited)
- lisinopril 4.13x — non-hepatic elimination + low Peff (renal CL, multiplier inappropriate)
- diclofenac 3.10x — Sprint 13 close-but-not-quite (literature midpoint multiplier 3.5)
```

Fill placeholders.

- [ ] **Step 6: Commit**

For Branches A/B (regenerated reports + narrative):

```bash
git add validation/reports/layer3_fih_dose.md \
        validation/reports/layer3_fih_dose.json \
        validation/reports/layer3_ivive_decomposition.md \
        validation/reports/layer3_ivive_decomposition.json
git commit -m "chore(sprint15): regenerated Layer 3 reports + narrative (CYP2D6 correction)"
```

For Branch C (only narrative + audit appended):

```bash
git add validation/reports/layer3_fih_dose.md
git commit -m "docs(sprint15): null-result narrative for propranolol audit"
```

(The decomposition report's §10 audit subsection was already committed in Task 1.)

---

## Task 5: Sprint 10 ticket update + full suite green

**Files:**
- Modify: `docs/superpowers/sprint10-ivive-bias-ticket.md`

- [ ] **Step 1: Append Sprint 15 status + Sprint 14 reconciliation**

Append to `docs/superpowers/sprint10-ivive-bias-ticket.md`. Use the template matching Task 2's outcome branch:

**Branch A (closed within 3x):**

```markdown

## Sprint 15 (CYP2D6 correction for propranolol) completed — 2026-04-24

Added `hepatic_clint_multiplier: <YOUR_M>` to propranolol.yaml (<YOUR_CITATIONS> midpoint for CYP2D6/high-extraction base HLM IVIVE gap). Reuses Sprint 12 infrastructure unchanged.

**Propranolol delta:**
- Sprint 14: MRSD 0.3502 mg, fold 28.55x (outside 3x)
- Sprint 15: MRSD <S15_PROPRA_MRSD> mg, fold <S15_PROPRA_FOLD>x (WITHIN 3x ✓)

**Layer 3 Tier A within-3x progression:**
- Sprint 9:  5/12 = 41.7% (FAILED)
- Sprint 11: 7/12 = 58.3% (FAILED by 2%)
- Sprint 12: 8/12 = 66.7% (§8 PASSED)
- Sprint 13: 8/12 = 66.7% (close-but-not-quite for diclofenac)
- Sprint 14: 8/12 = 66.7% (diazepam null result)
- Sprint 15: <M_COUNT>/12 = <PCT>% (propranolol closed)

Other 11 compounds: zero delta confirmed (multiplier propranolol-scoped).

### Sprint 14 ticket reconciliation

Sprint 14 ticket stated: *"propranolol 28.55x — ACAT oral F computation architectural gap (Sprint 11 oral migration barely improved from 36.3 → 28.55; ACAT likely gives F≈0.80 vs literature 0.26). Not a multiplier fix."*

This was a hypothesis based on quick mental model, not on F-decomposition data. Sprint 15 audit (Task 1) explicitly captured Fa/Fg/Fh from the Pipeline output:

- Fa = <FA_VAL> (vs hypothesized ACAT issue → false; ACAT functioning correctly)
- Fg = <FG_VAL> (non-CYP3A4 substrate, gut metabolism negligible)
- Fh = <FH_VAL> (the actual gap source — too high, traceable to CLint underprediction)
- F_oral = <F_VAL> (overpredicted ~3.4x vs literature 0.26)

The gap is the **same Sprint 12/13 multiplier-template pattern**, not ACAT architecture. Sprint 14's claim was unverified at the time and is corrected here per CLAUDE.md §6.5 honesty.

**Remaining unresolved gaps:**
- diazepam 4.91x — very low fu_p well-stirred sensitivity (Sprint 14 honest null; framework-limited)
- lisinopril 4.13x — non-hepatic elimination (Sprint 17 candidate; renal CL refinement)
- diclofenac 3.10x — Sprint 13 close-but-not-quite

§8 status: PASSED at <PCT>% — improved from 66.7% to <PCT>%.
```

**Branch B (close-but-not-quite):**

```markdown

## Sprint 15 (CYP2D6 correction for propranolol — partial closure) completed — 2026-04-24

Added `hepatic_clint_multiplier: <YOUR_M>` to propranolol.yaml (<YOUR_CITATIONS> midpoint for CYP2D6 IVIVE gap; cited range supports only <UPPER_BOUND>x — applying higher multiplier would violate §6.5 honesty).

**Propranolol delta:**
- Sprint 14: MRSD 0.3502 mg, fold 28.55x (outside 3x)
- Sprint 15: MRSD <S15_PROPRA_MRSD> mg, fold <S15_PROPRA_FOLD>x (still outside 3x by <margin>)

**Layer 3 Tier A within-3x:** 8/12 = 66.7% (unchanged; propranolol close but does not cross).

Honest interpretation: substantial improvement (28.55 → <S15_PROPRA_FOLD>x) but literature-cited multiplier range insufficient to fully close to 3x. Further closure requires Sprint 16+ architectural work (per-CYP2D6 ML-recalibration or extended-clearance modeling).

### Sprint 14 ticket reconciliation

Sprint 14 ticket stated: *"propranolol 28.55x — ACAT oral F computation architectural gap (Sprint 11 oral migration barely improved from 36.3 → 28.55; ACAT likely gives F≈0.80 vs literature 0.26). Not a multiplier fix."*

This was a hypothesis based on quick mental model, not on F-decomposition data. Sprint 15 audit (Task 1) explicitly captured Fa/Fg/Fh from the Pipeline output:

- Fa = <FA_VAL> (vs hypothesized ACAT issue → false; ACAT functioning correctly)
- Fg = <FG_VAL> (non-CYP3A4 substrate, gut metabolism negligible)
- Fh = <FH_VAL> (the actual gap source — too high, traceable to CLint underprediction)
- F_oral = <F_VAL> (overpredicted ~3.4x vs literature 0.26)

The gap is the **same Sprint 12/13 multiplier-template pattern**, not ACAT architecture. Sprint 14's claim was unverified at the time and is corrected here per CLAUDE.md §6.5 honesty.
```

**Branch C (null):**

```markdown

## Sprint 15 (propranolol audit, null result) — 2026-04-24

F-decomposition audit (scripts/sprint15_audit.py) examined empirical Fa/Fg/Fh for propranolol. <One of:>

- "Audit revealed Fa = <FA_VAL> < 0.85, indicating ACAT/absorption issue, not CLint. Multiplier inappropriate. Sprint 16 ACAT investigation flagged for propranolol oral first-pass mechanism."
- "WebSearch yielded no propranolol-specific in vivo/in vitro CLint ratio with ≥ 2 primary sources. Generic CYP2D6 bias is mentioned in reviews but no specific quantitative ratio cited. Per CLAUDE.md §6.5, no multiplier applied."
- "Literature support only ≤ 5x for propranolol IVIVE under-prediction, far below the empirically-needed `m_close_3x = <M_3X>` for §8 closure. Inflating beyond cited range would violate §6.5; null accepted."

**Layer 3 Tier A within-3x:** 8/12 = 66.7% (unchanged from Sprint 14).

### Sprint 14 ticket reconciliation

Sprint 14 ticket framed propranolol as "ACAT architectural gap, not a multiplier fix." Sprint 15 audit confirms the F-gap exists but determined that <multiplier inappropriate per audit findings | literature support insufficient>. Sprint 14's conclusion (no multiplier) was approximately correct, though the mechanistic framing (ACAT specifically) <was confirmed by audit | was not the primary issue per audit>.

**Remaining unresolved gaps (Sprint 16+ candidates):**
- propranolol 28.55x — Sprint 16 architectural (ACAT or CYP2D6 IVIVE deep dive)
- diazepam 4.91x — Sprint 14 framework-limited
- lisinopril 4.13x — non-hepatic elimination
```

Fill all placeholders with actual values from Tasks 1, 2, 4 outputs.

- [ ] **Step 2: Run full test suite**

```bash
pytest -q
```

Expected:
- Branches A/B: **941 passed** (938 Sprint 14 baseline + 3 new propranolol integration tests)
- Branch C: **938 passed** (no new tests, no code changes)

If any test fails: STOP. Report BLOCKED. Likely flaky-test patterns from prior sprints:
- `test_point_mode_aafe_unchanged` — known flaky, retry once with `pip install -e ".[dev]"` reinstall
- Schema test failures: re-check Task 2 Step 4 YAML syntax

- [ ] **Step 3: Commit ticket**

```bash
git add docs/superpowers/sprint10-ivive-bias-ticket.md
git commit -m "docs(sprint15): Sprint 10 ticket — Sprint 15 CYP2D6 propranolol results + Sprint 14 reconciliation"
```

Check `git status` for stray `validation/reports/layer2_human_pk.*` drift. **Do NOT commit those.**

- [ ] **Step 4: Summarize commits**

```bash
git log --oneline main..HEAD
```

Expected commits:
- Branch A or B: 4 (audit + data + test + reports + ticket = 5 commits, but reports + narrative may collapse to 1) — between 4 and 6 commits
- Branch C: 2-3 (audit + null narrative + ticket)

---

## Self-Review Notes (for the implementer)

- **Audit output guides everything.** Don't pre-decide the branch. Run Task 1 first. The empirical sweep curve will show which `m` values close the gap; literature search (Task 2) determines whether such `m` is justifiable.

- **Multiplier value 15 is illustrative only.** The Task 2 step 1 WebSearch may surface specific propranolol Obach 1999 Table 3 values (or related secondary references) that point to 12, 18, or 20. Use the verified midpoint, not a default.

- **Honesty discipline (§6.5) is paramount.** If literature only supports up to 8x and Charon's `m_close_3x = 12`, this is Branch B (close-but-not-quite). Do NOT inflate the multiplier to 12 to force §8 closure. The Sprint 13 diclofenac precedent (3.5 cited midpoint, fold landed 3.10x — 0.10x outside 3x) is the model: cite honestly and report the close-but-not-quite outcome.

- **Sprint 14 ticket reconciliation MUST be in the ticket text.** The Sprint 14 hypothesis ("ACAT architectural") was unverified. Sprint 15 should explicitly correct it (Branch A/B) or partially confirm it with caveats (Branch C). Not reconciling = leaving a contradictory framing in the repo.

- **`m_close_3x` from audit may exceed literature support.** This is the Sprint 13 scenario: empirically need m≥12 to close 3x, but literature only supports 5-8x. Branch B with honest reporting.

- **Test tolerance [4, 25] in Task 3.** Wider than Sprint 12 [4,12] or Sprint 13 [2,5] — reflects CYP2D6 literature uncertainty. Expected ratio for `m=15`: ~14-15 (in moderate-extraction regime). For `m=8`: ~7-8. For `m=20`: ~18-20. The wide band accommodates the entire plausible range.

- **No other compound YAML should be touched.** Verify `git status` shows only `propranolol.yaml` modified before each commit.

- **Commit count estimates:**
  - Branch A: 5 commits (audit, data, test, reports, ticket) + merge
  - Branch B: 5 commits (same as A; only narrative differs) + merge
  - Branch C: 3 commits (audit, null narrative, ticket) + merge

- **Estimated final test count:**
  - Branches A/B: 941 (938 + 3 integration tests)
  - Branch C: 938 (unchanged)

- **The audit script is preserved as committed code.** Future sprints may rerun it (e.g., to verify Sprint 16 architectural changes) without rewriting. Reusable infrastructure.
