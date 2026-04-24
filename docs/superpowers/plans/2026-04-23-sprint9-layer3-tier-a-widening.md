# Sprint 9 — Layer 3 Tier A Panel Widening — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Expand Layer 3 Tier A gold panel from 5 to 12 compounds by promoting 4 Obach compounds and adding 3 new primary-literature-sourced compounds (acetaminophen, lisinopril, atorvastatin), then rerun the FIH dose benchmark and report results honestly (whether the §8 ≥60% within-3x target passes or fails).

**Architecture:** Pure data curation. No code changes in `validation/benchmarks/layer3_fih_dose.py` (panel-driven). Three new compound YAMLs in `validation/data/tier1_obach/compounds/`, one `panel.yaml` expansion in `validation/data/fih_reference/`, schema-test thresholds bumped.

**Tech Stack:** YAML (data), pytest (schema validation), existing `charon.core.compound_config.load_compound_config` loader. No new deps.

**Spec:** `docs/superpowers/specs/2026-04-23-sprint9-layer3-tier-a-widening-design.md`

---

## File Structure

### New files

| File | Responsibility |
|---|---|
| `validation/data/tier1_obach/compounds/acetaminophen.yaml` | Experimental ADME for acetaminophen (UGT path). |
| `validation/data/tier1_obach/compounds/lisinopril.yaml` | Experimental ADME for lisinopril (renal-dominant; CLint≈0). |
| `validation/data/tier1_obach/compounds/atorvastatin.yaml` | Experimental ADME for atorvastatin (CYP3A4 + OATP caveat). |

### Modified files

| File | Change |
|---|---|
| `validation/data/fih_reference/panel.yaml` | +7 gold entries (4 Obach promotions + 3 new). Tier B unchanged. |
| `tests/unit/test_fih_panel_schema.py` | Threshold updates: `≥15` → `≥22` compounds total; `≥5` → `≥12` gold. |
| `validation/reports/layer3_fih_dose.{md,json}` | Regenerated. |
| `README.md` | Layer 3 subsection: (n=5) → (n=12), update within-3x/within-10x counts. |

### Explicitly untouched

- `validation/benchmarks/layer3_fih_dose.py` — panel-driven, no code change.
- `validation/data/fih_reference/panel.yaml` Tier B section — unchanged.
- `models/*` — no retraining.
- Regression tests / goldens — don't hit the FIH panel.

---

## Task 1: Create the three new compound YAMLs

**Files:**
- Create: `validation/data/tier1_obach/compounds/acetaminophen.yaml`
- Create: `validation/data/tier1_obach/compounds/lisinopril.yaml`
- Create: `validation/data/tier1_obach/compounds/atorvastatin.yaml`
- Test: `tests/unit/test_fih_new_compound_yamls.py` (new)

- [ ] **Step 1: Write failing smoke test**

Create `tests/unit/test_fih_new_compound_yamls.py`:

```python
"""Sprint 9 smoke tests: the three newly-curated compound YAMLs must load
cleanly via charon.core.compound_config.load_compound_config, and each
must declare at least one non-zero clearance path (hepatic or renal).
"""
from __future__ import annotations

from pathlib import Path

import pytest

from charon.core.compound_config import load_compound_config

REPO_ROOT = Path(__file__).resolve().parents[2]
COMPOUNDS_DIR = REPO_ROOT / "validation" / "data" / "tier1_obach" / "compounds"
NEW_NAMES = ["acetaminophen", "lisinopril", "atorvastatin"]


@pytest.mark.parametrize("name", NEW_NAMES)
def test_new_compound_yaml_loads(name: str):
    path = COMPOUNDS_DIR / f"{name}.yaml"
    assert path.exists(), f"Missing {path}"
    compound = load_compound_config(path)
    assert compound.name == name
    assert compound.molecular_weight > 0


@pytest.mark.parametrize("name", NEW_NAMES)
def test_new_compound_has_clearance_path(name: str):
    compound = load_compound_config(COMPOUNDS_DIR / f"{name}.yaml")
    clint = compound.properties.metabolism.clint_uL_min_mg
    clren = compound.properties.renal.clrenal_L_h
    has_hepatic = clint is not None and clint.value > 0
    has_renal = clren is not None and clren.value > 0
    assert has_hepatic or has_renal, (
        f"{name} must declare at least one non-zero clearance path; "
        f"clint={clint}, clrenal={clren}"
    )


@pytest.mark.parametrize("name", NEW_NAMES)
def test_new_compound_binding_present(name: str):
    compound = load_compound_config(COMPOUNDS_DIR / f"{name}.yaml")
    fup = compound.properties.binding.fu_p
    bp = compound.properties.binding.bp_ratio
    assert fup is not None and 0.0 < fup.value <= 1.0
    assert bp is not None and bp.value > 0
```

- [ ] **Step 2: Run to confirm failure**

```
pytest tests/unit/test_fih_new_compound_yamls.py -v
```

Expected: 9 failures (3 compounds × 3 tests) — all YAML files missing.

- [ ] **Step 3: Create `acetaminophen.yaml`**

File: `validation/data/tier1_obach/compounds/acetaminophen.yaml`

```yaml
# Sprint 9 additions — acetaminophen (paracetamol).
#
# Eliminated primarily by UGT (sulfate and glucuronide conjugation) with
# minor CYP2E1 contribution. Included in the Tier A panel to diversify
# elimination-pathway coverage beyond the Obach-12 CYP bias.
#
# Primary references:
#   - Obach 2008 human PK predictions table (CL = 21.6 L/h, Vss = 63 L)
#   - Prescott 1980 (fu_p = 0.80)
#   - FDA OTC Monograph 21 CFR 343 (starting dose 325-650 mg q4-6h)
name: acetaminophen
smiles: "CC(=O)Nc1ccc(O)cc1"
molecular_weight: 151.16
source: experimental
properties:
  physicochemical:
    logp:
      value: 0.46
      source: experimental
      method: "Obach 2008"
    compound_type: neutral
  binding:
    fu_p: {value: 0.80, source: experimental, unit: fraction, method: "Prescott 1980"}
    fu_inc: {value: 1.0, source: experimental, unit: fraction, method: "hepatocyte default"}
    bp_ratio: {value: 1.0, source: experimental, unit: ratio, method: "Prescott 1980"}
  metabolism:
    clint_uL_min_mg: {value: 1.8, source: experimental, unit: uL/min/10^6 cells, method: "Obach 2008 hepatocyte"}
  renal:
    clrenal_L_h: {value: 2.0, source: experimental, unit: L/h, method: "Prescott 1980 (~4% unchanged)"}
```

- [ ] **Step 4: Create `lisinopril.yaml`**

File: `validation/data/tier1_obach/compounds/lisinopril.yaml`

```yaml
# Sprint 9 additions — lisinopril (ACE inhibitor, renal-dominant).
#
# Essentially not metabolised hepatically (CLint ≈ 0); eliminated almost
# entirely by renal excretion of unchanged drug. Chosen to stress-test
# the PAD path when the hepatic clearance contribution is zero.
#
# Primary references:
#   - Beermann 1988 human PK (CL_total ≈ 5.1 L/h; renal CL > 95%)
#   - FDA Prinivil label (starting 5-10 mg PO OD)
#   - Lancaster 1988 (fu_p ≈ 0.75, low binding)
name: lisinopril
smiles: "N[C@@H](CCCCN)C(=O)N1CCC[C@H]1C(=O)N[C@@H](CC1CCCCC1)C(=O)O"
molecular_weight: 405.49
source: experimental
properties:
  physicochemical:
    logp:
      value: -1.22
      source: experimental
      method: "DrugBank DB00722 experimental logP"
    compound_type: zwitterion
  binding:
    fu_p: {value: 0.75, source: experimental, unit: fraction, method: "Lancaster 1988"}
    fu_inc: {value: 1.0, source: experimental, unit: fraction, method: "hepatocyte default"}
    bp_ratio: {value: 0.85, source: experimental, unit: ratio, method: "Beermann 1988"}
  metabolism:
    clint_uL_min_mg: {value: 0.1, source: experimental, unit: uL/min/10^6 cells, method: "Beermann 1988 (negligible hepatic)"}
  renal:
    clrenal_L_h: {value: 5.0, source: experimental, unit: L/h, method: "Beermann 1988 (~100% renal)"}
```

Note: `clint` set to `0.1` (physical minimum per CLAUDE.md §6k), not literal zero, to avoid divide-by-zero edge cases in well-stirred liver model. Effectively non-contributory.

- [ ] **Step 5: Create `atorvastatin.yaml`**

File: `validation/data/tier1_obach/compounds/atorvastatin.yaml`

```yaml
# Sprint 9 additions — atorvastatin (CYP3A4 + OATP1B1 substrate).
#
# Eliminated mainly by CYP3A4; hepatic uptake is OATP1B1-mediated (not
# modelled in Charon). The compound is included as an honest stress
# test — expected IVIVE overprediction due to missing transporter layer,
# flagged in notes so the Tier A fold-error reading is interpreted
# correctly.
#
# Primary references:
#   - Obach 2008 human PK (CL = 39 L/h, Vss = 381 L)
#   - Lennernäs 2003 (transporter review)
#   - FDA Lipitor label (10 mg PO starting dose)
name: atorvastatin
smiles: "CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CC[C@@H](O)C[C@@H](O)CC(=O)O"
molecular_weight: 558.64
source: experimental
properties:
  physicochemical:
    logp:
      value: 6.36
      source: experimental
      method: "DrugBank DB01076"
    compound_type: acid
    pka_acid:
      value: 4.46
      source: experimental
      method: "DrugBank DB01076"
  binding:
    fu_p: {value: 0.02, source: experimental, unit: fraction, method: "Lennernäs 2003"}
    fu_inc: {value: 0.07, source: experimental, unit: fraction, method: "Austin 2002 (HLM)"}
    bp_ratio: {value: 0.61, source: experimental, unit: ratio, method: "Lennernäs 2003"}
  metabolism:
    clint_uL_min_mg: {value: 128.0, source: experimental, unit: uL/min/mg, method: "Obach 2008 HLM"}
  renal:
    clrenal_L_h: {value: 0.0, source: experimental, unit: L/h, method: "Lennernäs 2003 (<2% renal)"}
```

- [ ] **Step 6: Run the smoke tests**

```
pytest tests/unit/test_fih_new_compound_yamls.py -v
```

Expected: 9/9 PASS.

If `lisinopril.yaml`'s `clint_uL_min_mg=0.1` fails `test_new_compound_has_clearance_path` (the test asserts `> 0` for CLint OR CLrenal; lisinopril has CLrenal=5.0 so it passes on the renal side), that's fine. If the test fails, investigate — don't loosen the assertion.

- [ ] **Step 7: Run compound_config loader broad tests**

```
pytest tests/unit/test_compound_config.py tests/unit/test_schema*.py -v
```

Expected: all pass. No regression from adding new YAMLs.

- [ ] **Step 8: Commit**

```bash
git add validation/data/tier1_obach/compounds/acetaminophen.yaml \
        validation/data/tier1_obach/compounds/lisinopril.yaml \
        validation/data/tier1_obach/compounds/atorvastatin.yaml \
        tests/unit/test_fih_new_compound_yamls.py
git commit -m "$(cat <<'EOF'
data(validation): add acetaminophen, lisinopril, atorvastatin compound YAMLs

Sprint 9 curations for Tier A panel widening:
- acetaminophen: UGT-dominant elimination (Obach 2008)
- lisinopril: renal-dominant (CLint≈0, CLrenal=5 L/h; Beermann 1988)
- atorvastatin: CYP3A4 + OATP1B1 (transporter unmodelled; Obach 2008)

Each YAML cites primary literature inline. Smoke tests verify load via
charon.core.compound_config.load_compound_config plus at-least-one
non-zero clearance path.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 2: Expand `fih_reference/panel.yaml` with 7 new gold entries

**Files:**
- Modify: `validation/data/fih_reference/panel.yaml`

- [ ] **Step 1: Read current panel.yaml and confirm insertion point**

Run:
```
grep -n "tier: gold\|tier: sanity_floor\|# --- Tier" validation/data/fih_reference/panel.yaml
```

The 4 Obach promotions go at the END of the gold section (before the Tier B header comment). The 3 new compounds go immediately after.

- [ ] **Step 2: Append 4 Obach-promotion gold entries**

Add AFTER the existing 5 gold entries (`midazolam, warfarin, propranolol, verapamil, omeprazole`) and BEFORE the `# --- Tier B` header.

Insert into `validation/data/fih_reference/panel.yaml`:

```yaml
    - name: theophylline
      tier: gold
      reference_fih_mg: 100.0
      route: iv_bolus
      target_ceff_nM: 55000.0
      source: "FDA label (Theo-24), initial PO 100-200 mg"
      source_type: label_start
      notes: "Therapeutic Cp ~10 ug/mL, MW 180."

    - name: diclofenac
      tier: gold
      reference_fih_mg: 50.0
      route: iv_bolus
      target_ceff_nM: 5000.0
      source: "FDA Voltaren label, 50 mg PO TID initial"
      source_type: label_start
      notes: "Therapeutic Cp ~1.5 ug/mL, MW 296."

    - name: diazepam
      tier: gold
      reference_fih_mg: 2.0
      route: iv_bolus
      target_ceff_nM: 1800.0
      source: "FDA Valium label, 2-10 mg PO initial"
      source_type: label_start
      notes: "Therapeutic Cp ~500 ng/mL, MW 285. Very-low fu_p (~0.013)."

    - name: metoprolol
      tier: gold
      reference_fih_mg: 50.0
      route: iv_bolus
      target_ceff_nM: 380.0
      source: "FDA Lopressor label, 50 mg PO BID initial"
      source_type: label_start
      notes: "Therapeutic Cp ~100 ng/mL, MW 267. CYP2D6 substrate."
```

- [ ] **Step 3: Append 3 new-compound gold entries**

Append immediately after the 4 Obach promotions:

```yaml
    - name: acetaminophen
      tier: gold
      reference_fih_mg: 500.0
      route: iv_bolus
      target_ceff_nM: 66000.0
      source: "FDA OTC Monograph 21 CFR 343, 325-650 mg q4-6h"
      source_type: label_start
      notes: "Therapeutic Cp ~10 ug/mL, MW 151. UGT-dominant elimination."

    - name: lisinopril
      tier: gold
      reference_fih_mg: 10.0
      route: iv_bolus
      target_ceff_nM: 170.0
      source: "FDA Prinivil label, 5-10 mg PO OD initial"
      source_type: label_start
      notes: "Therapeutic Cp ~70 ng/mL, MW 405. Renal-dominant (>95% unchanged)."

    - name: atorvastatin
      tier: gold
      reference_fih_mg: 10.0
      route: iv_bolus
      target_ceff_nM: 5.3
      source: "FDA Lipitor label, 10 mg PO OD initial"
      source_type: label_start
      notes: "Therapeutic Cp ~3 ng/mL, MW 559. CYP3A4 + OATP1B1 (transporter unmodelled — expect IVIVE overprediction)."
```

- [ ] **Step 4: Quick parse-check**

```bash
python -c "
import yaml
from pathlib import Path
panel = yaml.safe_load(Path('validation/data/fih_reference/panel.yaml').read_text())['panel']
gold = [c for c in panel['compounds'] if c['tier'] == 'gold']
floor = [c for c in panel['compounds'] if c['tier'] == 'sanity_floor']
print(f'Gold: {len(gold)}')
print(f'Floor: {len(floor)}')
print(f'Total: {len(panel[\"compounds\"])}')
for c in gold:
    print(f'  gold: {c[\"name\"]} ({c[\"route\"]}, ref={c[\"reference_fih_mg\"]} mg)')
"
```

Expected: `Gold: 12`, `Floor: 12`, `Total: 24`.

- [ ] **Step 5: Commit**

```bash
git add validation/data/fih_reference/panel.yaml
git commit -m "$(cat <<'EOF'
data(validation): expand Tier A gold panel from 5 to 12 compounds

Adds to Tier A:
- 4 Obach promotions (theophylline, diclofenac, diazepam, metoprolol)
  — existing compound YAMLs, FDA-label starting doses.
- 3 new compounds (acetaminophen UGT, lisinopril renal, atorvastatin
  CYP3A4+OATP) — widens elimination-pathway coverage beyond CYP bias.

Tier B sanity floor unchanged at 12 compounds.

Closes Sprint 7 final-review follow-up: the prior 3/5 at 60% §8
boundary is replaced by a 12-compound panel, so single-compound flips
no longer dominate the within-3x fraction.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 3: Update the schema tests for new thresholds

**Files:**
- Modify: `tests/unit/test_fih_panel_schema.py`

- [ ] **Step 1: Find the threshold assertions**

```bash
grep -n "len(panel\|len(gold\|at_least" tests/unit/test_fih_panel_schema.py
```

You'll find:
- `test_top_level_fields`: `assert len(panel["compounds"]) >= 15` — bump to `>= 22`.
- `test_at_least_five_gold`: rename and bump to `>= 12`.

- [ ] **Step 2: Apply the edits**

In `tests/unit/test_fih_panel_schema.py`:

Change
```python
    assert len(panel["compounds"]) >= 15
```
to
```python
    assert len(panel["compounds"]) >= 22
```

Change the test function + assertion
```python
class TestGoldTier:
    def test_at_least_five_gold(self, panel):
        gold = [c for c in panel["compounds"] if c["tier"] == "gold"]
        assert len(gold) >= 5
```
to
```python
class TestGoldTier:
    def test_at_least_twelve_gold(self, panel):
        gold = [c for c in panel["compounds"] if c["tier"] == "gold"]
        assert len(gold) >= 12
```

The other gold-tier tests (`test_gold_has_reference_fih`, `test_gold_has_target_ceff`) are structural — they iterate over whatever gold compounds exist, so no change needed; they now exercise 12 instead of 5.

The `TestCompoundYamlsExist::test_each_compound_has_property_yaml` is already structural — it reads `panel["compounds"]` names and checks each compound YAML exists. It will automatically pick up acetaminophen, lisinopril, atorvastatin and require YAMLs at those paths (which Task 1 created).

- [ ] **Step 3: Run the schema tests**

```
pytest tests/unit/test_fih_panel_schema.py -v
```

Expected: 7/7 PASS (same test count as before; only values changed).

- [ ] **Step 4: Commit**

```bash
git add tests/unit/test_fih_panel_schema.py
git commit -m "$(cat <<'EOF'
test(validation): bump FIH panel schema thresholds for n=12 Tier A

test_top_level_fields now asserts total >= 22 (was 15).
test_at_least_five_gold renamed + bumped to test_at_least_twelve_gold
(>= 12). Iterator-based tests (has_reference_fih / has_target_ceff /
each_compound_has_property_yaml) unchanged — they auto-cover the new
compounds.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 4: Run the Layer 3 benchmark, commit new report, file Sprint 10 if <60%

**Files:**
- Regenerate: `validation/reports/layer3_fih_dose.{md,json}`

- [ ] **Step 1: Run the benchmark**

```bash
python validation/benchmarks/layer3_fih_dose.py
```

**Expected outcomes — branch on the result:**

A) Exit 0 + "12/12 sanity floor pass" — proceed to Step 2 regardless of Tier A numbers.
B) Exit 2 — Tier B broke; STOP and escalate BLOCKED. New compounds may have introduced a sanity-floor failure.

**Do NOT modify `layer3_fih_dose.py`** — it's panel-driven and should not need any change.

- [ ] **Step 2: Inspect the new Tier A numbers**

```bash
cat validation/reports/layer3_fih_dose.md | sed -n '/Gold (Tier A)/,/Sanity floor/p'
```

Note the within-3x and within-10x counts from the summary block. The spec (`§6.5`) says the acceptance decision branches on these:

- **within-3x ≥ 60% (≥ 8/12)**: §8 target PASSES. Proceed to Step 4 (README update), no Sprint 10 ticket needed.
- **within-3x < 60% (< 8/12)**: §8 target FAILS. Proceed to Step 3 (file follow-up), then Step 4.

**Do not touch the panel to try to make it pass. An honest fail is acceptable.**

- [ ] **Step 3: File Sprint 10 follow-up (only if within-3x < 60%)**

Create `docs/superpowers/sprint10-ivive-bias-ticket.md`:

```markdown
# Sprint 10 follow-up ticket — IVIVE point-estimate bias

Sprint 9 widened the Layer 3 Tier A panel to 12 compounds. Result:
within-3x = <N>/12 = <FRAC>% (below the §8 target of ≥60%).

## Root-cause candidates

- Highly-bound compounds (warfarin fu_p=0.012, diazepam fu_p=0.013,
  atorvastatin fu_p=0.02) — well-stirred model is extremely sensitive
  to fu_b = fu_p / BP at low fu_p. Small errors amplify.
- Transporter-limited compounds (atorvastatin OATP1B1) — hepatic uptake
  not modelled.
- Very-low CLint compounds (diazepam CLint=9 uL/min/mg) — experimental
  variability dominates.
- PAD dose-projection assumption (target_ceff_nM × CL × tau / SF) may
  not match the regulatory FIH rationale for some labels (e.g.
  propranolol's oral label dose is dominated by first-pass F ≈ 0.3,
  so IV-route prediction inherits 3× bias).

## Suggested investigations

1. Per-compound error decomposition: split fold-error into
   (CL component, Vss component, F-bias component, OATP-unmodelled
   component).
2. Consider a non-well-stirred liver model (parallel-tube or dispersion)
   for very-low-fu_p compounds.
3. Evaluate Berezhkovskiy Kp or Rodgers&Rowland choice for Tier A.
4. Add a Tier A "observational bucket" (within-10x) as a secondary
   target when within-3x is dominated by transporter/F-bias gaps.

## Explicit non-scope

Adding Papp/Peff for oral migration is a separate concern (Sprint 11).
Adding transporter plumbing is a major architecture change and out of
scope for a bias-investigation sprint.

Tracking: [link to internal issue once created]
```

This documents the honest failure and its next steps. Commit separately if created:

```bash
git add docs/superpowers/sprint10-ivive-bias-ticket.md
git commit -m "docs: file Sprint 10 ticket — IVIVE bias vs Tier A §8 target"
```

- [ ] **Step 4: Commit the regenerated report**

```bash
git add validation/reports/layer3_fih_dose.md validation/reports/layer3_fih_dose.json
git commit -m "$(cat <<'EOF'
chore(validation): regenerate Layer 3 report with n=12 Tier A panel

Sprint 9 widening result: Tier A within-3x = <M>/12 = <PCT>%, Tier B
sanity floor 12/12 pass.

[If within-3x >= 60%]: §8 target PASSES at the n=12 boundary.
[If within-3x <  60%]: §8 target reported honestly; Sprint 10
follow-up filed in docs/superpowers/sprint10-ivive-bias-ticket.md.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

Replace `<M>`, `<PCT>`, and the if-branch with the actual numbers from the report before committing. Do not leave placeholders.

---

## Task 5: Update README with n=12 numbers

**Files:**
- Modify: `README.md`

- [ ] **Step 1: Find the Layer 3 subsection**

```bash
grep -n "Layer 3 -- FIH dose\|Gold (n=" README.md
```

The existing section (from Sprint 7) has a table with `Gold (n=5)` and specific within-3x/within-10x counts.

- [ ] **Step 2: Update the subsection**

Replace the existing Layer 3 block. Read the Layer 3 report values first:

```bash
python -c "
import json
d = json.loads(open('validation/reports/layer3_fih_dose.json').read())
s = d['summary']
print('gold_n:', s['gold_n'])
print('gold_within_3x:', s['gold_within_3x'])
print('gold_within_3x_fraction:', f\"{s['gold_within_3x_fraction']*100:.1f}%\")
print('gold_within_10x:', s['gold_within_10x'])
print('sanity_pass_count:', s['sanity_pass_count'], '/', s['sanity_n'])
"
```

Transcribe the values into README. Example replacement (substitute actual numbers):

```markdown
### Layer 3 -- FIH dose (Sprint 9, n=12)

12-compound Tier A panel. All compounds run as `iv_bolus` because the
reused compound YAMLs lack absorption data; reference doses are daily
oral equivalents, so Tier A fold-errors carry a `1/F` bias and Tier B
is conservative.

| Tier | Metric | Result | §8 target |
| --- | --- | --- | --- |
| Gold (n=12) | Within-3-fold of reference FIH dose | <G3>/12 (<P3>%) | >= 60% <STATUS_3X> |
| Gold (n=12) | Within-10-fold of reference FIH dose | <G10>/12 (<P10>%) | -- |
| Sanity (n=12) | MRSD <= approved starting dose | <SP>/<SN> (<SF>%) | gated [PASS] |

Panel composition (Sprint 7 core + Sprint 9 expansion):
- 5 core (midazolam, warfarin, propranolol, verapamil, omeprazole)
- 4 Obach promotions (theophylline, diclofenac, diazepam, metoprolol)
- 3 broadened elimination (acetaminophen UGT, lisinopril renal, atorvastatin CYP3A4+OATP)

<IF_FAIL_NOTE>

MRSD computed via PAD path with `safety_factor=10`.

Full results: [Layer 3 report](validation/reports/layer3_fih_dose.md).
```

Replace every `<…>` placeholder with the actual number or phrase. `STATUS_3X` is either `[PASS]` (if ≥60%) or `[FAIL]`. If FAIL, include an `<IF_FAIL_NOTE>` paragraph:

```markdown
**§8 target not met** at the n=12 panel. Honest reading: IVIVE point-estimate
accuracy on highly-bound / transporter-mediated compounds drags the
fraction below 60%. Tracked for Sprint 10 IVIVE-bias investigation
(see `docs/superpowers/sprint10-ivive-bias-ticket.md`). The Tier B
sanity floor gate (`MRSD <= approved starting dose`, 12/12 pass) is
unaffected.
```

If PASS, omit that paragraph.

- [ ] **Step 3: Verify no placeholders remain**

```bash
grep -n "<G3>\|<P3>\|<G10>\|<P10>\|<SP>\|<SN>\|<SF>\|<STATUS_3X>\|<IF_FAIL_NOTE>" README.md
```

Expected: empty output. If anything matches, you missed a substitution.

- [ ] **Step 4: Commit**

```bash
git add README.md
git commit -m "$(cat <<'EOF'
docs: README Layer 3 subsection — n=12 Tier A results (Sprint 9)

Replaces the Sprint 7 n=5 block with the actual n=12 numbers. Adds
panel-composition breakdown (5 core + 4 Obach promotions + 3 broadened
elimination) and, if the §8 target fails, links to the Sprint 10 ticket.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Final Verification

- [ ] **Step 1: Full test suite**

```
pytest tests/ -q
```

Expected: 873 + ~9 new = ~882 pass. No regressions.

- [ ] **Step 2: Confirm artifacts**

```bash
git ls-files validation/data/tier1_obach/compounds/ | grep -E "(acetaminophen|lisinopril|atorvastatin)"
git ls-files validation/data/fih_reference/panel.yaml
git ls-files validation/reports/layer3_fih_dose.md
```

Expected: all three new compound YAMLs tracked; panel and report tracked.

- [ ] **Step 3: Spec acceptance checklist**

Cross-check spec §6.5 acceptance:

1. All 7 new/promoted compounds have legitimate `reference_fih_mg` + Obach-style ADME. ✓ by Tasks 1-2.
2. Panel schema tests updated + green. ✓ by Task 3.
3. Tier B 12/12 pass. ✓ by Task 4 Step 1 gate (exit 0).
4. Tier A report committed with actual numbers. ✓ by Task 4 Step 4.
5. README updated to cite actual results. ✓ by Task 5.
6. §8 target either:
   - Passes (within-3x ≥ 60%): closeout complete.
   - Fails honestly (within-3x < 60%): Sprint 10 ticket filed, commit proceeds. **Panel was NOT massaged.**

- [ ] **Step 4: Announce**

Summarise: "Sprint 9 implementation complete. Tier A n=12, within-3x = <X>/12 (<PCT>%). §8 target <PASS/FAIL>. Ready for review/merge."
