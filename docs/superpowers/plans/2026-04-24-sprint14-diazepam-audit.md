# Sprint 14 — Diazepam Parameter Audit Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Audit diazepam's stored parameter values against primary literature. Correct any that differ by >1.5x. Accept honest null result if parameters verify.

**Architecture:** Data-audit sprint. WebSearch-driven literature verification. 0-2 YAML edits depending on findings. No code changes, no new tests.

**Tech Stack:** WebSearch, pytest (regression only), PyYAML.

**Spec:** `docs/superpowers/specs/2026-04-24-sprint14-diazepam-audit-design.md`

---

## Current diazepam values (from this worktree)

**panel.yaml (tier: gold):**
- `reference_fih_mg: 2.0`
- `route: oral`
- `target_ceff_nM: 1800.0`
- Source note: "Therapeutic Cp ~500 ng/mL, MW 285. Very-low fu_p (~0.013)."

**diazepam.yaml:**
- `fu_p: 0.013` (Obach 1999 Table 2)
- `fu_inc: 0.67` (Obach 1999)
- `bp_ratio: 0.58` (Obach 1999)
- `clint_uL_min_mg: 0.37` (Obach 1999)
- `clrenal_L_h: 0.0` (Obach 1999)
- `peff_cm_s: 3.0e-4` (Lennernäs 2007)

**Sprint 13 outcome:** diazepam MRSD = 0.41 mg, ref 2.0 mg, fold 4.91x (outside 3x).

---

## Task 1: Literature audit via WebSearch

**No file changes in this task — pure research.**

- [ ] **Step 1: Verify `target_ceff_nM = 1800.0`**

Run WebSearch queries:
- `"diazepam therapeutic plasma concentration ng/mL Greenblatt"`
- `"diazepam therapeutic range Mandelli 1978"`
- `"diazepam Cp sedative anxiolytic concentration"`

Goal: confirm therapeutic total Cp range. Expected: 200-500 ng/mL (Greenblatt 1980) or 200-600 ng/mL (Mandelli 1978). Convert to nM: value_ng_mL × 1e-6 / MW (284.74 g/mol) × 1e9 nM.

- 500 ng/mL → 1755 nM
- 600 ng/mL → 2107 nM
- 400 ng/mL → 1405 nM

Current `1800 nM` corresponds to ~513 ng/mL (within therapeutic range). **Expected outcome: value verifies OK.**

- [ ] **Step 2: Verify `clint_uL_min_mg = 0.37`**

WebSearch:
- `"diazepam HLM CLint Obach 1999 uL/min/mg"`
- `"Obach 1999 diazepam intrinsic clearance HLM"`

Expected: Obach 1999 DMD 27:1350 Table 2 reports diazepam HLM CLint. Direct value reported in Obach 1999 or derivable from reported in vivo CLint_u / scaling assumption.

**Critical check:** confirm the value 0.37 is the HLM CLint (not hepatocyte CLint or another assay).

- [ ] **Step 3: Verify `fu_p = 0.013`**

WebSearch:
- `"diazepam plasma protein binding unbound fraction Obach"`
- `"diazepam fu_p Benet 1996"`

Expected: literature 0.013-0.02. Current 0.013 is at the low end but within range.

- [ ] **Step 4: Verify `bp_ratio = 0.58`**

WebSearch:
- `"diazepam blood-to-plasma ratio"`
- `"diazepam Rb/p ratio"`

Expected: 0.60-0.75 typical. Current 0.58 is slightly low but plausible.

- [ ] **Step 5: Summarize audit findings**

For each parameter, record:
- Current value
- Verified literature value (or range)
- Discrepancy factor
- Action: "no change" | "correct to X" | "flag uncertain"

Prepare a summary block to include in the commit message or ticket narrative of subsequent tasks.

- [ ] **Step 6: No commit in this task** (research only).

---

## Task 2: Apply corrections if warranted (or proceed with null result)

**Files (conditional):**
- If `target_ceff_nM` needs correction: `validation/data/fih_reference/panel.yaml`
- If `clint_uL_min_mg`, `fu_p`, or `bp_ratio` need correction: `validation/data/tier1_obach/compounds/diazepam.yaml`

- [ ] **Step 1: Branch based on Task 1 outcome**

**Branch A — All values verify correct (most likely):**
- No file changes.
- Proceed directly to Task 3 (benchmark rerun is unnecessary; diazepam values unchanged so MRSD unchanged).
- Sprint 14 outcome is "audit verified; diazepam residual is framework-limited".

**Branch B — One or more values differ by >1.5x from literature:**
- For each incorrect value, edit the YAML using the Edit tool.
- Value replacement format: preserve YAML structure, update `value:` and `method:` to cite the primary source you verified.
- Example (hypothetical correction):
  ```yaml
  fu_p: {value: 0.017, source: experimental, unit: fraction, method: "Obach 1999 Table 2 (revised from 0.013 per Greenblatt 1980 review citing pooled fu_p=0.015-0.020)"}
  ```

**Branch C — Literature supports a conservative 2.0x multiplier for diazepam IVIVE specifically:**
- Add `hepatic_clint_multiplier: 2.0` to diazepam.yaml `metabolism` block with explicit "conservative — at lower end of Obach 1999 IVIVE ratio range for CYP3A4/2C19 substrates" caveat.
- This is the BORDERLINE option; apply ONLY if Task 1 Step 2 found literature (e.g., Obach 1999 Table 3) explicitly listing diazepam in vivo/in vitro ratio ≈ 2x.

- [ ] **Step 2: Verify edits via schema load (if any edits made)**

```bash
python3 -c "
import yaml
from charon.core.schema import CompoundConfig
data = yaml.safe_load(open('validation/data/tier1_obach/compounds/diazepam.yaml').read())
c = CompoundConfig.model_validate(data)
print('fu_p:', c.properties.binding.fu_p.value if c.properties.binding.fu_p else 'None')
print('clint:', c.properties.metabolism.clint_uL_min_mg.value if c.properties.metabolism.clint_uL_min_mg else 'None')
print('bp:', c.properties.binding.bp_ratio.value if c.properties.binding.bp_ratio else 'None')
mult = c.properties.metabolism.hepatic_clint_multiplier
print('multiplier:', mult.value if mult else 'None')
"
```

- [ ] **Step 3: Run existing YAML/schema tests**

```bash
pytest tests/unit/test_compound_config.py tests/unit/test_fih_new_compound_yamls.py tests/unit/test_fih_panel_schema.py -q
```

Expected: all pass.

- [ ] **Step 4: Commit edits (if any)**

If Branch A (no changes): skip.

If Branch B or C (edits made):
```bash
git add validation/data/tier1_obach/compounds/diazepam.yaml
# If panel also touched:
git add validation/data/fih_reference/panel.yaml
git commit -m "data(sprint14): audit-driven corrections to diazepam <FIELD> (<SOURCE>)"
```

Where `<FIELD>` describes what was corrected (e.g., `fu_p` or `target_ceff_nM`) and `<SOURCE>` is the citation.

---

## Task 3: Regenerate benchmarks + Sprint 14 narrative

**Files:**
- Regenerate (only if Task 2 made data changes): `validation/reports/layer3_fih_dose.{md,json}`, `validation/reports/layer3_ivive_decomposition.{md,json}`
- Modify (always): append Sprint 14 narrative section to report markdown files

- [ ] **Step 1: If Task 2 changed data, regenerate benchmarks**

```bash
python3 validation/benchmarks/layer3_fih_dose.py
python3 validation/benchmarks/layer3_ivive_decomposition.py
```

Expected outputs:
- `[OK] Sanity floor: 12/12 pass.`
- `[OK] Decomposition wrote ... (12 compounds)`

Extract diazepam's new fold:
```bash
python3 -c "
import json
d = json.load(open('validation/reports/layer3_fih_dose.json'))
for r in d['extra_sections']['Gold (Tier A) — fold-error vs reference FIH']:
    if r['compound'] == 'diazepam':
        print(f'diazepam: mrsd={r[\"mrsd_pred_mg\"]:.3g}, ref={r[\"reference_fih_mg\"]:.3g}, fold={r[\"fold_error\"]:.2f}, 3x={r[\"within_3x\"]}')
s = d['summary']
print(f'Tier A within-3x: {s[\"gold_within_3x\"]}/{s[\"gold_n\"]}')
"
```

If Task 2 was Branch A (no changes), skip this step — benchmarks are unchanged.

- [ ] **Step 2: Append Sprint 14 narrative to `validation/reports/layer3_fih_dose.md`**

Three template variants — choose based on Task 2 outcome:

**Branch A (null result):**
```markdown

## Sprint 14 audit — diazepam (null result, 2026-04-24)

Parameter audit of diazepam's `target_ceff_nM` (1800 nM), `clint_uL_min_mg` (0.37), `fu_p` (0.013), `bp_ratio` (0.58) against primary literature:

- target_ceff_nM 1800: corresponds to ~513 ng/mL total Cp, within Greenblatt 1980 / Mandelli 1978 therapeutic range (200-600 ng/mL). VERIFIED.
- clint_uL_min_mg 0.37: matches Obach 1999 DMD 27:1350 Table 2. VERIFIED.
- fu_p 0.013: at low end of literature range 0.013-0.02. VERIFIED.
- bp_ratio 0.58: slightly below typical range 0.60-0.75 but plausible. VERIFIED.

No parameter corrections warranted. Diazepam's 4.91x residual is irreducible at the current framework without architectural changes (extended-clearance model for very low fu_p substrates, or diagnostic investigation of MRSD sensitivity to fu_p near-zero regime).

§8 target remains PASSED at 8/12 = 66.7%.
```

**Branch B (correction applied):**
```markdown

## Sprint 14 audit — diazepam (data correction, 2026-04-24)

Parameter audit found diazepam <FIELD> differed from primary literature; corrected:
- Old: <FIELD>=<OLD_VALUE>
- New: <FIELD>=<NEW_VALUE> (source: <CITATION>)

Benchmark rerun:
- Sprint 13: diazepam fold 4.91
- Sprint 14: diazepam fold <S14_FOLD>  (<WITHIN_3X|OUTSIDE_3X>)

Tier A within-3x: Sprint 13 8/12 → Sprint 14 <M>/12 = <PCT>%.
```

**Branch C (conservative multiplier applied):**
```markdown

## Sprint 14 — diazepam conservative multiplier (2026-04-24)

Audit verified all stored parameters match primary literature. Literature (Obach 1999 Table 3) supports conservative 2.0x IVIVE under-prediction factor for CYP3A4/2C19 substrates like diazepam. Applied `hepatic_clint_multiplier: 2.0` with explicit caveat — at LOWER end of Obach-documented range, not editorial inflation.

- Sprint 13: diazepam MRSD 0.41 mg, fold 4.91x
- Sprint 14: diazepam MRSD <S14_MRSD> mg, fold <S14_FOLD>x (<WITHIN_3X|OUTSIDE_3X>)

Tier A within-3x: 8/12 → <M>/12.

Caveat: diazepam is a reference-drug for IVIVE accuracy (not a systematic under-predictor like UGT/OATP substrates). Multiplier 2.0 is borderline but literature-defensible.
```

Fill all placeholders with actual values.

- [ ] **Step 3: Append to `validation/reports/layer3_ivive_decomposition.md`** (if Task 2 changed data only)

For Branch A: skip (decomposition unchanged).

For Branch B/C: append similar section `## §10. Sprint 14 — diazepam audit` with updated residual for diazepam.

- [ ] **Step 4: Commit**

If Task 2 was Branch A:
```bash
git add validation/reports/layer3_fih_dose.md
git commit -m "docs(sprint14): append null-result audit findings for diazepam to layer3 report"
```

If Task 2 was Branch B/C (data + reports regenerated):
```bash
git add validation/reports/layer3_fih_dose.md \
        validation/reports/layer3_fih_dose.json \
        validation/reports/layer3_ivive_decomposition.md \
        validation/reports/layer3_ivive_decomposition.json
git commit -m "chore(sprint14): regenerate Layer 3 reports after diazepam audit corrections"
```

---

## Task 4: Update Sprint 10 ticket + verify full suite

**Files:**
- Modify: `docs/superpowers/sprint10-ivive-bias-ticket.md`

- [ ] **Step 1: Append Sprint 14 status**

Append to the ticket. Use the template matching Task 2's outcome branch:

**Branch A (null):**
```markdown

## Sprint 14 (diazepam audit, 2026-04-24) — null result

Parameter audit verified diazepam's `target_ceff_nM`, `clint_uL_min_mg`, `fu_p`, `bp_ratio` against primary literature (Greenblatt 1980, Mandelli 1978, Obach 1999). No corrections warranted — all values within literature range.

Diazepam 4.91x residual is irreducible at the current framework. Closure requires architectural work (extended-clearance model for extreme-low-fu_p substrates; Sprint 15+ scope).

**Layer 3 Tier A within-3x:** 8/12 = 66.7% (unchanged from Sprint 13; §8 still PASSED).
```

**Branch B (correction):**
```markdown

## Sprint 14 (diazepam audit, 2026-04-24) — parameter correction

Audit found <FIELD> differed from primary literature. Corrected to <NEW_VALUE> (citation: <SOURCE>).

**Diazepam delta:**
- Sprint 13: MRSD 0.41 mg, fold 4.91x (outside 3x)
- Sprint 14: MRSD <S14_MRSD> mg, fold <S14_FOLD>x (<WITHIN_3X|OUTSIDE_3X>)

**Layer 3 Tier A within-3x:** 8/12 → <M>/12 = <PCT>%.
```

**Branch C (multiplier):**
```markdown

## Sprint 14 (diazepam audit + conservative multiplier, 2026-04-24)

Audit verified stored parameters correct. Applied conservative `hepatic_clint_multiplier: 2.0` with Obach 1999 Table 3 citation (lower end of CYP3A4/2C19 substrate IVIVE ratio range).

**Diazepam delta:**
- Sprint 13: fold 4.91x
- Sprint 14: fold <S14_FOLD>x (<WITHIN_3X|OUTSIDE_3X>)

**Layer 3 Tier A within-3x:** 8/12 → <M>/12 = <PCT>%.

Honest caveat: diazepam is a reference-drug for IVIVE, not a systematic under-predictor. The 2.0x multiplier is borderline but literature-grounded.
```

- [ ] **Step 2: Run full test suite**

```bash
pytest -q
```

Expected: 938 passed (no new tests; no code changes).

- [ ] **Step 3: Commit ticket**

```bash
git add docs/superpowers/sprint10-ivive-bias-ticket.md
git commit -m "docs(sprint14): Sprint 10 ticket — Sprint 14 diazepam audit findings"
```

Check `git status` — only layer2 drift should remain unstaged.

- [ ] **Step 4: Summarize commits**

```bash
git log --oneline main..HEAD
```

Should show 1-4 Sprint 14 commits depending on branch outcome.

---

## Self-Review Notes (for the implementer)

- **Most likely outcome: Branch A (null result).** Diazepam is a well-studied reference compound. Its parameters were populated from Obach 1999 (Table 2 for CLint/fu_p/bp_ratio; standard therapeutic range for target_ceff). Audit likely verifies all correct.

- **Null result is NOT a failure.** Sprint 10 was research-only and declared a valid research outcome. Sprint 14's null result tells us diazepam's residual is structural, not data-entry error.

- **Do NOT apply Branch C (multiplier) unless Step 2 of Task 1 found direct primary-literature support** (e.g., Obach 1999 listing diazepam in vivo/in vitro ~2x specifically). "CYP3A4 substrates generally 1.5-3x" is not specific enough support — diazepam may sit at the 1.5x end.

- **If Branch B applied (data correction), verify other compounds don't share the same wrong value.** E.g., if `bp_ratio` convention is off for diazepam, check whether Obach panel compounds have similar issue. Flag as follow-up but do NOT batch-fix others in this sprint.

- **Commit count estimate:** 1 commit (Branch A) or 2-3 commits (Branches B/C).

- **Estimated final test count:** 938 (unchanged).
