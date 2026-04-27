# Charon Layer 3 FIH Dose Benchmark

**Generated:** 2026-04-27T17:52:23.934018+00:00
**Panel:** charon_sprint7_fih

## Summary

| key | value |
| --- | --- |
| gold_n | 12 |
| gold_within_3x | 8 |
| gold_within_3x_fraction | 0.6667 |
| gold_within_10x | 12 |
| sanity_n | 12 |
| sanity_pass_count | 12 |
| sanity_pass_fraction | 1 |
| sanity_failures | - |

## Results

*(no results)*

## Gold (Tier A) — fold-error vs reference FIH

| compound | route | mrsd_pred_mg | limiting_method | reference_fih_mg | source_type | fold_error | within_3x | within_10x | source |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| midazolam | oral | 1.456 | pad | 1 | briefing | 1.456 | [PASS] | [PASS] | NDA 20-942 briefing, single-dose IV Phase 1 |
| warfarin | oral | 2.03 | pad | 2 | label_start | 1.015 | [PASS] | [PASS] | FDA label (Coumadin), initial dose 2-5 mg |
| propranolol | oral | 2.075 | pad | 10 | label_start | 4.819 | [FAIL] | [PASS] | NDA 16-418 (Inderal), approved starting dose |
| verapamil | oral | 16.05 | pad | 40 | label_start | 2.492 | [PASS] | [PASS] | NDA 18-817 (Calan), lowest approved PO dose |
| omeprazole | oral | 31.95 | pad | 20 | label_start | 1.597 | [PASS] | [PASS] | NDA 19-810 (Prilosec), standard 20 mg PO OD |
| theophylline | oral | 97.88 | pad | 100 | label_start | 1.022 | [PASS] | [PASS] | FDA label (Theo-24), initial PO 100-200 mg |
| diclofenac | oral | 16.15 | pad | 50 | label_start | 3.096 | [FAIL] | [PASS] | FDA Voltaren label, 50 mg PO TID initial |
| diazepam | oral | 0.4074 | pad | 2 | label_start | 4.91 | [FAIL] | [PASS] | FDA Valium label, 2-10 mg PO initial |
| metoprolol | oral | 23.53 | pad | 50 | label_start | 2.125 | [PASS] | [PASS] | FDA Lopressor label, 50 mg PO BID initial |
| acetaminophen | oral | 412.7 | pad | 500 | label_start | 1.211 | [PASS] | [PASS] | FDA OTC Monograph 21 CFR 343, 325-650 mg q4-6h |
| lisinopril | oral | 3.241 | pad | 10 | label_start | 3.085 | [FAIL] | [PASS] | FDA Prinivil label, 5-10 mg PO OD initial |
| atorvastatin | oral | 16.98 | pad | 10 | label_start | 1.698 | [PASS] | [PASS] | FDA Lipitor label, 10 mg PO OD initial |

## Sanity floor (Tier B) — MRSD <= approved starting dose

| compound | route | mrsd_pred_mg | limiting_method | approved_starting_dose_mg | pass_floor | source |
| --- | --- | --- | --- | --- | --- | --- |
| theophylline | iv_bolus | 85.21 | pad | 100 | [PASS] | FDA label, initial oral 100-200 mg; therapeutic Cp ~10 μg/mL, MW 180 |
| antipyrine | iv_bolus | 63.57 | pad | 500 | [PASS] | Historical analgesic; commonly 500 mg PO adult; Cp ~9 μg/mL, MW 188 |
| caffeine | iv_bolus | 26.19 | pad | 100 | [PASS] | OTC monograph, 100-200 mg typical; Cp ~5 μg/mL, MW 194 |
| warfarin | iv_bolus | 1.914 | pad | 5 | [PASS] | FDA label, upper bound of starting range |
| diclofenac | iv_bolus | 12.52 | pad | 50 | [PASS] | FDA label (Voltaren), 50 mg PO TID; Cp ~1.5 μg/mL, MW 296 |
| midazolam | iv_bolus | 0.5784 | pad | 2.5 | [PASS] | Conservative anaesthetic pre-med IV |
| propranolol | iv_bolus | 0.8651 | pad | 40 | [PASS] | FDA label, upper initial PO |
| diazepam | iv_bolus | 0.3261 | pad | 2 | [PASS] | FDA label (Valium), 2-10 mg initial; therapeutic Cp ~500 ng/mL, MW 285 |
| metoprolol | iv_bolus | 6.895 | pad | 50 | [PASS] | FDA label (Lopressor), 50 mg PO BID initial; Cp ~100 ng/mL, MW 267 |
| verapamil | iv_bolus | 4.914 | pad | 80 | [PASS] | FDA label (Calan), 80 mg PO TID standard |
| omeprazole | iv_bolus | 13.36 | pad | 20 | [PASS] | FDA label (Prilosec), 20 mg PO OD standard |
| dextromethorphan | iv_bolus | 0.679 | pad | 30 | [PASS] | OTC antitussive, 30 mg PO q6-8h; Cp ~8 ng/mL, MW 271 |

## Notes

- Tier A (gold): fold-error vs published FIH or approved-start dose. Report-only.
- Tier B (sanity_floor): MRSD_pred <= approved_starting_dose_mg. GATED — any failure blocks release.
- MRSD computed via PAD path (target_ceff_nM) with safety_factor=10.
- Sources: inline per compound in panel.yaml.

<!-- BEGIN_PRESERVED_HISTORY -->

## Sprint 15 comparison (CYP2D6 enhancement for propranolol — partial closure)

Sprint 14 (8/12 = 66.7% within-3x, §8 PASSED): propranolol at 28.55x was the largest remaining residual.
Sprint 15 (propranolol `hepatic_clint_multiplier: 6.0`): Tier A within-3x = 8/12 = 66.7% (unchanged; propranolol close but stays outside 3x).

**Propranolol-only:**
- Sprint 14: MRSD = 0.3502 mg, fold 28.55x (outside 3x)
- Sprint 15: MRSD = 2.075 mg, fold 4.82x (outside 3x by 1.82x)

Multiplier applied: 6.0 (Hu 2020 / Hallifax & Houston 2010 / Wood 2017 / Chiba 2009 — CYP2D6/high-extraction base IVIVE gap; cited range only supports up to ~9x — applying higher multiplier would violate CLAUDE.md §6.5 honesty).

**Honest interpretation:** Substantial improvement (28.55 → 4.82x) but the literature-cited multiplier range is insufficient to fully close to 3x. Further closure requires either (a) finding additional primary literature with a higher cited ratio, or (b) architectural work (per-CYP2D6 ML-recalibration of CLint, or extended-clearance modeling of in vivo CYP2D6 turnover). Both are Sprint 16+ scope.

### Sprint 14 ticket reconciliation

Sprint 14 (diazepam audit) ticket characterized propranolol's residual as *"ACAT oral F computation architectural gap... Not a multiplier fix."* Sprint 15 audit (Task 1, see §10 of decomposition report) explicitly captured `Fa, Fg, Fh, F_oral` from the pipeline output and found:

- Fa ≈ 0.949 — ACAT functioning normally for high-permeability propranolol
- Fg ≈ 1.007 — non-CYP3A4 substrate, no gut metabolism issue (slightly above 1.0 is solver-tolerance noise around 1.0)
- Fh ≈ 0.923 — analytically traceable to CLint underprediction in HLM
- F_oral ≈ 0.882 (vs literature 0.26) — overprediction by ~3.4x

The gap is therefore the **same Sprint 12/13 multiplier-template pattern** (HLM CLint underprediction for IVIVE-known substrate classes), **not an ACAT architectural issue**. Sprint 14's claim was a hypothesis based on a quick mental model and is corrected here per CLAUDE.md §6.5 honesty.

## Sprint 17 Branch B — lisinopril Peff back-calibration (2026-04-27)

Audit-first sprint surfaced that Charon **over-predicts F_oral for lisinopril** (0.33 pred vs 0.25 obs, Beermann 1988). Sprint 16 closure had misframed this as "renal CL refinement"; Sprint 17 audit refuted that — CL_renal=5.0 is direct experimental input.

Root cause: `peff_cm_s: 0.3e-4` in YAML, which the file's existing comment self-flagged as "NO primary Peff measurement located". Lisinopril is BCS III + PEPT1 substrate (Knutter 2008); ACAT does not model PEPT1.

Correction: `peff_cm_s = 2.10e-5` (back-calibrated to Beermann 1988 F_obs=0.25 typical adult). lisinopril row update:
- Sprint 16: MRSD 2.424 mg, fold 4.126x
- Sprint 17: MRSD 3.241 mg, fold 3.085x (Branch B close-but-not-quite)

§8 status: 8/12 = 66.7% unchanged (per CLAUDE.md §6.5, Peff cannot honestly go below Beermann 1988 lower-bound). PEPT1 transporter modeling deferred as architectural sprint candidate.

Updated residual classifications:
- propranolol 4.82x — Sprint 15 Branch B (CYP2D6)
- diazepam 4.91x — Sprint 14 honest null (low fu_p)
- lisinopril 3.085x — Sprint 17 Branch B (Peff calibration; PEPT1 not modeled) ← updated from 4.13x
- diclofenac 3.10x — Sprint 13 Branch B (UGT/CYP2C9)
