# Charon Layer 3 FIH Dose Benchmark

**Generated:** 2026-04-24T15:40:51.303208+00:00
**Panel:** charon_sprint7_fih

## Summary

| key | value |
| --- | --- |
| gold_n | 12 |
| gold_within_3x | 8 |
| gold_within_3x_fraction | 0.6667 |
| gold_within_10x | 10 |
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
| propranolol | oral | 0.3502 | pad | 10 | label_start | 28.55 | [FAIL] | [FAIL] | NDA 16-418 (Inderal), approved starting dose |
| verapamil | oral | 16.05 | pad | 40 | label_start | 2.492 | [PASS] | [PASS] | NDA 18-817 (Calan), lowest approved PO dose |
| omeprazole | oral | 31.95 | pad | 20 | label_start | 1.597 | [PASS] | [PASS] | NDA 19-810 (Prilosec), standard 20 mg PO OD |
| theophylline | oral | 97.88 | pad | 100 | label_start | 1.022 | [PASS] | [PASS] | FDA label (Theo-24), initial PO 100-200 mg |
| diclofenac | oral | 4.89 | pad | 50 | label_start | 10.23 | [FAIL] | [FAIL] | FDA Voltaren label, 50 mg PO TID initial |
| diazepam | oral | 0.4074 | pad | 2 | label_start | 4.91 | [FAIL] | [PASS] | FDA Valium label, 2-10 mg PO initial |
| metoprolol | oral | 23.53 | pad | 50 | label_start | 2.125 | [PASS] | [PASS] | FDA Lopressor label, 50 mg PO BID initial |
| acetaminophen | oral | 412.7 | pad | 500 | label_start | 1.211 | [PASS] | [PASS] | FDA OTC Monograph 21 CFR 343, 325-650 mg q4-6h |
| lisinopril | oral | 2.424 | pad | 10 | label_start | 4.126 | [FAIL] | [PASS] | FDA Prinivil label, 5-10 mg PO OD initial |
| atorvastatin | oral | 16.98 | pad | 10 | label_start | 1.698 | [PASS] | [PASS] | FDA Lipitor label, 10 mg PO OD initial |

## Sanity floor (Tier B) — MRSD <= approved starting dose

| compound | route | mrsd_pred_mg | limiting_method | approved_starting_dose_mg | pass_floor | source |
| --- | --- | --- | --- | --- | --- | --- |
| theophylline | iv_bolus | 85.21 | pad | 100 | [PASS] | FDA label, initial oral 100-200 mg; therapeutic Cp ~10 μg/mL, MW 180 |
| antipyrine | iv_bolus | 63.57 | pad | 500 | [PASS] | Historical analgesic; commonly 500 mg PO adult; Cp ~9 μg/mL, MW 188 |
| caffeine | iv_bolus | 26.19 | pad | 100 | [PASS] | OTC monograph, 100-200 mg typical; Cp ~5 μg/mL, MW 194 |
| warfarin | iv_bolus | 1.914 | pad | 5 | [PASS] | FDA label, upper bound of starting range |
| diclofenac | iv_bolus | 4.265 | pad | 50 | [PASS] | FDA label (Voltaren), 50 mg PO TID; Cp ~1.5 μg/mL, MW 296 |
| midazolam | iv_bolus | 0.5784 | pad | 2.5 | [PASS] | Conservative anaesthetic pre-med IV |
| propranolol | iv_bolus | 0.2755 | pad | 40 | [PASS] | FDA label, upper initial PO |
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

## Sprint 12 comparison (OATP enhancement for atorvastatin)

Sprint 11 (oral, no OATP): Tier A within-3x = 7/12 = 58.3% (§8 FAILED by ~2%).
Sprint 12 (oral + atorvastatin OATP multiplier=8.0): Tier A within-3x = 8/12 = 66.7% (§8 PASSED).

Per-compound fold-error deltas (Sprint 11 → Sprint 12):

| Compound | Sprint 11 fold | Sprint 12 fold | Δ |
|---|---:|---:|:---:|
| midazolam | 1.46 | 1.46 | ~ |
| warfarin | 1.02 | 1.02 | ~ |
| propranolol | 28.55 | 28.55 | ~ |
| verapamil | 2.49 | 2.49 | ~ |
| omeprazole | 1.60 | 1.60 | ~ |
| theophylline | 1.02 | 1.02 | ~ |
| diclofenac | 10.23 | 10.23 | ~ |
| diazepam | 4.91 | 4.91 | ~ |
| metoprolol | 2.13 | 2.13 | ~ |
| acetaminophen | 1.21 | 1.21 | ~ |
| lisinopril | 4.13 | 4.13 | ~ |
| atorvastatin | 4.64 | 1.70 | ↓↓ |

**Interpretation:** The multiplier is atorvastatin-specific; other 11 compounds should show zero delta within numerical noise. Any non-zero delta on other compounds indicates accidental cross-contamination in the benchmark pipeline.

Confirmed: all 11 non-atorvastatin compounds show identical fold values to Sprint 11 (zero cross-contamination).

Atorvastatin-only:
- Sprint 11 MRSD = 2.15 mg (fold 4.64 outside 3x)
- Sprint 12 MRSD = 16.98 mg (fold 1.70 WITHIN 3x)

Multiplier applied: 8.0 (Izumi 2018 / Barton 2013 OATP1B1 IVIVE midpoint).
