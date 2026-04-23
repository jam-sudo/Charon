# Charon Layer 3 FIH Dose Benchmark

**Generated:** 2026-04-23T08:09:12.581482+00:00
**Panel:** charon_sprint7_fih

## Summary

| key | value |
| --- | --- |
| gold_n | 5 |
| gold_within_3x | 3 |
| gold_within_3x_fraction | 0.6 |
| gold_within_10x | 4 |
| sanity_n | 12 |
| sanity_pass_count | 12 |
| sanity_pass_fraction | 1 |
| sanity_failures | - |

## Results

*(no results)*

## Gold (Tier A) — fold-error vs reference FIH

| compound | route | mrsd_pred_mg | limiting_method | reference_fih_mg | source_type | fold_error | within_3x | within_10x | source |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| midazolam | iv_bolus | 0.5784 | pad | 1 | briefing | 1.729 | [PASS] | [PASS] | NDA 20-942 briefing, single-dose IV Phase 1 |
| warfarin | iv_bolus | 1.914 | pad | 2 | label_start | 1.045 | [PASS] | [PASS] | FDA label (Coumadin), initial dose 2-5 mg |
| propranolol | iv_bolus | 0.2755 | pad | 10 | label_start | 36.3 | [FAIL] | [FAIL] | NDA 16-418 (Inderal), approved starting dose |
| verapamil | iv_bolus | 4.914 | pad | 40 | label_start | 8.14 | [FAIL] | [PASS] | NDA 18-817 (Calan), lowest approved PO dose |
| omeprazole | iv_bolus | 13.36 | pad | 20 | label_start | 1.497 | [PASS] | [PASS] | NDA 19-810 (Prilosec), standard 20 mg PO OD |

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
