# Charon Layer 3 FIH Dose Benchmark

**Generated:** 2026-04-24T06:34:06.234099+00:00
**Panel:** charon_sprint7_fih

## Summary

| key | value |
| --- | --- |
| gold_n | 12 |
| gold_within_3x | 7 |
| gold_within_3x_fraction | 0.5833 |
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
| atorvastatin | oral | 2.154 | pad | 10 | label_start | 4.642 | [FAIL] | [PASS] | FDA Lipitor label, 10 mg PO OD initial |

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

## Sprint 11 comparison (oral migration)

Sprint 9 (iv_bolus) Tier A within-3x: 5/12 = 41.7% (§8 FAILED).
Sprint 11 (oral)     Tier A within-3x: 7/12 = 58.3% (§8 FAIL — target is >=60%).

Per-compound fold-error deltas (iv→oral):

| Compound | Sprint 9 fold | Sprint 11 fold | improvement |
|---|---:|---:|:---:|
| midazolam | 1.73 | 1.46 | ✓ |
| warfarin | 1.04 | 1.02 | ~ |
| propranolol | 36.30 | 28.55 | ✓ |
| verapamil | 8.14 | 2.49 | ✓ |
| omeprazole | 1.50 | 1.60 | ~ |
| theophylline | 1.17 | 1.02 | ✓ |
| diclofenac | 11.70 | 10.23 | ~ |
| diazepam | 6.13 | 4.91 | ✓ |
| metoprolol | 7.25 | 2.13 | ✓ |
| acetaminophen | 1.74 | 1.21 | ✓ |
| lisinopril | 13.36 | 4.13 | ✓ |
| atorvastatin | 70.90 | 4.64 | ✓ |

**Interpretation:**

- Compounds that moved INTO within-3x from outside: verapamil, metoprolol
- Compounds that moved OUT of within-3x: none (Peff curation is sound)
- Residuals for propranolol, diclofenac, diazepam, lisinopril, atorvastatin remain large — Sprint 12/13 target compounds.

Sprint 10's prediction that removing the 42% route_bias aggregate would shift the panel significantly toward within-3x compliance is supported by this result.
