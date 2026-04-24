# Charon Layer 3 FIH Dose Benchmark

**Generated:** 2026-04-24T17:23:58.750666+00:00
**Panel:** charon_sprint7_fih

## Summary

| key | value |
| --- | --- |
| gold_n | 12 |
| gold_within_3x | 8 |
| gold_within_3x_fraction | 0.6667 |
| gold_within_10x | 11 |
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
| diclofenac | oral | 16.15 | pad | 50 | label_start | 3.096 | [FAIL] | [PASS] | FDA Voltaren label, 50 mg PO TID initial |
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
| diclofenac | iv_bolus | 12.52 | pad | 50 | [PASS] | FDA label (Voltaren), 50 mg PO TID; Cp ~1.5 μg/mL, MW 296 |
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

## Sprint 13 comparison (UGT/CYP2C9 enhancement for diclofenac)

Sprint 12 (8/12 = 66.7% within-3x, §8 PASSED): diclofenac at 10.23x was the next-largest remaining residual.
Sprint 13 (diclofenac `hepatic_clint_multiplier: 3.5`): Tier A within-3x = 8/12 = 66.7%.

Per-compound deltas (Sprint 12 → Sprint 13):

| Compound | Sprint 12 fold | Sprint 13 fold | Δ |
|---|---:|---:|:---:|
| midazolam | 1.46 | 1.456 | ~ |
| warfarin | 1.02 | 1.015 | ~ |
| propranolol | 28.55 | 28.55 | ~ |
| verapamil | 2.49 | 2.492 | ~ |
| omeprazole | 1.60 | 1.597 | ~ |
| theophylline | 1.02 | 1.022 | ~ |
| diclofenac | 10.23 | 3.096 | ↓ |
| diazepam | 4.91 | 4.910 | ~ |
| metoprolol | 2.13 | 2.125 | ~ |
| acetaminophen | 1.21 | 1.211 | ~ |
| lisinopril | 4.13 | 4.126 | ~ |
| atorvastatin | 1.70 | 1.698 | ~ |

**Diclofenac-only:**
- Sprint 12: MRSD = 4.89 mg, fold 10.23 (outside 3x)
- Sprint 13: MRSD = 16.15 mg, fold 3.096 (OUTSIDE_3X — just above the 3.0 boundary)

Multiplier applied: 3.5 (Miners 2006 / Rowland 2013 / Obach 1999 midpoint for UGT2B7+CYP2C9 IVIVE gap).
Other 11 compounds show zero delta (multiplier is diclofenac-scoped; all ratios within ±0.3%).

**Interpretation:** Meaningful improvement from 10.23x, but diclofenac's fold sits at ~3.1x — just outside the 3x boundary. §8 target remains PASSED at 8/12. A more aggressive multiplier (4.0) would close the gap but exceeds the literature-supported midpoint (3.5). Honest reporting per spec §6.5.

## Sprint 14 audit — diazepam (null result, 2026-04-24)

Parameter audit of diazepam's stored values against primary literature (Greenblatt 1981 PMID:6790582; Jones & Larsson 2004 PMID:15257067; Obach 1999 DMD 27:1350 Table 2):

| Parameter | Stored | Literature | Verdict |
|---|---|---|---|
| `target_ceff_nM` | 1800 nM (≈513 ng/mL total Cp) | Therapeutic range 200-600 ng/mL (Greenblatt 1980 / Mandelli 1978); Greenblatt 1981 long-term SS mean ~329 ng/mL | WITHIN RANGE |
| `clint_uL_min_mg` | 0.37 | Obach 1999 Table 2 (authoritative) | EXACT MATCH |
| `fu_p` | 0.013 | Greenblatt 1981 pooled mean 0.0148 (range 0.0085-0.0230, n=62) | WITHIN RANGE (1.14x below mean) |
| `bp_ratio` | 0.58 | 0.51-0.59 (Jones & Larsson 2004 PMID:15257067; independent 2022 PBPK model uses exactly 0.58) | EXACT MATCH |

No parameter corrections warranted. No primary-literature source explicitly documents a diazepam-specific in vivo/in vitro CLint ratio ≈ 2x — Branch C (conservative multiplier) does not qualify under the honest-evidence standard.

**Conclusion:** Diazepam's 4.91x residual is IRREDUCIBLE at the current framework. Closure requires architectural work — most likely extended-clearance model handling for extreme-low-fu_p substrates (fu_p=0.013 sits at the well-stirred model's sensitivity boundary) or Vss overprediction investigation (documented in Sprint 3b-1.5 as a known R&R Kp issue for lipophilic neutral drugs).

§8 target remains PASSED at 8/12 = 66.7% (unchanged from Sprint 13).
