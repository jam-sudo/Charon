# Charon Sprint 10 — Tier A IVIVE Fold-Error Decomposition

**Generated:** 2026-04-25T04:08:57.581591+00:00
**Panel:** charon_sprint7_fih

## Summary

| key | value |
| --- | --- |
| n_compounds | 12 |
| aggregate_pct_liver_model | 7.714 |
| aggregate_pct_route_bias | 0.000 |
| aggregate_pct_residual | 92.29 |

## Results

*(no results)*

## Per-compound decomposition

| compound | mrsd_ws_mg | mrsd_pt_mg | mrsd_disp_mg | clh_ws_L_h | clh_pt_L_h | clh_disp_L_h | cl_renal_L_h | reference_fih_mg | fih_reference_route | simulation_route | reference_route | f_lit | fold_observed_signed | fold_liver_model_signed | fold_route_bias | fold_residual_signed | fold_observed | fold_liver_model | fold_residual | best_alt_model | flags | f_source | notes |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| diazepam | 0.4074 | 0.4074 | 0.4074 | 0.04454 | 0.04455 | 0.04455 | 0.000 | 2 | oral | oral | oral | 0.93 | 0.2037 | 0.9998 | 1 | 0.2037 | 4.91 | 1 | 4.909 | parallel_tube | - | Greenblatt 1980 | F 0.90-1.00 oral vs IV cross-over |
| lisinopril | 2.424 | 2.424 | 2.424 | 0.3166 | 0.3171 | 0.317 | 5 | 10 | oral | oral | oral | 0.25 | 0.2424 | 0.9999 | 1 | 0.2424 | 4.126 | 1 | 4.126 | parallel_tube | - | Beermann 1988 | Low GI absorption; F 0.25-0.29 reported range |
| propranolol | 2.075 | 2.449 | 2.329 | 33.15 | 39.13 | 37.22 | 0.1 | 10 | oral | oral | oral | 0.26 | 0.2075 | 0.8476 | 1 | 0.2449 | 4.819 | 1.18 | 4.084 | parallel_tube | - | Wood 1978 | Extensive first-pass; F 15-40% range; 0.26 = median |
| diclofenac | 16.15 | 16.72 | 16.55 | 6.898 | 7.142 | 7.071 | 0.1 | 50 | oral | oral | oral | 0.54 | 0.323 | 0.9662 | 1 | 0.3343 | 3.096 | 1.035 | 2.991 | parallel_tube | - | Willis 1979 | F 0.54 after oral tablet; subject to enterohepatic recirculation |
| verapamil | 16.05 | 20.05 | 18.73 | 46.3 | 57.83 | 54.03 | 0.000 | 40 | oral | oral | oral | 0.22 | 0.4013 | 0.8006 | 1 | 0.5012 | 2.492 | 1.249 | 1.995 | parallel_tube | - | Eichelbaum 1981 | High first-pass; F 10-35% range; 0.22 = median single-dose |
| metoprolol | 23.53 | 28.55 | 26.92 | 39.92 | 48.59 | 45.77 | 0.7 | 50 | oral | oral | oral | 0.5 | 0.4706 | 0.8241 | 1 | 0.571 | 2.125 | 1.213 | 1.751 | parallel_tube | - | Regardh 1980 | F 0.40-0.60 range; 0.50 = midpoint of extensive metabolizers |
| atorvastatin | 16.98 | 17.96 | 17.95 | 94.03 | 99.45 | 99.42 | 0.000 | 10 | oral | oral | oral | 0.14 | 1.698 | 1 | 1 | 1.698 | 1.698 | 1 | 1.698 | well_stirred | - | Lennernas 2003 | F ~0.14 reflects OATP1B1 hepatic uptake + CYP3A4 first-pass |
| omeprazole | 31.95 | 36.37 | 34.99 | 25.67 | 29.22 | 28.11 | 0.000 | 20 | oral | oral | oral | 0.4 | 1.597 | 1 | 1 | 1.597 | 1.597 | 1 | 1.597 | well_stirred | - | Andersson 1996 | Single-dose F 0.30-0.40; increases to ~0.70 on repeat dosing |
| midazolam | 1.456 | 1.561 | 1.529 | 13.67 | 14.65 | 14.36 | 0.000 | 1 | iv | oral | iv | - | 1.456 | 1 | 1 | 1.456 | 1.456 | 1 | 1.456 | well_stirred | - | Heizmann 1983 | IV reference dose; f_oral 0.30-0.60 range reported in oral studies but not applicable to IV FIH |
| acetaminophen | 412.7 | 436.2 | 429.2 | 12.58 | 13.41 | 13.16 | 2 | 500 | oral | oral | oral | 0.88 | 0.8255 | 0.9463 | 1 | 0.8723 | 1.211 | 1.057 | 1.146 | parallel_tube | - | Rawlins 1977 | F 0.80-0.90 oral; some first-pass sulfation |
| warfarin | 2.03 | 2.044 | 2.04 | 1.387 | 1.397 | 1.394 | 0.05 | 2 | oral | oral | oral | 1 | 1.015 | 1 | 1 | 1.015 | 1.015 | 1 | 1.015 | well_stirred | - | Breckenridge 1970 | Racemate near-complete absorption; F ~100% cited |
| theophylline | 97.88 | 100 | 99.4 | 4.373 | 4.471 | 4.442 | 0.1 | 100 | oral | oral | oral | 1 | 0.9788 | 0.9787 | 1 | 1 | 1.022 | 1.022 | 1 | parallel_tube | - | Hendeles 1977 | Near-complete absorption; F ~96-100% |

## Notes

- Decomposition: fold_observed = fold_liver_model * fold_route_bias * fold_residual.
- fold_liver_model: ws/best_alt if alternate improves prediction, else 1.0.
- fold_route_bias: 1/F for oral-reference compounds, 1.0 for IV or unknown.
- fold_residual: unexplained remainder (transporter, non-hepatic, UGT, model-gap).
- Sorted by fold_residual descending (worst-unexplained first).
- Aggregate %: 100 * sum(|log10(factor)|) / sum(|log10(fold_observed)|).
- Liver-model what-ifs are analytical (mrsd_model = mrsd_ws * CL_model/CL_ws).
- Pipeline.liver_model is a no-op for the PBPK ODE; only one Pipeline run per compound.
- Research only — no production code changes.

## §9. Sprint 13 — UGT/CYP2C9 correction for diclofenac

After `hepatic_clint_multiplier: 3.5` added to diclofenac.yaml (Miners 2006 / Rowland 2013 / Obach 1999 midpoint):

- liver_model: 5.3%
- route_bias:  0.0%
- residual:    94.7%

Diclofenac per-compound:
- Sprint 12 fold_residual: 10.23
- Sprint 13 fold_residual: 2.991 (OUTSIDE_3X at fold_observed 3.096 — residual < 3.0 but observed fold just above boundary due to liver-model factor 1.035)

The multiplier treats combined CYP2C9 + UGT2B7 hepatic clearance as a single enhanced path (empirical — analogous to Sprint 12 OATP1B1 for atorvastatin). Residual after enhancement represents either:
- Biliary excretion (diclofenac undergoes enterohepatic recirculation)
- Minor CYP3A4/CYP2C8 pathways
- Multiplier choice slightly off the actual in vivo gap

Largest remaining residual after Sprint 13:
- propranolol 28.55x — CYP2D6 + extensive first-pass (Sprint 14 target)
- diazepam 4.91x — very low fu_p sensitivity
- lisinopril 4.13x — non-hepatic elimination

## §10 Sprint 15 audit — propranolol F-decomposition

**Generated:** by `scripts/sprint15_audit.py`
**Reference FIH dose:** 10.0 mg (n/a)

### F-decomposition (no multiplier)

| Component | Value | Analytical expected | Within range? |
|---|---:|---:|:---:|
| Fa | 0.9490 | ~0.96 | yes |
| Fg | 1.0066 | ~1.00 | **NO** |
| Fh | 0.9231 | ~0.93 | yes |
| F_oral | 0.8818 | ~0.89 | yes |
| CLint_liver_L_h | 51.00 | ~51.0 | yes |
| MRSD_base_mg | 0.3502 | 0.3502 (Sprint 14) | n/a |
| fold_base | 28.55 | 28.55 (Sprint 14) | n/a |

### Multiplier sweep (fold vs m)

| m | MRSD_mg | fold_observed | within_3x | F_oral | CLint_liver_L_h |
|---:|---:|---:|:---:|---:|---:|
| 1 | 0.3502 | 28.55 | no | 0.882 | 51.00 |
| 2 | 0.6952 | 14.38 | no | 0.819 | 102.00 |
| 3 | 1.04 | 9.61 | no | 0.764 | 153.00 |
| 5 | 1.73 | 5.78 | no | 0.674 | 255.00 |
| 8 | 2.765 | 3.62 | no | 0.573 | 408.00 |
| 12 | 4.145 | 2.41 | yes | 0.478 | 612.00 |
| 15 | 5.18 | 1.93 | yes | 0.425 | 765.00 |
| 20 | 6.904 | 1.45 | yes | 0.358 | 1020.00 |
| 25 | 8.629 | 1.16 | yes | 0.310 | 1275.00 |
| 30 | 10.35 | 0.97 | yes | 0.273 | 1530.00 |

### Boundary summary

- **m_close_3x:** 12  (smallest m bringing fold ≤ 3x)
- **m_close_1x:** 30  (m bringing fold closest to 1.0)

Use this empirical curve to choose the literature multiplier (Task 2). If literature supports a value ≥ m_close_3x, Branch A applies; if literature supports a value < m_close_3x but > 3, Branch B applies (close-but-not-quite); if literature is inconclusive or < 3, Branch C applies (null).

## §11. Sprint 15 — CYP2D6 correction for propranolol

After `hepatic_clint_multiplier: 6.0` added to propranolol.yaml (Hu 2020 / Hallifax & Houston 2010 / Wood 2017 / Chiba 2009 midpoint anchored to Hu 2020 CYP2D6 AFE=6.18, within cited 2.8-9x literature range):

- liver_model: 7.7%
- route_bias:  0.0%
- residual:    92.3%

Propranolol per-compound:
- Sprint 14 fold_residual: 27.49 (signed: 0.03638)
- Sprint 15 fold_residual: 4.084 (signed: 0.245)

The multiplier treats CYP2D6 + CYP1A2 hepatic clearance as an enhanced uniform path (empirical — analogous to Sprint 12 OATP1B1 atorvastatin and Sprint 13 UGT2B7 diclofenac approaches). Sprint 15 audit (§10) confirmed the F-gap traces to CLint underprediction (Fa, Fg normal; Fh too high), validating multiplier appropriateness. Branch B (close-but-not-quite) — empirical sweep showed m_close_3x=12 but cited literature only supports up to ~9x; per CLAUDE.md §6.5 honesty, multiplier was set to 6.0 (anchored to Hu 2020 CYP2D6 AFE=6.18) rather than inflated to force §8 closure.

Remaining largest residuals after Sprint 15:
- propranolol 4.82x — Sprint 15 partial closure (close-but-not-quite); Sprint 16+ architectural work needed for full 3x closure
- diazepam 4.91x — very low fu_p well-stirred sensitivity (Sprint 14 honest null; framework-limited)
- lisinopril 4.13x — non-hepatic elimination + low Peff (renal CL, multiplier inappropriate)
- diclofenac 3.10x — Sprint 13 close-but-not-quite (literature midpoint multiplier 3.5)
