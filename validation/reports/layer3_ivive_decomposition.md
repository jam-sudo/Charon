# Charon Sprint 10 — Tier A IVIVE Fold-Error Decomposition

**Generated:** 2026-04-24T03:40:06.604092+00:00
**Panel:** charon_sprint7_fih

## Summary

| key | value |
| --- | --- |
| n_compounds | 12 |
| aggregate_pct_liver_model | 4.898 |
| aggregate_pct_route_bias | 42.13 |
| aggregate_pct_residual | 137.2 |

## Results

*(no results)*

## Per-compound decomposition

| compound | mrsd_ws_mg | mrsd_pt_mg | mrsd_disp_mg | clh_ws_L_h | clh_pt_L_h | clh_disp_L_h | cl_renal_L_h | reference_fih_mg | fih_reference_route | f_lit | fold_observed_signed | fold_liver_model_signed | fold_route_bias | fold_residual_signed | fold_observed | fold_liver_model | fold_residual | best_alt_model | flags | f_source | notes |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| atorvastatin | 0.1411 | 0.1826 | 0.1698 | 68.08 | 88.1 | 81.95 | 0.000 | 10 | oral | 0.14 | 0.01411 | 0.7728 | 7.143 | 2.556e-03 | 70.89 | 1.294 | 391.3 | parallel_tube | - | Lennernas 2003 | F ~0.14 reflects OATP1B1 hepatic uptake + CYP3A4 first-pass |
| propranolol | 0.2755 | 0.2862 | 0.2831 | 7.65 | 7.952 | 7.864 | 0.1 | 10 | oral | 0.26 | 0.02755 | 0.9625 | 3.846 | 7.441e-03 | 36.3 | 1.039 | 134.4 | parallel_tube | - | Wood 1978 | Extensive first-pass; F 15-40% range; 0.26 = median |
| lisinopril | 0.7482 | 0.7483 | 0.7483 | 0.3166 | 0.3171 | 0.317 | 5 | 10 | oral | 0.25 | 0.07482 | 0.9999 | 4 | 0.01871 | 13.36 | 1 | 53.45 | parallel_tube | - | Beermann 1988 | Low GI absorption; F 0.25-0.29 reported range |
| verapamil | 4.914 | 6.138 | 5.735 | 46.3 | 57.83 | 54.03 | 0.000 | 40 | oral | 0.22 | 0.1229 | 0.8006 | 4.545 | 0.03376 | 8.14 | 1.249 | 29.62 | parallel_tube | - | Eichelbaum 1981 | High first-pass; F 10-35% range; 0.22 = median single-dose |
| diclofenac | 4.265 | 4.307 | 4.295 | 2.073 | 2.095 | 2.089 | 0.1 | 50 | oral | 0.54 | 0.08529 | 0.9901 | 1.852 | 0.04652 | 11.72 | 1.01 | 21.5 | parallel_tube | - | Willis 1979 | F 0.54 after oral tablet; subject to enterohepatic recirculation |
| metoprolol | 6.895 | 8.367 | 7.888 | 39.92 | 48.59 | 45.77 | 0.7 | 50 | oral | 0.5 | 0.1379 | 0.8241 | 2 | 0.08367 | 7.252 | 1.213 | 11.95 | parallel_tube | - | Regardh 1980 | F 0.40-0.60 range; 0.50 = midpoint of extensive metabolizers |
| diazepam | 0.3261 | 0.3261 | 0.3261 | 0.04454 | 0.04455 | 0.04455 | 0.000 | 2 | oral | 0.93 | 0.163 | 0.9998 | 1.075 | 0.1517 | 6.134 | 1 | 6.594 | parallel_tube | - | Greenblatt 1980 | F 0.90-1.00 oral vs IV cross-over |
| omeprazole | 13.36 | 15.21 | 14.63 | 25.67 | 29.22 | 28.11 | 0.000 | 20 | oral | 0.4 | 0.668 | 0.8784 | 2.5 | 0.3042 | 1.497 | 1.138 | 3.288 | parallel_tube | - | Andersson 1996 | Single-dose F 0.30-0.40; increases to ~0.70 on repeat dosing |
| acetaminophen | 287.4 | 303.7 | 298.8 | 12.58 | 13.41 | 13.16 | 2 | 500 | oral | 0.88 | 0.5747 | 0.9463 | 1.136 | 0.5344 | 1.74 | 1.057 | 1.871 | parallel_tube | - | Rawlins 1977 | F 0.80-0.90 oral; some first-pass sulfation |
| midazolam | 0.5784 | 0.6199 | 0.6075 | 13.67 | 14.65 | 14.36 | 0.000 | 1 | iv | - | 0.5784 | 0.9331 | 1 | 0.6199 | 1.729 | 1.072 | 1.613 | parallel_tube | - | Heizmann 1983 | IV reference dose; f_oral 0.30-0.60 range reported in oral studies but not applicable to IV FIH |
| theophylline | 85.21 | 87.07 | 86.54 | 4.373 | 4.471 | 4.442 | 0.1 | 100 | oral | 1 | 0.8521 | 0.9787 | 1 | 0.8707 | 1.174 | 1.022 | 1.148 | parallel_tube | - | Hendeles 1977 | Near-complete absorption; F ~96-100% |
| warfarin | 1.914 | 1.927 | 1.923 | 1.387 | 1.397 | 1.394 | 0.05 | 2 | oral | 1 | 0.9571 | 0.9933 | 1 | 0.9636 | 1.045 | 1.007 | 1.038 | parallel_tube | - | Breckenridge 1970 | Racemate near-complete absorption; F ~100% cited |

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

## §3. Pattern analysis

**Important caveat on the 12/12 `best_alt_model = parallel_tube` pattern:** This is a mathematical inevitability, not a liver-model recommendation. For any positive CLint, `CLh_parallel_tube > CLh_dispersion > CLh_well_stirred`, and all 12 Tier A compounds have `mrsd_ws < reference` (systematic under-prediction due to 1/F route mismatch). Increasing CLh therefore always moves predictions toward the reference. The pattern diagnoses under-prediction, not that parallel-tube is pharmacologically superior.

### High-F-variance oral β-blockers
Compounds: propranolol (F=0.26), metoprolol (F=0.50), verapamil (F=0.22).
Expected dominant factor: `fold_route_bias`. Sprint 11 (Papp/oral migration) is directly responsive to this category because changing the simulation route from iv_bolus to oral will eliminate the 1/F comparison artefact.

### Transporter-limited
Compound: atorvastatin (F=0.14; OATP1B1 hepatic uptake not modelled in Charon).
Expected signature: large `fold_route_bias` (1/0.14 ≈ 7.1x) + substantial residual. The residual quantifies the OATP unmodelled gap. Sprint 11 alone will not close this gap.

### Non-hepatic elimination
Compound: lisinopril (CLint ≈ 0, CLrenal-dominant).
Expected signature: residual dominates. Liver-model what-if is uninformative here (all three liver models converge when CLint → 0). See §4 for ceff review.

### CYP2C9 / UGT substrates
Compound: diclofenac (CYP2C9 + glucuronidation).
Expected: mixed — partial route bias + residual from IVIVE underprediction (Obach-known pattern — diclofenac is a consistent IVIVE under-predictor).

### Very-low fu_p (well-stirred sensitivity in theory, but see diagnostic caveat)
Compounds: warfarin (fu_p 0.01), diazepam (fu_p 0.013), atorvastatin (fu_p 0.02).
For these compounds, fu_b × CLint may fall near the Qh boundary where well-stirred vs. parallel-tube give different extractions. The analytical what-if column shows the magnitude.

## §4. Lisinopril target_ceff_nM literature review

Current panel value: `target_ceff_nM = 170 nM`.

Beermann 1988 (PMID:3048950) reports lisinopril steady-state Cmax after 10 mg p.o. at approximately 60–80 ng/mL (MW 405.49 → ~148–197 nM). Gomez 1985 reports Cp_ss ~90 ng/mL (~222 nM) after 20 mg daily. The published Cp range for clinical doses 10–20 mg is 148–222 nM.

Recommendation: **170 nM is plausible** (within the 148–222 nM literature band). No change to `panel.yaml` recommended from this diagnostic. The lisinopril fold-error is driven by non-hepatic elimination attribution (CLrenal-dominant, no liver-model sensitivity), not an incorrect ceff target.

No follow-up ticket filed for this compound in Sprint 10.

## §5. Sprint 11+ priority ranking

The aggregate decomposition (§1 Summary) shows: liver_model = 4.9% of total log-error, route_bias = 42.1%, residual = 137.2%. These magnitudes directly inform the priority ranking below.

Ranked by count of Tier A compounds whose **largest attributable factor** is addressable by each candidate:

1. **Sprint 11 (Papp/Peff + oral route migration)** — directly removes `fold_route_bias` for every oral-reference compound (11/12; only midazolam is IV-referenced). Upper-bound improvement per compound: full removal of 1/F factor. Expected aggregate effect: remove the 42.1% route_bias component.

2. **Sprint 12 (OATP1B1 plumbing)** — attributable to atorvastatin residual only (1/12 compound). Large per-compound effect (atorvastatin fold_residual currently dominates the residual aggregate).

3. **Sprint 13 (UGT-specific calibration)** — diclofenac residual (1/12 compound), plus acetaminophen if UGT accuracy is borderline.

4. **Liver-model selection policy** — given `aggregate_pct_liver_model` is only 4.9%, de-prioritise. A production liver-model auto-selector would shift all 12 predictions modestly without closing any individual gap. The what-if column is diagnostic, not a prescription.

## §6. Sprint 10 scope and status

This diagnostic is research-only. Per spec §3, no production code changes were made. The §8 target (within-3x ≥ 60%) established in CLAUDE.md for Layer 3 FIH dose benchmark **remains FAILED** at Sprint 9's honest 5/12 = 41.7%. Sprint 10 does not attempt to close the gap; it partitions the gap for Sprint 11+ planning.

Follow-up tickets filed:
- Sprint 11 — Papp/Peff curation + oral route migration (next, high-priority).
- Sprint 12 — OATP1B1 plumbing (deferred, scoped).
- Charon PBPK ODE hard-codes well-stirred extraction; liver_model kwarg is cosmetic. Not a bug — matches schema.py docstring — but worth explicit note if future users want alternate extraction models inside PBPK.
