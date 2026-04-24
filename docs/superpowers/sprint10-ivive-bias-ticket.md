# Sprint 10 follow-up ticket — IVIVE point-estimate bias

Sprint 9 widened the Layer 3 Tier A panel to 12 compounds. Result: **within-3x = 5/12 = 41.7%** (below the §8 target of ≥60%). Within-10x = 8/12 = 66.7%.

## Per-compound fold-errors (Sprint 9 result)

| Compound | MRSD_pred (mg) | Ref FIH (mg) | Fold | 3x | 10x |
|---|---:|---:|---:|:---:|:---:|
| midazolam | 0.58 | 1 | 1.73 | ✓ | ✓ |
| warfarin | 1.91 | 2 | 1.04 | ✓ | ✓ |
| propranolol | 0.28 | 10 | 36.3 | | |
| verapamil | 4.91 | 40 | 8.14 | | ✓ |
| omeprazole | 13.36 | 20 | 1.50 | ✓ | ✓ |
| theophylline | 85.2 | 100 | 1.17 | ✓ | ✓ |
| diclofenac | 4.27 | 50 | 11.7 | | |
| diazepam | 0.33 | 2 | 6.13 | | ✓ |
| metoprolol | 6.90 | 50 | 7.25 | | ✓ |
| acetaminophen | 287.4 | 500 | 1.74 | ✓ | ✓ |
| lisinopril | 0.75 | 10 | 13.36 | | |
| atorvastatin | 0.14 | 10 | 70.9 | | |

## Root-cause candidates (from observed pattern)

- **High-F-variance oral β-blockers**: propranolol (fold 36x) and metoprolol (fold 7x) have extensive hepatic first-pass; IV-predicted MRSD inherits 1/F bias when compared to oral reference. Expected per Sprint 7's documented 1/F caveat.
- **Transporter-limited**: atorvastatin (fold 71x). OATP1B1 uptake is not modelled. Confirmed the honest stress test outcome — this is documented in the compound YAML.
- **Non-hepatic elimination gaps**: lisinopril (fold 13x) has CLint≈0 + CLrenal=5 L/h. Either CLrenal scaling in the PAD path is off or the target_ceff_nM choice (170 nM ≈ 70 ng/mL) is low; re-check against Beermann 1988.
- **CYP2C9/UGT compound**: diclofenac (fold 12x). Obach-known IVIVE under-predictor.
- **Very-low fu_p**: diazepam (fold 6x, fu_p=0.013). Well-stirred model sensitivity to fu_b near-zero is a systemic weakness.

## Suggested Sprint 10 investigations

1. **Per-compound error decomposition**: split fold-error into (CL_hepatic, CL_renal, F-bias, transporter-unmodelled, Kp-method).
2. **Alternate liver models for extreme fu_p**: parallel-tube or dispersion vs well-stirred for warfarin/diazepam/atorvastatin (fu_p < 0.03).
3. **Empirical F-correction lookup**: accept oral references by adjusting MRSD with a known-F factor for the 3-4 worst offenders; flag others honestly.
4. **Review target_ceff_nM values**: lisinopril 170 nM vs literature — verify Beermann 1988 Cp.

## Explicit non-scope

Papp/Peff curation + oral route migration (Sprint 11), OATP transporter plumbing (major architecture), Kp method choice beyond R&R (research-grade).

## Tracking

Filed at `docs/superpowers/sprint10-ivive-bias-ticket.md` during Sprint 9 closeout (2026-04-23, commit-TBD).

## Status — Sprint 10 diagnostic merged (2026-04-23)

Diagnostic produced at `validation/reports/layer3_ivive_decomposition.md`.

**Aggregate attribution (n=12 Tier A):**
- Liver-model choice: **4.9%** (analytical what-if; low because production PBPK ODE implicitly fixes well-stirred extraction — see §6 architectural note)
- 1/F route bias (oral-reference vs IV-simulation mismatch): **42.1%**
- Residual (unmodelled transporters, non-hepatic routes, UGT, model gaps): **137.2%** (exceeds 100% because signed factors partially cancel across compounds; atorvastatin alone contributes ~391x fold_residual)

**Largest residuals (worst-unexplained):** atorvastatin (391x, OATP1B1), propranolol (~135x), verapamil, metoprolol, lisinopril (renal).

**Investigation outcomes (vs original 5 candidates):**
1. Per-compound error decomposition: **DONE** (full table in §2 of report).
2. Alternate liver models for extreme fu_p: **DONE** (analytical what-if; aggregate effect modest at 4.9%, larger per-compound for high-extraction β-blockers and verapamil).
3. Empirical F-correction lookup: **not applied** (per spec §3 non-goals — research report only; Sprint 11 will address F-bias via oral route migration, not lookup).
4. Review target_ceff_nM values (lisinopril): **DONE** (§4 confirms 170 nM plausible, no change).
5. Architectural note: **NEW FINDING** — Pipeline.liver_model is a no-op for the PBPK ODE (matches schema.py documented behaviour). Sprint 10 orchestrator works around this via analytical MRSD scaling.

**Remediation tickets:**
- Sprint 11 — Papp/Peff curation + oral route migration (next, high-priority, directly addresses 42% route_bias aggregate).
- Sprint 12 — OATP1B1 plumbing (deferred, attributable to atorvastatin's dominant residual).
- Sprint 13 — UGT / CYP2C9 calibration refresh (diclofenac residual).
- §8 target remains FAILED at Sprint 9's honest 5/12 = 41.7%. Not revisited by Sprint 10.
