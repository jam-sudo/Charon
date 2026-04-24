# Sprint 10 — IVIVE point-estimate bias diagnostic (research, no code in prod path)

**Status:** Design approved (brainstorming), pending user spec review before plan.

## 1. Context

Sprint 9 merged (commit `ac61b72`) widened Layer 3 Tier A to n=12 and produced an honest §8 failure: **within-3x = 5/12 = 41.7%** (target ≥60%), within-10x = 8/12 = 66.7%. Tier B sanity floor remained 12/12. Per CLAUDE.md §6.5 honesty clause the panel was not massaged. The follow-up ticket (`docs/superpowers/sprint10-ivive-bias-ticket.md`) lists per-compound fold-errors and five root-cause categories (high-F-variance oral β-blockers, transporter-limited, non-hepatic elimination, CYP2C9/UGT, very-low fu_p).

Sprint 10 does **not** try to close the §8 gap. It attributes each Tier A fold-error to a small number of named factors so Sprint 11+ scope is defined by data rather than intuition.

## 2. Goal

Produce a diagnostic report that decomposes each of the 12 Tier A fold-errors into three attributable factors, ranks Sprint 11+ work by attributable-error contribution, and leaves the production pipeline untouched.

## 3. Non-goals

- Empirical F-correction injection into the prod path (Sprint 10-ticket item 3).
- OATP / transporter modelling (deferred — post-Sprint 12).
- Papp/Peff curation and oral route migration (Sprint 11).
- Edits to `validation/data/fih_reference/panel.yaml` (lisinopril `target_ceff_nM` change, if warranted, is filed as a separate follow-up).
- Changes to `src/charon/core/liver_models.py` behaviour or production default model selection.

## 4. Architecture — 3-factor log-additive decomposition

For each Tier A compound the observed fold-error (predicted MRSD over reference FIH dose) is split as:

```
log10(fold_observed) = log10(fold_liver_model)
                     + log10(fold_route_bias)
                     + log10(fold_residual)
```

where

- `fold_liver_model` = `well_stirred_MRSD / best_alternate_MRSD`. `best_alternate` is whichever of `parallel_tube` or `dispersion` brings the MRSD closest to the reference (minimum `|log10(pred/ref)|`). If neither moves the prediction closer, this factor is 1.0 (no explanatory power).
- `fold_route_bias` = 1 / F_literature for oral-only references; 1.0 for IV-bolus references or when F is unknown (then flagged). Intuition: the production pipeline predicts IV MRSD; when compared against an oral clinical reference we inherit a 1/F bias that is a route-comparison artefact, not an IVIVE error.
- `fold_residual` = `fold_observed / (fold_liver_model × fold_route_bias)`. This is the unexplained remainder — the candidate for transporter / non-hepatic / UGT / low-fu_p / model-gap attribution.

Invariant: `log10(fold_observed) == log10(fold_liver_model) + log10(fold_route_bias) + log10(fold_residual)` within `rtol=1e-9`. Tested.

### 4.1 Why three factors and not more

Charon has many knobs (fu_p, CLint, Kp, Fg, CLrenal, …). All Tier A compounds use Tier 1 experimental overrides on fu_p and CLint (panel.yaml), so their within-compound contribution is 0 by construction. Fg is computed inside ACAT and Tier A uses `iv_bolus` so Fg is not on the IV path. Renal CLint is from YAML literature. That leaves **liver-model choice** (a Charon-side modelling decision) and **1/F route mismatch** (a comparison artefact) as the two factors Sprint 10 can cleanly quantify; everything else collapses into `residual`. A fourth "Kp-method" factor was considered but rejected for Sprint 10 since the IV MRSD path uses CL-dominated Cp profiles where Kp mainly shifts Vss and t1/2 rather than AUC in the PAD window — expanding scope to include Kp would mean re-running R&R vs. alternate (Poulin-Theil, Schmitt) implementations, which is a research effort of its own.

### 4.2 Alternate liver models (ticket item 2)

`src/charon/core/liver_models.py` already exposes `well_stirred`, `parallel_tube`, `dispersion`. Sprint 10 computes all three MRSDs per compound for the decomposition. This is analytical: it does not change which model production uses. Compounds with `fu_p < 0.03` (warfarin, diazepam, atorvastatin) are highlighted in §3 of the report since well-stirred is most sensitive to fu_b near-zero.

### 4.3 target_ceff_nM review (ticket item 4)

Lisinopril reference FIH = 10 mg, predicted MRSD = 0.75 mg (fold 13.4x). CLrenal path works but `target_ceff_nM` = 170 nM in panel.yaml may be too low; a brief literature cross-check (Beermann 1988 and Gomez 1985 steady-state Cp data) is folded into the report as §4 and either confirms the current value or recommends an updated value. No YAML edit in this sprint — the recommendation is filed as a follow-up ticket.

## 5. Components

### 5.1 `src/charon/translational/decomposition.py` (new)

Pure-function library. No I/O except what the orchestrator injects.

Public API:

- `decompose_fold_error(mrsd_ws: float, mrsd_pt: float, mrsd_disp: float, f_lit: float | None, route_ref: str, fih_reference_mg: float) -> DecompositionResult` — returns a dataclass with `fold_observed`, `fold_liver_model`, `fold_route_bias`, `fold_residual`, `best_alt_model_name`, and `flags: list[str]`.
- `select_best_alternate_liver_model(ws: float, pt: float, disp: float, reference: float) -> tuple[str, float]` — returns `(name, mrsd)` of whichever of {parallel_tube, dispersion} minimises `|log10(mrsd/reference)|`. If both are farther from reference than well_stirred, returns `("well_stirred", ws)` (factor = 1.0).
- `compute_route_bias_factor(route_ref: str, f_lit: float | None) -> tuple[float, list[str]]` — returns `(factor, flags)`. IV reference → `(1.0, [])`. Oral + known F → `(1/F, [])`. Oral + unknown F → `(1.0, ["f_unknown"])`.

Dependencies: `math` only. No scipy, no pydantic, no RDKit. Small file (~150 LOC).

### 5.2 `validation/data/fih_reference/bioavailability.csv` (new)

Curated literature-F table for the 12 Tier A compounds. Columns:

| Column | Type | Example |
|---|---|---|
| `compound` | str | `midazolam` |
| `fih_reference_route` | str | `oral` or `iv` |
| `f_oral` | float \| empty | `0.40` |
| `f_source` | str | `Heizmann 1983` |
| `f_doi_or_pmid` | str | `PMID:6352420` |
| `notes` | str | `"p.o. 15 mg in healthy volunteers, F from IV cross-over"` |

Twelve rows, one per Tier A compound. The `fih_reference_route` field captures the route of the *clinical FIH reference dose* in `panel.yaml` (independent of panel.yaml's simulation `route`, which is uniformly `iv_bolus` since Sprint 7). Pure-IV references leave `f_oral` empty. All sources must be primary literature or a well-sourced secondary (Obach 1999 supplement is acceptable when it cites primary).

### 5.3 `validation/benchmarks/layer3_ivive_decomposition.py` (new)

Orchestrator script. Reads panel.yaml, bioavailability.csv, and Sprint 9 `layer3_fih_dose.json`. For each compound:

1. Load baseline MRSD from Sprint 9 JSON (well_stirred, the production default).
2. Recompute MRSD with `parallel_tube` and `dispersion` by re-running the translational path with the alternate model keyword passed to `liver_models.get_liver_model(name)`. Only the liver-model call site changes; everything upstream (Layer 1 predictions, fu_p, Kp) is held fixed.
3. Look up `f_oral` and `fih_reference_route` from bioavailability.csv.
4. Call `decompose_fold_error()`.
5. Aggregate and write `validation/reports/layer3_ivive_decomposition.{json, md}`.

Exit code 0 on success, 1 on CSV-schema violation, 1 on decomposition failure. No `sys.exit(2)` (that's reserved for Tier B sanity-floor failures in Sprint 7's benchmark and Sprint 10 does not gate).

### 5.4 Report structure: `validation/reports/layer3_ivive_decomposition.md`

Sections (numbered, markdown tables):

- **§1 Summary** — one-paragraph recap of Sprint 9 result + this report's role. Aggregate attribution percentages defined as `100 × Σ|log10(factor_k)| / Σ|log10(fold_observed)|` summed over the 12 compounds. Example rendering: "Of the cumulative log-fold-error across 12 compounds, liver-model choice accounts for N%, 1/F route bias accounts for N%, residual accounts for N%." The three percentages may sum to more than 100% when factors partially cancel within compounds — noted in the report.
- **§2 Per-compound decomposition table** — 12 rows. Columns: compound, `fold_obs`, `fold_liver`, `fold_route`, `fold_residual`, `best_alt_model`, `flags`. Sorted by `fold_residual` descending (worst-unexplained first).
- **§3 Pattern analysis** — grouped by the five root-cause categories from the Sprint 10 ticket. Each group: compounds in that group, typical decomposition signature, literature context.
- **§4 `target_ceff_nM` review (lisinopril)** — literature summary, current value 170 nM, recommended range or confirmation, citation.
- **§5 Sprint 11+ priority ranking** — ranked list: for each sprint candidate (Papp/oral migration; OATP plumbing; non-hepatic elimination tooling), estimated Tier A fold-error reduction if the attributable factor were zeroed, number of compounds affected.

Also emit `validation/reports/layer3_ivive_decomposition.json` with the per-compound structured breakdown for machine consumption.

## 6. Data flow

```
panel.yaml (12 compounds) ─┐
bioavailability.csv ──────┤
layer3_fih_dose.json ─────┤
                          ↓
              re-run MRSD path per compound × 3 liver models
                          ↓
              decompose_fold_error() per compound
                          ↓
              aggregate + pattern analysis
                          ↓
              report.md + report.json
```

## 7. Testing

### 7.1 Unit tests — `tests/unit/translational/test_decomposition.py`

1. `test_additivity_invariant_synthetic` — synthetic fold factors; `log10(fold_obs) == log10(liver) + log10(route) + log10(residual)` within `rtol=1e-9`.
2. `test_select_best_alternate_picks_closest_to_reference` — synthetic MRSDs; ensures correct model picked when parallel_tube is closer.
3. `test_select_best_alternate_falls_back_to_well_stirred_when_no_improvement` — when both alternates are further from reference, returns `well_stirred`.
4. `test_route_bias_iv_reference_returns_unity` — `route_ref="iv"`, any F → `(1.0, [])`.
5. `test_route_bias_oral_with_known_F` — `route_ref="oral"`, `f_lit=0.26` → factor ≈ 3.846, no flags.
6. `test_route_bias_oral_with_unknown_F` — `route_ref="oral"`, `f_lit=None` → `(1.0, ["f_unknown"])`.
7. `test_decompose_preserves_observed_fold_end_to_end` — using real Sprint 9 numbers for one compound (midazolam), verify that the returned `fold_observed` matches the input ratio.

### 7.2 Integration test — `tests/integration/test_layer3_decomposition_benchmark.py`

1. `test_decomposition_script_runs_all_12_compounds` — runs the orchestrator (may take a minute), asserts `report.json` has 12 entries, each with all required keys, all numeric.
2. `test_decomposition_log_sum_matches_observed` — for each of 12 compounds in `report.json`, asserts log-additivity within `rtol=1e-6`.
3. `test_bioavailability_csv_covers_tier_a` — reads `bioavailability.csv`, asserts every Tier A compound in panel.yaml has a row (F value may be empty, but the row must exist with a citation).

### 7.3 No test for

- Specific decomposition percentages (would be re-tuning the report on itself).
- Lisinopril `target_ceff_nM` literature recommendation (a literature claim, not a code property).

## 8. Success criteria

- 12/12 Tier A compounds have a complete decomposition record (numeric or explicitly flagged).
- Report identifies which factor (liver model, 1/F route, residual) dominates aggregate error in log-space.
- Sprint 11 candidate (Papp/oral migration) has a quantified expected Tier A improvement: count of compounds whose `fold_route_bias` is the largest factor, × the magnitude if zeroed.
- §8 target remains FAIL after Sprint 10 — explicitly stated in the report's §1.
- All unit + integration tests pass.
- No changes to files under `src/charon/core/liver_models.py`, `src/charon/translational/pad.py`, or `validation/data/fih_reference/panel.yaml`.

## 9. Risk register

| Risk | Mitigation |
|---|---|
| Decomposition explains <50% of aggregate log-fold-error (residual dominates) | Report this outcome honestly; it means Sprint 11 (1/F fix) is insufficient, escalate the finding. This is a valid research outcome, not a failure. |
| Literature F values disagree across sources for a compound (e.g. propranolol 0.25–0.40 range) | Pick median of cited range, note the range in `notes`, sensitivity-check in report §3. |
| Parallel-tube or dispersion numerically unstable at very low fu_b | Already validated in Sprint 3 tests; if any compound breaks, fall back to well_stirred with `"liver_model_numerical"` flag. |
| Reader interprets "fold_route_bias" as "Charon is wrong about F" when it is a reference-path comparison artefact | Report §1 explicitly distinguishes *route-comparison artefact* from *IVIVE error*. |
| Sprint 11 scope expectation inflation based on Sprint 10 numbers | Report §5 uses conservative "if factor → 1.0" bounds (upper bound on improvement), labels them upper bounds. |

## 10. Files touched (exhaustive)

**New:**
- `src/charon/translational/decomposition.py`
- `validation/data/fih_reference/bioavailability.csv`
- `validation/benchmarks/layer3_ivive_decomposition.py`
- `validation/reports/layer3_ivive_decomposition.md`
- `validation/reports/layer3_ivive_decomposition.json`
- `tests/unit/translational/test_decomposition.py`
- `tests/integration/test_layer3_decomposition_benchmark.py`

**Not touched (confirmed):**
- `src/charon/core/liver_models.py`
- `src/charon/translational/pad.py`, `dose_projector.py`
- `validation/data/fih_reference/panel.yaml`
- Any file under `src/charon/predict/`, `src/charon/pbpk/`, `src/charon/uncertainty/`

## 11. Follow-up tickets (filed during Sprint 10, not resolved in-sprint)

- `sprint10-followup-lisinopril-ceff.md` — proposed `target_ceff_nM` update with literature citation, pending user decision.
- Update `sprint10-ivive-bias-ticket.md` status to "diagnosed (Sprint 10 merged); remediation in Sprint 11 / 12".

## 12. Estimated test count

Sprint 9 ended at 882 tests. Sprint 10 adds ~10 tests (7 unit + 3 integration). Expected end state: ~892 tests.
