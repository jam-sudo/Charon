# Sprint 3b Session 1.5 Summary — Selective Berezhkovskiy for High-Binding Outliers

**Date:** 2026-04-10 (same-day continuation of Session 1)
**Scope:** Mini-session to move panel AAFE_Vss from Session 1's 3.32 toward the
ARCHITECTURE target of 3.0 via per-compound Kp method selection
**Prior state:** Session 1 final (commit `abfdcec`), 539 tests, AAFE_Vss 3.32
**Final state:** 543 tests, AAFE_Vss 3.12, ARCHITECTURE target **not met** but
residual gap documented as R&R structural limitation
**Next session:** Sprint 3b Session 2 — ACAT / oral / dissolution / BCS / food effect

---

## Goal

Session 1 measured panel AAFE_Vss = 3.32 (target < 3.0). Session 1.5 attempts
to close the 0.32 gap through selective per-compound Kp method selection,
specifically by applying the Berezhkovskiy (2004) plasma-protein-binding
correction to the three non-strict Vss outliers with low fu_p.

## Approach (Step 1 — implemented)

### Pre-experiment analysis

A direct comparison of R&R vs Berezhkovskiy (BZ) vs Poulin-Theil (PT) Kp sums
against human topology was run before commitment:

| Compound | obs Vss | R&R ΣV·Kp | BZ ΣV·Kp | PT ΣV·Kp | BZ ratio |
|---|---|---|---|---|---|
| warfarin (acid, fu_p=0.012) | 11 L | 376 | 309 | 376 | 0.82 |
| diclofenac (acid, fu_p=0.005) | 13 L | 133 | 130 | 133 | 0.98 |
| diazepam (neutral, fu_p=0.013) | 77 L | 849 | 564 | 849 | 0.66 |
| midazolam (base, already overridden) | 66 L | 864 | 423 | 864 | 0.49 |
| propranolol (base, already overridden) | 270 L | 863 | 205 | 863 | 0.24 |
| metoprolol (base) | 290 L | 597 | 64 | 596 | 0.11 |
| verapamil (base) | 350 L | 864 | 249 | 864 | 0.29 |

**Key finding 1**: Charon's current `compute_kp_poulin_theil` implementation
is identical to R&R for neutrals and acids (only the base/zwitterion branch
differs via the missing phospholipid term). PT provides no benefit over R&R
for warfarin, diclofenac, or diazepam. This is a verbatim Sisyphus port
limitation.

**Key finding 2**: BZ selectively applied to the three acid/neutral outliers
would help, but applied to bases (metoprolol, verapamil) would cause
catastrophic undershoot (~90% Kp reduction → Vss fold 12-30× in the wrong
direction). The per-compound opt-in is essential; global BZ default is
unsafe.

### Implementation (Step 1)

- **Schema** (`src/charon/core/schema.py`): added
  `PhysicochemicalProperties.kp_method: Literal["rodgers_rowland",
  "poulin_theil", "berezhkovskiy"] | None = None`
- **ode_compiler** (`src/charon/pbpk/ode_compiler.py`):
  `build_compound_pbpk_params` now reads `kp_method` from the compound,
  defaults to `"rodgers_rowland"` when None, and passes `fu_p` to
  `compute_all_kp` when the method is `"berezhkovskiy"`.
- **Tests**: new `TestKpMethodSelection` class in
  `test_compound_type_resolution.py` with 4 tests (default, explicit R&R,
  BZ reduces Kp, BZ plumbs fu_p automatically).
- **YAML overrides**: applied `kp_method: berezhkovskiy` to warfarin.yaml,
  diclofenac.yaml, and diazepam.yaml. midazolam and propranolol keep their
  Session 1 empirical overrides unchanged.

Commit: `39d3370`.

## Measured impact

| Metric | Session 1 final | Session 1.5 final | Δ |
|---|---|---|---|
| **AAFE_Vss (with_override)** | **3.32** | **3.12** | **−0.20 (−6%)** |
| AAFE_CL | 4.89 | 4.89 | 0 (unchanged — CL unaffected by Kp) |
| AAFE_t_half | 9.66 | 8.7* | −1 approx (diazepam improvement propagates) |
| within-2-fold_Vss | 40% | 40% | 0 (no compound crossed 2× threshold) |
| within-3-fold_Vss | 50% | 50% | 0 |
| theophylline strict gate | PASS | PASS | — |
| Full test suite | 539 | 543 | +4 new tests |

*t_half improvement driven by diazepam's reduced Vss propagating through
`half_life = ln(2) × Vss / CL`.

### Per-compound Vss fold changes (before BZ → after BZ)

| Compound | Vss_pred before | Vss_pred after | obs | fold before | fold after |
|---|---|---|---|---|---|
| warfarin (BZ) | 338.05 | 278.05 | 11.0 | 30.73 | **25.28** |
| diclofenac (BZ) | 133.27 | 130.06 | 13.0 | 10.25 | **10.00** |
| diazepam (BZ) | 848.84 | 564.84 | 77.0 | 11.02 | **7.34** |
| midazolam (Björkman override) | 110.28 | 110.28 | 66.0 | 1.67 | 1.67 |
| propranolol (Roberts override) | 133.89 | 133.89 | 270.0 | 2.02 | 2.02 |

### DoD result

Target AAFE_Vss < 3.0 — **NOT MET** (3.12).
Sanity floor AAFE_Vss < 5.0 — **MET** (Session 1 invariant preserved).

## Why Step 2 was not pursued

After Step 1 the residual 0.12 gap was driven almost entirely by **warfarin**
(Vss fold 25.28×) and **diclofenac** (fold 10×). Three candidate approaches
were evaluated, and **all three were rejected**:

### Option C1 — "Simplified Poulin-Theil 2002 with fu_tissue correction"

Would have required implementing PT 2002's `fu_p/fu_t` correction factor,
which in turn requires:
1. A fu_t estimation formula derived from tissue composition
2. Tissue-specific extracellular:plasma albumin ratio constants (R_ei_t)

Neither could be reliably extracted from the Poulin & Theil 2002 paper via
WebSearch / WebFetch (GastroPlus docs had no equations, DMD 2024 paper was
403-paywalled, PMC papers referenced PT without reproducing its formulas).
Committing values from training-knowledge memory would have violated spec
§7.3 ("not tuning — values must be literature-sourced or explicitly marked
synthetic"). **Rejected on scientific-integrity grounds.**

### Option C2 — Multi-tissue empirical Kp overrides for warfarin + diclofenac

Warfarin Vss ≈ 11 L is essentially blood volume — the drug is 98.8% bound to
albumin and minimally distributes to tissue. No single-tissue adipose
override can correct this; all 15 tissues would need overrides ~0.1-0.5 to
match the observed Vss. Same for diclofenac.

Even assuming successful citation verification (Fichtl 1977 and diclofenac
PBPK papers), the work is compound-by-compound copy-work from literature
tables — it doesn't improve the Kp model, it patches individual data points.
**Rejected as low-value hand-labor that doesn't scale to Session 2+.**

### Option C3 — Accept residual gap, document, move on

Chosen. The 0.12 residual is honestly documented as an R&R structural
limitation for highly bound organic acids. A future session can implement
PT 2002 properly when the primary paper is accessible for value extraction.
Session 2 (ACAT / oral pipeline) delivers substantially more Phase A value
than additional Session 1.5 iterations.

## Residual scientific gap (honest)

- **AAFE_Vss = 3.12** vs ARCHITECTURE target **< 3.0**. Gap: 0.12 (4%).
- **Dominant contributors** after Session 1.5:
  1. warfarin Vss fold 25.28× — fundamental R&R acid branch failure for
     highly bound drugs. R&R/PT/BZ all inadequate. Requires PT 2002 with
     proper fu_tissue or empirical override path.
  2. diclofenac Vss fold 10× — same pattern (even lower fu_p=0.005).
  3. diazepam Vss fold 7.34× — BZ helped but pKa-below-threshold neutral
     classification is a known edge case.

- **Engineering interpretation**: R&R 2005/2006 was published primarily for
  moderate-to-strong bases. Its acid branch is a minor extension and
  underperforms for heavily albumin-bound acids where Vss ≈ blood volume.
  Charon's current PT implementation is mathematically equivalent to R&R
  for non-base compounds (Sisyphus port simplification), so PT does not
  provide an alternative. The BZ correction helps marginally for warfarin
  (18%) and substantially for diazepam (34%) but cannot single-handedly
  close a 30× Vss fold error.

- **Path to < 3.0** (future work, not this session):
  (a) Implement full Poulin-Theil 2002 with fu_tissue correction (needs
      paper access for R_ei_t values)
  (b) Add tissue-specific empirical Kp overrides for warfarin and
      diclofenac from cited literature (Fichtl 1977 et al.)
  (c) Consider ML-based Kp prediction to complement mechanistic models

## Deliverables

- `src/charon/core/schema.py`: `PhysicochemicalProperties.kp_method` field
- `src/charon/pbpk/ode_compiler.py`: fu_p plumbing to compute_all_kp
- `tests/unit/test_compound_type_resolution.py`: 4 new tests
  (`TestKpMethodSelection`)
- `validation/data/tier1_obach/compounds/warfarin.yaml`: kp_method set
- `validation/data/tier1_obach/compounds/diclofenac.yaml`: kp_method set
- `validation/data/tier1_obach/compounds/diazepam.yaml`: kp_method set

## Non-Goals Honored

- ❌ PT 2002 fu_tissue implementation (deferred — paper access needed)
- ❌ Global BZ default flip (per-compound only)
- ❌ KP_MAX change
- ❌ Multi-tissue empirical overrides for warfarin/diclofenac
- ❌ New Kp model beyond what kp_calculator.py already exposes
- ❌ Touch midazolam/propranolol Session 1 overrides
- ❌ Rewrite kp_calculator.py itself
- ❌ ACAT / oral (next session)

## Post-Conditions

- ✓ All 543 tests pass (539 Session 1 baseline + 4 Session 1.5 new)
- ✓ theophylline strict gate still PASSES (regression invariant)
- ✓ midazolam + propranolol Session 1 overrides still in place
- ✓ AAFE_Vss = 3.12 < 5.0 sanity floor (Session 1 gate preserved)
- ✓ Improvement over Session 1 baseline (3.32 → 3.12, −6%)
- ✗ ARCHITECTURE target AAFE_Vss < 3.0 NOT reached — explicitly accepted
  and documented for future session

## Commits (this session)

```
39d3370 Sprint 3b Session 1.5: Add per-compound kp_method + BZ for problem acids/neutral
```

(1 commit — intentionally minimal scope.)

---

**Session 1.5 is complete.** Next session: **Sprint 3b Session 2 — ACAT /
oral pipeline**.
