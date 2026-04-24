# Sprint 11 — Tier A Oral Route Migration Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Migrate Layer 3 Tier A panel from `iv_bolus` to `oral` route, curating Peff/Papp for 11 compounds, to eliminate the 42% route_bias aggregate identified by Sprint 10.

**Architecture:** Data curation (11 compound YAMLs gain `permeability` block with literature-sourced Peff) + panel config flip (Tier A route: oral) + a small orchestrator update (`route_bias_override` kwarg on `decompose_fold_error`; simulation-route tracking in the decomposition orchestrator). No production code path changes beyond the one orchestrator tweak.

**Tech Stack:** PyYAML, pytest. WebSearch/WebFetch used for literature verification during curation.

**Spec:** `docs/superpowers/specs/2026-04-24-sprint11-oral-route-migration-design.md`

---

## File Structure

**Modified compound YAMLs (11 files):**

| File | Add block |
|---|---|
| `validation/data/tier1_obach/compounds/warfarin.yaml` | `permeability.peff_cm_s` |
| `validation/data/tier1_obach/compounds/propranolol.yaml` | `permeability.peff_cm_s` |
| `validation/data/tier1_obach/compounds/verapamil.yaml` | `permeability.peff_cm_s` |
| `validation/data/tier1_obach/compounds/omeprazole.yaml` | `permeability.peff_cm_s` |
| `validation/data/tier1_obach/compounds/theophylline.yaml` | `permeability.peff_cm_s` |
| `validation/data/tier1_obach/compounds/diclofenac.yaml` | `permeability.peff_cm_s` |
| `validation/data/tier1_obach/compounds/diazepam.yaml` | `permeability.peff_cm_s` |
| `validation/data/tier1_obach/compounds/metoprolol.yaml` | `permeability.peff_cm_s` |
| `validation/data/tier1_obach/compounds/acetaminophen.yaml` | `permeability.peff_cm_s` |
| `validation/data/tier1_obach/compounds/lisinopril.yaml` | `permeability.peff_cm_s` |
| `validation/data/tier1_obach/compounds/atorvastatin.yaml` | `permeability.peff_cm_s` |

(midazolam, felodipine, nifedipine already have permeability from Sprint 3b-2a.)

**Other modifications:**

| File | Change |
|---|---|
| `validation/data/fih_reference/panel.yaml` | Tier A 12 entries: `route: iv_bolus` → `route: oral`. Update `route_note`. |
| `src/charon/translational/decomposition.py` | Add `route_bias_override: float \| None = None` kwarg to `decompose_fold_error` |
| `validation/benchmarks/layer3_ivive_decomposition.py` | Track `simulation_route` per row, pass `route_bias_override=1.0` when routes match |
| `tests/unit/test_fih_panel_schema.py` | +2 guards (Tier A has permeability; Tier A panel route=oral) |
| `tests/unit/test_decomposition.py` | +1 test (route_bias_override bypasses F lookup) |
| `tests/integration/test_fih_pipeline.py` | +1 parametric smoke test (12 Tier A oral runs succeed) |
| `validation/reports/layer3_fih_dose.{md,json}` | Regenerated from oral panel |
| `validation/reports/layer3_ivive_decomposition.{md,json}` | Regenerated with simulation_route awareness |

---

## Task 1: Curate Peff for 11 Tier A compounds

**Files:** 11 YAMLs under `validation/data/tier1_obach/compounds/` (warfarin, propranolol, verapamil, omeprazole, theophylline, diclofenac, diazepam, metoprolol, acetaminophen, lisinopril, atorvastatin)

**Literature values to use** (Peff in cm/s, derived from published human-jejunum perfusion or Caco-2 Papp via `papp_to_peff()` correlation):

| Compound | Peff (cm/s) | Source (primary citation hint) | Confidence |
|---|---|---|---|
| warfarin | 3.8e-4 | Fagerholm 1996 / Winiwarter 1998 compilation | med |
| propranolol | 2.9e-4 | Winiwarter 1998 | high |
| verapamil | 6.8e-4 | Varma 2010 J Med Chem 53:1098 (Table 1) | med |
| omeprazole | 2.5e-4 | Sun 2002 (Caco-2 compilation) | med |
| theophylline | 3.3e-4 | Winiwarter 1998 | high |
| diclofenac | 2.2e-4 | Sun 2002 (Caco-2 Papp ~30e-6 → Peff) | med |
| diazepam | 3.0e-4 | Lennernäs 2007 Eur J Pharm Sci 29(3-4):278 | high |
| metoprolol | 1.3e-4 | Winiwarter 1998 / Amidon | high |
| acetaminophen | 3.1e-4 | Winiwarter 1998 | high |
| lisinopril | 0.3e-4 | Knutter 2008 (PEPT1 substrate; low passive Peff) | low |
| atorvastatin | 0.45e-4 | Wu 2000 / Lennernäs 2003 (low passive; OATP1B1 uptake) | low |

**Verification responsibility:** Before committing, the implementer should use WebSearch/WebFetch to confirm each value against its cited source when confidence is "med" or "low". If the literature disagrees with the above by more than 3-fold, prefer the verified primary-source value and note the discrepancy in the YAML `method` string. If the implementer cannot verify a value (offline, inaccessible source), use the value from the table above and annotate `source: compilation` (not `literature`).

- [ ] **Step 1: Verify each Peff value (11 compounds)**

For each of the 11 compounds, attempt `WebSearch "<compound> peff jejunum human"` or similar. Note the verified value. If the verified value is within 2-fold of the suggestion table, use the table value for consistency. If it disagrees by >2-fold, use the verified value and cite the primary paper.

No commit yet — this step is research only.

- [ ] **Step 2: Add `permeability` block to each of the 11 compound YAMLs**

The exact edit pattern (replace `<VALUE>`, `<CITATION_STRING>` with the verified values from Step 1):

Locate the `properties:` block in each YAML. Find the location between `binding:` and `metabolism:` (or wherever makes sense for alphabetical / already-established ordering — follow the pattern in `midazolam.yaml` lines 24-28).

Add:

```yaml
  permeability:
    peff_cm_s:
      value: <VALUE>
      source: literature
      unit: cm/s
      method: "<CITATION_STRING>"
```

Apply this pattern to ALL 11 compound YAMLs.

For compounds where only compilation (secondary) source is available, use `source: compilation` instead of `literature`.

- [ ] **Step 3: Verify YAMLs parse correctly**

```bash
python3 -c "
import yaml
from pathlib import Path
for p in sorted(Path('validation/data/tier1_obach/compounds').glob('*.yaml')):
    data = yaml.safe_load(p.read_text())
    perm = data.get('properties', {}).get('permeability', {})
    print(f\"{p.stem:15s}  peff={perm.get('peff_cm_s', {}).get('value', 'MISSING')}\")
"
```

Expected: 14 lines (all compounds in the folder). For the 11 Tier A compounds + midazolam + felodipine + nifedipine, Peff value printed. For antipyrine, caffeine, dextromethorphan (Tier B or unused), no permeability — print "MISSING".

- [ ] **Step 4: Run existing tests to confirm no regression**

```bash
pytest tests/unit/test_fih_new_compound_yamls.py tests/unit/test_compound_config.py -q
```

Expected: pass.

- [ ] **Step 5: Commit**

```bash
git add validation/data/tier1_obach/compounds/warfarin.yaml \
        validation/data/tier1_obach/compounds/propranolol.yaml \
        validation/data/tier1_obach/compounds/verapamil.yaml \
        validation/data/tier1_obach/compounds/omeprazole.yaml \
        validation/data/tier1_obach/compounds/theophylline.yaml \
        validation/data/tier1_obach/compounds/diclofenac.yaml \
        validation/data/tier1_obach/compounds/diazepam.yaml \
        validation/data/tier1_obach/compounds/metoprolol.yaml \
        validation/data/tier1_obach/compounds/acetaminophen.yaml \
        validation/data/tier1_obach/compounds/lisinopril.yaml \
        validation/data/tier1_obach/compounds/atorvastatin.yaml
git commit -m "data(sprint11): curate Peff for 11 Tier A compounds (Obach/Winiwarter/Lennernas)"
```

---

## Task 2: Schema guard — Tier A YAMLs have permeability, panel route = oral

**Files:**
- Modify: `tests/unit/test_fih_panel_schema.py`

- [ ] **Step 1: Read the existing schema test file**

```bash
cat tests/unit/test_fih_panel_schema.py
```

Familiarise with the existing pattern (`panel = yaml.safe_load(...)`, extraction of tier="gold" subset).

- [ ] **Step 2: Write the failing test**

Append these two tests to `tests/unit/test_fih_panel_schema.py`:

```python
def test_tier_a_compounds_have_permeability():
    """Sprint 11: Tier A oral migration requires Peff or Papp in each
    compound YAML."""
    panel_path = Path(__file__).resolve().parents[2] / "validation" / "data" / "fih_reference" / "panel.yaml"
    compounds_dir = Path(__file__).resolve().parents[2] / "validation" / "data" / "tier1_obach" / "compounds"
    panel = yaml.safe_load(panel_path.read_text())["panel"]
    tier_a_names = {c["name"] for c in panel["compounds"] if c["tier"] == "gold"}
    missing = []
    for name in sorted(tier_a_names):
        yaml_path = compounds_dir / f"{name}.yaml"
        data = yaml.safe_load(yaml_path.read_text())
        perm = data.get("properties", {}).get("permeability", {})
        has_peff = "peff_cm_s" in perm and perm["peff_cm_s"].get("value") is not None
        has_papp = "papp_nm_s" in perm and perm["papp_nm_s"].get("value") is not None
        if not (has_peff or has_papp):
            missing.append(name)
    assert not missing, (
        f"Tier A compounds missing Peff/Papp: {missing}. "
        f"Sprint 11 requires oral-route permeability data."
    )


def test_tier_a_panel_route_is_oral():
    """Sprint 11: Tier A panel entries simulate oral route."""
    panel_path = Path(__file__).resolve().parents[2] / "validation" / "data" / "fih_reference" / "panel.yaml"
    panel = yaml.safe_load(panel_path.read_text())["panel"]
    wrong_route = [
        c["name"] for c in panel["compounds"]
        if c["tier"] == "gold" and c["route"] != "oral"
    ]
    assert not wrong_route, (
        f"Tier A compounds with non-oral route: {wrong_route}. "
        f"Sprint 11 expects route=oral for all gold-tier entries."
    )
```

If `Path` is not already imported at the top of the test file, add `from pathlib import Path` and `import yaml`.

- [ ] **Step 3: Run tests to verify they fail (expected: the route test fails; the permeability test might pass after Task 1)**

```bash
pytest tests/unit/test_fih_panel_schema.py -v
```

Expected:
- `test_tier_a_panel_route_is_oral` FAILS: panel still has `iv_bolus`
- `test_tier_a_compounds_have_permeability` PASSES (Task 1 already committed the YAML edits)
- Pre-existing tests in the file PASS

- [ ] **Step 4: (Implementation comes in Task 3 — keep the failing test for now)**

No action. Proceed to Task 3 where the panel.yaml route flip happens.

- [ ] **Step 5: Commit (tests only)**

```bash
git add tests/unit/test_fih_panel_schema.py
git commit -m "test(sprint11): schema guards for Tier A permeability + oral route"
```

---

## Task 3: Flip panel.yaml Tier A route to oral

**Files:**
- Modify: `validation/data/fih_reference/panel.yaml`

- [ ] **Step 1: Read the current panel.yaml**

```bash
head -100 validation/data/fih_reference/panel.yaml
```

Confirm the structure: Tier A compounds (tier: gold) have `route: iv_bolus`.

- [ ] **Step 2: Apply the edit — change `route: iv_bolus` → `route: oral` for Tier A gold entries only**

**Important:** Tier B (sanity_floor) entries MUST keep `route: iv_bolus`. Only gold-tier entries flip to oral.

Method: use `sed` scoped to the Tier-A block, or use a Python script to be safe:

```bash
python3 <<'EOF'
import yaml
from pathlib import Path

path = Path("validation/data/fih_reference/panel.yaml")
data = yaml.safe_load(path.read_text())
count = 0
for compound in data["panel"]["compounds"]:
    if compound["tier"] == "gold" and compound.get("route") == "iv_bolus":
        compound["route"] = "oral"
        count += 1
print(f"Flipped {count} Tier A entries to oral")

# Write back preserving YAML structure
path.write_text(yaml.safe_dump(data, sort_keys=False, default_flow_style=False, width=100))
EOF
```

Expected stdout: `Flipped 12 Tier A entries to oral`.

After running, verify Tier B is untouched:
```bash
grep -E "(name|tier|route):" validation/data/fih_reference/panel.yaml | head -40
```

Tier A entries should show `route: oral`; Tier B entries should show `route: iv_bolus`.

- [ ] **Step 3: Update the `route_note` block**

Open `validation/data/fih_reference/panel.yaml` in an editor. Find the `route_note: |` block (which currently explains the IV simulation caveat) and replace its content with:

```yaml
  route_note: |
    Sprint 11 (2026-04-24): Tier A migrated to oral route after curating
    Peff for 11 compounds (Obach 2008 / Winiwarter 1998 / Lennernäs 2007
    primary sources). This eliminates the 1/F comparison artefact Sprint 10
    identified (route_bias aggregate 42%). Tier B remains iv_bolus for
    compatibility with Sprint 7 sanity-floor tests; migration is a separate
    future sprint.
```

- [ ] **Step 4: Run the schema tests from Task 2 (should now pass)**

```bash
pytest tests/unit/test_fih_panel_schema.py -v
```

Expected: ALL tests pass, including `test_tier_a_panel_route_is_oral`.

- [ ] **Step 5: Commit**

```bash
git add validation/data/fih_reference/panel.yaml
git commit -m "data(sprint11): flip Tier A panel route iv_bolus -> oral"
```

---

## Task 4: `route_bias_override` kwarg on `decompose_fold_error`

**Files:**
- Modify: `src/charon/translational/decomposition.py`
- Modify: `tests/unit/test_decomposition.py`

**Rationale:** When the simulation route equals the clinical reference route (e.g., both oral), there is no 1/F comparison artefact. The orchestrator needs to signal this; a `route_bias_override` kwarg lets it pass `1.0` directly.

- [ ] **Step 1: Write the failing test**

Append to `tests/unit/test_decomposition.py`:

```python
def test_decompose_with_route_bias_override_bypasses_f_lookup():
    """When simulation_route == reference_route, the orchestrator passes
    route_bias_override=1.0 and the decomposition returns fold_route_bias=1.0
    regardless of f_lit."""
    result = decompose_fold_error(
        mrsd_ws=10.0,
        mrsd_pt=5.0,
        mrsd_disp=15.0,
        f_lit=0.5,         # Would normally produce fold_route_bias=2.0
        route_ref="oral",
        fih_reference_mg=5.0,
        route_bias_override=1.0,
    )
    assert result.fold_route_bias == 1.0
    assert result.flags == ()


def test_decompose_with_route_bias_override_none_uses_default():
    """override=None (default) preserves the original F-lookup path."""
    result = decompose_fold_error(
        mrsd_ws=10.0,
        mrsd_pt=5.0,
        mrsd_disp=15.0,
        f_lit=0.5,
        route_ref="oral",
        fih_reference_mg=5.0,
        route_bias_override=None,
    )
    assert result.fold_route_bias == pytest.approx(2.0, rel=1e-12)
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
pytest tests/unit/test_decomposition.py::test_decompose_with_route_bias_override_bypasses_f_lookup -v
```

Expected: FAIL with `TypeError: decompose_fold_error() got an unexpected keyword argument 'route_bias_override'`.

- [ ] **Step 3: Implement the override kwarg**

Open `src/charon/translational/decomposition.py`. Find `decompose_fold_error` (around line 135). Extend its signature and body:

**Before:**
```python
def decompose_fold_error(
    mrsd_ws: float,
    mrsd_pt: float,
    mrsd_disp: float,
    f_lit: float | None,
    route_ref: Literal["iv", "oral"],
    fih_reference_mg: float,
) -> DecompositionResult:
```

Add the new kwarg at the end:

```python
def decompose_fold_error(
    mrsd_ws: float,
    mrsd_pt: float,
    mrsd_disp: float,
    f_lit: float | None,
    route_ref: Literal["iv", "oral"],
    fih_reference_mg: float,
    route_bias_override: float | None = None,
) -> DecompositionResult:
```

In the function body, replace the existing `compute_route_bias_factor` call:

**Before:**
```python
    fold_route, flags = compute_route_bias_factor(
        route_ref=route_ref, f_lit=f_lit
    )
```

**After:**
```python
    if route_bias_override is not None:
        if route_bias_override <= 0:
            raise ValueError(
                f"route_bias_override must be > 0, got {route_bias_override}"
            )
        fold_route, flags = route_bias_override, ()
    else:
        fold_route, flags = compute_route_bias_factor(
            route_ref=route_ref, f_lit=f_lit
        )
```

Update the function docstring to document the new kwarg:

Find the `Args:` block in the docstring and add:
```
        route_bias_override: If provided, bypasses compute_route_bias_factor
            and uses this value directly. Used by the orchestrator when
            simulation_route == reference_route (no 1/F artefact).
```

- [ ] **Step 4: Run the new tests**

```bash
pytest tests/unit/test_decomposition.py -v
```

Expected: all 17+ tests PASS (existing 16 + 2 new).

- [ ] **Step 5: Commit**

```bash
git add src/charon/translational/decomposition.py tests/unit/test_decomposition.py
git commit -m "feat(sprint11): route_bias_override kwarg on decompose_fold_error"
```

---

## Task 5: Orchestrator update — track simulation_route, use override when routes match

**Files:**
- Modify: `validation/benchmarks/layer3_ivive_decomposition.py`

- [ ] **Step 1: Read the current orchestrator**

```bash
sed -n '100,180p' validation/benchmarks/layer3_ivive_decomposition.py
```

Find the `run_panel` function and locate the `decompose_fold_error(...)` call inside the `for entry in tier_a:` loop.

- [ ] **Step 2: Edit the orchestrator**

In `validation/benchmarks/layer3_ivive_decomposition.py`, inside `run_panel`, change the per-compound loop to:

**Before:**
```python
        result = decompose_fold_error(
            mrsd_ws=mrsds["well_stirred"],
            mrsd_pt=mrsds["parallel_tube"],
            mrsd_disp=mrsds["dispersion"],
            f_lit=bioav_row["f_oral"],
            route_ref=bioav_row["fih_reference_route"],
            fih_reference_mg=float(entry["reference_fih_mg"]),
        )
```

**After:**
```python
        simulation_route = "oral" if entry["route"].startswith("oral") else "iv"
        reference_route = bioav_row["fih_reference_route"]
        route_override = 1.0 if simulation_route == reference_route else None
        result = decompose_fold_error(
            mrsd_ws=mrsds["well_stirred"],
            mrsd_pt=mrsds["parallel_tube"],
            mrsd_disp=mrsds["dispersion"],
            f_lit=bioav_row["f_oral"],
            route_ref=reference_route,
            fih_reference_mg=float(entry["reference_fih_mg"]),
            route_bias_override=route_override,
        )
```

Add two columns to the row dict appended to `rows`:

Find the `rows.append({` block and add:
```python
            "simulation_route": simulation_route,
            "reference_route": reference_route,
```

- [ ] **Step 3: Verify the decomposition integration tests still pass**

```bash
pytest tests/integration/test_layer3_decomposition_benchmark.py -v
```

Expected: 5 tests pass. The log-additivity test in particular should still hold — route_bias=1.0 satisfies the invariant trivially (log10(1.0)=0).

If `test_liver_model_whatif_produces_nonzero_attribution` fails because liver factor is also ~0 for some runs, that's OK for this step — don't change it yet. The invariant is what matters.

- [ ] **Step 4: Run the orchestrator manually to see new numbers**

```bash
python3 validation/benchmarks/layer3_ivive_decomposition.py --output-stem /tmp/sprint11_decomp_check
```

Expected: `[OK] Decomposition wrote ... (12 compounds)`.

Inspect:
```bash
python3 -c "
import json
d = json.load(open('/tmp/sprint11_decomp_check.json'))
print('Summary:', d['summary'])
"
```

Expected: `aggregate_pct_route_bias` should be near 0% (simulation matches reference for all 12 compounds now). `aggregate_pct_residual` should now dominate the attribution.

- [ ] **Step 5: Commit**

```bash
git add validation/benchmarks/layer3_ivive_decomposition.py
git commit -m "feat(sprint11): orchestrator tracks simulation_route + uses override when routes match"
```

---

## Task 6: Smoke test — 12 Tier A compounds run oral pipeline

**Files:**
- Modify: `tests/integration/test_fih_pipeline.py`

- [ ] **Step 1: Read the existing integration test**

```bash
head -40 tests/integration/test_fih_pipeline.py
```

Note the imports and existing test patterns.

- [ ] **Step 2: Write the failing smoke test**

Append to `tests/integration/test_fih_pipeline.py`:

```python
@pytest.mark.parametrize("compound_name", [
    "midazolam", "warfarin", "propranolol", "verapamil", "omeprazole",
    "theophylline", "diclofenac", "diazepam", "metoprolol", "acetaminophen",
    "lisinopril", "atorvastatin",
])
def test_tier_a_oral_pipeline_runs_to_mrsd(compound_name):
    """Sprint 11 smoke: every Tier A compound completes Pipeline(route=oral)
    with a positive finite MRSD."""
    import math
    import yaml
    from pathlib import Path
    from charon import Pipeline
    from charon.core.schema import CompoundConfig, DoseProjectionConfig

    repo_root = Path(__file__).resolve().parents[2]
    compound_yaml = (
        repo_root / "validation" / "data" / "tier1_obach" / "compounds"
        / f"{compound_name}.yaml"
    )
    panel_yaml = repo_root / "validation" / "data" / "fih_reference" / "panel.yaml"

    compound = CompoundConfig.model_validate(yaml.safe_load(compound_yaml.read_text()))
    panel = yaml.safe_load(panel_yaml.read_text())["panel"]
    entry = next(
        c for c in panel["compounds"]
        if c["name"] == compound_name and c["tier"] == "gold"
    )

    pipe = Pipeline(
        compound,
        route=entry["route"],
        dose_mg=1.0,
        dose_projection=DoseProjectionConfig(
            target_ceff_nM=float(entry["target_ceff_nM"]),
            safety_factor=10.0,
            tau_h=24.0,
        ),
    )
    result = pipe.run()
    assert result.dose_recommendation is not None, f"{compound_name}: no dose rec"
    mrsd = float(result.dose_recommendation.mrsd_mg)
    assert mrsd > 0, f"{compound_name}: mrsd={mrsd} (must be > 0)"
    assert not math.isnan(mrsd), f"{compound_name}: mrsd is NaN"
    assert not math.isinf(mrsd), f"{compound_name}: mrsd is inf"
```

- [ ] **Step 3: Run the smoke test**

```bash
pytest tests/integration/test_fih_pipeline.py::test_tier_a_oral_pipeline_runs_to_mrsd -v
```

Expected: 12/12 PASS. Wall-clock ~30-60s (one Pipeline.run per compound, including PBPK simulation).

If any compound fails (e.g. ACAT throws on lisinopril's very-low Peff): STOP. This is a real integration issue — either the Peff curation in Task 1 is too low for ACAT's numerical stability, or there's a Pipeline bug. Report as BLOCKED and let the controller investigate.

- [ ] **Step 4: Commit**

```bash
git add tests/integration/test_fih_pipeline.py
git commit -m "test(sprint11): Tier A oral pipeline smoke (12 compounds)"
```

---

## Task 7: Regenerate benchmarks + append report narrative

**Files:**
- Regenerate: `validation/reports/layer3_fih_dose.{md,json}`
- Regenerate: `validation/reports/layer3_ivive_decomposition.{md,json}`

- [ ] **Step 1: Run the Layer 3 FIH dose benchmark (oral panel)**

```bash
python3 validation/benchmarks/layer3_fih_dose.py
```

Expected stdout: `[OK] Sanity floor: 12/12 pass.` (Tier B is unchanged and still passes.)

Also inspect the Tier A fold-error counts:

```bash
python3 -c "
import json
d = json.load(open('validation/reports/layer3_fih_dose.json'))
s = d['summary']
print(f'Tier A within-3x: {s[\"gold_within_3x\"]}/{s[\"gold_n\"]} = {100*s[\"gold_within_3x_fraction\"]:.1f}%')
print(f'Tier A within-10x: {s[\"gold_within_10x\"]}/{s[\"gold_n\"]}')
print(f'Tier B pass: {s[\"sanity_pass_count\"]}/{s[\"sanity_n\"]}')
"
```

- [ ] **Step 2: Run the decomposition orchestrator (oral panel)**

```bash
python3 validation/benchmarks/layer3_ivive_decomposition.py
```

Expected: `[OK] Decomposition wrote ... (12 compounds)`.

- [ ] **Step 3: Append Sprint 11 comparison narrative to the FIH report**

Read `validation/reports/layer3_fih_dose.md` — it already has §1 Summary + §2 Tier A table + §3 Tier B table. Append a new section:

```markdown

## Sprint 11 comparison (oral migration)

Sprint 9 (iv_bolus) Tier A within-3x: 5/12 = 41.7% (§8 FAILED).
Sprint 11 (oral)     Tier A within-3x: <M>/12 = <PCT>% (§8 <PASS|FAIL>).

Per-compound fold-error deltas (iv→oral):

| Compound | Sprint 9 fold | Sprint 11 fold | delta |
|---|---:|---:|:---|
| midazolam | 1.73 | <v> | <d> |
| warfarin | 1.04 | <v> | <d> |
| propranolol | 36.30 | <v> | <d> |
| verapamil | 8.14 | <v> | <d> |
| omeprazole | 1.50 | <v> | <d> |
| theophylline | 1.17 | <v> | <d> |
| diclofenac | 11.70 | <v> | <d> |
| diazepam | 6.13 | <v> | <d> |
| metoprolol | 7.25 | <v> | <d> |
| acetaminophen | 1.74 | <v> | <d> |
| lisinopril | 13.36 | <v> | <d> |
| atorvastatin | 70.90 | <v> | <d> |

Fill <v> with the Sprint 11 per-compound `fold_error` from the §2 table. Fill <d> with the ratio Sprint11/Sprint9 (value <1 = improved).

Replace <M>, <PCT>, <PASS|FAIL> with the actual numbers from Step 1.
```

**Important:** You must manually fill in the `<v>` / `<d>` / `<M>` / `<PCT>` / `<PASS|FAIL>` values from the regenerated Sprint 11 report. Do NOT leave placeholders.

- [ ] **Step 4: Append Sprint 11 narrative to the decomposition report**

Read `validation/reports/layer3_ivive_decomposition.md` — already has §1–§6 from Sprint 10. Append:

```markdown

## §7. Sprint 11 — post-oral-migration attribution

After Sprint 11 flipped Tier A simulation route iv_bolus → oral, the decomposition re-runs with the following aggregate attribution:

- liver_model choice: <NEW_LIVER>%
- 1/F route bias: <NEW_ROUTE>% (expected near 0%; simulation route now matches reference)
- residual: <NEW_RESIDUAL>%

The route-bias collapse confirms Sprint 10's diagnostic: the 42% previously attributed to 1/F was a route-comparison artefact. What remains in the residual is the true IVIVE / transporter / mechanistic gap. Compounds with the largest remaining residual:

- atorvastatin: <ATORVA_RESIDUAL>x (OATP1B1 unmodelled — Sprint 12 target)
- lisinopril: <LISI_RESIDUAL>x (non-hepatic elimination + low Peff)
- <others>

Sprint 12 scope is now tighter: OATP1B1 plumbing remains the largest single-compound residual. Sprint 13 (UGT/CYP2C9) addresses diclofenac residual. Remaining residual across the panel defines the fundamental IVIVE prediction ceiling for Charon's current mechanistic coverage.
```

Replace placeholders with actual values from the regenerated report.

- [ ] **Step 5: Commit the regenerated reports + narrative**

```bash
git add validation/reports/layer3_fih_dose.md \
        validation/reports/layer3_fih_dose.json \
        validation/reports/layer3_ivive_decomposition.md \
        validation/reports/layer3_ivive_decomposition.json
git commit -m "chore(sprint11): regenerated Layer 3 reports + Sprint 11 narrative (oral migration)"
```

---

## Task 8: Update Sprint 10 ticket + full test suite green

**Files:**
- Modify: `docs/superpowers/sprint10-ivive-bias-ticket.md`

- [ ] **Step 1: Append Sprint 11 status to the ticket**

Open `docs/superpowers/sprint10-ivive-bias-ticket.md`. Append at the end:

```markdown

## Sprint 11 (oral migration) completed — 2026-04-24

Tier A simulation route migrated iv_bolus → oral after Peff curation for 11 compounds. Benchmark re-run at `validation/reports/layer3_fih_dose.md`.

**Layer 3 Tier A result:** within-3x = <M>/12 = <PCT>% (§8 <PASS|FAIL>).
- Sprint 9 baseline: 5/12 = 41.7%
- Sprint 11 result: <M>/12 = <PCT>%
- Delta: +<DELTA> compounds within 3x

**Post-oral decomposition attribution:**
- liver_model: <NEW_LIVER>%
- route_bias: <NEW_ROUTE>% (collapsed as expected from 42.1% to near 0)
- residual: <NEW_RESIDUAL>% (now the pure signal)

Dominant remaining residuals: atorvastatin (OATP1B1), lisinopril (non-hepatic), diclofenac (UGT/CYP2C9 under-prediction). Sprint 12/13 scope still applies.
```

Fill placeholders from the regenerated Sprint 11 reports.

- [ ] **Step 2: Run the full test suite**

```bash
pytest -q
```

Expected:
- All tests pass
- Total count ~913 (909 baseline + 2 schema tests + 2 decomposition override tests + 1 parametric smoke = 913; actual may vary slightly)
- No new warnings beyond pre-existing `pytest.mark.slow`

If any test fails: STOP. Sprint 11 should not break production paths. A failure is likely either a Peff curation issue (ACAT numerical failure) or a missed schema edge case.

- [ ] **Step 3: Commit the ticket update**

```bash
git add docs/superpowers/sprint10-ivive-bias-ticket.md
git commit -m "docs(sprint11): Sprint 10 ticket — Sprint 11 oral migration results"
```

---

## Self-Review Notes (for the implementer)

- **Peff curation quality gate:** Task 1 is the biggest risk. If a Peff value is wildly wrong (>3x off literature), the oral MRSD will be wrong and Sprint 11's benchmark result will reflect that noise. The implementer should use WebSearch to verify at least the "low confidence" entries (lisinopril, atorvastatin) before committing.

- **Low-Peff numerical stability:** lisinopril Peff ≈ 0.3e-4 cm/s. ACAT's absorption rate constant scales with Peff (`k_abs = 2 × Peff × 3600 / radius`). A low Peff = very slow absorption = potentially long-tail ODE integration. If the smoke test in Task 6 reports a timeout or NaN, the implementer should first check if the low Peff is wrong (literature error) and only after verifying it's correct, investigate the solver.

- **Honest reporting:** If Tier A within-3x does NOT improve, that is the correct outcome to report. Sprint 11's success criteria says ≥50% is strong progress; if we get 4/12 = 33% (worse than Sprint 9's 5/12), then the oral migration introduced more error than it removed, and that is important data. Do not revert. Report honestly per CLAUDE.md §6.5.

- **Tier B sanity floor:** MUST pass 12/12 after Sprint 11. Tier B route is untouched — any regression means the orchestrator accidentally affected shared state or the panel edit was malformed.

- **Integration test parametrisation:** `@pytest.mark.parametrize` generates 12 individual test cases. Each runs a full Pipeline.run() which is ~3-5 seconds including PBPK ODE. Total smoke runtime ~45-60s. This is acceptable for integration tests. Do not try to parallelize.

- **Commit hygiene:** 8 feature commits total (Tasks 1-8) + merge. All commits should be individually reviewable.

- **Estimated final test count:** ~913. If significantly different, investigate — Sprint 11 is not supposed to add or subtract tests beyond the 5 listed here.
