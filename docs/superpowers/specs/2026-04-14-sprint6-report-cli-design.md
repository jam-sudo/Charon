# Sprint 6 — Report Generation + CLI

**Date:** 2026-04-14
**Scope:** Collector → ReportData → Markdown/JSON export, and 5-command
CLI entry point (`predict`, `simulate`, `translate`, `recommend`, `report`).
**Prior state:** Sprint 5 complete (678 tests, full SMILES → oral/IV PK →
FIH dose → 90% CI pipeline). `report/` and `cli/` directories exist as
empty scaffolds.
**Expected scale:** ~12 tasks.

---

## 1. Goal

Turn the working Charon pipeline into a shippable tool:

1. **Auto-generated FIH dose rationale report** (regulatory-oriented
   Markdown) with a machine-readable JSON sibling.
2. **Five CLI commands** covering each pipeline entry point, wired via
   `pyproject.toml` `[project.scripts]`.
3. **Zero new runtime dependencies** — argparse + f-string only.
   matplotlib/jinja2 are explicitly out of scope for Phase A.

**Non-goals:**

- HTML or PDF rendering (Phase B)
- Interactive dashboards (Phase C)
- Cp-time plot images (Phase B — JSON carries the raw arrays so users
  can plot themselves)
- Jinja2 template engine (Phase B — existing empty YAML templates are
  kept in place for future migration)
- Brief / short-form report mode (Phase B)

---

## 2. Architecture

```
PipelineResult
      │
      ▼
[1] Collector (collector.py)        ← pure data flattening
      │  PipelineResult → ReportData (frozen dataclass)
      ▼
[2] Narrative (narrative.py)        ← pure rendering
      │  ReportData → Markdown string
      ▼
[3] Export (export.py)              ← file I/O
      │  writes <out>.md and <out>.json
      ▼
CLI (cli/main.py) orchestrates Pipeline → Collector → Narrative → Export
```

Each stage is a separate module with one responsibility, importable and
testable in isolation.

---

## 3. Report Data Model (`report/collector.py`)

### 3.1 `ReportData` dataclass

```python
from dataclasses import dataclass, field

@dataclass(frozen=True)
class ReportData:
    # Identity
    compound_name: str
    smiles: str
    molecular_weight: float | None
    source: str                              # "predicted" / "experimental" / "mixed"
    compound_type: str | None                # acid / base / neutral / zwitterion

    # Layer 1: ADME properties, flattened for tabular rendering
    # Each entry: {"value", "ci_lower", "ci_upper", "unit", "source", "flag"}
    properties: dict[str, dict]

    # Bridge: IVIVE audit (populated from pipeline metadata)
    ivive_summary: dict                      # {clh_L_h, clint_liver_L_h,
                                             #  extraction_ratio, liver_model,
                                             #  fu_b, steps: list[dict]}

    # Layer 2: PK parameters (one flat dict of floats / None)
    pk_params: dict[str, float | None]       # cmax, tmax, auc_0_inf, auc_0_24,
                                             # half_life, cl_apparent, vss,
                                             # bioavailability, fa, fg, fh
    pk_table: list[dict]                     # sampled Cp-time points (see 3.3)
    route: str
    dose_mg: float
    duration_h: float

    # Layer 3: FIH dose projection (None if not requested)
    dose_recommendation: dict | None         # flattened FIHDoseRecommendation

    # Layer 4: Uncertainty (None if not requested)
    uncertainty: dict | None                 # flattened UncertaintyResult

    # Layer 0 warnings + run metadata
    warnings: list[str] = field(default_factory=list)
    metadata: dict = field(default_factory=dict)
    timestamp: str = ""                      # ISO-8601
    charon_version: str = "0.1.0"
```

### 3.2 `collect()` function

```python
def collect(
    result: PipelineResult,
    *,
    warnings: list[str] | None = None,
    timestamp: str | None = None,
) -> ReportData:
```

Responsibilities:

1. **Flatten compound properties** — iterate through
   `result.compound.properties.{physicochemical, permeability, binding,
   metabolism, safety, renal}` and emit a flat `{name: {value, ci_lower,
   ci_upper, unit, source, flag}}` dict.  `PredictedProperty` fields
   map directly; `None` values are skipped.
2. **Build IVIVE summary** — pull `clint_liver_L_h`, `cl_renal_L_h`,
   `fu_b`, `liver_model`, `compound_type` directly from
   `result.metadata` (these are written by the existing pipeline).  The
   collector surfaces these fields verbatim and does **not** attempt to
   recompute CLh or extraction ratio — those would require re-running
   the liver model and risk drifting from the PBPK ODE's behavior.
   Narrative section 4 uses `cl_apparent` from `pk_params` as the in-vivo
   clearance estimate.
3. **Flatten PK parameters** — one pass through `PKParameters` fields
   into a plain dict.
4. **Build PK table** — sample `time_h` / `cp_plasma` / `cp_blood` at
   canonical timepoints (see 3.3).
5. **Flatten dose recommendation** — if `result.dose_recommendation` is
   not None, convert the Pydantic model + nested HED/MABEL/PAD results
   to plain dicts (use `model_dump()`).
6. **Flatten uncertainty** — if `result.uncertainty` is not None,
   convert the dataclass to a plain dict (use `asdict()`).
7. **Timestamp** — default to `datetime.now(timezone.utc).isoformat()`
   if the caller does not pass one.

The collector performs no numerical computation beyond indexing and
copying.  Any derived metric must be computed upstream in the pipeline.

### 3.3 PK table sampling

Canonical timepoints (hours):
`[0.0, 0.25, 0.5, 1.0, 2.0, 4.0, 6.0, 8.0, 12.0, 24.0, 48.0, 72.0]`

For each canonical time ≤ `max(time_h)`, find the nearest index via
`np.searchsorted` and emit:

```python
{"time_h": t, "cp_plasma_ug_L": float(cp_plasma[i]),
 "cp_blood_ug_L":  float(cp_blood[i])}
```

Rationale: fixed canonical points avoid accidental information loss
near Cmax while keeping the table short and human-readable.  The full
time/Cp arrays are preserved in the JSON sibling via `pk_table_full`
when the caller passes `include_full_profile=True`.

---

## 4. Narrative Rendering (`report/narrative.py`)

### 4.1 Public API

```python
def render_report(data: ReportData) -> str:
    """ReportData → complete Markdown document string."""
```

### 4.2 Section-level helpers (private but individually tested)

Each returns a Markdown fragment ending with a trailing newline.

```python
_render_header(data)               # title + metadata line
_render_executive_summary(data)    # Section 1
_render_compound_profile(data)     # Section 2
_render_adme_table(data)           # Section 3
_render_ivive_audit(data)          # Section 4
_render_pk_results(data)           # Section 5
_render_dose_projection(data)      # Section 6
_render_uncertainty(data)          # Section 7
_render_limitations(data)          # Section 8
_render_appendix(data)             # Section 9
```

`render_report` concatenates them in order with `"\n\n"` separators.

### 4.3 Section contents

**Section 1 — Executive Summary**
- One-line conclusion: `"Recommended FIH starting dose: {mrsd} mg ({route})"`.
- If uncertainty was run: append `"with 90% CI [{lo} – {hi}] mg, confidence: {HIGH|MEDIUM|LOW}"`.
- Limiting method (HED / MABEL / PAD).
- Two-sentence paragraph summarising compound identity and main caveat
  (if any critical warning was raised).

**Section 2 — Compound Profile**
- Bullet list: SMILES, MW, salt form (if any, with salt factor),
  source (predicted/experimental/mixed), compound type.

**Section 3 — ADME Predictions**
- Markdown table: `Property | Value | 90% CI | Unit | Source | Flag`.
- Rows rendered in a fixed order:
  `logp, pka_acid, pka_base, fu_p, fu_inc, bp_ratio, clint_uL_min_mg,
   papp_nm_s, peff_cm_s, herg_ic50_uM`.
- Missing entries are skipped (no blank rows).
- CI column shows `"-"` when `ci_lower`/`ci_upper` are None.
- Flag column shows `"-"` when empty.

**Section 4 — IVIVE & Hepatic Clearance**
- Narrative: "CLint was scaled via the {liver_model} model. Whole-liver
  CLint = {clint_liver_L_h} L/h. fu_b = {fu_b}. Renal CL =
  {cl_renal_L_h} L/h. In-vivo apparent CL from simulation = {cl_apparent} L/h."
- Bullet list of those values pulled from `ivive_summary` and `pk_params`.
- No ConversionStep audit dump in the Markdown; the JSON sibling carries
  the full metadata dict (which includes whatever audit fields the
  pipeline already persisted).

**Section 5 — PK Simulation Results**
- Preamble: route, dose, duration, species.
- PK parameters table: `Parameter | Value | Unit`
  (Cmax ug/L, Tmax h, AUC0-inf ug·h/L, AUC0-24 ug·h/L, t1/2 h,
   CL apparent L/h, Vss L, F, Fa, Fg, Fh).
- Concentration-time table from `data.pk_table`: `Time (h) | Cp_plasma (ug/L) | Cp_blood (ug/L)`.

**Section 6 — FIH Dose Projection**
- Table: `Method | MRSD (mg) | Inputs | Status`.
  HED / MABEL / PAD rows always appear; missing methods show
  `"insufficient inputs"` in Status.
- Limiting method highlighted with `**bold**` MRSD value.
- Footer: `"Safety factor: {sf}  |  Salt factor: {salt_factor}"`.
- Full rationale text (already produced by `project_fih_dose`) printed
  in a fenced code block for clarity.

**Section 7 — Uncertainty Analysis**
- Skipped entirely (header omitted) if `data.uncertainty is None`.
- Summary line: `"Dose: {point} mg [{lo} – {hi}] (90% CI)  —  {confidence}"`.
- Sensitivity table: `Parameter | Importance (%) | Rank`, sorted by
  descending importance, percentages rounded to 1 decimal.
- R² diagnostic: `"SRC regression R² = {r2:.3f}"` plus a one-line warning
  if `r2 < 0.7`.
- Recommendation string (from `UncertaintyResult.recommendation`) in a
  blockquote.
- Sample counts (n_samples, n_successful, convergence_met).

**Section 8 — Limitations & Caveats**
- Fixed boilerplate paragraph noting: well-stirred liver (unless other),
  R&R Kp caveat for weak bases, IR formulation assumption, conformal
  marginal coverage, MRSD uses apparent PK.
- Appended items pulled from `data.warnings` and any `flag` strings
  attached to ADME properties (e.g. `"clint_tier2_ml"`).

**Section 9 — Appendix**
- Charon version, timestamp, seed (if in metadata), liver model,
  solver method, solver nfev.
- Full metadata dict rendered as a fenced `yaml` block (stable key order).

### 4.4 Formatting conventions

- Numeric values rendered with `format_value()` helper:
  - `abs(v) >= 1000` or `abs(v) < 0.01` → `f"{v:.3g}"` (scientific).
  - otherwise → `f"{v:.4g}"` (fixed-ish).
  - `None` → `"-"`.
- All tables use pipe syntax with header separator row.
- Section headings use `##`, subsections `###`.
- Never embed raw Python `repr()` in the document.

---

## 5. Export (`report/export.py`)

### 5.1 API

```python
def export_markdown(data: ReportData, path: Path) -> Path:
    """Render and write Markdown.  Returns resolved absolute path."""

def export_json(
    data: ReportData,
    path: Path,
    *,
    include_full_profile: bool = False,
    full_profile: dict | None = None,
) -> Path:
    """Write ReportData as JSON.  Optionally embed the full time/Cp arrays."""

def export_report(
    data: ReportData,
    output: Path,
    *,
    full_profile: dict | None = None,
) -> tuple[Path, Path]:
    """Produce sibling .md and .json.  If *output* ends in .md, the JSON
    path replaces the suffix with .json; otherwise .md/.json are both
    appended.  Returns (md_path, json_path)."""
```

### 5.2 JSON schema

The JSON is a direct `dataclasses.asdict(data)` with two extensions:

- `numpy.ndarray` and `np.floating` are converted to Python lists / floats
  via a small helper.
- If `include_full_profile=True`, a `full_profile` key is attached with
  `{"time_h": [...], "cp_plasma_ug_L": [...], "cp_blood_ug_L": [...]}`.

Indentation: `json.dumps(..., indent=2, sort_keys=False, default=_json_default)`.

### 5.3 `figures.py` and `templates/` — unchanged

Both stay as empty scaffolds so their names remain reserved for
Phase B (matplotlib plots) without touching the import surface now.

---

## 6. CLI (`cli/main.py`)

### 6.1 Entry point

Add to `pyproject.toml`:

```toml
[project.scripts]
charon = "charon.cli.main:main"
```

`main()` parses `sys.argv`, dispatches to one of five subcommands, and
returns an int exit code (`0` on success, `1` on user error, `2` on
internal failure).

### 6.2 Subcommands

All commands share these top-level options:

- `--compound-name NAME` (default: SMILES)
- `--species {human,rat,dog,monkey}` (default: `human`)
- `--liver-model {well_stirred,parallel_tube,dispersion}` (default:
  `well_stirred`)
- `--json` — emit result as a JSON dict to stdout instead of a table
- `--quiet` / `-q` — suppress the banner line
- `-h` / `--help`

#### 6.2.1 `charon predict <SMILES>`

- Runs Layer 0 + Layer 1 only (no ODE).
- Calls `predict_properties(smiles)` directly, not `Pipeline.from_smiles`.
- Default stdout: a Markdown-style ADME table (re-uses
  `_render_adme_table` on a minimal `ReportData`).
- `--json`: dumps `CompoundProperties.model_dump()`.

#### 6.2.2 `charon simulate <SMILES>`

- Required: `--route {iv_bolus,iv_infusion,oral}`, `--dose MG`.
- Optional: `--duration HOURS` (default 72), `--infusion-duration HOURS`
  (only meaningful with `iv_infusion`).
- Runs Pipeline through Layer 2.  No dose projection, no uncertainty.
- Default stdout: PK parameters table + PK sampling table (from the
  collector).

#### 6.2.3 `charon translate <SMILES>`

- Required: one or more of
  `--noael MG_KG --noael-species SPECIES`,
  `--target-kd NM`,
  `--target-ceff NM`.
- Optional: `--safety-factor X` (default 10), `--tau HOURS` (default 24).
- Runs Pipeline + dose projection (Layer 3).  No uncertainty.
- Default stdout: the existing `rationale` text from
  `FIHDoseRecommendation`, followed by an HED/MABEL/PAD comparison table.

#### 6.2.4 `charon recommend <SMILES>`

- Same inputs as `translate`.
- Adds `--uncertainty` flag; when present, enables Layer 4 with
  `--n-samples N` (default 500) and `--seed N` (default 42).
- Default stdout: headline dose + CI line, then translate-style table.

#### 6.2.5 `charon report <SMILES>`

- Same inputs as `recommend`.
- Additional required: `--output PATH` (file base; `.md` added if absent).
- Optional: `--include-full-profile` to embed full time/Cp arrays in JSON.
- Runs full pipeline (UQ off by default; `--uncertainty` opts in).
- Writes `<out>.md` and `<out>.json`.  Prints the two paths to stdout
  and exits 0.

### 6.3 Shared plumbing

Private helpers in `cli/main.py`:

```python
def _build_pipeline(args, *, dose_projection=None, uncertainty=None) -> Pipeline:
    """Common factory: Pipeline.from_smiles(...) with shared options."""

def _print_table(rows: list[dict], headers: list[str]) -> None:
    """Render a list of dicts as a Markdown-style table on stdout."""

def _fail(msg: str, code: int = 1) -> NoReturn:
    """Print error to stderr, sys.exit(code).  No traceback."""
```

### 6.4 Error handling

- Unparseable SMILES → `_fail("Invalid SMILES: ...", code=1)`.
- Missing required inputs for a method → `_fail(...)` with a hint
  about which flag to pass.
- Unexpected internal exception → print the exception class + message
  (not the traceback) and exit code 2.  Traceback is shown only when
  `--debug` is passed.

---

## 7. File Changes

### New files

| File | Purpose |
|------|---------|
| `src/charon/report/collector.py` | `ReportData` dataclass + `collect()` |
| `src/charon/report/narrative.py` | `render_report()` + section renderers |
| `src/charon/report/export.py` | `export_markdown/json/report` |
| `src/charon/report/__init__.py` | Public symbol exports |
| `src/charon/cli/main.py` | argparse dispatcher + 5 subcommand handlers |
| `src/charon/cli/__init__.py` | Re-export `main` |
| `tests/unit/test_collector.py` | ≥8 tests: flattening, IVIVE, PK table |
| `tests/unit/test_narrative.py` | ≥10 tests: per-section rendering |
| `tests/unit/test_export.py` | ≥5 tests: file write + JSON schema |
| `tests/unit/test_cli.py` | ≥12 tests: one per subcommand + error paths |
| `tests/integration/test_report_e2e.py` | ≥3 tests: full SMILES → MD+JSON |

### Modified files

| File | Changes |
|------|---------|
| `pyproject.toml` | Add `[project.scripts] charon = "charon.cli.main:main"` |
| `src/charon/report/figures.py` | Stays empty (`__all__ = []` only) |
| `src/charon/report/templates/*.yaml` | Unchanged (empty scaffolds) |

---

## 8. Data Flow Example

```
$ charon report "Cc1ncc(n1C)[C@@H](c2ccccc2)OC(=O)N" \
    --route oral --dose 5 \
    --noael 50 --noael-species rat \
    --uncertainty --n-samples 200 \
    --output midazolam_fih

[1/4] Running pipeline...
[2/4] Running uncertainty propagation (200 samples)...
[3/4] Collecting report data...
[4/4] Writing report...

Wrote:
  midazolam_fih.md   (regulatory Markdown)
  midazolam_fih.json (machine-readable data)
```

---

## 9. Validation

### 9.1 Acceptance criteria

| Criterion | Target |
|-----------|--------|
| `charon --help` lists 5 subcommands | Yes |
| Each subcommand runs end-to-end on a valid SMILES | Smoke tested |
| Report MD renders without raw `None` / `repr()` text | Regex check |
| Report JSON is valid JSON and round-trips via `json.loads` | Parser test |
| Report JSON contains every ReportData field | Key set equality |
| PK table has at most 12 rows and is monotonically increasing in time | Assertion |
| CLI `--json` on every subcommand produces valid JSON to stdout | json.loads |
| Invalid SMILES exits 1 with a clean message (no traceback) | subprocess test |
| `charon report` writes both `.md` and `.json` when given `--output foo.md` | File existence |
| All 678 existing tests still pass | Regression |
| New tests: ≥35 | Coverage target |

### 9.2 End-to-end smoke test

`tests/integration/test_report_e2e.py` runs:

1. Midazolam oral 5 mg + NOAEL 50 mg/kg rat + `--uncertainty --n-samples 50`.
2. Verifies MD file contains sections 1–9, dose recommendation block,
   sensitivity table.
3. Verifies JSON parses and has `uncertainty.ci_90_lower_mg > 0`.

### 9.3 CLI regression gate

`tests/unit/test_cli.py` invokes `main([...])` directly (not via
subprocess) to stay fast and capture stdout/stderr via `capsys`.

---

## 10. Known Limitations

### 10.1 No plots in Phase A

The report carries raw time/Cp arrays in JSON so external tooling can
plot.  Adding matplotlib would pull in a large transitive stack for
marginal benefit.

### 10.2 Fixed template

Section layout is hard-coded in `narrative.py`.  Users cannot reorder
or customise sections without editing Python.  Jinja2 migration is
Phase B work and the empty `templates/*.yaml` files are the anchor
point.

### 10.3 Regulator language is author-reviewed

The boilerplate in Section 8 (Limitations) is static text written by
the developer.  It is not legal advice and must be reviewed by the
sponsor before IND submission.  The report is an engineering output,
not a regulatory filing.

### 10.4 No localisation

All text is in English.  Non-English sponsors would need to translate
the rendered Markdown manually.

---

## 11. Definition of Done

- [ ] `ReportData` dataclass and `collect()` implemented and tested.
- [ ] `render_report()` produces a well-formed Markdown document with all
      nine sections.
- [ ] `export_report()` writes `.md` and `.json` siblings.
- [ ] `pyproject.toml` has `[project.scripts] charon = ...`.
- [ ] All five subcommands run end-to-end on at least one reference
      compound each.
- [ ] `charon report` smoke test passes end-to-end with uncertainty.
- [ ] All existing tests still pass (no regression).
- [ ] ≥35 new tests added across collector / narrative / export / cli /
      e2e.
- [ ] `figures.py` and `templates/*.yaml` remain empty scaffolds.
