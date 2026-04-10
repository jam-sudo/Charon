# Obach 1999 Tier-1 Human IV PK Validation Panel

Loaded by `validation/benchmarks/layer2_human_pk.py`.

**Source:** Obach RS, *Drug Metab Dispos* 27(11):1350-1359, 1999.
"Prediction of human clearance of twenty-nine drugs from hepatic
microsomal intrinsic clearance data".

This directory contains:

- `panel.yaml` — panel metadata and observed PK values (Table 6)
- `compounds/*.yaml` — per-compound `CompoundConfig` files with
  literature-sourced physicochemical / binding / metabolism / renal
  values

Each numeric value in every compound YAML cites its source (primarily
Obach 1999 Table 2 for in vitro values; Obach 1999 Table 6 for in vivo
reference PK). Override values (in compound `distribution.empirical_kp_by_tissue`)
require verified literature citations recorded in the `method` field —
see `docs/superpowers/specs/2026-04-10-sprint3b-kp-override-validation-design.md`
Section 7.2 for the citation protocol.
