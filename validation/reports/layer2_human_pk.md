# Charon Layer 2 Human PBPK Benchmark — Obach 1999 Tier-1 Panel

**Generated:** 2026-04-16T10:26:51.169601+00:00
**Panel:** obach_1999_tier1

## Summary

| metric | AAFE | within_2_fold | within_3_fold |
| --- | --- | --- | --- |
| CL (L/h) | 3.992 | 0.3333 | 0.4167 |
| Vss (L) | 3.028 | 0.4167 | 0.5 |
| t_half (h) | 9.389 | 0.1667 | 0.1667 |

## Targets

| metric | target | met |
| --- | --- | --- |
| CL (L/h) | AAFE < 2.5 | [FAIL] |
| Vss (L) | AAFE < 3.0 | [FAIL] |

## Results

| compound | cl_pred | cl_obs | cl_fold | cl_pass_2x | vss_pred | vss_obs | vss_fold | vss_pass_2x | t_half_pred | t_half_obs | t_half_fold | t_half_pass_2x | strict_target |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| theophylline | 3.263 | 2.9 | 1.125 | [PASS] | 37.09 | 35 | 1.06 | [PASS] | 9.488 | 8 | 1.186 | [PASS] | [PASS] |
| antipyrine | 2.642 | 2.8 | 1.06 | [PASS] | 55.01 | 43 | 1.279 | [PASS] | 16.5 | 11 | 1.5 | [PASS] | [FAIL] |
| caffeine | 2.142 | 6 | 2.802 | [FAIL] | 43.21 | 42 | 1.029 | [PASS] | 15.35 | 4.9 | 3.132 | [FAIL] | [FAIL] |
| warfarin | 0.7808 | 0.19 | 4.109 | [FAIL] | 278 | 11 | 25.28 | [FAIL] | 274.8 | 37 | 7.426 | [FAIL] | [FAIL] |
| diclofenac | 1.142 | 16 | 14.02 | [FAIL] | 112.6 | 13 | 8.658 | [FAIL] | 79.52 | 1.2 | 66.27 | [FAIL] | [FAIL] |
| midazolam | 5.77 | 21 | 3.639 | [FAIL] | 110.3 | 66 | 1.671 | [PASS] | 24.9 | 3 | 8.301 | [FAIL] | [FAIL] |
| propranolol | 4.717 | 50 | 10.6 | [FAIL] | 133.9 | 270 | 2.017 | [FAIL] | 28.09 | 3.9 | 7.204 | [FAIL] | [FAIL] |
| diazepam | 0.02586 | 1.6 | 61.87 | [FAIL] | 564.8 | 77 | 7.336 | [FAIL] | 1.519e+04 | 43 | 353.2 | [FAIL] | [FAIL] |
| metoprolol | 17.54 | 63 | 3.592 | [FAIL] | 86.03 | 290 | 3.371 | [FAIL] | 22.65 | 4 | 5.663 | [FAIL] | [FAIL] |
| verapamil | 12.27 | 60 | 4.892 | [FAIL] | 101.8 | 350 | 3.44 | [FAIL] | 47.13 | 4 | 11.78 | [FAIL] | [FAIL] |
| omeprazole | 19.25 | 34.2 | 1.776 | [PASS] | 139.2 | 24.5 | 5.683 | [FAIL] | 12.37 | 0.9 | 13.75 | [FAIL] | [FAIL] |
| dextromethorphan | 43.79 | 50 | 1.142 | [PASS] | 298.8 | 250 | 1.195 | [PASS] | 30.9 | 3.5 | 8.828 | [FAIL] | [FAIL] |

## Notes

- Panel: n=12 compounds (Obach 1999 Tier-1)
- Mode: R&R + empirical Kp overrides (with_override)
- Strict-gate failures: 0
- AAFE targets: CL < 2.5, Vss < 3.0 (ARCHITECTURE.md §8)
- BDF solver (scipy.integrate.solve_ivp, method='BDF')
- fu_p applied only via fu_b = fu_p/BP in liver model (no double-application)
