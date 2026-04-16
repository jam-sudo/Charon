# Charon Layer 1 ADMET Benchmark — adme_reference.csv (n=153)

**Generated:** 2026-04-16T10:26:46.167828+00:00
**Panel:** adme_reference_v1

> Evaluated on the conformal calibration set (n=153). Point-estimate metrics are meaningful; 90% CI coverage is tautological by construction (P90 was fit on the same set) and does NOT reflect true out-of-sample performance. CLint excluded: unit mismatch between Charon (hepatocyte uL/min/10^6) and reference (recombinant uL/min/pmol CYP3A4).

## Summary

| property | n | metric_type | mae | rmse | r2 | within_0.5_logP_pct |
| --- | --- | --- | --- | --- | --- | --- |
| logP | 151 | linear | 0.9697 | 1.29 | 0.6567 | 34.44 |
| fu_p | 151 | log | - | - | - | - |
| bp_ratio | 151 | log | - | - | - | - |
| peff_cm_s | 0.000 | - | - | - | - | - |

## Targets

| property | target | met |
| --- | --- | --- |
| logP | MAE < 1.0 | [PASS] |
| fu_p | AAFE < 2.0 | [FAIL] |
| bp_ratio | AAFE < 2.0 | [PASS] |
| peff_cm_s | AAFE < 2.0 | - |

## Results

| name | smiles | logP_obs | logP_pred | logP_err | fup_obs | fup_pred | rbp_obs | rbp_pred | peff_obs | peff_pred |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| midazolam | Clc1ccc2c(c1)C(=NCc1nccn1C)c1ccccc1N2 | 3.89 | 3.809 | -0.0811 | 0.035 | 0.02732 | 0.55 | 0.5623 | 2.800e-04 | - |
| warfarin | CC(=O)CC(c1ccccc1)c1c(O)c2ccccc2oc1=O | 2.7 | 3.61 | 0.9096 | 5.000e-03 | 0.01779 | 0.72 | 0.558 | 1.500e-04 | - |
| propranolol | CC(C)NCC(O)COc1cccc2ccccc12 | 3.48 | 2.578 | -0.9025 | 0.13 | 0.1252 | 0.83 | 0.6232 | 3.200e-04 | - |
| metformin | CN(C)C(=N)NC(=N)N | -1.43 | -1.034 | 0.3958 | 1 | 0.4901 | 0.95 | 0.8367 | 6.000e-06 | - |
| caffeine | Cn1c(=O)c2c(ncn2C)n(C)c1=O | -0.07 | -1.029 | -0.9593 | 0.65 | 0.3039 | 0.72 | 0.6868 | 9.000e-05 | - |
| alprazolam | Cc1nnc2n1-c1ccc(Cl)cc1C(=NC2)c1ccccc1 | 2.12 | 3.58 | 1.46 | 0.2 | 0.03876 | 0.72 | 0.5674 | 2.200e-04 | - |
| diazepam | CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21 | 2.82 | 3.154 | 0.3338 | 0.015 | 0.03658 | 0.68 | 0.5665 | 2.500e-04 | - |
| amlodipine | CCOC(=O)C1=C(COCCN)NC(C)=C(C(=O)OC)C1c1ccccc1Cl | 3 | 2.266 | -0.7337 | 0.07 | 0.1409 | 0.92 | 0.6325 | 1.800e-04 | - |
| atorvastatin | CC(C)c1n(CC[C@@H](O)C[C@@H](O)CC(=O)O)c(-c2ccccc2)c(-c2ccc(F)cc2)c1C(=O)Nc1ccccc1 | 4.06 | 6.314 | 2.254 | 0.02 | 3.899e-03 | 0.85 | 0.5512 | 9.000e-05 | - |
| carbamazepine | NC(=O)N1c2ccccc2C=Cc2ccccc21 | 2.45 | 3.387 | 0.9372 | 0.25 | 0.07438 | 0.72 | 0.5835 | 1.700e-04 | - |
| celecoxib | Cc1ccc(-c2cc(C(F)(F)F)nn2-c2ccc(S(N)(=O)=O)cc2)cc1 | 3.51 | 3.514 | 3.920e-03 | 0.025 | 0.03928 | 0.72 | 0.5677 | 1.600e-04 | - |
| clopidogrel | COC(=O)[C@@H](c1ccccc1Cl)N1CCc2sccc2C1 | 2.21 | 3.674 | 1.464 | 0.02 | 0.07808 | 0.72 | 0.5957 | 1.800e-04 | - |
| codeine | COc1ccc2CC3N(C)CCC4=CC(O)C(c1c24)O3 | 1.19 | 1.728 | 0.5385 | 0.25 | 0.382 | 0.72 | 0.7735 | 1.000e-04 | - |
| dextromethorphan | COc1ccc2CC3N(C)CCC4(CCCC(c1c24)C3)O | 3.37 | 2.801 | -0.5694 | 0.17 | 0.2187 | 0.72 | 0.678 | 2.100e-04 | - |
| erythromycin | CCC1OC(=O)C(C)C(OC2CC(C)(OC)C(O)C(C)O2)C(C)C(OC2OC(C)CC(N(C)C)C2O)C(C)(O)CC(C)C(=O)C(C)C(O)C1(C)O | 3.06 | 1.786 | -1.274 | 0.27 | 0.04105 | 0.82 | 0.574 | 3.000e-05 | - |
| fluconazole | OC(Cn1cncn1)(Cn1cncn1)c1ccc(F)cc1F | 0.5 | 0.7358 | 0.2358 | 0.89 | 0.3316 | 0.75 | 0.6992 | 1.200e-04 | - |
| ibuprofen | CC(C)Cc1ccc(C(C)C(=O)O)cc1 | 3.97 | 3.073 | -0.8968 | 5.000e-03 | 0.1354 | 0.72 | 0.5926 | 2.500e-04 | - |
| imipramine | CN(C)CCCN1c2ccccc2CCc2ccccc21 | 4.28 | 3.875 | -0.405 | 0.11 | 0.1329 | 0.8 | 0.6277 | 2.900e-04 | - |
| ketoconazole | CC(=O)Oc1ccc(C2(Cn3ccnc3)OCC(O2)c2ccc(Cl)cc2Cl)cc1 | 4.35 | 4.756 | 0.4063 | 0.015 | 0.01455 | 0.8 | 0.5565 | 1.100e-04 | - |
| lansoprazole | Fc1cc2c([nH]c(S(=O)Cc3ncc(CC(F)(F)F)cn3)n2)cc1 | 2.72 | 2.905 | 0.1846 | 0.02 | 0.04955 | 0.72 | 0.5723 | 1.900e-04 | - |
| metoprolol | COCCc1ccc(OCC(O)CNC(C)C)cc1 | 1.88 | 1.613 | -0.2668 | 0.88 | 0.6062 | 0.72 | 0.9047 | 1.800e-04 | - |
| nifedipine | COC(=O)C1=C(C)NC(C)=C(C(=O)OC)C1c1ccccc1[N+](=O)[O-] | 2.2 | 2.176 | -0.0244 | 0.045 | 0.07757 | 0.72 | 0.5954 | 2.000e-04 | - |
| omeprazole | COc1ccc2[nH]c(S(=O)Cc3ncc(C)c(OC)c3C)nc2c1 | 2.23 | 2.9 | 0.6697 | 0.05 | 0.03741 | 0.72 | 0.5668 | 1.600e-04 | - |
| simvastatin | CCC(C)(C)C(=O)O[C@H]1C[C@@H](O)C=C2C=C[C@H](C)[C@H](CC[C@@H](O)CC(=O)O)[C@@H]21 | 4.68 | 3.08 | -1.6 | 0.05 | 0.1094 | 0.8 | 0.5844 | 1.200e-04 | - |
| verapamil | COc1ccc(CCN(C)CCCC(C#N)(c2ccc(OC)c(OC)c2)C(C)C)cc1OC | 3.79 | 5.093 | 1.303 | 0.1 | 0.02368 | 0.9 | 0.5639 | 2.500e-04 | - |
| amoxicillin | CC1(C)SC2C(NC(=O)C(N)c3ccc(O)cc3)C(=O)N2C1C(=O)O | -1.85 | 0.0237 | 1.874 | 0.8 | 0.2897 | 0.87 | 0.6804 | 2.000e-04 | - |
| ciprofloxacin | O=C(O)c1cn(C2CC2)c2cc(N3CCNCC3)c(F)cc2c1=O | 0.28 | 1.583 | 1.303 | 0.7 | 0.432 | 0.95 | 0.7444 | 5.000e-05 | - |
| doxycycline | OC1=C(C(N)=O)C(=O)C2(O)C(O)C3C(O)c4c(O)cccc4C(C)(O)C3CC2=C1N(C)C | -0.02 | -0.9294 | -0.9094 | 0.07 | 0.06409 | 0.8 | 0.5875 | 5.000e-05 | - |
| clarithromycin | CCC1OC(=O)C(C)C(OC2CC(C)(OC)C(O)C(C)O2)C(C)C(OC2OC(C)CC(N(C)C)C2O)C(C)(OC)CC(C)C(=O)C(C)C(O)C1(C)O | 1.83 | 2.44 | 0.6097 | 0.4 | 0.05119 | 0.8 | 0.5799 | 5.000e-05 | - |
| azithromycin | CCC1OC(=O)C(C)C(OC2CC(C)(OC)C(O)C(C)O2)C(C)C(OC2OC(C)CC(N(C)C)C2O)C(C)(O)CC(C)CN(C)C(C)C(O)C1(C)O | 2.47 | 1.901 | -0.5693 | 0.51 | 0.05132 | 0.6 | 0.58 | 6.000e-05 | - |
| tetracycline | OC1=C(C(N)=O)C(=O)C2(O)C(O)C3C(O)c4c(O)cccc4C(O)(C)C3CC2=C1N(C)C | -1.3 | -0.9294 | 0.3706 | 0.35 | 0.06409 | 0.8 | 0.5875 | 2.000e-05 | - |
| trimethoprim | COc1cc(Cc2cnc(N)nc2N)cc(OC)c1OC | 0.91 | 1.258 | 0.3476 | 0.55 | 0.2115 | 0.8 | 0.6737 | 1.000e-04 | - |
| sulfamethoxazole | Cc1cc(NS(=O)(=O)c2ccc(N)cc2)no1 | 0.89 | 1.366 | 0.476 | 0.35 | 0.1531 | 0.8 | 0.6396 | 1.000e-04 | - |
| nitrofurantoin | O=C1CN(/N=C/c2ccc(o2)[N+](=O)[O-])C(=O)N1 | -0.47 | 0.0735 | 0.5435 | 0.6 | 0.1542 | 0.8 | 0.6194 | 5.000e-05 | - |
| levofloxacin | CC1COc2c(N3CCN(C)CC3)c(F)cc3c(=O)c(C(=O)O)cn1c23 | -0.39 | 1.544 | 1.934 | 0.69 | 0.3154 | 0.9 | 0.6919 | 1.000e-04 | - |
| lovastatin | CCC(C)C(=O)OC1CC(O)C=C2C=CC(C)C(CCC3CC(O)CC(=O)O3)C21 | 4.26 | 2.92 | -1.34 | 0.05 | 0.08971 | 0.8 | 0.5904 | 1.500e-04 | - |
| pravastatin | CC[C@H](C)C(=O)O[C@H]1C[C@@H](O)C=C2C=C[C@H](C)[C@H](CC[C@@H](O)C[C@@H](O)CC(=O)O)[C@@H]21 | 0.59 | 2.44 | 1.85 | 0.5 | 0.2071 | 0.8 | 0.6152 | 4.000e-05 | - |
| fluvastatin | OC(CC(=O)O)CC=Cc1c(c2ccccc2)c2ccccc2n1C1CC1 | 3.24 | 4.882 | 1.642 | 0.015 | 0.01539 | 0.8 | 0.5548 | 1.200e-04 | - |
| pitavastatin | OC(CC(=O)O)CC=Cc1c(c2ccccc2)c2ccccc2n1C1CC1 | 2.69 | 4.882 | 2.192 | 0.04 | 0.01539 | 0.8 | 0.5548 | 8.000e-05 | - |
| lisinopril | N[C@@H](CCc1ccccc1)C(=O)N1C[C@@H](C(=O)O)C[C@H]1C(=O)O | -3.31 | 0.3328 | 3.643 | 0.75 | 0.3188 | 0.9 | 0.6934 | 5.000e-06 | - |
| enalapril | CCOC(=O)C(CCc1ccccc1)NC(C)C(=O)N1CCCC1C(=O)O | 0.07 | 1.605 | 1.535 | 0.5 | 0.1811 | 0.9 | 0.6315 | 8.000e-05 | - |
| losartan | CCCCc1nc(Cl)c(CO)n1Cc1ccc(-c2ccccc2-c2nnn[nH]2)cc1 | 3.15 | 4.267 | 1.117 | 0.015 | 0.01165 | 0.8 | 0.5552 | 1.200e-04 | - |
| captopril | CC(CS)C(=O)N1CCCC1C(=O)O | 0.34 | 0.6279 | 0.2879 | 0.7 | 0.4326 | 0.9 | 0.6863 | 6.000e-05 | - |
| ramipril | CCOC(=O)C(CCc1ccccc1)NC(C)C(=O)N1C2CCCC2CC1C(=O)O | 1.33 | 2.383 | 1.053 | 0.27 | 0.1832 | 0.8 | 0.6324 | 7.000e-05 | - |
| candesartan | CCOc1nc2cccc(C(=O)O)c2n1Cc1ccc(-c2ccccc2-c2nn[nH]n2)cc1 | 5.23 | 4.029 | -1.201 | 5.000e-03 | 8.095e-03 | 0.8 | 0.5525 | 8.000e-05 | - |
| irbesartan | CCCCC1=NC2(CCCC2)c2ccc(-c3ccc(-c4nn[nH]n4)cc3)cc2N1 | 3.61 | 5.317 | 1.707 | 0.04 | 0.01153 | 0.8 | 0.5567 | 1.000e-04 | - |
| olmesartan | CCCc1nc(C(C)(C)O)c(C(=O)O)n1Cc1ccc(-c2ccccc2-c2nn[nH]n2)cc1 | 3.83 | 3.657 | -0.1734 | 5.000e-03 | 0.01763 | 0.8 | 0.5556 | 7.000e-05 | - |
| telmisartan | CCCc1nc2c(C)cc(-c3ccc(-c4nn[nH]n4)cc3)cc2n1Cc1ccccc1C(=O)O | 7.2 | 4.891 | -2.309 | 5.000e-03 | 8.784e-03 | 0.8 | 0.5528 | 1.000e-04 | - |
| valsartan | CCCCC(=O)N(Cc1ccc(-c2ccccc2-c2nn[nH]n2)cc1)C(C(=O)O)C(C)C | 4.1 | 4.162 | 0.0617 | 0.05 | 0.01935 | 0.8 | 0.5561 | 8.000e-05 | - |
| atenolol | CC(C)NCC(O)COc1ccc(CC(N)=O)cc1 | 0.16 | 0.4521 | 0.2921 | 0.95 | 0.2232 | 0.9 | 0.6806 | 5.000e-05 | - |
| carvedilol | COc1ccccc1OCCNCC(O)COc1cccc2[nH]c3ccccc3c12 | 3.15 | 3.738 | 0.588 | 0.02 | 0.02674 | 0.8 | 0.5656 | 2.000e-04 | - |
| bisoprolol | CC(C)NCC(O)COc1ccc(COCCOC(C)C)cc1 | 1.87 | 2.366 | 0.4959 | 0.7 | 0.3304 | 0.8 | 0.7433 | 1.500e-04 | - |
| labetalol | CC(NCC(O)c1ccc(O)c(C(N)=O)c1)CCc1ccccc1 | 1.24 | 2.135 | 0.8954 | 0.5 | 0.2426 | 0.8 | 0.6919 | 1.500e-04 | - |
| nebivolol | OC(c1ccc2c(c1)C(F)(F)OC2)C1CCN(CC(O)c2ccc3c(c2)OC(F)(F)C3)C1 | 3.01 | 3.883 | 0.8727 | 0.02 | 0.0879 | 0.8 | 0.6014 | 1.800e-04 | - |
| timolol | CC(C)(C)NCC(O)COc1noc2cc(OC3CCCO3)ccc12 | 1.61 | 2.471 | 0.8609 | 0.4 | 0.1352 | 0.8 | 0.6291 | 1.500e-04 | - |
| sotalol | CC(NCC(O)c1ccc(NS(C)(=O)=O)cc1)C | -1.01 | 1.09 | 2.099 | 0.85 | 0.3897 | 0.9 | 0.778 | 5.000e-05 | - |
| nadolol | CC(C)(C)NCC(O)c1ccc2c(c1)C(O)C(O)CC2 | 0.71 | 1.449 | 0.7386 | 0.7 | 0.7219 | 0.9 | 0.9723 | 7.000e-05 | - |
| glipizide | Cc1cnc(C(=O)NCCc2ccc(S(=O)(=O)NC(=O)NC3CCCCC3)cc2)cn1 | 1.45 | 2.078 | 0.6281 | 0.02 | 0.05644 | 0.8 | 0.5754 | 9.000e-05 | - |
| glyburide | COc1ccc(Cl)cc1C(=O)NCCc1ccc(S(=O)(=O)NC(=O)NC2CCCCC2)cc1 | 4.79 | 3.642 | -1.148 | 4.000e-03 | 0.01361 | 0.8 | 0.5561 | 1.000e-04 | - |
| pioglitazone | O=C1NC(=O)SC1Cc1ccc(OCCc2ncccc2)cc1 | 2.95 | 2.597 | -0.3528 | 5.000e-03 | 0.03611 | 0.8 | 0.5663 | 1.500e-04 | - |
| rosiglitazone | CN(CCOc1ccc(CC2SC(=O)NC2=O)cc1)c1ccccn1 | 2.1 | 2.491 | 0.3909 | 5.000e-03 | 0.03754 | 0.8 | 0.572 | 1.200e-04 | - |
| sitagliptin | N#Cc1cc(F)c(F)cc1CC(N)CC(=O)N1CCn2c(nnc2C(F)(F)F)C1 | 0.43 | 1.749 | 1.319 | 0.62 | 0.2969 | 0.85 | 0.7237 | 8.000e-05 | - |
| empagliflozin | OCC1OC(c2cc(Cc3ccc4c(c3)OCC(c3ccc(Cl)cc3)O4)ccc2Cl)C(O)C(O)C1O | 1.29 | 3.612 | 2.322 | 0.17 | 0.03675 | 0.8 | 0.5665 | 1.000e-04 | - |
| canagliflozin | OCC1OC(c2cc(F)c(-c3ccc4c(c3)CCC(Cc3ccc(F)cc3)O4)cc2)C(O)C(O)C1O | 2.09 | 3.083 | 0.993 | 0.01 | 0.05651 | 0.8 | 0.5754 | 1.000e-04 | - |
| dapagliflozin | OCC1OC(c2cc(Cc3ccc(OCC4CCCCO4)cc3)ccc2Cl)C(O)C(O)C1O | 2.6 | 2.394 | -0.2064 | 0.09 | 0.1192 | 0.8 | 0.6036 | 1.000e-04 | - |
| saxagliptin | N#CC(C1CCCC1(O)O)N1C(=O)CC2(CC3CC4CC(C3)CC2C4)C1 | 0.3 | 2.424 | 2.124 | 0.75 | 0.1312 | 0.8 | 0.609 | 8.000e-05 | - |
| linagliptin | Cc1nc(=O)n(/C=C/c2ccccn2)c2c1ncn2C1CC(N2CCC(=O)CC2)C1 | 3.21 | 2.293 | -0.9173 | 0.21 | 0.1013 | 0.8 | 0.6092 | 8.000e-05 | - |
| fluoxetine | CNCCC(Oc1ccc(C(F)(F)F)cc1)c1ccccc1 | 4.05 | 4.435 | 0.385 | 0.06 | 0.02679 | 0.8 | 0.5657 | 2.500e-04 | - |
| sertraline | CNC1CCC(c2ccc(Cl)c(Cl)c2)c2ccccc21 | 4.97 | 5.18 | 0.2096 | 0.015 | 0.01008 | 0.8 | 0.5559 | 2.500e-04 | - |
| paroxetine | Fc1ccc(C2CCNCC2COc2ccc3c(c2)OCO3)cc1 | 3.95 | 3.327 | -0.6235 | 0.04 | 0.145 | 0.8 | 0.6349 | 2.200e-04 | - |
| escitalopram | N#CCC(c1ccc(F)cc1)c1ccc2c(c1)C(CCCN1CCCC1)CO2 | 2.83 | 5.223 | 2.393 | 0.44 | 0.03089 | 0.8 | 0.5681 | 2.000e-04 | - |
| venlafaxine | COc1ccc(C(CN(C)C)C2(O)CCCCC2)cc1 | 3.56 | 3.036 | -0.5244 | 0.73 | 0.3534 | 0.8 | 0.7567 | 2.200e-04 | - |
| duloxetine | CNCC(Oc1cccc2ccccc12)c1cccs1 | 4.07 | 4.241 | 0.1708 | 0.03 | 8.417e-03 | 0.8 | 0.5549 | 2.500e-04 | - |
| bupropion | CC(NC(C)(C)C)C(=O)c1cccc(Cl)c1 | 3.56 | 3.299 | -0.2607 | 0.15 | 0.1507 | 0.8 | 0.6381 | 2.500e-04 | - |
| mirtazapine | CN1CCN2c3ncccc3Cc3ccccc3C2C1 | 2.9 | 2.479 | -0.4211 | 0.15 | 0.295 | 0.8 | 0.7226 | 2.300e-04 | - |
| amitriptyline | CN(C)CCC=C1c2ccccc2CCc2ccccc21 | 4.92 | 4.169 | -0.7514 | 0.05 | 0.1107 | 0.8 | 0.6147 | 2.800e-04 | - |
| nortriptyline | CNCCC=C1c2ccccc2CCc2ccccc21 | 4.43 | 3.826 | -0.6036 | 0.08 | 0.09597 | 0.8 | 0.6061 | 2.600e-04 | - |
| rivaroxaban | O=C1OC[C@@H](n2cc(-c3ccc(N4CCOCC4=O)cc3)nn2)c2ccc(Cl)cc2N1 | 0.38 | 3.113 | 2.733 | 0.07 | 0.06587 | 0.8 | 0.5796 | 1.000e-04 | - |
| apixaban | COc1ccc(-n2nc(C(N)=O)c3c2C(=O)N(c2ccc(N4CCOCC4)cc2)CC3)cc1 | 1.26 | 2.019 | 0.7593 | 0.13 | 0.2287 | 0.8 | 0.6838 | 8.000e-05 | - |
| dabigatran | CN1C(=O)N(CCC(=O)OCC)c2cc(c3ccc(NC(=N)N)cc3)ccc21 | 0.37 | 2.266 | 1.896 | 0.65 | 0.1823 | 0.8 | 0.6566 | 5.000e-05 | - |
| aspirin | CC(=O)Oc1ccccc1C(=O)O | 1.19 | 1.31 | 0.1201 | 0.9 | 0.2091 | 0.9 | 0.6159 | 2.000e-04 | - |
| edoxaban | CC(C)N1CCC(NC(=O)C2CC(NC(=O)c3cnc(Cl)s3)CN2C(=O)C(NC(=O)C(O)C2CCCCC2)c2ccc(Cl)cc2)CC1 | 1.36 | 3.937 | 2.577 | 0.55 | 0.02565 | 0.8 | 0.565 | 7.000e-05 | - |
| naproxen | COc1ccc2cc(C(C)C(=O)O)ccc2c1 | 3.18 | 3.037 | -0.1435 | 3.000e-03 | 0.04802 | 0.72 | 0.5651 | 2.000e-04 | - |
| diclofenac | OC(=O)Cc1ccccc1Nc1c(Cl)cccc1Cl | 4.51 | 4.364 | -0.1459 | 5.000e-03 | 9.920e-03 | 0.72 | 0.5545 | 2.000e-04 | - |
| meloxicam | Cc1cnc(NC(=O)C2=C(O)c3ccccc3S(=O)(=O)N2C)s1 | 0.1 | 1.951 | 1.851 | 6.000e-03 | 0.01366 | 0.72 | 0.5561 | 1.500e-04 | - |
| indomethacin | COc1ccc2c(c1)c(CC(=O)O)c(C)n2C(=O)c1ccc(Cl)cc1 | 4.27 | 3.927 | -0.3427 | 4.000e-03 | 0.01612 | 0.72 | 0.5551 | 1.800e-04 | - |
| acetaminophen | CC(=O)Nc1ccc(O)cc1 | 0.46 | 1.351 | 0.8906 | 0.83 | 0.2719 | 0.9 | 0.6724 | 1.500e-04 | - |
| tramadol | COc1cccc(C2(O)CCCCC2CN(C)C)c1 | 2.77 | 2.635 | -0.1354 | 0.8 | 0.3745 | 0.8 | 0.7691 | 2.000e-04 | - |
| ketorolac | OC(=O)C1CCc2cccc(-c3ccncc3)c2N1 | 1.43 | 2.56 | 1.13 | 6.000e-03 | 0.102 | 0.72 | 0.5959 | 1.200e-04 | - |
| piroxicam | CN1C(C(=O)Nc2ccccn2)=C(O)c2ccccc2S1(=O)=O | 1.86 | 1.581 | -0.279 | 5.000e-03 | 0.09703 | 0.72 | 0.5937 | 1.300e-04 | - |
| phenytoin | O=C1NC(=O)C(c2ccccc2)(c2ccccc2)N1 | 2.47 | 1.77 | -0.7004 | 0.1 | 0.1324 | 0.78 | 0.6096 | 1.500e-04 | - |
| levetiracetam | CCC(C(N)=O)N1CCCC1=O | -0.64 | -0.1273 | 0.5127 | 0.9 | 0.5978 | 0.9 | 0.819 | 5.000e-05 | - |
| valproic_acid | CCCC(CCC)C(=O)O | 2.75 | 2.287 | -0.4626 | 0.07 | 0.5472 | 0.9 | 0.7224 | 1.500e-04 | - |
| lamotrigine | Nc1nnc(-c2cccc(Cl)c2Cl)c(N)n1 | 0.97 | 2.01 | 1.04 | 0.45 | 0.1051 | 0.9 | 0.6115 | 1.300e-04 | - |
| topiramate | OCC1(OS(N)(=O)=O)OC2COC3(CCCCC3C)OC2C1OC(C)(C)C | -0.22 | 0.7994 | 1.019 | 0.85 | 0.3869 | 0.9 | 0.7241 | 8.000e-05 | - |
| clonazepam | O=C1CN=C(c2ccccc2Cl)c2cc([N+](=O)[O-])ccc2N1 | 2.41 | 3.038 | 0.6277 | 0.15 | 0.04126 | 0.72 | 0.5686 | 1.800e-04 | - |
| gabapentin | NCC1(CC(=O)O)CCCCC1 | -1.1 | 1.37 | 2.47 | 0.97 | 0.3956 | 0.9 | 0.728 | 5.000e-05 | - |
| pregabalin | CC(C)C[C@H](CN)CC(=O)O | -1.35 | 1.082 | 2.432 | 1 | 0.6466 | 0.9 | 0.841 | 5.000e-05 | - |
| memantine | CC12CC3CC(N)(C1)CC(C)(C3)C2 | 2.09 | 2.694 | 0.6041 | 0.55 | 0.3112 | 0.9 | 0.7321 | 2.000e-04 | - |
| itraconazole | CC(c1ncn(c1)-c1ccc(N2CCN(c3ccc(OCC4COC(Cn5ccnc5)(c5ccc(Cl)cc5Cl)O4)cc3)CC2)cc1)OC | 5.66 | 6.757 | 1.097 | 5.000e-03 | 0.02077 | 0.8 | 0.5622 | 1.200e-04 | - |
| voriconazole | CC(c1ncncc1F)C(O)(Cn1cncn1)c1ccc(F)cc1 | 1.27 | 2.038 | 0.7678 | 0.42 | 0.3792 | 0.8 | 0.7206 | 1.500e-04 | - |
| posaconazole | OCC(Oc1ccc(-c2ccn(-c3ccc(N4CCN(c5ccc(OCC6COC(c7ccc(F)cc7F)(c7ccc(Cl)cc7Cl)O6)cc5)CC4)cc3)n2)cc1)CC | 3.72 | 9.296 | 5.576 | 0.02 | 6.955e-03 | 0.8 | 0.5541 | 8.000e-05 | - |
| cyclosporine | CCC1NC(=O)C(C(O)C(C)CC=CC)N(C)C(=O)C(C(C)C)N(C)C(=O)C(CC(C)C)N(C)C(=O)C(CC(C)C)N(C)C(=O)C(C)NC(=O)C(C)NC(=O)C(CC(C)C)N(C)C(=O)C(C(C)C)NC(=O)C(CC(C)C)N(C)C(=O)CN(C)C1=O | 2.92 | 3.269 | 0.349 | 0.07 | 0.057 | 2.4 | 0.5756 | 8.000e-05 | - |
| tacrolimus | COC1CC(=O)C=C(C)/C=C(\\C)C[C@@H](CC(=O)[C@H](CC(C=C1OC)=O)O[C@@H]1O[C@@H](C)[C@@H](O)[C@H](OC)C1)[C@@H]1C[C@@H](O)C(\\C=C\\C=C(/C)[C@@H](OC)C(=O)[C@@H](C)[C@@H](O)C(\\C)=C\\[C@@H](CC=O)C(=O)O1)=O | 3.96 | 3.698 | -0.2616 | 0.01 | 0.04015 | 12 | 0.5681 | 9.000e-05 | - |
| sirolimus | COC1CC(O)CC(C)CC(OC)C(=O)C(OC)CC(C=CC=CC=CC(OC2OC(C)CC(O)C2OC)C(C)CC2CC(=O)C(C=C(C)CC(O)CC(C=CC1C)OC)O2)C | 4.33 | 6.42 | 2.09 | 0.08 | 0.01998 | 6.5 | 0.559 | 7.000e-05 | - |
| imatinib | Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1 | 3.74 | 4.59 | 0.8503 | 0.07 | 0.02098 | 0.8 | 0.5623 | 1.500e-04 | - |
| gefitinib | COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1 | 3.16 | 4.276 | 1.116 | 0.1 | 0.03009 | 0.8 | 0.5676 | 1.500e-04 | - |
| erlotinib | COCCOc1cc2ncnc(Nc3cccc(C#C)c3)c2cc1OCCOC | 2.21 | 3.405 | 1.195 | 0.06 | 0.03833 | 0.8 | 0.5724 | 1.300e-04 | - |
| tamoxifen | CCC(=C(c1ccccc1)c1ccccc1)c1ccc(OCCN(C)C)cc1 | 6.3 | 5.996 | -0.3039 | 0.01 | 4.069e-03 | 0.8 | 0.5524 | 2.500e-04 | - |
| ritonavir | CC(C)[C@H](NC(=O)N(C)Cc1csc(C(C)C)n1)C(=O)N[C@@H](C[C@H](O)[C@H](Cc1ccccc1)NC(=O)OCc1cncs1)Cc1ccccc1 | 4.14 | 5.905 | 1.765 | 0.02 | 0.01244 | 0.8 | 0.5556 | 5.000e-05 | - |
| lopinavir | CC(C)c1nc(CN(C)C(=O)NC(CC(O)C(Cc2ccccc2)NC(=O)COc2c(C)cccc2C)Cc2ccccc2)cs1 | 4.89 | 6.194 | 1.304 | 0.02 | 6.747e-03 | 0.8 | 0.553 | 4.000e-05 | - |
| efavirenz | OC1(C#CC2CC2)NC(=O)Oc2ccc(Cl)cc21 | 4.6 | 2 | -2.599 | 5.000e-03 | 0.1447 | 0.8 | 0.6151 | 2.000e-04 | - |
| oseltamivir | CCOC(=O)C1=CC(OC(CC)CC)C(NC(C)=O)C(N)C1 | 0.36 | 1.285 | 0.9254 | 0.58 | 0.4529 | 0.8 | 0.8149 | 1.000e-04 | - |
| tenofovir | Nc1ncnc2c1ncn2[C@@H](CO)OCP(=O)(O)O | -1.02 | -0.9488 | 0.0712 | 0.99 | 0.2524 | 0.9 | 0.6977 | 2.000e-05 | - |
| emtricitabine | Nc1nc(=O)n(C2CSC(CO)O2)cc1F | -1.4 | -0.455 | 0.945 | 0.93 | 0.2993 | 0.9 | 0.7251 | 3.000e-05 | - |
| nirmatrelvir | CC1(C)C2CCC1(CC(=O)NC(CC(=O)C(F)(F)F)C#N)C(NC(=O)C1CC1)C2 | 1.24 | 2.627 | 1.387 | 0.31 | 0.1204 | 0.8 | 0.6042 | 8.000e-05 | - |
| loratadine | CCOC(=O)N1CCC(=C2c3ccc(Cl)cc3CCc3ncccc32)CC1 | 4.02 | 4.888 | 0.8678 | 0.025 | 0.01828 | 0.8 | 0.5582 | 2.200e-04 | - |
| cetirizine | OC(=O)COCCN1CCN(C(c2ccccc2)c2ccc(Cl)cc2)CC1 | 1.7 | 3.148 | 1.448 | 0.06 | 0.0352 | 0.8 | 0.5658 | 8.000e-05 | - |
| fexofenadine | CC(C)(c1ccc(C(O)CCCN2CCC(C(O)(c3ccccc3)c3ccccc3)CC2)cc1)C(=O)O | 2.18 | 5.511 | 3.331 | 0.27 | 0.01964 | 0.8 | 0.5588 | 5.000e-05 | - |
| diphenhydramine | CN(C)CCOC(c1ccccc1)c1ccccc1 | 3.27 | 3.354 | 0.0842 | 0.24 | 0.2352 | 0.8 | 0.6876 | 2.000e-04 | - |
| chlorpheniramine | CN(C)CCC(c1ccc(Cl)cc1)c1ccccn1 | 3.38 | 3.819 | 0.4386 | 0.28 | 0.1419 | 0.8 | 0.633 | 2.000e-04 | - |
| furosemide | NS(=O)(=O)c1cc(C(=O)O)c(NCc2ccco2)cc1Cl | 2.03 | 1.891 | -0.1393 | 0.01 | 0.08767 | 0.72 | 0.5894 | 5.000e-05 | - |
| hydrochlorothiazide | NS(=O)(=O)c1cc2c(cc1Cl)NCNS2(=O)=O | -0.07 | -0.3513 | -0.2813 | 0.58 | 0.2189 | 0.9 | 0.6781 | 5.000e-05 | - |
| spironolactone | CC12CCC(=O)C=C1CCC1C2C(=O)CC2(C)C1CCC2(O)C(=O)CSC | 2.76 | 3.361 | 0.6005 | 0.02 | 0.05065 | 0.8 | 0.5728 | 1.500e-04 | - |
| torsemide | CC(C)NC(=O)NS(=O)(=O)c1cc(C)c(NC2CCCC2)cc1 | 1.24 | 2.746 | 1.506 | 0.015 | 0.08642 | 0.8 | 0.6006 | 8.000e-05 | - |
| chlorthalidone | NS(=O)(=O)c1ccc(C2(O)NC(=O)c3ccccc32)cc1Cl | 0.8 | 0.9242 | 0.1242 | 0.25 | 0.09503 | 0.8 | 0.5928 | 7.000e-05 | - |
| olanzapine | Cc1cc2c(s1)Nc1ccccc1N2C1=NCCN1C | 2.84 | 3.553 | 0.713 | 0.07 | 0.04337 | 0.8 | 0.5754 | 2.300e-04 | - |
| risperidone | Cc1nc2n(CC3CCN(CCc4noc5cc(F)ccc45)CC3)c(=O)c3ccccc3c2[nH]1 | 3.27 | 4.421 | 1.151 | 0.1 | 0.04276 | 0.8 | 0.575 | 2.000e-04 | - |
| quetiapine | OCCOCCN1CCN(C2=Nc3ccccc3Sc3ccccc32)CC1 | 2.81 | 2.856 | 0.046 | 0.17 | 0.04907 | 0.8 | 0.5787 | 2.000e-04 | - |
| aripiprazole | O=C1CCc2ccc(OCCCCN3CCN(c4cccc5c4OCC5)CC3)cc2N1 | 3.56 | 3.487 | -0.0726 | 0.01 | 0.08056 | 0.8 | 0.5971 | 2.000e-04 | - |
| haloperidol | OC1(c2ccc(Cl)cc2)CCN(CCCC(=O)c2ccc(F)cc2)CC1 | 3.58 | 4.426 | 0.8456 | 0.08 | 0.09056 | 0.8 | 0.603 | 2.200e-04 | - |
| clozapine | CN1CCN(C2=Nc3ccccc3Nc3cc(Cl)ccc32)CC1 | 3.23 | 3.723 | 0.4927 | 0.05 | 0.06858 | 0.8 | 0.5901 | 2.300e-04 | - |
| ziprasidone | O=c1cc(-n2ccc(CCN3CCN(c4ccc(Cl)cc4)CC3)c2)[nH]c2ccccc12 | 3.24 | 4.337 | 1.097 | 0.02 | 0.02368 | 0.8 | 0.5639 | 2.200e-04 | - |
| paliperidone | OC1(c2noc3cc(F)ccc23)CCN(CCCC2=NNC(=O)c3ccccc32)CC1 | 1.12 | 3.12 | 2 | 0.26 | 0.09573 | 0.8 | 0.606 | 1.200e-04 | - |
| digoxin | C[C@@H]1O[C@@H](O[C@@H]2C[C@H](O)[C@@H](O[C@@H]3C[C@H](O)[C@@H](O[C@@H]4C[C@H](O)[C@@H](OC5CC(CO)=CC(=O)O5)C(C)O4)C(C)O3)C(C)O2)C[C@H](O)[C@H]1O | 1.26 | -1.306 | -2.566 | 0.75 | 0.2293 | 0.9 | 0.6532 | 8.000e-05 | - |
| colchicine | COc1cc2c(c(OC)c1OC)CC(NC(C)=O)C1CC(=O)C(OC)=CC=C12 | 1.85 | 2.276 | 0.4259 | 0.61 | 0.09814 | 0.8 | 0.5942 | 1.500e-04 | - |
| quinidine | COc1ccc2nccc([C@@H](O)[C@@H]3CC4CCN3C=C4)c2c1 | 3.44 | 2.885 | -0.5553 | 0.2 | 0.1546 | 0.8 | 0.6404 | 2.000e-04 | - |
| phenobarbital | O=C1NC(=O)C(CC)(c2ccccc2)C(=O)N1 | 1.47 | 0.7004 | -0.7696 | 0.55 | 0.2021 | 0.8 | 0.6409 | 1.000e-04 | - |
| dexamethasone | C[C@@H]1C[C@H]2[C@@H]3CCC4=CC(=O)C=C[C@]4(C)[C@@]3(F)[C@@H](O)C[C@]2(C)[C@@]1(O)C(=O)CO | 1.83 | 1.896 | 0.0657 | 0.26 | 0.161 | 0.8 | 0.6224 | 1.500e-04 | - |
| prednisolone | C[C@@]12C=CC(=O)C=C1CC[C@@H]1[C@@H]2[C@@H](O)C[C@@]2(C)[C@H]1CC[C@]2(O)C(=O)CO | 1.57 | 1.558 | -0.0124 | 0.3 | 0.1785 | 0.85 | 0.6303 | 1.500e-04 | - |
| methylprednisolone | CC1CC2C3CCC(O)(C(=O)CO)C3(C)CC(O)C2C2(C)C=CC(=O)CC12 | 1.82 | 1.884 | 0.0635 | 0.24 | 0.1371 | 0.8 | 0.6117 | 1.800e-04 | - |
| fentanyl | CCC(=O)N(c1ccccc1)C1CCN(CCc2ccccc2)CC1 | 4.05 | 4.137 | 0.0867 | 0.16 | 0.08669 | 0.8 | 0.6007 | 3.000e-04 | - |
| morphine | Oc1ccc2CC3N(C)CCC4=CC(O)C(c1c24)O3 | 0.89 | 1.425 | 0.5355 | 0.65 | 0.3484 | 0.8 | 0.7538 | 8.000e-05 | - |
| oxycodone | COC1=CC=C2CC3N(C)CCC4(CCCC(=O)C14)C3(O)O2 | 0.65 | 1.583 | 0.9328 | 0.55 | 0.3392 | 0.8 | 0.7484 | 1.200e-04 | - |
| hydrocodone | COc1ccc2CC3N(C)CCC4(CCCC(=O)c1c24)C3O | 1.43 | 1.921 | 0.4907 | 0.74 | 0.4333 | 0.8 | 0.8035 | 1.500e-04 | - |
| naloxone | OC1C2CC3=CC(=O)CCC3(CC1OC2=O)N1CC2CCC1C2 | 1.46 | 1.195 | -0.2649 | 0.54 | 0.2468 | 0.8 | 0.6944 | 1.500e-04 | - |
| lithium_carbonate | [Li+].[Li+].[O-]C([O-])=O | -5 | -8.439 | -3.439 | 1 | 0.2486 | 1 | 0.6619 | 1.000e-05 | - |
| methadone | CCC(=O)C(CC(C)N(C)C)(c1ccccc1)c1ccccc1 | 5.01 | 4.292 | -0.718 | 0.14 | 0.04086 | 0.8 | 0.5739 | 3.000e-04 | - |
| theophylline | Cn1c(=O)c2[nH]cnc2n(C)c1=O | -0.02 | -1.04 | -1.02 | 0.6 | 0.2845 | 0.72 | 0.678 | 1.000e-04 | - |
| allopurinol | O=c1[nH]cnc2[nH]ncc12 | -1.1 | -0.3538 | 0.7462 | 0.97 | 0.2509 | 0.9 | 0.6629 | 7.000e-05 | - |
| febuxostat | CC(C)COc1cc(C#N)c(-c2cc(C)on2)cc1C(=O)O | 1.69 | 3.255 | 1.565 | 6.000e-03 | 0.03736 | 0.8 | 0.5618 | 1.300e-04 | - |

## Excluded

| name | reason |
| --- | --- |
| rosuvastatin | invalid SMILES |
| buprenorphine | invalid SMILES |

## Notes

- CLint excluded: unit mismatch between Charon hepatocyte (µL/min/10^6) and reference recombinant CYP3A4 (µL/min/pmol).
- peff_cm_s: predicted only when permeability module is active (ACAT pipeline).
- fu_p min-positive guard: obs and pred both > 1e-4. bp_ratio guard: both > 0.1. peff guard: both > 1e-8.
- logP source: RDKit Crippen cLogP (derived). fu_p source: XGBoost ml_ensemble. bp_ratio source: empirical formula.
