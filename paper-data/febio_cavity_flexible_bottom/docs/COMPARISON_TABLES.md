# Comparison tables (generated from `results/*.csv` by `scripts/build_tables.py`; see `FINDINGS_REPORT.md` for interpretation)

All values are steady point-A top-surface u_y (m) and vertical interface force F_y (N per m depth, pressure contribution for FEBio; total pressure+viscous for the references where they report it). FEBio values are conditional estimates (all three runs error-terminated; M3 not steady to 1e-4). Reference classes are **not** merged into one ranking. The manuscript values live in the manuscript repository (main 863b8ba); PR #22 / PR #12 values are provisional (both PRs open).

## 1. FEBio vs same-dx solids4foam (StVK-confirmed class)

| level | dx | u_y FEBio | u_y s4f partitioned | rel | F_y FEBio (F1, latest accepted) | F_y s4f QM #618 | rel | u_x FEBio | u_x s4f QM | ratio |
|---|---|---|---|---|---|---|---|---|---|---|
| M1 | 0.1 | -0.211967 | -0.204568 | +3.6 % | -4.7176 | -5.1226 | -7.9 % | 3.57e-04 | 7.45e-05 | 4.8 |
| M2 | 0.05 | -0.226233 | -0.231445 | -2.3 % | -4.9339 | -5.2266 | -5.6 % | 3.35e-04 | 9.08e-05 | 3.7 |
| M3 | 0.025 | -0.230038 | -0.238404 | -3.5 % | -4.9835 | -5.2504 | -5.1 % | 3.26e-04 | 9.54e-05 | 3.4 |

Extrapolations (3 levels, indicative): FEBio u_y -0.2314 m, F_y about -5.00 N/m; solids4foam u_y -0.2407 (partitioned) / -0.2406 (QM), F_y about -5.26 (manuscript). Gaps -3.9 % and -4.9 %.

## 2. Reference values used (`results/reference_values.csv`)

| source | variant | dx | u_y(A) | u_x(A) | F_y | status |
|---|---|---|---|---|---|---|
| solids4foam | partitioned (Aitken), strict | 0.1 | -0.20457 | 7.46e-05 | -5.1528 | reference (StVK) |
| solids4foam | partitioned (Aitken), strict | 0.05 | -0.23145 | 9.04e-05 | -5.2388 | reference (StVK) |
| solids4foam | partitioned (Aitken), strict | 0.025 | -0.23840 | 9.47e-05 | -5.2566 | reference (StVK) |
| solids4foam | partitioned (Aitken), strict | 0.0125 | -0.24011 | 9.59e-05 | -5.2599 | reference (StVK) |
| solids4foam | quasi-monolithic, #618 rerun (manuscript data) | 0.1 | -0.20360 | 7.45e-05 | -5.1226 | reference (StVK) |
| solids4foam | quasi-monolithic, #618 rerun (manuscript data) | 0.05 | -0.23104 | 9.08e-05 | -5.2266 | reference (StVK) |
| solids4foam | quasi-monolithic, #618 rerun (manuscript data) | 0.025 | -0.23820 | 9.54e-05 | -5.2504 | reference (StVK) |
| solids4foam | quasi-monolithic, #618 rerun (manuscript data) | 0.0125 | -0.23999 | 9.65e-05 | -5.2560 | reference (StVK) |
| Tukovic et al. (2018) |  | 0.1 | -0.23324 |  | -5.1603 | reference, constitutive law unresolved / unconfirmed (not StVK-confirmed) |
| Tukovic et al. (2018) |  | 0.05 | -0.24571 |  | -5.2290 | reference, constitutive law unresolved / unconfirmed (not StVK-confirmed) |
| Tukovic et al. (2018) |  | 0.025 | -0.24955 |  | -5.2494 | reference, constitutive law unresolved / unconfirmed (not StVK-confirmed) |
| Tukovic et al. (2018) |  | 0.0125 | -0.25062 |  | -5.2548 | reference, constitutive law unresolved / unconfirmed (not StVK-confirmed) |
| ANSYS |  | 0.1 | -0.24727 | 7.98e-05 | -5.1440 | reference, constitutive law unresolved / unconfirmed (not StVK-confirmed) |
| ANSYS |  | 0.05 | -0.25054 | 8.08e-05 | -5.2180 | reference, constitutive law unresolved / unconfirmed (not StVK-confirmed) |
| ANSYS |  | 0.025 | -0.25197 | 8.08e-05 | -5.2500 | reference, constitutive law unresolved / unconfirmed (not StVK-confirmed) |
| ANSYS |  | 0.0125 | -0.25274 | 4.13e-05 | -5.2630 | reference, constitutive law unresolved / unconfirmed (not StVK-confirmed) |
| COMSOL |  | 0.1 | -0.23600 | 2.19e-04 | -5.1088 | reference, constitutive law unresolved / unconfirmed (not StVK-confirmed) |
| COMSOL |  | 0.05 | -0.24100 | 1.70e-04 | -5.2276 | reference, constitutive law unresolved / unconfirmed (not StVK-confirmed) |
| COMSOL |  | 0.025 | -0.24300 | 1.47e-04 | -5.2956 | reference, constitutive law unresolved / unconfirmed (not StVK-confirmed) |
| COMSOL |  | 0.0125 | -0.24300 | 1.17e-04 | -5.2952 | reference, constitutive law unresolved / unconfirmed (not StVK-confirmed) |

Tukovic et al. values were read from a published vector figure (force rescaled x20 for 0.05 m depth, confirmed by the author per PR #22). ANSYS used a 10 mm slab (x100). COMSOL: 3 significant figures for u_y, solid model and force-integration convention unconfirmed. Extrapolations (provisional, PR #22): Tukovic u_y -0.2510, F_y -5.257; ANSYS u_y about -0.2536 (order about 1), F_y -5.27; COMSOL u_y plateau -0.243 +/- 0.001, F_y -5.295 (+0.7 %), not extrapolable.

## 3. The two FEBio force definitions at the same stored times (`results/force_reconstruction_F1_F4.csv`)

| state | t (s) | flag | F1 | F2 (= F4) | F3 | F1-F2 | D (Pa) | F2 vs QM |
|---|---|---|---|---|---|---|---|---|
| M1 | 900 | 0 | -4.7176 | -5.4375 | -5.4372 | 0.7198 | 0.3647 | +6.15 % |
| M2 | 750 | 0 | -4.9339 | -5.2748 | -5.2738 | 0.3409 | 0.1730 | +0.92 % |
| M3 | 500 | 0 | -4.9809 | -5.2469 | -5.2463 | 0.2660 | 0.1344 | -0.07 % |
| M3_t717p5_failure_exit_state | 717.5 | 2 | -4.9835 | -5.2494 | -5.2488 | 0.2659 | 0.1344 | -0.02 % |

flag 2 = state written at the solver-failure exit (not a converged step; equal to the last accepted step t=715 to 3e-5). Pressure contribution only. F2 is not established as the physical force; F2 = F4 is an algebraic identity for planar facets.

## 4. Coupled summary (`results/coupled_febio_summary.csv`)

| level | run | job | last accepted t | u_y | u_x | F_y F1 | odd-even bracket on F1 | p mean (Pa) | mass imbalance | cores x wall | CPU-h | terminated normally |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| M1 | C2_fsi_pilot_f1 | 14016 | 972.5 | -0.211967 | 3.57e-04 | -4.7176 | [-4.575, -4.860] | 1.1784 | 1.2e-04 | 4 x 15:14.43 | 1.02 | False |
| M2 | C3_fsi_M2 | 14019 | 870 | -0.226233 | 3.35e-04 | -4.9339 | [-4.952, -4.915] | 1.2320 | -1.9e-05 | 4 x 1:06:16 | 4.42 | False |
| M3 | C3_fsi_M4 | 14021 | 715 | -0.230038 | 3.26e-04 | -4.9835 | [-5.033, -4.934] | 1.2447 | -2.3e-05 | 8 x 3:31:29 | 28.20 | False |

See also `results/mesh_convergence_table.csv` (differences, observed orders, Richardson), `results/numerical_sensitivity_C4.csv`, `results/row_parity_by_level.csv`, `results/solver_termination_and_cost.csv`.

