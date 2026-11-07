# Results: CC-S1 (static wall), CC-F1 (fixed-shape fluid), D1 (20-step coupled diagnostic)

Stock FEBio v4.13 (binary sha256 `f5ae1ed6f347f875c392f93d847f2369f437345a9b785296d46b0bd6fbb72f02`). Numbers below were recomputed from the stored raw outputs by an independent script (`scripts/reverify.py`, 141 quantities, all within stated tolerances; corrections to the original report are listed in `provenance/INCONSISTENCIES.md`). Kinds of quantity are kept apart: **[S]** direct solver output (logged node data); **[R]** reconstruction from solver output (F2 element mean, free-body load, volume balance); **[A]** analytical estimate; **[H]** hypothesis. The original report is preserved in `archive/PILOT_RESULTS_original.md`.

## 1. CC-S1 static clamped wall, solid only (jobs 14214-14217; 2 cores; 20 load steps; hex8 SVK, plane strain, follower 200 Pa)
Reference: static clamped midpoint -0.14362832 m, quarter points -0.09645925 m (oomph-lib Kirchhoff-Love beam, `reference/collapsibleChannel_verification_references.json`). Midpoint = interface node at x = 10, quarter = x = 7.5 [S]. Table: `results/s1_static_wall.csv`; figure `figures/fig1_s1_static_wall_convergence.png`.

| wall mesh (cols x rows) | aspect | midpoint (m) | midpoint error | quarter-point error | Newton total | gate | result |
|---|---|---|---|---|---|---|---|
| 80 x 8 (tutorial) | 20 | -0.126920 | -11.63 % | -19.46 % | 91 | 3 % | FAIL |
| 160 x 4 | 5 | -0.139242 | -3.054 % | -6.27 % | 92 | 3 % | FAIL by 0.054 pp (strict) |
| 160 x 8 | 5 | -0.139320 | -2.9996 % (reported as -3.000 %) | -6.17 % | 92 | informative | - |
| 320 x 8 | 5 | -0.142542 | -0.757 % | -1.71 % | 92 | 1 % | PASS |

- Error ratios 3.88 (80 -> 160, 8 rows) and 3.96 (160 -> 320, 8 rows): second order in the in-plane spacing. Richardson (p = 2) limit: -0.143616 m (-0.009 % of the beam value) from 160x8/320x8 (same rows); the report's +0.01 % came from the 160x4/320x8 pair (rows 4 vs 8): -0.143642 m, +0.009 %. Either way the extrapolated hex8 value is within 0.01 % of the beam value (INC-3). Thickness rows 4 vs 8 change the midpoint by 0.05 pp.
- **Reference check (my own re-derivation, `scripts/beam_ref.py`, results in `results/beam_reference_rederivation.json`).** A geometrically nonlinear plane-strain Kirchhoff-Love beam (E_eff = E/(1 - nu^2) = 5.0e8 Pa, h = 0.05 m, L = 10 m, clamped, follower pressure 200 Pa, cubic Hermite elements, Newton) gives midpoint -0.141914 / -0.143100 / -0.143444 m for 10 / 20 / 40 elements (observed order 1.78), extrapolated -0.143585 m (observed order) or -0.143559 m (p = 2): -0.03 % / -0.05 % relative to the oomph-lib value -0.143628 m, which is the size of the oomph-lib tutorial-resolution uncertainty (0.13 % static refinement 1 vs 3, 0.05 % continuum offset, per the reference README). It is the same model class as the oomph-lib beam (not an independent continuum solution), and the oomph-lib beam source was not read: it shows that the static value is reproducible from the stated beam equations and parameters to about 0.05 %, not that the beam equals the continuum wall. No hex20 or continuum FEBio reference run was made (no new simulations). The linear (small-deflection) result 200 L^4/(384 EI) = 1.0 m shows the wall is membrane dominated at the static load.
- Classification: tutorial wall 80x8: **DIAGNOSTIC DISCREPANCY** (PROVEN, hex8 too stiff by 11.6 %). 320x8: **QUALIFIED AGREEMENT** (not VERIFIED: one code, in-plane refinement only, reference is a beam model). 160x4: strict gate missed (see the A1 review in `PLAN_AMENDMENTS_AND_ACCEPTANCE.md`).

## 2. CC-F1 fixed collapsed shape, fluid only (job 14218; 4 cores; 31 s; 7,040 cells; K = 1e7, mu = 0.02, rho = 1; dt 2.5 s x 40 steps, rhoi 0)
Shape = S1 320x8 final interface (`inputs/shape_S1_320x8.json`, regenerated bit-identically by `scripts/make_shape.py`). Definitions: `PARITY_DEFINITIONS.md`. Table: `results/f1_fixed_shape_summary.csv`; wall rows `results/f1_wall_rows_final.csv`; figure `figures/fig2_f1_wall_normal_row_parity.png`.

| quantity | value | pre-registered gate | result |
|---|---|---|---|
| steps / Newton | 40/40, mean 1.575, no negative Jacobian | all steps | PASS |
| mass balance Q_out/Q_in - 1 [S] | 6.0e-7 (Q_in 0.8810 m^2/s per unit depth) | 1e-3 | PASS |
| steady drift of F1, F2 over last 20 % [R] | 1.2e-11, 4.7e-12 | 1e-3 | PASS |
| streamwise parity a_t, wall row [S] | 2.4e-4 Pa = 6.8e-5 of V (V = 3.51 Pa) | 5 % of V | PASS |
| F1 / F2 net F_y | 31.68 / 34.07 N/m (-6.99 %) | 3 % (3-8 % conditional) | conditional |
| F1 / F2 first-mode load G | 15.33 / 16.77 (-8.63 %) | stop above 8 % | FAIL |

- Wall-normal row parity (smooth in x) [S]: |p(wall) - p(row 2)| <= 4.1e-3 Pa (rms 2.4e-3 Pa; the report said 1e-3, INC-4); p(wall) - p(row 1) ranges from -1.34 to +0.55 Pa (rms 0.87 Pa); D = p0 - 2 p1 + p2 from -2.67 to +1.10 Pa (rms 1.75 Pa). So the wall row and row 2 carry one pressure family, row 1 the other. F1 [R of the applied nodal load] samples the wall row only, F2 [R] averages the two families.
- Junction rule (b) was literally true by 7e-4 / 5e-4 Pa (a monotone gradient), judged a non-event; no pin run (recorded).
- Classification: steady fluid mass balance, steadiness and streamwise parity (own gates, one mesh, no external reference): **QUALIFIED AGREEMENT**. Applied load (F1 vs F2 steady, first-mode): **DIAGNOSTIC DISCREPANCY**.
- Interpretation limit [A]: the steady force is small (F_y about 32 N/m, 1.6 % of p_ext L = 2,000 N/m), so the relative spread overstates the dynamic importance of a 2.4 N/m absolute spread (0.06 % of 4.2e2 Pa x 10 m); D1 tested the dynamic regime.

## 3. D1 coupled diagnostic, 20 steps through the first collapse (job 14219; 4 cores; 2 m 18 s; 9,600 cells; dt 0.025; rhoi 0.9; K = 1e7; E_mesh = 1)
Purpose: the single test of A2's hypothesis H_D (the spread is small in the dynamic regime). **It is a diagnostic, not a validation.** Tables: `results/d1_load_spread_by_step.csv`, `results/d1_freebody_identification.csv`, `results/d1_newton_per_step.csv`; figures `fig3`, `fig4`, `fig5`.
- Health [S]: 20/20 steps, Newton 5-6 per step (mean 5.25), no negative Jacobian, normal termination. Volume balance [R]: max |integral (Q_in - Q_out) dt - channel area change| = 0.24 % of the peak area change; the inlet flux reverses sign at about t = 0.17 s.
- Applied-load spread [R]: net wall force and first-mode load from the nodal wall row (F1, what the wall receives) against the element mean (F2) on the moving interface, force pushing the wall up positive:

| t (s) | F1_y | F2_y | (F1 - F2)/F2 | G spread | mean y-parity D (Pa) |
|---|---|---|---|---|---|
| 0.10 | 791 | 663 | +19.4 % | +19.4 % | 51 |
| 0.15 | 1,342 | 1,078 | +24.5 % | +24.6 % | 106 |
| 0.20 | 1,702 | 1,301 | +30.8 % | +31.2 % | 161 |
| 0.25 | 1,659 | 1,172 | +41.6 % | +42.4 % | 196 |
| 0.30 | 1,313 | 808 | +62.4 % | +65.1 % | 203 |
| 0.35 | 522 | 116 | F2 near zero (349 %) | 317 % | 164 |
| 0.40 | -717 | -865 | -17.0 % | -18.6 % | 60 |
| 0.45 | -2,093 | -1,843 | +13.6 % | +14.4 % | -101 |
| 0.50 | -2,842 | -2,234 | +27.2 % | +29.3 % | -245 |

  H_D criterion (maximum spread <= 3 % over t >= 0.10 s with |F2_y| >= 100 N/m): maximum 349 % (near the zero crossing of F2) and 19-65 % in the loading phase 0.10-0.30 s: **H_D REFUTED**. The streamwise a_t stays below 3e-3 Pa at all steps; the mode is wall-normal.
- Pressure at x = 10, t = 0.25 s [S] on rows (wall, 1, 2, channel floor): 193.4, 77.0, 192.9, 191.7 Pa: the wall row and row 2 agree (even family, about 193 Pa), row 1 (odd family) is 116 Pa lower; amplitude (193.4 - 77.0)/2 = 58 Pa about the two-family mean 135 Pa. (The report's "+-50 Pa on 130 Pa, 36 %" is a rounded version, INC-7.)
- Free-body identification [R], valid only while |u_mid| <= 0.1425 m (t <= 0.35 s): inertia is negligible (rho_s h = 0.05 kg/m^2), so the uniform-equivalent fluid pressure that the wall balance requires is p_f = p_ext - p_el(u_mid) with p_el from the S1 320x8 static curve, compared with the psi-weighted equivalents of F1, F2 and the odd row (Pa, positive pushes the wall up):

| t | identified | F1 equivalent | F2 equivalent | odd-row equivalent |
|---|---|---|---|---|
| 0.20 | 173.7 | 182.2 | 139.0 | 95.7 |
| 0.225 | 183.3 | 184.0 | 134.7 | 85.2 |
| 0.25 | 179.8 | 178.5 | 125.4 | 72.4 |
| 0.30 | 141.3 | 141.2 | 85.8 | 30.4 |
| 0.35 | 58.9 | 59.4 | 14.0 | -31.4 |

  The load that moves the wall equals F1 to within 0.1-1 % for t = 0.225-0.35 s (0.4, 0.7, 0.1, 0.9 %); F2 is 27-76 % below it. At 0.075-0.15 s the static-shape assumption is poor (see the figure). The stored JSON has the opposite sign for the identified value (convention p_el + inertia - p_ext, INC-5); magnitudes agree. This supports (it does not by itself prove) that the load applied to the wall is the parity-contaminated nodal wall-row pressure. It does **not** show which of F1 or F2 is the physically correct fluid load: the true pressure varies smoothly across the 0.0625 m cell (physical cross-row variation rho a_y h about 0.3 Pa [A]), so neither family is exact at the nodes; F2 is the cavity-validated better estimate only by a different, indirect argument (see PR #27).
- Wall response against oomph-lib (informational) [S]: FEBio midpoint -0.0523 / -0.0874 / -0.1248 / -0.1764 / -0.1975 / -0.2016 m at t = 0.25 / 0.30 / 0.35 / 0.425 / 0.475 / 0.50 s against the converged oomph-lib -0.0560 / -0.0979 / -0.1451 / -0.2041 / -0.2157 / -0.2101 m. Maximum difference 12.8 % of the 0.216 m peak at t = 0.425 s (rms 6.5 % over t <= 0.5 s; 11.3 % maximum against oomph at the same dt = 0.025 s); early collapse (t <= 0.2 s) within 0.7 % of the peak. The FEBio wall is still descending at t = 0.5 s (min at the last step: trough NOT bracketed), so its trough is later than 0.5 s (> +6.7 % against 0.4685 s) and the value at 0.5 s is 4 % of the peak short.
- Relation to the load-balance law ([A], `LOAD_BALANCE_ANALYSIS.md`): the lag has the sign and the size that the law gives for a one-signed amplification of the applied load by 25-60 % (0.23 eps trough, eps/2 period). This is **consistent with, not proof of**, the parity mechanism as the sole cause: single mesh and time step, rhoi 0.9, E_mesh 1, K = 1e7, the hex8 wall 0.76 % too stiff, and the unknown error of which of F1/F2 is correct are not separated.

## 4. Classification (five terms, as reported, with the corrections of this package)
| quantity | class | basis |
|---|---|---|
| static wall, tutorial mesh 80x8 | DIAGNOSTIC DISCREPANCY | S1 |
| static wall, 320x8 | QUALIFIED AGREEMENT | S1 (not verified) |
| steady fluid: mass balance, steadiness, streamwise parity (one mesh, own gates) | QUALIFIED AGREEMENT | F1 |
| applied nodal wall load, steady and dynamic | DIAGNOSTIC DISCREPANCY | F1, D1 (PROVEN as a data statement; cause SUPPORTED) |
| first trough value and time, period, rebound | NOT ESTABLISHED | no full 140-step run |
| coupled FEBio vs oomph-lib comparison | NOT ESTABLISHED | no full run; 20-step diagnostic only |
