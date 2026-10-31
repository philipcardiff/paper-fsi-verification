# Stage C: rigid cylinder + flag (CFD1-type) with the FSI1 inlet/outlet formulation

One level (8,330 hex8 cells), two runs: C1u (unpinned, diagnostic) and C1p (pinned flag-root junction). Level 2 and the sensitivity runs were not made: the pre-registered rule "the first failed gate stops stage C" applied after level 1 (C3, C4, C5 failed by one to two orders of magnitude; a second-order refinement could not cure a 4 Pa amplitude against a 0.048 Pa limit). Model and boundary conditions: `PROVENANCE_AND_CONFIGURATION.md`. Normalisation: the slab is one element thick (TZ = 0.015 m) with w_z = 0 and z-uniform fields; all forces here are line integrals of the z = 0 layer, i.e. **N per metre of spanwise depth**, directly comparable with the benchmark (a force summed over the 3-D slab would be TZ = 0.015 times larger). The xplt velocity per element equals the mean of the logged nodal velocities to 1.5e-8 m/s (max velocity 0.40 m/s): the fields are z-uniform, so the 2-D reduction is exact for the data used.

## 1. Reference and provenance
CFD1 (Re = 20, rigid flag) from the Featflow CFD page (fetched 2026-10-31), drag and lift on cylinder + flag: level 2+0 (2,304 elements) 14.2730 / 1.11692; 3+0 14.2880 / 1.11852; 4+0 14.2919 / 1.11896; 5+0 14.2927 / 1.11904; **6+0 (589,824 elements) 14.2929 / 1.11905 (used)**; 5+1, 5+2, 5+3 agree to the printed digits. The values 14.29 / 1.119 quoted in the Turek-Hron paper tables could not be checked (paper not retrievable): UNRESOLVED; they agree with the Featflow numbers to the quoted digits only if recalled correctly. The coupled FSI1 numbers (drag 14.29426, lift 0.7637460) are not used for the rigid flag: the lift differs by 32 % because the flag deflects. Force definition (Featflow page: none given; solids4foam `forces` function object: pressure + viscous over cylinder and flag): F = integral over the wetted surface of (p n_f - 2 mu D n_f) dS with n_f out of the fluid.

## 2. Mass conservation, steady state, Newton
| quantity | C1u | C1p | criterion |
|---|---|---|---|
| inlet flux Q_in - 0.082 | 2.8e-15 | 2.8e-15 | exact discrete flux |
| Q_out/Q_in - 1 | -1.22e-04 | -1.22e-04 | <= 1e-3: met |
| Q(x = 1.0)/Q_in - 1 | 2.56e-04 | 2.56e-04 | <= 1e-3: met |
| drag spread, last 20 % (t = 32-40 s) | 1.2e-02 | 1.3e-11 | <= 1e-3: C1u not steady, C1p met |
| lift spread, last 20 % | 5.8e-03 | 5.6e-11 | <= 3e-3 |
| equilibrium iterations (80 steps) | 266 (3.3 per step) | 262 (3.3 per step) | mean <= 8: met |
Both runs: normal termination (log inspected), no retries or cut-backs, no "Max nr of reformations" message. Per-step maxima were not tabulated.
Unpinned run C1u: the two flag-root junction nodes (concave wall-wall corners) drift: nodal pressures 470.8 and 587.9 Pa against -131.8 .. 195.5 Pa on all other obstacle nodes, drift over the last 20 % of the run 0.51 and 0.63 of max|p|. The pre-declared pin rule fired; C1u is a diagnostic. The pin does not change the parity mode (cylinder odd-even amplitude 14.8651 Pa in both runs; F3 identical to 6 digits).

## 3. Pressure parity: measurement definitions and results (C1p)
Definitions (all in `scripts/parity_diagnostics.py` and `scripts/hronturek_cfd_analyze.py`; the two flag faces are measured separately; nodes are the z = 0 layer wall nodes sorted by x; the 3-4 nodes nearest the root and tip corners are excluded):
- **a_t (P1)**: alternating along-wall amplitude, rms over inner wall nodes of (p_k - (p_(k-1) + p_(k+1))/2)/2 (a pure odd-even mode of amplitude a gives a). The gate quantity (<= 0.048 Pa).
- **D_row**: mean of p_wall - p_row1 over inner wall nodes (row 1 = the node across the adjacent quad); for a smooth field one cell of wall-normal pressure variation, for a wall-row parity mode about +-2a.
- **pbar**: two-row average (p_wall + p_row1)/2: cancels an alternation between adjacent rows. **dp_raw** = mean (p_top - p_bottom) at equal x on the raw wall nodes; **dp_2row** = the same with pbar. The physical cross-flag pressure difference is estimated by dp_2row (an estimate: it assumes the mode alternates between adjacent rows) and is kept separate from the spurious part.
- P2 of the original analysis (wall-row offset from a quadratic extrapolation through the first three rows, top -81 Pa, bottom +3.8 Pa) amplifies an alternating pattern (the extrapolation weights 3, -3, 1 magnify it about 7 times); its magnitude should NOT be read as a pressure offset. The gate failure rests on a_t, which does not have this defect.

| measure | flag top face (y = 0.21) | flag bottom face (y = 0.19) |
|---|---|---|
| a_t (Pa) | 4.109 | 2.630 |
| mean wall pressure (Pa) | 23.90 | 34.54 |
| D_row mean, rms (Pa) | -19.99, 22.95 | -3.03, 13.86 |
| two-row mean pressure (Pa) | 33.90 | 36.05 |
| dp_raw (top - bottom) | -10.64 Pa | |
| dp_2row (top - bottom) | -2.15 Pa; pressure-only transverse force on the flag from dp_2row 0.609 N/m | |
Along the vertical line x = 1.0 m the nodal pressure alternates node by node: [39.01, 7.637, 39.034, 7.65, 39.06, 7.667, 39.091, 7.689] ... Pa (amplitude 15.8 Pa on a mean of 23.5 Pa), while ux is smooth (`figures/fig3_stageC_pressure_parity.png`). The same alternation exists along the grid lines from x = 0.6 to 2.4 m (amplitude 20-23 Pa near the flag tip, 1.4 Pa at x = 2.38; `results/stageC_raw/C1p_L1_f1.25_K1e6/stageC_diag.json`), and in x (column means alternate), and on the cylinder wall (14.9 Pa). Gate C5 (<= 0.048 Pa): FAILED by a factor 55-85 on the wall chains. Context: the uniform transverse load equivalent to the FSI1 tip deflection is 0.48 Pa (analytical estimate), the parity-cleaned cross-flag difference on the rigid flag is 2.15 Pa, the raw nodal difference 10.6 Pa. A small parity amplitude would not have proved a correct integrated force; here the amplitude is large and the force test fails independently (section 4).

## 4. Force reconstructions (C1p, level 1; N/m)
Each line is a different numerical reconstruction from the same stored nodal fields (p = -K e_f, nodal velocities). None is a measured FEBio coupled force and none is presented as the true force: F1 and F-actual are reconstructions of a formula (FEBio's FSI traction was not exercised: no coupled run); F2 is not an independent solver reference; F3b depends on the control-volume reconstruction and on the momentum balance; F4 depends on my reading of FEBio's weak form.

| reconstruction | drag F_x | error vs 14.2929 | lift F_y | error vs 1.11905 | remark |
|---|---|---|---|---|---|
| F1  nodal-pressure facet integral + viscous at the facet Gauss points (one-sided) | 16.54131 | +15.73 % | 2.87588 | +156.99 % | wall-integrated |
| F-actual (reconstruction of the FEBio fluid-FSI traction formula on the rigid-flag fields) | 16.47683 | +15.28 % | 2.87427 | +156.85 % | nodal pressure + element-mean viscous stress (FEBio form) |
| F2  element-mean pressure + viscous at the element centre | 13.99036 | -2.12 % | 1.10655 | -1.12 % | wall-integrated |
| F3  control volume, nodal pressure on the outer boundary | 11.74136 | -17.85 % | 47.01128 | +4101.00 % | contaminated by the parity mode on the channel walls |
| F3b control volume, element-mean pressure on the outer boundary | 14.66320 | +2.59 % | 1.11667 | -0.21 % | momentum balance |
| F4  weak-form nodal reactions on the obstacle nodes (integration by parts) | 14.33785 | +0.31 % | 1.11222 | -0.61 % | uses the stored nodal fields and the solver's weak form |
Parts (drag): F1 pressure 10.039, viscous 6.503; F2 pressure 7.562, viscous 6.428; F4 pressure 7.755, viscous 6.583, convective 0.0003. Full term table: `results/force_budget_C1p_L1_f1.25_K1e6.csv`.
Pre-registered gates: C3 (drag within 3 % at level 1; F1) FAILED (+15.73 %); C4 (F1 vs F2 vs F3 within 1 % of the drag, 3 % of the lift) FAILED (F1-F2 15.4 %, F1-F3 29.0 % of F1; lift F1-F2 61.5 %). Signs and units: n_f out of the fluid; the force on the obstacle is positive in +x for drag and +y for lift (the Featflow lift of +1.119 is upward); the sum includes both flag faces, the tip and the cylinder; forces per metre of depth.
**Independent reproduction.** `scripts/hronturek_cfd1_independent.py` was written separately from the first analysis script (geometry re-read from the .feb input, boundary re-derived from element topology, vectorised kernels with closed-form bilinear derivatives). It reproduces F1, F2, F3, F3b and the mass balance of both runs to 1e-13 (`reproduce.sh`; same definitions, different code); it also adds the weak-form reaction F4 and checks the data path: the logged dilatation equals the xplt dilatation to 7.2e-12 (scale 2e-4) and the xplt element pressure equals the mean of the four nodal pressures to 7.6e-6 Pa. The F4 identity closes: sum of nodal weak-form forces over all nodes minus the convective volume integral = 5e-15 (algebraic identity); the residual of the steady discrete momentum equations at the 8040 interior nodes sums to -0.0147 N/m in x (0.1 % of the drag), maximum single-node residual 0.0083 N/m at x = 0.6, y = 0.183 (near the flag tip) against a typical obstacle nodal force of 0.63 N/m: the stored field satisfies the solver's interior momentum equations at the 1 % level of my reconstruction.

## 5. Why F2 and F3b differ by about 4.8 % (quantified budget, C1p, drag)
F3b - F2 = 0.673 N/m = 4.8 % of F2 (4.6 % of F3b). Using the weak-form force F4 as a common yardstick (not as truth):
- F4 - F2 = +0.347: obstacle side, F2 is low through its pressure (+0.192, element-mean pressure against the weak-form pressure) and its viscous part (+0.155: the one-sided element-centre velocity gradient underestimates the consistent wall shear by 2.4 %).
- F3b - F4 = +0.325: outer-boundary side, made of +0.600 (facet sigma.n on the outer boundary minus the weak nodal sum: inlet pressure +0.003; the channel walls plus the outlet corners together +0.597, i.e. the one-sided wall shear on the 2.5 m top and bottom walls is 20.09 N/m against 20.69 N/m in the weak form, 3 % smaller), -0.290 (momentum flux through inlet and outlet -0.190 against the convective volume integral -0.480: the difference is rho w div(w), the discrete continuity/compressibility term) and +0.015 (interior residual).
- The inlet pressure term is not the cause: element-mean 34.563 against nodal weak 34.560 N/m. The nodal-pressure version of the same control volume (F3) differs by 2.9 N/m at the inlet and by a factor of 40 in lift on the walls because of the parity mode.
- Not separated: the influence of the pins (it changed F3 by < 1e-5), parity on the outer walls beyond the element-mean treatment, the effect of the curved O-grid on the velocity-gradient estimate. The budget is therefore a decomposition by term, not an attribution of truth: F2 and F3b err in opposite directions by roughly 2.4 % (F2 low) and 2.3 % (F3b high) around F4, all from first-order wall-gradient and element-averaging effects. **Closeness of F2, F3b or F4 to 14.2929 N/m (-2.1 %, +2.6 %, +0.3 %) is not evidence that any of them is correct**: these are single-level numbers, the reconstructions disagree with each other by more than the pre-registered 1 %, and the nodal field underneath is parity-contaminated.

## 6. Classification
Pressure parity and nodal-pressure force F1 / F-actual: DIAGNOSTIC DISCREPANCY. F2, F3b, F4 vs the CFD1 reference: NOT ESTABLISHED. Mass balance: criterion met (internal). Level 2, K and pin sensitivities, one-way prediction: NOT TESTED.
