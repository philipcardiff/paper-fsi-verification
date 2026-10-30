# FEBio v4.13 cavityFlexibleBottom: consolidated findings

Status: **archival diagnostic**. FEBio v4.13's monolithic fluid-FSI solver did **not** reproduce the benchmark force or displacement at the three mesh levels run; none of the three coupled runs terminated normally; the cause of the force gap is **not established**. Nothing here is a verification, a correction, or a proposal to change FEBio.
Evidence labels: PROVEN / STRONGLY SUPPORTED / HYPOTHESIS / FAILED / NOT TESTED. Numerical residual convergence, physical accuracy and benchmark verification are three different things; a PASS in this report is a stated numerical criterion, never benchmark verification. Forces are per metre of depth (N/m). `results/TEST_MATRIX.csv` has one row per test; discrepancies between source reports and raw data are in `provenance/INCONSISTENCIES.md`.

## 1. Question and set-up
Does FEBio v4.13 (stock tag, no source change) solve the cavityFlexibleBottom benchmark (Re 100, StVK plate E=500 Pa, nu=0.4, plane strain, clamped ends, p=0 at the outlet) to the steady values of solids4foam (u_y at A = (4,0): -0.2407 m extrapolated; F_y about -5.26 N/m) when the problem is run as a monolithic, weakly compressible (K=1e6 Pa), generalised-alpha pseudo-transient fluid-FSI model on conformal hex8 meshes (M1/M2/M3: 1688/6752/27008 fluid + 160/640/2560 plate cells)? The FEBio model matches the benchmark in geometry, material law, fluid parameters, inlet flux, outlet datum, depth normalisation and output definitions (tests A1-A8), with two documented differences: four floor-end corner nodes are pinned to p=1.15 Pa (a numerical constraint, not a benchmark condition) and the fluid is the equal-order hex8 weakly compressible scheme.

## 2. Coupled results (tests D1-D6, E4)
| level | dx | last accepted t | u_y(A) (m) | F_y raw, F1 (N/m) | termination | steady spread (last 10 %) |
|---|---|---|---|---|---|---|
| M1 (C2_fsi_pilot_f1, job 14016) | 0.1 | 972.5 s | -0.211967 | -4.7176 | ERROR step 390 (t=975), max reformations on the dilatation floor | 7.6e-5 (met, window [875,972.5]) |
| M2 (C3_fsi_M2, job 14019) | 0.05 | 870 s | -0.226233 | -4.9339 | ERROR step 349 (t=872.5), same | 1.8e-4 (not met) |
| M3 (C3_fsi_M4, job 14021, 8 cores, 28.2 CPU-h) | 0.025 | 715 s | -0.230038 | -4.9835 | ERROR step 287 (t=717.5), same | 1.05e-3 (not met) |
- **All three runs terminated abnormally** and are conditional steady estimates. The error is the relative dilatation-increment criterion (1e-12 of the step's own increment, about 1e-22) sitting on the arithmetic floor at steady state while every other norm is at 1e-23..1e-27: a stopping-rule artefact, not a divergence (0 cut-backs, no negative Jacobian, 7.1 iterations/step). A changed dilatation tolerance was not tried. M3 is not steady to the stated criterion (decaying 44 s oscillation, about +/-2e-4 m).
- Mass imbalance Q_out/Q_in - 1: +1.2e-4, -1.9e-5, -2.3e-5 (numerical conservation, passes 1e-3).
- Same-dx comparison with solids4foam (u_y vs partitioned, F_y vs quasi-monolithic #618): u_y +3.6/-2.2/-3.5 %, **F_y -7.9/-5.6/-5.1 %**. Indicative three-level extrapolations: u_y -0.2314 (order 1.91; -3.9 % vs -0.2407), F_y about -5.00 (-4.9 % vs -5.26; the F_y order is not reliable because the M2->M3 change equals the odd-even bracket). The gap **flattens, it does not close**: a mesh-resolution error alone is not supported. (STATUS: PROVEN for the numbers, INDICATIVE for extrapolation, FAILED to match the reference.)
- The plate deflects as expected for the load it receives (u_y/F_y 0.0449/0.0459/0.0462 m per N/m; solids4foam 0.0458 extrapolated), so the displacement gap is the force gap (STRONGLY SUPPORTED). u_x(A) is 3.4-4.8 times the solids4foam value (floor pressure tilts left to right).

## 3. What was excluded
| candidate | result | label |
|---|---|---|
| bulk modulus (K x10) | u_y -1.1e-6, F_y -1.4e-6 relative (M1) | PROVEN (M1) |
| pseudo-time step (dt/2) | +4.1e-5 / +2.2e-5, same size as the steady-window spread | PROVEN (M1) |
| tolerances 100x tighter | identical to 8 digits (criteria not binding) | PROVEN (M1) |
| corner pins | bound 0.41/0.19/0.10 % of F_y, wrong sign to close the gap | PROVEN (linear attribution) |
| viscous traction | 0.03 % of F_y (2-point estimate; the solver-applied viscous traction is NOT reconstructable from stored output) | STRONGLY SUPPORTED |
| shared pressure BCs (outlet p=0, free underside, rho, flux, no-slip) | equal in the two models | PROVEN |
| plate / hex8 stiffness | explains part of the displacement gap at M1/M2 only (4.5 %, 1.3 %, 0.3 %), not the force | PROVEN |
| pressure transfer (S2) | 0.025 % self-consistency, pressure part only | PROVEN (pressure part); cannot detect a wrong pressure |
| deformation difference | <= 0.04 % of F_y | STRONGLY SUPPORTED |
| force integration / reconstruction | F1 reconstructed from stored fields to 1e-7 | PROVEN |

## 4. Force definitions at the same stored times (tests E6, F10)
`results/force_reconstruction_F1_F4.csv` (pressure contribution; **viscous contribution not independently reconstructed**: a 2-point estimate bounds it at <=0.0015 N/m, but the solver-applied value uses Gauss-point gradients that are not stored).
| state | F1 applied (nodal wall-row p; the load that moves the plate) | F2 FEBio element-mean surface-pressure definition (= F4 for planar facets) | F1-F2 | solids4foam QM, same dx |
|---|---|---|---|---|
| M1 t=900 | -4.7176 | -5.4375 | +0.720 | -5.1226 |
| M2 t=750 | -4.9339 | -5.2748 | +0.341 | -5.2266 |
| M3 t=500 | -4.9809 | -5.2469 | +0.266 | -5.2504 |
| M3 t=717.5 (failure-exit state, flag 2; separate row) | -4.9835 | -5.2494 | +0.266 | -5.2504 |
- F1 is **fully exposed** to a row-to-row (wall-normal) alternation of the nodal pressure; F2 averages the wall row and the first row above it. The parity part of F1 (the difference to a 1-2-1 smoothed field) accounts for F1-F2 to within 0.0002-0.001 N/m (PROVEN arithmetic decomposition). Row-parity amplitude D = p(row1) - p(row0) = 0.363/0.171/0.134 Pa (M1/M2/M3; windows in `INCONSISTENCIES.md` I-5; 0.1335 Pa in the M3 failure-exit state); the alternation is flat through at least six rows and deterministic (unchanged by dt/2, K x10, tolerances; present with the plate fixed).
- **F2 is not established as the physical interface force.** It agrees with the reference at M3 (-0.07 % at t=500) but is +6.1 % off at M1 and +0.9 % at M2 (opposite sign to F1), and it is by construction FEBio's own averaging, not an independent estimator. This table is a re-evaluation of stored solver fields, **not an independent solver verification**.

## 5. Where the alternation comes from (tests F4-F12; `SOURCE_CODE_INVESTIGATION.md`)
FEBio v4.13's fluid-FSI domain is an equal-order trilinear hex8 element for relative velocity and dilatation, with the pressure p=-K e_f, a strong-form pressure gradient in the momentum row, and **no interior pressure stabilisation** (PROVEN by source search). The operator [A G; G^T -(c/K)M] has soft parity pressure modes whose eigenvalues scale as h^2 while smooth modes do not (PROVEN on two levels of one geometry); they are not null modes (the compressibility term is 2000 times smaller, and the linear emulation relaxes a perturbation in one step) and are excited by the flow's own viscous (42 %) and convective (58 %) force terms in the separated cavity flow (exact linear attribution on the C1b field: fixed values -1e-4 Pa). Class assessment: **discrete numerical deficiency of the unstabilised equal-order pair: STRONGLY SUPPORTED; setup error: not supported; source defect: not supported** (the independently re-assembled equations reproduce FEBio's converged C1b state to 9.4e-7 / 2.2e-7); whether the parity-free value is the true force and whether it would restore agreement is **NOT ESTABLISHED**. The four corner pins are necessary (without them the corner pressure DOF is undetermined: Schur eigenvalue 3.5e-10) and are not the cause.

## 6. Two research stabilisations failed (tests G1-H4; `FAILED_STABILISATION_HISTORY.md`)
- Round 1 (polynomial pressure projection, alpha = 3/4, fixed-plate C1b): the mode is removed (D 0.4207 -> 9.4e-5 Pa) and F1 ~ F2 (0.105 %), but **mass conservation fails** (Q_out/Q_in - 1 = -1.55 %, stock +9e-5), the smooth pressure shifts 3-7 %, reactions +12 %. The agreement F1 = F2 is a consequence of parity removal and carries no information about the true force (F1 moved +10.1 %, F2 -6.3 %).
- Round 2 (flux-preserving redistribution, alpha = 1e-3, exact Poiseuille channel): global row balance to machine precision, but three of six pre-registered criteria fail (nodal column flux 4.9e-4 > 1e-4, dp/dx -0.29 % > 0.2 %, y-alternation only about 2x lower near the outlet). Cavity test not justified.
- Neither is a correction or a solution; neither is proposed for FEBio; no coupled steady run was made with either; the FSI-domain implementation was exercised for 5 steps only.

## 7. Reference classes (see `COMPARISON_TABLES.md`)
StVK-confirmed: solids4foam partitioned and quasi-monolithic (shared code lineage, a consistency pair). Constitutive law unresolved: Tukovic et al. (u_y -0.2506 finest), ANSYS (-0.2527), COMSOL (-0.243): their displacement spread (4-5 %) is a plate-law effect (plate-only Q9: neo-Hookean +4.8 %, Hencky +5.7 % vs StVK) and is a separate question from the FEBio numerical problem discussed in manuscript-repository PRs #12 (ANSYS/COMSOL results) and #22 (COMSOL/Tukovic references and Richardson extrapolation; both open when this was prepared). The **forces** of all references agree to 0.2-0.7 %; FEBio's raw F1 is about 5 % below them at every level. No source of the solids4foam floor pressure or viscous traction field exists in reach, so no local comparison was possible.

## 8. Conclusion
FEBio v4.13 with this model is stable, mass-conserving to 1e-4 and self-consistent (pressure part), converges in dx at order about 2 to values about 4-5 % below the solids4foam StVK benchmark in force and displacement, ends every long run on an iteration-floor error, and carries a diagnosed numerical row-parity pressure mode that contaminates the applied wall load. Two experimental remedies failed. The benchmark is **not** verified for FEBio and the cause of the gap is not identified. Minimal criteria for reopening are in the package `README.md`.
