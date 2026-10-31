# Coupled torsion: the Newton failures, the stall analysis and the failed attempts

**Coupled torsion with the proposal's Mooney-Rivlin law in stock FEBio v4.13 is NOT ESTABLISHED.** No coupled field exists to compare with the reference. This file documents what was run, what failed, what is supported and what is unresolved. Figure: `../figures/fig5_coupled_newton.png`; statistics: `../results/run_summary.csv`, `newton_iterations_per_step.csv`, `rtc_m1_step2_residuals.csv` (from the solver logs).

## 1. Smoke tests on the tiny M0 mesh (16 x 2 x 4: 32 solid, 64 fluid cells; 20 steps of 0.01 s; not scientific results)
Coupled `fluid-FSI` model: Mooney-Rivlin three-field hex8 solid, `math` rotation BC on the inner circle, `fluid-FSI traction` on the interface, conformal no-slip. All with the stock binary.

| Job | Variation | Outcome (solver log) |
|---|---|---|
| 14197 | fluid module alone, attempt 1 | input rejected (unrecognised `dtol` tag); 14199 (attempt 2) completes, 20 steps, 1 iteration per step |
| 14198 | coupled MR, rho_inf = 1, min_residual 1e-30 | step 1 (t = 0.01) fails to converge: the relative norms cannot be met at drive amplitudes of 1e-46 rad; leads to amendment A1 (min_residual 1e-18) |
| 14200 | coupled MR, rho_inf = 1, E_mesh = 1 Pa, A1 | steps 1-4 accepted, **step 5 (t = 0.05, rotation 8.8e-11 rad) does not converge** (squared residual 5.6e-12 at the start of the step, 1.5e-8 after the first update (about 2,700 x), then stagnating near 2-5e-9 for the remaining iterations until "Max nr of reformations reached") |
| 14201 | same, E_mesh = 100 Pa | same stall at step 5 |
| 14203 | `use_three_field_hex` tag | input rejected (unrecognised tag): not a test |
| 14204 | MR with the domain attribute `three_field="off"` | same stall at step 5 |
| 14205 | outer fluid-mesh circle clamped | same stall at step 5 |
| 14208 | fluid viscosity 0.3 instead of 300 Pa s | same stall at step 5 |
| 14202 | St Venant-Kirchhoff solid instead of MR, rho_inf = 1 | **converges, 20/20 steps**, 1-10 iterations per step (mean 6.85) |
| 14206 / 14207 | fluid-FSI mixture only, interface circle prescribed to rotate (no solid) | 14206: input rejected (material id); 14207 converges, 1-3 iterations per step |
| 14212 | coupled MR, **rho_inf = 0.9** (single pre-registered diagnostic) | **converges, 20/20 steps**, iterations 1,1,1,1,4,6,7,8,8,9,9,9,9,9,9,9,9,10,10,10 (mean 6.95) |

Reading: the difficulty sits in the coupling of the Mooney-Rivlin solid with the fluid-FSI mixture in this configuration, not in the rotating ALE mesh (14207) nor in viscous traction (14208). The stall is independent of E_mesh, the three-field switch, the outer mesh clamp and the viscosity.

## 2. Source analysis (no run) and hypothesis H-N
`FEElasticSolidDomain::UpdateElementStresses` (`FEBioMech/FEElasticSolidDomain.cpp`, lines 648-671; the three-field domain has the same condition, `FE3FieldElasticSolidDomain.cpp` from line 567) adjusts the stress for strain-energy conservation, `s += D ((W_t - W_p)/(dt J) - s:D)/(D:D)`, **only when alpha_f == 0.5 exactly, i.e. rho_inf = 1, and only when D:D exceeds machine epsilon**. For the proposal's law, W = U(J) + dev W with U = K/4 (J^2 - 2 ln J - 1), the increment W_t - W_p at strains of 1e-10 is cancellation noise (about K/4 x 1e-16), so the correction is noise divided by D:D, and the Newton tangent does not contain it. The stall appeared at step 5 where the rotation (8.8e-11 rad) makes D:D reach machine epsilon, is independent of everything varied above, is absent for the mixture alone (no solid), and absent for St Venant-Kirchhoff (energy from E without cancellation): all as observed.

**H-N (the stall is caused by that correction acting on round-off-level energy): SUPPORTED, not PROVEN.** One diagnostic (14212, rho_inf = 0.9, alpha_f = 0.526, correction inactive) removed the stall at M0. Its pre-registered acceptance also required a mean of at most 6 iterations per step; the observed 6.95 does not meet it (the iteration pattern equals the SVK run: a separate slow-convergence effect, UNRESOLVED, consistent with an inexact tangent term; not tested). The coordinator's rule keyed on the stall and the RT-C gate itself is at most 12 iterations, so RT-C M1 proceeded with rho_inf = 0.9 (declared deviation, amendment A3). Other mechanisms were examined on paper and not supported (solid velocity/displacement mismatch of the prescribed rotation, bulk/shear ratio K/mu = 2.17, skeleton block).

## 3. RT-C M1 (job 14213, 4 cores, rho_inf = 0.9, reformation cap 25): FAILED numerically
M1 = 512 solid + 1,024 fluid cells, 200 steps planned. Step 1 (t = 0.01) accepted (squared residual 1.5e-21). **Step 2 (t = 0.02) diverged**: successive squared residual norms 1.5e-18, 1.0e-14, 5.3e-15, 3.6e-15, 7.8e-12, 1.7e-8, 4.0e-5, 9.1e-2, 2.1e2, 5.2e5, "Problem is diverging", then "32 negative jacobians detected", failure to converge at t = 0.02 (3 s wall, 12 CPU-s). The drive amplitude at t = 0.02 is 0.1 S(t) sin 3t = 3.2e-24 rad (S = 5.4e-22), i.e. round-off level, so the growth is independent of any physical amplitude. The M0 remedy (rho_inf = 0.9, no energy-conserving correction) therefore did **not** carry over to M1: a different or additional instability takes over. **Cause UNRESOLVED**; remaining candidates, untested: inexact interface/mixture tangent, near-singular massless-skeleton block at the M1 geometry. By the coordinator's rule (one diagnostic, no sweep) RT-C was stopped; RT-C M2 was not run.

Consequence: **alpha_B, T_B, sigma_rt and H_f of the coupled problem are NOT ESTABLISHED.** The analytical expectation (first-order traction alone: alpha_B 8.0 / 4.0 %, H_f 8.2 / 4.1 %, T_B 6.7 / 3.5 % at M1 / M2) stays an untested estimate; the declared classification rule (a miss within +-25 % of it is a DIAGNOSTIC DISCREPANCY with confirmed mechanism) was never applied because no coupled result exists. The coupled ALE-with-rotating-mesh path beyond the smoke tests is NOT TESTED.

## 4. Preserved failed or superseded attempts
Static: 14186/14188 (law test design error, not evidence), 14190/14191 (duplicate node-set name in the input), 14192 (valid, outer reaction not logged). Smoke: 14197, 14198, 14203, 14206 (inputs rejected or min_residual, listed above). Their inputs (`../inputs/failed_attempts/`) and logs (`../results/raw/*/failed/`) are kept; they count towards the 0.1350 CPU-h (`ACCOUNTING.md`).
