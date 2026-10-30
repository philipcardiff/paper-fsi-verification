# PHASE D FINAL COMPARISON: ringAddedMass, independent FEBio v4.13 evidence

Scope: FEBio v4.13 fluid-FSI, **inviscid-limit formulation** (slip outer wall, slip fluid-ring interface with w.n = 0 enforced, E_mesh = 1e-3 Pa regulariser, K = 1e8 Pa, rhoi = 1, 100 steps per period, 4 periods, hex8 ring and fluid, 1 core). Pristine v4.13 binary + per-node initial-condition workaround for the v4.13 fluid-IC defect (drafted locally, not posted). All numbers below are read from `results/phaseD_*.csv`, which are built by `tools/phaseD_results.py` from solver output and the official solids4foam reference JSON (commit 1576edcc), not typed in.

## 0. This is not the same boundary-value problem as solids4foam
solids4foam solves a nearly inviscid **finite-viscosity** problem (nu_f = 1e-7 m^2/s, no-slip interface and wall, partitioned IQN-ILS, high-order solid). FEBio here solves the **free-slip inviscid limit** that the analytical reference assumes. The two agree on the analytic limit, not on the discrete problem; FEBio's own no-slip attempts (preserved, Phase B/C: -6.0 % M1, -2.7 % M2 ratio error) did **not** reproduce the benchmark because the strong nodal no-slip constraint on an unresolved boundary layer perturbs the added mass. That is a model-form difference, not shown to be a FEBio defect.

## 1. Corrections disclosed
1. **Ratio tolerance.** I earlier called 0.3 % the "ratio screening tolerance". The official JSON gives ratio **0.1 %** (omega 0.3 %, omegaExtrapolated 0.1 %, ratioThinRing 0.25 %). All Phase D figures/tables use 0.1 %. The Phase C documents' 0.3 % wording is superseded.
2. **Reference rounding.** Earlier analysis used 8-digit rounded omega (1.5e-7 relative). Phase D uses the exact JSON values (dry 2.5550421958, weak 2.4356148807 / 0.9532581829, moderate 1.8045131724 / 0.7062557226, strong 0.7686544910 / 0.3008382767). The generator's dt still uses the rounded value (effect ~1e-10 on phase error). The Phase C CSV error columns still use the rounded reference (1.5e-5 % difference).

## 2. Dry/wet convention (reuse of dry runs is legitimate)
solids4foam: 100 steps per period of the exact frequency of each configuration (dt differs dry vs wet), ratio = wet/dry on the same mesh, QoI D(t) = (u_r(0) - u_r(90 deg))/2, frequency from zero crossings / damped fit. FEBio: the same, dt = T_ref/100 with the dry frequency's period for dry and the wet frequency's period for wet, so omega*dt is matched (≈ 0.0628) and the generalized-alpha phase error (-(omega dt)^2/12, rhoi = 1) cancels in the ratio. The dry ring contains no fluid or regulariser, so one dry run per mesh serves all three densities. Dry reuse: B1 (M1), D1 (M2), C1 M3 dry. Dry M1/M2/M3 omega: 2.62106, 2.57108, 2.55843 rad/s.

## 3. Main comparison (E_mesh = 1e-3 Pa, FEBio vs exact vs solids4foam README values)
| regime | mesh | FEBio omega_wet | omega err | FEBio ratio err | s4f ratio err | n=2 force err | leakage rms | amp. loss/period |
|---|---|---|---|---|---|---|---|---|
| weak | M1 | 2.498643 | +2.588 % | +0.0037 % | 0.007 % | -2.13 % | 1.03 % | 5e-5 % |
| weak | M2 | 2.451004 | +0.632 % | +0.0041 % | 0.002 % | -1.06 % | 0.34 % | 4e-5 % |
| moderate | M1 | 1.851253 | +2.590 % | +0.0061 % | 0.038 % | -2.12 % | 1.05 % | 3e-4 % |
| moderate | M2 | 1.815971 | +0.635 % | +0.0072 % | 0.011 % | -1.06 % | 0.40 % | 3e-4 % |
| moderate | **M3** | 1.807043 | +0.140 % | +0.0075 % | 0.003 % | -0.53 % | 0.32 % | 3e-4 % |
| strong | M1 | 0.788795 | +2.620 % | +0.0355 % | 0.068 % | -2.12 % | 1.23 % | 1.7e-3 % |
| strong | M2 | 0.773787 | +0.668 % | +0.0399 % | 0.021 % | -1.07 % | 0.76 % | 1.8e-3 % |

Observations (strictly what the data support):
- **Ratio**: all seven cases are within the official 0.1 % band (largest +0.040 %, strong M2). The FEBio moderate ratio error is +0.0061, +0.0072, +0.0075 %: **flat/slightly rising, not decreasing**. The official criterion "ratio error decreases with refinement" is therefore **not met** for FEBio (weak and strong too). The ratio is mesh-converged to ~0.001 pp at a small non-zero offset; the offset is below the official tolerance but I have not explained it (strong plateau ~+0.04 %: not the regulariser, not K (~1e-5 %), not leakage, which falls while the error does not). Candidates (untested): temporal error not cancelling exactly at the high added-mass level, mesh-dependent n=2 coupling, dilatation noise (strong margins ~0.99).
- **Absolute frequency** (moderate): +2.590, +0.635, +0.140 % for M1-M3. Differences 0.03528 and 0.00893 rad/s, ratio 3.95: observed order **1.98** (three levels, consistency check only). M3 passes the 0.3 % criterion; M1/M2 and weak/strong M2 do not. FEBio uses standard hex8 solids, solids4foam a high-order solid, so absolute-frequency errors at equal "mesh level" are not like-for-like. Wet omega change M2->M3 is 0.49 %, above the 0.1 % "last two meshes" criterion: **not mesh-converged to 0.1 % in absolute frequency**.
- **Regulariser**: E_mesh 1 -> 1e-3 Pa removed an offset of +0.17..+0.18 % in the moderate ratio (+0.170/+0.178/+0.180 % -> +0.006/+0.007/+0.0075 %); mesh-independent effect (-0.164/-0.171/-0.172 pp). Earlier E_mesh = 1 results are retained as diagnostics.
- **Time order**: not tested. Only a paired dt-halving at M1 (E_mesh = 1, moderate): ratio shift +9e-5 %. The official "time order >= 1.5" criterion is untested.
- **Added-mass effective density** versus the exact continuum: weak -0.09 %, moderate -0.03 %, strong -0.09 %.

## 4. Pressure, interface force, leakage
- n=2 pressure harmonic follows the exact radial profile (A r^2 + B r^-2) to 0.04 % at all levels (layer 0 / exact profile 0.9993..0.9997); fitted n=2 phase 178.1-178.2 deg at every level (identical across levels; its reference convention is that of `tools/pressure_modes.py`, and I have not independently checked the ~1.8 deg departure from 180 deg). Amplitudes: weak 0.044, moderate 0.32, strong 1.38 Pa (scales with rho_f).
- Integrated n=2 interface force error vs m_a u_r: -2.12 % (M1), -1.06 % (M2), -0.53 % (M3): halves per refinement and **equals the exact-profile sampling offset at the first node layer** (black ticks in fig D3), i.e. it is a sampling/measure effect of nodal tractions (first order), not a force imbalance. The frequency error converges at second order, so the two must not be quoted as the same order.
- Leakage (w.n weakly enforced): rms relative to ring normal velocity 1.0-1.2 % (M1), 0.34-0.76 % (M2), 0.32 % (M3, moderate). Strong is the worst (1.23 / 0.76 %). Net flux 0.3-1.0 %. Leakage decreases with refinement but does not vanish; volume is conserved (mean J-1 ≈ 5e-10 weak/moderate, 2.9e-9 strong).
- Numerical damping and energy: amplitude loss per period 4e-5 to 1.8e-3 %, no growth; kinetic energy at integer periods within 0.1 % (0.9991..1.0000); half-period KE excursion 0.60 % (strong M1) down to 0.05 % (weak M2) (excursion shrinks with refinement).
- Mesh quality: no Jacobian warnings, no negative Jacobians, no max_refs events in any run; J deviation ≤ 1.7e-8.

## 5. Nonlinear convergence evidence
All 7 runs: 400 steps, ~4.0-4.6 Newton iterations per step (strong M2: 4.55), no reformation-limit events. Worst-step achieved/tolerance (M3): residual 1.2e-2, displacement 1.8e-3, velocity 8e-5, **dilatation 0.96**. Dilatation margins at worst step: weak M1 0.37, weak M2 0.87, moderate M2 0.76, moderate M3 0.96, strong M1 0.989, strong M2 0.992. The dilatation increment sits on a flat arithmetic floor (|inc|^2 ~ 1e-30..1e-29), so the stopping rule cannot be tightened below ~5e-7 at K = 1e8.

## 6. Iterative-error test (D6, failed; honest status)
D6 tightened dtol/vtol 100x and etol/rtol 1e4x on moderate M2. Steps 1-74 gave bit-identical D(t) and iteration counts to the baseline (those criteria were not binding), but step 75 hit the displacement/residual floors and the run failed, so no frequency was extracted. **The iterative error in omega was not isolated.** Supported: d/v/r/e criteria not binding over steps 1-74; dilatation increment at its floor. Not supported: any measured delta-omega from tolerance. A bounded follow-up is proposed in PHASE_D_ITERATIVE_ERROR.md §4 (dtol 1e-7, rtol 1e-12; ~27 min), not run (six-run budget exhausted). Failure preserved in `runs/phaseD/D6_iter_moderate_M2`.

## 7. Official-criteria scorecard (FEBio, `results/phaseD_acceptance_checks.csv`)
Passed (all three regimes): finest ratio error ≤ 0.1 %; ratio vs thin-ring ≤ 0.25 %; amplitude loss < 2 %/period; growth ≤ 0.2 %/period. Moderate also passes finest omega ≤ 0.3 % (M3). **Failed:** omega ≤ 0.3 % at finest level for weak/strong (only M2 run); ratio error decreases with refinement (all); omega change between last two meshes ≤ 0.1 % (all). **Not tested:** time order.

## 8. Publication-readiness decision
Good numerical agreement: yes (ratio within 0.04 % of exact in all regimes; moderate absolute frequency +0.14 % at M3). Mesh-convergence evidence: yes for moderate omega (order ≈ 1.98, 3 levels) and for the interface force; **no** for the ratio (flat offset) and weak/strong have two levels only. Independent solver consistency: partial (different formulation: inviscid-limit vs finite viscosity; the FEBio ratio errors are of the same size as, or slightly above, solids4foam's at the finest level). Remaining uncertainty: iterative error unquantified, time order untested, ratio plateau unexplained (strong +0.04 %), dilatation margins ~0.99. Model-form limitation: not the no-slip finite-viscosity problem.
**Decision: suitable for a DRAFT evidence PR labelled as independent FEBio evidence with these limitations; NOT fully verified, not a Category A claim.**
