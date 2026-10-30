# Failed experimental stabilisation history (REJECTED RESEARCH PROTOTYPES)

Neither stabilisation is a correction or a solution; neither is proposed for FEBio. The diffs (`../experimental/`) are research prototypes of a Dohrmann-Bochev polynomial pressure projection added to the continuity row, built in a separate directory from the pristine tag (the stock source and binary were never touched; stock `febio4` sha256 f5ae1ed6... unchanged after all work). Switch: environment variable `FEBIO_PPP_ALPHA`. Derivation and critical review of earlier claims: `reports/FEBIO_PRESSURE_PROJECTION_DERIVATION.md`. Budget: round 1 0.68 CPU-h (limit 3.0), round 2 0.046 CPU-h (limit 0.75); sacct totals in `../provenance/slurm_jobs.csv`. Every acceptance threshold was fixed before the corresponding run; nothing was tuned to -5.25 N/m or any solids4foam value.

## Round 1: alpha = 3/4 (predeclared from the Schur-complement bound 3/(4 mu))
| stage | test | result | verdict |
|---|---|---|---|
| 2 | element residual/tangent, 12 cases | alpha=0 bitwise no-op; FD 2.5e-13 / 3.5e-10; parity energy exact | PASS (G1) |
| 2 | tiny channel, alpha=3/4 as first built | max reformations at step 1 (Newton contraction 0.55 per iteration) | FAILED (G5) |
| 2 | diagnosis | lagged inlet rim pressure; remedy (rim once per step when alpha>0) gives quadratic convergence and the same converged state | PASS (G5) |
| 2 | tiny exact Poiseuille channel | y-alternation 6.3e-4 -> 1e-8 Pa; profile error 1.01 -> 0.26 %; dp/dx error -0.04 -> -0.79 %; flux at x=2 -0.74 % (predicted alpha (h/H)^2 Q = 0.75 %) | mass conservation deteriorates (G3) |
| 3 | alpha=0 control on C1b | field differences 6.6e-12 / 1.9e-12 / 8.3e-12 (<= pre-registered 1e-8), 1179 vs 1180 iterations | PASS 5/5 (G2) |
| 4 | alpha=3/4 on C1b | D 0.4207 -> 9.4e-5 Pa (S2 pass); momentum residual 8.4e-7 (S4 pass); F1-F2 0.8196 -> 0.0054 N/m (S7 pass); **Q_out/Q_in - 1 = -1.55 % (S5 FAIL; stock +9e-5; flux deficit 0.78 % already at x=1.03)**; smooth pressure -3.4/-5.3/+3.0/+6.6 % by region, velocity L2 3.7 % (S3 FAIL); reaction sum +12 % (S6 FAIL) | **FAILED** (G4) |
| 5 | coupled M1 with the stabilisation | not run (gate) | NOT TESTED (G7) |
Why it failed (derived, matched numerically): the zero-sum-per-element stabilisation rows do not cancel at inlet nodes with prescribed flux and in the unenforced rows (outlet, pins, rim); the consistency error is first order in alpha and O((h/H)^2) with O(1) constants, not small at h/H = 0.1. The straight-channel test shows the interior profile and gradient move closer to exact (1.0 -> 0.26 %, 0.98 -> 0.21 %) while the flux is 0.7 % worse. That F1 and F2 agree after the stabilisation is a consequence of parity removal; both moved (F1 +10.1 %, F2 -6.3 %), so it is not evidence about the true force. Scaling of the defect with alpha and (h/H)^2 is a prediction, **NOT TESTED** (no second alpha was run).

## Round 2: flux-preserving redistribution, alpha = 1e-3
Form: stabilisation rows of constrained-dilatation nodes (outlet, rim, pins) redistributed uniformly over the free rows of the element, and the constrained rows set to zero; alpha = 1e-3 declared from mu x the lowest smooth Schur eigenvalue (not from any force).
| test | result | verdict |
|---|---|---|
| assembled random mesh, 16 constrained nodes | free-row sum: old -9.7e-6, new -2.9e-20; FD 1.2e-12 / 3.2e-9 | PASS (H1) |
| P0 rim-once at alpha=0 | 5.9e-12 relative to stock | PASS (H2) |
| P1 max abs nodal column flux error < 1e-4 | 4.86e-4 (x2 -4e-6, x7 1.5e-5, x10 7.5e-5, x12 1.9e-4, x13 2.9e-4, outlet 4.9e-4) | **FAIL** |
| P2 dp/dx within 0.2 % | -0.29 % | **FAIL** |
| P3 profile <= 0.5 % of centreline | 2.0e-4 | PASS |
| P4 y-alternation <= 6.3e-5 Pa | 3.0e-4 at x=12 (stock 6.2e-4: only about 2x lower) | **FAIL** |
| P5 re-assembled residuals | momentum 1.5e-9, continuity 1.4e-7 (7.6e-5 without the term) | PASS |
| P6 convergence | 274 iterations (1.37/step) vs 492 | PASS |
Stopped at the first failure by the pre-registered rule; no second alpha. Diagnosis (login-node post-processing of the same run): the stabilisation behaves as derived (cumulative redistributed sum 1e-5..2e-5 as predicted; pair-sum flux error 2.0e-5), but the nodal column flux, which the continuity rows constrain only pairwise, develops an x-sawtooth growing from 2e-5 at the inlet to 4.9e-4 at the outlet (stock 2e-6): about 300-fold amplification. Cause of the amplification: inferred mechanism (the stock x-sawtooth is held small by coupling to the pressure x-sawtooth that the term now damps), **UNRESOLVED**, not tested. A cavity fluid-only test was judged not justified. The FSI domain of round 2 is compiled and element-tested only.

## Not done
M2 validation, any coupled steady run, a second alpha, jump-based or macro-element projection variants, Taylor-Hood / higher-order fluid elements. Possible next design directions are in `reports/FEBIO_STABILIZATION_FINAL_DECISION.md` (proposals, not executed).
