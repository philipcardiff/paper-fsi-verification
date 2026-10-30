# RING FINAL COMPARISON: independent FEBio v4.13 evidence for ringAddedMass (final bounded campaign)

Status: independent evidence, **draft; NOT fully verified, not a Category A claim.** All FEBio numbers are read from `runs/**/report.json` and the CSVs `RING_SPATIAL_STUDY.csv`, `RING_TEMPORAL_STUDY.csv`, `RING_COMPRESSIBILITY_STUDY.csv`, `RING_COMPUTATIONAL_COST.csv` (built by `tools/ring_build_csvs.py`); exact references from the official solids4foam JSON (commit 1576edcc). Pristine FEBio v4.13 (sha256 f5ae1ed6...) for every run, Slurm, host xenosim (shared).

## 0. Three different things are compared; do not conflate them
1. **Exact mathematical reference**: incompressible, inviscid added-mass solution (official JSON). 
2. **Closest physical-model comparison**: FEBio here solves the **free-slip inviscid limit** (slip outer wall and slip fluid-ring interface with w.n = 0, weak compressibility K = 1e8 Pa, regulariser E_mesh = 1e-3 Pa, hex8 solid, generalized-alpha rhoi = 1, monolithic). solids4foam solves a finite-viscosity (nu_f = 1e-7 m^2/s) **no-slip** problem with partitioned IQN-ILS (or Robin-Neumann, or quasi-monolithic) coupling and BDF2. Both approach the exact reference; they are not the same discrete problem and no strict discrete equivalence is claimed. FEBio's original no-slip runs (-6.0 % M1, -2.7 % M2 ratio error) are kept as model-form diagnostics.
3. **Solid discretisation differs**: FEBio hex8 vs solids4foam high-order (README tables) or "standard solid" (PR #13 raw CSVs). Absolute frequencies mostly measure the solid discretisation; the wet/dry ratio is the comparable quantity.

## 1. What was run in this campaign (12 of 13 authorised simulations; no duplicates; no M4)
| case | purpose | Slurm | cores | wall | allocated CPU-h |
|---|---|---|---|---|---|
| R1_strong_M1_K1e9 | Philip's K test (only K changed vs D2; `diff` of the .feb shows only `<k>` and the title comment) | 14001 | 2 | 184 s | 0.10 |
| R_MESH1_weak_M3, R_MESH2_strong_M3 | complete the three-level spatial study at K = 1e8 | 14002, 14003 | 8 + 8 | 4684 s, 4668 s | 10.4 + 10.4 |
| R_TIME_dry_M2_s50 / s25 | time study, dry | 14004, 14008 | 2 | 17 s, 10 s | 0.01 |
| R_TIME_{weak,moderate,strong}_M2_s50 | time study, wet | 14005-14007 | 4 each | 374 s (moderate) typical | 0.4 each |
| R_TIME_{weak,moderate,strong}_M2_s25 | time study, wet | 14009-14011 | 4 each | ~210-230 s | 0.25 each |
| R4_strong_M2_loosetol | iterative-error bracket (all stopping tolerances 100x looser, ftol 10x) | 14014 | 4 | 786 s | 0.87 |
Not run (by design): strong M2 K = 1e9 (the K effect was immaterial), any M4, any viscosity sweep. One authorised slot remains unused. Dry M1/M2/M3 were reused after regenerating their inputs with the current generator and diffing them (identical apart from comments). All 12 runs: all steps converged, normal termination, no Jacobian or max_refs events; every case is in `MULTIBENCHMARK_JOB_LEDGER.csv`.

## 2. Spatial study (K = 1e8 Pa, E_mesh = 1e-3 Pa, 100 steps/period): complete, three levels in every regime
| regime | mesh | FEBio omega_wet (rad/s) | FEBio omega err % | dry omega err % (same mesh) | **FEBio ratio err %** | s4f IQN-ILS high-order, README (unsigned) % | s4f IQN-ILS standard solid, PR #13 % | s4f quasi-monolithic standard solid, PR #13 % |
|---|---|---|---|---|---|---|---|---|
| dry | M1 | 2.621063 | +2.584 | - | - | -0.098 (omega err) | +20.66 (omega err) | same run (omega err) |
| dry | M2 | 2.571081 | +0.628 | - | - | -0.110 (omega err) | +6.17 (omega err) | same run (omega err) |
| dry | M3 | 2.558432 | +0.133 | - | - | -0.129 (omega err) | +1.58 (omega err) | same run (omega err) |
| weak | M1 | 2.498644 | +2.588 | +2.584 | **+0.0037** | 0.007 | n/a | n/a |
| weak | M2 | 2.451005 | +0.632 | +0.628 | **+0.0041** | 0.002 | n/a | n/a |
| weak | M3 | 2.438950 | +0.137 | +0.133 | **+0.0042** | 0.001 | n/a | n/a |
| moderate | M1 | 1.851254 | +2.590 | +2.584 | **+0.0061** | 0.038 | +0.0689 | -0.0158 |
| moderate | M2 | 1.815972 | +0.635 | +0.628 | **+0.0072** | 0.011 | +0.0181 | -0.0016 |
| moderate | M3 | 1.807043 | +0.140 | +0.133 | **+0.0075** | 0.003 | +0.0051 | -0.0001 |
| strong | M1 | 0.788796 | +2.620 | +2.584 | **+0.0355** | 0.068 | +0.1219 | -0.0312 |
| strong | M2 | 0.773788 | +0.668 | +0.628 | **+0.0399** | 0.021 | +0.0331 | -0.0035 |
| strong | M3 | 0.769990 | +0.174 | +0.133 | **+0.0410** | 0.006 | +0.0088 | -0.0004 |

Provenance: FEBio = this work. s4f "high-order IQN-ILS" = official README tables (printed values, **unsigned ratio errors**, high-order solid; the raw CSV is not in the repository). s4f "standard solid" = **raw PR #13 CSVs** (`paper-data/ringAddedMass_monolithic_{runs,mesh}.csv`, branch `ring-monolithic`; signed; ratio error given as a fraction and converted). The weak regime has no PR #13 mesh-study row. s4f Robin-Neumann: README gives omega error (mesh 2) only (-0.111 / -0.116 / -0.121 %), no ratio error. PR #13 standard-solid absolute errors (+20.7 / +6.2 / +1.6 % for the dry ring) are much larger than FEBio's hex8 (+2.6 / +0.63 / +0.13 %): FEBio's solid is not like-for-like with either s4f solid.

Observed behaviour (strictly from these data):
- **omega**: all three regimes and the dry ring converge at observed order 1.98 (three levels, hex8). M3 absolute errors: +0.137 / +0.140 / +0.174 % (weak/moderate/strong). The wet omega change between the last two meshes is 0.49 % (official 0.1 % criterion not met; absolute frequency is not mesh-converged to 0.1 %). Per Philip, the 0.3 % absolute criterion is reported but not used as pass/fail.
- **ratio error**: within the 0.1 % band at every level (largest +0.041 %). It **rises slowly and converges to a mesh-independent floor** (observed order 1.91-1.92 in all three regimes):
| regime | observed order of the ratio error M1-M3 | Richardson floor of the ratio error % | equivalent additive shift at the floor (1e-4 rad^2/s^2) |
|---|---|---|---|
| weak | 1.91 | +0.00427 | 5.07 |
| moderate | 1.92 | +0.00763 | 4.97 |
| strong | 1.91 | +0.04147 | 4.90 |
  So "decreases with refinement" is not met, and it is not a plateau of unknown extent: it is a converged offset. This is the answer to the earlier "floor of X %" question: X = 0.0043 / 0.0076 / 0.041 % (extrapolated, three-level, order about 1.9; an extrapolation, not a measurement).

## 3. Does the small ratio error come from cancellation of larger absolute errors? Yes at coarse meshes, no at the floor
| regime | mesh | wet omega err % | dry omega err % | wet - dry % | ratio err % | cancellation factor (wet err / ratio err) | equivalent additive shift of omega_wet^2 (2 omega_exact^2 eps, 1e-4 rad^2/s^2) |
|---|---|---|---|---|---|---|---|
| weak | M1 | +2.588 | +2.584 | +0.0038 | +0.0037 | 692.0 | 4.44 |
| weak | M2 | +0.632 | +0.628 | +0.0042 | +0.0041 | 153.0 | 4.90 |
| weak | M3 | +0.137 | +0.133 | +0.0042 | +0.0042 | 32.3 | 5.02 |
| moderate | M1 | +2.590 | +2.584 | +0.0063 | +0.0061 | 425.1 | 3.97 |
| moderate | M2 | +0.635 | +0.628 | +0.0073 | +0.0072 | 87.9 | 4.71 |
| moderate | M3 | +0.140 | +0.133 | +0.0075 | +0.0075 | 18.6 | 4.90 |
| strong | M1 | +2.620 | +2.584 | +0.0364 | +0.0355 | 73.9 | 4.19 |
| strong | M2 | +0.668 | +0.628 | +0.0401 | +0.0399 | 16.7 | 4.71 |
| strong | M3 | +0.174 | +0.133 | +0.0411 | +0.0410 | 4.2 | 4.85 |
- At M1/M2 the wet and dry absolute errors (0.6 to 2.6 %) are almost equal because they come from the same hex8 ring discretisation, so the ratio error is 17 to 690 times smaller than the absolute error. Small ratio errors on coarse meshes therefore partly **are** cancellation, and say nothing about the absolute accuracy of the solid.
- As the solid error falls (M3), the cancellation factor drops to 4 (strong) and the offset becomes visible as the wet-minus-dry difference.
- Equivalent additive shift: the ratio errors in all three regimes, at every mesh, correspond to a shift of omega_wet^2 of about **(4.0-5.1) x 1e-4 rad^2/s^2, independent of rho_f** (4.90-5.07e-4 at the extrapolated floors, spread +/-2 %). That is why the error is 10x larger in the strong case (omega_wet^2 is small there). This is a descriptive regularity of the data (`figures/ring_fig3_plateau_diagnostic.png`); **its physical or numerical cause has not been identified.**

## 4. Compressibility (Philip's request): strong M1, K = 1e8 -> 1e9 Pa, nothing else changed
| | K = 1e8 (D2) | K = 1e9 (R1) |
|---|---|---|
| omega_wet (rad/s) | 0.788795780 | 0.788796881 |
| wet/dry ratio error | +0.035455 % | **+0.035595 %** |
| n=2 pressure amplitude (Pa) | 1.38193 | 1.38194 |
| interface leakage rms | 1.232 % | 1.236 % |
| amplitude loss per period | 1.696e-5 | 1.696e-5 |
| Newton iterations per step | 4.06 | 4.01 |

- **The strong-case plateau does not drop tenfold; it does not move** (+0.00014 pp, i.e. a relative omega change of +1.40e-6, equal to the fit-window spread of 1.44e-6). The hypothesis that weak compressibility at K = 1e8 causes the plateau is **not supported**.
- Theory (linear acoustics in the actual annulus: Helmholtz solution, rigid wall at 1.5 m, n = 2, actual c = sqrt(K/rho_f) = 267 m/s at K = 1e8 for the strong case; `tools/ring_compress_theory.py`) gives an added-mass change of order (omega b/c)^2 x 0.15 = +2.8e-6 at K = 1e8 and +2.8e-7 at K = 1e9, i.e. a predicted ratio shift of -1.3e-6 (-0.00013 %) at K = 1e8. That is about 280 times smaller than the observed +0.0355 % and **of the opposite sign** (below the first acoustic resonance compressibility raises, not lowers, the effective added mass in this calculation). Measured vs predicted omega shifts agree in size: strong K 1e8 -> 1e9: +1.40e-6 measured (floor 1.4e-6) vs +1.14e-6 predicted; moderate M1 K 1e8 -> 1e7 (Phase C, E_mesh = 1): -4.01e-6 measured (floor 2e-7) vs -3.47e-6 predicted (`figures/ring_fig4_compressibility.png`). K is therefore a ~1e-6 effect in omega here and cannot be the plateau. A strong M2 K = 1e9 run would add nothing and was not run.

## 5. Time-step study (M2; 25 / 50 / 100 steps per period; dt from each case's exact period)
| regime | steps/period | FEBio omega err % | **FEBio ratio err %** | FEBio loss %/period | s4f IQN-ILS high-order ratio err % (README, unsigned) | s4f IQN-ILS std ratio err % (PR #13) | s4f quasi-mono std ratio err % (PR #13) |
|---|---|---|---|---|---|---|---|
| dry | 25 | +0.129 | - | -9.5e-08 | - | - | - |
| dry | 50 | +0.527 | - | 1.0e-07 | - | - | - |
| dry | 100 | +0.628 | - | 8.6e-08 | - | - | - |
| weak | 25 | +0.134 | **+0.0040** | 4.0e-05 | 0.013 | +0.0151 | -0.0003 |
| weak | 50 | +0.531 | **+0.0041** | 4.1e-05 | 0.009 | +0.0103 | -0.0004 |
| weak | 100 | +0.632 | **+0.0041** | 4.1e-05 | 0.002 | +0.0036 | -0.0002 |
| moderate | 25 | +0.137 | **+0.0071** | 3.2e-04 | 0.069 | +0.0753 | -0.0034 |
| moderate | 50 | +0.535 | **+0.0072** | 3.2e-04 | 0.048 | +0.0565 | -0.0025 |
| moderate | 100 | +0.635 | **+0.0072** | 3.2e-04 | 0.011 | +0.0181 | -0.0016 |
| strong | 25 | +0.169 | **+0.0391** | 1.8e-03 | 0.125 | +0.1262 | -0.0070 |
| strong | 50 | +0.567 | **+0.0397** | 1.8e-03 | 0.087 | +0.1010 | -0.0052 |
| strong | 100 | +0.668 | **+0.0399** | 1.8e-03 | 0.021 | +0.0332 | -0.0035 |
- **omega converges at observed order 1.987** in all four cases (dry, weak, moderate, strong; three levels, outside the noise floor of 1e-6), i.e. second order as expected for generalized-alpha with rhoi = 1 (not BDF2; no claim that FEBio must reproduce BDF2's phase error). Richardson dt -> 0 for the dry ring: omega = 2.57195 (+0.662 % vs exact), so at 100 steps the temporal error is -0.034 % and at 25 steps -0.53 %. The low-looking dry error at 25 steps (+0.13 %) is an accidental cancellation of the positive spatial and the negative temporal error, not accuracy.
- **The ratio error is almost independent of dt** (strong: +0.0391 -> +0.0399 %, a change of +0.0008 pp from 25 to 100 steps; Richardson dt -> 0: +0.0399 %). The temporal error cancels in the ratio, as designed (omega dt matched), and **the plateau is not temporal**. In contrast the solids4foam IQN-ILS ratio error falls strongly with steps per period (strong standard solid: 0.126 -> 0.101 -> 0.033 %), and the quasi-monolithic one is about -0.0035 to -0.007 %, ten times below FEBio's +0.040 % in the strong case.
- Amplitude loss in FEBio: 4e-5 %, 3e-4 %, 1.8e-3 % per period (weak/moderate/strong), independent of dt; solids4foam IQN-ILS (README, high-order solid) loses 2.3 / 2.9 / 3.4 % per period at 25 steps and 0.15 / 0.64 / 1.14 % at 100 steps (mesh 2).

## 6. Iterative solver error (R4)
Earlier tightening (100x) failed on arithmetic floors (D6). Here the perturbation goes the other way and cannot hit a floor: strong M2 with dtol = vtol = 1e-4, etol = 1e-6, rtol = 1e-8 (100x looser), ftol = 1e-5 (10x looser). Newton iterations drop from 4.55 to 4.00 per step (1819 -> 1601 residual evaluations), yet **omega_wet is identical to D3 in all printed digits (0.773787885179203; relative difference 0.0), the n=2 force and pressure agree to 1e-12**. So the accepted-residual (stopping-rule) error is orders of magnitude below any quantity discussed here: under 1e-12 relative on omega even after loosening, with the usual caveat that this is one regime at one mesh and the dilatation increment floor (~1e-15) cannot be probed. It does not address spatial or temporal error.

## 7. Cost (descriptive; machines, cores and solvers differ)
Wall/CPU/RSS for every FEBio run are in `RING_COMPUTATIONAL_COST.csv`. Moderate M2/M3 on 1 core: 1597 s / 11538 s (3.2 core-h), 357 / 1405 MB. Weak/strong M3 on **8 cores: 4684 s / 4668 s, 10.4 allocated core-h each at 40 % CPU utilisation** (the 1-core moderate M3 used 3.2 h): 8 cores bought a 2.5x speed-up for 3.3x the allocated CPU-hours; 4 cores would have been the better choice (consistent with the Phase A scaling, efficiency 0.3-0.7 beyond 4 cores). R1 on 2 cores: 184 s for 298 s on 1 core (1.6x). s4f README times (serial, OpenFOAM v2412, same host): mesh 2, 100 steps: IQN-ILS 125 / 185 / 406 s, Robin-Neumann 156 / 213 / 575 s (weak/moderate/strong); PR #13 raw execution times: IQN-ILS 51 / 73 / 120 s, quasi-monolithic 86 / 89 / 91 s (m2, n100). FEBio at M2, 1 core: 1756 / 1597 / 2033 s: 5-14x slower than the README IQN-ILS times and 20-40x slower than the PR #13 raw times (different loads and hardware state), with 1.4 GB peak memory at M3 (s4f memory not recorded). Not a fair speed comparison (hex8 vs high-order, monolithic with direct Pardiso vs segregated).

## 8. R5: like-for-like free-slip formulation in solids4foam (read-only feasibility audit; nothing was modified or run)
solids4foam ringAddedMass (development, 0/fluid/U): interface `newMovingWallVelocity` (no-slip moving wall), `outerWall` `fixedValue (0 0 0)` (no-slip), symmetry on the cut planes, `pimpleFluid`, p `zeroGradient`. A matched free-slip variant looks **straightforward for the outer wall** (`slip`, standard) and **plausible for the interface** (`movingWallSlipVelocity`, present in OpenFOAM 10; its flux consistency with solids4foam's mesh-velocity handling and with the traction transfer in IQN-ILS/Robin-Neumann must be checked, and solids4foam's own `newMovingWallVelocity` would have to be replaced or extended). It requires a working solids4foam + OpenFOAM build and a new verification run set (a solver change on Philip's side), so it was **not attempted**. Value: it would separate "slip vs no-slip / finite viscosity" from "solver" in the +0.04 % offset. Recommendation: only if Philip wants it; the accepted comparison remains against the exact inviscid reference.

## 9. Official-criteria scorecard for FEBio (`results/phaseD_acceptance_checks.csv`)
Pass (all regimes): finest ratio error <= 0.1 %; ratio vs thin-ring <= 0.25 %; amplitude loss and growth criteria; **time order >= 1.5 (1.987, new)**; finest omega <= 0.3 % (M3: +0.137 / +0.140 / +0.174 %, reported not decisive). Fail / not meaningful: ratio error decreases with refinement (converged floor instead, section 2); omega change between last two meshes <= 0.1 % (0.49 %; absolute frequency not mesh-converged to that level for a hex8 solid).

## 10. Uncertainties and what is NOT resolved
- **The cause of the +0.0043 / +0.0076 / +0.041 % floors is unknown.** Excluded by measurement: bulk modulus (section 4), time step (section 5), mesh (section 2), stopping tolerances (section 6), regulariser (E_mesh 1 -> 1e-3 Pa removed 0.164 pp; linear scaling leaves about 2 % of the moderate floor). Not excluded: free-slip model-form (weakly enforced w.n = 0; leakage 0.24-1.2 % rms, falling with refinement but not to zero), nodal-traction interface force (-0.53 % at M3, first order), initial-condition/start-up content, amplitude-dependent (nonlinear) effects, FSI formulation details. Two cheap discriminating runs are proposed but were **not run** (outside the authorised matrix): halve the initial amplitude (a nonlinear origin would scale the offset by 1/4) and E_mesh 1e-3 -> 1e-5.
- Weak M3 had one step needing 16 Newton iterations (still converged; dilatation margin 0.999 at the worst step), and the dilatation increment sits at its arithmetic floor (margins up to 0.97-0.999 at M3).
- Three mesh levels give an observed order only for omega (1.98) and the ratio-error offset (1.91): consistency checks, not asymptotic proofs; the extrapolated floors rely on them.
- 8-core runs (M3 weak/strong) used a different thread count from the 1-core runs (2- and 4-core for others): round-off level differences (1e-13) only, as the R4 and D3 comparison (4 vs 1 core) shows.
- solids4foam columns: high-order values are printed README values (unsigned ratio error); standard-solid values are raw PR #13 data. No solids4foam run was repeated; no Ansys/COMSOL comparison exists for this benchmark.
- Not fully verified. Remaining model-form difference (slip vs no-slip) and the unexplained floor keep this at "independent evidence with stated limitations".

## 11. Figures (all inspected visually)
`figures/ring_fig1_spatial.png` (absolute errors and ratio error vs mesh, with s4f raw data), `ring_fig2_temporal.png`, `ring_fig3_plateau_diagnostic.png`, `ring_fig4_compressibility.png`, `ring_fig5_cost.png`; updated `fig_D1`-`fig_D5` (now with weak/strong M3).
