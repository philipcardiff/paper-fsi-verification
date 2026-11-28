# solids4foam learnings register

The verification campaign keeps finding subtle choices in solids4foam, each one a fix, an option or a set-up detail, that change results. This register records each finding while it is fresh.

**Why it exists.** After the papers are submitted, each entry has to reach solids4foam in a robust form, so that a user who has not read the papers cannot pick the wrong option by accident. That means one of:

- a corrected default;
- a warning or a fatal error on the unsafe combination;
- a fixed tutorial;
- as a last resort, documentation.

The solver paper (paper-monolithic-fsi) keeps a matching register for the quasi-monolithic coupling. Entries that affect both couplings are cross-referenced.

**How to add an entry.** Give it the next free number, and never renumber entries. Fill in every field, citing evidence by issue, PR or paper-data path. Update the *State* field when solids4foam changes.

**Status values:**

- **Open:** nothing in solids4foam yet.
- **Tutorial only:** fixed in the affected tutorials; an uninformed user can still hit it in a new case.
- **Issue filed:** a solids4foam issue exists.
- **Robust:** a default, check or error in the code prevents it.
- **Done:** robust and documented.

## Summary

| # | Topic | Area | Effect if wrong | Status |
|---:|---|---|---|---|
| 1 | Continuity-flux interpolation in `pimpleFluid` | partitioned fluid | 0.05% F_y bias, slowly converging (cavity) | Tutorial only (#743, open) |
| 2 | `fixedFluxPressure` on a moving FSI wall | partitioned fluid | silently zeroGradient | Open |
| 3 | Interface pressure condition | partitioned fluid | O(h) wall pressure; transients | Issue filed (#721) |
| 4 | Interface viscous traction form | both couplings | O(1) on a deforming wall | Open |
| 5 | Wall-normal velocity gradient at moving walls | partitioned fluid | 1.5–3.8% in FSI3 | Robust (#671) |
| 6 | Default coupling convergence test | coupling | accepts unconverged steps | Issue filed (#552) |
| 7 | Strict test at a true steady state | coupling | runs to nOuterCorr, false warning | Issue filed (#710) |
| 8 | First coupling iterate without a predictor | coupling | late added-mass crashes | Robust (#492) |
| 9 | Unrelaxed PIMPLE outer loop at large νΔt/h² | partitioned fluid | divergence (Mok 128²) | Open (#700) |
| 10 | Mesh families with cell-size jumps | meshing | spoils convergence studies | Open |
| 11 | Wall arcs on curved two-region meshes | meshing | ~1 µm AMI wall mismatch | Tutorial only |
| 12 | cfMesh octree aligned with interface faces | meshing | bulged flap, bad cells | Open |
| 13 | `#codeStream` static cache on macOS | tutorials | NaN or FPE, non-deterministic | Issue filed (#711) |
| 14 | Units of the interface-displacement log line | output | misread as metres | Open |
| 15 | Benchmark definitions in tutorials | tutorials | wrong reference comparison | Open |
| 16 | Steady state depends on the pseudo-time step | partitioned fluid | up to 0.2 % on coarse meshes | Open |
| 17 | Boundary non-orthogonal correction is zero on OpenFOAM.com | moving-wall and other BCs | viscous traction error on skewed wall cells | Issue filed (#748) |
| 18 | First-order pressure in `pimpleFluid` on a manufactured FSI solution | partitioned fluid | p order ~1, outlet–wall corner error | Issue filed (#750) |
| 19 | Robin–Neumann restart is not exact | restart / coupling | ≈130 Pa pressure jump; 0.4 % displacement perturbation (3dTube) | Issue filed (#782) |

## Entries

### 1. Continuity-flux interpolation in `pimpleFluid`

- **Found:** 2026-11-20, cavityFlexibleBottom. Dylan Armfield's SU2 + CalculiX comparison (pfv#102) showed our partitioned F_y limit 0.054% from SU2, FEniCSx and the quasi-monolithic solution, outside every GCI.
- **What goes wrong:** `pimpleFluid` builds the continuity flux from `fvc::flux(HbyA)`. That looks up `flux(HbyA)` in `interpolationSchemes` and falls back to `default linear`, without skew correction.
  - The `interpolate(U) skewCorrected linear` entry in the tutorials reaches only `Uf`, not the continuity flux.
  - On a strongly deformed mesh, here the plate ends of the cavity, the error converges slowly (about 1.5× per halving at the finest levels). Four levels then give a falsely tight GCI.
  - The quasi-monolithic solver has interpolated its continuity flux with skew correction since #597.
- **Evidence:** /Users/philipc/claude-runs/cavgap, build 206437983, strict Aitken.
  - Adding only `flux(HbyA) skewCorrected linear;` moves the levels 1–3 F_y limit from −5.26050 (+0.052% from SU2) to −5.25695 (−0.016%; GCI 0.14%).
  - The quasi-monolithic levels 1–3 limit is −5.25742.
  - Agent B's catalogue (paper-monolithic-fsi `paper-data/a2_discretisation_differences.md`, D4) shows the quasi-monolithic and partitioned results agree to below 0.001% when both use linear.
  - Level 4 confirms it. The levels 2–4 limits are −0.240583 m / −5.257947 N/m (u_y / F_y; orders 1.98 / 1.86, GCI 0.31% / 0.045%), against −0.240671 / −5.260619 with `linear`. They are 0.010% / 0.003% from SU2, 0.016% / 0.003% from FEniCSx and 0.001% / 0.004% from the quasi-monolithic #618 solution.
  - Level 4 used `nOuterCorr 40` and fluid linear-solver tolerances of 1e-8, as the earlier level 4 did. Data: `paper-data/cavity_partitioned_skew/`.
- **Right choice:** skew-corrected interpolation of the continuity flux on any mesh that moves or deforms.
- **State:** Tutorial only (solids4foam #743, open). The paper's cavity sequence uses the skew-corrected flux.
- **Robust fix to consider:** make `pimpleFluid` use the `interpolate(U)` scheme for `flux(HbyA)` when no explicit entry exists, or warn when a moving mesh runs on the `default linear` fallback. Check the other moving-mesh tutorials, such as FSI3 and the 3-D tube.

### 2. `fixedFluxPressure` on a moving FSI wall

- **Found:** 2026-11-20, Seevani Bali's preCICE FSI3 (pfv#104).
- **What goes wrong:**
  - On a wall where velocity is fixed, `constrainHbyA` sets HbyA_b = U_b, and the ddt correction is zero there. `constrainPressure` therefore receives a zero numerator, and `fixedFluxPressure` gives dp/dn = 0.
  - A user expecting the normal-momentum gradient −ρ n·a gets zeroGradient without being told. Her "base" and "zerograd" runs agree to 2e-5.
  - This applies to `pimpleFoam` and to `pimpleFluid`, which calls `constrainPressure`.
- **Right choice:** `movingWallPressure` for the acceleration term (see entry 3), or an explicit zeroGradient where that is the intent.
- **State:** Open.
- **Robust fix to consider:** a warning from the FSI interface set-up when the fluid interface pressure is `fixedFluxPressure`; documentation in the FSI tutorials' README.

### 3. Interface pressure condition

- **Found:** October–November 2026, ringAddedMass, Hessenthaler, womersleyTube and the solver-paper comparison.
- **What goes wrong:**
  - `zeroGradient` gives an O(h) wall-pressure error, d·dp/dn. On a moving wall it misses −ρ n·a, and on any wall the viscous part μ ∂²u_n/∂n².
  - `movingWallPressure` adds the acceleration only.
  - The quasi-monolithic solver uses the full normal-momentum gradient, which includes the viscous terms.
  - `movingWallPressure` fails in enclosed fluid domains (ringAddedMass with IQN-ILS: "Continuity error cannot be removed"). Under IQN-ILS it can also stall temporal convergence at small Δt.
- **Right choice:** Philip's standard since 2026-11-17 is `movingWallPressure` for every partitioned run, with the omitted viscous term stated in the text.
- **State:** Issue filed: solids4foam#721 (viscous term, enclosed domains, removal of the hard `refCast`).
- **Robust fix to consider:** the #721 proposal, then make it the FSI tutorials' default.

### 4. Interface viscous traction form

- **Found:** solver-paper comparison of the couplings (paper-monolithic-fsi `a2_discretisation_differences.md`, D6), womersleyTube.
- **What goes wrong:**
  - Partitioned passes `n & (-devReff)`; quasi-monolithic passes `rho*nuEff*snGrad(U)`.
  - On a deforming wall they differ by μ grad(u·n), a physical term that does not vanish under refinement. In the Womersley tube it is about 3e-4 relative.
  - At a stationary wall the difference is O(h).
  - A related rule (#583, rejected): the traction must be the fluid's own discrete wall flux. Quadratic wall-shear corrections degrade the partitioned convergence.
- **Right choice:** not settled; the two couplings should use the same, correct, conservative form.
- **State:** Open. Owned by the solver paper.
- **Robust fix to consider:** a shared traction routine for both couplings, tested on the Womersley tube.

### 5. Wall-normal velocity gradient at moving walls

- **Found:** November 2026, FSI3. The corrected sequence moved by 1.5–3.8% (u_x −3.1%, u_y amplitude −1.5%, drag amplitude −3.8%, lift amplitude −2.2% at 2x′).
- **What goes wrong:** before #671, `newMovingWallVelocity` and `elasticWallVelocity` set the wall-normal gradient of u_n inconsistently with continuity.
- **State:** Robust. #671 made the continuity-based gradient the default. Every paper result computed with these conditions before #671 is superseded or caveated (pfv#93, pfv#97).
- **Remaining:** none, except checking that no tutorial pins the old behaviour.

### 6. Default coupling convergence test

- **Found:** October 2026, IQN-ILS steady and quasi-steady phases.
- **What goes wrong:** the default test accepts a step when min(r1, r2) is below the tolerance. Near steady state, r2 is tiny from the first iteration, so a one-iteration step is accepted and the solution creeps.
- **Right choice:** `requireAllResidualMeasures yes`, the strict test.
- **State:** Issue filed: solids4foam#552.
- **Robust fix to consider:** make the strict test the default, once entry 7 is fixed.

### 7. Strict test at a true steady state

- **Found:** 2026-11-14, partitioned Hessenthaler rF2.
- **What goes wrong:** at round-off steadiness, r1 (relative to the step's largest residual) cannot fall below 1e-4. Every step then runs to `nOuterCorr` and prints "The resulting solution may be invalid", although it is converged. This costs 4–6× the coupling iterations.
- **State:** Issue filed: solids4foam#710, which proposes an absolute floor for r1.
- **Robust fix to consider:** the #710 floor; needed before entry 6's default can change.

### 8. First coupling iterate without a predictor

- **Found:** October 2026, HronTurekFsi3 late crashes (#383, #489).
- **What goes wrong:** the relaxed first iterate without a predictor causes an added-mass spike.
- **State:** Robust. `predictor` defaults to `yes` since #492.

### 9. Unrelaxed PIMPLE outer loop at large νΔt/h²

- **Found:** 2026-11-10, Mok lid-driven cavity at 128².
- **What goes wrong:** with no PIMPLE relaxation, the fluid outer loop diverges inside a frozen-interface coupling iteration (νΔt/h² = 4.1 at 128², against 2.3 at 96²). It looks like a coupling failure, but it is not one.
- **Right choice:** p 0.3 / U 0.7 relaxation with residual control, as in run `m128_relax2`. It reaches t = 70 within 0.5% of Kratos.
- **State:** Open. Documented in solids4foam#700, but the tutorial `fvSolution` is unchanged until it is tested on the coarse levels.
- **Robust fix to consider:** relaxed defaults in the tutorial, and a log hint when the p residual grows across outer correctors.

### 10. Mesh families with cell-size jumps

- **Found:** October 2026, the Hron–Turek ring+box family.
- **What goes wrong:**
  - Ring+box: neighbour size ratios of 2–2.8:1 through the wake on every level, so refining does not remove the jump error and the observed order is unreliable.
  - FSI3 on the old family also showed a first-order issue.
- **Right choice:** audit the neighbour size ratios before any convergence study. FSI3 now uses agent B's Gmsh family.
- **State:** Open (practice, not code).
- **Robust fix to consider:** a small utility, or a `checkMesh` function object, that reports the distribution of neighbour size ratios, used in the verification drivers.

### 11. Wall arcs on curved two-region meshes

- **Found:** November 2026, 3-D pressure-pulse tube.
- **What goes wrong:** blockMeshDict arc points written to too few digits put the fluid and solid walls about 1 µm apart, so the AMI interface sees a gap.
- **State:** Tutorial only. solids4foam 9b52c2303 writes the arcs at full precision, and the paper uses the corrected meshes.
- **Robust fix to consider:** a check in the FSI interface set-up that reports the largest fluid–solid interface face-centre gap at start-up.

### 12. cfMesh octree aligned with interface faces

- **Found:** 2026-11-10, Hessenthaler feF1.
- **What goes wrong:** cfMesh centres its octree on the surface bounding box. When the flap faces land on octree planes, they bulge by 0.02 mm, and the strict `checkMesh` finds concave cells and bad tets.
- **Right choice:** shift the octree phase by adding one unused point to the `.fms`, then check with `checkMesh -allGeometry -allTopology`.
- **State:** Open (meshing practice).

### 13. `#codeStream` static cache on macOS

- **Found:** 2026-11-14, womersleyTube.
- **What goes wrong:** a function-local static inside an inline function is shared across `#codeStream` libraries. The result is NaN or an FPE with identical inputs, varying from run to run.
- **State:** Issue filed: solids4foam#711.

### 14. Units of the interface-displacement log line

- **What goes wrong:** "Maximal accumulated displacement" is normal motion × deltaCoeffs, in cell widths, not metres. Values of 3–6 are normal; it was misread once (#687, retracted).
- **State:** Open.
- **Robust fix to consider:** state the unit in the log line.

### 15. Benchmark definitions in tutorials

Set-up details that cost time because the tutorial or its reference did not say them:

- **Cavity force depth:** the published forces are for a 0.05 m depth, and the tutorial is 1 m deep.
- **Cavity inflow profile:** parabolic, with a peak of 1.5 m/s (mean 1 m/s). Neither the tutorial README nor, until 2026-11-23, the paper stated it, so an independent COMSOL solution used a uniform inflow of the same mean and came out about 2 % high in both u_y and F_y (pfv#88).
- **Cavity plate law:** the Tuković et al. plate law is not confirmed; it may be neo-Hookean rather than StVK.
- **collapsibleChannel:** needs the Poiseuille initial flow (`setPoiseuilleFlow`). The oomph-lib reference is first order in `--refine`.
- **FSI3:** ring+box meshes as in entry 10. `linearUpwind grad(U)` convection and the coupling tolerances of the shared set-up.

**State:** Open. Each tutorial README should state these.

### 16. Steady state depends on the pseudo-time step

- **Found:** 2026-11-21, cavityFlexibleBottom. This was the follow-up to entry 1.
- **What goes wrong:** a pseudo-transient steady solution from `pimpleFluid` depends on Δt. The Rhie–Chow (pressure–velocity coupling) term carries the coefficient rAU, which depends on Δt.
  - On level 1, reducing Δt from 0.5 to 0.0625 s changes F_y by 0.18 % with skew correction and by 0.07 % without.
  - The dependence vanishes rapidly under refinement: 0.03 % on level 2 and 0.001 % on level 3.
  - The solver tolerances do not matter (≤ 0.02 %).
  - A study that scales Δt with h mixes this into the mesh convergence. Here it is harmless: a fixed-Δt sequence extrapolates within 0.003 %.
- **Evidence:** `paper-data/cavity_partitioned_skew/README.md`, section "Δt check".
- **Right choice:** in a steady study, state Δt per level and check one level at a second Δt. On coarse meshes, a "steady" result is not Δt-independent.
- **State:** Open (practice).
- **Robust fix to consider:** a note in the steady FSI tutorials' READMEs. Longer term, a Δt-independent Rhie–Chow coefficient for pseudo-transient steady runs (as in steady SIMPLE), or a steady fluid model for steady FSI.

### 17. Boundary non-orthogonal correction is zero on OpenFOAM.com

- **Found:** 2026-11-21, by solver-paper agent B on the manufactured-solution FSI case (solver-paper register entries 28–32).
- **What goes wrong:** several conditions build their correction vector as k = (I − nn)·`patch().delta()`:
  - the moving-wall and interface velocity conditions, `newMovingWallVelocity` and `elasticWallVelocity`;
  - the inlet profiles;
  - `solidSymmetry`, `normalDisplacement`, `tractionPressureDisplacement` and `fixedTemperature`.

  On OpenFOAM.com, `fvPatch::delta()` is purely normal for non-coupled patches, so k ≡ 0 and the correction does nothing. On foam-extend it does correct, so the forks differ on skewed wall cells.
- **Effect:** an uncorrected wall-normal velocity gradient gives a viscous-traction error of about ν tanθ |∂U/∂τ|, which does not vanish under refinement if the wall angle θ does not. Where it matters:
  - it is zero where the wall velocity does not vary along the wall, e.g. at a steady state (cavity, Hessenthaler);
  - it is small where the near-wall cells are orthogonal;
  - it is largest on skewed near-wall cells of a deforming wall: FSI3 around the flag, the collapsible channel.
- **Right choice:** `patchCorrectionVectors(patch())`, as `fixedDisplacement`, `solidTraction` and `fixedGradientCorrected` already use.
- **State:** Issue filed: solids4foam#748. Once a fix exists, run an A/B on the transient partitioned cases with non-orthogonal wall cells, FSI3 and the collapsible channel first.
- **Robust fix to consider:** the #748 fix, plus a unit test that the correction is non-zero on a skewed patch on both forks.

### 18. First-order pressure in `pimpleFluid` on a manufactured FSI solution

- **Found:** 2026-11-21, manufactured-solution FSI case (solver-paper agent B).
- **What goes wrong:** with the partitioned couplings, fluid p converges at about first order over the whole domain. Fluid U converges at about 1.3, and there is an outlet–bottom-wall corner error with an L∞ order of 0.4–0.8. Meanwhile the solid is second order, and the quasi-monolithic fluid gives a pressure order of 1.7–2.1 on the same meshes.
- **Ruled out (2026-11-21, `/Users/philipc/claude-runs/mms-pfv`):** skew-corrected continuity flux (entry 1), `movingWallPressure` on the interface (entry 3), and the exact pressure gradient on the bottom wall. Each changes the fluid orders by at most 0.1.
- **Lead (solver-paper agent B, 2026-11-21, solids4foam #750 comment 7063555399 and #752):** `pimpleFluid`'s `ddtCorr`, on by default, accounts for much of the first-order pressure. With Aitken on n16–n64, Δt refined with h:
  - default: p L2 orders 0.88 and 0.94;
  - unlimited `backward 1`: 1.13 and 1.40;
  - `ddtCorr no`: 1.73 and 1.49, with p L2 at n64 2.4× lower and U L∞ 3.4× lower.

  This is not yet a recommendation, because Δt and h were refined together; separate fixed-Δt and fixed-mesh studies are needed. Related: entry 16, which is also a Rhie–Chow/Δt effect, and solver-paper register entry 33 (time-history-consistent Rhie–Chow for QM, #752).
- **Also open:** the partitioned solid and interface displacements are first order at t = 0.2 for both couplings, while QM is second order. From n64 to n128 at t = 0.4 their order drops to 0.7–0.8 (paper-data/mms_fsi). The interface pressure condition does not change this.
- **State:** Issue filed: solids4foam#750. Diagnostic variants are requested: fluid only with prescribed interface motion, and an outlet without the corner clash.
- **Robust fix to consider:** depends on the cause.

### 19. Robin–Neumann restart is not exact

- **Found:** 2026-11-27, 3dTube, restart test for the L4 run (`/Users/philipc/claude-runs/drafts/tube_L4_and_stvk.md`, section A).
- **What goes wrong:** a Robin–Neumann run restarted from a written time does not continue the uninterrupted solution.
  - L1 was written every 1 ms, stopped at 3 ms and resubmitted to 6 ms, with the L4 controls and job script (`restart yes` in the solid coefficients), and compared with a continuous L1 run.
  - Up to the restart the two runs agree bit for bit, even on a different node type.
  - In the first step after the restart the axis pressure at z = 25 mm jumps by ≈ 130 Pa (102.7 → −28.2 Pa, against 105.1 Pa continuous). Up to 6 ms the axis-pressure probes differ by 10–70 % of their running peak.
  - Point-A u_r differs by up to 0.4 % of its running peak (u_z by 0.3 %); at 6 ms u_r is 0.38 % low. That is comparable to the L2→L3 change of u_r,max (0.59 %).
- **Suspected cause (not confirmed in the code):** the Robin wall condition (`elasticWallPressure`) rebuilds its previous-step interface acceleration at every step from the solid's face acceleration. The solid writes `U` but not its old-time fields (`U_0`), so after a restart that acceleration starts from an incomplete old-time state. The Robin secant history is not written either, but that affects only convergence.
- **Consequence:** run verification cases in one segment. Restart directories are a fallback only; a restarted run carries a perturbation of the order of the finest mesh changes and must be flagged. A restart at a quieter time than the inlet switch-off may perturb less, but it will not be exact.
- **Evidence:** `sonic:/scratch/pcardiff/tube_stvk/L1rst/` (`prepare_rst.py`, `compare_restart.py`).
- **State:** Issue filed: solids4foam#782.
- **Robust fix to consider:** write and read the old-time solid fields the Robin condition needs (or the Robin acceleration itself, and the secant history), and add a restart regression test that compares a stop-and-restart run with a continuous one; until then, a warning when a Robin–Neumann run starts from a non-zero time.

### 20. 3dTube driver: c_p should be measured on the averaged axis pressure; law effects need resolved meshes

- **Found:** 2026-11-28, 3dTube, axisymmetric (1° wedge) sequence and c_p recomputation (`paper-data/3dtube_axisymmetric/`).
- **What goes wrong (c_p):** the verification driver of `tutorials/fluidSolidInteraction/3dTube` takes c_p from the first crossing of half the inlet pressure by the raw axis-pressure probe signal. The axis pressure rises slowly and carries a spurious oscillation of ≈ 0.2 ms period, whose amplitude near the front grows as numerical damping falls with h and Δt (2–24 Pa on L1, 10–40 Pa on 3-D L3). The crossing therefore drifts: the raw c_p changes by −0.39, −0.54, −0.81 % per level on the wedge and by −0.86 % when Δt is halved on 3-D L3.
- **Fix used in the paper:** apply a centred moving average over 0.4 ms (two periods) to the probe signal before the crossing. The averaged c_p changes by ≤ 0.06 % when Δt is halved on 3-D L2/L3, by ≤ 0.05 % between 0.2 and 0.4 ms windows on L3/L4, and converges to about 0.1–0.2 %. Script: `paper-data/3dtube_axisymmetric/scripts/cp_averaged.py`.
- **State:** follow-up solids4foam change after submission: make the 0.4 ms averaged c_p the driver's c_p (keep the raw value as a diagnostic column), and update the tutorial's reference JSON and README.
- **Law-effect lesson:** a comparison of two material laws on one coarse mesh can overstate the law effect. On the 3dTube L1 the Hookean–St Venant–Kirchhoff difference in u_r,max was −1.1 %; on resolved meshes it is −0.4 % (wedge L4) and −0.5 % (3-D L3), because on L1 the mesh error and the law effect mix. Measure model-form effects on at least two resolved levels before quoting them.
- **Robust fix to consider:** the driver change above, plus a regression check that c_p changes by less than ≈ 0.1 % when Δt is halved on L2.
