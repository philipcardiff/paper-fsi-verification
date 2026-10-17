# Hessenthaler Phase I verification findings

Evidence report for the Hessenthaler et al. (2017) Phase I benchmark
(`sec:hessenthaler` in `sectionBenchmarkSuite.tex`) before any manuscript
change. It addresses item 7 of `manuscript_update_2026-10.md`: strict-test
verification runs from t = 0, the PR #546 A/B check and the justification
of μ*. No `.tex` file is touched.

Compact data:
- `hessenthaler_phaseI_verification_steady.csv`: steady tip of every strict
  run, with tail window and drift.
- `hessenthaler_phaseI_verification_final.json`: the convergence triplets,
  the additivity check and the wall diagnostic.
- `hessenthaler_phaseI_verification.json`: the restart tests.

They are generated from the tip histories of the runs listed in Section 2.

## 1. Executive scientific conclusion

- **What is verified.** Phase I is used as a solution-verification case: one
  fixed mathematical model, solved on refined meshes and time steps. The
  MRI data enter only the validation (Section 13).
- **Fixed model.** μ* = 64.2 kPa (ν = 0.45) is held fixed for every run. It
  comes from the continuum limit of a mesh-converged zero-flow calibration
  (Section 4) and is never recalibrated per mesh.
- **Main method finding (Section 9).** The default solids4foam FSI
  convergence test let an iterative error accumulate that looked like
  physics.
  - The test accepts a step when the smaller of two residual measures is
    below the tolerance. One of them is normalised by the total interface
    displacement (tip about 17 mm).
  - After the first second, every step was accepted with the in-step
    residual r1 above the tolerance (median 6.8e-2).
  - The resulting "creep" (+0.003 to +0.005 mm/s, still present at 80 s)
    exceeded every spatial difference and distorted the transient.
  - Requiring both measures (the strict test) removes it. All results below
    are strict-test runs from t = 0, run to a steady state.
- **Solid refinement.** At fixed F1 and material parameters, the tip shows
  monotone solid-mesh convergence consistent with approximately second
  order (p = 2.21). The conditional GCI estimate on S3 is 0.034 mm
  (0.20 %).
- **Fluid refinement.** Asymptotic convergence has not been established.
  - F1 → F2 → F3 shows a monotone but slow decrease of the differences
    (ratio 0.81, p_obs = 0.62, formal order 2).
  - A provisional Richardson-based GCI indicator (p_obs, Fs = 1.25) is
    0.16 mm. It dominates the budget.
  - Upstream flow quantities converge at about second order. Refining only
    the walls explains about 20 % of the fluid effect on the tip.
- **Temporal, iterative and end-state effects are small.** The Δt
  sensitivity of the steady tip is ≤ 0.004 mm (0.02 %). The coupling
  tolerance changes it by ≤ 0.0004 mm, both at the end state and along the
  path. After about 40 s the remaining drift is < 1e-4 mm/s.
- **Best-resolved solution.** F3S3 gives a tip y of 16.8133 mm (16.819 mm
  after the PR #546 correction, Section 15).
  - Provisional linear aggregate allowance: about 0.20 mm (1.2 %), dominated
    by the fluid discretisation.
  - This is a conditional estimate, not a guaranteed bound or a confidence
    interval.
- **Validation.** F3S3 lies 0.40 mm above the MRI tip, and every refinement
  moves the result further away (S1 −0.21 mm → F3S3 +0.40 mm).

## 2. Provenance

- **solids4foam.** Branch `verification/hessenthaler-phaseI`, draft PR #585,
  final commit `f7f1489e6`. It is based on the PR #535 head `4f549b500`,
  and `development` was merged in after #535 merged. `src/` is unchanged.
  The compact results are in `verification/results/`. Verification driver
  (`tutorials/fluidSolidInteraction/hessenthalerFsi/verification/`):
  - spec keys `fsitest=strict` (`a2ccaf64c`), `reuse=0` (`735e525c8`) and
    `wall=<level>` (`c68379d08`);
  - extension write interval (`4defd8148`) and restart monitor files
    (`2f20dbf8f`);
  - README: the convergence-test finding (`9f33c42a2`) and the μ*
    justification (`797042dbe`).
- **Build.** OpenFOAM v2412 on MeluXina (EasyBuild foss-2024a, GCC 13.3),
  solids4foam `development` at `aee0c8e35`. Every compile line uses
  `-O2 -fno-tree-vectorize`.
- **Runs.** The strict runs are in `/project/home/p201136/philip-hv/hv_S`
  and the restart tests in `.../philip-hv/RT`. The default-test runs are
  archived (tars) in `.../philip-hv/archive-p201287`.

## 3. The mathematical model

- **Geometry.** The published fluid-domain surface plus a 250 mm outlet
  extension. The flap is 11 x 2 x 65 mm, clamped at z = 0, with the flow
  along +z.
- **Fluid.**
  - Incompressible Newtonian, ρ = 1163.3 kg/m³, ν = 1.0745e-5 m²/s,
    laminar.
  - Inlets: parabolic, with peaks of 630 (upper) and 615 (lower) mm/s,
    ramped as 3s² − 2s³ over 0.5 s. Outlet: p = 0 with `inletOutlet`
    velocity. No-slip walls; the flap is a moving wall.
  - Numerics: PIMPLE (1 outer and 3 pressure correctors), BDF2,
    `linearUpwind` convection with a cell-limited gradient, limited
    corrected Laplacians, ALE `velocityLaplacian` mesh motion.
- **Solid.**
  - Total-Lagrangian neo-Hookean, μ* = 64.2 kPa, ν = 0.45,
    ρ_s = 1058.3 kg/m³. The net buoyancy body force,
    g_eff = (ρ_s − ρ_f)/ρ_s g = 0.973 m/s² in +y, is ramped with the
    inflow.
  - Numerics: cell-centred FV, BDF2, PETSc SNES. Momentum stabilisation
    `diffStencilLaplacian` with scale factor 0.01, a consistent O(h²) term
    that is part of the discretisation.
- **Coupling.**
  - Partitioned IQN-ILS over a non-conforming AMI interface with a
    predictor.
  - Strict test (`requireAllResidualMeasures yes`), `nOuterCorr 60`,
    `couplingReuse 0`, `outerCorrTolerance 1e-4` (F3S3: 1e-3, justified in
    Section 7).
  - `allowUnconvergedCoupling yes`, and every step accepted above the
    tolerance is counted (Section 9).
- **End state.** Steady. The QoI is the tip displacement averaged over the
  last 4 s of a settled run.
- **Parameters.** Only the meshes, the time step and the iterative settings
  vary between runs.

## 4. Fixed shear modulus

The silicone is calibrated against the measured zero-flow tip deflection.
The flap alone under net buoyancy is damped to rest on the nested meshes
S1-S4 (ratio 2):

| Level | Cells | μ_cal, sf 0.01 (kPa) | μ_cal, sf 0.001 (kPa) |
| --- | --- | --- | --- |
| S1 (6x4x33) | 792 | 44.81 | 50.07 |
| S2 (12x8x65) | 6 240 | 58.42 | 60.68 |
| S3 (24x16x130) | 49 920 | 62.58 | 63.23 |
| S4 (48x32x260) | 399 360 | 63.61 | 63.79 |
| Richardson (S2-S4) | | 63.95 (p = 2.01, GCI 0.42) | 63.94 (p = 2.21, GCI 0.19) |

- **Choice.** μ* = 64.2 kPa was fixed from S1-S3 before S4 existed. It lies
  within the S2-S4 GCI of the 63.95 kPa limit.
- **Effect on the coupled tip.** About −0.007 mm per kPa, so the 0.25 kPa
  offset is worth about +0.002 mm. Any continuum statement refers to the
  model with μ* = 64.2 kPa.
- **The tutorial value.** 58.45 kPa is the S2 calibration. It compensates
  the S2 discretisation error.
- **Original authors.** With ν = 0.45, about 64 kPa is consistent with the
  ≈ 61 kPa of Hessenthaler et al. for a nearly incompressible model.
- **Calibration uncertainty.** The calibration uncertainty from the measured
  deflection is a model (validation) uncertainty.

## 5. Mesh families

- **Solid.** S1-S3 as above, nested, with S2 → S3 refining by exactly 2.
- **Fluid.** cfMesh `cartesianMesh`. Level Fk scales every cell size of the
  coarse dictionary by 2^(−(k−1)/2): a √2 refinement per level, not nested.

| Level | Cells | Interface faces | Interface / wall cell size (mm) |
| --- | --- | --- | --- |
| F1 | 231 034 | 7 184 | 0.505 / 2.02 |
| F2 | 605 497 | 14 252 | 0.357 / 1.43 |
| F3 | 1 531 752 | 28 062 | 0.253 / 1.01 |
| F1 + F3 walls | 459 198 | 7 184 | 0.505 / 1.01 |

## 6. Steady state and end time

- **Approach.** In the strict runs the tip approaches its steady value
  geometrically, with τ ≈ 4 s on F1S2. The 2-s means from 26 to 45 s are
  16.6495, 16.6516, 16.6529, 16.6538, 16.6544, 16.6547, 16.6549, 16.6550,
  16.6550 and 16.6551 mm.
- **Control restart.** An unchanged restart from 45 to 55 s stays at
  16.6551 mm, with a drift of −5e-7 mm/s.
- **End-of-run drift.** Every steady value in Section 8 is a 4-s mean with a
  final drift of at most about 1e-4 mm/s.
- **Settling time.** The coarse levels need 35-40 s and the finer ones up to
  40-45 s. The tutorial's 10 s end time is far from steady: at 10 s the F1S2
  tip is still 0.25 mm below its steady value.

## 7. Iterative, Δt and stabilisation sensitivity (F1S2 steady state)

Each restart from the converged state at 45 s changes one setting. The tip
is compared over 52-55 s with the unchanged control:

| Change | Steady tip y (mm) | Shift (mm) |
| --- | --- | --- |
| control | 16.6551 | 0 |
| Δt 0.002 → 0.001 s | 16.6515 | −0.0036 |
| Δt 0.002 → 0.004 s | 16.6579 | +0.0028 |
| tolerance 1e-4 → 1e-5, tight fluid solvers | 16.6551 | < 1e-4 |
| tolerance 1e-4 → 1e-3 | 16.6551 | < 1e-4 |
| stabilisation 0.01 → 0.001 | 16.7048 | +0.0497 |

- **Δt.** At a true steady state the BDF2 time derivatives vanish. The small
  monotone Δt dependence of the discrete steady state comes from the
  Δt-dependent face-flux (momentum) interpolation of the PIMPLE
  pressure-velocity coupling. Part of it may be residual settling: the
  Δt-variant restarts were still moving by about ±1e-4 mm/s. It is
  reported as a Δt sensitivity of 0.004 mm, not as a truncation-error
  estimate.
- **Coupling tolerance.**
  - At the end state: < 1e-4 mm.
  - Along the path: a strict run at tolerance 1e-3 from t = 0 is 0.015 mm
    ahead at 20 s but reaches 16.6554 mm (40-45 s), against 16.6550 mm at
    1e-4. The coupling (iterative) contribution is therefore ≤ 0.0004 mm.
    This justifies tolerance 1e-3 for the F3S3 run.
- **Stabilisation.** +0.050 mm on S2. The term is a consistent part of the
  spatial discretisation, so its effect is contained in the solid
  refinement (Section 8). The sensitivity is a robustness check, not an
  additional uncertainty.

## 8. Spatial convergence

Steady tip y (mm), strict test, Δt = 0.002 s:

| | S1 | S2 | S3 |
| --- | --- | --- | --- |
| F1 | 16.2006 | 16.6550 | 16.7530 |
| F2 | | 16.6920 | |
| F3 | | 16.7219 | 16.8133 (tolerance 1e-3) |

**Solid refinement (F1; S1 → S2 → S3, ratio 2).**
- Differences: +0.4544 and +0.0979 mm. The convergence ratio is 0.215, so
  p = 2.21.
- Richardson: 16.7799 mm. GCI on S3 (Fs = 1.25): 0.034 mm (0.20 %).
- Statement: at fixed F1 and material parameters, the tip shows monotone
  solid-mesh convergence consistent with approximately second-order
  accuracy, with a conditional S3 GCI estimate of 0.034 mm.

**Fluid refinement (S2; F1 → F2 → F3, ratio √2).**
- Differences: +0.0370 and +0.0299 mm. The convergence ratio is 0.81, so
  p_obs = 0.62 against a formal order of 2.
- Richardson(p_obs): 16.847 mm.
- Provisional GCI indicator (Fs = 1.25): 0.16 mm.
- Asymptotic convergence has not been established. Alternative rules give
  smaller values that rest on unsupported assumptions: 3|Δ32| = 0.090 mm,
  or 0.072-0.090 mm with p = 1 assumed. The p_obs-based 0.16 mm is
  retained as the indicator.
- Supporting evidence (it shows where the fluid error is not, but does not
  establish tip convergence):
  - Upstream quantities converge at about second order. Lower-jet peak at
    z = 10 mm: 572.1, 576.9, 579.4 mm/s (p ≈ 2.0). Upper inlet pressure:
    3.50, 2.59, 2.29 Pa (p ≈ 3.1). The lower inlet pressure is the slowest:
    3.33, 2.54, 1.96 Pa (p ≈ 0.9).
  - Refining only the F1 walls to F3 size shifts the tip by +0.013 mm,
    about 20 % of F1 → F3 (0.067 mm). Wall resolution alone does not drive
    the fluid error. This does not allocate the remaining error uniquely.

**Additivity.**
- F3S3 = 16.8133 mm, against F3S2 + (F1S3 − F1S2) = 16.8198 mm. The
  interaction is −0.0065 mm.
- The solid and fluid corrections are approximately additive over the
  tested rectangle. This does not bound the interaction in the continuum
  limit.

**Matched sequence (F1S1, F2S2, F3S3).** 16.2006, 16.6920 and 16.8133 mm.
The refinement ratios are mixed (2 for the solid, √2 for the fluid), so no
single Richardson order is derived from it.

## 9. Finding: the default FSI convergence test accumulates an iterative error

- **The test.** A step is accepted when `min(r1, r2)` is below
  `outerCorrTolerance`. r1 is the interface residual relative to its
  largest value in the step. r2 is normalised by the TOTAL interface
  displacement.
- **Per-step histories.** From the default-test F2S2 run, 40,000 steps over
  0-80 s. Values are taken at the accepting iteration:

  | Window (s) | Iterations per step | r1 (median) | Steps with r1 > 1e-2 | r2 (median) |
  | --- | --- | --- | --- | --- |
  | 0-1 | 4.65 | 9.3e-4 | 4 % | 5.4e-5 |
  | 1-5 | 3.19 | 4.9e-2 | 87 % | 3.9e-5 |
  | 5-80 (each window) | 3.10-3.12 | 6.8e-2 | 90-91 % | 3.5e-5 |

- **Consequences in the default-test runs:**
  - The tip rose 0.3-0.4 mm between 8 and 30 s and was still rising at
    80 s (16.84/16.90/16.98 mm at 30/60/80 s on F2S2).
  - Rank replicates differed by 0.03-0.06 mm, and Δt = 0.001 sat
    0.17-0.3 mm above Δt = 0.002.
  - The fluid interface lagged the deformed solid by about 0.03 mm in −y
    (mean gap 0.047-0.066 mm).
  - The early transient was distorted (tip at 0.5 s: 3.66 mm vs 1.91 mm
    strict).
- **Strict test.** 10-15 iterations per step, median r1 5e-5. On F1S2 the
  mean interface gap is 0.0045 mm with no systematic offset. The 32- and
  48-rank replicates agree to 0.001 mm.
- **Diagnosis.** Neither component drifts on its own:
  - with the geometry frozen, the fluid settles in about 2 s;
  - the solid alone is constant to 1e-15 mm.
- **Restarts keep history.** The residual compares displacement increments
  (`incrementalResiduals` true), so a mismatch built up under the default
  test is never corrected afterwards. Strict restarts from default-test
  states stayed at 17.0195 and 16.90 mm. Only strict runs from t = 0
  define one trajectory and one end state.
- **Robustness and cost of the strict test.**
  - With `couplingReuse 5`, stale secant modes made the first update of a
    step diverge three times (interface updates of 32, 72 and 992 m).
    `couplingReuse 0` is robust.
  - Near steady state, a step that starts at round-off cannot reduce r1 by
    a further 1e-4, because of the noise floor. Such steps run to 60
    iterations and are accepted with r1 ≈ 1.2-3e-4 and r2 ≤ 6e-11 (F1S2:
    5 % of steps over 30-45 s; F1S3: 25 %). At tolerance 1e-3 this does
    not occur (F3S3: 0 % over 30-40 s). The restart tests show that such
    steps do not affect the steady value.
  - Recommendation: use the strict test with a tolerance above the noise
    floor (1e-3 here).
- **Scope.** This is solids4foam issue #552. Any steady or slowly varying
  FSI case with a large total interface displacement is exposed.

## 10. Default-test results (superseded, kept as the record)

The default-test campaign (F1-F4 x S1-S3 at 15-80 s, with Δt, tolerance
and stabilisation variants) is in
`verification/postProcessing/verification_runs.csv` on the solids4foam
branch. Its spatial series is non-monotone inside the iterative drift. None
of its values is used for verification.

## 11. High-order solid

Not used. Earlier high-order zero-flow calibrations (cubic, 62-64 kPa) agree
with the standard-solid continuum limit of Section 4.

## 12. Numerical uncertainty budget (best-resolved tip, F3S3)

Deterministic estimates and sensitivities, combined as a linear aggregate
allowance. There is no basis for treating them as independent random
errors, so no RSS is used.

| Contribution | Value (mm) | Basis |
| --- | --- | --- |
| Fluid discretisation | 0.16 | provisional Richardson-based GCI indicator, p_obs = 0.62, Fs = 1.25; asymptotic range not established |
| Solid discretisation | 0.034 | GCI on S3 (F1), p = 2.21; transfer to F3 supported by the −0.0065 mm interaction |
| Δt | 0.004 | steady-state Δt sensitivity (Section 7) |
| Coupling / iterative | 0.0004 | tolerance 1e-3 vs 1e-4 from t = 0 |
| End state | ≤ 0.001 | final drift < 1e-4 mm/s; control restart |
| Aggregate allowance | ≈ 0.20 (1.2 %) | conditional and provisional |

- **Bias correction.** The PR #546 transfer-weight defect is a known bias,
  not an uncertainty. It is corrected: +0.006 mm on S3 levels, so F3S3
  becomes 16.819 mm (Section 15).
- **Extrapolation scenario (not an established interval).** F3S3 + the
  solid correction (+0.027) + the fluid correction (+0.030 to +0.125)
  gives a continuum tip of about 16.88-16.97 mm for μ* = 64.2 kPa.
- **Wording.** The tip is numerically resolved to about 1.2 % (conditional),
  with the fluid discretisation dominant. Δt and iterative contributions
  are below 0.03 % and 0.003 %.

## 13. Validation against MRI

Driver `validation_metrics`, applied to the time-averaged centreline and the
voxel velocities of each strict run. The MRI tip is 16.41 mm.

| Strict run | Tip y (mm) | Tip − MRI (mm) | Centreline RMS / max (mm) | Velocity d_mean / d_max (/630 mm/s) |
| --- | --- | --- | --- | --- |
| F1S1 | 16.199 | −0.21 | 0.27 / 0.59 | 0.051 / 0.348 |
| F1S2 | 16.654 | +0.24 | 0.29 / 0.62 | 0.050 / 0.351 |
| F2S2 | 16.691 | +0.28 | 0.30 / 0.64 | 0.050 / 0.332 |
| F1S3 | 16.753 | +0.34 | 0.37 / 0.75 | 0.050 / 0.351 |
| F3S3 | 16.813 | +0.40 | 0.40 / 0.78 | 0.049 / 0.332 |

- **Tip and centreline.** Every refinement moves the tip, and the
  centreline, further from the measurement. With the extrapolation
  scenario of Section 12, the discrepancy of the continuum model would be
  about 0.47-0.56 mm (2.8-3.4 %).
- **Velocities.** The velocity metrics barely change with refinement.
- **What this shows.** Agreement with the MRI does not indicate numerical
  convergence, and the coarse-mesh agreement is partly a cancellation of
  discretisation error.
- **What would be needed to judge the discrepancy.** An experimental
  uncertainty (segmentation, registration, geometry, loading and material
  calibration), which is not available here. The voxel size (1.302 x 1.302
  x 6 mm) is not a measurement uncertainty.

## 14. Independent numerical comparison

No COMSOL solution is available.

## 15. PR #546 exposure

- **Defect (a): -O3 aliasing.** Not applicable: the build is -O2.
- **Defect (b): barycentric weights.** It applies only to solid interface
  fan triangles with an area below 0.0158 mm²: the S3 side and tip faces,
  and S4. S1 and S2 are unaffected.
- **A/B, old vs fixed build, 750 steps from the same state:**
  - F1S2: bit-identical.
  - F1S3: +0.006 mm (mean tip over the last second).
- **Correction.** This bias is applied as a correction to the S3 levels
  (Section 12).

## 16. Limitations

- **Fluid convergence.** The fluid meshes are not nested, and the tip's
  fluid convergence is not demonstrated to be asymptotic. A fourth fluid
  level (4.2 M cells, about 800 node-hours with the strict test) was not
  affordable.
- **Additivity.** It is checked on the F1-F3 x S2-S3 rectangle only.
- **Cost.** Hessenthaler Phase I needs runs far longer than the tutorial's
  10 s, and the strict convergence test (or an equivalent absolute
  interface criterion). It costs 3-5x the default test.
- **Validation.** It lacks an experimental uncertainty.

## 17. Recommended manuscript claims

May claim:
- Phase I is a solution-verification case with a fixed model (μ* = 64.2 kPa
  from a mesh-converged calibration). It is validated against MRI
  separately.
- The default convergence test of the partitioned solver produced an
  accumulating iterative error that mimicked a slow physical creep. A
  strict test removes it. This is an example of a converged-looking QoI
  that is not converged.
- At fixed F1, solid refinement converges monotonically at about second
  order (S3 GCI 0.034 mm). Under fluid refinement the tip converges
  monotonically but slowly (p_obs ≈ 0.6). Its provisional GCI indicator
  (0.16 mm) dominates the budget.
- Δt and coupling contributions are negligible (≤ 0.004 and ≤ 0.0004 mm).
- The best-resolved tip is 16.81 mm (16.82 mm with the PR #546
  correction), numerically resolved to about 1.2 % (conditional).
- Refinement moves the solution away from the MRI tip. The best-resolved
  model lies about 0.4 mm (2.5 %) above it.

Must not claim:
- A demonstrated total numerical accuracy of 1.2 %, or a statistically
  justified (RSS) uncertainty.
- Asymptotic fluid convergence, or an established continuum interval. The
  continuum extrapolation (16.88-16.97 mm) is a scenario only.
- Validation on the grounds that the discrepancy is "within one voxel".
- A single Richardson order from the mixed-ratio matched sequence.
- That the coarse-mesh agreement with the MRI (16.5-16.7 mm in the
  tutorial) reflects numerical accuracy.
