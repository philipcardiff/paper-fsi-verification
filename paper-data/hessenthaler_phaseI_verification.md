# Hessenthaler Phase I verification findings

Evidence report for the Hessenthaler et al. (2017) Phase I benchmark
(`sec:hessenthaler` in `sectionBenchmarkSuite.tex`) before any manuscript
change. It addresses item 7 of `manuscript_update_2026-10.md` (strict-test
verification runs from t = 0, the PR #546 A/B check, the justification of
μ*). No `.tex` file is touched.

Compact data: `hessenthaler_phaseI_verification_steady.csv` (steady tip
displacement of every strict run and restart test) and
`hessenthaler_phaseI_verification.json` (the same, with tail windows and
drift rates). They are generated from the tip histories of the runs listed
in Section 2.

**STATUS: DRAFT.** Sections 5-8 (spatial convergence) are incomplete. The
F3S2, F1S3 and the F2S2 extension to 45 s are still running. Every other
section is final.

## 1. Executive scientific conclusion (provisional)

- **Verification vs validation.** Phase I is used here as a numerical
  verification case: one fixed mathematical model, solved on refined meshes
  and time steps. The MRI data are used only for validation, after the
  numerical picture is established (Section 13).
- **Fixed model.** The model has one fixed shear modulus, μ* = 64.2 kPa
  (ν = 0.45), taken from the continuum limit of a mesh-converged zero-flow
  calibration (Section 4). It is not recalibrated per mesh.
- **The default FSI convergence test produced an iterative error that
  looked like physics** (Section 9). This is the main finding. With
  `outerCorrTolerance 1e-4`, the solids4foam default accepts a step when
  the smaller of two residual measures is below the tolerance. One measure
  is normalised by the total interface displacement, which is about 17 mm
  here. After the first second, every step was accepted with the in-step
  residual r1 above the tolerance (median 6.8e-2). The interface lag
  accumulated into a tip "creep" of +0.003 to +0.005 mm/s that persisted to
  80 s. That creep exceeded every spatial difference, and it distorted the
  whole transient (tip at t = 0.5 s: 3.66 mm vs 1.91 mm).
- **With both measures required (strict test), the solution is
  well-defined.** On F1S2 the tip settles with a time constant of about 4 s
  to 16.6550 mm, flat to 3e-5 mm/s over 40-45 s. Runs on two rank counts
  agree to 0.001 mm. The fluid interface sits on the deformed solid to
  0.0045 mm on average (default test: 0.066 mm, with a systematic 0.03 mm
  lag).
- **Temporal and iterative errors are negligible.** At the F1S2 steady
  state, halving or doubling the time step shifts the tip by -0.0036 and
  +0.0028 mm (0.02 %). Coupling tolerances of 1e-3 and 1e-5 change it by
  less than 1e-4 mm.
- **Spatial convergence (pending; Section 8).** The solid refinement
  dominates: F1 with S1 → S2 gives +0.46 mm. The fluid effect is smaller:
  F1 → F2 on S2 gives about +0.04 mm. The stabilisation term contributes
  0.05 mm on S2 and vanishes as O(h²).

## 2. Provenance

- **solids4foam.** Branch `verification/hessenthaler-phaseI`, based on the
  PR #535 head `4f549b500`. `src/` is unchanged. The verification driver
  (`tutorials/fluidSolidInteraction/hessenthalerFsi/verification/`) adds
  the spec keys `fsitest=strict` (`a2ccaf64c`) and `reuse=0` (`735e525c8`).
  The README documents the convergence-test finding (`9f33c42a2`) and the
  μ* justification (`797042dbe`).
- **Build.** OpenFOAM v2412, MeluXina (EasyBuild foss-2024a, GCC 13.3),
  solids4foam `development` at `aee0c8e35`. Compiled at
  `-O2 -fno-tree-vectorize` (every compile line in `log.Allwmake`), so the
  build is not affected by the -O3 aliasing defect of PR #546. The PR #546
  A/B is reported in Section 15.
- **Runs.** MeluXina, accounts p201287, p201136 and p201141. The strict runs
  are under `/project/home/p201136/philip-hv/hv_S`, the restart tests under
  `.../philip-hv/RT`, and the default-test runs under
  `/project/home/p201287/philip/s4f-verif-v2412/runs/hv_*`.

## 3. The mathematical model

- **Geometry.** The published fluid-domain surface with a 250 mm outlet
  extension. The flap is 11 x 2 x 65 mm, clamped at z = 0, with the flow
  along +z.
- **Fluid.** Incompressible Newtonian, ρ = 1163.3 kg/m³,
  ν = 1.0745e-5 m²/s, laminar. Inlets are parabolic, with peaks of 630
  (upper) and 615 (lower) mm/s, ramped as 3s² − 2s³ over 0.5 s. The outlet
  has p = 0 and `inletOutlet` velocity. Walls are no-slip; the flap is a
  moving wall.
  - Discretisation: PIMPLE (1 outer and 3 pressure correctors), BDF2,
    `linearUpwind` convection with a cell-limited gradient, limited
    corrected Laplacians, ALE `velocityLaplacian` mesh motion (quadratic
    inverse-distance diffusivity).
- **Solid.** Total-Lagrangian neo-Hookean, μ* = 64.2 kPa, ν = 0.45,
  ρ_s = 1058.3 kg/m³. Buoyancy enters as a net body force,
  g_eff = (ρ_s − ρ_f)/ρ_s g = 0.973 m/s² in +y, ramped with the inflow.
  - Discretisation: cell-centred FV, BDF2 d2dt2, PETSc SNES, momentum
    stabilisation `diffStencilLaplacian` with scale factor 0.01. The
    stabilisation is a consistent O(h²) term, so it is part of the
    discretisation error.
- **Coupling.** Partitioned IQN-ILS over a non-conforming AMI interface,
  with a predictor and `outerCorrTolerance 1e-4`.
  - Verification runs: strict test (`requireAllResidualMeasures yes`,
    `nOuterCorr 60`, `allowUnconvergedCoupling yes`) and `couplingReuse 0`.
  - Steps still above the tolerance after 60 iterations are accepted and
    counted. There are two kinds:
    - In the start-up transient (t < 10 s): 0.1-0.5 % of steps, with r1
      up to ~7e-3.
    - Near the steady state (F1S2, 30-45 s: 384 of 7,500 steps, 5 %):
      r1 stops at about 1.25e-4, just above 1e-4, while the absolute
      measure r2 is at most 3e-12. A step that starts at round-off cannot
      reduce its residual by a further 1e-4 (a noise floor). The restart
      tests of Section 7 show that such steps do not affect the steady
      value.
- **End state.** Steady. The QoIs are the tip displacement averaged over
  the last seconds of a run that has settled (drift < 1e-4 mm/s).
- **Fixed parameters.** Only the meshes, the time step and the iterative
  settings change between runs. No physical or model parameter changes.

## 4. Fixed shear modulus

The silicone is calibrated against the measured zero-flow tip deflection.
The flap alone under net buoyancy, damped to rest, is solved on the nested
meshes S1-S4 (ratio 2 in every direction):

| Level | Cells | μ_cal, sf 0.01 (kPa) | μ_cal, sf 0.001 (kPa) |
| --- | --- | --- | --- |
| S1 (6x4x33) | 792 | 44.81 | 50.07 |
| S2 (12x8x65) | 6 240 | 58.42 | 60.68 |
| S3 (24x16x130) | 49 920 | 62.58 | 63.23 |
| S4 (48x32x260) | 399 360 | 63.61 | 63.79 |
| Richardson (S2-S4) | | 63.95 (p = 2.01, GCI 0.42) | 63.94 (p = 2.21, GCI 0.19) |

- **Choice of μ\*.** Both stabilisation levels converge at second order to
  63.9-64.0 kPa. μ* = 64.2 kPa was fixed from S1-S3 before S4 existed, and
  it lies within the S2-S4 GCI.
- **Effect on the coupled tip.** About −0.007 mm per kPa (F1S2, legacy
  runs at 58.45 and 64.2 kPa). The 0.25 kPa offset is therefore worth about
  0.002 mm.
- **The tutorial value.** 58.45 kPa is exactly the S2 calibration
  (58.42 kPa). It compensates the S2 discretisation error and is not a
  model parameter.
- **Comparison with the original authors.** With ν = 0.45 the bending
  stiffness needs a higher μ than with ν ≈ 0.5 (by 3-10 % depending on the
  beam/plate limit). 64 kPa is therefore consistent with the ≈ 61 kPa
  calibration of Hessenthaler et al.
- **What the calibration uncertainty is.** The measurement uncertainty of
  the calibration is a validation (model) uncertainty, not a numerical one.

## 5. Mesh families

- **Solid.** S1-S4 above. They are nested, and S2 → S3 → S4 refines by
  exactly 2.
- **Fluid.** cfMesh `cartesianMesh`. Level Fk scales every cell size of the
  tutorial's coarse dictionary by 2^(−(k−1)/2), so each level refines by a
  nominal √2. The meshes are not nested.

| Level | Cells | Interface faces | Interface cell size (mm) |
| --- | --- | --- | --- |
| F1 | 231 034 | 7 184 | 0.505 |
| F2 | 605 497 | 14 252 | 0.357 |
| F3 | 1 531 752 | 28 062 | 0.253 |
| F4 | 4 177 457 | 56 704 | 0.179 |

## 6. Steady state and end time

- **F1S2 (strict test).** The 2-s means from 26 to 45 s are 16.6495,
  16.6516, 16.6529, 16.6538, 16.6544, 16.6547, 16.6549, 16.6550, 16.6550
  and 16.6551 mm. The increments shrink geometrically (τ ≈ 4 s), and the
  drift over 40-45 s is 3e-5 mm/s. A restart with unchanged settings
  (45-55 s) stays at 16.6551 mm with a drift of −5e-7 mm/s.
- **End time.** The end-state error after 30 s is below 0.003 mm on F1S2.
  F2S2 settles more slowly: still +0.001 mm/s at 32 s, with a 45 s
  extension running.
- **The tutorial's 10 s end time is not sufficient** for any level: at 10 s
  the F1S2 tip is still 0.25 mm below its steady value.

## 7. Iterative, temporal and stabilisation sensitivity (F1S2 steady state)

Each restart starts from the converged state at 45 s and changes one
setting. The steady tip is compared over 52-55 s with the unchanged
control:

| Change | Steady tip y (mm) | Shift (mm) |
| --- | --- | --- |
| control | 16.6551 | 0 |
| Δt 0.002 → 0.001 s | 16.6515 | −0.0036 |
| Δt 0.002 → 0.004 s | 16.6579 | +0.0028 |
| tolerance 1e-4 → 1e-5, tight fluid solvers | 16.6551 | < 1e-4 |
| tolerance 1e-4 → 1e-3 | 16.6551 | < 1e-4 |
| stabilisation 0.01 → 0.001 | 16.7048 | +0.0497 |

- **Time step.** The steady value depends weakly on Δt, through the
  Δt-dependent face-flux interpolation of the pressure-velocity coupling.
  The shift is monotone, at 0.003-0.004 mm per factor 2 (0.02 %). The Δt
  = 0.001 and 0.004 restarts were still converging slowly (±1e-4 mm/s).
- **Coupling tolerance.** Its end-state effect is below 1e-4 mm. Path
  independence is checked separately: the from-zero tolerance 1e-3 run is
  0.015 mm ahead at 20 s, and is being extended to 45 s.
- **Stabilisation.** +0.05 mm (0.3 %) on S2. It is part of the solid
  discretisation error and is expected to fall by about 4x per solid level.

## 8. Spatial convergence (PENDING)

Steady tip y (mm), strict test, Δt = 0.002 s:

| | S1 | S2 | S3 |
| --- | --- | --- | --- |
| F1 | 16.196 (30 s, drift 1e-3 mm/s) | 16.6550 (45 s) | running |
| F2 | | 16.688 (32 s, drift 1e-3 mm/s; 45 s running) | |
| F3 | | running | |

The orders, Richardson estimates and uncertainties will be added when
F1S3, F3S2 and the F2S2 extension complete.

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

  Strict test: 10-11 iterations per step, median r1 5e-5 and r2 < 1e-8.
- **Consequences in the default-test runs:**
  - The tip rose by 0.3-0.4 mm between 8 and 30 s and was still rising at
    80 s (16.84/16.90/16.98 mm at 30/60/80 s on F2S2).
  - Rank replicates ended up 0.03-0.06 mm apart, and Δt = 0.001 sat
    0.17-0.3 mm above Δt = 0.002.
  - The fluid interface lagged the solid by about 0.03 mm in −y (mean gap
    0.047-0.066 mm).
- **Diagnosis.** Neither component drifts on its own:
  - with the geometry frozen, the fluid settles in about 2 s;
  - the solid alone under buoyancy is constant to 1e-15 mm.
- **Restarts keep history.** The residual compares displacement increments
  (`incrementalResiduals` true), so a position mismatch accumulated under
  the default test is never corrected afterwards. Strict restarts from two
  default-test states stay at 17.0195 and 16.90 mm. Only runs that are
  strict from t = 0 define one trajectory and one end state.
- **Robustness of the strict test.** With `couplingReuse 5`, stale secant
  modes made the first IQN-ILS update of a step diverge three times (32 m,
  72 m, 992 m), each after accepted unconverged steps. `couplingReuse 0` is
  robust.
- **Scope.** This is solids4foam issue #552. Any steady or slowly varying
  FSI case with a large total interface displacement is exposed.

## 10. Default-test results (superseded, kept as the record)

The default-test campaign (F1-F4 x S1-S3 at 15, 30, 45, 60 and 80 s; Δt;
tolerances; stabilisation) is in
`verification/postProcessing/verification_runs.csv` on the solids4foam
branch. Its spatial series is non-monotone inside the iterative drift.
None of its values is used for verification.

## 11. High-order solid

Not used. The earlier high-order zero-flow calibrations (cubic, 62-64 kPa)
agree with the standard-solid continuum limit of Section 4 and support the
same μ*.

## 12. Numerical uncertainty budget (PENDING)

## 13. Validation against MRI (PRELIMINARY)

These metrics are computed with the driver's `validation_metrics` from the
time-averaged centreline and voxel velocities of each run (MRI tip
16.41 mm). They are not used anywhere above.

| Strict run | Tip y (mm) | Tip − MRI (mm) | Centreline RMS / max (mm) | Velocity d_mean / d_max (/630 mm/s) |
| --- | --- | --- | --- | --- |
| F1S1 (30 s) | 16.189 | −0.22 | 0.278 / 0.597 | 0.051 / 0.347 |
| F1S2 (45 s) | 16.654 | +0.24 | 0.291 / 0.618 | 0.050 / 0.351 |
| F2S2 (20 s, not settled) | 16.610 | +0.20 | 0.261 / 0.591 | 0.049 / 0.326 |

- Solid refinement moves the tip from below the measurement (S1) to above
  it (S2). The finer levels move it further above.
- The converged model will lie about 0.3-0.4 mm (about 2 %) above the MRI
  tip. That is well inside one MRI voxel (1.302 x 1.302 x 6 mm).
- The velocity metrics barely change with refinement.
- Final numbers will be taken on the best-resolved run.

## 14. Independent numerical comparison

No COMSOL solution is available.

## 15. PR #546 exposure

- **Defect (a): -O3 aliasing.** Not applicable: the build is -O2.
- **Defect (b): barycentric weights.** It applies only where solid
  interface fan triangles have an area below 0.0158 mm². That excludes
  S1/S2 entirely. It affects the S3 side and tip faces, and S4.
- **A/B, old vs fixed build, 750 steps from the same state:**
  - F1S2: bit-identical tip and pressure/viscous forces at every step.
  - F1S3: the mean tip over the last second changes by +0.006 mm (0.04 %),
    15x smaller than S2 → S3.
- **For the remaining S3 runs.** The strict S3 runs use the pre-#546 build.
  This exposure (≤ 0.01 mm) is included in the S3 uncertainty.

## 16. Limitations

- The fluid meshes are not nested. The interface sizes are tied to the
  level.
- Hessenthaler Phase I needs runs much longer than the tutorial's 10 s, and
  it needs the strict convergence test.
- The strict test costs about 3.5x the default one.

## 17. Recommended manuscript claims (PROVISIONAL)

To be finalised with Sections 8, 12 and 13.
