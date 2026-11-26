# Final partitioned results of the manufactured FSI solution (26 Nov 2026)

These files are the source of Table `tab:mms_fsi` (Section `sec:mms_fsi`), the
partitioned columns of Table `tab:refvals_mms` (Appendix `app:reference_values`)
and the partitioned numbers in Appendix `app:lessons`. They supersede the
partitioned series in `../../summary_errors_t0.2_t0.4.csv`,
`../../coupling_iterations.csv` and `../../order_tables/`, which were computed
before the fixes and are kept for the record.

## Provenance

- **Build:** solids4foam `development` at 7a2f0a63e, merged locally with
  - solids4foam PR #756 (`tractionPatchSurfaceGradient`: surface-gradient
    evaluation of the deformed traction area of the total Lagrangian solid), and
  - solids4foam PR #776 (`pimpleFluid` `fluxConsistentPatches` also on
    fixed-pressure patches with a Neumann velocity: the outlet face flux is
    built from the boundary velocity),
  - plus a one-line merge fix in `nonLinGeomTotalLagVelocitySolid` (compact
    layout passed to the quadrature stress update; quasi-monolithic solid only).
- **Platform:** OpenFOAM v2512, PETSc 3.22.2, macOS (arm64), Apple Clang, `Opt`.
- **Environment:** `/Users/philipc/claude-runs/env-mmspart.sh` (isolated user
  build in `~/OpenFOAM/mmspart-v2512`).
- **Case generator:** `case/makeCase.py` of this directory tree with the changes
  in `makeCase.patch` (new options `--nOuter`, `--pimpleExtra` and
  `--exactOldTime`). Options used for every run:

  ```bash
  --outerTol 1e-6 --pimpleExtra "fluxConsistentPatches (outlet);" --nOuter 2 \
    --ddtCorr no --solidSurfGrad --exactOldTime
  ```

  with Δt = 0.32/n, run to t = 0.4.
- **Runs:** `/Users/philipc/claude-runs/mms-partitioned/runs/cp-*` (not copied).

## Case settings, in words

- Fluid (`pimpleFluid`): two outer correctors, four pressure correctors, one
  non-orthogonal corrector; time-derivative correction of the face flux
  (`ddtCorr`) off; outlet face flux from the boundary velocity.
- Solid (displacement-form total Lagrangian): deformed traction area from a
  surface gradient over the boundary face; old-time displacements `D_0`,
  `D_0_0`, `D_0_0_0` set to the exact solution at -Δt, -2Δt, -3Δt.
- Couplings as in `../../README.md` (Aitken: `newMovingWallVelocity` +
  `fixedFluxPressure`; Robin–Neumann: `elasticWallVelocity` +
  `elasticWallPressure`), tolerance 1e-6, legacy test, at most 200 iterations.

## Files

| File | Content |
|---|---|
| `aitken_full_*`, `robin_full_visc_*` | full variant, n = 16, 32, 64, 128 (Table `tab:mms_fsi`) |
| `robin_full_*` | full variant, Robin–Neumann with the acceleration-only Robin condition (superseded by `robin_full_visc_*`; quoted in the text for comparison) |
| `aitken_flat_*`, `robin_flat_*` | flat variant, n = 16, 32, 64 |
| `aitken_full_noSSG_*` | A/B without the surface-gradient traction area, n = 64, 128 |
| `before_aitken_full_*` | the earlier case settings on this build, n = 16, 32, 64 |
| `*_errors.csv` | errors at t = 0.1, 0.2, 0.3 and 0.4 for each mesh (columns as in `../../README.md`) |
| `*_orders.md` | order tables at t = 0.2 and 0.4, from `case/orders.py` |
| `coupling_iterations.csv` | coupling iterations per step (mean, max), from `postProcessing/fsiResiduals.dat`, as in `../../make_csv.py` |
| `makeCase.patch` | the changes to the case generator |

## Robin–Neumann with the viscous and source terms (`robin_full_visc_*`)

- **Build:** the tree above (development e38f3b06f, which contains #756 and
  #776) merged with the interface-pressure branch `fix/movingWallPressure-721`
  (merge beb712f4e: `elasticWallPressure` with the viscous term of the
  wall-normal momentum equation and the fvOptions momentum source,
  `viscousTerm yes`). OpenFOAM v2512, macOS (arm64), Apple Clang, `Opt`.
- **Case options:** as above (`robin`, `--pimpleExtra "fluxConsistentPatches (outlet);"
  --nOuter 2 --ddtCorr no --solidSurfGrad --exactOldTime`), Δt = 0.32/n, t = 0.4.
- **Runs:** `/Users/philipc/claude-runs/mms-robin-viscous/runs/cp-robin-visc-src/n{16,32,64}`
  and `cp-robin-visc-src128/n128` (not copied). Files written by
  `mms_robin_visc.py` (in this directory) from `postProcessing/mmsFsiErrors/errors.dat`
  (same columns as `robin_full_errors.csv`) and of `fsiResiduals.dat`.
- **Effect:** p L2 orders 1.81, 1.96, 1.76 (were 1.61, 1.77, 1.66); n128 p L2
  error 1.29e-3 (was 2.45e-3); p L∞ on n64 2.1e-2 (was 6.1e-2). The viscous term
  alone, without the source, made the errors worse.
- **Flat variant:** not rerun; `robin_flat_*` use the acceleration-only condition.

## Remaining limitations

- Robin–Neumann, acceleration-only condition (`robin_full_*`): pressure L2
  orders 1.61, 1.77, 1.66. Predicted interface-row error (h/2)·4.3 = 0.067 (n64),
  0.034 (n128); measured 0.061, 0.031. Removed by the full condition above.
- Aitken interface-face pressure (`intP_L2`): orders 1.69, 1.60, 1.28 at t = 0.4;
  cause not identified.
