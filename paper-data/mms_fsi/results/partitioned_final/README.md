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
| `aitken_full_*`, `robin_full_*` | full variant, n = 16, 32, 64, 128 (Table `tab:mms_fsi`) |
| `aitken_flat_*`, `robin_flat_*` | flat variant, n = 16, 32, 64 |
| `aitken_full_noSSG_*` | A/B without the surface-gradient traction area, n = 64, 128 |
| `before_aitken_full_*` | the earlier case settings on this build, n = 16, 32, 64 |
| `*_errors.csv` | errors at t = 0.1, 0.2, 0.3 and 0.4 for each mesh (columns as in `../../README.md`) |
| `*_orders.md` | order tables at t = 0.2 and 0.4, from `case/orders.py` |
| `coupling_iterations.csv` | coupling iterations per step (mean, max), from `postProcessing/fsiResiduals.dat`, as in `../../make_csv.py` |
| `makeCase.patch` | the changes to the case generator |

## Remaining limitations

- Robin–Neumann pressure: L2 orders 1.61, 1.77, 1.66; the interface pressure
  condition imposes only the acceleration term of the normal momentum balance.
  Predicted interface-row error (h/2)·4.3 = 0.067 (n64), 0.034 (n128); measured
  0.061, 0.031. Fix in progress.
- Aitken interface-face pressure (`intP_L2`): orders 1.69, 1.60, 1.28 at t = 0.4;
  cause not identified.
