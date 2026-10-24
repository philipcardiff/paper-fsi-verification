# newtonIcoFluid continuity fix: monolithic results need rerunning

Note for the agents working on the quasi-monolithic verification cases
(2026-10-19). No `.tex` file is touched.

## Summary

The weaker monolithic FSI3 oscillation found in the Hron–Turek comparison
(monolithic `u_y` amplitude about 10% below IQN-ILS on the same 2x mesh) was a
defect in the fluid model `newtonIcoFluid`, not in the coupling. It is fixed in
solids4foam draft PR
[#597](https://github.com/solids4foam/solids4foam/pull/597) (branch
`fix/newtonicofluid-continuity-flux`, based on
`feature-petsc-snes-quasi-monolithic`).

**Every monolithic result recorded before this fix is affected and should be
rerun with it**:

- Hron–Turek (FSI1, FSI3);
- womersleyTube;
- 3dTube;
- Mok cavity;
- beamInCrossFlow;
- any other case using `newtonIcoFluid` with a non-`linear` `interpolate(U)`
  scheme.

All monolithic tutorial variants use `interpolate(U) skewCorrected linear`, so
all are affected.

## Cause

`newtonIcoFluid` convects momentum with `phi = interpolate(U) & Sf`, but its
continuity residual was `-fvc::div(U)`, which interpolates `U` with the
separate `div(U)` scheme (`Gauss linear`). With `skewCorrected linear` for
`interpolate(U)`, the flux that convects momentum is therefore not
divergence-free at convergence, and this dissipates transient flow. The fix
builds continuity from the same flux, `fvc::div(fvc::interpolate(U) & Sf)`.

## Evidence

**Rigid plate.** Fluid only, the Hron–Turek geometry with the plate held
rigid (CFD3-like, Re = 200, Δt = 1e-3). `newtonIcoFluid` was restarted from a
fully shedding `pimpleFluid` state, so late onset of shedding is not mistaken
for damping.

| Mesh | `pimpleFluid` | Before fix | After fix |
| --- | --- | --- | --- |
| 1x | grows to a limit cycle | decays ~5% per cycle | grows |
| 2x lift amplitude | 6.49 | 5.66 (−13%) | 6.27, rising (−3.4%) |

**FSI3, coupled.** 1x tutorial meshes for both couplings, with an AMI
interface for the monolithic coupling:

| Coupling | `u_y` amplitude (mm) | `u_x` mean (mm) | Frequency (Hz) |
| --- | --- | --- | --- |
| IQN-ILS | 29.83 | −2.23 | 5.587 |
| Monolithic, before | 24.22 (−18.8%) | −1.53 | 5.556 |
| Monolithic, after | 28.83 (−3.4%) | −2.14 | 5.587 |

**Ruled out while looking for the cause.** Each gave identical or worse
results:

- the gradient limiter in `div(phi,U)`;
- mesh motion and the coupling (fluid-only static mesh identical);
- the plate pressure condition;
- the Rhie–Chow scale;
- solver tolerances;
- PIMPLE splitting error.

**Unaffected cases:**

- `ringAddedMass` monolithic and non-conformal verification: wet/dry ratio
  errors still at most 0.004%;
- `decayingTaylorGreenVortex`, which uses `linear` interpolation.

**`blobInTreacle`.** Its monolithic mesh study was rerun; the extrapolated
steady displacement is now within 0.08% of the partitioned value (previously
0.12%).

## Practical notes

- **Onset of shedding.** It is noise-seeded, and a tightly converged Newton
  solve is quieter than PIMPLE, so comparisons over a fixed window from rest
  can confuse a late start with damping. Compare limit cycles, or restart from
  a common shedding state.
- **Comparing couplings.** The monolithic tutorial variant uses coarser
  conformal solid meshes than the partitioned variant, so its plate is
  stiffer. Use the same meshes (with `interfaceTransferMethod AMI`) when
  comparing couplings.
- **Full investigation.** On xenosim, see
  `/home/philipc/OpenFOAM/philipc-v2512/newtonIcoFluid-shedding-check/README.md`.
