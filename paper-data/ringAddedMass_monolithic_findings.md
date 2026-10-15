# ringAddedMass: quasi-monolithic coupling

Evidence for the quasi-monolithic Newton coupling on the ring added-mass
problem (Section `sec:ring`), the only case in which the experimental
quasi-monolithic solver is used.

## Provenance

- solids4foam branch `feature-ring-monolithic-verification` (commit
  `dcbe04505`), on `feature-petsc-snes-quasi-monolithic-update` (PR #580), which
  ports the quasi-monolithic solver of PR #233 onto `development`. The
  solver library was built from `feature-petsc-snes-components` (PR #581),
  whose serial results are bit-identical to PR #580.
- OpenFOAM v2512, macOS (Apple M-series), serial runs, 14 at a time; the
  execution times are therefore indicative only.
- Command: `tutorials/fluidSolidInteraction/ringAddedMass/verification/Allverify
  --study monolithic --cores 14`.
- `ringAddedMass_monolithic_runs.csv` is the driver's
  `ringAddedMass_results.csv` (21 runs: the dry ring and the IQN-ILS and
  monolithic couplings at three added-mass levels, 25, 50 and 100 steps per
  period, refinement-factor-2 mesh, standard solid).
- `ringAddedMass_monolithic_mesh.csv`: separate moderate-level runs on the
  refinement-factor-1, 2 and 4 meshes at 100 steps per period (the factor-2
  rows repeat the study runs).

## Set-up

The monolithic velocity solid has the standard second-order residual only, so
every run uses the standard solid and the comparison is through the wet/dry
frequency ratio, which removes the solid error (standard and high-order
ratios agree to 0.003% on the finest mesh, Section `sec:ring`). The
monolithic system converges in two Newton iterations per step (linear
problem).

## Findings

| Level | Steps | Ratio error, IQN-ILS | Ratio error, monolithic | Loss/period, IQN-ILS | Loss/period, monolithic |
|---|---:|---:|---:|---:|---:|
| weak | 25 | +0.015% | +1.022% | 2.70% | 7.01% |
| weak | 50 | +0.010% | +0.275% | 0.47% | 3.84% |
| weak | 100 | +0.004% | +0.052% | 0.16% | 2.01% |
| moderate | 25 | +0.075% | +1.428% | 3.27% | 34.28% |
| moderate | 50 | +0.056% | +0.312% | 0.98% | 19.57% |
| moderate | 100 | +0.018% | -0.021% | 0.65% | 10.59% |
| strong | 25 | +0.126% | +0.305% | 3.77% | 54.12% |
| strong | 50 | +0.101% | -0.044% | 1.50% | 32.71% |
| strong | 100 | +0.033% | -0.207% | 1.15% | 18.25% |

- Frequency: monolithic time orders 2.04, 2.34, 2.22; agreement with IQN-ILS
  at 100 steps per period +0.05%, -0.04%, -0.24%. At the strong level the
  monolithic ratio error moves away from zero under time-step refinement.
- Damping: the monolithic amplitude loss grows with the added mass, falls at
  first order in the time step (observed orders 0.80, 0.71, 0.57, successive
  ratios approaching 2) and hardly depends on the mesh (moderate, 100 steps:
  12.45%, 10.59%, 9.90% on factors 1, 2, 4; IQN-ILS 1.22%, 0.65%, 0.36%).
- Interpretation: a first-order, added-mass-proportional numerical
  dissipation in the time discretisation of the monolithic interface
  coupling; the fluid and solid discretisations are shared with the
  partitioned runs. Root cause not yet identified. A previous review of the
  port noted that the solid residual's interface force averages a traction
  from an older time level, a candidate mechanism.
