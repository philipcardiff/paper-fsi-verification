# ringAddedMass: quasi-monolithic coupling

Evidence for the quasi-monolithic Newton coupling on the ring added-mass
problem (Section `sec:ring`), the only case in which the experimental
quasi-monolithic solver is used. The study found two defects in the solver,
both corrected; the data of every stage are kept.

## Provenance

- solids4foam branch `feature-ring-monolithic-verification` (PR #582), on
  `feature-petsc-snes-quasi-monolithic-update` (PR #580): study driver
  `dcbe04505`, BDF2 solid `c008e8dd1`, interface pressure `3e93daef6`.
- OpenFOAM v2512, macOS (Apple M-series), serial runs, 10-14 at a time; the
  execution times are therefore indicative only.
- Command: `tutorials/fluidSolidInteraction/ringAddedMass/verification/Allverify
  --study monolithic --cores 10`.
- Data (refinement-factor-2 mesh, standard solid, 25/50/100 steps per period,
  three added-mass levels, unless stated):
  - `ringAddedMass_monolithic_runs.csv`: corrected solver (final).
  - `ringAddedMass_monolithic_mesh.csv`: corrected solver, moderate and strong
    levels, meshes 1/2/4, 100 steps per period, with IQN-ILS and dry runs.
  - `*_bdf2.csv`: after the first correction only (BDF2 solid).
  - `*_lagged.csv`: first version (built from PR #581, serial results
    bit-identical to PR #580).

## Set-up

The monolithic velocity solid has the standard second-order residual only, so
every run uses the standard solid and the comparison is through the wet/dry
frequency ratio, which removes the solid error. The dry ring with the
velocity solid gives the same frequencies as the total-displacement solid on
all meshes (3.08290, 2.71258, 2.59528 rad/s). The monolithic system converges
in two Newton iterations per step (linear problem).

## Defect 1: one-step lag of the interface force (time)

- First version: the solid residual averaged the stress at the old and an
  extrapolated future time; on the interface the old-time half was the
  traction from two steps before, so the interface force lagged by one step.
  Damper of coefficient m_a ω² Δt, loss π (m_a/M) ω Δt per period:
  1.8/9.9/17.9% predicted, 2.0/10.6/18.3% measured at 100 steps per period.
- Removing only the lag gave a centred scheme with |amplification| = 1 at
  every ω Δt; a mode with a period of about 4.5 steps then grew in
  blobInTreacle.
- Correction: backward (BDF2) displacement, D_k = (4 D_{k-1} − D_{k-2} +
  2 Δt U_k)/3, stress and interface traction once at the new time.
- Diagnosis: Codex Astra and Fable, independently.

## Defect 2: zero-gradient interface pressure (space)

- After defect 1: loss as IQN-ILS, but ratio error -0.194/-0.117/-0.063%
  (moderate) and -0.355/-0.213/-0.115% (strong) on meshes 1/2/4, about first
  order; excess added mass of about a quarter of a fluid cell.
- Cause: the zeroGradient pressure halves the normal pressure gradient that
  the least-squares gradient gives in the wall cells, and makes the cell
  pressure the wall pressure; two first-order errors that partly cancel
  (Fable's model of the discrete fluid on the actual mesh).
  Ruled out by experiment: interface Rhie-Chow flux (no change), solid model
  (dry velocity solid identical). Sensitivity to the Rhie-Chow scale factor
  (-0.26/-0.18/-0.08% at 0.5/1/2, mesh 1) pointed to the wall-cell closure.
- Correction: fixedGradient interface pressure with
  dp/dn = n . (nu lap U − a), a the interface acceleration, set by the
  coupling; interface velocity extrapolated from the solid cell centres.

## Findings (corrected solver)

| Level | Steps | Ratio error, IQN-ILS | Ratio error, monolithic | Loss/period, IQN-ILS | Loss/period, monolithic | Loss/period, first version |
|---|---:|---:|---:|---:|---:|---:|
| weak | 100 | +0.004% | -0.000% | 0.16% | 0.045% | 2.01% |
| moderate | 100 | +0.018% | -0.002% | 0.65% | 0.040% | 10.59% |
| strong | 25 | +0.126% | -0.007% | 3.77% | 2.57% | 54.12% |
| strong | 50 | +0.101% | -0.005% | 1.50% | 0.36% | 32.71% |
| strong | 100 | +0.033% | -0.004% | 1.15% | 0.049% | 18.25% |

- The monolithic loss equals the dry ring's (2.57/0.36/0.046%): the coupling
  adds no dissipation.
- Mesh (100 steps): ratio error -0.016/-0.002/-0.0001% (moderate),
  -0.031/-0.004/-0.0004% (strong) on meshes 1/2/4; IQN-ILS
  +0.069/+0.018/+0.005% and +0.122/+0.033/+0.009%. The apparent order of
  about 3 exceeds the formal second order because the finest values reach the
  precision of the frequency measurement.
