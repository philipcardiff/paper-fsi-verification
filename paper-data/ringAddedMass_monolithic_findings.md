# ringAddedMass: quasi-monolithic coupling

Evidence for the quasi-monolithic Newton coupling on the ring added-mass
problem (Section `sec:ring`), the only case in which the experimental
quasi-monolithic solver is used. The study found a defect in the solver, which
was then corrected; both sets of results are kept.

## Provenance

- solids4foam branch `feature-ring-monolithic-verification` (PR #582; corrected
  solver at commit `c008e8dd1`, study driver at `dcbe04505`), on
  `feature-petsc-snes-quasi-monolithic-update` (PR #580), which ports the
  quasi-monolithic solver of PR #233 onto `development`.
- OpenFOAM v2512, macOS (Apple M-series), serial runs, 14 at a time; the
  execution times are therefore indicative only.
- Command: `tutorials/fluidSolidInteraction/ringAddedMass/verification/Allverify
  --study monolithic --cores 14`.
- `ringAddedMass_monolithic_runs.csv`: corrected solver (backward solid), the
  driver's `ringAddedMass_results.csv` (21 runs: the dry ring and the IQN-ILS
  and monolithic couplings at three added-mass levels, 25, 50 and 100 steps per
  period, refinement-factor-2 mesh, standard solid).
- `ringAddedMass_monolithic_mesh.csv`: corrected solver, moderate and strong
  levels on the refinement-factor-1, 2 and 4 meshes at 100 steps per period.
- `ringAddedMass_monolithic_runs_lagged.csv` and
  `ringAddedMass_monolithic_mesh_lagged.csv`: the first version of the solver,
  built from `feature-petsc-snes-components` (PR #581, serial results
  bit-identical to PR #580).

## Set-up

The monolithic velocity solid has the standard second-order residual only, so
every run uses the standard solid and the comparison is through the wet/dry
frequency ratio, which removes the solid error (standard and high-order
ratios agree to 0.003% on the finest mesh, Section `sec:ring`). The
monolithic system converges in two Newton iterations per step (linear
problem).

## Defect and correction

- First version: the solid residual averaged the stress at the old and an
  extrapolated future time, with the BDF2 acceleration; on the interface the
  old-time half was the traction from two steps before, so the interface force
  lagged by one step. For an added-mass force this is a damper of coefficient
  m_a ω² Δt, giving a loss of π (m_a/M) ω Δt per period: 1.8/9.9/17.9%
  predicted, 2.0/10.6/18.3% measured at 100 steps per period. Diagnosis
  confirmed independently by two reviewers (Codex Astra and Fable).
- Removing the lag gave a centred scheme with |amplification| = 1 at every
  ω Δt (no dissipation; stiff limit λ → ±i, a period of about 4 steps). The ring
  then matched IQN-ILS, but in blobInTreacle a mode with a period of about
  4.5 steps grew after the steady state was reached.
- Correction: the solid integrates its displacement with the backward (BDF2)
  scheme of ddt(U), D_k = (4 D_{k-1} − D_{k-2} + 2 Δt U_k)/3, and evaluates the
  stress and interface traction once, at the new time; the Adams-Bashforth
  displacement only moves the fluid mesh before the solve. blobInTreacle is
  then stable with the same steady state.

## Findings (corrected solver)

| Level | Steps | Ratio error, IQN-ILS | Ratio error, monolithic | Loss/period, IQN-ILS | Loss/period, monolithic (lagged) | Loss/period, monolithic (corrected) |
|---|---:|---:|---:|---:|---:|---:|
| weak | 100 | +0.004% | -0.021% | 0.16% | 2.01% | 0.17% |
| moderate | 100 | +0.018% | -0.117% | 0.65% | 10.59% | 0.71% |
| strong | 25 | +0.126% | -0.016% | 3.77% | 54.12% | 4.28% |
| strong | 50 | +0.101% | -0.122% | 1.50% | 32.71% | 1.85% |
| strong | 100 | +0.033% | -0.213% | 1.15% | 18.25% | 1.26% |

- Time: frequency orders 1.91, 2.01, 2.13; amplitude-loss orders 2.71, 2.33,
  2.05 (weak, moderate, strong).
- Mesh (100 steps per period): loss 1.37/0.71/0.39% (moderate) and
  2.41/1.26/0.69% (strong) on meshes 1/2/4, as IQN-ILS; ratio error
  -0.194/-0.117/-0.063% (moderate) and -0.355/-0.213/-0.115% (strong), about
  first order, against +0.069/+0.018/+0.005% and +0.122/+0.033/+0.009% for
  IQN-ILS (about second order). A first-order spatial error in the monolithic
  interface treatment, growing with the added mass, remains; not yet
  diagnosed.
