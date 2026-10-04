# beamInCrossFlow graded-mesh verification findings

## 1. Executive scientific conclusion

The graded family materially improves the allocation of fluid cells and gives
three-level monotone trends for the original-case displacements and `F_x` at a
fixed time step. Their local observed orders are about 1.0 for `u_x(A)` and
`u_y(A)` and 1.5 for `F_x`, but the 8x calculation reaches the unchanged
100-iteration coupling limit, so the available evidence does not establish an
asymptotic range. `F_y` remains non-monotone. The hypothesis is therefore only
partly confirmed: near-body grading improves the principal trends efficiently,
but force cancellation and coupling/solid-interface sensitivities remain.

The study does establish two conclusions independently of the final numerical
trend. First, the old global-uniform sequence is a poor allocation of cells:
two thirds of its fluid cells lie in the downstream blocks while the first
cell normal to the beam remains the global spacing. Second, the commonly used
Richter numbers are not a fully matching independent reference for the present
solids4foam problem. The source problem differs in beam position, monitoring
point interpretation and inlet peak speed.

## 2. Provenance

- solids4foam branch: `verification/beam-crossflow-graded-mesh`
- solids4foam commit: `7d1012ada277bff935a89a6e7472155937c0d933`
- starting commit: `aee0c8e3502675749ae89185cc957633ebd6d854`
  (`development`, fetched 2026-10-04)
- OpenFOAM: OpenCFD v2512
- solids4foam solver: `solids4Foam`, PETSc-enabled build
- platform: `xenosim`, Linux 5.15, two AMD EPYC 9684X 96-core sockets
- date: 2026-10-04
- coupling: partitioned IQN-ILS with direct interface mapping, predictor,
  relative mode filter `0.01`, outer tolerance `1e-6`, at most 100 iterations
- time integration: BDF2 (`backward`) in fluid and solid
- production commands:

  ```bash
  source ~/bin/load-openfoam v2512
  ./verification/Allverify --case original --study mesh --family graded
  ./verification/Allverify --case original --study mesh --family graded \
      --levels 4,8 --cores 64
  ./verification/Allverify --case original --study temporal --family graded \
      --cores 64
  ./verification/Allverify --case modified --study mesh --family graded \
      --levels 2,4 --cores auto
  ```

The first original sweep was stopped after L1 and resumed from L2 to avoid
repeating completed levels. The first modified attempt used 64 ranks and
failed on L1; the listed retry uses the driver's validated factor-dependent
rank counts. Both failures are part of the evidence rather than discarded
solver tuning.

The commands were run from
`tutorials/fluidSolidInteraction/beamInCrossFlow`. Compact source results are
versioned beside the verification driver in solids4foam and copied to this
paper repository. No raw field directories are committed.

## 3. Existing uniform-mesh evidence

The old family refines every fluid and solid block uniformly. Its base mesh has
14,592 fluid and 256 solid cells with `0.025 m` spacing: four cells through the
beam thickness, eight cells in each half-interface direction, and one
`0.025 m`-wide first cell at the beam (centre distance `0.0125 m`). At 8x it
has 7,602,176 total cells and
`0.003125 m` spacing. Its time step was reduced with its mesh spacing, from
`0.05 s` to `0.00625 s`, so it is a combined space--time path.

The versioned plots and recovered function-object records reproduce the old
trends. Original `u_x(A)` is approximately
`3.922e-5, 5.028e-5, 5.69e-5, 5.954e-5 m`; its final change is about 4.6%.
Modified `u_x(A)` is
`0.0097815, 0.0121414, 0.013145 (3x), 0.0136386, 0.0143007 m`; its final
change is about 4.9%. Original `F_y` is non-monotone. The old run did not
version its CSV, and not every final force component is recoverable; this is a
data-retention discrepancy, recorded explicitly in the reconstructed
historical CSV rather than silently filled.

Force-component records explain why force convergence need not mirror
displacement. For the modified 1x--4x sequence, pressure `F_y` decreases from
0.1714 to 0.1281 N while viscous `F_y` increases from 0.0555 to 0.1121 N. Their
sum therefore rises and then falls. `F_x`, by contrast, is pressure-dominated
and rises monotonically over the recovered levels. This cancellation is direct
evidence for near-body sensitivity, although it does not prove that fluid
resolution is the only limiting error.

A field diagnostic on the completed original L2 result reinforces that
interpretation. The maxima of `|grad(U)|=34.13 1/s` and kinematic
`|grad(p)|=6.412 m/s2` (equivalent to `6.412 kPa/m` at the specified density)
occur in the same fluid cell centred at approximately
`(0.4487,0.1969,-0.1969) m`, immediately upstream of the beam's free outer
corner. The controlling feature is therefore demonstrably local to the body,
but its location also implicates the sharp corner/free-end geometry rather
than proving that a smooth-wall boundary layer alone controls the error.

## 4. New graded mesh family

The new family keeps the original conformal 11-block topology. Its fixed total
block expansion ratios are `xUp=0.125`, `xDown=8`, `yOuter=6` and
`zOuter=1/6`. They cluster cells on the two flow-normal beam faces, around the
free end and side, and through the near wake, with coarsening toward external
boundaries. All block counts are multiplied by 1, 2, 4 and 8; the regions and
grading are not retuned between levels. The solid and interface counts use the
same factor.

| Level | Factor | Fluid | Solid | Near (m) | Far (m) | Thick. | Face y x z |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| L0 | 1 | 7,584 | 256 | 0.0108071 | 0.0927085 | 4 | 8 x 8 |
| L1 | 2 | 60,672 | 2,048 | 0.00548932 | 0.0466999 | 8 | 16 x 16 |
| L2 | 4 | 485,376 | 16,384 | 0.00276512 | 0.0234344 | 16 | 32 x 32 |
| L3 | 8 | 3,883,008 | 131,072 | 0.00138756 | 0.0117380 | 32 | 64 x 64 |

The representative controlling spacing is close to halved at every level.
The corresponding first-cell-centre distances from a locally planar beam face
are approximately 0.005404, 0.002745, 0.001383 and 0.000694 m. “Near” in the
table is the full interface-adjacent cell width, not a wall-function `y+`
quantity.
`checkMesh -allTopology -allGeometry` reports `Mesh OK` in both regions at L0
and L3. Fluid maximum aspect ratio changes from 8.159 to 7.904; maximum
non-orthogonality is zero and skewness remains at round-off level. Thus a
changing quality defect is not a plausible explanation for the convergence
trend. L3 decomposition over 64 ranks is also balanced: per-rank fluid cell
counts vary by only 2.02% and solid counts by 1.97%. The L3 coupling failure
cannot be attributed to a gross cell-count imbalance.

### Exact problem definition for independent reproduction

The fluid computational half-domain is
`[0,1.5] x [0,0.4] x [-0.4,0] m`; the physical domain is its reflection about
`z=0`. The undeformed solid is
`[0.45,0.55] x [0,0.2] x [-0.2,0] m`. The `y=0` solid face is clamped and
`z=0` is a solid/fluid symmetry plane. The inlet is at `x=0`, zero-gauge
pressure outlet at `x=1.5`, and all other channel walls are no-slip. The inlet
is parabolic in `y` and `z`, with the cosine ramp implemented by
`paraboloidInletVelocity`.

Fluid properties are `rho=1000 kg/m3` and `nu=0.001 m2/s`. The solid has
`rho=1000 kg/m3`, Poisson ratio 0.4 and St Venant--Kirchhoff elasticity. The
original form uses `E=1.4 MPa`, peak inlet speed `0.2 m/s`, and a ramp ending
at `4 s`; the modified form uses `E=10 kPa`, `0.3 m/s`, and a ramp ending at
`1 s`. Both are evaluated at `8 s`.

Point A is `(0.45,0.15,-0.15) m`. Its three displacement components are
interpolated by `solidPointDisplacement`. Forces are the pressure-plus-viscous
integrals exerted by the fluid on the half-beam interface: positive `x` is
downstream, positive `y` upward. `2u_z(A)` is used only when comparing with the
published symmetry-paired transverse difference. A future COMSOL calculation
should reproduce these definitions literally and report both pressure and
viscous force components. The expected steady check is the relative change of
each principal QoI over `t=7...8 s`, in addition to the `t=8 s` value; this
study treats a change below 0.5% as steady for the reported precision.

## 5. Original case results

| L | Cells | `u_x(A)` (m) | change | `u_y(A)` (m) | `u_z(A)` (m) |
| --- | ---: | ---: | ---: | ---: | ---: |
| 0 | 7,840 | 4.08348e-5 | -- | 1.56846e-5 | -3.18602e-7 |
| 1 | 62,720 | 5.21672e-5 | 27.75% | 2.10399e-5 | -7.38201e-7 |
| 2 | 501,760 | 5.79343e-5 | 11.06% | 2.36525e-5 | -9.26181e-7 |
| 3 | 4,014,080 | not converged | -- | not converged | not converged |

| L | `F_x` (N) | change | `F_y` (N) | `F_z` (N) |
| --- | ---: | ---: | ---: | ---: |
| 0 | 1.197706 | -- | 0.1077488 | -0.0417130 |
| 1 | 1.279401 | 6.82% | 0.1065742 | -0.0463250 |
| 2 | 1.308220 | 2.25% | 0.1082129 | -0.0504448 |
| 3 | not converged | -- | not converged | not converged |

L0--L2 are monotone for every listed quantity except `F_y`. The three-point
local orders are 0.975 for `u_x(A)`, 1.035 for `u_y(A)`, 1.158 for `u_z(A)`
and 1.503 for `F_x`. The `F_z` order is only 0.163 and its changes shrink too
slowly to support extrapolation. `F_y` reverses direction and has no observed
order. L2 is 2.20% below Richter's Table 6 extrapolated displacement and 1.42%
below its extrapolated `F_x`, but those agreements cannot be treated as error
estimates because the benchmark definitions differ. No Richardson value is
reported: the last displacement changes are still 11--25%, and L3 did not
complete the first time step within the unchanged coupling limit.

At identical solid/interface factors, replacing the old uniform fluid mesh by
the graded fluid mesh raises `u_x(A)` by 4.11%, 3.76% and 1.82% on L0--L2;
the corresponding `F_x` changes are -0.53%, 1.87% and 1.10%. This decreasing
but non-negligible family difference is evidence that near-body fluid
resolution contributed to the old displacement trend. It is not evidence that
fluid resolution was the only limitation: completed graded L2 has a slightly
finer near-body cell than old uniform 8x, but only half its solid/interface
factor, and its displacement remains 2.70% lower.

## 6. Modified large-deformation case results

| L | Cells | `u_x(A)` (m) | `u_y(A)` (m) | `u_z(A)` (m) |
| --- | ---: | ---: | ---: | ---: |
| 0 | 7,840 | 0.0100906 | 0.00350698 | -0.000103346 |
| 1 | 62,720 | 0.0126639 | 0.00444752 | -0.000190036 |
| 2 | 501,760 | 0.0139699 | 0.00487255 | -0.000219652 |

| L | `F_x` (N) | `F_y` (N) | `F_z` (N) |
| --- | ---: | ---: | ---: |
| 0 | 2.164837 | 0.274350 | -0.261011 |
| 1 | 2.296244 | 0.239080 | -0.272527 |
| 2 | 2.344948 | 0.229769 | -0.281750 |

Every listed modified QoI is monotone over L0--L2. Local orders are 0.978,
1.146, 1.549, 1.432, 1.921 and 0.320 for `u_x`, `u_y`, `u_z`, `F_x`, `F_y`
and `F_z`, respectively. The finest changes remain 9.6--15.6% for
displacements but only 2.1--3.9% for forces; these are promising trends, not
an asymptotic demonstration. L2 is 4.51% below Tukovic `u_x`, 2.55% below
`u_y`, and its symmetry-paired `2u_z` is 1.72% below the Figure 28 value.
Thus the sequence approaches the old same-lineage values, but an independent
limit is still absent. At L2 all point components and `F_x`/`F_y` change by at
most 0.42% over `t=7...8 s`; `F_z` still changes by 0.93%.

## 7. Displacement versus force convergence

The graded family makes `u_x(A)`, `u_y(A)` and pressure-dominated `F_x`
monotone in the original case, but not at the same rate: `F_x` changes by
6.82% then 2.25% and has local order 1.50, whereas `u_x(A)` changes by 27.75%
then 11.06% and has order 0.97. Thus agreement of one force does not imply
displacement convergence. `F_y` is more difficult because decreasing pressure
and increasing viscous contributions partially cancel; its total falls by
1.09% and then rises by 1.54%. This remains unresolved and no order is fitted.

Observed orders are reported only for three consecutive values whose two
successive differences have the same sign. Blank orders mean that the sequence
is non-monotone or that fewer than three comparable levels exist. A Richardson
extrapolation is not reported unless the finest three levels have a stable,
credible order and demonstrably shrinking changes.

## 8. Temporal-error assessment

On original L2, halving `deltaT` from 0.0125 to 0.00625 s changes `u_x(A)`,
`u_y(A)`, `F_x` and `F_y` by 0.035%, 0.039%, 0.020% and 0.071%, respectively.
The changes in the smaller `u_z(A)` and `F_z` components are 0.203% and 0.099%.
These are all much smaller than the L1--L2 spatial changes, so temporal error
is demonstrably subordinate for the reported original-case sequence. This
one-mesh check does not prove the same statement independently for every
modified-case level.

The mesh family itself holds `deltaT=0.00625 s` on all four levels. It is
therefore a spatial sequence, not the combined space--time path used by the old
uniform family.

## 9. Literature/reference provenance

Richter's open 2012 source, *Goal-oriented error estimation for
fluid--structure interaction problems* (doi:10.1016/j.cma.2012.02.014), gives
five globally refined linear finite-element results in Table 5. The finest
level has 7,600,775 unknowns, `F_x=1.3380 N` and
`u_x(A)=5.9202e-5 m`. Table 6 extrapolates the finest three levels to
`F_x=1.327 +/- 0.01 N` and `u_x(A)=5.924e-5 +/- 1e-7 m`, with a claimed
relative accuracy of at most 1%. Thus the commonly quoted `1.33 N` and
`5.95e-5 m` are rounded benchmark values, not one refinement level.

The configuration is not exact. Richter places the solid at
`x=0.4...0.5 m`, samples `A=(0.45,0.15,+0.15) m` within that thickness, and
states peak inlet speed `0.3 m/s`. Solids4foam places the solid at
`x=0.45...0.55 m`, samples its upstream face at
`A=(0.45,0.15,-0.15) m`, and sets peak speed `0.2 m/s` in the original form.
The sign of `z` is only the symmetry-half convention; the consequential point
difference is `x`: Richter's point is at mid-thickness, whereas the current
solids4foam point lies on the upstream face. Its transient ramp is also absent
from Richter's stationary definition. The reference is independent FE evidence
for a closely related benchmark, not an exact-reference validation of this
case. This study preserves the current point so that the mesh effect is not
confounded with a benchmark-definition change; the point should be reconciled
in separate follow-up work.

Tukovic et al. (2018, doi:10.21278/TOF.42301) report the solids4foam-lineage
values used for both variants. They used a locally concentrated unstructured
mesh (273,539 fluid and 6,661 solid cells) but did not publish a full spatial
convergence sequence; the modified values are read from Figure 28.

Gillebaart's 2016 dissertation (doi:10.4233/uuid:c078909a-a39c-47b5-8ca6-
2cbc9e04486e), section 2.3.4, has the same domain dimensions, beam dimensions,
modified material properties, peak speed and matching-interface IQN-ILS
setup. It is nevertheless not the same transient problem: the solid is held
fixed until `5 s`, loads are ramped from `5` to `6 s`, and only then are full
fluid loads applied. It studies temporal field errors and forces but does not
tabulate matching steady point-A QoIs. It cannot replace an independent
reference for the present modified case.

## 10. Cost comparison

| Run | Cells | MPI ranks | Wall time | Steps | Mean/max coupling iterations |
| --- | ---: | ---: | ---: | ---: | ---: |
| Original L0 | 7,840 | 1 | 595 s | 1,280 | 6.06 / 10 |
| Original L1 | 62,720 | 4 | 1,477 s | 1,280 | 5.77 / 10 |
| Original L2 | 501,760 | 64 | 3,561 s | 1,280 | 5.98 / 14 |
| Original L3 | 4,014,080 | 64 | failed, step 1 | 0 | 100 at failure |
| Original L2, `dt=0.0125 s` | 501,760 | 64 | 3,236 s | 640 | 8.06 / 20 |
| Modified L0 | 7,840 | 4 | 315 s | 1,280 | 6.66 / 12 |
| Modified L1 | 62,720 | 4 | 1,770 s | 1,280 | 6.92 / 13 |
| Modified L2 | 501,760 | 32 | 5,095 s | 1,280 | 8.38 / 39 |

Completed graded L2 is the clearest like-for-like cost result: its `0.002765 m`
near-body spacing is 13% finer than old uniform 8x (`0.003125 m`), yet it has
501,760 rather than 7,602,176 cells, a factor of 15.2 reduction. On the same
64-rank machine it took about 0.99 h, versus the recorded approximately 11 h
for old uniform 8x. These are not identical mesh topologies, but their time
steps are both `0.00625 s` and the comparison directly measures the intended
reallocation of resolution.

New L3 has 4,014,080 total cells, 47.2% fewer than old uniform 8x, while its
representative near-body spacing is 2.25 times smaller. A uniform 16x mesh
would have spacing `0.0015625 m` (slightly coarser than new L3) and about 60.8
million total cells. The graded mesh therefore uses about 15 times fewer cells
than a uniform mesh with comparable controlling resolution. Cell-count and
runtime savings do not by themselves establish a better convergence rate.

## 11. Remaining uncertainties

The fourth original level and the first over-decomposed modified L1 attempt
failed the unchanged 100-iteration coupling limit. These failures are retained
as evidence and no tolerance was relaxed to manufacture a trend. Three
completed levels give only one local-order estimate, not proof of a stable
asymptotic range. The original transverse force remains non-monotone, and the
small out-of-plane components are not convincingly converged.

The most fundamental uncertainty is reference definition, not solver
repeatability: neither form currently has an independent result for the exact
solids4foam geometry, point and loading history. The clamp-edge singularity,
solid/interface order and force quadrature can also limit force convergence
after fluid near-body resolution is improved.

## 12. Consequences for the manuscript

The evidence strengthens the claim that global uniform refinement is an
inefficient route and that displacement and integrated force must be assessed
separately. It weakens any claim that the commonly quoted Richter values are
an exact independent validation, because the beam location, point-A meaning,
inlet speed and solution history differ. The original form now has useful
three-level spatial evidence but remains provisional rather than a strong
validation case; the modified form deserves still weaker status because its
reference is same-lineage. Until the exact benchmark definition has an
independent reference and the finest-level coupling sensitivity is resolved,
`beamInCrossFlow` should not carry a central main-text verification claim.

No manuscript `.tex` file was changed in this work.

## 13. Recommended remaining work

### Essential

- Resolve the L3 coupling failure without changing physical or spatial
  discretisation settings: first test decomposition sensitivity and then, if
  needed, diagnose the IQN-ILS history/initialisation. Do not merely relax the
  convergence tolerance.
- Repeat enough fine levels to decide whether the L0--L2 local orders persist;
  do not use a Richardson estimate until they do.
- Reconcile the exact monitoring-point and geometry definition with the
  intended benchmark. The present study deliberately retains the current
  solids4foam point, `(0.45,0.15,-0.15) m`, on the upstream face.

### High-value

- Generate an independent COMSOL result for the *exact* solids4foam geometry
  and definitions above, prioritising the modified case and including a mesh
  and time-step study.
- Separate solid/interface resolution from fluid resolution by holding the
  fine graded fluid mesh fixed and refining only the solid/interface where
  feasible.
- Retain pressure and viscous force components at every level; their
  cancellation is needed to interpret `F_y`.

### Optional

- Test a second, smoothly graded topology or goal-oriented local adaptation to
  distinguish topology bias from resolution.
- Compare force evaluation on the moving interface with an equivalent clamp
  reaction when available.
- Archive compact histories and selected final fields in the future Zenodo
  evidence package.

## 14. Suggested manuscript wording

> A parametrically graded beam-in-cross-flow family reduced cell count by a
> factor of about 15 at comparable near-body spacing. Three completed levels
> gave monotone displacement and streamwise-force trends in both forms, with
> local orders near one for streamwise displacement and 1.4--1.5 for
> streamwise force. Original transverse force remained non-monotone, and the
> finest original level reached the coupling-iteration limit. The case is
> therefore retained as provisional numerical evidence, not exact validation;
> the commonly quoted Richter values correspond to a different beam position,
> mid-thickness sample point, inlet speed and stationary formulation.

## Explicit answers to the scientific questions

1. **Was the old uniform family the reason convergence was poor?**
   It was an important efficiency problem, not a demonstrated sole cause.
   Grading changes the solution at matched solid resolution and delivers
   comparable near-body spacing with 15.2 times fewer cells, but residual
   solid/interface and coupling sensitivities remain.
2. **Does the graded family converge monotonically?**
   Original `u_x`, `u_y`, `u_z`, `F_x` and `F_z` do over L0--L2; original
   `F_y` does not. All modified L0--L2 QoIs are monotone. Original L3 has no
   solution because the coupling limit is reached.
3. **Which QoIs show credible observed orders?**
   Original `u_x(A)`, `u_y(A)` and `F_x` have useful local orders of 0.975,
   1.035 and 1.503. They are provisional three-point estimates, not confirmed
   asymptotic orders. The small `u_z` component gives 1.158 but is less robust.
4. **Which QoIs remain unresolved?**
   `F_y`, `F_z`, the small out-of-plane displacement, and every L3 trend;
   modified-case conclusions are additionally limited by their weak reference.
5. **Does `F_x` converge differently from `u_x(A)`?**
   Yes. Original L1--L2 changes are 2.25% for `F_x` versus 11.06% for
   `u_x(A)`, with local orders 1.50 and 0.97 respectively.
6. **What happens to `F_y`?**
   It falls and then rises in the original graded sequence, consistent with
   cancellation between oppositely evolving pressure and viscous components;
   no order or extrapolation is defensible.
7. **Is original Richter agreement robust or partly coincidental?**
   It is partly coincidental as validation evidence. L2 is close numerically,
   but Richter's beam position, mid-thickness point, inlet speed and stationary
   formulation do not match the present upstream-face transient case.
8. **Does the modified case converge to the old Tukovic value?**
   It trends toward it: L2 errors are 4.51% (`u_x`), 2.55% (`u_y`) and 1.72%
   (`2u_z`). The remaining 9.6--15.6% mesh changes and same-code provenance
   prevent treating that destination as independently established.
9. **Is temporal error negligible?**
   It is subordinate on original L2: principal-QoI changes on halving the time
   step are 0.020--0.071%, far below spatial changes. This has not been proved
   independently for the modified form.
10. **How much cheaper is the graded family?**
    At slightly finer near-body spacing, completed graded L2 uses 15.2 times
    fewer cells and about 11 times less 64-rank wall time than old uniform 8x.
11. **Is this now a strong main-paper verification case?**
    No. It is useful numerical evidence, but the reference-definition mismatch,
    non-monotone force, and failed finest coupling level keep it provisional.
12. **Is independent COMSOL still valuable?**
    Yes. An exact-definition COMSOL mesh/time study remains high-value,
    especially for the modified case, and should report the same point and
    separated pressure/viscous forces.
