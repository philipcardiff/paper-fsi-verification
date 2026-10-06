# beamInCrossFlow graded-mesh verification findings

## 1. Executive scientific conclusion

The graded family materially improves the allocation of fluid cells and gives
four completed fixed-time-step levels for both case forms. On the finest
refinement, original `u_x(A)` and `F_x` change by 3.83% and 0.69%; modified
`u_x(A)` and `F_x` change by 3.60% and 0.64%. The mixed-platform fine-level
apparent orders are about 1.4 for `u_x(A)` and 1.7 for `F_x` in both forms.
Dedicated L2 controls show that version/platform differences are at most 0.29%
for `u_x(A)` and 0.062% for `F_x`, supporting those principal trends. The
small `u_z(A)` is strongly platform-sensitive, original `F_y` remains
non-monotone over the full family, and no Richardson extrapolation is
justified. The hypothesis is therefore confirmed as an efficiency diagnosis,
but not as proof that near-body fluid resolution was the sole error source.

The study does establish two conclusions independently of the final numerical
trend. First, the old global-uniform sequence is a poor allocation of cells:
two thirds of its fluid cells lie in the downstream blocks while the first
cell normal to the beam remains the global spacing. Second, the commonly used
Richter numbers are not a fully matching independent reference for the present
solids4foam problem. The source problem differs in beam position, monitoring
point interpretation and inlet peak speed.

## 2. Provenance

- solids4foam branch: `verification/beam-crossflow-graded-mesh`
- solids4foam commit: `f6ff138cfe624679b9e1b0648b45af143d43cdd4`
  (results update; branch history contains
  mesh implementation commit `7d1012ada277bff935a89a6e7472155937c0d933`)
- starting commit: `aee0c8e3502675749ae89185cc957633ebd6d854`
  (`development`, fetched 2026-10-04)
- OpenFOAM: OpenCFD v2512 for L0--L2 and temporal data; OpenCFD v2412 for L3
  and platform controls
- solids4foam solver: `solids4Foam`, PETSc-enabled build
- platforms: `xenosim`, Linux 5.15, two AMD EPYC 9684X 96-core sockets;
  MeluXina Slurm CPU nodes for L3 and controls
- dates: 2026-10-04 to 2026-10-06
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

L3 and platform controls were subsequently run on MeluXina on 2026-10-04 to
2026-10-06 using OpenFOAM v2412, PETSc 3.22 and OpenMPI 5.0.3. Both production
L3 cases used 128 ranks, `deltaT=0.00625 s`, and otherwise identical case,
scheme, tolerance and coupling settings. The Slurm production jobs were
`5293133` and `5293134`; the L2 controls were `5303748` and `5307350`.

The first original L3 attempt on v2512/64 ranks reached the 100-iteration cap
at its first step. On v2412, one-step diagnostics converged in 36 iterations
on 32 ranks and 23 on 128 ranks; the unchanged 128-rank production case then
completed. A modified v2412 L2 control stalled at `t=1.68125 s` on 32 ranks,
whereas its unchanged 128-rank repeat completed. These failures are retained
as decomposition/platform sensitivity evidence rather than hidden by solver
tuning.

The commands were run from
`tutorials/fluidSolidInteraction/beamInCrossFlow`. Compact source results are
versioned beside the verification driver in solids4foam and copied to this
paper repository as `beam_crossflow_graded_mesh_results.csv` and
`beam_crossflow_graded_platform_controls.csv`. No raw field directories are
committed.

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
| 3 | 4,014,080 | 6.01533e-5 | 3.83% | 2.44087e-5 | -1.73881e-6 |

| L | `F_x` (N) | change | `F_y` (N) | `F_z` (N) |
| --- | ---: | ---: | ---: | ---: |
| 0 | 1.197706 | -- | 0.1077488 | -0.0417130 |
| 1 | 1.279401 | 6.82% | 0.1065742 | -0.0463250 |
| 2 | 1.308220 | 2.25% | 0.1082129 | -0.0504448 |
| 3 | 1.317197 | 0.69% | 0.1087008 | -0.0516540 |

L0--L3 are monotone for every listed quantity except `F_y`. The L1--L3
apparent orders are 1.378 for `u_x(A)`, 1.789 for `u_y(A)`, 1.683 for `F_x`
and 1.768 for `F_z`. They mix v2512 L1--L2 with v2412 L3, so they are
platform-qualified rather than strict single-executable orders. The v2412 L2
control differs from v2512 by 0.26% in `u_x(A)` and only 0.0015% in `F_x`;
using that control, the final changes are 3.56% and 0.685%. Thus the principal
fine-level trends are robust. In contrast, L2 `u_z(A)` differs by 77% between
controls and its mixed-platform apparent order is negative. `F_y` reverses
direction over the full family, although its finest three values increase
with shrinking differences. L3 is 1.10% above the rounded Richter
displacement and 0.96% below its `F_x`, but those agreements are not error
estimates because the benchmark definitions differ. No Richardson value is
reported.

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
| 3 | 4,014,080 | 0.0144735 | 0.00496883 | -0.000431879 |

| L | `F_x` (N) | `F_y` (N) | `F_z` (N) |
| --- | ---: | ---: | ---: |
| 0 | 2.164837 | 0.274350 | -0.261011 |
| 1 | 2.296244 | 0.239080 | -0.272527 |
| 2 | 2.344948 | 0.229769 | -0.281750 |
| 3 | 2.359966 | 0.2280017 | -0.2837244 |

Every listed modified QoI is monotone over L0--L3. L2--L3 changes are 3.60%,
1.98%, 96.6%, 0.64%, 0.77% and 0.70% for `u_x`, `u_y`, `u_z`, `F_x`, `F_y`
and `F_z`. Mixed-platform L1--L3 apparent orders for all but `u_z` are 1.375,
2.142, 1.697, 2.397 and 2.224; `u_z` has negative order. The v2412 L2 control
differs from v2512 by 0.29% in `u_x`, 0.74% in `u_y`, 0.062% in `F_x`, 0.19%
in `F_y` and 0.27% in `F_z`, but by 87% in the small `u_z`. Thus the principal
displacement and force trends are robust, while the out-of-plane component is
not. L3 is 1.07% below Tukovic `u_x` and 0.62% below `u_y`; its
symmetry-paired `2u_z` is 93% larger in magnitude than the Figure 28 value.
The old apparent `u_z` agreement was therefore not robust. Over `t=7...8 s`,
L3 principal displacements and `F_x` change by at most 0.39%, while `F_y` and
`F_z` still change by 0.91% and 0.97%.

## 7. Displacement versus force convergence

The graded family makes `u_x(A)`, `u_y(A)` and pressure-dominated `F_x`
monotone in the original case, but not at the same rate: the final `F_x`
change is 0.69%, whereas the final `u_x(A)` change is 3.83%. The corresponding
modified changes are 0.64% and 3.60%. Thus integrated drag is substantially
closer to mesh independence than point displacement, and agreement of one
force does not imply displacement convergence. Original `F_y` is more
difficult because decreasing pressure and increasing viscous contributions
partially cancel; it falls on L0--L1 and then rises. Its full sequence is
non-monotone, so the positive finest-three apparent order is not promoted to a
Richardson estimate.

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
| Original L3 | 4,014,080 | 128 | 59,018 s | 1,280 | 6.82 / 37 |
| Original L2, `dt=0.0125 s` | 501,760 | 64 | 3,236 s | 640 | 8.06 / 20 |
| Modified L0 | 7,840 | 4 | 315 s | 1,280 | 6.66 / 12 |
| Modified L1 | 62,720 | 4 | 1,770 s | 1,280 | 6.92 / 13 |
| Modified L2 | 501,760 | 32 | 5,095 s | 1,280 | 8.38 / 39 |
| Modified L3 | 4,014,080 | 128 | 36,504 s | 1,280 | 10.37 / 75 |

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
than a uniform mesh with comparable controlling resolution. It completed both
forms without relaxing any coupling criterion. Cell-count and runtime savings
do not by themselves establish a better convergence rate.

## 11. Remaining uncertainties

All four production levels completed, but L3 required a different OpenFOAM
version/platform and decomposition after the original v2512/64-rank attempt
failed the unchanged 100-iteration limit. L2 controls bound the resulting
principal-QoI difference below 0.75%, but the orders remain mixed-platform
estimates rather than a strict single-executable sequence. Original transverse
force remains non-monotone over the full family, modified transverse forces
are still changing by about 1% over the final second, and the small
out-of-plane displacement is strongly platform-sensitive and unconverged.

The most fundamental uncertainty is reference definition, not solver
repeatability: neither form currently has an independent result for the exact
solids4foam geometry, point and loading history. The clamp-edge singularity,
solid/interface order and force quadrature can also limit force convergence
after fluid near-body resolution is improved.

## 12. Consequences for the manuscript

The evidence strongly supports the claim that global uniform refinement is an
inefficient route and that displacement and integrated force must be assessed
separately. It weakens any claim that the commonly quoted Richter values are
an exact independent validation, because the beam location, point-A meaning,
inlet speed and solution history differ. The original form now has useful
four-level numerical verification evidence, but remains provisional as a
validation case. The modified form has similarly strong internal mesh evidence
but an even weaker same-lineage reference. Until the exact benchmark definition
has an independent reference, `beamInCrossFlow` should not carry a central
main-text validation claim.

No manuscript `.tex` file was changed in this work.

## 13. Recommended remaining work

### Essential

- Reproduce L2--L3 with one OpenFOAM executable/platform before using the
  apparent orders for Richardson extrapolation.
- Diagnose why the fine cases and modified L2 control are decomposition-
  sensitive despite unchanged coupling settings; do not merely relax the
  convergence tolerance.
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
> factor of about 15 at comparable near-body spacing. Four levels completed in
> both forms. On the final refinement, streamwise displacement changed by
> 3.6--3.8%, whereas streamwise force changed by only 0.6--0.7%, demonstrating
> materially different displacement and integrated-force convergence. L2
> platform controls support these principal trends, although the finest-level
> orders mix OpenFOAM versions and the small out-of-plane displacement is not
> robust. Original transverse force remains non-monotone. The case is retained
> as provisional numerical evidence, not exact validation; the commonly quoted
> Richter values correspond to a different beam position, mid-thickness sample
> point, inlet speed and stationary formulation.

## Explicit answers to the scientific questions

1. **Was the old uniform family the reason convergence was poor?**
   It was an important efficiency problem, not a demonstrated sole cause.
   Grading changes the solution at matched solid resolution and delivers
   comparable near-body spacing with 15.2 times fewer cells, but residual
   solid/interface and coupling sensitivities remain.
2. **Does the graded family converge monotonically?**
   Original `u_x`, `u_y`, `u_z`, `F_x` and `F_z` do over L0--L3; original
   `F_y` does not. All modified L0--L3 QoIs are monotone, although monotonicity
   of the platform-sensitive `u_z` is not evidence of convergence.
3. **Which QoIs show credible observed orders?**
   The finest three values give useful apparent orders of 1.38 (`u_x`), 1.79
   (`u_y`) and 1.68 (`F_x`) for the original form, and 1.37, 2.14 and 1.70 for
   the modified form. They are mixed-platform estimates, not confirmed
   asymptotic orders. Force orders other than original `F_y` are also positive;
   `u_z` has negative order and is not credible.
4. **Which QoIs remain unresolved?**
   The small out-of-plane displacement, original full-family `F_y`, strict
   single-platform orders, and an exact independent reference. Modified-case
   conclusions are additionally limited by their weak reference and residual
   final-second transverse-force drift.
5. **Does `F_x` converge differently from `u_x(A)`?**
   Yes. Original L2--L3 changes are 0.69% for `F_x` versus 3.83% for
   `u_x(A)`; modified changes are 0.64% versus 3.60%. Drag is substantially
   closer to mesh independence than displacement.
6. **What happens to `F_y`?**
   It falls and then rises in the original graded sequence, consistent with
   cancellation between oppositely evolving pressure and viscous components.
   The finest three points have shrinking changes, but the full sequence is
   non-monotone, so no Richardson extrapolation is defensible. Modified `F_y`
   is monotone and changes by 0.77% on L3.
7. **Is original Richter agreement robust or partly coincidental?**
   It is partly coincidental as validation evidence. L2 is close numerically,
   but Richter's beam position, mid-thickness point, inlet speed and stationary
   formulation do not match the present upstream-face transient case.
8. **Does the modified case converge to the old Tukovic value?**
   Its principal components trend toward it: L3 errors are 1.07% (`u_x`) and
   0.62% (`u_y`). The symmetry-paired `2u_z` instead moves far from the old
   value and is platform-sensitive. Same-code provenance prevents treating the
   apparent destination as independently established.
9. **Is temporal error negligible?**
   It is subordinate on original L2: principal-QoI changes on halving the time
   step are 0.020--0.071%, far below spatial changes. This has not been proved
   independently for the modified form.
10. **How much cheaper is the graded family?**
    At slightly finer near-body spacing, completed graded L2 uses 15.2 times
    fewer cells and about 11 times less 64-rank wall time than old uniform 8x.
11. **Is this now a strong main-paper verification case?**
    Not yet as an independent validation case. It is now strong internal
    numerical evidence, but the reference-definition mismatch, mixed-platform
    finest order, original non-monotone force and unresolved `u_z` keep it
    provisional.
12. **Is independent COMSOL still valuable?**
    Yes. An exact-definition COMSOL mesh/time study remains high-value,
    especially for the modified case, and should report the same point and
    separated pressure/viscous forces.
