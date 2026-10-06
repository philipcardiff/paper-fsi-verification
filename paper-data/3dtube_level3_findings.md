# 3dTube level-3 verification findings

Evidence report for the `3dTube` benchmark (Section `sec:3dtube`). It records
what the third mesh level and the supporting runs show, for review before the
manuscript is changed. No `.tex` file is modified by this report.

> **Status (6 October 2026): Sections 1-11 are superseded where they give
> numbers or convergence verdicts.** The investigation of the cross-platform
> `u_z,min` discrepancy (Section 12) found two defects in solids4foam that
> affected the runs behind Sections 1-11: an aliasing miscompilation of the
> fluid wall viscous force (wrong axial wall shear on `-O3` builds, i.e. every
> xenosim run) and degenerate barycentric weights in the FSI point transfer
> (a smeared displacement transfer on mesh levels 2 and 3 on every
> platform, but not on level 1). Both are fixed. The corrected three-level
> results, convergence verdicts and consequences are in Sections 12.6-12.7
> and replace the corresponding parts of Sections 1, 4, 5, 8, 9 and 11. The
> time-discretisation (implicit-Euler) conclusion of Section 6 is confirmed
> with the corrected runs; the reference provenance (Section 3) is unchanged.
> Sections 1-11 are kept as the record of what the uncorrected runs showed.

## 1. Executive scientific conclusion

Level 3 (1.02 M fluid and 0.41 M solid cells, `Δt = 6.25e-6 s`) ran to
completion with every coupling criterion met, and changes the scientific
reading of the case in three ways:

1. **It does not establish a converging sequence for the headline QoI.** The
   peak radial displacement `u_r,max(A)` changes by +0.07% from level 1 to 2
   but by -0.81% from level 2 to 3. The apparent two-level convergence of
   `u_r,max` (0.05% in the manuscript) was a coincidence. The level-2-to-3
   change is spatial (it persists at a fixed time step and on a second
   platform) and is not iterative error (tightening every tolerance changes
   it by < 3e-6 relative). No observed order can be given for `u_r,max`,
   arrival time or pulse speed; only `u_z,min` decreases monotonically, with
   an order of about 0.7-1.1, well below the nominal 2.
2. **The time-discretisation explanation of the literature discrepancy is
   strengthened and sharpened.** Implicit Euler at the published
   `Δt = 1e-4 s` lowers this code's peak by 0.036-0.037 mm (levels 1 and 2),
   i.e. 92-114% of the 0.033-0.039 mm gap between the published implicit-Euler
   peaks and the backward solution (86-126% within the figure-reading
   accuracy); it also brings the axial trough and the reflected radial
   trough from 11-27% to 2-5% of the published values. The split is between
   time integrators, not between finite elements and finite volumes: Eken's
   fluid is a finite volume discretisation.
3. **A new, unresolved platform dependence affects the axial displacement.**
   Two platforms (OpenFOAM v2512/x86 here, v2412 EasyBuild on MeluXina)
   agree on `u_r,max`, arrival time and pulse speed to 0.1-0.2%, but `u_z,min`
   differs by 1.5% on level 1, growing to 3.0% on level 3. It is
   deterministic on each platform and not caused by decomposition, solver,
   tolerance, OpenFOAM release, mesh or uninitialised memory; its cause is
   not known. Until resolved, `u_z,min` should not carry a convergence or
   accuracy claim.

The temporal error at the level-2 and level-3 time steps (0.1-0.2%) is
small against the level-2-to-3 spatial changes (0.5-0.8%), so the scaled
refinement path is primarily a spatial study. Robin-Neumann and IQN-ILS
agree to 0.07% (QoIs) and 0.26% (history). 3dTube remains valuable as the
only 3-D transient case and as the case that explains the literature split,
but its self-convergence evidence is weaker than the manuscript currently
implies, and it should not be described as converged.

## 2. Provenance

- **Code:** solids4foam, https://github.com/solids4foam/solids4foam, branch
  `verification/3dtube-level3` (from `development` at `aee0c8e35`), commit
  `15a652f09` (full: 15a652f098d94d3c028332187c33663e559ee04a). All solver runs used commits whose solver and tutorial are
  identical to it (later commits change only the verification driver's
  evaluation and the README).
- **Main platform:** OpenFOAM v2512 (Ubuntu package, build
  `_87ed40d2-20251219`), PETSc 3.24 development build with Hypre, Open MPI,
  GCC; 192-core AMD node `xenosim` through Slurm; solids4foam built into a
  private library directory.
- **Replicate platform:** MeluXina (LuxProvide), AMD EPYC 7H12, OpenFOAM
  v2412 (EasyBuild `foss-2024a`), PETSc 3.22.0 with Hypre 2.31.0, Open MPI
  5.0.3, GCC 13.3.0.
- **Dates:** 4-5 October 2026.
- **Commands** (from `tutorials/fluidSolidInteraction/3dTube/verification`):

  ```bash
  ./Allverify --levels 1,2,3                         # scaled mesh study (auto ranks 4/8/32)
  ./Allverify --levels 1,2,3 --probe-z-shift 7.8125e-5   # c_p with off-face probes
  ./Allverify --levels 1,2,3 --delta-t 2.5e-5        # fixed-dt mesh study
  ./Allverify --levels 3 --delta-t 1.25e-5           # level-3 temporal check
  ./Allverify --levels 1,2,3 --tight-tolerances      # iterative-error check
  ./Allverify --study timestep                       # backward, level 1
  ./Allverify --study timestep --time-scheme Euler   # Euler, level 1
  ./Allverify --time-scheme Euler --levels 1,2 --delta-t 1e-4
  ./Allverify --study literature
  ./Allverify --study coupling                       # Robin vs IQN-ILS
  python3 scripts/analyse_3dTube_evidence.py         # compact evidence
  ```

- **Cost:** level 1 504 s on 4 ranks; level 2 3 641 s on 8; level 3
  22 502 s (6.25 h) on 32 ranks (720 core-hours); on MeluXina level 3 took
  17 512 s on 128 ranks.
- **Machine-readable results** (solids4foam repository,
  `tutorials/fluidSolidInteraction/3dTube/verification/results/`):
  `3dTube_evidence.json` (all values, changes, orders and their
  classification), `3dTube_mesh_study.csv`, `3dTube_runs.csv`,
  `3dTube_meluxina_runs.csv`. A copy of the mesh-study table is
  `paper-data/3dtube_level3_mesh_study.csv`.

## 3. Existing benchmark and reference status

Fernández-Moubachir form of the Formaggia et al. pressure-pulse tube: length
5 cm, inner radius 0.5 cm, wall 0.1 cm, clamped ends; fluid ρ = 1000 kg/m³,
ν = 3e-6 m²/s; wall ρ = 1200 kg/m³, E = 3e5 Pa, ν = 0.3 (Hookean, small
strain); inlet pressure 1333.2 Pa for 3 ms, ramped to zero by 3.1 ms; zero
outlet pressure; quarter tube with symmetry planes; run to 20 ms. Point A is
the inner wall at mid-length. Verified in the scripts: the mesh family keeps
the block topology (fluid O-grid with a fixed square core, solid 4 x 20 x 80
on level 1) and multiplies every block division by 2 per level, so the levels
are geometrically similar with ratio 2 (blockMesh arcs; the meshes of two
platforms agree to 2e-17 m); the time step halves with the mesh
(2.5e-5, 1.25e-5, 6.25e-6 s); Robin-Neumann coupling with displacement,
pressure and leakage-flux tolerances 1e-6, 1e-5, 1e-5; fluid relTol 1e-3,
solid Newton rtol 1e-6; backward (BDF2) in fluid and solid.

QoI definitions (unchanged): `u_r,max(A)` peak radial displacement, parabolic
fit through the largest sample and its neighbours; `u_z,min(A)` incident
axial trough (minimum for t < 6.5 ms, at ≈ 4.7 ms); `u_r,min(A), t > 14 ms`
reflected radial trough; `t_arr(A)` first time `u_r(A)` reaches half its first
peak; `c_p` least-squares slope of the half-inlet-pressure arrival on the
axis, z = 1-4 cm; `c_wall` the same from the wall displacement.

| Reference | Method | Time integration | Mesh | Extraction | Own convergence evidence | Independence |
|---|---|---|---|---|---|---|
| Tuković et al. (2018) | partitioned cell-centred FV (fluid and solid), IQN-ILS, StVK wall, full tube | backward, Δt 2.5e-5 (also 5e-5, 1e-4) | one mesh, 449 600 fluid / 288 000 solid hex | vector graphics of the PDF (exact) | time step only (peak 0.15203 → 0.15631 → 0.15694 mm); no mesh study | same code lineage: weak, not independent |
| Lozovskiy et al. (2019) | monolithic FE, P2-P1 Taylor-Hood fluid, P2 displacement, StVK, Ani3D, MUMPS; geometry and advection extrapolated (one linear solve per step) | implicit Euler, Δt 1e-4 | three tetrahedral meshes, finest 89 232 fluid / 38 016 solid | raster figure, ±3e-6 m, ±0.05 ms | three meshes: largest fine-to-finer history difference 0.7% (axial), 2.3% (radial); no time-step study | independent code |
| Eken (2016); Eken & Sahin (2016) | monolithic ALE; side-centred (staggered) unstructured FV fluid + Galerkin FE StVK wall; Newton-Krylov | first-order implicit Euler in the fluid, Δt 1e-4 ("to be consistent with Gee et al.") | hex meshes M1 (324 122 DOF), M2 (2 557 571 DOF) | raster figure (thesis Fig. 5.6, M2), ±3e-6 m | two meshes shown, no quantified study; no time-step study | independent code |

Provenance was checked against the original Lozovskiy et al. (2019) and
Eken & Sahin (2016, IJNMF) papers in this work; the Tuković entries are from
the repository's earlier extraction (the PDF could not be retrieved in this
session) and the Eken thesis figure was not re-read. Corrections to the
current text: Eken is **not** a finite element fluid; Lozovskiy et al. switch
the pressure off instantaneously and lag the geometry; neither FE paper
performed a time-step study, so neither published peak is converged in time.
Tuković stays a same-lineage, weak reference; Lozovskiy and Eken are
independent but unconverged in time and raster-digitised, so all three are
category C-weak.

## 4. Spatial / scaled-refinement results

Scaled path (mesh and Δt refined together), v2512. `c_p` from the
probe-shifted replicate (identical solution, see §8). Displacements mm, times
ms, speeds m/s.

| Quantity | L1 | L2 | L3 | Δ L1→L2 | Δ L2→L3 | Observed order | Richardson | Verdict |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| u_r,max(A) | 0.15964 | 0.15976 | 0.15847 | +0.07% | -0.81% | undefined | none | non-monotone, change grows: **unreliable** |
| t(u_r,max) | 7.202 | 7.207 | 7.254 | +0.07% | +0.65% | undefined | none | flat peak; extraction-limited |
| u_z,min(A) | -0.08954 | -0.08880 | -0.08835 | 0.83% | 0.51% | 0.73 | not justified (p far from 2) | monotone; platform uncertainty 1.5-3% exceeds the change (§8) |
| u_r,min(A), t>14 ms | -0.09998 | -0.09352 | -0.09566 | 6.7% | 2.3% | undefined | none | non-monotone |
| t_arr(A) | 5.872 | 5.857 | 5.830 | -0.26% | -0.46% | undefined | none | change grows |
| c_p | 4.651 | 4.632 | 4.591 | -0.41% | -0.88% | undefined | none | change grows; at the probe sampling noise (≈0.5-1%) |
| c_wall | 5.065 | 5.143 | 5.064 | +1.5% | -1.5% | undefined | none | non-monotone |

Fixed time step `Δt = 2.5e-5 s` on all levels (pure mesh refinement), v2512:
`u_r,max` 0.15964, 0.15961, 0.15833 (0.02%, 0.80%; undefined); `u_z,min`
-0.08954, -0.08863, -0.08818 (1.02%, 0.51%; order 1.01; formal Richardson
-0.0877 mm); `t_arr` 5.872, 5.853, 5.820 (0.33%, 0.56%; undefined). On
MeluXina: `u_r,max` 0.15978, 0.15927, 0.15820 (0.32%, 0.67%; undefined);
`u_z,min` -0.08819, -0.08638, -0.08551 (2.05%, 1.00%; order 1.06; formal
Richardson -0.0847 mm); `c_p` (shifted probes) 4.654, 4.639, 4.633.

Uncertainty bands (Roache GCI, `Fs = 3`, nominal `p = 2`, i.e. equal to the
last change): `u_r,max(A)` ±0.0013 mm (0.8%), `u_z,min(A)` ±0.00045 mm
(0.5%, before the platform term), `t_arr(A)` ±0.027 ms, `c_p` ±0.04 m/s.

Comparison with the references on level 3 (backward): `u_r,max` +1.0% and its
time +0.3% vs Tuković (largest history difference 0-10 ms: 2.1% of the peak);
`c_p` 4.6% below the 4.81 m/s thick-wall estimate (Tuković Eq. 36) and 1.1%
above their simulated 4.54 m/s (different method). Against Lozovskiy and Eken
the backward results differ by 28-33% (`u_r,max`), 11-19% (`u_z,min`) and
25-27% (late trough); these are not like-for-like (see §6).

**Answers.** Level 3 does not establish a credible converging sequence for
`u_r,max`, `t_arr` or `c_p`: their changes grow from L2→L3. Only `u_z,min`
behaves like a sequence entering an asymptotic range (monotone on both
paths and platforms), with a sub-nominal order of about 0.7-1.1. No order
other than the `u_z,min` one is defensible, and that one is qualified by the
platform dependence. `u_r,max` changes by 0.81% from level 2 to 3, `u_z,min`
by 0.51%. The pulse speed does not stabilise in the strict sense (changes
0.41% → 0.88%), but all `c_p` values lie in 4.59-4.65 m/s, and the remaining
changes are of the size of the probe sampling noise.

## 5. Temporal-error assessment

The default sweep is a combined space-time refinement, so the temporal error
was measured directly on each mesh (backward):

| Mesh | Δt change (s) | Δu_r,max | Δu_z,min | Δt_arr |
|---|---|---:|---:|---:|
| L1 | 2.5e-5 → 1.25e-5 | +0.07% | -0.01% | -0.12% |
| L2 | 2.5e-5 → 1.25e-5 | +0.09% | +0.19% | +0.07% |
| L3 | 2.5e-5 → 1.25e-5 → 6.25e-6 | +0.08%, +0.01% | +0.12%, +0.07% | +0.14%, +0.03% |
| L3 (MeluXina) | 2.5e-5 → 1.25e-5 → 6.25e-6 | +0.09%, +0.01% | +0.12%, +0.07% | +0.12%, +0.03% |

Temporal error at the level-2 and level-3 time steps is 0.1-0.2%, against
level-2-to-3 spatial changes at fixed Δt of 0.5-0.8%: the scaled path is
mesh-dominated (by a factor of about four or more), and the fixed-step
sequence leads to the same conclusions, so no larger fixed-Δt campaign is
needed. On level 1 the backward time-step study (Δt 1e-4 → 1.25e-5) gives
`u_r,max` changes of 1.42%, 0.78%, 0.07% (irregular; formal order 3.5, not
reported as meaningful).

Caveat: the backward scheme with a time step small relative to the mesh
develops a growing spurious mode. On level 1 at `Δt = 6.25e-6 s` it grows from
≈ 9 ms and diverges at 13.9 ms (every coupling iteration converged); at
`1.25e-5 s` it appears after 17 ms but stays bounded; implicit Euler shows
none. The scaled levels 2 and 3 show no trace of it to 20 ms, and the
principal QoIs (≤ 9 ms) are unaffected; the cause (the Robin wall conditions
with BDF2 are a candidate) was not investigated. A level-1 run at the level-3
time step is therefore not a valid temporal reference.

## 6. Time-discretisation explanation of the literature discrepancy

| Run (level 1 unless stated) | u_r,max (mm) | t(u_r,max) (ms) | u_z,min (mm) | u_r,min t>14 ms (mm) |
|---|---:|---:|---:|---:|
| Euler, Δt 1e-4 | 0.12233 | 7.005 | -0.07601 | -0.07756 |
| Euler, Δt 1e-4, level 2 | 0.12107 | 7.037 | -0.07515 | -0.07611 |
| Euler, Δt 5e-5 | 0.13638 | 7.078 | -0.08097 | -0.08573 |
| Euler, Δt 2.5e-5 | 0.14639 | 7.121 | -0.08430 | -0.09036 |
| Euler, Δt 1.25e-5 | 0.15259 | 7.155 | -0.08650 | -0.09367 |
| Backward, level 3 | 0.15847 | 7.254 | -0.08835 | -0.09566 |
| Lozovskiy et al. (2019), Euler 1e-4 | 0.1194 | 7.00 | -0.0740 | -0.0755 |
| Eken (2016), Euler 1e-4 | 0.1242 | 6.69 | -0.0794 | -0.0763 |
| Tuković et al. (2018), backward 2.5e-5 | 0.15694 | 7.23 | – | – |

Fraction of the published `u_r,max` gap reproduced by switching this code
from backward to Euler at `Δt = 1e-4 s` (damping / gap; range from the
±3e-6 m reading accuracy):

| Euler run | Reference | Gap to backward L3 | Gap to Tuković |
|---|---|---:|---:|
| level 1 | Lozovskiy | 92% (86-100%) | 96% (89-105%) |
| level 1 | Eken | 105% (97-116%) | 110% (101-122%) |
| level 2 | Lozovskiy | 96% (89-104%) | 100% (92-108%) |
| level 2 | Eken | 109% (100-120%) | 114% (105-126%) |

With the published discretisation this code lands within +2.5%/-1.5%
(level 1) and +1.4%/-2.5% (level 2) of the two published peaks, and within
1.7-4.3% of their axial and reflected troughs. The Euler damping is
essentially mesh-independent (1% between levels 1 and 2), and Euler converges
slowly in time (order 0.5-0.8; still 4.4% low at Δt = 1.25e-5 s), so neither
published peak is converged in time. Within the reading accuracy, the time
integration reproduces essentially all of the `u_r,max` difference and much
of the trough differences. It does not explain the 0.3 ms earlier peak of
Eken, and the references also differ in wall model (StVK vs Hookean, peak
hoop strain ≈ 3%), pressure switch-off and meshes, so the restrained claim
remains appropriate: *much of the discrepancy is reproduced by the time
integration*. An additional Euler time step was run (the whole Euler
sequence, 1e-4 to 1.25e-5 s); it strengthens the conclusion that the
published FE values are time-discretisation-limited.

## 7. Coupling/iterative error

- Robin-Neumann vs IQN-ILS (level 1, backward, Δt 2.5e-5, v2512): QoIs differ
  by ≤ 0.066%, radial histories by ≤ 0.26% of the peak; Robin 4.78 mean / 9
  max iterations, IQN-ILS 15.52 / 21. This is an order of magnitude below the
  level-2-to-3 changes. The same comparison on MeluXina gives 0.06% (QoIs) and
  0.26% (history), with 15.47 mean IQN-ILS iterations.
- Level 3 met all three Robin criteria at every step (worst final residuals
  4.8e-7, 1.0e-5, 5.4e-6 vs 1e-6, 1e-5, 1e-5), with no stalled steps; mean
  iterations fall with refinement (4.78, 4.19, 4.06; maximum 9 on every
  level): no pathological growth.
- Targeted tight-tolerance diagnostic (run because the L2→L3 change in
  `u_r,max` was unexpected): fluid relTol 1e-3 → 1e-6, pressure tolerance
  1e-9, PIMPLE control 1e-7, solid Newton 1e-9 with linear 1e-8, FSI 1e-8.
  Point-A QoIs change by ≤ 3e-6 relative on levels 1, 2 and 3 (level 3 to
  8.9 ms). Iterative error does not contaminate any principal QoI.

## 8. Remaining uncertainties

1. **Platform dependence of `u_z,min` (unresolved).** MeluXina (v2412
   EasyBuild) and the earlier Apple M1 runs give a 1.5% (L1) to 3.0% (L3)
   shallower axial trough than v2412 and v2512 on xenosim. `u_r,max`, `t_arr`,
   `c_p` agree to 0.1-0.2%. Each platform is invariant to ranks (1-8), solid
   preconditioner, tolerances and NaN-initialised memory; meshes are
   identical; the platforms differ already in the axial wall force of the
   first fluid solve (-2.60e-6 vs -6.53e-6 N, radial equal to six digits),
   with and without limiters or least-squares gradients. IQN-ILS
   (Dirichlet-Neumann, no Robin boundary conditions) shows the same split
   (-0.08950 vs -0.08815 mm), so the coupling is not the cause. Root cause unknown;
   the two platforms converge to different `u_z,min` limits, so this is not a
   discretisation error and does not vanish with the mesh.
2. **No asymptotic range for `u_r,max`, `t_arr`, `c_p`.** A fourth level
   (8.2 M fluid cells, about 16x level 3) would be needed to say more; the
   sub-nominal and non-monotone behaviour may reflect the near-singular wave
   front (3 ms pulse with a 0.1 ms ramp-off) and the coupled wall-fluid
   boundary layer.
3. **`c_p` sampling.** Cell-value probes on cell faces: fixed by an axial
   shift for the tie, but the remaining sampling noise (≈ Δz/c per station)
   is of the size of the L2→L3 change; interpolating probes would remove it.
4. **References.** All three are digitised; two are unconverged in time; one
   is same-lineage. Model differences (StVK vs Hookean, instantaneous vs
   ramped switch-off) are not quantified here.
5. **Small-Δt BDF2 mode** on coarse meshes (§5); not a problem on the scaled
   path but unexplained.

## 9. Consequences for the manuscript

Strengthened:

- The literature-discrepancy explanation: now quantified (92-114% of the
  `u_r,max` gap, troughs also reproduced), shown mesh-independent and backed
  by an Euler time-step sequence showing the published FE results are not
  converged in time.
- Coupling claims: Robin-Neumann vs IQN-ILS agreement is an order of
  magnitude below the discretisation changes; iterative error is shown
  negligible on all three levels.
- Temporal error: now measured on every level, not only on level 1.

Need weakening or correction:

- "`u_r,max` converged" (Lesson paragraph: "`u_r,max` converged, `u_z,min`
  not") must be reversed or removed: `u_r,max` changes by 0.81% at level 3 and
  is not in an asymptotic range; `u_z,min` is the only monotone QoI.
- "The working time step is converged to ≈ 0.1%" holds for the time step, but
  should not be read as overall convergence.
- "Lozovskiy2019 and Eken2016 (independent monolithic finite element codes)":
  Eken's fluid is finite volume; say "independent monolithic codes (finite
  element; finite volume-finite element)". Replace "finite element/finite
  volume difference" by "implicit-Euler/BDF2 difference".
- "largely a time-discretisation effect": supported for `u_r,max`; keep
  "much of the discrepancy is reproduced", and do not extend to arrival time.
- Any statement that the case is "verified" or "converged" in
  `sectionRecommendedSuite.tex` / `sectionBenchmarkSuite.tex`.

Tables and numbers to update:

- `tab:3dtube`: add level 3, replace levels 1-2 with the v2512 values
  (`u_z,min` -0.08954 / -0.08880; `c_p` 4.651 / 4.632 with shifted probes),
  add the L2→L3 change column and an order/verdict column.
- Text: "Two mesh levels are available ... a third (≈ 1.4e6 cells) has not
  been run" → three levels; the Tuković comparison (+1.0% at level 3, not
  +1.8%); `c_p` (4.59 m/s on level 3, 4.6% below 4.81 m/s).
- `sectionBenchmarkSuite.tex` summary table row: `3dTube ... 2` → 3 mesh
  levels; `sectionRecommendedSuite.tex` "needs level 3" → done.
- Methodology `\todo` on versions: 3dTube now recorded on OpenFOAM v2512
  (plus a v2412 replicate), not v2412 only.
- `benchmark_inventory.md`, `benchmarks.yaml`, `manuscript_todos.md` entries
  for 3dTube.

Reference category: stays **C-weak** (no change). Tuković remains a
same-lineage weak reference; Lozovskiy and Eken are independent but unconverged
in time.

Main text: **keep 3dTube in the main text**, re-framed: its value is (i) the
only 3-D transient FSI case, (ii) the quantified demonstration that a 20-25%
code-to-code discrepancy in the literature is a time-integration artefact,
and (iii) a cautionary example in which two levels suggested convergence that
a third did not confirm. It should not be presented as a converged
benchmark solution.

## 10. Recommended remaining work

Essential before publication:

- Resolve, or at least bound and disclose, the `u_z,min` platform dependence
  (start from the first fluid solve: compare the momentum matrix and wall
  boundary values of the two builds on the level-1 case; test a third
  platform/compiler). Until then report `u_z,min` with a ±3% platform
  uncertainty or drop it as a convergence QoI.
- Update `tab:3dtube` and the text from `results/3dTube_evidence.json`
  (generate the table from the CSV, manuscript todo F8).

High-value strengthening:

- **COMSOL (or other independent code) remains worthwhile**, now mainly for
  two reasons: an independent BDF2-type solution at converged Δt would give
  an external converged reference for `u_r,max` (where self-convergence is not
  demonstrated), and an independent `u_z,min` would arbitrate the platform
  dependence. It should report `u_r,max`, `u_z,min`, `t_arr` at A for two
  time integrators (BDF2 and implicit Euler at 1e-4 s) on at least two meshes.
- Interpolating pressure probes (`interpolationScheme cellPoint`) for `c_p`.
- A level 4 (≈ 8 M fluid cells; ≈ 16x level 3, ≈ 10-20 k core-hours) only if
  the paper wants a convergence statement for `u_r,max`; otherwise not needed.

Optional:

- Diagnose the small-Δt BDF2 mode (Robin boundary conditions).
- St Venant-Kirchhoff wall and instantaneous switch-off runs to quantify the
  remaining model differences with Lozovskiy/Eken.

## 11. Suggested manuscript wording

> Three mesh levels were computed (16 000 to 1 024 000 fluid cells, Δt halved
> with the mesh); the time-step error at the two finest levels is 0.1-0.2%,
> so the sequence is dominated by the mesh. The levels do not show
> convergence of the peak radial displacement: it changes by +0.07% and then
> by -0.81%, so the agreement of the first two levels was fortuitous, and its
> level-3 value (0.1585 mm, 1.0% above the earlier finite volume result of
> Tuković et al.) carries an estimated discretisation uncertainty of about 1%.
> The axial trough decreases monotonically (0.83%, 0.51%), with an observed
> order of about 0.7, well below the nominal order, and is sensitive at the
> 1.5-3% level to the software platform; the arrival time and pulse speed
> change by less than 1%. Iterative and coupling errors are below 0.1%.
>
> The published monolithic results (implicit Euler, Δt = 10⁻⁴ s) lie
> 21-25% below the BDF2 peak. Run with the same time integration, the
> present code reproduces the published peaks to within 2.5% and lowers its
> own peak by an amount equal, within the accuracy of the digitised figures,
> to the published difference; the axial and reflected troughs follow. Much
> of the apparent code-to-code discrepancy, and within the reading accuracy
> all of it for the peak radial displacement, is therefore a consequence of
> first-order time integration in the references, which are not converged in
> time; it is not a difference between finite element and finite volume
> discretisations.

## 12. Cross-platform u_z discrepancy investigation

The question: on two platforms the benchmark was deterministic but
`u_z,min(A)` differed by 1.5% (level 1) to 3% (level 3), while `u_r,max`,
arrival time and pulse speed agreed to 0.1-0.2%. Platforms: **xenosim**
(AMD EPYC 9684X; OpenFOAM v2512 and v2412 Ubuntu packages, GCC 11.4; the
solids4foam library compiled with the wmake default `-O3`; PETSc 3.24
development build) and **MeluXina** (AMD EPYC 7H12; OpenFOAM v2412 EasyBuild
foss-2024a, GCC 13.3, whose wmake rules compile with `-O2 -fno-tree-vectorize
-march=znver2`; PETSc 3.22). The investigation found **two independent
defects**, both in solids4foam, both now fixed. The full evidence, the
standalone reproducer and the diagnostics are on the solids4foam branch in
`tutorials/fluidSolidInteraction/3dTube/verification/platform/`.

### 12.1 Minimal reproducer

- **Defect 1** is visible on level 1, first time step, first FSI iteration,
  first PIMPLE corrector (`nOuterCorr 1`, `allowUnconvergedCoupling yes`,
  `nOuterCorrectors 1`): the axial component of the fluid force on the wall
  is -2.597e-6 N (xenosim) against -6.519e-6 N (MeluXina), the radial
  components equal to 16 digits.
- **Defect 2** is visible on level 2, first time step, second FSI
  iteration: the first solid solve is identical on both platforms (SNES
  norms to 1e-11), the second starts from a residual 0.6% different.
- Inputs: bit-identical case directories copied between the machines,
  meshes identical to 2e-17 m (points, owners, neighbours, patch ordering),
  dictionaries and fields identical. Case generation is excluded.

### 12.2 First-divergence evidence

| Stage (level, step, iteration) | Quantity | xenosim vs MeluXina |
|---|---|---|
| L1, plain OpenFOAM `pimpleFoam`, rigid wall, 1 step | all fields, wall force | equal to 1e-15 (axial viscous force 5.4626698912391e-06 N on both) |
| L1, solids4foam, step 1, it. 1, after the fluid solve | U, p, phi, Uf, meshPhi, motion fields (every cell and face) | max relative difference 2e-15 |
| same | wall `snGrad(U)`, cell `grad(U)`, wall `devReff` | identical (to 1e-16) |
| same | wall traction `rho*(nf & -devReff_b)`, axial | **differs: face 0 -3.9e-5 vs -0.659 Pa; total -2.597e-6 vs -6.519e-6 N** |
| L2, step 1, it. 1 (after fix 1) | fluid fields | 1e-8 relative (ordinary iterate-path round-off) |
| L2, step 1, it. 2 | solid interface displacement increment, input to the transfer | identical to 1e-10 |
| same | the same field after `transferPointsZoneToZone` (solid to fluid) | **differs 0.15-0.2% in the sums, asymmetric in x and y on both platforms** |
| same | fluid wall `pointMotionU` | differs 12-43% |

### 12.3 Wall-force decomposition

On the first fluid solve the axial wall force is entirely viscous: the
pressure contribution `sum(p Sf_z)` is 7e-19 N on both platforms. The
viscous traction is `rho*(nf & -devReff)` on the wall. Its inputs are
identical on both platforms: the wall normal (axial component < 3e-15), the
wall `snGrad(U)` (the wall cells are orthogonal, `|k| < 1e-19`, so the
non-orthogonal correction and the registered `grad(U)` that
`elasticWallVelocity::snGrad` looks up do not contribute), the boundary
gradient and the wall `devReff`. On xenosim, the same expression evaluated
on a *named* copy of the operands, or by a hand-written loop, gives
MeluXina's value; only the one-line expression on tmp operands is wrong. The
x and y components are right and the z component is wrong:

`z = (n&(-D))_x (-D_xz) + (n&(-D))_y (-D_yz)`, i.e. the z component is formed
from the already overwritten x and y of the result, which reproduces the
xenosim face values exactly.

### 12.4 Version/toolchain tests

| Test | Result | Rules in / out |
|---|---|---|
| OpenFOAM v2412 vs v2512 on xenosim (solids4foam, level 1) | identical | OpenFOAM release out |
| old solids4foam commit (`d35d59e11`) vs current | identical | solids4foam history out |
| MPI ranks 1, 2, 4, 8 (xenosim, level 1); serial on both | identical per platform | decomposition and reductions out (defect 1) |
| solid preconditioner hypre vs LU; tight tolerances; IQN-ILS vs Robin | identical per platform | linear solvers, iterative error, coupling scheme out |
| `FOAM_SETNAN`, `FOAM_SIGFPE` on both | unchanged | uninitialised memory out |
| convection limiter off; Gauss instead of least-squares gradients | split persists | operator/scheme choice out |
| plain `pimpleFoam` (no solids4foam code) | platforms agree to 1e-15 | OpenFOAM's own operators out |
| **standalone OpenFOAM program: `tmp<vectorField> & (-tmp<symmTensorField>)`** | wrong at `-O3` with GCC 11.4 and GCC 13.3 (error 7e-4 on values of 7e-4); exact at `-O2`, `-O1`, `-O0`, `-O3 -march=native` and `-O3 -D__restrict__=`; FMA contraction irrelevant | **defect 1 = aliasing undefined behaviour exposed by `-O3`** |
| solids4foam built `-O2` vs `-O3` on xenosim, level 2 | identical (after fix 1) | defect 2 not compiler-dependent |
| dump of the point-transfer weights, level 2 | all 6601 points have weights (1/3, 1/3, 1/3) | **defect 2 = degenerate barycentric fallback** |

### 12.5 Root cause

**Two root causes, both identified with reproducible discriminating tests
and both fixed.**

1. **Defect 1 (classification: compiler/toolchain sensitivity of an
   undefined-behaviour operator use, C/D).** solids4foam computed the fluid
   wall viscous force as `rho*(nf() & (-devReff().boundaryField()[patch]))`.
   OpenFOAM's field operators reuse the storage of the first tmp operand
   (the normal) for the vector result, while their loops access result and
   operands through `__restrict__` pointers; the result therefore aliases an
   input the compiler is told is not aliased. At `-O3` GCC computes and
   stores x and y, then forms z from the overwritten normal. Builds at `-O2`
   (the EasyBuild OpenFOAM rules) are correct. Effect: the axial wall shear
   was wrong (by up to 2.5x on the first step; 1.5-3% in `u_z,min`) on
   xenosim at every level. Fix: named normal field in every fluid model, and
   the same pattern removed from `newtonIcoFluid` and the `solidTractions`
   function object (solids4foam commit `635ff464f`). The hazard is in
   OpenFOAM's field algebra itself and should be reported to OpenCFD.
2. **Defect 2 (classification: FSI mapping defect, F; refinement-dependent).**
   The solid-to-fluid point transfer (`amiInterfaceToInterfaceMapping`, also
   `amiZoneInterpolation`) takes barycentric weights from
   `triangle::pointToBarycentric`, which returns the "degenerate" weights
   (1/3, 1/3, 1/3) when `d00*d11 - d01^2 < SMALL`. That quantity is
   `4 A^2` in m^4: about 1.5e-14 on level 1, 9.4e-16 on level 2 and 6e-17 on
   level 3, against `SMALL = 1e-15`. Every interface point of levels 2 and 3
   was therefore interpolated as a three-point average, level 1 correctly.
   The smeared transfer broke the x-y symmetry of the quarter tube (solid
   wall force x and y 0.3% apart) and, because the average depends on which
   triangle the addressing search picks, results depended on decomposition
   (xenosim level 2: 0.16000 mm on 8 ranks, 0.16029 on 32) and build. Fix: a
   scale-invariant weight computation (`triangleWeights.H`, solids4foam
   commit `0ff7f462b`).

**Defect 2 also contaminated the original mesh study on every platform**:
level 1 and levels 2-3 were coupled through different interpolation
operators, so the earlier non-monotone `u_r,max` sequence (Section 4) was not
evidence about the discretisation.

### 12.6 Quantified platform uncertainty and corrected results

After both fixes:

| | u_r,max (mm) | u_z,min (mm) | t_arr (ms) | c_p (m/s) |
|---|---:|---:|---:|---:|
| L1 xenosim (4 ranks) / MeluXina (8) | 0.15978 / 0.15978 | -0.08819 / -0.08819 | 5.873 / 5.873 | 4.654 / 4.654 |
| L2 xenosim (8) / MeluXina (32) | 0.15946 / 0.15946 | -0.08659 / -0.08659 | 5.843 / 5.843 | 4.635 / 4.634 |

**Platform uncertainty in `u_z,min(A)`: below 0.01% (no difference at the
resolution of the reported digits on levels 1 and 2; on the level-1 and
level-2 first-step reproducers the wall forces agree to 12 significant
digits).** It is two orders of magnitude below
the level-2-to-3 change. Level 3 was computed once (MeluXina, 128 ranks); a
second level-3 run on xenosim was not possible within the queue (estimated
start a week later), so the level-3 platform agreement is inferred from
levels 1 and 2 and from the now rank-independent formulation.

Corrected three-level study (both fixes; L1-L2 xenosim, L3 MeluXina; Δt
halved with the mesh; `c_p` with probes off the cell faces):

| Quantity | L1 | L2 | L3 | Δ1→2 | Δ2→3 | Order | Verdict |
|---|---:|---:|---:|---:|---:|---:|---|
| u_r,max(A) (mm) | 0.15978 | 0.15946 | 0.15851 | -0.20% | -0.60% | undefined | monotone, change grows: **not asymptotic** |
| t(u_r,max) (ms) | 7.201 | 7.215 | 7.252 | +0.19% | +0.50% | undefined | flat peak, extraction-limited |
| u_z,min(A) (mm) | -0.08819 | -0.08659 | -0.08571 | 1.87% | 1.02% | 0.87 | monotone, ~first order (nominal 2) |
| t_arr(A) (ms) | 5.873 | 5.843 | 5.825 | -0.51% | -0.31% | 0.72 | monotone, sub-nominal |
| c_p (m/s) | 4.654 | 4.635 | 4.594 | -0.41% | -0.88% | undefined | change grows; at the probe sampling noise |
| u_r,min(A), t>14 ms (mm) | -0.10014 | -0.09409 | -0.09603 | 6.3% | 2.0% | undefined | non-monotone |
| mean FSI iterations | 4.78 | 3.87 | 3.18 | | | | all steps converged, max 9 |

Fixed `Δt = 2.5e-5 s` (pure mesh refinement): `u_r,max` 0.15978, 0.15921,
0.15835 (0.36%, 0.54%; undefined); `u_z,min` -0.08819, -0.08643, -0.08554
(2.0%, 1.03%; order 0.98); `t_arr` 5.873, 5.837, 5.815 (0.61%, 0.38%; order
0.65). Temporal error at the fine levels: 2.5e-5 → 1.25e-5 s on L2 changes
`u_r,max` +0.16%, `u_z,min` +0.18%, `t_arr` +0.10%; 2.5e-5 → 6.25e-6 s on L3
+0.10%, +0.20%, +0.17%. The scaled path is mesh-dominated (0.4-1.0% mesh
against 0.1-0.2% time step).

Literature (unchanged in substance): level-3 `u_r,max` is 1.0% above
Tuković et al. (2018); implicit Euler at `Δt = 1e-4 s` lands within
+2.5%/-1.4% (L1) and +1.3%/-2.6% (L2) of Lozovskiy/Eken and its damping
(0.036-0.038 mm) is 92-115% of the published gap (86-126% within the reading
accuracy); Robin-Neumann and IQN-ILS agree to 0.07% (history 0.26%).

### 12.7 Consequences for the manuscript

- **`u_z,min` can remain a QoI.** The platform dependence was not a property
  of the benchmark or of the axial displacement; it was two code defects,
  and after the fixes the platform uncertainty is below 0.01%. The axial
  displacement is not ill-conditioned (classification G is excluded: the
  corrected platforms agree to the reported digits, and the sequence is
  monotone).
- **The convergence statements change.** Report: `u_z,min` and `t_arr`
  converge monotonically at about first order (0.87-0.98 and 0.65-0.72);
  `u_r,max` decreases monotonically but its changes are not yet decreasing
  (0.20% → 0.60%), so no order and no asymptotic claim; its level-3
  uncertainty is about 0.6% or more. `u_z,min` is now the *better-behaved*
  displacement QoI; the earlier statement "`u_r,max` converged, `u_z,min` not"
  is wrong in both halves and must not be used.
- **Every 3dTube number recorded before commit `0ff7f462b` must be
  replaced**, including the two-level `tab:3dtube` values in the current
  manuscript (recorded on an Apple M1 with the defective point transfer on
  level 2): use the table in 12.6. The FE/FV time-integration conclusion and
  its numbers are confirmed (Euler L1 `u_r,max` 0.12241 mm, against 0.12233 mm
  on xenosim before the fixes).
- **`u_r,max` and the time-integration conclusion are not affected by the
  platform question** (u_r,max agreed to 0.1% even before the fixes; the
  Euler damping is 0.036-0.038 mm in every variant).
- **A methodological point worth a sentence in the paper:** the verification
  campaign found two solver defects that no single-platform run would have
  shown: a compiler-dependent miscompilation and a refinement-dependent
  mapping degeneracy that switched on between levels 1 and 2. Cross-platform
  replication and symmetry checks are cheap and were decisive.
- **COMSOL** remains useful for `u_r,max` (not self-converged) but is no
  longer needed to arbitrate `u_z,min`.

### 12.8 Remaining work

- Not needed for `u_z,min`: the cause is identified, fixed and verified.
- Optional: a second corrected level-3 run on another platform (the queue
  did not allow it), and an upstream report to OpenCFD about tmp reuse with
  `__restrict__` in compound inner products.
- Small-time-step behaviour: PENDING_TS5_PAPER
- The earlier sections' general recommendations still hold for `u_r,max`:
  a fourth level or an independent converged solution would be needed for a
  convergence statement on the peak.
