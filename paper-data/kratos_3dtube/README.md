# 3dTube with Kratos Multiphysics (steps 1 and 2)

An independent-code solution of the 3dTube pressure pulse for `sec:3dtube`. Kratos shares no code with
solids4foam. This folder holds the input generator, the run scripts, every run's history at the nine wall
stations with its coupling record, the QoI tables and the figures. No `.tex` file is changed. The discussion is in
issue #36.

**Main results.**

- **Step 1.** With implicit Euler at Δt = 1e-4 s (level 1), Kratos reproduces the published finite element
  values. Its peak u_r,max(A) is 0.12334 mm, 22.2% below solids4foam's backward level 3, as solids4foam's
  own Euler run (22.8%), Lozovskiy et al. (24.7%) and Eken (21.6%) are. The 20–25% gap between the
  published values and solids4foam is first-order time integration, not a code difference.
- **Step 2.** With BDF2 and Δt halved with the mesh (levels 1 and 2), the Kratos peak rises by 30%, as in
  solids4foam (30.5%). At level 2 Kratos is within 0.2% of solids4foam's level 3 in u_z,min(A), c_p and
  the peak time, and 1.5% above it in u_r,max(A).
- **Arrival time.** t_arr(A) is 1.4% later in Kratos at matched levels. The definition is the same; the
  difference is a wall front that travels 1.4–1.5% slower between z = 5 and 25 mm at both levels.

## Code and platform

- Kratos Multiphysics 10.3.0 from PyPI (KratosMultiphysics, KratosFSIApplication,
  KratosFluidDynamicsApplication, KratosStructuralMechanicsApplication, KratosMeshMovingApplication,
  KratosMappingApplication and KratosLinearSolversApplication, all 10.3.0) with CPython 3.11.18, as for
  `../kratos_mok_lid_driven_cavity/`.
- The same UCD workstation, an Intel i9-10900X (10 cores) with 62 GB and Ubuntu 18.04.6, running OpenMP
  with 2 to 10 threads. Every run's `run_info.json` records the version string, host, CPU and thread count, and
  `time.txt` the wall time and peak memory.
- The QoI rules (`peak`, `first_peak`, `arrival_time`, `front_arrival`, `fit_speed`), the stations and the
  digitised published histories are those of the solids4foam driver
  `tutorials/fluidSolidInteraction/3dTube/verification/scripts/tube_3d_verification.py` and its
  `reference/` directory, solids4foam commit `a421f063d`. The scripts import them and do not copy them.

## Problem and discretisation

The case is that of the solids4foam tutorial `tutorials/fluidSolidInteraction/3dTube`. The tube is 5 cm
long, with inner radius 0.5 cm and wall thickness 0.1 cm, and is clamped at both ends. The fluid has
ρ_f = 1000 kg/m³ and ν_f = 3e-6 m²/s, and the wall ρ_s = 1200 kg/m³, E = 3e5 Pa and ν_s = 0.3. The inlet
pressure is 1333.2 Pa to 3 ms, falling linearly to 0 at 3.1 ms, and the outlet pressure is 0. A quarter of
the tube is modelled. Point A is on the inner wall at mid-length, (0, 0.005, 0.025) m.

- **Mesh.** The solids4foam block layout, with a fluid O-grid (core 10 × 10, outer blocks 5 radial × 20
  around, 80 along) and a wall of 4 × 20 × 80, all times the level factor. Each fluid hexahedron is split into 6
  tetrahedra (two prisms of three tetrahedra, minimum-index rule, conforming). Level 1 has 96 000
  tetrahedra and 18 306 nodes (16 000 hexahedra in solids4foam) and 6 400 wall hexahedra; level 2 has
  768 000 tetrahedra, 137 011 nodes and 51 200 wall hexahedra. Fluid and wall nodes coincide on the
  interface. The circular arcs use exact coordinates; the solids4foam master tutorial rounds them to four
  digits, which moves the wall by about 1e-4 of the radius.
- **Fluid.** QSVMS P1/P1 tetrahedra with `dynamic_tau = 1`, Newton tolerance 1e-6 (relative), AMGCL 1e-8.
  Kratos has no 3D quadrilateral wall condition for hexahedral fluid elements, and `TimeIntegratedQSVMS`
  cannot compute the right-hand side alone, which the FSI reactions need.
- **Wall.** `SmallDisplacementElement3D8N` with `LinearElastic3DLaw` (small strain, as in the tutorial).
- **Mesh motion.** Structural similarity, fixed at both ends.
- **Boundary conditions.** Inlet pressure applied as a fixed PRESSURE and as EXTERNAL_PRESSURE through
  `ApplyOutletProcess`; outlet p = 0; inlet and outlet velocity natural; normal velocity, mesh and wall
  displacement fixed on the symmetry planes; wall clamped at both ends.
- **Coupling.** Dirichlet–Neumann with MVQN. A step is converged when the RMS interface displacement
  residual |r|/√n_DOF is below 1e-9 m (about 1e-5 of the peak), with at most 50 iterations. Kratos accepts
  a step at the limit with only a warning, so `MainKratos_tube.py` records for each step whether it met
  the tolerance and its residual at acceptance.
- **Time integration.** Step 1 uses implicit Euler. Kratos's QSVMS scheme is BDF2-based, so
  `MainKratos_tube.py` overwrites the scheme's three BDF coefficients with [1/Δt, −1/Δt, 0] after each
  step's initialisation, replaces the mesh solver's velocity helper by `TimeDiscretization.BDF1()`, and the
  wall uses `backward_euler`. Step 2 uses BDF2 in the fluid, the mesh velocity and the wall.
- **Histories.** At every step, `history.csv` records the wall displacement (u_r = u_y and u_z) of the
  wall node at (0, 0.005, z) and the fluid pressure at (0.0005, 0.0005, z), interpolated in the containing
  tetrahedron, at the nine stations z = 5, 10, …, 45 mm, with the coupling solves, the converged flag and
  the residual at acceptance.

## Checks

- **Wall.** Static internal pressure 1333.2 Pa, clamped and symmetric (`scripts/wall_static_check.py`).
  The mid-length inner u_r is 0.12131 mm (level 1) and 0.12138 mm (factor 2), against 0.1208 mm from Lamé
  with zero axial strain (+0.4%, the finite length and the clamps). No locking.
- **Openings and implicit Euler.** A rigid tube under a constant 1333.2 Pa (`scripts/fluid_rigid_check.py`).
  The axis velocity at mid-length follows Δp/(ρL) t to 0.04% after one step and below 0.01% at 1 ms with
  the BDF1 override. Inlet PRESSURE and EXTERNAL_PRESSURE are 1333.2 Pa on every node.
- **BDF2 start.** Kratos's BDF2 starts from a two-step rest history, which delays a step load by Δt/2 for
  the whole run (rigid tube, −33% after one step, −5% at 1 ms). BDF2 runs therefore take a BDF1 first step
  in the fluid and the mesh velocity, as OpenFOAM's `backward` scheme does. The wall keeps Kratos's start,
  because the structural scheme holds its own copy of the coefficients. Its effect should be small, since
  the wall inertia is about a quarter of the added mass (ρ_s h / (ρ_f R) = 0.24).
- **Fluid stabilisation.** With `dynamic_tau = 0`, the QSVMS default, τ is about 1 ms here (ten Euler
  steps, comparable to the pulse) and the level-1 Euler peak came out 7% low (0.1145 mm,
  `runs/f1_dt1e-4_euler_nobfs`). That run then diverged at 10.3 ms, when backflow started at the inlet;
  with Kratos's backflow stabilisation from t = 0 it ran to the end, with a peak of 0.1116 mm
  (`runs/f1_dt1e-4_euler_bfs`). With `dynamic_tau = 1` (τ about Δt) the coupling does not diverge and no
  backflow stabilisation is needed, and all reported runs use it. Tables of the two superseded runs are in
  `results/step1_euler_dtau0_superseded/`.
- **Wall linear solver.** Kratos's default for the wall is the single-threaded direct `sparse_lu`, which
  at factor 2 (about 178 000 DOFs) takes minutes per coupling iteration. The step-2 runs use CG with AMG
  (`--structure-solver amgcl`, tolerance 1e-12). Over 10 level-1 steps it agrees with the direct solve to
  5e-15 m with the same coupling iterations (`runs/_debug/f1_rec_amg_10steps` against
  `runs/_superseded/f1_dt2.5e-5_bdf2_rec`).
- **Coupling accelerator.** The level-1 production run (full MVQN, AMG wall) agrees with the full MVQN and
  direct wall solve to 4e-9 m over the first 50 steps (`runs/_crashed/f1_dt2.5e-5_bdf2`, which stopped at
  step 51 during a period of faulty memory on the workstation), and `MVQN_recursive` agrees with it to
  5e-9 m (`runs/_superseded/f1_dt2.5e-5_bdf2_rec`), both at the 1e-9 m RMS tolerance. `MVQN_recursive`
  needs about 11.5 coupling solves per step against 7.1, so the full MVQN is used.
  `MVQN_randomized_SVD` crashes reproducibly at step 3 in Kratos 10.3.0 (null pointer in
  `KratosFSIApplication.so`) and is not used.
- **Axis pressure.** The pressure differs by a few Pa between accelerators at the same tolerance, because
  it follows the wall acceleration. It enters only c_p, through the arrival of half the inlet pressure.

## Step 1, implicit Euler, level 1, Δt = 1e-4 s

`runs/f1_dt1e-4_euler_dtau1_nobfs`: 201 steps, 56 min on 4 threads, 1.4 GB, 4.5 coupling solves per step
(max 13), every step met the tolerance. Kratos's backflow stabilisation at the openings is not needed with
`dynamic_tau = 1`; a run with it from t = 0 (`runs/f1_dt1e-4_euler_dtau1_bfs`) differs by 2.4–2.5% in
the three displacement extrema.

| QoI | Kratos | Lozovskiy (2019) | Eken (2016) | solids4foam Euler L1 |
|---|---|---|---|---|
| u_r,max(A) (mm) | 0.12334 | 0.1194 | 0.1242 | 0.12241 |
| t(u_r,max) (ms) | 7.12 | 7.00 | 6.69 | 7.005 |
| u_z,min(A) (mm) | −0.07154 | −0.0740 | −0.0794 | −0.07601 (before the fixes of solids4foam PR #546) |
| u_r,min(A), t > 14 ms (mm) | −0.07689 | −0.0755 | −0.0763 | −0.07756 (before the fixes) |

The full table is in `results/step1_euler/summary.md`, and the histories against the published curves
are in `figures/step1_point_A_histories.png`.

## Step 2, BDF2, levels 1 and 2, Δt halved with the mesh

| Run | Level | Δt (s) | Steps | Wall time | Threads | Memory | Coupling solves mean / max |
|---|---|---|---|---|---|---|---|
| `f1_dt2.5e-5_bdf2_amg` | 1 | 2.5e-5 | 800 | 2.6 h | 6 | 1.2 GB | 4.62 / 14 |
| `f2_dt1.25e-5_bdf2_amg` | 2 | 1.25e-5 | 1601 (to 20.0125 ms) | 33.2 h | 10 | 12.4 GB | 5.31 / 50 |

Every step met the tolerance except the first step of level 2 (below).

| QoI | Kratos L1 | solids4foam L1 | Kratos L2 | solids4foam L2 | solids4foam L3 |
|---|---|---|---|---|---|
| u_r,max(A) (mm) | 0.16035 | 0.15978 | 0.16096 | 0.15946 | 0.15851 |
| t(u_r,max) (ms) | 7.175 | 7.201 | 7.262 | 7.215 | 7.252 |
| u_z,min(A) (mm) | −0.08198 | −0.08819 | −0.08555 | −0.08659 | −0.08571 |
| t_arr(A) (ms) | 5.954 | 5.873 | 5.926 | 5.843 | 5.825 |
| c_p (m/s) | 4.578 | 4.654 | 4.583 | 4.635 | 4.594 |

The solids4foam values are the backward mesh study (`../3dtube_level3_findings.md`, Section 12.6, and
`../3dtube_level3_ecore/`).

- **u_z,min.** The 7% gap at level 1 is a mesh effect. The two codes approach from opposite sides, and
  Kratos level 2 is within 0.2% of solids4foam's level 3.
- **u_r,max.** Kratos is 0.4% above solids4foam at level 1 and 0.9% at level 2, and 1.5% above its level 3.
  The Kratos peak rises with refinement (+0.4%) while solids4foam's falls (−0.2%, then −0.6%), so neither
  sequence is asymptotic in the peak.
- **Timing.** The peak time agrees to 0.1% at the finest levels, and c_p to 0.2% against solids4foam's
  level 3. t_arr is discussed below.
- **Time integration.** The BDF2 peak is 30.0% above the Kratos Euler peak, as in solids4foam (30.5%), so
  step 2 confirms the step-1 conclusion with second-order time integration.

The tables are in `results/step2_bdf2/summary_level1.md` and `summary_levels12.md`, and the histories
in `figures/step2_point_A.png`.

### First step of level 2

The first step stalled in MVQN. The residual fell to a floor of 1.66e-9 m at iteration 15, drifted up,
and the step was accepted at 4.3e-8 m after 50 iterations. `results/step2_bdf2/level2_step1_coupling_residuals.csv`
has the step-1 residuals of four runs with the same mesh and time step.

- With the fluid Newton tolerance at 1e-8 and AMGCL at 1e-10 (`runs/_debug/f2_amg_tight_4steps`, stopped
  in step 1, residuals in its `log.txt`), the residuals follow the production run to a relative 3e-7 up
  to the floor and to 6e-4 at iteration 21. A run with the production tolerances on 5 threads
  (`runs/_debug/f2_full_amg_3steps`) differs from production by more than that, and converged in 28
  iterations. So the inner tolerances are not the cause. The stall is in the MVQN update, and thread-count
  round-off decides whether the step slips under the tolerance.
- Against that converged run, the production displacements differ by up to 1.4e-7 m (8e-4 of the level-2
  peak) at the end of step 1 and by at most 1.8e-8 m (1.1e-4 of the peak) after three steps.
  `runs/_debug/f2_rec_amg_3steps` (`MVQN_recursive`) agrees with the converged full-MVQN run to 4e-9 m.
- The production run is kept, with its first step reported as over the tolerance.

## Arrival time

`scripts/arrival_time.py` writes `results/arrival_time/arrival_time.md`.

- `tube_qoi.py` takes t_arr from the driver's `front_arrival`, which is the first upward crossing of 50% of the
  incident peak, the incident peak being the first local maximum above 50% of the history maximum.
- To compare full histories, the solids4foam level 1 was rerun to 10 ms
  (`./Allverify --quick --cores 1`, OpenFOAM v2606, serial, 590 s). It reproduces the published level 1
  (t_arr 5.872 against 5.873 ms, identical wall arrivals at every station).
  `results/arrival_time/solids4foam_level1_rerun/` has its monitors, settings and driver output.
- The offset is not a fixed time shift. Kratos minus solids4foam wall arrival grows with distance from the
  inlet, from +0.013 ms at z = 5 mm to +0.081 ms at A (level 1), and from +0.012 to +0.083 ms (level 2).
  The least-squares time shift that aligns the incident fronts at level 1 grows the same way, from
  −0.025 ms at 5 mm to +0.068 ms at A, so the front itself is late, not only its 50% crossing.
- Between 5 and 25 mm the wall front travels at 4.203 and 4.238 m/s in Kratos (levels 1 and 2) and at
  4.263, 4.302 and 4.326 m/s in solids4foam (levels 1 to 3). Both codes speed up by about 0.8% per level,
  and Kratos stays 1.4–1.5% slower at matched levels. The same gap shows in c_p at matched levels (−1.6%,
  −1.1%).
- The inputs are the same in both codes (small-strain linear elastic wall, material data, inlet pulse),
  and solids4foam's Robin–Neumann and IQN-ILS couplings agree to 0.07%. What remains is the
  discretisation (hexahedral finite elements and P1 QSVMS against cell-centred finite volumes).
- Tuković et al. (2018), with the same rule, give t_arr(A) = 5.889 ms, between solids4foam level 3
  (5.825 ms) and Kratos level 2 (5.926 ms).
- The solids4foam axis-pressure probes take the value of the containing cell. At level 1 the per-station
  pressure arrivals of this rerun (unshifted probes) are about 0.13 ms (one cell) earlier than those of the
  published run with `--probe-z-shift 7.8125e-5`, while c_p differs by 0.2%.

## Runs

| Run | Step | Settings | Status |
|---|---|---|---|
| `f1_dt1e-4_euler_dtau1_nobfs` | 1 | Euler, level 1, Δt 1e-4 | reported |
| `f1_dt1e-4_euler_dtau1_bfs` | 1 | the same with backflow stabilisation from t = 0 | sensitivity |
| `f1_dt1e-4_euler_bfs` | 1 | `dynamic_tau = 0`, backflow stabilisation from t = 0 | superseded |
| `f1_dt1e-4_euler_nobfs` | 1 | `dynamic_tau = 0`, no backflow stabilisation | superseded (peak 7% low), diverged at 10.3 ms (history valid to 10.2 ms, `STOPPED.txt`) |
| `f1_dt2.5e-5_bdf2_amg` | 2 | BDF2, level 1, Δt 2.5e-5 | reported |
| `f2_dt1.25e-5_bdf2_amg` | 2 | BDF2, level 2, Δt 1.25e-5 | reported |
| `_crashed/f1_dt2.5e-5_bdf2` | 2 | level 1, full MVQN, direct wall solve, 50 steps | check |
| `_superseded/f1_dt2.5e-5_bdf2_rec` | 2 | level 1, `MVQN_recursive`, direct wall solve, 63 steps | check |
| `_debug/f1_rec_amg_10steps` | 2 | level 1, `MVQN_recursive`, AMG wall solve, 10 steps | check |
| `_debug/f2_full_amg_3steps` | 2 | level 2, production inputs, 5 threads, 3 steps | check |
| `_debug/f2_rec_amg_3steps` | 2 | level 2, `MVQN_recursive`, 3 steps | check |
| `_debug/f2_amg_tight_4steps` | 2 | level 2, tight inner tolerances, stopped in step 1 | check (`log.txt`) |

Attempts that ended for other reasons are not included. These are runs with the backflow stabilisation
switched on during the run, which made MVQN diverge, `MVQN_randomized_SVD` runs that crashed, and
level-2 runs with the direct wall solve that did not finish step 1. Their notes stay on the workstation
and are available on request.

## Files

| Path | Contents |
|---|---|
| `scripts/tube_case.py` | Writes an input set (for example `tube_case.py <out> --factor 1 --dt 2.5e-5 --scheme bdf2 --structure-solver amgcl`) and records its options in `case.json` |
| `scripts/MainKratos_tube.py` | Runs a case with the inlet pulse, the implicit-Euler override and the BDF1 first step, and writes the station histories, the coupling record and `run_info.json` |
| `scripts/run_case.sh` | `MAIN=MainKratos_tube.py scripts/run_case.sh <input dir> <run dir> <threads>` |
| `scripts/tube_qoi.py` | The QoI tables in `results/step*/` |
| `scripts/arrival_time.py` | The arrival-time comparison in `results/arrival_time/` |
| `scripts/plot_tube_history.py`, `scripts/plot_tube_qoi.py`, `scripts/plot_step2.py` | The figures (need matplotlib) |
| `scripts/wall_static_check.py`, `scripts/fluid_rigid_check.py` | The wall and rigid-tube checks |
| `runs/<name>/` | `history.csv`, `case.json`, the `ProjectParameters.json` and materials files that ran, the `MainKratos_tube.py` that ran, `inputs.sha256` (including the meshes), `run_info.json`, `run_meta.txt`, `time.txt`; `log.txt` for the check runs |
| `results/step1_euler/`, `results/step1_euler_dtau0_superseded/`, `results/step2_bdf2/` | `summary*.md` and the matching `*_histories.csv`; the level-2 step-1 residuals |
| `results/arrival_time/` | `arrival_time.md`, `pointA_level1.csv` (point-A histories of both codes at level 1), and the solids4foam rerun |
| `figures/` | `step1_point_A_histories.png`, `step1_qoi_comparison.png`, `step2_point_A.png` |
| `SHA256SUMS` | Checksums of every file in this folder |

## Reproducing

From this folder, with Kratos 10.3.0 installed (`KRATOS_ENV` may name a script that activates it):

```bash
python scripts/tube_case.py gen/f1_dt2.5e-5_bdf2_amg --factor 1 --dt 2.5e-5 --scheme bdf2 --structure-solver amgcl
MAIN=MainKratos_tube.py scripts/run_case.sh gen/f1_dt2.5e-5_bdf2_amg myruns/f1_dt2.5e-5_bdf2_amg 6
```

The current `tube_case.py`, given the options in each run's `case.json`, regenerates the meshes and input
files of all six main runs with the same sha256 as in their `inputs.sha256`. Older runs used earlier
versions of `MainKratos_tube.py`; each run folder has the one that ran.

The tables need the solids4foam driver (commit `a421f063d`, or the current master, whose QoI functions
are the same):

```bash
D=<solids4foam>/tutorials/fluidSolidInteraction/3dTube/verification/scripts/tube_3d_verification.py
python scripts/tube_qoi.py $D runs/f2_dt1.25e-5_bdf2_amg runs/f1_dt2.5e-5_bdf2_amg --out results/step2_bdf2/summary_levels12.md
python scripts/arrival_time.py $D ../3dtube_level3_ecore/mesh_robin_backward_pz7.8125e-05.csv \
    results/arrival_time/solids4foam_level1_rerun runs/f1_dt2.5e-5_bdf2_amg runs/f2_dt1.25e-5_bdf2_amg \
    --out results/arrival_time/arrival_time.md --csv results/arrival_time/pointA_level1.csv
```

Rerunning `tube_qoi.py` from this layout reproduces all four `summary*.md` files and their history CSVs
exactly.

## Limitations

- Two levels only. A third level (factor 3, ratio 1.5) would take one to two weeks on the workstation and
  would need `MVQN_recursive` to fit in memory; the factor-4 level needs an MPI build. Philip judged two
  levels enough for the paper (#36), because the remaining u_r,max difference is in a quantity the
  solids4foam sequence has not converged in either.
- The fluid stabilisation (`dynamic_tau = 1`) and the tetrahedral fluid mesh are numerical choices that
  differ from solids4foam. The cause of the 1.4% slower front (wall or fluid discretisation) is not
  separated; a Kratos run with only the wall or only the fluid refined would do it.
- The wall keeps Kratos's BDF2 start from a two-step rest history (see Checks).
- The first step of level 2 is over the coupling tolerance (see above).
- The VTK fields (written every 1 ms) and the production `log.txt` files are not included; they are
  available on request.
