# cavityFlexibleBottom with preCICE: SU2 fluid + CalculiX solid, steady state (November 2026)

This is an independent-code solution of the steady cavityFlexibleBottom problem (issue #36). The fluid code (SU2, a vertex-based finite volume solver), the solid code (CalculiX, finite elements) and the coupling library (preCICE) are all independent of solids4foam. The platform, build and coupling stack are the same as in `../fsi1_su2_calculix`.

## Blind provenance

These results were computed **blind** and posted to #36 before any comparison (comment of 16 November 2026, https://github.com/philipcardiff/paper-fsi-verification/issues/36#issuecomment-6928621052).
- The case definition comes only from the solids4foam tutorial `tutorials/fluidSolidInteraction/cavityFlexibleBottom` at `development` commit `1b5ab09` (`1b5ab0963c847cd9a6d64937f0397434d7b054ef`): the description and geometry in its README, and its case files (`system/`, `constant/`, `0/`).
- Not consulted before the numbers were posted: the tutorial's `verification/` directory, its result plots and gnuplot scripts, the README's results/literature/verification sections, the paper's cavity section, the cavity folders in `paper-data/`, `simulationData/`, and the ANSYS/COMSOL results.
- The references were opened only after the solids4foam sequence was shared on #36 (19 November 2026, https://github.com/philipcardiff/paper-fsi-verification/issues/36#issuecomment-7015729752). No run, setting or value in this package was changed afterwards; the comparison section and `scripts/compare_references.py` are the only additions.

## Platform

As `../fsi1_su2_calculix` (README, section Platform): AMD Ryzen Threadripper PRO 9985WX workstation, Ubuntu 24.04.5, GCC 13.3.0, Open MPI 4.1.6; SU2 7.5.1 (`INC_NAVIER_STOKES`) with the preCICE SU2 adapter (`fde3106`); CalculiX 2.20 (`ccx_preCICE`, adapter 2.20.2); preCICE 3.5.0, pyprecice 3.5.0. Environment: `scripts/env-su2-precice.sh` (the FSI1 file). Jobs ran through Slurm with dedicated cores.

## Case (from the solids4foam tutorial)

- Channel x = 0–14 m, y = 1–2 m, over a cavity x = 2–6 m down to y = 0 (rounded corners, as the tutorial blockMesh). The cavity floor is a flexible plate x = 2–6 m, y = −0.4–0 m, clamped at both ends, bottom traction-free.
- Fluid: ρ = 1 kg/m³, ν = 0.01 m²/s. Parabolic inlet with a 1.5 m/s peak (mean 1 m/s), outlet p = 0, no-slip walls.
- Solid: St Venant–Kirchhoff, E = 500 Pa, ν = 0.4, plane strain.
- Point A = (4, 0): the mid-span of the plate's top (wetted) surface.

## Set-up

- **Fluid mesh:** the tutorial's fluid blockMeshDict (`case/fluid-blockMeshDict`, copied unchanged from `system/fluid/blockMeshDict` at `1b5ab09`), with every in-plane block division multiplied by the level L, meshed by OpenFOAM v2412 `blockMesh` and converted to 2-D `.su2` (`scripts/make_fluid_mesh.py`, the FSI1 converter). The tutorial mesh is level 1 (dx = 0.1 m); levels L = 1, 2, 4, 8, 16 halve dx each time.
- **Solid mesh:** (40 L) × (4 L) C3D8I elements, one layer of 1 m depth with u_z = 0 (plane strain), refined together with the fluid (`scripts/solid_mesh.py`). Fluid and solid interface nodes coincide at every level (|Δx| < 2e-15 m).

| Level | Fluid cells | Fluid nodes | Plate edges | Solid elements (x × y) |
|---|---:|---:|---:|---:|
| 1 | 1 688 | 1 847 | 40 | 40 × 4 |
| 2 | 6 752 | 7 069 | 80 | 80 × 8 |
| 4 | 27 008 | 27 641 | 160 | 160 × 16 |
| 8 | 108 032 | 109 297 | 320 | 320 × 32 |
| 16 | 432 128 | 434 657 | 640 | 640 × 64 |

- **Fluid settings** (`case/fluid-su2/*.cfg`):
  - Dimensional incompressible solver, constant density, ρ = 1, μ = 0.01 Pa s. Inlet profile at the inlet nodes (`scripts/make_inlet_profile.py`).
  - FDS with MUSCL and no limiter, weighted-least-squares gradients, default `BETA_FACTOR`.
  - BDF2 dual time stepping; implicit Euler in pseudo-time, CFL 1000, FGMRES + ILU (1e-8). Inner iterations stop when the RMS residuals of pressure and both velocity components are below 1e-10 (absolute), at most 1 000.
  - Mesh motion: linear elasticity with wall-distance stiffness; only the plate moves.
  - Forces: SU2's force coefficients with force reference ½ρU²A = 1, so CD/CL are F_x/F_y on the plate in N per metre depth (gauge pressure plus viscous stress).
- **Solid** (`case/solid-calculix/plate.inp`): CalculiX `*STATIC` with `NLGEOM` and `*ELASTIC` (exactly St Venant–Kirchhoff). The adapter solves each coupling window as one full-load static increment from the previous converged state. The preCICE load is applied with a constant amplitude, because a `*STATIC` step otherwise ramps loads from zero within the increment. The slab is 1 m deep because SU2 forces are per metre depth.
- **Coupling** (`case/precice-config.base.xml`): parallel-implicit with IQN-ILS (residual-sum preconditioner, QR2 filter 1e-3, 10 reused windows, initial relaxation 0.1), relative convergence limit 1e-6 on force and displacement, at most 50 iterations: as `../fsi1_su2_calculix`. One difference: the mapping is nearest-neighbour both ways (conservative for force, consistent for displacement), which is exact because the interface nodes coincide. (The tutorial couples with Aitken relaxation and a 1e-5 tolerance.)
- **Steady-state protocol:** only the steady state is computed, not the tutorial's damped transient.
  - Stage 1: steady SU2 with a rigid plate, residuals below 1e-12 (`steady_rigid.cfg`).
  - Stage 2: coupled from it, with fluid time step = coupling window = **100 s** for 40 windows (t = 4 000 s). The large step damps a slow cavity-flow mode (time scale about 100 s, of the order of depth²/ν) that small steps leave drifting for hundreds of seconds at the finest level. The steady values do not depend on the step (next section).
- **Evaluation** (`scripts/evaluate_coupled.py`, `scripts/study_table.py`):
  - u_x(A), u_y(A) from a preCICE watch-point at (4, 0) on the solid mesh, cross-checked against CalculiX's `*NODE PRINT` of the two nodes at A.
  - F_x, F_y: the converged (last) SU2 solve of each window.
  - Steady value = the final window; steadiness from the late-time range (next section).

## Mesh study (fluid and solid refined together, Δt = 100 s)

| Level | u_y(A) (m) | F_y (N/m) | u_x(A) (m) | F_x (N/m) |
|---|---:|---:|---:|---:|
| 1 | −0.24078516 | −5.3656746 | 9.096707e-5 | −0.0088271917 |
| 2 | −0.23977776 | −5.2717124 | 9.3850868e-5 | −0.013758867 |
| 4 | −0.24012430 | −5.2594312 | 9.6377986e-5 | −0.015372309 |
| 8 | −0.24036175 | −5.2579794 | 9.6937733e-5 | −0.016026244 |
| 16 | −0.24046943 | −5.2578071 | 9.69987e-5 | −0.016325298 |

Observed order p = ln((f_L − f_2L)/(f_2L − f_4L))/ln 2, Richardson value f_4L + (f_4L − f_2L)/(2^p − 1), GCI on the finest level of the triple with a safety factor of 1.25 (`cavity_su2_study.md`):

| Quantity | Levels | Observed order | Richardson | GCI (finest) |
|---|---|---:|---:|---:|
| u_y(A) | 1-2-4 | not monotone | – | – |
| u_y(A) | 2-4-8 | 0.55 | −0.24087862 | 0.27 % |
| u_y(A) | 4-8-16 | 1.14 | −0.24055878 | 0.046 % |
| F_y | 1-2-4 | 2.94 | −5.2575846 | 0.044 % |
| F_y | 2-4-8 | 3.08 | −5.2577848 | 0.005 % |
| F_y | 4-8-16 | 3.07 | −5.2577839 | 0.001 % |
| u_x(A) | 1-2-4 | 0.19 | 1.1428e-4 | 23 % |
| u_x(A) | 2-4-8 | 2.17 | 9.7097e-5 | 0.21 % |
| u_x(A) | 4-8-16 | 3.20 | 9.7006e-5 | 0.010 % |
| F_x | 1-2-4 | 1.61 | −0.016157 | 6.4 % |
| F_x | 2-4-8 | 1.30 | −0.016472 | 3.5 % |
| F_x | 4-8-16 | 1.13 | −0.016577 | 1.9 % |

- **F_y** converges at about third order, consistently over all three triples: **F_y = −5.25778 N/m**.
- **u_y(A)** is not monotone at 1-2-4 (L2 is the least deflected) and its order is still rising (0.55, then 1.14), so its Richardson value is approximate. Best estimate **u_y(A) = −0.2405 m** (L16 −0.24047, extrapolated −0.24056).
- A likely reason for the u_y behaviour: the magnitude of the fluid load decreases with refinement (at about third order), while a separate CalculiX check of the plate alone under a fixed load (not included) deflected more with refinement, at about order 1.3. The two contributions to u_y then have opposite signs and different orders.
- **F_x** (a small force, about 0.3 % of |F_y|) converges at first order; u_x(A) is of order 1e-4 m.

## Time-step independence and steadiness

Runs at Δt = 5, 25 and 50 s (L1, L4) and 5 s (L16), each compared with Δt = 100 s (`dt-check/dt_table.md`, from `dt-check/dt_check.csv` and `histories/` with `scripts/dt_table.py`):

| Level | Δt (s) | Windows | Δu_y(A) (m) | Δu_x(A) (m) | ΔF_y (N/m) | ΔF_x (N/m) | max rel. Δ |
|---|---:|---:|---:|---:|---:|---:|---:|
| 1 | 5 | 50 | +3.70e-08 | −4.59e-11 | −1.40e-08 | +3.32e-09 | 5.0e-07 |
| 1 | 25 | 40 | 0 | 0 | 0 | 0 | 0 |
| 1 | 50 | 40 | 0 | 0 | 0 | 0 | 0 |
| 4 | 5 | 50 | 0 | −5.00e-13 | +5.00e-09 | +2.00e-11 | 5.2e-09 |
| 4 | 25 | 40 | 0 | 0 | −1.00e-09 | −1.00e-11 | 6.5e-10 |
| 4 | 50 | 40 | 0 | 0 | +1.00e-09 | 0 | 1.9e-10 |
| 16 | 5 | 50 | +1.13e-07 | +1.57e-09 | +5.90e-06 | +1.79e-07 | 1.6e-05 |

Δ = value at that Δt minus the Δt = 100 s value; 0 means equal to the 9–10 digits of the watch-point and history output.

- At L1 and L4 the steady values are independent of the step to 1.5e-7 (L1, u_y) and 5e-9 (L4) relative, i.e. to the output precision and the coupling tolerance. The Δt = 5 s runs end at t = 250 s, where L1 is still settling at the 1e-7 level.
- The L16 Δt = 5 s run differs from Δt = 100 s by 1.1e-7 m in u_y and 5.9e-6 N/m in F_y, far below the L8–L16 differences (1.1e-4 m and 1.7e-4 N/m). That run is not fully steady: see Notes.

Late-time range of each quantity over the last 25 Δt = 100 s windows (t = 1 600–4 000 s):

| Level | u_y(A) (m) | u_x(A) (m) | F_y (N/m) | F_x (N/m) |
|---|---:|---:|---:|---:|
| 1 | 0 | 1.6e-12 | 1.0e-09 | 2.4e-11 |
| 2 | 1.0e-09 | 9.0e-13 | 8.0e-09 | 3.0e-11 |
| 4 | 5.0e-09 | 4.4e-12 | 1.6e-07 | 6.8e-10 |
| 8 | 1.3e-07 | 2.2e-10 | 3.5e-06 | 1.3e-08 |
| 16 | 4.1e-07 (1.7e-6 rel.) | 2.2e-09 | 1.3e-05 (2.5e-6 rel.) | 8.5e-08 |

The L8 and L16 values wander at the level of the 1e-6 coupling tolerance (L16: SU2's inner stop at 1e-10 absolute leaves about 1e-5 relative noise in the nodal forces), two to three orders of magnitude below the L8–L16 differences. u_y(A) is within 1e-6 (relative) of its final value from t = 600 s (L1–L4), 900 s (L8) and 2 900 s (L16) on (`cavity_su2_evaluate.txt`).

## Rigid-plate check

Stage 1 alone (steady SU2, rigid flat plate, residuals below 1e-12; `checks/rigid_plate_steady.csv`, `.txt`): F_y = −5.4254 / −5.3210 / −5.3071 / −5.3054 / −5.3051 N/m (L1–L16), observed order 2.9–3.0. The coupled F_y at L16 is 0.9 % smaller in magnitude than the rigid one, so the plate deflection relieves part of the load. The rigid-plate field data are not included; `scripts/evaluate_rigid.py` and `scripts/run-steady.sh` are the scripts that produced the table.

## Cost and coupling

From `scripts/cost_table.py logs/<run> ...`:

| Run | Slurm job | Cores (SU2 + CalculiX) | Stage 1 wall (s) | Coupled wall (s) | Core-hours (coupled) | Windows | Mean / max coupling its | Windows at cap | SU2 inner its per solve |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| L1 | 873 | 3 (2 + 1) | 1 | 20 | 0.02 | 40 | 1.70 / 15 | 0 | 20.7 |
| L2 | 877 | 3 (2 + 1) | 3 | 70 | 0.06 | 40 | 1.65 / 15 | 0 | 25.2 |
| L4 | 876 | 5 (4 + 1) | 36 | 251 | 0.35 | 40 | 1.65 / 15 | 0 | 34.2 |
| L8 | 878 | 9 (8 + 1) | 240 | 1163 | 2.91 | 40 | 1.70 / 15 | 0 | 54.8 |
| L16 | 879 | 32 (31 + 1) | 1685 † | 3583 | 31.85 | 40 | 1.70 / 15 | 0 | 124.9 |

† L16 stage 1 is the rigid-plate check (job 860, 16 cores; `logs/L16-stage1/`), whose restart was reused: same mesh, inlet and configuration except the iteration cap (40 000 instead of 20 000; it converged in 2 370 iterations).

Core-hours count all allocated cores; CalculiX is serial and mostly waits for SU2. The cost is almost all SU2 inner iterations. The first window takes 15 coupling iterations; from the seventh window on every window converges in one (one window with three at L8).

## Comparison with the other solutions

Made after the blind values above were posted (`scripts/compare_references.py`; `comparison/reference_limits.csv`, `comparison/reference_comparison.csv`). Sources:
- **solids4foam partitioned:** St Venant–Kirchhoff, strict coupling test, Δt halved per level, as posted on #36 (`comparison/s4f_partitioned_issue36.csv`);
- **solids4foam quasi-monolithic:** the #618 and #597 builds (`../cavity_monolithic_618/`);
- **FEniCSx monolithic:** `../fenicsx_cavityflexiblebottom/data/coupled_levels.csv`. Its README describes it as a second implementation by the same group, not an independent one;
- **Tuković et al. 2018, ANSYS, COMSOL:** the paper's reference levels (`../cavity_plate_material_law/reference_levels.csv`).

Each limit is the Richardson value of the sequence's three finest levels (ours: L4–L8–L16). Differences are reference / ours − 1:

| Solution | u_y(A) limit (m) | order | vs ours | F_y limit (N/m) | order | vs ours |
|---|---:|---:|---:|---:|---:|---:|
| SU2 + CalculiX (this package) | −0.240559 | 1.14 | | −5.257784 | 3.07 | |
| solids4foam quasi-monolithic #618 | −0.240580 | 2.00 | +0.009 % | −5.257741 | 2.08 | −0.001 % |
| solids4foam quasi-monolithic #597 | −0.240581 | 2.00 | +0.009 % | −5.257746 | 2.11 | −0.001 % |
| FEniCSx monolithic | −0.240544 | 1.28 | −0.006 % | −5.257801 | 2.18 | +0.000 % |
| solids4foam partitioned (#36) | −0.240671 | 2.02 | +0.046 % | −5.260619 | 2.33 | +0.054 % |
| Tuković et al. 2018 | −0.251033 | 1.84 | +4.35 % | −5.256744 | 1.92 | −0.020 % |
| ANSYS | −0.253638 | 0.89 | +5.44 % | −5.271895 | 1.30 | +0.27 % |
| COMSOL | not monotone | | | not monotone | | |

- **Monolithic solutions:** this partitioned solution agrees with both solids4foam quasi-monolithic builds and with FEniCSx to within 0.01 % in u_y and 0.001 % in F_y. That is within our u_y GCI (0.046 %) and close to our F_y GCI (0.001 %).
- **solids4foam partitioned:** the difference is 0.05 % in both quantities. In u_y this is within either sequence's uncertainty. In F_y it is larger than our GCI (0.001 %) and than the sequence's own GCI on its finest level (about 0.02 %, with a safety factor of 1.25), while the other three limits agree with each other to 0.001 %.
- **Approach to the limit:** at the same Δx the solids4foam sequences converge from below at second order. Their u_y(A) at Δx = 0.1 m is 15 % smaller than ours, and 0.10–0.16 % smaller at Δx = 0.0125 m. Our u_y is already within 0.1 % of its limit at L1, but it is not monotone (see the mesh study). FEniCSx's u_y is within 0.75 % of ours at every level. A likely reason is that both solid discretisations represent bending well on coarse meshes: C3D8I incompatible modes here, quadratic (P2) elements in FEniCSx.
- **COMSOL:** its finest values (Δx = 0.0125 m) are 1.1 % (u_y) and 0.7 % (F_y) above ours in magnitude, and its sequence is not monotone. Five solutions (two solids4foam quasi-monolithic builds, FEniCSx, this package and the solids4foam partitioned sequence) agree to within about 0.05 %, which places the COMSOL difference with COMSOL.
- **Tuković et al. 2018 and ANSYS:** they differ by 4–5 % in u_y and agree with ours in F_y to 0.02 % (Tuković) and 0.27 % (ANSYS). This repository already discusses those differences and they are not examined further here.

## Notes

- **Why Δt = 100 s (the L16 Δt = 5 s run):** the first L16 sequence used Δt = 5 s for 50 windows (job 869, `logs/L16-dt5/`). Window 21 reached the 50-iteration coupling cap: the displacement residual converged, but the force residual stalled at about 1e-5, the nodal-force noise left by SU2's absolute 1e-10 inner stop at L16, above the 1e-6 force tolerance. The accepted window excited the slow cavity mode, and from window 24 every window converged in one coupling iteration while the state crept monotonically (u_y by up to about 1e-7 m per window) to t = 250 s (`logs/L16-dt5/solid-calculix/precice-Solid-*.log`, `histories/dt-check/L16-dt5.csv`). A rerun with preCICE `<min-iterations value="2"/>` gave the same creep (not included). Larger steps damp that mode within a few windows (Δt = 100 s: u_y settles to 1e-6 in 6–9 windows at L1–L8), so the whole sequence was run at Δt = 100 s.
- **Time-step check runs** (`scripts/cost_table.py`): L1 at Δt = 5/25/50 s: jobs 865/871/872, 41/20/20 s coupled on 2 + 1 cores; L4: jobs 867/874/875, 462/310/281 s on 4 + 1 cores; L16 at Δt = 5 s: job 869, 7 517 s on 31 + 1 cores (66.8 core-hours, one window at the cap). The L1 and L4 Δt = 100 s runs of the main table are the Δt-scan runs (jobs 873 and 876).
- **Adapter marker name:** the SU2 adapter hard-codes the coupled marker name `interface`; the tutorial blockMesh already names the plate patch `interface`.
- **Unsteady inlet profile:** unsteady SU2 reads the inlet profile with a time-step suffix (`inlet_profile_00000.dat`); `run-case.sh` copies it.
- **Watch-point and CalculiX output:** the watch-point log's time starts at the coupled stage. CalculiX's `*NODE PRINT` at A (`solid-calculix/plate.dat`, 7 digits) agrees with the watch-point at every level (`cavity_su2_evaluate.txt`).

## Files

- `cavity_su2_steady.csv`: one row per level (Δt = 100 s): steady values, late-time changes, settling times, coupled wall clock, windows, coupling and inner iterations (`evaluate_coupled.py --csv`).
- `cavity_su2_evaluate.txt`, `cavity_su2_study.md`: the evaluation output and the mesh-study tables (`study_table.py`).
- `histories/L*.csv`: u_x(A), u_y(A), F_x, F_y per window for the Δt = 100 s runs; `histories/dt-check/*.csv`: the same for the time-step check runs.
- `dt-check/`: `dt_check.csv` and `evaluate.txt` (all runs), `dt_table.md` (the two tables above).
- `checks/rigid_plate_steady.csv`, `.txt`: the rigid-plate steady sequence.
- `comparison/`: `reference_limits.csv` and `reference_comparison.csv` (`scripts/compare_references.py`), and `s4f_partitioned_issue36.csv`, the solids4foam partitioned sequence as posted on #36.
- `case/`: the template: SU2 `steady_rigid.cfg` and `coupled.cfg`, CalculiX `plate.inp` and `config.yml`, `precice-config.base.xml`, and the tutorial fluid blockMeshDict. The `@…@` placeholders, meshes and inlet profile are filled in by `scripts/setup_case.py`.
- `scripts/`:
  - `setup_case.py`: generates a case for a level (fluid mesh with `make_fluid_mesh.py`, inlet profile with `make_inlet_profile.py`, solid mesh with `solid_mesh.py`). It regenerates the L1 and L2 cases used here byte for byte (meshes, configuration files, inlet profile, solid mesh).
  - `run-case.sh`: Slurm batch script, the FSI1 one with only names and comments changed: stage 1 (steady rigid plate), stage 2 (coupled), timing, hang guard and progress watchdog. `run-steady.sh`: stage 1 alone (the rigid-plate check).
  - `evaluate_coupled.py`, `study_table.py`, `dt_table.py`, `cost_table.py`: the evaluation; `compare_references.py`: the comparison with the other solutions; `evaluate_rigid.py`: the rigid-plate table (needs the full stage-1 output); `check_su2_mesh.py`: mesh check.
  - `env-su2-precice.sh`: the environment (as `../fsi1_su2_calculix`; the SU2 build fix is `../fsi1_su2_calculix/scripts/su2-7.5.1-cstdint.patch`).
- `logs/<run>/`: for each run the Slurm output (`slurm.out`), case summary, exit statuses, wall-clock times, the SU2 coupled history (`fluid-su2/history_00002.csv`), the preCICE coupling-iteration, convergence and point-A watch-point logs, CalculiX's `plate.dat` (`*NODE PRINT` at A), and the last lines of the solver logs. The layout mirrors a case folder, so `evaluate_coupled.py logs/<run>` reproduces the tables. The repository `.gitignore` excludes `*.out` and `*.log`, so these files were added with `git add -f`. Full solver logs and field data are available on request.

To reproduce one level (Slurm; from a study directory holding `scripts/` and `case/` copied from here, with OpenFOAM for `blockMesh` and the paths in `scripts/env-su2-precice.sh` edited):

```bash
python3 scripts/setup_case.py 1 cases/L1 --dt 100 --t-end 4000 --write-every 10
sbatch --ntasks=3 --mem=8G --time=01:00:00 --export=ALL,CASE=cases/L1 scripts/run-case.sh
python3 scripts/evaluate_coupled.py cases/L1 --csv results.csv --hist-dir histories --tail 10
```

and for the tables from the files here:

```bash
python3 scripts/evaluate_coupled.py logs/L1 logs/L2 logs/L4 logs/L8 logs/L16 --csv cavity_su2_steady.csv --hist-dir histories --tail 10
python3 scripts/study_table.py cavity_su2_steady.csv
python3 scripts/dt_table.py dt-check/dt_check.csv histories
python3 scripts/compare_references.py
```
