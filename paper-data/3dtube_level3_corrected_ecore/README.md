# 3dTube level 3 on the corrected meshes, desktop (2026-11-20)

Level 3 of the `3dTube` mesh study on the corrected meshes (full-precision arc and vertex points), as requested on issue #57 (17 Nov). Same machine as `../3dtube_level3_ecore/` (#64).

## Platform

- Host: ECORE, AMD Ryzen Threadripper 9970X, 32 cores / 64 threads, 125 GiB RAM, Ubuntu 24.04.4 LTS.
- GCC 13.3.0 (`linux64GccDPInt32Opt`), Open MPI 4.1.6, OpenFOAM v2512, PETSc 3.21.2 with hypre 2.31.0 (as in #64).
- solids4foam: `development` at `206437983`, with `system/fluid/blockMeshDict` and `system/solid/blockMeshDict` replaced by the corrected files from `9b52c2303`. It was built into its own platforms directory, so the #64 build stays untouched (`scripts/env-tube-l3.sh`).
- Environment: `OMP_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1` (recorded in `logs/log.environment`).

## Command

```bash
source scripts/env-tube-l3.sh
cd tutorials/fluidSolidInteraction/3dTube/verification
./Allverify --study mesh --levels 3 --cores 32
```

`scripts/run-tube-l3.sh` is the exact launcher (run in tmux). The run went from 16:29 to 23:27 on 20 November 2026. The machine was otherwise idle (load average about 0.3 before the start).

## Driver fix for the post-processing

The solver completed all 3 200 steps, but the driver then stopped with

```
... solidPointDisplacement_verificationWallZ050.dat has 3200 time levels, expected 3200 at dt = 6.25e-06
```

The monitors write the time with 6 significant digits, which cannot resolve Δt = 6.25e-6 s beyond t = 10 ms (for example, 0.01000625 is written as 0.0100062). The spacing check (0.1 %) then fails. Levels 1 and 2 are not affected, because their times are exact in 6 digits.

This is the problem fixed in `4bf07e18d` ("snap monitor times to the time-step grid"). That commit is on `verification/3dtube-level3` but not in `development` at `206437983`. I applied only its time-snapping part (`require_complete` and the iteration lookup in `coupling_summary`) as `scripts/tube-l3-driver-time-snap.patch`. I then reran the post-processing on the finished case with `--reuse` (`logs/log.Allverify.reuse`). The stored run settings matched, so the solver was not rerun. `logs/log.Allverify` is the original driver log with the error.

## Results

| Level | Fluid / solid cells | Δt (s) | u_r,max(A) (mm) | t(u_r,max) (ms) | u_z,min(A) (mm) | t_arr(A) (ms) | c_p (m/s) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 (corrected, Philip) | 16 000 / 6 400 | 2.5e-5 | 0.15958 | | | | 4.664 |
| 2 (corrected, Philip) | 128 000 / 51 200 | 1.25e-5 | 0.15924 | | | | 4.651 |
| **3 (corrected, this run)** | 1 024 000 / 409 600 | 6.25e-6 | **0.15831** | 7.251 | -0.08570 | 5.823 | **4.625** |

c_p is unshifted (no `--probe-z-shift`), as for the corrected L1/L2.

All 3 200 time steps met the robin convergence criteria. The driver's checks pass except "at least two levels are required", which is expected for a single-level run (see `verification_summary.md`).

## Cost (wall clock, 32 ranks)

| Level | Steps | Mean FSI iter. | Wall clock (s) |
|---:|---:|---:|---:|
| 3 | 3 200 | 3.18 | 25 060 (6 h 58 min) |

## Files

- `verification_summary.md`: driver checks.
- `mesh_robin_backward.csv`: driver output (one row, level 3).
- `mesh_robin_backward.png`: point-A history against the published curves.
- `histories/robin_backward_mesh3.csv`: point-A histories (radial and axial displacement, axis pressure), one row per time step.
- `postProcessing/`: the case's raw monitor files. These are the wall displacements at the nine stations (`0/solidPointDisplacement_*.dat`), the axis pressure probes (`verificationAxisPressure/fluid/0/p`), the solid and fluid forces, and the FSI convergence data and residuals.
- `verification_settings.json`: the run settings the driver recorded for the case.
- `logs/log.solids4Foam.xz`: the solver log (53 MB uncompressed).
- `logs/log.blockMesh.{fluid,solid}`, `logs/log.decomposePar.{fluid,solid}`, `logs/log.Allrun.case`: meshing, decomposition and the case's Allrun log.
- `logs/log.Allverify`, `logs/log.Allverify.reuse`, `logs/log.environment`: the driver logs (original and post-processing rerun) and the environment.
- `scripts/env-tube-l3.sh`, `scripts/run-tube-l3.sh`: environment and launcher.
- `scripts/tube-l3-driver-time-snap.patch`: the driver patch described above.

## Point-A history figure (`figure/`, Fig. `fig:3dtube_histories`)

`figures/3dtube_pointA_histories.pdf`: radial and axial wall displacement at point A for levels
1–3 on the corrected meshes (Robin–Neumann, BDF2), with the published curves. It replaces the
gnuplot draft `mesh_robin_backward.png` for the paper.

| File | Content |
|---|---|
| `figure/tube_histories_figure.py` | Draws the figure in the shared style (`../figstyle`). Inputs by argument. |
| `figure/data/robin_backward_mesh{1,2}.csv` | Levels 1 and 2: the driver's point-A histories (`postProcessing/histories/` of the corrected L1/L2 study in `~/claude-runs/tube3d/post671/3dTube/verification`), every time step. Checked identical to `solidPointDisplacement_verificationWallZ250.dat` (columns Dy, Dz) of the cases in `work/robin_backward_mesh{1,2}`. |
| `histories/robin_backward_mesh3.csv` | Level 3 (this run), used as is. |
| `figure/data/reference/` | The digitised published curves and the reference JSON of the solids4foam `3dTube/verification/reference` directory (same copy as above): Tuković et al. (2018) Fig. 25(a), Δt = 2.5e-5 s, radial only, to 10 ms, from the vector graphics; Lozovskiy et al. (2019) Fig. 3 (finer mesh) and Eken (2016, PhD thesis, Fig. 5.6, mesh M2; companion of Eken & Sahin 2016), implicit Euler, Δt = 1e-4 s, raster-digitised to about ±3e-6 m. Provenance in each file's header. |

```bash
cd paper-data/3dtube_level3_corrected_ecore/figure
python3 tube_histories_figure.py \
  --level "L1=data/robin_backward_mesh1.csv" \
  --level "L2=data/robin_backward_mesh2.csv" \
  --level "L3=../histories/robin_backward_mesh3.csv" \
  --reference-dir data/reference --out ../../../figures/3dtube_pointA_histories.pdf
```

The script prints u_r,max and u_z,min (t < 6.5 ms) per level: 0.15958, 0.15924, 0.15831 mm and
−0.08817, −0.08657, −0.08569 mm, as in Table `tab:3dtube` (level 3 u_z,min there is −0.08570,
from the driver's interpolated extremum).
