# blobInTreacle: second implementation (FEniCSx) and interface-formulation audit

Status: **draft supplementary evidence, same research group, not an independent verification by a different group.**
Text: `../../supplement_fenicsx_blobintreacle.tex` (builds alone with `latexmk -pdf supplement_fenicsx_blobintreacle.tex`). Related: issue #50 (this is *not* the FEBio run that issue asked for).

The wording follows an internal review of the evidence package by a further agent of the same model family, followed by a second internal quality gate (neither is an external or human review; reviewed state of the code repository: master ee9379f plus uncommitted run directories). The first review accepted the steady-state sequences, the M2 transmission loss and the limited BDF2 result, and rejected the claims listed in "What is not claimed" below.

## Contents
| file | content |
|---|---|
| `data/steady_levels.csv` | per-level steady values, four models (FEniCSx M2-BE, FEniCSx M1-BE, solids4foam M2 original, solids4foam E-AREA diagnostic), with source document per row (hand-transcribed from the signed-off run records) |
| `data/steady_convergence.csv` | observed order, Richardson limit, GCI (finest and coarsest triple) |
| `data/steady_crosscode.csv` | FEniCSx M2 against solids4foam M2: same level and limits |
| `data/modelform_M1_M2.csv`, `data/fourway_M1_M2.csv` | M1/M2 sensitivity; the four-way table |
| `data/time_values_t1_F1_M2.csv`, `time_orders_t1_F1_M2.csv`, `time_errors_vs_finer_reference_t1_F1_M2.csv`, `time_rms_history_t0_1_F1_M2.csv` | temporal study at t = 1, F = 1, M2 |
| `data/bdf2_transient_spatial_t1_M2.csv`, `bdf2_transient_spatial_t1_orders.csv`, `bdf2_transient_spatial_t1_crossF_ux.csv` | BDF2 transient spatial study at t = 1 (dt 0.005, F = 0.5 to 4; solver repo commit 7b4e4e1, code 1484d94), orders and limits, cross-F coincidence table |
| `data/earea_ratio_reconciliation.csv` | E-AREA / M2 ratios, like-for-like (same triple) against mixed-level |
| `data/s4f_time_sensitivity_t1.csv` | solids4foam time-step sensitivity at t = 1 (no order identified) |
| `data/interface_force_bins_M2_F4.csv` | M2 load on the solid per 15-degree bin, F = 4 / x4 (coupled and fixed-shape fluid) |
| `data/force_conservation_evidence.csv`, `data/uy_sensitivity_chain.csv` | closure, transmission and power evidence; the chain behind the 0.7% uy gap |
| `scripts/convergence.py` | the campaign's convergence routine, copied unchanged (67 lines of standard-library Python) |
| `scripts/analyze_steady.py` | regenerates `steady_convergence`, `steady_crosscode`, `modelform`, `fourway` CSVs from `steady_levels.csv` |
| `scripts/make_figures.py` | figures from the CSVs (numpy, matplotlib) |
| `figures/` | `fenicsx_blob_spatial`, `fenicsx_blob_temporal`, `fenicsx_blob_transient_spatial` (pdf, png) |

Not included: raw solver histories, state dumps, meshes (regenerated and hash-checked, see below), the solver source.

## Reproduce the tables and figures
Use `python3 -B` (or `PYTHONDONTWRITEBYTECODE=1`) so that no `__pycache__` is left behind (also ignored by `.gitignore`).
```
cd paper-data/fenicsx_blobintreacle
python3 -B -I scripts/analyze_steady.py   # steady_convergence, crosscode, modelform, ratio table, base-reaction ratios, fourway
python3 -B -I scripts/analyze_time.py     # bdf2_pairwise_vs_prediction.csv
python3 -B -I scripts/make_figures.py     # needs numpy and matplotlib
python3 -B scripts/transfer_ratios_from_raw.py SOLVER_REPO data/interface_transfer_ratios.csv   # needs the solver repository
```
The hand-transcribed CSVs and the report each comes from (with sha256) are listed in `data/SOURCES.md`.
The derived CSVs reproduce the sign-off values (e.g. M2 u_x p 1.414, limit 0.120286, GCI 0.0477 %; M1 u_x p 1.519, 0.129372, 0.0319 %). The transcribed inputs were spot-checked against raw histories by the reviewer (10+ values, no discrepancy).

## Review lineage and process detail (moved from the supplement)
Evidence produced, run and verified by AI agents of one model family for the same research group; two internal reviews (a referee-style review and a quality gate; neither external nor human). Preregistrations (26, sha256) were verified by the first review. These details are not repeated in the supplement beyond one disclosure paragraph.

## Provenance
- Solver: monolithic ALE FSI, DOLFINx 0.10.0, unknowns (u Q2, p Q1, d Q2) on one mesh; PETSc SNES. The solver repository (`fenicsx-fsi-verification`) is not part of this PR; commits below refer to it.
- Frozen code commits: **c8a4b0a** (M2-BE, F = 0.5, 1, 2, dt 0.01, t = 10; BE time study); **1484d94** (M2-BE F = 4, dt 0.04, t = 6; BDF2 levels dt 0.04 to 0.005; differs from c8a4b0a only in options that default to the old behaviour); **6f85a48** (M1 via `run_blob_b1.py`, B1 study; consistent CN); **d17df8b** (BDF2 dt 0.0025 and 0.00125 references); **9f4b6e8** (CN with `--theta-pressure-geometry`); **722429d** (restart tests); solids4foam **1576edcc4**, OpenFOAM v2512 (E-AREA is a user library, upstream untouched: 26 md5 sums re-checked).
- Mesh sha256 (quad meshes `meshes/generated/blob_f*_quad.msh`):
  - F0.5 `39eebdbb51f85406dc6e252d5e452a5fddd947f130a0432b544b3129aa521a80`
  - F1 `2c8ce4d9fdc712a472c79562a2c6c518e96b321582735c8d1196f66e732d5c26`
  - F2 `599a81e64ca630c6356b1fe170570debcdde811b2b860d84cd8f490bb7147796`
  - F4 `50e9ac2d4e693dd1f1b5564cce2b60561895af3fcddeb6124fd10b7c30b6507e`
- Environment: micromamba env `fenicsx010` from an explicit lock (python 3.11.15, dolfinx 0.10.0, basix 0.10.0, ufl 2025.2.1, PETSc/petsc4py 3.25.6 real, MUMPS 5.8.2, hypre 3.2.0, numpy 2.4.6, MPICH 5.0.1, mpi4py 4.1.2); lock file sha256 prefix `72764221...b33`, wrapper `environment/run_in_env.sh` sha256 prefix `321dc817...006`.
- Launcher / Slurm: `slurm/submit_run.sh` -> `slurm/job.sbatch` (partition `main`, `-n N -c 1`, N = 2, 4, 8 for F = 1, 2, 4), `run_in_env.sh mpiexec -n N python3 -u slurm/launch_driver.py ...`, node xenosim (shared). Each run directory has `SUBMIT_CMD`, `HOST.txt`, `MANIFEST.json` (code commit, input hashes) and `EXIT` written by the driver (the driver status, not Slurm's, is authoritative: a Slurm "COMPLETED" job with solver rc 1 occurred once on the solids4foam side and is recorded).
- Example command (M2, F = 1, backward Euler): `run_blob.py --mesh blob_f1_quad.msh --dt 0.01 --t-end 10 --theta 1.0 --interface-measure M2 --diag`; M1 uses `run_blob_b1.py ... --interface-measure M1 --dt 0.04 --t-end 6`; BDF2 runs `--theta 0.5 --fluid-scheme bdf2`.
- Historical solids4foam results (`steady_levels.csv`, `steady_crosscode.csv`, ...) use the legacy wall BC (`newMovingWallVelocity`, zero wall-normal gradient of the wall-normal velocity), build 1576edcc4. They remain valid for the legacy BC.
- Cost: FEniCSx campaign 46.9 CPU-h (all jobs including failed and cancelled; M1 mesh study 3.2 CPU-h, M2 backward-Euler family 13.2 CPU-h); solids4foam side runs about 2.2 CPU-h for the t = 1 time series plus the regenerated mesh series.
- Preregistrations (26, sha256 files all re-verified) are in the solver repository.
- Reproduction check: fresh archives of c8a4b0a and 6f85a48 reproduce the F = 1 runs to 4e-15 relative at every step.

## Retained failed or excluded runs
- A1 reference R1 (BDF2, dt 0.0025): the first launch, Slurm job 16419, run directory `runs/a1_R1_bdf2_dt0.0025_n2_20261106T144042Z_ecc75a`, ended `EXIT FAIL SystemExit(2)` and the solver never started; the valid R1 used throughout is job 16450, `..._20261106T145125Z_975b2a`.
- A1 T1 restart runs at dt 0.04 (BDF2 and CN) end at t = 0.98 and are excluded from the t = 1 analysis (t = 1 is not on their time grid).
- Also retained in the solver repository: cancelled and failed M2 and M1-CN runs (M1 CN failed or was killed at F = 2, at F = 1 with dt 0.02 and 0.005; only F = 0.5, 1 with dt 0.01 are valid), cancelled pinned M2 runs, and a mislaunched solids4foam run (Slurm COMPLETED with solver rc 1).

## solids4foam PR #671 rerun (proposed branch, NOT an official release)
- Source: solver repository `s4f_side/WALL671_RESULTS.md` (commits 6452eae, 5b004a2), raw tables `s4f_side/w671/WALL671_TABLES.md`, `WALL671_DIST.md`, `WALL671_TUT_TABLES.md`, `data/provenance.csv`, `BUILD_INFO_head.txt` (sha256 in `data/SOURCES.md`). Preregistrations PREREG_WALL671*.md were committed before the runs.
- Build: HEAD da1044a4f2e1af7bc0c645cffad9336fc3aa5c4c (pull/671/head, unmerged); BASE 6a49bf0a914b4ec184bd0f717efd0eccb004db0a (parent of the first #671 commit). OpenFOAM v2512, PETSc 3.24, g++ 11.4.0, same Allwmake recipe as the 1576edcc4 build (only -j2). Executable `solids4Foam` md5 48ce4dd3b27f0c34bf07655e89f1bcf7 (identical to the 1576edcc4 build); `libsolids4FoamModels.so` md5: HEAD dcd1f8fa32e9119f69740d8fd85d3ed3, BASE 4aa620aa3728655e89bf69b53000a330, 1576edcc4 build 4cc785b78122a8be56d43a98f07ad0ac. The head is not 1576edcc4 plus the BC: 33 upstream commits lie between (the commit list file in the solver repository was not recounted here); the #671 change proper is 6 commits touching `newMovingWallVelocity*` and `elasticWallVelocity*` in src/.
- Arms: LEGACY = `continuityNormalGradient no;`, CORRECTED = key absent (the #671 default), same executable. Gate G0: the legacy arm reproduces the accepted legacy results bitwise (four runs; also the BASE build). 46 blob-matrix run directories and 11 tutorial runs, validity from solver logs, not Slurm codes.
- Data here: `data/s4f_w671_levels.csv` and `s4f_w671_force_ratio.csv` are extracted from `WALL671_TABLES.md` by `scripts/extract_w671.py SOLVER_REPO out.csv` (not hand-transcribed); `scripts/analyze_w671.py` writes `s4f_w671_convergence.csv`, `s4f_w671_effect.csv`, `s4f_w671_crosscode.csv`, `s4f_w671_ratios.csv`. All reproduce the numbers of `WALL671_RESULTS.md` (limits, gaps, ratios, 1.24 and 5.1 and 14 times the combined GCI). The `gap` columns use s4f/FEniCSx - 1 (campaign convention), opposite in sign to FEniCSx/s4f - 1 used in the supplement's cross-code table.
- Failed preregistered predictions (reported as failed): the BC effect on blob displacement mesh independent (spread 0.54 pp in u_x); the u_y gap to FEniCSx closing to < 0.2 % (0.70 % to 0.53 %, below the 2u resolution of 0.21 %); identical mesh quality (differences <= 1.5 %); collapsibleChannel effect 0.5 to 3 % (0.08 %); Hron-Turek FSI1 steady 2 to 10 % reduction (<= 0.2 % at 2x). Held: legacy bitwise; forces < 1 %, dp < 0.5 %, transmission and power ratios unchanged; E-AREA/M2 ratios unchanged; dt-independence; Hron-Turek FSI3 t = 2.1 tip displacement -19 % reproduced.
- Not blob, not in the supplement: collapsibleChannel IQN-ILS dt 0.025 (full run, 0.08 %), Robin dt 0.00625 (start-up window t <= 0.325 only, cancelled at the CPU-h cap), Hron-Turek FSI3/FSI1; these are in `WALL671_RESULTS.md` sections 3 and 4 only.
- Open: E-AREA x4 (the preregistered trigger was met by |Delta u_y| = 0.56 % at x2; about 10 CPU-h; not run); origin of the remaining about 0.5 % u_y gap to FEniCSx; same-shape and static-wall tests not repeated with the #671 arm.
- Internal inconsistencies found in `WALL671_RESULTS.md` (not edited there): section 2.6 says force-distribution bins were not generated while section 5 reports them (WALL671_DIST.md exists; the section 5 statement is used); the legacy-arm x1 steady run differs from the accepted x1 run by 3e-5 relative in u_x (partition noise), so legacy limits are 0.1203491 against the accepted 0.1203488 ('exactly' in the summary is loose); F_x change at x1 is -0.52 % (summary says <= 0.5 %); the transient limit comparison legacy 0.19 % is 1.41 times the combined GCI (the summary calls it consistent, which holds only under the 2u rule).

## Mixed composition of the M2 sequence (disclosed)
F = 0.5, 1, 2: dt 0.01, t = 10, commit c8a4b0a; F = 4: dt 0.04, t = 6, commit 1484d94. The steady states are time-step independent to <= 4e-10 (8 printed digits), so the mixture does not enter the GCI at the printed level. Pinned and unpinned MPI wrappers were used (speed only). `git_dirty=true` in some manifests is from deleted environment symlink files only.

## What is not claimed
No second-order temporal convergence for blobInTreacle; no independent verification (different group); no statement that M1 or M2 is the correct model, or that solids4foam is defective; no uy agreement within a propagated uncertainty (it agrees to 0.7 %, not resolved below about 1 %); no cross-code transient agreement (the t = 1 u_x limits differ by 1.41 x combined GCI; u_y and F_x not accepted; the time error at dt = 0.005 is not isolated) and no interface-shape comparison; the x-force statement is a difference between the fluid-reported force on the current surface and the solid-applied force on the reference areas, not an unconditional physical force loss; E-AREA is a same-lineage diagnostic, not a verification of the coupling; no level-3 claim.

## Known open points (from the review)
Pointwise stress closure misses by 1.9 to 2.9 % (explained as a corner strong-form defect, about h^0.6; the fluid-side weak-vs-facet difference and the M2 correction-term share of the base reaction are not decomposed); x4 E-AREA level not run; Liu steady/transient values not obtained from the journal article; run directories `runs/a1_*` of the solver repository untracked at the review date and two R1 reference run directories (one used) not yet recorded in its failed-runs register; no review by a human outside the lineage.
