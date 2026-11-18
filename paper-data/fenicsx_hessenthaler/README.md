# Hessenthaler Phase I: second implementation (FEniCSx) and matched early-time comparison

Status: **provisional, evidence only. Not a completed Hessenthaler verification, not an MRI validation, not a second-group or independent verification (same research group, same AI model lineage).**
Text: `../../supplement_fenicsx_hessenthaler.tex` (builds alone with `latexmk -pdf supplement_fenicsx_hessenthaler.tex`).

What this directory is: the data, scripts, case files and provenance behind the supplement. The solids4foam results are **our own reconstruction** of the protocol that the main text (`sec:hessenthaler`, repo HEAD 4dfee64) describes for the feature-edge campaign; they are **not** the corrected feF1-feF3 runs of the author, which were not available. Where the two could differ is listed in the supplement (Table S1) and in `s4f_case/`.

## Contents
| path | content |
|---|---|
| `data/s4f_feF1_ours_0_0.75.csv` | per-step data of our solids4foam run (375 steps, t = 0.002 ... 0.75 s): tip displacement, fluid- and solid-interface force totals, coupling iterations, last FSI residual |
| `data/s4f_flapTip.dat`, `data/s4f_log_extract.txt`, `data/s4f_full_log_gz.sha256` | function-object output; extract of the solver log (times, force totals, FSI residuals, coupling iterations); sha256 of the full compressed log (6.3 MB, not included) |
| `data/Gmi_n32b_history.csv`, `data/Gmi_dt01_history.csv`, `data/Gc_T3b_history.csv` | FEniCSx step histories (Gmi: t <= 0.69 s at dt 0.02 s, then 0.01 s to 0.75 s; Gc: development level, to 1.21 s) |
| `data/centreline_0p69_*.csv` | deformed-centreline data of both codes at t = 0.69 s |
| `data/matched_comparison_Gmi_Gc_vs_ours_s4f.csv`, `data/common_times_tip_forces.csv`, `data/summary_stats.json`, `data/centreline_0p69_stats.json` | matched-time tables and statistics (as produced by the campaign scripts) |
| `data/checks/` | standalone checks: material tests, zero-flow solid levels (S1-S3, mu = 64.2 and 60 kPa), geometry gates (Gc, Gmi, G0), fluid-only gate C, force-decomposition states, Jacobian check, solver tests (preconditioner failures preserved) |
| `figures/` | figures of the supplement (pdf) |
| `scripts/regenerate_from_packaged_data.py` | regenerates `tables_matched.tex`, statistics and the three figures from `data/` only (numpy, matplotlib) |
| `scripts/extract_s4f.py`, `compare_matched.py` | campaign scripts that produced the CSVs from the solids4foam case and the FEniCSx histories (absolute paths of the original machine) |
| `fenicsx_source/` | the FEniCSx code that produced the histories (3-D neo-Hookean, monolithic ALE FSI, drivers, mesh scripts) |
| `s4f_case/` | case files (deltas to the pinned tutorial), build, mesh and run scripts, mesh checks and hashes of the solids4foam reconstruction |
| `prereg/` | preregistration, addenda A1-A7 and their sha256 (each written before the corresponding runs) |
| `ledger/` | CPU-hour ledger reconciled with sacct |
| `SOURCES.md`, `SHA256SUMS` | versions, source commits, input and output hashes |

## Reproduce the tables and figures
From this directory: `python3 -B scripts/regenerate_from_packaged_data.py regenerated` (needs numpy and matplotlib); compare `regenerated/tables_matched.tex` with `tables_matched.tex` (identical apart from the first comment line) and the figures with `figures/`. Then `latexmk -pdf supplement_fenicsx_hessenthaler.tex` in the repository root.
Not reproducible from this directory without the original machine: the solver runs themselves (cost, see supplement; the FEniCSx checkpoints are rank-local) and mesh generation (scripts and hashes are included, the surface file is the published one from the solids4foam tutorial).

## Path adaptation of shipped scripts
The scripts, case files and logs were adapted for shipping: absolute paths of the authoring machine were replaced by `HESS_ROOT` (campaign root, FEniCSx scripts; `os.environ`, default `.`) and `S4F_ROOT` (solids4foam work tree, shell scripts and `s4f_case/make_case.py`), `#SBATCH -o` paths made relative, and the log headers and the A6 sentence that named the home directory anonymised. `MANIFEST_path_adaptation.csv` lists, for every adapted file, the sha256 of the as-run version, the sha256 of the shipped version and the number of replacements. The FEniCSx files in `fenicsx_source/` are shipped flat (the as-run layout was `src/hessenthaler/`, `benchmarks/hessenthaler/`, `geometry/scripts/`), so the scripts are a record of the code that produced the results, not a ready-to-run package; the regeneration of tables and figures from `data/` is independent of them.

## Not included
Meshes and checkpoints, the full solids4foam log, the solids4foam and OpenFOAM sources (commit and version in `SOURCES.md`), MRI data (not used here).
