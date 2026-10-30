# Reproducibility audit

Three levels, stated separately and without softening. "Package alone" = python3 (numpy, scipy, matplotlib) and the files in this directory, run from the package root. Test record: `provenance/REPRODUCIBILITY_TEST_RECORD.md` (every script below was run in a scratch copy of the package).

## 1. Reproducible from the package alone
| what | how | result of the scratch run |
|---|---|---|
| all FEBio inputs (20 inputs incl. the 6 MB M3 input) regenerate from the generators and match the recorded sha256 | `scripts/check_inputs.py` | 18 identical; 2 differ (early failed attempts `B0_fluid_f1`, `B0_fluid_straight`; monitor node sets only, INCONSISTENCIES I-9); shipped files match their recorded file hashes |
| every headline number (u_y, u_x, F_y, mass imbalance, F1..F4, D, S1/S2, CPU-h) re-derived from the raw stored per-run outputs | `scripts/verify_against_raw.py` -> `results/VERIFICATION_LOG.csv` | 68/68 within tolerance |
| all summary tables in `results/*.csv` (coupled summary, solver termination and cost from the logs, mesh convergence, sensitivity, parity, comparison with references, stabilisation tables, plate-only) | `scripts/build_tables.py` | rebuilt |
| force reconstruction F1..F4 and the row-parity decomposition at the stored states | `scripts/reconstruct_forces_from_facets.py` on `results/raw/facets/*.npz` | equals `results/raw/pressure_mode/b3_forces_summary.csv` to <=4e-7 |
| figures `figures/fig_*.png` | `scripts/make_figures.py` | rebuilt |
| discrete-operator analyses (Schur spectra and h^2 scaling, backward-Euler Stokes emulation, straight-channel emulation, projection matrix checks) | `scripts/b4_operator.py`, `b4_modal_scaling.py` (part 1), `b4_straight_channel.py`, `b4_straight_profile.py`, `b4_transient_stokes.py`, `b7_ppp_check.py`, `b7_stabilisation_check.py` | identical to the stored JSON in `results/raw/pressure_mode/` (relative 1e-8) |
| plate-only Q9 response to the stored pair-mean / raw wall pressure | `scripts/run_plate_q9_pairmean.py` | identical to `results/raw/b2/plate_q9_pairmean_M3.json` |
| Poiseuille inlet pre-accounting (analytic) | `scripts/poiseuille_check.py pre` | runs |
`bash reproduce.sh` runs all of the above (about 4 minutes on one login-node core).
**Caveat that applies to all of it:** the "raw" files are the stored outputs of the original runs. Recomputing tables from them checks bookkeeping and arithmetic; it is **not** an independent solver verification and does not re-solve any equation of the benchmark. The facet tables are a re-evaluation of stored FEBio fields.

## 2. Needs the retained HPC outputs (not shipped; 0.51 GB; names, sizes and sha256 in `provenance/hpc_retained_outputs.csv`)
| script | needs | largest files |
|---|---|---|
| `b3_force_reconstruction.py` (creates the facet tables), `cavity_xplt_audit.py`, `b2_ybranch.py`, `b3_momentum_balance.py`, `b2_febio_history.py` | `*.xplt`, `.feb`, `model_meta.json`, `Mon_iface/inlet/outlet.txt` of C2, C3_fsi_M2, C3_fsi_M4 | M3 xplt 12.1 MB; M2 xplt 3.7 MB; C2 xplt 1.8 MB (Mon_iface/inlet/outlet are shipped as `.gz`) |
| `b4_residual_check.py`, `b4_pressure_attribution.py`, `b3_reaction_force_C1b.py`, `b4_modal_scaling.py` part 2, `stab_analysis.py`, `compare_runs.py`, `compare_stage4.py` | `Mon_all.txt` of the fixed-plate C1b runs (stock and experimental) | 60.6 MB each |
| `t5_analysis.py`, `t2_analysis.py`, `t3_analysis.py`, `poiseuille_check.py post` | `Mon_*.txt` of the tiny-channel and Poiseuille runs | 17-25 MB `Mon_all` each |
| `cavity_fsi_check.py`, `cavity_c1_check.py`, `cavity_plate_analyze.py`, `orig_cavity_figures_*` | run directories with Mon / point files | as above |
Set `CAVITY_RUNS=<directory containing the run folders>` for these scripts (default `runs/cavityNext`). The derived outputs of these steps that the findings rely on are shipped in `results/raw/` (JSON, npz, csv), so the numbers can be inspected without the HPC files.
`scripts/b2_profile_compare.py` was not run: it needs a reference p(x) table that does not exist (no solids4foam pressure field was ever available).

## 3. Needs a FEBio rerun or other unavailable resources
| item | needs | cost |
|---|---|---|
| any coupled / fluid-only result | FEBio v4.13 built from the tag with MKL Pardiso, the inputs here, a Slurm-like node | M1 1.0, M2 4.4, M3 28.2 allocated CPU-h (load-dependent); C4 0.9-1.9 each; fluid-only 0.19-0.37 |
| bit-exact stock binary | a clean rebuild and binary comparison (not done) | - |
| experimental-prototype runs | applying the diffs in `experimental/` to the v4.13 tag, building, running with `FEBIO_PPP_ALPHA` | round 1 0.68, round 2 0.046 CPU-h (the compiled element test harnesses are not shipped) |
| solids4foam values | the maintainers' runs (not committed; this package cites them generically) and the manuscript / PR data | - |
| M4 (dx 0.0125) | about 4-5 x M3, more than 120 CPU-h, about 14 GB | not run |
Reruns are **not bitwise reproducible** even with the identical binary (thread-order noise about 1e-11 relative on a tiny model).

## 4. What cannot be reproduced from anything retained
The solids4foam floor pressure and viscous traction (never available); the plate reaction forces (FEBio returned exact zeros); the solver-applied viscous traction (Gauss-point gradients not stored); the compiled tiny-model element test harness outputs beyond the retained text.
