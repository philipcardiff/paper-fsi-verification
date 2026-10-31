# Reproduction record

## A. Binary and source identity (checked when this package was prepared)
- FEBio binary sha256 `f5ae1ed6f347f875c392f93d847f2369f437345a9b785296d46b0bd6fbb72f02`: re-computed from the build, and recorded by all 25 runs that wrote `run_meta_start.json` (`run_meta_sanitized.csv`: one distinct value).
- Source: FEBio tag v4.13, commit `32ae206ff4881dfb54f62296cd1558e58ed9fcc6`, working tree clean (no local changes, no patch, no plugin); every solver log reports `version 4.13.0.32ae206ff`. The binary is linked against MKL (user-space libraries) for PARDISO; the build tree itself is not part of the package.
- Vendored reference: `reference_hashes.csv`; the files fetched read-only from the two draft PR heads (solids4foam PR #637 head 30ae00eb, manuscript PR #25 head 97f9efef) have the recorded hashes.
- Pre-registration hash: `../docs/PREREGISTRATION.md`.

## B. How the simulations were run (not re-run here)
Per job (Slurm, shared partition, `--cpus-per-task` 2, 4 for RT-C): in the run directory, with the MKL libraries on the library path, `OMP_NUM_THREADS = MKL_NUM_THREADS = cores`, `MKL_DYNAMIC=FALSE`, `OMP_PROC_BIND=close`, `OMP_PLACES=cores`:
`/usr/bin/time -v -o time.txt febio4 -i <TAG>.feb -silent`.
Inputs are produced by `scripts/make_rtr_feb.py` with the arguments recorded in `results/raw/<TAG>/model_meta.json` (`scripts/regen_inputs.py` maps them to command-line options). Examples:
```
python3 -I scripts/make_rtr_feb.py --mode fluid --mesh M1 --dt 0.01 --t-end 2 --rhoi 1 --out RTB_fluid_M1 --tag RTB_fluid_M1
python3 -I scripts/make_rtr_feb.py --mode twist --mesh M3 --dt 0.1 --t-end 2 --amp 0.5 --out RTA_twist_M3 --tag RTA_twist_M3
python3 -I scripts/make_rtr_feb.py --mode fsi --mesh M1 --dt 0.01 --t-end 2 --rhoi 0.9 --max-refs 25 --out RTC_fsi_M1 --tag RTC_fsi_M1
```
Each run's exit code was never taken as evidence: normal termination, step count and fields were checked from the log (`../results/run_summary.csv`).

## C. Non-simulation reproduction check run for this package
`./reproduce.sh` run in a scratch copy of the package on a one-core login node (Python 3.10, numpy 1.21.5, scipy 1.8.0, matplotlib 3.5.1); wall time 4 min 10 s; result **REPRODUCTION CHECKS PASSED**:
1. Vendored reference hashes (3 files): match.
2. Inputs regenerated from the shipped generator: 13 of 19 identical to the executed inputs by sha256 (all RT-A, RT-B, RT-C inputs, and the smoke inputs SMK_fluid_M0, SMK_ale_M0, SMK_fsi_M0_mu0p3, SMK_fsi_M0_rhoi0p9); 4 smoke inputs (SMK_fsi_M0, _cm, _svk, _tf0) identical modulo the order of NodeSet lines; 2 superseded first designs (RTA_law1, RTA_law1b) are not regenerable by the final generator and are shipped as run. Failed attempts (`inputs/failed_attempts/`) are shipped as run.
3. Campaign analysis scripts re-run on the stored raw outputs: 10 result JSON files identical to `results/campaign_analysis/` (relative difference below 2e-19).
4. Independent re-derivation (`scripts/independent_check.py`): `results/independent_check.json` identical (5,004 entries). Its values against the campaign values: law identity 1.2246712e-12 (both), rotation 7.014e-13 m (both), T_A errors identical (same reference function), RT-B as-run and aligned errors agree to 3-4 digits (differences only in T3, see INCONSISTENCIES.md item 11), wall velocity equals the table at t - dt/2 to 2.5e-13 rad/s.
5. Newton and termination statistics recomputed from the solver logs: three CSV files identical.
6. Tables identical, five figures regenerated.
7. Spectral fixed-geometry reference regenerated from the omega table (Chebyshev N = 40, Radau 1e-12): arrays T_B, H_f, L, p identical to the cached `fluid_exact_N40_dt0.01*.npz` (difference 0).
8. Analytical coupled estimate re-run: identical to `results/predicted_coupled_element_mean_AE.json`.
9. `SHA256SUMS`: all files match.

## D. What a reader cannot check from this package
The FEBio runs themselves, the solver binary, the source-code claims against the FEBio tag (the quoted lines are listed in `../docs/BC_TIME_CONVENTION.md` and `../docs/EQUIVALENCE.md` and can be checked against tag v4.13), and the draft PRs' own content beyond the hashes recorded in `reference_hashes.csv`.
