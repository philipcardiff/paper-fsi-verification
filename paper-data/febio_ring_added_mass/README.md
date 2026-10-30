# Independent FEBio evidence for ringAddedMass: reproducibility package

Status: **independent FEBio v4.13 evidence, draft; not a Category A verification claim.** Companion to the solids4foam `ringAddedMass` verification (commit 1576edcc, see `reference/`).

## Layout
- `tools/` generator (`make_ring_feb.py`), run script (`run_febio_job.sh`, Slurm), post-processing (`phaseD_postprocess.sh`, `run_report.py`, `pressure_modes.py`, `interface_leakage.py`, `convergence_margins.py`, `dilatation_floor.py`, `iter_compare.py`) and table/figure builders (`phaseC_results.py`, `phaseD_results.py`).
- `runs/<phase>/<case>/` for every case used: the `.feb` input (only where < 0.8 MB; larger mesh inputs are regenerated and hash-checked), `model_meta.json` (all parameters, generator args, input sha256), `report.json`, `pressure_modes.json`, `interface_leakage.json`, convergence margins, `probes.txt`, `time.txt` (wall time, RSS), `run_meta_*.json` (FEBio binary sha256, host, exit code), logs (gz). No `.xplt` or nodal dumps (large).
- `results/` summary CSVs (`phaseD_febio_summary.csv`, `phaseD_s4f_vs_febio.csv`, `phaseD_acceptance_checks.csv`, `phaseC_*.csv`, verbatim solids4foam README tables).
- `figures/` fig_D1-D5 (final) and fig_C1-C6 (supporting history). `reference/` official solids4foam reference JSON + README. `docs/` reports. `SHA256SUMS`.

## Final selected results vs supporting history
- **Final (Phase D)**: `runs/phaseD/D1_weak_M1, D4_weak_M2, D5_moderate_M3, D2_strong_M1, D3_strong_M2` with `runs/phaseC/C2_Emesh1e-3_wet_f1/f2` (moderate M1/M2); dry runs `B1_dry_hex8_f1` (M1), `D1_dry_f2` (M2), `C1_M3_dry_f4` (M3). Inviscid-limit formulation: slip outer wall, slip fluid-ring interface with w.n = 0 enforced, E_mesh = 1e-3 Pa, K = 1e8 Pa, rhoi = 1, 100 steps/period, 4 periods.
- **Diagnostic history (kept, not hidden)**: original no-slip runs (Phase B: -6.0 % ratio error M1, -2.7 % M2); wall-slip-only (-1.1 %); interface/outer-wall investigation (`docs/PHASE_C_INTERFACE_DIAGNOSTIC.md`); E_mesh = 1 Pa baseline (+0.17..+0.18 % ratio error, mesh-independent offset traced to the regulariser: `docs/PHASE_C_REGULARISER_DIAGNOSTIC.md`); interface leakage and force convergence (`docs/PHASE_D_FINAL_COMPARISON.md` §4); failed tightened-tolerance runs (Phase C tolerance run, D6: `docs/PHASE_D_ITERATIVE_ERROR.md`).
- **FEBio initial-condition defect**: in v4.13, a `math`-expression `initial fluid velocity` IC crashes (segmentation fault). Workaround used for every wet run: one constant-vector IC per fluid node (`--ic-mode pernode`), with the compatible potential-flow field. The pristine v4.13 binary is used throughout; a source patch exists locally and was never used for results. A report for FEBio is being drafted separately.

## Build and solver provenance
FEBio v4.13 built in user space from the tagged source with Intel MKL (Pardiso, non-symmetric); binary sha256 is in each `run_meta_start.json`. Fluid-FSI module, generalized-alpha rhoi = 1, full Newton, `reform_each_time_step`, tolerances dtol = vtol = ftol = 1e-6, etol = 1e-8, rtol = 1e-10, 1 core per run, Slurm partition `main`, host xenosim (shared node).

## Reproduce
1. `./reproduce.sh` regenerates the final inputs from the generator and compares them with the shipped inputs (hash), then rebuilds all CSVs and figures D1-D5 from the stored run outputs. This needs only Python 3 with numpy/matplotlib.
2. Re-running a simulation needs FEBio 4.13 + MKL (`tools/febio_env.sh` shows our paths; edit them), then `tools/run_febio_job.sh <run_dir> <tag>.feb` and `tools/phaseD_postprocess.sh <wet> <dry> <NT> <NR>` (NT, NR = 48f, 12f angular/radial cells for mesh factor f = 1, 2, 4).
Case commands: `python3 tools/make_ring_feb.py --out <dir> --tag <tag> --rho-f <14|140|1400> --k 1e8 --e-mesh 0.001 --rhoi 1 --mesh-factor <1|2|4> --steps-per-period 100 --periods 4 --ic potential --ic-field compatible --ic-mode pernode --wall slip --interface slip`; dry cases omit `--rho-f` (exact arguments are recorded in each `model_meta.json`).

## Known gaps
Iterative error not isolated; time order untested; two mesh levels only for weak/strong; ratio error flat (not decreasing) with refinement; finite-viscosity no-slip problem not reproduced. See `docs/PHASE_D_FINAL_COMPARISON.md`.
