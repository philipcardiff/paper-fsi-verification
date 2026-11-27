# Turek–Hron FSI3, Gmsh level 1: coupling-tolerance and time-step checks

Evidence for the "FSI3 time step and coupling" paragraph of Section 5.3.1 (`sec:fsi23`). These runs
replace the coupling-tolerance and time-step studies of the earlier uniform block-structured
(ring+box) family, so every FSI3 result of the paper is on the Gmsh family.

## Runs

All three runs use the level-1 set-up of `../README.md` (96 712 fluid and 10 320 solid cells;
IQN-ILS, `outerCorrTolerance 1e-6`, `nOuterCorr 50`, `couplingReuse 5`,
`requireAllResidualMeasures yes`; BDF2; rigid flag to t = 2 s, coupled to t = 6.2 s), with one
change each:

| Run | Change from the baseline | Job (UCD Sonic, 24 ranks) | Wall time |
|---|---|---|---|
| `g1` (baseline) | – (the L1 production run, `../g123.csv` row `gmsh1`) | – | 9.7 h |
| `g1_tol1e-7` | `outerCorrTolerance 1e-7` | 1191796, t = 0–6.2 s, exit 0 | 12.7 h |
| `g1_dt1e-4` | Δt = 1e-4 s (`writeInterval 2500`, 0.25-s checkpoints kept) | 1192165, restart at t = 2.0 s to 6.2 s, exit 0 | 20.0 h + 0.8 h (rigid phase) |

- **Build:** solids4foam `9e6d1bfb2` (`9e6d1bfb270345d08dfac2a2ecc65ae84e4b0a25`), the build of
  levels 1–2, in `/scratch/pcardiff/OpenFOAM/gmshpart-v2512` (OpenFOAM v2512,
  `linux64GccDPInt32Opt`), UCD Sonic (`runs/*/run_provenance.txt`).
- **Start-up of `g1_dt1e-4`** (`runs/g1_dt1e-4/sensitivity.txt`): the first two attempts (jobs
  1191797 and 1192060) aborted at t = 2.0001 s, the first coupled step, in coupling iteration 2.
  The second coupling iterate loads the flag with about 5.6 kN/m (0.49 kN/m in the baseline at the
  same stage), the added-mass response to the first interface update, which scales with 1/Δt².
  The solid Newton solve failed under it (DIVERGED_DTOL with the full step; DIVERGED_LINEAR_SOLVE
  with a backtracking line search). Job 1192165 restarted from the t = 2.0 checkpoint of attempt 2
  (rigid phase, identical set-up) with the solid `snes_linesearch_type bt` (was `basic`) and
  `stopOnPetscError no`; tolerances unchanged. The first coupled step then took a slightly
  different path and converged in 24 iterations.
  - The line search alone did not pass the step (attempt 2 had it); the restart, which perturbed
    the first iterates, did. A softer start of the coupling (e.g. a smaller initial relaxation
    factor) would be the robust fix.
  - Audit of the whole log (`runs/g1_dt1e-4/audit_solid_failures.{py,txt}`): 0 solid SNES/KSP
    failures in the 42 000 coupled steps, and no PETSc solve of any kind reported a failure (0
    "did not converge" lines against 518 704 converged nonlinear solves). `stopOnPetscError no`
    therefore never acted; every coupled step was accepted below 1e-6 and none reached the limit.
- Case data: Sonic `/scratch/pcardiff/fsi3gmsh_l1sens/cases`; postProcessing and logs on the Mac
  in `~/claude-runs/fsi3gmsh/results/g1_{tol1e-7,dt1e-4}`.

## Analysis

```bash
python3 ~/claude-runs/fsi3ho/scripts/ht_coupled_analysis.py \
  --case g1=~/claude-runs/fsi3gmsh/l1sens_base_g1 \
  --case g1_tol1e-7=~/claude-runs/fsi3gmsh/results/g1_tol1e-7 \
  --case g1_dt1e-4=~/claude-runs/fsi3gmsh/results/g1_dt1e-4 \
  --t-end 6 --json g1_sens.json --csv g1_sens.csv
python3 compute_changes.py g1_sens.csv ../g123.csv > changes.csv
```

The baseline is `results/g1` with its log decompressed (`l1sens_base_g1`); its QoIs equal the
`gmsh1` row of `../g123.csv` exactly (checked by `compute_changes.py`). QoIs as in `../README.md`
(last five full u_y periods before t = 6 s). The baseline's `fsiResiduals.dat` holds only
t > 6 s, so its iteration counts in `coupling.csv` are from the log, as in `../coupling.csv`.

| File | Content |
|---|---|
| `g1_sens.csv`, `g1_sens.json` | Output of `ht_coupled_analysis.py` on the three runs (m, N/m, Hz). |
| `changes.csv` | QoIs, relative changes from the baseline (to its magnitude, as in the paper tables) and the L1→L2 change. |
| `coupling.csv` | Coupling iterations per coupled step (t = 2–6 s), steps at the limit, r1 at acceptance. |
| `compute_changes.py` | Writes `changes.csv` (standard library only). |
| `runs/` | `run_provenance.txt`, `sensitivity.txt`, the Slurm output of job 1192165 and the solid-failure audit. |

## Results

| QoI | baseline (L1) | tol 1e-7 | Δ tol | Δt 1e-4 | Δ Δt | L1→L2 |
|---|---|---|---|---|---|---|
| u_x (max+min)/2 (mm) | −2.8898 | −2.8898 | −0.0005 % | −2.8897 | +0.0054 % | −1.77 % |
| u_x amplitude (mm) | 2.7370 | 2.7370 | +0.0004 % | 2.7370 | +0.0012 % | +1.65 % |
| u_y amplitude (mm) | 34.982 | 34.982 | +0.0002 % | 34.979 | −0.0066 % | +0.74 % |
| frequency (Hz) | 5.4714 | 5.4714 | −0.0000 % | 5.4714 | +0.0003 % | +0.04 % |
| drag (max+min)/2 (N/m) | 460.99 | 460.99 | −0.0001 % | 461.01 | +0.0052 % | +0.04 % |
| drag amplitude (N/m) | 28.427 | 28.427 | +0.0019 % | 28.440 | +0.047 % | −0.50 % |
| lift amplitude (N/m) | 163.69 | 163.69 | −0.0022 % | 163.70 | +0.0076 % | −1.70 % |
| coupling its/step (max) | 5.05 (20) | 6.39 (50) | | 5.69 (24) | | |
| median r1 at acceptance | 5.7e-7 | 5.1e-8 | | 5.9e-7 | | |

- **Coupling tolerance:** every QoI changes by ≤ 0.0022 % (lift amplitude), 18 times less than the
  smallest L1→L2 change and 200–800 times less than the others. The tip histories differ by at
  most 9.5e-7 m in u_y and 1.6e-7 m in u_x over t = 2–6 s. Two of the 20 000 coupled steps
  reached the iteration limit and were accepted at r1 = 1.2e-7 and 1.4e-7.
- **Time step:** the drag amplitude changes by 0.047 %, a tenth of its (non-monotone) L1→L2
  change; every other QoI by ≤ 0.008 %, ≥ 90 times less than the L1→L2 changes of the
  displacements and the lift amplitude and 7–140 times less than those of the frequency and drag
  mean. Δt falls by only √2 per level, so the temporal part of a level-to-level change is smaller
  than the change on halving Δt for any order: the level-to-level changes are essentially
  spatial. The tip histories differ by up to 4e-4 m in u_y, a phase offset of about 0.2 ms
  carried over from the different start of the coupled phase; it does not enter the per-period
  QoIs.
