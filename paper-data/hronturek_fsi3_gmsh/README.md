# Turek–Hron FSI3 on the Gmsh unstructured fluid family (levels 1–3)

These are the FSI3 results of the paper (Section 5.3.1, Table `tab:hronturek_gmsh`). They replace
the ring+box block-structured family (`hronturek_fsi3_4x_findings.md`, `hronturek_production_*`),
which stays in paper-data but is no longer cited as the result. The same runs are the IQN-ILS
rows of the companion solver paper (Section 6.1), which uses a different window (t = 5–6 s) and
its own statistics driver.

## Files

| File | Content |
|---|---|
| `g123.csv`, `g123.json` | Output of `ht_coupled_analysis.py` on the three levels (QoIs, per-period spreads, flag-only forces, coupling statistics). Lengths in m, forces in N/m. |
| `orders.csv` | Successive changes, observed orders, provisional GCIs and differences from Failer & Richter (2020) and Featflow, from `compute_orders.py`. |
| `coupling.csv` | Coupling iterations per step, steps at the iteration limit and the residual r1 at acceptance per level. |
| `window_sensitivity.csv` | The five-period window slid back period by period (window end 4.55–5.84 s), per level. |
| `compute_orders.py` | `python3 compute_orders.py g123.csv > orders.csv` (standard library only). |
| `sensitivity_L1/` | Level-1 coupling-tolerance (1e-7) and time-step (1e-4 s) checks; own README. |

## Analysis command

```bash
python3 ~/claude-runs/fsi3ho/scripts/ht_coupled_analysis.py \
  --case gmsh1=/Users/philipc/claude-runs/fsi3gmsh/results/g1 \
  --case gmsh2=/Users/philipc/claude-runs/fsi3gmsh/results/g2 \
  --case gmsh3=/Volumes/OpenFOAM/philipc-v2512/claude-runs/fsi3-l3-shared \
  --t-end 6 --json g123.json --csv g123.csv
# window_sensitivity.csv: the same command with --drift-from 3.5
```

Rerun on 25 Nov 2026. The CSV is identical, apart from the path column, to the `g123.csv` the
solver-paper agent wrote next to the level-3 data, and the level-1 and level-2 rows are identical
to the earlier `fsi3gmsh/results/g12.csv`.

QoIs: the benchmark definitions over the last five full u_y periods before t = 6 s (level 3:
4.919–5.833 s), mean (max+min)/2 and amplitude (max−min)/2 per period, averaged; frequency =
periods/duration; drag and lift on cylinder + flag per unit depth.

## Results

| QoI | L1 | L2 | L3 | L1→L2 | L2→L3 | order | GCI (prov.) | L3 vs F&R | L3 vs Featflow |
|---|---|---|---|---|---|---|---|---|---|
| u_x (max+min)/2 (mm) | −2.8898 | −2.9411 | −2.9470 | −1.77 % | −0.20 % | 6.2 | 0.03 % | −1.97 % | −2.33 % |
| u_x amplitude (mm) | 2.7370 | 2.7821 | 2.7860 | +1.65 % | +0.14 % | 7.0 | 0.02 % | +2.09 % | +2.43 % |
| u_y amplitude (mm) | 34.982 | 35.239 | 35.327 | +0.74 % | +0.25 % | 3.1 | 0.16 % | +0.85 % | +0.96 % |
| frequency (Hz) | 5.4714 | 5.4737 | 5.4746 | +0.04 % | +0.02 % | 3.0 | 0.01 % | – | +0.27 % |
| drag (max+min)/2 (N/m) | 460.99 | 461.17 | 461.28 | +0.04 % | +0.02 % | 1.5 | 0.04 % | +0.17 % | +0.17 % |
| drag amplitude (N/m) | 28.427 | 28.285 | 28.335 | −0.50 % | +0.18 % | – | – | +1.52 % | +2.14 % |
| lift amplitude (N/m) | 163.69 | 160.90 | 159.09 | −1.70 % | −1.13 % | 1.2 | 2.64 % | +1.59 % | +3.37 % |

- Effective refinement ratio r = (387 219 / 96 712)^(1/4) = 1.415 (per level 1.417 and 1.412);
  Δt is reduced by √2 per level, so the study is scaled (combined space–time).
- Observed order only where monotone with a decreasing change. None is close to 2, so per the
  paper's Section 4 rules no Richardson value is given and the GCI (F_s = 1.25, observed order)
  is a provisional indicator. The drag amplitude is not monotone.
- References: Failer & Richter (2020) level 6, Δt = 0.001 s: u_x −2.890 ± 2.729 mm, u_y amplitude
  35.03 mm, drag 460.5 ± 27.91 N/m, lift amplitude 156.6 N/m (no frequency given). Featflow level 4
  (Δt = 0.00025 s): u_x −2.88 ± 2.72 mm, u_y amplitude 34.99 mm at 5.46 Hz, drag 460.5 ± 27.74 N/m,
  lift amplitude 153.91 N/m.
- The displacements converge from below and cross the Failer & Richter values between levels 1
  and 2 (level 1: −0.14 % u_y amplitude, +0.01 % u_x mean); the lift amplitude converges towards
  them from above.
- Stationarity: sliding the window back by up to three periods changes the level-3 values by at
  most 0.03 % (u_x mean 0.017 %, u_x amplitude 0.018 %, u_y amplitude 0.022 %, frequency 0.009 %,
  drag mean 0.007 %, drag amplitude 0.027 %, lift amplitude 0.014 %). This is a seventh or less of
  the L2→L3 changes, except for the frequency and the drag mean (L2→L3 0.016 % and 0.024 %, about
  two and four times the variation), whose orders are therefore poorly determined. For windows ending between 5.1 and 5.8 s the
  L2→L3 differences of the u_y and lift amplitudes vary by about 5 %, that of the u_x mean by
  about 15 %.
- Force noise: a 4 ms moving average lowers the extrema-based lift amplitude by 0.13–0.14 % on
  every level (`lift_amplitude_smoothed`), so the noise does not grow as Δt falls.

## Set-up (verified against the case files)

- **Build:** solids4foam branch `verification/fsi3-gmsh-partitioned`, commit `9e6d1bfb2`
  (`run_provenance.txt` of level 3: `9e6d1bfb270345d08dfac2a2ecc65ae84e4b0a25` for all five
  segments). It is the quasi-monolithic branch `9dd35dce1` with development `206437983` (after
  #671; contains #546) merged, and is needed for the PETSc SNES path of
  `meshMotionTotalDisplacement`. OpenFOAM v2512 (`_87ed40d256-20251219`), `linux64GccDPInt32Opt`.
  Levels 1–2: UCD Sonic (GCC, `Opt`). Level 3: MeluXina (GCC 13.3 toolchain; the optimisation
  flags of this build were not recorded).
- **Fluid meshes:** Gmsh 4.11.1, quad-dominant (Blossom recombination), solver-paper generator
  `ht_mesh_gmsh2.py --level L --bl-y1 3e-4 --tip-y1 3e-4` (paper-monolithic-fsi `e99e2c3`,
  `simulationData/fsi3-gmsh-meshes`; copy in `~/claude-runs/fsi3gmsh/gen`). The size field
  h = s (h0 + 0.06 d), capped at 3.5 mm in the wake band and 8 mm elsewhere, is scaled by
  s = 2^(−(L−1)/2); structured layer collars around the cylinder and flag (first cell 0.300, 0.212,
  0.150 mm) and along the channel walls.
  - The first cell of the flag-tip corner blocks is capped at 0.300, 0.300 and 0.286 mm, so the
    family is not geometrically similar there.
  - Cells: 96 712, 194 183, 387 219. Flag faces along each long side: 430, 608, 860.
  - Max non-orthogonality 28.5°, 38.8°, 44.6°.
  - md5 of `points`: 13890c5f…, a2eac580…, 8a74309d….
- **Solid meshes:** two-block plate, 62+368 × 24, 88+520 × 34, 124+736 × 48 cells (430 × 24,
  608 × 34, 860 × 48; 10 320, 20 672, 41 280 cells), matching the fluid spacing along the flag.
  St Venant–Kirchhoff, E = 5.6 MPa, ν = 0.4, ρ = 1000 kg/m³, plane strain, total Lagrangian, PETSc
  SNES (rtol 1e-9, stol 1e-12).
- **Time step:** 2e-4, 1.41421e-4, 1e-4 s, constant from t = 0 (rigid flag to t = 2 s, coupled
  from t = 2 s), BDF2, to t = 6.2 s.
- **Fluid:** `div(phi,U) Gauss linearUpwind grad(U)` (unlimited, `grad(U)` least squares),
  `interpolate(U) skewCorrected linear`, PIMPLE 3/3/1, absolute tolerances 1e-9 (p, U),
  `movingWallPressure` on the flag (acceleration term only), corrected wall-normal velocity
  derivative (#671, continuity-based), impulsive parabolic inflow, 2 m/s mean.
- **Mesh motion:** `solids4foamSolidModel` with `meshMotionTotalDisplacement` (total displacement
  from the initial points, so the motion is a function of the current interface position);
  `diffusionElastic` (mu 1, K 3), stiffness scaled by `inverseDistance 2(plate cylinder)`;
  PETSc SNES.
- **Coupling:** IQN-ILS, `outerCorrTolerance 1e-6`, `nOuterCorr 50`, `couplingReuse 5`,
  `requireAllResidualMeasures yes` (strict test), predictor on, relaxation 0.05, AMI transfer.
- Case generator: `~/claude-runs/fsi3gmsh/mk_corrected_case.py` + `apply_shared.py` +
  `apply_shared_fsi.py`; deviations from the earlier corrected set-up in
  `~/claude-runs/fsi3gmsh/README.md`.

## Coupling (`coupling.csv`)

Iterations per coupled step (t = 2–6 s) 5.05, 5.44, 6.01; over the QoI window 4.87, 5.23, 5.78;
maximum 20, 21, 24; no step reached the limit of 50; r1 at acceptance median 5.7e-7, 6.0e-7,
6.3e-7, maximum below 1e-6 on every level (88 284 coupled steps). On level 1 `fsiResiduals.dat`
holds only t > 6 s, so its counts are from the log:

```bash
LC_ALL=C awk '/^Time = /{t=$3+0} /^FSI relative residual1 norm/{n[t]++}
  END{for(k in n) if(k>2 && k<=6.0000001){s+=n[k]-1;c++} print c, s/c}' log.solids4Foam
```

(the same count gives 5.44 on level 2, equal to the `fsiResiduals.dat` value). Level 3 was
restarted four times; `fsiResiduals.dat.upto_5.2424` was rebuilt from the logs for t < 5.2424
(iterations = residual lines − 1, checked against the current segment).

## Runs and cost (t = 0–6.2 s)

| Level | Platform | Ranks | Wall time |
|---|---|---|---|
| 1 | UCD Sonic, one node | 24 | 9.7 h (ClockTime 34 959 s) |
| 2 | UCD Sonic, one node | 32 | 25.4 h (ClockTime 91 428 s) |
| 3 | MeluXina, one exclusive node | 128 | ≈190 node-hours, five restart-safe jobs; ≈16 s per step in the late coupled phase |

## Paths

- Level 1–2 results: `~/claude-runs/fsi3gmsh/results/g1`, `g2` (logs, postProcessing); cases in
  `~/claude-runs/fsi3gmsh/cases/g1_np24`, `g2_np32` (Sonic `/scratch/pcardiff/fsi3gmsh/cases`).
- Level 3: `/Volumes/OpenFOAM/philipc-v2512/claude-runs/fsi3-l3-shared/` (postProcessing, logs,
  system, fsiProperties, `run_provenance.txt`); MeluXina
  `/project/home/p201287/philip/fsi3-l3/g3_np128` (case files, postProcessing, logs).
- Solver-paper report on level 3: paper-monolithic-fsi issue 13, comment 7176633587.

## Coupling-tolerance and time-step checks (`sensitivity_L1/`)

On level 1, a coupling tolerance of 1e-7 (job 1191796) and a halved time step Δt = 1e-4 s (job
1192165) against the baseline: every QoI changes by ≤ 0.0022 % with the tolerance, and by
0.047 % (drag amplitude) and ≤ 0.008 % (all others) with the time step, against L1→L2 changes of
0.04–1.8 %. These replace, in the paper, the coupling studies of the uniform block-structured
family. Details in `sensitivity_L1/README.md`.

## Not done on this family

No time-step study on the finer levels (Δt is scaled with h), no coupling-algorithm comparison,
no solid-only refinement and no second-code solution. The second-solid-code run on the ring+box
1x′ mesh (`fsi3_precice_calculix/`) is on another mesh.

## History figure (`figure/`, Fig. `fig:fsi3_histories`)

`figures/fsi3_gmsh_histories.pdf`: tip u_y, u_x, lift and drag (cylinder + flag, per unit depth)
of L1–L3 over the five-period QoI window, with enlargements of one u_y maximum, one u_x minimum
and one lift maximum, against Failer & Richter (2020) level 6 and Featflow level 4 (mean ±
amplitude; the u_y and lift reference lines are centred on the L3 (max+min)/2, since only the
amplitudes are compared). Each level is shifted so that its last upward u_y crossing before
t = 6 s (the window end: 5.8399, 5.8329, 5.8325 s) is at t − t0 = 0.

| File | Content |
|---|---|
| `figure/extract_histories.py` | Stitches the restart segments of postProcessing (last row of each time, as `ht_coupled_analysis.py`), divides the forces by the depth 0.015 m and writes t = 4.5–6 s. With `--check`, recomputes the window statistics from the written CSVs and compares them with `g123.csv`. |
| `figure/data/fsi3_gmsh_L{1,2,3}_histories.csv` | Columns `t_s, ux_m, uy_m, drag_N_per_m, lift_N_per_m`, every time step (not decimated), t = 4.5–6 s; 0.47, 0.66 and 0.93 MB. |
| `figure/fsi3_histories_figure.py` | Draws the figure from the CSVs, in the shared style (`../figstyle`). |

```bash
cd paper-data/hronturek_fsi3_gmsh/figure
python3 extract_histories.py --out-dir data --check ../g123.csv \
  --case L1=<results>/g1 --case L2=<results>/g2 --case L3=<fsi3-l3-shared>
python3 fsi3_histories_figure.py --data-dir data --out ../../../figures/fsi3_gmsh_histories.pdf
```

The inputs were `~/claude-runs/fsi3gmsh/results/g1`, `g2` and
`/Volumes/OpenFOAM/philipc-v2512/claude-runs/fsi3-l3-shared`. The check reproduces every QoI of
`g123.csv` (frequency, (max+min)/2 and amplitude of u_x, u_y, drag and lift) to a relative
difference of 5e-7 or less (the lift (max+min)/2, a near-zero mean; all others ≤ 5e-8), so the
CSVs carry the same data as the table.
