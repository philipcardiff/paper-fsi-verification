# 3-D pressure-pulse tube: axisymmetric (1° wedge) mesh sequence and the averaged c_p

Evidence for paper Section `sec:3dtube` (the *Axisymmetric sequence* bullet, the
*Problem* paragraph on the wall law, the c_p definition and Table `tab:3dtube`),
Appendix E (`tab:refvals_3dtube`) and the tube entries of Section 6.

Two things are here:

1. the axisymmetric model of the tube, computed on four levels (L1–L4) with a
   Hookean wall, as a mesh-only series at fixed Δt, and on four levels (S1–S4)
   with a St Venant–Kirchhoff wall;
2. the pressure-front speed c_p recomputed, for every tube result the paper
   reports, with the averaged definition (below), from the stored axis-pressure
   probe histories.

## Set-up

The physical problem, the solver settings and the QoI definitions are those of
the 3-D case: the `3dTube` tutorial at solids4foam `206437983` with the
blockMeshDicts of `9b52c2303` (`template/3dTube`, without its images), Robin–Neumann
coupling (`fixedRelaxation`, relaxation 1, `elasticWallPressure`, nOuterCorr 100,
tolerance 1e-6), BDF2 in fluid and solid, `pimpleFluid` with linearUpwind
convection, Hookean wall (`linearElastic`, `linearGeometryTotalDisplacement`;
E 0.3 MPa, ν 0.3, ρ 1200 kg/m³), fluid ρ 1000 kg/m³, ν 3e-6 m²/s, inlet
1333.2 Pa to 3 ms, ramped to 0 at 3.1 ms, end time 0.02 s.

What differs from the 3-D quarter tube (`scripts/make_case.py`):

- **Geometry:** a 1° wedge (±0.5°), one cell in θ; axis x, radial y, wedge
  across z (the `womersleyTube` layout; solids4foam treats z as the
  empty/wedge direction of a 2-D solid). The symmetry planes become `wedge`
  patches (`solidWedge` for D), and the fluid has an `empty` axis patch.
- **Interface:** conformal, `interfaceTransferMethod directMap` (3-D: AMI).
- **Solid Jacobian:** assembled compact-stencil Jacobian instead of the 3-D
  matrix-free `snes_mf_operator`, which stalls at ≈1e-4 of the initial
  residual on the wedge. Tightening `snes_stol` from 1e-6 to 1e-8 changes no
  QoI by more than 2e-6 relative (run `L1_stol8`).
- **Monitors:** the nine wall stations of the 3-D driver (z = 5–45 mm), with
  u_r = sign(D_y)·hypot(D_y, D_z) and u_z = D_x; axis-pressure probes at
  r = 0.707 mm, the radius of the 3-D probes at (0.5, 0.5) mm.
- **Wedge-angle control** (`L1_w2`, 2°): every QoI within 1.7e-4 of the 1° run.

Case set-ups as generated: `case/L1` (Hookean) and `case/S1` (StVK); the other
levels differ only in the cell counts, Δt and decomposition
(`runs/<case>/axi_settings.json`).

### Levels and jobs (UCD Sonic)

| Case | Fluid cells (r × z) | Solid cells | Δt (s) | Ranks | Job | Wall clock |
|---|---:|---:|---:|---:|---:|---:|
| L1 | 1 200 (15 × 80) | 320 | 2.5e-5 | 1 | 1195942 | 23 min |
| L2 | 4 800 (30 × 160) | 1 280 | 1.25e-5 | 8 | 1196114 | 32 min |
| L3 | 19 200 (60 × 320) | 5 120 | 6.25e-6 | 6 | 1196117 | 2.9 h |
| L4 | 76 800 (120 × 640) | 20 480 | 3.125e-6 | 12 | 1196118 | 11.7 h |
| L1m, L2m, L4m | as L1, L2, L4 | | 6.25e-6 | 1, 8, 12 | 1200218, 1200219, 1200220 | |
| S1–S4 (StVK) | as L1–L4 | | as L1–L4 | 1, 8, 6, 12 | 1196115, 1196116, 1196120, 1196121 | |
| L1_stol8, L1_w2 (controls) | as L1 | | 2.5e-5 | 1 | 1196112, 1196113 | |

- StVK: `StVenantKirchhoffElastic` with
  `nonLinearGeometryTotalLagrangianTotalDisplacement`, all else identical.
- Every step of every reported run met the Robin convergence criteria; no SNES
  solve failed. Mean coupling iterations per step 10.1, 6.7, 6.0, 5.3 (L1–L4).
- **L1m** develops the BDF2 instability of the 3-D L1 at Δt = 6.25e-6 s:
  exponential growth from ≈9.5–10 ms, by a factor of ≈7 per 0.5 ms, failure at
  12.3 ms. **L2m** stops at 17.5 ms (a coupling step reached nOuterCorr after
  episodes of 70–90 iterations per step from 11.5 ms). Both are clean to
  ≥ 9.5 ms; their QoIs are from histories truncated at 10.3 ms (all point-A
  QoIs and c_p are fixed by ≈8.6 ms; on L1–L3 the truncated extraction
  reproduces the full-run values exactly).
- An L5 (≈ 5 days on the cores available) was submitted as job 1196119 and
  cancelled before it started.
- Build: Sonic `s4f-671bso`, solids4foam `64896b0a5` = `206437983` plus one
  local commit (an outlet condition this case does not use), OpenFOAM v2512,
  GCC 11.5, OpenMPI 5.0.5, PETSc 3.24.1 with hypre (`runs/<case>/log.environment`).
- Raw runs: `sonic:/scratch/pcardiff/tube_axi/runs/<case>` (`log.solids4Foam`,
  `postProcessing/`). Earlier submissions 1195350–1195385 and 1195506–1195518
  were invalid (a job-script error, then a wedge with the axis along z, in which
  solids4foam suppresses the axial wall displacement); they are not used.

## The averaged pressure-front speed c_p

The paper's c_p is the least-squares slope dz/dt, over the stations
z = 10–40 mm, of the first time at which the axis-pressure probe signal rises
through half the inlet pressure, **after a centred moving average of the probe
signal over 0.4 ms** (all samples within ±0.2 ms; two periods of the ≈0.2 ms
spurious axis-pressure oscillation). Otherwise the definition is that of the
3dTube verification driver, whose unaveraged value is kept as "raw".

- `scripts/cp_averaged.py` (standard library only) computes it from an OpenFOAM
  probes file or a Kratos `history.csv`. With window 0 it reproduces the
  driver's raw c_p of every run to all printed digits (e.g. 3-D L1–L3
  4.6639, 4.6511, 4.6253 m/s; Kratos 4.578, 4.583 m/s).
- Probe times are written to 6 significant digits; the script snaps them to the
  time-step grid.
- `results/cp_averaged.csv`: raw and averaged (0.1, 0.2, 0.4 ms) c_p and the
  station arrivals for every case below. `results/osc_amplitude.csv`: the
  oscillation amplitude max|p − p̄₀.₄ₘₛ| within ±1 ms of the arrival.
- Probe histories: `probes/*.p.xz` (wedge runs, 3-D L1/L2 of Table `tab:3dtube`,
  the 3-D Δt runs and the 3-D StVK runs); the 3-D L3 probes are in
  `../3dtube_level3_corrected_ecore/postProcessing/`, the Kratos histories in
  `../kratos_3dtube/runs/`.

Reproduce (from this directory):

```bash
K=../kratos_3dtube/runs
L3=../3dtube_level3_corrected_ecore/postProcessing/verificationAxisPressure/fluid/0/p
args=""; for f in probes/*.p.xz; do args="$args $(basename $f .p.xz)=$f"; done
python3 scripts/cp_averaged.py --window 0 --window 1e-4 --window 2e-4 --window 4e-4 \
    --truncate 0.0103 $args 3d_L3=$L3 \
    kratos_L1=$K/f1_dt2.5e-5_bdf2_amg/history.csv \
    kratos_L2=$K/f2_dt1.25e-5_bdf2_amg/history.csv > results/cp_averaged.csv
```

### Why average

On the axis the pressure rises slowly (≈0.25 Pa/µs) and carries an oscillation
of ≈0.2 ms period. Its amplitude near the front grows with refinement, as the
numerical damping falls, and mainly with Δt: max|p − p̄| at z = 10–40 mm is
2–24 Pa on L1, 10–40 Pa on 3-D L3 and 13–37 Pa on wedge L4 (wedge L4m, the L4 mesh at the
L3 Δt: 8–37 Pa, as L3). The oscillation shifts the half-inlet crossing, so the
raw c_p drifts with the mesh (wedge −0.39, −0.54, −0.81 % per level) and with
Δt (3-D L3, Δt halved: −0.86 %). The 0.4 ms average removes it.

### Averaged c_p (m/s), all tube results

| Result | raw | 0.2 ms | **0.4 ms** |
|---|---:|---:|---:|
| 3-D L1, L2, L3 (Table `tab:3dtube`) | 4.6639, 4.6511, 4.6253 | 4.6635, 4.6450, 4.6385 | **4.6616, 4.6418, 4.6362** |
| 3-D L2, Δt 1.25e-5 (Sonic control) → 6.25e-6 | 4.6524 → 4.6490 | | **4.6435 → 4.6429** (−0.013 %) |
| 3-D L3, Δt 6.25e-6 → 3.125e-6 | 4.6253 → 4.5854 (−0.86 %) | | **4.6362 → 4.6334** (−0.061 %) |
| 3-D StVK L1, L2, L3 | 4.6954, 4.6838, 4.6683 | | **4.6944, 4.6765, 4.6689** |
| 3-D Hookean L1 (Sonic control of the StVK runs) | 4.6639 | | 4.6615 |
| Wedge L1–L4 | 4.6700, 4.6517, 4.6267, 4.5894 | 4.6696, 4.6458, 4.6351, 4.6352 | **4.6677, 4.6427, 4.6328, 4.6357** |
| Wedge L1m, L2m, L3, L4m (Δt 6.25e-6) | 4.6982, 4.6512, 4.6267, 4.5918 | | **4.6683, 4.6452, 4.6328, 4.6280** |
| Wedge S1–S4 (StVK) | 4.7014, 4.6817, 4.6659, 4.6259 | | **4.7004, 4.6743, 4.6648, 4.6593** |
| Kratos L1, L2 (BDF2) | 4.5783, 4.5834 | | **4.5770, 4.5812** |

- **3-D:** −0.42 %, then −0.12 % (observed order 1.8 from L1–L3).
- **Wedge:** −0.54 %, −0.21 %, +0.06 % (L3→L4); at fixed Δt −0.49 %, −0.27 %,
  −0.10 %. Halving Δt on the wedge L4 mesh changes it by +0.17 %; on 3-D L2 and
  L3 by −0.01 % and −0.06 %; the platform control on 3-D L2 (Sonic against
  macOS) differs by 0.04 %.
- **Window:** 0.2 against 0.4 ms changes it by ≤ 0.05 % on L3 and L4 (3-D L3
  0.05 %, wedge L3 0.05 %, wedge L4 0.01 %, 3-D L3 at halved Δt 0.03 %); any
  window from 0.35 to 0.6 ms stays within 0.1 % of the 0.4 ms value.
- **Conclusion:** the averaged c_p is converged to about 0.1–0.2 %, with no
  Richardson estimate: the wedge sequence is not monotone at L4, and Δt still
  moves it by 0.17 % there. Its time-step sensitivity on the 3-D levels
  (≤ 0.06 %) replaces the ≈1 % of the raw definition.

## Results (Hookean wall, Δt halved with h)

Point A at z = 25 mm; definitions of Table `tab:3dtube`. Changes are
(f_{k+1} − f_k)/|f_k|; p = log2((f_k − f_{k+1})/(f_{k+1} − f_{k+2})).
Full precision in `results/qoi_hooke_full.csv`, `analysis_hooke.md`.

| Quantity | L1 | L2 | L3 | L4 | L1→L2 | L2→L3 | L3→L4 | p (L2–L4) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| u_r,max (mm) | 0.159721 | 0.159509 | 0.158366 | 0.158264 | −0.13 % | −0.72 % | −0.06 % | 3.5 (above nominal; no extrapolation) |
| u_z,min (mm) | −0.088482 | −0.086885 | −0.085788 | −0.085542 | +1.81 % | +1.26 % | +0.29 % | **2.15**, f∞ = −0.085470 |
| t_arr (ms) | 5.8705 | 5.8416 | 5.8231 | 5.8319 | −0.49 % | −0.32 % | +0.15 % | – |
| c_p, averaged (m/s) | 4.6677 | 4.6427 | 4.6328 | 4.6357 | −0.54 % | −0.21 % | +0.06 % | – |
| c_p, raw (m/s) | 4.6700 | 4.6517 | 4.6267 | 4.5894 | −0.39 % | −0.54 % | −0.81 % | – |
| c_wall (m/s) | 5.0693 | 5.1394 | 5.0681 | 5.0642 | +1.38 % | −1.39 % | −0.08 % | – |

- The L2–L4 triple of u_r,max gives p = 3.5, far above the formal order, so by
  the paper's rules it is pre-asymptotic and no Richardson value is used; its
  triple limit (0.158254 mm) lies 0.006 % from L4.
- Fixed Δt (`qoi_fixeddt_t103.csv`): u_r,max −0.44, −0.61, −0.03 %; u_z,min
  +1.89, +1.24, +0.28 % (p = 2.17). The pre-asymptotic behaviour of u_r,max is
  therefore spatial.

### Wedge against the 3-D quarter tube (matched levels)

3-D − wedge, relative to the wedge: u_r,max −0.09, −0.17, −0.04 %; u_z,min +0.35,
+0.36, +0.11 %; t_arr +0.02, −0.00, +0.01 %; averaged c_p −0.13, −0.02, +0.07 %;
c_wall −0.04, −0.02, −0.02 % (L1, L2, L3). The 3-D L3 u_r,max lies 0.03 % above
the wedge L4 value, and the 3-D L3 u_z,min 0.26 % below (deeper than) the wedge
Richardson limit −0.08547 mm.

## St Venant–Kirchhoff wall and the law effect

`results/qoi_stvk_full.csv` (S1–S4 to 20 ms), `qoi_stvk_t103.csv` (truncated, as
first reported; identical point-A QoIs), and `results/3d_stvk/` (3-D StVK L1–L3
and the Sonic Hookean L1 control: driver CSVs, jobs 1192830, 1195242, 1195243,
1192829; raw runs in `sonic:/scratch/pcardiff/tube_stvk/`).

Law effect (StVK − Hookean)/|Hookean| on the same level:

| Quantity | 3-D L1 | 3-D L3 | Wedge L1 | L2 | L3 | **L4** |
|---|---:|---:|---:|---:|---:|---:|
| u_r,max | −1.11 % | −0.50 % | −1.11 % | −0.56 % | −0.52 % | **−0.42 %** |
| u_z,min | −3.19 % | −3.40 % | −3.20 % | −3.29 % | −3.40 % | **−3.47 %** |
| t_arr | −0.81 % | −0.63 % | −0.81 % | −0.64 % | −0.64 % | **−0.77 %** |
| c_p, averaged | +0.70 % | +0.70 % | +0.70 % | +0.68 % | +0.69 % | **+0.51 %** |
| c_p, raw | +0.68 % | +0.93 % | +0.67 % | +0.65 % | +0.85 % | +0.80 % |
| c_wall | +0.04 % | +0.30 % | +0.07 % | −1.21 % | +0.30 % | +0.32 % |

- On L1 the u_r,max effect (−1.1 %) mixes the law effect with the larger L1
  mesh error; on the finest levels it is ≈ −0.4 to −0.5 %.
- Peak against Tuković et al. (2018, StVK, 0.15694 mm, 7.23 ms): 3-D StVK L3
  0.15751 mm (+0.36 %, at 7.20 ms), against +0.87 % for the Hookean 3-D L3;
  wedge S4 0.157602 mm (+0.42 %).

## Files

- `scripts/`: `make_case.py` (case generator from `template/3dTube`),
  `run_case.sh` (Slurm job), `post_axi.py` (QoIs with the 3-D driver's
  functions; `--truncate`), `analyse.py` (orders, limits, 3-D comparison),
  `plot_histories.py`, `cp_averaged.py`, `osc_amplitude.py`, and `driver/` (the
  3dTube verification driver at `206437983` with the monitor-time snapping
  patch, and its reference JSON).
- `case/L1`, `case/S1`: generated case set-ups (no mesh).
- `runs/<case>/`: `axi_settings.json`, `log.environment`, `log.run`, `qoi.csv`.
- `results/`: QoI CSVs, `analysis_hooke.md`, `cp_averaged.csv`,
  `osc_amplitude.csv`, point-A histories (`histories/`), and
  `tube_axi_vs_3d_pointA.png` (wedge L1–L4 against 3-D L1–L3).
- `probes/`: axis-pressure probe histories (xz).

The working notes behind this package used an earlier smoothing script with a
window whose end samples depended on rounding; its averaged values differ from
`cp_averaged.py` by up to 0.1 % (wedge L1). The values above and in the paper
are those of `cp_averaged.py`.
