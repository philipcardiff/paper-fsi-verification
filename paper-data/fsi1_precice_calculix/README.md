# Hron–Turek FSI1 with preCICE: OpenFOAM fluid + CalculiX solid (November 2026)

This is an independent-solid solution of FSI1 (issue #36). The solid code and the coupling library are independent of solids4foam. The fluid is OpenFOAM, as in solids4foam, so the route is partly independent.

## Platform

- **Machine:** AMD Ryzen Threadripper 9970X desktop (32 cores, 125 GiB), Ubuntu 24.04.4, GCC 13.3.0, Open MPI 4.1.6.
- **Fluid:** OpenFOAM v2512 `pimpleFoam`, built from source, with the preCICE OpenFOAM adapter 1.4.0.
- **Solid:** CalculiX 2.20 (`ccx_preCICE`, adapter 2.20.2, Ubuntu noble package).
- **Coupling:** preCICE 3.4.1 (Ubuntu noble package).
- **Environment:** `OMP_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1` (`scripts/env-precice.sh`). The system FFTW is used instead of OpenFOAM's ThirdParty copy, because the system PETSc that preCICE links needs it.
- **Load:** no other jobs ran during the study. The runs were chained one after another, so the load average of 1–2 at each start (`logs/*/log.environment`) is the chain's own fluid and solid processes. The wall-clock times are clean.

## Set-up

The fluid side follows the solids4foam FSI1 verification (`TurekHron2006/verification`, `--benchmark fsi1`), so the comparison isolates the solid and the coupling.

- **Fluid mesh:** the `TurekHron2006` fluid blockMesh, with every in-plane block division multiplied by the level. This gives 5 336, 21 344 and 85 376 cells at 1x, 2x and 4x.
- **Fluid settings:**
  - ν = 1e-3 m²/s, ρ = 1000 kg/m³.
  - Parabolic inlet with a 0.3 m/s maximum (mean 0.2 m/s), ramped by 0.5(1 − cos(πt/2)) over 2 s (`codedFixedValue`).
  - Backward time scheme, `linearUpwind`, PIMPLE 3 outer / 3 pressure / 1 non-orthogonal corrector, solver tolerances 1e-9.
  - Mesh motion by `displacementLaplacian` with quadratic inverse-distance diffusivity.
- **Solid:**
  - CalculiX `*DYNAMIC` with `NLGEOM` and `*ELASTIC`, which is exactly St Venant–Kirchhoff. E = 1.4 MPa, ν = 0.4, ρ = 1000 kg/m³.
  - The plate is (105·L) × (6·L) × 1 C3D8I elements, 0.015 m deep (the fluid depth), with u_z = 0 (plane strain).
  - It is clamped on the cylinder face, and quasi-2D coupled through the adapter.
- **Pseudo-transient protocol, as in solids4foam:**
  - Δt = 0.025 s at 1x and 2x, 0.0125 s at 4x, to t = 30 s.
  - 0–2 s: the inlet ramp with a rigid, uncoupled plate.
  - 2–30 s: coupled.
- **Coupling:**
  - preCICE parallel-implicit with IQN-ILS (residual-sum preconditioner, QR2 filter 1e-3, 10 reused windows).
  - Relative convergence limit 1e-6 on both force and displacement, at most 50 iterations.
  - Force is mapped from fluid face centres to solid nodes, and displacement from solid nodes to fluid face nodes. Both use RBF mapping (compact C6 basis, support radius 0.03 m, i.e. 9–36 solid node spacings); the force mapping is conservative and the displacement mapping consistent.
- **Evaluation:**
  - Point A comes from a preCICE watch-point at (0.6, 0.2).
  - Drag and lift come from `forces` on `plate` + `cylinder`, divided by 0.015 m. Only the converged (last) entry of each time step is used.
  - Each quantity is the mean over the closing 20 % of the run (t = 24.4–30 s). The spread there must be below 1e-3, as in the solids4foam test.

## Mesh study (fluid and solid refined together)

| Quantity | 1x | 2x | 4x | Featflow 7+0 |
|---|---:|---:|---:|---:|
| u_x(A) (m) | 2.22770e-5 | 2.25482e-5 | 2.26847e-5 | 2.270493e-5 |
| u_y(A) (m) | 8.09851e-4 | 8.15650e-4 | 8.17387e-4 | 8.208773e-4 |
| drag (N/m) | 14.21665 | 14.27168 | 14.29062 | 14.29426 |
| lift (N/m) | 0.760195 | 0.762540 | 0.763115 | 0.7637460 |

| Quantity | Error 1x / 2x / 4x | Observed order (1x–2x–4x) | Richardson value | vs Featflow |
|---|---:|---:|---:|---:|
| u_x | −1.88 / −0.69 / −0.09 % | 0.99 | 2.2823e-5 | +0.52 % |
| u_y | −1.34 / −0.64 / −0.43 % | 1.74 | 8.1813e-4 | −0.33 % |
| drag | −0.54 / −0.16 / −0.03 % | 1.54 | 14.3006 | +0.04 % |
| lift | −0.46 / −0.16 / −0.08 % | 2.03 | 0.76330 | −0.06 % |

Every error decreases monotonically. The Richardson values agree with the solids4foam four-level study (`../fsi1_four_level`, 2x–4x–8x) to within 0.02 % in every quantity: 2.2824e-5, 8.1825e-4, 14.3026 and 0.76318.

At 4x, the errors are close to those of the solids4foam high-order (p = 3) solid (`../fsi1_high_order`): −0.27 / −0.48 / −0.04 / −0.10 %.

As in that study, u_x approaches Featflow from below, while the Featflow sequence approaches its limit from above.

## Solid-only refinement on the 1x fluid mesh

| Quantity | Solid 105×6 | 210×12 | 420×24 | Change 105×6 → 420×24 |
|---|---:|---:|---:|---:|
| u_x(A) (m) | 2.22770e-5 | 2.22303e-5 | 2.22309e-5 | −0.21 % |
| u_y(A) (m) | 8.09851e-4 | 8.10682e-4 | 8.11035e-4 | +0.15 % |
| drag (N/m) | 14.21665 | 14.21664 | 14.21664 | −0.000 % |
| lift (N/m) | 0.760195 | 0.759840 | 0.759690 | −0.07 % |

A 16-fold solid refinement moves u_y(A) by 0.15 %. The C3D8I plate is therefore converged at 105 × 6, and the remaining 1x error (−1.34 % in u_y) comes from the fluid mesh. The solids4foam standard (p = 1) solid on the same fluid mesh gives −14.21 % at 1x and −4.63 % at 2x.

## Cost and coupling

| Run | Wall clock (s) | Cores | Coupled windows | Mean / max coupling iterations | Unconverged windows |
|---|---:|---:|---:|---:|---:|
| 1x | 222 | 2 (1 fluid + 1 solid) | 1120 | 2.11 / 17 | 0 |
| 2x | 1 598 | 2 | 1120 | 2.13 / 16 | 0 |
| 4x | 8 090 † | 9 (8 fluid + 1 solid) | 2239 | 2.30 / 20 | 0 |
| 1x fluid, 210×12 solid | 759 | 2 | 1120 | 2.10 / 16 | 0 |
| 1x fluid, 420×24 solid | 3 346 | 2 | 1120 | 2.11 / 16 | 0 |

The wall-clock times include the uncoupled 0–2 s ramp (15 s at 1x, 101 s at 2x, 227 s at 4x).

† 4x: the time from start to the last completed window (t = 29.9875 s): 227 s ramp + 7 857 s coupled (`ClockTime` in the fluid log). The final window ended in the CalculiX abort described under Notes. The total in `logs/level4/log.run` (30 150 s) also includes about 6 h in which the aborted solver process hung before it was killed, and is not the cost.

The largest spread over the closing 20 % is in the lift: about 5e-4 on all 1x and 2x runs, and 1.9e-3 at 4x. At 4x this is a ±0.1 % step-to-step coupling noise at the 1e-6 tolerance, not drift. The 1 s means of the lift over t = 20–30 s agree to 7e-6, so the tabulated 4x mean is converged in time to the digits shown. The solids4foam p = 3 4x run shows the same kind of lift fluctuation. The displacement spreads are at most 1.1e-5.

## Notes

- **4x final window:** CalculiX aborted inside its SPOOLES solver (`ChvManager_newObjectOfSizeNbytes`: bad input, a negative size) in the last of 2240 windows, at t = 29.99 s. The 2239 completed windows all converged, and the run had been steady since about t = 10 s, so the results are unaffected. The run's `log.run` shows the resulting non-zero return codes, and its total time includes the hang (see the cost table).
- **Mapping choice:**
  - A global RBF support (0.35 m, the plate length, as in the preCICE `turek-hron-fsi3` tutorial) made the displacement mapping ill-conditioned on the 2x and 4x solids. The force residual then stalled at about 1e-2 and the coupling diverged. A local 0.03 m support fixed this.
  - The preliminary level-1 result posted on #36 used nearest-neighbour force mapping (u_y −0.55 %). With the final RBF mapping, level 1 gives −1.34 %. The mapping therefore affects u_y by about 0.8 % at 1x. All results here use the final mapping.
- **Coupling tolerance:** at level 1, a relative tolerance of 1e-5 left ±0.2 % step-to-step noise in the lift. 1e-7 hit the 50-iteration cap in 7 transient windows. 1e-6 gave the same steady values to 4–5 digits with no unconverged windows.
- **Independence:** the fluid is OpenFOAM, as in solids4foam, so this is the "partly independent" route. The solid (CalculiX, finite elements) and the coupling (preCICE IQN-ILS) are independent.

## Files

- `fsi1_precice_sweep.csv`: one row per run, with steady values, closing-window spreads, coupling statistics and cost.
- `histories/*.csv`: time histories of u_x(A), u_y(A), drag and lift (converged value per step) for the five runs.
- `fsi1_precice_histories.png`: the 1x/2x/4x histories with the Featflow values (`scripts/plot-histories.gnuplot`).
- `case/`: the template, i.e. the preCICE config, CalculiX `plate.inp` and `config.yml`, and the fluid `0/`, `constant/` and `system/`. The `@DT@` placeholders and the meshes are filled in by `scripts/setup-case.py`.
- `scripts/`:
  - `setup-case.py`: generates a case for a fluid and solid level.
  - `run-case.sh`: runs the ramp, then the coupled run, with timing.
  - `run-study.sh`: runs all five cases.
  - `evaluate.py`, `study-table.py`: the evaluation.
  - `env-precice.sh`: the environment.
- `logs/<run>/`: `log.run` (stages and timing) and `log.environment`. The solver logs and field data are not included (available on request).

To reproduce:

```bash
source scripts/env-precice.sh
scripts/run-study.sh /path/to/workdir
```
