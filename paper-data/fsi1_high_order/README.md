# Hron–Turek flap and FSI1 with the high-order (cubic) solid

solids4foam development `0a51b0f0c` (2026-11-05), OpenFOAM v2512, Apple M-series (Clang), 2026-11-05/06.

## `flap/`: solid alone

The FSI3 flap, 0.35 × 0.02 m, St Venant–Kirchhoff, E = 1.4 MPa, ν = 0.4, plane strain.
It is clamped at x = 0 and loaded on the upper face by a pressure ramped to 20 Pa in 10 static steps.
The tip deflection is about 10% of the length, as in FSI3, and has no polynomial exact solution.

- Meshes have 6, 12, 24 and 48 square cells through the thickness.
- p = 1 is the standard solid (`diffStencilLaplacian` 0.5).
- p = 2 and 3 use `highOrderCoeffs` (moving least squares, `faceStencilExtraCells` 15, `alpha` stabilisation 0.1, `highOrderJacobian true`).
- p = 4 is not supported (orders 1–3 only).
- `flap_study.sh` generates and runs the cases, and `flap_table.py` writes `flap_results.md`.
- The error is measured against the Richardson value of the p = 3 sequence, −33.70365 mm.

| Cells thick | p = 1 | p = 2 | p = 3 |
|---:|---:|---:|---:|
| 6 | 21.6% | 5.9% | 0.17% |
| 12 | 7.0% | 0.91% | 0.066% |
| 24 | 1.9% | 0.16% | 0.030% |
| 48 | 0.52% | 0.033% | 0.014% |

The p = 3 sequence has a floor near 0.01–0.03% (observed order 1.1 on the finest three meshes).
The clamped corner and the residual stabilisation are the likely causes; this is not investigated further.

## `fsi1/`: Hron–Turek FSI1 with the tutorial verification driver

Both runs use the tutorial driver: `Allverify --benchmark fsi1 --study mesh --levels 1,2 --cores 4`.

- **p = 1:** `TurekHron2006` unchanged; results in `p1_sweep.csv`.
- **p = 3:** a copy whose `constant/solid/solidProperties` differs only as in `p3_solidProperties.diff`; results in `p3_sweep.csv`.

The driver's pseudo-time step (0.025 s) and coupling test differ from the four-level study (`../fsi1_four_level`, Δt = 0.00625 s, strict test).
At 2x the p = 1 driver run agrees with the four-level 2x row to 0.01 percentage points in every QoI; at 1x it differs by 0.4% in u_x only.

The p = 3 2x run shared the machine with a 16-rank job, so its wall time (5 308 s against 1 639 s for p = 1) overstates the cost.
The 1x pair (1 860 s against 1 030 s, also partly contended) is a better guide.

Errors are measured against the Featflow finest level:

| Level | Solid | u_x | u_y | drag | lift |
|---|---|---:|---:|---:|---:|
| 1x | p = 1 | −2.81% | −14.21% | −0.80% | +5.19% |
| 1x | p = 3 | −3.92% | −1.74% | −0.80% | −0.46% |
| 2x | p = 1 | −0.97% | −4.63% | −0.24% | +1.56% |
| 2x | p = 3 | −1.33% | −0.83% | −0.24% | −0.17% |
| 4x | p = 3 | −0.27% | −0.48% | −0.04% | −0.10% |

The 4x p = 3 level ran on MeluXina (job 5426960, 32 cores, 11.8 h, same commit), results in `p3_4x_sweep.csv`.
The driver halves the step with each mesh level, so 4x uses Δt = 0.0125 s.
The driver marks the run as failed because the lift spread over the closing window is 3.5 × 10⁻³ against its 10⁻³ tolerance.
The tip displacement is steady over t = 10–30 s (u_x to 1 × 10⁻⁵ and u_y to 6 × 10⁻⁵ relative), and the lift fluctuates by about ±0.2% about its mean, so the tabulated values are converged in time to the digits shown.
