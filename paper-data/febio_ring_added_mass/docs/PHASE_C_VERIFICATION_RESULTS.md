# PHASE C VERIFICATION RESULTS (ringAddedMass, FEBio v4.13 fluid-FSI)

Scope: Stage C0 (interface diagnostic, see `PHASE_C_INTERFACE_DIAGNOSTIC.md`) and Stage C1 (focused verification of the **inviscid-limit formulation**: free-slip fluid-ring interface + free-slip outer wall). This is an analytical-limit formulation, **not** a replica of the finite-viscosity (nu_f = 1e-7) no-slip benchmark. Pristine v4.13 binary (sha256 `f5ae1ed6…`), per-node fluid IC workaround, MKL Pardiso, 1 core per run, K = 1e8, rho_f = 140, rhoi = 1, 100 steps/period unless stated, compatible IC (v2).

All 8 budgeted runs were used: C0, M2 wet, K, tolerance (failed), dt dry, dt wet, M3 dry, M3 wet. Reference: dry 2.55504257, wet 1.80451344, ratio 0.70625572 rad/s. The 0.3 % wet/dry tolerance is a screening criterion only.

## 1. Mesh study (matched dry/wet; `results/phaseC_mesh_convergence.csv`)

| level | ring x fluid cells (radial) | wet nodes | omega_dry | err_dry | omega_wet | err_wet | ratio | ratio err | mu_eff |
|---|---|---|---|---|---|---|---|---|---|
| M1 | 4 / 12 | 1666 | 2.621064 | +2.584 % | 1.854296 | +2.759 % | 0.707459 | **+0.170 %** | 0.9980 |
| M2 | 8 / 24 | 6402 | 2.571081 | +0.628 % | 1.819071 | +0.807 % | 0.707512 | **+0.178 %** | 0.9977 |
| M3 | 16 / 48 | 25090 | 2.558432 | +0.133 % | 1.810158 | +0.313 % | 0.707526 | **+0.180 %** | 0.9976 |

- Dry and wet frequencies converge at an observed difference-based order of **1.98** each (M1-M2-M3, factor-2 refinement; three levels only, so this is a consistency check, not a formal order claim). Richardson with p = 2 assumed: dry 2.554216 (-0.032 %), wet 1.807187 (+0.148 %).
- The **ratio** error is flat: +0.170, +0.178, +0.180 %. Increments 0.0075 then 0.0020 (consistent with a small second-order tail; extrapolated roughly +0.181 %). The ratio is therefore mesh-converged to within ~0.002 percentage points and sits at **+0.18 %**, inside the 0.3 % screening tolerance, but not at zero. The persistent offset is **not** a spatial-resolution effect.
- The wet absolute frequency retains an extrapolated +0.15 % error versus a dry -0.03 %; the difference between them is the 0.18 % in the ratio.

## 2. Pressure and traction (`results/phaseC_pressure_traction.csv`, `figures/fig_C4_pressure_traction.png`)

| level | n=2 added mass layer0 / exact profile | n=2 integrated force error | leakage rms / max (relative to ring normal velocity / V0) | n=2 pressure phase | KE at integer periods | mean(J-1) max |
|---|---|---|---|---|---|---|
| M1 | 0.9993 | -2.12 % | 1.04 % / 1.25 % | ~178 deg | 0.9995-0.9999 | 5.1e-10 |
| M2 | 0.9996 | -1.06 % | 0.40 % / 0.49 % | ~178 deg | 0.9994-0.9998 | 5.3e-10 |
| M3 | 0.9997 | -0.53 % | 0.32 % / 0.43 % | 178.2 deg | 0.9995-0.9999 | 5.3e-10 |

- The measured n=2 pressure amplitude follows the exact radial profile (`A r^2 + B r^-2`) to 0.03 % at M3 (profile-corrected criterion P2' of C0). The integrated n=2 force error halves with each refinement (first-order behaviour of the nodal-traction measure, -2.1 -> -1.06 -> -0.53 %); the frequency error is second order, so the two measures converge differently and neither should be quoted as the other's order.
- Interface normal leakage (a weakly enforced condition) falls with h (1.04 -> 0.40 -> 0.32 % rms), mean leaked volume ~5e-7 m^2 per unit thickness at M3, negligible. No unphysical leakage, no drift in volume (mean J-1 ~5e-10), energy at integer periods conserved to <0.06 %, amplitude loss ~3e-6 per period.
- Pressure noise relative to the analytic amplitude is ~6 % at every level (`p_noise_rel`), mesh-independent, residual power above 20 rad/s: a checkerboard/dilatation-noise level, not a physical mode. Its origin is untested; it does not bias the frequency fit (window spread 3e-5 %).

## 3. Cost (`results/phaseC_costs.csv`, `figures/fig_C6_cost.png`)

| level | dry wall / wet wall | wet nodes | wet RSS |
|---|---|---|---|
| M1 | 12 s / 256 s | 1666 | 114 MB |
| M2 | 48 s / 1493 s | 6402 | 342 MB |
| M3 | 244 s / 9300 s (2 h 35 min) | 25090 | 1386 MB |

Wet cost scales roughly 5.8x then 6.2x per refinement on 1 core, close to the earlier extrapolation (~3.5 h predicted, 2.6 h actual, node load dependent). The wall time is dominated by the sparse factorisation; Phase A measured poor thread scaling.

## 4. Sensitivities at M1 (see `PHASE_C_NUMERICAL_UNCERTAINTY.md`)

| change | d omega_wet | d omega_dry | d ratio | ratio error after |
|---|---|---|---|---|
| dt halved (paired) | +0.0261 % | +0.0260 % | +9e-5 % | 0.17050 % |
| K 1e8 -> 1e7 | -0.0004 % | (none) | -0.0004 % | 0.17000 % |
| mesh M1 -> M2 | -1.900 % | -1.907 % | +0.0075 % | 0.1779 % |

## 5. Comparison with the Phase B no-slip family (matched dry)

| configuration | ratio error |
|---|---|
| B2 no-slip, original IC | -6.60 % |
| R1 no-slip, compatible IC | -6.01 % |
| D4 E_mesh 1e-3 (no-slip) | -6.16 % |
| D2 no-slip, M2 | -2.67 % |
| D5 wall slip, interface no-slip (M1) | -1.11 % |
| **C0 wall slip + interface slip (M1)** | **+0.170 %** |
| C1 M2 / M3 | +0.178 / +0.180 % |

The original no-slip runs are preserved and not reclassified. The slip results answer a different (inviscid, free-slip) model than the no-slip benchmark: the residual of the no-slip family is attributed to the strong nodal no-slip constraints (viscous boundary layer unresolved at nu_f = 1e-7 relative to the mesh), not to a FEBio defect. Whether native FEBio can reproduce the **finite-viscosity no-slip** case at the 0.3 % level is **untested** (only M1/M2 no-slip meshes were run).

## 6. Figures and data
`figures/fig_C1_frequency_vs_mesh.png`, `fig_C2_ratio_vs_mesh.png`, `fig_C3_errors.png`, `fig_C4_pressure_traction.png`, `fig_C5_temporal_iterative_sensitivity.png`, `fig_C6_cost.png`. CSVs in `results/`: `phaseC_mesh_convergence.csv`, `phaseC_pressure_traction.csv`, `phaseC_costs.csv`, `phaseC_time_iter_sensitivity.csv`. Per-run reports: `runs/phaseC/*/report.json`, `pressure_modes.json`, `interface_leakage.json`.

## 7. Addendum: figure review and regulariser run
Figures C1-C6 were visually reviewed and corrected for presentation (no numbers changed): mesh axes labelled "finer to the right"; C4 now marks leakage as n/a where w=0 is enforced (previously drawn as zero), shows the exact n=2 profile at every layer and the layer-0 sampling offset for the force panel, and uses explicit family labels; C5 now has ratio and wet-omega panels, an "iterative error not quantified" note, and the E_mesh row; C6 annotates mesh levels. The E_mesh 1 -> 1e-3 run is in `results/phaseC_regulariser_diagnostic.csv` and `PHASE_C_REGULARISER_DIAGNOSTIC.md`; the table in section 1 is unchanged (baseline E_mesh = 1 Pa), and the +0.18 % plateau there is now attributed at M1 to the regulariser (confirmation at M2/M3 pending).

## 8. Correction (added in Phase D)
The "0.3 % screening tolerance" for the wet/dry ratio used in this document is **wrong**: the official solids4foam JSON gives the ratio tolerance as 0.1 % (omega 0.3 %). Phase D (`PHASE_D_FINAL_COMPARISON.md`) uses 0.1 % and the exact JSON references. The E_mesh = 1 Pa +0.18 % ratio offsets in this document exceed the official 0.1 % band; the regulariser diagnostic removed them (+0.006/+0.007/+0.008 % at E_mesh = 1e-3 Pa).
