> **ARCHIVED WORKING REPORT (verbatim except path/host sanitisation).** Written during the investigation; file paths such as `runs/`, `tools/`, `cavity_agent/` refer to the original working tree, not to this package (see `docs/INDEX.md` for the package equivalents). Later reports supersede earlier ones where they disagree; known numerical/wording inconsistencies and their resolution are listed in `provenance/INCONSISTENCIES.md`. Nothing here claims a verified FEBio result.

# CAVITY MESH CONVERGENCE: FEBio v4.13, three levels vs the current manuscript values (coordinator, 2026-10-30)
Same generator, physics and numerics at all levels (K 1e6, dt 2.5 s, ramp 25 s, E_mesh 1e-3, corner pin 1.15 Pa, hex8 StVK plate E = 500, nu = 0.4, plane strain, clamped ends, no-slip conformal interface, parabolic `fluid normal velocity` inlet, outlet p = 0). Point A = top wetted surface (4, 0); per 1 m depth.

**Validity of every row: NONE of M1, M2, M3 terminated normally.** M1 (14016, t = 975 s) and M2 (14019, t = 872.5 s) ended on `Max nr of reformations reached` (dilatation-increment floor); **M3 (14021, 8 cores, 3 h 31 min, exit 1) ended the same way after the step converged at t = 715 s.** 0 cut-backs, 0 negative Jacobians in all three. The values are *conditional steady estimates*: the fields were steady (or nearly) before the failure, but a failed termination is not a clean completion.

| | M1 (dx 0.1) | M2 (dx 0.05) | M3 (dx 0.025) | s4f partitioned (M1/M2/M3, finest 0.0125) | s4f QM (same) | extrapolated s4f |
|---|---|---|---|---|---|---|
| fluid / plate cells | 1688 / 160 | 6752 / 640 | 27008 / 2560 | | | |
| last accepted t (s) | 972.5 | 870 | 715 | | | |
| **u_y(A) (m)** | -0.211967 | -0.226233 | **-0.230038** (window [605,715] mean -0.230068) | -0.2046 / -0.2314 / -0.2384 (-0.2401) | -0.2036 / -0.2310 / -0.2382 (-0.2400) | -0.2407 |
| relative to s4f same dx | +3.6 % | -2.2 % | **-3.5 %** | | | |
| steady spread u_y over last 10 % | 7.6e-5 | 1.8e-4 | **1.05e-3** (criterion 1e-4: NOT met; decaying ~45 s oscillation of +/-1.2e-4 about -0.23007) | | | |
| u_x(A) (m) | 3.571e-4 | 3.354e-4 | **3.257e-4** | s4f QM 7.4e-5 .. 9.7e-5, ANSYS 4.1e-5 | | |
| **F_y pressure raw (N/m)** | -4.7176 | -4.9339 | **-4.9835** | | -5.1226 / -5.2266 / -5.2504 (-5.2560) | -5.26 |
| relative to s4f QM same dx | -7.9 % | -5.6 % | **-5.1 %** | | | |
| F_x pressure raw (N/m) | +0.02145 | +0.02198 | +0.02194 | | | |
| viscous F_y / F_x (2-pt, xplt audit) | -0.0007 / -0.0197 | -0.0012 / -0.0244 | not computed | | | |
| total F_x (pressure + viscous) | +0.0017 | -0.0025 | not computed | s4f total not in the record | | |
| odd-even amplitude on the interface (Pa) | -0.0708 | +0.0076 | **+0.0237** | | | |
| odd-even bracket on F_y | +/-3.0 % | +/-0.4 % | **+/-1.0 % (-5.0330 / -4.9340)** | | | |
| interface p min/max/mean (Pa) | 1.052/1.344/1.178 | 1.150/1.389/1.232 | 1.150/1.400/1.245 | README floor 1.28-1.31 | | |
| mass imbalance Q_out/Q_in - 1 | 1.2e-4 | -1.9e-5 | -2.3e-5 | | | |
| Newton iterations per step (mean / max) | 7.1 / 25 | 7.1 / 25 | 7.1 / 25 | | | |
| cost (alloc. cores x wall) | 4 x 15 min = 1.0 CPU-h | 4 x 66 min = 4.4 CPU-h | 8 x 211 min = **28.2 CPU-h** | | | |
Force comparisons use the raw pressure-only `F_y` (M2 total incl. viscous differs by 0.025 %). Not filtered. FEBio reaction forces are identically zero in the xplt (no independent total).

## Convergence assessment (observed orders only where the change exceeds the noise)
- **u_y(A):** differences -0.014266 (M1->M2) and -0.003805 (M2->M3): ratio 3.75, observed order **1.9**. The M2->M3 change (3.8e-3) is 16x the M3 oscillation amplitude (+/-1.2e-4) and 3x the M3 window uncertainty, so the order is outside the noise. Extrapolated (p = 2 / p = 1.9): **-0.2313 / -0.2314 m**, i.e. **-3.9 %** from the solids4foam extrapolation -0.2407. Part of that gap is the hex8 plate (S1: the 160 x 16 hex8 plate alone is 0.3 % too stiff), the rest follows the load deficit (below).
- **F_y:** differences -0.2162 and -0.0496: ratio 4.36, nominal order 2.1, **but the M2->M3 change (0.0496) equals the M3 odd-even bracket (+/-0.0495) and the bracket does not decrease monotonically (3.0 -> 0.4 -> 1.0 %)**: the observed order of F_y is therefore *not reliable* and I do not claim it. Extrapolating anyway as a bound: about **-5.00 N/m**, i.e. **-4.9 % from -5.26**. The deficit against the solids4foam same-dx value goes -7.9, -5.6, -5.1 %: it **shrinks and then flattens; it does not close** (the +/-1 % bracket is five times smaller than the deficit).
- **u_x(A):** 3.57, 3.35, 3.26e-4: converges to about 3.2e-4, **3.4-4.4x the solids4foam value and not approaching it**. Not explained; consistent with a left/right asymmetry of the floor pressure that differs from solids4foam (M2 pinned-corner neighbours: left 1.17-1.20, right 1.36-1.39 Pa).
- **Plate response vs load:** u_y/F_y = 0.0449 (M1), 0.0459 (M2), 0.0462 (M3), extrapolated 0.0463 m per N/m vs 0.0458 for solids4foam: the plate deflects as expected for the load it receives (S2: 0.025 %); the displacement deficit is a force (pressure-field) deficit, not a plate or transfer defect.
- **Conclusion:** u_y(A) converges (order about 1.9, conditional) to about -0.2314 m; F_y converges to about -5.0 N/m with +/-1 % odd-even uncertainty. **Both are about 4-5 % below the solids4foam StVK benchmark and the gap is not shrinking to zero with mesh refinement**, so this is not a mesh-resolution problem; the displacement and the force agree with each other (same pressure deficit). Not claimed: convergence to the benchmark, verification, or a clean solver completion.
