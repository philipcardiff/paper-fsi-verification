# Stage B: structural check of the flag (CSM1, CSM2)

Reference (Featflow CSM page, fetched 2026-10-31; level 4+3, steady state, E_s = 1.4e6 (CSM1) / 5.6e6 (CSM2), nu_s = 0.4, rho_s = 1000, g = 2 m/s^2): CSM1 u_x = -7.18777e-3 m, u_y = -66.1029e-3 m; CSM2 u_x = -0.469006e-3 m, u_y = -16.9740e-3 m (the 5+1 values differ in the fifth digit: -66.1023, -7.18767, -16.9739, -0.469000 mm). The page does not state the geometry, clamp, point-A coordinates, the plane-strain assumption or the constitutive law (the image name `csm1b_svk` suggests SVK); these were taken from the tutorial and the Featflow definition page (`PROVENANCE_AND_CONFIGURATION.md` section 2). The numbers supplied for the task were identical to the page.

## Results (final state, load factor 1; `results/stageB_summary.csv`)
| case | elements | u_x (mm) | u_y (mm) | u_x error | u_y error | Newton iterations per load step |
|---|---|---|---|---|---|---|
| B_CSM1_h20_nt2 | 70 | -7.14358 | -65.86394 | -0.615 % | -0.362 % | 7.0 |
| B_CSM1_h20_nt4 | 280 | -7.19859 | -66.15398 | +0.151 % | +0.077 % | 7.0 |
| B_CSM1_h8_nt16 | 4496 | -7.15400 | -65.96222 | -0.470 % | -0.213 % | 7.0 |
| B_CSM1_h8_nt16_root025 | 4480 | -7.01643 | -65.23724 | -2.384 % | -1.310 % | 7.0 |
| B_CSM1_h8_nt2 | 70 | -4.45778 | -52.14953 | -37.981 % | -21.109 % | 6.7 |
| B_CSM1_h8_nt4 | 280 | -6.30316 | -61.92822 | -12.307 % | -6.315 % | 7.0 |
| B_CSM1_h8_nt8 | 1120 | -6.96476 | -65.08289 | -3.103 % | -1.543 % | 7.0 |
| B_CSM2_h20_nt2 | 70 | -0.46603 | -16.91050 | -0.635 % | -0.374 % | 4.0 |
| B_CSM2_h20_nt4 | 280 | -0.46974 | -16.98785 | +0.157 % | +0.082 % | 4.0 |
| B_CSM2_h8_nt16 | 4496 | -0.46670 | -16.93560 | -0.492 % | -0.226 % | 4.0 |
| B_CSM2_h8_nt16_root025 | 4480 | -0.45740 | -16.74206 | -2.475 % | -1.366 % | 4.0 |
| B_CSM2_h8_nt2 | 70 | -0.28614 | -13.25532 | -38.990 % | -21.908 % | 4.0 |
| B_CSM2_h8_nt4 | 280 | -0.40910 | -15.84965 | -12.772 % | -6.624 % | 4.0 |
| B_CSM2_h8_nt8 | 1120 | -0.45384 | -16.69814 | -3.233 % | -1.625 % | 4.0 |
Linear-regime runs (g = -0.02, E = 1.4e6; informational): | B_LIN_h8_nt16 | 4496 | -0.00075 | -0.67867 |  |  | 3.0 | ; | B_LIN_h8_nt8 | 1120 | -0.00073 | -0.66912 |  |  | 3.0 |. Plane-strain Euler-Bernoulli: -0.68311 mm, with Timoshenko shear (kappa = 5/6): -0.68607 mm (the hex8 nt 16 value is 0.65 % / 1.08 % below them; beam theory is only accurate to about 1 % at h/L = 0.057 with a full-face clamp). Nonlinearity of CSM1: u_y(g = 2)/(100 u_y(g = 0.02)) = 0.972.

## Mesh convergence (`results/stageB_convergence.csv`, `figures/fig1_stageB_convergence.png`)
| benchmark | quantity | monotone nt 2-4-8-16 | observed order (nt 4-8-16) | Richardson error vs reference | finest hex8 error | hex20 nt 4 error | hex20 nt 2 error | root 0.25, nt 16 error |
|---|---|---|---|---|---|---|---|---|
| CSM1 | uy_m | True | 1.84 | +0.30 % | -0.213 % | +0.077 % | -0.362 % | -1.31 % |
| CSM1 | ux_m | True | 1.81 | +0.58 % | -0.470 % | +0.151 % | -0.615 % | -2.38 % |
| CSM2 | uy_m | True | 1.84 | +0.32 % | -0.226 % | +0.082 % | -0.374 % | -1.37 % |
| CSM2 | ux_m | True | 1.80 | +0.61 % | -0.492 % | +0.157 % | -0.635 % | -2.47 % |
Observed order 1.8 from nt 4-8-16; hex8 converges from below and its Richardson limit is 0.30 % (u_y) and 0.59-0.61 % (u_x) BEYOND the reference; hex20 nt 4 is already +0.08 % / +0.15 %. The converged FEBio value therefore differs from Featflow by 0.1-0.3 % in u_y: UNRESOLVED (candidates: Featflow's own discretisation, the root geometry, the 0.35 vs 0.351 length, the clamp definition; none was tuned or tested further). Locking: hex8 with nu = 0.4 suffers bending (shear) locking, -21 %, -6.3 %, -1.5 %, -0.2 % in u_y for nt = 2, 4, 8, 16 at aspect ratio 1; hex20 nt 2 is already -0.36 % (PROVEN by comparison); volumetric locking is not the dominant mechanism because lambda/mu = 4 (not near incompressible). hex8 therefore needs at least 8 elements across the 0.02 m thickness for 1 % in displacement.

## Nonlinear convergence
All 16 jobs: 10 load steps, every step converged with the declared tolerances, normal termination (log inspected, not the exit code); 6.7-7.0 equilibrium iterations per step for CSM1 (full Newton with line search, 70 iterations in 10 steps at nt 16), 4.0 for CSM2 and 3.0 for the linear case; no cut-backs; max |u_z| = 0.

## Body-force sign and units (diagnostic job 14165)
FEBio `body_load type="const"` multiplies the stated vector by the density and integrates over the reference volume (`FEBioMech/FEElasticSolidDomain.cpp` `BodyForce`: `fa = -H density f J0`), so the input is an acceleration (rho f is the force density, rho|g| = 2000 N/m^3). The first job (generator v1, y written as -2) deflected the flag UP (u_y = +61.93 mm at nt = 4); the written y was therefore flipped once to +2 (generator v2) and all inputs regenerated; the v1 job is kept (`results/stageB_raw/DIAG_v1_signtest_CSM1_h8_nt4`, ledger row marked diagnostic). Final check against physics: point A ends at y = 0.134 m, i.e. 66 mm BELOW its initial y = 0.2 m (every B_* run: `point_A_deformed_y` column, u_y < 0): the tip deflects down, and the magnitude agrees with Featflow to 0.2-0.5 %. PROVEN that, in this build, the sign of the written y component is opposite to the acceleration it represents.

## Support reactions and equilibrium evidence
The logged reaction forces are zero (`results/stageB_raw/*/react_root.txt`). Cause (CONFIRMED FROM SOURCE): the support is a `zero displacement` boundary condition, which makes the DOFs fixed (no equation numbers); `FESolidSolver2::Residual` writes nodal reaction loads only for nodes with an equation number or a prescribed DOF (`FEBioMech/FESolidSolver2.cpp` l.1630-1642), so a fixed DOF never receives a reaction. This is a limitation of how the support was specified (a prescribed zero displacement or a rigid body would have reported reactions), not evidence about equilibrium; the plan's check (sum R_y = weight, 1e-6) was therefore NOT performed. The generator-v2 inputs also request the `reaction forces` plot variable: zeros as well.
Replacement from stored output only (`scripts/equilibrium_from_xplt.py`, `results/stageB_equilibrium_from_xplt.csv`, `figures/fig2_stageB_shear_equilibrium.png`): the stored xplt holds one Cauchy stress tensor per element; the free-body shear equilibrium of the part of the flag beyond a vertical section, sum of the section-averaged sigma_xy dA against minus the weight rho|g|(L-s) h TZ, is 1.0002 / 1.0004 / 0.9964 / 0.9942 / 0.9940 (CSM2, hex8 nt 16, five sections from the root to the tip), 0.996-1.001 (CSM2, hex20 nt 4) and 1.003-1.005 near the root falling to 0.91 at the tip for CSM1 (the section rotates by up to 0.28 rad in CSM1; the check uses x-normal sections and the small-slope approximation, so deviations of order 1 - cos(theta) are expected); the moment of sigma_xx about the section centre agrees with the beam value to 1-2 % (CSM2) away from the tip. SUPPORTED, approximate (element-average stress, rotation neglected, about 0.5 % scatter): not a 1e-6 equilibrium check and not a substitute for reaction forces.

## Flag-root geometry
Root on the cylinder circle (x = 0.2489898, length 0.351010) was chosen before the runs because (i) the Featflow definition states that the bar's left end is "fully attached to the fixed cylinder" and the nominal x = 0.25 end lies 1.0 mm outside the circle at y = 0.19/0.21 (it would not touch the cylinder), and (ii) the tutorial's blockMeshDict uses exactly these vertices. The x = 0.25 variant (length 0.35) was run as a diagnostic at nt = 16 and is 1.3 % (u_y) / 2.4-2.5 % (u_x) stiffer: it was not selected. The better agreement of the circle root supports, but is not an independent proof of, the Featflow geometry (UNRESOLVED without the Featflow mesh).

## Classification
Tip displacement, both components, both moduli: QUALIFIED AGREEMENT (strong, mesh-convergent, like-for-like model, 0.2-0.5 % hex8, 0.08-0.16 % hex20). Not labelled VERIFIED AGREEMENT: equilibrium evidence incomplete, a 0.1-0.3 % unexplained offset of the converged value, geometry of the reference not independently documented. Evidence label SUPPORTED.
