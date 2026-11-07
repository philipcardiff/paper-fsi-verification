# Cavity FSI xplt audit (M2 = C3_fsi_M2, dx 0.05; C2 = C2_fsi_pilot_f1, dx 0.1)

Tool: tools/cavity_xplt_audit.py RUNDIR OUT.json. Data: results/cavity_xplt_audit_{M2,C2}.json. Read-only; nothing tuned.
Forces are per unit depth (z = 0 layer line integrals, same convention as tools/cavity_fsi_check.py). Plate outward normal on the wetted face is +y (plate -> fluid);
force of fluid on plate = int(-p n + mu (grad w + grad w^T).n) dl, so a downward load is negative Fy.

## Reader and validation
xplt format 0x35 (upstream FEBioPlot/FEBioPlotFile.h): magic "BEF\0", then chunks (u32 tag, u32 size, payload), uncompressed (HDR_COMPRESSION = 0). Node data are (id, floats) sub-chunks of PLT_STATE_VAR_DATA
and node/element lists use 0-based node indices. Dictionary: nodal = displacement, fluid dilatation, reaction forces; domain (FMT_ITEM, one value per element) = fluid pressure, relative fluid velocity.
Note: there is NO nodal pressure/velocity in the xplt; nodal p = -K ef from "fluid dilatation"; velocity only as per-element averages.
Validation (numbers in JSON "validation"):
* u_y of point A (node 876 / 249) from xplt vs pointA.txt at every common state: max |diff| 5e-9 (M2: t=250,500,750; C2: t=100..900). M2 xplt t=750 uy = -0.2261943 (= pointA.txt at 750); the quoted final -0.2262325 is t=870, which is not an xplt state. u_z = 0 everywhere.
* xplt dilatation at the 81 interface nodes vs Mon_iface.txt ef: max diff 6e-14. Element "fluid pressure" vs -K * mean of 8 nodal ef: max diff 1e-7 Pa (max p 1.5 Pa) -> p = -K ef holds.
* C2: xplt element-average relative velocity vs mean of the 4 z=0 nodes of the same element from Mon_all.txt: max diff 1.8e-9 (so the per-element value is the centroid value, reader correct for velocity).
* Mon_iface wx, wy at the interface nodes are exactly 0 (no-slip, wall velocity 0 in steady state) -> w(wall)=0 used in the viscous estimate.

## States
xplt states written only every plot_stride steps. M2: t = 0, 250, 500, 750 (flag 0) and 872.5 (flag 2 = written on the solver-failure exit, not a converged step; the log t_last is 870). No state near 800 exists.
C2: t = 0,100,...,900 (flag 0) and 975 (flag 2). "Last complete" = M2 t=750, C2 t=900. Steadiness: M2 t=500 vs 750 and the failed t=872.5 state differ by <0.003 in Fy (<0.06 %); C2 800 vs 900 by 1e-4.

## (a) Reaction forces: NOT available
The "reaction forces" nodal variable is identically 0 (all nodes, all states, both runs; clamp-node sums 0). The FSI solver sets node loads only for dofs with equation index <= -2 (FEFluidFSISolver.cpp ~L1396) and apparently none are set here; the cause was not established without running FEBio. So the reaction-force total could not be used for a check, and "residual of the force balance" is not obtainable from this output. Not a result for or against any force value.

## (b) Pressure-only force from nodal p on the deformed interface (trapezoid on the z=0 polyline, exact for the trilinear field)
M2 t=750: Fy = -4.93389, Fx = +0.02198 (tool value at t=870: -4.93387 / +0.02198); t=500: -4.93098; t=872.5: -4.93387. C2 t=900: Fy = -4.71760, Fx = +0.02145.
Max interface slope 0.167 (M2), 0.157 (C2). Interface p: M2 min 1.150 (pin) max 1.3895 mean 1.2320; C2 min 1.052 max 1.3435 mean 1.178. Odd-even (even minus odd nodes) 0.0076 Pa (M2), -0.0708 Pa (C2).
Nodal p along x: JSON last_complete.interface_p_along_x (x_ref, x_def, y_def, p). M2 first nodes 1.150(pin) 1.2045 1.1719 1.1946 1.1696 ...; last 1.3616 1.3895 1.150(pin).

## (c) Viscous force
Method: wall facet (z=0 edge a-b, normal n = (-dy,dx)/L, outward from plate). Fluid hex E1 adjacent to the facet; 2nd-row hex E2 across the opposite face. With w = 0 at the wall, grad w = a (x) n, a = dw/dn, traction t = mu (grad w + grad w^T).n = mu (a + n (a.n)) (a.n ~ 0 for incompressible flow). a from the element-average velocities (= centroid values):
1pt: a = w1/s1 (first order, s1 = centroid distance along n, 0.026-0.031 m for M2); 2pt: quadratic w = a s + b s^2 through the E1, E2 centroids (second order on smooth fields).
M2 t=750: Fx_visc = -0.02380 (1pt), -0.02442 (2pt); Fy_visc = -0.00119 (1pt), -0.00124 (2pt).
C2 t=900: Fx_visc = -0.01862 (1pt), -0.01972 (2pt); nodal 3-point quadratic from Mon_all nodal w (independent of the element averaging): Fx = -0.02192, Fy = -0.00074; Fy_visc 2pt = -0.00067. So method spread in Fx_visc is ~ +-0.003 (10-15 %), in Fy_visc ~ 0.0001. Viscous Fx is negative (floor flow is reversed, cavity vortex); viscous Fy is tiny (only via the slope, since the wall shear is tangential).
M2 viscous has no nodal cross-check (Mon_all.txt was not logged for M2).

## (d) Total (pressure + viscous, 2pt), per unit depth
M2 t=750: Fy = -4.93513, Fx = -0.00245 (1pt: -4.93508, -0.00183). Pressure-only -4.93389 / +0.02198. Viscous changes Fy by -0.0012 (0.025 %) but flips and cancels Fx (+0.0220 -> -0.0024).
C2 t=900: Fy = -4.71827, Fx = +0.00172 (pressure-only -4.71760 / +0.02145).
Reaction-based total: unavailable (see a). Reference solids4foam total -5.23 at dx 0.05 differs from the FEBio field-based total by 0.295 N/m (5.9 %) in Fy; the viscous term and the pin bound below (-0.0096) cannot close this gap. Only the xplt fields were used; no statement is made here on the cause.

## 3. Pin bound (not a correction)
Pinned interface nodes have p = 1.1500000028 exactly (ef = -1.15e-6) at x=2 and 6 (all 4 nodes). Neighbours: M2 left 1.2045, 1.1719; right 1.3895, 1.3616 (interface values at x=5.95, 5.9). Replacing the two pinned z=0 nodal values by linear extrapolation of the two nearest nodes (left 1.2373, right 1.4174) changes Fy by -0.00957 (M2; 0.19 %) and Fx by +0.0005; C2: pinned 1.15, extrapolated 1.2405 / 1.4232, dFy = -0.0194 (0.41 %), dFx = +0.0008. The right pin node sits 0.24 Pa below its neighbour, the left 0.05 Pa; the effect is two half-facets only. It scales with dx (C2 about 2x M2).

## 4. Channel pressure (nodal -K ef, last complete state)
Outlet (x=14): ef = 0 exactly at all 42 outlet nodes (both z layers) in M2 and C2: p = 0. Inlet (x=0): mean 1.4298 (M2; trapezoid over y; min 1.212 max 1.635, i.e. strong odd-even in y), 1.4040 (C2); mean inlet flux int wx dy = 1.0000, outlet 0.99999.
The channel (y in [1,2], height 1) pressure carries a strong y checkerboard (columns alternate by ~0.085 Pa at x=13, ~0.4 Pa at x=1 for M2, e.g. 1.192/1.613 at x=0.98), the y-average (slab-average over 2 neighbouring columns) is clean. Slab-averaged mean pressure vs x (JSON channel_sections_all), linear fits (expected plane Poiseuille for H=1, U=1, mu=0.01: -0.12 Pa/m):
upstream x 0.3-1.7: dp/dx = -0.1237 (M2), -0.1248 (C2), mean p ~1.34/1.30 at x=1; above the cavity x 2.3-5.7: dp/dx = +0.0015 (M2), +0.0036 (C2), p ~1.28/1.30 (flat, i.e. NOT Poiseuille: the channel over the open cavity has no pressure drop); downstream 6.3-13.7: -0.1249 (M2) -0.1216 (C2), with local slope -0.139 at x 7-8 relaxing to -0.119 at x 12-13 (outlet development); p(13.5)=0.0517 vs 0.0600 expected from the outlet.
Extrapolating the downstream plane-Poiseuille line to x=2 gives 1.44 Pa, vs 1.28 measured above the cavity and 1.43 at the inlet. Inlet-outlet mean pressure drop 1.43 Pa over 14 m (0.102 Pa/m average) vs 1.68 Pa for plane Poiseuille all the way; the drop is concentrated outside the cavity span (about 0.12-0.125 Pa/m there) and is ~0 above it (dead-end cavity, pressure nearly uniform).
Cavity-side columns (x=2..6, y<1) were not averaged separately: the mesh is non-conforming in x between blocks.

## Caveats
* Only the z = 0 layer is used (flow plane, u_z = 0), per unit depth; force in N = value * TZ (0.1).
* Viscous estimate is first/second-order from element averages with w(wall)=0; spread ~10-15 % of Fx_visc (0.002-0.003 N/m); no effect on Fy at the 0.001 level.
* The xplt states are not at the final solver step (870 / 972.5); differences are < 0.06 % in Fy.
* The reaction-force check could not be done; no independent total force from FEBio's own solid residual is available from this output.
