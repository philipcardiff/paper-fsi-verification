# FEBio v4.13 source-code investigation (read-only; no source change)

Scope: `FEBioFluid` (fluid-FSI domain, FSI traction load, inlet load, plot routines) and `FECore/FESurface.cpp` of the **unmodified tag v4.13** (tree clean, 0 changed files; binary provenance limitation in `../provenance/BUILD_PROVENANCE.md`). Line numbers refer to that tag; the FEBio source is not copied here. Labels as in the findings report. Full audits: `reports/FEBIO_V413_PRESSURE_SOURCE_AUDIT.md`, `reports/FEBIO_FSI_TRACTION_ASSEMBLY_AUDIT.md`, `reports/FEBIO_PRESSURE_MODE_ANALYSIS.md`, `reports/FEBIO_PRESSURE_PROJECTION_DERIVATION.md`.

| item | finding | label | test id |
|---|---|---|---|
| element | `FEFluidFSIDomain3D` ("fluid-FSI-3D") for every solid-element shape; 7 unknowns per node (u, relative fluid velocity w, dilatation e_f); equal-order trilinear hex8, 8-point Gauss; one fluid-FSI domain class only (no Taylor-Hood-type option from the input file) | PROVEN | F4, F5 |
| pressure | not an unknown: p = -K e_f (`FELinearElasticFluid`); Gauss-point pressure is the trilinear interpolant of nodal e_f; plotted element pressure = mean of the 8 nodal values (xplt = -K mean e_f to 1e-7 Pa) | PROVEN | F4 |
| momentum row | pressure enters as the strong-form gradient `int N_i grad p dV`; no boundary pressure integral | PROVEN | F6 |
| continuity row | `int N_i e_f-dot/J + int gradN_i . w` (divergence integrated by parts); boundary flux term only via `FEFluidNormalVelocity::LoadVector` at the inlet | PROVEN | F6 |
| operator | kJw = G^T (max difference 1.7e-17); system [A G; G^T -(c/K) M] with c/K = 5.3e-7 at dt 2.5 s, rhoi 0 (alpha_f 1, alpha_m 2, gamma 1.5) | PROVEN | F6 |
| stabilisation | `grep stabiliz|pspg|brezzi|supg|upwind`: only outflow-boundary backflow / tangential-flow stabilisations; no interior pressure stabilisation in `FEFluidFSIDomain3D` or `FEFluidDomain3D` | PROVEN (v4.13) | F11 |
| row-parity mode | y-parity, checkerboard and x-parity families of the Schur complement S = G^T A^-1 G have eigenvalues proportional to h^2 (ratios 3.97/3.51/3.63 between h = 0.2 and 0.1) while the smooth mode is h-independent (0.97); the interior x-uniform y-alternation is an exact null of the interior gradient, made a soft mode by the outlet / wall rows | PROVEN (two levels) | F7, F8 |
| not a null mode | lowest parity eigenvalue 1e-3, 2000x the compressibility term; backward-Euler Stokes emulation removes a 0.1 Pa perturbation in one step, independent of K and dt: a forced steady response, not a stored history | PROVEN (linear emulation) | F7 |
| forcing | exact linear attribution on C1b, D = 0.417 Pa: viscous 0.177 + convective 0.240 + fixed values (outlet, rim, pins) -0.0001 | PROVEN | F7 |
| amplification chain | cell-scale content of the separated-flow forces projected on the soft modes' left singular vectors, divided by sigma ~ h | HYPOTHESIS | F7 |
| assembled residual | independent python re-assembly reproduces FEBio's converged fixed-plate C1b state: momentum 9.4e-7, continuity 2.2e-7 | PROVEN (fluid-only domain) | F9 |
| interface traction (F1) | `FEFluidFSITraction::LoadVector`: pressure from the 4 facet nodes' own e_f (bilinear, 2x2 Gauss), `f = gt p - sv gt`, viscous stress `sv` the Gauss mean over the adjacent element | PROVEN | F10 |
| plot definition (F2) | `FEPlotFluidSurfaceForce`: adjacent-element Gauss-mean stress x area vector; the pressure part equals the pair-row mean F4 for planar facets | PROVEN (planar facets) | F10 |
| consequence | F1 sees the wall-row nodal values (and therefore the row parity) fully; F2 averages rows 0 and 1; F1-F2 equals the parity part to <=0.001 N/m | PROVEN decomposition; STRONGLY SUPPORTED mechanism | E6, F10 |
| reaction check | because the momentum row carries no boundary pressure integral, the Dirichlet-node reaction does not contain the pressure traction (floor reaction -0.82 N/m vs -4.66 F1): there is no internal identity that picks F1 or F2 | PROVEN (algebra and data) | F10 |
| pins | corner pressure DOF of the two floor-end corners is almost undetermined without a constraint (Schur eigenvalue 3.5e-10 vs 2.5e-3 for the next mode); reproduces the C1 drift; pins contribute 1e-4 Pa to D | PROVEN | A7, C7 |
| inlet | `fluid normal velocity` parabolic map normalises the linear interpolant to mean 1; nodal profile = continuum parabola of flux 1/(1-h^2); `fluid velocity` + `math` evaluates the profile at the origin for nodal velocities in v4.13 | PROVEN (source + straight channel) | A5, C0, C1 |
| observation, no effect shown | `FESurface::LoadVector` sets `Hs = el.Gr(...)` (copy-paste); unused by `FEFluidFSITraction` | NOT YET TESTED for other loads | - |

## Classification
Class B, discrete numerical deficiency of the unstabilised equal-order Q1/Q1 pair: **STRONGLY SUPPORTED** for the existence and size of the parity. Class A (setup error): not supported. Class C (source defect): not supported (the documented element choice explains the behaviour; no implementation error found in the paths read). Class D (inadequate evidence): applies to the 5 % force gap: that it is entirely the parity offset is **not established** (F2 is +6.1 % off at M1 and +0.9 % at M2, no solids4foam pressure or traction field was available, the viscous part of F1 is unreconstructable, all coupled runs ended abnormally and M3 is not steady).

## What this investigation does not show
That FEBio is "wrong" or "incapable"; that F2 or the pair mean is the correct force; that any source change would restore the benchmark value. Source options (polynomial pressure projection, Brezzi-Pitkaranta; `reports/CAVITY_SOURCE_CORRECTION_OPTIONS.md`) were proposals; the projection was tested experimentally and failed (`FAILED_STABILISATION_HISTORY.md`).
