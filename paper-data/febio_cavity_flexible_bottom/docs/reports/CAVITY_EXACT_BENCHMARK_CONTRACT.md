> **ARCHIVED WORKING REPORT (verbatim except path/host sanitisation).** Written during the investigation; file paths such as `runs/`, `tools/`, `cavity_agent/` refer to the original working tree, not to this package (see `docs/INDEX.md` for the package equivalents). Later reports supersede earlier ones where they disagree; known numerical/wording inconsistencies and their resolution are listed in `provenance/INCONSISTENCIES.md`. Nothing here claims a verified FEBio result.

# CAVITY EXACT BENCHMARK CONTRACT (cavityFlexibleBottom, FEBio v4.13 campaign)

Status 2026-10-30, written by the CAVITY agent from the *actual* sources (not earlier summaries). No FSI run was made for this document. Sources read:
- solids4foam **development** `1576edc` (fetched; the local `upstream/s4f` checkout is `master` v2.4 and is stale for this case): tutorial `tutorials/fluidSolidInteraction/cavityFlexibleBottom` (README, `system/fluid/blockMeshDict`, `system/solid/blockMeshDict`, `constant/*`, `0/*`, `system/controlDict`) and its `verification/` (README, reference JSON, driver). Copy in `cavity_agent/src_dev/`.
- Manuscript `upstream/pc-paper/sectionBenchmarkSuite.tex` (section "Flexible-bottom cavity [PROVISIONAL]"), `paper-data/cavity_monolithic_618/quasi_monolithic_618.csv`, `paper-data/cavity_plate_material_law*`.
- Paper-repo PR #12 (ANSYS/COMSOL; read through the PR #12 comments), PR #22 (`reference_levels.csv`, `reference_richardson.csv`, copies in `cavity_agent/refs/`).
- Existing FEBio cavity reports: `CAVITY_BENCHMARK_DEFINITION.md`, `CAVITY_PLATE_ONLY_RESULTS.md`, `CAVITY_FLUID_ONLY_RESULTS.md`, `CAVITY_INLET_BOUNDARY_CONDITION_AUDIT.md`, `CAVITY_POISEUILLE_DIAGNOSTIC.md`.

## 1. Locked physical problem
| Item | Value (source) |
|---|---|
| Fluid | rho = 1 kg/m^3, nu = 0.01 m^2/s (mu = 0.01 Pa s), incompressible, laminar (`transportProperties`) |
| Channel | height H = 1 m, y in [1,2], inlet x = 0, outlet x = 14 (`blockMeshDict` vertices 0-25) |
| Inlet | parabolic `parabolicVelocity`, maxValue 1.5 m/s over [1,2] -> mean 1.0 m/s, Re = U H / nu = 100 |
| Outlet | p = 0 (fixedValue), U zeroGradient |
| Fluid walls | no slip (fixedValue 0) incl. the cavity side walls and the plate-end corners |
| Cavity | 4 m wide (x in [2,6]), depth about 0.75-1 m, with rounded upper corners (8 arc block edges, vertices 6-11, 14-15, 18-19); fluid-solid interface = the cavity floor y = 0, x in [2,6] (40 faces at dx = 0.1) |
| Solid | plate x in [2,6], y in [-0.4,0], 40 x 4 cells at dx = 0.1 (`system/solid/blockMeshDict`) |
| Solid law | `StVenantKirchhoffElastic`, E = 500 Pa, nu = 0.4, rho = 1000 kg/m^3, `planeStress no` (plane strain); total-Lagrangian, large displacement |
| Solid supports | `left` (x = 2) and `right` (x = 6) fixedDisplacement 0 over the full end faces; `bottom` (y = -0.4) and `interface` traction (the interface carries the fluid load; bottom is traction free) |
| Out-of-plane | one cell, 1 m thick, front/back `empty` = plane strain; all forces are **per 1 m depth = N/m** |
| Monitor point A | `solidPointDisplacement point (4 0 0.5)`: x = 4, **y = 0 (top, wetted surface)**, z = 0.5 (mid-depth; the 2-D solution is z-independent) |
| Interface force | `forces` function object on the *fluid* `interface` patch: pressure + viscous force of the fluid on the plate, total (Fx, Fy), rhoInf = 1 (kinematic p) |
| Steady target | steady state of a pseudo-transient run (s4f: first-order Euler, dt = 0.5, endTime 200, solid `dampingCoeff` 0.1, coupling starts at t = 10) |

Equivalence statement: the FEBio model must solve this same steady problem. FEBio is monolithic, dynamic (generalised-alpha), weakly compressible (p = -K e_f), so the **steady state is the only comparable quantity**; the route (time step, damping, ramp) is a numerical choice that must be shown not to change the steady state (campaign step C4).

## 2. Verified geometry equivalence (programmatic, `tools/make_cavity_fluid_feb.py` vs the development blockMeshDict)
- 26 vertices (z = 0) identical to 1e-9; z = 1 copies consistent; **15 blocks identical** (vertex order and cell counts); 8 arcs (z = 0) identical (mid-points). Cells at mesh factor 1: 1688 = the tutorial's level-1 fluid count; interface 40 cells (2 + 36 + 2), inlet 10 and outlet 10 faces. Generator-reported fluid area 18.0279 m^2.
- Solid: the tutorial mesh is 160 / 640 / 2560 / 10240 cells at dx = 0.1 / 0.05 / 0.025 / 0.0125; fluid 1688 / 6752 / 27008 / 108032 (`verification/README.md`). The generator's `--mesh-factor f` doubles every block division, which reproduces these fluid counts (x4 per level).
- Plate-only FEBio StVK hex20 (`CAVITY_PLATE_ONLY_RESULTS.md`): u_y(A) within 0.2-0.65 % of the plane-strain Q9 StVK reference at the same load; this supports StVK equivalence but is a structural diagnostic, not an FSI result.

## 3. Reference values (confirmed from the manuscript text and data, not from the brief)
Top-surface point A, steady (manuscript `sectionBenchmarkSuite.tex`, section cavity):
| dx (m) | solids4foam partitioned (Aitken) u_y | quasi-monolithic u_y (rerun #618) | QM F_y (N/m) |
|---|---|---|---|
| 0.1 | -0.2046 | -0.203601 | -5.122626 |
| 0.05 | -0.2314 | -0.231044 | -5.226559 |
| 0.025 | -0.2384 | -0.238200 | -5.250356 |
| 0.0125 | -0.2401 | -0.239986 | -5.255992 |
Extrapolated (p = 2): partitioned -0.2407, QM -0.2406 m. F_y converges to about -5.26 N/m (orders 2.1-2.4). u_x(A) is below 1.1e-4 m (QM 7.4e-5 to 9.7e-5 m). **The brief's numbers are confirmed**; the partitioned level-1 value -0.2046 differs from the -0.2036 listed for "solids4foam" in PR #22 (that row is the QM value), so the two solids4foam variants are kept separate.

Provenance cautions found:
1. The tutorial `verification/README.md` table (u_y -0.20246 / -0.22883 / -0.23563, F_y -5.153 / -5.239 / -5.256) is about 1 % smaller in |u_y| than the manuscript values. Its own note says values recorded before October 2026 were *underside* values (the old monitor point (4,-1,0.5) lies outside the plate, so the nearest underside node was used). Those README u_y values are therefore **not** top-surface references and are not used; the manuscript/PR #22 top-surface values are.
2. Force depth: tutorial = 1 m, so F_y in N is N/m. Tuković et al. quote per 0.05 m depth (x20, confirmed by Ž. Tuković per PR #22). ANSYS used a 10 mm depth (x100). FEBio slab thickness TZ = 0.1 m: **every FEBio force must be divided by TZ** (the generator already writes `TZ`; the post-processing must do it explicitly).
3. The old underside displacement (-0.2025 / -0.2288 at (4,-1)) must not be used (brief, and PR #12 comments).

## 4. Cross-solver reference classes (not merged; used in the later comparison)
- **StVK-confirmed (same constitutive law):** solids4foam partitioned and quasi-monolithic. They share code lineage, so they are a consistency pair, not independent references.
- **Unresolved constitutive difference:** Tuković et al. (2018): u_y(0.0125) = -0.2506, extrapolated about -0.2510 m, +4.1 % vs solids4foam; displacement matches a neo-Hookean-like plate while the force matches StVK to 0.01 % (PR #22). ANSYS (PR #12): -0.2527 (finest), extrapolated about -0.2536 (order about 1), +5.3 %, consistent with a Hencky/log-strain law; not StVK (no built-in StVK). COMSOL (PR #22): plateau u_y = -0.243 +/- 0.001, F_y = -5.295 N/m (+0.7 % load), solid model unconfirmed, not extrapolable.
- **Historical/supplementary:** the tutorial README comparison (-0.232 at dx = 0.1 from a published diagram) and the old underside values.
Plate-only evidence (paper-data, `plate_laws.py`, Q9, uniform 1.315 Pa follower pressure): StVK -0.23887 / -0.24047 (20x2 / 40x4); neo-Hookean +4.8 %, Hencky +5.7 %. These identify the plate law as the main cause of the displacement spread; the force agrees on all references to 0.2 % (COMSOL 0.7 %).

## 5. Equivalence gaps and decisions for FEBio (stated, not hidden)
| Gap | Decision / how it is controlled |
|---|---|
| Incompressible fluid (s4f) vs weakly compressible FEBio | K = 1e6 Pa baseline (Mach ~1.5e-3; Poiseuille test: flux imbalance 1.7e-6); K sensitivity is a planned C4 run |
| Fluid inlet | validated built-in parabolic `fluid normal velocity` + rim pressure (`CAVITY_POISEUILLE_DIAGNOSTIC.md`); the discrete inlet profile is the nodal parabola of flux 1/(1-h^2) at dx = 0.1 (1.0101 x the nominal), so the inlet flux is exact (1 m^2/s) and the nominal-profile mismatch is O(h^2) and shrinks with refinement |
| Odd-even pressure mode (0.0124 Pa in the straight channel, frozen) | measured and reported on raw forces in C1; no filtering is applied in any reported number; any filtered value is labelled as such |
| Plane strain / depth | one-element slab, `wz` and `uz` fixed on all nodes, forces / TZ |
| Solid element | FEBio fluid-FSI uses a conforming shared-node interface: the hex8 solid must match the hex8 fluid interface nodes (40 x 4 hex8 at dx = 0.1). The plate-only test used hex20 (a different element family from the coupled model); a hex8 structural check precedes the pilot because hex8 is stiffer in bending |
| Pseudo-transient route and damping | generalised-alpha dynamics with numerical dissipation; the plate's lowest dry period is about 50 s and its fluid damping is weak, so a start-up route that settles is needed; the steady state must be shown independent of it (C4) |
| Material law | StVK only (`St.Venant-Kirchhoff`); any other law would be a separately labelled experiment and is **not** planned |
| Reference status | none of the available references is both StVK-confirmed and independent of solids4foam; this limits what "agreement" can mean and is carried into the final comparison |

## 6. Gate C0
**PASS (benchmark definition locked; fluid geometry identical to the development tutorial).** Open, carried forward: (i) hex8 vs hex20 solid; (ii) start-up/damping independence; (iii) reaction-force logging gave zeros in the solid module, so the plate equilibrium check needs another route (plotfile reaction forces or an end-face stress integral); (iv) references are not independent of solids4foam for StVK. Coupled FSI is not run until C1 passes.

## 7. Status update at the end of the campaign (2026-10-30)
- Open item (i) hex8 vs hex20 solid: resolved by S1 (hex8 plates 40x4/40x8/80x8/160x16: u_y(A) -0.23026/-0.23467/-0.23789/-0.24023 vs Q9 extrapolation -0.2410, order about 1.7): the conformal hex8 plate is 4.5 %/1.3 %/0.3 % too stiff at the three coupled levels.
- (ii) start-up/damping independence: partly tested (dt/2: 4e-5 in u_y, 2e-5 in F_y; K, tolerance: nil). The ramp length, rhoi and E_mesh were not varied.
- (iii) reaction forces: still zero in the solid-module log; the plate equilibrium check was replaced by S2 (plate-only under the FSI interface pressure reproduces the coupled u_y to 0.03 %).
- (iv) references: still no independent StVK reference. The equivalence gate for the force is **not met**: FEBio's converged F_y is about -5.00 N/m versus -5.26 N/m (see `CAVITY_MESH_STUDY.md`, `CAVITY_CROSS_SOLVER_COMPARISON.md`).
