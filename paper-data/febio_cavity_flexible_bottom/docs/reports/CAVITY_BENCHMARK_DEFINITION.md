> **ARCHIVED WORKING REPORT (verbatim except path/host sanitisation).** Written during the investigation; file paths such as `runs/`, `tools/`, `cavity_agent/` refer to the original working tree, not to this package (see `docs/INDEX.md` for the package equivalents). Later reports supersede earlier ones where they disagree; known numerical/wording inconsistencies and their resolution are listed in `provenance/INCONSISTENCIES.md`. Nothing here claims a verified FEBio result.

# CAVITY BENCHMARK DEFINITION (cavityFlexibleBottom): read-only audit

Status: read-only audit, no FEBio case written, no cavity simulation run. Sources inspected 2026-10-29:
- Manuscript `philipcardiff/paper-fsi-verification` main `863b8ba` (`git fetch` shows no newer commit): `sectionTestCases.tex` (Case 1), `sectionBenchmarkSuite.tex` (§ Flexible-bottom cavity [PROVISIONAL]), `paper-data/cavity_plate_material_law_findings.md`, `paper-data/cavity_plate_material_law/plate_laws.py`, `paper-data/cavity_monolithic_618/README.md`.
- solids4foam `development` `1341312` tutorial `tutorials/fluidSolidInteraction/cavityFlexibleBottom` (README, `system/*/blockMeshDict`, `0/*`, `constant/*`, `controlDict`).
- Manuscript PRs (state as of today): **#12** OPEN (Ansys + COMSOL results, E. Tandis), **#22** OPEN (COMSOL and Tuković references, Richardson extrapolation, `evidence/cavity-comsol-richardson`), #14/#18/#20/#23 merged. **#25** is the user's draft PR (radialTorsionalRing). #22's `reference_levels.csv` and `reference_richardson.csv` were fetched read-only from its branch into `cavity/reference/` (provisional until #22 merges).
- FEBio v4.13 source (`upstream/febio-v4.13`).

## 1. Definition (as the current manuscript/tutorial state it)

| item | value | source |
|---|---|---|
| Dimensionality | 2-D plane strain (`planeStress no`), 1 m depth (z = 0..1, `empty` front/back) | tutorial `mechanicalProperties`, solid/fluid blockMeshDict |
| Fluid domain | channel x in [0, 14], y in [1, 2] (H = 1 m) with a cavity below: x in [2, 6], y in [0, 1]; top edges of the cavity are rounded (radius ~0.25 m arcs at x ~ 2 and 6, y 0.75-1.0); upper wall y = 2 | fluid blockMeshDict; geometry figure |
| Solid | rectangle x in [2, 6], y in [-0.4, 0]: 4 m x 0.4 m (L/h = 10), forms the cavity floor | solid blockMeshDict |
| Supports | `left` (x = 2) and `right` (x = 6) faces: `fixedDisplacement` 0 over the full thickness (clamped both ends); `bottom` (y = -0.4) traction free; interface = top face y = 0 | `0/solid/D` |
| Constitutive law | **Saint Venant-Kirchhoff**, E = 500 Pa, nu = 0.4, rho_s = 1000 kg/m^3 (lambda = 714.29, mu = 178.57, K = 833.33 Pa) | tutorial, manuscript Table |
| Solid damping | mass-proportional `dampingCoeff` 0.1 /s (pseudo-transient acceleration only) | `solidProperties` |
| Fluid | incompressible Newtonian, rho_f = 1 kg/m^3, nu_f = 0.01 m^2/s (mu = 0.01 Pa s), Re = 100 on H and mean inlet speed | tutorial |
| Inlet | parabolic u_x = 1.5 - 6 (y - 1.5)^2 on y in [1, 2] (mean 1 m/s, max 1.5), u_y = 0 | `inlet.fixedProfile` / `parabolicVelocity maxValue 1.5` |
| Outlet | x = 14: p = 0 (fixed), velocity zero gradient | `0/fluid/p`, `U` |
| Walls | no-slip on channel walls and cavity side walls | `U` |
| Interface | no-slip (fluid velocity = plate velocity); fluid pressure and viscous traction act on the plate; solids4foam moves the fluid mesh to the interface | `fsiProperties` |
| Time integration / stopping | transient pseudo-time to steady state: Euler (first order) both solvers, dt = 0.5 s (tutorial; the paper's sensitivity study varies dt 0.5-20 s), tutorial end 200 s, paper runs to the automatic steady criterion (relative change of u_y and F_y), steady at ~120 s; finest level needed a smaller dt | tutorial README, manuscript |
| Coupling tolerance | partitioned: Aitken, outer tolerance 1e-5, up to 20 outer iterations (strict test used for paper results) | `fsiProperties`, manuscript |
| Quantities of interest | steady u_y at **point A on the top (wetted) surface at (4, 0)**; integrated interface force F_y per unit depth | manuscript |
| Mesh levels | uniform refinement dx = 0.1, 0.05, 0.025, 0.0125 m (solid 40x4 at dx = 0.1) | manuscript / solid blockMeshDict |

## 2. Reference values and their uncertainty (`cavity/reference/reference_levels.csv`, `reference_richardson.csv`; provisional, from PR #22 and PR #12)

| source | u_y at A, finest dx = 0.0125 (m) | extrapolated u_y (m) | F_y finest (N/m) | extrapolated F_y | note |
|---|---|---|---|---|---|
| solids4foam (StVK, partitioned) | -0.2401 | -0.2406 (p ~ 2.15) | -5.256 | -5.2574 | monotone, observed order 1.97-2.15 |
| Tuković et al. 2018 | -0.25062 | -0.2510 (p ~ 1.8) | -5.2548 | -5.2567 | read from published figure vector paths; force published per 0.05 m depth, rescaled; law apparently neo-Hookean-like |
| Ansys | -0.25274 | ~ -0.2536 (p ~ 0.9-1.2, not asymptotic) | -5.263 | -5.272 | material law unconfirmed (behaves like Hencky/neo-Hookean) |
| COMSOL | -0.243 | plateau -0.243 +/- 0.001 | -5.2952 | not extrapolable | F_y offset +0.7 %, solid model unconfirmed |

- Integrated force agrees across sources to 0.01-0.2 % (solids4foam vs Tuković) and to 0.2-0.7 % (Ansys, COMSOL), so the fluid load is common.
- Displacement differs by up to 5 %: the plate-only check in the manuscript data (uniform follower pressure 1.315 Pa, 20x2/40x4 Q9) gives u_y ratios to StVK of 1.048 (neo-Hookean), 1.057 (Hencky), 1.05 (solids4foam-split neo-Hookean), consistent with the reference differences. **The accepted reading in the manuscript is: the benchmark specifies StVK; Tuković and Ansys displacements are consistent with other large-strain laws.** This was not independently verified here.
- There is no single settled reference value. The appropriate comparison target for a StVK FEBio run is the solids4foam value (-0.2406 m, -5.257 N/m) and, for the force, all four sources; displacement against Tuković/Ansys is expected to differ by ~4-5 % by law, not by solver.

## 3. Ambiguities and risks found in the audit (do not resolve by tuning)

1. **Point A**: the tutorial monitors (4, -1, 0.5), outside the solid (plate bottom is at y = -0.4), so it records the nearest underside node; the paper corrects to the top surface (4, 0). The top/mid/underside values differ by ~0.5-1.2 %. FEBio must sample the **top-surface** node at x = 4.
2. **Depth normalisation**: solids4foam F_y is per 1 m depth; Tuković per 0.05 m (rescaled by 1/0.05); Ansys per 10 mm slab (x100). FEBio slab of thickness t must be reported as F_y / t. The manuscript lists the force as "N/m".
3. **Plane strain vs plane stress** appears in the manuscript text as a candidate explanation for the displacement discrepancy (factor 1/(1-nu^2) = 1.19), but the tutorial states `planeStress no`. The reference differences (4-5 %) are smaller than that factor and the plate-only study attributes them to the law; plane strain is the definition to use.
4. **Fluid incompressibility**: FEBio fluid-FSI is weakly compressible (p = -K e_f). Steady Re = 100 with ~1 Pa dynamic pressure needs K >> rho U^2; K will have to be demonstrated converged (a K sensitivity step as in ringAddedMass), not assumed.
5. **Outlet and inlet**: solids4foam outlet is p = 0, zero-gradient velocity; in FEBio the corresponding natural condition (zero dilatation = p 0 with traction-free normal velocity) is to be checked. A parabolic inlet needs a velocity profile in FEBio (math or per-node load curve); the per-node workaround from the ring campaign avoids the v4.13 fluid-IC math defect, but the velocity BC (not IC) path has not been tested for spatial expressions.
6. **Steady state target**: s4f accelerates via solid mass damping (0.1 /s). FEBio has no equivalent in the fluid-FSI solid in v4.13 as far as inspected; the steady state itself is independent of damping, but the pseudo-transient path (physical density 1000, dt, rhoi < 1 numerical dissipation) must be shown to reach the same steady state with a defined stopping criterion. **Untested.**
7. **Rounded cavity corners (arcs)**: need curved/graded meshing in the FEBio generator; the conformal fluid-solid interface at x = 2 and 6 meets sharp fluid corners at the clamped ends (singular stress at the clamped corners; same in all codes).
8. **Same-lineage warning**: Tuković (2018) and the solids4foam tutorial share authors, so agreement between them is not independent; the independent references are Ansys and COMSOL, both with unconfirmed solid models and not yet converged.

## 4. Constitutive law in FEBio 4.13

- Material type `St.Venant-Kirchhoff` exists (`FEBioMech/FEStVenantKirchhoff.cpp`, registered in `FEBioMechModule.cpp`), parameters `E`, `v`.
- Stress in the source (Cauchy): sigma = [lambda trE b + mu (b^2 - b)] / J with trE = (tr b - 3)/2. This equals F S F^T / J for S = lambda trE I + 2 mu E, E = (C - I)/2, i.e. the standard StVK law. Same lambda, mu as the tutorial, so the law is equivalent in form; finite-deformation kinematics are total Lagrangian with large strain in FEBio's nonlinear solid (not a small-strain code).
- Inside FEBio fluid-FSI the solid is a separate standard solid domain (as in the ring), so StVK can be used for the plate; the regulariser applies only to the fluid domain. Neo-Hookean and a Hencky-like law are also available for the law-sensitivity question, but **the principal baseline is StVK, no material tuning**.
- Equivalence of the implemented stress to the target law has been read from source, not yet verified numerically; the plate-only test (section 5) is the numerical check.

## 5. Proposed B0 test (not run)

Plate-only, FEBio, plane strain (one-element-thick slab, u_z = 0), x in [2, 6], y in [-0.4, 0], both ends clamped, StVK E = 500, nu = 0.4, uniform follower pressure 1.315 Pa on y = 0 (FEBio `pressure` surface load follows the deformed surface), static, load-stepped. Compare u_y at (4, 0) with the documented table (`plate_laws.py`: StVK 20x2 Q9 -0.23887, 40x4 -0.24047; neo-Hookean 1.048 ratio; linear -0.28816) and with the linear beam-with-shear estimate. Acceptance in terms of the ratio between laws and the StVK value within the plate-only mesh sensitivity (~0.7 % absolute); recorded before running.

## 6. Feasibility-relevant FEBio gaps found so far
- None found that blocks starting B0.
- Open for B1: curved-geometry meshing, spatial inlet velocity BC, steady-state strategy and K convergence (items 4-7 above).
- No FEBio native fluid-FSI cavity with flexible floor exists in the TestSuite per the Phase 0 audit (`PHASE0_AUDIT.md` line 163 lists only Couette, wedge, lid-driven cavity, tube cases).
