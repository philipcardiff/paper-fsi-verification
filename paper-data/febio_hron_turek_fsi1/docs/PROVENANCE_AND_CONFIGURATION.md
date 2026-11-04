# Provenance, exact configurations, mesh generation, acceptance criteria

## 1. Binaries and sources
| item | value |
|---|---|
| FEBio | v4.13, upstream git tag `v4.13`, commit `32ae206ff4881dfb54f62296cd1558e58ed9fcc6`, working tree clean (no local modification); built in user space with Intel MKL (Pardiso) |
| stock binary used for all 19 jobs | `febio4`, sha256 `f5ae1ed6f347f875c392f93d847f2369f437345a9b785296d46b0bd6fbb72f02` (recorded by `run_meta_start.json` of every job: all 19 identical) |
| not used | any plugin, the experimental pressure-projection build, any stabilisation or source change |
| solids4foam tutorial audited | `development`, commit `b200718f4f8390c90f9c7e3ded8b74f3e13941cc` (2026-10-30), `tutorials/fluidSolidInteraction/TurekHron2006`; the fluid blockMeshDict copy used for the mesh topology is `scripts/hronturek_data/TurekHron2006_fluid_blockMeshDict.s4f_b200718`, sha256 `88359ea97ad4b9167f92d6d40b450dae969470b848e2b2bb9cd6ae203db4e637` |
| references | Featflow FSI benchmark pages fetched 2026-10-31: `fsi_definitions.html` (geometry), `fsi_csm_tests.html` (CSM1/CSM2), `fsi_cfd_tests.html` (CFD1), `fsi_fsi_tests.html` (FSI1); repo copy `HronTurek_verification_references.json` of the tutorial agrees digit for digit with the FSI1 table. The Turek-Hron 2006 paper text could not be retrieved: any statement attributed to it is UNRESOLVED |
| hardware/software | shared Slurm node, 2 or 4 cores per job, `OMP_NUM_THREADS = MKL_NUM_THREADS = cores`; generators and analysis: Python 3 with numpy (matplotlib for figures) |
Script hashes (sha256, first 16 hex digits) at the time of packaging:

| script | sha256 |
|---|---|
| `scripts/check_inputs.py` | 51deb99d5d083790 |
| `scripts/compare_json.py` | 3371b8a16ee1fb42 |
| `scripts/equilibrium_from_xplt.py` | 860049319362e758 |
| `scripts/hronturek_bc_driver.py` | c07ade9a1c3709ff |
| `scripts/hronturek_cfd1_independent.py` | 518b146c71448662 |
| `scripts/hronturek_cfd_analyze.py` | a47d54215fa1ea88 |
| `scripts/hronturek_cfd_diag.py` | 705450c2ff510d0f |
| `scripts/hronturek_csm_analyze.py` | 67e6a5f0916322da |
| `scripts/hronturek_fsi_traction_recon.py` | d247edfd4e75d53e |
| `scripts/make_figures.py` | a166b434df39c2ed |
| `scripts/make_force_budget.py` | b5f64544281b0d1f |
| `scripts/make_hronturek_cfd_feb.py` | 29239afc7b0c1d1d |
| `scripts/make_hronturek_csm_feb.py` | 1e2a44d18128c974 |
| `scripts/make_stageB_tables.py` | 9dce40e1726ce77e |
| `scripts/parity_diagnostics.py` | 3fef90e647364758 |
| `scripts/subset_json.py` | 50f3adc6ed1f92d2 |
| `scripts/xplt_reader.py` | 403978df5d3a9fe2 |
`make_hronturek_cfd_feb.py` v1 generated the unpinned input C1u, v2 added `--pin-root-values` (the v2 no-pin output is bit-identical to the v1 file: checked by `reproduce.sh`). `make_hronturek_csm_feb.py` v2 writes the body-force sign convention and the `reaction forces` plot variable (see `STAGE_B_STRUCTURAL.md`); the v1 input of the diagnostic job is `inputs/stageB/DIAG_v1_signtest_B_CSM1_h8_nt4.feb`.

## 2. Exact model configuration
### Benchmark geometry (all stages)
Channel 2.5 x 0.41 m; cylinder centre (0.2, 0.2), r = 0.05; flag thickness 0.02 m, tip at x = 0.6 m, flag root on the cylinder circle at x = 0.2 + sqrt(0.05^2 - 0.01^2) = 0.2489898 m (y = 0.19 and 0.21), free length 0.351010 m; point A = (0.6, 0.2) (tip, mid-thickness). The Featflow definition page (fetched 2026-10-31) states: bar length l = 0.35, height h = 0.02, right-bottom corner at (0.6, 0.19), "the left end is fully attached to the fixed cylinder", A(0) = (0.6, 0.2). The nominal l = 0.35 puts the left end at x = 0.25, which at y = 0.19/0.21 is 1.0 mm outside the circle (x = 0.24899): an end at x = 0.25 would not touch the cylinder, so "attached to the cylinder" is read as a root on the circle, which is also where the tutorial's blockMeshDict puts the flag (vertices (248.9898, 190) and (248.9898, 210) mm). This reading was fixed before any run; the x = 0.25 variant was run only as a labelled diagnostic. Which root Featflow used numerically is UNRESOLVED (the agreement of the G1 variant, 0.2 % against 1.3 %, supports but does not prove it).
### Stage B (flag only)
FEBio `solid` module, `STATIC`, 10 linear load steps, full Newton (dtol 1e-8, etol 1e-12, rtol 0, lstol 0.9, max_refs 25, symmetric stiffness). Material `St.Venant-Kirchhoff`, E = 1.4e6 (CSM1) or 5.6e6 (CSM2), v = 0.4, density 1000. One-element-thick slab (TZ = 0.015 m), u_z = 0 on all nodes (plane strain), root face fully clamped (u_x = u_y = 0 on all root nodes). Body force: `<body_load type="const">` with y written as +2 (see the sign note: FEBio applies this with the opposite sign to the stated acceleration; the physical g is -2 m/s^2, the flag deflects down). Mesh: nt elements across the thickness, nx = round(L / (0.02/nt)) along the length (aspect ratio 1.00): hex8 nt = 2, 4, 8, 16 (nx = 35, 70, 140, 281; 70, 280, 1120, 4496 elements) and hex20 on the nt = 2, 4 grids; informational: root at 0.25 (nt 16), linear load g = -0.02 (nt 8, 16). nu = 0.4 is not near-incompressible (lambda/mu = 4); no mixed or reduced-integration formulation was used.
### Stage C (rigid flag and cylinder)
FEBio `fluid` module (fixed mesh), `DYNAMIC`, hex8, one-element-thick slab TZ = 0.015 with w_z = 0 on all nodes. Fluid material: density 1000, Newtonian mu = 1.0 Pa s (nu = 1e-3), kappa 0, bulk modulus K = 1e6 Pa (p = -K e_f), rhoi = 0. Domain: the tutorial's 24-block topology (cylinder O-grid, rigid flag excluded from the fluid), counts x 1.25 = 8,330 cells, 17,240 nodes per layer; area within 2.1e-5 of exact, minimum corner sine 0.892, maximum edge aspect ratio 17.4, first-cell wall height 1.6 / 4.5 / 9.1 mm (min / median / max). Boundary conditions: no-slip (`zero fluid velocity`) on cylinder, flag, top, bottom and the inlet/outlet wall-corner nodes; inlet `fluid normal velocity`, velocity -0.2 (mean speed), `parabolic` map, `prescribe_nodal_velocities`, `prescribe_rim_pressure`: the discrete inlet flux is exactly 0.2 x 0.41 g(t) (verified 3e-15), the nodal profile is 1.0006 x the exact parabola (peak 0.3002 m/s); outlet `zero fluid dilatation` on every outlet node (p = 0). Pseudo-transient to steady state: g(t) = (1 - cos(pi t/2))/2 for t <= 2 s then 1, dt = 0.5 s, 80 steps (t = 40 s, about 3 flow-through times); the CFD1 reference is steady, the transient is only FEBio's route to it. Solver: full Newton, reform each step, vtol = ftol = 1e-6, etol = rtol = 0, max_refs 15, lstol 0.9 (the declared relaxed-test fallback was not needed). C1u: unpinned. C1p: the two flag-root junction nodes (x = 0.2489898, y = 0.19 and 0.21, both z layers) carry `prescribed fluid dilatation` of p = -81.6042 Pa (y = 0.19) and -123.8709 Pa (y = 0.21), ramped with g(t); the values are the means of the two nearest wall neighbours of C1u's final state (pre-declared rule), uncertain by 10-40 Pa because those neighbours are themselves affected by the parity mode.
### Mesh generation and re-running
```
python3 scripts/make_hronturek_csm_feb.py --nt 16 --elem hex8 --E 1.4e6 --out RUN_DIR --tag B_CSM1_h8_nt16        # stage B (CSM1); --E 5.6e6 for CSM2; --root 0.25; --g -0.02
python3 scripts/make_hronturek_cfd_feb.py --out RUN_DIR --tag C1p_L1_f1.25_K1e6 --factor 1.25 --pin-root-values=-81.6042241147,-123.87086089350001   # stage C; omit the pin option for C1u
febio4 -i <tag>.feb -silent   # FEBio 4.13 + MKL, 2 (B) or 4 (C) threads; the log <tag>.log and the node_data files are written next to the input
```
`model_meta.json` written by each generator holds all arguments and the sha256 of the .feb. `scripts/hronturek_bc_driver.py` is the Slurm/ledger driver of the original workspace (site-specific, not needed to reproduce the numbers).

## 3. Acceptance criteria defined before the runs (summary; full text in `docs/archive/HRON_TUREK_FSI1_STAGE_BC_RESULTS.md` sections 1.4 and 3.4)
Stage B (finest hex8, root on the circle, both CSM1 and CSM2): PB1 u_y within 1 % of the reference; PB2 u_x within 2 % (second-order quantity: a 1 % u_y error implies about 2 % in u_x; written before running); PB3 all load steps converged, mean iterations <= 8, direction correct, max |u_z| = 0, and a reaction-force equilibrium check (NOT ACHIEVABLE, see `STAGE_B_STRUCTURAL.md`; replaced by independent checks and, in this package, a stress-based free-body check); PB4 monotone mesh behaviour with the observed order reported; PB5 hex20 nt 4 within 0.5 % (diagnostic).
Stage C: C1 mass balance <= 1e-3; C2 steady-state spread over the last 20 % of the run <= 1e-3 (drag), 3e-3 (lift); C3 drag within 3 % at level 1 and 1 % at level 2 with decreasing error (gate on F1, the solver-applied-equivalent force); C4 force consistency F1 vs F2 vs F3 <= 1 % of the drag (3 % of the lift); C5 pressure odd-even amplitude <= 0.048 Pa on both flag faces and a reduction >= 1.5 between levels; C6 normal termination, mean Newton iterations <= 8; C7 K x 10 and pin sensitivities (not run); C8 one-way structural prediction (not run). The first failed gate stops stage C (declared before the first run); level 1 failed C3, C4 and C5, so level 2 and C7/C8 were not run.
Corrections to the written criteria made before the affected results were interpreted: (i) the sign of the F3 formula in the pre-registration text was reversed (the code was corrected before the numbers were interpreted); (ii) the stage-B reaction check could not be performed (see `STAGE_B_STRUCTURAL.md`).
