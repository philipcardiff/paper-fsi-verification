> **ARCHIVED WORKING REPORT (verbatim except path/host sanitisation).** Written during the investigation; file paths such as `runs/`, `tools/`, `cavity_agent/` refer to the original working tree, not to this package (see `docs/INDEX.md` for the package equivalents). Later reports supersede earlier ones where they disagree; known numerical/wording inconsistencies and their resolution are listed in `provenance/INCONSISTENCIES.md`. Nothing here claims a verified FEBio result.

# CAVITY BENCHMARK EQUIVALENCE REVIEW (coordinator, 2026-10-30)
Question: does the FEBio model solve the same boundary-value problem as the manuscript's solids4foam case (`CAVITY_EXACT_BENCHMARK_CONTRACT.md` is the baseline)? Evidence: the generated `.feb` files (`tools/make_cavity_fsi_feb.py`), the contract, `CAVITY_PRESSURE_FORCE_AUDIT.md`, `CAVITY_MESH_CONVERGENCE.md`.

| Item | Manuscript / solids4foam | FEBio | Verdict |
|---|---|---|---|
| Fluid geometry | 15 blockMesh blocks, 26 vertices, 8 arcs | programmatically identical (agent-checked to 1e-9); cell counts 1688/6752/27008 match | equivalent |
| Fluid properties | rho 1, nu 0.01, Re 100 | rho 1, mu 0.01 | equivalent |
| Compressibility | incompressible | K = 1e6 (K x10: 1e-6 effect) | equivalent to 1e-6 |
| Inlet | parabolic, mean 1 | validated `fluid normal velocity`; discrete flux 1.0000 | equivalent (profile 1.0101 x nominal at dx 0.1, O(h^2)) |
| Outlet | p = 0, U zero-gradient | zero dilatation (p = 0) all 42 nodes, velocity natural | same datum; natural-vs-zeroGradient treatment untested |
| Pressure constraints | none other | **4-node corner pin 1.15 Pa** (value from FEBio's own C1) | **not equivalent** (added constraint; effect 0.2-0.4 % of F_y) |
| Plate | StVK E 500, nu 0.4, rho 1000, plane strain, clamped ends | same law; hex8 (not Q9) 3-D slab with u_z = 0; same clamps | equivalent in law, hex8 stiffer (0.3 % at 160 x 16) |
| Underside | traction free | traction free | equivalent |
| Interface | no-slip, conformal, total fluid traction | no-slip shared nodes, `fluid-FSI traction` | equivalent; pressure part verified (S2); viscous/total transfer not tested |
| Dynamics | pseudo-transient Euler, damping 0.1, coupling from t = 10 | generalised-alpha (rhoi 0), 25 s ramp | route differs; dt halving and tolerances change the steady state by 4e-5 / 0 (C4, M1) |
| Output definition | `forces` on interface, N/m (1 m depth), A on top at (4,0) | pressure + viscous, / TZ, A at (4,0) top | equivalent (verified from the fields) |
| Solver completion | normal | all three levels error-terminated on the dilatation-increment floor | **not equivalent in quality of completion** |

**Verdict:** the FEBio model matches the manuscript benchmark in geometry, material law, fluid parameters, inlet, outlet datum and output definitions, with two documented differences (the 4-node pin, the equal-order weakly compressible hex8 discretisation) and one quality issue (error terminations). The remaining ~5 % force/displacement difference is **not explained by any BC that can be found different**; it is therefore not yet attributable to a defect in either code. Not established: that solids4foam's own floor pressure is correct to 5 % in this configuration (its three variants agree to 0.2 % but share a lineage), that the outlet natural condition equals zeroGradient U, that FEBio's equal-order scheme converges to the incompressible solution at this Re without a persistent offset.
