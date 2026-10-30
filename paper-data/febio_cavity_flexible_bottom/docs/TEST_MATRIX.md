# Test matrix (one block per test)

Generated from `results/TEST_MATRIX.csv` (the machine-readable master; same content). Evidence labels: PROVEN / STRONGLY SUPPORTED / HYPOTHESIS / FAILED / NOT TESTED. 'not recorded' / 'not applicable' mean what they say; no value is invented. A PASS means a pre-registered or stated numerical criterion was met, never that the benchmark is verified.

| id | group | test | pass/fail | evidence |
|---|---|---|---|---|
| A1 | A physical model | Fluid and plate geometry equals the manuscript / solids4foam tutorial | PASS | PROVEN |
| A2 | A physical model | Solid material | PASS | STRONGLY SUPPORTED |
| A3 | A physical model | Depth normalisation (plane strain, forces per metre) | PASS (FEBio side); COMSOL convention unresolved | PROVEN (FEBio side) |
| A4 | A physical model | Fluid properties and weak compressibility | PASS (M1 only) | PROVEN (M1) |
| A5 | A physical model | Inlet and outlet conditions | PASS | PROVEN (inlet, straight channel C1); STRONGLY SUPPORTED (outlet equivalence) |
| A6 | A physical model | No-slip conformal FSI interface | PASS (pressure part) | PROVEN (pressure part, vertical) |
| A7 | A physical model | Pressure constraints (added numerical constraints) | PASS (bound) | PROVEN (operator and linear attribution); bound is an estimate |
| A8 | A physical model | Monitoring point A on the wetted top surface | PASS | PROVEN |
| B1 | B structural | Plate-only StVK check, hex20 20x2 | PASS | STRONGLY SUPPORTED |
| B2 | B structural | Plate-only StVK check, hex20 40x4 | PASS | STRONGLY SUPPORTED |
| B3 | B structural | Constitutive-law checks: why displacement references spread by 4-5 % | n/a (context) | STRONGLY SUPPORTED |
| B4 | B structural | hex8 plate refinement (element family of the coupled plate) | PASS (converges); hex8 stiffer at coupled levels (4.5 %, 1.3 %, 0.3 % at 40x4, 80x8, 160x16) | PROVEN |
| B5 | B structural | S2: pressure-transfer / displacement compliance check | PASS | PROVEN for the pressure-part self-consistency |
| B6 | B structural | Displacement/force compliance ratio | consistent at the finest levels | STRONGLY SUPPORTED |
| B7 | B structural | Plate reaction-force equilibrium | NOT EVALUABLE | NOT TESTED |
| B8 | B structural | Remaining structural differences | PASS (cannot explain the deficit) | STRONGLY SUPPORTED |
| C0 | C fluid-only | Early fluid-only attempts with nodal Dirichlet or 'fluid velocity' math inlets | FAIL (preserved) | STRONGLY SUPPORTED (source audit) / HYPOTHESIS (nodal-inlet omits flux term; no K-matched run) |
| C1 | C fluid-only | Poiseuille reference: straight 14x1 channel with the built-in parabolic inlet | PASS (inlet treatment); sawtooth caveat | PROVEN (straight channel) |
| C2 | C fluid-only | C1: fixed-plate fluid-only cavity without corner pin | FAIL (A3, A4) | PROVEN (drift) ; corner mechanism PROVEN by the later operator analysis |
| C3 | C fluid-only | C1b: corrected corner constraints (4 corner nodes pinned to 1.15 Pa) | PASS A3,A5-A8; FAIL A4 as written | PROVEN (steady, flux); the A4 explanation is STRONGLY SUPPORTED (derived from the same data) |
| C4 | C fluid-only | Flux conservation | PASS | PROVEN |
| C5 | C fluid-only | Pressure modes in the fluid-only fields | n/a | PROVEN (existence, in the solver unknown e_f); mechanism STRONGLY SUPPORTED (see F7, F8) |
| C6 | C fluid-only | Momentum and continuity residuals of the converged C1b field (independent python re-assembly of the source element equations) | PASS | PROVEN (for the fluid-only domain of C1b) |
| C7 | C fluid-only | Outlet and corner conditions | PASS | PROVEN (spectrum) |
| D1 | D coupled FEBio | M1 coupled FSI (C2 pilot) | FAIL (abnormal termination); spread 7.6e-5 passes on [875,972.5] | PROVEN (values as stored); conditional steady estimate |
| D2 | D coupled FEBio | M2 coupled FSI | FAIL (abnormal termination); steady spread 1.79e-4 NOT met | PROVEN (values as stored); conditional estimate |
| D3 | D coupled FEBio | M3 coupled FSI | FAIL (abnormal termination); steady spread 1.05e-3 NOT met (criterion 1e-4) | PROVEN (values as stored); NOT steady to 1e-4 |
| D4 | D coupled FEBio | Forces and displacements vs same-dx solids4foam | FAIL (not in agreement) | PROVEN (numbers); cause of the gap NOT established |
| D5 | D coupled FEBio | Solver convergence history and termination reason | FAIL | PROVEN |
| D6 | D coupled FEBio | Steady-state acceptance | M1 met; M2, M3 not met | PROVEN |
| E1 | E numerical sensitivity | Compressibility K x10 | PASS | PROVEN (M1) |
| E2 | E numerical sensitivity | Time step | PASS (independent within about 1e-4) | PROVEN (M1) |
| E3 | E numerical sensitivity | Stopping tolerances | PASS | PROVEN (M1) |
| E4 | E numerical sensitivity | Mesh refinement of the coupled model | FAIL (converges to a value about 4-5 % below the reference) | PROVEN (numbers); INDICATIVE (extrapolation) |
| E5 | E numerical sensitivity | Corner-pin effect on the force | PASS (cannot explain the deficit) | PROVEN (linear attribution); bound is an estimate |
| E6 | E numerical sensitivity | Pressure-force reconstruction F1..F4 at stored states (never mixing times) | PASS (reconstruction); F2 is NOT established as the physical force | PROVEN (arithmetic decomposition); which value is physical is NOT ESTABLISHED |
| E7 | E numerical sensitivity | Viscous contribution (2-point wall-gradient estimate) | PASS (bound) | STRONGLY SUPPORTED (order of magnitude only) |
| E8 | E numerical sensitivity | Parameters never varied | NOT TESTED | NOT TESTED |
| F1 | F source analysis | Source version | PASS | PROVEN (tag); STRONGLY SUPPORTED (binary built from it) |
| F2 | F source analysis | Source-tree cleanliness | PASS | PROVEN (at check time) |
| F3 | F source analysis | Binary provenance | NOT PROVEN bit-exactly | STRONGLY SUPPORTED |
| F4 | F source analysis | Pressure interpolation | n/a | PROVEN |
| F5 | F source analysis | Dilatation degrees of freedom | n/a | PROVEN |
| F6 | F source analysis | Momentum operator | PASS | PROVEN |
| F7 | F source analysis | Row-parity pressure mode (Schur analysis) | n/a | STRONGLY SUPPORTED (soft inf-sup-deficient mode of the unstabilised equal-order pair); mechanism chain HYPOTHESIS |
| F8 | F source analysis | Mesh-dependent eigenvalue scaling | n/a | PROVEN for the two levels studied |
| F9 | F source analysis | Assembled residual reconstruction | PASS | PROVEN (fluid-only domain) |
| F10 | F source analysis | Interface traction assembly (two pathways) | PASS | PROVEN (what each routine does); STRONGLY SUPPORTED (F1 contains the parity offset, F2 does not) |
| F11 | F source analysis | Search for interior pressure stabilisation in v4.13 | n/a | PROVEN (v4.13 tree searched) |
| F12 | F source analysis | Row parity independent of mesh motion and of coupling | n/a | STRONGLY SUPPORTED |
| G1 | G stabilisation round 1 (RESEARCH PROTOTYPE) | Element-level tests of the polynomial-pressure-projection (Dohrmann-Bochev) term | PASS | PROVEN (element level) |
| G2 | G stabilisation round 1 (RESEARCH PROTOTYPE) | alpha=0 control on C1b (does the experimental build reproduce stock?) | PASS (5/5 pre-registered criteria) | PROVEN (fluid-only C1b, 4 threads); SUPPORTED for the FSI domain (5-step start-up only) |
| G3 | G stabilisation round 1 (RESEARCH PROTOTYPE) | Tiny exact Poiseuille channel, alpha=3/4 | mass conservation deteriorates (as predicted) | PROVEN (this case) |
| G4 | G stabilisation round 1 (RESEARCH PROTOTYPE) | C1b fixed-plate fluid-only with alpha=3/4 (one predeclared alpha) | FAIL (S3, S5, S6); S2 mode suppression, S4, S7, S1, S8 pass | FAILED (as a correction); mode suppression PROVEN for this case |
| G5 | G stabilisation round 1 (RESEARCH PROTOTYPE) | Diagnosis of the Newton stall at alpha=3/4 (lagged rim pressure) | initial build FAILED; remedy PASS | PROVEN (tiny channel) |
| G6 | G stabilisation round 1 (RESEARCH PROTOTYPE) | FSI domain (coupled M1 start-up, 5 steps) | PASS (alpha=0); alpha>0 untested beyond 5 steps | SUPPORTED (5 steps) |
| G7 | G stabilisation round 1 (RESEARCH PROTOTYPE) | Stage 5: one coupled M1 run with the stabilisation; M2 fluid pair | NOT RUN | NOT TESTED |
| H1 | H stabilisation round 2 (RESEARCH PROTOTYPE) | Assembled tests of the flux-preserving (redistribution) form | PASS | PROVEN (machine precision) |
| H2 | H stabilisation round 2 (RESEARCH PROTOTYPE) | P0: rim-once at alpha=0 reproduces stock | PASS | PROVEN (this case) |
| H3 | H stabilisation round 2 (RESEARCH PROTOTYPE) | Tiny exact Poiseuille channel, alpha=1e-3 (flux-preserving form; alpha declared from mu x lowest smooth Schur eigenvalue, not from forces) | FAIL (3 of 6 + P0 control passes) | FAILED |
| H4 | H stabilisation round 2 (RESEARCH PROTOTYPE) | Fluid-only cavity test with the round-2 form | NOT RUN | NOT TESTED |

## A1: Fluid and plate geometry equals the manuscript / solids4foam tutorial

- **input configuration**: Channel x 0-14 m, y 1-2 m; cavity x 2-6 m, depth about 1 m with rounded upper corners; plate x 2-6, y -0.4-0 m
- **source binary provenance**: not applicable (generator check)
- **mesh timestep**: level-1 fluid 1688 cells = tutorial level 1; counts 6752/27008 at 2x/4x
- **cpu allocation**: login node, seconds
- **observed result**: 26 vertices identical to 1e-9, 15 blocks and 8 arcs identical; fluid cell counts 1688/6752/27008 and plate 160/640/2560 equal the tutorial levels
- **numerical acceptance criterion**: programmatic equality of vertices, blocks, arcs, cell counts
- **pass fail**: PASS
- **strength of evidence**: PROVEN
- **remaining limitations**: Compared with the solids4foam development tutorial blockMeshDict (commit 1576edc, maintainers' repository), not with a solids4foam run
- **source file or csv**: docs/reports/CAVITY_EXACT_BENCHMARK_CONTRACT.md sec. 2; scripts/make_cavity_fluid_feb.py
- **reproducible from retained data**: partly: generator and hashes retained; tutorial blockMeshDict not copied (needs the solids4foam repository)

## A2: Solid material

- **input configuration**: St.Venant-Kirchhoff E=500 Pa, nu=0.4, rho=1000 kg/m3, large-displacement total-Lagrangian solid, plane strain by u_z=0 on all nodes
- **source binary provenance**: stock febio4 sha256 f5ae1ed6...72f02 (v4.13, source tree clean; no bit-exact binary proof)
- **mesh timestep**: hex8 plate 40x4 / 80x8 / 160x16 cells (M1/M2/M3)
- **cpu allocation**: not applicable
- **observed result**: Same law and parameters as the tutorial (planeStress no); stress expression read from FEStVenantKirchhoff.cpp equals F S F^T/J
- **numerical acceptance criterion**: parameter and source-formula equality; numerical check by plate-only tests B1-B4
- **pass fail**: PASS
- **strength of evidence**: STRONGLY SUPPORTED
- **remaining limitations**: No independent StVK reference exists (solids4foam only; ANSYS/COMSOL/Tukovic laws unresolved)
- **source file or csv**: docs/reports/CAVITY_BENCHMARK_DEFINITION.md sec. 4; inputs/C2_fsi_pilot_f1.feb
- **reproducible from retained data**: yes (inputs shipped); source read needs FEBio v4.13 source

## A3: Depth normalisation (plane strain, forces per metre)

- **input configuration**: One-element slab TZ=0.1 m, wz and uz fixed on all nodes, all forces divided by TZ
- **source binary provenance**: stock febio4 sha256 f5ae1ed6...72f02 (v4.13, source tree clean; no bit-exact binary proof)
- **mesh timestep**: not applicable
- **cpu allocation**: not applicable
- **observed result**: Forces reported in N per m depth; ANSYS x100 (10 mm), Tukovic x20 (0.05 m) per their authors; COMSOL convention unanswered
- **numerical acceptance criterion**: unit consistency
- **pass fail**: PASS (FEBio side); COMSOL convention unresolved
- **strength of evidence**: PROVEN (FEBio side)
- **remaining limitations**: COMSOL depth convention not confirmed
- **source file or csv**: docs/reports/CAVITY_EXACT_BENCHMARK_CONTRACT.md sec. 3
- **reproducible from retained data**: yes

## A4: Fluid properties and weak compressibility

- **input configuration**: rho=1, mu=0.01 (Re=100), p=-K e_f with K=1e6 Pa (Mach about 1.5e-3)
- **source binary provenance**: stock febio4 sha256 f5ae1ed6...72f02 (v4.13, source tree clean; no bit-exact binary proof)
- **mesh timestep**: not applicable
- **cpu allocation**: not applicable
- **observed result**: K x10 changes u_y by -1.1e-6 and F_y by -1.4e-6 relative at M1 (E1)
- **numerical acceptance criterion**: insensitivity of steady values to K
- **pass fail**: PASS (M1 only)
- **strength of evidence**: PROVEN (M1)
- **remaining limitations**: not tested at M2/M3; incompressible solids4foam vs weakly compressible FEBio
- **source file or csv**: results/numerical_sensitivity_C4.csv
- **reproducible from retained data**: yes (table); rerun needs FEBio

## A5: Inlet and outlet conditions

- **input configuration**: Built-in parabolic 'fluid normal velocity' inlet (mean 1 m/s, discrete flux exactly 1) with rim pressure; outlet zero fluid dilatation (p=0), velocity natural
- **source binary provenance**: stock febio4 sha256 f5ae1ed6...72f02 (v4.13, source tree clean; no bit-exact binary proof)
- **mesh timestep**: not applicable
- **cpu allocation**: not applicable
- **observed result**: Inlet nodal profile is the continuum parabola of flux 1.0101 at h=0.1 (discrete flux exact); outlet p=0 at all 42 nodes
- **numerical acceptance criterion**: flux 1.0000+/-1e-4; outlet datum equal to the reference
- **pass fail**: PASS
- **strength of evidence**: PROVEN (inlet, straight channel C1); STRONGLY SUPPORTED (outlet equivalence)
- **remaining limitations**: natural outlet vs OpenFOAM zeroGradient U: cavity-exit loss equivalence not tested; solids4foam discrete inlet flux not seen
- **source file or csv**: docs/reports/CAVITY_INLET_BOUNDARY_CONDITION_AUDIT.md; docs/reports/CAVITY_POISEUILLE_DIAGNOSTIC.md
- **reproducible from retained data**: yes (inputs, logs gz)

## A6: No-slip conformal FSI interface

- **input configuration**: Shared nodes between fluid interface row and plate top row; w=0 on interface nodes; 'fluid-FSI traction' load on 40/80/160 facets
- **source binary provenance**: stock febio4 sha256 f5ae1ed6...72f02 (v4.13, source tree clean; no bit-exact binary proof)
- **mesh timestep**: M1/M2/M3
- **cpu allocation**: not applicable
- **observed result**: Interface orientation verified through the sign of F_y and downward deflection; traction pathway read from source (F1)
- **numerical acceptance criterion**: orientation and sign consistent
- **pass fail**: PASS (pressure part)
- **strength of evidence**: PROVEN (pressure part, vertical)
- **remaining limitations**: viscous / total traction transfer not tested; FEBio reaction forces are identically zero so no plate force balance
- **source file or csv**: docs/reports/FEBIO_FSI_TRACTION_ASSEMBLY_AUDIT.md
- **reproducible from retained data**: partly (facet tables shipped)

## A7: Pressure constraints (added numerical constraints)

- **input configuration**: 4 floor-end corner nodes pinned to p=1.15 Pa (value from the clean floor mean of C1); inlet rim pressure
- **source binary provenance**: stock febio4 sha256 f5ae1ed6...72f02 (v4.13, source tree clean; no bit-exact binary proof)
- **mesh timestep**: not applicable
- **cpu allocation**: not applicable
- **observed result**: Pins are required (undetermined pressure DOF, Schur eigenvalue 3.5e-10) and contribute -1e-4 Pa to the row parity; force effect bound 0.41/0.19/0.10 % (M1/M2/M3), wrong sign to close the deficit
- **numerical acceptance criterion**: pin effect small against the 5 % deficit
- **pass fail**: PASS (bound)
- **strength of evidence**: PROVEN (operator and linear attribution); bound is an estimate
- **remaining limitations**: a pin is NOT a benchmark condition; alternative pin values / algebraic corner constraints not tested
- **source file or csv**: docs/reports/FEBIO_PRESSURE_MODE_ANALYSIS.md sec. 3; docs/reports/CAVITY_PRESSURE_FORCE_AUDIT.md sec. 2, 7
- **reproducible from retained data**: partly (b4 operator scripts shipped; pin bound from xplt audit JSON)

## A8: Monitoring point A on the wetted top surface

- **input configuration**: Node at (4,0), z layer 0; u_y(A), u_x(A) from pointA.txt
- **source binary provenance**: stock febio4 sha256 f5ae1ed6...72f02 (v4.13, source tree clean; no bit-exact binary proof)
- **mesh timestep**: not applicable
- **cpu allocation**: not applicable
- **observed result**: u_y(A) read at the top surface node; xplt u_y equals pointA.txt to 5e-9
- **numerical acceptance criterion**: top-surface node, not the underside value
- **pass fail**: PASS
- **strength of evidence**: PROVEN
- **remaining limitations**: old underside values (-0.2025/-0.2288) must not be used
- **source file or csv**: results/raw/coupled/*/pointA.txt.gz; scripts/verify_against_raw.py
- **reproducible from retained data**: yes

## B1: Plate-only StVK check, hex20 20x2

- **input configuration**: uniform follower p=1.315 Pa, clamped ends, plane strain, 10 load steps
- **source binary provenance**: stock febio4 sha256 f5ae1ed6...72f02 (v4.13, source tree clean; no bit-exact binary proof)
- **mesh timestep**: 20x2 hex20, 393 nodes, 10 steps
- **cpu allocation**: 2 cores, 1 s (job 13809)
- **observed result**: u_y(A) -0.237313 m vs Q9 reference -0.23887 (-0.65 %)
- **numerical acceptance criterion**: within about 1 % of the Q9 reference of the same mesh label
- **pass fail**: PASS
- **strength of evidence**: STRONGLY SUPPORTED
- **remaining limitations**: reference is another numerical solution; reaction forces reported as zero (equilibrium not verified)
- **source file or csv**: results/plate_only_checks.csv
- **reproducible from retained data**: yes (table); rerun needs FEBio

## B2: Plate-only StVK check, hex20 40x4

- **input configuration**: as B1
- **source binary provenance**: stock febio4 sha256 f5ae1ed6...72f02 (v4.13, source tree clean; no bit-exact binary proof)
- **mesh timestep**: 40x4 hex20, 1343 nodes
- **cpu allocation**: 2 cores, 2 s (job 13810)
- **observed result**: -0.239985 m vs Q9 -0.24047 (-0.20 %); two-level extrapolations agree to 0.05 % (assuming p=2)
- **numerical acceptance criterion**: as B1
- **pass fail**: PASS
- **strength of evidence**: STRONGLY SUPPORTED
- **remaining limitations**: as B1; extrapolation assumes p=2 from two levels
- **source file or csv**: results/plate_only_checks.csv
- **reproducible from retained data**: yes (table)

## B3: Constitutive-law checks: why displacement references spread by 4-5 %

- **input configuration**: manuscript plate-only study (Q9, uniform 1.315 Pa): StVK, neo-Hookean, Hencky
- **source binary provenance**: not applicable (manuscript data, not rerun here)
- **mesh timestep**: 20x2, 40x4 Q9
- **cpu allocation**: not applicable
- **observed result**: StVK -0.23887/-0.24047; neo-Hookean +4.8 %, Hencky +5.7 % relative to StVK
- **numerical acceptance criterion**: law spread explains reference displacement spread
- **pass fail**: n/a (context)
- **strength of evidence**: STRONGLY SUPPORTED
- **remaining limitations**: results are from the manuscript repository (PR #20), not recomputed; not about the FEBio force deficit
- **source file or csv**: docs/reports/CAVITY_EXACT_BENCHMARK_CONTRACT.md sec. 4
- **reproducible from retained data**: no (manuscript data not copied)

## B4: hex8 plate refinement (element family of the coupled plate)

- **input configuration**: 4 independent hex8 StVK plates in one input, uniform 1.315 Pa follower pressure, 20 steps
- **source binary provenance**: stock febio4 sha256 f5ae1ed6...72f02 (v4.13, source tree clean; no bit-exact binary proof)
- **mesh timestep**: 40x4, 40x8, 80x8, 160x16 hex8
- **cpu allocation**: 2 cores, 6 s (job 14013)
- **observed result**: u_y(A) -0.23026/-0.23467/-0.23789/-0.24023 m: 4.5 %/2.6 %/1.3 %/0.3 % smaller than the Q9 extrapolation -0.2410 (observed order about 1.7)
- **numerical acceptance criterion**: convergence to the Q9 value
- **pass fail**: PASS (converges); hex8 stiffer at coupled levels (4.5 %, 1.3 %, 0.3 % at 40x4, 80x8, 160x16)
- **strength of evidence**: PROVEN
- **remaining limitations**: partial cause of the displacement gap at M1/M2 only; Q9 reference extrapolated from two levels
- **source file or csv**: results/plate_only_checks.csv
- **reproducible from retained data**: yes (table); rerun needs FEBio

## B5: S2: pressure-transfer / displacement compliance check

- **input configuration**: hex8 80x8 plate loaded only by the M2 interface pressure (face-mean of the two nodal values per facet), follower load
- **source binary provenance**: stock febio4 sha256 f5ae1ed6...72f02 (v4.13, source tree clean; no bit-exact binary proof)
- **mesh timestep**: 80x8 hex8, 20 steps
- **cpu allocation**: 2 cores, 2 s (job 14025)
- **observed result**: u_y(A) -0.226176 vs coupled -0.226233 (0.025 %); u_x(A) 3.92e-4 vs 3.35e-4
- **numerical acceptance criterion**: within 2 % of the coupled u_y
- **pass fail**: PASS
- **strength of evidence**: PROVEN for the pressure-part self-consistency
- **remaining limitations**: a self-consistency test: it feeds FEBio's own wall-row pressure to the same plate and cannot detect an error in the pressure or distinguish F1 from F2 (the face-mean of nodal p equals the F1 vertical force); viscous traction not tested
- **source file or csv**: results/plate_only_checks.csv; docs/reports/CAVITY_VALIDATION_GATE_AUDIT.md sec. 3
- **reproducible from retained data**: yes (input regenerated from stored M2 interface log)

## B6: Displacement/force compliance ratio

- **input configuration**: u_y(A)/F_y of the coupled results vs solids4foam QM
- **source binary provenance**: stock febio4 sha256 f5ae1ed6...72f02 (v4.13, source tree clean; no bit-exact binary proof)
- **mesh timestep**: M1/M2/M3
- **cpu allocation**: not applicable
- **observed result**: FEBio 0.0449/0.0459/0.0462 m per N/m (extrapolated 0.0463) vs solids4foam 0.0397/0.0442/0.0454 at the same dx (0.0458 extrapolated)
- **numerical acceptance criterion**: same compliance if the plate and load transfer are equivalent
- **pass fail**: consistent at the finest levels
- **strength of evidence**: STRONGLY SUPPORTED
- **remaining limitations**: load distribution (tilt: u_x 3.4-4.8x larger than solids4foam) and StVK nonlinearity also enter; reference compared by ratio only
- **source file or csv**: results/comparison_febio_vs_references.csv
- **reproducible from retained data**: yes

## B7: Plate reaction-force equilibrium

- **input configuration**: reaction forces of the clamped nodes (solid module and xplt)
- **source binary provenance**: stock febio4 sha256 f5ae1ed6...72f02 (v4.13, source tree clean; no bit-exact binary proof)
- **mesh timestep**: not applicable
- **cpu allocation**: not applicable
- **observed result**: reaction forces are identically zero in the node log and xplt for plate-only and coupled runs
- **numerical acceptance criterion**: sum of reactions = applied load
- **pass fail**: NOT EVALUABLE
- **strength of evidence**: NOT TESTED
- **remaining limitations**: cause not established (solver code assigns node loads to equation index <= -2 only; not tested without a run); says nothing for or against any force value
- **source file or csv**: docs/reports/CAVITY_PLATE_ONLY_RESULTS.md; docs/reports/CAVITY_PRESSURE_FORCE_AUDIT.md sec. 1
- **reproducible from retained data**: yes (documented negative result)

## B8: Remaining structural differences

- **input configuration**: hex8 3-D slab with u_z=0 vs Q9/FV plane strain
- **source binary provenance**: stock febio4 sha256 f5ae1ed6...72f02 (v4.13, source tree clean; no bit-exact binary proof)
- **mesh timestep**: not applicable
- **cpu allocation**: not applicable
- **observed result**: hex8 stiffness (B4) and load asymmetry (u_x) are the known differences; deformation difference changes F_y by <= 0.04 % (bounded from C1b vs M1)
- **numerical acceptance criterion**: small against the 5 % force deficit
- **pass fail**: PASS (cannot explain the deficit)
- **strength of evidence**: STRONGLY SUPPORTED
- **remaining limitations**: bound from two fixed/deformed M1 states only
- **source file or csv**: docs/reports/CAVITY_SOLVER_FIELD_COMPARISON.md sec. 8
- **reproducible from retained data**: partly

## C0: Early fluid-only attempts with nodal Dirichlet or 'fluid velocity' math inlets

- **input configuration**: 5 runs (K 1e4/1e8, mu 0.01/0.1, straight and cavity geometry)
- **source binary provenance**: stock febio4 sha256 f5ae1ed6...72f02 (v4.13, source tree clean; no bit-exact binary proof)
- **mesh timestep**: dt 0.2-0.5 s, 1688 or 1400 hex8
- **cpu allocation**: 2 cores, 20-96 s each (jobs 13856, 13864, 13875, 13881, 13888); total 0.14 CPU-h
- **observed result**: all failed (negative Jacobian or max reformations) within 1.8-17 s of model time; inlet dilatation spikes next to wall corners; 'math' profile evaluated at the origin in v4.13 FEFluidVelocity::Update
- **numerical acceptance criterion**: steady solution
- **pass fail**: FAIL (preserved)
- **strength of evidence**: STRONGLY SUPPORTED (source audit) / HYPOTHESIS (nodal-inlet omits flux term; no K-matched run)
- **remaining limitations**: two of the five shipped inputs differ from the current generator output only in diagnostic monitor node sets
- **source file or csv**: docs/reports/CAVITY_FLUID_ONLY_RESULTS.md; docs/reports/CAVITY_INLET_BOUNDARY_CONDITION_AUDIT.md; inputs/early_failed_attempts/
- **reproducible from retained data**: yes (inputs, logs gz); rerun needs FEBio

## C1: Poiseuille reference: straight 14x1 channel with the built-in parabolic inlet

- **input configuration**: K=1e6, mu=0.01, ramp 10 s, t=150 s, outlet zero dilatation
- **source binary provenance**: stock febio4 sha256 f5ae1ed6...72f02 (v4.13, source tree clean; no bit-exact binary proof)
- **mesh timestep**: 140x10 hex8 (1400), dt 0.5 s, 300 steps
- **cpu allocation**: 4 cores, 57 s (job 13977), 0.063 CPU-h
- **observed result**: inlet flux 1.0000000000, outlet +1.7e-6, dp/dx -0.1212124 vs -0.1212122 (actual-flux value; +1.01 % vs nominal -0.12 explained by the nodal inlet parabola), profile error 1.1e-6; frozen x-sawtooth 0.0124 Pa and y-alternation 0.0013 Pa
- **numerical acceptance criterion**: mass imbalance 1e-3; gradient vs actual-flux value
- **pass fail**: PASS (inlet treatment); sawtooth caveat
- **strength of evidence**: PROVEN (straight channel)
- **remaining limitations**: does not test cavity geometry, convection, force or coupling; steady state evidenced by solution, Newton tables degenerate after about step 46
- **source file or csv**: docs/reports/CAVITY_POISEUILLE_DIAGNOSTIC.md; results/raw/fluid_only/B0_poiseuille_normalvel/poiseuille_check.json
- **reproducible from retained data**: partly (summary JSON; Mon_all needs HPC output)

## C2: C1: fixed-plate fluid-only cavity without corner pin

- **input configuration**: Re=100, K=1e6, dt 0.5, ramp 10, t_end 300
- **source binary provenance**: stock febio4 sha256 f5ae1ed6...72f02 (v4.13, source tree clean; no bit-exact binary proof)
- **mesh timestep**: 1688 hex8, dt 0.5, 600 steps
- **cpu allocation**: 4 cores, 336 s (job 14012), 0.373 CPU-h
- **observed result**: floor pressure drifts linearly (2.1 -> 44.9 Pa mean over 300 s) at 43 nodes of the two floor-end corners; F_y raw -70.3 N/m
- **numerical acceptance criterion**: A3 steady (<1e-4), A4 upstream dp/dx within 2 %
- **pass fail**: FAIL (A3, A4)
- **strength of evidence**: PROVEN (drift) ; corner mechanism PROVEN by the later operator analysis
- **remaining limitations**: failed run preserved
- **source file or csv**: docs/reports/CAVITY_FLUID_ONLY_RESULTS.md; docs/reports/CAVITY_VALIDATION_GATE_AUDIT.md
- **reproducible from retained data**: yes (logs, summary JSON); Mon_all HPC

## C3: C1b: corrected corner constraints (4 corner nodes pinned to 1.15 Pa)

- **input configuration**: as C1 plus 'prescribed fluid dilatation' p=1.15 Pa at the 4 corner nodes
- **source binary provenance**: stock febio4 sha256 f5ae1ed6...72f02 (v4.13, source tree clean; no bit-exact binary proof)
- **mesh timestep**: 1688 hex8, dt 0.5, 600 steps
- **cpu allocation**: 4 cores, 173 s (job 14015), 0.19 CPU-h
- **observed result**: steady (F_y spread 1e-12), floor p 1.029-1.348 Pa, mass imbalance +9.1e-5, F_y raw -4.6637 N/m, odd-even bracket +/-2.8 %; A4 upstream dp/dx raw +60 % (FAIL as written; pair-averaged -0.3 %, criterion defect diagnosed)
- **numerical acceptance criterion**: A3, A5-A8 pass; A4 within 2 %
- **pass fail**: PASS A3,A5-A8; FAIL A4 as written
- **strength of evidence**: PROVEN (steady, flux); the A4 explanation is STRONGLY SUPPORTED (derived from the same data)
- **remaining limitations**: strict fluid-only gate not passed; the coupled runs proceeded as a diagnostic by documented judgement
- **source file or csv**: docs/reports/CAVITY_FLUID_ONLY_RESULTS.md; results/raw/fluid_only/C1b_fluid_fixedplate_f1_pin/c1_check.json
- **reproducible from retained data**: partly (summary JSON, logs, Mon_iface gz; Mon_all HPC 60 MB)

## C4: Flux conservation

- **input configuration**: inlet vs outlet and x=1.03, 12.99 columns
- **source binary provenance**: stock febio4 sha256 f5ae1ed6...72f02 (v4.13, source tree clean; no bit-exact binary proof)
- **mesh timestep**: C1b / coupled M1-M3
- **cpu allocation**: not applicable
- **observed result**: C1b +9.1e-5; coupled +1.2e-4 / -1.9e-5 / -2.3e-5 (Q_out/Q_in - 1)
- **numerical acceptance criterion**: 1e-3
- **pass fail**: PASS
- **strength of evidence**: PROVEN
- **remaining limitations**: nodal column flux; continuity rows constrain only pair sums
- **source file or csv**: results/coupled_febio_summary.csv; results/ppp_round1_C1b.csv
- **reproducible from retained data**: yes (Mon_inlet/outlet gz)

## C5: Pressure modes in the fluid-only fields

- **input configuration**: C1b (M1 mesh) and straight-channel field analysis
- **source binary provenance**: stock febio4 sha256 f5ae1ed6...72f02 (v4.13, source tree clean; no bit-exact binary proof)
- **mesh timestep**: not applicable
- **cpu allocation**: not applicable
- **observed result**: cavity: row parity D (row1-row0) 0.415-0.421 Pa (window dependent), present everywhere in the domain; straight channel: x-sawtooth 0.0124 Pa, y-alternation 0.0013 Pa
- **numerical acceptance criterion**: descriptive
- **pass fail**: n/a
- **strength of evidence**: PROVEN (existence, in the solver unknown e_f); mechanism STRONGLY SUPPORTED (see F7, F8)
- **remaining limitations**: amplification chain from forcing to soft modes is HYPOTHESIS
- **source file or csv**: results/row_parity_by_level.csv; docs/reports/FEBIO_PRESSURE_MODE_ANALYSIS.md
- **reproducible from retained data**: partly

## C6: Momentum and continuity residuals of the converged C1b field (independent python re-assembly of the source element equations)

- **input configuration**: nodal w, e_f of C1b at t=300 and the validated operator
- **source binary provenance**: stock febio4 sha256 f5ae1ed6...72f02 (v4.13, source tree clean; no bit-exact binary proof)
- **mesh timestep**: M1 fixed plate
- **cpu allocation**: login node
- **observed result**: free-row momentum residual 9.4e-7 relative, continuity 2.2e-7
- **numerical acceptance criterion**: <= 1e-5
- **pass fail**: PASS
- **strength of evidence**: PROVEN (for the fluid-only domain of C1b)
- **remaining limitations**: FSI domain (mesh/solid terms) not re-assembled
- **source file or csv**: results/raw/pressure_mode/b4_residual_check.json
- **reproducible from retained data**: no (needs C1b Mon_all.txt, 60 MB, HPC)

## C7: Outlet and corner conditions

- **input configuration**: outlet zero dilatation; corner pins; inlet rim
- **source binary provenance**: stock febio4 sha256 f5ae1ed6...72f02 (v4.13, source tree clean; no bit-exact binary proof)
- **mesh timestep**: not applicable
- **cpu allocation**: not applicable
- **observed result**: outlet ef=0 at all outlet nodes (max abs 0); constraint sets are disjoint; with pins no Schur eigenvalue below 2.5e-3
- **numerical acceptance criterion**: no undetermined pressure DOF
- **pass fail**: PASS
- **strength of evidence**: PROVEN (spectrum)
- **remaining limitations**: see A7
- **source file or csv**: docs/reports/CAVITY_PRESSURE_FORCE_AUDIT.md sec. 7; results/raw/pressure_mode/b4_operator_spectra.json
- **reproducible from retained data**: yes (b4_operator.py is package-only)

## D1: M1 coupled FSI (C2 pilot)

- **input configuration**: generalised-alpha rhoi=0, dt 2.5 s, inlet ramp 25 s, K=1e6, E_mesh=1e-3 Pa, pin 1.15 Pa, hex8 StVK plate, full Newton; t_end 1250 s requested
- **source binary provenance**: stock febio4 sha256 f5ae1ed6...72f02 (v4.13, source tree clean; no bit-exact binary proof)
- **mesh timestep**: 1688 fluid + 160 plate hex8, dx 0.1, dt 2.5 s
- **cpu allocation**: 4 cores, 914 s FEBio wall (15:15 Slurm), 1.02 CPU-h (job 14016, 3660 CPU-s), shared loaded node
- **observed result**: latest accepted t=972.5 s (step 389): u_y(A) -0.211967 m, u_x(A) 3.57e-4, F_y raw (F1) -4.7176 N/m, mass imbalance +1.2e-4; 389 steps, 2782 equilibrium iterations (7.15/step); ERROR termination at step 390 (t=975): 'Max nr of reformations reached', dilatation increment on the arithmetic floor (1.8e-33 vs 1.7e-34 required), other norms 1e-24..1e-27
- **numerical acceptance criterion**: normal termination; steady spread <1e-4 over the last 10 %; mass 1e-3
- **pass fail**: FAIL (abnormal termination); spread 7.6e-5 passes on [875,972.5]
- **strength of evidence**: PROVEN (values as stored); conditional steady estimate
- **remaining limitations**: failed exit; planned window [1125,1250] never reached; the force is 7.9 % below solids4foam QM at the same dx
- **source file or csv**: results/coupled_febio_summary.csv; results/solver_termination_and_cost.csv; results/raw/coupled/C2_fsi_pilot_f1/
- **reproducible from retained data**: yes (table, logs gz, pointA gz); rerun needs FEBio

## D2: M2 coupled FSI

- **input configuration**: generalised-alpha rhoi=0, dt 2.5 s, inlet ramp 25 s, K=1e6, E_mesh=1e-3 Pa, pin 1.15 Pa, hex8 StVK plate, full Newton; t_end 900 s
- **source binary provenance**: stock febio4 sha256 f5ae1ed6...72f02 (v4.13, source tree clean; no bit-exact binary proof)
- **mesh timestep**: 6752 + 640 hex8, dx 0.05, dt 2.5 s
- **cpu allocation**: 4 cores, 1:06:17 (job 14019), 4.42 CPU-h
- **observed result**: latest accepted t=870 s (step 348): u_y(A) -0.226233, u_x 3.35e-4, F_y raw -4.9339, mass -1.9e-5; 348 steps, 2461 iterations (7.07/step); ERROR termination at step 349 (t=872.5), same mechanism (dilatation 1.9e-32 vs 5.2e-33 required)
- **numerical acceptance criterion**: as D1
- **pass fail**: FAIL (abnormal termination); steady spread 1.79e-4 NOT met
- **strength of evidence**: PROVEN (values as stored); conditional estimate
- **remaining limitations**: as D1; F_y 5.6 % below QM
- **source file or csv**: results/coupled_febio_summary.csv
- **reproducible from retained data**: yes (table, logs); rerun needs FEBio

## D3: M3 coupled FSI

- **input configuration**: generalised-alpha rhoi=0, dt 2.5 s, inlet ramp 25 s, K=1e6, E_mesh=1e-3 Pa, pin 1.15 Pa, hex8 StVK plate, full Newton; t_end 900 s
- **source binary provenance**: stock febio4 sha256 f5ae1ed6...72f02 (v4.13, source tree clean; no bit-exact binary proof)
- **mesh timestep**: 27008 + 2560 hex8, dx 0.025, dt 2.5 s
- **cpu allocation**: 8 cores, 3:31:29 (job 14021), 28.2 CPU-h; 423 % CPU use on 8 allocated cores (shared loaded node)
- **observed result**: latest accepted t=715 s (step 286): u_y(A) -0.230038, u_x 3.26e-4, F_y raw -4.9835, mass -2.3e-5; 286 steps, 2031 iterations (7.10/step); ERROR termination at step 287 (t=717.5): dilatation 6.2e-30 vs 3.9e-31 required; tail is a decaying ~44 s oscillation (e-fold about 86 s, +/-1.2e-4 m)
- **numerical acceptance criterion**: as D1
- **pass fail**: FAIL (abnormal termination); steady spread 1.05e-3 NOT met (criterion 1e-4)
- **strength of evidence**: PROVEN (values as stored); NOT steady to 1e-4
- **remaining limitations**: values carry about +/-2e-4 m (u_y) and +/-1e-3 relative; M4 not run; xplt state at 717.5 is a failure-exit state (flag 2)
- **source file or csv**: results/coupled_febio_summary.csv; results/raw/b2/history_C3_fsi_M4.json
- **reproducible from retained data**: yes (table, logs); rerun needs FEBio (28 CPU-h)

## D4: Forces and displacements vs same-dx solids4foam

- **input configuration**: FEBio raw F1 and u_y vs solids4foam QM #618 (F_y) and partitioned (u_y)
- **source binary provenance**: stock febio4 sha256 f5ae1ed6...72f02 (v4.13, source tree clean; no bit-exact binary proof)
- **mesh timestep**: M1/M2/M3
- **cpu allocation**: not applicable
- **observed result**: F_y -7.9/-5.6/-5.1 % vs QM; u_y +3.6/-2.2/-3.5 % vs partitioned; extrapolated u_y -0.2314 (-3.9 % vs -0.2407) and F_y about -5.00 (-4.9 % vs -5.26): INDICATIVE
- **numerical acceptance criterion**: agreement with the benchmark values (reference, not target)
- **pass fail**: FAIL (not in agreement)
- **strength of evidence**: PROVEN (numbers); cause of the gap NOT established
- **remaining limitations**: three levels, all error-terminated, M3 not steady, F_y order not reliable (M2->M3 change equals the odd-even bracket)
- **source file or csv**: results/comparison_febio_vs_references.csv; results/mesh_convergence_table.csv
- **reproducible from retained data**: yes

## D5: Solver convergence history and termination reason

- **input configuration**: FEBio logs
- **source binary provenance**: stock febio4 sha256 f5ae1ed6...72f02 (v4.13, source tree clean; no bit-exact binary proof)
- **mesh timestep**: M1/M2/M3
- **cpu allocation**: not applicable
- **observed result**: 0 cut-backs, no negative Jacobian, 7.1 iterations/step; all three end with max reformations (25) on the dilatation-increment norm at steady state; every other norm at its 1e-23..1e-27 floor
- **numerical acceptance criterion**: normal termination
- **pass fail**: FAIL
- **strength of evidence**: PROVEN
- **remaining limitations**: the relative dilatation criterion (1e-12 of the step's own increment ~1e-22) sits at the arithmetic floor: a stopping-rule artefact at steady state, not a divergence; not tested with a changed dilatation tolerance
- **source file or csv**: results/solver_termination_and_cost.csv
- **reproducible from retained data**: yes (logs gz)

## D6: Steady-state acceptance

- **input configuration**: relative spread of u_y(A) over the last 10 % of the available history
- **source binary provenance**: stock febio4 sha256 f5ae1ed6...72f02 (v4.13, source tree clean; no bit-exact binary proof)
- **mesh timestep**: M1/M2/M3
- **cpu allocation**: not applicable
- **observed result**: 7.6e-5 (met, window [875,972.5]) / 1.8e-4 (not met) / 1.05e-3 (not met)
- **numerical acceptance criterion**: < 1e-4
- **pass fail**: M1 met; M2, M3 not met
- **strength of evidence**: PROVEN
- **remaining limitations**: M3 residual is a decaying plate-fluid mode (period 44 s) rather than drift
- **source file or csv**: results/coupled_febio_summary.csv
- **reproducible from retained data**: yes

## E1: Compressibility K x10

- **input configuration**: K=1e7 (pin dilatation scaled), t_end 900
- **source binary provenance**: stock febio4 sha256 f5ae1ed6...72f02 (v4.13, source tree clean; no bit-exact binary proof)
- **mesh timestep**: M1, dt 2.5
- **cpu allocation**: 4 cores, 13:04 (job 14023), 0.87 CPU-h
- **observed result**: u_y -1.1e-6, F_y -1.4e-6 relative
- **numerical acceptance criterion**: changes below the steady-window spread (1e-4)
- **pass fail**: PASS
- **strength of evidence**: PROVEN (M1)
- **remaining limitations**: M1 only
- **source file or csv**: results/numerical_sensitivity_C4.csv
- **reproducible from retained data**: yes (table)

## E2: Time step

- **input configuration**: dt 1.25 s
- **source binary provenance**: stock febio4 sha256 f5ae1ed6...72f02 (v4.13, source tree clean; no bit-exact binary proof)
- **mesh timestep**: M1, 720 steps
- **cpu allocation**: 4 cores, 28:50 (job 14022), 1.92 CPU-h
- **observed result**: u_y +4.1e-5, F_y +2.2e-5 relative (|value| larger)
- **numerical acceptance criterion**: as E1
- **pass fail**: PASS (independent within about 1e-4)
- **strength of evidence**: PROVEN (M1)
- **remaining limitations**: two dt values give no order; same size as the window spread (3.8e-4)
- **source file or csv**: results/numerical_sensitivity_C4.csv
- **reproducible from retained data**: yes

## E3: Stopping tolerances

- **input configuration**: dtol 1e-5, etol 1e-6, rtol 1e-8 (vtol, ftol unchanged: on the arithmetic floor)
- **source binary provenance**: stock febio4 sha256 f5ae1ed6...72f02 (v4.13, source tree clean; no bit-exact binary proof)
- **mesh timestep**: M1
- **cpu allocation**: 4 cores, 13:17 (job 14024), 0.89 CPU-h
- **observed result**: u_y identical to 8 digits, F_y 5e-9 relative
- **numerical acceptance criterion**: as E1
- **pass fail**: PASS
- **strength of evidence**: PROVEN (M1)
- **remaining limitations**: the tightened criteria were not binding; vtol/ftol/dilatation tolerance not varied
- **source file or csv**: results/numerical_sensitivity_C4.csv
- **reproducible from retained data**: yes

## E4: Mesh refinement of the coupled model

- **input configuration**: M1-M3
- **source binary provenance**: stock febio4 sha256 f5ae1ed6...72f02 (v4.13, source tree clean; no bit-exact binary proof)
- **mesh timestep**: dx 0.1/0.05/0.025
- **cpu allocation**: 1.02 + 4.42 + 28.2 CPU-h
- **observed result**: u_y order 1.91, extrapolation -0.2314; F_y nominal order 2.12 (unreliable), extrapolation about -5.00; deficit vs QM -7.9/-5.6/-5.1 % shrinks then flattens
- **numerical acceptance criterion**: convergence toward the reference
- **pass fail**: FAIL (converges to a value about 4-5 % below the reference)
- **strength of evidence**: PROVEN (numbers); INDICATIVE (extrapolation)
- **remaining limitations**: see D4
- **source file or csv**: results/mesh_convergence_table.csv
- **reproducible from retained data**: yes

## E5: Corner-pin effect on the force

- **input configuration**: replace pinned values by linear extrapolation of the two nearest nodes (bound)
- **source binary provenance**: stock febio4 sha256 f5ae1ed6...72f02 (v4.13, source tree clean; no bit-exact binary proof)
- **mesh timestep**: M1/M2/M3
- **cpu allocation**: login node
- **observed result**: -0.0194 / -0.0096 / -0.0050 N/m (0.41 / 0.19 / 0.10 %), larger |F_y|
- **numerical acceptance criterion**: small vs 5 %
- **pass fail**: PASS (cannot explain the deficit)
- **strength of evidence**: PROVEN (linear attribution); bound is an estimate
- **remaining limitations**: no run with a different pin value
- **source file or csv**: results/raw/xplt_audit/cavity_xplt_audit_{C2,M2,M3}.json
- **reproducible from retained data**: yes (JSON)

## E6: Pressure-force reconstruction F1..F4 at stored states (never mixing times)

- **input configuration**: facet tables extracted from the xplt: M1 t=900, M2 t=750, M3 t=500 (and M3 t=717.5 failure-exit state)
- **source binary provenance**: stock febio4 sha256 f5ae1ed6...72f02 (v4.13, source tree clean; no bit-exact binary proof)
- **mesh timestep**: M1/M2/M3
- **cpu allocation**: 0.6 s on the login node
- **observed result**: F1 -4.7176/-4.9339/-4.9809; F2 -5.4375/-5.2748/-5.2469; M3 t=717.5: F1 -4.9835, F2 -5.2494; F1-F2 = 0.720/0.341/0.266 N/m equals the row-parity part to <=0.001; F2 = F4 (<=5e-8)
- **numerical acceptance criterion**: reconstruction reproduces the tool F_y to 1e-7
- **pass fail**: PASS (reconstruction); F2 is NOT established as the physical force
- **strength of evidence**: PROVEN (arithmetic decomposition); which value is physical is NOT ESTABLISHED
- **remaining limitations**: a re-evaluation of stored fields, not independent solver verification; viscous part of F1 not reconstructable; pair-mean agrees with solids4foam at M3 only (+6.1 % M1, +0.9 % M2)
- **source file or csv**: results/force_reconstruction_F1_F4.csv; results/raw/facets/
- **reproducible from retained data**: yes from facet tables; the facet tables themselves need the retained xplt (HPC)

## E7: Viscous contribution (2-point wall-gradient estimate)

- **input configuration**: element-average velocities of the first fluid rows, one-sided estimate
- **source binary provenance**: stock febio4 sha256 f5ae1ed6...72f02 (v4.13, source tree clean; no bit-exact binary proof)
- **mesh timestep**: M1/M2/M3
- **cpu allocation**: login node
- **observed result**: F_y,visc -0.0007/-0.0012/-0.0015 N/m (<=0.03 % of F_y); F_x,visc -0.0197/-0.0244/-0.0259
- **numerical acceptance criterion**: small vs 5 %
- **pass fail**: PASS (bound)
- **strength of evidence**: STRONGLY SUPPORTED (order of magnitude only)
- **remaining limitations**: NOT independently reconstructed: the solver-applied viscous traction (Gauss-point mean over the adjacent element) is not stored; method spread +/-0.003 in F_x
- **source file or csv**: results/raw/xplt_audit/
- **reproducible from retained data**: yes (JSON); recomputation needs xplt

## E8: Parameters never varied

- **input configuration**: rhoi, ramp length, E_mesh, pin value, dilatation tolerance, M2/M3 sensitivities, M4 mesh (dx 0.0125, about 120 CPU-h)
- **source binary provenance**: not applicable
- **mesh timestep**: not applicable
- **cpu allocation**: not applicable
- **observed result**: not run
- **numerical acceptance criterion**: not applicable
- **pass fail**: NOT TESTED
- **strength of evidence**: NOT TESTED
- **remaining limitations**: start-up route independence is shown only for dt/2 at M1
- **source file or csv**: docs/reports/CAVITY_NEXT_DECISION.md
- **reproducible from retained data**: no

## F1: Source version

- **input configuration**: FEBio source tag v4.13
- **source binary provenance**: git describe v4.13, HEAD 32ae206ff4881dfb54f62296cd1558e58ed9fcc6 (shallow grafted checkout)
- **mesh timestep**: not applicable
- **cpu allocation**: not applicable
- **observed result**: run logs print 'version 4.13.0.32ae206ff'
- **numerical acceptance criterion**: tag equals the version string
- **pass fail**: PASS
- **strength of evidence**: PROVEN (tag); STRONGLY SUPPORTED (binary built from it)
- **remaining limitations**: see F3
- **source file or csv**: provenance/BUILD_PROVENANCE.md
- **reproducible from retained data**: yes (log gz)

## F2: Source-tree cleanliness

- **input configuration**: git status / git diff of the build tree, checked before and after all experimental work
- **source binary provenance**: stock tree
- **mesh timestep**: not applicable
- **cpu allocation**: not applicable
- **observed result**: 0 changed files, empty diff; the unrelated bug-report patch for initial fluid velocity is not applied
- **numerical acceptance criterion**: clean
- **pass fail**: PASS
- **strength of evidence**: PROVEN (at check time)
- **remaining limitations**: a statement about the working tree, not about the compiled binary
- **source file or csv**: provenance/BUILD_PROVENANCE.md
- **reproducible from retained data**: no (needs the maintainers' source tree)

## F3: Binary provenance

- **input configuration**: CMake cache, timestamps, hashes
- **source binary provenance**: stock febio4 sha256 f5ae1ed6...72f02 (v4.13, source tree clean; no bit-exact binary proof)
- **mesh timestep**: not applicable
- **cpu allocation**: not applicable
- **observed result**: Release -O3, MKL, OpenMP; source files mtime 18:21:07, libfebiofluid 18:23:58, febio4 18:24:02 (2026-10-29): binary linked after the last source write
- **numerical acceptance criterion**: binary built from the clean tag
- **pass fail**: NOT PROVEN bit-exactly
- **strength of evidence**: STRONGLY SUPPORTED
- **remaining limitations**: no clean rebuild and binary comparison was done: LIMITATION, not disguised
- **source file or csv**: provenance/BUILD_PROVENANCE.md
- **reproducible from retained data**: no

## F4: Pressure interpolation

- **input configuration**: FEFluidFSIDomain3D (registered fluid-FSI-3D)
- **source binary provenance**: stock febio4 sha256 f5ae1ed6...72f02 (v4.13, source tree clean; no bit-exact binary proof); source read at v4.13
- **mesh timestep**: not applicable
- **cpu allocation**: not applicable
- **observed result**: p=-K e_f; w and e_f share the same trilinear nodal space (equal order, 8-point Gauss); element 'fluid pressure' output = mean of the 8 nodal values
- **numerical acceptance criterion**: source inspection + data (xplt element pressure = -K mean nodal e_f to 1e-7 Pa)
- **pass fail**: n/a
- **strength of evidence**: PROVEN
- **remaining limitations**: source lines quoted from the maintainers' tree (FEBio source not copied)
- **source file or csv**: docs/reports/FEBIO_V413_PRESSURE_SOURCE_AUDIT.md
- **reproducible from retained data**: partly (needs FEBio source)

## F5: Dilatation degrees of freedom

- **input configuration**: Activate() of the FSI domain
- **source binary provenance**: stock febio4 sha256 f5ae1ed6...72f02 (v4.13, source tree clean; no bit-exact binary proof)
- **mesh timestep**: not applicable
- **cpu allocation**: not applicable
- **observed result**: every non-excluded fluid node has w and e_f active; only one fluid-FSI domain class exists, so a compatible (Taylor-Hood-type) space cannot be requested from the input
- **numerical acceptance criterion**: source inspection
- **pass fail**: n/a
- **strength of evidence**: PROVEN
- **remaining limitations**: hex20/hex27 fluid elements share the class (equal order again): NOT TESTED
- **source file or csv**: docs/reports/FEBIO_V413_PRESSURE_SOURCE_AUDIT.md sec. 2
- **reproducible from retained data**: partly

## F6: Momentum operator

- **input configuration**: ElementInternalForce / ElementStiffness blocks, rest-state linearisation
- **source binary provenance**: stock febio4 sha256 f5ae1ed6...72f02 (v4.13, source tree clean; no bit-exact binary proof)
- **mesh timestep**: straight and cavity M1 meshes
- **cpu allocation**: login node
- **observed result**: pressure enters the momentum row as the strong-form gradient (no boundary integral); continuity block kJw = G^T (max difference 1.7e-17); system [A G; G^T -(c/K) M] with c/K = 5.3e-7
- **numerical acceptance criterion**: numerical identity
- **pass fail**: PASS
- **strength of evidence**: PROVEN
- **remaining limitations**: linearised about rest on a fixed mesh, Stokes + mass term; convection not in the operator
- **source file or csv**: results/raw/pressure_mode/b4_operator_spectra.json; scripts/b4_operator.py
- **reproducible from retained data**: yes (package-only script; run recorded in docs/REPRODUCIBILITY_AUDIT.md)

## F7: Row-parity pressure mode (Schur analysis)

- **input configuration**: S=G^T A^-1 G eigenproblem; backward-Euler Stokes emulation; exact attribution on C1b
- **source binary provenance**: stock febio4 sha256 f5ae1ed6...72f02 (v4.13, source tree clean; no bit-exact binary proof)
- **mesh timestep**: straight channel h=0.2/0.1; cavity M1
- **cpu allocation**: login node
- **observed result**: y-parity eigenvalue 1.0e-3 at h=0.1 (about 2000 x the compressibility term, 100 x softer than the softest smooth mode); not a null mode with memory (relaxes in one step in the emulation); attribution of D=0.417 Pa: viscous 0.177 + convective 0.240 + fixed values -0.0001
- **numerical acceptance criterion**: classification of the mode
- **pass fail**: n/a
- **strength of evidence**: STRONGLY SUPPORTED (soft inf-sup-deficient mode of the unstabilised equal-order pair); mechanism chain HYPOTHESIS
- **remaining limitations**: emulation is linear Stokes; the amplification of the real forcing (58 % convective) is not isolated; ratios from two mesh levels of one geometry
- **source file or csv**: docs/reports/FEBIO_PRESSURE_MODE_ANALYSIS.md; results/raw/pressure_mode/
- **reproducible from retained data**: partly (spectra scripts package-only; attribution needs C1b Mon_all)

## F8: Mesh-dependent eigenvalue scaling

- **input configuration**: straight channel h = 0.2, 0.1 (0.05 noisy)
- **source binary provenance**: stock febio4 sha256 f5ae1ed6...72f02 (v4.13, source tree clean; no bit-exact binary proof)
- **mesh timestep**: not applicable
- **cpu allocation**: login node
- **observed result**: parity families scale as h^2 (ratios 3.97 / 3.51 / 3.63), smooth mode h-independent (0.97)
- **numerical acceptance criterion**: scaling
- **pass fail**: n/a
- **strength of evidence**: PROVEN for the two levels studied
- **remaining limitations**: two levels, one geometry; the observed row-parity amplitude D falls much more slowly than h^2 (0.365/0.173/0.134 Pa)
- **source file or csv**: results/raw/pressure_mode/b4_modal_scaling.json
- **reproducible from retained data**: yes (package-only script)

## F9: Assembled residual reconstruction

- **input configuration**: python re-assembly from the source blocks on the FEBio C1b field
- **source binary provenance**: stock febio4 sha256 f5ae1ed6...72f02 (v4.13, source tree clean; no bit-exact binary proof)
- **mesh timestep**: M1 fixed plate
- **cpu allocation**: login node
- **observed result**: momentum 9.4e-7, continuity 2.2e-7 relative
- **numerical acceptance criterion**: <=1e-5
- **pass fail**: PASS
- **strength of evidence**: PROVEN (fluid-only domain)
- **remaining limitations**: see C6
- **source file or csv**: results/raw/pressure_mode/b4_residual_check.json
- **reproducible from retained data**: no (needs C1b Mon_all)

## F10: Interface traction assembly (two pathways)

- **input configuration**: FEFluidFSITraction::LoadVector vs FEPlotFluidSurfaceForce
- **source binary provenance**: stock febio4 sha256 f5ae1ed6...72f02 (v4.13, source tree clean; no bit-exact binary proof)
- **mesh timestep**: M1/M2/M3
- **cpu allocation**: login node
- **observed result**: F1 (applied load): nodal e_f of the 4 facet nodes, 2x2 Gauss, plus viscous stress from the adjacent element Gauss mean; F2 (plot definition): element Gauss-mean pressure; F2 equals the pair-row mean F4 for planar facets (algebraic identity)
- **numerical acceptance criterion**: reconstruction of the stored F_y
- **pass fail**: PASS
- **strength of evidence**: PROVEN (what each routine does); STRONGLY SUPPORTED (F1 contains the parity offset, F2 does not)
- **remaining limitations**: F2 is not established as the physical interface force; the built-in surface-force plot was not requested in the cavity inputs
- **source file or csv**: docs/reports/FEBIO_FSI_TRACTION_ASSEMBLY_AUDIT.md; results/force_reconstruction_F1_F4.csv
- **reproducible from retained data**: yes (facet tables)

## F11: Search for interior pressure stabilisation in v4.13

- **input configuration**: grep stabiliz|pspg|brezzi|supg|upwind over FEBioFluid
- **source binary provenance**: not applicable (source search)
- **mesh timestep**: not applicable
- **cpu allocation**: not applicable
- **observed result**: only boundary outflow stabilisations (backflow, tangential flow) and solute/PF variants; none in FEFluidFSIDomain3D / FEFluidDomain3D
- **numerical acceptance criterion**: search
- **pass fail**: n/a
- **strength of evidence**: PROVEN (v4.13 tree searched)
- **remaining limitations**: other FEBio versions not searched
- **source file or csv**: docs/reports/FEBIO_V413_PRESSURE_SOURCE_AUDIT.md sec. 2
- **reproducible from retained data**: no (needs source)

## F12: Row parity independent of mesh motion and of coupling

- **input configuration**: C1b fixed mesh vs C2 moving mesh; dt/2, K x10, tolerance (E1-E3)
- **source binary provenance**: stock febio4 sha256 f5ae1ed6...72f02 (v4.13, source tree clean; no bit-exact binary proof)
- **mesh timestep**: M1
- **cpu allocation**: not applicable
- **observed result**: D = 0.415 (C1b fixed) vs 0.363 (C2); identical to 4-5 digits under dt/2, K x10, tolerance; absent in the straight channel
- **numerical acceptance criterion**: descriptive
- **pass fail**: n/a
- **strength of evidence**: STRONGLY SUPPORTED
- **remaining limitations**: deterministic discretisation property; solids4foam pressure field never seen
- **source file or csv**: results/row_parity_by_level.csv
- **reproducible from retained data**: yes (ybranch JSON)

## G1: Element-level tests of the polynomial-pressure-projection (Dohrmann-Bochev) term

- **input configuration**: 12 random cases (rotation, 15 % distortion, aspect, 2^3/3^3 Gauss, alphaf 1/0.5) against the compiled library
- **source binary provenance**: experimental build febio4 85c344a7...; libfebiofluid 1bff35b8... (round 1)
- **mesh timestep**: not applicable
- **cpu allocation**: login node (seconds)
- **observed result**: alpha=0 bitwise no-op 12/12; dS/de vs FD 2.5e-13; dS/dx vs FD 3.5e-10; parity energy exact to 2e-16
- **numerical acceptance criterion**: FD agreement; no-op at alpha=0
- **pass fail**: PASS
- **strength of evidence**: PROVEN (element level)
- **remaining limitations**: element tests only
- **source file or csv**: results/raw/ppp_experiment/stage2_element_tests.txt
- **reproducible from retained data**: partly (test output retained; the test harness needs a compiled experimental library, source not shipped)

## G2: alpha=0 control on C1b (does the experimental build reproduce stock?)

- **input configuration**: C1b input + one XML comment line, FEBIO_PPP_ALPHA=0
- **source binary provenance**: experimental build febio4 85c344a7...; libfebiofluid 1bff35b8... (round 1)
- **mesh timestep**: 1688 hex8, dt 0.5, 600 steps
- **cpu allocation**: 4 cores, 150 s (job 14107), 0.167 CPU-h
- **observed result**: nodal fields differ by 6.6e-12 (w_x), 1.9e-12 (w_y), 8.3e-12 (e_f) over 601 blocks; 1179 vs 1180 iterations; F1, F2, D, flux and residuals equal to 1e-9
- **numerical acceptance criterion**: pre-registered: field difference <= 1e-8; iterations within 1 %
- **pass fail**: PASS (5/5 pre-registered criteria)
- **strength of evidence**: PROVEN (fluid-only C1b, 4 threads); SUPPORTED for the FSI domain (5-step start-up only)
- **remaining limitations**: not bitwise (thread-order noise floor 7.6e-12); alpha=0 call path includes the rim code of source change #2
- **source file or csv**: results/ppp_alpha0_control.csv; docs/reports/FEBIO_ALPHA_ZERO_CONTROL.md
- **reproducible from retained data**: partly (summary JSON; field comparison needs the Mon files)

## G3: Tiny exact Poiseuille channel, alpha=3/4

- **input configuration**: T3 straight channel, t=100 s
- **source binary provenance**: experimental build febio4 85c344a7...; libfebiofluid 1bff35b8... (round 1)
- **mesh timestep**: 1400 hex8, dt 0.5
- **cpu allocation**: 2 cores, 26 s (job 14106) vs stock 62 s (job 14105)
- **observed result**: y-alternation 6.3e-4 -> 1e-8 Pa, profile error 1.01 % -> 0.26 %, dp/dx error -0.04 % -> -0.79 %; flux at x=2 1.0000 -> 0.99264 (-0.74 %), outlet 0.99708
- **numerical acceptance criterion**: derivation prediction: flux deficit alpha (h/H)^2 Q = 0.75 %
- **pass fail**: mass conservation deteriorates (as predicted)
- **strength of evidence**: PROVEN (this case)
- **remaining limitations**: one alpha, one mesh
- **source file or csv**: results/ppp_tiny_channel_T3.csv
- **reproducible from retained data**: yes (table)

## G4: C1b fixed-plate fluid-only with alpha=3/4 (one predeclared alpha)

- **input configuration**: C1b + comment, FEBIO_PPP_ALPHA=0.75, rim update once per step
- **source binary provenance**: experimental build febio4 85c344a7...; libfebiofluid 1bff35b8... (round 1)
- **mesh timestep**: 1688 hex8, dt 0.5, 600 steps
- **cpu allocation**: 4 cores, 48 s (job 14108), 0.053 CPU-h
- **observed result**: D 0.4207 -> 9.4e-5 Pa (x4500); F1-F2 0.8196 -> 0.0054 N/m (F1 -4.6637 -> -5.1335, F2 -5.4833 -> -5.1388); BUT Q_out/Q_in-1 +9e-5 -> -1.55 %; smooth pressure regional shifts -5.3..+6.6 %; velocity L2 3.7 %; reaction sum +12 %; momentum residual 8.4e-7 (satisfied); 776 vs 1180 iterations
- **numerical acceptance criterion**: pre-registered S1-S8: S3 smooth fields <=2-3 %, S5 mass <=1e-3, S6 reactions within 2 %
- **pass fail**: FAIL (S3, S5, S6); S2 mode suppression, S4, S7, S1, S8 pass
- **strength of evidence**: FAILED (as a correction); mode suppression PROVEN for this case
- **remaining limitations**: agreement F1~F2 is a consequence of parity removal, not evidence about the true force; no value claimed correct; scaling of the defect with alpha and (h/H)^2 NOT TESTED
- **source file or csv**: results/ppp_round1_C1b.csv; results/ppp_round1_comparison_stock_vs_alpha075.json; docs/reports/FEBIO_STABILIZED_FLUID_TEST.md
- **reproducible from retained data**: partly (summary JSON/npz; re-analysis needs Mon_all, HPC)

## G5: Diagnosis of the Newton stall at alpha=3/4 (lagged rim pressure)

- **input configuration**: T1 tiny channel with / without rim, max_refs 60
- **source binary provenance**: experimental build febio4 85c344a7...; libfebiofluid 1bff35b8... (round 1)
- **mesh timestep**: 1400 hex8, 10 steps
- **cpu allocation**: 2 cores, 2-30 s (jobs 14096, 14098, 14099)
- **observed result**: alpha 3/4 as first built: max reformations at step 1 (contraction 0.55); without rim quadratic (4 it/step); max_refs 60 converges at 43 it/step; remedy (rim once per step) gives quadratic convergence and the same converged state
- **numerical acceptance criterion**: convergence
- **pass fail**: initial build FAILED; remedy PASS
- **strength of evidence**: PROVEN (tiny channel)
- **remaining limitations**: source change #2 is a second modification
- **source file or csv**: docs/reports/FEBIO_EXPERIMENTAL_BUILD.md sec. 3
- **reproducible from retained data**: partly (logs gz)

## G6: FSI domain (coupled M1 start-up, 5 steps)

- **input configuration**: T4 stock / alpha 0 / alpha 3/4, t_end 12.5 s
- **source binary provenance**: experimental build febio4 85c344a7...; libfebiofluid 1bff35b8... (round 1)
- **mesh timestep**: 1688 + 160, dt 2.5
- **cpu allocation**: 2 cores, 9-15 s (jobs 14109-14111)
- **observed result**: alpha=0 equals stock to 9e-12 (34 = 34 it); alpha 3/4 converges (27 it), transient differs (start-up)
- **numerical acceptance criterion**: alpha=0 equals stock
- **pass fail**: PASS (alpha=0); alpha>0 untested beyond 5 steps
- **strength of evidence**: SUPPORTED (5 steps)
- **remaining limitations**: no steady coupled run
- **source file or csv**: results/raw/ppp_experiment/T4_fsi_stock_vs_a0.json
- **reproducible from retained data**: yes (JSON)

## G7: Stage 5: one coupled M1 run with the stabilisation; M2 fluid pair

- **input configuration**: planned, gated on Stage 4
- **source binary provenance**: not applicable
- **mesh timestep**: not applicable
- **cpu allocation**: 0 (not run)
- **observed result**: not run: Stage 4 failed
- **numerical acceptance criterion**: Stage 4 gate
- **pass fail**: NOT RUN
- **strength of evidence**: NOT TESTED
- **remaining limitations**: coupled benefit of any stabilisation unknown; the 1.5 % conservation defect exceeds the effect being studied (0.1-5 % on the force)
- **source file or csv**: docs/reports/FEBIO_STABILIZATION_FINAL_DECISION.md
- **reproducible from retained data**: no

## H1: Assembled tests of the flux-preserving (redistribution) form

- **input configuration**: random 7x4x1 perturbed mesh, 16 constrained nodes (outlet column, rim, pins), alpha 1e-3; 40 random elements/masks
- **source binary provenance**: experimental build febio4 85c344a7...; libfebiofluid 6627dadd... (round 2)
- **mesh timestep**: not applicable
- **cpu allocation**: login node
- **observed result**: free-row sum of the stabilisation residual: old form -9.7e-6, new form -2.9e-20 (2e-15 of max|S|); constrained rows 0; dS/de FD 1.2e-12; dS/dx FD 3.2e-9; alpha=0 bitwise unchanged
- **numerical acceptance criterion**: machine-precision balance; FD <= 1e-8
- **pass fail**: PASS
- **strength of evidence**: PROVEN (machine precision)
- **remaining limitations**: global row balance only; says nothing about nodal column flux
- **source file or csv**: results/raw/ppp_experiment/round2_conservation_tests.txt
- **reproducible from retained data**: partly (output retained; harness needs the compiled experimental library)

## H2: P0: rim-once at alpha=0 reproduces stock

- **input configuration**: T5 tiny channel, t=100 s
- **source binary provenance**: experimental build febio4 85c344a7...; libfebiofluid 6627dadd... (round 2)
- **mesh timestep**: 1400 hex8, dt 0.5
- **cpu allocation**: 2 cores, 40 s (job 14118)
- **observed result**: 5.9e-12 relative on all Mon files; 274 vs 492 iterations
- **numerical acceptance criterion**: <= 1e-8
- **pass fail**: PASS
- **strength of evidence**: PROVEN (this case)
- **remaining limitations**: rim-once alone changes the iteration path, not the converged state
- **source file or csv**: results/ppp_round2_criteria.csv
- **reproducible from retained data**: yes (JSON)

## H3: Tiny exact Poiseuille channel, alpha=1e-3 (flux-preserving form; alpha declared from mu x lowest smooth Schur eigenvalue, not from forces)

- **input configuration**: T5 tiny channel, t=100 s, pre-registered P1-P6, stop at first failure, no second alpha
- **source binary provenance**: experimental build febio4 85c344a7...; libfebiofluid 6627dadd... (round 2)
- **mesh timestep**: 1400 hex8, dt 0.5
- **cpu allocation**: 2 cores, 36 s (job 14119), 0.020 CPU-h
- **observed result**: P1 nodal column flux 4.86e-4 (> 1e-4) FAIL; P2 dp/dx -0.29 % (> 0.2 %) FAIL; P4 y-alternation 3.0e-4 Pa at x=12 (stock 6.2e-4: only about 2x lower; criterion 6.3e-5) FAIL; P3 profile 2e-4 PASS; P5 momentum 1.5e-9 / continuity 1.4e-7 PASS; P6 274 it (1.37/step) PASS
- **numerical acceptance criterion**: pre-registered P1-P6
- **pass fail**: FAIL (3 of 6 + P0 control passes)
- **strength of evidence**: FAILED
- **remaining limitations**: cause (x-sawtooth of nodal flux amplified about 300x) is an inferred mechanism, UNRESOLVED; the pair-sum flux (what the continuity rows constrain) is conserved to 2e-5
- **source file or csv**: results/ppp_round2_criteria.csv; results/ppp_round2_T5_poiseuille.csv; docs/reports/FEBIO_FLUX_PRESERVING_PROJECTION.md
- **reproducible from retained data**: yes (JSON; Mon files HPC)

## H4: Fluid-only cavity test with the round-2 form

- **input configuration**: not run
- **source binary provenance**: not applicable
- **mesh timestep**: not applicable
- **cpu allocation**: 0 (not run)
- **observed result**: not justified after H3
- **numerical acceptance criterion**: gated
- **pass fail**: NOT RUN
- **strength of evidence**: NOT TESTED
- **remaining limitations**: FSI domain of round 2 compiled and element-tested only
- **source file or csv**: docs/reports/FEBIO_FLUX_PRESERVING_PROJECTION.md sec. D
- **reproducible from retained data**: no

