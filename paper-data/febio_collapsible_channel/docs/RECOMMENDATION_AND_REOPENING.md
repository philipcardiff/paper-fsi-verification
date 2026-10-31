# Recommendation and reopening criteria

**Recommendation: do not continue the collapsibleChannel coupled comparison with stock FEBio v4.13.** Archive this package as diagnostic evidence. A full coupled run would only quantify the discrepancy; it could not verify the solver by itself.

**Why.** The load that FEBio applies to the wall (nodal wall-row pressure) differs from the element-mean pressure by -8.6 % (steady first-mode load) and by +19 to +65 % in the dynamic loading phase, through a wall-normal row-parity mode of the equal-order hex8 fluid operator; the wall load actually applied matches the nodal one (free-body, 0.1-1 %). A meaningful coupled comparison needs a fluid operator without the wall-normal pressure parity (a stabilised or mixed-order fluid element). **Stock v4.13 does not provide one** (earlier campaign tests of experimental stabilisations failed their criteria and are not proposed; see PRs #27/#28 packages). A plugin or a source change is outside this pilot and is not proposed here.

**Costs (NOT run; from measured D1 speed 6.9 s per step on 4 cores, pilot report section 8):**
| item | allocated CPU-h |
|---|---|
| full 140-step coupled run at the pilot mesh (to record trough, period and decay as the parity-biased FEBio answer, labelled DIAGNOSTIC) | about 1.1 |
| one 2x finer-in-y fluid mesh F1 (does the y-parity fall under refinement) | about 0.3 |
| time-step study (dt 0.0125 and 0.00625) | about 2.2 and 4.5 |
| K sweep | about 1.1 each |
| level 2 (4x cells, dt/2) | about 9-11 |
The first two are the cheapest remaining information within stock FEBio (about 1.4 CPU-h) and are recommended only if someone wants the parity-biased answer on record.

**Reopening criteria (any of):**
1. A fluid formulation without wall-normal pressure row parity is available (stock or reviewed plugin) and passes a fixed-shape F1 test with |F1 - F2|/|F2| <= 3 % for F_y and G on this collapsed shape and a dynamic-step spread test (the A2 criterion, <= 3 % where |F2_y| >= 100 N/m).
2. A reference for the physically correct wall pressure on the discrete interface (e.g. a converged solution) becomes available, so F1 vs F2 can be judged against it.
3. Maintainers state that a parity-biased FEBio trough/period/decay record is wanted; then run the 1.1 CPU-h full run and the 0.3 CPU-h finer-y F1, labelled DIAGNOSTIC.
When reopening, use the 320x8 hex8 wall (or a locking-free wall) with conformal fluid columns, K = 1e7, and a second mesh level and time step before using the word verification.

**Relation to PRs #27 and #28 (same pressure-parity mechanism class; not shown to be identical).** PR #27 (cavityFlexibleBottom) and PR #28 (Hron-Turek FSI1) record the same kind of finding: the nodal wall-row pressure that FEBio's fluid-FSI traction applies carries a row-parity component that the element-mean pressure does not (cavity: F1 about 5-8 % below the solids4foam force and F1 - F2 equal to the parity part; FSI1 rigid-flag test: parity 16-23 Pa on 25-30 Pa, drag +15.7 % by F1). **Demonstrated here:** in this geometry, a wall-normal row parity in the fixed-shape and dynamic fluid data, the F1 vs F2 spread, and that the wall load equals F1. **Hypothesis:** that the cause is the same unstabilised equal-order Q1/Q1 operator in all three geometries; the mechanism-level argument is in the earlier packages, and this package does not show that the cause is identical in each geometry (for example the dependence of the amplitude on the smooth load vs viscous forcing, and the refinement behaviour, were not measured here). PR #26 (ring) and PR #29 (radialTorsionalRing torsion) are not part of this relation: the ring flow is smooth and parity-free, and the torsion QoI is parity-immune at leading order. No existing PR or branch covers collapsibleChannel (checked with `gh`, read-only, 2026-10-31).
