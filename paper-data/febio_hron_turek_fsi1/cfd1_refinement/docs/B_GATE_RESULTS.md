# Campaign B, gate 1: Hron-Turek CFD1 rigid-flag refinement (preregistration `O2_PREREGISTRATION_B.md`, hash in `PREREG_B.sha256`)

Question: does the nodal-pressure wall force F1 (what FEBio's `fluid-FSI traction` applies) converge under refinement towards the Featflow converged CFD1 reference (drag 14.2929, lift 1.11905 N/m per metre depth)? Stage-C level 1 had F1 +15.73 % in drag and +157 % in lift.
Configuration: stage-C model unchanged (FEBio `fluid`, hex8 slab, rho 1000, mu 1, K 1e6, 80 steps to 40 s), factors 1.25 (L1, existing), 1.875 (L1.5, 18,606 cells), 2.5 (L2, 33,368 cells); junction pins by the stage-C rule (neighbour mean from the unpinned level run). All runs normal termination, 80/80 steps, mass imbalance <= 3.3e-4, steady to 1e-11.

## Results
| run | cells | F1 drag err % | F1 lift err % | F2 drag/lift err % | F3 drag/lift err % | a_t top/bottom/cyl (Pa) | mass imb | drag spread | junction p (Pa) |
|---|---|---|---|---|---|---|---|---|---|
| L1 (1.25) pinned [stage C] | 8330 | +15.73 | +157.0 | -2.12 / -1.12 | -17.85 / +4101.0 | 4.11 / 2.63 / 14.87 | -1.2e-04 | 1.3e-11 | [-81.6042241147, -123.87086089400002] |
| L1.5 (1.875) pinned (rule) | 18606 | -0.83 | +13.7 | -1.62 / +0.06 | +6.88 / -11.8 | 1.09 / 1.23 / 6.54 | 1.2e-04 | 1.8e-12 | [9.89986964937, 4.45979062992] |
| L2 (2.5) pinned (rule) | 33368 | +22.89 | -26.9 | -1.19 / +0.08 | +56.68 / -93.5 | 13.91 / 13.86 / 69.08 | 3.3e-04 | 1.8e-12 | [50.7441321248, 46.2058930776] |
| L1.5 pins = L1 values | 18606 | +0.02 | +15.1 | -1.30 / +0.60 | +6.88 / -11.8 | 1.08 / 1.22 / 6.54 | 1.2e-04 | 1.2e-12 | [-81.6042241147, -123.87086089400002] |
| L1.5 pins = 0 | 18606 | -0.77 | +13.4 | -1.60 / -0.02 | +6.88 / -11.8 | 1.09 / 1.23 / 6.54 | 1.2e-04 | 9.2e-13 | [0.0, 0.0] |
| L1.5 pins = 2x rule | 18606 | -0.88 | +13.9 | -1.64 / +0.14 | +6.88 / -11.8 | 1.09 / 1.23 / 6.54 | 1.2e-04 | 1.7e-12 | [19.799739298699997, 8.91958125984] |
| L1.5 unpinned | 18606 | -1.50 | +13.3 | -1.88 / -0.07 | +6.88 / -11.8 | 1.09 / 1.24 / 6.54 | 1.2e-04 | 1.5e-03 | [93.8753338653, 97.0685387156] |
| L2 unpinned | 33368 | +23.02 | -27.0 | -1.14 / +0.03 | +56.68 / -93.5 | 13.91 / 13.86 / 69.08 | 3.3e-04 | 2.7e-04 | [26.2262979177, 26.590079743300002] |

## Preregistered criteria
| Criterion | Result |
|---|---|
| B1 monotone decrease of the F1 drag and lift errors | **FAILED**: drag +15.73, -0.83, +22.89 %; lift +157, +13.7, -26.9 % |
| B2 F1 drag error at L2 at most 0.7 of L1 | **FAILED**: 22.89 / 15.73 = 1.46 |
| B3 flag odd-even amplitude falls by 1.5 from L1 to L2 | **FAILED**: top 4.11 -> 13.91 Pa, bottom 2.63 -> 13.86 Pa (cylinder 14.9 -> 69.1 Pa); at L1.5 it did fall (1.09 / 1.23 / 6.5 Pa) |
| B4 soundness at each level | met |
No coupled FSI1 run was made, as registered (B1-B4 had to hold first).

## Pin sensitivity (exploratory addendum B2, `O2_PREREGISTRATION_B2.md`)
At L1.5 the F1 drag error is -0.83 % (rule pins), -0.77 % (pins 0), -0.88 % (twice the rule values) and +0.02 % with the level-1 pin values (-81.6, -123.9 Pa), i.e. up to 0.85 pp across a pin range of 200 Pa; the lift error changes by at most 1.8 pp. The stage-C criterion (<= 0.3 % of the drag) is not met for the level-1 values, so the pins matter at the 1 % level, but they do not explain the level dependence: L2 pinned and unpinned agree (+22.89 / +23.02 %), and L1.5 unpinned gives -1.50 %.

## Reading
- The smooth reconstruction F2 (element-mean pressure, an O(h)-biased estimate) converges monotonically (drag error -2.12, -1.62, -1.19 %; lift -1.12, +0.06, +0.08 %), so the physical pressure field underneath is converging towards the reference.
- The nodal-pressure force F1 does not: its error changes sign and is not reduced by refinement, and the wall odd-even pressure amplitude is not reduced either (cylinder 14.9, 6.5, 69.1 Pa). The level-dependence is therefore not a resolution effect in the sense of the ring (whose nodal n=2 force converged at first order, -2.1, -1.06, -0.53 %), and a coupled FSI1 with FEBio's nodal traction cannot be verified on these meshes.
- F1 at L1.5 (-0.8 % drag) is an accident of the mode amplitude at that level, not evidence of convergence.
Closeness of F2 to the reference is not evidence that F2 is the correct interface load (it has an intrinsic O(h) bias even for exact pressure).

---
# Addendum B3 (exploratory, `O2_PREREGISTRATION_B3.md`): parity of the block cell counts, and the consistent reaction F4 on eight meshes

All eight meshes: normal termination, 80/80 steps, steady to 1e-11, mass imbalance <= 7e-4 (3.3e-4 at the worst before B3; B3 meshes 2.5e-5 to 6.9e-4). Pins: stage-C rule (levels 1.25, 1.875, 2.5) or 0 Pa (B3 meshes; the pin sensitivity at level 1.875 is <= 0.9 pp in drag).
| factor | cells | odd-count blocks | F1 drag / lift error % | F2 drag / lift error % | **F4 (consistent reaction) drag / lift error %** | cylinder odd-even a_t (Pa) |
|---|---|---|---|---|---|---|
| 1.25 (stage C, pins by rule) | 8330 | 19 | +15.73 / +157 | -2.12 / -1.12 | **+0.31 / -0.61** | 14.9 |
| 1.5 (pins 0) | 11996 | 6 | +23.70 / -71 | -1.64 / +0.15 | **+0.27 / +0.59** | 65.6 |
| 1.625 (pins 0) | 14191 | 29 | +12.50 / +661 | -1.71 / +0.01 | **+0.23 / +0.25** | 19.0 |
| 1.875 (pins by rule) | 18606 | 18 | -0.83 / +14 | -1.62 / +0.06 | **+0.16 / +0.29** | 6.5 |
| 2.0 (pins 0) | 21344 | 0 | +18.85 / -1 | -1.75 / +0.07 | **+0.21 / +0.37** | 29.9 |
| 2.375 (pins 0) | 29951 | 29 | +12.72 / +689 | -1.28 / +0.06 | **+0.13 / +0.21** | 16.5 |
| 2.5 (pins by rule) | 33368 | 6 | +22.89 / -27 | -1.19 / +0.08 | **+0.12 / +0.28** | 69.1 |
| 3.0 (pins 0) | 48024 | 26 | +3.68 / -292 | -1.27 / +0.14 | **+0.08 / +0.36** | 12.2 |

Preregistered predictions: **P1 met** (a_t cylinder at 1.5, 6 odd, 65.6 Pa, is 3.5 times that at 1.625, 29 odd, 19.0 Pa); **P2 failed** (a_t at 2.0, 0 odd, 29.9 Pa, is 4.6 times that at 1.875 but only 1.8 times that at 2.375, 16.5 Pa; both were required); **P3 met** (16.5 Pa at 2.375 against 69.1 Pa at 2.5); **P4 failed** (|F1 drag error| on the odd-rich meshes is 12.5, 12.7 and 3.7 %, not at most 5 % and not decreasing). Parity therefore correlates with the odd-even amplitude (the two meshes with 6 odd counts have 66-69 Pa, the others 6-30 Pa) but does not control the force error: the nodal-pressure force F1 is +12 to +24 % in drag on six of the eight meshes and its lift error is unbounded (the lift is a small difference of large terms).
**The consistent nodal reaction F4** (the fluid's own variationally consistent reaction at the obstacle nodes, an independent implementation) is within +0.08 to +0.31 % in drag and +0.2 to +0.6 % in lift on every mesh and its drag error falls with refinement (+0.31, +0.27, +0.23, +0.16, +0.21, +0.13, +0.12, +0.08 %); the element-mean F2 sits at -1.2 to -2.1 % (its known O(h) bias).

## Mathematical identification (derivation, established; arithmetic check in `reaction_check.py` on stock ring output)
FEBio's fluid momentum row is ff = sv:grad N + (grad p) N with strong pressure gradient (`FEFluidFSIDomain3D.cpp`), and at a node with w fixed that row is discarded. Per element, boundary integral of N_a p n dS = integral (N_a grad p + p grad N_a) dV. The force the traction load applies is the boundary integral; the variationally consistent force on the obstacle is the weak form, whose pressure part is integral p grad N_a. Their difference is exactly the strong-form residual R_strong,a = integral N_a grad p dV (plus viscous and inertial terms) at the constrained nodes: F1 - F4 = sum of R_strong over the obstacle nodes. The momentum equation is not enforced at the constrained nodes, so R_strong is the discretisation residual there; it is small where the layer is in balance (ring: the strong reaction's work is 0.12 % of the ring energy exchange) and O(10 %) of the drag on CFD1 (flow with a stagnation point and wake), where the equal-order nodal pressure carries alternating modes that the strong gradient integral does not see but the boundary integral does. The null-mode property is shown by the Q1/Q1 operator audit (`q1q1_parity_operator.json`: the continuity/momentum coupling has exact alternating null modes, 8 zero eigenvalues of the pressure Schur complement on 8x8 and 9x9 meshes alike, controlled only by the compressibility term, of relative size gamma/(K dt)).
**Established:** the CFD1 force discrepancy of the stock traction is not a resolution effect and not a pin or parity effect; it is the inconsistency between the boundary-integral traction and the weak-form reaction, carried by the nodal-pressure alternating mode. **Not established:** that the formulation as a whole is wrong (the underlying pressure field converges: F2 and F4 do), and a stock-solver measurement of the applied load (it was reconstructed from FEBio's own formula; FEBio's coupled load was not measured on CFD1).
