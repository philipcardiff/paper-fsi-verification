# Equivalence of the FEBio model with the proposed problem (torsion member)

The complete pre-run audit (written before any FEBio run) is archived in `archive/RTR_EQUIVALENCE_AUDIT.md`. This file keeps what the classifications depend on. Labels: CS = code-source reading, CI = input read, AE = analytical estimate, MD = documentation mismatch.

## 1. The proposal as specified
Plane-strain annulus, one cell in depth (dz = 0.1 m); solid ring 0.5 <= R <= 1 m, neo-Hookean (E = 2e4 Pa, nu = 0.3, rho_s = 1000, mu = 7692.31, K = 16666.7 Pa; W = mu/2 (Ibar1 - 3) + K/4 (J^2 - 1 - 2 ln J)); fluid 1 <= r <= 2 m (rho = 1000, mu = 300 Pa s, convection on); the inner solid circle is rotated by alpha_A(t) = 0.1 S(t) sin(3t) (C-infinity ramp over 1 s); outer fluid boundary open, p = 0 and zero swirl shear stress; rest start; 0-2 s. Meshes M1/M2/M3 = 64x8x16 / 128x16x32 / 256x32x64 (azimuthal x solid radial x fluid radial), uniform radial spacing. Quantities: alpha_B, T_B = b^2 sigma_rt, sigma_rt, H_f per radian and per metre of depth.

## 2. The inner-boundary discrepancy ("free of radial displacement" versus "radially clamped"): RESOLVED
Three statements exist in the draft: the tutorial README (PR #637) says "the solid is free of radial displacement there"; the reference script docstring says "R = A: clamped radially (u = 0), prescribed rotation"; manuscript PR #25 says "radially clamped". The executed condition is defined by the input and the solver, not by prose:
- Input `0/solid/D`, patch `inner`: `type fixedRotation`, rotation axis z, origin (0,0,0), angle from `constant/timeVsAngle` (CI, PR #637 head 30ae00eb).
- Solver `fixedRotationFvPatchVectorField::updateCoeffs` (CS, solids4foam source): the patch is a `fixedValue` vector field; the prescribed value is the Rodrigues-rotated original face-centre position minus the original position, i.e. **both displacement components** of every patch face are prescribed as a rigid rotation about the origin. The rotated point stays on the circle r = A, so the radial displacement is zero.

**Resolution: the executed condition is radially clamped and rigidly rotated.** The README sentence is wrong or at best ambiguous (a free radial displacement would be an unconstrained component). FEBio equivalent: both in-plane displacement components of every inner-circle node are prescribed by the exact rotation expression; this is the same condition. This is a documentation defect of the draft (README) and is recorded as an inconsistency of the source material in `../provenance/INCONSISTENCIES.md`.

## 3. FEBio representation
| Element | FEBio v4.13 | Status |
|---|---|---|
| Geometry | full 360-degree one-layer hex8 slab, u_z = 0 (fluid w_z = 0), nodes exactly on the circles (polygonal faces, 64/128/256 azimuthal cells), conformal interface | SUPPORTED |
| Solid law | `Mooney-Rivlin` c1 = mu/2, c2 = 0, k = K, NIKE3D pressure model (hex8 uncoupled materials use FEBio's three-field domain). Continuum identity with the proposal's law by algebra (CS `FEUncoupledMaterial.h`, `FEMooneyRivlin`); discrete identity tested numerically (RT-A law) | PROVEN (algebra and numerics, 1.2e-12) |
| Rotation BC | `prescribed displacement` with a `math` expression of X, Y, t (reads the end-of-step time) | PROVEN (7e-13 m) |
| Fluid | `fluid` (RT-B) or `fluid-FSI` (RT-C) mixture, weakly compressible (K_f = 1e7 Pa, p = -K_f e_f, p/K_f about 5e-6), bulk viscosity zero, convection on | DIFFERENCE (K_f not varied) |
| Viscous stress | Newtonian, from the absolute fluid velocity (CS `FEFluidFSIDomain3D`) | PROVEN (CS) |
| Interface traction | `fluid-FSI traction`: f = gt p - sv gt with gt the covariant facet vector, p from the nodal dilatation, sv the **Gauss-point mean of the viscous stress of the whole attached element** (CS `FEFluidFSITraction::GetFluidStress`, lines 136-154 and 181-189) | PROVEN (CS); consequence first-order in h |
| Outer boundary | natural velocity (zero viscous traction) plus `zero fluid dilatation`: identical to the proposal's open condition for the swirl | SUPPORTED |
| Mesh skeleton | fluid-FSI mixture carries a neo-Hookean skeleton E_mesh = 1 Pa, nu = 0, density 0 | DIFFERENCE (stated, bounded AE below 1e-4 of the shear stress, not varied) |
| Time integration | generalized-alpha, rho_inf = 1 (RT-B, M0 smoke) or 0.9 (M0 diagnostic, RT-C); RT-A is static, monolithic Newton, dt = 0.01 s | DIFFERENCE (proposal: BDF2, partitioned IQN-ILS, dt = 0.0025 s) |
| Fixed-geometry RT-B drive | wall circle velocity = tabulated d(alpha_B)/dt of the proposal's coupled reference (1 ms table, linear load curve), omega_w(t) (-y, x) | the fixed-geometry reference sees the same table |
| alpha_B location | nodes on the circle r = B (the proposal uses the polygon face-centre radius B cos(pi/n_theta), a 6e-4-level correction at M1, not applied) | DIFFERENCE (only relevant to the coupled run) |

## 4. Verdict
For the quantities that depend on the swirl, a like-for-like comparison is achievable with stock FEBio. Differences that enter every classification: K_f = 1e7 (declared, not swept), hex8 equal-order fluid with element-mean interface viscous traction (first-order), monolithic generalized-alpha at dt = 0.01 s with no time-step study, three-field solid domain (tested by RT-A), the mesh skeleton, no s4f data reused.
