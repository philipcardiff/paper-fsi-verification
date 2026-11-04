# Limitations, unresolved questions and reopening criteria

## Limitations
1. **The reference is our own draft.** The static torsion energy, the coupled reduced model and the fixed-geometry fluid solver are by the same author as the proposal; no second implementation of the coupled model exists. The static solid reference was exercised by an independent code lineage (agreement at second order), and the fixed-geometry fluid reference was cross-checked here by a finite-difference BDF2 solver (1.3e-5 of max |T_B|); the coupled reference was not compared with anything.
2. **One time step only (dt = 0.01 s) and no K_f study.** The temporal part of the aligned velocity error (0.19 % at M1) is not separable from the spatial part; acoustic and compressibility effects (estimated below 1e-5) were not measured. Hence no VERIFIED AGREEMENT for any fluid quantity.
3. **The dt/2 alignment is post hoc.** The gated, as-run RT-B criteria failed. The aligned numbers are the correct comparison only because the offset was proven from source and from a diagnostic run; they are not a re-definition of the gates.
4. **Reconstructed torques.** T1-T3 are own algebra on logged nodal fields (FEBio logs no fluid-FSI wall reaction force), validated by the exact reference but not solver output. T3 is a consistency torque, not an applied load.
5. **RT-B parity result is flow-specific.** An axisymmetric swirl has no odd-even pressure mode to excite; the passed parity gate does not speak to FEBio's equal-order hex8 parity behaviour in other flows (documented for the cavity and FSI1 packages).
6. **No coupled result.** RT-C M1 diverged; M2, M3, dt and K runs were not made; the coupled ALE-with-rotating-mesh path is tested only by M0 smoke cases (64 fluid cells) and one failed M1 run.
7. **hex8 only, stock settings only.** Tolerances vtol = ftol = 1e-6, etol = rtol = 0, min_residual = 1e-18 (amendment A1), E_mesh = 1 Pa not varied.
8. **Source-reading claims** (dt/2 convention, element-mean traction, energy-conserving correction) are code reading of the v4.13 tag with the quoted line numbers; the dt/2 claim is also demonstrated by data, the other two are not isolated by a dedicated experiment (traction: confirmed by an advance quantitative prediction; correction: one diagnostic).
9. The generator was edited after the plan was registered (provenance: `../provenance/INCONSISTENCIES.md`); all gated inputs nevertheless regenerate bit-identically from the shipped generator.

## Unresolved
- Why RT-C M1 diverges at rho_inf = 0.9 from a round-off-level drive (inexact tangent of the interface or mixture block, near-singular massless skeleton block at M1: untested).
- The 4-10 Newton iterations per step of the coupled smoke runs (SVK at rho_inf = 1 and MR at rho_inf = 0.9), against 1-3 for the mixture alone.
- Whether the coupled rotation is affected by the first-order interface traction as the analytical expectation says (alpha_B 8 % and 4 % at M1 and M2): untested.
- Whether the time-label offset of the fluid module is intended (a question for the FEBio developers).

## Reopening criteria
Reopen the coupled test (RT-C) if any of the following becomes available, at the pre-registered gates of `PREREGISTRATION.md` (alpha_B, H_f at most 5 % / 2 %, T_B and sigma_rt at most 10 % / 4 % at M1 / M2, orders at least 0.8 and 1.5, at most 12 Newton iterations on average):
1. A FEBio version or setting in which the coupled Mooney-Rivlin smoke case converges at rho_inf = 1 and the coupled M1 run does not diverge from round-off-level drives (for example a documented fix of the energy-conserving correction or of the tangent), or a documented, justified starting amplitude (a ramp that does not begin at round-off) agreed before the run.
2. A second independent implementation of the proposal's coupled reference (so that the draft's reference is itself tested), which would also make a coupled FEBio comparison more informative.
3. A developer answer on the intended load-curve time convention, or a decision to compensate it in inputs by default.
Reopen the fixed-geometry fluid test for a full study if budget allows (cost estimate in `ACCOUNTING.md`): M3 (16,384 cells), a dt halving at M1 and M2, and a K_f pair, to convert the QUALIFIED AGREEMENTS into a time-and-mesh convergence statement and to test the predicted 3.1 % T2 error at M3. A claim that the first-order traction is removed would need a source change or plugin, which is out of scope here and not proposed.
