# PHASE C NUMERICAL UNCERTAINTY (slip-wall + slip-interface formulation)

Error sources are kept separate. Quantities are signed errors relative to the exact continuum reference; the screened metric is the wet/dry ratio (ratio error +0.18 % at M3). No formal convergence order is claimed from only three levels.

| source | evidence | effect on omega (wet) | effect on ratio | status |
|---|---|---|---|---|
| **Spatial** | M1/M2/M3 matched dry/wet; difference-based order 1.98 for both dry and wet | 2.759 -> 0.807 -> 0.313 % (extrapolated ~ +0.15 %) | +0.170 -> +0.178 -> +0.180 % (increments 0.0075, 0.0020) | quantified; ratio is mesh-converged to ~0.002 pp; absolute omega still has ~+0.15 % wet bias |
| **Temporal** | one paired dry/wet dt-halving at M1 | +0.026 % wet, +0.026 % dry (rhoi = 1 midpoint, phase error ~ -(omega dt)^2/12) | +9e-5 pp (cancels) | quantified at M1 only; not repeated at M3; the ratio is insensitive, the absolute frequency is not |
| **Compressibility (K)** | K 1e8 -> 1e7 at M1 | -0.0004 % | -0.0004 pp | negligible; K=1e8 adequate |
| **Iterative** | one tolerance-tightening attempt **failed**; margins at C0, M2, M3 | not quantified | not quantified | **open** (see below) |
| **Viscous / model form** | not varied in the slip family; no-slip family differs by ~6 pp at M1 | n/a | n/a | the slip configuration is an inviscid-limit formulation, not the benchmark's finite-viscosity case |
| **Mesh regulariser (E_mesh = 1 Pa)** | no-slip D4 (E_mesh 1e-3) moved omega by -0.160 % | not tested in slip family | candidate for the +0.18 % plateau (inference only) | **untested** here |
| **Initial condition** | compatible IC v2 vs v1 in no-slip family; v2 adds non-modal KE | not tested in slip family | n/a | partial |

## Iterative error: what is and is not known
- Convergence margins (achieved / required): residual, energy, displacement and velocity norms are met 3-10 orders below tolerance at every level. The **dilatation** norm is the weakest: worst case 0.48 (M1), 0.44 (M2) and **0.99 (M3)** of tolerance, median 0.045 at M3 (4 iterations per step, max 7). The M3 worst case is at the edge; a step that stops on dilatation at that level has an uncontrolled dilatation increment of order the tolerance.
- The attempted tolerance-tightening run (`runs/phaseC/C1_tol_wet_f1`: dtol = vtol = ftol = 1e-8, etol = 1e-10, rtol = 1e-14) **failed to converge at step 1**. The dilatation increment has a noise floor of about 5e-8 in this formulation, so tolerances below roughly ftol ~ 2e-7 (relative form) cannot be met. The failure is preserved, counted against the budget, and is a floor finding, not a defect claim. Consequently the iterative error is bounded only indirectly: the frequency is insensitive to anything that shifted between M1/M2/M3 despite the tolerance margin varying (0.48 -> 0.99), the 6 % pressure noise is mesh-independent, and the window spread of the frequency fit is 3e-5 %. This bounds the iterative contribution to omega at well below the other errors but **is not a measurement** of it.

## Interpretation limits
- The temporal check used M1 only; the cancellation in the ratio (9e-5 pp) can be expected to persist but was not verified at M3.
- The +0.18 % plateau has no measured cause. The M1 value (0.170 %) is already within 0.01 pp of the M3 value; the offset is a systematic, mesh-independent term (regulariser stiffness, residual viscous or leakage effects, or the IC), not discretisation of the continuum equations.
- Third-level extrapolations use p = 2 assumed; they are indicative (wet +0.148 %, dry -0.032 %), not Richardson-validated.

## Combined statement
For this formulation the numerical uncertainty on the ratio from spatial, temporal and K effects is below 0.01 pp, while the systematic offset of +0.18 % remains unexplained. Iterative error is bounded but not quantified. Viscous/model-form effects are not covered by these numbers.
