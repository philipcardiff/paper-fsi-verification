# PHASE D ITERATIVE-ERROR INVESTIGATION (moderate, M2, E_mesh = 1e-3 Pa)

Script: `tools/dilatation_floor.py`, `tools/iter_compare.py`; run `runs/phaseD/D6_iter_moderate_M2` (Slurm 13717, 1 core, **failed at step 75 of 400, preserved**). One simulation of the six authorised.

## 1. Dilatation floor, investigated from logs before running (no new simulation)
FEBio prints SQUARED norms: CURRENT = |increment|^2, REQUIRED = tol^2 x |accumulated step increment|^2.

| log | tol | final-iteration \|inc\|^2 (min / median / max) | accumulated step \|inc\|^2 (median) | final relative increment (median / max) |
|---|---|---|---|---|
| C2 M2 (E_mesh 1e-3) | 1e-6 | 1.2e-31 / 3.9e-29 / 2.7e-28 | 3.0e-15 | 1.2e-7 / 8.7e-7 |
| C1 M2 (E_mesh 1) | 1e-6 | 4.0e-32 / 4.2e-29 / 3.3e-28 | 3.0e-15 | 1.2e-7 / 6.6e-7 |
| C1 M3 (E_mesh 1) | 1e-6 | 1.1e-30 / 1.9e-28 / 1.7e-27 | 4.1e-15 | 2.1e-7 / 1.0e-6 |
| failed Phase C tolerance run, step 1 (25 iterations) | 1e-8 | plateau 9.3e-31 (1e-31 to 1e-29 over iterations 3-25) | 5.2e-17 | 1.3e-7 |

- The dilatation increment does not decay below |inc|^2 of about 1e-30 to 1e-29 (|inc| ~ 1e-15) however many Newton iterations are taken (25 in the failed run): it is a flat arithmetic floor, not slow convergence. The other norms keep falling by 5-10 more orders in the same iterations.
- Relative to the accumulated step increment the floor is about 1e-8 to 1e-7 (largest at step 1 where the pressure change is smallest, e_f ~ 4e-9). The Phase C tolerance of 1e-8 was therefore 10x below the floor even on typical steps, and the baseline 1e-6 is only 3-10x above it on the worst steps (margin 0.4-0.99 of the tolerance, M3 0.99, strong-added-mass M1 0.99).
- A plausible origin (consistent in magnitude, not proven): e_f = J - 1 ~ 1e-9 stored near J = 1 in double precision has resolution ~2e-16 per node, so the increment norm over the nodes is ~1e-15. If that is the cause, the relative floor scales as 1/|e_f|, i.e. it would fall ~10x at K = 1e7 (larger e_f for the same pressure) and the dilatation tolerance could then be tightened ~10x. Untested.
- Conclusion: the dilatation increment criterion cannot be meaningfully tightened below about 5e-7 at K = 1e8. ftol was therefore kept at 1e-6 in D6.

## 2. D6 control strategy and outcome
Changed relative to the C2 M2 baseline (one `.feb` diff: 4 numbers on two lines): `dtol` 1e-6 -> 1e-8, `vtol` 1e-6 -> 1e-8, `etol` 1e-8 -> 1e-12, `rtol` 1e-10 -> 1e-14; `ftol` kept at 1e-6. Everything else byte-identical.

- Steps 1-74 (0.745 periods) converged in 3.03 iterations/step, **identical to the baseline**, and the probe histories (D(t), velocity) and the first zero crossing are identical to the printed digits: `max|D_tight - D_base| / max|D_base| = 0`. Tolerances that are 100x (d, v) and 10^4x (e, r) tighter were already satisfied at the iteration where the baseline stopped, because those increments collapse by 5-10 orders in the final Newton iteration: they were not the binding stopping criteria.
- Step 75 (t = 2.611 s) failed ("Max nr of reformations reached", 25 iterations): at the plateau the residual norm^2 sat at 2.8e-18 to 3.7e-18 against a required 3.6e-19 (initial 3.6e-5, rtol 1e-14), and the displacement increment^2 at 1.5e-26 to 3.3e-26 against a required 1.44e-27 (accumulated step increment^2 1.44e-11, dtol 1e-8). Both are floors: displacement increment floor ~4e-8 of the step increment (dtol must be >= ~5e-8 at this step) and residual floor ~8e-14 of the initial (rtol must be >= ~2e-13).
- Mistake in my control design (disclosed): I judged the feasible tightening from the minima of the final-iteration values in the baseline log instead of the plateau at the worst step. At step 75 the step increment is small (small ring motion near a velocity turning point), so the absolute floor is a larger fraction of it, and 100x on dtol was beyond the floor (about 27x was feasible there). The step-1 plateau used in my estimate was not representative.
- Result: no frequency can be extracted from 75 steps (one zero crossing at 1.730 s, identical in both runs), so the tightened run does not give a Delta-omega. **The iterative error in the ring frequency was not isolated within the six-run budget.**

## 3. What can and cannot be said about the iterative error
- Supported: (i) the d/v/r/e stopping criteria were not binding over steps 1-74: 100x-10^4x tighter tolerances give bit-identical probe histories and iteration counts. (ii) The dilatation increment is at its arithmetic floor (flat, independent of iteration count) at the baseline stop, so further Newton iterations do not change the solution beyond ~1e-15 in e_f (~1e-7 Pa at K = 1e8, i.e. ~3e-7 of the peak pressure 0.4 Pa). (iii) The frequency-fit window spread is ~3e-5 %, amplitude drift per period ~3e-6 (moderate M2).
- Not supported: a measured Delta-omega from tolerance changes over the full run; any statement on the dilatation criterion below 1e-6; K = 1e7 combined tolerance sensitivity.
- Residual risk: the worst-step dilatation margins at the baseline tolerance are 0.76 (M2, E_mesh 1e-3), 0.99 (M3 E_mesh 1, strong M1). The solution is at the arithmetic floor there, so the error attributable to the stopping rule is bounded by the floor, but this is an argument, not a measurement.

## 4. Proposed single follow-up (not run; needs your approval and counts as a seventh run)
Moderate M2, E_mesh = 1e-3, `dtol` = `vtol` = 5e-8 (the largest feasible tightening of the displacement floor, 20x), `rtol` = 1e-12, `etol` = 1e-10, `ftol` = 1e-6 unchanged. Predicted: converges (rtol floor 2e-13 < 1e-12; dtol floor 4e-8 < 5e-8, thin margin, so a retry could still fail on a small-motion step; a safer choice is 1e-7). Expected cost 27 min on 1 core. Alternative that attacks the dilatation floor itself: K = 1e7 with ftol = 2e-7, but this confounds with the K effect (measured -4e-4 % in omega at M1 moderate E_mesh 1) and needs a K = 1e7 baseline.
