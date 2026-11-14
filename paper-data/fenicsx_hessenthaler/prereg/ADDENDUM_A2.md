# Addendum A2 (2026-11-11), registered BEFORE any run with the changed mesh motion
Numerical (not physical) change: the fluid-mesh motion equation. Variant "total" (used in all runs so far, kept and reported):
  int k grad_X d : grad_X psi dX = 0 on the fluid cells (harmonic extension of the TOTAL displacement in the reference configuration).
Variant "incr" (new): harmonic extension of the displacement INCREMENT dd = d1 - d0 on the previous-step configuration x0 = X + d0 (F0 = I + grad d0):
  int k inner(grad_X(dd) F0^-1, grad_X(psi) F0^-1) dX = 0   (k = mean volume / reference cell volume as before; the current-configuration J0 cancels).
Same boundary conditions (d = 0 on inlets, outlet, wall, clamp; interface d rows from the kinematic condition). The fluid, solid, interface and
boundary-condition physics are unchanged; the choice of the arbitrary ALE map changes only discretisation error. Option mmotion=incr in benchmarks/hessenthaler/run_fsi.py.
Pilot protocol (fixed now): Gc dev level, 16 ranks, dt = 0.02 s, t = 0 -> 1.0 s. Pilot passes iff, for every step, Newton <= 12 iterations,
min fluid detF >= 0.2, min solid detF > 0.3, |Qout - Qin|/Qin <= 1e-3, no step failure that needs dt reduction. Any breach stops the run and is reported.
Production (Gmi) only if the pilot passes, within the remaining CPU-h budget (cap ~400 CPU-h total; ~160 remained at registration), with the same gates and stop rules;
dt <= 0.01 and dt/2 control are requested by the coordinator but are limited by budget: the achieved t_end and dt are reported as they are, no extrapolation.
Thresholds of PREREG_HESSENTHALER_PHASE_I.md are unchanged.
