# PREREG_COL1: collapsibleChannel solver development, structural preflight (Agent COL-1), written 2026-11-06 BEFORE any run

Model (unchanged, from COLLAPSIBLE_AUDIT.md 769d038): plane-strain TL St Venant-Kirchhoff continuum wall, E = 4.55e8, nu = 0.3, rho_s = 1, length L = 10 (channel 5<x<15), thickness h = 0.05, whole end faces clamped, external pressure on the OUTER face.

## Load definition (T1)
Static external pressure p = 200 Pa (the plateau of the cosine ramp), TRUE pressure acting on the DEFORMED outer surface (follower, Nanson: residual term + p * int_top J F^{-T} N . v dS_0, force per unit depth, pointing into the channel i.e. -y). Measure: force per unit deformed outer-face length; the oomph beam applies per unit deformed MIDPLANE length (differs by (1 + kappa h/2), mean curvature of a clamped-clamped wall integrates to zero, so total force agrees to O(h/L)^2; counted as model form, not corrected).
Monitored quantity: u_y of the wall bottom face (y = 1) at x_local = 5 (midpoint, channel x=10) and 2.5, 7.5 (channel 7.5, 12.5). Reference (oomph CSV): -0.14362832 (mid), -0.09645925 (quarters).

## T1 ladder (all runs listed here; nothing else may be quoted as a gate)
- Q2 (quadrilateral, degree 2), ny = 8, nx = 80, 160, 320 (primary ladder, as in audit P1).
- Q3, ny = 8, nx = 40, 80, 160.
- Thickness: Q2, nx = 160, ny = 2, 4, 8; Q3, nx = 80, ny = 2, 4, 8.
- Linear check: p = 0.01 Pa, Q2 160x8, vs clamped Timoshenko value p L^4/(384 D) + p L^2/(8 kappa G h) (D = E h^3/(12(1-nu^2)), kappa = 5/6). Gate G0: within 0.5%.
- Load applied in 20 equal continuation steps, SNES newtonls basic line search, rtol 1e-12; a run is VALID only if every step has SNES converged reason > 0.
- Analysis with verification/convergence.py analyse (>=3 levels, monotone, p in [0.5,3] or not accepted).
GATES (primary: Q2 ny=8 nx=320 finest): G1 midpoint within 0.3% of -0.14362832; G2 quarter points within 0.5%; G3 nx ladder monotone and (observed order >= 1.5 OR the successive differences are < 1e-6 relative, i.e. converged to noise). Thickness ladder and Q3 are reported diagnostics (locking check: thickness-ladder spread and Q2-Q3 difference at matched mesh; "no locking" = ny=4 and ny=8 agree within 0.1% and Q2/Q3 within 0.1%).
Model form vs discretisation: model-form = (converged continuum limit, Richardson on the nx ladder at ny=8 and Q3 ladder) minus beam, reported separately from discretisation error = (finest minus limit). Expected 0.01-0.05% (audit); not a gate.
STOP RULE: if G1 fails for Q2, switch element family ONCE (Q3 ladder, already in the ladder); if Q3 finest also fails G1, STOP and report (no T2/T3 claims dependent on it; T2 may still be reported as a module test since it does not depend on the gate but T3 is blocked).

## T2 transient solid (src/solid/svk_dynamic.py)
Schemes implemented: "newmark" (beta = 1/4, gamma = 1/2, average acceleration; the oomph wall scheme) and "bdf2" = second derivative by the one-sided second-order backward formula (2 d^{n+1} - 5 d^n + 4 d^{n-1} - d^{n-2})/dt^2 (OpenFOAM `backward` d2dt2, the solids4foam tutorial solid scheme; first step uses (d^{n+1} - 2 d^n + d^{n-1})/dt^2 with d^{-1} = d^0 - dt v0 + dt^2 a0/2; exact s4f start-up not verified, U), velocity by BDF2 first derivative. Static mode (no inertia) shares the residual. Inertia rho_s a . v dX_0 (TL, reference density).
Tests (Q2, 160x4 unless stated; 2 MPI ranks; dt in units of the linear first-mode period T1 = 2 pi / omega1, omega1 = 22.373 sqrt(D/(rho h L^4)) clamped-clamped Euler-Bernoulli, ~72 rad/s):
- F1 free vibration, tiny amplitude (static shape at p = 0.01 Pa released at t=0), Newmark dt = T1/200: measured period (zero-crossing/phase fit of mid displacement over 5 periods) within 0.5% of the Euler-Bernoulli estimate; period error converges with dt (dt, dt/2, dt/4) with observed order in [1.5, 2.5] for Newmark.
- F2 energy: large amplitude (static shape at 200 Pa released, p = 0), t <= 1.5 s: Newmark total energy (kinetic + int psi_SVK dX_0) drift < 1% of E0 at dt = T1/100; BDF2 energy is reported (expected non-increasing; dissipation quantified).
- F3 small-amplitude agreement: constant p = 0.01 Pa suddenly applied: time mean of the midpoint displacement over an integer number of periods equals the static value within 1% (Newmark); slow ramp (10 T1) tracks the static value within 0.5%.
- F4 time convergence on the actual 0-0.25 s ramp + hold to t = 0.5 s at 200 Pa (follower), mid displacement at t = 0.5, dt = 0.0125, 0.00625, 0.003125, both schemes, analyse rules; reported with observed order; gate p >= 1.5 for Newmark only if the differences are monotone.
All jobs: Slurm via slurm/submit_run.sh, >= 2 ranks, converged flags checked from outputs (EXIT file AND the JSON `all_converged` flag), ledgered by job.sbatch. Budget <= 1.0 CPU-h total.

## T3
Conditional: build coupled solver only if G1/G2 pass AND a collapsible section exists in s4f_side/WALL671_RESULTS.md (not present at the time of this writing, checked 2026-11-06). Otherwise hand back without coding it.

## Amendment A1 (2026-11-06, written before the T2 dynamic runs; T1 job 16739 already submitted unchanged)
To stay in budget: F1 uses Q2 80x4 with dt = T1/100, T1/200, T1/400 and 3 periods (period from the mean spacing of 3 same-direction zero crossings of the mid displacement), plus one spatial check at 160x4, dt=T1/200; gate unchanged (Newmark period within 0.5% of the Euler-Bernoulli estimate at dt=T1/200, observed order in [1.5,2.5] on the three dts if the successive differences exceed 1e-7 relative). F2/F3/F4 use Q2 80x4 (F3 sudden load: mean over 3 measured periods; ramp: 10 T1). F4 uses four dts 0.0125 ... 0.0015625 (finest three analysed; all triples reported). Development smoke runs (jobs 16729 smoke, 16730s dbg) were run before this to debug the point-evaluation helper (a bug returning u_x instead of u_y was found and fixed); they carry no gate information.

## Amendment A2 (2026-11-06, after F2 first result, before F2b/F4 reruns)
F2 at dt = T1/100 gave Newmark max energy drift 1.97% (> 1% gate: FAIL as registered, kept as the registered result) and BDF2 77% (dissipation). Adding F2b (Newmark dt = T1/200, T1/400, diagnostic, no new gate: reported order of the drift decay) to test whether the drift is a time-step effect. Also: the F2F4 job (16741) executed only F2 (F4 absent in its output; cause being checked); F4 is re-run in its own job.

## Amendment A3 (2026-11-06, after F4 result)
F4 (job 16743; the earlier F2F4 job did not execute F4, cause not found, F4 rerun alone) at dt >= 0.0016 s shows ringing of the DRY wall (period ~0.09 s linear, shorter when stiffened by the 200 Pa membrane state; dt 0.0125..0.0016 = T1/7..T1/55): values at t=0.25/0.5 are non-monotone in dt, so the registered F4 gate cannot be evaluated (not in the asymptotic range); registered as NOT MET, not extrapolated. Added F4b (diagnostic): 0.5 s cosine ramp with no hold (smooth loading, little residual ringing), dt = 0.00625 ... 0.00078125, w at t=0.4 and 0.5, analysed with convergence.analyse, both schemes.

## Amendment A4 (2026-11-06, after F1 first result)
F1 first result (job 16740, cancelled after the 80x4 levels; values retained as PRE-FIX): period rel. error vs Euler-Bernoulli -0.228%, -0.200%, -0.174% at dt = T1/100, /200, /400, i.e. NOT converging at O(dt^2): cause identified: release from a deformed state started with a0 = 0 (inconsistent initial acceleration). Fix: `SVKDynamicSolid.consistent_acceleration` (solve M a0 = -R_int(d0) + load at t=0) now default in `set_initial`. F1, F2, F3 (F3 unaffected: rest start, but rerun because its job was cancelled) and F2b (dt=T1/200 only, budget) are rerun with the fix; the pre-fix F1/F2/F2b values are reported in COL1_RESULTS.md as superseded, F4/F4b (rest start, p0 = 0, a0 = 0 consistent) are unaffected. Gates unchanged.

## Amendment A5 (2026-11-06, after the fixed F1 rerun)
Fixed F1: period errors vs Euler-Bernoulli (80x4) -0.231%, -0.2007%, -0.1743% at dt = T1/100,/200,/400; 160x4 at T1/200: -0.138%. The EB-gate (0.5%) is MET. The registered order test on the 3-crossing period estimate is contaminated (successive differences 0.030%, 0.026% do not decay: the mid-point signal contains higher odd modes whose numerical phase shifts with dt move the zero crossings; not an O(dt^2) estimate) and is therefore NOT MET as registered. Replacement diagnostic F1c (registered here before it is run): solution-level time-convergence, max_t |w_dt - w_{T1/800}| / max|w| over 3 periods for dt = T1/50, /100, /200 (Newmark, 80x4), observed order expected in [1.5,2.5]. (F3 job failed on `np.trapz` removed from numpy 2: replaced by `np.trapezoid`, rerun.)
