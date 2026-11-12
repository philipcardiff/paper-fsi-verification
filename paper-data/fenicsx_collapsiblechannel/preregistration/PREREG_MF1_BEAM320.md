# PREREG_MF1_BEAM320: open-loop beam replay with 320 beam elements (written 2026-11-11, before the run)

Request: review of PR #71 asks for the MF1 open-loop replay with 320 beam elements, same load (S2, dt = 0.00625 s) and controls as before. Nothing is changed to improve agreement.

## 1. Existing setup (recorded; unchanged)
- Load history: `colP_data/data_MF1_S2_dt0.00625_loads.npz` (arrays x (321), gx, gy (240 x 321), pext (240), t (240), dt = 0.00625): the fluid residual vector at the 321 P2 interface nodes of the coupled S2 run (ny = 16) at every step; traction density per reference length t = M^-1 g with the undeformed P2 interface mass matrix (spacing 1/16); force on the wall = -t; external pressure as follower at the end-of-step time.
- Beam: `colP_beam.py` (copy of the COL-2 Hermite Kirchhoff-Love beam: E_eff = E/(1-nu^2) = 5.0e8 (E = 4.55e8, nu = 0.3), h = 0.05, rho = 1 (consistent translational mass rho h), clamped ends, complex-step Jacobian, Newton tolerance 1e-5 (assembly floor), 30 iterations), dead load interpolated from the P2 traction at the 5 Gauss points per element, time operator d2dt2 (u = BDF2 derivative of q, acceleration BDF2 of u, backward Euler on step 1; the coupled wall operator), observation: beam axis displacement Y at x = 7.5, 10, 12.5.
- Continuum reference (unchanged, not rerun): `replay_MF1.json` continuum replay of the same load on the wall-only submesh (160 x 4 P2 triangles), midplane (y = 1.025) and bottom (y = 1) displacements at the same points, controls C1 1.7e-14 and C2 1.0e-7 passed.
- Definitions: normalisation P = 0.216 m, percent of P; gap G = max_t |Y_beam - d_mid| / P; V1 = first trough on a cubic spline of the 241-point history (V123 of the production analysis); signed trough difference dV1 = V1(beam) - V1(continuum midplane) (negative: beam deeper); continuum minus beam = -dV1.
- Previous results (`data/mf1_gap.csv`, `mf1_resolution.csv`), x = 7.5 / 10 / 12.5: ne80 G 0.020 / 0.016 / 0.022 and dV1 -0.0205 / +0.0163 / -0.0222; ne160 G 0.053 / 0.033 / 0.056 and dV1 -0.0533 / -0.0329 / -0.0555 (continuum-minus-beam +0.053 / +0.033 / +0.056); ne160 vs ne80 max |dY|/P 0.033 / 0.049 / 0.033. AMEND2 resolution rule: |G(160) - G(80)| > 0.05 % P means unresolved; H6' (G <= 0.15 % P) passed.

## 2. The run
One Slurm job, 2 ranks (beam on rank 0), env fenicsx010 via `environment/run_in_env.sh`, HYDRA_BOOTSTRAP=fork, driver `scripts/colP_mf1_beam_refine.py` (the beam part of `colP_mf1_replay.py`, unchanged; element list as argument): `ne=320` only. A preceding small reproduction job (`ne=160,80`, about 10 s) checks the existing results from the saved load file and scripts. No broader refinement (no 640, no other dt, no other load).
Expected cost: 160 elements took 7.4 s and 80 elements 2.7 s wall; the dense complex-step Jacobian makes the cost grow faster than linearly, so 320 elements are expected at 1-5 min wall: <= 0.1 CPU-h (2 ranks; time limit 30 min). Stop and report if the projected cost exceeds 5 CPU-h. No duplicate run is active (squeue checked).

## 3. Acceptance tests (before looking at the 320 result)
A1 reproduction: the recomputed ne160 and ne80 Y histories equal the stored ones to <= 1e-9 (relative to max |Y|); otherwise report the discrepancy and do not interpret.
A2 run validity: status OK from the driver's own flags: all 240 steps done (reached_end), every Newton converged within 30 iterations (ok true), no launcher-only EXIT.
A3 data hygiene: same time grid as the continuum (t = 0, 0.00625, ..., 1.5), same observation points.
## 4. Reported quantities (no target, no extrapolation)
For ne = 80, 160, 320 at x = 7.5, 10, 12.5: beam first trough V1 and its change per refinement; signed dV1 against the continuum midplane and its change; G (max over t) and its change; beam-versus-beam differences max |dY|/P (160 vs 80, 320 vs 160); time of the beam and continuum V1 and the time difference; Newton counts; sign of dV1 for each ne (in particular x = 10, where the sign changed between 80 and 160).
Classification, fixed now: the gap is called *settled at the 0.05 level* only if |G(320) - G(160)| <= 0.05 at all three points AND the sign of dV1 is the same for ne160 and ne320 at all three points AND the two beam differences (160 vs 80, 320 vs 160) are not both larger than the gap itself; otherwise *not settled*. No convergence order, no limit and no model-form estimate is derived from three levels; the changes are reported as they are.
Source data, scripts and CSVs are shipped; the supplement text is updated either way.
