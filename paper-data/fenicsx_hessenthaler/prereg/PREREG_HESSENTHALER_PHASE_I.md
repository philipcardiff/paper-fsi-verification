# Preregistration: Hessenthaler Phase I, FEniCSx cross-code check (written BEFORE any production run)

Date: 2026-11-11 (UTC), coordinating agent. Worktree branch `hessenthaler-phaseI` off master 7eab623.
Status of this document: its sha256 is recorded in prereg/PREREG.sha256 and committed before any production
(coupled or fluid-only/solid-only) run. Thresholds are NOT changed afterwards; amendments are appended as
dated, separately hashed addenda and never replace a threshold.

## 0. Framing (binding)
* Primary: cross-code solution verification FEniCSx vs solids4foam (s4f). The FEniCSx work is a same-group,
  same-AI-lineage second implementation and cross-code check, NOT independent verification.
* Secondary, separate: comparison with the Hessenthaler et al. (2017) MRI data (figshare 10.6084/m9.figshare.4141836.v1),
  conditional on the zero-flow calibration. Not a new experimental validation. Agreement is not evidence of convergence.
* Corrected s4f reference (feF1-feF3, solids4foam 517a44921, OpenFOAM v2512, UCD Sonic, started 2026-11-07): status at
  survey time = IN PROGRESS, no corrected results published (PR #585 head b65d327f1, PR #666 open, issues #552/#665 open).
  Chamfered-mesh results are SUPERSEDED and are never mixed with the corrected family. Until corrected feF results exist
  every cross-code statement is PRELIMINARY and limited to model/definition/geometry/solid equivalence and to
  independently checkable analytic quantities.

## 1. Model (fixed; no tuning)
* Geometry: published fluid surface (tutorial fluidDomain.stl.gz sha256 0be13c04..., metres; inlets z=-29.5 mm,
  outlet z=339.5 mm = 250 mm extension), flap 11 x 2 x 65 mm, base clamped at z=0, flow +z, buoyancy +y.
  Fluid and solid meshes share one interface surface (conforming, monolithic single mesh).
* Fluid: incompressible Newtonian, rho_f=1163.3, nu=1.074529e-5 (mu_f=12.50 mPa s), laminar.
  Inlets: parabolic, peaks 630 mm/s (upper, y>0) and 615 mm/s (lower), ramp 3s^2-2s^3, s=min(t/0.5,1).
  Outlet: p=0 with the backflow-stabilised traction sigma n = beta rho min(u.n,0) u, beta=0.2 (Esmaily-Moghadam), the
  benchmark definition. No-slip walls. Fluid solved without gravity.
* Solid: total-Lagrangian compressible neo-Hookean, W = mu/2 (Ibar1-3) + K/4 (J^2-1-2 ln J), mu=64.2 kPa, nu=0.45
  => K = 2 mu (1+nu)/(3(1-2nu)) = 620.6 kPa, rho_s=1058.3. Verified equal to s4f neoHookeanElastic law at 517a44921
  (sigma = mu/J dev(bbar) + K/2 (J^2-1)/J I; line 213-222 of the law source). Body force rho_s*g_eff with
  g_eff=(0,0.973306,0) m/s^2 ramped by the same ramp. The tutorial value 58.45 kPa is the s4f S2 calibration and is NOT used.
* mu is never changed to improve any agreement.

## 2. Acceptance criteria (gates, in order; a failed gate blocks later simulations until diagnosed)
GATE A (geometry/mesh, no solver): flap interface area within 0.1 % of 1712 mm^2; flap thickness 2.000 mm, width 11.000 mm,
 length 65.000 mm, tip z=65.000 mm within 1e-6 m; the same interface facets belong to the fluid and the solid cell sets
 (shared nodes, zero gap); outward solid normals consistent (sum of oriented areas of closed solid < 1e-12 m^2); no
 chamfer (interface facets coplanar with the 5 flap faces within 1e-9 m); inlet areas within 0.5 % of pi 10.95^2 mm^2
 and outlet within 0.5 % of pi 38.1^2 mm^2 (polygonal STL: gap reported, not hidden); all cells positive volume,
 minimum scaled Jacobian (gmsh SICN) > 0.05 for the production mesh.
GATE B (solid only): (B1) analytic P vs UFL dW/dF rel. diff < 1e-10; consistent tangent vs central FD < 1e-6 rel.;
 frame indifference of W < 1e-12 and sigma(QF)=Q sigma Q^T < 1e-10; pure-rotation stress = 0 < 1e-9 mu.
 (B2) zero-flow static flap: Newton converges (|r|/|r0|<1e-9) with J>0 everywhere; asymptotic Newton convergence
 (final step reduces residual by >= 1e2); three nested solid levels. PREDICTION (registered now): continuum-limit tip y
 at mu=64.2 kPa lies in [29.20, 29.60] mm (s4f continuum calibration mu_cal=63.95 kPa for 29.50 mm gives
 about 29.38 mm). Acceptance: finest-level tip within 0.5 mm of 29.50 mm AND the 3-level sequence monotone. A miss is
 reported as a failed prediction, mu is NOT retuned.
GATE C (fluid only, flap rigid and undeformed, full ramp value): |Q_in,computed - Q_in,exact|/Q_in,exact <= 2 % per inlet
 (polygon-inscribed STL inlet vs analytic 21.9 mm circle; the exact value is computed from the mesh inlet geometry and the
 discrepancy to the nominal is reported separately); |Q_out + backflow - sum Q_in|/sum Q_in <= 1e-3; Newton converges to
 rel 1e-8 in <= 12 iterations from the ramped state; no reversed flow of net inflow at the outlet face
 (backflow fraction of outlet flux < 1 %); pressure/traction conventions documented.
GATE D (first coupled pilot): per step Newton converged (rel 1e-8, <= 12 its); min(J_ale) > 0 and min(detF_ale/detF_ale0)
 >= 0.2 (ALE quality) and J_solid > 0.3; kinematic condition u=dd/dt at the interface satisfied to < 1e-6 relative of
 the interface speed (checked as post-processing, not by construction alone); action-reaction: |F_fluid_on_solid +
 [solid inertia + internal + body + clamp reaction]| <= max(2 % of |F_fluid|, 5 uN); M1 (exact Nanson) interface measure;
 mass: |Q_out - sum Q_in| / sum Q_in <= 1e-3 (net flap volume change is O(1e-6) of Q dt); tip y finite and within
 [0, 30] mm and monotone-ramp behaviour: no sign inversion. Stop-rules as in the mandate.
GATE E (scaling, once): same restart state, same 10 steps at 16/32/64 cores if memory/CPU budget allows; report only.
SETTLED-STATE test (any claim "steady"): tip y drift |dy/dt| < 5e-4 mm/s averaged over a window of >= 4 s, flux and
 force drifts below 0.1 %/s, three consecutive windows. Otherwise the result is labelled "transient t = ..." only.
CONVERGENCE claims: only with >= 3 valid monotone levels, observed order p in [0.5, 3]; fluid, solid and time refinement
 kept separate; never a continuum limit from fewer levels.
CROSS-CODE (Campaign F): defined comparison metrics = tip y (steady/at t), deformed centreline at x=0 vs z, net
 fluid force on flap (y,z), inlet/outlet flux, pressure on the flap mid-line; compared only to corrected-family s4f
 data at the same level/time; difference classes: matched-mesh / matched-time / extrapolated-limit / model-choice.
 No numerical tolerance for "agreement" is set before s4f data exist; a small discrepancy is not validation.
MRI (Campaign H): no tolerance (no experimental uncertainty available; voxel 1.302x1.302x6 mm is not an uncertainty).
 Metrics: tip y vs 16.41 mm, centreline RMS and max deviation, velocity component RMS/L2 over the 278 points. Only reported.
 Not used for mesh/dt selection nor mu.

## 3. Resource caps
<= 64 concurrent cores; <= ~400 allocated CPU-h total; wall target <= 9 h per job; >= 2 cores per substantive job; Slurm
only; HYDRA_BOOTSTRAP=fork (run_in_env.sh); never use cores allocated to others; check squeue/sinfo before each submission;
EXIT OK of the launcher is never trusted: check converged flags, reached_end, logs. Ledger: ledger/HESSENTHALER_LEDGER.csv.

## 4. Run matrix (planned; cost estimates, will be refined from measured step costs)
See docs/PHASE0_PLAN.md.
