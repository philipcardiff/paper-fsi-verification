# Hessenthaler Phase I: exact model and numerical method of the FEniCSx implementation

All quantities SI. Reference configuration Omega_0 = fluid cells (tag 1) U solid cells (tag 2) of ONE conforming tetrahedral mesh.
Geometry: published fluid surface (s4f tutorial fluidDomain.stl.gz sha256 0be13c049ddf84eb97c7d48af0b6c7680455e9b03342846f0373ad5be3ed815f;
inlets z = -29.5 mm, outlet z = 339.5 mm); flap 11 x 2 x 65 mm (x in [-5.5, 5.5], y in [-1, 1], z in [0, 65] mm), clamped at z = 0, flow +z, buoyancy +y.
The wall patch (122 610 triangles, highly anisotropic) is quadric-decimated to 1 838 triangles (max deviation from the published surface 0.12 mm,
area +0.004 %); inlet/outlet discs, the flap interface (exact 5 faces) and the clamp rectangle are exact; remeshed with gmsh (Delaunay/Frontal, size fields), micro-edges (< 0.15 mm)
collapsed. The same interface triangles belong to a fluid and a solid tet (Gate A).

Unknowns: u (P2 vector, fluid velocity / solid velocity, continuous), p (P1, fluid; pinned to 0 strictly inside the solid), d (P2 vector, solid displacement / harmonic mesh displacement).
Time discretisation: BDF2 (first step backward Euler), constant dt, all terms implicit.

Fluid (ALE in the reference configuration, F = I + grad d, J = det F, grad_x f = grad f F^-1, w = d_t d):
  rho_f [ d_t(J u) - u d_t J ] . phi  (discrete GCL, conservative form)  + rho_f J (grad_x u)(u - w) . phi + 2 mu_f J sym(grad_x u) : sym(grad_x phi)
  - p J div_x phi + q J div_x u = 0,   rho_f = 1163.3, nu_f = 1.074529e-5 (mu_f = 12.50 mPa s). No gravity (hydrostatic part absorbed in p).
  Inlets: u = (0, 0, V_peak (1 - r^2/R^2) s(t)), R = 10.95 mm centred at y = +-27.15 mm, V_peak = 0.630 (upper) / 0.615 (lower) m/s; s(t) = 3 s^2 - 2 s^3, s = min(t/0.5, 1).
  Walls: u = 0. Outlet: p = 0 natural, viscous gradient-form correction (-mu_f (grad u)^T n . phi), backflow-stabilised traction t = sigma n = +beta rho_f min(u.n, 0) u (opposes u at backflow, dissipative; an earlier draft/first runs had the opposite sign, corrected after the independent math review), beta = 0.2
  (the |u.n| is smoothed with sqrt(x^2 + 1e-12)).
Mesh motion: int k grad d : grad psi dX = 0 on the fluid cells, k = mean cell volume / cell volume (small cells stiffer); d = 0 on inlets, outlet, wall, clamp; the rows of the interface d dofs are replaced by the kinematic condition.
Solid (total Lagrangian, compressible neo-Hookean; verified identical to solids4foam neoHookeanElastic at 517a44921):
  W = mu/2 (J^(-2/3) tr(F F^T) - 3) + K/4 (J^2 - 1 - 2 ln J),  P = dW/dF = mu J^(-2/3) (F - tr(F F^T)/3 F^-T) + K/2 (J^2 - 1) F^-T,
  mu = 64.2 kPa (fixed), nu = 0.45, K = 2 mu (1 + nu) / (3 (1 - 2 nu)) = 620.6 kPa, rho_s = 1058.3.
  rho_s d_t u . phi + P : grad phi - s(t) rho_s g_eff e_y . phi = 0,  g_eff = (rho_s - rho_f)/rho_s * 9.81 = 0.973306 m/s^2 magnitude, acting in +y (net buoyancy),
  (d_t d - u) . psi mu_s = 0 on the solid cells (kinematic condition; the factor mu_s only scales the rows).
Interface: u, d continuous (single mesh) => kinematic and no-slip conditions hold exactly in the discrete space; the traction balance is the sum of the fluid and solid momentum residuals
  (natural condition): the force on the solid is int J sigma_f F^-T N dA0 (exact Nanson, "M1").
Nonlinear solver: PETSc SNES Newton (backtracking), UFL-derived exact Jacobian, rtol 1e-8 / atol 1e-12. Linear solver: GMRES (rtol 1e-6) preconditioned by a parallel MUMPS LU factorisation of the
  Jacobian that is REUSED (SNES lag, persistent) and refactorised when the Krylov iterations per Newton step exceed a threshold (kthr). MUMPS: ICNTL(14)=80, parallel analysis (ICNTL(28)=2, ICNTL(29)=2).
  A global LU of a 1.1M-dof P2/P1/P2 3D problem needs ~157 GB (measured, 32 ranks, fluid-only) and ~7 min for the first factorisation; it is therefore used only as a lagged preconditioner.
Environment: dolfinx 0.10.0, PETSc 3.25 (MUMPS), MPICH 5 (OFI/shm, HYDRA_BOOTSTRAP=fork), gmsh 4.15.2, numpy 2.4.6, Slurm, node xenosim.
