"""3-D monolithic ALE FSI for the Hessenthaler Phase I flap on ONE conforming mesh (fluid + solid tets).
Unknowns (u, p, d): u P2 vector (fluid velocity / solid velocity, continuous -> no-slip & kinematic condition),
p P1 (fluid pressure; pinned to 0 strictly inside the solid), d P2 vector (solid displacement / harmonic mesh displacement).
Fluid cells: ALE Navier-Stokes in the reference configuration, BDF2 with discrete-GCL conservative form (src/ale/forms.py).
Solid cells: total-Lagrangian neo-Hookean (nh3d), rho_s u_t . phi + P(F(d)):grad phi - rho_s g_eff(t) e_y . phi = 0, d_t = u.
The interface traction is the natural sum of the two momentum residuals: exact Nanson (M1) J sigma_f F^-T N dA0.
Outlet: p=0 natural + backflow-stabilised traction sigma n = +beta rho min(u.n,0) u (beta = 0.2; sign corrected after math review), 'gradient' viscous form.
Ramp s(t) = 3s^2-2s^3, s = min(t/0.5,1) multiplies the inlet velocity and the buoyancy.
mode 'coupled' or 'rigid' (solid and interface held fixed: fluid-only gate C)."""
from __future__ import annotations
import numpy as np
import ufl
import basix.ufl
from dolfinx import fem
from petsc4py import PETSc
from mpi4py import MPI
from ..solvers.snes import SNESSolver
from ..ale.forms import fluid_ale_residual_bdf2, harmonic_residual, kin
from ..ale.mesh_motion import stiffness_weight
from . import nh3d

IDS = dict(fluid=1, solid=2, upper=1, lower=2, outlet=3, wall=4, interface=5, clamp=6)
RHO_F, NU_F = 1163.3, 1.074529e-5
MU_F = RHO_F * NU_F
RHO_S, G_EFF = 1058.3, 0.973306
YC = {"upper": +27.15e-3, "lower": -27.15e-3}
RIN = 10.95e-3
VPEAK = {"upper": 0.630, "lower": 0.615}
BETA = 0.2

def ramp(t, t_ramp=0.5):
    s = min(max(t / t_ramp, 0.0), 1.0)
    return 3 * s ** 2 - 2 * s ** 3


class HessenthalerFSI:
    def __init__(self, msh, cell_tags, facet_tags, mu_s=nh3d.MU, nu_s=nh3d.NU, mode="coupled", prefix="fsi_", snes_options=None, quad_degree=4, mesh_exp=1.0, mmotion="total"):
        import time as _t; _t0 = _t.time()
        self.msh, self.mode, self.mu_s, self.nu_s = msh, mode, mu_s, nu_s
        self.K_s = nh3d.bulk(mu_s, nu_s)
        gdim, tdim = 3, 3
        msh.topology.create_connectivity(2, 3); msh.topology.create_connectivity(3, 3)
        cell = msh.basix_cell()
        Pu = basix.ufl.element("Lagrange", cell, 2, shape=(3,)); Pp = basix.ufl.element("Lagrange", cell, 1); Pd = basix.ufl.element("Lagrange", cell, 2, shape=(3,))
        self.W = W = fem.functionspace(msh, basix.ufl.mixed_element([Pu, Pp, Pd]))
        self.Vu, self.mu_map = W.sub(0).collapse(); self.Vp, self.mp_map = W.sub(1).collapse(); self.Vd, self.md_map = W.sub(2).collapse()
        self.w, self.w0, self.wm = (fem.Function(W, name=n) for n in ("w", "w0", "wm"))
        self.bdf = [fem.Constant(msh, PETSc.ScalarType(c)) for c in (1.0, -1.0, 0.0)]
        self.dt = fem.Constant(msh, PETSc.ScalarType(1.0)); self.t = 0.0; self.nsteps = 0; self.dt_prev = None
        self.gfac = fem.Constant(msh, PETSc.ScalarType(0.0))        # buoyancy ramp factor
        md = {"quadrature_degree": quad_degree}
        self.cell_tags, self.facet_tags = cell_tags, facet_tags
        self.dx = ufl.Measure("dx", domain=msh, subdomain_data=cell_tags, metadata=md)
        self.ds = ufl.Measure("ds", domain=msh, subdomain_data=facet_tags)
        dxf, dxs = self.dx(IDS["fluid"]), self.dx(IDS["solid"]); self.dxf, self.dxs = dxf, dxs
        u1, p1, d1 = ufl.split(self.w); u0, _, d0 = ufl.split(self.w0); um, _, dm = ufl.split(self.wm)
        phi, q, psi = ufl.TestFunctions(W)
        a0, a1, a2 = self.bdf; dt = self.dt
        if mode == "rigid":     # d == 0 identically: plain NS on the fixed mesh (no ALE terms), d/solid rows are identity equations
            z3 = fem.Constant(msh, np.zeros(3))
            F_f = fluid_ale_residual_bdf2(u1, p1, z3, u0, z3, um, z3, dt, a0, a1, a2, RHO_F, MU_F, phi, q, dxf)
        else:
            F_f = fluid_ale_residual_bdf2(u1, p1, d1, u0, d0, um, dm, dt, a0, a1, a2, RHO_F, MU_F, phi, q, dxf)
        nn = ufl.FacetNormal(msh)
        # outlet: natural p=0, gradient-form viscous correction, backflow-stabilised traction
        F_f += -MU_F * ufl.dot(ufl.dot(ufl.grad(u1).T, nn), phi) * self.ds(IDS["outlet"])
        un = ufl.dot(u1, nn); F_f += -BETA * RHO_F * 0.5 * (un - ufl.sqrt(un ** 2 + 1e-12)) * ufl.dot(u1, phi) * self.ds(IDS["outlet"])   # traction t = +beta rho min(u.n,0) u (dissipative; sign fixed after math review, see analysis/math_review)
        P = nh3d.P_analytic(ufl.Identity(3) + ufl.grad(d1), mu_s, self.K_s)
        body = ufl.as_vector((0.0, RHO_S * G_EFF, 0.0))
        F_s = (RHO_S * ufl.dot((a0 * u1 + a1 * u0 + a2 * um) / dt, phi) * dxs + ufl.inner(P, ufl.grad(phi)) * dxs
               - self.gfac * ufl.dot(body, phi) * dxs
               + ufl.dot((a0 * d1 + a1 * d0 + a2 * dm) / dt - u1, psi) * mu_s * dxs
               + 1e-30 * p1 * q * dxs)       # sparsity/diagonal for the pinned solid-interior pressure rows (value ~1e-30, no effect)
        self.kmesh = stiffness_weight(msh, "invvol")
        if mesh_exp != 1.0: self.kmesh.x.array[:] = self.kmesh.x.array ** mesh_exp      # k = (mean cell volume / cell volume)^mesh_exp
        if mmotion == "incr":     # Addendum A2: harmonic extension of the increment on the previous configuration
            F0 = ufl.Identity(3) + ufl.grad(d0); F0i = ufl.inv(F0); dd = d1 - d0
            F_m = self.kmesh * ufl.inner(ufl.dot(ufl.grad(dd), F0i), ufl.dot(ufl.grad(psi), F0i)) * dxf
        else:
            F_m = harmonic_residual(d1, psi, dxf, k=self.kmesh)
        self.F_fluid, self.F_solid, self.F_mesh = F_f, F_s, F_m
        # ---- boundary conditions
        fdim = 2; self.g_in = fem.Function(self.Vu); bcs = []
        def dofs(sub, V, tags):
            f = np.unique(np.concatenate([facet_tags.find(t) for t in tags])).astype(np.int32)
            return fem.locate_dofs_topological((W.sub(sub), V), fdim, f)
        zero_u, zero_d = fem.Function(self.Vu), fem.Function(self.Vd)
        bcs.append(fem.dirichletbc(self.g_in, dofs(0, self.Vu, [IDS["upper"], IDS["lower"]]), W.sub(0)))
        bcs.append(fem.dirichletbc(zero_u, dofs(0, self.Vu, [IDS["wall"], IDS["clamp"]]), W.sub(0)))
        bcs.append(fem.dirichletbc(zero_d, dofs(2, self.Vd, [IDS["upper"], IDS["lower"], IDS["outlet"], IDS["wall"], IDS["clamp"]]), W.sub(2)))
        fcells, scells = cell_tags.find(IDS["fluid"]), cell_tags.find(IDS["solid"])
        pf = fem.locate_dofs_topological((W.sub(1), self.Vp), tdim, fcells)[1]
        ps_par, ps_col = fem.locate_dofs_topological((W.sub(1), self.Vp), tdim, scells)
        keep = ~np.isin(ps_col, pf)
        zero_p = fem.Function(self.Vp)          # collective: created on every rank
        bcs.append(fem.dirichletbc(zero_p, [ps_par[keep], ps_col[keep]], W.sub(1)))
        if mode == "rigid":
            # solid and interface held fixed: u = 0, d = 0 on all solid cells and interface; mesh d = 0 everywhere
            so_u = fem.locate_dofs_topological((W.sub(0), self.Vu), tdim, scells)
            so_d = fem.locate_dofs_topological((W.sub(2), self.Vd), tdim, scells)
            fl_d = fem.locate_dofs_topological((W.sub(2), self.Vd), tdim, fcells)
            bcs.append(fem.dirichletbc(zero_u, so_u, W.sub(0))); bcs.append(fem.dirichletbc(zero_d, so_d, W.sub(2))); bcs.append(fem.dirichletbc(zero_d, fl_d, W.sub(2)))
        self.bcs = bcs
        if msh.comm.rank == 0: print(f'[setup] forms/BCs built {_t.time()-_t0:.1f}s', flush=True)
        idof = dofs(2, self.Vd, [IDS["interface"]])[0]
        if mode == "rigid":
            Fmain = F_f + ufl.dot(u1, phi) * dxs + ufl.dot(d1, psi) * self.dx + 1e-30 * p1 * q * dxs      # identity rows for the pinned dofs (sparsity/diagonal)
        else:
            Fmain = F_f + F_s
        self.snes = SNESSolver(Fmain, self.w, bcs, prefix=prefix, F_masked=(None if mode == "rigid" else F_m),
                               masked_dofs=(None if mode == "rigid" else idof), options=snes_options)
        if msh.comm.rank == 0: print(f'[setup] SNES/forms compiled {_t.time()-_t0:.1f}s', flush=True)
        # fixed coordinates of inlet dofs for the profile
        self._inlet_dofs = {}
        for name in ("upper", "lower"):
            f = facet_tags.find(IDS[name]); self._inlet_dofs[name] = fem.locate_dofs_topological(self.Vu, fdim, f)
        self._Xu = self.Vu.tabulate_dof_coordinates()

    def set_time(self, t):
        self.t = t
        r = ramp(t); self.gfac.value = r
        g = self.g_in.x.array.reshape(-1, 3); g[:] = 0.0
        for name, dofs in self._inlet_dofs.items():
            x = self._Xu[dofs]
            rr2 = (x[:, 0] ** 2 + (x[:, 1] - YC[name]) ** 2) / RIN ** 2
            g[dofs, 2] = VPEAK[name] * np.maximum(1.0 - rr2, 0.0) * r
        self.g_in.x.scatter_forward()

    def step(self, dt, predictor=False):
        self.dt.value = dt
        t_new = self.t + dt
        self.set_time(t_new)
        if self.nsteps >= 1:      # variable-step BDF2 (w = dt_n/dt_{n-1}); equals (3/2,-2,1/2) for constant dt
            w_ = dt / self.dt_prev; c = ((1 + 2 * w_) / (1 + w_), -(1 + w_), w_ ** 2 / (1 + w_))
        else:
            c = (1.0, -1.0, 0.0)
        for k, v in zip(self.bdf, c): k.value = v
        if predictor and self.nsteps >= 1:      # Addendum A3b: linear extrapolation of the initial guess
            om = dt / self.dt_prev
            self.w.x.array[:] = self.w0.x.array + om * (self.w0.x.array - self.wm.x.array); self.w.x.scatter_forward()
        res = self.snes.solve()
        return res

    def accept(self):
        self.wm.x.array[:] = self.w0.x.array; self.w0.x.array[:] = self.w.x.array; self.nsteps += 1; self.dt_prev = float(self.dt.value)

    def revert(self, t_old):
        self.w.x.array[:] = self.w0.x.array; self.w.x.scatter_forward(); self.t = t_old

    # ---------------- diagnostics
    def sub(self, i): return self.w.sub(i).collapse()
    def flux(self, tag):
        u = ufl.split(self.w)[0]; n = ufl.FacetNormal(self.msh)
        return self.msh.comm.allreduce(fem.assemble_scalar(fem.form(ufl.dot(u, n) * self.ds(tag))), op=MPI.SUM)
    def backflow(self, tag):
        u = ufl.split(self.w)[0]; n = ufl.FacetNormal(self.msh); un = ufl.dot(u, n)
        return self.msh.comm.allreduce(fem.assemble_scalar(fem.form(0.5 * (un - abs(un)) * self.ds(tag))), op=MPI.SUM)
    def detF_stats(self):
        d = ufl.split(self.w)[2]
        Q = fem.functionspace(self.msh, ("DG", 0)); f = fem.Function(Q)
        out = []
        for cell_id, nm in ((IDS["fluid"], "fluid"), (IDS["solid"], "solid")):
            f.interpolate(fem.Expression(ufl.det(ufl.Identity(3) + ufl.grad(d)), Q.element.interpolation_points))
            idx = self.cell_tags.find(cell_id)
            vals = f.x.array[idx] if len(idx) else np.array([np.inf])
            out.append((self.msh.comm.allreduce(float(vals.min()), op=MPI.MIN), self.msh.comm.allreduce(float(vals.max() if len(idx) else -np.inf), op=MPI.MAX)))
        return dict(fluid=out[0], solid=out[1])
    def point_value(self, field, pts):
        from dolfinx import geometry
        tree = geometry.bb_tree(self.msh, 3); P = np.asarray(pts, float).reshape(-1, 3)
        cand = geometry.compute_collisions_points(tree, P); col = geometry.compute_colliding_cells(self.msh, cand, P)
        out = np.full((len(P), 3), np.nan)
        for i in range(len(P)):
            l = col.links(i)
            if len(l): out[i] = np.asarray(field.eval(P[i:i + 1], l[:1])).reshape(-1)[:3]
        return out
    def tip(self):
        d = self.w.sub(2).collapse()
        v = self.point_value(d, [[0.0, 0.0, 65e-3 - 1e-9]])
        g = self.msh.comm.allgather(v[0])
        g = [x for x in g if not np.isnan(x[0])]
        return np.array(g[0]) if g else np.full(3, np.nan)
    def centreline(self, n=66):
        d = self.w.sub(2).collapse(); z = np.linspace(1e-9, 65e-3 - 1e-9, n)
        v = self.point_value(d, np.c_[np.zeros(n), np.zeros(n), z])
        allv = self.msh.comm.allgather(v); out = np.full((n, 3), np.nan)
        for a in allv:
            m = ~np.isnan(a[:, 0]); out[m] = a[m]
        return z, out
    def _arg(self, F): return [a for a in ufl.algorithms.extract_arguments(F) if a.number() == 0][0]
    def _lift(self, tag, comp):
        Vsc, _ = self.W.sub(0).sub(comp).collapse()
        d = fem.locate_dofs_topological((self.W.sub(0).sub(comp), Vsc), 2, self.facet_tags.find(tag))[0]
        lift = fem.Function(self.W); lift.x.array[np.unique(d)] = 1.0; return lift
    def interface_force_weak(self):
        """Force on the solid from the fluid momentum residual at interface u-dofs, and the solid-side residual at the same dofs."""
        of, os_ = np.zeros(3), np.zeros(3)
        for i in range(3):
            lift = self._lift(IDS["interface"], i)
            rf = fem.assemble_scalar(fem.form(ufl.replace(self.F_fluid, {self._arg(self.F_fluid): lift})))
            rs = fem.assemble_scalar(fem.form(ufl.replace(self.F_solid, {self._arg(self.F_solid): lift})))
            of[i] = self.msh.comm.allreduce(rf, op=MPI.SUM); os_[i] = self.msh.comm.allreduce(rs, op=MPI.SUM)
        return -of, os_
    def clamp_reaction(self):
        out = np.zeros(3)
        for i in range(3):
            lift = self._lift(IDS["clamp"], i)
            r = fem.assemble_scalar(fem.form(ufl.replace(self.F_solid, {self._arg(self.F_solid): lift})))
            out[i] = self.msh.comm.allreduce(r, op=MPI.SUM)
        return out
    def interface_force_direct(self):
        """int_Gamma0 J sigma_f F^-T N_s dA0 with N_s the solid-outward reference normal, using fluid-side quantities (interior facets)."""
        Q0 = fem.functionspace(self.msh, ("DG", 0)); chi = fem.Function(Q0)
        chi.x.array[self.cell_tags.indices[self.cell_tags.values == IDS["solid"]]] = 1.0; chi.x.scatter_forward()
        u, p, d = ufl.split(self.w); n = ufl.FacetNormal(self.msh); ns = n("+") * (2 * chi("+") - 1)
        F, J = kin(d); cfp, cfm = 1 - chi("+"), 1 - chi("-")
        gu = ufl.dot(ufl.grad(u), ufl.inv(F)); gu_f = cfp * gu("+") + cfm * gu("-")
        F_f = cfp * F("+") + cfm * F("-"); J_f = cfp * J("+") + cfm * J("-")
        sig = 2 * MU_F * ufl.sym(gu_f) - p("+") * ufl.Identity(3)
        dS = ufl.Measure("dS", domain=self.msh, subdomain_data=self.facet_tags)
        t = J_f * ufl.dot(sig, ufl.dot(ufl.inv(F_f).T, ns))
        pr = -p("+") * J_f * ufl.dot(ufl.inv(F_f).T, ns)
        out, outp = np.zeros(3), np.zeros(3)
        for i in range(3):
            out[i] = self.msh.comm.allreduce(fem.assemble_scalar(fem.form(t[i] * dS(IDS["interface"]))), op=MPI.SUM)
            outp[i] = self.msh.comm.allreduce(fem.assemble_scalar(fem.form(pr[i] * dS(IDS["interface"]))), op=MPI.SUM)
        return out, outp

    def balance(self):
        """Registered action-reaction ingredients for the solid at the just-solved step (call BEFORE accept(), so w0/wm are the n, n-1 states):
        I = int rho_s a dV (exact BDF2 acceleration), f_b = ramp rho_s g_eff V e_y, kinematic residual (L2, relative) of (BDF2 d_t d - u) on the solid cells."""
        u1, p1, d1 = ufl.split(self.w); u0, _, d0 = ufl.split(self.w0); um, _, dm = ufl.split(self.wm)
        a0, a1, a2 = self.bdf; dt = self.dt
        acc = (a0 * u1 + a1 * u0 + a2 * um) / dt
        comm = self.msh.comm
        I = np.array([comm.allreduce(fem.assemble_scalar(fem.form(RHO_S * acc[i] * self.dxs)), op=MPI.SUM) for i in range(3)])
        V = comm.allreduce(fem.assemble_scalar(fem.form(1.0 * self.dxs)), op=MPI.SUM)
        fb = np.array([0.0, float(self.gfac.value) * RHO_S * G_EFF * V, 0.0])
        kin_ = (a0 * d1 + a1 * d0 + a2 * dm) / dt - u1
        num = comm.allreduce(fem.assemble_scalar(fem.form(ufl.dot(kin_, kin_) * self.dxs)), op=MPI.SUM)
        den = comm.allreduce(fem.assemble_scalar(fem.form(ufl.dot(u1, u1) * self.dxs)), op=MPI.SUM)
        return I, fb, float(np.sqrt(num / max(den, 1e-300))), float(np.sqrt(den / V))

    # ------------------------------------------------------------------ Phase 2: independent force functionals on a fixed state
    def _lift_cells(self, comp, cells):
        Vsc, _ = self.W.sub(0).sub(comp).collapse()
        d = fem.locate_dofs_topological((self.W.sub(0).sub(comp), Vsc), 3, cells)[0]
        lift = fem.Function(self.W); lift.x.array[np.unique(d)] = 1.0; return lift

    def force_decomp(self):
        """Call BEFORE accept(). Reconstructs, component by component (x,y,z), the fluid and structural interface functionals of the ASSEMBLED equations:
        W_all = -F_f(lift_iface) (nodal consistent force on the solid, all interface nodes incl. those shared with the clamp), split into time/convection, viscous, pressure terms;
        W_free = same without the clamp-constrained nodes; F_dir = facet traction integral (cell-wise gradients); R_BC = (F_f+F_s)(lift_clamp) reaction rows;
        Fs_all = F_s(lift_solid_cells) (must equal I - f_b); I = int rho_s a; f_b buoyancy. Identity: I - f_b = R_BC + W_all."""
        comm = self.msh.comm
        u1, p1, d1 = ufl.split(self.w); u0, _, d0 = ufl.split(self.w0); um, _, dm = ufl.split(self.wm)
        a0, a1, a2 = self.bdf; dt = self.dt
        F1, J1 = kin(d1); _, J0 = kin(d0); _, Jm = kin(dm)
        w = (a0 * d1 + a1 * d0 + a2 * dm) / dt
        dJu = (a0 * J1 * u1 + a1 * J0 * u0 + a2 * Jm * um) / dt; dJ = (a0 * J1 + a1 * J0 + a2 * Jm) / dt
        gu = ufl.dot(ufl.grad(u1), ufl.inv(F1))
        out = {}
        sc = lambda form: comm.allreduce(fem.assemble_scalar(fem.form(form)), op=MPI.SUM)
        ifc = self.facet_tags.find(IDS["interface"]); clc = self.facet_tags.find(IDS["clamp"]); scells = self.cell_tags.find(IDS["solid"])
        terms = {k: np.zeros(3) for k in ("time", "conv", "visc", "pres", "outlet")}; W_free = np.zeros(3)
        Fs_all = np.zeros(3); Ff_cl = np.zeros(3); Fs_cl = np.zeros(3); Ff_all = np.zeros(3)
        for i in range(3):
            li = self._lift(IDS["interface"], i); lc = self._lift(IDS["clamp"], i); ls = self._lift_cells(i, scells)
            li_free = fem.Function(self.W); li_free.x.array[:] = li.x.array; li_free.x.array[lc.x.array > 0.5] = 0.0
            phi = ufl.split(li)[0]; gp = ufl.dot(ufl.grad(phi), ufl.inv(F1))
            terms["time"][i] = sc(RHO_F * ufl.dot(dJu - u1 * dJ, phi) * self.dxf)
            terms["conv"][i] = sc(RHO_F * J1 * ufl.dot(ufl.dot(gu, u1 - w), phi) * self.dxf)
            terms["visc"][i] = sc(2 * MU_F * J1 * ufl.inner(ufl.sym(gu), ufl.sym(gp)) * self.dxf)
            terms["pres"][i] = sc(-p1 * J1 * ufl.tr(gp) * self.dxf)
            Ff_all[i] = sc(ufl.replace(self.F_fluid, {self._arg(self.F_fluid): li}))
            W_free[i] = -sc(ufl.replace(self.F_fluid, {self._arg(self.F_fluid): li_free}))
            Ff_cl[i] = sc(ufl.replace(self.F_fluid, {self._arg(self.F_fluid): lc}))
            Fs_cl[i] = sc(ufl.replace(self.F_solid, {self._arg(self.F_solid): lc}))
            Fs_all[i] = sc(ufl.replace(self.F_solid, {self._arg(self.F_solid): ls}))
        I, fb, kin_rel, urms = self.balance()
        Fd, Fp = self.interface_force_direct()
        out.update(F_dir=Fd, F_dir_pressure=Fp, F_dir_viscous=Fd - Fp, W_all=-Ff_all, W_free=W_free, Ff_clamp=Ff_cl, Fs_clamp=Fs_cl, R_BC=Ff_cl + Fs_cl, Fs_all=Fs_all, I=I, f_b=fb,
                   weak_time=-terms["time"], weak_conv=-terms["conv"], weak_visc=-terms["visc"], weak_pres=-terms["pres"])
        out["identity_residual"] = (I - fb) - (out["R_BC"] + out["W_all"])          # must vanish (assembled equations)
        out["Fs_all_minus_(I-fb)"] = Fs_all - (I - fb)                                  # must vanish (assembly consistency)
        out["W_all_minus_F_dir"] = out["W_all"] - Fd
        out["terms_sum_minus_W_all"] = (-terms["time"] - terms["conv"] - terms["visc"] - terms["pres"]) - out["W_all"]
        return {k: (v.tolist() if hasattr(v, "tolist") else v) for k, v in out.items()}

    # ------------------------------------------------------------------ solver-architecture PoC: physics-informed block index sets
    def block_index_sets(self):
        """Global-index ISs of the three physical blocks of the monolithic Jacobian (owned dofs):
        S = solid + interface (u, d on solid cells incl. interface nodes), A = ALE interior (d on fluid-only nodes), F = fluid (u fluid-only, p, incl. pinned solid pressures)."""
        from petsc4py import PETSc
        W = self.W; im = W.dofmap.index_map; lo = im.local_range[0]; nown = im.size_local
        fc, sc = self.cell_tags.find(IDS["fluid"]), self.cell_tags.find(IDS["solid"])
        def cell_dofs(i, comp, cells):
            sp = W.sub(i) if comp is None else W.sub(i).sub(comp)
            Vc, _ = sp.collapse(); d = np.asarray(fem.locate_dofs_topological((sp, Vc), 3, cells)[0]); return set(d[d < nown].tolist())
        reg = {}
        for i, nm in ((0, "u"), (1, "p"), (2, "d")):
            inF, inS = set(), set()
            for c in (range(3) if i != 1 else [None]): inF |= cell_dofs(i, c, fc); inS |= cell_dofs(i, c, sc)
            sub = np.asarray(W.sub(i).collapse()[1]); owned = set(sub[sub < nown].tolist())
            reg[nm + "F"] = (inF - inS) & owned; reg[nm + "I"] = inF & inS & owned; reg[nm + "S"] = (inS - inF) & owned
        mk = lambda sset: PETSc.IS().createGeneral((np.array(sorted(sset), dtype=np.int32) + lo).astype(np.int32), comm=self.msh.comm)
        return mk(reg["uI"] | reg["uS"] | reg["dI"] | reg["dS"]), mk(reg["dF"]), mk(reg["uF"] | reg["pF"] | reg["pI"] | reg["pS"])
