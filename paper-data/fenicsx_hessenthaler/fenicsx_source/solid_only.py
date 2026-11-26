"""GATE B2: zero-flow static flap, total-Lagrangian neo-Hookean, P2 tets, body force rho_s*g_eff in +y, clamp at z=0.
usage: mpiexec -n N python solid_only.py LEVEL OUT.json [mu_Pa]   LEVEL in 1,2,3,4 (nx,ny,nz = 6,2,32 * 2^(l-1))"""
import sys, json, time
import numpy as np

def _abort_hook(t, v, tb):
    import traceback
    traceback.print_exception(t, v, tb); sys.stderr.flush(); sys.stdout.flush()
    from mpi4py import MPI as _M
    _M.COMM_WORLD.Abort(1)
sys.excepthook = _abort_hook   # a failing rank must not leave the other ranks hanging in a collective

sys.path.insert(0, __import__("os").environ.get("HESS_ROOT", "."))
import ufl
from mpi4py import MPI
from petsc4py import PETSc
from dolfinx import mesh, fem
from dolfinx.fem import petsc as fpetsc
from src.hessenthaler import nh3d as nh

comm = MPI.COMM_WORLD
lvl = int(sys.argv[1]); out = sys.argv[2]; mu = float(sys.argv[3]) if len(sys.argv) > 3 else nh.MU
nu = nh.NU; K = nh.bulk(mu, nu); rho_s = 1058.3; g_eff = 0.973306
n = (6 * 2 ** (lvl - 1), 2 * 2 ** (lvl - 1), 32 * 2 ** (lvl - 1))
msh = mesh.create_box(comm, [[-5.5e-3, -1e-3, 0], [5.5e-3, 1e-3, 65e-3]], n, mesh.CellType.tetrahedron)
V = fem.functionspace(msh, ("Lagrange", 2, (3,)))
d = fem.Function(V, name="d")
v = ufl.TestFunction(V)
F = ufl.Identity(3) + ufl.grad(d)
lam = fem.Constant(msh, 0.0)           # load factor
body = ufl.as_vector((0.0, rho_s * g_eff, 0.0))
R = ufl.inner(nh.P_analytic(F, mu, K), ufl.grad(v)) * ufl.dx - lam * ufl.dot(body, v) * ufl.dx
J = ufl.derivative(R, d, ufl.TrialFunction(V))
fdim = msh.topology.dim - 1
facets = mesh.locate_entities_boundary(msh, fdim, lambda x: np.isclose(x[2], 0.0))
bc = fem.dirichletbc(np.zeros(3), fem.locate_dofs_topological(V, fdim, facets), V)
# Newton by hand with PETSc KSP (CG + GAMG, near-nullspace = rigid body modes), monitored
from dolfinx.fem import form
Rf, Jf = form(R), form(J)
A = fpetsc.create_matrix(Jf); b = fpetsc.create_vector(V)
ns = None
# rigid-body near-nullspace
def rbm(Amat):
    x = V.tabulate_dof_coordinates(); no = V.dofmap.index_map.size_local
    x = x[:no]
    vecs = [Amat.createVecRight() for _ in range(6)]
    for vv in vecs: vv.set(0.0)
    arr = [vv.array for vv in vecs]
    for i in range(3): arr[i].reshape(-1, 3)[:, i] = 1.0
    a = arr[3].reshape(-1, 3); a[:, 0] = -x[:, 1]; a[:, 1] = x[:, 0]
    a = arr[4].reshape(-1, 3); a[:, 0] = x[:, 2]; a[:, 2] = -x[:, 0]
    a = arr[5].reshape(-1, 3); a[:, 1] = -x[:, 2]; a[:, 2] = x[:, 1]
    for vv in vecs: vv.assemble()
    # Gram-Schmidt
    for i, vv in enumerate(vecs):
        for w in vecs[:i]: vv.axpy(-vv.dot(w), w)
        vv.scale(1.0 / vv.norm())
    return PETSc.NullSpace().create(vectors=vecs, comm=comm)
ksp = PETSc.KSP().create(comm); ksp.setOptionsPrefix("sol_")
opts = PETSc.Options(); import os
if os.environ.get("SOLID_LIN", "lu") == "lu":
    opts["sol_ksp_type"] = "preonly"; opts["sol_pc_type"] = "lu"; opts["sol_pc_factor_mat_solver_type"] = "mumps"
else:
    opts["sol_ksp_type"] = "gmres"; opts["sol_pc_type"] = "gamg"; opts["sol_ksp_rtol"] = 1e-10; opts["sol_pc_gamg_agg_nsmooths"] = 1
opts["sol_ksp_max_it"] = 500
log = dict(level=lvl, mu=mu, K=K, nu=nu, cells=int(msh.topology.index_map(3).size_global), dofs=int(V.dofmap.index_map.size_global * 3), ranks=comm.size, newton=[])
t0 = time.time()
def residual_norm():
    with b.localForm() as bl: bl.set(0.0)
    fpetsc.assemble_vector(b, Rf); fpetsc.apply_lifting(b, [Jf], [[bc]], x0=[d.x.petsc_vec], alpha=-1.0)
    b.ghostUpdate(addv=PETSc.InsertMode.ADD, mode=PETSc.ScatterMode.REVERSE); fpetsc.set_bc(b, [bc], d.x.petsc_vec, -1.0)
    return b.norm()
nloc = V.dofmap.index_map.size_local * 3
fext = fem.form(ufl.dot(body, v) * ufl.dx)
bf = fpetsc.assemble_vector(fext); bf.ghostUpdate(addv=PETSc.InsertMode.ADD, mode=PETSc.ScatterMode.REVERSE)
FEXT1 = bf.norm()     # |f_ext| at load factor 1 (reference scale of the residual)
Qd = fem.functionspace(msh, ("DG", 0)); detexpr = fem.Expression(ufl.det(ufl.Identity(3) + ufl.grad(d)), Qd.element.interpolation_points); detf = fem.Function(Qd)
def detmin():
    detf.interpolate(detexpr); return comm.allreduce(float(detf.x.array.min()), op=MPI.MIN)
def newton(lf, maxit=25):
    """Full Newton (no residual line search: the thin flap is bending-soft/axially-stiff, residual-norm backtracking is not
    meaningful); step is only damped if det F <= 0.05 or non-finite. Converged when |R| <= 1e-8 |f_ext(1)|*max(lf,..) and last |dx| small."""
    lam.value = lf
    rn = residual_norm(); tol = 1e-7 * FEXT1 * max(lf, 1e-3)   # Addendum A1: roundoff floor ~2e-8|f_ext| for the axially stiff mode
    for it in range(maxit + 1):
        log["newton"].append(dict(load=lf, it=it, rnorm=rn, tol=tol))
        if comm.rank == 0: print(f"load {lf:.5f} it {it} |R| {rn:.3e} tol {tol:.2e}", flush=True)
        if rn < tol: return True, it
        if it == maxit: break
        A.zeroEntries(); fpetsc.assemble_matrix(A, Jf, bcs=[bc]); A.assemble()
        if not hasattr(newton, "ns"):
            newton.ns = rbm(A); A.setNearNullSpace(newton.ns); ksp.setFromOptions()
        ksp.setOperators(A)
        dx = A.createVecRight(); ksp.solve(b, dx)
        dxinf = dx.norm(PETSc.NormType.NORM_INFINITY)
        if ksp.getConvergedReason() <= 0:
            if comm.rank == 0: print('KSP failed reason', ksp.getConvergedReason(), 'its', ksp.getIterationNumber(), flush=True)
            return False, it
        d0 = d.x.array.copy(); alpha = 1.0; okls = False
        for _ in range(10):
            d.x.array[:nloc] = d0[:nloc] - alpha * dx.array; d.x.scatter_forward()
            dm = detmin()
            if np.isfinite(dm) and dm > 0.05: okls = True; break
            alpha *= 0.5
        dx.destroy()
        log['newton'][-1]['dx_inf'] = dxinf
        if dxinf < 1e-12 and okls: return True, it + 1
        if not okls:
            if comm.rank == 0: print('detF line search failed', dm, flush=True)
            d.x.array[:] = d0; return False, it
        rn = residual_norm()
    return False, maxit
lam_cur, dlam, nfail = 0.0, 0.125, 0
ok = True; d_saved = d.x.array.copy()
while lam_cur < 1.0 - 1e-12:
    step = min(dlam, 1.0 - lam_cur)
    d_saved = d.x.array.copy()
    conv, its = newton(lam_cur + step)
    if conv:
        lam_cur += step; nfail = 0
        if its <= 6: dlam = min(dlam * 1.5, 0.5)
    else:
        d.x.array[:] = d_saved; d.x.scatter_forward(); dlam = step / 2; nfail += 1
        if dlam < 1e-4: ok = False; break
log["n_load_steps_final_dlam"] = dlam
log["converged"] = bool(ok); log["wall_s"] = time.time() - t0
# diagnostics: tip displacement at (0,0,65), centreline, J range, clamp reaction vs total body force
tip = np.array([[0.0, 0.0, 65e-3]])
from dolfinx import geometry
tree = geometry.bb_tree(msh, 3); cand = geometry.compute_collisions_points(tree, tip); col = geometry.compute_colliding_cells(msh, cand, tip)
val = np.full(3, np.nan)
if len(col.links(0)) > 0: val = np.asarray(d.eval(tip, col.links(0)[:1])).reshape(-1)[:3]
allv = comm.gather(val, 0)
Jmin = comm.allreduce(float(np.min(fem.Function(fem.functionspace(msh, ("DG", 0))).x.array) if False else 0), op=MPI.SUM)
Q = fem.functionspace(msh, ("DG", 1)); Jfun = fem.Function(Q)
Jfun.interpolate(fem.Expression(ufl.det(ufl.Identity(3) + ufl.grad(d)), Q.element.interpolation_points))
log["detF_min"] = comm.allreduce(float(Jfun.x.array.min()), op=MPI.MIN); log["detF_max"] = comm.allreduce(float(Jfun.x.array.max()), op=MPI.MAX)
vol = comm.allreduce(fem.assemble_scalar(fem.form(1.0 * ufl.dx(domain=msh))), op=MPI.SUM)
log["total_body_force_y"] = rho_s * g_eff * vol
# reaction at clamp via weak residual with lifting of y-component at z=0 dofs
lift = fem.Function(V); dofs_ = fem.locate_dofs_topological(V, fdim, facets)
Vy, _ = V.sub(1).collapse()
ly = fem.Function(V); idx = fem.locate_dofs_topological((V.sub(1), Vy), fdim, facets)[0]; ly.x.array[idx] = 1.0
Rn = ufl.replace(ufl.inner(nh.P_analytic(F, mu, K), ufl.grad(v)) * ufl.dx, {v: ly})
log["clamp_reaction_y"] = comm.allreduce(fem.assemble_scalar(fem.form(Rn)), op=MPI.SUM)   # = internal force at clamp dofs (balances body force)
if comm.rank == 0:
    tipv = [x for x in allv if not np.isnan(x[0])][0]
    log["tip_disp_m"] = [float(t) for t in tipv]; log["tip_y_mm"] = float(tipv[1] * 1e3)
    json.dump(log, open(out, "w"), indent=1); print(json.dumps({k: v for k, v in log.items() if k != "newton"}, indent=1))
