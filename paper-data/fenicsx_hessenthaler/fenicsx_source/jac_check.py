"""Phase 2 solver diagnostic: exact Jacobian vs central finite differences of the assembled residual at one state (Gc, mmotion=incr, t = 0.69 -> step 0.71, first Newton iterate),
for random directions restricted to each field (u, p, d) and to the interface/solid/fluid dofs. usage: mpiexec -n 16 python jac_check.py CKPTDIR"""
import sys, os, json
import numpy as np
sys.path.insert(0, __import__("os").environ.get("HESS_ROOT", "."))
def _h(t, v, tb):
    import traceback; traceback.print_exception(t, v, tb); sys.stderr.flush()
    from mpi4py import MPI as M; M.COMM_WORLD.Abort(1)
sys.excepthook = _h
from mpi4py import MPI
from petsc4py import PETSc
from dolfinx.fem import petsc as fpetsc
from src.hessenthaler import mesh_io
from src.hessenthaler.fsi3d import HessenthalerFSI
comm = MPI.COMM_WORLD; out = sys.argv[1]
m = mesh_io.load(__import__("os").environ.get("HESS_ROOT", ".") + "/geometry/meshes/Gc/hess_Gc_clean.npz", comm)
fsi = HessenthalerFSI(m.mesh, m.cell_tags, m.facet_tags, mode="coupled", prefix="jc_", snes_options={}, mmotion=os.environ.get("MM", "incr"))
z = np.load(os.path.join(out, f"ckpt_rank{comm.rank}_of{comm.size}.npz"))
fsi.w.x.array[:] = z["w"]; fsi.w0.x.array[:] = z["w0"]; fsi.wm.x.array[:] = z["wm"]; fsi.t = float(z["t"]); fsi.nsteps = int(z["nsteps"]); fsi.dt_prev = 0.02
fsi.dt.value = 0.02; fsi.set_time(fsi.t + 0.02)
for k, v in zip(fsi.bdf, (1.5, -2.0, 0.5)): k.value = v
s = fsi.snes
x0 = fsi.w.x.petsc_vec.copy()
# initial guess = linear predictor (A3b) so that the test is at a representative iterate
x0.array[:] = fsi.w0.x.array[:x0.getLocalSize()] + (fsi.w0.x.array[:x0.getLocalSize()] - fsi.wm.x.array[:x0.getLocalSize()])
fpetsc.set_bc(x0, fsi.bcs)
def F(x):
    b = s.b.duplicate(); s._residual(None, x, b); return b
s._jacobian(None, x0, s.A, s.A); A = s.A
W = fsi.W; im = W.dofmap.index_map; lo = im.local_range[0]; nown = im.size_local
res = {}
rng = np.random.default_rng(comm.rank + 1)
def sub_idx(i):
    mp = np.asarray(W.sub(i).collapse()[1]); return mp[mp < nown]
# non-BC mask
bcd = np.concatenate([np.asarray(bc.dof_indices()[0]) for bc in fsi.bcs]); bcd = bcd[bcd < nown]
sets = {"all": np.arange(nown), "u": sub_idx(0), "p": sub_idx(1), "d": sub_idx(2)}
scells = fsi.cell_tags.find(2)
for name, idx in sets.items():
    idx = np.setdiff1d(idx, bcd)
    for eps in (1e-4, 1e-6, 1e-8):
        v = x0.duplicate(); v.set(0.0); v.array[idx] = rng.standard_normal(len(idx)); v.assemble()
        # scale: eps relative to |x| (dof scales differ): use v scaled to 1e-3 * typical state scale
        sc = 1e-3 * max(x0.norm(), 1e-12) / max(v.norm(), 1e-12); v.scale(sc * eps / 1e-4)
        xp = x0 + v; xm = x0 - v
        fp, fm = F(xp), F(xm)
        fd = fp - fm; fd.scale(0.5)
        jv = A.createVecLeft(); A.mult(v, jv)
        # BC rows are identity in A; zero them in both
        for vec in (fd, jv): vec.array[bcd] = 0.0
        e = fd - jv
        rel = e.norm() / max(jv.norm(), 1e-300)
        res[f"{name}_eps{eps:g}"] = rel
        if comm.rank == 0: print(name, eps, "rel err J v vs FD:", rel, flush=True)
if comm.rank == 0: json.dump(res, open(os.path.join(out, "jac_check.json"), "w"), indent=1)
