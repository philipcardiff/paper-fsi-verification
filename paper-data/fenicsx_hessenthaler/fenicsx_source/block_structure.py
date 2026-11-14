"""Inspect the ACTUAL unknown ordering and block structure of our monolithic Jacobian (Gc, mmotion=incr, state t=0.69 -> step 0.71, predictor guess), 16 ranks.
Index sets by physical region: u/p/d x {fluid-only, interface, solid-only}. Prints nnz and Frobenius norm of every block, ordering info, symmetry measure and diagonal info."""
import sys, os, json
import numpy as np
sys.path.insert(0, __import__("os").environ.get("HESS_ROOT", "."))
def _h(t, v, tb):
    import traceback; traceback.print_exception(t, v, tb); sys.stderr.flush()
    from mpi4py import MPI as M; M.COMM_WORLD.Abort(1)
sys.excepthook = _h
from mpi4py import MPI
from petsc4py import PETSc
from dolfinx import fem
from dolfinx.fem import petsc as fpetsc
from src.hessenthaler import mesh_io
from src.hessenthaler.fsi3d import HessenthalerFSI
comm = MPI.COMM_WORLD; out = sys.argv[1]
m = mesh_io.load(__import__("os").environ.get("HESS_ROOT", ".") + "/geometry/meshes/Gc/hess_Gc_clean.npz", comm)
fsi = HessenthalerFSI(m.mesh, m.cell_tags, m.facet_tags, mode="coupled", prefix="bs_", snes_options={}, mmotion="incr")
z = np.load(os.path.join(out, f"ckpt_rank{comm.rank}_of{comm.size}.npz"))
fsi.w.x.array[:] = z["w"]; fsi.w0.x.array[:] = z["w0"]; fsi.wm.x.array[:] = z["wm"]; fsi.t = float(z["t"]); fsi.nsteps = int(z["nsteps"]); fsi.dt_prev = 0.02
fsi.dt.value = 0.02; fsi.set_time(fsi.t + 0.02)
for k, v in zip(fsi.bdf, (1.5, -2.0, 0.5)): k.value = v
s = fsi.snes; x = fsi.w.x.petsc_vec.copy()
n = x.getLocalSize(); x.array[:] = fsi.w0.x.array[:n] + (fsi.w0.x.array[:n] - fsi.wm.x.array[:n]); fpetsc.set_bc(x, fsi.bcs)
s._residual(None, x, s.b); s._jacobian(None, x, s.A, s.A); A = s.A
W = fsi.W; im = W.dofmap.index_map; lo = im.local_range[0]; nown = im.size_local
fc, sc = fsi.cell_tags.find(1), fsi.cell_tags.find(2)
def cell_dofs(i, comp, cells):
    sp = W.sub(i) if comp is None else W.sub(i).sub(comp)
    Vc, _ = sp.collapse()
    d = fem.locate_dofs_topological((sp, Vc), 3, cells)[0]; d = np.asarray(d); return set(d[d < nown].tolist())
regions = {}
ncomp = {0: 3, 1: None, 2: 3}
for i, nm in ((0, "u"), (1, "p"), (2, "d")):
    inF, inS = set(), set()
    for c in (range(3) if ncomp[i] else [None]):
        inF |= cell_dofs(i, c, fc); inS |= cell_dofs(i, c, sc)
    sub = np.asarray(W.sub(i).collapse()[1]); owned = set(sub[sub < nown].tolist())
    regions[nm + "_fluid"] = sorted((inF - inS) & owned); regions[nm + "_iface"] = sorted(inF & inS & owned); regions[nm + "_solid"] = sorted((inS - inF) & owned)
names = [k for k in regions if len(regions[k]) or True]
IS = {k: PETSc.IS().createGeneral((np.array(v, dtype=np.int32) + lo).astype(np.int32), comm=comm) for k, v in regions.items()}
sizes = {k: comm.allreduce(len(v), op=MPI.SUM) for k, v in regions.items()}
res = dict(sizes=sizes, nnz={}, fro={}, ordering={})
tot = A.norm(PETSc.NormType.FROBENIUS)
for a in names:
    for b in names:
        if sizes[a] == 0 or sizes[b] == 0: continue
        B = A.createSubMatrix(IS[a], IS[b])
        nz = comm.allreduce(int(B.getInfo()["nz_used"]), op=MPI.SUM); nr = B.norm(PETSc.NormType.FROBENIUS)
        if nz: res["nnz"][f"{a}|{b}"] = nz; res["fro"][f"{a}|{b}"] = nr / tot
        B.destroy()
# ordering of unknowns inside W (first 12 local dofs of rank 0 -> field label)
lab = np.full(nown, -1)
for i in range(3):
    sub = np.asarray(W.sub(i).collapse()[1]); lab[sub[sub < nown]] = i
res["ordering"]["rank0_first_30_field_labels(0=u,1=p,2=d)"] = lab[:30].tolist()
res["ordering"]["bs_of_W_dofmap"] = int(W.dofmap.index_map_bs)
res["matrix"] = dict(rows=int(A.getSize()[0]), nnz_total=int(A.getInfo()["nz_used"]), fro=tot)
if comm.rank == 0:
    json.dump(res, open(os.path.join(out, "block_structure.json"), "w"), indent=1)
    print(json.dumps(res["sizes"]))
    keys = [k for k in names if sizes[k]]
    print("Frobenius norm of block / total (x1e3); rows=test (residual), cols=unknown")
    print(" " * 10 + " ".join("%9s" % k[:9] for k in keys))
    for a in keys: print("%-10s" % a[:10] + " ".join("%9.3f" % (1e3 * res["fro"].get(f"{a}|{b}", 0.0)) for b in keys))
