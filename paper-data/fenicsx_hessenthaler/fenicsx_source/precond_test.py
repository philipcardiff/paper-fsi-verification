"""Addendum A3 (a): preconditioner comparison on ONE Newton Jacobian/residual (Gc, mmotion=incr, ckpt t=0.69 -> step to 0.71), 16 ranks.
usage: mpiexec -n 16 python precond_test.py OUTDIR  (OUTDIR contains ckpt_rank*_of16.npz)"""
import sys, os, time, json, resource
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
fsi = HessenthalerFSI(m.mesh, m.cell_tags, m.facet_tags, mode="coupled", prefix="pt_", snes_options={}, mmotion="incr")
z = np.load(os.path.join(out, f"ckpt_rank{comm.rank}_of{comm.size}.npz"))
fsi.w.x.array[:] = z["w"]; fsi.w0.x.array[:] = z["w0"]; fsi.wm.x.array[:] = z["wm"]; fsi.t = float(z["t"]); fsi.nsteps = int(z["nsteps"]); fsi.dt_prev = 0.02
fsi.dt.value = 0.02; fsi.set_time(fsi.t + 0.02)
c = (1.5, -2.0, 0.5)
for k, v in zip(fsi.bdf, c): k.value = v
s = fsi.snes
x = fsi.w.x.petsc_vec.copy(); fpetsc.set_bc(x, fsi.bcs)
s._residual(None, x, s.b); s._jacobian(None, x, s.A, s.A)
A, b = s.A, s.b
W = fsi.W; im = W.dofmap.index_map; lo = im.local_range[0]; nown = im.size_local
def IS(i):
    mp = np.asarray(W.sub(i).collapse()[1]); mp = mp[mp < nown]
    return PETSc.IS().createGeneral((mp + lo).astype(np.int32), comm=comm)
isu, isp, isd = IS(0), IS(1), IS(2)
N = A.getSize()[0]
res = dict(ranks=comm.size, dofs=int(N), nnz=int(A.getInfo()["nz_used"]), b_norm=b.norm())
def run(name, opts, fs=False):
    pre = f"c{name}_"; o = PETSc.Options()
    for k, v in opts.items(): o[pre + k] = v
    ksp = PETSc.KSP().create(comm); ksp.setOptionsPrefix(pre); ksp.setOperators(A)
    if fs: ksp.setType("gmres"); ksp.getPC().setType("fieldsplit"); ksp.getPC().setFieldSplitIS(("u", isu), ("p", isp), ("d", isd))
    ksp.setFromOptions()
    comm.Barrier(); t0 = time.time(); ksp.setUp(); comm.Barrier(); tsetup = time.time() - t0
    xx = A.createVecRight(); t0 = time.time(); ksp.solve(b, xx); comm.Barrier(); tsolve = time.time() - t0
    r = A.createVecLeft(); A.mult(xx, r); r.axpy(-1.0, b)
    rss = comm.allreduce(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, op=MPI.SUM) / 1e6
    d = dict(setup_s=tsetup, solve_s=tsolve, its=ksp.getIterationNumber(), reason=int(ksp.getConvergedReason()), rel_res=r.norm() / b.norm(), peak_rss_gb_sum=rss,
             cpu_h_setup_plus_solve=comm.size * (tsetup + tsolve) / 3600)
    res[name] = d
    if comm.rank == 0: print(name, d, flush=True)
    ksp.destroy()
ONLY = os.environ.get("CFG", "all")
if ONLY in ("all", "A"):
    run("A_fieldsplit_gamg", {"ksp_type": "gmres", "ksp_rtol": 1e-6, "ksp_max_it": 100, "ksp_gmres_restart": 100, "pc_fieldsplit_type": "multiplicative",
        "fieldsplit_u_ksp_type": "preonly", "fieldsplit_u_pc_type": "gamg", "fieldsplit_p_ksp_type": "preonly", "fieldsplit_p_pc_type": "jacobi",
        "fieldsplit_d_ksp_type": "preonly", "fieldsplit_d_pc_type": "gamg"}, fs=True)
if ONLY in ("all", "Aright"):
    run("A_right_pc", {"ksp_type": "gmres", "ksp_pc_side": "right", "ksp_rtol": 1e-6, "ksp_max_it": 200, "ksp_gmres_restart": 200, "pc_fieldsplit_type": "multiplicative",
        "fieldsplit_u_ksp_type": "preonly", "fieldsplit_u_pc_type": "gamg", "fieldsplit_p_ksp_type": "preonly", "fieldsplit_p_pc_type": "jacobi",
        "fieldsplit_d_ksp_type": "preonly", "fieldsplit_d_pc_type": "gamg"}, fs=True)
if ONLY in ("all", "C"):
    run("C_asm_ilu1", {"ksp_type": "gmres", "ksp_rtol": 1e-6, "ksp_max_it": 100, "ksp_gmres_restart": 100, "pc_type": "asm", "pc_asm_overlap": 1, "sub_pc_type": "ilu", "sub_pc_factor_levels": 1})
if ONLY in ("all", "R"):
    run("R_lu_mumps", {"ksp_type": "preonly", "pc_type": "lu", "pc_factor_mat_solver_type": "mumps", "mat_mumps_icntl_14": 80, "mat_mumps_icntl_28": 2, "mat_mumps_icntl_29": 2})
if comm.rank == 0: json.dump(res, open(os.path.join(out, "precond_test_" + ONLY + ".json"), "w"), indent=1)
