"""Phase 2 scaling: ONE rank-independent saved state (zero state, first Newton iterate of step 1: t = 0.02, BE), Gc mesh. Times assembly of F and J, MUMPS LU setup (A4 options),
repeated triangular solves, and peak RSS for the given number of ranks. usage: mpiexec -n N python scaling_test.py OUT.json [mesh]"""
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
comm = MPI.COMM_WORLD
mesh = sys.argv[2] if len(sys.argv) > 2 else "Gc"
m = mesh_io.load(f"{__import__('os').environ.get('HESS_ROOT', '.')}/geometry/meshes/{mesh}/hess_{mesh}_clean.npz", comm)
t0 = time.time(); fsi = HessenthalerFSI(m.mesh, m.cell_tags, m.facet_tags, mode="coupled", prefix="sc_", snes_options={}, mmotion="incr"); tsetup = time.time() - t0
fsi.dt.value = 0.02; fsi.set_time(0.02)
for k, v in zip(fsi.bdf, (1.0, -1.0, 0.0)): k.value = v
s = fsi.snes; x = fsi.w.x.petsc_vec.copy(); fpetsc.set_bc(x, fsi.bcs)
def timed(f, n=3):
    comm.Barrier(); t = time.time()
    for _ in range(n): f()
    comm.Barrier(); return (time.time() - t) / n
tF = timed(lambda: s._residual(None, x, s.b)); tJ = timed(lambda: s._jacobian(None, x, s.A, s.A))
A, b = s.A, s.b
o = PETSc.Options()
for k, v in {"ksp_type": "preonly", "pc_type": "lu", "pc_factor_mat_solver_type": "mumps", "mat_mumps_icntl_14": 80, "mat_mumps_icntl_28": 2, "mat_mumps_icntl_29": 2}.items(): o["scl_" + k] = v
ksp = PETSc.KSP().create(comm); ksp.setOptionsPrefix("scl_"); ksp.setOperators(A); ksp.setFromOptions()
comm.Barrier(); t0 = time.time(); ksp.setUp(); comm.Barrier(); tfac = time.time() - t0
xx = A.createVecRight(); tsol = timed(lambda: ksp.solve(b, xx), 5)
r = A.createVecLeft(); A.mult(xx, r); r.axpy(-1.0, b)
rss = comm.allreduce(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, op=MPI.SUM) / 1e6
res = dict(ranks=comm.size, mesh=mesh, dofs=int(A.getSize()[0]), t_form_setup=tsetup, t_assemble_F=tF, t_assemble_J=tJ, t_lu_setup=tfac, t_solve=tsol, rel_res=r.norm() / b.norm(), peak_rss_gb_sum=rss)
if comm.rank == 0: json.dump(res, open(sys.argv[1], "w"), indent=1); print(res, flush=True)
