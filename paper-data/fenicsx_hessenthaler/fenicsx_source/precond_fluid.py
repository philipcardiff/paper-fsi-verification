"""Addendum A5 experiment (derived from precond_poc.py). Bounded PoC: physics-informed 3-block (solid+interface S | ALE interior A | fluid F) block-Gauss-Seidel fieldsplit preconditioner on OUR monolithic Jacobian.
Same discrete equations: GMRES on the identical (equilibrated) matrix, true unpreconditioned residual, compared with the exact MUMPS LU solution.
State: Gc, mmotion=incr, t = 0.69 -> step 0.71, predictor guess; 16 ranks. usage: mpiexec -n 16 python precond_poc.py CKPTDIR VARIANTS(comma: ref,V1,V2)"""
import sys, os, time, json, resource
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
comm = MPI.COMM_WORLD; out = sys.argv[1]; variants = sys.argv[2].split(",")
MESH = os.environ.get("MESH", "Gc"); m = mesh_io.load(f"{__import__('os').environ.get('HESS_ROOT', '.')}/geometry/meshes/{MESH}/hess_{MESH}_clean.npz", comm)
fsi = HessenthalerFSI(m.mesh, m.cell_tags, m.facet_tags, mode="coupled", prefix="poc0_", snes_options={}, mmotion=os.environ.get("MM", "incr"))
z = np.load(os.path.join(out, f"ckpt_rank{comm.rank}_of{comm.size}.npz"))
fsi.w.x.array[:] = z["w"]; fsi.w0.x.array[:] = z["w0"]; fsi.wm.x.array[:] = z["wm"]; fsi.t = float(z["t"]); fsi.nsteps = int(z["nsteps"]); fsi.dt_prev = 0.02
fsi.dt.value = 0.02; fsi.set_time(fsi.t + 0.02)
for k, v in zip(fsi.bdf, (1.5, -2.0, 0.5)): k.value = v
s = fsi.snes; x0 = fsi.w.x.petsc_vec.copy(); n = x0.getLocalSize()
x0.array[:] = fsi.w0.x.array[:n] + (fsi.w0.x.array[:n] - fsi.wm.x.array[:n]); fpetsc.set_bc(x0, fsi.bcs)
s._residual(None, x0, s.b); s._jacobian(None, x0, s.A, s.A); A0, b0 = s.A, s.b
W = fsi.W; im = W.dofmap.index_map; lo = im.local_range[0]; nown = im.size_local
fc, sc = fsi.cell_tags.find(1), fsi.cell_tags.find(2)
def cell_dofs(i, comp, cells):
    sp = W.sub(i) if comp is None else W.sub(i).sub(comp)
    Vc, _ = sp.collapse(); d = np.asarray(fem.locate_dofs_topological((sp, Vc), 3, cells)[0]); return set(d[d < nown].tolist())
reg = {}
for i, nm in ((0, "u"), (1, "p"), (2, "d")):
    inF, inS = set(), set()
    for c in (range(3) if i != 1 else [None]): inF |= cell_dofs(i, c, fc); inS |= cell_dofs(i, c, sc)
    sub = np.asarray(W.sub(i).collapse()[1]); owned = set(sub[sub < nown].tolist())
    reg[nm + "F"] = (inF - inS) & owned; reg[nm + "I"] = inF & inS & owned; reg[nm + "S"] = (inS - inF) & owned
mk = lambda sset: PETSc.IS().createGeneral((np.array(sorted(sset), dtype=np.int32) + lo).astype(np.int32), comm=comm)
ISS = mk(reg["uI"] | reg["uS"] | reg["dI"] | reg["dS"]); ISA = mk(reg["dF"]); ISF = mk(reg["uF"] | reg["pF"] | reg["pI"] | reg["pS"])
if comm.rank == 0: print('IS built', flush=True)
ISuF = mk(reg['uF']); ISp = mk(reg['pF'] | reg['pI'])
sz = lambda I: comm.allreduce(I.getLocalSize(), op=MPI.SUM)
info = dict(ranks=comm.size, dofs=int(A0.getSize()[0]), S=sz(ISS), A=sz(ISA), F=sz(ISF))
# ---- equilibration: row then column 2-norm scaling (same equations; reported residuals are for the ORIGINAL system)
def row_norms(M):
    r0, r1 = M.getOwnershipRange(); v = M.createVecLeft(); a = v.array
    for i in range(r1 - r0):
        _, vals = M.getRow(r0 + i); a[i] = np.sqrt(np.sum(vals ** 2))
    return v
rn = row_norms(A0); rn.reciprocal()
if comm.rank == 0: print('row norms done', flush=True)

A1 = A0.duplicate(copy=True)
if comm.rank == 0: print('dup done', flush=True)
A1.diagonalScale(rn, None)
if comm.rank == 0: print('rowscale done', flush=True)
N = A1.getSize()[0]; colsq = np.zeros(N); r0, r1 = A1.getOwnershipRange()
for i in range(r0, r1):
    cols, vals = A1.getRow(i); colsq[cols] += vals ** 2
comm.Allreduce(MPI.IN_PLACE, colsq, op=MPI.SUM)
cn = A1.createVecRight(); c0, c1 = cn.getOwnershipRange(); cn.array[:] = 1.0 / np.sqrt(colsq[c0:c1])
A1.diagonalScale(None, cn)
if comm.rank == 0: print('equilibrated', flush=True)
b1 = b0.copy(); b1.pointwiseMult(b1, rn)        # D_r b ;  x = D_c y
res = dict(info=info)
def rss(): return comm.allreduce(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, op=MPI.SUM) / 1e6
def unscaled_rel_res(y):
    xx = y.copy(); xx.pointwiseMult(xx, cn); r = A0.createVecLeft(); A0.mult(xx, r); r.axpy(-1.0, b0); return r.norm() / b0.norm(), xx
def run_ref():
    ksp = PETSc.KSP().create(comm); o = PETSc.Options()
    for k, v in {"ksp_type": "preonly", "pc_type": "lu", "pc_factor_mat_solver_type": "mumps", "mat_mumps_icntl_14": 80, "mat_mumps_icntl_28": 2, "mat_mumps_icntl_29": 2}.items(): o["ref_" + k] = v
    ksp.setOptionsPrefix("ref_"); ksp.setOperators(A0); ksp.setFromOptions()
    comm.Barrier(); t = time.time(); ksp.setUp(); comm.Barrier(); ts = time.time() - t
    x = A0.createVecRight(); t = time.time(); ksp.solve(b0, x); comm.Barrier(); tv = time.time() - t
    r = A0.createVecLeft(); A0.mult(x, r); r.axpy(-1.0, b0)
    d = dict(setup_s=ts, solve_s=tv, rel_res=r.norm() / b0.norm(), peak_rss_gb=rss(), cpu_h=comm.size * (ts + tv) / 3600)
    if comm.rank == 0: print("ref", d, flush=True)
    return x, d

def mkksp(name, opts):
    pre = f"a5{name}_"; o = PETSc.Options()
    for k, v in opts.items(): o[pre + k] = v
    ksp = PETSc.KSP().create(comm); ksp.setOptionsPrefix(pre); ksp.setOperators(A1)
    ksp.setType("gmres"); pc = ksp.getPC(); pc.setType("fieldsplit"); pc.setFieldSplitIS(("S", ISS), ("A", ISA), ("F", ISF)); ksp.setFromOptions()
    return ksp
def solve(ksp, label, xref=None):
    hist = []; ksp.setMonitor(lambda k, it, r_: hist.append(float(r_)))
    comm.Barrier(); t = time.time(); ksp.setUp(); comm.Barrier(); ts = time.time() - t
    y = A1.createVecRight(); y.set(0.0); t = time.time(); ksp.solve(b1, y); comm.Barrier(); tv = time.time() - t
    try: fr = int(ksp.getPC().getFailedReason())
    except Exception: fr = None
    rr, xx = unscaled_rel_res(y)
    r = A0.createVecLeft(); A0.mult(xx, r); r.axpy(-1.0, b0); bn = b0.norm()
    blk = {}
    for nm, I in (("S", ISS), ("A", ISA), ("F_u", ISuF), ("F_p", ISp)):
        blk[nm] = r.getSubVector(I).norm() / bn
    d = dict(label=label, pc_failed_reason=fr, setup_s=ts, solve_s=tv, total_s=ts + tv, its=ksp.getIterationNumber(), reason=int(ksp.getConvergedReason()), true_rel_res=rr, per_block_rel_residual=blk,
             cpu_h=comm.size * (ts + tv) / 3600, hist=hist[:: max(1, len(hist) // 15)])
    if xref is not None:
        e = xx.copy(); e.axpy(-1.0, xref); d["rel_diff_vs_reference"] = e.norm() / xref.norm()
    if comm.rank == 0: print(label, {k: v for k, v in d.items() if k != "hist"}, flush=True); print("  hist", d["hist"], flush=True)
    return d, xx
base = {"ksp_rtol": 1e-6, "ksp_max_it": 60, "ksp_gmres_restart": 100, "ksp_pc_side": "right", "pc_fieldsplit_type": "multiplicative",
        "fieldsplit_S_ksp_type": "preonly", "fieldsplit_S_pc_type": "lu", "fieldsplit_S_pc_factor_mat_solver_type": "mumps",
        "fieldsplit_A_ksp_type": "preonly", "fieldsplit_A_pc_type": "gamg", "fieldsplit_A_pc_gamg_threshold": 0.0}
V1 = dict(base); V1.update({"fieldsplit_F_ksp_type": "preonly", "fieldsplit_F_pc_type": "lu", "fieldsplit_F_pc_factor_mat_solver_type": "mumps"})
L1 = dict(base); L1.update({"fieldsplit_F_ksp_type": "preonly", "fieldsplit_F_pc_type": "fieldsplit", "fieldsplit_F_pc_fieldsplit_detect_saddle_point": "",
    "fieldsplit_F_pc_fieldsplit_type": "schur", "fieldsplit_F_pc_fieldsplit_schur_fact_type": "upper", "fieldsplit_F_pc_fieldsplit_schur_precondition": "self",
    "fieldsplit_F_fieldsplit_0_ksp_type": "preonly", "fieldsplit_F_fieldsplit_0_pc_type": "asm", "fieldsplit_F_fieldsplit_0_pc_asm_overlap": 1,
    "fieldsplit_F_fieldsplit_0_sub_pc_type": "ilu", "fieldsplit_F_fieldsplit_0_sub_pc_factor_levels": 1,
    "fieldsplit_F_fieldsplit_1_ksp_type": "preonly", "fieldsplit_F_fieldsplit_1_pc_type": "lsc", "fieldsplit_F_fieldsplit_1_pc_lsc_scale_diag": "",
    "fieldsplit_F_fieldsplit_1_lsc_ksp_type": "preonly", "fieldsplit_F_fieldsplit_1_lsc_pc_type": "gamg"})
res = dict(info=info); xref = None
for step in variants:
    if step == "V1":
        k = mkksp("V1", V1); d, xref = solve(k, "V1_reference"); res["V1_reference"] = d; k.destroy()
    elif step == "V1rebuild":      # defect-fix validation: brand-new KSP/PC (as after a refactor), same index sets, same matrix
        k = mkksp("V1b", V1); d, _ = solve(k, "V1_rebuilt_KSP", xref); res["V1_rebuilt_KSP"] = d; k.destroy()
    elif step == "L1a":      # troubleshooting: same LSC Schur, A_uu by exact LU (isolates the ILU/ASM sub-solver from the LSC part)
        o = dict(L1); o.update({"fieldsplit_F_fieldsplit_0_pc_type": "lu", "fieldsplit_F_fieldsplit_0_pc_factor_mat_solver_type": "mumps"})
        for kk in [k_ for k_ in o if k_.startswith("fieldsplit_F_fieldsplit_0_pc_asm") or k_.startswith("fieldsplit_F_fieldsplit_0_sub_")]: o.pop(kk)
        k = mkksp("L1a", o); d, _ = solve(k, "L1a_LSC_with_LU_Auu", xref); res["L1a"] = d
    elif step == "L1b":      # diagnostic (hypothesis: the flap load is a small difference of large face pressures, so a single LSC application is too inexact): F solved by an inner GMRES (rtol 1e-3) preconditioned by the L1a block Schur/LSC
        o = dict(L1); o.update({"fieldsplit_F_fieldsplit_0_pc_type": "lu", "fieldsplit_F_fieldsplit_0_pc_factor_mat_solver_type": "mumps", "fieldsplit_F_ksp_type": "gmres", "fieldsplit_F_ksp_rtol": 1e-3, "fieldsplit_F_ksp_max_it": 30, "ksp_max_it": 40})
        for kk in [k_ for k_ in o if k_.startswith("fieldsplit_F_fieldsplit_0_pc_asm") or k_.startswith("fieldsplit_F_fieldsplit_0_sub_")]: o.pop(kk)
        k = mkksp("L1b", o); d, _ = solve(k, "L1b_inner_GMRES_on_F", xref); res["L1b"] = d
    elif step == "L1":
        k = mkksp("L1", L1); d, _ = solve(k, "L1_fluid_Schur_LSC", xref); res["L1"] = d
if comm.rank == 0: json.dump(res, open(os.path.join(out, "a5_" + "_".join(variants) + ".json"), "w"), indent=1, default=float)
