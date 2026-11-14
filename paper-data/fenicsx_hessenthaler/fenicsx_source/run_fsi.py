"""Drive HessenthalerFSI. usage: mpiexec -n N python run_fsi.py MESH.msh OUTDIR mode(rigid|coupled) dt t_end [restart_dir] [key=val ...]
Writes history.csv (append), diag.json, checkpoints (w, w0, wm as .npy per rank not used: bp4 via adios2 if available) -> here simple numpy gather on rank 0 is avoided;
uses dolfinx.io.VTXWriter-free HDF5 checkpoints through numpy rank-local arrays (restart valid for the SAME mpi size and mesh partition)."""
import sys, os, time, json
import numpy as np

def _abort_hook(t, v, tb):
    import traceback
    traceback.print_exception(t, v, tb); sys.stderr.flush(); sys.stdout.flush()
    from mpi4py import MPI as _M
    _M.COMM_WORLD.Abort(1)
sys.excepthook = _abort_hook   # a failing rank must not leave the other ranks hanging in a collective

sys.path.insert(0, __import__("os").environ.get("HESS_ROOT", "."))
from mpi4py import MPI
from petsc4py import PETSc
from src.hessenthaler import mesh_io
from src.hessenthaler.fsi3d import HessenthalerFSI, IDS
from src.hessenthaler import nh3d

import faulthandler; faulthandler.dump_traceback_later(int(os.environ.get('FH_S', '240')), repeat=True)
comm = MPI.COMM_WORLD
mesh_path, outdir, mode, dt, t_end = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4]), float(sys.argv[5])
kv = dict(a.split("=") for a in sys.argv[6:])
os.makedirs(outdir, exist_ok=True)
t0 = time.time()
m = mesh_io.load(mesh_path, comm, 0)
msh, ct, ft = m.mesh, m.cell_tags, m.facet_tags
ncell = msh.topology.index_map(3).size_global
if comm.rank == 0: print(f'[setup] mesh loaded {time.time()-t0:.1f}s cells {ncell}', flush=True)
snes_opts = {"snes_type": "newtonls", "snes_linesearch_type": "bt", "snes_rtol": 1e-8, "snes_atol": 1e-12, "snes_stol": 0.0, "snes_max_it": 14,
             "ksp_type": "gmres", "ksp_rtol": 1e-6, "ksp_max_it": 100, "pc_type": "lu", "pc_factor_mat_solver_type": "mumps",
             "snes_lag_preconditioner": int(kv.get("lagpc", 1)), "snes_lag_preconditioner_persists": "true" if int(kv.get("lagpc", 1)) != 1 else "false",
             "mat_mumps_icntl_14": 80, "mat_mumps_icntl_28": 2, "mat_mumps_icntl_29": 2}
if int(kv.get('monitor', 0)):
    snes_opts.update({'snes_monitor': '', 'snes_converged_reason': '', 'ksp_converged_reason': '', 'ksp_monitor': '', 'snes_linesearch_monitor': ''})
for k_, v_ in kv.items():
    if k_.startswith('petsc_'): snes_opts[k_[6:]] = v_
fsi = HessenthalerFSI(msh, ct, ft, mu_s=float(kv.get("mu", nh3d.MU)), mode=mode, snes_options=snes_opts, quad_degree=int(kv.get("quad", 4)), mesh_exp=float(kv.get("mexp", 1.0)), mmotion=kv.get("mmotion", "total"))
if int(kv.get("blockpc", 0)):      # PoC: physics-informed block-Gauss-Seidel fieldsplit (S | A | F), sub-solves: S: MUMPS LU, A: GAMG, F: MUMPS LU (same discrete equations)
    from petsc4py import PETSc as _P
    o_ = _P.Options()
    for k_, v_ in {"pc_type": "fieldsplit", "ksp_pc_side": "right", "pc_fieldsplit_type": "multiplicative", "fieldsplit_S_ksp_type": "preonly", "fieldsplit_S_pc_type": "lu", "fieldsplit_S_pc_factor_mat_solver_type": "mumps",
                   "fieldsplit_A_ksp_type": "preonly", "fieldsplit_A_pc_type": "gamg", "fieldsplit_A_pc_gamg_threshold": 0.0,
                   "fieldsplit_F_ksp_type": "preonly", "fieldsplit_F_pc_type": "lu", "fieldsplit_F_pc_factor_mat_solver_type": "mumps", "ksp_gmres_restart": 100, "ksp_max_it": 200}.items(): o_["fsi_" + k_] = v_
    kspb = fsi.snes.snes.getKSP(); kspb.setType("gmres"); pcb = kspb.getPC(); pcb.setType("fieldsplit"); iS_, iA_, iF_ = fsi.block_index_sets()
    pcb.setFieldSplitIS(("S", iS_), ("A", iA_), ("F", iF_)); kspb.setFromOptions()
ndof = fsi.W.dofmap.index_map.size_global * fsi.W.dofmap.index_map_bs
if comm.rank == 0: print(f"mesh {mesh_path} cells {ncell} W dofs {ndof} ranks {comm.size} setup {time.time()-t0:.1f}s", flush=True)
hist = os.path.join(outdir, "history.csv")
if comm.rank == 0 and not os.path.exists(hist):
    open(hist, "w").write("step,t,dt,converged,reason,newton,fnorm,wall_s,ksp_its,refactor,solve_s,peak_rss_gb,tip_x,tip_y,tip_z,Qin_up,Qin_lo,Qout,backflow_out,detF_f_min,detF_s_min,Fy_fluid,Fz_fluid,Fx_fluid,clamp_y,Ix,Iy,Iz,fbY,kin_rel,u_rms,ar_x,ar_y,ar_z\n")
t, step = 0.0, 0
nsteps = int(round((t_end - 0.0) / dt)) if "tstart" not in kv else int(round((t_end - float(kv["tstart"])) / dt))
ck = os.path.join(outdir, f"ckpt_rank{comm.rank}_of{comm.size}.npz")
if "restart" in kv and os.path.exists(ck):
    z = np.load(ck); fsi.w.x.array[:] = z["w"]; fsi.w0.x.array[:] = z["w0"]; fsi.wm.x.array[:] = z["wm"]; t = float(z["t"]); fsi.nsteps = int(z["nsteps"]); fsi.dt_prev = float(z["dt_prev"]) if "dt_prev" in z.files and float(z["dt_prev"]) > 0 else dt; step = int(z["step"]); fsi.t = t
    fsi.w.x.scatter_forward(); fsi.w0.x.scatter_forward(); fsi.wm.x.scatter_forward()
    if comm.rank == 0: print("restarted at t", t, flush=True)
if int(kv.get("post_only", 0)):      # diagnostics of the restored state (no solve)
    out = {}
    out["t"] = t; out["tip_m"] = fsi.tip().tolist()
    Ff, Fp = fsi.interface_force_direct(); out["F_on_solid_direct_N"] = Ff.tolist(); out["F_pressure_part_N"] = Fp.tolist()
    fw, sw = fsi.interface_force_weak(); out["F_on_solid_weak_fluidresidual_N"] = fw.tolist(); out["solid_residual_at_interface_dofs_N"] = sw.tolist()
    out["clamp_residual_N"] = fsi.clamp_reaction().tolist(); out["detF"] = fsi.detF_stats()
    out["Q_up"] = fsi.flux(IDS["upper"]); out["Q_lo"] = fsi.flux(IDS["lower"]); out["Q_out"] = fsi.flux(IDS["outlet"])
    z_, cl_ = fsi.centreline(66); out["centreline_z_m"] = z_.tolist(); out["centreline_d_m"] = cl_.tolist()
    # kinematic condition u = (a0 d1 + a1 d0 + a2 dm)/dt at all dofs: discrete residual norm relative to |u| (solid cells)
    if comm.rank == 0:
        json.dump(out, open(os.path.join(outdir, "post_t%.3f.json" % t), "w"), indent=1); print(json.dumps({k: v for k, v in out.items() if "centreline" not in k}), flush=True)
    sys.exit(0)
ckpt_every = int(kv.get("ckpt", 10)); failed = False; kits_prev = 0; kcum = 0
if int(kv.get("diag_only", 0)):      # exercise all diagnostics on the initial state (no solve)
    fsi.set_time(0.0)
    print("tip", fsi.tip(), "flux", fsi.flux(IDS["upper"]), fsi.flux(IDS["outlet"]), fsi.backflow(IDS["outlet"]), "detF", fsi.detF_stats(), flush=True)
    print("force", fsi.interface_force_direct(), "clamp", fsi.clamp_reaction(), "weak", fsi.interface_force_weak(), flush=True)
    z_, cl_ = fsi.centreline(5); print("centreline", cl_, flush=True)
    sys.exit(0)
while t < t_end - 1e-12:
    ts = time.time(); t_old = t
    dt_now = float(kv['dt2']) if ('dt2' in kv and t >= float(kv.get('tsw', 0.5)) - 1e-9) else dt
    _ts = time.time(); res = fsi.step(dt_now, predictor=bool(int(kv.get('predictor', 0)))); nretry = 0; solve_s = time.time() - _ts
    while not res.converged and nretry < (0 if int(kv.get('noretry', 0)) else 3):      # recovery: (1) refactor the preconditioner, (2,3) halve dt; every attempt is logged
        nretry += 1
        if comm.rank == 0: print(f"retry {nretry} at t={t_old+dt_now:.4f} reason {res.reason} its {res.iterations} fnorm {res.fnorm:.2e}", flush=True)
        fsi.revert(t_old); fsi.snes.snes.getKSP().getPC().reset(); kcum = 0
        if nretry >= 2: dt_now *= 0.5
        res = fsi.step(dt_now, predictor=bool(int(kv.get('predictor', 0))))
    kits = fsi.snes.snes.getLinearSolveIterations(); dk = kits - kits_prev; kits_prev = kits
    refactor = 0
    kcum += dk
    ctrl = os.path.join(outdir, "refactor_now")
    do_ref = (float(kv.get('kcum', 0)) > 0 and kcum > float(kv['kcum'])) or os.path.exists(ctrl)
    do_ref = comm.bcast(do_ref, root=0)
    if do_ref:
        refactor = 1; kcum = 0
        if int(kv.get("blockpc", 0)):      # A5 fix: PCReset drops the fieldsplit splits -> rebuild by creating a NEW KSP with the same prefix/options/index sets and installing it in the SNES
            _old = fsi.snes.snes.getKSP(); _new = PETSc.KSP().create(comm); _new.setOptionsPrefix(_old.getOptionsPrefix()); _new.setType("gmres")
            _pc = _new.getPC(); _pc.setType("fieldsplit"); _iS, _iA, _iF = fsi.block_index_sets(); _pc.setFieldSplitIS(("S", _iS), ("A", _iA), ("F", _iF)); _new.setFromOptions()
            fsi.snes.snes.setKSP(_new)
        else:
            fsi.snes.snes.getKSP().getPC().reset()
        if comm.rank == 0 and os.path.exists(ctrl): os.remove(ctrl)
    ok = bool(res.converged)
    flag_reason = res.reason
    if not ok:
        fsi.revert(t_old); failed = True
        if comm.rank == 0: print(f"STEP FAILED at t={t_old+dt:.4f} reason {flag_reason} its {res.iterations} fnorm {res.fnorm:.3e}", flush=True)
        break
    if int(kv.get('forcediag', 0)) and mode == 'coupled':
        fd = fsi.force_decomp()
        if comm.rank == 0: json.dump(fd, open(os.path.join(outdir, 'forcediag_t%.4f.json' % fsi.t), 'w'), indent=1)
    bal = fsi.balance() if mode == "coupled" else (np.zeros(3), np.zeros(3), 0.0, 0.0)
    fsi.accept(); t = fsi.t; step += 1
    tip = fsi.tip(); qu, ql = fsi.flux(IDS["upper"]), fsi.flux(IDS["lower"]); qo = fsi.flux(IDS["outlet"]); bfo = fsi.backflow(IDS["outlet"])
    dj = fsi.detF_stats()
    if mode == "coupled":
        Ff, _ = fsi.interface_force_direct(); cr = fsi.clamp_reaction()
    else:
        Ff, _ = fsi.interface_force_direct(); cr = np.zeros(3)
    row = [step, t, dt_now, int(ok), flag_reason, res.iterations, res.fnorm, time.time() - ts, dk, refactor, solve_s, comm.allreduce(__import__('resource').getrusage(__import__('resource').RUSAGE_SELF).ru_maxrss, op=MPI.SUM) / 1e6, *tip, qu, ql, qo, bfo, dj["fluid"][0], dj["solid"][0], Ff[1], Ff[2], Ff[0], cr[1], *bal[0], bal[1][1], bal[2], bal[3], *(Ff + bal[1] + cr - bal[0])]
    if comm.rank == 0:
        open(hist, "a").write(",".join(f"{x:.10g}" if isinstance(x, float) else str(x) for x in row) + "\n")
        print(f"step {step} t {t:.4f} newton {res.iterations} |F| {res.fnorm:.2e} wall {time.time()-ts:.1f}s tip_y {tip[1]*1e3:.4f}mm Qin {-(qu+ql):.4e} Qout {qo:.4e} detFmin f/s {dj['fluid'][0]:.3f}/{dj['solid'][0]:.3f}", flush=True)
    if step % ckpt_every == 0:
        np.savez(ck, w=fsi.w.x.array, w0=fsi.w0.x.array, wm=fsi.wm.x.array, t=t, nsteps=fsi.nsteps, step=step, dt_prev=(fsi.dt_prev or 0.0))
np.savez(ck, w=fsi.w.x.array, w0=fsi.w0.x.array, wm=fsi.wm.x.array, t=t, nsteps=fsi.nsteps, step=step, dt_prev=(fsi.dt_prev or 0.0))
if comm.rank == 0:
    json.dump(dict(reached_end=(not failed and t >= t_end - 1e-12), failed=failed, t_final=t, steps=step, cells=int(ncell), dofs=int(ndof), ranks=comm.size, wall_s=time.time() - t0, snes_options=snes_opts,
                   mesh=mesh_path), open(os.path.join(outdir, "diag.json"), "w"), indent=1)
    open(os.path.join(outdir, "EXIT_STATUS"), "w").write("reached_end" if (not failed and t >= t_end - 1e-12) else "FAILED")
