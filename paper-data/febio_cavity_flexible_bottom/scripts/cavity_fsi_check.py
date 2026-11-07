#!/usr/bin/env python3
import os, os as _os; HERE = _os.path.dirname(_os.path.abspath(__file__)); PKG = _os.path.dirname(HERE)   # package root; run directories of retained HPC outputs: env CAVITY_RUNS (default runs/cavityNext)
"""Post-processing of a coupled cavity FSI run. usage: cavity_fsi_check.py RUNDIR [t_window_frac=0.1] -> RUNDIR/fsi_check.json
Quantities (per unit depth): u_y(A), u_x(A) (top surface node at (4,0), z = 0 layer); steady spread over the last window; vertical/horizontal pressure force on the
DEFORMED interface polyline (p = -K ef, linear between nodes: exact for the trilinear field; viscous part not included, estimated from the fixed-plate run: |Fx_visc| ~ 0.02 N/m);
odd-even amplitude along the interface; inlet/outlet flux; Newton statistics; cost. No filtering is applied to any reported number."""
import sys, os, re, json, numpy as np
sys.path.insert(0, os.path.dirname(__file__)); import cavity_traction as ct
rd = sys.argv[1]; frac = float(sys.argv[2]) if len(sys.argv) > 2 else 0.1
meta = json.load(open(os.path.join(rd, "model_meta.json"))); tag = meta["tag"]; K, mu = meta["k"], meta["mu"]; nodes = ct.read_feb_nodes(os.path.join(rd, tag + ".feb"))
def read_cols(path):   # {t: {id: tuple(values)}}
    out, t = {}, None
    for ln in open(path):
        if ln.startswith("*Time"): t = float(ln.split("=")[1]); out[t] = {}
        elif not ln.startswith("*") and ln.strip(): f = ln.split(); out[t][int(f[0])] = tuple(map(float, f[1:]))
    return out
A = read_cols(os.path.join(rd, "pointA.txt")); aid = meta["point_A_node"]; ts = sorted(A)
uxA = np.array([A[t][aid][2] for t in ts]); uyA = np.array([A[t][aid][3] for t in ts])   # columns after id: x y ux uy
w = np.array(ts) >= ts[-1] * (1 - frac); spread = lambda a: float((a[w].max() - a[w].min()) / max(abs(a[-1]), 1e-30))
res = {"tag": tag, "t_last": ts[-1], "n_times": len(ts), "uy_A_final": float(uyA[-1]), "ux_A_final": float(uxA[-1]), "uy_A_spread_last_window_rel": spread(uyA), "ux_A_spread_last_window_abs": float(uxA[w].max() - uxA[w].min()),
       "uy_A_history": {str(t): float(uyA[k]) for k, t in enumerate(ts) if k % max(1, len(ts) // 25) == 0}}
I = read_cols(os.path.join(rd, "Mon_iface.txt")); iface_ids = sorted(next(iter(I.values())).keys(), key=lambda i: nodes[i][0])
def forces(t):
    X = np.array([nodes[i][0] + I[t][i][0] for i in iface_ids]); Y = np.array([nodes[i][1] + I[t][i][1] for i in iface_ids]); P = np.array([-K * I[t][i][4] for i in iface_ids])
    dx, dy = np.diff(X), np.diff(Y); pb = 0.5 * (P[1:] + P[:-1])
    return dict(Fy_pressure=float(-np.sum(pb * dx)), Fx_pressure=float(np.sum(pb * dy)), p=P, X=X, Y=Y)
Fh = [forces(t) for t in ts if t in I]; tI = [t for t in ts if t in I]
Fy = np.array([f["Fy_pressure"] for f in Fh]); Fx = np.array([f["Fx_pressure"] for f in Fh]); wI = np.array(tI) >= tI[-1] * (1 - frac)
sp = lambda a: float((a[wI].max() - a[wI].min()) / max(abs(a[-1]), 1e-30))
P = Fh[-1]["p"]; idx = np.arange(len(P)); ev, od = P[idx % 2 == 0], P[idx % 2 == 1]
res["interface"] = {"Fy_pressure_raw": float(Fy[-1]), "Fx_pressure_raw": float(Fx[-1]), "Fy_spread_rel": sp(Fy), "p_min": float(P.min()), "p_max": float(P.max()), "p_mean": float(P.mean()),
                    "oddeven_amp_even_minus_odd_Pa": float(ev.mean() - od.mean()), "n_nodes": len(P)}
# bracket from even-only / odd-only interpolation of p on the deformed interface (reported as uncertainty, not applied)
X, Y = Fh[-1]["X"], Fh[-1]["Y"]
def fy_sub(sel):
    Ps = np.interp(np.arange(len(P)), idx[sel], P[sel]); pb = 0.5 * (Ps[1:] + Ps[:-1]); return float(-np.sum(pb * np.diff(X)))
res["interface"]["Fy_even_only"] = fy_sub(idx % 2 == 0); res["interface"]["Fy_odd_only"] = fy_sub(idx % 2 == 1)
res["interface"]["deflection_profile_uy"] = [float(I[tI[-1]][i][1]) for i in iface_ids][::4]
# fluxes
def flux(name):
    d = read_cols(os.path.join(rd, name)); t = max(d); ids = sorted(d[t], key=lambda i: nodes[i][1]); y = [nodes[i][1] for i in ids]; return float(np.trapz([d[t][i][0] for i in ids], y))
res["Q_in"] = flux("Mon_inlet.txt"); res["Q_out"] = flux("Mon_outlet.txt"); res["mass_imbalance_rel"] = res["Q_out"] / res["Q_in"] - 1
log = open(os.path.join(rd, tag + ".log"), errors="ignore").read(); steps = re.split(r"===== beginning time step", log)[1:]
nt = [len(s.split("convergence summary")[0].split("convergence norms")[1:]) for s in steps]
res["newton"] = {"steps": len(steps), "iters_mean": float(np.mean(nt)), "iters_max": int(max(nt)), "steps_degenerate_single_table": int(sum(n < 2 for n in nt)), "first_degenerate_step": next((k for k, n in enumerate(nt, 1) if n < 2), None)}
res["log"] = {"normal_termination": "N O R M A L   T E R M I N A T I O N" in log, "errors": len(re.findall(r"(?i)error", log)), "negative_jacobian": bool(re.search(r"(?i)negative jacobian", log)),
              "cutbacks": len(re.findall(r"(?i)retrying|cutback|reducing time step", log)), "warnings_distinct": sorted(set(re.findall(r"\*\s+([A-Z][^*\n]{8,70})\s+\*", log)))[:6]}
tm = open(os.path.join(rd, "time.txt")).read(); g = lambda p: re.search(p, tm)
res["cost"] = {"wall": g(r"Elapsed \(wall clock\) time.*?: ([\d:.]+)").group(1), "user_s": float(g(r"User time \(seconds\): ([\d.]+)").group(1)), "sys_s": float(g(r"System time \(seconds\): ([\d.]+)").group(1)), "max_rss_MB": int(g(r"Maximum resident set size \(kbytes\): (\d+)").group(1)) / 1024, "cores": json.load(open(os.path.join(rd, "run_meta_start.json")))["cpus_per_task"]}
json.dump(res, open(os.path.join(rd, "fsi_check.json"), "w"), indent=1, default=float)
print(json.dumps({k: v for k, v in res.items() if k != "uy_A_history"}, indent=1, default=float)); print("uy(A) history", res["uy_A_history"])
