#!/usr/bin/env python3
import os, os as _os; HERE = _os.path.dirname(_os.path.abspath(__file__)); PKG = _os.path.dirname(HERE)   # package root; run directories of retained HPC outputs: env CAVITY_RUNS (default runs/cavityNext)
"""Interface traction and conservation diagnostics for the FEBio cavity fluid-only / straight-channel runs (fixed wall).
Reads Mon_all.txt (z = 0 layer: id x y wx wy ef; x,y printed as 0 by FEBio, coordinates come from the .feb) and model_meta.json (K, mu, TZ).
Quantities are per unit depth (the z = 0 layer line integrals; no TZ factor needed because the slab is plane (wz = 0) and the line integral
along x of a nodal field is already per unit depth). Raw nodal pressure p = -K ef; trapezoid over wall nodes = exact integral of the trilinear interpolant.
Wall shear: mu * du/dy at the wall from a 3-point one-sided (quadratic) fit through the first three nodes of the column (second order); the 1st-order
two-node value is also returned. Fy = -int p dx (pressure; the viscous normal stress 2 mu dv/dy = 0 on a flat no-slip wall), Fx = int mu du/dy dx.
usage: cavity_traction.py RUNDIR WALL(floor|bottom) [x0 x1]      floor = cavity floor y = 0 ; bottom = straight channel lower wall y = 1
"""
import sys, os, re, json, numpy as np

def read_feb_nodes(feb):
    return {int(m.group(1)): tuple(map(float, m.group(2).split(","))) for m in re.finditer(r'<node id="(\d+)">([^<]+)</node>', open(feb, encoding="latin-1").read())}

def read_mon(path):
    out, t = {}, None
    for ln in open(path):
        if ln.startswith("*Time"): t = float(ln.split("=")[1]); out[t] = {}
        elif not ln.startswith("*") and ln.strip():
            f = ln.split(); out[t][int(f[0])] = tuple(map(float, f[3:6]))
    return out

def load(rd):
    meta = json.load(open(os.path.join(rd, "model_meta.json"))); nodes = read_feb_nodes(os.path.join(rd, meta["tag"] + ".feb"))
    mon = read_mon(os.path.join(rd, "Mon_all.txt")); return meta, nodes, mon

def wall_columns(nodes, mon, ywall, x0, x1, tol=1e-6):
    t0 = next(iter(mon)); ids = set(mon[t0])
    z0 = {i: nodes[i] for i in ids if abs(nodes[i][2]) < tol}
    wall = sorted([(x, i) for i, (x, y, z) in z0.items() if abs(y - ywall) < tol and x0 - tol <= x <= x1 + tol])
    cols = []
    for x, iw in wall:
        up = sorted([(y, i) for i, (xx, y, z) in z0.items() if abs(xx - x) < tol and y > ywall + tol])[:2]
        cols.append((x, iw, [iw] + [i for y, i in up], [ywall] + [y for y, i in up]))
    return cols

def quantities(meta, nodes, mon, cols, t):
    K, mu = meta["k"], meta["mu"]; xs = np.array([c[0] for c in cols])
    p = np.array([-K * mon[t][c[1]][2] for c in cols])
    tau2 = []; tau1 = []
    for x, iw, ids, ys in cols:
        u = np.array([mon[t][i][0] for i in ids]); y = np.array(ys)
        tau2.append(mu * np.polyfit(y - y[0], u, 2)[1])   # du/dy at the wall of the quadratic through the first 3 nodes
        tau1.append(mu * (u[1] - u[0]) / (y[1] - y[0]))
    tau2, tau1 = np.array(tau2), np.array(tau1)
    Fy_raw = -np.trapz(p, xs)
    idx = np.round((xs - xs[0]) / 0.1).astype(int); ev, od = p[idx % 2 == 0], p[idx % 2 == 1]
    ps = p.copy(); ps[1:-1] = 0.25 * p[:-2] + 0.5 * p[1:-1] + 0.25 * p[2:]    # 1-2-1 smoothing: DIAGNOSTIC ONLY (never a reported result)
    return dict(t=t, Fy_raw=Fy_raw, Fy_smoothed121_DIAGNOSTIC=-np.trapz(ps, xs), Fx_viscous_2nd=np.trapz(tau2, xs), Fx_viscous_1st=np.trapz(tau1, xs),
                oddeven_amp_Pa=float(ev.mean() - od.mean()), p_mean=float(p.mean()), xs=xs, p=p, tau2=tau2, tau1=tau1)

if __name__ == "__main__":
    rd, wall = sys.argv[1], sys.argv[2]; x0, x1 = (float(sys.argv[3]), float(sys.argv[4])) if len(sys.argv) > 4 else (2.0, 6.0)
    meta, nodes, mon = load(rd); ts = sorted(mon); cols = wall_columns(nodes, mon, 0.0 if wall == "floor" else 1.0, x0, x1)
    q = quantities(meta, nodes, mon, cols, ts[-1])
    print(f"n wall nodes {len(cols)} t_last {ts[-1]}")
    for k in ("Fy_raw", "Fy_smoothed121_DIAGNOSTIC", "Fx_viscous_2nd", "Fx_viscous_1st", "oddeven_amp_Pa", "p_mean"): print(k, f"{q[k]:.8g}")
    print("tau2 mid-wall", np.median(q["tau2"][5:-5]) if len(q["tau2"]) > 12 else q["tau2"])
