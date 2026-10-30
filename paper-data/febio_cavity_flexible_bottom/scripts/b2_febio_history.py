#!/usr/bin/env python3
import os, os as _os; HERE = _os.path.dirname(_os.path.abspath(__file__)); PKG = _os.path.dirname(HERE)   # package root; run directories of retained HPC outputs: env CAVITY_RUNS (default runs/cavityNext)
"""B2: FEBio time histories from Mon_iface.txt (nodal ef, ux, uy of the z=0 interface nodes) and pointA.txt: F_y(t), F_x(t) (pressure only, deformed polyline, per unit depth), even/odd-branch forces, u_y(A,t); window statistics and oscillation decay.
usage: b2_febio_history.py RUNDIR OUT.json   (reads only)"""
import sys, os, re, json, numpy as np
rd, outp = sys.argv[1], sys.argv[2]
meta = json.load(open(os.path.join(rd, "model_meta.json"))); tag, K = meta["tag"], meta["k"]
f = open(os.path.join(rd, tag + ".feb"), encoding="latin-1").read()
ids = np.array(re.search(r'<NodeSet name="InterfaceZ0">([^<]+)</NodeSet>', f).group(1).split(","), int)
X0 = {}
want = set(ids.tolist())
for m in re.finditer(r'<node id="(\d+)">([^<]+)</node>', f):
    n = int(m.group(1))
    if n in want: X0[n] = [float(v) for v in m.group(2).split(",")]
xr = np.array([X0[n][0] for n in ids]); o = np.argsort(xr); ids = ids[o]; xr = xr[o]; pos = {n: i for i, n in enumerate(ids)}
T, FY, FX, FYe, FYo, PM = [], [], [], [], [], []
cur = None
with open(os.path.join(rd, "Mon_iface.txt")) as fh:
    t = None; ef = np.zeros(len(ids)); ux = np.zeros(len(ids)); uy = np.zeros(len(ids))
    def flush():
        P = -K * ef; X = xr + ux; Y = uy; dx, dy = np.diff(X), np.diff(Y)
        def F(Pn): pb = 0.5 * (Pn[1:] + Pn[:-1]); return float(-np.sum(pb * dx)), float(np.sum(pb * dy))
        i = np.arange(len(P)); Pe = np.where(i % 2 == 0, P, np.roll(P, -1)); Po = np.where(i % 2 == 1, P, np.roll(P, 1)); Pe[-1] = P[-1] if (len(P) - 1) % 2 == 0 else P[-2]; Po[0] = P[1]
        fy, fx = F(P); T.append(t); FY.append(fy); FX.append(fx); FYe.append(F(Pe)[0]); FYo.append(F(Po)[0]); PM.append(float(P.mean()))
    for ln in fh:
        if ln.startswith("*Time"):
            if t is not None: flush()
            t = float(ln.split("=")[1])
        elif ln.startswith("*"): continue
        elif ln.strip():
            v = ln.split(); n = int(v[0])
            if n in pos: j = pos[n]; ux[j], uy[j], ef[j] = float(v[1]), float(v[2]), float(v[5])
    flush()
T = np.array(T); FY = np.array(FY); FX = np.array(FX); FYe = np.array(FYe); FYo = np.array(FYo)
A = {}; tt = None
for ln in open(os.path.join(rd, "pointA.txt")):
    if ln.startswith("*Time"): tt = float(ln.split("=")[1])
    elif not ln.startswith("*") and ln.strip(): A[tt] = list(map(float, ln.split()[1:]))
tA = np.array(sorted(A)); uyA = np.array([A[x][3] for x in tA]); uxA = np.array([A[x][2] for x in tA])
res = {"tag": tag, "n_times": int(len(T)), "t_last": float(T[-1]), "note": "F in N/m (pressure only, nodal -K ef, trapezoid on the deformed z=0 polyline); even/odd = replace every other node by its neighbour (branch bound); u_y(A) from pointA.txt"}
def win(a, b, name, y, tx):
    m = (tx >= a) & (tx <= b); v = y[m]; return {"window": [a, b], "n": int(m.sum()), "mean": float(v.mean()), "min": float(v.min()), "max": float(v.max()), "spread_rel": float((v.max() - v.min()) / abs(v.mean()))}
tl = T[-1]
res["windows"] = {}
for w in (50, 100, 200):
    res["windows"][f"last_{w}s"] = {"Fy": win(tl - w, tl, "Fy", FY, T), "Fy_even": win(tl - w, tl, "", FYe, T), "Fy_odd": win(tl - w, tl, "", FYo, T), "Fx": win(tl - w, tl, "", FX, T), "uyA": win(tl - w, tl, "", uyA, tA)}
# oscillation: local extrema of uyA and Fy about the last-100 s mean, envelope decay
def extrema(tx, y):
    d = np.diff(y); idx = [i for i in range(1, len(y) - 1) if (y[i] - y[i - 1]) * (y[i + 1] - y[i]) < 0]; return [(float(tx[i]), float(y[i])) for i in idx]
for nm, tx, y in (("uyA", tA, uyA), ("Fy", T, FY)):
    ref = float(y[tx >= tl - 100].mean()); ex = [(a, b - ref) for a, b in extrema(tx, y) if a > 0.3 * tl]
    res[f"extrema_{nm}"] = {"reference_mean_last100": ref, "extrema_t_dev": ex}
    if len(ex) >= 4:
        tt_ = np.array([e[0] for e in ex]); aa = np.log(np.abs([e[1] for e in ex]) + 1e-30); c = np.polyfit(tt_, aa, 1)
        per = float(2 * np.mean(np.diff(tt_))); res[f"extrema_{nm}"].update({"decay_rate_per_s_logfit": float(c[0]), "efold_time_s": float(-1 / c[0]) if c[0] < 0 else None, "period_s": per})
res["final"] = {"Fy": float(FY[-1]), "Fy_even": float(FYe[-1]), "Fy_odd": float(FYo[-1]), "Fx": float(FX[-1]), "uyA": float(uyA[-1]), "uxA": float(uxA[-1]), "p_mean": float(PM[-1])}
res["series_every"] = {"t": T[::4].tolist(), "Fy": FY[::4].tolist(), "Fy_even": FYe[::4].tolist(), "Fy_odd": FYo[::4].tolist(), "Fx": FX[::4].tolist(), "uyA": [float(A[x][3]) for x in T[::4] if x in A]}
json.dump(res, open(outp, "w"), indent=1)
print(tag, "final", res["final"]); 
for k, v in res["windows"].items(): print(k, "Fy mean %.5f spread %.4f%%  uyA mean %.6f spread %.4f%%" % (v["Fy"]["mean"], 100 * v["Fy"]["spread_rel"], v["uyA"]["mean"], 100 * v["uyA"]["spread_rel"]))
for nm in ("uyA", "Fy"):
    e = res[f"extrema_{nm}"]; print(nm, {k: e.get(k) for k in ("reference_mean_last100", "efold_time_s", "period_s")}, [(round(a), round(b, 5)) for a, b in e["extrema_t_dev"][-8:]])
