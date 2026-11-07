#!/usr/bin/env python3
"""Coupled collapsibleChannel analysis (tools/make_cc_feb.py --mode fsi): wall-midpoint history vs the oomph-lib references, first trough (parabolic refinement), rebound/second trough/period,
history error E_hist, per-step applied-load spread F1 vs F2 (net vertical force and first-mode load; same definitions as tools/cc_f1_analyze.py but on the MOVING interface), y-parity
D = p0 - 2 p1 + p2, streamwise a_t, inlet/outlet flux and channel-volume conservation, Newton iterations per step.  usage: cc_fsi_analyze.py RUN_DIR [OUT.json]"""
import sys, os, re, json, numpy as np
rd = os.path.realpath(sys.argv[1]); OUT = sys.argv[2] if len(sys.argv) > 2 else rd + "/fsi_result.json"
RF = os.environ.get("CC_REFERENCE_DIR") or os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "reference")   # PACKAGE EDIT: reference CSV directory (original: ../reference_data/cc)
meta = json.load(open(rd + "/model_meta.json")); M = json.load(open(rd + "/mesh2d.json")); K = meta["k"]
nodes2 = np.array(M["nodes2"]); xs = np.array(M["xs"]); nr = M["nr"]; i0, i1 = M["i0w"], M["i1w"]; F = lambda i, j: i * (nr + 1) + j; nw = i1 - i0
def read(fn):
    out = {}; t = None
    for l in open(fn):
        if l.startswith("*Time"): t = float(l.split("=")[1]); out[t] = {}
        elif not l.startswith("*") and l.strip():
            f = l.split(); out[t][int(f[0])] = [float(v) for v in f[1:]]
    return out
mif = read(rd + "/Mon_iface.txt"); mr = read(rd + "/Mon_rows.txt"); mi = read(rd + "/Mon_inlet.txt"); mo = read(rd + "/Mon_outlet.txt")
times = sorted(t for t in mif if t > 0); iface = M["iface"]
def loadref(name):
    rows = [l.split(",") for l in open(f"{RF}/collapsibleChannel_oomph_{name}.csv") if not l.startswith("#") and not l.startswith("time")]
    d = np.array(rows, float); return d[:, 0], d[:, 1], d[:, 2], d[:, 3]
phi = lambda x: 0.5 * (1 - np.cos(2 * np.pi * (x - 5.0) / 10.0))
xi = xs[i0:i1 + 1]; kmid = nw // 2; kq1 = nw // 4; kq3 = 3 * nw // 4
U = np.array([[mif[t][iface[k] + 1][1] for k in range(nw + 1)] for t in times]); UX = np.array([[mif[t][iface[k] + 1][0] for k in range(nw + 1)] for t in times])
mid = U[:, kmid]; q1 = U[:, kq1]; q3 = U[:, kq3]; T = np.array(times)
def trough(t, u, lo=0.0, hi=0.9):
    m = (t >= lo) & (t <= hi); k = int(np.argmin(np.where(m, u, 9)))
    if k >= len(t) - 1 or k < 1 or not m[k + 1]: return None, None, int(k)     # extremum not bracketed inside the record: no trough reached
    c = np.polyfit(t[k - 1:k + 2], u[k - 1:k + 2], 2); tt = -c[1] / (2 * c[0]); return float(np.polyval(c, tt)), float(tt), int(k)
def extrema(t, u):
    k = np.where((u[1:-1] < u[:-2]) & (u[1:-1] <= u[2:]))[0] + 1; kM = np.where((u[1:-1] > u[:-2]) & (u[1:-1] >= u[2:]))[0] + 1; out = []
    for k_ in k:
        if 1 <= k_ < len(t) - 1: c = np.polyfit(t[k_ - 1:k_ + 2], u[k_ - 1:k_ + 2], 2); tt = -c[1] / (2 * c[0]); out.append(("min", float(tt), float(np.polyval(c, tt))))
    for k_ in kM:
        c = np.polyfit(t[k_ - 1:k_ + 2], u[k_ - 1:k_ + 2], 2); tt = -c[1] / (2 * c[0]); out.append(("max", float(tt), float(np.polyval(c, tt))))
    return sorted(out, key=lambda e: e[1])
res = dict(run=rd, K=K, n_steps=len(times), t_last=times[-1], wall_mid_last=float(mid[-1]))
u1, t1, k1 = trough(T, mid); res.update(trough_u=u1, trough_t=t1, trough_reached=u1 is not None)
res["extrema"] = extrema(T, mid)
refs = {}
for nm in ("reference", "dt40"):
    t_, a_, m_, c_ = loadref(nm); refs[nm] = (t_, a_, m_, c_)
    ru, rt, _ = trough(t_, m_, 0.0, 0.9) if False else (None, None, None)
RU = {"reference": (-0.2160, 0.4685), "dt40": (-0.2129, 0.4735)}
for nm, (ru, rt) in RU.items():
    if u1 is not None: res[f"trough_err_vs_{nm}"] = u1 / ru - 1; res[f"trough_time_err_vs_{nm}"] = t1 / rt - 1
for nm in ("reference", "dt40"):
    t_, a_, m_, c_ = refs[nm]
    for lab, sim, col in (("mid", mid, m_), ("q1", q1, a_), ("q3", q3, c_)):
        for tmax in (1.5, T[-1]):
            sel = T <= tmax + 1e-9
            ri = np.interp(T[sel], t_, col); fin = t_[-1] >= T[sel][-1] - 1e-9
            if not fin: continue
            res[f"E_hist_{lab}_{nm}_to{tmax:g}"] = float(np.max(np.abs(sim[sel] - ri)) / 0.216); res[f"rms_hist_{lab}_{nm}_to{tmax:g}"] = float(np.sqrt(np.mean((sim[sel] - ri) ** 2)) / 0.216)
# applied-load spread per step
rows = []
for j, t in enumerate(times):
    p = lambda i, jj: -K * mr[t][F(i, jj) + 1][2]
    p0 = np.array([p(i, nr) for i in range(i0, i1 + 1)]); p1 = np.array([p(i, nr - 1) for i in range(i0, i1 + 1)]); p2 = np.array([p(i, nr - 2) for i in range(i0, i1 + 1)])
    X = xi + UX[j]; Y = 1.0 + U[j]; dx = np.diff(X); dy = np.diff(Y)
    pa = 0.5 * (p0[:-1] + p0[1:]); pe = 0.25 * (p0[:-1] + p0[1:] + p1[:-1] + p1[1:]); xm = 0.5 * (xi[:-1] + xi[1:])
    F1 = np.sum(pa * dx); F2 = np.sum(pe * dx); G1 = np.sum(pa * phi(xm) * dx); G2 = np.sum(pe * phi(xm) * dx)
    k_ = np.arange(4, nw - 3); at = float(np.sqrt(np.mean(((p0[k_] - 0.5 * (p0[k_ - 1] + p0[k_ + 1])) / 2) ** 2)))
    D = p0 - 2 * p1 + p2; Dm = D[4:-4]
    qin = float(np.trapz([mi[t][F(0, jj) + 1][0] for jj in range(nr + 1)], [nodes2[F(0, jj)][1] for jj in range(nr + 1)]))
    qout = float(np.trapz([mo[t][F(len(xs) - 1, jj) + 1][0] for jj in range(nr + 1)], [nodes2[F(len(xs) - 1, jj)][1] for jj in range(nr + 1)]))
    vol = float(np.trapz(U[j], xi))   # signed area change of the channel (negative = collapse)
    rows.append(dict(t=t, F1y=F1, F2y=F2, G1=G1, G2=G2, spread_F=(F1 - F2) / F2 if abs(F2) > 1e-12 else None, spread_G=(G1 - G2) / G2 if abs(G2) > 1e-12 else None, a_t=at, D_rms=float(np.sqrt(np.mean(Dm ** 2))), D_mean=float(Dm.mean()), D_maxabs=float(np.abs(Dm).max()),
                     p_wall_range=float(p0[4:-4].max() - p0[4:-4].min()), Qin=qin, Qout=qout, vol=vol, mid=float(mid[j])))
res["rows"] = rows
sel = [r for r in rows if r["t"] >= 0.10 and abs(r["F2y"]) >= 100.0]
res["load_spread_steps"] = len(sel)
if sel:
    res["max_spread_F"] = max(abs(r["spread_F"]) for r in sel); res["max_spread_G"] = max(abs(r["spread_G"]) for r in sel)
    res["max_D_rms"] = max(r["D_rms"] for r in rows); res["max_abs_F2y"] = max(abs(r["F2y"]) for r in rows)
# solid free-body (quasi-static + inertia) identification of the uniform-equivalent fluid pressure from the wall midpoint history and the S1 static curve (S1_320x8: all 20 load steps)
S1 = os.path.join(os.path.dirname(rd), "S1_320x8", "Mon_bottom.txt")
if os.path.exists(S1) and nw == 320:
    steps = read(S1); ts = sorted(steps); umid = np.array([steps[tt][iface_id][3] for tt in ts]) if False else None
    ids = sorted(steps[ts[0]]); mm = {tt: np.array([steps[tt][i] for i in ids]) for tt in ts}
    pe = np.array([200.0 * tt for tt in ts]); um = np.array([mm[tt][len(ids) // 2][3] for tt in ts])     # columns: x y ux uy -> uy index 3
    psi = mm[ts[-1]][:, 3] / mm[ts[-1]][len(ids) // 2][3]; psim = 0.5 * (psi[:-1] + psi[1:]); xm0 = 0.5 * (xi[:-1] + xi[1:]); dx0 = np.diff(xi)
    Ipsi = float(np.sum(psim * dx0)); cin = float(np.sum(psim ** 2 * dx0) / Ipsi)
    pel = lambda u: float(np.interp(-u, -um[::1] if False else -um, pe))   # uniform-equivalent elastic restoring pressure at midpoint displacement u (monotone curve, valid for |u| <= 0.1425)
    pext = lambda t: 200.0 * (0.5 * (1 - np.cos(np.pi * t / 0.25)) if t < 0.25 else 1.0)
    acc = np.gradient(np.gradient(mid, T), T)
    ident = []
    for j, t in enumerate(times):
        if -um[-1] < -mid[j] or j < 2: continue
        pid = pel(mid[j]) + 0.05 * cin * acc[j] - pext(t)     # p_fluid,eq = p_el(u) + rho_s h c u'' - p_ext   (Pa, uniform-equivalent, + pushes the wall up)
        r = rows[j]; pa = None
        p = lambda i, jj, tt=t: -K * mr[tt][F(i, jj) + 1][2]
        p0 = np.array([p(i, nr) for i in range(i0, i1 + 1)]); p1 = np.array([p(i, nr - 1) for i in range(i0, i1 + 1)])
        pf1 = 0.5 * (p0[:-1] + p0[1:]); pf2 = 0.25 * (p0[:-1] + p0[1:] + p1[:-1] + p1[1:])
        ident.append(dict(t=t, p_identified=pid, p_F1_equiv=float(np.sum(pf1 * psim * dx0) / Ipsi), p_F2_equiv=float(np.sum(pf2 * psim * dx0) / Ipsi), p_odd_equiv=float(np.sum(0.5 * (p1[:-1] + p1[1:]) * psim * dx0) / Ipsi)))
    res["freebody_ident"] = ident; res["psi_Ipsi"] = Ipsi; res["psi_cin"] = cin
# mass/volume conservation: integral of (Qin - Qout) dt vs the channel area change (trapezoid in time)
t_ = np.array([0.0] + [r["t"] for r in rows]); dq = np.array([0.0] + [r["Qin"] - r["Qout"] for r in rows]); cum = np.concatenate([[0.0], np.cumsum(0.5 * (dq[1:] + dq[:-1]) * np.diff(t_))])
vol = np.array([0.0] + [r["vol"] for r in rows]); res["volume_residual_max_over_peak"] = float(np.max(np.abs(cum - vol)) / max(np.max(np.abs(vol)), 1e-12))   # d(area change)/dt = Qin - Qout
# Newton
logf = [f for f in os.listdir(rd) if f.endswith(".log") and not f.startswith("run")][0]; txt = open(rd + "/" + logf).read()
its = [int(m.group(1)) for m in re.finditer(r"number of iterations\s*:\s*(\d+)", txt)]
res["newton_per_step"] = its; res["newton_mean"] = float(np.mean(its)) if its else None; res["newton_max"] = max(its) if its else None
res["normal_termination"] = "N O R M A L   T E R M I N A T I O N" in txt; res["negative_jacobian"] = ("negative jacobian" in txt.lower()); res["error_termination"] = ("E R R O R   T E R M I N A T I O N" in txt)
json.dump(res, open(OUT, "w"), indent=1, default=float)
for k, v in res.items():
    if k not in ("rows", "newton_per_step", "extrema"): print(f"{k:34s} {v}")
print("extrema", [(a, round(b, 4), round(c, 4)) for a, b, c in res["extrema"][:6]]); print("newton", its)
