#!/usr/bin/env python3
"""CC-F1 analysis (frozen-geometry fluid-only run of tools/make_cc_feb.py --mode fluid): mass balance, steady drift, wall-pressure forces F1/F2 (net vertical force and first-mode generalised load),
parity a_t / D_n, junction-pin rule diagnostics.  Definitions (declared in COLLAPSIBLE_CHANNEL_PILOT_PLAN.md):
 F1 = sum over wall facets of (p_a + p_b)/2 * (-dy, dx)  (facet-linear nodal wall pressure; applied-load equivalent of the nodal FSI traction);
 F2 = sum of the element-mean pressure (mean of the 4 nodal pressures of the adjacent fluid cell) * (-dy, dx);  (-dy, dx) = facet vector pointing out of the fluid (into the wall), force ON the wall = p * that vector.
 G = sum p * phi(x_ref) * dx_facet with phi = (1 - cos(2 pi (x - 5)/10))/2.   a_t = rms over interior wall nodes (3 facets excluded at each end) of (p_k - (p_{k-1} + p_{k+1})/2)/2.
 D_n = p_wall - (2 p_row1 - p_row2) (linear extrapolation of rows 1, 2 to the wall).  p = -K e_f (kinematic = Pa at rho = 1).
usage: cc_f1_analyze.py RUN_DIR [OUT.json]"""
import sys, os, json, numpy as np
rd = os.path.realpath(sys.argv[1]); OUT = sys.argv[2] if len(sys.argv) > 2 else rd + "/f1_result.json"
meta = json.load(open(rd + "/model_meta.json")); M = json.load(open(rd + "/mesh2d.json")); K = meta["k"]
nodes2 = np.array(M["nodes2"]); xs = np.array(M["xs"]); nr = M["nr"]; i0, i1 = M["i0w"], M["i1w"]; F = lambda i, j: i * (nr + 1) + j
def read(fn):
    out = {}; t = None
    for l in open(fn):
        if l.startswith("*Time"): t = float(l.split("=")[1]); out[t] = {}
        elif not l.startswith("*") and l.strip():
            f = l.split(); out[t][int(f[0])] = [float(v) for v in f[1:]]
    return out
mr = read(rd + "/Mon_rows.txt"); mi = read(rd + "/Mon_inlet.txt"); mo = read(rd + "/Mon_outlet.txt")
times = sorted(t for t in mr if t > 0)
nw = i1 - i0; xi = xs[i0:i1 + 1]
phi = lambda x: 0.5 * (1 - np.cos(2 * np.pi * (x - 5.0) / 10.0))
def pressure(t, i, j): return -K * mr[t][F(i, j) + 1][2]
def forces(t):
    p0 = np.array([pressure(t, i, nr) for i in range(i0, i1 + 1)]); p1 = np.array([pressure(t, i, nr - 1) for i in range(i0, i1 + 1)]); p2 = np.array([pressure(t, i, nr - 2) for i in range(i0, i1 + 1)])
    P = nodes2[[F(i, nr) for i in range(i0, i1 + 1)]]; dx = np.diff(P[:, 0]); dy = np.diff(P[:, 1])
    F1 = np.array([np.sum(0.5 * (p0[:-1] + p0[1:]) * (-dy)), np.sum(0.5 * (p0[:-1] + p0[1:]) * dx)])
    pe = 0.25 * (p0[:-1] + p0[1:] + p1[:-1] + p1[1:]); F2 = np.array([np.sum(pe * (-dy)), np.sum(pe * dx)])
    xm = 0.5 * (xi[:-1] + xi[1:]); G1 = np.sum(0.5 * (p0[:-1] + p0[1:]) * phi(xm) * dx); G2 = np.sum(pe * phi(xm) * dx)
    return p0, p1, p2, F1, F2, G1, G2
def parity(p, skip=3):
    k = np.arange(1 + skip, len(p) - 1 - skip); return float(np.sqrt(np.mean(((p[k] - 0.5 * (p[k - 1] + p[k + 1])) / 2) ** 2)))
def flux(m, ids_y):
    ids = [F(0, j) + 1 for j in range(nr + 1)] if m is mi else [F(len(xs) - 1, j) + 1 for j in range(nr + 1)]
    wx = np.array([m[t][i][0] for i in ids]); y = np.array([nodes2[F(0 if m is mi else len(xs) - 1, j)][1] for j in range(nr + 1)]); return float(np.trapz(wx, y))
hist = []
for t in times:
    p0, p1, p2, F1, F2, G1, G2 = forces(t)
    hist.append(dict(t=t, F1x=F1[0], F1y=F1[1], F2x=F2[0], F2y=F2[1], G1=G1, G2=G2, Qin=flux(mi, None), Qout=flux(mo, None), pj_lo=p0[0], pj_hi=p0[-1]))
t = times[-1]; p0, p1, p2, F1, F2, G1, G2 = forces(t)
n80 = int(round(0.8 * len(times))); h80 = hist[n80 - 1]; hl = hist[-1]
Dn = p0 - (2 * p1 - p2); inner = slice(4, len(p0) - 4)
V = float(p0[inner].max() - p0[inner].min()); pmaxwall = float(np.abs(p0).max())
nb = lambda k: p0[k]
res = dict(run=rd, K=K, t_last=t, n_steps=len(times), mesh=dict(nw=nw, n_fluid=meta.get("n_fluid_elem")),
    Qin=hl["Qin"], Qout=hl["Qout"], mass_balance=hl["Qout"] / hl["Qin"] - 1,
    F1y=F1[1], F2y=F2[1], F1x=F1[0], F2x=F2[0], rel_F1_F2_y=(F1[1] - F2[1]) / F2[1], G1=G1, G2=G2, rel_G=(G1 - G2) / G2,
    eps_abs_equiv=abs(F1[1] - F2[1]) / (420.0 * 10.0),
    steady_dF2y_last20=abs(hl["F2y"] - h80["F2y"]) / abs(hl["F2y"]), steady_dF1y_last20=abs(hl["F1y"] - h80["F1y"]) / abs(hl["F1y"]),
    junction_drift_lo=abs(hl["pj_lo"] - h80["pj_lo"]) / pmaxwall, junction_drift_hi=abs(hl["pj_hi"] - h80["pj_hi"]) / pmaxwall,
    junction_p=[float(p0[0]), float(p0[-1])], neighbour_range_lo=[float(p0[1:4].min()), float(p0[1:4].max())], neighbour_range_hi=[float(p0[-4:-1].min()), float(p0[-4:-1].max())],
    a_t_wall=parity(p0), a_t_row1=parity(p1), a_t_row2=parity(p2), a_t_over_V=parity(p0) / V, V=V, p_wall_max=float(p0.max()), p_wall_min=float(p0.min()),
    D_n_mean=float(Dn[inner].mean()), D_n_rms=float(np.sqrt(np.mean(Dn[inner] ** 2))), p_wall_mean=float(p0[inner].mean()),
    pin_value_lo=float(p0[1:3].mean()), pin_value_hi=float(p0[-3:-1].mean()))
res["pin_rule_b_lo"] = not (res["neighbour_range_lo"][0] <= p0[0] <= res["neighbour_range_lo"][1]); res["pin_rule_b_hi"] = not (res["neighbour_range_hi"][0] <= p0[-1] <= res["neighbour_range_hi"][1])
res["pin_rule_c"] = bool(max(res["junction_drift_lo"], res["junction_drift_hi"]) > 1e-3)
log = [f for f in os.listdir(rd) if f.endswith(".log") and not f.startswith("run")][0]; txt = open(rd + "/" + log).read()
import re
m = re.search(r"Average number of equilibrium iterations \.+ : ([\d.]+)", txt); res["newton_avg"] = float(m.group(1)) if m else None
res["normal_termination"] = "N O R M A L   T E R M I N A T I O N" in txt; res["negative_jacobian"] = ("negative jacobian" in txt.lower())
json.dump(res, open(OUT, "w"), indent=1, default=float)
json.dump(hist, open(OUT.replace(".json", "_hist.json"), "w"), indent=0)
for k, v in res.items(): print(f"{k:24s} {v}")
