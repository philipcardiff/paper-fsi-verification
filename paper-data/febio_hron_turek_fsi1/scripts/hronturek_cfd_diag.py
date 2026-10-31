#!/usr/bin/env python3
"""Stage C diagnostics from EXISTING output only (no job): checkerboard/stripe structure, element-mean control-volume balance F3b (diagnostic, not a pre-registered gate),
flag transverse load from nodal vs element-mean pressure, and an ANALYTICAL (Euler-Bernoulli plane-strain) tip-deflection estimate from those loads.  usage: hronturek_cfd_diag.py RUN_DIR"""
import sys, os, json, math
import numpy as np
RUN = sys.argv[1]
sys.argv = ['x', RUN]
exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "hronturek_cfd_analyze.py")).read().split("# ---- time history")[0])
P = -K * A[-1, :, 2]; W = A[-1, :, :2]
out = {}
# 1. stripe/checkerboard amplitude along vertical grid lines (nodes sorted by y): a = rms of second difference/2 along the line, and fit of (-1)^j
lines = {}
for i in range(n2): lines.setdefault(round(XY[i, 0], 7), []).append(i)
res = []
for x, ids in sorted(lines.items()):
    if len(ids) < 40: continue
    ids = sorted(ids, key=lambda i: XY[i, 1]); p = P[ids]
    d = [(p[k] - 0.5 * (p[k - 1] + p[k + 1])) / 2 for k in range(1, len(p) - 1)]
    sgn = np.array([(-1) ** k for k in range(len(p))]); amp = float(np.dot(p - np.convolve(p, np.ones(5) / 5, mode="same"), sgn) / len(p))
    res.append((x, len(ids), float(np.sqrt(np.mean(np.square(d)))), float(p.mean()), float(np.ptp(p[::2]) if len(p) > 3 else 0)))
out["vertical_lines_x_n_a_t_pmean"] = [[round(a, 3), b, round(c, 3), round(d, 3)] for a, b, c, d, e in res]
# 2. element-mean pressure on outer boundaries -> F3b
def elem_p(qi): return P[Q[qi]].mean()
F = np.zeros(2); tops = np.zeros(2); bots = np.zeros(2)
for (e, qi, k) in bed:
    i, j = e; kd = kind(i, j)
    if kd == "obst": continue
    d = XY[j] - XY[i]; ln = math.hypot(*d); n = np.array([d[1], -d[0]]) / ln
    Fp = elem_p(qi) * ln * (-n); Fv = np.zeros(2); Fm = np.zeros(2)
    for s in G2:
        D, _ = gradD(qi, *EDGE_PAR[k](s), W); Fv += 2 * MU * (D @ n) * ln / 2
        t = 0.5 * (s + 1); u = W[i] * (1 - t) + W[j] * t; Fm += RHO * u * (u @ n) * ln / 2
    if kd in ("inlet", "outlet"): F -= Fm
    F += Fp + Fv
out["F3b_element_mean_pressure"] = F.tolist()
# 3. flag facets: transverse/axial force and equivalent tip deflection (EB, plane strain, E=1.4e6, h=0.02) from nodal (F1-type) and element-mean (F2-type) pressure; pressure only
E_, nu_, h = 1.4e6, 0.4, 0.02; EI = E_ / (1 - nu_**2) * h**3 / 12; L = 0.6 - XROOT
def flag_load(mode):
    xs, qs = [], []
    for (e, qi, k) in bed:
        i, j = e
        if kind(i, j) != "obst": continue
        xm, ym = 0.5 * (XY[i] + XY[j])
        if xm < XROOT + 1e-6 or xm > 0.6 - 1e-6 or not (abs(ym - 0.19) < 1e-9 or abs(ym - 0.21) < 1e-9): continue
        d = XY[j] - XY[i]; ln = abs(d[0]); n = np.array([d[1], -d[0]]) / math.hypot(*d)
        p = 0.5 * (P[i] + P[j]) if mode == "nodal" else elem_p(qi)
        xs.append(xm); qs.append(p * ln * n[1])      # force on the flag in y from this facet (pressure only; n_f points into the flag)
    xs = np.array(xs); qs = np.array(qs)
    s = xs - XROOT; fy = float(qs.sum())
    # tip deflection of a cantilever from point loads q_k at s_k: delta = sum q s^2 (3L - s)/(6 EI)
    delta = float(np.sum(qs * s**2 * (3 * L - s) / (6 * EI)))
    return fy, delta, float(np.sum(qs[np.array([True if abs(0)==0 else False for _ in qs])]))
for mode in ("nodal", "element"):
    fy, delta, _ = flag_load(mode); out[f"flag_pressure_Fy_{mode}"] = fy; out[f"flag_EB_tip_uy_{mode}_m"] = delta
out["reference_uy_FSI1_m"] = 8.208773e-4; out["note"] = "AE: Euler-Bernoulli plane-strain cantilever, pressure-only loads from the rigid-flag field; no FEBio run"
json.dump(out, open(RUN + "/stageC_diag.json", "w"), indent=1)
print(json.dumps(out, indent=1))
