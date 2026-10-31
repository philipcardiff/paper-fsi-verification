#!/usr/bin/env python3
"""Pressure-parity diagnostics of a CFD1 run (final state), definitions fixed here and documented in docs/STAGE_C_CFD1.md.  Independent of hronturek_cfd_analyze.py (connectivity from the .feb).
Wall chains: obstacle nodes on the flag top (y = 0.21) and bottom (y = 0.19) faces, x in [x_root, 0.6], sorted by x; inner nodes = excluding the 4 nodes nearest each end.
 a_t   alternating along-wall amplitude  rms_k [ (p_k - (p_{k-1}+p_{k+1})/2) / 2 ]                      (a pure odd-even mode of amplitude a gives exactly a; Pa)
 D_row mean_k [ p_wall,k - p_row1,k ]   (row-1 node = the node across the adjacent quad; for a smooth field this is the wall-normal pressure variation over one cell, for a wall-row parity mode it is +-2a)
 pbar  two-row average (p_wall + p_row1)/2 per node (alternation between adjacent rows cancels); dp_raw = mean(p_top,wall - p_bottom,wall) at equal x; dp_2row = mean(pbar_top - pbar_bottom) at equal x
 line  amplitude of the node-by-node alternation along the vertical grid line x = 1.0 (rms second difference / 2) and the line-mean pressure.
usage: parity_diagnostics.py RUN_DIR OUT.json"""
import sys, os, re, json, math, glob
import numpy as np
RUN = os.path.realpath(sys.argv[1]); OUT = sys.argv[2]
feb = glob.glob(RUN + "/*.feb")[0]; txt = open(feb, encoding="latin-1").read()
K = float(re.search(r"<k>([^<]+)</k>", txt).group(1))
nodes = np.array([[float(v) for v in m.group(2).split(",")] for m in re.finditer(r'<node id="(\d+)">([^<]+)</node>', txt)])
els = np.array([[int(v) - 1 for v in m.group(1).split(",")] for m in re.finditer(r'<elem id="\d+">([^<]+)</elem>', txt)])
z0 = np.where(np.abs(nodes[:, 2]) < 1e-12)[0]; remap = -np.ones(len(nodes), int); remap[z0] = np.arange(len(z0)); XY = nodes[z0, :2]; Q = remap[els[:, :4]]
for ln in open(RUN + "/Mon_all.txt"):
    if ln.startswith("*Time"): rows = []
    elif not ln.startswith("*") and ln.strip(): rows.append(ln.split())
P = -K * np.array([float(r[3]) for r in rows])
XR = 0.2 + math.sqrt(0.05**2 - 0.01**2)
cnt = {}
for qi, q in enumerate(Q):
    for k in range(4): a, b = int(q[k]), int(q[(k + 1) % 4]); cnt.setdefault((min(a, b), max(a, b)), []).append((a, b, qi))
bed = [v[0] for v in cnt.values() if len(v) == 1]
above = {}
for a, b, qi in bed:
    q = list(Q[qi]); k = q.index(a); assert q[(k + 1) % 4] == b
    above.setdefault(a, set()).add(q[(k + 3) % 4]); above.setdefault(b, set()).add(q[(k + 2) % 4])
def chain(y):
    ids = [i for a, b, _ in bed for i in (a, b) if abs(XY[i, 1] - y) < 1e-9 and XR - 1e-9 <= XY[i, 0] <= 0.6 + 1e-9]
    return sorted(set(ids), key=lambda i: XY[i, 0])
res = {"run": os.path.basename(RUN), "definitions": __doc__}
pb = {}
for nm, y in (("top", 0.21), ("bottom", 0.19)):
    c = chain(y); p = P[c]; inner = list(range(4, len(c) - 4))
    at = [(p[k] - 0.5 * (p[k - 1] + p[k + 1])) / 2 for k in inner]
    r1 = []
    for k in inner:
        nb = [m for m in above.get(c[k], []) if m not in c]
        # the wall-normal neighbour: the single node of the adjacent quads that is not on the chain and is closest in x
        m = min(nb, key=lambda m: abs(XY[m, 0] - XY[c[k], 0])) if nb else None; r1.append(m)
    drow = [p[k] - P[m] for k, m in zip(inner, r1) if m is not None]
    pbar = {c[k]: 0.5 * (p[k] + P[m]) for k, m in zip(inner, r1) if m is not None}
    pb[nm] = dict(chain=c, pbar=pbar, p=p)
    res[nm] = dict(n_nodes=len(c), n_inner=len(inner), a_t=float(np.sqrt(np.mean(np.square(at)))), D_row_mean=float(np.mean(drow)), D_row_rms=float(np.sqrt(np.mean(np.square(drow)))),
                   p_wall_mean=float(p[inner].mean()), pbar_mean=float(np.mean(list(pbar.values()))))
xt = {round(XY[i, 0], 9): i for i in pb["top"]["chain"]}; xb = {round(XY[i, 0], 9): i for i in pb["bottom"]["chain"]}
com = sorted(set(xt) & set(xb)); com_in = com[4:-4]
res["dp_raw_mean_top_minus_bottom"] = float(np.mean([P[xt[x]] - P[xb[x]] for x in com_in]))
ok = [x for x in com_in if xt[x] in pb["top"]["pbar"] and xb[x] in pb["bottom"]["pbar"]]
res["dp_2row_mean_top_minus_bottom"] = float(np.mean([pb["top"]["pbar"][xt[x]] - pb["bottom"]["pbar"][xb[x]] for x in ok]))
res["dp_2row_Fy_flag_pressure_only"] = float(-np.trapz([pb["top"]["pbar"][xt[x]] - pb["bottom"]["pbar"][xb[x]] for x in ok], ok))
ids = sorted([i for i in range(len(XY)) if abs(XY[i, 0] - 1.0) < 1e-9], key=lambda i: XY[i, 1]); p = P[ids]
res["line_x1.0"] = dict(n=len(ids), amplitude=float(np.sqrt(np.mean([((p[k] - 0.5 * (p[k - 1] + p[k + 1])) / 2) ** 2 for k in range(1, len(p) - 1)]))), mean=float(p.mean()),
                        first_values=[round(float(v), 3) for v in p[:8]])
json.dump(res, open(OUT, "w"), indent=1); print(json.dumps({k: v for k, v in res.items() if k != "definitions"}, indent=1))
