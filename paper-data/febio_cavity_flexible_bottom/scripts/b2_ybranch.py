#!/usr/bin/env python3
import os, os as _os; HERE = _os.path.dirname(_os.path.abspath(__file__)); PKG = _os.path.dirname(HERE)   # package root; run directories of retained HPC outputs: env CAVITY_RUNS (default runs/cavityNext)
"""DIAGNOSTIC (not a correction): y-alternation of the nodal pressure in the cavity near the floor.  For every floor node (reference y=0) the nodes of the same column at y_ref = h_y, 2 h_y, ... are read (z = 0 layer), giving
the row-to-row branch difference D = p(row1) - p(row0) and the pair mean (p0+p1)/2.  F_y with the raw floor row (the reported number) and with the pair mean are both given; the second is NOT a result, it measures how much of the force is carried by the alias.
usage: b2_ybranch.py RUNDIR TIME OUT.json"""
import sys, os, json, importlib.util, numpy as np
spec = importlib.util.spec_from_file_location("aud", HERE + "/cavity_xplt_audit.py"); aud = importlib.util.module_from_spec(spec); spec.loader.exec_module(aud)
rd, tsel, outp = sys.argv[1], float(sys.argv[2]), sys.argv[3]
meta = json.load(open(os.path.join(rd, "model_meta.json"))); tag, K = meta["tag"], meta["k"]
o = aud.read_xplt(os.path.join(rd, tag + ".xplt")); import re; _f = open(os.path.join(rd, tag + ".feb"), encoding="latin-1").read(); ns = {m.group(1): np.array(m.group(2).split(","), int) for m in re.finditer(r'<NodeSet name="(\w+)">([^<]+)</NodeSet>', _f)}
S = o["states"]; k = [i for i, s in enumerate(S) if abs(s["time"] - tsel) < 1e-6][0]; V = S[k]["vars"]
X0 = o["X0"]; Xd = X0 + V["displacement"] if "displacement" in V else X0.copy(); P = -K * V["fluid dilatation"][:, 0]
z0 = np.abs(X0[:, 2]) < 1e-9; ids = ns.get("InterfaceZ0", ns.get("Mon_iface")) - 1; ids = ids[np.argsort(X0[ids, 0])]
fl = np.zeros(len(X0), bool); fl[np.unique(np.array([e[1:] for e in o["domains"][-1]["elems"]]))] = True
rows = {}   # per floor node: p of the first 6 rows above
for n in ids:
    m = np.where(z0 & fl & (np.abs(X0[:, 0] - X0[n, 0]) < 1e-7) & (X0[:, 1] >= -1e-9) & (X0[:, 1] <= 0.75 + 1e-9))[0]; m = m[np.argsort(X0[m, 1])]; rows[int(n)] = (X0[m, 1], P[m])
R = np.array([r[1][:6] for r in rows.values()]); Y = np.array([r[0][:6] for r in rows.values()])
x = X0[ids, 0]; Xi, Yi = Xd[ids, 0], Xd[ids, 1]
def F(Pn): dx, dy = np.diff(Xi), np.diff(Yi); pb = 0.5 * (Pn[1:] + Pn[:-1]); return float(-np.sum(pb * dx)), float(np.sum(pb * dy))
p0, p1, p2 = R[:, 0], R[:, 1], R[:, 2]
# interior floor nodes only (pins at both ends have a constrained value, excluded from the statistics of D)
mid = slice(2, -2)
res = {"tag": tag, "time": tsel, "h_y_first_row": float(Y[mid, 1].mean()), "n_floor_nodes": int(len(ids)),
       "D_row1_minus_row0_mean_interior": float((p1 - p0)[mid].mean()), "D_row1_minus_row0_min_max": [float((p1 - p0)[mid].min()), float((p1 - p0)[mid].max())],
       "D_row2_minus_row1_mean_interior": float((p2 - p1)[mid].mean()), "p_row0_mean": float(p0.mean()), "p_row1_mean": float(p1.mean()), "p_row2_mean": float(p2.mean()),
       "rows_mean_interior_first6": [float(R[mid, j].mean()) for j in range(R.shape[1])],
       "Fy_floor_row_raw": F(p0)[0], "Fy_pairmean_row0_row1_DIAGNOSTIC_ONLY": F(0.5 * (p0 + p1))[0], "Fy_pairmean_endpins_kept": None,
       "Fy_row1_only_DIAGNOSTIC_ONLY": F(p1)[0], "Fy_row0_row1_row2_121_DIAGNOSTIC_ONLY": F(0.25 * p0 + 0.5 * p1 + 0.25 * p2)[0]}
json.dump(dict(res, per_node={"x_ref": x.tolist(), "p_row0": p0.tolist(), "p_row1": p1.tolist(), "p_row2": p2.tolist(), "y_rows": Y.tolist()}), open(outp, "w"), indent=1)
print(json.dumps({k_: v for k_, v in res.items() if k_ != "Fy_pairmean_endpins_kept"}, indent=1))
