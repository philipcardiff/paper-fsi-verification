#!/usr/bin/env python3
import os, os as _os; HERE = _os.path.dirname(_os.path.abspath(__file__)); PKG = _os.path.dirname(HERE)   # package root; run directories of retained HPC outputs: env CAVITY_RUNS (default runs/cavityNext)
"""B4: (1) h-scaling of the soft (parity) eigenvalues of the pressure Schur operator S = G^T A^-1 G (S q = lam M q) for the straight channel at h = 0.2, 0.1, 0.05 (steady A);
(2) modal content of FEBio's own C1b pressure (Mon_all t = 300) in the eigenbasis of the cavity-M1 Schur operator: how much of the field (and of its parity part) sits in the soft modes.  usage: OUT.json [fmax]"""
import sys, json, numpy as np, scipy.sparse.linalg as sla, scipy.linalg as sl
sys.path.insert(0, HERE); import b4_operator as B
sys.path.insert(0, HERE); import cavity_xplt_audit as aud
def schur(g, pins=None, c=0.0):
    X = g["X"]; A, G, M = B.assemble(X, g["quads"], c=c); wfree, pfree, wfix, pfix = B.bcs(g, pins); Gf = G[wfree][:, pfree]; Mp = M[pfree][:, pfree]
    lu = sla.splu(A[wfree][:, wfree].tocsc()); S = Gf.T @ lu.solve(Gf.toarray()); S = 0.5 * (S + S.T); return S, Mp.toarray(), pfree
out = {}; fmax = float(sys.argv[2]) if len(sys.argv) > 2 else 2
for f in (0.5, 1, 2):
    if f > fmax: continue
    import make_cavity_fluid_feb as GEO
    GEO.V = [(0, 1), (14, 1), (14, 2), (0, 2)]; GEO.BLOCKS = [((0, 1, 2, 3), int(140 * f), int(10 * f))]; GEO.ARCS = {}
    nodes, quads = GEO.mesh(1); B.GEO.V, B.GEO.BLOCKS, B.GEO.ARCS = GEO.V, GEO.BLOCKS, GEO.ARCS
    X = np.array(nodes); eps = 1e-6; cnt = {}
    for q in quads:
        for e in ((q[0], q[1]), (q[1], q[2]), (q[2], q[3]), (q[3], q[0])): cnt.setdefault(frozenset(e), []).append(e)
    bn = {i for e in [v[0] for v in cnt.values() if len(v) == 1] for i in e}
    inlet = sorted(i for i in bn if abs(X[i, 0]) < eps); outlet = sorted(i for i in bn if abs(X[i, 0] - 14) < eps); walls = sorted(bn - set(inlet) - set(outlet))
    g = dict(X=X, quads=quads, inlet=inlet, outlet=outlet, iface=[], walls=walls, in_corner=[i for i in inlet if abs(X[i, 1] - 1) < eps or abs(X[i, 1] - 2) < eps], out_corner=[i for i in outlet if abs(X[i, 1] - 1) < eps or abs(X[i, 1] - 2) < eps])
    S, Mp, pfree = schur(g); ev, V = sl.eigh(S, Mp, subset_by_index=[0, 40]); info = []
    for k in range(40):
        full = np.zeros(len(X)); full[pfree] = V[:, k]; fr = set(pfree.tolist()); cx = cy = tx = ty = 0
        for q in quads:
            for a, b in ((q[0], q[1]), (q[1], q[2]), (q[2], q[3]), (q[3], q[0])):
                if a in fr and b in fr:
                    if abs(X[a, 1] - X[b, 1]) < 1e-9: tx += 1; cx += full[a] * full[b] < 0
                    else: ty += 1; cy += full[a] * full[b] < 0
        info.append((float(ev[k]), cx / max(tx, 1), cy / max(ty, 1)))
    cls = lambda ax, ay: "smooth" if ax < 0.2 and ay < 0.2 else ("y-parity" if ay > 0.8 and ax < 0.2 else ("x-parity" if ax > 0.8 and ay < 0.2 else ("checker" if ax > 0.8 and ay > 0.8 else "other/mixed")))
    first = {}
    for lam, ax, ay in info:
        first.setdefault(cls(ax, ay), lam)
    out[f"h={0.1 / f:g}"] = dict(n_p=len(pfree), lowest_eig=[i[0] for i in info[:10]], first_eig_by_class=first, classes_lowest10=[cls(i[1], i[2]) for i in info[:10]], eps_CK=B.CC / B.K)
    print(f"h={0.1 / f:g}", first, [round(i[0], 6) for i in info[:6]], [cls(i[1], i[2]) for i in info[:6]]); sys.stdout.flush()
json.dump(out, open(sys.argv[1], "w"), indent=1)   # part (1) is package-only
if not os.path.exists(os.environ.get("CAVITY_RUNS", "runs/cavityNext") + "/C1b_fluid_fixedplate_f1_pin/Mon_all.txt"): print("part (2) skipped: needs the retained C1b Mon_all.txt (HPC output, see docs/REPRODUCIBILITY_AUDIT.md)"); sys.exit(0)
# (2) cavity modal content of FEBio's C1b pressure
B.GEO.V, B.GEO.BLOCKS, B.GEO.ARCS = B._ORIG
g = B.geometry("cavity", 1); X = g["X"]; n = len(X); S, Mp, pfree = schur(g, 1.15); ev, V = sl.eigh(S, Mp)
rd = os.environ.get("CAVITY_RUNS", "runs/cavityNext") + "/C1b_fluid_fixedplate_f1_pin"; Mo = aud.read_mon_block(rd + "/Mon_all.txt", 300.0)
p = np.zeros(n)
for i, v in Mo.items(): p[i - 1] = -1e6 * v[4]
pf = p[pfree]; coef = V.T @ (Mp @ pf); e2 = coef ** 2; tot = e2.sum()
fr = set(pfree.tolist()); alt = np.zeros(len(pfree))
for k in range(len(pfree)):
    full = np.zeros(n); full[pfree] = V[:, k]; c = t = 0
    for q in g["quads"]:
        for a, b in ((q[0], q[1]), (q[1], q[2]), (q[2], q[3]), (q[3], q[0])):
            if a in fr and b in fr: t += 1; c += full[a] * full[b] < 0
    alt[k] = c / max(t, 1)
band = lambda lo, hi: float(e2[(ev >= lo) & (ev < hi)].sum() / tot)
out["C1b_pressure_modal_content_cavityM1_schur"] = dict(total_M_norm2=float(tot), frac_lowest_eig_lt_0p01=band(0, 0.01), frac_0p01_to_0p1=band(0.01, 0.1), frac_0p1_to_1=band(0.1, 1), frac_gt_1=band(1, 1e9),
      frac_in_modes_edgealt_ge_0p45=float(e2[alt >= 0.45].sum() / tot), frac_smooth_modes_edgealt_lt_0p2=float(e2[alt < 0.2].sum() / tot), n_modes_edgealt_ge_0p45=int((alt >= 0.45).sum()),
      median_lambda_of_modes_edgealt_ge_0p45=float(np.median(ev[alt >= 0.45])), median_lambda_of_modes_edgealt_lt_0p2=float(np.median(ev[alt < 0.2])), n_modes=len(ev),
      lowest5_lambda=[float(x) for x in ev[:5]], lowest5_coef_over_total=[float(e2[k] / tot) for k in range(5)])
print(json.dumps(out["C1b_pressure_modal_content_cavityM1_schur"], indent=1))
json.dump(out, open(sys.argv[1], "w"), indent=1)
