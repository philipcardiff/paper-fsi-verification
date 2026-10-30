#!/usr/bin/env python3
import os, os as _os; HERE = _os.path.dirname(_os.path.abspath(__file__)); PKG = _os.path.dirname(HERE)   # package root; run directories of retained HPC outputs: env CAVITY_RUNS (default runs/cavityNext)
"""B7 (matrix check only): Dohrmann-Bochev polynomial pressure projection (PPP) on the validated element operator: pressure block -(eps M + gamma (M - M_Pi)) with M_Pi = element-wise (int N_i)(int N_j)/V, gamma = alpha/mu.
Effects on (1) parity/smooth lowest eigenvalues of S + gamma (M - M_Pi) (straight channel h = 0.2, 0.1), (2) the cavity M1 corner DOFs WITHOUT pins, (3) global conservation: column sums of (M - M_Pi) vanish.  usage: OUT.json"""
import sys, json, numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla, scipy.linalg as sl
sys.path.insert(0, HERE); import b4_operator as B
import make_cavity_fluid_feb as GEO
def mpi(X, quads, n):
    rows, cols, vals = [], [], []
    for q in quads:
        x = X[list(q)]; V = 0; I = np.zeros(4)
        for (r, s) in B.GP:
            dr, ds = B.dN(r, s); dJ = np.linalg.det(np.array([[dr @ x[:, 0], dr @ x[:, 1]], [ds @ x[:, 0], ds @ x[:, 1]]])); V += dJ; I += B.N(r, s) * dJ
        for a in range(4):
            for b in range(4): rows.append(q[a]); cols.append(q[b]); vals.append(I[a] * I[b] / V)
    return sp.csr_matrix((vals, (rows, cols)), shape=(n, n))
def classes(g, pfree, V, ev, k=40):
    X, quads = g["X"], g["quads"]; fr = set(pfree.tolist()); out = {}
    for kk in range(k):
        full = np.zeros(len(X)); full[pfree] = V[:, kk]; cx = cy = tx = ty = 0
        for q in quads:
            for a, b in ((q[0], q[1]), (q[1], q[2]), (q[2], q[3]), (q[3], q[0])):
                if a in fr and b in fr:
                    if abs(X[a, 1] - X[b, 1]) < 1e-9: tx += 1; cx += full[a] * full[b] < 0
                    else: ty += 1; cy += full[a] * full[b] < 0
        ax, ay = cx / max(tx, 1), cy / max(ty, 1)
        c = "smooth" if ax < 0.2 and ay < 0.2 else ("y-parity" if ay > 0.8 and ax < 0.2 else ("x-parity" if ax > 0.8 and ay < 0.2 else ("checker" if ax > 0.8 and ay > 0.8 else "mixed")))
        out.setdefault(c, float(ev[kk]))
    return out
def straight(f):
    GEO.V = [(0, 1), (14, 1), (14, 2), (0, 2)]; GEO.BLOCKS = [((0, 1, 2, 3), int(140 * f), int(10 * f))]; GEO.ARCS = {}; B.GEO.V, B.GEO.BLOCKS, B.GEO.ARCS = GEO.V, GEO.BLOCKS, GEO.ARCS
    nodes, quads = GEO.mesh(1); X = np.array(nodes); eps = 1e-6; cnt = {}
    for q in quads:
        for e in ((q[0], q[1]), (q[1], q[2]), (q[2], q[3]), (q[3], q[0])): cnt.setdefault(frozenset(e), []).append(e)
    bn = {i for e in [v[0] for v in cnt.values() if len(v) == 1] for i in e}; inlet = sorted(i for i in bn if abs(X[i, 0]) < eps); outlet = sorted(i for i in bn if abs(X[i, 0] - 14) < eps)
    return dict(X=X, quads=quads, inlet=inlet, outlet=outlet, iface=[], walls=sorted(bn - set(inlet) - set(outlet)), in_corner=[i for i in inlet if abs(X[i, 1] - 1) < eps or abs(X[i, 1] - 2) < eps], out_corner=[i for i in outlet if abs(X[i, 1] - 1) < eps or abs(X[i, 1] - 2) < eps])
out = {}
for f in (0.5, 1):
    g = straight(f); X = g["X"]; n = len(X); A, G, M = B.assemble(X, g["quads"], c=0.0); wfree, pfree, _, _ = B.bcs(g, None); Gf = G[wfree][:, pfree]; Mp = M[pfree][:, pfree].toarray(); MP = mpi(X, g["quads"], n)
    S = Gf.T @ sla.splu(A[wfree][:, wfree].tocsc()).solve(Gf.toarray()); S = 0.5 * (S + S.T); D = (M - MP)[pfree][:, pfree].toarray(); res = {}
    for alpha in (0.0, 0.25, 1.0, 4.0):
        ev, V = sl.eigh(S + alpha / B.MU * D, Mp, subset_by_index=[0, 40]); res[f"alpha={alpha}"] = classes(g, pfree, V, ev); print(f"h={0.1/f:g} alpha={alpha}", {k: round(v, 5) for k, v in res[f'alpha={alpha}'].items()}); sys.stdout.flush()
    res["colsum_max_abs(M-MPi)_all_nodes"] = float(abs((M - MP).sum(0)).max()); out[f"straight_h={0.1/f}"] = res
B.GEO.V, B.GEO.BLOCKS, B.GEO.ARCS = B._ORIG; g = B.geometry("cavity", 1); X = g["X"]; n = len(X); A, G, M = B.assemble(X, g["quads"], c=0.0); wfree, pfree, _, _ = B.bcs(g, None); Gf = G[wfree][:, pfree]; Mp = M[pfree][:, pfree].toarray(); MP = mpi(X, g["quads"], n)
S = Gf.T @ sla.splu(A[wfree][:, wfree].tocsc()).solve(Gf.toarray()); S = 0.5 * (S + S.T); D = (M - MP)[pfree][:, pfree].toarray(); cav = {}
for alpha in (0.0, 1.0):
    ev, V = sl.eigh(S + alpha / B.MU * D, Mp, subset_by_index=[0, 6]); cav[f"alpha={alpha}"] = dict(lowest7=[float(x) for x in ev], classes=classes(g, pfree, V[:, :], ev, 7))
    print("cavity M1 NO PINS", alpha, [f"{x:.3g}" for x in ev])
out["cavity_M1_nopin"] = cav; json.dump(out, open(sys.argv[1], "w"), indent=1)
