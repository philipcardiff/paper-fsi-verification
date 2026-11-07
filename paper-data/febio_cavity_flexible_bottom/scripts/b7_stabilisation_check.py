#!/usr/bin/env python3
import os, os as _os; HERE = _os.path.dirname(_os.path.abspath(__file__)); PKG = _os.path.dirname(HERE)   # package root; run directories of retained HPC outputs: env CAVITY_RUNS (default runs/cavityNext)
"""B7 (paper-and-matrix check only, no FEBio): effect of a Brezzi-Pitkaranta term tau L (L = int gradN_i.gradN_j, tau = alpha h^2/mu) added to the pressure block of the validated FEBio element operator
[A G; G^T -(eps M + tau L)] on the soft parity eigenvalues of the Schur operator S + tau L (S q = lam M q), straight channel h = 0.2, 0.1, and on the smooth mode.  usage: OUT.json"""
import sys, json, numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla, scipy.linalg as sl
sys.path.insert(0, HERE); import b4_operator as B
import make_cavity_fluid_feb as GEO
def lap(X, quads, n):
    rows, cols, vals = [], [], []
    for q in quads:
        x = X[list(q)]
        for (r, s) in B.GP:
            dr, ds = B.dN(r, s); J = np.array([[dr @ x[:, 0], dr @ x[:, 1]], [ds @ x[:, 0], ds @ x[:, 1]]]); dJ = np.linalg.det(J); gr = (np.linalg.inv(J) @ np.vstack([dr, ds])).T
            for a in range(4):
                for b in range(4): rows.append(q[a]); cols.append(q[b]); vals.append(gr[a] @ gr[b] * dJ)
    return sp.csr_matrix((vals, (rows, cols)), shape=(n, n))
out = {}
for f in (0.5, 1):
    GEO.V = [(0, 1), (14, 1), (14, 2), (0, 2)]; GEO.BLOCKS = [((0, 1, 2, 3), int(140 * f), int(10 * f))]; GEO.ARCS = {}; B.GEO.V, B.GEO.BLOCKS, B.GEO.ARCS = GEO.V, GEO.BLOCKS, GEO.ARCS
    nodes, quads = GEO.mesh(1); X = np.array(nodes); n = len(X); eps = 1e-6; cnt = {}
    for q in quads:
        for e in ((q[0], q[1]), (q[1], q[2]), (q[2], q[3]), (q[3], q[0])): cnt.setdefault(frozenset(e), []).append(e)
    bn = {i for e in [v[0] for v in cnt.values() if len(v) == 1] for i in e}; inlet = sorted(i for i in bn if abs(X[i, 0]) < eps); outlet = sorted(i for i in bn if abs(X[i, 0] - 14) < eps)
    g = dict(X=X, quads=quads, inlet=inlet, outlet=outlet, iface=[], walls=sorted(bn - set(inlet) - set(outlet)), in_corner=[i for i in inlet if abs(X[i, 1] - 1) < eps or abs(X[i, 1] - 2) < eps], out_corner=[i for i in outlet if abs(X[i, 1] - 1) < eps or abs(X[i, 1] - 2) < eps])
    A, G, M = B.assemble(X, quads, c=0.0); wfree, pfree, _, _ = B.bcs(g, None); Gf = G[wfree][:, pfree]; Mp = M[pfree][:, pfree].toarray(); L = lap(X, quads, n)[pfree][:, pfree].toarray()
    S = Gf.T @ sla.splu(A[wfree][:, wfree].tocsc()).solve(Gf.toarray()); S = 0.5 * (S + S.T); h = 0.1 / f; res = {}
    for alpha in (0.0, 1 / 48, 1 / 12, 1 / 6):
        tau = alpha * h * h / B.MU; ev, V = sl.eigh(S + tau * L, Mp, subset_by_index=[0, 40]); alt = []
        for k in range(40):
            full = np.zeros(n); full[pfree] = V[:, k]; fr = set(pfree.tolist()); cx = cy = tx = ty = 0
            for q in quads:
                for a, b in ((q[0], q[1]), (q[1], q[2]), (q[2], q[3]), (q[3], q[0])):
                    if a in fr and b in fr:
                        if abs(X[a, 1] - X[b, 1]) < 1e-9: tx += 1; cx += full[a] * full[b] < 0
                        else: ty += 1; cy += full[a] * full[b] < 0
            alt.append((cx / max(tx, 1), cy / max(ty, 1)))
        cls = lambda ax, ay: "smooth" if ax < 0.2 and ay < 0.2 else ("y-parity" if ay > 0.8 and ax < 0.2 else ("x-parity" if ax > 0.8 and ay < 0.2 else ("checker" if ax > 0.8 and ay > 0.8 else "mixed")))
        first = {}
        for k in range(40): first.setdefault(cls(*alt[k]), float(ev[k]))
        res[f"alpha={alpha:.4f}"] = dict(tau=tau, first_eig_by_class=first, lowest=float(ev[0]))
        print(f"h={h} alpha={alpha:.4f}", {k: round(v, 5) for k, v in first.items()}); sys.stdout.flush()
    out[f"h={h}"] = res
json.dump(out, open(sys.argv[1], "w"), indent=1)
