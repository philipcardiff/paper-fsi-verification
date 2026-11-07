#!/usr/bin/env python3
import os, os as _os; HERE = _os.path.dirname(_os.path.abspath(__file__)); PKG = _os.path.dirname(HERE)   # package root; run directories of retained HPC outputs: env CAVITY_RUNS (default runs/cavityNext)
"""B4: spatial structure and h-scaling of the discrete steady x-sawtooth in the straight channel (linear operator, steady limit via one big backward-Euler step ensemble). usage: OUT.json"""
import sys, json, numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla
sys.path.insert(0, HERE); import b4_operator as B
def steady(f, K=1e6, dt=5.0, nsteps=60, nbc="rim"):
    g = B.geometry("straight", f); X = g["X"]; n = len(X); c = 1.0 / dt
    A0, G, M = B.assemble(X, g["quads"], c=0.0); At, _, _ = B.assemble(X, g["quads"], c=c); Mv = (At - A0) / (B.RHO * c)
    wfree, pfree, wfix, pfix = B.bcs(g, None); wfix_l = np.array(sorted(set(range(2 * n)) - set(wfree.tolist()))); prof = B.inlet_profile(g); bin_, _ = B.inlet_load(g)
    wp = np.zeros(2 * n)
    for i, ph in prof.items(): wp[2 * i] = ph
    GT = G.T.tocsr(); Gp = G[wfree][:, pfree]; GTw = GT[pfree][:, wfree]; Mpp = M[pfree][:, pfree]; pfix_l = sorted(pfix); Gpfix = G[wfree][:, pfix_l]; Mpfix = M[pfree][:, pfix_l]
    lu = sla.splu(sp.bmat([[At[wfree][:, wfree], Gp], [GTw, -(1.0 / (K * dt)) * Mpp]]).tocsc()); Awp = At[wfree][:, wfix_l]; GTwp = GT[pfree][:, wfix_l]; Mvff = Mv[wfree][:, wfree]
    inl_free = [i for i in g["inlet"] if i not in g["in_corner"]]; w = np.zeros(2 * n); p = np.zeros(n)
    for s in range(nsteps):
        pf = np.array([np.mean(p[inl_free]) if nd in g["in_corner"] else 0.0 for nd in pfix_l])
        rw = -Awp @ wp[wfix_l] - Gpfix @ pf + (1.0 / dt) * (Mvff @ w[wfree]); rp = -bin_[pfree] - GTwp @ wp[wfix_l] + (1.0 / (K * dt)) * (Mpp @ p[pfree]) + (1.0 / (K * dt)) * (Mpfix @ (pf - p[pfix_l]))
        x = lu.solve(np.concatenate([rw, rp])); w[wfree] = x[:len(wfree)]; w[wfix_l] = wp[wfix_l]; p[pfree] = x[len(wfree):]; p[pfix_l] = pf
    return X, p, g
out = {}
for f in (1, 2):
    X, p, g = steady(f); h = 0.1 / f; ix = np.round(X[:, 0] / h).astype(int); jy = np.round((X[:, 1] - 1) / h).astype(int); jm = jy.max() // 2
    row = (jy == jm)
    pts = [];
    for xs in (0.5, 1, 2, 4, 7, 10, 12, 13, 13.5):
        i0 = int(round(xs / h)); a = p[row & (ix == i0)][0]; b = p[row & (ix == i0 + 1)][0]; pts.append((xs, float(a - b)))
    mid = (X[:, 0] > 2) & (X[:, 0] < 12); ev = p[row & mid & (ix % 2 == 0)].mean() - p[row & mid & (ix % 2 == 1)].mean()
    out[f"f{f}"] = dict(h=h, even_minus_odd_mid=float(ev), stations_even_minus_next_odd=pts, pressure_drop_per_cell=float(0.12121 * h))
    print(f, "h", h, "even-odd mid", ev, "per-cell drop", 0.12121 * h, pts); sys.stdout.flush()
json.dump(out, open(sys.argv[1], "w"), indent=1)
