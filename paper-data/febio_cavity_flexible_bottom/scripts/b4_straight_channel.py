#!/usr/bin/env python3
import os, os as _os; HERE = _os.path.dirname(_os.path.abspath(__file__)); PKG = _os.path.dirname(HERE)   # package root; run directories of retained HPC outputs: env CAVITY_RUNS (default runs/cavityNext)
"""B4: straight-channel (B0 run: 140x10 cells of 0.1, dt 0.5 s, ramp 10 s, K 1e6, mu 0.01) linear Stokes+inertia emulation with the FEBio element blocks and BCs of the model (backward Euler).
Compares the x-sawtooth that FEBio showed (0.01239 Pa frozen from t=20 to 150, appears in the ramp) and the y-alternation (0.00126 Pa) with the linear operator's own response.  usage: OUT.json"""
import sys, json, numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla
sys.path.insert(0, HERE); import b4_operator as B
def go(K=1e6, dt=0.5, tend=150.0, ramp=10.0, scheme="BE"):
    g = B.geometry("straight", 1); X = g["X"]; n = len(X); c = 1.0 / dt
    A0, G, M = B.assemble(X, g["quads"], c=0.0); At, _, _ = B.assemble(X, g["quads"], c=c); Mv = (At - A0) / (B.RHO * c)
    wfree, pfree, wfix, pfix = B.bcs(g, None); wfix_l = np.array(sorted(set(range(2 * n)) - set(wfree.tolist()))); prof = B.inlet_profile(g); bin_, _ = B.inlet_load(g)
    wp = np.zeros(2 * n)
    for i, ph in prof.items(): wp[2 * i] = ph
    GT = G.T.tocsr(); Gp = G[wfree][:, pfree]; GTw = GT[pfree][:, wfree]; Mpp = M[pfree][:, pfree]; pfix_l = sorted(pfix); Gpfix = G[wfree][:, pfix_l]; Mpfix = M[pfree][:, pfix_l]
    Kmat = sp.bmat([[At[wfree][:, wfree], Gp], [GTw, -(1.0 / (K * dt)) * Mpp]]).tocsc(); lu = sla.splu(Kmat); Awp = At[wfree][:, wfix_l]; GTwp = GT[pfree][:, wfix_l]; Mvff = Mv[wfree][:, wfree]
    inl_free = [i for i in g["inlet"] if i not in g["in_corner"]]; pidx = {nd: k for k, nd in enumerate(pfree)}
    w = np.zeros(2 * n); p = np.zeros(n); hist = []
    jy = np.round((X[:, 1] - 1) / 0.1).astype(int); ix = np.round(X[:, 0] / 0.1).astype(int)
    for s in range(1, int(round(tend / dt)) + 1):
        t = s * dt; lc = min(1.0, t / ramp); wpt = lc * wp; pf = np.zeros(len(pfix_l))
        for k, nd in enumerate(pfix_l):
            if nd in g["in_corner"]: pf[k] = np.mean(p[inl_free])
        rw = -Awp @ wpt[wfix_l] - Gpfix @ pf + (1.0 / dt) * (Mvff @ w[wfree])
        rp = -bin_[pfree] * lc - GTwp @ wpt[wfix_l] + (1.0 / (K * dt)) * (Mpp @ p[pfree]) + (1.0 / (K * dt)) * (Mpfix @ (pf - p[pfix_l]))
        x = lu.solve(np.concatenate([rw, rp])); w[wfree] = x[:len(wfree)]; w[wfix_l] = wpt[wfix_l]; p[pfree] = x[len(wfree):]; p[pfix_l] = pf
        if s % 10 == 0:
            m = (X[:, 0] > 2) & (X[:, 0] < 12); row = m & (jy == 5)   # y = 1.5
            ev = p[row & (ix % 2 == 0)].mean() - p[row & (ix % 2 == 1)].mean()
            col = (ix == 70); cj = jy[col]; py = p[col]; ya = py[cj % 2 == 0].mean() - py[cj % 2 == 1].mean()
            hist.append(dict(t=t, x_sawtooth_even_minus_odd=float(ev), y_alt_section_x7=float(ya), dpdx_pairavg=float(np.polyfit(X[row & (ix % 2 == 0), 0], p[row & (ix % 2 == 0)], 1)[0])))
    return hist
if __name__ == "__main__":
    out = {}
    for nm, kw in (("BE_dt0.5_K1e6", {}), ("BE_dt0.25", dict(dt=0.25)), ("BE_K1e7", dict(K=1e7))):
        h = go(**kw); out[nm] = h; print(nm, [(x["t"], round(x["x_sawtooth_even_minus_odd"], 5), round(x["y_alt_section_x7"], 5), round(x["dpdx_pairavg"], 5)) for x in (h[0], h[len(h) // 5], h[len(h) // 3], h[-1])]); sys.stdout.flush()
    json.dump(out, open(sys.argv[1], "w"), indent=1)
