#!/usr/bin/env python3
import os, os as _os; HERE = _os.path.dirname(_os.path.abspath(__file__)); PKG = _os.path.dirname(HERE)   # package root; run directories of retained HPC outputs: env CAVITY_RUNS (default runs/cavityNext)
"""B3 extra (F5): discrete reaction of the fluid momentum equations on the floor (wall) nodes, from FEBio's own fixed-plate solution C1b (Mon_all nodal w, ef at t=300), using the assembly validated in b4_residual_check.py.
The force the fluid exerts on the wall in the discrete (Galerkin) sense is minus the momentum-equation residual contributions at the Dirichlet wall nodes: F5 = - sum_i (A w + G p + conv)_i  (per unit depth).  Compared with F1 (nodal p, LoadVector-style), F2/F4 (pair mean).
Not a physical claim: it is what a discretely consistent coupling would transmit for THIS element; labelled a diagnostic.  usage: OUT.json"""
import sys, json, numpy as np
sys.path.insert(0, HERE); import b4_operator as B
sys.path.insert(0, HERE); import cavity_xplt_audit as aud
exec(open(""+HERE+"/b4_residual_check.py").read().split("Rv = A0 @ w")[0].split("g = B.geometry")[0])  # imports only
g = B.geometry("cavity", 1); X = g["X"]; n = len(X); K = 1e6
rd = os.environ.get("CAVITY_RUNS", "runs/cavityNext") + "/C1b_fluid_fixedplate_f1_pin"; Mo = aud.read_mon_block(rd + "/Mon_all.txt", 300.0)
w = np.zeros(2 * n); p = np.zeros(n)
for i, v in Mo.items(): k = i - 1; w[2 * k], w[2 * k + 1], p[k] = v[2], v[3], -K * v[4]
A0, G, M = B.assemble(X, g["quads"], c=0.0)
def conv(wv):
    r = np.zeros(2 * n)
    for q in g["quads"]:
        x = X[list(q)]; wl = wv[np.array([[2 * i, 2 * i + 1] for i in q]).ravel()].reshape(4, 2)
        for (rr, ss) in B.GP:
            dr, ds = B.dN(rr, ss); J = np.array([[dr @ x[:, 0], dr @ x[:, 1]], [ds @ x[:, 0], ds @ x[:, 1]]]); dJ = np.linalg.det(J); grad = (np.linalg.inv(J) @ np.vstack([dr, ds])).T; H = B.N(rr, ss)
            a = (wl.T @ grad) @ (H @ wl)
            for aa in range(4): r[2 * q[aa]:2 * q[aa] + 2] += B.RHO * H[aa] * a * dJ
    return r
Rv, Rp, Rc = A0 @ w, G @ p, conv(w); ifs = g["iface"]; xs = X[ifs, 0]; o = np.argsort(xs); ifs = [ifs[i] for i in o]; xs = X[ifs, 0]
sel = lambda R_: np.array([R_[2 * i:2 * i + 2] for i in ifs])
tot = {nm: (-sel(R_).sum(0)).tolist() for nm, R_ in (("viscous", Rv), ("pressure_gradient", Rp), ("convective", Rc), ("total", Rv + Rp + Rc))}
chains = B.rows_above_floor(g); p0 = np.array([p[a] for a in ifs]); p1 = np.array([p[chains[a][1]] for a in ifs])
F1 = -np.sum(0.5 * (p0[1:] + p0[:-1]) * np.diff(xs)); F2 = -np.sum(0.5 * (0.5 * (p0[1:] + p0[:-1]) + 0.5 * (p1[1:] + p1[:-1])) * np.diff(xs) * 1.0) * 1.0
F2 = -np.sum(0.5 * ((p0[1:] + p0[:-1]) / 2 + (p1[1:] + p1[:-1]) / 2) / 1.0 * np.diff(xs) * 1.0) / 1.0 * 1.0
F2 = -np.sum(np.diff(xs) * ((p0[1:] + p0[:-1] + p1[1:] + p1[:-1]) / 4))
sub = (xs > 2.0 + 1e-9) & (xs < 6.0 - 1e-9)
out = dict(F_on_wall_discrete_reaction_Fx_Fy_per_m=tot, F1_nodal_wallrow_Fy=float(F1), F2_pairmean_Fy=float(F2), D_mean_row1_minus_row0=float(np.mean((p1 - p0)[4:-4])), n_floor_nodes=len(ifs),
           note="F5 = -sum over floor nodes of (A w + G p + conv); fixed plate (C1b, mesh M1, pinned corners). Corner rows include the pinned-dilatation nodes. Per unit depth.")
print(json.dumps(out, indent=1)); json.dump(out, open(sys.argv[1], "w"), indent=1)
