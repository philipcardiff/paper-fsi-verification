#!/usr/bin/env python3
import os, os as _os; HERE = _os.path.dirname(_os.path.abspath(__file__)); PKG = _os.path.dirname(HERE)   # package root; run directories of retained HPC outputs: env CAVITY_RUNS (default runs/cavityNext)
"""B4: attribution of the near-floor row-parity pressure in FEBio's own fixed-plate solution C1b (t = 300 s, Mon_all) to the forces that the pressure must balance.
The discrete momentum equation (validated in b4_residual_check.py to 9e-7) is   G_f p_free = -(A w + conv(w))_free - G_fixed p_fixed ,  G has full column rank (corner nodes pinned), so p is a LINEAR function of the force terms:
p = p[viscous] + p[convective] + p[boundary-fixed values].  Each part is computed by least squares with the SAME operator, then the row-parity D = mean(p_row1 - p_row0) is evaluated on each.  usage: OUT.json"""
import sys, json, numpy as np, scipy.linalg as sl
sys.path.insert(0, HERE); import b4_operator as B
sys.path.insert(0, HERE); import cavity_xplt_audit as aud
g = B.geometry("cavity", 1); X = g["X"]; n = len(X); K = 1e6
rd = os.environ.get("CAVITY_RUNS", "runs/cavityNext") + "/C1b_fluid_fixedplate_f1_pin"; Mo = aud.read_mon_block(rd + "/Mon_all.txt", 300.0)
w = np.zeros(2 * n); p = np.zeros(n)
for i, v in Mo.items(): k = i - 1; w[2 * k], w[2 * k + 1], p[k] = v[2], v[3], -K * v[4]
A0, G, M = B.assemble(X, g["quads"], c=0.0)
exec(open(""+HERE+"/b3_reaction_force_C1b.py").read().split("def conv(wv):")[1].join(["def conv(wv):", ""]).split("Rv, Rp, Rc")[0]) if False else None
def conv(wv):
    r = np.zeros(2 * n)
    for q in g["quads"]:
        x = X[list(q)]; wl = wv[np.array([[2 * i, 2 * i + 1] for i in q]).ravel()].reshape(4, 2)
        for (rr, ss) in B.GP:
            dr, ds = B.dN(rr, ss); J = np.array([[dr @ x[:, 0], dr @ x[:, 1]], [ds @ x[:, 0], ds @ x[:, 1]]]); dJ = np.linalg.det(J); grad = (np.linalg.inv(J) @ np.vstack([dr, ds])).T; H = B.N(rr, ss)
            a = (wl.T @ grad) @ (H @ wl)
            for aa in range(4): r[2 * q[aa]:2 * q[aa] + 2] += B.RHO * H[aa] * a * dJ
    return r
wfree, pfree, wfix, pfix = B.bcs(g, 1.15); pfix_l = sorted(pfix); Gf = G[wfree][:, pfree].toarray(); Gx = G[wfree][:, pfix_l].toarray()
fv = -(A0 @ w)[wfree]; fc = -conv(w)[wfree]; fb = -Gx @ p[pfix_l]
sol = {}
for nm, f in (("viscous", fv), ("convective", fc), ("fixed_boundary_values", fb)):
    x, *_ = np.linalg.lstsq(Gf, f, rcond=1e-10); sol[nm] = x
tot = sum(sol.values()); err = np.linalg.norm(tot - p[pfree]) / np.linalg.norm(p[pfree])
chains = B.rows_above_floor(g); ifs = sorted(g["iface"], key=lambda i: X[i, 0]); sel = [a for a in ifs if 2.4 < X[a, 0] < 5.6]
pos = {nd: k for k, nd in enumerate(pfree)}
def D(vec):
    return float(np.mean([vec[pos[chains[a][1]]] - vec[pos[chains[a][0]]] for a in sel if chains[a][0] in pos and chains[a][1] in pos]))
res = dict(reconstruction_rel_error_of_p=float(err), D_total_field=D(p[pfree]), D_from_sum=D(tot), parts={nm: dict(D=D(v), mean_p_row0_interior=float(np.mean([v[pos[a]] for a in sel if a in pos])), norm_force=float(np.linalg.norm(f))) for (nm, v), f in zip(sol.items(), (fv, fc, fb))},
           cond_note="lstsq rcond 1e-10; G full column rank iff the 4 corner pressures are pinned (b4_operator spectra)", n_free_p=len(pfree), n_free_w=len(wfree))
print(json.dumps(res, indent=1)); json.dump(res, open(sys.argv[1], "w"), indent=1)
