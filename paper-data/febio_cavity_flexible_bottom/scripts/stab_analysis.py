#!/usr/bin/env python3
import os, os as _os; HERE = _os.path.dirname(_os.path.abspath(__file__)); PKG = _os.path.dirname(HERE)   # package root; run directories of retained HPC outputs: env CAVITY_RUNS (default runs/cavityNext)
"""Diagnostics for a C1b-type fixed-plate fluid-only run (stock or experimental PPP). Independent python re-assembly (cavity_agent/pressure_mode/b4_operator.py).
usage: stab_analysis.py RUN_DIR ALPHA OUT.json [T=300]   (ALPHA is the value used by the run; used ONLY to evaluate the stabilisation term in the python residual check)
Outputs: momentum residual (stock momentum equations; unchanged by the PPP term), continuity residual WITHOUT and WITH the PPP term, mass flux balance, row-parity D (Pa), x-parity, F1/F2, pair-mean (smooth) pressure,
discrete wall/inlet reactions, plus raw fields for pair-wise comparison (npz next to OUT)."""
import sys, os, json, numpy as np, scipy.sparse as sp
rd, ALPHA, outp = sys.argv[1], float(sys.argv[2]), sys.argv[3]; T = float(sys.argv[4]) if len(sys.argv) > 4 else 300.0
os.environ["RD_PPP"] = rd
src = open(""+HERE+"/b4_residual_check.py").read()
src = src.replace('rd = R + "/runs/cavityNext/C1b_fluid_fixedplate_f1_pin"', 'rd = __import__("os").environ["RD_PPP"]').replace("M = aud.read_mon_block(rd + \"/Mon_all.txt\", 300.0)", f"M = aud.read_mon_block(rd + \"/Mon_all.txt\", {T})")
src = src.replace('json.dump(res, open(sys.argv[1], "w"), indent=1); print(json.dumps(res, indent=1))', "")
sys.argv = [sys.argv[0], "/dev/null"]; ns = {"__name__": "ns"}; exec(compile(src, "b4_residual_check_adapted", "exec"), ns)
B, g, X, n, p, w, K, Rv, Rp, Rc, Rtot, G, A0, Mm, res = (ns[k] for k in ("B", "g", "X", "n", "p", "w", "K", "Rv", "Rp", "Rc", "Rtot", "G", "A0", "Mm", "res"))
wfree, pfree, wfix, pfix = ns["wfree"], ns["pfree"], ns["wfix"], ns["pfix"]; MU = B.MU
def mpi(X, quads, n):
    rows, cols, vals = [], [], []
    for q in quads:
        x = X[list(q)]; V = 0; I = np.zeros(4)
        for (r, s) in B.GP:
            dr, ds = B.dN(r, s); dJ = np.linalg.det(np.array([[dr @ x[:, 0], dr @ x[:, 1]], [ds @ x[:, 0], ds @ x[:, 1]]])); V += dJ; I += B.N(r, s) * dJ
        for a in range(4):
            for b in range(4): rows.append(q[a]); cols.append(q[b]); vals.append(I[a] * I[b] / V)
    return sp.csr_matrix((vals, (rows, cols)), shape=(n, n))
MP = mpi(X, g["quads"], n); S = -(ALPHA / MU) * ((Mm - MP) @ p)      # PPP term in the continuity rows (p variable, z-uniform plane-strain => per unit depth same factor as G^T w)
b, prof = B.inlet_load(g); GTw = G.T @ w; Rc0 = GTw + b; Rc1 = Rc0 + S; pf = np.array(sorted(set(pfree.tolist())))
res["continuity_check"] = dict(alpha_used=ALPHA, norm_GTw=float(np.linalg.norm(GTw[pf])), resid_without_PPP_term=float(np.linalg.norm(Rc0[pf])), resid_with_PPP_term=float(np.linalg.norm(Rc1[pf])), norm_PPP_term=float(np.linalg.norm(S[pf])),
    rel_resid_with_PPP=float(np.linalg.norm(Rc1[pf]) / np.linalg.norm(GTw[pf])), sum_PPP_term_all_nodes=float(S.sum()), sum_PPP_term_constrained_e_rows=float(S[[i for i in range(n) if i not in set(pf.tolist())]].sum()), sum_abs_PPP_term=float(np.abs(S).sum()))
# parity amplitudes and pair-mean fields along columns above the floor
chains = B.rows_above_floor(g); ifs = sorted(g["iface"], key=lambda i: X[i, 0]); xs = X[ifs, 0]
p0 = np.array([p[a] for a in ifs]); p1 = np.array([p[chains[a][1]] for a in ifs]); p2 = np.array([p[chains[a][2]] for a in ifs])
inn = (xs > 2.4 - 1e-9) & (xs < 5.6 + 1e-9)
res["parity"] = dict(D_row1_minus_row0_mean_Pa=float(np.mean((p1 - p0)[inn])), D_rms_Pa=float(np.sqrt(np.mean((p1 - p0)[inn] ** 2))), wallrow_mean_Pa=float(p0[inn].mean()), row1_mean_Pa=float(p1[inn].mean()),
    a_parity_mean_Pa=float(np.mean(((p0 - 2 * p1 + p2) / 4)[inn])), a_parity_rms_Pa=float(np.sqrt(np.mean(((p0 - 2 * p1 + p2) / 4)[inn] ** 2))))
d = p0[inn]; i0 = np.arange(len(d)); lap = np.array([(-1.0) ** k * (d[k - 1] - 2 * d[k] + d[k + 1]) / 4 for k in range(1, len(d) - 1)]); res["parity"]["xparity_wallrow_rms_Pa"] = float(np.sqrt(np.mean(lap ** 2))); res["parity"]["xparity_wallrow_signed_mean_Pa"] = float(lap.mean())
F1 = -np.sum(0.5 * (p0[1:] + p0[:-1]) * np.diff(xs)); F2 = -np.sum(np.diff(xs) * ((p0[1:] + p0[:-1] + p1[1:] + p1[:-1]) / 4))
res["forces_N_per_m_pressure_Fy"] = dict(F1_applied_nodal_wallrow=float(F1), F2_element_mean_pairrow=float(F2), F1_minus_F2=float(F1 - F2), rel_F1_minus_F2=float((F1 - F2) / abs(F2)))
# whole-domain pair-mean pressure (row pairs along every chain) for stock-vs-PPP comparison; save the raw fields
smooth = p.copy()
for a, ch in chains.items():
    for j in range(len(ch) - 1): pass
# discrete reactions: sum of momentum residual over fixed-w nodes (wall + inlet) and over free (should vanish) ; flux balance from Mon
fixedw = sorted(set(i for i in range(n) if i in wfix)); fx = np.array([[Rtot[2 * i], Rtot[2 * i + 1]] for i in fixedw]); fr = np.array([[Rtot[2 * i], Rtot[2 * i + 1]] for i in range(n) if i not in wfix])
iface_react = np.array([[Rtot[2 * i], Rtot[2 * i + 1]] for i in ifs]).sum(0)
res["momentum_balance"] = dict(sum_reaction_fixedw_nodes_xy=fx.sum(0).tolist(), sum_residual_freew_nodes_xy=fr.sum(0).tolist(), sum_all_nodes_xy=(fx.sum(0) + fr.sum(0)).tolist(), floor_reaction_xy=(-iface_react).tolist(), momentum_flux_scale_rho_U2_H=1.0)
# fluxes from Mon columns (trapezoid) at x = 1, 13, 14 and inlet
def flux(xt):
    ids = [i for i in range(n) if abs(X[i, 0] - xt) < 1e-6 and 1 - 1e-9 <= X[i, 1] <= 2 + 1e-9]; ids.sort(key=lambda i: X[i, 1]); yy = X[ids, 1]; ww = np.array([w[2 * i] for i in ids]); return float(np.trapz(ww, yy))
res["flux"] = {k: flux(v) for k, v in (("x0", 0.0), ("x1", X[np.argmin(abs(X[:, 0] - 1.0)), 0]), ("x13", X[np.argmin(abs(X[:, 0] - 13.0)), 0]), ("x14", 14.0))}
res["flux"]["Qout_over_Qin_minus_1"] = res["flux"]["x14"] / res["flux"]["x0"] - 1
json.dump(res, open(outp, "w"), indent=1)
np.savez(outp.replace(".json", ".npz"), X=X, p=p, w=w, p0=p0, p1=p1, p2=p2, xs=xs)
print(json.dumps({k: res[k] for k in ("continuity_check", "parity", "forces_N_per_m_pressure_Fy", "momentum_balance", "flux", "momentum_residual_free_rows")}, indent=1))
