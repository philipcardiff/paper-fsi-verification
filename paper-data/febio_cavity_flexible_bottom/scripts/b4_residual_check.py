#!/usr/bin/env python3
import os, os as _os; HERE = _os.path.dirname(_os.path.abspath(__file__)); PKG = _os.path.dirname(HERE)   # package root; run directories of retained HPC outputs: env CAVITY_RUNS (default runs/cavityNext)
"""B4 (A)(B): evaluate the FEBio v4.13 discrete steady equations (reassembled from the source blocks: viscous sv.gradN_i, strong-form pressure gradient N_i grad p (dp = -K), convective rho N_i (Lf w), continuity gradN_i.w with the inlet flux load)
on the FEBio fixed-plate fluid-only solution C1b (nodal w, ef from Mon_all.txt at t = 300 s) and report (i) the residual of FEBio's own converged field under my assembly (validates the assembly),
(ii) how visible the near-floor row-parity pressure component is to the discrete momentum equation (|G p_parity| vs |G p_smooth| vs |viscous|), (iii) the continuity residual.  usage: OUT.json"""
import sys, json, numpy as np, scipy.sparse as sp
sys.path.insert(0, HERE); import b4_operator as B
sys.path.insert(0, HERE); import cavity_xplt_audit as aud
R = PKG; rd = os.environ.get("CAVITY_RUNS", "runs/cavityNext") + "/C1b_fluid_fixedplate_f1_pin"; K = 1e6
g = B.geometry("cavity", 1); X = g["X"]; n = len(X); M = aud.read_mon_block(rd + "/Mon_all.txt", 300.0); assert len(M) == n
w = np.zeros(2 * n); ef = np.zeros(n)
for i, v in M.items():
    k = i - 1; w[2 * k], w[2 * k + 1], ef[k] = v[2], v[3], v[4]
p = -K * ef; A0, G, Mm = B.assemble(X, g["quads"], c=0.0)
# convective term  int rho N_i (grad w) w
def conv(wv):
    r = np.zeros(2 * n)
    for q in g["quads"]:
        x = X[list(q)]; wl = wv[np.array([[2 * i, 2 * i + 1] for i in q]).ravel()].reshape(4, 2)
        for (rr, ss) in B.GP:
            dr, ds = B.dN(rr, ss); J = np.array([[dr @ x[:, 0], dr @ x[:, 1]], [ds @ x[:, 0], ds @ x[:, 1]]]); dJ = np.linalg.det(J); grad = (np.linalg.inv(J) @ np.vstack([dr, ds])).T; H = B.N(rr, ss)
            Lf = wl.T @ grad; wg = H @ wl; a = Lf @ wg
            for aa in range(4): r[2 * q[aa]:2 * q[aa] + 2] += B.RHO * H[aa] * a * dJ
    return r
Rv = A0 @ w; Rp = G @ p; Rc = conv(w); Rtot = Rv + Rp + Rc
wfree, pfree, wfix, pfix = B.bcs(g, 1.15)
res = {}
nrm = lambda v: float(np.linalg.norm(v[wfree]))
res["momentum_residual_free_rows"] = dict(viscous=nrm(Rv), pressure_gradient=nrm(Rp), convective=nrm(Rc), total=nrm(Rtot), total_over_pressure_gradient=nrm(Rtot) / nrm(Rp), max_abs_total=float(np.abs(Rtot[wfree]).max()))
# continuity: G^T w - inlet load  (free p rows) ; ef-dot = 0 steady ; (gradef.w term O(1/K) neglected, check below)
b, prof = B.inlet_load(g); GTw = G.T @ w; Rcont = GTw - np.where(np.isin(np.arange(n), g["inlet"]), 0, 0)
# inlet flux load sign: continuity row i: int gradN_i.w + load_i = 0 with load = -int N_i phi (VN = -phi): residual = GTw + b; residual = GTw + load = GTw - b
Rcont = GTw + b
pf = np.array(sorted(set(pfree.tolist()))); res["continuity_residual_free_p_rows"] = dict(norm=float(np.linalg.norm(Rcont[pf])), max_abs=float(np.abs(Rcont[pf]).max()), norm_GTw=float(np.linalg.norm(GTw[pf])), note="O(compress) (1/K) terms ignored: steady, gradef.w/Jf ~1e-6")
# parity decomposition on the near-floor patch (rows 0..7 above the floor nodes, x in [2.4,5.6])
chains = B.rows_above_floor(g); pa = np.zeros(n); ps = p.copy(); amps = []; nodes_patch = []
for a, ch in chains.items():
    if 2.4 < X[a, 0] < 5.6 and len(ch) >= 7:
        pv = np.array([p[m] for m in ch[:6]]); jj = np.arange(6); Amat = np.vstack([np.ones(6), (-1.0) ** jj]).T; (s, al), *_ = np.linalg.lstsq(Amat, pv, rcond=None)
        amps.append(al)
        for j in range(8):
            if j < len(ch): pa[ch[j]] = al * (-1.0) ** j; nodes_patch.append(ch[j])
Gpa = G @ pa; Gp = G @ p
# momentum rows of rows 1..4 of the patch, interior
rows_eval = []
for a, ch in chains.items():
    if 2.8 < X[a, 0] < 5.2 and len(ch) >= 7:
        for j in range(1, 5): rows_eval += [2 * ch[j], 2 * ch[j] + 1]
rows_eval = np.array([r for r in rows_eval if r in set(wfree.tolist())])
e = lambda v: float(np.linalg.norm(v[rows_eval]))
res["parity_visibility_rows1to4_x2.8to5.2"] = dict(n_rows=int(len(rows_eval)), mean_parity_amp_a=float(np.mean(amps)), D_row1_minus_row0=float(-2 * np.mean(amps)),
          norm_G_p_total=e(Gp), norm_G_p_parity_part=e(Gpa), norm_viscous=e(Rv), norm_convective=e(Rc), norm_total_residual=e(Rtot),
          ratio_Gp_parity_over_Gp_total=e(Gpa) / e(Gp), ratio_Gp_parity_over_viscous=e(Gpa) / e(Rv),
          note="parity patch truncated at row 7 (edge effects excluded by evaluating rows 1-4 only); G p_parity at interior rows of a uniform patch is exactly 0 for an ideal 2-D checkerboard in y")
# pressure parity content on the whole solution: local 5-point roughness
nb = {}
for q in g["quads"]:
    for k in range(4): nb.setdefault(q[k], set()).update([q[(k + 1) % 4], q[(k - 1) % 4]])
rough = np.zeros(n)
for i in range(n):
    if len(nb[i]) == 4 and all(len(nb[j]) >= 3 for j in nb[i]): rough[i] = p[i] - np.mean([p[j] for j in nb[i]])
reg = {"inlet_x<1": (X[:, 0] < 1), "channel_above_cavity_y>=1_x2..6": (X[:, 1] >= 1) & (X[:, 0] > 2) & (X[:, 0] < 6), "cavity_y<0.75_x2..6": (X[:, 1] < 0.75) & (X[:, 0] > 2) & (X[:, 0] < 6), "downstream_x>6.5": (X[:, 0] > 6.5) & (X[:, 1] >= 1)}
res["pressure_roughness_p_minus_mean_of_4_edge_neighbours_Pa_rms"] = {k: float(np.sqrt(np.mean(rough[m & (rough != 0)] ** 2))) for k, m in reg.items()}
json.dump(res, open(sys.argv[1], "w"), indent=1); print(json.dumps(res, indent=1))
