#!/usr/bin/env python3
"""Stage C post-processing for the Hron-Turek CFD1 rigid-flag runs, exactly as pre-registered in HRON_TUREK_FSI1_STAGE_BC_RESULTS.md section 3.3-3.4.
usage: hronturek_cfd_analyze.py RUN_DIR [--all-steps]      (RUN_DIR contains model_meta.json, mesh2d.npz, Mon_all.txt, <tag>.log)
Writes RUN_DIR/stageC_analysis.json and prints a summary.  Forces per unit depth (z = 0 layer line integrals), n_f = outward fluid normal (into the obstacle)."""
import sys, os, json, re, math
import numpy as np

RUN = os.path.realpath(sys.argv[1]); ALLSTEPS = "--all-steps" in sys.argv
meta = json.load(open(RUN + "/model_meta.json")); M = np.load(RUN + "/mesh2d.npz")
K, MU, RHO, UBAR, H = meta["k"], 1.0, 1000.0, 0.2, 0.41
XY, Q = M["xy"], M["quads"]; n2 = int(M["n2"])
XROOT = 0.2 + math.sqrt(0.05**2 - 0.01**2)

def read_steps(path):
    out, cur = [], None
    for ln in open(path):
        if ln.startswith("*Time"): cur = [float(ln.split("=")[1]), []]; out.append(cur)
        elif ln.startswith("*") or not ln.strip(): continue
        else: cur[1].append(ln.split())
    T = np.array([c[0] for c in out]); A = np.array([[[float(x) for x in r[1:4]] for r in c[1]] for c in out])  # (steps, nodes, 3): wx wy ef
    ids = np.array([int(r[0]) for r in out[0][1]]); assert np.all(ids == np.arange(1, n2 + 1)), "node order"
    return T, A

T, A = read_steps(RUN + "/Mon_all.txt")
# ---- boundary edges
cnt = {}
for qi, q in enumerate(Q):
    for k in range(4):
        e = (int(q[k]), int(q[(k + 1) % 4])); cnt.setdefault(frozenset(e), []).append((e, qi, k))
bed = [v[0] for v in cnt.values() if len(v) == 1]            # (edge(i,j), quad, local edge index)
inl, outl, bot, top, obs = [set(M[k].tolist()) for k in ("inlet", "outlet", "bottom", "top", "obst")]
def kind(i, j):
    for nm, s in (("inlet", inl), ("outlet", outl), ("bottom", bot), ("top", top), ("obst", obs)):
        if i in s and j in s: return nm
    raise RuntimeError("edge kind")
G2 = np.array([-1 / math.sqrt(3), 1 / math.sqrt(3)])        # 2-point Gauss on [-1,1]
def shape(xi, eta):
    N = 0.25 * np.array([(1 - xi) * (1 - eta), (1 + xi) * (1 - eta), (1 + xi) * (1 + eta), (1 - xi) * (1 + eta)])
    dN = 0.25 * np.array([[-(1 - eta), -(1 - xi)], [(1 - eta), -(1 + xi)], [(1 + eta), (1 + xi)], [-(1 + eta), (1 - xi)]])
    return N, dN
def gradD(qi, xi, eta, W):
    """symmetric velocity gradient (2x2) at (xi,eta) of quad qi from nodal velocities W (n2,2)"""
    q = Q[qi]; X = XY[q]; N, dN = shape(xi, eta); J = X.T @ dN; g = dN @ np.linalg.inv(J)       # (4,2) dN/dx
    L = W[q].T @ g                                                                             # L[i,j] = du_i/dx_j
    return 0.5 * (L + L.T), N
# local edge k of the quad (corner k -> k+1) in parent coordinates
EDGE_PAR = [lambda s: (s, -1.0), lambda s: (1.0, s), lambda s: (-s, 1.0), lambda s: (-1.0, -s)]   # s in [-1,1] along the edge direction (corner k -> k+1)
def facet_data(i, j, qi, k, P, W, mode):
    """mode 'F1': pressure trapezoid + viscous at facet Gauss points; 'F2': element-mean p, D at element centre; returns pressure force vec, viscous force vec (force ON obstacle, n_f)"""
    d = XY[j] - XY[i]; ln = math.hypot(*d); n = np.array([d[1], -d[0]]) / ln
    q = Q[qi]
    if mode == "F1":
        Fp = 0.5 * (P[i] + P[j]) * ln * n
        Fv = np.zeros(2)
        for s in G2:
            D, _ = gradD(qi, *EDGE_PAR[k](s), W); Fv += -2 * MU * (D @ n) * ln / 2
    else:
        Fp = P[q].mean() * ln * n
        D, _ = gradD(qi, 0.0, 0.0, W); Fv = -2 * MU * (D @ n) * ln
    return Fp, Fv
def obstacle_forces(P, W):
    out = {}
    for mode in ("F1", "F2"):
        acc = {"cyl": [np.zeros(2), np.zeros(2)], "flag": [np.zeros(2), np.zeros(2)], "flag_top": [np.zeros(2), np.zeros(2)], "flag_bot": [np.zeros(2), np.zeros(2)], "flag_tip": [np.zeros(2), np.zeros(2)]}
        for (e, qi, k) in bed:
            i, j = e
            if kind(i, j) != "obst": continue
            xm, ym = 0.5 * (XY[i] + XY[j])
            if (xm - 0.2) ** 2 + (ym - 0.2) ** 2 < 0.0505 ** 2 and xm < XROOT + 1e-6 and not (abs(ym - 0.19) < 1e-9 or abs(ym - 0.21) < 1e-9): grp = "cyl"
            elif (xm - 0.2) ** 2 + (ym - 0.2) ** 2 < 0.0505 ** 2 and xm <= XROOT + 1e-4: grp = "cyl"
            else: grp = "flag"
            Fp, Fv = facet_data(i, j, qi, k, P, W, mode)
            acc[grp][0] += Fp; acc[grp][1] += Fv
            if grp == "flag":
                sub = "flag_tip" if xm > 0.6 - 1e-6 else ("flag_top" if ym > 0.2 else "flag_bot")
                acc[sub][0] += Fp; acc[sub][1] += Fv
        out[mode] = acc
    return out
def f3(P, W):
    F = np.zeros(2); Qin = Qout = 0.0
    for (e, qi, k) in bed:
        i, j = e; kd = kind(i, j)
        if kd == "obst": continue
        d = XY[j] - XY[i]; ln = math.hypot(*d); n = np.array([d[1], -d[0]]) / ln
        Fp = 0.5 * (P[i] + P[j]) * ln * (-n)          # -p n
        Fv = np.zeros(2); Fm = np.zeros(2)
        for s in G2:
            D, _ = gradD(qi, *EDGE_PAR[k](s), W); Fv += 2 * MU * (D @ n) * ln / 2             # + 2 mu D n
            t = 0.5 * (s + 1); u = W[i] * (1 - t) + W[j] * t; Fm += RHO * u * (u @ n) * ln / 2
        if kd in ("inlet", "outlet"): F -= Fm
        F += (Fp + Fv)                                  # F_obs = int sigma n - int rho u (u.n) (outer boundaries, n outward); Fp = -p n .. so sigma n = Fp + Fv with Fp = -p n
        flux = 0.5 * ((W[i] + W[j]) @ n) * ln
        if kd == "inlet": Qin -= flux
        if kd == "outlet": Qout += flux
    return F, Qin, Qout
# vertical section x = 1.0
sec = sorted([i for i in range(n2) if abs(XY[i, 0] - 1.0) < 1e-9], key=lambda i: XY[i, 1])
def q_section(W): return float(sum(0.5 * (W[sec[a], 0] + W[sec[a + 1], 0]) * (XY[sec[a + 1], 1] - XY[sec[a], 1]) for a in range(len(sec) - 1)))
# ---- flag wall chains for parity
def wall_chain(ytarget, xlo=XROOT, xhi=0.6):
    idx = [i for i in M["obst"].tolist() if abs(XY[i, 1] - ytarget) < 1e-9 and xlo - 1e-9 <= XY[i, 0] <= xhi + 1e-9]
    return sorted(idx, key=lambda i: XY[i, 0])
adj = {}
for q in Q:
    for k in range(4): adj.setdefault(int(q[k]), set()).update([int(q[(k + 1) % 4]), int(q[(k - 1) % 4])])
def inward_next(i, prev):
    """next node along the normal chain: the neighbour across the quad that is not on the boundary line; found via quad connectivity"""
    # choose the quad containing i (and the boundary neighbours) and take the diagonal-opposite's edge neighbour: use quads of i on the fluid side
    best = None
    for qi in np.where((Q == i).any(axis=1))[0]:
        q = list(Q[qi]); k = q.index(i); nxt, prv, opp = q[(k + 1) % 4], q[(k - 1) % 4], q[(k + 2) % 4]
        for cand in (nxt, prv):
            if cand not in bnd_nodes_obst and cand != prev: best = cand if best is None else best
    return best
bnd_nodes_obst = set(M["obst"].tolist())
def normal_chain(i, n=3):
    """three nodes above wall node i along the structured normal direction (through quad opposite corners)"""
    chain = [i]; cur = i; prevq = None
    for _ in range(n):
        # quads containing cur on the fluid side not already used
        cands = [qi for qi in np.where((Q == cur).any(axis=1))[0] if qi != prevq]
        # the wall quad(s) for the first step: take the quad where cur is on a boundary edge (obstacle); later steps: the quad that contains cur and does not contain the previous chain node
        sel = None
        for qi in cands:
            q = list(Q[qi])
            if chain[-2:-1] and chain[-2] in q: continue
            sel = qi; break
        if sel is None: break
        q = list(Q[sel]); k = q.index(cur)
        # nodes adjacent to cur in this quad: next/prev; the one not in the wall set (for first step) / not previous
        nbrs = [q[(k + 1) % 4], q[(k - 1) % 4]]
        if len(chain) == 1: nbrs = [x for x in nbrs if x not in bnd_nodes_obst]
        else: nbrs = [x for x in nbrs if x != chain[-2]]
        if not nbrs: break
        nxt = nbrs[0]; chain.append(nxt); cur = nxt; prevq = sel
    return chain
def parity(P, wall):
    ch = wall_chain(wall)
    xs = XY[ch, 0]; p = P[ch]; skip = 3                                  # exclude 3 elements at the root and tip corners
    inner = list(range(1 + skip, len(ch) - 1 - skip))
    at = [(p[k] - 0.5 * (p[k - 1] + p[k + 1])) / 2 for k in inner]
    rows = {}
    for r in (1, 2):
        rr = []
        for k in inner:
            c = normal_chain(ch[k], 3)
            rr.append(None if len(c) <= r else c[r])
        rows[r] = rr
    atr = {}
    for r in (1, 2):
        # alternating amplitude along the row-r chain through nodes above inner wall nodes (needs neighbours k-1,k+1)
        v = []
        for a in range(1, len(inner) - 1):
            ids = [rows[r][a - 1], rows[r][a], rows[r][a + 1]]
            if None in ids: continue
            v.append((P[ids[1]] - 0.5 * (P[ids[0]] + P[ids[2]])) / 2)
        atr[r] = float(np.sqrt(np.mean(np.square(v)))) if v else None
    dn = []
    for k in inner:
        c = normal_chain(ch[k], 3)
        if len(c) < 4: continue
        d = [np.linalg.norm(XY[c[m]] - XY[c[0]]) for m in range(4)]
        # Lagrange quadratic extrapolation to d=0 from nodes 1,2,3
        L = 0.0
        for m in (1, 2, 3):
            w = 1.0
            for o in (1, 2, 3):
                if o != m: w *= (0 - d[o]) / (d[m] - d[o])
            L += w * P[c[m]]
        dn.append(P[ch[k]] - L)
    return dict(n_nodes=len(ch), a_t=float(np.sqrt(np.mean(np.square(at)))), a_t_row1=atr[1], a_t_row2=atr[2],
                D_n_mean=float(np.mean(dn)) if dn else None, D_n_rms=float(np.sqrt(np.mean(np.square(dn)))) if dn else None, p_mean=float(p[inner].mean()),
                x0=float(xs[inner[0]]), x1=float(xs[inner[-1]]))
def cross_flag(P):
    top_c, bot_c = wall_chain(0.21), wall_chain(0.19)
    xt = {round(XY[i, 0], 9): i for i in top_c}; xb = {round(XY[i, 0], 9): i for i in bot_c}
    common = sorted(set(xt) & set(xb)); dp = np.array([P[xt[x]] - P[xb[x]] for x in common]); xs = np.array(common)
    # transverse force on the flag from dp (pressure only): F_y = int (p_bot - p_top) dx (trapezoid)
    Fy = float(np.trapz(-dp, xs)); skip = 3
    inner = slice(1 + skip, len(common) - 1 - skip)
    return dict(n_common=len(common), dp_mean=float(dp[inner].mean()), dp_rms=float(np.sqrt(np.mean(dp[inner] ** 2))), Fy_pressure_flagTopBottom=Fy, dp_first=float(dp[0]), dp_last=float(dp[-1]))
def cyl_parity(P):
    # contour order of cylinder wall nodes (angle about the centre), nodes with x < XROOT region and r~0.05
    ids = [i for i in M["obst"].tolist() if abs(math.hypot(XY[i, 0] - 0.2, XY[i, 1] - 0.2) - 0.05) < 1e-7]
    ang = sorted(ids, key=lambda i: math.atan2(XY[i, 1] - 0.2, XY[i, 0] - 0.2)); p = P[ang]
    # exclude the flag attachment arc (angles near 0)
    th = np.array([math.atan2(XY[i, 1] - 0.2, XY[i, 0] - 0.2) for i in ang]); sel = np.where(np.abs(th) > 0.5)[0]
    v = [(p[k] - 0.5 * (p[k - 1] + p[k + 1])) / 2 for k in sel if 0 < k < len(ang) - 1 and (k - 1) in set(sel) and (k + 1) in set(sel)]
    return float(np.sqrt(np.mean(np.square(v))))
# ---- time history (F1 drag/lift at every step) and final-step analysis
def fields(s):
    W = A[s, :, :2]; P = -K * A[s, :, 2]; return P, W
hist = []
steps = range(len(T)) if ALLSTEPS else [len(T) - 1]
for s in (range(len(T)) if True else steps):
    P, W = fields(s)
    ob = obstacle_forces(P, W) if (ALLSTEPS or s >= int(0.75 * len(T)) or s == len(T) - 1) else None
    if ob is None: continue
    F1 = ob["F1"]["cyl"][0] + ob["F1"]["cyl"][1] + ob["F1"]["flag"][0] + ob["F1"]["flag"][1]
    hist.append((float(T[s]), float(F1[0]), float(F1[1])))
hist = np.array(hist)
# steady-state spread over the last 20 % of the run t in [0.8 t_end, t_end]
t_end = T[-1]; w = hist[hist[:, 0] >= 0.8 * t_end - 1e-9]
spread_drag = float((w[:, 1].max() - w[:, 1].min()) / abs(w[:, 1].mean())); spread_lift = float((w[:, 2].max() - w[:, 2].min()) / abs(w[:, 2].mean()))
P, W = fields(len(T) - 1)
ob = obstacle_forces(P, W)
def tot(m): a = ob[m]; return (a["cyl"][0] + a["cyl"][1] + a["flag"][0] + a["flag"][1])
F1t, F2t = tot("F1"), tot("F2"); F3t, Qin, Qout = f3(P, W)
res = dict(run=os.path.basename(RUN), n_cells=meta["n_elem"], K=K, t_final=float(t_end), n_steps_logged=len(T),
           F1=F1t.tolist(), F2=F2t.tolist(), F3=F3t.tolist(),
           parts={m: {g: {"pressure": ob[m][g][0].tolist(), "viscous": ob[m][g][1].tolist()} for g in ob[m]} for m in ob},
           drag_err_F1=float((F1t[0] - 14.2929) / 14.2929), drag_err_F2=float((F2t[0] - 14.2929) / 14.2929), drag_err_F3=float((F3t[0] - 14.2929) / 14.2929),
           lift_err_F1=float((F1t[1] - 1.11905) / 1.11905), lift_err_F2=float((F2t[1] - 1.11905) / 1.11905), lift_err_F3=float((F3t[1] - 1.11905) / 1.11905),
           consist_drag_F1F2=float((F1t[0] - F2t[0]) / F1t[0]), consist_drag_F1F3=float((F1t[0] - F3t[0]) / F1t[0]),
           consist_lift_F1F2=float((F1t[1] - F2t[1]) / abs(F1t[1])), consist_lift_F1F3=float((F1t[1] - F3t[1]) / abs(F1t[1])),
           Qin=Qin, Qout=Qout, mass_imbalance=float(Qout / Qin - 1), Qsec_x1_over_Qin=float(q_section(W) / Qin - 1), Qin_minus_0p082=float(Qin - UBAR * H),
           spread_drag=spread_drag, spread_lift=spread_lift,
           parity_top=parity(P, 0.21), parity_bottom=parity(P, 0.19), cross_flag=cross_flag(P), parity_cyl_at=cyl_parity(P),
           p_max=float(P.max()), p_min=float(P.min()), hist_tail=hist[-5:].tolist())
# pin-rule diagnostics: junction nodes vs wall range, drift
root = M["root"].tolist(); wall_ids = [i for i in M["obst"].tolist() if i not in root]
pw = P[wall_ids]; res["junction"] = dict(p_junction=[float(P[i]) for i in root], wall_range=[float(pw.min()), float(pw.max())],
    drift_last20pct=[float((lambda ser: ser.max() - ser.min())(np.array([-K * A[s, i, 2] for s in range(len(T)) if T[s] >= 0.8 * t_end - 1e-9])) / max(abs(pw).max(), 1e-30)) for i in root])
nb = []
for i in root:
    nbrs = [m for m in adj[i] if m in bnd_nodes_obst]; res["junction"].setdefault("neighbour_p", []).append([float(P[m]) for m in nbrs]); nb.append(float(np.mean([P[m] for m in nbrs])))
res["junction"]["pin_value_rule_mean_of_neighbours"] = nb
json.dump(res, open(RUN + "/stageC_analysis.json", "w"), indent=1)
np.savetxt(RUN + "/F1_history.txt", hist, header="t drag_F1 lift_F1")
def pr(k, v): print(f"{k:34s} {v}")
pr("run", res["run"]); pr("cells", res["n_cells"]); pr("steps logged / t_final", f"{len(T)} / {t_end}")
pr("F1 drag, lift", f"{F1t[0]:.5f} {F1t[1]:.5f}   err {res['drag_err_F1']*100:+.3f}% {res['lift_err_F1']*100:+.3f}%")
pr("F2 drag, lift", f"{F2t[0]:.5f} {F2t[1]:.5f}   err {res['drag_err_F2']*100:+.3f}% {res['lift_err_F2']*100:+.3f}%")
pr("F3 drag, lift", f"{F3t[0]:.5f} {F3t[1]:.5f}   err {res['drag_err_F3']*100:+.3f}% {res['lift_err_F3']*100:+.3f}%")
pr("consistency drag F1-F2, F1-F3", f"{res['consist_drag_F1F2']*100:+.3f}% {res['consist_drag_F1F3']*100:+.3f}%")
pr("consistency lift F1-F2, F1-F3", f"{res['consist_lift_F1F2']*100:+.3f}% {res['consist_lift_F1F3']*100:+.3f}%")
pr("mass imbalance out/in, x=1/in", f"{res['mass_imbalance']:.2e} {res['Qsec_x1_over_Qin']:.2e}   Qin-0.082 = {res['Qin_minus_0p082']:.2e}")
pr("steady spread drag, lift (last 20%)", f"{spread_drag:.2e} {spread_lift:.2e}")
for nm in ("parity_top", "parity_bottom"):
    p_ = res[nm]; pr(nm, f"a_t {p_['a_t']:.4f} Pa (row1 {p_['a_t_row1']}, row2 {p_['a_t_row2']}) D_n {p_['D_n_mean']} rms {p_['D_n_rms']} p_mean {p_['p_mean']:.3f} n={p_['n_nodes']}")
pr("cross-flag dp mean, rms", f"{res['cross_flag']['dp_mean']:.4f} {res['cross_flag']['dp_rms']:.4f}")
pr("cylinder a_t", f"{res['parity_cyl_at']:.4f} Pa"); pr("p range", f"{res['p_min']:.2f} .. {res['p_max']:.2f} Pa")
pr("junction", json.dumps(res["junction"]))
