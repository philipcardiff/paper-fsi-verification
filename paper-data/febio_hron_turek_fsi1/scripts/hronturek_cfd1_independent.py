#!/usr/bin/env python3
"""INDEPENDENT second implementation of the Stage C force/parity reductions for the Hron-Turek CFD1 rigid-flag runs (written separately from tools/hronturek_cfd_analyze.py:
geometry and connectivity are re-read from the .feb input (not mesh2d.npz), the boundary is re-derived from element topology, all kernels are vectorised and the velocity
gradients use explicit closed-form bilinear shape-function derivatives).  Adds a variationally consistent nodal-reaction force F4 and a term-by-term budget of the F2/F3b difference.
usage: hronturek_cfd1_independent.py RUN_DIR [OUT.json]       (RUN_DIR: <tag>.feb, Mon_all.txt, <tag>.xplt optional)
Per-unit-depth forces (slab 2-D, z = 0 layer).  Pressure p = -K e_f (nodal e_f from the log).  Convention: force ON the obstacle = - int sigma n_f dS, n_f out of the fluid."""
import sys, os, re, json, math, glob
import numpy as np

RUN = os.path.realpath(sys.argv[1]); OUT = sys.argv[2] if len(sys.argv) > 2 else os.path.join(RUN, "independent_analysis.json")
feb = glob.glob(RUN + "/*.feb")[0]; tag = os.path.basename(feb)[:-4]
txt = open(feb, encoding="latin-1").read()
K = float(re.search(r"<k>([^<]+)</k>", txt).group(1)); MU = float(re.search(r"<mu>([^<]+)</mu>", txt).group(1)); RHO = float(re.search(r"<density>([^<]+)</density>", txt).group(1))
nodes = np.array([[float(v) for v in m.group(2).split(",")] for m in re.finditer(r'<node id="(\d+)">([^<]+)</node>', txt)])
els = np.array([[int(v) - 1 for v in m.group(1).split(",")] for m in re.finditer(r'<elem id="\d+">([^<]+)</elem>', txt)])
z0 = np.where(np.abs(nodes[:, 2]) < 1e-12)[0]; remap = -np.ones(len(nodes), int); remap[z0] = np.arange(len(z0))
XY = nodes[z0, :2]; Q = remap[els[:, :4]]; assert (Q >= 0).all(); n2 = len(z0)
# --- final-state nodal data from the log: columns wx wy ef; last time block
blk, T = None, None
for ln in open(RUN + "/Mon_all.txt"):
    if ln.startswith("*Time"): T = float(ln.split("=")[1]); rows = []; blk = rows
    elif not ln.startswith("*") and ln.strip(): blk.append(ln.split())
A = np.array([[float(x) for x in r[1:4]] for r in blk]); ids = np.array([int(r[0]) for r in blk]); assert (ids - 1 == z0).all()
Wv, EF = A[:, :2], A[:, 2]; P = -K * EF
# --- shape functions at 2x2 Gauss points
g = 1 / math.sqrt(3); GP = [(-g, -g), (g, -g), (g, g), (-g, g)]
def N_dN(xi, eta):
    N = 0.25 * np.array([(1 - xi) * (1 - eta), (1 + xi) * (1 - eta), (1 + xi) * (1 + eta), (1 - xi) * (1 + eta)])
    dN = 0.25 * np.array([[-(1 - eta), -(1 - xi)], [(1 - eta), -(1 + xi)], [(1 + eta), (1 + xi)], [-(1 + eta), (1 - xi)]])
    return N, dN
def elem_geom(xi, eta):
    """vectorised over elements: returns N (4,), dNdx (ne,4,2), detJ (ne,)"""
    N, dN = N_dN(xi, eta); X = XY[Q]                                   # (ne,4,2)
    J = np.einsum("nka,kb->nab", X, dN)                                  # J[a,b] = dx_a/dxi_b
    det = J[:, 0, 0] * J[:, 1, 1] - J[:, 0, 1] * J[:, 1, 0]
    inv = np.empty_like(J); inv[:, 0, 0] = J[:, 1, 1] / det; inv[:, 1, 1] = J[:, 0, 0] / det; inv[:, 0, 1] = -J[:, 0, 1] / det; inv[:, 1, 0] = -J[:, 1, 0] / det
    dNdx = np.einsum("kb,nba->nka", dN, inv)                              # inv[b,a] = dxi_b/dx_a
    return N, dNdx, det
# --- (1) weak-form nodal forces f_i = int [ sigma : grad N_i + rho (w.grad) w N_i ] dV (steady), sigma = -p I + 2 mu D
F4node = np.zeros((n2, 2)); Pgp = {}
for (xi, eta) in GP:
    N, dNdx, det = elem_geom(xi, eta)
    wq = Wv[Q]; pq = P[Q]
    L = np.einsum("nka,nkb->nab", wq, dNdx)                               # L[i,j] = dw_i/dx_j
    D = 0.5 * (L + np.transpose(L, (0, 2, 1))); p = pq @ N
    sig = 2 * MU * D; sig[:, 0, 0] -= p; sig[:, 1, 1] -= p
    wgp = np.einsum("k,nka->na", N, wq); conv = RHO * np.einsum("nab,nb->na", L, wgp)
    contrib = np.einsum("nab,nkb->nka", sig, dNdx) + conv[:, None, :] * N[None, :, None]   # (ne,4,2)
    contrib *= det[:, None, None]                                          # Gauss weight 1
    np.add.at(F4node, Q.reshape(-1), contrib.reshape(-1, 2))
# --- boundary edges from topology
edge_count = {}
for qi, q in enumerate(Q):
    for k in range(4):
        a, b = int(q[k]), int(q[(k + 1) % 4]); edge_count.setdefault((min(a, b), max(a, b)), []).append((a, b, qi, k))
bed = [v[0] for v in edge_count.values() if len(v) == 1]
LCH, HCH = 2.5, 0.41
def cls(a, b):
    xa, ya = XY[a]; xb, yb = XY[b]
    if abs(xa) < 1e-9 and abs(xb) < 1e-9: return "inlet"
    if abs(xa - LCH) < 1e-9 and abs(xb - LCH) < 1e-9: return "outlet"
    if abs(ya) < 1e-9 and abs(yb) < 1e-9: return "bottom"
    if abs(ya - HCH) < 1e-9 and abs(yb - HCH) < 1e-9: return "top"
    return "obstacle"
kinds = np.array([cls(a, b) for a, b, qi, k in bed])
obst_nodes = sorted({a for (a, b, qi, k), kd in zip(bed, kinds) if kd == "obstacle"} | {b for (a, b, qi, k), kd in zip(bed, kinds) if kd == "obstacle"})
all_bnd = {a for a, b, qi, k in bed} | {b for a, b, qi, k in bed}
interior = np.array(sorted(set(range(n2)) - all_bnd))
# interior momentum residual (steady-state check of the stored field): relative to the typical nodal force magnitude
res_int = np.linalg.norm(F4node[interior], axis=1); scale = np.linalg.norm(F4node[obst_nodes], axis=1).max()
# --- obstacle force from the consistent nodal reactions: force of fluid on obstacle = - sum over obstacle nodes of f_i (f_i = force the wall exerts on the fluid)
F4 = -F4node[obst_nodes].sum(axis=0)
# --- facet kernels
def facet(a, b, qi, k):
    d = XY[b] - XY[a]; ln = math.hypot(*d); n = np.array([d[1], -d[0]]) / ln; return ln, n
def grad_at(qi, xi, eta):
    N, dN = N_dN(xi, eta); X = XY[Q[qi]]; J = X.T @ dN; Ji = np.linalg.inv(J); dNdx = dN @ Ji; Lm = Wv[Q[qi]].T @ dNdx; return 0.5 * (Lm + Lm.T), N
par = [lambda s: (s, -1.0), lambda s: (1.0, s), lambda s: (-s, 1.0), lambda s: (-1.0, -s)]
gs = [-g, g]
def obstacle_force(mode):
    F = {"cyl": np.zeros(2), "flag": np.zeros(2)}; Fp_ = {"cyl": np.zeros(2), "flag": np.zeros(2)}; Fv_ = {"cyl": np.zeros(2), "flag": np.zeros(2)}
    for (a, b, qi, k), kd in zip(bed, kinds):
        if kd != "obstacle": continue
        ln, n = facet(a, b, qi, k)
        xm = 0.5 * (XY[a, 0] + XY[b, 0]); ym = 0.5 * (XY[a, 1] + XY[b, 1])
        on_circle = abs(math.hypot(XY[a, 0] - 0.2, XY[a, 1] - 0.2) - 0.05) < 1e-6 and abs(math.hypot(XY[b, 0] - 0.2, XY[b, 1] - 0.2) - 0.05) < 1e-6
        grp = "cyl" if on_circle else "flag"
        if mode == 1:
            pf = 0.5 * (P[a] + P[b]) * ln * n
            fv = np.zeros(2)
            for s in gs:
                D, _ = grad_at(qi, *par[k](s)); fv += -2 * MU * D @ n * ln / 2
        else:
            pf = P[Q[qi]].mean() * ln * n
            D, _ = grad_at(qi, 0.0, 0.0); fv = -2 * MU * D @ n * ln
        Fp_[grp] += pf; Fv_[grp] += fv
    return Fp_, Fv_
F1p, F1v = obstacle_force(1); F2p, F2v = obstacle_force(2)
tot = lambda p, v: sum(p[g_] + v[g_] for g_ in p)
F1, F2 = tot(F1p, F1v), tot(F2p, F2v)
# --- control volume with a selectable pressure and term budget
def control_volume(pmode):
    terms = {k_: np.zeros(2) for k_ in ("p_inlet", "p_outlet", "p_top", "p_bottom", "v_inlet", "v_outlet", "v_top", "v_bottom", "mom_inlet", "mom_outlet")}
    Qin = Qout = 0.0
    for (a, b, qi, k), kd in zip(bed, kinds):
        if kd == "obstacle": continue
        ln, n = facet(a, b, qi, k)
        pf = (0.5 * (P[a] + P[b]) if pmode == "nodal" else P[Q[qi]].mean()) * ln * (-n)
        terms["p_" + kd] += pf
        fv = np.zeros(2); fm = np.zeros(2)
        for s in gs:
            D, _ = grad_at(qi, *par[k](s)); fv += 2 * MU * D @ n * ln / 2
            t = 0.5 * (s + 1); u = Wv[a] * (1 - t) + Wv[b] * t; fm += RHO * u * (u @ n) * ln / 2
        terms["v_" + kd] += fv
        if kd in ("inlet", "outlet"): terms["mom_" + kd] += fm
        fl = 0.5 * ((Wv[a] + Wv[b]) @ n) * ln
        if kd == "inlet": Qin -= fl
        if kd == "outlet": Qout += fl
    F = sum(v for kk, v in terms.items() if kk[0] in "pv") - terms["mom_inlet"] - terms["mom_outlet"]
    return F, terms, Qin, Qout
F3n, T3n, Qin, Qout = control_volume("nodal"); F3b, T3b, _, _ = control_volume("element")
# --- parity (independent definition: second difference along the wall chain, cavity-style) on the flag walls
def chain(yv):
    ids_ = [i for i in obst_nodes if abs(XY[i, 1] - yv) < 1e-9 and XY[i, 0] > 0.24898 and XY[i, 0] < 0.6000001]
    return sorted(ids_, key=lambda i: XY[i, 0])
def at(yv):
    c = chain(yv); p = P[c]; s = 4; v = [(p[k] - 0.5 * (p[k - 1] + p[k + 1])) / 2 for k in range(1 + s, len(c) - 1 - s)]; return float(np.sqrt(np.mean(np.square(v)))), len(c)
res = dict(run=tag, K=K, mu=MU, rho=RHO, t_final=T, n_cells=len(Q), n_nodes_layer=n2,
           F1=F1.tolist(), F2=F2.tolist(), F3nodal=F3n.tolist(), F3b=F3b.tolist(), F4_weak_reaction=F4.tolist(),
           F1_parts={g_: dict(p=F1p[g_].tolist(), v=F1v[g_].tolist()) for g_ in F1p}, F2_parts={g_: dict(p=F2p[g_].tolist(), v=F2v[g_].tolist()) for g_ in F2p},
           F3b_terms={k_: v.tolist() for k_, v in T3b.items()}, F3nodal_terms={k_: v.tolist() for k_, v in T3n.items()},
           Qin=Qin, Qout=Qout, mass_imbalance=Qout / Qin - 1,
           interior_residual_max_over_obstacle_nodal_force_max=float(res_int.max() / scale), interior_residual_rms_over_scale=float(np.sqrt(np.mean(res_int**2)) / scale),
           a_t_top=at(0.21), a_t_bottom=at(0.19),
           refs=dict(drag=14.2929, lift=1.11905),
           n_obstacle_nodes=len(obst_nodes), n_interior=len(interior))
# --- F2 - F3b budget: pressure on obstacle vs outer; difference of drag by terms
res["budget_drag"] = dict(F2_drag=float(F2[0]), F3b_drag=float(F3b[0]), diff=float(F3b[0] - F2[0]), diff_pct_of_F2=float(100 * (F3b[0] - F2[0]) / F2[0]),
    F3b_terms_x={k_: float(v[0]) for k_, v in T3b.items()}, F2_pressure_x=float(sum(F2p[g_][0] for g_ in F2p)), F2_viscous_x=float(sum(F2v[g_][0] for g_ in F2v)))
# --- xplt cross-check of the data path (element pressure / dilatation)
try:
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import importlib.util
    spec = importlib.util.spec_from_file_location("xa", os.path.join(os.path.dirname(os.path.abspath(__file__)), "xplt_reader.py")); xa = importlib.util.module_from_spec(spec); sys.argv = ["x"]
    try: spec.loader.exec_module(xa)
    except SystemExit: pass
    r = xa.read_xplt(glob.glob(RUN + "/*.xplt")[0]); S = r["states"][-1]; efx = S["vars"]["fluid dilatation"][:, 0][z0]; pe = S["vars"]["fluid pressure"][0][:, 0]
    res["xplt_check"] = dict(time=float(S["time"]), max_abs_diff_ef_log_vs_xplt=float(np.abs(efx - EF).max()), ef_scale=float(np.abs(EF).max()),
                             elem_pressure_vs_mean_nodal_pressure_max_abs_diff_Pa=float(np.abs(pe - P[Q].mean(axis=1)).max()),
                             elem_pressure_vs_GP_average_note="xplt element pressure is the element average of the Gauss-point pressure p = -K e_f(GP)")
except Exception as exc:
    res["xplt_check"] = f"not available: {exc}"
json.dump(res, open(OUT, "w"), indent=1)
f = lambda v: f"({v[0]:.5f}, {v[1]:.5f})"
print(tag, "t", T); print(" F1", f(F1), " F2", f(F2), " F3nodal", f(F3n), " F3b", f(F3b)); print(" F4 (weak-form nodal reactions)", f(F4))
print(" mass imbalance", Qout / Qin - 1, " interior residual max/scale", res["interior_residual_max_over_obstacle_nodal_force_max"])
print(" budget F3b terms x:", {k_: round(float(v[0]), 4) for k_, v in T3b.items()}); print(" xplt check", res["xplt_check"])

# ---- supplementary budget (added after the first independent reproduction; definitions fixed here): weak-form identity and where the interior residual sits
outer_nodes = sorted(all_bnd - set(obst_nodes))
conv_total = np.zeros(2)
for (xi, eta) in GP:
    N, dNdx, det = elem_geom(xi, eta); wq = Wv[Q]; L = np.einsum("nka,nkb->nab", wq, dNdx); wgp = np.einsum("k,nka->na", N, wq)
    conv_total += (RHO * np.einsum("nab,nb->na", L, wgp) * det[:, None]).sum(axis=0)
S_outer = F4node[outer_nodes].sum(axis=0); S_int = F4node[interior].sum(axis=0); S_obst = F4node[obst_nodes].sum(axis=0)
imax = interior[np.argmax(res_int)]
res["weak_identity"] = dict(sum_obstacle_f=S_obst.tolist(), sum_outer_f=S_outer.tolist(), sum_interior_f=S_int.tolist(), convective_volume_integral=conv_total.tolist(),
    closure_sum_all_minus_conv=(S_obst + S_outer + S_int - conv_total).tolist(), interior_residual_max_at_xy=XY[imax].tolist(),
    interior_residual_max=float(res_int.max()), obstacle_nodal_force_scale=float(scale), n_outer=len(outer_nodes))
# F4 split into pressure / viscous / convective parts on the obstacle nodes
parts = {"pressure": np.zeros((n2, 2)), "viscous": np.zeros((n2, 2)), "convective": np.zeros((n2, 2))}
for (xi, eta) in GP:
    N, dNdx, det = elem_geom(xi, eta); wq = Wv[Q]; pq = P[Q]; L = np.einsum("nka,nkb->nab", wq, dNdx); D = 0.5 * (L + np.transpose(L, (0, 2, 1))); p = pq @ N
    wgp = np.einsum("k,nka->na", N, wq); conv = RHO * np.einsum("nab,nb->na", L, wgp)
    cp = -p[:, None, None] * dNdx * det[:, None, None]                                  # -p grad N_i
    cv = 2 * MU * np.einsum("nab,nkb->nka", D, dNdx) * det[:, None, None]
    cc = conv[:, None, :] * N[None, :, None] * det[:, None, None]
    for nm, c in (("pressure", cp), ("viscous", cv), ("convective", cc)): np.add.at(parts[nm], Q.reshape(-1), c.reshape(-1, 2))
res["F4_parts"] = {nm: (-parts[nm][obst_nodes].sum(axis=0)).tolist() for nm in parts}
json.dump(res, open(OUT, "w"), indent=1)
print(" weak identity:", json.dumps(res["weak_identity"])); print(" F4 parts:", res["F4_parts"])

# ---- weak-form nodal sums on the outer boundary by group (corner nodes: inlet/outlet take precedence over top/bottom) for the F3b-F4 budget
def grp(i):
    if abs(XY[i, 0]) < 1e-9: return "inlet"
    if abs(XY[i, 0] - LCH) < 1e-9: return "outlet"
    if abs(XY[i, 1]) < 1e-9: return "bottom"
    return "top"
by = {k_: F4node[[i for i in outer_nodes if grp(i) == k_]].sum(axis=0).tolist() for k_ in ("inlet", "outlet", "top", "bottom")}
res["weak_outer_by_group"] = by
res["facet_sigma_n_by_group"] = {k_: (np.array(T3b["p_" + k_]) + np.array(T3b["v_" + k_])).tolist() for k_ in ("inlet", "outlet", "top", "bottom")}
json.dump(res, open(OUT, "w"), indent=1)
print(" outer groups weak:", by); print(" outer groups facet:", res["facet_sigma_n_by_group"])
