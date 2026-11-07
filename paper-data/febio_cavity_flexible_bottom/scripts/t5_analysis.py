#!/usr/bin/env python3
import os, os as _os; HERE = _os.path.dirname(_os.path.abspath(__file__)); PKG = _os.path.dirname(HERE)   # package root; run directories of retained HPC outputs: env CAVITY_RUNS (default runs/cavityNext)
"""Round-2 tiny Poiseuille channel analysis (criteria P0-P6 of FEBIO_FLUX_PRESERVING_PROJECTION.md). usage: t5_analysis.py OUT.json"""
import sys, re, os, json, numpy as np, scipy.sparse as sp
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); from compare_runs import read_blocks, iters
R = os.path.join(HERE, "../../runs/cavityExp/"); K = 1e6; MU = 0.01; ALPHA = 1e-3
NODES = {int(m.group(1)): tuple(float(v) for v in m.group(2).split(",")) for m in re.finditer(r'<node id="(\d+)">([^<]+)</node>', open(R + "T3_stock/T3_straight.feb").read())}
def col(run, f, t=None):
    b = read_blocks(R + run + f"/Mon_{f}.txt"); t = t or max(b); a = b[t]; y = np.array([NODES[int(i)][1] for i in a[:, 0]]); x = np.array([NODES[int(i)][0] for i in a[:, 0]]); o = np.argsort(y); return t, x[o], y[o], a[o]
out = {}
# P0 rim-once alpha=0 vs stock steady state, all Mon files, last block
p0 = {}
for f in sorted(x[4:-4] for x in os.listdir(R + "T3_stock") if x.startswith("Mon_")):
    a = read_blocks(R + f"T3_stock/Mon_{f}.txt"); b = read_blocks(R + f"T5_rimonce_a0/Mon_{f}.txt"); t = max(a); A, B_ = a[t], b[t]
    vs = max(np.abs(A[:, 3]).max(), np.abs(A[:, 4]).max()); es = max(np.abs(A[:, 5]).max(), 1e-300)
    p0[f] = dict(wx=float(np.abs(A[:, 3] - B_[:, 3]).max() / vs), wy=float(np.abs(A[:, 4] - B_[:, 4]).max() / vs), ef=float(np.abs(A[:, 5] - B_[:, 5]).max() / es))
out["P0_rim_once_alpha0_vs_stock_last_block"] = p0
out["P0_worst"] = max(max(v.values()) for v in p0.values())
def analyse(run):
    r = {}; ps = {}; Qin = None
    t, x, y, a = col(run, "inlet"); Qin = float(np.trapz(a[:, 3], y)); inlet_wx = a[:, 3].copy(); r["Q_inlet"] = Qin
    for f, xx in (("x2", 2), ("x7", 7), ("x10", 10), ("x12", 12), ("x13", 13), ("outlet", 14)):
        t, x, y, a = col(run, f); p = -K * a[:, 5]; Q = float(np.trapz(a[:, 3], y)); ps[xx] = p.mean()
        r[f] = dict(Q_over_Qin_minus_1=Q / Qin - 1, profile_maxdiff_over_center=float(np.abs(a[:, 3] - inlet_wx).max() / inlet_wx.max()), p_mean=float(p.mean()),
                    y_alt_rms_Pa=float(np.sqrt(np.mean(((p[:-2] - 2 * p[1:-1] + p[2:]) / 4) ** 2))) if f != "outlet" else 0.0)
    xs = np.array([2, 7, 10, 12, 13]); r["dpdx_fit"] = float(np.polyfit(xs, [ps[k] for k in xs], 1)[0]); r["dpdx_discrete_poiseuille"] = -0.12 / (1 - 0.01); r["dpdx_rel_err"] = r["dpdx_fit"] / r["dpdx_discrete_poiseuille"] - 1
    r["max_abs_flux_err_all_sections"] = max(abs(r[f]["Q_over_Qin_minus_1"]) for f in ("x2", "x7", "x10", "x12", "x13", "outlet"))
    r["y_alt_rms_max_Pa"] = max(r[f]["y_alt_rms_Pa"] for f in ("x2", "x7", "x10", "x12", "x13"))
    r["max_profile_err_x7_x10_x12"] = max(r[f]["profile_maxdiff_over_center"] for f in ("x7", "x10", "x12"))
    log = [f for f in os.listdir(R + run) if f.endswith(".log") and f != "run.log"][0]; its = iters(R + run + "/" + log); txt = open(R + run + "/" + log, errors="ignore").read()
    r["iterations_total"] = sum(its); r["iterations_per_step_mean"] = float(np.mean(its)); r["normal_termination"] = "N O R M A L   T E R M I N A T I O N" in txt; r["max_reformations_event"] = "Max nr of reformations" in txt
    m = re.search(r"Total elapsed time[ .]*:\s*([\d:]+)\s*\(([\d.]+) sec\)", txt); r["wall_s"] = float(m.group(2)) if m else None
    return r
for run in ("T3_stock", "T3_ppp_a075", "T5_rimonce_a0", "T5_ppp_a1e-3"):
    try: out[run] = analyse(run)
    except Exception as e: out[run] = {"error": str(e)}
# P5 independent re-assembly on the straight channel
import importlib; sys.path.insert(0, HERE + "/../pressure_mode"); import b4_operator as B
src = open(HERE + "/../pressure_mode/b4_residual_check.py").read().split("# parity decomposition")[0]
src = src.replace('rd = R + "/runs/cavityNext/C1b_fluid_fixedplate_f1_pin"', 'rd = __import__("os").environ["RD_PPP"]').replace('g = B.geometry("cavity", 1)', 'g = B.geometry("straight", 1)').replace("wfree, pfree, wfix, pfix = B.bcs(g, 1.15)", "wfree, pfree, wfix, pfix = B.bcs(g, None)").replace('aud.read_mon_block(rd + "/Mon_all.txt", 300.0)', 'aud.read_mon_block(rd + "/Mon_all.txt", 100.0)')
def pp_check(run, alpha):
    os.environ["RD_PPP"] = os.path.abspath(R + run); sys.argv = ["x", "/dev/null"]; ns = {"__name__": "ns"}; exec(compile(src, "adapted", "exec"), ns)
    g, X, n, p, w, G, Mm, wfree, pfree, pfix, res = (ns[k] for k in ("g", "X", "n", "p", "w", "G", "Mm", "wfree", "pfree", "pfix", "res"))
    quads = g["quads"]; nbrs = None
    # PPP~ term with redistribution (constrained = pfix, equal weights)
    S = np.zeros(n)
    for q in quads:
        x = X[list(q)]; V = 0; I = np.zeros(4); pq = []; Wq = []; Hq = []
        for (r, s) in B.GP:
            dr, ds = B.dN(r, s); dJ = np.linalg.det(np.array([[dr @ x[:, 0], dr @ x[:, 1]], [ds @ x[:, 0], ds @ x[:, 1]]])); Hh = B.N(r, s); Wq.append(dJ); Hq.append(Hh); pq.append(Hh @ p[list(q)]); V += dJ; I += Hh * dJ
        pbar = sum(Wq[k] * pq[k] for k in range(4)) / V; s_e = np.array([-(alpha / MU) * sum(Wq[k] * Hq[k][a] * (pq[k] - pbar) for k in range(4)) for a in range(4)])
        con = np.array([1 if i in pfix else 0 for i in q]); nF = 4 - con.sum()
        if nF == 0: continue
        C = s_e[con == 1].sum(); s_t = np.where(con == 1, 0.0, s_e + C / nF)
        for a in range(4): S[q[a]] += s_t[a]
    b, prof = B.inlet_load(g); Rc = G.T @ w + b + S; pf = np.array(sorted(set(pfree.tolist())))
    return dict(momentum_rel=res["momentum_residual_free_rows"]["total_over_pressure_gradient"], continuity_rel_with_PPP=float(np.linalg.norm(Rc[pf]) / np.linalg.norm((G.T @ w)[pf])), continuity_rel_without_PPP=float(np.linalg.norm((G.T @ w + b)[pf]) / np.linalg.norm((G.T @ w)[pf])), sum_S_free=float(S[pf].sum()), sum_S_all=float(S.sum()))
out["P5_reassembly"] = {"T3_stock": pp_check("T3_stock", 0.0), "T5_ppp_a1e-3": pp_check("T5_ppp_a1e-3", ALPHA)}
json.dump(out, open(sys.argv[1] if len(sys.argv) > 1 and sys.argv[1] != "/dev/null" else HERE + "/results/T5_analysis.json", "w"), indent=1)
print(json.dumps(out, indent=1))
