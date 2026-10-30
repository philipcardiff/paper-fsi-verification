import os, os as _os; HERE = _os.path.dirname(_os.path.abspath(__file__)); PKG = _os.path.dirname(HERE)   # package root; run directories of retained HPC outputs: env CAVITY_RUNS (default runs/cavityNext)
import sys, json, numpy as np
sys.path.insert(0, HERE); import b4_operator as B
A = np.load("results/stock_C1b_analysis.npz"); C = np.load("results/C1b_ppp_alpha075_analysis.npz"); Z = np.load("results/C1b_ppp_alpha0_analysis.npz")
g = B.geometry("cavity", 1); X = g["X"]; n = len(X); K = 1e6
nb = {}
for q in g["quads"]:
    for k in range(4): nb.setdefault(q[k], set()).update([q[(k + 1) % 4], q[(k - 1) % 4]])
def smooth(p):   # parity-removing 5-point average (diagnostic filter, identical for both runs; NOT used in any solver)
    s = p.copy()
    for i in range(n):
        if len(nb[i]) == 4: s[i] = 0.5 * p[i] + 0.125 * sum(p[j] for j in nb[i])
    return s
pa, pc, wa, wc = A["p"], C["p"], A["w"].reshape(n, 2), C["w"].reshape(n, 2)
sa, sc = smooth(pa), smooth(pc)
res = {}
res["pressure_rel_L2_diff_raw"] = float(np.linalg.norm(pc - pa) / np.linalg.norm(pa)); res["pressure_rel_L2_diff_smooth5pt"] = float(np.linalg.norm(sc - sa) / np.linalg.norm(sa))
res["pressure_rel_max_diff_smooth_over_max_p"] = float(np.abs(sc - sa).max() / np.abs(sa).max())
res["velocity_rel_L2_diff"] = float(np.linalg.norm(wc - wa) / np.linalg.norm(wa)); res["velocity_max_abs_diff_over_U"] = float(np.abs(wc - wa).max())
res["pressure_rms_roughness_p_minus_smooth_Pa"] = dict(stock=float(np.sqrt(np.mean((pa - sa) ** 2))), ppp075=float(np.sqrt(np.mean((pc - sc) ** 2))))
for nm, (xl, xh, yl, yh) in {"inlet_x<1": (0, 1, 1, 2), "channel_above_cavity": (2, 6, 1, 2), "cavity": (2, 6, 0, 0.75), "downstream": (6.5, 14, 1, 2)}.items():
    m = (X[:, 0] >= xl) & (X[:, 0] <= xh) & (X[:, 1] >= yl) & (X[:, 1] <= yh)
    res.setdefault("region_smooth_p_mean_Pa", {})[nm] = dict(stock=float(sa[m].mean()), ppp075=float(sc[m].mean()), rel=float(sc[m].mean() / sa[m].mean() - 1))
# flux at every node column along the channel (y in [1,2]) for the straight part
cols = {}
for xv in (1.0, 2.0, 7.0, 9.0, 12.0, 13.0, 14.0):
    xc = X[np.argmin(abs(X[:, 0] - xv)), 0]; ids = [i for i in range(n) if abs(X[i, 0] - xc) < 1e-6 and X[i, 1] >= 1 - 1e-9]; ids.sort(key=lambda i: X[i, 1]); yy = X[ids, 1]
    cols[f"x={xc:.3f}"] = dict(Q_stock=float(np.trapz(wa[ids, 0], yy)), Q_ppp075=float(np.trapz(wc[ids, 0], yy)))
res["flux_columns_channel"] = cols
# discrete reaction on the floor: pair-mean and pressure force from the smooth part
res["Fy_pressure_stock_F1_F2"] = [-4.6637167409, -5.483270450046501]
json.dump(res, open("results/stage4_comparison.json", "w"), indent=1); print(json.dumps(res, indent=1))
