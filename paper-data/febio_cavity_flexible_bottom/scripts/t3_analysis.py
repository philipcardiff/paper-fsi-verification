import os, os as _os; HERE = _os.path.dirname(_os.path.abspath(__file__)); PKG = _os.path.dirname(HERE)   # package root; run directories of retained HPC outputs: env CAVITY_RUNS (default runs/cavityNext)
import sys, re, numpy as np, json
sys.path.insert(0, "."); from compare_runs import read_blocks
R = "../../runs/cavityExp/"; K = 1e6; mu = 0.01
NODES = {int(m.group(1)): tuple(float(v) for v in m.group(2).split(",")) for m in re.finditer(r'<node id="(\d+)">([^<]+)</node>', open(R + "T3_stock/T3_straight.feb").read())}
def col(run, f, t=None):
    b = read_blocks(R + run + f"/Mon_{f}.txt"); t = t or max(b); a = b[t]; y = np.array([NODES[int(i)][1] for i in a[:, 0]]); x = np.array([NODES[int(i)][0] for i in a[:, 0]]); o = np.argsort(y); return t, x[o], y[o], a[o]
out = {}
for run in ("T3_stock", "T3_ppp_a075"):
    r = {}
    ps = {}
    for f, xx in (("x2", 2), ("x7", 7), ("x10", 10), ("x12", 12), ("x13", 13)):
        t, x, y, a = col(run, f); p = -K * a[:, 5]; wx = a[:, 3]; yy = y; Q = np.trapz(wx, yy); exact = 6 * (yy - 1) * (2 - yy)   # unit-mean Poiseuille
        ps[xx] = p.mean()
        r[f] = dict(t=t, Q=float(Q), wx_err_max_vs_exact_poiseuille=float(np.abs(wx - exact).max()), wx_err_rel_vs_exact_poiseuille=float(np.abs(wx - exact).max() / 1.5), p_mean=float(p.mean()), y_alt_rms=float(np.sqrt(np.mean(((p[:-2] - 2 * p[1:-1] + p[2:]) / 4) ** 2))), wy_max=float(np.abs(a[:, 4]).max()))
    xs = np.array(sorted(ps)); r["dpdx_fit_x2_to_x13"] = float(np.polyfit(xs, [ps[k] for k in xs], 1)[0]); r["dpdx_exact"] = -12 * mu * 1.0 / 1.0
    r["p_x2_minus_x13_over_11"] = float((ps[2] - ps[13]) / 11)
    # time history of p mean at x7 (last 5 blocks) for steadiness
    b = read_blocks(R + run + "/Mon_x7.txt"); ts = sorted(b); r["p_x7_last3"] = [float((-K * b[t][:, 5]).mean()) for t in ts[-3:]]; r["p_x7_t=50"] = float((-K * b[ts[100]][:, 5]).mean()); r["t_last"] = ts[-1]
    out[run] = r
print(json.dumps(out, indent=1)); json.dump(out, open("results/T3_straight_vs_exact.json", "w"), indent=1)
