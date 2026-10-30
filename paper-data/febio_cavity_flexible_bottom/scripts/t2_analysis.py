import os, os as _os; HERE = _os.path.dirname(_os.path.abspath(__file__)); PKG = _os.path.dirname(HERE)   # package root; run directories of retained HPC outputs: env CAVITY_RUNS (default runs/cavityNext)
import sys, re, numpy as np, json
sys.path.insert(0, "."); from compare_runs import read_blocks
R = "../../runs/cavityExp/"; K = 1e6
NODES = {int(m.group(1)): tuple(float(v) for v in m.group(2).split(",")) for m in re.finditer(r'<node id="(\d+)">([^<]+)</node>', open(R + "T1_stock/T1_straight.feb").read())}
def load(run):
    out = {}
    for f in ("x2", "x7", "x10", "x12", "x13", "outlet", "inlet"):
        b = read_blocks(R + run + f"/Mon_{f}.txt"); t = max(b); a = b[t]; a = np.hstack([a, np.array([[NODES[int(i)][0], NODES[int(i)][1]] for i in a[:, 0]])]); a = a[np.argsort(a[:, 7])]; out[f] = (t, a)
    return out
A, B = load("T1_stock"), load("T2_ppp_a075"); res = {}
for f in A:
    t, a = A[f]; _, b = B[f]; Qa = np.trapz(a[:, 3], a[:, 7]); Qb = np.trapz(b[:, 3], b[:, 7]); pa = -K * a[:, 5]; pb = -K * b[:, 5]
    # alternation in y: amplitude of (-1)^j component after removing smooth part (1-2-1 filter on the column)
    def alt(p): return float(np.sqrt(np.mean(((p[:-2] - 2 * p[1:-1] + p[2:]) / 4) ** 2)))
    res[f] = dict(t=t, x=float(a[0, 6]), Q_stock=float(Qa), Q_ppp=float(Qb), Q_ppp_over_stock_minus1=float(Qb / Qa - 1), pmean_stock=float(pa.mean()), pmean_ppp=float(pb.mean()), pmean_diff=float(pb.mean() - pa.mean()), p_colwise_maxdiff=float(np.abs(pa - pb).max()), p_scale=float(np.abs(pa).max()),
                  yalt_rms_stock_Pa=alt(pa), yalt_rms_ppp_Pa=alt(pb), wx_maxdiff_over_U=float(np.abs(a[:, 3] - b[:, 3]).max()))
print(json.dumps(res, indent=1)); json.dump(res, open("results/T2_ppp_a075_vs_stock_columns.json", "w"), indent=1)
