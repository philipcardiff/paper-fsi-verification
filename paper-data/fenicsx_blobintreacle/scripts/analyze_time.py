"""Pairwise BDF2 orders at t = 1 (M2, F = 1) against the finer same-scheme run dt_ref = 0.00125,
and the orders a PURE second-order error would give for the same levels (reference carries its own
error C*dt_ref^2).  Usage: python3 -B scripts/analyze_time.py   (standard library only)
Reads data/time_values_t1_F1_M2.csv, writes data/bdf2_pairwise_vs_prediction.csv"""
import csv, math, os
here = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
V = {}
for r in csv.DictReader(open(os.path.join(here, "time_values_t1_F1_M2.csv"))):
    if r["scheme"] == "BDF2": V[float(r["dt"])] = r
ref = V[0.00125]; dts = [0.01, 0.005, 0.0025]; dr = 0.00125
pred = [math.log((dts[i] ** 2 - dr ** 2) / (dts[i + 1] ** 2 - dr ** 2)) / math.log(2) for i in range(2)]
cols = {"apex_ux": "apex_ux_t1", "apex_uy": "apex_uy_t1", "Fx": "Fx_fluid_on_solid_t1", "Fy": "Fy_fluid_on_solid_t1"}
with open(os.path.join(here, "bdf2_pairwise_vs_prediction.csv"), "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["observable", "pair", "observed_order", "pure_second_order_prediction", "deviation", "within_0.05"])
    for name, c in cols.items():
        e = [abs(float(V[d][c]) - float(ref[c])) for d in dts]
        for i in range(2):
            o = math.log(e[i] / e[i + 1]) / math.log(2)
            w.writerow([name, "dt %g/%g" % (dts[i], dts[i + 1]), "%.2f" % o, "%.2f" % pred[i], "%+.2f" % (o - pred[i]), "yes" if abs(o - pred[i]) <= 0.05 else "no"])
