#!/usr/bin/env python3
"""Rigid-plate steady SU2: F_x, F_y on the plate (N/m; CD/CL with force reference 1),
p at A = (4, 0), mean p on the plate, residuals, observed orders (refinement ratio 2)."""
import csv, math
from pathlib import Path

rows = []
for L in (1, 2, 4, 8, 16):
    d = Path(f"L{L}")
    h = d / "fluid-su2/history_steady.csv"
    if not (d / "exit_status.txt").exists() or not h.exists():
        print(f"L{L}: not finished"); continue
    r = [dict((k.strip().strip('"'), v.strip()) for k, v in x.items()) for x in csv.DictReader(open(h))][-1]
    s = list(csv.DictReader(open(d / "fluid-su2/surface_steady.csv")))
    pts = sorted((float(x["x"]), float(x["Pressure"])) for x in s if abs(float(x["y"])) < 1e-9)
    pA = min(pts, key=lambda q: abs(q[0] - 4.0))
    pmean = sum(0.5 * (pts[k][1] + pts[k+1][1]) * (pts[k+1][0] - pts[k][0]) for k in range(len(pts) - 1)) / 4.0
    rows.append(dict(level=L, Fx=float(r["CD"]), Fy=float(r["CL"]), pA=pA[1], xA=pA[0], pmean=pmean, iters=int(r["Inner_Iter"]),
                     res=max(float(r["rms[P]"]), float(r["rms[U]"]), float(r["rms[V]"])),
                     exit=(d / "exit_status.txt").read_text().strip(), wall=(d / "walltime.txt").read_text().strip()))
print(f"{'L':>3} {'F_x (N/m)':>14} {'F_y (N/m)':>14} {'p(A) (Pa)':>12} {'mean p (Pa)':>12} {'iters':>6} {'max rms':>8} {'exit':>4} {'wall s':>6}")
for q in rows:
    print(f"{q['level']:>3} {q['Fx']:14.8f} {q['Fy']:14.8f} {q['pA']:12.8f} {q['pmean']:12.8f} {q['iters']:6d} {q['res']:8.2f} {q['exit']:>4} {q['wall']:>6}")
for key in ("Fx", "Fy", "pA"):
    v = [q[key] for q in rows]
    out = []
    for k in range(len(v) - 2):
        r = (v[k] - v[k+1]) / (v[k+1] - v[k+2]) if v[k+1] != v[k+2] else float("nan")
        out.append(f"L{rows[k]['level']}-{rows[k+2]['level']}: " + (f"{math.log(r)/math.log(2):.2f}" if r > 0 else f"non-monotone ({r:.2f})"))
    print(f"observed order {key}: " + ", ".join(out))
with open("rigid_fluid.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
