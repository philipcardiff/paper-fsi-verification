#!/usr/bin/env python3
"""FSI1 study table: closing-20 % means, errors vs Featflow 7+0, observed orders and
Richardson extrapolation from levels 1/2/4 (as in the solids4foam FSI1 study), and the
solid-only refinement on the 1x fluid mesh.     study-table.py STUDYDIR"""
import glob, math, sys
from pathlib import Path
S = Path(sys.argv[1]); depth = 0.015; t0 = 2.0; tEnd = 30.0; tWin = tEnd - 0.2 * (tEnd - t0)
ref = {"ux": 2.270493e-5, "uy": 8.208773e-4, "drag": 14.29426, "lift": 0.7637460}
s4f = {"ux": (2.20672e-5, 2.24843e-5), "uy": (7.04258e-4, 7.82820e-4),
       "drag": (14.18027, 14.26063), "lift": (0.803445, 0.775620)}   # solids4foam 1x, 2x

def means(case):
    f = {}
    for p in sorted(glob.glob(str(S / case / "fluid-openfoam/postProcessing/forces/*/force.dat"))):
        for l in open(p):
            if l[0] != "#" and l.strip():
                c = l.split(); f[round(float(c[0]), 8)] = (float(c[1]) / depth, float(c[2]) / depth)
    d = {}
    for l in open(S / case / "solid-calculix/precice-Solid-watchpoint-PointA.log"):
        c = l.split()
        if c[0] != "Time":
            d[round(float(c[0]) + t0, 8)] = (float(c[5]), float(c[6]))
    m = lambda h, i: (lambda v: sum(v) / len(v))([x[i] for t, x in h.items() if t >= tWin - 1e-9])
    return {"ux": m(d, 0), "uy": m(d, 1), "drag": m(f, 0), "lift": m(f, 1)}

L = {k: means(k) for k in ("level1-rbf", "level2", "level4")}
F = {k: means(k) for k in ("level1-rbf", "fluid1-solid2", "fluid1-solid4")}
e = lambda v, q: 100 * (v - ref[q]) / ref[q]
print("Mesh study (fluid and solid refined together); values = mean over t = 24.4-30 s\n")
print("| Quantity | 1x | 2x | 4x | Featflow 7+0 | Error 1x / 2x / 4x | p (3 levels) | Extrapolated | p vs Featflow (2x->4x) |")
print("|---|---:|---:|---:|---:|---:|---:|---:|---:|")
for q in ref:
    a, b, c = (L[k][q] for k in L)
    r = (b - a) / (c - b)
    p = math.log(r) / math.log(2) if r > 0 else float("nan")
    ext = c + (c - b) / (2 ** p - 1) if r > 0 else float("nan")
    pf = math.log(abs(b - ref[q]) / abs(c - ref[q])) / math.log(2)
    print(f"| {q} | {a:.6g} | {b:.6g} | {c:.6g} | {ref[q]:.7g} | {e(a,q):+.2f} / {e(b,q):+.2f} / {e(c,q):+.2f} % "
          f"| {p:.2f} | {ext:.6g} ({e(ext,q):+.2f} %) | {pf:.2f} |")
print("\nSolid-only refinement on the 1x fluid mesh (5 336 cells)\n")
print("| Quantity | solid 105x6 | 210x12 | 420x24 | change 105x6 -> 420x24 | solids4foam 1x / 2x error |")
print("|---|---:|---:|---:|---:|---:|")
for q in ref:
    a, b, c = (F[k][q] for k in F)
    print(f"| {q} | {a:.6g} | {b:.6g} | {c:.6g} | {100*(c-a)/a:+.3f} % | {e(s4f[q][0],q):+.2f} / {e(s4f[q][1],q):+.2f} % |")
