#!/usr/bin/env python3
"""Observed orders of accuracy from postProcessing/mmsFsiErrors/errors.dat.

Usage: orders.py run1 run2 run3 ... [--times 0.1 0.2 0.3 0.4] [--cols U_L2 p_L2 ...]
            [--ratio 2] [--max] [--md out.md]

Runs are ordered coarse to fine with a constant refinement ratio. For each
time, the error of each quantity and the observed order
log(e_coarse/e_fine)/log(ratio) between successive runs are printed. With
--max, the maximum over 0 < t <= time is used instead of the value at time.
"""
import argparse
import math
import os

ap = argparse.ArgumentParser()
ap.add_argument("runs", nargs="+")
ap.add_argument("--times", nargs="+", type=float, default=[0.1, 0.2, 0.3, 0.4])
ap.add_argument("--cols", nargs="+", default=["U_L2", "U_Linf", "p_L2", "p_Linf", "D_L2", "D_Linf",
                                              "Us_L2", "intD_L2", "intTrac_L2", "intTrac_Linf", "intP_L2"])
ap.add_argument("--ratio", type=float, default=2.0)
ap.add_argument("--max", action="store_true")
ap.add_argument("--md", default=None)
a = ap.parse_args()


def load(run):
    f = os.path.join(run, "postProcessing", "mmsFsiErrors", "errors.dat")
    with open(f) as fh:
        lines = fh.read().splitlines()
    names = lines[0].lstrip("#").split()
    data = [list(map(float, l.split())) for l in lines[1:] if l.strip()]
    return names, data


tabs = []
out = []
loaded = [load(r) for r in a.runs]
names = loaded[0][0]
for T in a.times:
    rows = []
    for (nm, data) in loaded:
        if a.max:
            sel = [d for d in data if d[0] <= T + 1e-9]
            if not sel or abs(sel[-1][0] - T) > 1e-6:
                rows.append(None)
                continue
            rows.append({c: max(d[nm.index(c)] for d in sel) for c in a.cols})
        else:
            m = [d for d in data if abs(d[0] - T) < 1e-6]
            rows.append({c: m[0][nm.index(c)] for c in a.cols} if m else None)
    hdr = f"t = {T}" + (" (max over time)" if a.max else "")
    out.append(f"\n### {hdr}\n")
    out.append("| quantity | " + " | ".join(os.path.basename(r.rstrip('/')) for r in a.runs) + " | orders |")
    out.append("|---|" + "---|" * len(a.runs) + "---|")
    for c in a.cols:
        vals = [r[c] if r else None for r in rows]
        ords = []
        for i in range(len(vals) - 1):
            if vals[i] and vals[i + 1] and vals[i] > 0 and vals[i + 1] > 0:
                ords.append(f"{math.log(vals[i] / vals[i + 1]) / math.log(a.ratio):.2f}")
            else:
                ords.append("-")
        out.append(f"| {c} | " + " | ".join(f"{v:.3e}" if v else "-" for v in vals) + " | " + ", ".join(ords) + " |")
text = "\n".join(out)
print(text)
if a.md:
    with open(a.md, "w") as fh:
        fh.write(text + "\n")
