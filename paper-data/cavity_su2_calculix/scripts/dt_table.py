#!/usr/bin/env python3
"""Time-step independence of the steady values and late-time wander of the dt = 100 s runs.

    dt_table.py DT_CHECK.csv HIST_DIR [--last 25]

DT_CHECK.csv: evaluate_coupled.py --csv output with the dt = 5/25/50 s runs (case names L<n>-dt<dt>)
and the dt = 100 s runs (case names L<n>). For every level and quantity, the difference of each
dt run from the dt = 100 s value (absolute and relative). HIST_DIR: the dt = 100 s histories
(evaluate_coupled.py --hist-dir); for each level, max - min of each quantity over the last
--last windows.
"""
import argparse
import csv
from pathlib import Path

Q = (("uy", "u_y(A) (m)"), ("ux", "u_x(A) (m)"), ("Fy", "F_y (N/m)"), ("Fx", "F_x (N/m)"))
H = {"uy": "uy_m", "ux": "ux_m", "Fy": "Fy_N_per_m", "Fx": "Fx_N_per_m"}

ap = argparse.ArgumentParser()
ap.add_argument("csv")
ap.add_argument("hist_dir")
ap.add_argument("--last", type=int, default=25)
a = ap.parse_args()

rows = list(csv.DictReader(open(a.csv)))
ref = {r["level"]: r for r in rows if r["case"] == f"L{r['level']}"}
print("| Level | dt (s) | windows | " + " | ".join(f"Δ{n}" for _, n in Q) + " | max rel. Δ |")
print("|---|---:|---:|" + "---:|" * (len(Q) + 1))
for r in rows:
    if r["case"] == f"L{r['level']}":
        continue
    b = ref[r["level"]]
    d = {q: float(r[q]) - float(b[q]) for q, _ in Q}
    rel = max(abs(d[q]) / abs(float(b[q])) for q, _ in Q)
    print(f"| {r['level']} | {float(r['dt']):g} | {r['windows']} | " + " | ".join(f"{d[q]:+.2e}" for q, _ in Q)
          + f" | {rel:.1e} |")
print(f"\nΔ = value at that dt minus the dt = 100 s value.\n")
print(f"| Level | " + " | ".join(f"range of {n}, last {a.last} windows" for _, n in Q) + " |")
print("|---|" + "---:|" * len(Q))
for L in sorted(ref, key=int):
    h = list(csv.DictReader(open(Path(a.hist_dir) / f"L{L}.csv")))[-a.last:]
    rng = {q: max(float(x[H[q]]) for x in h) - min(float(x[H[q]]) for x in h) for q, _ in Q}
    print(f"| {L} | " + " | ".join(f"{rng[q]:.1e} ({rng[q] / abs(float(ref[L][q])):.1e} rel.)" for q, _ in Q) + " |")
