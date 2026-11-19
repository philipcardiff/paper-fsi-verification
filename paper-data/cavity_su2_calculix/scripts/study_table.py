#!/usr/bin/env python3
"""cavityFlexibleBottom mesh study (SU2 + CalculiX, fluid and solid refined together), from
evaluate_coupled.py --csv output (one row per level; levels 1, 2, 4, 8[, 16], dx = 0.1/L m).

    study_table.py RESULTS.csv

For every consecutive triple (L, 2L, 4L): observed order p = ln((f_L - f_2L)/(f_2L - f_4L)) / ln 2,
Richardson f_R = f_4L + (f_4L - f_2L)/(2^p - 1), GCI(4L) = Fs |(f_4L - f_2L)/f_4L| / (2^p - 1), Fs = 1.25.
Non-monotone triples (ratio <= 0) get no order. No reference values: blind comparison.
"""
import csv
import math
import sys

Q = ("uy", "ux", "Fy", "Fx")
rows = {int(r["level"]): r for r in csv.DictReader(open(sys.argv[1]))}
Ls = sorted(rows)
f = {q: {L: float(rows[L][q]) for L in Ls} for q in Q}

print("| Level | dx (m) | u_y(A) (m) | u_x(A) (m) | F_y (N/m) | F_x (N/m) |")
print("|---|---:|---:|---:|---:|---:|")
for L in Ls:
    print(f"| {L} | {0.1/L:g} | " + " | ".join(f"{f[q][L]:.8g}" for q in Q) + " |")
print()
print("| Quantity | Levels | Ratio | Observed order | Richardson | GCI (finest) |")
print("|---|---|---:|---:|---:|---:|")
for q in Q:
    for a in range(len(Ls) - 2):
        L1, L2, L4 = Ls[a:a + 3]
        if not (L2 == 2 * L1 and L4 == 2 * L2):
            continue
        f1, f2, f4 = f[q][L1], f[q][L2], f[q][L4]
        r = (f1 - f2) / (f2 - f4) if f2 != f4 else float("nan")
        if r > 0 and not math.isnan(r):
            p = math.log(r) / math.log(2)
            fr = f4 + (f4 - f2) / (2 ** p - 1)
            gci = 1.25 * abs((f4 - f2) / f4) / (2 ** p - 1)
            print(f"| {q} | {L1}-{L2}-{L4} | {r:.3f} | {p:.2f} | {fr:.8g} | {100*gci:.3f} % |")
        else:
            print(f"| {q} | {L1}-{L2}-{L4} | {r:.3f} | not monotone | – | – |")
print()
print("| Level | Windows | Coupling its mean | SU2 inner its per solve | Coupled wall (s) | Settled (u_y within 1e-6) from t (s) |")
print("|---|---:|---:|---:|---:|---:|")
for L in Ls:
    r = rows[L]
    print(f"| {L} | {r['windows']} | {float(r['coupling_its_mean']):.2f} | {float(r['inner_its_mean']):.1f} | "
          f"{r.get('wall_coupled_s', '')} | {r['uy_settled_t']} |")
