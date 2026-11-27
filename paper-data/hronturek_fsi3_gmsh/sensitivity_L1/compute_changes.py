#!/usr/bin/env python3
"""Relative changes (to the magnitude of the baseline, as in the paper tables) of the L1 sensitivity runs from the L1 baseline, next to the L1->L2 changes.

Usage: python3 compute_changes.py g1_sens.csv ../g123.csv > changes.csv   (standard library only)
"""
import csv
import sys

QOIS = [
    ("ux_midrange", "u_x (max+min)/2", "m"),
    ("ux_amplitude", "u_x amplitude", "m"),
    ("uy_amplitude", "u_y amplitude", "m"),
    ("frequency", "frequency", "Hz"),
    ("drag_midrange", "drag (max+min)/2", "N/m"),
    ("drag_amplitude", "drag amplitude", "N/m"),
    ("lift_amplitude", "lift amplitude", "N/m"),
]


def rows(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


sens = {r["name"]: r for r in rows(sys.argv[1])}
lev = rows(sys.argv[2])
base, tol, dt = sens["g1"], sens["g1_tol1e-7"], sens["g1_dt1e-4"]
l1, l2 = lev[0], lev[1]
w = csv.writer(sys.stdout, lineterminator="\n")
w.writerow(["qoi", "unit", "baseline_L1", "tol1e-7", "change_tol_pct", "dt1e-4",
            "change_dt_pct", "L2", "change_L1_L2_pct"])
for key, name, unit in QOIS:
    b = float(base[key])
    if abs(b - float(l1[key])) > 1e-12 * abs(b):
        sys.exit(f"baseline {key} differs from the L1 row of {sys.argv[2]}")
    t, d, v2 = float(tol[key]), float(dt[key]), float(l2[key])
    w.writerow([name, unit, f"{b:.8g}", f"{t:.8g}", f"{100*(t-b)/abs(b):+.4f}", f"{d:.8g}",
                f"{100*(d-b)/abs(b):+.4f}", f"{v2:.8g}", f"{100*(v2-b)/abs(b):+.3f}"])
