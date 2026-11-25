#!/usr/bin/env python3
"""Successive changes, observed orders, provisional GCIs and reference differences
for the Turek-Hron FSI3 Gmsh family (levels 1-3), from g123.csv.

    python3 compute_orders.py g123.csv > orders.csv

Rules (paper, Section 4): an observed order is computed only where the three
values are monotone with a decreasing change; no QoI has an order compatible
with the formal order 2, so no Richardson value is given and the GCI
(F_s = 1.25, observed order, centred on level 3) is a provisional indicator.
The refinement ratio is the effective one of the family,
r = (N3/N1)^(1/4) from the fluid cell counts (2-D, so sqrt of the cell ratio
per level); dt is reduced by sqrt(2) per level, so the study is "scaled".
"""
import csv
import math
import sys

CELLS = (96712, 194183, 387219)
R = (CELLS[2] / CELLS[0]) ** 0.25
# Failer & Richter (2020), level 6, dt = 0.001 s; Featflow level 4, dt = 0.00025 s
FR = {"ux_midrange": -2.890, "ux_amplitude": 2.729, "uy_amplitude": 35.03,
      "drag_midrange": 460.5, "drag_amplitude": 27.91, "lift_amplitude": 156.6}
FF = {"ux_midrange": -2.88, "ux_amplitude": 2.72, "uy_amplitude": 34.99, "frequency": 5.46,
      "drag_midrange": 460.5, "drag_amplitude": 27.74, "lift_amplitude": 153.91}
QOI = (("ux_midrange", 1e3, "mm"), ("ux_amplitude", 1e3, "mm"), ("uy_amplitude", 1e3, "mm"),
       ("frequency", 1.0, "Hz"), ("drag_midrange", 1.0, "N/m"), ("drag_amplitude", 1.0, "N/m"),
       ("lift_amplitude", 1.0, "N/m"))


def main() -> int:
    rows = list(csv.DictReader(open(sys.argv[1])))
    w = csv.writer(sys.stdout)
    w.writerow(["qoi", "unit", "L1", "L2", "L3", "change_L1_L2_pct", "change_L2_L3_pct",
                "r_eff", "observed_order", "gci_provisional_pct", "label",
                "L3_vs_FailerRichter_pct", "L3_vs_Featflow_pct"])
    for key, scale, unit in QOI:
        f1, f2, f3 = (float(r[key]) * scale for r in rows)   # coarse to fine
        d12, d23 = f2 - f1, f3 - f2
        c12, c23 = 100 * d12 / abs(f1), 100 * d23 / abs(f2)
        if d12 * d23 > 0 and abs(d23) < abs(d12):
            p = math.log(d12 / d23) / math.log(R)
            gci = 100 * 1.25 * abs(d23 / f3) / (R ** p - 1)
            label = "monotone; order not close to 2: provisional GCI"
            ps, gs = f"{p:.2f}", f"{gci:.3f}"
        else:
            ps, gs, label = "", "", "not monotone: no order"
        fr = f"{100 * (f3 - FR[key]) / abs(FR[key]):+.2f}" if key in FR else ""
        ff = f"{100 * (f3 - FF[key]) / abs(FF[key]):+.2f}" if key in FF else ""
        w.writerow([key, unit, f"{f1:.5g}", f"{f2:.5g}", f"{f3:.5g}", f"{c12:+.3f}", f"{c23:+.3f}",
                    f"{R:.4f}", ps, gs, label, fr, ff])
    return 0


if __name__ == "__main__":
    sys.exit(main())
