#!/usr/bin/env python3
"""Convergence analysis of the axisymmetric 3dTube levels, and comparison
with the 3-D quarter-tube levels L1-L3 of Table tab:3dtube.

Usage: analyse.py qoi_L1.csv qoi_L2.csv ... [--md out.md]
(each CSV as written by post_axi.py; rows are sorted by level)
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

QOIS = [
    ("ur_max_A_m", "u_r,max (mm)", 1e3, 6),
    ("t_ur_max_A_s", "t(u_r,max) (ms)", 1e3, 4),
    ("uz_min_A_m", "u_z,min (mm)", 1e3, 6),
    ("t_arrival_A_s", "t_arr (ms)", 1e3, 4),
    ("wave_speed_pressure_m_s", "c_p (m/s)", 1.0, 4),
    ("wave_speed_wall_m_s", "c_wall (m/s)", 1.0, 4),
    ("ur_min_late_A_m", "u_r,min late (mm)", 1e3, 6),
]

# 3-D quarter tube, corrected meshes (paper Table tab:3dtube; driver CSVs of
# tube3d/post671 for L1/L2 and paper-data/3dtube_level3_corrected_ecore for L3)
THREE_D = {
    1: dict(ur_max_A_m=0.00015958314062500001, t_ur_max_A_s=0.007201171875,
            uz_min_A_m=-8.817142461612284e-05, t_arrival_A_s=0.0058715997190624085,
            wave_speed_pressure_m_s=4.663896331674805,
            wave_speed_wall_m_s=5.067093004807056,
            ur_min_late_A_m=-9.9972344687017e-05),
    2: dict(ur_max_A_m=0.0001592429642857143, t_ur_max_A_s=0.007215178571428573,
            uz_min_A_m=-8.656945698398577e-05, t_arrival_A_s=0.005841378282759041,
            wave_speed_pressure_m_s=4.651108565828932,
            wave_speed_wall_m_s=5.138673154479264,
            ur_min_late_A_m=-9.397018775641026e-05),
    3: dict(ur_max_A_m=0.00015830720000000001, t_ur_max_A_s=0.007251250000000004,
            uz_min_A_m=-8.569550065789474e-05, t_arrival_A_s=0.005823331463825014,
            wave_speed_pressure_m_s=4.625250577609168,
            wave_speed_wall_m_s=5.0669531190767065,
            ur_min_late_A_m=-9.5913825e-05),
}


def order(f1, f2, f3):
    a, b = f1 - f2, f2 - f3
    if a == 0 or b == 0 or a * b < 0:
        return math.nan
    return math.log(abs(a / b), 2.0)


def richardson(f2, f3, p):
    if not math.isfinite(p) or p <= 0:
        return math.nan
    return f3 + (f3 - f2) / (2.0 ** p - 1.0)


def fmt(x, d):
    return "–" if not math.isfinite(x) else f"{x:.{d}f}"


def pct(x):
    return "–" if not math.isfinite(x) else f"{100 * x:+.3f} %"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("csvs", nargs="+", type=Path)
    ap.add_argument("--md", type=Path)
    ap.add_argument("--label", default="Hookean")
    args = ap.parse_args()
    rows = []
    for p in args.csvs:
        rows += list(csv.DictReader(p.open()))
    rows.sort(key=lambda r: int(r["level"]))
    levels = [int(r["level"]) for r in rows]
    out = []
    w = out.append
    w(f"### Axisymmetric levels ({args.label})\n")
    w("| Level | Fluid / solid cells | Δt (s) | mean coupl. it. | clock (s) | cores |")
    w("|---:|---:|---:|---:|---:|---:|")
    for r in rows:
        w(f"| {r['level']} | {int(r['fluid_cells']):,} / {int(r['solid_cells']):,} | "
          f"{float(r['delta_t_s']):.4g} | {float(r.get('mean_fsi_iterations') or 'nan'):.2f} | "
          f"{float(r.get('clock_time_s') or 'nan'):.0f} | {r['cores']} |")
    w("")
    hdr = "| Quantity | " + " | ".join(f"L{l}" for l in levels) + " |"
    w(hdr)
    w("|---|" + "---:|" * len(levels))
    vals = {}
    for key, name, scale, d in QOIS:
        v = [float(r[key]) for r in rows]
        vals[key] = v
        w(f"| {name} | " + " | ".join(fmt(x * scale, d) for x in v) + " |")
    w("")
    w("Successive relative changes (f_{k+1} − f_k)/|f_k|:\n")
    w("| Quantity | " + " | ".join(f"L{a}→L{b}" for a, b in zip(levels, levels[1:])) + " |")
    w("|---|" + "---:|" * (len(levels) - 1))
    for key, name, scale, d in QOIS:
        v = vals[key]
        w(f"| {name} | " + " | ".join(pct((b - a) / abs(a)) for a, b in zip(v, v[1:])) + " |")
    w("")
    w("Observed orders p = log2((f_{k}−f_{k+1})/(f_{k+1}−f_{k+2})) and Richardson "
      "limits f_∞ = f_{k+2} + (f_{k+2}−f_{k+1})/(2^p−1) from each triple:\n")
    trip = list(zip(levels, levels[1:], levels[2:]))
    w("| Quantity | " + " | ".join(f"p (L{a}–L{c})" for a, b, c in trip) + " | "
      + " | ".join(f"f∞ (L{a}–L{c})" for a, b, c in trip) + " |")
    w("|---|" + "---:|" * (2 * len(trip)))
    for key, name, scale, d in QOIS:
        v = vals[key]
        ps, lims = [], []
        for i in range(len(v) - 2):
            p = order(*v[i:i + 3])
            ps.append(p)
            lims.append(richardson(v[i + 1], v[i + 2], p))
        w(f"| {name} | " + " | ".join(fmt(p, 2) for p in ps) + " | "
          + " | ".join(fmt(x * scale, d) for x in lims) + " |")
    w("")
    # Error of each level against the finest-triple Richardson limit
    w("Relative difference from the finest-triple Richardson limit (where it exists) "
      "and, for comparison, the 3-D values at the same level:\n")
    w("| Quantity | " + " | ".join(f"axi L{l}" for l in levels)
      + " | " + " | ".join(f"3-D L{l}" for l in sorted(THREE_D)) + " |")
    w("|---|" + "---:|" * (len(levels) + len(THREE_D)))
    for key, name, scale, d in QOIS:
        v = vals[key]
        if len(v) >= 3:
            p = order(*v[-3:])
            lim = richardson(v[-2], v[-1], p)
        else:
            lim = math.nan
        w(f"| {name} | " + " | ".join(pct((x - lim) / abs(lim)) for x in v) + " | "
          + " | ".join(pct((THREE_D[l][key] - lim) / abs(lim)) for l in sorted(THREE_D)) + " |")
    w("")
    w("3-D minus axisymmetric at matching levels, relative to the axisymmetric value:\n")
    common = [l for l in sorted(THREE_D) if l in levels]
    w("| Quantity | " + " | ".join(f"L{l}" for l in common) + " |")
    w("|---|" + "---:|" * len(common))
    for key, name, scale, d in QOIS:
        v = dict(zip(levels, vals[key]))
        w(f"| {name} | " + " | ".join(pct((THREE_D[l][key] - v[l]) / abs(v[l])) for l in common) + " |")
    text = "\n".join(out) + "\n"
    print(text)
    if args.md:
        args.md.write_text(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
