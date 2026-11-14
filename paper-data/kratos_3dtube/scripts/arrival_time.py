#!/usr/bin/env python3
"""Arrival-time comparison at point A and along the tube: Kratos BDF2 levels 1 and 2 against the
solids4foam backward mesh study, and Kratos level 1 against a solids4foam level-1 rerun with full
histories.

All front rules are the solids4foam driver's (front_arrival, first_peak, arrival_time, peak), imported
from tube_3d_verification.py. The reference JSON is found next to the driver (reference/) or one level
up (../reference/, the layout of the solids4foam repository).

Usage:
  python arrival_time.py <tube_3d_verification.py> <s4f mesh-study csv> <s4f level-1 rerun dir>
                         <Kratos level-1 run dir> <Kratos level-2 run dir>
                         [--out report.md] [--csv pointA_level1.csv]
"""
import csv
import importlib.util
import json
import os
import sys

Z_MM = [5, 10, 15, 20, 25, 30, 35, 40, 45]
FRACTIONS = [0.1, 0.25, 0.5, 0.75, 0.9]


def load_driver(path):
    spec = importlib.util.spec_from_file_location("tubedriver", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    here = os.path.dirname(os.path.abspath(path))
    for ref_dir in (os.path.join(here, "reference"), os.path.join(here, "..", "reference")):
        f = os.path.join(ref_dir, "3dTube_verification_references.json")
        if os.path.exists(f):
            return mod, json.load(open(f)), ref_dir
    raise SystemExit("3dTube_verification_references.json not found next to the driver")


def kratos(run_dir, end=None):
    rows = list(csv.DictReader(open(os.path.join(run_dir, "history.csv"))))
    if end is not None:
        rows = [r for r in rows if float(r["time"]) <= end + 1e-12]
    t = [float(r["time"]) for r in rows]
    wall = {z: [(0.0, 0.0)] + list(zip(t, [float(r[f"ur_z{10 * z:03d}"]) for r in rows])) for z in Z_MM}
    pres = {z: [(0.0, 0.0)] + list(zip(t, [float(r[f"p_z{10 * z:03d}"]) for r in rows])) for z in Z_MM}
    return t, wall, pres


def interp(h, t):
    for (t0, v0), (t1, v1) in zip(h, h[1:]):
        if t0 <= t <= t1:
            return v0 + (v1 - v0) * (t - t0) / (t1 - t0)
    return float("nan")


def best_shift(a, b, t0, t1):
    """Shift s (s) minimising the sum of (a(t) - b(t - s))^2 over [t0, t1], on a 2.5e-6 s grid."""
    ts = [t0 + i * (t1 - t0) / 400 for i in range(401)]
    best = None
    for k in range(-80, 81):
        s = k * 2.5e-6
        e = sum((interp(a, t) - interp(b, t - s)) ** 2 for t in ts)
        if best is None or e < best[1]:
            best = (s, e)
    return best[0]


def speed(z0, a0, z1, a1):
    return (z1 - z0) * 1e-3 / (a1 - a0)


def main():
    argv = sys.argv[1:]
    opts = {}
    for key in ("--out", "--csv"):
        if key in argv:
            i = argv.index(key)
            opts[key] = argv[i + 1]
            del argv[i:i + 2]
    drv, ref, ref_dir = load_driver(argv[0])
    study_csv, rerun, k1_dir, k2_dir = argv[1:5]
    pin = float(ref["inlet_pressure_Pa"])
    s4f = {int(r["level"]): r for r in csv.DictReader(open(study_csv)) if r["time_scheme"] == "backward"}
    L = ["# 3dTube arrival times, Kratos against solids4foam\n",
         "Front rules of the solids4foam driver: wall arrival is the first upward crossing of 50% of the "
         "station's incident peak (the first local maximum above 50% of the history maximum); pressure "
         "arrival is the first crossing of 50% of the inlet pressure.\n"]

    # 1. Wall arrivals per station against the solids4foam mesh study
    arrivals = {}
    for lev, run in ((1, k1_dir), (2, k2_dir)):
        _, wall, _ = kratos(run)
        arrivals[lev] = {z: drv.front_arrival(wall[z], ref) for z in Z_MM}
    L.append("## Wall arrival per station (ms), matched levels\n")
    L.append("| z (mm) | Kratos L1 | s4f L1 | diff | Kratos L2 | s4f L2 | diff | s4f L3 |")
    L.append("|---|---|---|---|---|---|---|---|")
    for z in Z_MM:
        k = f"wall_arrival_z{z:02d}mm_s"
        a1, a2 = arrivals[1][z][0], arrivals[2][z][0]
        s1, s2, s3 = (float(s4f[l][k]) for l in (1, 2, 3))
        L.append(f"| {z} | {1e3 * a1:.3f} | {1e3 * s1:.3f} | {1e3 * (a1 - s1):+.3f} | {1e3 * a2:.3f} | "
                 f"{1e3 * s2:.3f} | {1e3 * (a2 - s2):+.3f} | {1e3 * s3:.3f} |")
    L.append("\nBeyond 30 mm the incident front merges with the reflection from the clamped outlet and the "
             "first-peak rule can pick different features in the two codes.\n")

    # 2. Wall-front speed between 5 and 25 mm
    L.append("## Wall-front speed between z = 5 and 25 mm (m/s)\n")
    L.append("| | level 1 | level 2 | level 3 |")
    L.append("|---|---|---|---|")
    ks = [speed(5, arrivals[l][5][0], 25, arrivals[l][25][0]) for l in (1, 2)]
    ss = [speed(5, float(s4f[l]["wall_arrival_z05mm_s"]), 25, float(s4f[l]["wall_arrival_z25mm_s"]))
          for l in (1, 2, 3)]
    L.append(f"| Kratos | {ks[0]:.3f} | {ks[1]:.3f} | – |")
    L.append(f"| solids4foam | {ss[0]:.3f} | {ss[1]:.3f} | {ss[2]:.3f} |")
    L.append(f"| Kratos / solids4foam − 1 | {100 * (ks[0] / ss[0] - 1):+.2f}% | "
             f"{100 * (ks[1] / ss[1] - 1):+.2f}% | – |\n")

    # 3. Level 1 against the solids4foam rerun, full histories
    end = 0.01
    t, wall, pres = kratos(k1_dir, end)
    s_wall = {z: [(tt, r) for tt, r, _ in drv.displacement_history(_Path(rerun), drv.station_name(z * 1e-3), end)]
              for z in Z_MM[:6]}
    s_p = drv.pressure_histories(_Path(rerun), len(ref["stations_z_m"]), end, float(ref["fluid"]["rho_kg_m3"]))
    L.append("## Level 1 against the solids4foam level-1 rerun (full histories to 10 ms)\n")
    L.append("Crossing-time difference Kratos − solids4foam (ms) at fractions of each code's own incident "
             "peak; the least-squares time shift that aligns the two incident fronts over the 3 ms before "
             "the solids4foam peak; and the pressure arrival (the solids4foam axis probes take the value of "
             "the containing cell, so they carry up to half a cell of position error).\n")
    L.append("| z (mm) | " + " | ".join(f"f = {f}" for f in FRACTIONS)
             + " | peak K / s4f (µm) | front shift (ms) | p arrival K / s4f (ms) |")
    L.append("|---" * (len(FRACTIONS) + 4) + "|")
    for i, z in enumerate(Z_MM[:6]):
        diffs = []
        for f in FRACTIONS:
            c = []
            for h in (wall[z], s_wall[z]):
                pt, pv = drv.first_peak(h, 0.5)
                c.append(drv.arrival_time([(a, b) for a, b in h if a <= pt], f * pv))
            diffs.append(1e3 * (c[0] - c[1]))
        pk = [drv.first_peak(h, 0.5) for h in (wall[z], s_wall[z])]
        sh = best_shift(wall[z], s_wall[z], max(0.0, pk[1][0] - 3e-3), pk[1][0])
        pa = [drv.arrival_time(h, 0.5 * pin) for h in (pres[z], s_p[i])]
        L.append(f"| {z} | " + " | ".join(f"{d:+.3f}" for d in diffs)
                 + f" | {1e6 * pk[0][1]:.2f} / {1e6 * pk[1][1]:.2f} | {1e3 * sh:+.3f} | "
                 f"{1e3 * pa[0]:.3f} / {1e3 * pa[1]:.3f} |")

    # 4. Published curve with the same rule
    tk = []
    for line in open(os.path.join(ref_dir, "Tukovic2018_fig25a_pointA.csv")):
        x = line.strip().split(",")
        if line.startswith("#") or x[0] == "time_s" or x[1] == "nan":
            continue
        tk.append((float(x[0]), float(x[1])))
    L.append(f"\nTuković et al. (2018), Fig. 25(a), same rule: t_arr(A) = "
             f"{1e3 * drv.front_arrival(tk, ref)[0]:.3f} ms.\n")

    text = "\n".join(L) + "\n"
    if "--out" in opts:
        open(opts["--out"], "w").write(text)
    print(text)
    if "--csv" in opts:
        with open(opts["--csv"], "w") as fo:
            fo.write("time_s,ur_kratos_m,ur_s4f_m,p_kratos_Pa,p_s4f_Pa\n")
            for tt in t:
                fo.write(f"{tt:.8g},{interp(wall[25], tt):.8e},{interp(s_wall[25], tt):.8e},"
                         f"{interp(pres[25], tt):.6g},{interp(s_p[4], tt):.6g}\n")


def _Path(p):
    from pathlib import Path
    return Path(p)


if __name__ == "__main__":
    main()
