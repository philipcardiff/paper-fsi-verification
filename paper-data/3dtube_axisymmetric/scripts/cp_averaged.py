#!/usr/bin/env python3
"""Pressure-front speed c_p of the pressure-pulse tube from axis-pressure probes.

c_p is the least-squares slope dz/dt of the first time at which the axis
pressure rises through half the inlet pressure (1333.2 Pa / 2), at the stations
z = 10-40 mm.  With --window W > 0 the probe signal is first replaced by its
centred moving average over W seconds (all samples with |t' - t| <= W/2; the
window is truncated at the ends of the record).  The paper's c_p uses
W = 0.4 ms, two periods of the ~0.2 ms axis-pressure oscillation; W = 0 is the
raw definition of the 3dTube verification driver.

Inputs (one or more, any mix):
  * an OpenFOAM probes file  postProcessing/verificationAxisPressure/fluid/0/p
    (optionally xz-compressed, *.xz)
    (kinematic pressure; multiplied by rho = 1000 kg/m^3), nine stations
    z = 5, 10, ..., 45 mm;
  * a Kratos history.csv (columns p_z050 ... p_z450, in Pa; *.csv or *.csv.xz).

Usage: cp_averaged.py [--window W ...] [--truncate T] LABEL=PATH ...
Prints a CSV: label, window_s, c_p_m_s, then the nine station arrivals (s).
Only the Python standard library is used.
"""
import argparse
import bisect
import csv
import lzma
import math
import sys

STATIONS = [0.005 * k for k in range(1, 10)]
FIT = [i for i, z in enumerate(STATIONS) if 0.01 - 1e-9 <= z <= 0.04 + 1e-9]
P_IN = 1333.2
RHO = 1000.0


def _open(path):
    return lzma.open(path, "rt") if path.endswith(".xz") else open(path)


def read_openfoam(path):
    rows = {}
    with _open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            v = [float(x) for x in line.split()]
            if len(v) != len(STATIONS) + 1:
                raise SystemExit(f"{path}: unexpected column count {len(v)}")
            rows[v[0]] = v[1:]          # last row per time wins
    t = sorted(rows)
    # Probe times are written to 6 significant digits; snap them to the
    # time-step grid (dt from the early, exactly written times).
    dt = sorted(b - a for a, b in zip(t[:50], t[1:51]))[24]
    t_snap = [round(x / dt) * dt for x in t]
    return [[(ts, RHO * rows[x][i]) for ts, x in zip(t_snap, t)]
            for i in range(len(STATIONS))]


def read_kratos(path):
    with _open(path) as f:
        rows = list(csv.DictReader(f))
    t = [float(r["time"]) for r in rows]
    out = []
    for z in STATIONS:
        name = f"p_z{round(z * 1e4):03d}"
        out.append(list(zip(t, [float(r[name]) for r in rows])))
    return out


def moving_average(h, w):
    if w <= 0.0:
        return h
    t = [a for a, _ in h]
    cs = [0.0]
    for _, v in h:
        cs.append(cs[-1] + v)
    eps = 1e-9 * w
    out = []
    for ti in t:
        lo = bisect.bisect_left(t, ti - 0.5 * w - eps)
        hi = bisect.bisect_right(t, ti + 0.5 * w + eps)
        out.append((ti, (cs[hi] - cs[lo]) / (hi - lo)))
    return out


def arrival_time(h, level):
    """First upward crossing of level, linearly interpolated (as the driver)."""
    for (t0, v0), (t1, v1) in zip(h, h[1:]):
        if v0 < level <= v1:
            return t0 + (level - v0) / (v1 - v0) * (t1 - t0)
    return math.nan


def fit_speed(z, t):
    if len(t) < 2 or not all(math.isfinite(x) for x in t):
        return math.nan
    mz, mt = sum(z) / len(z), sum(t) / len(t)
    stt = sum((x - mt) ** 2 for x in t)
    stz = sum((x - mt) * (y - mz) for x, y in zip(t, z))
    return stz / stt


def c_p(hist, window, truncate=None):
    arrivals = []
    for h in hist:
        if truncate is not None:
            h = [(a, b) for a, b in h if a <= truncate + 1e-12]
        arrivals.append(arrival_time(moving_average(h, window), 0.5 * P_IN))
    return fit_speed([STATIONS[i] for i in FIT],
                     [arrivals[i] for i in FIT]), arrivals


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--window", type=float, action="append",
                    help="averaging window in s (repeatable; default 0 and 4e-4)")
    ap.add_argument("--truncate", type=float, default=None,
                    help="drop samples after this time (s) before averaging")
    ap.add_argument("cases", nargs="+", metavar="LABEL=PATH")
    a = ap.parse_args()
    windows = a.window or [0.0, 4e-4]
    w = csv.writer(sys.stdout)
    w.writerow(["label", "window_s", "c_p_m_s"]
               + [f"arrival_z{round(z * 1e3):02d}mm_s" for z in STATIONS])
    for spec in a.cases:
        label, path = spec.split("=", 1)
        hist = read_kratos(path) if ".csv" in path else read_openfoam(path)
        for win in windows:
            cp, arr = c_p(hist, win, a.truncate)
            w.writerow([label, f"{win:g}", f"{cp:.6f}"] + [f"{x:.8e}" for x in arr])


if __name__ == "__main__":
    main()
