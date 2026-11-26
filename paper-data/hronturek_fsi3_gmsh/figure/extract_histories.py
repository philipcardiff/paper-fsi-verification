#!/usr/bin/env python3
"""Extract the late FSI3 histories of the Gmsh levels into small CSV files.

For each case, the tip displacement (solidPointDisplacement_pointDisp.dat) and
the cylinder + flag force (fluid/forces/<t0>/force.dat) are read from every
restart directory of postProcessing, stitched (the last row of every time is
kept, as in ht_coupled_analysis.py), divided by the depth 0.015 m, and cropped
to t_from <= t <= t_end. The CSV columns are

    t_s, ux_m, uy_m, drag_N_per_m, lift_N_per_m

With --check, the benchmark statistics (last five full u_y periods before
t_end, the definitions of ht_coupled_analysis.py) are recomputed from the
cropped data and compared with g123.csv.

    python3 extract_histories.py --case L1=<g1> --case L2=<g2> --case L3=<l3> \
        --out-dir data --check ../g123.csv

The window statistics depend on the data from t_end - 1.5 s onwards (the
crossing detector uses the mean and range of u_y over t >= 4.5 s), so keep
--t-from at 4.5 for the check to reproduce g123.csv.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

DEPTH = 0.015


def rows(path: Path) -> np.ndarray:
    out = []
    for line in path.read_text().splitlines():
        line = line.replace("(", " ").replace(")", " ")
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        try:
            out.append([float(v) for v in line.split()])
        except ValueError:
            continue
    n = min(len(r) for r in out)
    return np.array([r[:n] for r in out])


def unique_time(a: np.ndarray) -> np.ndarray:
    """Keep the last row of every time (restarts append or rewrite)."""
    _, idx = np.unique(a[::-1, 0], return_index=True)
    return a[::-1][idx]


def stitch(files: list[Path]) -> np.ndarray:
    # restart directories are named by their start time: read them in time order
    files = sorted(files, key=lambda p: float(p.parent.name) if p.parent.name.replace(".", "").isdigit() else 0.0)
    return unique_time(np.vstack([rows(p) for p in files]))


def history(case: Path) -> dict:
    pp = case / "postProcessing"
    disp = stitch(list(pp.glob("*/solidPointDisplacement_pointDisp.dat")))
    force = stitch(list(pp.glob("fluid/forces/*/force.dat")))
    return {"t": disp[:, 0], "ux": disp[:, 1], "uy": disp[:, 2],
            "tf": force[:, 0], "drag": force[:, 1] / DEPTH, "lift": force[:, 2] / DEPTH}


def crossings(t: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Upward crossings of the mean (hysteresis 0.25 of the range)."""
    mean = y.mean()
    band = 0.25 * (y.max() - y.min())
    out, armed = [], y[0] < mean - band
    for k in range(1, len(y)):
        if y[k] < mean - band:
            armed = True
        elif armed and y[k] >= mean:
            a, b = y[k - 1] - mean, y[k] - mean
            out.append(t[k - 1] + a / (a - b) * (t[k] - t[k - 1]))
            armed = False
    return np.array(out)


def window_stats(h: dict, periods: int = 5) -> dict:
    cr = crossings(h["t"], h["uy"])
    bounds = cr[-(periods + 1):]
    out = {"window_start": bounds[0], "window_end": bounds[-1],
           "frequency": (len(bounds) - 1) / (bounds[-1] - bounds[0])}
    for q, tkey in (("ux", "t"), ("uy", "t"), ("drag", "tf"), ("lift", "tf")):
        t, y = h[tkey], h[q]
        mid, amp = [], []
        for a, b in zip(bounds[:-1], bounds[1:]):
            seg = y[(t >= a) & (t <= b)]
            mid.append(0.5 * (seg.max() + seg.min()))
            amp.append(0.5 * (seg.max() - seg.min()))
        out[f"{q}_midrange"] = float(np.mean(mid))
        out[f"{q}_amplitude"] = float(np.mean(amp))
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--case", action="append", required=True, help="LABEL=case directory")
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--t-from", type=float, default=4.5)
    ap.add_argument("--t-end", type=float, default=6.0)
    ap.add_argument("--check", type=Path, help="g123.csv to compare the window statistics with")
    a = ap.parse_args()
    a.out_dir.mkdir(parents=True, exist_ok=True)
    ref = list(csv.DictReader(a.check.open())) if a.check else []
    worst = 0.0
    for k, spec in enumerate(a.case):
        label, _, path = spec.partition("=")
        h = history(Path(path).expanduser())
        m = (h["t"] >= a.t_from - 1e-9) & (h["t"] <= a.t_end + 1e-9)
        mf = (h["tf"] >= a.t_from - 1e-9) & (h["tf"] <= a.t_end + 1e-9)
        t, tf = h["t"][m], h["tf"][mf]
        if len(t) != len(tf) or np.max(np.abs(t - tf)) > 1e-7:
            raise SystemExit(f"{label}: displacement and force times differ")
        out = a.out_dir / f"fsi3_gmsh_{label}_histories.csv"
        with out.open("w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["t_s", "ux_m", "uy_m", "drag_N_per_m", "lift_N_per_m"])
            for row in zip(t, h["ux"][m], h["uy"][m], h["drag"][mf], h["lift"][mf]):
                w.writerow([f"{row[0]:.7f}"] + [f"{v:.9g}" for v in row[1:]])
        print(f"{label}: {len(t)} rows, t = {t[0]:.4f}-{t[-1]:.4f} s -> {out}")
        if ref:
            # recompute from the cropped data, as written
            d = np.loadtxt(out, delimiter=",", skiprows=1)
            hh = {"t": d[:, 0], "tf": d[:, 0], "ux": d[:, 1], "uy": d[:, 2], "drag": d[:, 3], "lift": d[:, 4]}
            s = window_stats(hh)
            r = ref[k]
            print(f"  window {s['window_start']:.4f}-{s['window_end']:.4f} s (vs {r['name']})")
            for key, v in s.items():
                if key.startswith("window"):
                    continue
                rv = float(r[key])
                rel = abs(v - rv) / abs(rv)
                worst = max(worst, rel)
                print(f"  {key:16s} {v: .7g}  g123 {rv: .7g}  rel. diff {rel:.1e}")
    if ref:
        print(f"largest relative difference from g123.csv: {worst:.1e}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
