#!/usr/bin/env python3
"""Amplitude of the axis-pressure oscillation near the front: max |p - pbar|,
pbar the centred 0.4 ms moving average, within +-1 ms of the averaged
half-inlet arrival, at the stations z = 10-40 mm.
Usage: osc_amplitude.py [--truncate T] LABEL=PATH ...  (uses cp_averaged.py)"""
import argparse
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cp_averaged as C  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--truncate", type=float, default=None)
ap.add_argument("cases", nargs="+")
a = ap.parse_args()
print("label," + ",".join(f"z{round(C.STATIONS[i]*1e3):02d}mm_Pa" for i in C.FIT))
for spec in a.cases:
    label, path = spec.split("=", 1)
    hist = C.read_kratos(path) if ".csv" in path else C.read_openfoam(path)
    out = []
    for i in C.FIT:
        h = hist[i]
        if a.truncate is not None:
            h = [(t, v) for t, v in h if t <= a.truncate + 1e-12]
        s = C.moving_average(h, 4e-4)
        ta = C.arrival_time(s, 0.5 * C.P_IN)
        out.append(max(abs(p - q) for (t, p), (_, q) in zip(h, s) if abs(t - ta) <= 1e-3))
    print(label + "," + ",".join(f"{x:.1f}" for x in out))
