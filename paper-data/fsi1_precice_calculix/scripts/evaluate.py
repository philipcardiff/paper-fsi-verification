#!/usr/bin/env python3
"""Steady FSI1 values of one case against the Featflow 7+0 reference.

    evaluate.py CASEDIR

Forces: last (converged) line per time of postProcessing/forces/*/force.dat,
divided by the 0.015 m depth. Point A: preCICE watch-point on the solid mesh.
Steady state: relative spread over the closing 20 % of the run (same test as the
solids4foam FSI1 verification, tolerance 1e-3).
"""
import glob, sys
from pathlib import Path

case = Path(sys.argv[1]); depth = 0.015; tEnd = 30.0; t0 = 2.0
ref = {"ux": 2.270493e-5, "uy": 8.208773e-4, "drag": 14.29426, "lift": 0.7637460}

forces = {}
for f in sorted(glob.glob(str(case / "fluid-openfoam/postProcessing/forces/*/force.dat"))):
    for line in open(f):
        if line.startswith("#") or not line.strip():
            continue
        c = line.split()
        forces[round(float(c[0]), 8)] = (float(c[1]) / depth, float(c[2]) / depth)
disp = {}
for line in open(case / "solid-calculix/precice-Solid-watchpoint-PointA.log"):
    c = line.split()
    if c[0] == "Time":
        continue
    disp[round(float(c[0]) + t0, 8)] = (float(c[5]), float(c[6]))

hist = {"ux": {t: v[0] for t, v in disp.items()}, "uy": {t: v[1] for t, v in disp.items()},
        "drag": {t: v[0] for t, v in forces.items()}, "lift": {t: v[1] for t, v in forces.items()}}
tStart = tEnd - 0.2 * (tEnd - t0)
print(f"case {case}, last time: fluid {max(forces)}, solid {max(disp)}")
print(f"{'qty':5} {'value':>13} {'Featflow 7+0':>13} {'error':>8} {'spread(last 20%)':>17}")
out = {}
for q, h in hist.items():
    v = h[max(h)]
    win = [x for t, x in h.items() if t >= tStart - 1e-9]
    spread = (max(win) - min(win)) / abs(v)
    out[q] = v
    print(f"{q:5} {v:13.6e} {ref[q]:13.6e} {100*(v-ref[q])/ref[q]:7.2f}% {spread:17.1e}"
          + ("" if spread < 1e-3 else "  NOT STEADY"))

its = [int(l.split()[2]) for l in open(case / "solid-calculix/precice-Solid-iterations.log")
       if l.split()[0].isdigit()]
print(f"coupled windows {len(its)}, coupling iterations mean {sum(its)/len(its):.2f}, max {max(its)}")
