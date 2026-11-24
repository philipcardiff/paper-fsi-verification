#!/usr/bin/env python3
"""Temporal orders on a fixed mesh from differences of successive time-step solutions.

Usage: timeOrders.py run_dt1 run_dt2 run_dt3 ... [--time 0.4]   (dt halving)

For successive runs the differences d_k = || f(dt_k) - f(dt_{k+1}) || are
computed cell-wise (volume-weighted RMS):
  - solid D and U on the fixed reference solid mesh (exact comparison);
  - solid interface traction on the solid interface faces;
  - fluid U and p cell-wise. The fluid meshes of the runs differ by the mesh
    motion, so with --exact the exact solution is first subtracted at each
    run's own cell centres (written by postProcess -func writeCellCentres
    with the environment --env or $MMS_FSI_ENV); without it the fluid
    differences are only approximate.
Observed order = log2(d_k / d_{k+1}). The differences remove the spatial
error, which is common to all runs on the same mesh.
"""
import argparse
import math
import os
import re

import numpy as np

import foamio

ap = argparse.ArgumentParser()
ap.add_argument("runs", nargs="+")
ap.add_argument("--time", default="0.4")
ap.add_argument("--env", default=os.environ.get("MMS_FSI_ENV"),
                help="OpenFOAM environment for writeCellCentres (default $MMS_FSI_ENV)")
ap.add_argument("--exact", default=None,
                help="evalExact executable: subtract the exact fluid solution at each run's own cell "
                     "centres before differencing (the fluid meshes of the runs differ)")
a = ap.parse_args()


def entry(path, patch, key):
    txt = open(path).read()
    b = txt.index("boundaryField")
    k = re.compile(r"\b" + patch + r"\b\s*\{").search(txt, b).end()
    seg = txt[k:]
    m = re.search(r"\b" + key + r"\s+nonuniform\s+List<\w+>\s*", seg)
    v, _ = foamio._parse_list(seg, m.end())
    return v


def load(run):
    t = os.path.join(run, a.time)
    out = {}
    out["solid D"] = foamio.read_field(os.path.join(t, "solid", "D"))[:, :2]
    out["solid U"] = foamio.read_field(os.path.join(t, "solid", "U"))[:, :2]
    out["interface traction"] = entry(os.path.join(t, "solid", "D"), "interface", "traction")[:, :2]
    out["fluid U"] = foamio.read_field(os.path.join(t, "fluid", "U"))[:, :2]
    out["fluid p"] = foamio.read_field(os.path.join(t, "fluid", "p"))
    if a.exact:
        import subprocess
        cpath = os.path.join(t, "fluid", "C")
        if not os.path.exists(cpath):
            subprocess.run("bash -lc 'cd %s && source %s >/dev/null 2>&1; postProcess -func writeCellCentres "
                           "-region fluid -time %s > log.cc 2>&1'" % (run, a.env, a.time), shell=True)
        C = foamio.read_field(cpath)
        ex = foamio.exact("fluid", [(c[0], c[1], float(a.time)) for c in C], a.exact)
        out["fluid U"] = out["fluid U"] - ex[:, :2]
        out["fluid p"] = out["fluid p"] - ex[:, 2]
    return out


data = [load(r) for r in a.runs]
names = list(data[0].keys())
print("| quantity | " + " | ".join(f"‖{os.path.basename(a.runs[i])}−{os.path.basename(a.runs[i+1])}‖"
                                    for i in range(len(a.runs) - 1)) + " | orders |")
print("|---|" + "---|" * (len(a.runs) - 1) + "---|")
for nm in names:
    diffs = []
    for i in range(len(data) - 1):
        d = data[i][nm] - data[i + 1][nm]
        d = np.sqrt(np.mean(d ** 2 if d.ndim == 1 else np.sum(d ** 2, axis=1)))
        diffs.append(d)
    ords = [f"{math.log(diffs[i] / diffs[i + 1]) / math.log(2):.2f}" for i in range(len(diffs) - 1)]
    print(f"| {nm} | " + " | ".join(f"{d:.2e}" for d in diffs) + " | " + ", ".join(ords) + " |")
