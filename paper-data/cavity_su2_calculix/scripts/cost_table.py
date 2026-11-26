#!/usr/bin/env python3
"""Cost table from the run logs: Slurm job, SU2 ranks, stage-1 (steady rigid plate) and coupled wall
clock, core-hours of the coupled stage (all allocated cores), coupling and SU2 inner iterations.

    cost_table.py LOGDIR [LOGDIR ...]     (each LOGDIR as logs/<run>: slurm.out, walltime_*.txt, solid-calculix/...)
"""
import csv
import re
import sys
from pathlib import Path

print("| Run | Slurm job | Cores (SU2 + CalculiX) | Stage 1 wall (s) | Coupled wall (s) | Core-hours (coupled) | Windows | Mean / max coupling its | Windows at cap | SU2 inner its per solve |")
print("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
for d in map(Path, sys.argv[1:]):
    s = (d / "slurm.out").read_text()
    job = re.search(r"^Job (\d+)", s, re.M).group(1)
    ntasks = int(re.search(r"tasks (\d+)", s).group(1))
    nf = int(re.search(r"SU2 ranks (\d+)", s).group(1))
    st = (d / "walltime_steady.txt").read_text().strip() if (d / "walltime_steady.txt").exists() else "†"
    wc = int((d / "walltime_coupled.txt").read_text())
    its = [int(l.split()[2]) for l in open(d / "solid-calculix/precice-Solid-iterations.log") if l.split() and l.split()[0].isdigit()]
    h = sorted((d / "fluid-su2").glob("history_0*.csv"))[-1]
    inner = [int(r[next(k for k in r if k.strip().strip('"') == "Inner_Iter")]) + 1 for r in csv.DictReader(open(h))]
    print(f"| {d.name} | {job} | {ntasks} ({nf} + 1) | {st} | {wc} | {wc * ntasks / 3600:.2f} | {len(its)} | "
          f"{sum(its) / len(its):.2f} / {max(its)} | {sum(i >= 50 for i in its)} | {sum(inner) / len(inner):.1f} |")
