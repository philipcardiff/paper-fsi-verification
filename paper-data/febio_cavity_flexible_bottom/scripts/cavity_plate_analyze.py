#!/usr/bin/env python3
import os, os as _os; HERE = _os.path.dirname(_os.path.abspath(__file__)); PKG = _os.path.dirname(HERE)   # package root; run directories of retained HPC outputs: env CAVITY_RUNS (default runs/cavityNext)
"""B0 plate-only checks. usage: cavity_plate_analyze.py <run_dir> <ref_uy>. Prints and writes plate_report.json."""
import sys, json, os, re
d, ref = sys.argv[1], float(sys.argv[2])
def last_block(path):
    rows = []
    for l in open(path):
        if l.startswith('*Step') or l.startswith('*Time') or l.startswith('*Data') or not l.strip(): continue
        rows.append([float(x) for x in l.split()])
    return rows
A = last_block(f'{d}/pointA.txt')           # columns: id x y ux uy per step
uyA = A[-1][-1]; hist = [r[-1] for r in A]
def sumR(f):
    rows = last_block(f); n = len(rows) // max(1, len(A)); last = rows[-n:] if n else rows
    return sum(r[2] for r in last)         # id Rx Ry Rz -> Ry
Ry = sumR(f'{d}/react_left.txt') + sumR(f'{d}/react_right.txt')
meta = json.load(open(f'{d}/model_meta.json'))
Fapplied = meta['p'] * 4.0 * meta['thickness_z']   # projected-area upper bound
log = open(f'{d}/run.log', errors='ignore').read()
rep = dict(uy_A=uyA, ref=ref, err_pct=100 * (uyA - ref) / ref, uy_history=hist, sum_Ry=Ry, applied_p_times_area_undeformed=Fapplied,
           Ry_over_applied=Ry / Fapplied, normal_termination='N O R M A L   T E R M I N A T I O N' in log,
           iters=len(re.findall(r'convergence summary', log)), elapsed=re.findall(r'Total elapsed time[ .:]*([\d:]+)', log))
json.dump(rep, open(f'{d}/plate_report.json', 'w'), indent=1); print(json.dumps({k: v for k, v in rep.items() if k != 'uy_history'}, indent=1))
