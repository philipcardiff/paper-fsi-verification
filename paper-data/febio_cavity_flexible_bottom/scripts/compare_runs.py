#!/usr/bin/env python3
import os, os as _os; HERE = _os.path.dirname(_os.path.abspath(__file__)); PKG = _os.path.dirname(HERE)   # package root; run directories of retained HPC outputs: env CAVITY_RUNS (default runs/cavityNext)
"""Compare two FEBio fluid-only runs (stock vs experimental) from their Mon_*.txt node logs and .log files.
usage: compare_runs.py RUN_A RUN_B OUT.json   (A = reference). Relative differences are max|dA|/max|A| per field per time block (A's max over that block)."""
import sys, re, json, numpy as np, os
def read_blocks(path):
    out = {}; t = None; rows = []
    for ln in open(path):
        if ln.startswith("*Time"):
            if t is not None: out[t] = np.array(rows)
            t = float(ln.split("=")[1]); rows = []
        elif ln.startswith("*"): continue
        else:
            v = ln.split()
            if len(v) >= 6: rows.append([float(x) for x in v])
    if t is not None: out[t] = np.array(rows)
    return out
def iters(log):
    its = [int(m.group(1)) for m in re.finditer(r"number of iterations\s*:\s*(\d+)", open(log, errors="ignore").read())]
    return its
def total_time(log):
    m = re.search(r"Total elapsed time[ .]*:\s*([\d:]+)\s*\(([\d.]+) sec\)", open(log, errors="ignore").read()); return float(m.group(2)) if m else None
def term(log):
    s = open(log, errors="ignore").read(); return ("N O R M A L   T E R M I N A T I O N" in s, "E R R O R   T E R M I N A T I O N" in s)
if __name__ == "__main__":
    A, B, out = sys.argv[1:4]; res = {}
    la = [f for f in os.listdir(A) if f.endswith(".log") and f != "run.log"][0]; lb = [f for f in os.listdir(B) if f.endswith(".log") and f != "run.log"][0]
    ia, ib = iters(os.path.join(A, la)), iters(os.path.join(B, lb))
    res["iterations"] = dict(A_total=sum(ia), B_total=sum(ib), A_steps=len(ia), B_steps=len(ib), identical_sequence=(ia == ib), max_abs_step_diff=int(max(abs(np.array(ia[:min(len(ia),len(ib))])-np.array(ib[:min(len(ia),len(ib))])))) if ia and ib else None)
    res["wall_s"] = dict(A=total_time(os.path.join(A, la)), B=total_time(os.path.join(B, lb))); res["termination_A_normal_error"] = term(os.path.join(A, la)); res["termination_B_normal_error"] = term(os.path.join(B, lb))
    fields = {}
    for mon in sorted(f for f in os.listdir(A) if f.startswith("Mon_") and os.path.exists(os.path.join(B, f))):
        ba, bb = read_blocks(os.path.join(A, mon)), read_blocks(os.path.join(B, mon)); worst = {"wx": 0, "wy": 0, "ef": 0}; bit = True; nb = 0
        for t in ba:
            if t not in bb: continue
            a, b = ba[t], bb[t]
            if a.ndim < 2 or b.ndim < 2 or a.shape[1] < 6: continue
            nb += 1
            if a.shape != b.shape: bit = False; continue
            bit &= bool(np.array_equal(a, b))
            vsc = max(np.abs(a[:, 3]).max(), np.abs(a[:, 4]).max(), 1e-300)   # velocity scale of the block (both components), so that tiny wy values are not divided by themselves
            for k, name in ((3, "wx"), (4, "wy"), (5, "ef")):
                sc = vsc if name != "ef" else max(np.abs(a[:, k]).max(), 1e-300); worst[name] = max(worst[name], float(np.abs(a[:, k] - b[:, k]).max() / sc))
        fields[mon] = dict(blocks=nb, bitwise_identical=bit, worst_rel_diff=worst)
    res["fields"] = fields; json.dump(res, open(out, "w"), indent=1); print(json.dumps(res, indent=1))
