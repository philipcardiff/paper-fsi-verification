#!/usr/bin/env python3
"""compare_json.py A.json B.json [rtol]: recursive comparison of two JSON files (numbers to relative tolerance rtol with a 1e-300 floor-free absolute floor of rtol*1e-3 of the largest magnitude in A)."""
import json, sys
def flat(x, p=""):
    if isinstance(x, dict):
        for k, v in x.items(): yield from flat(v, p + "/" + k)
    elif isinstance(x, list):
        for i, v in enumerate(x): yield from flat(v, p + f"[{i}]")
    else: yield p, x
a, b = dict(flat(json.load(open(sys.argv[1])))), dict(flat(json.load(open(sys.argv[2])))); rtol = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-9
if set(a) != set(b): print("KEY MISMATCH", sorted(set(a) ^ set(b))[:5]); sys.exit(1)
num = [abs(v) for v in a.values() if isinstance(v, (int, float)) and not isinstance(v, bool)]
floor = rtol * 1e-3 * (max(num) if num else 1.0); worst = 0.0
for k, v in a.items():
    w = b[k]
    if isinstance(v, (int, float)) and not isinstance(v, bool):
        d = abs(v - w) / max(abs(v), abs(w), floor / rtol if rtol else 1e-300); worst = max(worst, d)
        if d > rtol: print("DIFF", k, v, w); sys.exit(1)
    elif v != w: print("DIFF", k, v, w); sys.exit(1)
print(f"identical within rtol {rtol:g} (worst relative difference {worst:.2e}, {len(a)} entries)")
