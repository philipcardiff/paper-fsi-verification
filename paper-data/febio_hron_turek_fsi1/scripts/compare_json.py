#!/usr/bin/env python3
"""compare_json.py A.json B.json [rel_tol]: recursively compares numeric leaves (keys present in both; strings ignored); exit 1 on mismatch."""
import sys, json, math
a, b = json.load(open(sys.argv[1])), json.load(open(sys.argv[2])); tol = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-9; bad = []; n = 0
def walk(x, y, path):
    global n
    if isinstance(x, dict) and isinstance(y, dict):
        for k in x:
            if k in y: walk(x[k], y[k], path + "/" + k)
    elif isinstance(x, list) and isinstance(y, list) and len(x) == len(y):
        for i, (u, v) in enumerate(zip(x, y)): walk(u, v, f"{path}[{i}]")
    elif isinstance(x, (int, float)) and isinstance(y, (int, float)) and not isinstance(x, bool):
        n += 1
        if abs(x - y) > tol * max(1.0, abs(x), abs(y)): bad.append((path, x, y))
walk(a, b, "")
print(f"{sys.argv[1].split('/')[-1]}: {n} numeric values compared, {len(bad)} mismatches"); [print("  ", m) for m in bad[:5]]; sys.exit(1 if bad else 0)
