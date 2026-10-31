#!/usr/bin/env python3
"""Compare two JSON files recursively (numbers with relative tolerance).  usage: compare_json.py A.json B.json REL_TOL [ignored_key ...]   exit 1 if different."""
import sys, json, math
a, b, tol, ign = json.load(open(sys.argv[1])), json.load(open(sys.argv[2])), float(sys.argv[3]), set(sys.argv[4:]); bad = []
def cmp(x, y, p):
    if isinstance(x, dict) and isinstance(y, dict):
        for k in sorted(set(x) | set(y)):
            if k in ign: continue
            if k not in x or k not in y: bad.append(f"{p}/{k} missing"); continue
            cmp(x[k], y[k], f"{p}/{k}")
    elif isinstance(x, (list, tuple)) and isinstance(y, (list, tuple)):
        if len(x) != len(y): bad.append(f"{p} length {len(x)} vs {len(y)}"); return
        for i, (u, v) in enumerate(zip(x, y)): cmp(u, v, f"{p}[{i}]")
    elif isinstance(x, (int, float)) and isinstance(y, (int, float)) and not isinstance(x, bool):
        if not (x == y or abs(x - y) <= tol * max(abs(x), abs(y), 1e-300) or (math.isnan(x) and math.isnan(y))): bad.append(f"{p}: {x} vs {y}")
    elif x != y: bad.append(f"{p}: {x!r} vs {y!r}")
cmp(a, b, "")
print(f"compare_json {sys.argv[1].split('/')[-1]}: " + ("IDENTICAL within tolerance" if not bad else f"{len(bad)} DIFFERENCES, first: {bad[0]}")); sys.exit(1 if bad else 0)
