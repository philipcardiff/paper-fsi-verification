import sys
from pathlib import Path
import numpy as np
def rows(p):
    out = {}
    for l in Path(p).read_text().splitlines():
        if not l.strip() or l.lstrip().startswith("#"): continue
        v = [float(x) for x in l.replace("(", " ").replace(")", " ").split()]
        out[v[0]] = np.array(v[1:])
    return out
a, b = sys.argv[1], sys.argv[2]
for f in ["postProcessing/0/solidPointDisplacement_displacement.dat", "postProcessing/fluid/forces/0/force.dat"]:
    pick = lambda d: Path(d)/f if (Path(d)/f).is_file() else Path(d)/Path(f).name
    A, B = rows(pick(a)), rows(pick(b))
    ts = sorted(set(A) & set(B))
    d = np.array([np.abs(A[t]-B[t]) for t in ts]); s = np.array([np.abs(A[t]) for t in ts])
    rel = d / np.maximum(s, 1e-300)
    print(f, "common t", ts[0], "-", ts[-1], "n", len(ts), "identical rows", int((d.max(axis=1)==0).sum()),
          "max abs diff per column", d.max(axis=0), "max rel (col 0)", rel[:,0].max())
