import os, os as _os; HERE = _os.path.dirname(_os.path.abspath(__file__)); PKG = _os.path.dirname(HERE)   # package root; run directories of retained HPC outputs: env CAVITY_RUNS (default runs/cavityNext)
"""DIAGNOSTIC ONLY: plate-only Q9 StVK response to the pair-mean ((row0+row1)/2) floor pressure of b2_ybranch.py (what the plate would deflect if loaded by the branch mean). usage: run_plate_q9_pairmean.py YBRANCH.json NX NY OUT.json"""
import sys, json, numpy as np, importlib.util
sp_ = importlib.util.spec_from_file_location("pq", HERE + "/plate_q9_nonuniform.py"); pq = importlib.util.module_from_spec(sp_); sp_.loader.exec_module(pq)
d = json.load(open(sys.argv[1])); nx, ny = int(sys.argv[2]), int(sys.argv[3]); x = np.array(d["per_node"]["x_ref"]); p0 = np.array(d["per_node"]["p_row0"]); p1 = np.array(d["per_node"]["p_row1"])
res = {"source": sys.argv[1], "nx": nx, "ny": ny}
for nm, p in (("raw_row0", p0), ("pairmean", 0.5 * (p0 + p1))):
    uy, ux, um = pq.run(nx, ny, "stvk", 1.0, nsteps=10, pfun=lambda xx, p=p: float(np.interp(xx, x, p))); res[nm] = {"uy_A": uy, "ux_A": ux, "uy_mid": um}; print(nm, uy, ux, flush=True)
json.dump(res, open(sys.argv[4], "w"), indent=1)
