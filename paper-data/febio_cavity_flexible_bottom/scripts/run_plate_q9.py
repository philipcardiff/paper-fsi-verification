import os, os as _os; HERE = _os.path.dirname(_os.path.abspath(__file__)); PKG = _os.path.dirname(HERE)   # package root; run directories of retained HPC outputs: env CAVITY_RUNS (default runs/cavityNext)
"""Plate-only Q9 StVK response to (a) FEBio interface p(x) from the xplt, (b) uniform pressure with the same F_y, (c) uniform 5.2504/4 (s4f QM level-3 F_y). usage: run_plate_q9.py AUDIT.json STATE_KEY NX NY OUT.json"""
import sys, json, numpy as np, importlib.util
sp_ = importlib.util.spec_from_file_location("pq", HERE + "/plate_q9_nonuniform.py"); pq = importlib.util.module_from_spec(sp_); sp_.loader.exec_module(pq)
a, key, nx, ny, out = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
d = json.load(open(a))[key]; xr = np.array(d["interface_p_along_x"]["x_ref"]); p = np.array(d["interface_p_along_x"]["p"])
pf = lambda x: float(np.interp(x, xr, p))
Fy = -d["pressure_only"]["Fy"]; pm = Fy / 4.0
res = {"source": a, "state": key, "time": d["time"], "Fy_febio_pressure": -Fy, "nx": nx, "ny": ny, "p_mean_force_equiv": pm}
for nm, q, f in (("febio_p_of_x", 1.0, pf), ("uniform_same_Fy", pm, None), ("uniform_s4f_QM_L3_Fy_5.2504", 5.2504 / 4, None)):
    uy, ux, um = pq.run(nx, ny, "stvk", q, nsteps=10, pfun=f); res[nm] = {"uy_A": uy, "ux_A": ux, "uy_mid": um}; print(nm, uy, ux, um, flush=True)
json.dump(res, open(out, "w"), indent=1)
