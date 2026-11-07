#!/usr/bin/env python3
import os, os as _os; HERE = _os.path.dirname(_os.path.abspath(__file__)); PKG = _os.path.dirname(HERE)   # package root; run directories of retained HPC outputs: env CAVITY_RUNS (default runs/cavityNext)
"""B2: compare interface distributions p(x) (and, if present, tangential traction) of FEBio against a reference table WITHOUT smoothing or shifting.
usage: b2_profile_compare.py FEBIO_AUDIT.json STATE_KEY REF.csv OUT.json [--coord xref|arc]
REF.csv columns (header required): x_ref, p [, x_def, y_def, tx, ty]   (p in Pa = kinematic p * rho, rho = 1; x_ref = undeformed x in m, increasing; tx,ty = viscous traction of fluid on plate, Pa)
Common coordinate: 'xref' (default; undeformed floor position, the same material point in both solids) or 'arc' (arc length of the deformed polyline, normalised to the length of each curve: use only if x_def,y_def are given in REF).
Comparison: the reference is linearly interpolated onto the FEBio nodes (raw nodal p of FEBio, odd-even mode untouched, reported separately as even/odd means); differences are reported node by node and as integrals
 (pressure force F_y, F_x of each curve with its OWN deformed geometry, plus FEBio p integrated on the reference geometry as a consequence-vs-cause split)."""
import sys, json, csv, numpy as np
a, key, ref, outp = sys.argv[1:5]; coord = sys.argv[sys.argv.index("--coord") + 1] if "--coord" in sys.argv else "xref"
d = json.load(open(a))[key]; I = d["interface_p_along_x"]; xr, xd, yd, p = (np.array(I[k]) for k in ("x_ref", "x_def", "y_def", "p"))
R = list(csv.DictReader(open(ref))); col = lambda k: np.array([float(r[k]) for r in R]) if k in R[0] else None
rx, rp, rxd, ryd = col("x_ref"), col("p"), col("x_def"), col("y_def"); o = np.argsort(rx); rx, rp = rx[o], rp[o]
if rxd is not None: rxd, ryd = rxd[o], ryd[o]
def arcs(X, Y): s = np.r_[0, np.cumsum(np.hypot(np.diff(X), np.diff(Y)))]; return s / s[-1]
if coord == "arc":
    assert rxd is not None; cf, cr = arcs(xd, yd), arcs(rxd, ryd)
else: cf, cr = xr, rx
pr = np.interp(cf, cr, rp)                                   # reference p at the FEBio nodes (linear, no smoothing)
diff = p - pr; i = np.arange(len(p)); mask = (i > 0) & (i < len(p) - 1)   # nodes 0 and N are the corner pins in FEBio (1.15 Pa constraint): reported but flagged
def F(X, Y, P): dx, dy = np.diff(X), np.diff(Y); pb = 0.5 * (P[1:] + P[:-1]); return float(-np.sum(pb * dx)), float(np.sum(pb * dy))
res = {"febio": a, "state": key, "time": d["time"], "coordinate": coord, "n_nodes": int(len(p)), "nodes_flagged_pinned": [0, int(len(p) - 1)],
       "p_mean_febio": float(p.mean()), "p_mean_ref_on_febio_nodes": float(pr.mean()), "diff_mean_all": float(diff.mean()), "diff_mean_unpinned": float(diff[mask].mean()), "diff_rms_unpinned": float(np.sqrt((diff[mask] ** 2).mean())),
       "diff_even_nodes_mean": float(diff[(i % 2 == 0) & mask].mean()), "diff_odd_nodes_mean": float(diff[(i % 2 == 1) & mask].mean()),
       "febio_oddeven_amp_even_minus_odd": float(p[i % 2 == 0].mean() - p[i % 2 == 1].mean())}
thirds = lambda v: [float(v[(cf >= lo) & (cf <= hi)].mean()) for lo, hi in ((cf.min(), cf.min() + (cf.max() - cf.min()) / 3), (cf.min() + (cf.max() - cf.min()) / 3, cf.min() + 2 * (cf.max() - cf.min()) / 3), (cf.min() + 2 * (cf.max() - cf.min()) / 3, cf.max()))]
res["thirds_febio"] = thirds(p); res["thirds_ref"] = thirds(pr)
res["F_pressure_febio_own_geometry"] = F(xd, yd, p)
if rxd is not None:
    res["F_pressure_ref_own_geometry_on_its_nodes"] = F(rxd, ryd, rp)
    res["F_pressure_febio_p_on_ref_geometry"] = F(np.interp(cf, cr, rxd) if coord == "arc" else np.interp(xr, rx, rxd), np.interp(cf, cr, ryd) if coord == "arc" else np.interp(xr, rx, ryd), p)
    res["F_pressure_ref_p_on_febio_geometry"] = F(xd, yd, pr)
res["F_y_from_projection_unit_check"] = float(-np.trapz(p, xd))
json.dump(dict(res, node_table={"x_ref": xr.tolist(), "p_febio": p.tolist(), "p_ref_interp": pr.tolist(), "diff": diff.tolist()}), open(outp, "w"), indent=1)
print(json.dumps(res, indent=1))
