"""Matched-physical-time comparison FEniCSx (ours, Gmi and Gc dev) vs OUR solids4foam feF1-protocol reconstruction (not Philip's runs). usage: compare_matched.py s4f.csv out_prefix
Only identical physical times are compared. Forces: s4f 'solid interface' = force received by the solid; FE F_dir = facet traction integral on the solid (history.csv Fy_fluid/Fz_fluid); FE nodal forces exist only at the force-decomposition states."""
import csv, sys, json
import numpy as np
R = __import__("os").environ.get("HESS_ROOT", ".") + "/"
def rd(p): return list(csv.DictReader(open(R + p)))
def series(paths):
    d = {}
    for p in paths:
        for x in rd(p): d[round(float(x["t"]), 6)] = x
    return d
gmi = series(["runs/D/Gmi_coupled_n32b/history.csv", "runs/P/Gmi_dt01/history.csv"])
gc = series(["runs/S/T3b/history.csv"])
s4f = {round(float(x["t"]), 6): x for x in csv.DictReader(open(sys.argv[1]))}
times = [t for t in (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75) if t in s4f]
out = []
for t in times:
    s = s4f[t]; row = dict(t=t, s4f_tip_mm=float(s["tip_dy_m"]) * 1e3, s4f_Fy_on_solid_uN=float(s["Fs_y"]) * 1e6, s4f_Fz_on_solid_uN=float(s["Fs_z"]) * 1e6,
                           s4f_Fy_fluidside_signflipped_uN=-float(s["Ff_y"]) * 1e6, s4f_Fz_fluidside_signflipped_uN=-float(s["Ff_z"]) * 1e6, s4f_coupling_its=s["its"])
    for nm, d in (("Gmi", gmi), ("Gc", gc)):
        x = d.get(t)
        if x is None: continue
        row[nm + "_tip_mm"] = float(x["tip_y"]) * 1e3; row[nm + "_Fy_Fdir_uN"] = float(x["Fy_fluid"]) * 1e6; row[nm + "_Fz_Fdir_uN"] = float(x["Fz_fluid"]) * 1e6
        row[nm + "_tip_rel_diff_vs_s4f"] = row[nm + "_tip_mm"] / row["s4f_tip_mm"] - 1
        row[nm + "_Fy_rel_diff_vs_s4f"] = row[nm + "_Fy_Fdir_uN"] / row["s4f_Fy_on_solid_uN"] - 1
        row[nm + "_Fz_rel_diff_vs_s4f"] = row[nm + "_Fz_Fdir_uN"] / row["s4f_Fz_on_solid_uN"] - 1
    out.append(row)
keys = sorted({k for r in out for k in r}, key=lambda k: (k != "t", k))
with open(sys.argv[2] + ".csv", "w") as f:
    w = csv.DictWriter(f, fieldnames=keys); w.writeheader(); [w.writerow(r) for r in out]
for r in out: print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()})
