"""Interface transfer ratios from raw solver output (not shipped; paths are arguments).
Usage: python3 -B scripts/transfer_ratios_from_raw.py SOLVER_REPO OUT.csv
SOLVER_REPO = checkout of the FEniCSx campaign repository (verification/audit/steady_obs2_f*.json,
s4f_side/data/dist_orig_x*.csv, s4f_side/series_results.json).
Definitions (x components):
  F_fluid  = fluid-reported force on the CURRENT interface
             FEniCSx: Fx_dir = F_M1exact_x (exact interface integral of sigma_f J F^-T N)
             solids4foam: forces function object (series_results.json)
  F_applied = force applied to the SOLID on reference areas
             FEniCSx: F_M2applied_x; solids4foam: sum of written traction x |S_ref| (dist_orig_x*.csv)
  R_base   = support (base) reaction: FEniCSx Rbase_x (nodal residual); solids4foam sum of boundary stress (estimator)
  deficit  = 1 - F_applied / F_fluid        (basis: fluid-reported force)
  excess   = F_fluid / F_applied - 1        (basis: applied force; this is the '13-14 %' normalisation)
  base_deficit = 1 - |R_base| / F_fluid     (basis: fluid-reported force, stress-estimator reaction for solids4foam)"""
import csv, json, os, sys
repo, out = sys.argv[1], sys.argv[2]
rows = []
for f in ["0.5", "1", "2", "4"]:
    r = json.load(open(os.path.join(repo, "verification/audit/steady_obs2_f%s.json" % f)))["records"][-1]
    for c, (ff, fa, rb) in {"x": (r["Fx_dir"], r["F_M2applied_x"], r["Rbase_x"]), "y": (r["Fy_dir"], r["F_M2applied_y"], r["Rbase_y"])}.items():
        rows.append(["FEniCSx M2", "F" + f, c, ff, fa, rb])
S = json.load(open(os.path.join(repo, "s4f_side/series_results.json")))
def dist_sum(x):
    fx = fy = 0.0
    for line in open(os.path.join(repo, "s4f_side/data/dist_orig_x%s.csv" % x)):
        if line.startswith("#") or not line.strip(): continue
        p = line.split(","); fx += float(p[2]); fy += float(p[3])
    return fx, fy
for x in ["0.5", "1", "2", "4"]:
    d = S["x" + x]; app = dist_sum(x)
    for k, c in enumerate("xy"):
        rows.append(["solids4foam M2 (legacy BC)", "x" + x, c, d["fluid_force_tA"]["total"][k], app[k], d["base_reaction"]["t4"][k]])
with open(out, "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["code", "level", "component", "F_fluid_reported_current", "F_applied_to_solid_reference", "R_base", "deficit_1_minus_applied_over_fluid", "excess_fluid_over_applied_minus_1", "base_deficit_1_minus_abs_Rbase_over_fluid", "base_over_applied"])
    for code, lev, c, ff, fa, rb in rows:
        w.writerow([code, lev, c, "%.7g" % ff, "%.7g" % fa, "%.7g" % rb, "%.4f" % (1 - fa / ff), "%.4f" % (ff / fa - 1), "%.4f" % (1 - abs(rb) / abs(ff)), "%.4f" % (abs(rb) / abs(fa))])
