"""Extract the inputs of the local-comparison and corner sections from the preflight (C1) raw results; also derives the local floor-pressure and plate-shape tables.
usage: python3 -B extract_local_and_corner.py C1_DIR DUMP_DAT OUT_DATA_DIR
  C1_DIR  : .../fenicsx-cavity-C1 (runs/c2_dump_*, runs/c1_ladder_*, runs/c1_corner*_)
  DUMP_DAT: paper-data/cavity_strict_M2_floor/cavity_strict_M2_floor_t200.dat (PR #67)
Writes plate_only_nt32_control.csv, graded_mesh_attempts.csv, floor_pressure_local_comparison.csv, plate_shape_bins.csv.  Requires numpy."""
import sys, os, glob, json, csv
import numpy as np
c1, dump, out = sys.argv[1:4]
W = lambda n, h, r: csv.writer(open(os.path.join(out, n), "w", newline=""), lineterminator="\n").writerows([h] + r)
J = lambda f: json.load(open(f))
# nt = 4, 8, 16 from the preflight ladder, nt = 32 from the corner-experiment control
lad = {}
for f in glob.glob(c1 + "/runs/c1_ladder_*/results/*.json"):
    d = J(f); lad[(d["load"], d["degree"], d["nt"])] = d
rows = []
for ld in ("uniform", "follower"):
    for nt in (4, 8, 16):
        d = lad[(ld, 2, nt)]; rows.append([ld, 2, nt, d["ndofs"], d["uy_A"], d["converged"], "c1_ladder (preflight)"])
    for run in ("c1_corner_n2_*", "c1_corner2_n2_*"):
        f = glob.glob(c1 + "/runs/%s/results/c1plate_unif32_%s_q2_nt32__*.json" % (run, ld))[0]; d = J(f)
        rows.append([ld, 2, 32, d["ndofs"], d["uy_A"], d["converged"], "c1_corner control (%s)" % ("job 16716" if run == "c1_corner_n2_*" else "job 16722, amendment 1")])
W("plate_only_nt32_control.csv", ["load", "degree", "nt", "ndofs", "uy_A_m", "converged", "source"], rows)
# graded-mesh attempts: values are written ONLY for converged cases (diverged states are discarded)
rows = []
for run, job in (("c1_corner_n2_*", 16716), ("c1_corner2_n2_*", 16722)):
    for f in sorted(glob.glob(c1 + "/runs/%s/results/c1plate_graded_*.json" % run)):
        d = J(f); ok = bool(d["converged"])
        rows.append([job, d["load"], d["degree"], d["nt"], d["ndofs"], d["steps"], ok, d["snes_reason"], d["uy_A"] if ok else ""])
W("graded_mesh_attempts.csv", ["slurm_job", "load", "degree", "nt", "ndofs", "load_steps", "converged", "snes_reason", "uy_A_m_(converged_cases_only)"], rows)
# local floor pressure: FEniCSx rigid floor on the s4f level-2 shape against the PR #67 dump
a = np.loadtxt(dump); x, y, Sx, Sy, p = a[:, 0], a[:, 1], a[:, 2], a[:, 3], a[:, 5]
def fe_p(f):
    d = J(f); fx = np.array(d["floor_x"]); fp = np.array(d["floor_p"]); o = np.argsort(fx); return np.interp(x, fx[o], fp[o])
files = sorted(glob.glob(c1 + "/runs/c2_dump_*/results/*.json"), key=lambda f: -float(J(f)["mesh"].split("_h")[-1][:-4]))
pf = fe_p(files[-1]); Fs, Ff = (p * Sy).sum(), (pf * Sy).sum(); dF = (pf - p) * Sy; r = pf / p
cs, cf = (p * Sy * x).sum() / Fs, (pf * Sy * x).sum() / Ff
rows = [["sum p*Sf_y, s4f level 2 (PR 67 dump)", Fs], ["sum p*Sf_y, FEniCSx rigid floor h=0.0125 on the same faces", Ff], ["relative difference (percent)", 100 * (Ff / Fs - 1)],
        ["mean of ratio p_FE/p_s4f", r.mean()], ["std of ratio (percent)", 100 * r.std()], ["sum|dF| / |sum dF| (1 = all faces same sign)", np.abs(dF).sum() / abs(dF.sum())],
        ["centre of pressure s4f (m)", cs], ["centre of pressure FEniCSx (m)", cf], ["difference of centres (m)", cs - cf],
        ["s4f level-2 F_y against its own strict limit -5.260691 (percent)", 100 * (-5.23881 / -5.260691 - 1)],
        ["integrated viscous F_x, FEniCSx -0.02717 against s4f -0.02680735 (percent; s4f does not export shear: not like-for-like)", 100 * (-0.02717 / -0.02680735 - 1)]]
for f in files[:-1]:
    d = J(f); pp = fe_p(f); rows.append(["RMS of local p difference to the finest FEniCSx floor, " + os.path.basename(d["mesh"]) + " (fraction of mean s4f p)", float(np.sqrt(np.mean((pp - pf) ** 2)) / p.mean())])
W("floor_pressure_local_comparison.csv", ["quantity", "value"], rows)
# plate shape: FEniCSx plate-only follower Q3 nt16 under the recorded s4f level-2 load against the s4f level-2 deformed floor y(x)
d = lad[("follower", 3, 16)]; fe = np.interp(x, d["uy_top_profile_x"], d["uy_top_profile"]); rows = []
for lo, hi in [(2, 2.4), (2.4, 3), (3, 3.5), (3.5, 4.5), (4.5, 5.1), (5.1, 5.6), (5.6, 6.0)]:
    m = (x >= lo) & (x < hi); rows.append([lo, hi, fe[m].mean(), y[m].mean(), y[m].mean() - fe[m].mean(), 100 * (y[m].mean() - fe[m].mean()) / abs(d["uy_A"]), y[m].mean() / fe[m].mean()])
rows.append(["all faces", "", fe.mean(), y.mean(), (y - fe).mean(), 100 * (y - fe).mean() / abs(d["uy_A"]), ""])
W("plate_shape_bins.csv", ["x_lo_m", "x_hi_m", "FEniCSx_plate_uy_m", "s4f_level2_y_m", "s4f_minus_FE_m", "difference_percent_of_|uy_A|", "ratio_s4f_over_FE"], rows)
print("ok")
