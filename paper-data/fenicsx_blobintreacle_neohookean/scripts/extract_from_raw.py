#!/usr/bin/env python3
"""Regenerate the FEniCSx input CSVs of this directory from the raw run records (evidence tarball).

usage: python3 -B -I scripts/extract_from_raw.py EVIDENCE_DIR [OUT_DIR]
EVIDENCE_DIR contains runs/nhblob_<id>_.../RUNDIR/{history.csv,summary.json,interface_final.csv} and runs/nhblob_<id>_.../EXIT.
Writes fe_steady_levels.csv, fe_time_study_t1.csv, fe_run_validity.csv, fe_levels_t1_t3.csv into OUT_DIR (default ../data).
Only the standard library is used."""
import csv, glob, json, math, os, re, sys

ev = sys.argv[1]
out = sys.argv[2] if len(sys.argv) > 2 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
runs = {}
for d in sorted(glob.glob(os.path.join(ev, "runs", "nhblob_*"))):
    m = re.match(r"nhblob_(.+?)_f[0-9.]+_dt", os.path.basename(d))
    if m: runs[m.group(1)] = d


def hist(rid):
    return list(csv.DictReader(open(os.path.join(runs[rid], "RUNDIR", "history.csv"))))


def at(h, t):
    r = [x for x in h if abs(float(x["t"]) - t) < 1e-6]
    return r[0]


def f(x):
    return float(x) if x not in ("", None) else float("nan")


COLS = (("ux", "apex_ux"), ("uy", "apex_uy"), ("Fx", "F_M1x"), ("Fy", "F_M1y"), ("dp", "dp"), ("Fx_weak", "Fwex_x"), ("Fy_weak", "Fwex_y"))
CELLS = {"0.5": 150, "1": 600, "2": 2400, "4": 9600}
# steady levels at t = 6 (V2 drift over the last 0.5 time units), and t = 1, t = 3 values
with open(os.path.join(out, "fe_steady_levels.csv"), "w", newline="") as fh, open(os.path.join(out, "fe_levels_t1_t3.csv"), "w", newline="") as fh2:
    w = csv.writer(fh); w.writerow(["level", "cells", "ux", "uy", "Fx", "Fy", "dp", "Fx_weak", "Fy_weak", "drift_ux", "drift_uy", "drift_Fx", "drift_Fy", "drift_dp"])
    w2 = csv.writer(fh2); w2.writerow(["level", "cells", "t", "ux", "uy", "Fx", "Fy", "dp", "Fx_weak", "Fy_weak"])
    for lv in CELLS:
        h = hist("S" + lv); last = h[-1]; win = [r for r in h if float(r["t"]) >= float(last["t"]) - 0.5 - 1e-9]
        vals = [f(last[c]) for _, c in COLS]
        drift = []                                   # V2: (max - min) / |last| over the last 0.5 time units
        for c in ("apex_ux", "apex_uy", "F_M1x", "F_M1y", "dp"):
            v = [f(r[c]) for r in win]; drift.append((max(v) - min(v)) / abs(f(last[c])))
        w.writerow([lv, CELLS[lv]] + ["%.17g" % v for v in vals] + ["%.17g" % v for v in drift])
        for t in (1.0, 3.0):
            x = at(h, t); w2.writerow([lv, CELLS[lv], "%g" % t] + ["%.17g" % f(x[c]) for _, c in COLS])
# time study at t = 1 (x1): values and interface-vertex change between successive dt
dts = [0.25, 0.125, 0.0625, 0.03125, 0.015625, 0.0078125]
vals, ifc = {}, {}
for dt in dts:
    rid = "T%g" % dt
    x = at(hist(rid), 1.0); vals[dt] = [f(x[c]) for _, c in COLS[:5]]
    ifc[dt] = [(float(r["dx"]), float(r["dy"])) for r in csv.DictReader(open(os.path.join(runs[rid], "RUNDIR", "interface_final.csv")))]
with open(os.path.join(out, "fe_time_study_t1.csv"), "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["dt", "ux", "uy", "Fx", "Fy", "dp", "dux", "iface_l2_all", "iface_l2_vertices"])
    for i, dt in enumerate(dts):
        row = ["%.17g" % dt] + ["%.17g" % v for v in vals[dt]]
        if i:
            p = dts[i - 1]
            d_all = math.sqrt(sum((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 for a, b in zip(ifc[dt], ifc[p])))
            d_vert = math.sqrt(sum((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 for a, b in list(zip(ifc[dt], ifc[p]))[::2]))
            row += ["%.17g" % abs(vals[dt][0] - vals[p][0]), "%.17g" % d_all, "%.17g" % d_vert]
        else:
            row += ["", "", ""]
        w.writerow(row)
# time study repeated with a tighter SNES tolerance (A1T) and with a 16-substep first step (A2T), values at t = 1
with open(os.path.join(out, "fe_addenda_t1_values.csv"), "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["series", "dt", "ux", "uy", "Fx", "Fy", "dp"])
    for ser, pre in (("main", "T"), ("A1_tight_SNES", "A1T"), ("A2_16_substeps", "A2T")):
        for dt in dts:
            x = at(hist("%s%g" % (pre, dt)), 1.0); w.writerow([ser, "%g" % dt] + ["%.17g" % f(x[c]) for _, c in COLS[:5]])
# validity, closure and mesh-quality summary for every run (V1, V6, V7)
with open(os.path.join(out, "fe_run_validity.csv"), "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["run", "steps", "reached_end", "failed", "exit_ok", "snes_reason3_steps", "snes_min_reason", "closure_iface_max", "closure_global_max", "P_M1_last", "Jmin_solid", "Jmin_fluid", "sj_fluid_min", "sj_solid_min", "area_ratio_fluid_min", "area_ratio_solid_min", "wall_s", "cores", "cpu_h"])
    for rid in runs:
        d = runs[rid]; s = json.load(open(os.path.join(d, "RUNDIR", "summary.json"))); h = hist(rid)
        ex = open(os.path.join(d, "EXIT")).read().split()[0] if os.path.exists(os.path.join(d, "EXIT")) else "missing"
        skip = 1 if rid.startswith("A2T") else 0          # step-1 diagnostic uses the last-substep history (not meaningful)
        ci = max(max(abs(f(r["Fwex_x"]) - f(r["Rsex_x"])), abs(f(r["Fwex_y"]) - f(r["Rsex_y"]))) for r in h[skip:])
        cg = max(max(abs(f(r["Rsex_x"]) + f(r["root_Rx"]) - f(r["inert_x"])), abs(f(r["Rsex_y"]) + f(r["root_Ry"]) - f(r["inert_y"]))) for r in h[skip:])
        reasons = [int(r["snes_reason"]) for r in h]
        mn = lambda c: min((f(r[c]) for r in h if r.get(c) not in ("", None)), default=float("nan"))
        w.writerow([rid, len(h), s["reached_end"], s["failed"], ex, sum(1 for x in reasons if x == 3), min(reasons), "%.3e" % ci, "%.3e" % cg, "%.3e" % abs(f(h[-1]["P_M1"])),
                    "%.6g" % mn("Jmin_solid"), "%.6g" % mn("Jmin_fluid"), "%.6g" % mn("sj_fluid_min"), "%.6g" % mn("sj_solid_min"), "%.6g" % mn("area_ratio_fluid_min"), "%.6g" % mn("area_ratio_solid_min"),
                    "%.1f" % s["wall_s"], s["cores"], "%.4f" % s["cpu_h"]])
# sampled mesh-quality history of the x4 repeat (the extended logging exists in the A1 runs only; sampled every 10 steps)
with open(os.path.join(out, "fe_ale_x4_samples.csv"), "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["t", "Jmin_fluid", "Jmin_solid", "sj_fluid_min", "sj_solid_min"])
    for r in hist("A1S4"):
        if r.get("Jmin_fluid"): w.writerow(["%.4g" % float(r["t"]), r["Jmin_fluid"], r["Jmin_solid"], r["sj_fluid_min"], r["sj_solid_min"]])
print("wrote", out)
