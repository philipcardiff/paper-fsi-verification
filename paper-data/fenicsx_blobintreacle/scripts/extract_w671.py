"""Extract the solids4foam #671 A/B results from the solver repository's raw tables.
Usage: python3 -B scripts/extract_w671.py SOLVER_REPO data/s4f_w671_levels.csv
Reads SOLVER_REPO/s4f_side/w671/WALL671_TABLES.md (generated there by analyze_w671.py from the run logs).
LEGACY = continuityNormalGradient no (old behaviour) and CORRECTED = #671 default, same executable
(proposed solids4foam #671 branch, head da1044a4f, UNMERGED, not an official release).
Writes steady (t_a = 4) and transient (t = 1, dt = 0.005) values per model, arm and level, and the
applied/fluid-reported force ratios of the mesh-quality table."""
import csv, re, sys, os
repo, out = sys.argv[1], sys.argv[2]
text = open(os.path.join(repo, "s4f_side/w671/WALL671_TABLES.md")).read().split("\n")
cells = {"x0.5": 150, "x1": 600, "x2": 2400, "x4": 9600}
rows = []
def cellsplit(l): return [c.strip() for c in l.strip().strip("|").split("|")]
mode = None
for l in text:
    if l.startswith("## Steady (t_a = 4), model"): mode = ("steady", l.split("model")[1].strip()); continue
    if l.startswith("## Transient t = 1"): mode = ("transient", None); continue
    if l.startswith("#"): mode = None if not l.startswith("###") else mode
    if mode is None or not l.startswith("|") or "---" in l or l.startswith("| level") or l.startswith("| model | level"): continue
    c = cellsplit(l)
    if mode[0] == "steady" and len(c) >= 10:
        lev, arm, valid, ux, uy, Fx, Fy, dp = c[:8]
        base = [float(v) for v in c[9].split(",")]
        model, dt = mode[1], (0.005 if lev == "x4" else 0.01)
        vals = {"apex_ux": ux, "apex_uy": uy, "Fx_fluid": Fx, "Fy_fluid": Fy, "dp": dp, "Rx_base": base[0], "Ry_base": base[1]}
        for o, v in vals.items(): rows.append([model, arm, lev, cells[lev], dt, "steady_t4", o, v])
    elif mode[0] == "transient" and len(c) >= 13 and c[3] in ("True", "False"):
        model, lev, arm, valid, ux, uy, Fx, Fy, dp, P, W, Fa, Ff = c[:13]
        for o, v in {"apex_ux": ux, "apex_uy": uy, "Fx_fluid": Fx, "Fy_fluid": Fy, "dp": dp}.items():
            rows.append([model, arm, lev, cells[lev], 0.005, "transient_t1", o, v])
        pa, pf = [float(v) for v in P.split(",")]
        rows.append([model, arm, lev, cells[lev], 0.005, "transient_t1", "P_applied_over_P_fluid", "%.6f" % (pa / pf)])
with open(out, "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["model", "arm", "level", "cells_total", "dt", "state", "observable", "value"])
    w.writerows(rows)
# applied/fluid force ratios (mesh-quality table)
rat = []
for l in text:
    m = re.match(r"\| ([ST])_(orig|earea)_(x[0-9.]+)(_dt[0-9.]+)?_(cor|leg) \|.*\| ([0-9.]+), ([0-9.]+) \| [-0-9.e]+ \|", l)
    if m: rat.append([m.group(1), m.group(2), m.group(3), m.group(4) or "", "CORRECTED" if m.group(5) == "cor" else "LEGACY", m.group(6), m.group(7)])
with open(out.replace("levels", "force_ratio"), "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["state(S steady,T t=1)", "model", "level", "dt", "arm", "F_applied_over_F_fluid_x", "F_applied_over_F_fluid_y"]); w.writerows(rat)
