#!/usr/bin/env python3
"""Builds results/stageB_summary.csv, stageB_convergence.csv from results/stageB_raw/*  (pointA.txt, <case>.log, model_meta.json) using the pre-registered reference values.
usage: make_stageB_tables.py RAW_DIR OUT_DIR"""
import sys, os, re, json, glob, csv, math
raw, out = os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2])
REF = {"CSM1": (-7.18777e-3, -66.1029e-3), "CSM2": (-0.469006e-3, -16.9740e-3)}
def last_block(path):
    cur = None
    for ln in open(path):
        if ln.startswith("*Time"): cur = []
        elif not ln.startswith("*") and ln.strip() and cur is not None: cur.append([float(x) for x in ln.split()])
    return cur
rows = {}
for d in sorted(glob.glob(raw + "/B_*")):
    c = os.path.basename(d); m = json.load(open(d + "/model_meta.json")); a = last_block(d + "/pointA.txt")[0]; log = open(f"{d}/{c}.log", errors="replace").read()
    steps = len(re.findall(r"converged at time", log)); its = int(re.search(r"Total number of equilibrium iterations\s*\.+\s*:\s*(\d+)", log).group(1))
    rhs = int(re.search(r"Total number of right hand evaluations\s*\.+\s*:\s*(\d+)", log).group(1))
    k = "CSM1" if "CSM1" in c else "CSM2" if "CSM2" in c else "LIN"
    r = dict(case=c, benchmark=k, element=m["elem"], nt=m["nt"], nx=m["nx"], elements=m["n_elem"], nodes=m["n_nodes"], root_x=round(m["x_root"], 9), length=round(m["length"], 6), E=m["E"], g=m["g"],
             ux_m=a[3], uy_m=a[4], uz_m=a[5], point_A_deformed_y=a[2], load_steps_converged=steps, equilibrium_iterations=its, iterations_per_step=its / max(steps, 1), rhs_evaluations=rhs,
             normal_termination=("N O R M A L   T E R M I N A T I O N" in log))
    if k in REF:
        r["ux_ref_m"], r["uy_ref_m"] = REF[k]; r["ux_err_pct"] = 100 * (a[3] - REF[k][0]) / REF[k][0]; r["uy_err_pct"] = 100 * (a[4] - REF[k][1]) / REF[k][1]
    rows[c] = r
os.makedirs(out, exist_ok=True)
keys = []
for r in rows.values():
    for kk in r:
        if kk not in keys: keys.append(kk)
with open(out + "/stageB_summary.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, keys); w.writeheader(); w.writerows(rows.values())
conv = []
for k in ("CSM1", "CSM2"):
    ys = [rows[f"B_{k}_h8_nt{n}"] for n in (2, 4, 8, 16)]
    for q in ("uy_m", "ux_m"):
        v = [abs(r[q]) for r in ys]; mono = all(v[i] < v[i + 1] for i in range(3))
        d1, d2 = v[2] - v[1], v[3] - v[2]; p = math.log(d1 / d2) / math.log(2); ex = v[3] + d2 / (2 ** p - 1)
        ref = abs(REF[k][1 if q == "uy_m" else 0])
        conv.append(dict(benchmark=k, quantity=q, monotone_increasing_nt_2_4_8_16=mono, observed_order_nt_4_8_16=p, richardson_abs=ex, richardson_err_pct=100 * (ex - ref) / ref,
                         finest_hex8_err_pct=100 * (v[3] - ref) / ref, hex20_nt4_err_pct=100 * (abs(rows[f"B_{k}_h20_nt4"][q]) - ref) / ref, hex20_nt2_err_pct=100 * (abs(rows[f"B_{k}_h20_nt2"][q]) - ref) / ref,
                         root025_nt16_err_pct=100 * (abs(rows[f"B_{k}_h8_nt16_root025"][q]) - ref) / ref))
with open(out + "/stageB_convergence.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, list(conv[0])); w.writeheader(); w.writerows(conv)
for c in conv: print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in c.items()})
