#!/usr/bin/env python3
"""Stage B analysis: reads runs/hronturekBC/B_*/ (pointA.txt, react_root.txt, <case>.log, model_meta.json) and applies the pre-registered rules (HRON_TUREK_FSI1_STAGE_BC_RESULTS.md section 1.4)."""
import os, re, json, glob, math, sys
R = os.path.abspath(sys.argv[1]) if len(sys.argv) > 1 else os.getcwd()   # directory containing the B_* run directories
REF = {"CSM1": (-7.18777e-3, -66.1029e-3), "CSM2": (-0.469006e-3, -16.9740e-3)}
def last_block(path):
    cur = None
    for ln in open(path):
        if ln.startswith("*Time"): cur = []
        elif not ln.startswith("*") and ln.strip() and cur is not None: cur.append([float(x) for x in ln.split()])
    return cur
def newton(log):
    its = [int(m) for m in re.findall(r"Nonlinear solution status: time= [0-9.e+-]+\s*\n\s*stiffness updates\s*=\s*\d+\s*\n\s*right hand side evaluations\s*=\s*(\d+)", log)]
    return its
def one(d):
    c = os.path.basename(d); m = json.load(open(d + "/model_meta.json"))
    a = last_block(d + "/pointA.txt")[0]; x, y, ux, uy, uz = a[1], a[2], a[3], a[4], a[5]
    rb = last_block(d + "/react_root.txt"); Rx = sum(r[1] for r in rb); Ry = sum(r[2] for r in rb); Rz = sum(r[3] for r in rb)
    log = open(f"{d}/{c}.log", errors="replace").read()
    steps = len(re.findall(r"converged at time", log)); cut = len(re.findall(r"[Rr]etrying|cut|Reducing time step", log))
    rhs = [int(v) for v in re.findall(r"right hand side evaluations\s*=\s*(\d+)", log)]
    W = m["weight_total"]
    return dict(case=c, nt=m["nt"], elem=m["elem"], E=m["E"], g=m["g"], root=m["x_root"], nel=m["n_elem"], x=x, y=y, ux=ux, uy=uy, uz=uz, Ry_over_W=Ry / W, Rx_over_W=Rx / W,
                steps=steps, max_rhs=max(rhs) if rhs else None, normal="N O R M A L   T E R M I N A T I O N" in log)
rows = [one(d) for d in sorted(glob.glob(R + "/B_*")) if os.path.exists(d + "/pointA.txt")]
json.dump(rows, open(R + "/stageB_results.json", "w"), indent=1)
print(f"{'case':28s} {'nel':>5s} {'ux[mm]':>10s} {'uy[mm]':>10s} {'ux err%':>8s} {'uy err%':>8s} {'Ry/W':>12s} {'Rx/W':>10s} {'uz':>8s} steps")
for r in rows:
    k = "CSM1" if "CSM1" in r["case"] else "CSM2" if "CSM2" in r["case"] else None
    if k:
        ex, ey = (r["ux"] - REF[k][0]) / REF[k][0] * 100, (r["uy"] - REF[k][1]) / REF[k][1] * 100
    else: ex = ey = float("nan")
    print(f"{r['case']:28s} {r['nel']:5d} {r['ux']*1e3:10.5f} {r['uy']*1e3:10.5f} {ex:8.3f} {ey:8.3f} {r['Ry_over_W']:12.8f} {r['Rx_over_W']:10.2e} {r['uz']:8.1e} {r['steps']} rhs<= {r['max_rhs']} normal={r['normal']}")
