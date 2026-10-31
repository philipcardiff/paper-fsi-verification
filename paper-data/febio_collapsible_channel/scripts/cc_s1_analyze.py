#!/usr/bin/env python3
"""CC-S1 analysis: static clamped wall midpoint/quarter displacement vs the oomph-lib beam values (-0.14362832, -0.09645925).  usage: cc_s1_analyze.py RUN_DIR [RUN_DIR...]"""
import sys, json, re, numpy as np
REF = {"mid": -0.14362832, "q": -0.09645925}
for rd in sys.argv[1:]:
    rows = []; T = None
    for l in open(f"{rd}/Mon_bottom.txt"):
        if l.startswith("*Time"): T = float(l.split("=")[1]); rows = []
        elif not l.startswith("*") and l.strip(): rows.append(l.split())
    A = np.array(rows, float); x = A[:, 1]; ux = A[:, 3]; uy = A[:, 4]
    f = lambda xx: float(uy[int(np.argmin(abs(x - xx)))])
    mid, q1, q3 = f(10.0), f(7.5), f(12.5)
    log = open(glob := [p for p in __import__("glob").glob(f"{rd}/*.log") if "run" not in p.split("/")[-1]][0]).read()
    it = re.search(r"Total number of equilibrium iterations \.+ : (\d+)", log); avg = re.search(r"Average number of equilibrium iterations \.+ : ([\d.]+)", log)
    nt = "N O R M A L   T E R M I N A T I O N" in log
    res = dict(run=rd, T=T, mid=mid, e_mid=mid / REF["mid"] - 1, q1=q1, e_q1=q1 / REF["q"] - 1, q3=q3, e_q3=q3 / REF["q"] - 1, max_ux=float(np.abs(ux).max()), newton_total=int(it.group(1)) if it else None, newton_avg=float(avg.group(1)) if avg else None, normal_termination=nt)
    json.dump(res, open(f"{rd}/s1_result.json", "w"), indent=1)
    print(f"{rd.split('/')[-1]:12s} T={T} mid {mid:.6f} ({100*res['e_mid']:+.3f} %)  q1 {q1:.6f} ({100*res['e_q1']:+.3f} %)  q3 {q3:.6f} ({100*res['e_q3']:+.3f} %)  max|ux| {res['max_ux']:.2e}  Newton {res['newton_total']} (avg {res['newton_avg']}) normal {nt}")
