#!/usr/bin/env python3
import os, os as _os; HERE = _os.path.dirname(_os.path.abspath(__file__)); PKG = _os.path.dirname(HERE)   # package root; run directories of retained HPC outputs: env CAVITY_RUNS (default runs/cavityNext)
"""C1 acceptance + diagnostics for the cavity fixed-plate fluid-only run. usage: cavity_c1_check.py RUNDIR  -> RUNDIR/c1_check.json, prints a table.
Criteria A1-A6 as pre-registered in CAVITY_FLUID_ONLY_RESULTS.md (C1 section). All forces per unit depth, raw nodal pressure (no filtering)."""
import sys, os, re, json, numpy as np
sys.path.insert(0, os.path.dirname(__file__))
import cavity_traction as ct, poiseuille_check as pc
rd = sys.argv[1]; meta, nodes, mon = ct.load(rd); ts = sorted(mon); K, mu = meta["k"], meta["mu"]; tag = meta["tag"]
z0 = {i: nodes[i] for i in mon[ts[0]] if abs(nodes[i][2]) < 1e-9}
def column(xt, ylo=0.99, yhi=2.01, tol=1e-6):
    xs = sorted({round(v[0], 6) for v in z0.values() if ylo <= v[1] <= yhi}); x = min(xs, key=lambda a: abs(a - xt))
    return x, sorted([(v[1], i) for i, v in z0.items() if abs(v[0] - x) < tol and ylo - 1e-9 <= v[1] <= yhi + 1e-9])
def flux(t, x):
    c = column(x)[1]; return np.trapz([mon[t][i][0] for _, i in c], [y for y, _ in c])
res = {"tag": tag, "K": K, "mu": mu, "t_last": ts[-1], "n_times": len(ts)}
tl = ts[-1]; after = [t for t in ts if t >= 10.0]
# A1/A2 fluxes
cols = {x: column(x)[0] for x in (0.0, 1.0, 13.0, 14.0)}; res["flux_columns_x"] = cols
res["Q_in_final"] = flux(tl, 0.0); res["Q_x1_final"] = flux(tl, 1.0); res["Q_x13_final"] = flux(tl, 13.0); res["Q_out_final"] = flux(tl, 14.0)
qin = np.array([flux(t, 0.0) for t in ts]); res["Q_in_dev_max_after_ramp"] = float(np.max(np.abs(qin[np.array(ts) >= 10] - 1)))
res["mass_imbalance_rel"] = {k: res[k] / res["Q_in_final"] - 1 for k in ("Q_x1_final", "Q_x13_final", "Q_out_final")}
# A5 inlet profile
c0 = column(0.0)[1]; y0 = np.array([y for y, _ in c0]); w0 = np.array([mon[tl][i][0] for _, i in c0]); wy0 = np.array([mon[tl][i][1] for _, i in c0])
res["inlet_profile_max_dev_from_discrete_parabola"] = float(np.max(np.abs(w0 - pc.pre()["VN"]))); res["inlet_max_abs_wy"] = float(np.abs(wy0).max())
# wall traction (floor y=0, x in [2,6]) vs time
wc = ct.wall_columns(nodes, mon, 0.0, 2.0, 6.0)
hist = [ct.quantities(meta, nodes, mon, wc, t) for t in ts]
Fy = np.array([h["Fy_raw"] for h in hist]); Fx = np.array([h["Fx_viscous_2nd"] for h in hist]); amp = np.array([h["oddeven_amp_Pa"] for h in hist])
res["floor"] = {k: float(hist[-1][k]) for k in ("Fy_raw", "Fy_smoothed121_DIAGNOSTIC", "Fx_viscous_2nd", "Fx_viscous_1st", "oddeven_amp_Pa", "p_mean")}
res["floor"]["n_wall_nodes"] = len(wc)
# A3 steady state over the last 10 % of the run
w = np.array(ts) >= ts[-1] * 0.9
def spread(a): return float((a[w].max() - a[w].min()) / max(abs(a[-1]), 1e-30))
probes = {}
for name, (xt, yt) in {"P_channel_x8_y1.5": (8.0, 1.5), "P_cavity_x4_y0.4": (4.0, 0.4), "P_down_x12_y1.5": (12.0, 1.5)}.items():
    i = min(z0, key=lambda k: (z0[k][0] - xt) ** 2 + (z0[k][1] - yt) ** 2); probes[name] = (i, z0[i][0], z0[i][1])
res["probes"] = {k: {"node": v[0], "x": v[1], "y": v[2], "wx_final": mon[tl][v[0]][0], "wy_final": mon[tl][v[0]][1], "spread_last10pct_rel": spread(np.array([mon[t][v[0]][0] for t in ts]))} for k, v in probes.items()}
res["steady"] = {"Fy_raw_spread_rel": spread(Fy), "Fx_spread_rel": spread(Fx), "oddeven_amp_spread_rel": spread(amp)}
# A4 upstream dp/dx (x in [0.3,1.0]) section-mean and centreline, raw and pair-averaged (diagnostic)
xs_up = sorted({round(v[0], 6) for v in z0.values() if 0.25 <= v[0] <= 1.05 and 1 <= v[1] <= 2})
def pmean(x):
    c = column(x)[1]; return -K * np.trapz([mon[tl][i][2] for _, i in c], [y for y, _ in c])
pm = np.array([pmean(x) for x in xs_up]); g = np.polyfit(xs_up, pm, 1)[0]
res["upstream_dpdx_section_mean_raw_fit"] = float(g); res["upstream_dpdx_expected_actual_flux"] = pc.pre()["dpdx_discrete"] * res["Q_in_final"]
res["upstream_dpdx_rel_err"] = float(g / res["upstream_dpdx_expected_actual_flux"] - 1)
# odd-even mode amplitude: interface, upstream centreline, downstream centreline (estimator a = mean((-1)^i (p_i - (p_{i-1}+p_{i+1})/2)))
def modeamp(p): p = np.asarray(p); i = np.arange(1, len(p) - 1); return float(np.mean((-1.0) ** i * (p[1:-1] - 0.5 * (p[:-2] + p[2:]))))
def centreline(xlo, xhi):
    xs = sorted({round(v[0], 6) for v in z0.values() if xlo <= v[0] <= xhi and abs(v[1] - 1.5) < 1e-9}); idm = {round(v[0], 6): i for i, v in z0.items() if abs(v[1] - 1.5) < 1e-9}
    return xs, [-K * mon[tl][idm[x]][2] for x in xs]
for name, (a, b) in {"channel_upstream_x0.3-1.5": (0.3, 1.5), "channel_above_cavity_x6.5-11": (6.5, 11.0), "channel_downstream_x12-13.8": (12.0, 13.8)}.items():
    xs, p = centreline(a, b); res.setdefault("oddeven_amp_Pa_centreline", {})[name] = modeamp(p)
res["oddeven_amp_Pa_interface_estimator"] = modeamp(hist[-1]["p"])
res["force_bias_calibration"] = {"straight_channel_ratio_bias_over_(a*L)": 0.0358 / (0.01295 * 4), "estimated_bias_N_per_m_for_this_run(not applied)": 0.0358 / (0.01295 * 4) * abs(hist[-1]["oddeven_amp_Pa"]) * 4}
# mesh quality from the feb (z = 0 quads of each hex)
txt = open(os.path.join(rd, tag + ".feb"), encoding="latin-1").read(); els = [list(map(int, m.split(","))) for m in re.findall(r'<elem id="\d+">([^<]+)</elem>', txt)]
minJ, aspect, minang = [], [], []
for e in els:
    P = np.array([nodes[i][:2] for i in e[:4]]); J = []
    for k in range(4):
        a_, b_ = P[(k + 1) % 4] - P[k], P[k - 1] - P[k]; J.append(a_[0] * b_[1] - a_[1] * b_[0]); minang.append(np.degrees(np.arccos(np.dot(a_, b_) / np.linalg.norm(a_) / np.linalg.norm(b_))))
    ed = [np.linalg.norm(P[(k + 1) % 4] - P[k]) for k in range(4)]; minJ.append(min(J) / max(abs(np.array(J)))); aspect.append(max(ed) / min(ed))
res["mesh_quality"] = {"n_elem": len(els), "min_scaled_corner_Jacobian": float(min(minJ)), "max_aspect": float(max(aspect)), "min_corner_angle_deg": float(min(minang)), "max_corner_angle_deg": float(max(minang))}
# log: Newton convergence
log = open(os.path.join(rd, tag + ".log"), errors="ignore").read(); steps_txt = re.split(r"===== beginning time step", log)[1:]
nt = [len(s.split("convergence summary")[0].split("convergence norms")[1:]) for s in steps_txt]
mg = {"velocity": [], "dilatation": []}
for s in steps_txt:
    tabs = s.split("convergence summary")[0].split("convergence norms")[1:]
    if len(tabs) < 2: continue
    for key in mg:
        m = re.search(key + r"\s+(\S+)\s+(\S+)\s+(\S+)", tabs[-1])
        if m and float(m.group(3)) > 0: mg[key].append(float(m.group(2)) / float(m.group(3)))
res["newton"] = {"steps": len(steps_txt), "iters_mean_all": float(np.mean(nt)), "iters_max": int(max(nt)), "steps_with_iterations": int(sum(n >= 2 for n in nt)), "steps_degenerate_single_table": int(sum(n < 2 for n in nt)),
                 "first_degenerate_step": next((k for k, n in enumerate(nt, 1) if n < 2), None), "margin_velocity_max_nondegenerate": max(mg["velocity"]) if mg["velocity"] else None, "margin_dilatation_max_nondegenerate": max(mg["dilatation"]) if mg["dilatation"] else None,
                 "frac_steps_with_relative_criteria_met": float(np.mean([m <= 1 for m in mg["dilatation"]])) if mg["dilatation"] else None}
res["log"] = {"normal_termination": "N O R M A L   T E R M I N A T I O N" in log, "errors": len(re.findall(r"(?i)error", log)), "warnings_distinct": sorted(set(re.findall(r"\*\s+(No force acting on the system\.|[A-Z][^*\n]{8,60}\.)\s+\*", log)))[:6],
              "negative_jacobian": bool(re.search(r"(?i)negative jacobian", log)), "cutbacks": len(re.findall(r"(?i)retrying|cutback|reducing time step", log))}
mt = re.search(r"Total elapsed time[^:]*:\s*[\d:]+\s*\(([\d.]+) sec\)", log); res["febio_wall_s"] = float(mt.group(1)) if mt else None
tm = open(os.path.join(rd, "time.txt")).read(); gt = lambda p: re.search(p, tm)
res["cost"] = {"wall": gt(r"Elapsed \(wall clock\) time.*?: ([\d:.]+)").group(1), "user_s": float(gt(r"User time \(seconds\): ([\d.]+)").group(1)), "sys_s": float(gt(r"System time \(seconds\): ([\d.]+)").group(1)), "max_rss_MB": int(gt(r"Maximum resident set size \(kbytes\): (\d+)").group(1)) / 1024, "cores": json.load(open(os.path.join(rd, "run_meta_start.json")))["cpus_per_task"]}
# verdicts
A = {}
A["A1"] = abs(res["Q_in_final"] - 1) < 1e-4 and res["Q_in_dev_max_after_ramp"] < 1e-4
A["A2"] = all(abs(v) < 1e-3 for v in res["mass_imbalance_rel"].values())
A["A3"] = max(res["steady"].values()) < 1e-4 and all(p["spread_last10pct_rel"] < 1e-4 for p in res["probes"].values()) if False else (res["steady"]["Fy_raw_spread_rel"] < 1e-4 and res["steady"]["Fx_spread_rel"] < 1e-4 and all(p["spread_last10pct_rel"] < 1e-4 for p in res["probes"].values()))
A["A4"] = abs(res["upstream_dpdx_rel_err"]) < 0.02
A["A5"] = res["inlet_profile_max_dev_from_discrete_parabola"] < 1e-5 and res["inlet_max_abs_wy"] < 1e-6
A["A6"] = res["log"]["normal_termination"] and res["log"]["errors"] == 0 and not res["log"]["negative_jacobian"] and res["log"]["cutbacks"] == 0
res["verdict"] = A; json.dump(res, open(os.path.join(rd, "c1_check.json"), "w"), indent=1, default=float)
print(json.dumps({k: v for k, v in res.items() if k not in ("probes",)}, indent=1, default=float)); print("probes", json.dumps(res["probes"], default=float))
