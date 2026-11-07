#!/usr/bin/env python3
import os, os as _os; HERE = _os.path.dirname(_os.path.abspath(__file__)); PKG = _os.path.dirname(HERE)   # package root; run directories of retained HPC outputs: env CAVITY_RUNS (default runs/cavityNext)
"""Straight-channel Poiseuille check for the FEBio cavity fluid inlet (normal velocity, parabolic option).
  pre  : replicate FEParabolicMap (1D P1 Laplace across the inlet, natural BC in z) and give the discrete inlet profile/flux
         and the pressure gradient that the discrete solution must have.
  post RUNDIR : acceptance measurements from the logs.
H = 1 (y in [1,2]), TZ = 0.1 (hex8, 140 x 10 x 1), mu = 0.01 Pa s, rho = 1, K from model_meta.json.
"""
import sys, json, numpy as np
H, TZ, MU, NY = 1.0, 0.1, 0.01, 10
h = H / NY; y = 1 + h * np.arange(NY + 1); s = (y - 1) / H

def pre(ubar=1.0):
    # FEParabolicMap: -phi'' = 1, phi = 0 on rim (wall) nodes; P1 nodal values are exact for a quadratic
    phi = 0.5 * H * H * s * (1 - s)
    vbar = h * phi.sum()                     # Q/A of the linear interpolant (trapezoid)
    VN = phi / vbar * ubar                   # nodal inlet speeds (mean of the interpolant = ubar)
    Qh = TZ * h * VN.sum()                   # discrete flux with the trilinear shape functions (per TZ)
    exact = 1.5 * ubar * 4 * s * (1 - s)     # nominal parabola 1.5 ubar (1-(2s-1)^2)
    # the nodal profile is a parabola through the nodes => continuum flux of that parabola
    Qc = (VN.max() * 2 / 3) * H              # per unit depth, parabola of peak VN.max()
    return dict(VN=VN, Q_discrete_per_depth=Qh / TZ, peak=VN.max(), peak_over_nominal=VN.max() / (1.5 * ubar),
                Q_continuum_of_parabola=Qc, dpdx_discrete=-12 * MU * Qc / H**3, dpdx_nominal=-12 * MU * ubar * H / H**3,
                max_rel_dev_from_nominal=np.max(np.abs(VN - exact)) / (1.5 * ubar))

if __name__ == "__main__":
    if sys.argv[1] == "pre":
        r = pre()
        for k, v in r.items():
            print(k, np.array2string(v, precision=8) if k == "VN" else f"{v:.8g}")


# ---------------------------------------------------------------- post
import re, os

def read_feb_nodes(feb):
    nodes = {}
    for m in re.finditer(r'<node id="(\d+)">([^<]+)</node>', open(feb, encoding="latin-1").read()):
        x, yy, z = map(float, m.group(2).split(","))
        nodes[int(m.group(1))] = (x, yy, z)
    return nodes

def read_mon(path):
    """-> {time: {id: (wx, wy, ef)}}"""
    out, t = {}, None
    for ln in open(path):
        if ln.startswith("*Time"):
            t = float(ln.split("=")[1]); out[t] = {}
        elif not ln.startswith("*") and ln.strip():
            f = ln.split(); out[t][int(f[0])] = tuple(map(float, f[3:6]))
    return out

def post(rd):
    meta = json.load(open(os.path.join(rd, "model_meta.json"))); K = meta["k"]; mu = meta["mu"]
    tag = meta["tag"]; nodes = read_feb_nodes(os.path.join(rd, tag + ".feb"))
    mon = read_mon(os.path.join(rd, "Mon_all.txt")); ts = sorted(mon)
    res = {"tag": tag, "K": K, "mu": mu, "n_times": len(ts), "t_last": ts[-1]}
    # column structure (z = 0 layer): x -> sorted y
    cols = {}
    for i, (x, yv, z) in nodes.items():
        if z == 0.0 and i in mon[ts[-1]]:
            cols.setdefault(round(x, 6), []).append((yv, i))
    xs = sorted(cols)
    for c in cols.values(): c.sort()
    def flux(t, x):
        c = cols[x]; yy = np.array([a for a, _ in c]); w = np.array([mon[t][i][0] for _, i in c])
        return np.trapz(w, yy)           # trapezoid == integral of the trilinear interpolant (per unit depth)
    def pcl(t, x):
        c = cols[x]; i = [i for a, i in c if abs(a - 1.5) < 1e-9][0]
        return -K * mon[t][i][2]
    def pmean(t, x):
        c = cols[x]; yy = np.array([a for a, _ in c]); e = np.array([mon[t][i][2] for _, i in c])
        return -K * np.trapz(e, yy) / H
    # time histories
    Qin = np.array([flux(t, xs[0]) for t in ts]); Qout = np.array([flux(t, xs[-1]) for t in ts])
    res["Qin_t_last"] = Qin[-1]; res["Qout_t_last"] = Qout[-1]
    res["mass_imbalance_rel"] = (Qout[-1] - Qin[-1]) / Qin[-1]
    ramp = np.array([min(t / 10.0, 1.0) for t in ts])
    res["Qin_over_ramp_max_dev"] = float(np.max(np.abs(Qin[ramp >= 1] - 1.0))) if (ramp >= 1).any() else None
    res["max_abs_imbalance_after_ramp"] = float(np.max(np.abs((Qout - Qin)[ramp >= 1]))) if (ramp >= 1).any() else None
    # flux along the channel at the last time
    Qx = np.array([flux(ts[-1], x) for x in xs]); res["Q_along_x_min"] = Qx.min(); res["Q_along_x_max"] = Qx.max()
    # profile at x = 7
    c = cols[min(xs, key=lambda a: abs(a - 7.0))]; yy = np.array([a for a, _ in c]); w7 = np.array([mon[ts[-1]][i][0] for _, i in c]); wy7 = np.array([mon[ts[-1]][i][1] for _, i in c])
    p = np.polyfit(yy, w7, 2); Qc_fit = (p[0] * (1 + 2 ** 0 * 0) * 0)  # placeholder, replaced below
    # continuum flux of the fitted parabola over [1,2]
    a, b, c0 = p; Qc_fit = a * (8 - 1) / 3 + b * (4 - 1) / 2 + c0 * 1
    s = (yy - 1) / H; nominal = 1.5 * 4 * s * (1 - s)
    disc = pre()["VN"]
    res["profile_x7"] = dict(y=yy.tolist(), wx=w7.tolist(), peak=float(w7.max()), peak_over_nominal=float(w7.max() / 1.5),
                             max_abs_wy=float(np.abs(wy7).max()), max_dev_from_nominal=float(np.max(np.abs(w7 - nominal))),
                             max_dev_from_discrete_inlet_profile=float(np.max(np.abs(w7 - disc))),
                             parabola_fit_residual_max=float(np.max(np.abs(np.polyval(p, yy) - w7))), Qc_of_fitted_parabola=float(Qc_fit),
                             trapezoid_Q=float(np.trapz(w7, yy)))
    # pressure gradient from the centre line and the section-mean pressure, windows away from ends
    t = ts[-1]
    for name, fn in (("centreline", pcl), ("section_mean", pmean)):
        xs_a = np.array(xs); pp = np.array([fn(t, x) for x in xs])
        m = (xs_a >= 3) & (xs_a <= 11); g = np.polyfit(xs_a[m], pp[m], 1)
        res[f"dpdx_{name}_3_11"] = float(g[0]); res[f"p_inlet_{name}"] = float(pp[0]); res[f"p_outlet_{name}"] = float(pp[-1])
        res[f"fit_resid_max_{name}"] = float(np.max(np.abs(np.polyval(g, xs_a[m]) - pp[m])))
        loc = np.diff(pp) / np.diff(xs_a)
        res[f"local_gradient_range_{name}"] = [float(loc[(xs_a[:-1] >= 3) & (xs_a[1:] <= 11)].min()), float(loc[(xs_a[:-1] >= 3) & (xs_a[1:] <= 11)].max())]
        res[f"entry_len_{name}"] = float(next((xs_a[j] for j in range(len(loc)) if abs(loc[j] - g[0]) < 0.01 * abs(g[0])), np.nan))
        res[f"exit_len_{name}"] = float(14 - next((xs_a[j] for j in range(len(loc) - 1, -1, -1) if abs(loc[j] - g[0]) < 0.01 * abs(g[0])), np.nan))
        # steady-state change
        t2 = min(ts, key=lambda a: abs(a - (t - 10.0)))
        pp2 = np.array([fn(t2, x) for x in xs]); g2 = np.polyfit(xs_a[m], pp2[m], 1)
        res[f"dpdx_{name}_change_last10s_rel"] = float((g[0] - g2[0]) / g[0])
    # predictions
    Qc = res["profile_x7"]["Qc_of_fitted_parabola"]
    res["dpdx_pred_from_fitted_profile"] = -12 * mu * Qc / H**3
    res["dpdx_pred_discrete_inlet"] = pre()["dpdx_discrete"]; res["dpdx_pred_nominal"] = -0.12
    res["dpdx_rel_err_vs_fitted_profile_centreline"] = res["dpdx_centreline_3_11"] / res["dpdx_pred_from_fitted_profile"] - 1
    res["dpdx_rel_err_vs_fitted_profile_section_mean"] = res["dpdx_section_mean_3_11"] / res["dpdx_pred_from_fitted_profile"] - 1
    res["dpdx_rel_err_vs_discrete_inlet"] = res["dpdx_section_mean_3_11"] / res["dpdx_pred_discrete_inlet"] - 1
    res["dpdx_rel_err_vs_nominal_Q1"] = res["dpdx_section_mean_3_11"] / -0.12 - 1
    # dilatation field (z = 0 layer, all nodes, all times)
    allef = np.array([[mon[tt][i][2] for i in mon[tt]] for tt in ts])
    res["ef_min"] = float(allef.min()); res["ef_max"] = float(allef.max()); res["ef_final_min"] = float(allef[-1].min()); res["ef_final_max"] = float(allef[-1].max())
    res["ef_expected_inlet_abs"] = abs(res["dpdx_pred_from_fitted_profile"]) * 14 / K
    # check that the dilatation is monotone along x in the final state (centreline) -> no checkerboard
    pcl_last = np.array([pcl(t, x) for x in xs]); res["centreline_p_monotone_decreasing"] = bool(np.all(np.diff(pcl_last) < 0))
    # log analysis
    log = open(os.path.join(rd, tag + ".log"), errors="ignore").read()
    its = [int(x) for x in re.findall(r"number of iterations\s*:\s*(\d+)", log)]
    refs = [int(x) for x in re.findall(r"number of reformations\s*:\s*(\d+)", log)]
    res["steps_converged"] = len(its); res["iters_mean_all_steps"] = float(np.mean(its)) if its else None; res["iters_max"] = max(its) if its else None
    res["iters_first_step"] = its[0] if its else None; res["reformations_total"] = int(sum(refs))
    res["retries_or_cutbacks"] = len(re.findall(r"(?i)retrying|cutback|reducing time step", log))
    res["warnings"] = len(re.findall(r"(?i)warning", log)); res["errors"] = len(re.findall(r"(?i)error", log))
    res["negative_jacobian"] = bool(re.search(r"(?i)negative jacobian", log))
    res["normal_termination"] = "N O R M A L   T E R M I N A T I O N" in log
    # per-step final convergence tables. Steps whose log has a single table are "converged" at iteration 0 against a reference that is the
    # step's own (vanishing) increment ("No force acting on the system"): their cur/req margins are degenerate and are reported separately.
    steps_txt = re.split(r"===== beginning time step", log)[1:]
    mg = {"velocity": [], "dilatation": []}; resid = []; n_tabs = []
    for st in steps_txt:
        tabs = st.split("convergence summary")[0].split("convergence norms")[1:]; n_tabs.append(len(tabs))
        if len(tabs) < 2: continue
        last = tabs[-1]
        for key in mg:
            mm = re.search(key + r"\s+(\S+)\s+(\S+)\s+(\S+)", last)
            if mm and float(mm.group(3)) > 0: mg[key].append(float(mm.group(2)) / float(mm.group(3)))
        mr = re.search(r"residual\s+(\S+)\s+(\S+)", last)
        if mr and float(mr.group(1)) > 0: resid.append(float(mr.group(2)) / float(mr.group(1)))
    res["steps_with_newton_iterations"] = int(sum(1 for n in n_tabs if n >= 2)); res["steps_degenerate_zero_iteration"] = int(sum(1 for n in n_tabs if n < 2))
    res["first_degenerate_step"] = next((k for k, n in enumerate(n_tabs, 1) if n < 2), None)
    res["iters_mean_nondegenerate"] = float(np.mean([n for n in n_tabs if n >= 2])) if any(n >= 2 for n in n_tabs) else None
    res["margin_velocity_max_nondegenerate"] = max(mg["velocity"]) if mg["velocity"] else None
    res["margin_dilatation_max_nondegenerate"] = max(mg["dilatation"]) if mg["dilatation"] else None
    res["final_residual_over_initial_max_nondegenerate"] = max(resid) if resid else None
    res["warning_text_distinct"] = sorted(set(re.findall(r"\*\s+(No force acting on the system\.)", log)))
    # checkerboard (odd-even) pressure modes in the final state (z = 0 layer)
    t = ts[-1]
    def pnode(x, yv): return -K * [mon[t][i][2] for a_, i in cols[x] if abs(a_ - yv) < 1e-9][0]
    xs_a = np.array(xs); m = (xs_a >= 2) & (xs_a <= 12); idx = np.round(xs_a[m] / h).astype(int)
    def saw_x(yv):
        pp = np.array([pnode(x, yv) for x in xs_a[m]]); return float(pp[idx % 2 == 0].mean() - pp[idx % 2 == 1].mean())
    res["sawtooth_x_amp_even_minus_odd_Pa"] = {str(yv): saw_x(yv) for yv in (1.2, 1.5, 1.9)}
    res["sawtooth_x_amp_over_inlet_outlet_drop"] = abs(res["sawtooth_x_amp_even_minus_odd_Pa"]["1.5"]) / abs(res["p_inlet_section_mean"] - res["p_outlet_section_mean"])
    res["sawtooth_x_amp_over_per_cell_drop"] = abs(res["sawtooth_x_amp_even_minus_odd_Pa"]["1.5"]) / (abs(res["dpdx_pred_from_fitted_profile"]) * h)
    xc = xs[int(np.argmin([abs(a_ - 7.0) for a_ in xs]))]
    col = [pnode(xc, 1 + h * j) for j in range(NY + 1)]
    res["sawtooth_y_amp_x7_Pa"] = float(np.mean(col[0::2]) - np.mean(col[1::2]))
    sh = {}
    for tt in (5.0, 10.0, 15.0, 20.0, 50.0, 150.0):
        t_ = min(ts, key=lambda a_: abs(a_ - tt)); pp = np.array([-K * [mon[t_][i][2] for a_, i in cols[x] if abs(a_ - 1.5) < 1e-9][0] for x in xs_a[m]])
        sh[str(t_)] = float(pp[idx % 2 == 0].mean() - pp[idx % 2 == 1].mean())
    res["sawtooth_x_amp_history_Pa"] = sh
    # smoothed (even/odd-averaged) pressure gradient: mean of neighbouring pairs removes the odd-even mode exactly for a linear field
    pc_all = np.array([pcl(t, x) for x in xs_a]); sm = 0.5 * (pc_all[1:] + pc_all[:-1]); xm = 0.5 * (xs_a[1:] + xs_a[:-1]); mm_ = (xm >= 3) & (xm <= 11)
    res["dpdx_centreline_pair_averaged_3_11"] = float(np.polyfit(xm[mm_], sm[mm_], 1)[0]); res["pair_avg_fit_resid_max_Pa"] = float(np.max(np.abs(np.polyval(np.polyfit(xm[mm_], sm[mm_], 1), xm[mm_]) - sm[mm_])))
    # end-of-channel disturbance: deviation of the pair-averaged pressure from the fitted line
    gline = np.polyfit(xm[mm_], sm[mm_], 1); dev = sm - np.polyval(gline, xm)
    res["pair_avg_dev_max_x_lt_1_Pa"] = float(np.max(np.abs(dev[xm < 1]))); res["pair_avg_dev_max_x_gt_13_Pa"] = float(np.max(np.abs(dev[xm > 13])))
    res["pair_avg_dev_max_3_11_Pa"] = float(np.max(np.abs(dev[mm_])))
    res["p_inlet_pair_avg_line_Pa"] = float(np.polyval(gline, 0.0)); res["p_outlet_pair_avg_line_Pa"] = float(np.polyval(gline, 14.0))
    mt = re.search(r"Total elapsed time[^:]*:\s*[\d:]+\s*\(([\d.]+) sec\)", log); res["febio_wall_s"] = float(mt.group(1)) if mt else None
    tm = open(os.path.join(rd, "time.txt")).read()
    mw = re.search(r"Elapsed \(wall clock\) time.*?: ([\d:.]+)", tm); mu_ = re.search(r"User time \(seconds\): ([\d.]+)", tm); ms = re.search(r"System time \(seconds\): ([\d.]+)", tm); mr_ = re.search(r"Maximum resident set size \(kbytes\): (\d+)", tm)
    res["time_wall"] = mw.group(1) if mw else None; res["user_s"] = float(mu_.group(1)); res["sys_s"] = float(ms.group(1)); res["max_rss_MB"] = int(mr_.group(1)) / 1024
    res["mesh"] = dict(n_nodes=meta["n_nodes"], n_elem=meta["n_elem"])
    # mesh quality: hex8 uniform 0.1 x 0.1 x 0.1 cells
    return res

if __name__ == "__main__" and sys.argv[1] == "post":
    r = post(sys.argv[2]); json.dump(r, open(os.path.join(sys.argv[2], "poiseuille_check.json"), "w"), indent=1, default=float)
    for k, v in r.items():
        if k != "profile_x7": print(k, v)
    print("profile_x7", {k: v for k, v in r["profile_x7"].items() if k not in ("y", "wx")})
