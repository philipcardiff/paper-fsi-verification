"""Regenerate the matched-comparison tables, statistics and figures from the PACKAGED data only (relative paths).
usage (from the package root paper-data/fenicsx_hessenthaler/): python3 -B scripts/regenerate_from_packaged_data.py [OUTDIR]   (numpy, matplotlib)
Writes: OUTDIR/{matched_comparison.csv, summary_stats.json, centreline_0p69_stats.json, tables_matched.tex, fig_*.pdf}"""
import csv, json, sys, os
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
D = "data/"; OUT = (sys.argv[1] if len(sys.argv) > 1 else "regenerated").rstrip("/") + "/"; os.makedirs(OUT, exist_ok=True)
rd = lambda p: list(csv.DictReader(open(p)))
def series(paths):
    d = {}
    for p in paths:
        for x in rd(D + p): d[round(float(x["t"]), 6)] = x
    return d
gmi = series(["Gmi_n32b_history.csv", "Gmi_dt01_history.csv"]); gc = series(["Gc_T3b_history.csv"])
s4 = {round(float(x["t"]), 6): x for x in rd(D + "s4f_feF1_ours_0_0.75.csv")}
# ---- matched table (0.1 ... 0.75)
mt = []
for t in (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75):
    s = s4[t]; r = dict(t=t, s4f_tip_mm=float(s["tip_dy_m"]) * 1e3, s4f_Fy_on_solid_uN=float(s["Fs_y"]) * 1e6, s4f_Fz_on_solid_uN=float(s["Fs_z"]) * 1e6,
                        s4f_Fy_fluidside_signflipped_uN=-float(s["Ff_y"]) * 1e6, s4f_Fz_fluidside_signflipped_uN=-float(s["Ff_z"]) * 1e6, s4f_coupling_its=s["its"])
    for nm, d in (("Gmi", gmi), ("Gc", gc)):
        if t not in d: continue
        x = d[t]; r[nm + "_tip_mm"] = float(x["tip_y"]) * 1e3; r[nm + "_Fy_Fdir_uN"] = float(x["Fy_fluid"]) * 1e6; r[nm + "_Fz_Fdir_uN"] = float(x["Fz_fluid"]) * 1e6
        r[nm + "_tip_rel_diff_vs_s4f"] = r[nm + "_tip_mm"] / r["s4f_tip_mm"] - 1
        r[nm + "_Fy_rel_diff_vs_s4f"] = r[nm + "_Fy_Fdir_uN"] / r["s4f_Fy_on_solid_uN"] - 1; r[nm + "_Fz_rel_diff_vs_s4f"] = r[nm + "_Fz_Fdir_uN"] / r["s4f_Fz_on_solid_uN"] - 1
    mt.append(r)
keys = sorted({k for r in mt for k in r}, key=lambda k: (k != "t", k))
with open(OUT + "matched_comparison.csv", "w") as f:
    w = csv.DictWriter(f, fieldnames=keys); w.writeheader(); [w.writerow(r) for r in mt]
# ---- all common times
tt = sorted(t for t in gmi if t in s4 and 0.02 <= t <= 0.75 + 1e-9); rows = []
for t in tt:
    g, s = gmi[t], s4[t]; ts = float(s["tip_dy_m"]) * 1e3; tg = float(g["tip_y"]) * 1e3
    rows.append(dict(t=t, s4f_tip_mm=ts, Gmi_tip_mm=tg, rel=(tg / ts - 1 if ts > 0.01 else None)))
rel = np.array([r["rel"] for r in rows if r["rel"] is not None and r["t"] >= 0.1])
st = dict(mean=float(rel.mean()), rms=float(np.sqrt((rel ** 2).mean())), max_abs=float(np.abs(rel).max()), n=int(len(rel)))
json.dump({"tip_rel_diff_Gmi_vs_s4f_t_ge_0p1": st, "n_common_times": len(rows)}, open(OUT + "summary_stats.json", "w"), indent=1)
# ---- centreline
cs = np.array([[float(v) for v in r.values()] for r in rd(D + "centreline_0p69_s4f.csv")])
cf = np.array([[float(v) for v in r.values()] for r in rd(D + "centreline_0p69_FEniCSx_Gmi.csv")])
z0, dy, dz = cs[:, 0], cs[:, 2] * 1e3, cs[:, 3] * 1e3; zf, dyf, dzf = cf[:, 0] * 1e3, cf[:, 2] * 1e3, cf[:, 3] * 1e3
yi = np.interp(z0, zf, dyf)
cl = dict(t=0.69, n_points=len(z0), rms_mm=float(np.sqrt(np.mean((yi - dy) ** 2))), max_abs_mm=float(np.abs(yi - dy).max()), tip_y_s4f_mm=float(dy[-1]), tip_y_FE_mm=float(dyf[-1]))
json.dump(cl, open(OUT + "centreline_0p69_stats.json", "w"), indent=1)
# ---- tex tables
L = ["% regenerated; do not edit", "\\newcommand{\\matchedtable}{%", "\\begin{tabular}{rrrrrrrrr}\\toprule",
     "$t$ (s) & tip s4f & tip FE Gmi & diff (\\%) & $F_y$ s4f & $F_y$ FE & $F_z$ s4f & $F_z$ FE & s4f its\\\\",
     "& (mm) & (mm) & & ($\\mu$N) & ($\\mu$N) & ($\\mu$N) & ($\\mu$N) & \\\\ \\midrule"]
for r in mt:
    L.append("%.2f & %.4f & %.4f & %+.1f & %.1f & %.1f & %.1f & %.1f & %s\\\\" % (r["t"], r["s4f_tip_mm"], r["Gmi_tip_mm"], 100 * r["Gmi_tip_rel_diff_vs_s4f"],
             r["s4f_Fy_on_solid_uN"], r["Gmi_Fy_Fdir_uN"], r["s4f_Fz_on_solid_uN"], r["Gmi_Fz_Fdir_uN"], r["s4f_coupling_its"]))
L.append("\\bottomrule\\end{tabular}}")
L.append("\\newcommand{\\tipstats}{mean %+.1f\\,\\%%, rms %.1f\\,\\%%, maximum %.1f\\,\\%% (%d common times, $0.1\\le t\\le0.75$\\,s)}" % (100 * st["mean"], 100 * st["rms"], 100 * st["max_abs"], st["n"]))
L.append("\\newcommand{\\clstats}{rms %.3f\\,mm, maximum %.3f\\,mm (%d points)}" % (cl["rms_mm"], cl["max_abs_mm"], cl["n_points"]))
open(OUT + "tables_matched.tex", "w").write("\n".join(L) + "\n")
# ---- figures
t4 = np.array([float(x["t"]) for x in s4.values()]); y4 = np.array([float(x["tip_dy_m"]) for x in s4.values()]) * 1e3
fig, ax = plt.subplots(1, 2, figsize=(11, 3.8))
ax[0].plot(t4, y4, "-", color="k", lw=1.5, label="solids4foam feF1 reconstruction (OURS, dt 0.002 s)")
ax[0].plot([r["t"] for r in rows], [r["Gmi_tip_mm"] for r in rows], "o", ms=3, label="FEniCSx Gmi (118k tets, dt 0.02/0.01 s)")
gct = [t for t in sorted(gc) if t <= 0.75]; ax[0].plot(gct, [float(gc[t]["tip_y"]) * 1e3 for t in gct], ".", ms=3, label="FEniCSx Gc (dev level)")
ax[0].set_xlabel("t (s)"); ax[0].set_ylabel("tip y displacement (mm)"); ax[0].legend(fontsize=7); ax[0].set_title("matched physical times; transient, no steady value", fontsize=8)
ax[1].plot([r["t"] for r in rows if r["rel"] is not None and r["t"] >= 0.1], [100 * r["rel"] for r in rows if r["rel"] is not None and r["t"] >= 0.1], "o-", ms=3)
ax[1].axhline(0, color="k", lw=0.5); ax[1].set_xlabel("t (s)"); ax[1].set_ylabel("(FE Gmi / s4f - 1) tip (%)"); ax[1].set_title("relative difference of the tip", fontsize=8)
fig.tight_layout(); fig.savefig(OUT + "fig_tip_matched.pdf")
fig, ax = plt.subplots(1, 2, figsize=(11, 3.8))
for a, k, lab in ((ax[0], "y", "Fy"), (ax[1], "z", "Fz")):
    a.plot(t4, np.array([float(x["Fs_" + k]) for x in s4.values()]) * 1e6, "-", color="k", label="s4f 'solid interface' force (OURS)")
    a.plot(tt, [float(gmi[t][lab + "_fluid"]) * 1e6 for t in tt], "o", ms=3, label="FEniCSx Gmi facet traction F_dir")
    a.set_xlabel("t (s)"); a.set_ylabel(lab + " on the flap (uN)"); a.legend(fontsize=7)
ax[0].set_title("force on the solid; sign: fluid acting on the flap", fontsize=8); fig.tight_layout(); fig.savefig(OUT + "fig_forces_matched.pdf")
fig, ax = plt.subplots(figsize=(6, 3.8)); ax.plot(z0 + dz, dy, "k.-", label="solids4foam (OURS) t = 0.69 s"); ax.plot(zf + dzf, dyf, "-", label="FEniCSx Gmi t = 0.69 s")
ax.set_xlabel("z (mm)"); ax.set_ylabel("y (mm)"); ax.legend(fontsize=8); fig.tight_layout(); fig.savefig(OUT + "fig_centreline_0p69.pdf")
print(json.dumps(st), json.dumps(cl))

# ---- temporal control (Addendum A7): Gc dt 0.005 vs dt 0.02 vs s4f
if os.path.exists(D + "Gc_dt005_history.csv"):
    g5 = series(["Gc_dt005_history.csv"]); com = [t for t in sorted(g5) if t in gc and 0.1 <= t <= 0.75]
    tr = {}
    r1 = [g5[t]["tip_y"] for t in com]
    rel_tip = [(t, float(g5[t]["tip_y"]) / float(gc[t]["tip_y"]) - 1) for t in com]
    rel_tip2 = [x for x in rel_tip if x[0] >= 0.2]
    def peak(d, key, scale=1.0):
        c = [(float(d[t][key]) * 1e6 * scale, t) for t in d if 0.6 <= t <= 0.75]; return max(c)
    pk_s4 = max((float(s4[t]["Fs_y"]) * 1e6, t) for t in s4 if 0.6 <= t <= 0.75)
    pk_02 = peak(gc, "Fy_fluid"); pk_005 = peak(g5, "Fy_fluid"); pk_gmi = peak(gmi, "Fy_fluid")
    moved = (pk_005[0] - pk_02[0]) / (pk_s4[0] - pk_02[0])
    tr = dict(n_common_times=len(com), tip_rel_diff_dt005_vs_dt02_max_all=max(abs(x[1]) for x in rel_tip), tip_rel_diff_dt005_vs_dt02_max_t_ge_0p2=max(abs(x[1]) for x in rel_tip2),
              tip_rel_diff_dt005_vs_dt02_rms_t_ge_0p2=float(np.sqrt(np.mean([x[1] ** 2 for x in rel_tip2]))),
              Fy_peak_uN=dict(s4f=pk_s4, FE_Gc_dt002=pk_02, FE_Gc_dt0005=pk_005, FE_Gmi_dt002_001=pk_gmi), fraction_of_gap_closed_by_dt0005=float(moved))
    json.dump(tr, open(OUT + "temporal_control_stats.json", "w"), indent=1); print(json.dumps(tr))
    fig, ax = plt.subplots(1, 2, figsize=(11, 3.8))
    ts_ = [t for t in sorted(g5) if t <= 0.75]
    ax[0].plot(t4, np.array([float(x["Fs_y"]) for x in s4.values()]) * 1e6, "k-", label="s4f (OURS) dt 0.002 s")
    ax[0].plot(ts_, [float(g5[t]["Fy_fluid"]) * 1e6 for t in ts_], "-", label="FEniCSx Gc dt 0.005 s"); ax[0].plot([t for t in sorted(gc) if t <= 0.75], [float(gc[t]["Fy_fluid"]) * 1e6 for t in sorted(gc) if t <= 0.75], "o", ms=3, label="FEniCSx Gc dt 0.02 s")
    ax[0].set_xlim(0.4, 0.75); ax[0].set_xlabel("t (s)"); ax[0].set_ylabel("Fy on the flap (uN)"); ax[0].legend(fontsize=7); ax[0].set_title("temporal control (Gc dev level)", fontsize=8)
    ax[1].plot(t4, y4, "k-", label="s4f (OURS)"); ax[1].plot(ts_, [float(g5[t]["tip_y"]) * 1e3 for t in ts_], "-", label="FEniCSx Gc dt 0.005 s"); ax[1].plot([t for t in sorted(gc) if t <= 0.75], [float(gc[t]["tip_y"]) * 1e3 for t in sorted(gc) if t <= 0.75], "o", ms=3, label="FEniCSx Gc dt 0.02 s")
    ax[1].set_xlabel("t (s)"); ax[1].set_ylabel("tip y (mm)"); ax[1].legend(fontsize=7)
    fig.tight_layout(); fig.savefig(OUT + "fig_temporal_control.pdf")

    # ---- LaTeX: temporal table and paragraph (numbers from the statistics above)
    T = ["% regenerated; do not edit", "\\newcommand{\\temporaltable}{%", "\\begin{tabular}{rrrrrrr}\\toprule",
         "$t$ (s) & tip s4f & tip Gc $\\Delta t{=}0.005$ & tip Gc $\\Delta t{=}0.02$ & $F_y$ s4f & $F_y$ Gc 0.005 & $F_y$ Gc 0.02\\\\", "& (mm) & (mm) & (mm) & ($\\mu$N) & ($\\mu$N) & ($\\mu$N) \\\\ \\midrule"]
    for t in (0.1, 0.3, 0.5, 0.6, 0.62, 0.64, 0.66, 0.7, 0.75):
        f = lambda d, k, sc: ("%.1f" % (float(d[t][k]) * sc) if t in d else "--")
        g = lambda d, k, sc: ("%.4f" % (float(d[t][k]) * sc) if t in d else "--")
        T.append("%.2f & %s & %s & %s & %s & %s & %s\\\\" % (t, g(s4, "tip_dy_m", 1e3), g(g5, "tip_y", 1e3), g(gc, "tip_y", 1e3), f(s4, "Fs_y", 1e6), f(g5, "Fy_fluid", 1e6), f(gc, "Fy_fluid", 1e6)))
    T.append("\\bottomrule\\end{tabular}}")
    open(OUT + "tables_temporal.tex", "w").write("\n".join(T) + "\n")
    para = ("\\paragraph{Temporal control (development level).} Because the FEniCSx step ($0.02\\,\\mathrm{s}$) is ten times the solids4foam step, the Gc level was rerun from $t=0$ with $\\Delta t=0.005\\,\\mathrm{s}$ (same mesh and numerics, preregistered question and criterion in \\texttt{prereg/ADDENDUM\\_A7.md}). The tip changes by at most %.1f\\,\\%% (rms %.1f\\,\\%%) for $t\\ge0.2\\,\\mathrm{s}$ (and by %.0f\\,\\%% at $t=0.1\\,\\mathrm{s}$, where the tip is micrometres), and the peak of $F_y$ near $0.67\\,\\mathrm{s}$ rises from %.0f to %.0f\\,$\\mu$N, against %.0f\\,$\\mu$N in solids4foam (%.0f\\,\\%% of the gap closed): (sample maxima: the $\\Delta t=0.02\\,\\mathrm{s}$ peak is under-sampled, so this fraction is uncertain); on the accepted level Gmi the sample maximum of $F_y$ is only %.0f\\,$\\mu$N at $t=%.2f\\,\\mathrm{s}$, about half of the solids4foam peak, and $F_y$ has the opposite sign at $0.7\\,\\mathrm{s}$ (Table~\\ref{tab:matched}). The height of the $F_y$ peak is therefore largely a time-step effect in the FEniCSx data at $\\Delta t=0.02\\,\\mathrm{s}$ (shown on the development level only), and the tip agreement of Table~\\ref{tab:matched} contains temporal errors of the same order as the differences it reports (at $\\Delta t=0.005\\,\\mathrm{s}$ the development-level tip is $-6\\,\\%%$ below solids4foam at $0.6\\,\\mathrm{s}$ and $+2\\,\\%%$ above at $0.75\\,\\mathrm{s}$). This is one mesh and two time steps: it is neither a temporal convergence study nor an error bound. The step used in Table~\\ref{tab:matched} was not changed to improve agreement.\\par\\begin{table}[htbp]\\centering\\caption{Temporal control on the development level Gc (not an accepted level) against solids4foam (ours).}\\label{tab:temporal}\\scriptsize\\setlength{\\tabcolsep}{3pt}\\temporaltable\\end{table}"
            % (100 * tr["tip_rel_diff_dt005_vs_dt02_max_t_ge_0p2"], 100 * tr["tip_rel_diff_dt005_vs_dt02_rms_t_ge_0p2"], 100 * tr["tip_rel_diff_dt005_vs_dt02_max_all"], tr["Fy_peak_uN"]["FE_Gc_dt002"][0], tr["Fy_peak_uN"]["FE_Gc_dt0005"][0], tr["Fy_peak_uN"]["s4f"][0], 100 * tr["fraction_of_gap_closed_by_dt0005"], tr["Fy_peak_uN"]["FE_Gmi_dt002_001"][0], tr["Fy_peak_uN"]["FE_Gmi_dt002_001"][1]))
    open(OUT + "temporal_paragraph.tex", "w").write(para + "\n")
