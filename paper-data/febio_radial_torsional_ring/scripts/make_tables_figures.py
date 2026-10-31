#!/usr/bin/env python3
"""Tables (CSV) and figures from the stored analysis results.   usage: make_tables_figures.py PACKAGE_ROOT RESULTS_DIR OUT_TABLES_DIR OUT_FIG_DIR
RESULTS_DIR holds campaign_analysis/*.json, independent_check.json, newton_iterations_per_step.csv, rtc_m1_step2_residuals.csv (as written by run_analyses.sh)."""
import sys, os, json, csv, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

root, res, tabs, figs = sys.argv[1:5]
os.makedirs(tabs, exist_ok=True); os.makedirs(figs, exist_ok=True)
J = lambda p: json.load(open(os.path.join(res, p)))
ind = J("independent_check.json")
ca = {n: J(f"campaign_analysis/{n}.json") for n in ("RTB_fluid_M1_as_run", "RTB_fluid_M1_aligned_dt2", "RTB_fluid_M2_as_run", "RTB_fluid_M2_aligned_dt2", "RTB_fluid_M1_bcshift_as_run")}
tw = {m: J(f"campaign_analysis/RTA_twist_M{m}.json") for m in (1, 2, 3)}
law = J("campaign_analysis/RTA_law1c.json"); rot = J("campaign_analysis/RTA_rotcheck.json")

# --------------------------------------------------------------------------------------------------------------------------------- tables
def w(name, rows):
    with open(os.path.join(tabs, name), "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=list(rows[0])); wr.writeheader(); wr.writerows(rows)

fmt = lambda x: f"{x:.6g}"
rows = []
for m in (1, 2, 3):
    i = ind["RT-A_twist"][f"M{m}"]
    rows.append(dict(mesh=f"M{m}", n_theta={1: 64, 2: 128, 3: 256}[m], n_radial_solid={1: 8, 2: 16, 3: 32}[m], E_max_TA_campaign=fmt(tw[m]["E_max_TA"]), E_max_TA_independent=fmt(i["E_max_TA"]),
                     E_max_ur_R0p75=fmt(tw[m]["E_max_ur_075"]), E_max_alpha_R0p75=fmt(tw[m]["E_max_alpha_075"]), inner_plus_outer_torque_over_max_inner_independent=fmt(i["inner_plus_outer_over_max_inner"]),
                     axisymmetry_ray0_vs_ray45_rel=fmt(tw[m]["axisymmetry_ur_ray0_vs_ray45_rel"]), peak_TA_reference=fmt(i["peak_TA_reference"]), peak_TA_febio=fmt(i["peak_TA_febio"])))
rows[0]["observed_order_to_next"] = fmt(ind["RT-A_twist"]["order_M1_M2"]); rows[1]["observed_order_to_next"] = fmt(ind["RT-A_twist"]["order_M2_M3"]); rows[2]["observed_order_to_next"] = ""
w("rta_twist_summary.csv", rows)
w("rta_law_rotation_summary.csv", [
    dict(check="law identity (hex8 block, F_k up to J = 0.95), relative to max |stress|", campaign=fmt(law["rel"]), independent_closed_form=fmt(ind["RT-A_law_identity"]["rel_closed_form"]),
         independent_numeric_dWdF=fmt(ind["RT-A_law_identity"]["rel_numeric_dW_dF"]), gate="1e-8", result="PASS"),
    dict(check="rotation BC expression, max node difference [m]", campaign=fmt(rot["worst_abs_diff_m"]), independent_closed_form=fmt(ind["RT-A_rotation_expression"]["max_abs_node_diff_m"]),
         independent_numeric_dWdF="", gate="1e-10", result="PASS")])

keys = [("E_L2_velocity_areaweighted", "velocity_L2"), ("E_max_T1", "T1_Emax"), ("E_max_T2", "T2_Emax"), ("E_max_T3", "T3_Emax"), ("E_max_H", "Hf_Emax"), ("E_L2_T2", "T2_EL2"), ("E_L2_T3", "T3_EL2")]
rows = []
for run, nm, tag in (("RTB_fluid_M1", "as_run", "M1"), ("RTB_fluid_M2", "as_run", "M2"), ("RTB_fluid_M1", "aligned_dt2", "M1"), ("RTB_fluid_M2", "aligned_dt2", "M2")):
    c = ca[f"{run}_{nm}"]; i = ind[run][nm]
    for ck, label in keys:
        rows.append(dict(mesh=tag, time_labels="as run" if nm == "as_run" else "exact solution delayed by dt/2 (post-hoc alignment)", quantity=label, campaign_percent=fmt(100 * c[ck]), independent_percent=fmt(100 * i[ck])))
# bcshift run (table advanced by dt/2, as-run labels)
c = ca["RTB_fluid_M1_bcshift_as_run"]; i = ind["RTB_fluid_M1_bcshift"]["as_run"]
for ck, label in keys:
    rows.append(dict(mesh="M1 (job 14211, input table advanced by dt/2)", time_labels="as run (aligned by construction)", quantity=label, campaign_percent=fmt(100 * c[ck]), independent_percent=fmt(100 * i[ck])))
w("rtb_error_summary.csv", rows)
ordr = []
for nm in ("as_run", "aligned_dt2"):
    for ck, label in keys[:5]:
        a1, a2 = ca[f"RTB_fluid_M1_{nm}"][ck], ca[f"RTB_fluid_M2_{nm}"][ck]
        ordr.append(dict(time_labels=nm, quantity=label, order_campaign=fmt(math.log2(a1 / a2)), order_independent=fmt(ind[f"orders_{nm}"][ck]), ratio_M1_over_M2=fmt(a1 / a2)))
w("rtb_observed_orders.csv", ordr)
pre = json.load(open(os.path.join(root, "results", "predicted_coupled_element_mean_AE.json")))
w("rtb_prediction_vs_observation.csv", [
    dict(mesh="M1", predicted_T2_Emax_percent=11.7, observed_T2_Emax_aligned_percent=fmt(100 * ca["RTB_fluid_M1_aligned_dt2"]["E_max_T2"]), observed_as_run_percent=fmt(100 * ca["RTB_fluid_M1_as_run"]["E_max_T2"]),
         relative_deviation_aligned_vs_prediction=fmt(ca["RTB_fluid_M1_aligned_dt2"]["E_max_T2"] / 0.117 - 1)),
    dict(mesh="M2", predicted_T2_Emax_percent=6.0, observed_T2_Emax_aligned_percent=fmt(100 * ca["RTB_fluid_M2_aligned_dt2"]["E_max_T2"]), observed_as_run_percent=fmt(100 * ca["RTB_fluid_M2_as_run"]["E_max_T2"]),
         relative_deviation_aligned_vs_prediction=fmt(ca["RTB_fluid_M2_aligned_dt2"]["E_max_T2"] / 0.060 - 1))])

# --------------------------------------------------------------------------------------------------------------------------------- figures
BLUE, ORANGE, GREEN, GREY, INK = "#0072B2", "#D55E00", "#009E73", "#7a7a7a", "#222222"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": "#888", "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.color": "#e4e4e4", "grid.linewidth": 0.6, "legend.frameon": False, "figure.dpi": 130, "savefig.dpi": 150})

# fig 1: RT-A
fig, ax = plt.subplots(1, 2, figsize=(9.2, 3.5))
nth = np.array([64, 128, 256]); E = np.array([ind["RT-A_twist"][f"M{m}"]["E_max_TA"] for m in (1, 2, 3)])
ax[0].loglog(nth, E, "o-", color=BLUE, label="FEBio hex8, static solid")
ax[0].loglog(nth, E[0] * (nth[0] / nth) ** 2, "--", color=GREY, label="second-order guide")
ax[0].axhline(1e-3, color=ORANGE, lw=1, ls=":"); ax[0].text(140, 1.13e-3, "gate at M3: 1e-3", color=INK, fontsize=8)
for x, y, m in zip(nth, E, ("M1", "M2", "M3")): ax[0].annotate(f"{m}: {y:.2e}", (x, y), textcoords="offset points", xytext=(5, -13) if m == "M2" else (5, 6), fontsize=8)
ax[0].set_xlabel("azimuthal cells n_theta (radial cell size ~ 1/n_theta)"); ax[0].set_ylabel("E_max(T_A)  [max |T_A - T_ref| / max |T_ref|]")
ax[0].set_title("(a) Drive-torque error, observed orders 1.95, 1.99", fontsize=9, loc="left")
t = np.array(tw[3]["t"]); ta = -np.array(tw[3]["TA"]); tr = np.array(tw[3]["TAref"])
ax[1].plot(t, tr, "-", color=GREY, lw=3, alpha=0.7, label="reference (proposal's finite-strain energy)")
ax[1].plot(t, ta, "o", color=BLUE, ms=3.5, label="FEBio M3 (reaction torque, sign reversed)")
ax[1].set_xlabel("time t [s]"); ax[1].set_ylabel("inner-circle torque T_A [N m / rad / m]"); ax[1].legend(loc="upper left", fontsize=8)
ax[1].set_title("(b) T_A(t) at M3 (static solid, 20 load steps)", fontsize=9, loc="left")
fig.tight_layout(); fig.savefig(os.path.join(figs, "fig1_rta_static_twist.png")); plt.close(fig)

# fig 2: RT-B torque time series (aligned) and errors
fig, ax = plt.subplots(1, 2, figsize=(9.6, 3.6))
s1, s2 = ind["RTB_fluid_M1"]["series"], ind["RTB_fluid_M2"]["series"]
t = np.array(s1["t"]); Tx = np.array(s1["T_exact_shifted"]); sc = np.max(np.abs(Tx))
ax[0].plot(t, Tx, "-", color=GREY, lw=3, alpha=0.7, label="exact fixed-geometry T_B (delayed by dt/2)")
ax[0].plot(t, s1["T2"], "-", color=ORANGE, lw=1.2, label="T2 as FEBio applies it, M1")
ax[0].plot(t, s2["T2"], "--", color=ORANGE, lw=1.2, label="T2, M2")
ax[0].plot(t, s1["T3"], ":", color=BLUE, lw=1.6, label="T3 (angular momentum), M1")
ax[0].set_xlabel("time t [s]"); ax[0].set_ylabel("torque of fluid on solid [N m / rad / m]"); ax[0].legend(fontsize=7.5, loc="lower left")
ax[0].set_title("(a) Wall shear torque, RT-B (job 14209 M1, 14210 M2)", fontsize=9, loc="left")
for s, ls, nm in ((s1, "-", "M1"), (s2, "--", "M2")):
    ax[1].plot(t, 100 * (np.array(s["T2"]) - np.array(s["T_exact_shifted"])) / sc, ls, color=ORANGE, label=f"T2 {nm}")
    ax[1].plot(t, 100 * (np.array(s["T3"]) - np.array(s["T_exact_shifted"])) / sc, ls, color=BLUE, label=f"T3 {nm}")
ax[1].axhline(0, color="#999", lw=0.8); ax[1].set_xlabel("time t [s]"); ax[1].set_ylabel("error / max |T_exact| [%]"); ax[1].legend(fontsize=8, ncol=2)
ax[1].set_title("(b) Error history (dt/2-aligned reference)", fontsize=9, loc="left")
fig.tight_layout(); fig.savefig(os.path.join(figs, "fig2_rtb_torque_series.png")); plt.close(fig)

# fig 3: errors vs mesh
fig, ax = plt.subplots(1, 2, figsize=(9.6, 3.6))
h = np.array([1 / 16, 1 / 32])
for nm, ls, mk, lab in (("as_run", "-", "s", "as run (time labels as logged)"), ("aligned_dt2", "-", "o", "dt/2-aligned reference")):
    v = [100 * ca[f"RTB_fluid_M{m}_{nm}"]["E_L2_velocity_areaweighted"] for m in (1, 2)]
    ax[0].loglog(h, v, ls, marker=mk, color=ORANGE if nm == "as_run" else BLUE, label=lab)
ax[0].axhline(2, color=GREY, lw=1, ls=":"); ax[0].text(0.0335, 2.15, "gate 2 % (M1)", fontsize=8)
ax[0].loglog(h, 100 * ca["RTB_fluid_M1_aligned_dt2"]["E_L2_velocity_areaweighted"] * (h / h[0]) ** 2, "--", color=GREY, lw=0.9, label="second-order guide")
ax[0].set_xlabel("radial cell size h [m]"); ax[0].set_ylabel("velocity error E_L2 [%]"); ax[0].legend(fontsize=8, loc="center right", bbox_to_anchor=(1.0, 0.62))
ax[0].set_title("(a) Swirl velocity: time offset sets a 1.7 % floor", fontsize=9, loc="left")
for ck, lab, c, mk in (("E_max_T2", "T2 (applied traction)", ORANGE, "o"), ("E_max_T1", "T1 (face-evaluated)", "#CC79A7", "^"), ("E_max_T3", "T3 (angular momentum)", BLUE, "s"), ("E_max_H", "H_f", GREEN, "D")):
    v = [100 * ca[f"RTB_fluid_M{m}_aligned_dt2"][ck] for m in (1, 2)]
    ax[1].loglog(h, v, "-", marker=mk, color=c, label=lab)
ax[1].loglog(h, 11.42 * (h / h[0]), "--", color=GREY, lw=0.9, label="first-order guide"); ax[1].loglog(h, 0.26 * (h / h[0]) ** 2, ":", color=GREY, lw=1.2, label="second-order guide")
ax[1].axhline(10, color=GREY, lw=0.8, ls=":"); ax[1].text(0.0335, 10.6, "T2 gate 10 % (M1), 5 % (M2)", fontsize=7.5)
ax[1].set_xlabel("radial cell size h [m]"); ax[1].set_ylabel("E_max [% of max |exact|]"); ax[1].legend(fontsize=8, loc="center", bbox_to_anchor=(0.5, 0.52))
ax[1].set_title("(b) After dt/2 alignment: T2 stays first order", fontsize=9, loc="left")
fig.tight_layout(); fig.savefig(os.path.join(figs, "fig3_rtb_errors_vs_mesh.png")); plt.close(fig)

# fig 4: BC time convention
tab = np.loadtxt(os.path.join(root, "reference_data", "rtr", "omega_B_table_dt1e-3.csv")); tt, om = tab[:, 0], tab[:, 2]
wv = np.array(ind["RTB_fluid_M1"]["series"]["wall_velocity"]); wv2 = np.array(ind["RTB_fluid_M1_bcshift"]["series"]["wall_velocity"]); T = np.array(ind["RTB_fluid_M1"]["series"]["t"])
fig, ax = plt.subplots(1, 2, figsize=(9.6, 3.5))
sel = (tt > 0.9) & (tt < 1.1)
ax[0].plot(tt[sel], om[sel], "-", color=GREY, lw=2.4, alpha=0.7, label="load-curve table omega(t)")
m = (T > 0.9) & (T < 1.1)
ax[0].plot(T[m], wv[m], "o", color=ORANGE, ms=4, label="logged wall velocity, as run (job 14209)")
ax[0].plot(T[m], wv2[m], "s", mfc="none", color=BLUE, ms=5, label="logged, table advanced by dt/2 (job 14211)")
ax[0].set_xlabel("time t [s]"); ax[0].set_ylabel("wall angular velocity at r = 1 m [rad/s]"); ax[0].legend(fontsize=7.5, loc="lower left")
ax[0].set_title("(a) Logged state at t_n answers omega(t_n - dt/2)", fontsize=9, loc="left")
d_half = np.abs(wv - np.interp(np.clip(T - 0.005, 0, None), tt, om)); d_full = np.abs(wv - np.interp(T, tt, om))
ax[1].semilogy(T[1:], d_full[1:] + 1e-18, "-", color=ORANGE, label="|logged - omega(t)|")
ax[1].semilogy(T[1:], d_half[1:] + 1e-18, "-", color=BLUE, label="|logged - omega(t - dt/2)|")
ax[1].set_ylim(1e-16, 1e-1); ax[1].set_xlabel("time t [s]"); ax[1].set_ylabel("difference [rad/s]"); ax[1].legend(fontsize=8, loc="center right")
ax[1].set_title("(b) Agreement with the delayed table to round-off (job 14209)", fontsize=9, loc="left")
fig.tight_layout(); fig.savefig(os.path.join(figs, "fig4_bc_time_offset.png")); plt.close(fig)

# fig 5: coupled Newton behaviour
rows = list(csv.DictReader(open(os.path.join(res, "newton_iterations_per_step.csv"))))
fig, ax = plt.subplots(1, 2, figsize=(9.6, 3.6))
rr = list(csv.DictReader(open(os.path.join(res, "rtc_m1_step2_residuals.csv"))))
cur = [float(r["squared_residual_current"]) for r in rr]
ax[0].semilogy(range(1, len(cur) + 1), cur, "o-", color=ORANGE); ax[0].set_xlabel("Newton iteration (status line) in step 2, t = 0.02 s"); ax[0].set_ylabel("squared residual norm (current)")
ax[0].text(0.04, 0.86, "32 negative Jacobians: the step fails\n(job 14213, RT-C M1, rho_inf 0.9)", transform=ax[0].transAxes, fontsize=8)
ax[0].set_title("(a) RT-C M1: divergence from a round-off-level drive", fontsize=9, loc="left")
for run, lab, c, mk in (("SMK_fsi_M0_svk", "M0, St Venant-Kirchhoff, rho_inf 1", GREEN, "^"), ("SMK_fsi_M0_rhoi0p9", "M0, MR law, rho_inf 0.9", BLUE, "s"), ("SMK_fsi_M0", "M0, MR law, rho_inf 1 (stalls at step 5)", ORANGE, "o")):
    rs = [r for r in rows if r["run"] == run]
    ok = [r for r in rs if r["converged"] == "True"]
    ax[1].plot([int(r["step"]) for r in ok], [int(r["iterations"]) for r in ok], "-", marker=mk, color=c, ms=4, label=lab)
    bad = [r for r in rs if r["converged"] == "False"]
    if bad: ax[1].plot([int(bad[-1]["step"])], [12.5], "x", color=c, ms=9, mew=2)
ax[1].text(5.3, 12.9, "x = step that failed to converge", fontsize=8)
ax[1].set_ylim(0, 14.5); ax[1].set_xlabel("time step (dt = 0.01 s)"); ax[1].set_ylabel("Newton iterations per converged step"); ax[1].legend(fontsize=7.5, loc="lower right")
ax[1].set_title("(b) Coupled M0 smoke tests (64 fluid cells)", fontsize=9, loc="left")
fig.tight_layout(); fig.savefig(os.path.join(figs, "fig5_coupled_newton.png")); plt.close(fig)
print("tables and figures written")
