#!/usr/bin/env python3
"""Cavity figures: mesh convergence (u_y, F_y vs dx with solids4foam/other references), time sensitivity, cross-solver, cost. Uses fsi_check.json of available runs."""
import os, json, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
R = "runs/cavityNext"; L = lambda n: json.load(open(f"{R}/{n}/fsi_check.json")) if os.path.exists(f"{R}/{n}/fsi_check.json") else None
feb = {0.1: L("C2_fsi_pilot_f1"), 0.05: L("C3_fsi_M2"), 0.025: L("C3_fsi_M4")}; feb = {k: v for k, v in feb.items() if v}
ref = {  # solids4foam top-surface u_y (manuscript), QM F_y (paper-data quasi_monolithic_618.csv), Tukovic/Ansys/COMSOL (PR #22 reference_levels.csv)
 "solids4foam partitioned (Aitken)": ([0.1, 0.05, 0.025, 0.0125], [-0.2046, -0.2314, -0.2384, -0.2401], None),
 "solids4foam quasi-monolithic": ([0.1, 0.05, 0.025, 0.0125], [-0.203601, -0.231044, -0.2382, -0.239986], [-5.122626, -5.226559, -5.250356, -5.255992]),
 "Tukovic et al. 2018 (constitutive law unresolved)": ([0.1, 0.05, 0.025, 0.0125], [-0.23324, -0.24571, -0.24955, -0.25062], [-5.1603, -5.229, -5.2494, -5.2548]),
 "ANSYS (constitutive law unresolved)": ([0.1, 0.05, 0.025, 0.0125], [-0.24727, -0.25054, -0.25197, -0.25274], [-5.144, -5.218, -5.25, -5.263]),
 "COMSOL (solid model unconfirmed)": ([0.1, 0.05, 0.025, 0.0125], [-0.236, -0.241, -0.243, -0.243], [-5.1088, -5.2276, -5.2956, -5.2952])}
sty = {"solids4foam partitioned (Aitken)": ("C0", "o", "-"), "solids4foam quasi-monolithic": ("C1", "s", "-"), "Tukovic et al. 2018 (constitutive law unresolved)": ("0.5", "^", ":"), "ANSYS (constitutive law unresolved)": ("0.3", "v", ":"), "COMSOL (solid model unconfirmed)": ("0.65", "D", ":")}
fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.2))
for k, (x, uy, fy) in ref.items():
    c, m, ls = sty[k]; ax[0].plot(x, uy, marker=m, ls=ls, color=c, ms=4, lw=1, label=k)
    if fy: ax[1].plot(x, fy, marker=m, ls=ls, color=c, ms=4, lw=1, label=k)
if feb:
    dx = sorted(feb, reverse=True); ax[0].plot(dx, [feb[k]["uy_A_final"] for k in dx], "o-", color="C3", lw=2, ms=6, label="FEBio hex8 fluid-FSI (this work)")
    ax[1].plot(dx, [feb[k]["interface"]["Fy_pressure_raw"] for k in dx], "o-", color="C3", lw=2, ms=6, label="FEBio (raw pressure force)")
    for k in dx: ax[1].fill_between([k * 0.93, k * 1.07], feb[k]["interface"]["Fy_odd_only"], feb[k]["interface"]["Fy_even_only"], color="C3", alpha=0.25, lw=0)
for a_, yl in zip(ax, ("u_y(A) (m)", "F_y (N/m)")): a_.set_xscale("log"); a_.set_xlabel("dx (m)"); a_.set_ylabel(yl); a_.invert_xaxis(); a_.grid(alpha=0.25, which="both")
h, l = ax[0].get_legend_handles_labels(); fig.legend(h, l, loc="lower center", ncol=3, fontsize=7, frameon=False); ax[1].set_title("(b) steady vertical interface force per unit depth (shaded: FEBio odd-even bracket)", fontsize=8)
ax[0].set_title("(a) steady top-surface displacement at A", fontsize=9)
plt.tight_layout(rect=(0, 0.12, 1, 1)); plt.savefig("figures/cavity_mesh_convergence_cross_solver.png", dpi=130); plt.close()
# time/numerical sensitivity (M1): relative change against the baseline at t = 900 s (C2 history, identical trajectory)
base_uy, base_fy = -0.21195916, -4.717604
sens = [("tight d/e/r tol", "C4_tol_M1"), ("K 1e7", "C4_K1e7_M1"), ("dt 1.25", "C4_dt1p25_M1")]
fig, ax = plt.subplots(1, 2, figsize=(10, 3.4))
vals = [(n, L(r)) for n, r in sens if L(r)]
for j, (key, base, ttl) in enumerate((("uy", base_uy, "u_y(A)"), ("Fy", base_fy, "F_y (raw pressure force)"))):
    ch = [((d["uy_A_final"] if key == "uy" else d["interface"]["Fy_pressure_raw"]) / base - 1) * 1e6 for n, d in vals]
    ax[j].bar(range(len(vals)), ch, color="C3", width=0.5); ax[j].axhline(100, color="0.5", ls="--", lw=0.8); ax[j].axhline(-100, color="0.5", ls="--", lw=0.8)
    ax[j].set_xticks(range(len(vals))); ax[j].set_xticklabels([n for n, _ in vals], fontsize=8); ax[j].set_ylabel("relative change vs baseline (1e-6)"); ax[j].set_title(f"{ttl}, M1 steady, one change at a time", fontsize=9); ax[j].grid(alpha=0.25, axis="y")
    for i, c in enumerate(ch): ax[j].text(i, c, f"{c:+.1f}", ha="center", va="bottom" if c >= 0 else "top", fontsize=8)
ax[0].text(0.02, 0.95, "dashed: +/-1e-4 (steady-state window spread)", transform=ax[0].transAxes, fontsize=7, va="top")
plt.tight_layout(); plt.savefig("figures/cavity_numerical_sensitivity.png", dpi=130); plt.close()
# cost
rows = []
for n, lab in [("C1_fluid_fixedplate_f1", "C1 fluid-only M1"), ("C1b_fluid_fixedplate_f1_pin", "C1b fluid-only M1"), ("C2_fsi_pilot_f1", "FSI M1 (C2)"), ("C3_fsi_M2", "FSI M2"), ("C3_fsi_M4", "FSI M3"), ("C4_dt1p25_M1", "FSI M1 dt/2")]:
    p = f"{R}/{n}/time.txt"; m = f"{R}/{n}/run_meta_start.json"
    if os.path.exists(p) and os.path.exists(m):
        t = open(p).read(); import re; w = re.search(r"Elapsed \(wall clock\) time.*?: ([\d:.]+)", t)
        if w and os.path.exists(f"{R}/{n}/run_meta_end.json"):
            parts = [float(x) for x in w.group(1).split(":")]; sec = parts[-1] + 60 * parts[-2] + (3600 * parts[-3] if len(parts) == 3 else 0); cores = json.load(open(m))["cpus_per_task"]; rows.append((lab, sec / 3600 * cores, sec / 60))
fig, a_ = plt.subplots(figsize=(6.5, 3.4)); a_.barh([r[0] for r in rows], [r[1] for r in rows], color="C3"); a_.set_xlabel("allocated core-hours"); a_.invert_yaxis()
for i, r in enumerate(rows): a_.text(r[1], i, f"  {r[2]:.0f} min", va="center", fontsize=7)
a_.set_title("cavity runs: allocated core-hours (wall minutes labelled)", fontsize=9); plt.tight_layout(); plt.savefig("figures/cavity_cost_summary.png", dpi=130)
