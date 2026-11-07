#!/usr/bin/env python3
"""Regenerate the package figures from the retained data (results/*.csv, results/raw/*).  Package-only (numpy, matplotlib).  Usage from the package root: python3 -I scripts/make_figures.py
Figures: fig_mesh_convergence_vs_references, fig_force_definitions, fig_row_pressure_alternation, fig_numerical_sensitivity, fig_stabilisation_outcomes.
(The original figures in figures/cavity_*.png were produced from HPC-resident outputs by scripts/orig_cavity_figures_*.py and are not regenerated here.)"""
import csv, json, os
import numpy as np, matplotlib, matplotlib.ticker; matplotlib.use("Agg"); import matplotlib.pyplot as plt
C = {"febio": "#D55E00", "f2": "#CC79A7", "pt": "#0072B2", "qm": "#009E73", "tuk": "#7F7F7F", "ans": "#444444", "com": "#999999"}
plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.alpha": 0.25, "axes.spines.top": False, "axes.spines.right": False})
rd = lambda p: list(csv.DictReader(open(p))); f = float
ref = rd("results/reference_values.csv"); cs = rd("results/coupled_febio_summary.csv")[:3]; fr = {r["run_state"]: r for r in rd("results/force_reconstruction_F1_F4.csv")}
def series(src, var, key):
    rows = sorted([r for r in ref if r["source"] == src and var in r["variant"] and r[key] != ""], key=lambda r: -f(r["dx_m"])); return [f(r["dx_m"]) for r in rows], [f(r[key]) for r in rows]
dx = [f(r["dx_m"]) for r in cs]
NOTE = "FEBio: all three runs error-terminated (dilatation-increment floor); M3 not steady to 1e-4: conditional estimates, not a verification"
# 1 mesh convergence vs references
fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.6))
for a, key, yl in ((ax[0], "uy_A_top_m", "u_y(A), top surface (m)"), (ax[1], "Fy_Npm", "F_y (N/m per m depth)")):
    for src, var, col, mk, lab in (("solids4foam", "partitioned", C["pt"], "o", "solids4foam partitioned (StVK)"), ("solids4foam", "#618", C["qm"], "s", "solids4foam quasi-monolithic #618 (StVK)"),
                                   ("Tukovic et al. (2018)", "", C["tuk"], "^", "Tukovic et al. 2018 (law unresolved)"), ("ANSYS", "", C["ans"], "v", "ANSYS (law unconfirmed)"), ("COMSOL", "", C["com"], "D", "COMSOL (model unconfirmed)")):
        x, y = series(src, var, key)
        if x: a.plot(x, y, marker=mk, ms=4, lw=1, color=col, label=lab, ls="-" if "solids4foam" in src else ":")
    a.plot(dx, [f(r["uy_A_m" if key.startswith("uy") else "Fy_pressure_raw_Npm"]) for r in cs], "o-", color=C["febio"], lw=2, ms=6, label="FEBio v4.13, raw (F1: applied wall-row load)")
    if key == "Fy_Npm":
        a.plot(dx, [f(fr[k]["F2_Fy"]) for k in ("M1", "M2", "M3")], "d--", color=C["f2"], lw=1.2, ms=6, mfc="none", label="FEBio element-mean F2 (diagnostic; not shown to be the physical force); states t=900/750/500")
        for r, d in zip(cs, dx): a.fill_between([d * 0.93, d * 1.07], f(r["Fy_odd_only"]), f(r["Fy_even_only"]), color=C["febio"], alpha=0.2, lw=0)
    a.set_xscale("log"); a.invert_xaxis(); a.set_xlabel("dx (m)"); a.set_ylabel(yl); a.set_xticks([0.1, 0.05, 0.025, 0.0125]); a.set_xticklabels(["0.1", "0.05", "0.025", "0.0125"]); a.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax[0].set_title("(a) steady u_y at A", fontsize=9); ax[1].set_title("(b) steady F_y (shaded: FEBio odd-even bracket on F1)", fontsize=9)
h, l = ax[1].get_legend_handles_labels(); fig.legend(h, l, loc="lower center", ncol=2, fontsize=7.5, frameon=False); fig.suptitle(NOTE, fontsize=8, color="#7a0000")
plt.tight_layout(rect=(0, 0.17, 1, 0.95)); plt.savefig("figures/fig_mesh_convergence_vs_references.png", dpi=140); plt.close()
# 2 force definitions
qm = {r["dx_m"]: f(r["Fy_Npm"]) for r in ref if r["source"] == "solids4foam" and "#618" in r["variant"]}
fig, ax = plt.subplots(figsize=(7.2, 4.4)); keys = [("M1", "0.1", "t=900"), ("M2", "0.05", "t=750"), ("M3", "0.025", "t=500")]; x0 = np.arange(3)
for j, (k, lab, col, mk) in enumerate((("F1_Fy", "F1 applied (nodal wall-row p; moves the plate)", C["febio"], "o"), ("F2_Fy", "F2 element-mean = F4 (diagnostic)", C["f2"], "d"), (None, "solids4foam QM #618 (reference)", C["qm"], "s"))):
    v = [qm[d] if k is None else f(fr[m][k]) for m, d, _ in keys]; ax.plot(x0 + (j - 1) * 0.12, v, mk, color=col, ms=9, label=lab)
    for xx, vv in zip(x0 + (j - 1) * 0.12, v): ax.text(xx + 0.05, vv, f"{vv:.3f}", fontsize=7, va="center")
ax.set_xticks(x0); ax.set_xticklabels([f"{m} (dx {d})\nstate {t} s" for m, d, t in keys]); ax.set_ylabel("F_y, pressure contribution (N/m)"); ax.set_xlim(-0.4, 2.6)
ax.set_title("Pressure force on the plate, three definitions (stored states, same time per level)\nF2 vs reference: +6.1 % (M1), +0.9 % (M2), -0.07 % (M3, t=500): not a general identity", fontsize=8.5)
ax.legend(fontsize=7.5, loc="lower right", frameon=False); plt.tight_layout(); plt.savefig("figures/fig_force_definitions.png", dpi=140); plt.close()
# 3 row alternation
fig, ax = plt.subplots(1, 2, figsize=(10.5, 3.9))
for key, lab, col in (("M1", "M1 t=900", "#0072B2"), ("M2", "M2 t=750", "#009E73"), ("M3", "M3 t=500", C["febio"]), ("C1b_fixed_M1", "C1b fixed plate", "#7F7F7F")):
    d = json.load(open(f"results/raw/b2/ybranch_{key}.json")); ax[0].plot(range(6), d["rows_mean_interior_first6"], "o-", color=col, label=lab, ms=4)
ax[0].set_xlabel("node row above the floor (0 = interface row)"); ax[0].set_ylabel("mean nodal p over interior floor columns (Pa)"); ax[0].set_title("(a) alternation does not decay through six rows", fontsize=9); ax[0].legend(fontsize=7.5, frameon=False)
d = json.load(open("results/raw/b2/ybranch_M3_t717p5.json"))["per_node"]; x = np.array(d["x_ref"]); p0 = np.array(d["p_row0"]); p1 = np.array(d["p_row1"])
ax[1].plot(x, p0, "-", color=C["febio"], lw=1.2, label="row 0 (wall row, F1)"); ax[1].plot(x, p1, "-", color="#0072B2", lw=1.2, label="row 1"); ax[1].plot(x, 0.5 * (p0 + p1), "--", color=C["f2"], lw=1.2, label="pair mean (diagnostic only)")
ax[1].set_xlabel("x (m)"); ax[1].set_ylabel("p (Pa)"); ax[1].set_title("(b) M3, failure-exit state t=717.5 s: floor pressure along x", fontsize=9); ax[1].legend(fontsize=7.5, frameon=False)
plt.tight_layout(); plt.savefig("figures/fig_row_pressure_alternation.png", dpi=140); plt.close()
# 4 sensitivity
sens = rd("results/numerical_sensitivity_C4.csv")[1:]; fig, ax = plt.subplots(1, 2, figsize=(9, 3.4))
for a, key, ttl in ((ax[0], "rel_change_uy", "u_y(A)"), (ax[1], "rel_change_Fy", "F_y raw (F1)")):
    v = [f(r[key]) * 1e6 for r in sens]; a.bar(range(3), v, color=C["febio"], width=0.5); a.set_xticks(range(3)); a.set_xticklabels([r["change"] for r in sens], fontsize=7)
    a.set_ylabel("relative change vs baseline (1e-6)"); a.set_title(f"{ttl}: M1, one change at a time", fontsize=9); a.axhline(100, color="0.5", ls="--", lw=0.8); a.axhline(-100, color="0.5", ls="--", lw=0.8)
    for i, c in enumerate(v): a.text(i, c, f"{c:+.1f}", ha="center", va="bottom" if c >= 0 else "top", fontsize=8)
ax[0].text(0.02, 0.95, "dashed: +/-1e-4 (window spread)", transform=ax[0].transAxes, fontsize=7, va="top"); plt.tight_layout(); plt.savefig("figures/fig_numerical_sensitivity.png", dpi=140); plt.close()
# 5 stabilisation outcomes
r1 = rd("results/ppp_round1_C1b.csv"); r2 = rd("results/ppp_round2_T5_poiseuille.csv"); fig, ax = plt.subplots(1, 3, figsize=(12, 3.8))
ax[0].bar([0, 1], [abs(f(r1[0]["D_row1_minus_row0_x2p4_to_5p6_Pa"])), abs(f(r1[2]["D_row1_minus_row0_x2p4_to_5p6_Pa"]))], color=[C["febio"], "#0072B2"], width=0.5); ax[0].set_yscale("log"); ax[0].set_xticks([0, 1]); ax[0].set_xticklabels(["stock", "round 1\nalpha=3/4"])
ax[0].set_ylabel("|D| row parity (Pa), C1b M1"); ax[0].set_title("(a) mode suppressed 4500x", fontsize=9)
q = [abs(f(r1[0]["Qout_over_Qin_minus_1"])), abs(f(r1[2]["Qout_over_Qin_minus_1"]))]; ax[1].bar([0, 1], q, color=[C["febio"], "#0072B2"], width=0.5); ax[1].set_yscale("log"); ax[1].axhline(1e-3, color="k", ls="--", lw=0.8); ax[1].text(0.5, 1.3e-3, "gate 1e-3", fontsize=7, ha="center")
ax[1].set_xticks([0, 1]); ax[1].set_xticklabels(["stock", "round 1 alpha=3/4"]); ax[1].set_ylabel("|Qout/Qin - 1|, C1b M1"); ax[1].set_title("(b) round 1 FAILED mass conservation (-1.55 %)", fontsize=9)
v = [f(r["max_abs_flux_err_all_sections"]) for r in r2]; ax[2].bar(range(3), v, color=[C["febio"], "#999999", "#0072B2"], width=0.5); ax[2].set_yscale("log"); ax[2].axhline(1e-4, color="k", ls="--", lw=0.8); ax[2].text(1, 1.25e-4, "criterion 1e-4", fontsize=7, ha="center")
ax[2].set_xticks(range(3)); ax[2].set_xticklabels(["stock", "alpha=0\nrim once", "alpha=1e-3\nround 2"], fontsize=8); ax[2].set_ylabel("max |Q/Qin - 1|, tiny Poiseuille channel"); ax[2].set_title("(c) round 2 FAILED nodal-flux criterion (4.9e-4)", fontsize=9)
plt.tight_layout(); plt.savefig("figures/fig_stabilisation_outcomes.png", dpi=140); plt.close(); print("figures written")
