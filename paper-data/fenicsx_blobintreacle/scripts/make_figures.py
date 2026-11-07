"""Figures S1 and S2 of the FEniCSx blobInTreacle supplement, from data/*.csv.
Usage: python3 scripts/make_figures.py   (needs numpy, matplotlib; writes figures/*.pdf and *.png)"""
import csv, os
from collections import defaultdict
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({"font.size": 8, "mathtext.fontset": "cm", "font.family": "serif", "pdf.fonttype": 42})
here = os.path.dirname(os.path.abspath(__file__))
data, figs = os.path.join(here, "..", "data"), os.path.join(here, "..", "figures")
rd = lambda n: list(csv.DictReader(open(os.path.join(data, n))))
# Okabe-Ito colours; marker and line style carry identity as well
sty = {"fenicsx_M2_BE": ("#0072B2", "o", "-", "FEniCSx M2 (BE)"), "fenicsx_M1_BE": ("#D55E00", "s", "-", "FEniCSx M1 (BE)"),
       "s4f_M2_original": ("#009E73", "^", "--", "solids4foam M2 (original)"), "s4f_EAREA_diagnostic": ("#CC79A7", "v", ":", "solids4foam E-AREA (diagnostic)")}
S = defaultdict(list)
for r in rd("steady_levels.csv"):
    S[(r["model"], r["observable"])].append((int(r["cells_total"]), float(r["value"])))
lim = {(r["model"], r["observable"]): float(r["richardson_limit"]) for r in rd("steady_convergence.csv") if r["triple"] == "finest3" and r["richardson_limit"]}

fig, ax = plt.subplots(1, 2, figsize=(7.0, 2.8))
for a, obs, lab in zip(ax, ["apex_ux", "apex_uy"], ["apex $u_x$ (m)", "apex $u_y$ (m)"]):
    for m, (c, mk, ls, name) in sty.items():
        s = sorted(S[(m, obs)]); h = [n ** -0.5 for n, _ in s]; v = [x for _, x in s]
        a.plot(h, v, color=c, marker=mk, ls=ls, lw=1, ms=4, label=name)
        if (m, obs) in lim:
            a.plot([0], [lim[(m, obs)]], color=c, marker=mk, mfc="none", ms=6, ls="none")
    a.set_xlabel(r"$h=N_\mathrm{cells}^{-1/2}$"); a.set_ylabel(lab); a.set_xlim(-0.005, 0.085); a.grid(alpha=0.25, lw=0.4)
ax[0].legend(frameon=False, fontsize=6.5, loc="lower left")
ax[0].set_title("(a) steady apex $u_x$; open markers: Richardson limits", fontsize=7.5)
ax[1].set_title("(b) steady apex $u_y$", fontsize=7.5)
fig.tight_layout()
for e in ("pdf", "png"): fig.savefig(os.path.join(figs, "fenicsx_blob_spatial." + e), dpi=200)

# temporal errors at t=1 against the independent BDF2 dt=0.00125 reference (F=1, M2)
V = {}
for r in rd("time_values_t1_F1_M2.csv"):
    V[(r["scheme"], float(r["dt"]))] = r
ref = V[("BDF2", 0.00125)]
fig, ax = plt.subplots(1, 2, figsize=(7.0, 2.8))
for a, key, lab in zip(ax, ["apex_ux_t1", "Fy_fluid_on_solid_t1"], ["apex $u_x(t=1)$", "$F_y(t=1)$"]):
    for sch, (c, mk, ls) in {"BE": ("#D55E00", "s", "-"), "BDF2": ("#0072B2", "o", "-"), "CN_consistent": ("#009E73", "^", "--")}.items():
        pts = sorted((dt, abs(float(r[key]) - float(ref[key]))) for (s_, dt), r in V.items() if s_ == sch and dt >= 0.005)
        a.loglog([p[0] for p in pts], [p[1] for p in pts], color=c, marker=mk, ls=ls, lw=1, ms=4, label=sch.replace("_", " "))
    x = np.array([0.005, 0.04]); y0 = abs(float(V[("BDF2", 0.01)][key]) - float(ref[key]))
    a.loglog(x, y0 * (x / 0.01) ** 2, "k:", lw=0.7, label="slope 2"); a.loglog(x, y0 * 8 * (x / 0.01), "k-.", lw=0.7, label="slope 1")
    a.set_xticks([0.005, 0.01, 0.02, 0.04]); a.set_xticklabels(["0.005", "0.01", "0.02", "0.04"]); a.minorticks_off(); a.set_xlabel(r"$\Delta t$ (s)"); a.set_ylabel("|error| vs. BDF2, $\\Delta t=0.00125$"); a.grid(alpha=0.25, lw=0.4, which="both")
ax[0].legend(frameon=False, fontsize=6.5)
ax[0].set_title("(a) apex $u_x$", fontsize=7.5); ax[1].set_title("(b) interface force $F_y$", fontsize=7.5)
fig.tight_layout()
for e in ("pdf", "png"): fig.savefig(os.path.join(figs, "fenicsx_blob_temporal." + e), dpi=200)

# transient spatial study at t=1 (BDF2, dt=0.005): apex ux against h, with solids4foam at the same cell counts
T = rd("bdf2_transient_spatial_t1_M2.csv")
h = [int(r["cells_total"]) ** -0.5 for r in T]
fig, ax = plt.subplots(figsize=(3.6, 2.8))
ax.plot(h, [float(r["apex_ux_t1"]) for r in T], color="#0072B2", marker="o", lw=1, ms=4, label="FEniCSx M2, BDF2, $\\Delta t=0.005$")
ax.plot(h, [float(r["s4f_M2_same_level_ux_t1"]) for r in T], color="#009E73", marker="^", ls="--", lw=1, ms=4, label="solids4foam M2")
ax.plot([0], [0.0576543], color="#0072B2", marker="o", mfc="none", ms=6, ls="none")
ax.plot([0], [0.0577633], color="#009E73", marker="^", mfc="none", ms=6, ls="none")
ax.set_xlabel(r"$h=N_\mathrm{cells}^{-1/2}$"); ax.set_ylabel("apex $u_x(t=1)$ (m)"); ax.grid(alpha=0.25, lw=0.4); ax.legend(frameon=False, fontsize=6.5, loc="lower left")
fig.tight_layout()
for e in ("pdf", "png"): fig.savefig(os.path.join(figs, "fenicsx_blob_transient_spatial." + e), dpi=200)
