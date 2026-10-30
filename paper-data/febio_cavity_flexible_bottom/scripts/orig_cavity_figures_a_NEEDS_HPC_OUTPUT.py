#!/usr/bin/env python3
"""Cavity figures: (1) geometry + BCs, (2) representative result (C1b velocity magnitude, interface pressure with odd-even, C2 deflection)."""
import sys, os, json, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
sys.path.insert(0, os.path.dirname(__file__)); import make_cavity_fluid_feb as G, cavity_traction as ct
nodes2, quads = G.mesh(1)
fig, ax = plt.subplots(2, 1, figsize=(11, 6.2), gridspec_kw={"height_ratios": [1.0, 1.25]})
a = ax[0]; polys = [[nodes2[i] for i in q] for q in quads]
a.add_collection(PolyCollection(polys, facecolors="#d7e6f5", edgecolors="#8aa9c6", linewidths=0.25)); a.add_patch(plt.Rectangle((2, -0.4), 4, 0.4, fc="#e8c9a0", ec="#7a5a2a", lw=1.0))
a.set_xlim(-0.3, 14.3); a.set_ylim(-0.7, 2.35); a.set_aspect("equal"); a.set_xlabel("x (m)"); a.set_ylabel("y (m)", labelpad=2)
a.annotate("inlet: parabolic, mean 1 m/s\n(built-in fluid normal velocity)", (0, 1.5), (0.4, 2.12), fontsize=7, arrowprops=dict(arrowstyle="->", lw=0.6))
a.annotate("outlet: p = 0", (14, 1.5), (12.4, 2.12), fontsize=7, arrowprops=dict(arrowstyle="->", lw=0.6)); a.text(7, 2.12, "no-slip walls (top, cavity sides)", fontsize=7, ha="center")
a.text(4, -0.2, "StVK plate  E = 500 Pa, nu = 0.4, rho = 1000\n4 m x 0.4 m, clamped at both ends", fontsize=7, ha="center", va="center")
a.text(4, 0.4, "fluid-solid interface y = 0 (A at x = 4)", fontsize=7, ha="center"); a.plot([4], [0], "ko", ms=3); a.set_title("(a) cavityFlexibleBottom: fluid rho = 1, nu = 0.01 (Re = 100), plane strain; mesh = level 1 (dx = 0.1, 1688 cells)", fontsize=9)
b = ax[1]
rd = "runs/cavityNext/C1b_fluid_fixedplate_f1_pin"; meta, nd, mon = ct.load(rd); t = max(mon); z0 = {i: nd[i] for i in mon[t] if abs(nd[i][2]) < 1e-9}
ids = sorted(z0); X = np.array([z0[i][0] for i in ids]); Y = np.array([z0[i][1] for i in ids]); U = np.array([np.hypot(mon[t][i][0], mon[t][i][1]) for i in ids])
idmap = {}
for i in ids: idmap[(round(z0[i][0], 6), round(z0[i][1], 6))] = i
xy2i = lambda P: idmap[(round(P[0], 6), round(P[1], 6))]
vals = [np.mean([np.hypot(mon[t][xy2i(nodes2[k])][0], mon[t][xy2i(nodes2[k])][1]) for k in q]) for q in quads]
pc = PolyCollection(polys, array=np.array(vals), cmap="viridis", edgecolors="face", linewidths=0.2); pc.set_clim(0, 1.6); b.add_collection(pc)
plt.colorbar(pc, ax=b, pad=0.01, shrink=0.9, label="|w| (m/s), cell mean"); b.set_xlim(0, 14); b.set_ylim(-0.05, 2.05); b.set_aspect("equal")
b.set_title("(b) fixed-plate fluid-only (C1b): steady velocity magnitude", fontsize=9); b.set_xlabel("x (m)")
plt.subplots_adjust(left=0.07, right=0.97, top=0.95, bottom=0.07, hspace=0.25); plt.savefig("figures/cavity_geometry_and_flow.png", dpi=130); plt.close()
# representative interface result
fig, ax = plt.subplots(1, 2, figsize=(11, 3.8))
wc = ct.wall_columns(nd, mon, 0.0, 2.0, 6.0); q = ct.quantities(meta, nd, mon, wc, t); xs, p = q["xs"], q["p"]
ax[0].plot(xs, p, "o-", ms=3, lw=0.8, color="C0", label="raw nodal p (fixed plate, C1b)"); ax[0].plot(xs[::2], p[::2], "s", ms=4, mfc="none", color="C1", label="even-index nodes"); ax[0].plot(xs[1::2], p[1::2], "^", ms=4, mfc="none", color="C2", label="odd-index nodes")
ax[0].set_xlabel("x (m)"); ax[0].set_ylabel("p on the cavity floor (Pa)"); ax[0].legend(fontsize=7); ax[0].set_title("(a) interface pressure: odd-even mode (corner nodes pinned to 1.15 Pa)", fontsize=9)
if os.path.exists("runs/cavityNext/C2_fsi_pilot_f1/fsi_check.json"):
    d = json.load(open("runs/cavityNext/C2_fsi_pilot_f1/fsi_check.json")); y = np.array(d["interface"]["deflection_profile_uy"]); ax[1].plot(np.linspace(2, 6, len(y)), y, "o-", color="C3", label="FEBio M1, steady (every 4th node)")
    ax[1].axhline(-0.2407, color="0.4", ls="--", lw=0.8); ax[1].text(2.1, -0.2407 + 0.004, "solids4foam extrapolated u_y(A) = -0.2407 (reference, centre value)", fontsize=7)
    ax[1].set_xlabel("x (m)"); ax[1].set_ylabel("u_y of the plate top surface (m)"); ax[1].legend(fontsize=7); ax[1].set_title("(b) coupled pilot C2: plate top-surface deflection", fontsize=9)
plt.tight_layout(); plt.savefig("figures/cavity_interface_pressure_and_deflection.png", dpi=130)
