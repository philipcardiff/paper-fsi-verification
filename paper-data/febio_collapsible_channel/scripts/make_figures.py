#!/usr/bin/env python3
"""Regenerate the figures from the stored result tables (no simulation).  usage: make_figures.py RESULTS_DIR REFERENCE_DIR OUT_DIR  (RESULTS_DIR = results/ written by make_tables.py)"""
import sys, os, csv, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
RES, REF, OUT = sys.argv[1:4]; os.makedirs(OUT, exist_ok=True)
BL, VM, GR, BK = "#0072B2", "#D55E00", "#009E73", "#444444"
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.alpha": 0.25, "figure.dpi": 130})
def rd(n):
    r = list(csv.reader(open(f"{RES}/{n}"))); return r[0], np.array([[float(v) if v not in ("", "n/a") else np.nan for v in row] for row in r[1:]])
# fig1: S1 convergence
sp = {"80x8": 10 / 80, "160x4": 10 / 160, "160x8": 10 / 160, "320x8": 10 / 320}
rows = list(csv.reader(open(f"{RES}/s1_static_wall.csv")))[1:]
fig, ax = plt.subplots(figsize=(6.6, 4.4))
off = {"80x8": (-30, 10), "160x4": (10, 10), "160x8": (10, -16), "320x8": (10, 8)}
for r in rows:
    e = abs(float(r[3])); sq = r[0] == "160x4"
    ax.loglog(sp[r[0]], e, "s" if sq else "o", color=VM if sq else BL, ms=10 if sq else 6, mfc="none" if sq else BL, mew=1.6)
    ax.annotate(f"{r[0]}: {float(r[3]):+.2f} %", (sp[r[0]], e), textcoords="offset points", xytext=off[r[0]], fontsize=9)
xx = np.array([0.03, 0.14]); ax.loglog(xx, 0.757 * (xx / (10 / 320)) ** 2, "--", color=GR, lw=1.2, label="slope 2 (second order)")
ax.axhline(3.0, color=BK, lw=0.9, ls=":"); ax.text(0.0305, 3.25, "pre-registered 3 % gate", fontsize=8.5, color=BK); ax.axhline(1.0, color=BK, lw=0.9, ls=":"); ax.text(0.095, 1.08, "1 % gate (320x8)", fontsize=8.5, color=BK)
ax.legend(frameon=False, fontsize=8.5, loc="lower right")
ax.set_xlabel("in-plane wall element size dx (m)"); ax.set_ylabel("|midpoint error| vs beam value -0.143628 m (%)"); ax.set_title("CC-S1: static clamped wall, stock hex8 (open square: 160x4, 4 rows)", fontsize=10); ax.set_ylim(0.5, 20)
fig.tight_layout(); fig.savefig(f"{OUT}/fig1_s1_static_wall_convergence.png"); plt.close(fig)
# fig2: F1 wall rows
h, a = rd("f1_wall_rows_final.csv"); x, p0, p1, p2 = a[:, 0], a[:, 2], a[:, 3], a[:, 4]
fig, axs = plt.subplots(1, 2, figsize=(10, 4)); ax = axs[0]
ax.plot(x, p0, color=BL, lw=1.6, label="wall row (applied nodal pressure)"); ax.plot(x, p1, color=VM, lw=1.6, label="row 1 (first row off the wall)"); ax.plot(x, p2, color=GR, lw=1.0, ls="--", label="row 2 (coincides with the wall row to 4e-3 Pa)")
ax.set_xlabel("x (m)"); ax.set_ylabel("nodal pressure (Pa)"); ax.set_title("CC-F1 fixed shape, steady: pressure rows next to the wall", fontsize=10); ax.legend(frameon=False, fontsize=8.5)
ax = axs[1]; ax.plot(x, p0 - p1, color=VM, lw=1.4, label="p(wall) - p(row 1)"); ax.plot(x, p0 - 2 * p1 + p2, color=BK, lw=1.4, label="D = p0 - 2 p1 + p2"); ax.plot(x, (p0 - p2) * 100, color=GR, lw=1, label="100 x (p(wall) - p(row 2))")
ax.axhline(0, color=BK, lw=0.6); ax.set_xlabel("x (m)"); ax.set_ylabel("Pa"); ax.set_title("Wall-normal row parity (smooth in x)", fontsize=10); ax.legend(frameon=False, fontsize=8.5)
fig.tight_layout(); fig.savefig(f"{OUT}/fig2_f1_wall_normal_row_parity.png"); plt.close(fig)
# fig3: D1 loads and spread
h, a = rd("d1_load_spread_by_step.csv"); t = a[:, 0]; F1, F2, sF, D = a[:, 1], a[:, 2], a[:, 3], a[:, 8]
fig, axs = plt.subplots(1, 2, figsize=(10, 4.6)); ax = axs[0]
ax.plot(t, F1, "-o", ms=3, color=BL, label="F1: nodal wall pressure (what FEBio applies)"); ax.plot(t, F2, "-o", ms=3, color=VM, label="F2: element-mean pressure (reconstruction)"); ax.axhline(0, color=BK, lw=0.6)
ax.set_xlabel("t (s)"); ax.set_ylabel("net vertical wall force F_y (N/m)"); ax.set_title("D1 coupled diagnostic, 20 steps", fontsize=10); ax.legend(frameon=False, fontsize=8.5)
ax = axs[1]; m = (np.abs(F2) >= 100) & (t >= 0.10)
ax.plot(t[m], 100 * sF[m], "-o", ms=3, color=BL, label="(F1 - F2)/F2 (%), steps with |F2| >= 100 N/m"); ax.axhline(3, color=BK, ls=":", lw=1, label="3 % criterion of hypothesis H_D")
ax2 = ax.twinx(); ax2.plot(t, D, "-s", ms=3, color=GR, lw=1, label="mean y-parity D (right axis; all steps)"); ax2.set_ylabel("mean y-parity D = p0 - 2 p1 + p2 (Pa)", color=BK); ax2.grid(False); ax2.spines["right"].set_visible(True)
ax.set_xlabel("t (s)"); ax.set_ylabel("relative spread (%)"); ax.set_title("Spread and y-parity (green, right axis)", fontsize=10); hh, ll = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels(); ax.legend(hh + h2, ll + l2, frameon=False, fontsize=8.5, loc="upper center", bbox_to_anchor=(0.5, -0.17))
fig.tight_layout(); fig.savefig(f"{OUT}/fig3_d1_applied_load_spread.png"); plt.close(fig)
# fig4: D1 wall midpoint vs oomph
mid, om = a[:, 12], a[:, 13]
def load(n):
    q = np.array([l.split(",") for l in open(f"{REF}/collapsibleChannel_oomph_{n}.csv") if l[0] not in "#t"], float); return q[:, 0], q[:, 2]
tr, ur = load("reference"); t4, u4 = load("dt40")
fig, axs = plt.subplots(2, 1, figsize=(6.6, 5.6), sharex=True, gridspec_kw={"height_ratios": [3, 1.3]}); ax = axs[0]
ax.plot(tr[tr <= 0.55], ur[tr <= 0.55], color=BK, lw=1.5, label="oomph-lib converged reference"); ax.plot(t4[t4 <= 0.55], u4[t4 <= 0.55], "o", ms=3, color=GR, label="oomph-lib, dt = 0.025 s")
ax.plot(t, mid, "-s", ms=3.5, color=VM, label="FEBio v4.13, 20-step diagnostic (not a trough: still descending)"); ax.set_ylabel("wall midpoint u_y (m)"); ax.legend(frameon=False, fontsize=8.5, loc="lower left"); ax.set_title("D1 vs oomph-lib: informational, NOT a verification", fontsize=10)
axs[1].plot(t, 100 * a[:, 14], "-s", ms=3.5, color=VM); axs[1].axhline(0, color=BK, lw=0.6); axs[1].set_xlabel("t (s)"); axs[1].set_ylabel("FEBio - oomph\n(% of 0.216 m)")
fig.tight_layout(); fig.savefig(f"{OUT}/fig4_d1_wall_midpoint_vs_oomph.png"); plt.close(fig)
# fig5: free body
h, a = rd("d1_freebody_identification.csv"); t = a[:, 0]
fig, ax = plt.subplots(figsize=(6.2, 4.2))
ax.plot(t, -a[:, 1], "-o", ms=5, color=BK, label="identified from wall balance (free body)"); ax.plot(t, a[:, 2], "-s", ms=4, color=BL, label="F1 equivalent (nodal wall pressure)"); ax.plot(t, a[:, 3], "-^", ms=4, color=VM, label="F2 equivalent (element mean)"); ax.plot(t, a[:, 4], ":", color=GR, label="odd-row (row 1) equivalent")
ax.set_xlabel("t (s)"); ax.set_ylabel("uniform-equivalent fluid pressure on the wall (Pa)"); ax.set_title("Which load moves the wall? (valid only for |u_mid| <= 0.1425 m)", fontsize=10); ax.legend(frameon=False, fontsize=8.5)
fig.tight_layout(); fig.savefig(f"{OUT}/fig5_d1_freebody_identification.png"); plt.close(fig)
