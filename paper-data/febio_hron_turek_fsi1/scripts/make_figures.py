#!/usr/bin/env python3
"""Regenerates figures/*.png from the stored results (results/*.csv|json, results/stageC_raw/*/Mon_all.txt, inputs/stageC/*.feb).  usage: make_figures.py PACKAGE_ROOT [OUT_DIR]"""
import sys, os, re, json, csv, math, glob
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
R = os.path.abspath(sys.argv[1]); F = os.path.abspath(sys.argv[2]) if len(sys.argv) > 2 else R + "/figures"; os.makedirs(F, exist_ok=True)
plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.alpha": 0.3, "figure.dpi": 130})
# ---- fig 1: CSM convergence
rows = list(csv.DictReader(open(R + "/results/stageB_summary.csv")))
fig, ax = plt.subplots(1, 2, figsize=(9, 3.6), sharey=True)
for a, k in zip(ax, ("CSM1", "CSM2")):
    for el, mk, ls in (("hex8", "o", "-"), ("hex20", "s", "--")):
        r = sorted([x for x in rows if x["benchmark"] == k and x["element"] == el and "root025" not in x["case"]], key=lambda x: int(x["nt"]))
        a.plot([int(x["nt"]) for x in r], [float(x["uy_err_pct"]) for x in r], mk + ls, label=f"{el}, u_y")
        a.plot([int(x["nt"]) for x in r], [float(x["ux_err_pct"]) for x in r], mk + ls, mfc="none", label=f"{el}, u_x")
    g2 = [x for x in rows if x["benchmark"] == k and "root025" in x["case"]][0]
    a.plot([16], [float(g2["uy_err_pct"])], "kx", ms=8, label="hex8, root x = 0.25 (u_y), informational")
    a.axhline(0, color="k", lw=0.8); a.axhspan(-1, 1, color="g", alpha=0.08, label="+-1 % (u_y criterion)")
    a.set_xscale("log", base=2); a.set_xlabel("elements across the thickness, nt"); a.set_title(f"{k}: error of FEBio v4.13 vs Featflow reference at A")
ax[0].set_ylabel("error at point A (%)"); ax[0].set_ylim(-25, 5); ax[1].legend(fontsize=7, loc="lower right"); fig.tight_layout(); fig.savefig(F + "/fig1_stageB_convergence.png"); plt.close(fig)
# ---- fig 2: shear equilibrium from stored element stress
eq = list(csv.DictReader(open(R + "/results/stageB_equilibrium_from_xplt.csv")))
fig, ax = plt.subplots(figsize=(5.2, 3.6))
for case, mk in (("B_CSM2_h8_nt16", "o-"), ("B_CSM2_h20_nt4", "s--"), ("B_CSM1_h8_nt16", "^-")):
    r = [x for x in eq if x["case"] == case]; ax.plot([float(x["x_section_from_root"]) for x in r], [float(x["shear_ratio"]) for x in r], mk, label=case)
ax.axhline(1, color="k", lw=0.8); ax.set_xlabel("section position from the root, s (m)"); ax.set_ylabel("section shear resultant / (-weight of free part)")
ax.set_title("Free-body shear equilibrium from stored element stresses"); ax.legend(fontsize=7); fig.tight_layout(); fig.savefig(F + "/fig2_stageB_shear_equilibrium.png"); plt.close(fig)
# ---- stage C data
def load(case):
    t = open(f"{R}/inputs/stageC/{case}.feb", encoding="latin-1").read()
    nodes = np.array([[float(v) for v in m.group(2).split(",")] for m in re.finditer(r'<node id="(\d+)">([^<]+)</node>', t)])
    z0 = np.where(np.abs(nodes[:, 2]) < 1e-12)[0]; XY = nodes[z0, :2]
    for ln in open(f"{R}/results/stageC_raw/{case}/Mon_all.txt"):
        if ln.startswith("*Time"): rows_ = []
        elif not ln.startswith("*") and ln.strip(): rows_.append(ln.split())
    K = float(re.search(r"<k>([^<]+)</k>", t).group(1)); return XY, -K * np.array([float(r[3]) for r in rows_])
XY, P = load("C1p_L1_f1.25_K1e6")
ids = sorted([i for i in range(len(XY)) if abs(XY[i, 0] - 1.0) < 1e-9], key=lambda i: XY[i, 1])
fig, ax = plt.subplots(1, 2, figsize=(9, 3.6))
ax[0].plot(XY[ids, 1], P[ids], "o-", ms=3); ax[0].set_xlabel("y (m) along the vertical grid line x = 1.0 m"); ax[0].set_ylabel("nodal pressure p = -K e_f (Pa)"); ax[0].set_title("C1p level 1: node-by-node alternation across the channel")
XR = 0.2 + math.sqrt(0.05**2 - 0.01**2)
for yv, nm, c in ((0.21, "flag top face (y = 0.21)", "C0"), (0.19, "flag bottom face (y = 0.19)", "C3")):
    w = sorted([i for i in range(len(XY)) if abs(XY[i, 1] - yv) < 1e-9 and XR - 1e-9 <= XY[i, 0] <= 0.6 + 1e-9], key=lambda i: XY[i, 0])
    ax[1].plot(XY[w, 0], P[w], ".-", color=c, label=nm)
ax[1].set_xlabel("x (m)"); ax[1].set_ylabel("nodal wall pressure (Pa)"); ax[1].set_title("C1p: wall-node pressure on the two flag faces"); ax[1].legend(fontsize=7)
fig.tight_layout(); fig.savefig(F + "/fig3_stageC_pressure_parity.png"); plt.close(fig)
# ---- fig 4: forces
ind = json.load(open(R + "/results/independent_forces_C1p_L1_f1.25_K1e6.json")); rec = json.load(open(R + "/results/fsi_traction_reconstruction_C1p_L1_f1.25_K1e6.json"))
lab = ["F1", "F-actual\n(recon.)", "F2", "F4", "F3b"]
dr = [ind["F1"][0], rec["F_actual_recon"][0], ind["F2"][0], ind["F4_weak_reaction"][0], ind["F3b"][0]]; li = [ind["F1"][1], rec["F_actual_recon"][1], ind["F2"][1], ind["F4_weak_reaction"][1], ind["F3b"][1]]
fig, ax = plt.subplots(1, 2, figsize=(9.5, 3.8))
for a, v, ref, nm in ((ax[0], dr, 14.2929, "drag F_x (N/m)"), (ax[1], li, 1.11905, "lift F_y (N/m)")):
    a.bar(range(5), v, color=["C3", "C3", "C0", "C2", "C1"]); a.axhline(ref, color="k", ls="--", label=f"Featflow CFD1 level 6+0: {ref}"); a.set_xticks(range(5)); a.set_xticklabels(lab, fontsize=8)
    a.set_ylabel(nm); a.legend(fontsize=7, loc="upper right" if nm.startswith("lift") else "lower right"); [a.text(i, x / 2, f"{x:.2f}", ha="center", va="center", fontsize=8, color="w") for i, x in enumerate(v)]
ax[0].set_title("C1p level 1: drag, five reconstructions", fontsize=9); ax[1].set_title("lift (F1, F-actual: nodal p; F2, F3b: element-mean p; F4: weak form)", fontsize=8)
fig.tight_layout(); fig.savefig(F + "/fig4_stageC_force_reconstructions.png"); plt.close(fig)
print("figures written:", sorted(os.listdir(F)))
