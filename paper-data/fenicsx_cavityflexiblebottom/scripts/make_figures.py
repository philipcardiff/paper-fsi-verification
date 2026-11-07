"""Figures from the CSVs in ../data (numpy and matplotlib).  usage: python3 -B -I scripts/make_figures.py"""
import csv, os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); D = os.path.join(HERE, "..", "data"); F = os.path.join(HERE, "..", "figures")
rd = lambda n: list(csv.DictReader(open(os.path.join(D, n))))
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False, "pdf.compression": 0 if False else 6})
lv = rd("coupled_levels.csv"); o = {(r["quantity"], r["triple"]): r for r in rd("coupled_orders.csv")}
ref = {r["reference_variant"]: r for r in rd("reference_comparison.csv")}
dx = [float(r["dx_m"]) for r in lv]
fig, ax = plt.subplots(1, 2, figsize=(6.6, 3.0))
for a, k, q, lab in ((ax[0], "uy_A_m", "uy_A", r"$u_y(A)$ (m)"), (ax[1], "Fy_N_per_m", "Fy_total", r"$F_y$ (N/m)")):
    a.plot(dx, [float(r[k]) for r in lv], "o-", color="#1f4e79", label="FEniCSx M1, levels L1-L4")
    lim = o[(q, "L2-L4")]["richardson_limit"]
    a.axhline(float(lim), color="#1f4e79", ls="--", lw=0.8, label="FEniCSx L2-L4 limit")
    key = "ref_uy_limit" if q == "uy_A" else "ref_Fy_limit"
    for v, c, mk in (("partitioned_strict", "#b05a00", "s"), ("monolithic_618", "#4d7d2f", "^")):
        a.axhline(float(ref[v][key]), color=c, ls=":", lw=1.0, label="s4f limit: " + v.replace("_", " "))
    a.set_xscale("log"); a.set_xticks(dx); a.set_xticklabels(["0.1", "0.05", "0.025", "0.0125"]); a.minorticks_off(); a.invert_xaxis(); a.set_xlabel(r"$\Delta x$ (m)"); a.set_ylabel(lab)
h_, l_ = ax[0].get_legend_handles_labels(); fig.legend(h_, l_, fontsize=7, frameon=False, loc="lower center", ncol=2)
fig.tight_layout(rect=(0, 0.17, 1, 1)); fig.savefig(os.path.join(F, "fenicsx_cavity_spatial.pdf"), metadata={"CreationDate": None, "Creator": None, "Producer": None}); fig.savefig(os.path.join(F, "fenicsx_cavity_spatial.png"), dpi=150)
fig, ax = plt.subplots(figsize=(5.6, 2.4))
names = ["partitioned_strict", "partitioned_legacy_test", "monolithic_pre597", "monolithic_597", "monolithic_618"]
for i, n in enumerate(names):
    ax.plot(float(ref[n]["FEniCSx_M1_limit_minus_ref_limit_uy_percent"]), i, "o", color="#1f4e79")
ax.axvspan(-0.29, 0.29, color="0.9", zorder=0, label="s4f GCI (about 0.29 %)")
ax.set_yticks(range(len(names))); ax.set_yticklabels([n.replace("_", " ") for n in names]); ax.invert_yaxis()
ax.set_xlabel(r"FEniCSx M1 limit minus s4f variant limit, $u_y(A)$ (%)"); ax.legend(fontsize=7, frameon=False, loc="lower right")
fig.tight_layout(); fig.savefig(os.path.join(F, "fenicsx_cavity_reference_variants.pdf"), metadata={"CreationDate": None, "Creator": None, "Producer": None}); fig.savefig(os.path.join(F, "fenicsx_cavity_reference_variants.png"), dpi=150)
print("ok")
