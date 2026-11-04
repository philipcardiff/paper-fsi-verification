#!/usr/bin/env python3
"""FSI1 figures from stored results: coupled tip displacement errors (stock vs reaction-corrected candidate) and CFD1 force errors. Run from OVERNIGHT2/."""
import json
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
tip = json.load(open("coupled_tip_results.json"))
levels = [("1.25", 2, "FS1r_%s"), ("1.875", 4, "FS1c_f1.875%s"), ("2.0", 4, "FS1c_f2.0%s"), ("3.0", 6, "FS1c_f3.0%s"), ("3.75", 8, "FS1c_f3.75%s")]
fig, ax = plt.subplots(1, 2, figsize=(10.5, 4))
x = np.arange(len(levels))
for k, (comp, ttl) in enumerate((("uy_err_pct", "tip $u_y(A)$"), ("ux_err_pct", "tip $u_x(A)$"))):
    s = [tip[l[2] % "s"][comp] for l in levels]; p = [tip[l[2] % "p"][comp] for l in levels]
    ax[k].bar(x - 0.18, s, 0.36, color="#d95f02", label="stock FEBio v4.13"); ax[k].bar(x + 0.18, p, 0.36, color="#1b9e77", label="with reaction correction (candidate build)")
    ax[k].set_yscale("symlog", linthresh=1.0); ax[k].axhline(0, color="k", lw=0.6)
    ax[k].set_xticks(x); ax[k].set_xticklabels([f"f={l[0]}\nflag nt={l[1]}" for l in levels], fontsize=8); ax[k].set_ylabel("error against Featflow FSI1 (%), symlog axis"); ax[k].set_title(ttl + " (steady)", fontsize=10)
    for xi, v in zip(x + 0.18, p): ax[k].text(xi, v + (0.3 if v >= 0 else -0.3), f"{v:+.2f}", ha="center", va="bottom" if v >= 0 else "top", fontsize=6.5)
    for xi, v in zip(x - 0.18, s): ax[k].text(xi, v + (0.3 if v >= 0 else -0.3), f"{v:+.0f}" if abs(v) > 20 else f"{v:+.1f}", ha="center", va="bottom" if v >= 0 else "top", fontsize=6.5)
ax[0].legend(fontsize=7, loc="lower left")
fig.tight_layout(); fig.savefig("figures/fsi1_coupled_tip_errors.png", dpi=190); plt.close(fig)
rows = [("1.25", 8330, 15.73, -2.12, 0.31), ("1.5", 11996, 23.70, -1.64, 0.27), ("1.625", 14191, 12.50, -1.71, 0.23), ("1.875", 18606, -0.83, -1.62, 0.16), ("2.0", 21344, 18.85, -1.75, 0.21), ("2.375", 29951, 12.72, -1.28, 0.13), ("2.5", 33368, 22.89, -1.19, 0.12), ("3.0", 48024, 3.68, -1.27, 0.08)]
fig, ax = plt.subplots(figsize=(6.2, 4))
c = [r[1] for r in rows]
ax.plot(c, [r[2] for r in rows], "o-", color="#d95f02", label="F1: nodal-pressure boundary integral (the load FEBio applies)")
ax.plot(c, [r[3] for r in rows], "s-", color="#7570b3", label="F2: element-mean reconstruction")
ax.plot(c, [r[4] for r in rows], "^-", color="#1b9e77", label="F4: consistent nodal reaction")
ax.axhline(0, color="k", lw=0.6); ax.set_xscale("log"); ax.set_xlabel("fluid cells (CFD1 rigid flag)"); ax.set_ylabel("drag error against Featflow CFD1 (%)"); ax.legend(fontsize=6.5)
ax.set_title("CFD1 drag: three force definitions on eight meshes", fontsize=9)
fig.tight_layout(); fig.savefig("figures/fsi1_cfd1_force_definitions.png", dpi=190); plt.close(fig)
print("ok")
