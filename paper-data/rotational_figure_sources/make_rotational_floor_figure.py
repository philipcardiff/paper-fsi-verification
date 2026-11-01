"""Figure rotational_floor.pdf: reference-error convergence of the radial--torsional ring.

Input : ../rotational_floor_formulation.csv  (max-norm error / max|ref|, 0-2 s, dt = 2.5 ms)
        formulation A = tutorial as stored (cell leastSquares grad(p) and grad(U), stock newMovingWallVelocity)
        formulation B = A with only grad(p) changed to Gauss linear
        The A values equal the entries of rotational_convergence.csv.
Output: ../../figures/rotational_floor.pdf
Run   : python3 make_rotational_floor_figure.py
"""
import csv, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

here = os.path.dirname(os.path.abspath(__file__))
rows = list(csv.DictReader(open(os.path.join(here, "..", "rotational_floor_formulation.csv"))))

def series(form, member, qoi):
    d = {r["mesh"]: r for r in rows if (r["formulation"], r["member"], r["qoi"]) == (form, member, qoi)}
    e = [float(d[m]["E_max_over_peak"]) for m in ("M1", "M2", "M3")]
    n = [float(d[m]["noise_over_peak"]) if d[m]["noise_over_peak"] else 0.0 for m in ("M1", "M2", "M3")]
    return e, n

h = [4.0, 2.0, 1.0]                      # mesh spacing relative to M3
panels = [("alpha_B", r"interface rotation $\alpha_B$"),
          ("T_B", r"wall torque $T_B$"),
          ("H_f", r"fluid angular momentum $H_f$")]
styles = [("A", "torsion", "torsion, A", "#0072B2", "o", "-"),
          ("A", "combined", "combined, A (tutorial)", "#D55E00", "s", "-"),
          ("B", "combined", r"combined, B (Gauss $\nabla p$)", "#009E73", "^", "--")]

plt.rcParams.update({"font.size": 8, "axes.linewidth": 0.6, "pdf.fonttype": 42})
fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.7), sharey=True)
for ax, (qoi, title) in zip(axs, panels):
    for form, member, label, col, mk, ls in styles:
        e, n = series(form, member, qoi)
        if form != "A":  # the noise estimate (partition, coupling tolerance) was made for formulation A only
            n = [0.0] * 3
        if any(n):
            lo = [min(ni, 0.8 * ei) for ni, ei in zip(n, e)]
            ax.errorbar(h, e, yerr=[lo, n], fmt="none", ecolor=col, elinewidth=0.8, capsize=1.5, zorder=2)
        filled = [ei > 2 * ni for ei, ni in zip(e, n)]
        ax.plot(h, e, ls, color=col, lw=1.0, zorder=3)
        for hi, ei, fi in zip(h, e, filled):
            ax.plot([hi], [ei], mk, ms=4.5, mfc=col if fi else "white", mec=col, mew=1.0, zorder=4)
        ax.plot([], [], ls, color=col, marker=mk, ms=4.5, label=label)
    e0, _ = series("A", "torsion", qoi)
    ax.plot([4, 1], [e0[0], e0[0] / 16], ":", color="0.45", lw=0.9, zorder=1)
    ax.text(3.9, e0[0] / 3.2, "slope 2 (guide)", color="0.35", fontsize=6.5, ha="left", va="top", rotation=-34)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xticks([1, 2, 4]); ax.set_xticklabels(["M3", "M2", "M1"])
    ax.minorticks_off()
    ax.invert_xaxis()
    ax.set_title(title, fontsize=8)
    ax.grid(True, which="major", lw=0.3, color="0.85")
axs[0].set_ylabel(r"$\max_t|q-q_\mathrm{ref}|\,/\,\max_t|q_\mathrm{ref}|$")
axs[0].set_ylim(4e-4, 4e-2)
axs[1].set_xlabel("mesh (cell size halves from M1 to M3)", labelpad=2)
h_, l_ = axs[0].get_legend_handles_labels()
fig.legend(h_, l_, loc="lower center", ncol=3, fontsize=7, frameon=False)
fig.tight_layout(rect=(0, 0.09, 1, 1))
out = os.path.join(here, "..", "..", "figures", "rotational_floor.pdf")
fig.savefig(out)
print("wrote", os.path.normpath(out))
