"""3dTube level 1, implicit Euler at dt = 1e-4 s: peak radial displacement, axial trough and late radial
trough at point A, Kratos against the published Euler solutions and solids4foam. The second-order
(backward) solutions are shown in the u_r,max panel for the size of the Euler damping.

Values from results/level1_euler/summary.md (Kratos), the solids4foam driver's reference JSON
(Lozovskiy, Eken, Tukovic) and paper-data/3dtube_level3_findings.md (solids4foam).

Usage: pvpython plot_tube_qoi.py <out.png>
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

INK, INK2, GRID, BAND, BLUE = "#1f1f1e", "#5f5e5a", "#e4e3df", "#efeeea", "#2a78d6"

ROWS = ["Kratos 10.3.0 (this work)", "solids4foam Euler L1", "Lozovskiy et al. (2019)", "Eken (2016)"]
EULER = {  # magnitudes in mm
    "u_r,max(A)": [0.12334, 0.12241, 0.1194, 0.1242],
    "|u_z,min(A)|": [0.07154, 0.07601, 0.0740, 0.0794],
    "|u_r,min(A)|, t > 14 ms": [0.07689, 0.07756, 0.0755, 0.0763],
}
NOTE = {"|u_z,min(A)|": "solids4foam value from before its solver fixes",
        "|u_r,min(A)|, t > 14 ms": "solids4foam value from before its solver fixes"}
BACKWARD = [("solids4foam backward L3", 0.15851), ("Tuković et al. (2018), backward", 0.15694)]


def main():
    out = sys.argv[1]
    plt.rcParams.update({"font.size": 9.5, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                         "ytick.color": INK2, "axes.linewidth": 0.6, "font.family": "sans-serif"})
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 3.9), dpi=160, gridspec_kw={"width_ratios": [1.45, 1, 1]})
    fig.patch.set_facecolor("white")
    for col, (ax, (name, vals)) in enumerate(zip(axes, EULER.items())):
        rows = list(ROWS)
        values = list(vals)
        if col == 0:
            rows += [r for r, _ in BACKWARD]
            values += [v for _, v in BACKWARD]
        y = list(range(len(rows)))[::-1]
        ax.grid(True, axis="x", color=GRID, linewidth=0.6)
        ax.set_axisbelow(True)
        for side in ("top", "right", "left"):
            ax.spines[side].set_visible(False)
        lo, hi = min(vals[2], vals[3]), max(vals[2], vals[3])
        ax.axvspan(lo, hi, color=BAND, zorder=0)
        for yi, r, v in zip(y, rows, values):
            kratos = r.startswith("Kratos")
            backward = "backward" in r
            ax.plot([v], [yi], "o", markersize=8 if kratos else 6.5,
                    color=BLUE if kratos else (INK2 if not backward else "white"),
                    markeredgecolor=BLUE if kratos else INK2, markeredgewidth=1.2, zorder=3)
            ax.annotate(f"{v:.4f}", xy=(v, yi), xytext=(7, 0), textcoords="offset points", va="center",
                        fontsize=8, color=INK if kratos else INK2, fontweight="bold" if kratos else "normal")
        ax.set_yticks(y)
        ax.set_yticklabels(rows if col == 0 else [""] * len(rows), fontsize=8.5)
        if col == 0:
            for lab in ax.get_yticklabels():
                if lab.get_text().startswith("Kratos"):
                    lab.set_color(BLUE)
                    lab.set_fontweight("bold")
        ax.tick_params(axis="y", length=0)
        # The same relative scale (about +-12% of the Kratos value) for the Euler comparison in every
        # panel, so the spreads can be compared by eye; the first panel extends to the backward values
        k = vals[0]
        ax.set_xlim(0.88 * k, max(1.12 * k, 1.07 * max(values)))
        ax.set_ylim(-0.6, len(rows) - 0.4)
        ax.set_xlabel(f"{name} (mm)")
        title = {"u_r,max(A)": "Peak radial displacement", "|u_z,min(A)|": "Axial trough",
                 "|u_r,min(A)|, t > 14 ms": "Late radial trough"}[name]
        ax.set_title(title, loc="left", color=INK, fontsize=10)
        if name in NOTE:
            ax.text(0, -0.25, NOTE[name], transform=ax.transAxes, fontsize=7.5, color=INK2)
    # align the shorter panels on the Euler rows of the first panel
    for ax in axes[1:]:
        ax.set_ylim(-0.6 - 2, len(ROWS) - 0.4)
        ax.set_yticks([])
    fig.tight_layout()
    fig.text(0.01, 0.0, "Implicit Euler, Δt = 1e-4 s, level-1 meshes. Shaded: range of the two published Euler "
             "solutions. Hollow: second-order (backward) solutions, for the size of the Euler damping. "
             "Each panel spans about ±12% of the Kratos value.",
             fontsize=8, color=INK2)
    fig.savefig(out, facecolor="white", bbox_inches="tight", pad_inches=0.15)
    print("wrote", out)


if __name__ == "__main__":
    main()
