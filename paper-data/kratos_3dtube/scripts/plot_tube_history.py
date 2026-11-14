"""Radial and axial displacement at point A: Kratos runs against the published implicit-Euler
solutions (Lozovskiy et al. 2019, Eken 2016) and the BDF2-type solution of Tukovic et al. (2018).

Usage: pvpython plot_tube_history.py <reference dir> <out.png> <label>=<history.csv> [...]
<history.csv> is a Kratos run's history.csv (point A is the z = 0.025 m station).
"""
import csv
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

INK, INK2, GRID = "#1f1f1e", "#5f5e5a", "#e4e3df"
BLUE, ORANGE, AQUA, YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"   # validated slots 1-4


def read_ref(path):
    rows = [r for r in csv.DictReader(l for l in open(path) if not l.startswith("#"))]
    t = [float(r["time_s"]) * 1e3 for r in rows]
    rad = [float(r["radial_displacement_m"]) * 1e3 for r in rows]
    ax = [float(r["axial_displacement_m"]) * 1e3 for r in rows] if "axial_displacement_m" in rows[0] else None
    return t, rad, ax


def read_kratos(path):
    rows = list(csv.DictReader(open(path)))
    return ([float(r["time"]) * 1e3 for r in rows], [float(r["ur_z250"]) * 1e3 for r in rows],
            [float(r["uz_z250"]) * 1e3 for r in rows])


def main():
    ref_dir, out = sys.argv[1], sys.argv[2]
    runs = [a.split("=", 1) for a in sys.argv[3:]]
    plt.rcParams.update({"font.size": 9.5, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                         "ytick.color": INK2, "axes.linewidth": 0.6, "font.family": "sans-serif"})
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.5, 4.2), dpi=160)
    fig.patch.set_facecolor("white")
    series = []
    for label, path in runs:
        series.append((label, read_kratos(path), BLUE, "-", 2.0))
    lz = read_ref(os.path.join(ref_dir, "Lozovskiy2019_fig3_pointA.csv"))
    ek = read_ref(os.path.join(ref_dir, "Eken2016_fig5.6_pointA.csv"))
    tk = read_ref(os.path.join(ref_dir, "Tukovic2018_fig25a_pointA.csv"))
    series += [("Lozovskiy et al. (2019), Euler, Δt 1e-4", lz, ORANGE, (0, (4, 2)), 1.6),
               ("Eken (2016), Euler, Δt 1e-4", ek, AQUA, (0, (1, 1.5)), 1.8),
               ("Tuković et al. (2018), backward, Δt 2.5e-5", tk, YELLOW, (0, (6, 2, 1, 2)), 1.6)]
    for ax in (a1, a2):
        ax.grid(True, color=GRID, linewidth=0.6)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        ax.set_xlabel("time (ms)")
        ax.set_xlim(0, 20)
    for label, (t, rad, axl), color, ls, lw in series:
        a1.plot(t, rad, color=color, linestyle=ls, linewidth=lw, label=label)
        if axl is not None:
            a2.plot(t, axl, color=color, linestyle=ls, linewidth=lw)
    a1.set_ylabel("radial displacement at A (mm)")
    a2.set_ylabel("axial displacement at A (mm)")
    a1.set_title("Radial displacement, point A (inner wall, mid-length)", loc="left", color=INK, fontsize=10)
    a2.set_title("Axial displacement, point A", loc="left", color=INK, fontsize=10)
    fig.legend(*a1.get_legend_handles_labels(), loc="lower center", bbox_to_anchor=(0.5, -0.1), ncol=2,
               frameon=False, fontsize=9, labelcolor=INK)
    fig.tight_layout()
    fig.savefig(out, facecolor="white", bbox_inches="tight", pad_inches=0.15)
    print("wrote", out)


if __name__ == "__main__":
    main()
