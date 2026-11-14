"""Step-2 figure: point-A histories of the Kratos BDF2 levels 1 and 2 and the level-1 Euler run, with
Tukovic et al. (2018); and the incident front at A, Kratos BDF2 level 1 against a solids4foam backward
level-1 run.

Usage: python plot_step2.py <reference dir> <out.png> <euler L1 history> <bdf2 L1 history>
       <bdf2 L2 history> <pointA_L1.csv>
"""
import csv
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

INK, INK2, GRID = "#1f1f1e", "#5f5e5a", "#e4e3df"
BLUE, ORANGE, AQUA, YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"   # validated slots 1-4


def read_kratos(path):
    rows = list(csv.DictReader(open(path)))
    return ([float(r["time"]) * 1e3 for r in rows], [float(r["ur_z250"]) * 1e3 for r in rows],
            [float(r["uz_z250"]) * 1e3 for r in rows])


def read_ref(path):
    rows = [r for r in csv.DictReader(l for l in open(path) if not l.startswith("#"))]
    return [float(r["time_s"]) * 1e3 for r in rows], [float(r["radial_displacement_m"]) * 1e3 for r in rows]


def style(ax, xlabel="time (ms)"):
    ax.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.set_xlabel(xlabel)


def main():
    ref_dir, out, euler, l1, l2, front = sys.argv[1:7]
    plt.rcParams.update({"font.size": 9.5, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                         "ytick.color": INK2, "axes.linewidth": 0.6, "font.family": "sans-serif"})
    fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(15, 4.3), dpi=160,
                                     gridspec_kw={"width_ratios": [1.15, 1.15, 1.0]})
    fig.patch.set_facecolor("white")
    tk = read_ref(os.path.join(ref_dir, "Tukovic2018_fig25a_pointA.csv"))
    series = [("Kratos BDF2, level 2, Δt 1.25e-5", read_kratos(l2), BLUE, "-", 2.0),
              ("Kratos BDF2, level 1, Δt 2.5e-5", read_kratos(l1), ORANGE, (0, (4, 2)), 1.6),
              ("Kratos implicit Euler, level 1, Δt 1e-4", read_kratos(euler), AQUA, (0, (1, 1.5)), 1.8)]
    for label, (t, rad, axl), color, ls, lw in series:
        a1.plot(t, rad, color=color, linestyle=ls, linewidth=lw, label=label)
        a2.plot(t, axl, color=color, linestyle=ls, linewidth=lw)
    a1.plot(*tk, color=YELLOW, linestyle=(0, (6, 2, 1, 2)), linewidth=1.6,
            label="Tuković et al. (2018), backward, Δt 2.5e-5")
    for ax in (a1, a2):
        style(ax)
        ax.set_xlim(0, 20)
    a1.set_ylabel("radial displacement at A (mm)")
    a2.set_ylabel("axial displacement at A (mm)")
    a1.set_title("Radial displacement, point A", loc="left", color=INK, fontsize=10)
    a2.set_title("Axial displacement, point A", loc="left", color=INK, fontsize=10)

    rows = list(csv.DictReader(open(front)))
    t = [float(r["time_s"]) * 1e3 for r in rows]
    k = [float(r["ur_kratos_m"]) * 1e3 for r in rows]
    s = [float(r["ur_s4f_m"]) * 1e3 for r in rows]
    a3.plot(t, k, color=ORANGE, linestyle=(0, (4, 2)), linewidth=1.6, label="Kratos BDF2, level 1")
    a3.plot(t, s, color=INK2, linestyle="-", linewidth=1.4, label="solids4foam backward, level 1")
    style(a3)
    a3.set_xlim(3, 8)
    a3.set_ylim(0, 0.17)
    a3.axhline(0.5 * max(k), color=GRID, linewidth=1.0, zorder=0)
    a3.text(3.1, 0.5 * max(k) + 0.003, "50% of the incident peak", color=INK2, fontsize=8.5)
    a3.set_ylabel("radial displacement at A (mm)")
    a3.set_title("Incident front at A, level 1", loc="left", color=INK, fontsize=10)
    a3.legend(loc="upper left", frameon=False, fontsize=8.5, labelcolor=INK)
    fig.legend(*a1.get_legend_handles_labels(), loc="lower center", bbox_to_anchor=(0.37, -0.1), ncol=2,
               frameon=False, fontsize=9, labelcolor=INK)
    fig.tight_layout()
    fig.savefig(out, facecolor="white", bbox_inches="tight", pad_inches=0.15)
    print("wrote", out)


if __name__ == "__main__":
    main()
