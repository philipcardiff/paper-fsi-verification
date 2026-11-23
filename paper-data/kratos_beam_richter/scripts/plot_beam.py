#!/usr/bin/env python3
"""Figures for the Kratos beamInCrossFlow runs (need matplotlib).

levels.png:    steady u_x(A) and F_x against the level for both inflow peaks, Kratos levels 0 to 2 (F_x
               from the reactions, Fc, and the direct traction integral, Ff) with solids4foam levels 0 to 3
               and Richter (2012).
histories.png: u_x(A) and F_x against time for the six main runs, the march to the steady state.

Usage (from this folder): python scripts/plot_beam.py --runs runs --s4f ../beam_richter --out figures
"""
import argparse
import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

INK, INK2, GRID = "#1f1f1e", "#5f5e5a", "#e4e3df"
BLUE, ORANGE, AQUA, YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"   # validated slots 1-4
RICHTER = {"ux": 5.95e-5, "Fx": 1.33}
RUNS = {(0, 0.2): "f1_U0.2_dt0.05", (1, 0.2): "f2_U0.2_dt0.05", (2, 0.2): "f4_U0.2_dt0.05_amg",
        (0, 0.3): "f1_U0.3_dt0.05", (1, 0.3): "f2_U0.3_dt0.05", (2, 0.3): "f4_U0.3_dt0.05_amg"}


def history(run):
    rows = list(csv.DictReader(open(os.path.join(run, "history.csv"))))
    col = lambda k: [float(r[k]) for r in rows]
    return {"t": col("time"), "ux": col("ux_A"), "Fcx": col("Fcx"), "Ffx": col("Ffx")}


def s4f(d):
    out = {0.2: {}, 0.3: {}}
    for r in csv.DictReader(open(os.path.join(d, "r02_sequence.csv"))):
        out[0.2][int(r["level"][1:])] = (float(r["ux"]), float(r["Fx"]))
    for r in csv.DictReader(open(os.path.join(d, "production_sequence.csv"))):
        out[0.3][int(r["level"])] = (float(r["ux"]), float(r["fx"]))
    return out


def style(ax, xlabel):
    ax.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.set_xlabel(xlabel)


def levels_figure(runs, ref, out):
    fig, axes = plt.subplots(2, 2, figsize=(10, 7.2), dpi=160)
    for j, U in enumerate((0.2, 0.3)):
        kl = [0, 1, 2]
        k = [history(os.path.join(runs, RUNS[(lev, U)])) for lev in kl]
        sl = sorted(ref[U])
        a_u, a_f = axes[0, j], axes[1, j]
        a_u.plot(sl, [1e5 * ref[U][l][0] for l in sl], color=INK2, linewidth=1.4, marker="o", markersize=5,
                 label="solids4foam, levels 0–3")
        a_u.plot(kl, [1e5 * h["ux"][-1] for h in k], color=BLUE, linewidth=2.0, marker="s", markersize=6,
                 label="Kratos, levels 0–2")
        a_f.plot(sl, [ref[U][l][1] for l in sl], color=INK2, linewidth=1.4, marker="o", markersize=5,
                 label="solids4foam, levels 0–3")
        a_f.plot(kl, [h["Fcx"][-1] for h in k], color=BLUE, linewidth=2.0, marker="s", markersize=6,
                 label="Kratos Fc (interface reactions)")
        a_f.plot(kl, [h["Ffx"][-1] for h in k], color=ORANGE, linewidth=2.0, linestyle=(0, (4, 2)),
                 marker="D", markersize=5.5, label="Kratos Ff (direct traction integral)")
        if U == 0.2:
            a_u.axhline(1e5 * RICHTER["ux"], color=YELLOW, linewidth=1.6, linestyle=(0, (1.5, 2)),
                        label="Richter (2012)")
            a_f.axhline(RICHTER["Fx"], color=YELLOW, linewidth=1.6, linestyle=(0, (1.5, 2)), label="Richter (2012)")
        for ax in (a_u, a_f):
            style(ax, "level (solids4foam numbering)")
            ax.set_xticks([0, 1, 2, 3])
        a_u.set_ylabel("u_x(A) (1e-5 m)")
        a_f.set_ylabel("F_x on the half plate (N)")
        name = "corrected" if U == 0.2 else "as printed"
        a_u.set_title(f"u_x(A), peak {U} m/s ({name})", loc="left", color=INK, fontsize=10)
        a_f.set_title(f"F_x, peak {U} m/s ({name})", loc="left", color=INK, fontsize=10)
    axes[0, 0].legend(loc="lower right", frameon=False, fontsize=8.5, labelcolor=INK)
    axes[1, 0].legend(loc="lower right", frameon=False, fontsize=8.5, labelcolor=INK)
    fig.tight_layout()
    fig.savefig(out, facecolor="white", bbox_inches="tight", pad_inches=0.15)
    print("wrote", out)


def histories_figure(runs, out):
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.2), dpi=160)
    styles = {2: ("-", 2.0), 1: ((0, (4, 2)), 1.6), 0: ((0, (1, 1.5)), 1.5)}
    series = [(f"level {lev}, peak {U}", RUNS[(lev, U)], color, *styles[lev])
              for U, color in ((0.2, BLUE), (0.3, ORANGE)) for lev in (2, 1, 0)]
    for label, run, color, ls, lw in series:
        h = history(os.path.join(runs, run))
        a1.plot(h["t"], [1e5 * v for v in h["ux"]], color=color, linestyle=ls, linewidth=lw, label=label)
        a2.plot(h["t"], h["Fcx"], color=color, linestyle=ls, linewidth=lw)
    for ax in (a1, a2):
        style(ax, "time (s)")
        ax.set_xlim(0, 19)
    a1.set_ylabel("u_x(A) (1e-5 m)")
    a2.set_ylabel("Fc_x on the half plate (N)")
    a1.set_title("u_x(A), implicit Euler, Δt 0.05 s", loc="left", color=INK, fontsize=10)
    a2.set_title("F_x from the interface reactions", loc="left", color=INK, fontsize=10)
    a1.legend(loc="upper right", ncol=2, frameon=False, fontsize=8.5, labelcolor=INK)
    fig.tight_layout()
    fig.savefig(out, facecolor="white", bbox_inches="tight", pad_inches=0.15)
    print("wrote", out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", default="runs")
    ap.add_argument("--s4f", default="../beam_richter")
    ap.add_argument("--out", default="figures")
    a = ap.parse_args()
    plt.rcParams.update({"font.size": 9.5, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                         "ytick.color": INK2, "axes.linewidth": 0.6, "font.family": "sans-serif"})
    levels_figure(a.runs, s4f(a.s4f), os.path.join(a.out, "levels.png"))
    histories_figure(a.runs, os.path.join(a.out, "histories.png"))


if __name__ == "__main__":
    main()
