"""Geometry sketches of the four 2-D/axisymmetric benchmarks without one in
the paper (fig:geometry_panel): (a) the ring added-mass problem, (b) the Mok
lid-driven cavity, (c) the Womersley tube and (d) the collapsible channel
(dimensions from the solids4foam tutorial READMEs). Paper style
(../figstyle), drawn at the text width.

Usage: geometry_panel.py OUT.pdf"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "figstyle"))
import figstyle  # noqa: E402  (paper style: CM fonts, 8 pt at final size)
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Rectangle, Wedge, Circle, FancyArrowPatch  # noqa: E402

figstyle.use()
FLUID, SOLID, EDGE = "#dbe8f5", "#d62728", "black"
GREY = "0.35"
SMALL = 7


def dim(ax, p, q, text, off=(0, 0), **kw):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle="<->", mutation_scale=7, lw=0.6))
    if text:
        ax.text((p[0] + q[0]) / 2 + off[0], (p[1] + q[1]) / 2 + off[1], text,
                ha="center", va="center", **kw)


def ring(ax):
    ax.add_patch(Circle((0, 0), 1.5, fc=FLUID, ec=EDGE, lw=1.2))
    ax.add_patch(Wedge((0, 0), 1.05, 0, 360, width=0.1, fc=SOLID, ec="none"))
    ax.add_patch(Circle((0, 0), 0.95, fc="white", ec="none"))
    t = np.linspace(0, 2 * np.pi, 300)
    r = 1.0 + 0.08 * np.cos(2 * t)
    ax.plot(r * np.cos(t), r * np.sin(t), "--", color=GREY, lw=0.8)
    ax.plot([0, 1.5], [0, 0], color="0.55", lw=0.5)
    ax.plot([0, 0], [0, 1.5], color="0.55", lw=0.5)
    ax.text(-0.3, 0.3, "vacuum", ha="center", va="center")
    ax.text(-0.28, -0.45, "$n=2$ mode\n(dashed)", ha="center", va="center", color=GREY,
            fontsize=SMALL, linespacing=1.0)
    ax.text(-0.9, -0.86, "fluid", ha="center", va="center")
    ax.annotate("ring, $a=1$ m, $h=0.1$ m", (-0.74, 0.74), (-1.35, 1.78), ha="left", va="top",
                fontsize=SMALL, arrowprops=dict(arrowstyle="-", lw=0.5))
    ax.annotate("rigid wall", (-1.3, -0.75), (-1.78, -1.62), ha="left", va="bottom", fontsize=SMALL,
                arrowprops=dict(arrowstyle="-", lw=0.5))
    ang = -np.pi / 4
    dim(ax, (0, 0), (1.5 * np.cos(ang), 1.5 * np.sin(ang)), "")
    ax.text(1.14, -1.14, "$b=1.5$ m", fontsize=SMALL, ha="left", va="top")
    ax.annotate("modelled\nquadrant", (0.95, 0.95), (1.8, 1.5), ha="right", va="bottom",
                fontsize=SMALL, color=GREY, linespacing=1.0,
                arrowprops=dict(arrowstyle="-", lw=0.5, color=GREY))
    ax.set_xlim(-1.8, 1.8); ax.set_ylim(-1.7, 1.8)


def mok(ax):
    ax.add_patch(Rectangle((0, 0), 1, 1, fc=FLUID, ec="none"))
    ax.plot([0, 0], [0, 0.875], color=EDGE, lw=1.4)
    ax.plot([1, 1], [0, 0.875], color=EDGE, lw=1.4)
    ax.plot([0, 1], [1, 1], color=figstyle.OKABE_ITO[0], lw=2.2)
    x = np.linspace(0, 1, 100)
    ax.fill_between(x, -0.02, 0.0, color=SOLID, lw=0)
    ax.plot(x, -0.25 * np.sin(np.pi * x) ** 2 - 0.01, ":", color=SOLID, lw=0.8)
    ax.plot(0.5, 0, "ko", ms=3); ax.text(0.54, 0.04, "A", fontsize=SMALL)
    for y in (0.9, 0.94, 0.98):
        ax.add_patch(FancyArrowPatch((-0.13, y), (-0.01, y), arrowstyle="->",
                                     mutation_scale=5, lw=0.6))
    ax.add_patch(FancyArrowPatch((1.01, 0.94), (1.13, 0.94), arrowstyle="->",
                                 mutation_scale=5, lw=0.6))
    ax.text(0.5, 1.05, r"lid $v_x=1-\cos(2\pi t/5)$", ha="center", va="bottom", fontsize=SMALL)
    ax.text(-0.04, 0.8, "inflow", ha="right", va="top", fontsize=SMALL)
    ax.text(1.04, 0.8, "$p=0$", ha="left", va="top", fontsize=SMALL)
    ax.text(0.5, 0.55, "fluid", ha="center")
    ax.text(0.5, -0.33, "membrane, $h=0.002$ m,\n$E_S=250$ Pa (clamped)",
            ha="center", va="top", fontsize=SMALL, linespacing=1.0)
    dim(ax, (1.2, 0.875), (1.2, 1.0), "")
    ax.text(1.25, 0.9375, "0.125", fontsize=SMALL, va="center")
    dim(ax, (0, -0.72), (1, -0.72), "1 m", off=(0, -0.08), fontsize=SMALL)
    ax.set_xlim(-0.5, 1.55); ax.set_ylim(-0.86, 1.2)


def tube(ax):
    L, R, h = 15.0, 1.0, 0.1
    s = 3.0  # radial display stretch
    ax.add_patch(Rectangle((0, 0), L, R * s, fc=FLUID, ec="none"))
    ax.add_patch(Rectangle((0, R * s), L, h * s, fc=SOLID, ec="none"))
    ax.plot([0, L], [0, 0], "-.", color=GREY, lw=0.7)
    ax.text(L / 2, -0.2, "axis of symmetry (wedge)", ha="center", va="top", fontsize=SMALL)
    ax.text(L / 2, R * s / 2, "fluid", ha="center", va="center")
    ax.text(L / 2, R * s + h * s + 0.12, "elastic wall, $h=0.1$ m, traction-free outer surface",
            ha="center", va="bottom", fontsize=SMALL)
    ax.text(-0.25, R * s / 2, "inlet:\nexact $p$", ha="right", va="center", fontsize=SMALL,
            linespacing=1.0)
    ax.text(L + 0.25, R * s / 2, "outlet:\nexact $p$", ha="left", va="center", fontsize=SMALL,
            linespacing=1.0)
    dim(ax, (0, -1.15), (L, -1.15), "$L=15$ m", off=(0, -0.38), fontsize=SMALL)
    dim(ax, (L + 2.9, 0), (L + 2.9, R * s), "")
    ax.text(L + 3.1, R * s / 2, "$R=1$ m", rotation=90, ha="left", va="center", fontsize=SMALL)
    ax.set_xlim(-3.0, L + 4.0); ax.set_ylim(-1.9, 3.9)


def channel(ax):
    Lu, L, Ld, a, h = 5.0, 10.0, 10.0, 1.0, 0.05
    sy = 3.0  # transverse display stretch
    ax.add_patch(Rectangle((0, 0), Lu + L + Ld, a * sy, fc=FLUID, ec="none"))
    ax.plot([0, Lu + L + Ld], [0, 0], color=EDGE, lw=1.4)
    ax.plot([0, Lu], [a * sy] * 2, color=EDGE, lw=1.4)
    ax.plot([Lu + L, Lu + L + Ld], [a * sy] * 2, color=EDGE, lw=1.4)
    x = np.linspace(Lu, Lu + L, 200)
    w = -0.28 * sy * np.sin(np.pi * (x - Lu) / L) ** 2
    ax.fill_between(x, a * sy + w, a * sy + w + 0.12, color=SOLID, lw=0)
    ax.plot([Lu, Lu + L], [a * sy + 0.06] * 2, ":", color=GREY, lw=0.6)
    for xi in np.linspace(Lu + 2, Lu + L - 2, 4):
        ax.add_patch(FancyArrowPatch((xi, a * sy + 0.9), (xi, a * sy + 0.3),
                                     arrowstyle="->", mutation_scale=6, lw=0.6))
    ax.text(Lu + L / 2, a * sy + 1.0, r"$p_\mathrm{ext}$ (ramped)", ha="center", va="bottom",
            fontsize=SMALL)
    ax.text(-0.3, a * sy / 2, "$p_\\mathrm{in}$,\nparallel", ha="right", va="center", fontsize=SMALL,
            linespacing=1.0)
    ax.text(Lu + L + Ld + 0.3, a * sy / 2, "$p=0$,\nparallel", ha="left", va="center",
            fontsize=SMALL, linespacing=1.0)
    ax.text(Lu + L + Ld / 2 - 1.0, a * sy / 2, "fluid", ha="center", va="center")
    ax.text(Lu + L / 2, a * sy * 0.3, "elastic wall, $h=0.05$ m (clamped ends)", ha="center",
            va="center", fontsize=SMALL)
    dim(ax, (0, -0.6), (Lu, -0.6), r"$L_\mathrm{up}=5$", off=(0, -0.42), fontsize=SMALL)
    dim(ax, (Lu, -0.6), (Lu + L, -0.6), "$L=10$", off=(0, -0.42), fontsize=SMALL)
    dim(ax, (Lu + L, -0.6), (Lu + L + Ld, -0.6), r"$L_\mathrm{down}=10$", off=(0, -0.42),
        fontsize=SMALL)
    dim(ax, (Lu + L + Ld - 0.8, 0), (Lu + L + Ld - 0.8, a * sy), "")
    ax.text(Lu + L + Ld - 1.05, a * sy / 2, "$a=1$ m", rotation=90, ha="right", va="center",
            fontsize=SMALL)
    ax.set_xlim(-3.2, Lu + L + Ld + 3.2); ax.set_ylim(-1.5, a * sy + 1.5)


fig = plt.figure(figsize=(figstyle.TEXTWIDTH, 5.6))
gs = fig.add_gridspec(3, 2, height_ratios=[1.45, 0.72, 0.72], hspace=0.08, wspace=0.05,
                      left=0.01, right=0.99, top=0.97, bottom=0.01)
for f, spec, lab in ((ring, gs[0, 0], "(a)"), (mok, gs[0, 1], "(b)"), (tube, gs[1, :], "(c)"),
                     (channel, gs[2, :], "(d)")):
    ax = fig.add_subplot(spec); f(ax); ax.set_aspect("equal"); ax.axis("off")
    ax.text(0.0, 1.0, lab, transform=ax.transAxes, ha="left", va="top", fontsize=8)
fig.savefig(sys.argv[1] if len(sys.argv) > 1 else "benchmarks_geometry_panel.pdf")
