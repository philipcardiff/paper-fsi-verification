#!/usr/bin/env python3
"""Geometry sketch of the radial-torsional ring (vector PDF; matplotlib + numpy only; no data, no CFD).
Every number in the sketch is a documented setting of the case (see the section text).
Usage: python3 make_rotational_geometry.py [OUT_DIR]      (default OUT_DIR = ../../figures)"""
import math
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge, Circle, FancyArrowPatch, Rectangle

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..")
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "..", "figures")
plt.rcParams.update({"font.size": 8, "axes.linewidth": 0.6, "font.family": "serif", "mathtext.fontset": "cm",
                     "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42})
A, B, C = 0.5, 1.0, 2.0
SOLID, FLUID = "#E8DCC8", "#D6E8F5"
MEM = {"radial": "#009E73", "torsion": "#0072B2", "combined": "#D55E00"}
MESHSTY = {"ref": dict(color="k", lw=1.5, ls="-", label="reference"),
           "M1": dict(color="#CC79A7", lw=1.0, ls=":", label="M1"),
           "M2": dict(color="#E69F00", lw=1.0, ls="--", label="M2"),
           "M3": dict(color="#0072B2", lw=0.9, ls="-", label="M3")}


def S_cinf(tau):
    tau = np.clip(np.asarray(tau, float), 0.0, 1.0)
    with np.errstate(divide="ignore", over="ignore", invalid="ignore"):
        f = np.where(tau > 0, np.exp(-1.0 / np.where(tau > 0, tau, 1.0)), 0.0)
        g = np.where(tau < 1, np.exp(-1.0 / np.where(tau < 1, 1.0 - tau, 1.0)), 0.0)
    return f / (f + g)


def arrow(ax, p0, p1, **kw):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle=kw.pop("style", "-|>"), mutation_scale=kw.pop("ms", 7),
                                 lw=kw.pop("lw", 0.8), color=kw.pop("color", "k"), shrinkA=0, shrinkB=0, **kw))


def pol(r, deg):
    a = math.radians(deg)
    return r * math.cos(a), r * math.sin(a)


# ------------------------------------------------------------------ geometry and loading
def geometry():
    fig = plt.figure(figsize=(6.4, 5.3))
    ax = fig.add_axes([0.0, 0.33, 1.0, 0.66])
    ax.set_aspect("equal"); ax.axis("off")
    ax.set_xlim(-3.75, 3.75); ax.set_ylim(-2.3, 2.45)
    ax.add_patch(Wedge((0, 0), C, 0, 360, width=C - B, fc=FLUID, ec="none"))
    ax.add_patch(Wedge((0, 0), B, 0, 360, width=B - A, fc=SOLID, ec="#9C8B6E", hatch="////", lw=0))
    ax.add_patch(Circle((0, 0), A, fc="white", ec="none"))
    ax.add_patch(Circle((0, 0), A, fc="none", ec="k", lw=1.8))                       # inner boundary (prescribed rotation)
    ax.add_patch(Circle((0, 0), B, fc="none", ec="k", lw=1.0))                       # interface at t = 0
    ax.add_patch(Circle((0, 0), 1.10, fc="none", ec="k", lw=0.7, ls=(0, (3, 2))))    # deformed interface b(t) (schematic offset)
    ax.add_patch(Circle((0, 0), C, fc="none", ec="k", lw=1.0, ls=(0, (6, 2))))       # outer open boundary
    # prescribed rotation of the inner boundary
    arrow(ax, pol(0.30, 25), pol(0.30, 155), connectionstyle="arc3,rad=0.6", lw=1.1, ms=7, color=MEM["torsion"])
    ax.text(0, -0.17, r"$\alpha_A(t)$", ha="center", va="center", fontsize=8.5, color=MEM["torsion"])
    # region labels
    ax.text(*pol(0.77, 120), "solid", ha="center", va="center", fontsize=8,
            bbox=dict(boxstyle="round,pad=0.10", fc="white", ec="none", alpha=0.9))
    ax.text(*pol(1.52, 215), "fluid", ha="center", va="center", fontsize=8,
            bbox=dict(boxstyle="round,pad=0.10", fc="white", ec="none", alpha=0.9))
    # coordinates: radial and azimuthal unit vectors at a point in the fluid
    th = 28; P = pol(1.55, th)
    er = (math.cos(math.radians(th)), math.sin(math.radians(th))); et = (-er[1], er[0])
    arrow(ax, P, (P[0] + 0.42 * er[0], P[1] + 0.42 * er[1]), lw=1.0)
    arrow(ax, P, (P[0] + 0.42 * et[0], P[1] + 0.42 * et[1]), lw=1.0)
    ax.text(P[0] + 0.56 * er[0], P[1] + 0.56 * er[1] - 0.03, r"$\mathbf{e}_r$", fontsize=8.5, ha="left", va="center")
    ax.text(P[0] + 0.56 * et[0] - 0.05, P[1] + 0.56 * et[1] + 0.02, r"$\mathbf{e}_\theta$", fontsize=8.5, ha="center", va="bottom")
    # radii
    d = -62
    ax.plot(*zip((0, 0), pol(C, d)), color="0.35", lw=0.5)
    for r, lab in ((A, r"$A=0.5$ m"), (B, r"$B=1$ m"), (C, r"$c=2$ m")):
        x, y = pol(r, d)
        ax.plot([x], [y], marker="o", color="0.2", ms=2.5)
        ax.text(x + 0.07, y - 0.02, lab, fontsize=7, ha="left", va="center", color="0.1", bbox=dict(boxstyle="round,pad=0.08", fc="white", ec="none", alpha=0.85))
    ax.text(0.03, 0.0, r"$O$", fontsize=7, ha="left", va="bottom", color="0.3")
    # callouts
    kw = dict(arrowstyle="-", lw=0.5, color="0.3", shrinkA=0, shrinkB=1)
    ax.annotate("outer boundary $r=c$ (open):\n$p=p_c(t)$,  $\\sigma_{r\\theta}=0$ (tangentially\ntraction-free),\n$\\partial_r v_r=-v_r/c$ (exact)",
                xy=pol(C, 143), xytext=(-3.7, 1.95), fontsize=7, ha="left", va="top", arrowprops=kw)
    ax.annotate("inner boundary $R=A$:\nclamped radially ($u_r=0$),\nprescribed rotation $\\alpha_A(t)$\n(fixedRotation)",
                xy=pol(A, 205), xytext=(-3.7, -1.15), fontsize=7, ha="left", va="top", arrowprops=kw)
    ax.annotate("interface $r=b(t)=\\varphi(B,t)$\n($b(0)=B$; dashed: deformed):\n$\\mathbf{v}_f=\\partial_t\\mathbf{x}_s$ (no slip),\n$\\sigma_s\\mathbf{n}=\\sigma_f\\mathbf{n}$ (traction\ncontinuity)",
                xy=pol(1.10, 8), xytext=(2.12, 0.62), fontsize=7, ha="left", va="top", arrowprops=kw)
    ax.text(3.7, 2.3, "plane strain, per unit depth;\nrotationally invariant ($n=0$) branch;\nneo-Hookean solid, incompressible\nNewtonian fluid, ALE fluid mesh", fontsize=6.8, ha="right", va="top", color="0.25")
    ax.text(-3.72, 2.42, "(a)", fontsize=9, ha="left", va="top")
    # loading members ---------------------------------------------------
    t = np.linspace(0, 2.0, 801)
    S = S_cinf(t / 1.0)
    alpha = 0.1 * S * np.sin(3.0 * t)
    pc = -3.0 * S
    bx = fig.add_axes([0.095, 0.115, 0.38, 0.17]); cx = fig.add_axes([0.575, 0.115, 0.38, 0.17])
    bx.plot(t, alpha, color=MEM["torsion"], lw=1.8, label="torsion")
    bx.plot(t, alpha, color=MEM["combined"], lw=0.9, ls=(0, (3, 2)), label="combined")
    bx.plot(t, 0 * t, color=MEM["radial"], lw=1.2, label="radial")
    bx.set_ylabel(r"$\alpha_A$ (rad)"); bx.set_xlim(0, 2); bx.set_ylim(-0.125, 0.125); bx.set_xlabel(r"$t$ (s)")
    bx.set_title(r"(b) inner rotation $\alpha_A=\theta_0 S(t)\sin\Omega t$", fontsize=7.5, loc="left")
    bx.legend(frameon=False, fontsize=6.5, ncol=3, loc="upper right", handlelength=1.6, columnspacing=0.8, borderpad=0.1)
    cx.plot(t, pc, color=MEM["radial"], lw=1.8, label="radial")
    cx.plot(t, pc, color=MEM["combined"], lw=0.9, ls=(0, (3, 2)), label="combined")
    cx.plot(t, 0 * t, color=MEM["torsion"], lw=1.2, label="torsion")
    cx.set_ylabel(r"$p_c$ (kPa)"); cx.set_xlabel(r"$t$ (s)"); cx.set_xlim(0, 2); cx.set_ylim(-3.6, 0.6)
    cx.set_title(r"(c) outer pressure $p_c=P_0S(t)$", fontsize=7.5, loc="left")
    cx.legend(frameon=False, fontsize=6.5, ncol=3, loc="center right", handlelength=1.6, columnspacing=0.8, borderpad=0.1)
    for a_ in (bx, cx):
        a_.grid(True, lw=0.3, color="0.88"); a_.tick_params(labelsize=7)
    fig.text(0.5, 0.005, r"members: radial ($\theta_0=0$, $P_0=-3$ kPa), torsion ($\theta_0=0.1$ rad, $P_0=0$), combined (both); "
             r"$\Omega=3$ rad/s; $S$: $C^\infty$ ramp, $T_\mathrm{ramp}=1$ s", fontsize=6.5, va="bottom", ha="center")
    save(fig, "rotational_geometry.pdf")


def save(fig, name):
    fig.savefig(os.path.join(OUT, name), metadata={"CreationDate": None, "ModDate": None}); plt.close(fig)


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    geometry()
