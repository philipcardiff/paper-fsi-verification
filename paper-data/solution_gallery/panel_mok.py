"""Gallery panel: Mok lid-driven cavity at a peak of the membrane deflection.

Data: a serial rerun of the tutorial's 64^2 member of the verification mesh
study (64 x 64 fluid cells, 64 x 8 high-order solid, dt = 0.05 s, IQN-ILS;
tutorial at the commit of the paper's provenance table), written every
0.25 s and stopped in the periodic state; the snapshot is the peak at
t = 24.5 s (midpoint deflection 0.2871 m, against the 64^2 peak of 0.2870 m
over t = 20-70 s in the paper's table).

The fluid is coloured by the speed, with streamlines; the membrane is drawn
deformed, its undeformed position dashed.

Usage: panel_mok.py CASE OUT.pdf [--time T]"""
import argparse
import sys

import common
import fsi2d

TITLE = "Mok lid-driven cavity"
XLIM, YLIM = (-0.02, 1.02), (-0.03, 1.07)
ASPECT = (YLIM[1] - YLIM[0]) / (XLIM[1] - XLIM[0])


def load(case, time=24.5):
    return fsi2d.load(case, time)


def draw(fig, ax, cax, d, letter="f"):
    sp = fsi2d.speed(d)
    vmax = common.nice_ceil(sp.max())
    fsi2d.draw(fig, ax, cax, d, sp, common.Normalize(0, vmax), "viridis",
               r"$|\mathbf{u}|$ (m/s)", XLIM, YLIM, stream_density=1.0,
               cb_ticks=[0, 0.1, 0.2, 0.3, 0.4])
    common.scalebar(ax, 0.3, 0.05, 0.4, "0.4 m", above=True)
    # openings (0.875 < y < 1) on both sides, and the lid
    for x in (0.0, 1.0):
        ax.plot([x, x], [0.875, 1.0], color="#E69F00", lw=2.0,
                solid_capstyle="butt", zorder=9)
    ax.annotate("", (0.7, 1.0), (0.3, 1.0), zorder=12,
                arrowprops=dict(arrowstyle="-|>,head_length=0.3,head_width=0.15",
                                lw=0.8, color="black", shrinkA=0, shrinkB=0))
    ax.text(0.5, 1.003, "lid", ha="center", va="bottom", fontsize=8, zorder=12)
    ax.text(0.5, 0.255, "membrane", ha="center", va="top", fontsize=8, zorder=12)
    common.title(ax, letter, f"{TITLE}, $t = {d['t']:g}$ s")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("case")
    ap.add_argument("out")
    ap.add_argument("--time", type=float, default=24.5)
    a = ap.parse_args()
    d = load(a.case, a.time)
    print(f"t = {d['t']}, max speed {fsi2d.speed(d).max():.3f}, "
          f"D_y range {d['D'][:, 1].min():.4f} .. {d['D'][:, 1].max():.4f}")
    common.use_style()
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(2.0, 2.0))
    ax = fig.add_axes([0.01, 0.03, 0.76, 0.82])
    cax = fig.add_axes([0.80, 0.15, 0.04, 0.6])
    draw(fig, ax, cax, d)
    fig.savefig(a.out, dpi=400)
    print("wrote", a.out)


if __name__ == "__main__":
    sys.exit(main())
