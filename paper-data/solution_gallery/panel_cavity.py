"""Gallery panel: flexible-bottom cavity, steady state.

Data: a serial rerun of this paper's skew-corrected level-2 case (Aitken
coupling) to the final time t = 200 s, written at that time only.

The fluid is coloured by the speed, with streamlines; the elastic plate is
drawn deformed at true scale, its undeformed outline dashed.

Usage: panel_cavity.py CASE OUT.pdf [--time T] [--xmax X]"""
import argparse
import sys

import numpy as np

import common
import fsi2d

TITLE = "Flexible-bottom cavity"
XLIM, YLIM = (-0.05, 8.05), (-0.45, 2.05)
ASPECT = (YLIM[1] - YLIM[0]) / (XLIM[1] - XLIM[0])


def load(case, time=None):
    return fsi2d.load(case, time)


def draw(fig, ax, cax, d, letter="g"):
    sp = fsi2d.speed(d)
    vmax = common.nice_ceil(sp.max())
    fsi2d.draw(fig, ax, cax, d, sp, common.Normalize(0, vmax), "viridis",
               r"$|\mathbf{u}|$ (m/s)", XLIM, YLIM, stream_density=1.1)
    common.scalebar(ax, 6.6, -0.32, 1.0, "1 m", above=True)
    common.title(ax, letter, f"{TITLE}, $t = {d['t']:g}$ s")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("case")
    ap.add_argument("out")
    ap.add_argument("--time", type=float)
    a = ap.parse_args()
    d = load(a.case, a.time)
    print(f"t = {d['t']}, max speed {fsi2d.speed(d).max():.4f}, "
          f"min D_y {d['D'][:, 1].min():.5f}")
    common.use_style()
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(3.0, 1.2))
    ax = fig.add_axes([0.01, 0.03, 0.80, 0.82])
    cax = fig.add_axes([0.84, 0.15, 0.03, 0.6])
    draw(fig, ax, cax, d)
    fig.savefig(a.out, dpi=400)
    print("wrote", a.out)


if __name__ == "__main__":
    sys.exit(main())
