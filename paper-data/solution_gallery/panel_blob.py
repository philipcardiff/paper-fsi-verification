"""Gallery panel: blob in treacle, steady state.

Data: a serial rerun of the paper's third mesh level (mesh factor 2 of the
tutorial's verification mesh study, 1 920 fluid cells, dt = 0.01 s, the
finite-strain neo-Hookean solid of the paper's definition) to t = 3 s, when
the flow is steady.

The fluid is coloured by the speed, with streamlines; the half-cylinder is
drawn deformed (true scale), its undeformed outline dashed.

Usage: panel_blob.py CASE OUT.pdf [--time T]"""
import argparse
import sys

import common
import fsi2d

TITLE = "Blob in treacle"
XLIM, YLIM = (-0.03, 5.0), (-0.72, 1.03)
ASPECT = (YLIM[1] - YLIM[0]) / (XLIM[1] - XLIM[0])


def load(case, time=3.0):
    return fsi2d.load(case, time)


def draw(fig, ax, cax, d, letter="i"):
    sp = fsi2d.speed(d)
    vmax = common.nice_ceil(sp.max())
    fsi2d.draw(fig, ax, cax, d, sp, common.Normalize(0, vmax), "viridis",
               r"$|\mathbf{u}|$ (m/s)", XLIM, YLIM, stream_density=0.9)
    common.scalebar(ax, 3.6, -0.64, 1.0, "1 m", side="left")
    common.title(ax, letter, f"{TITLE}, $t = {d['t']:g}$ s")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("case")
    ap.add_argument("out")
    ap.add_argument("--time", type=float, default=3.0)
    a = ap.parse_args()
    d = load(a.case, a.time)
    print(f"t = {d['t']}, max speed {fsi2d.speed(d).max():.3f}, "
          f"max D_x {d['D'][:, 0].max():.5f}")
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
