"""Gallery panel: collapsible channel at the first trough.

Data: a rerun of this paper's fluid level F2 (12 800 fluid cells, solid
160 x 8 high-order, Robin-Neumann coupling, dt = 0.00625 s; the case of the
fluid-mesh study, configured with the tutorial's verification driver) to the
first trough, t = 0.46875 s, written at that time only. The wall-midpoint
displacement reproduces the study's F2 history (paper-data/
collapsible_mesh_study/F2.dat) to 1e-7 m.

The fluid is coloured by the speed, with streamlines; the elastic wall is
drawn deformed (true scale), its undeformed position dashed. Only the part
of the channel around the elastic section is shown.

Usage: panel_collapsible.py CASE OUT.pdf [--time T]"""
import argparse
import sys

import common
import fsi2d

TITLE = "Collapsible channel"
XLIM, YLIM = (3.0, 19.0), (-0.12, 1.1)
ASPECT = (YLIM[1] - YLIM[0]) / (XLIM[1] - XLIM[0])


def load(case, time=None):
    return fsi2d.load(case, time)


def draw(fig, ax, cax, d, letter="c"):
    sp = fsi2d.speed(d)
    vmax = common.nice_ceil(sp.max())
    fsi2d.draw(fig, ax, cax, d, sp, common.Normalize(0, vmax), "viridis",
               r"$|\mathbf{u}|$ (m/s)", XLIM, YLIM, stream_density=(2.2, 0.35),
               stream_n=400, cb_ticks=[0, 1, 2])
    common.scalebar(ax, 17.0, -0.1, 2.0, "2 m", above=False)
    common.title(ax, letter, f"{TITLE}, $t = {d['t']:.3f}$ s (first trough)")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("case")
    ap.add_argument("out")
    ap.add_argument("--time", type=float)
    a = ap.parse_args()
    d = load(a.case, a.time)
    print(f"t = {d['t']}, max speed {fsi2d.speed(d).max():.3f}, "
          f"min D_y {d['D'][:, 1].min():.5f}")
    common.use_style()
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(5.0, 0.75))
    ax = fig.add_axes([0.01, 0.12, 0.86, 0.7])
    cax = fig.add_axes([0.89, 0.15, 0.02, 0.6])
    draw(fig, ax, cax, d)
    fig.savefig(a.out, dpi=400)
    print("wrote", a.out)


if __name__ == "__main__":
    sys.exit(main())
