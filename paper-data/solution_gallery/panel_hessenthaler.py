"""Gallery panel: Hessenthaler flap (Phase I), velocity near the mid-plane.

Data: this paper's partitioned (IQN-ILS) Phase-I run on the rF1 fluid mesh,
time t = 45 s (end of the averaging window), reconstructed with
reconstructPar -region fluid -time 45 -fields '(U p)' and
reconstructPar -region solid -time 45 -fields '(D pointD)'.

The fluid and the flap are cut by the plane x = X0 (default 0.3 mm, the plane
of the mesh figure, which avoids cutting exactly along cell faces); the panel
shows the in-plane speed with streamlines of the in-plane velocity, and the
flap at true scale with its undeformed outline dashed. Coordinates in mm:
z (flow direction) horizontal, y vertical.

Usage: panel_hessenthaler.py CASE OUT.pdf [--time T] [--x0 X]"""
import argparse
import sys

import common
import fsi2d

TITLE = "Hessenthaler flap"
MM = 1e3
XLIM, YLIM = (-31.0, 131.0), (-45.5, 39.5)
ASPECT = (YLIM[1] - YLIM[0]) / (XLIM[1] - XLIM[0])


def load(case, time=45.0, x0=3e-4):
    return fsi2d.load(case, time, origin=(x0, 0.0, 0.1), normal=(1, 0, 0),
                      axes=(2, 1))


def draw(fig, ax, cax, d, letter="h"):
    sp = fsi2d.speed(d)
    vmax = common.nice_ceil(sp.max())
    fsi2d.draw(fig, ax, cax, d, sp, common.Normalize(0, vmax), "viridis",
               r"$|\mathbf{u}_\parallel|$ (m/s)", XLIM, YLIM, stream_density=0.8,
               scale=MM)
    common.scalebar(ax, 110.0, -43.0, 20.0, "20 mm", side="left")
    common.title(ax, letter, f"{TITLE}, $t = {d['t']:g}$ s")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("case")
    ap.add_argument("out")
    ap.add_argument("--time", type=float, default=45.0)
    ap.add_argument("--x0", type=float, default=3e-4)
    a = ap.parse_args()
    d = load(a.case, a.time, a.x0)
    print(f"t = {d['t']}, max in-plane speed {fsi2d.speed(d).max():.3f} m/s, "
          f"max D_y {d['D'][:, 1].max() * MM:.2f} mm")
    common.use_style()
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(2.8, 1.6))
    ax = fig.add_axes([0.01, 0.03, 0.80, 0.82])
    cax = fig.add_axes([0.84, 0.15, 0.03, 0.6])
    draw(fig, ax, cax, d)
    fig.savefig(a.out, dpi=400)
    print("wrote", a.out)


if __name__ == "__main__":
    sys.exit(main())
