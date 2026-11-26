"""Gallery panel: three-dimensional pressure-pulse tube.

Data: a rerun of this paper's level-1 case (corrected meshes, Robin-Neumann
coupling, backward time scheme, 8 ranks; case robin_backward_mesh1 of the
3dTube mesh study) to t = 5 ms, with fields written every 1 ms and
reconstructed (reconstructPar -region fluid/solid).

The quarter model is cut by its symmetry plane y = 0 (patch symmetry-y) and
mirrored about the axis, so the panel shows a full meridional section: the
fluid pressure (rho p, rho = 1000 kg/m^3) and the wall displaced with
EXAGGERATION times its displacement, with the undeformed wall dashed.

Usage: panel_tube3d.py CASE OUT.pdf [--time T]"""
import argparse
import sys

import numpy as np
from matplotlib.tri import Triangulation

import common
import foamio

TITLE = "3-D pressure-pulse tube"
RHO = 1000.0
EXAGGERATION = 10.0
MM = 1e3
XLIM, YLIM = (-0.4, 50.4), (-10.2, 8.0)
ASPECT = (YLIM[1] - YLIM[0]) / (XLIM[1] - XLIM[0])


def load(case, time=0.005):
    _, tf, df, t = foamio.patch(case, "fluid", "symmetry-y", time, ("p", "U"),
                                axes=(2, 0))
    _, ts, ds, _ = foamio.patch(case, "solid", "symmetry-y", t, ("pointD",),
                                axes=(2, 0))
    return dict(t=t, fluid=tf, p=RHO * df["p"], solid=ts, D=ds["pointD"])


def draw(fig, ax, cax, d, letter="e"):
    vmax = common.nice_ceil(d["p"].max() / 1e3)
    norm = common.Normalize(0.0, vmax)
    tf, ts, D = d["fluid"], d["solid"], d["D"]
    ex = EXAGGERATION
    for sgn in (1, -1):
        # the fluid mesh has moved with the (true-scale) wall: stretch it
        # radially so that it meets the exaggerated wall
        r = tf.y
        rw = 0.005
        # wall displacement at the fluid points, from the nearest wall point
        zf = tf.x
        inner = ts.y < rw + 1e-6 + np.abs(D[:, 0]).max()
        zi = ts.x[inner]
        o = np.argsort(zi)
        ur = np.interp(zf, zi[o], D[inner, 0][o])
        rr = r * (1 + (ex - 1) * ur / (rw + ur))
        tri = Triangulation(zf * MM, sgn * rr * MM, tf.triangles)
        pc = ax.tripcolor(tri, d["p"] / 1e3, cmap="viridis", norm=norm,
                          shading="gouraud", rasterized=True, zorder=1)
        trs = Triangulation((ts.x + ex * D[:, 2]) * MM,
                            sgn * (ts.y + ex * D[:, 0]) * MM, ts.triangles)
        common.solid_patch(ax, trs)
        for rad in (0.005, 0.006):
            ax.plot([0, 50], [sgn * rad * MM] * 2, **common.UNDEFORMED)
    ax.plot([0, 50], [0, 0], color="0.4", lw=0.4, ls=":", zorder=3)
    ax.set_xlim(*XLIM)
    ax.set_ylim(*YLIM)
    common.bare(ax)
    common.scalebar(ax, 0.0, -8.6, 10.0, "10 mm", above=False)
    common.title(ax, letter, f"3-D tube, $t = {1e3 * d['t']:.0f}$ ms, "
                 rf"wall displ. $\times{EXAGGERATION:g}$")
    common.colorbar(fig, pc, cax, r"$p$ (kPa)")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("case")
    ap.add_argument("out")
    ap.add_argument("--time", type=float, default=0.005)
    a = ap.parse_args()
    d = load(a.case, a.time)
    print(f"t = {d['t']}, p {d['p'].min():.1f} .. {d['p'].max():.1f} Pa, "
          f"max |D| {np.abs(d['D']).max(0)}")
    common.use_style()
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(5.0, 1.6))
    ax = fig.add_axes([0.03, 0.05, 0.84, 0.82])
    cax = fig.add_axes([0.89, 0.2, 0.02, 0.6])
    draw(fig, ax, cax, d)
    fig.savefig(a.out, dpi=400)
    print("wrote", a.out)


if __name__ == "__main__":
    sys.exit(main())
