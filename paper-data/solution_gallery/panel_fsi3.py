"""Gallery panel: Turek-Hron FSI3, vorticity and the deformed flag.

Data: this paper's partitioned run on the Gmsh level-2 mesh (case g2_np32 of
paper-data/hronturek_fsi3_gmsh), one written time in the periodic state,
reconstructed with reconstructPar -region fluid -fields '(U p)' and
reconstructPar -region solid -fields '(D pointD)'.

Usage: panel_fsi3.py CASE OUT.pdf [--time T] [--vmax W]"""
import argparse
import sys

import numpy as np

import common
import fsi2d

TITLE = r"Turek$-$Hron FSI3"
XLIM, YLIM = (-0.01, 2.51), (-0.03, 0.42)
ASPECT = (YLIM[1] - YLIM[0]) / (XLIM[1] - XLIM[0])


def load(case, time=None):
    return fsi2d.load(case, time, vorticity=True)


def draw(fig, ax, cax, d, letter="d", vmax=120.0):
    norm = common.sym_norm(vmax)
    fsi2d.draw(fig, ax, cax, d, d["vort"], norm, "RdBu_r",
               r"$\omega_z$ (1/s)", XLIM, YLIM, stream=False,
               cb_ticks=[-vmax, 0, vmax], extend="both")
    # cylinder
    th = np.linspace(0, 2 * np.pi, 200)
    ax.fill(0.2 + 0.05 * np.cos(th), 0.2 + 0.05 * np.sin(th), color="white",
            zorder=6)
    ax.plot(0.2 + 0.05 * np.cos(th), 0.2 + 0.05 * np.sin(th), color="0.1",
            lw=0.6, zorder=7)
    common.scalebar(ax, 2.1, 0.47, 0.4, "0.4 m", side="left")  # in the title row
    common.title(ax, letter, f"{TITLE}, $t = {d['t']:.2f}$ s")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("case")
    ap.add_argument("out")
    ap.add_argument("--time", type=float)
    ap.add_argument("--vmax", type=float, default=120.0)
    a = ap.parse_args()
    d = load(a.case, a.time)
    w = d["vort"]
    print(f"t = {d['t']}, vorticity range {w.min():.1f} .. {w.max():.1f}, "
          f"tip-region max D_y {d['D'][:, 1].max():.4f} m")
    common.use_style()
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(5.0, 1.15))
    ax = fig.add_axes([0.01, 0.12, 0.86, 0.75])
    cax = fig.add_axes([0.89, 0.2, 0.02, 0.6])
    draw(fig, ax, cax, d, vmax=a.vmax)
    fig.savefig(a.out, dpi=400)
    print("wrote", a.out)


if __name__ == "__main__":
    sys.exit(main())
