"""Hessenthaler flap geometry and mesh figure (fig:hessenthaler_geom):
(a) the flow domain on the mid-plane x = 0 (section of the published fluid
surface, section_x0.npz); (b) the rF1 fluid mesh cut 0.3 mm from the
mid-plane (rF1_x0.vtk), with the undeformed flap and the measured Phase I
centreline (phaseI_centreline.csv). See README.md for the data.

Usage: fig.py [--data DIR] [--out OUT.pdf] [--png OUT.png]
  DIR: directory with the three data files (default: this script's directory)
  OUT: default hessenthaler-geometry.pdf in the current directory
"""
import argparse
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--data", default=str(HERE))
ap.add_argument("--out", default="hessenthaler-geometry.pdf")
ap.add_argument("--png")
args = ap.parse_args()
D = Path(args.data)
sys.path.insert(0, str(HERE.parent / "figstyle"))
import figstyle  # noqa: E402  (paper style: CM fonts, 8 pt at final size)
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.collections import LineCollection, PolyCollection  # noqa: E402

figstyle.use()
sec = np.load(D / "section_x0.npz")

# Mid-plane slice of the rF1 fluid mesh (legacy VTK polydata)
L = open(D / "rF1_x0.vtk").read().split("\n")
i = 0
while not L[i].startswith("POINTS"):
    i += 1
n = int(L[i].split()[1])
vals = []
i += 1
while len(vals) < 3*n:
    vals += L[i].split()
    i += 1
P = np.array(vals, float).reshape(-1, 3)*1e3
while not L[i].startswith("POLYGONS"):
    i += 1
m, tot = map(int, L[i].split()[1:3])
tok = []
i += 1
while len(tok) < tot:
    tok += L[i].split()
    i += 1
tok = list(map(int, tok))
polys = []
k = 0
for j in range(m):
    nv = tok[k]
    polys.append(P[tok[k + 1:k + 1 + nv]][:, [2, 1]])   # (z, y)
    k += 1 + nv


def flap(ax):
    ax.add_patch(plt.Rectangle((0, -1), 65, 2, fc="#e08a2c", ec="k",
                               lw=0.5, zorder=5))


INLET, OUTLET = "#0072B2", "#D55E00"      # Okabe-Ito blue, vermillion
cols = {"wall": "k", "upperInlet": INLET, "lowerInlet": INLET,
        "outlet": OUTLET, "interface": "k"}

# equal axes widths: height ratios from the aspect ratios of the two views
fig, (a, b) = plt.subplots(2, 1, figsize=(figstyle.TEXTWIDTH, 5.5),
                           gridspec_kw={"height_ratios": [108 / 400, 79 / 131],
                                        "hspace": 0.2},
                           subplot_kw={"anchor": "W"})
fig.subplots_adjust(left=0.1, right=0.97, top=0.97, bottom=0.08)
for k, s in sec.items():
    if len(s) == 0:
        continue
    a.add_collection(LineCollection(s[:, :, ::-1], colors=cols[k],
                                    linewidths=1.4 if k in ("upperInlet", "lowerInlet", "outlet") else 0.7))
flap(a)
a.set_xlim(-48, 352)
a.set_ylim(-54, 54)
a.set_aspect("equal")
a.set_xlabel(r"$z$ (mm)", labelpad=1)
a.set_ylabel(r"$y$ (mm)")
for y0 in (27.15, -27.15):
    a.annotate("", xy=(-28, y0), xytext=(-47, y0),
               arrowprops=dict(arrowstyle="-|>", lw=0.8, color=INLET))
a.text(-46, 40.5, "upper inlet", color=INLET, fontsize=7, va="bottom")
a.text(-46, -40.5, "lower inlet", color=INLET, fontsize=7, va="top")
a.text(344, 0, "outlet", color=OUTLET, rotation=90, va="center",
       ha="left", fontsize=7)
a.annotate("flap", xy=(45, 1.2), xytext=(62, 18), fontsize=7,
           arrowprops=dict(arrowstyle="-", lw=0.5))
a.annotate("", xy=(120, -32), xytext=(120, -12),
           arrowprops=dict(arrowstyle="-|>", lw=0.8))
a.text(124, -26, r"$\mathbf{g}$", fontsize=8)
a.annotate("", xy=(95, -38.1), xytext=(95, 38.1),
           arrowprops=dict(arrowstyle="<->", lw=0.5))
a.text(98, 12, r"$\varnothing\,76.2$", fontsize=7)
a.text(165, 18, "outlet pipe extended by 250 mm", fontsize=7)
a.text(-14, 27.15, r"$\varnothing\,21.9$", fontsize=6, ha="center",
       va="center")
figstyle.panel_label(a, "(a)")

b.add_collection(PolyCollection(polys, facecolors="none", edgecolors="0.35",
                                linewidths=0.15, rasterized=True))
for k, s in sec.items():
    if len(s):
        b.add_collection(LineCollection(s[:, :, ::-1], colors=cols[k],
                                        linewidths=0.9))
flap(b)
cl = np.loadtxt(D / "phaseI_centreline.csv", delimiter=",", comments="#",
                skiprows=4)
b.plot(cl[:, 2], cl[:, 1], "--", color=OUTLET, lw=1.2, zorder=6,
       label="measured centreline, Phase I")
b.legend(loc="lower right", frameon=True, framealpha=1.0, edgecolor="0.7", fancybox=False)
b.set_xlim(-31, 100)
b.set_ylim(-39.5, 39.5)
b.set_aspect("equal")
b.set_xlabel(r"$z$ (mm)", labelpad=1)
b.set_ylabel(r"$y$ (mm)")
figstyle.panel_label(b, "(b)")
for ax in (a, b):
    ax.tick_params(top=False, right=False)
fig.savefig(args.out, dpi=600)
if args.png:
    fig.savefig(args.png, dpi=200)
