"""3-D pressure-pulse tube, perspective view of the quarter model at t = 5 ms
(fig:3dtube_3d): the image of render.py with labels, axis triad and the two
colour bars (p, u_r) in the paper style.

Usage: fig.py [--render R.png] [--out OUT.pdf] [--png PREVIEW.png] [--time-label TXT]
  R.png: image written by render.py; the anchor file R.json must sit next to it
         (default: render.png in this script's directory)
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--render", default=str(HERE / "render.png"))
ap.add_argument("--out", default="3dtube-3d.pdf")
ap.add_argument("--png")
ap.add_argument("--time-label", default=r"$t = 5$ ms")
args = ap.parse_args()
sys.path.insert(0, str(HERE.parent / "figstyle"))
import figstyle  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib as mpl  # noqa: E402

figstyle.use()
img = plt.imread(args.render)[..., :3]
meta = json.load(open(Path(args.render).with_suffix(".json")))
px = {k: np.array(p) for k, p in meta["px"].items()}

nonwhite = np.where(img.min(axis=2) < 0.985)
pad = 30
r0, r1 = max(nonwhite[0].min() - pad, 0), min(nonwhite[0].max() + pad, img.shape[0])
c0, c1 = max(nonwhite[1].min() - pad, 0), min(nonwhite[1].max() + pad, img.shape[1])
img = img[r0:r1, c0:c1]
for k in px:
    px[k] = px[k] - [c0, r0]
h, w = img.shape[:2]

W = figstyle.TEXTWIDTH
fx0, fx1 = 0.13, 0.74                      # image extent (figure fractions)
H = W*(fx1 - fx0)*h/w
fig = plt.figure(figsize=(W, H))
ax = fig.add_axes([fx0, 0.0, fx1 - fx0, 1.0])
ax.imshow(img, interpolation="lanczos")
ax.set_xlim(0, w)
ax.set_ylim(h, 0)
ax.set_axis_off()
ax.set_clip_on(False)

arrow = dict(arrowstyle="-", lw=0.5, color="k", shrinkA=1, shrinkB=1)


def lab(text, at, fx, fy, **kw):
    """Label at image fractions (fx, fy) pointing at anchor 'at'."""
    kw.setdefault("fontsize", 7)
    ax.annotate(text, xy=px[at], xytext=(fx*w, fy*h), arrowprops=arrow, annotation_clip=False,
                ha=kw.pop("ha", "center"), va=kw.pop("va", "center"), **kw)


lab("inlet", "inlet", -0.08, 0.80, ha="right")
lab("outlet", "outlet", 1.04, 0.10, ha="left")
lab(r"fluid pressure $p$ on the symmetry plane", "plane", 0.62, 0.92, ha="left")
lab(f"wall, radial displacement $u_r$ (shape $\\times{meta['scale']:g}$)", "urPeak", 0.33, 0.07, ha="center")
ax.text(-0.13*w, 0.55*h, args.time_label, fontsize=8, ha="left", va="center")

# axis triad: projected directions of 3 mm along x, y, z, to a common scale
o = np.array([-0.10*w, 0.30*h])
ds = {k: px["triad_" + k] - px["triad_o"] for k in ("ex", "ey", "ez")}
sc = 0.06*w/max(np.linalg.norm(d) for d in ds.values())
for k, s in (("ex", "x"), ("ey", "y"), ("ez", "z")):
    d = ds[k]*sc
    ax.annotate("", xy=o + d, xytext=o, annotation_clip=False,
                arrowprops=dict(arrowstyle="-|>", lw=0.7, color="k", shrinkA=0, shrinkB=0,
                                mutation_scale=6))
    ha, va = {"x": ("center", "bottom"), "y": ("right", "center"), "z": ("left", "center")}[s]
    ax.text(*(o + d), f" ${s}$ " if s != "x" else f"${s}$", fontsize=7, ha=ha, va=va, clip_on=False)


def cbar(x, cmap, lo, hi, label, ticks=None):
    cax = fig.add_axes([x, 0.18, 0.014, 0.64])
    cb = fig.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(lo, hi), cmap=cmap),
                      cax=cax, extend="min", ticks=ticks)
    cb.set_label(label)
    cb.ax.tick_params(direction="out", right=True)
    return cb


cbar(0.80, "viridis", meta["pmin"], meta["pmax"], r"$p$ (Pa)")
cbar(0.92, "Oranges", 0.0, meta["urmax"], r"$u_r$ (mm)")
fig.savefig(args.out, dpi=400)
if args.png:
    fig.savefig(args.png, dpi=200)
