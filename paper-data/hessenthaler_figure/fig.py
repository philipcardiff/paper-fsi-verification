"""Hessenthaler flap, 3-D perspective view (fig:hessenthaler_3d): the rendered
image of render.py with labels, axis triad and |U| colour bar in the paper style.

Usage: fig.py [--render R.png] [--out OUT.pdf] [--png PREVIEW.png]
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
ap.add_argument("--out", default="hessenthaler-3d.pdf")
ap.add_argument("--png")
args = ap.parse_args()
sys.path.insert(0, str(HERE.parent / "figstyle"))
import figstyle  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib as mpl  # noqa: E402

figstyle.use()
img = plt.imread(args.render)[..., :3]
meta = json.load(open(Path(args.render).with_suffix(".json")))
px = {k: np.array(p) for k, p in meta["px"].items()}

# crop the white margin (keep a small border)
nonwhite = np.where(img.min(axis=2) < 0.985)
r0, r1 = nonwhite[0].min(), nonwhite[0].max()
c0, c1 = nonwhite[1].min(), nonwhite[1].max()
pad = 20
r0, c0 = max(r0 - pad, 0), max(c0 - pad, 0)
r1, c1 = min(r1 + pad, img.shape[0]), min(c1 + pad, img.shape[1])
img = img[r0:r1, c0:c1]
for k in px:
    px[k] = px[k] - [c0, r0]
h, w = img.shape[:2]

W = figstyle.TEXTWIDTH
left, right_img = 0.11, 0.80        # image spans these fractions of the width
H = W*(right_img - left)*h/w + 0.2
fig = plt.figure(figsize=(W, H))
ax = fig.add_axes([left, 0.1/H, right_img - left, (H - 0.2)/H])
ax.imshow(img, interpolation="lanczos")
ax.set_xlim(-0.08*w, w)
ax.set_ylim(h, 0)
ax.set_position([left - 0.08*(right_img - left), 0.1/H, 1.08*(right_img - left), (H - 0.2)/H])
ax.set_axis_off()

INLET, OUTLET = "#0072B2", "#D55E00"
arrow = dict(arrowstyle="-", lw=0.5, color="k", shrinkA=1, shrinkB=1)


def lab(text, at, dx, dy, absolute=False, **kw):
    """Label pointing at anchor 'at'; text offset by (dx, dy) image fractions,
    or placed at (dx, dy) image fractions if absolute."""
    x, y = px[at]
    kw.setdefault("fontsize", 7)
    xt = (dx*w, dy*h) if absolute else (x + dx*w, y + dy*h)
    ax.annotate(text, xy=(x, y), xytext=xt, arrowprops=arrow,
                ha=kw.pop("ha", "center"), va=kw.pop("va", "center"), **kw)


lab("upper inlet", "upperInlet", -0.06, -0.16, color=INLET)
lab("lower inlet", "lowerInlet", 0.07, 0.12, color=INLET, ha="left")
lab("deformed flap, $t = 45$ s", "flapTip", 0.04, -0.33, ha="left")
lab("undeformed flap", "flap0Tip", 0.29, 0.95, absolute=True, va="top")
lab(r"mid-plane $x = 0$, coloured by $|\mathbf{u}|$", "sliceMid", 0.62, 0.95, absolute=True, va="top")
lab(f"outlet pipe, cut at $z = {meta['zcut']:.0f}$ mm\n(continues to $z = 339.5$ mm)", "cutTop",
    -0.04, -0.12, ha="right", va="bottom")

# axis triad: projected directions of 15 mm along x, y, z, drawn to a common scale
o = np.array([0.86*w, 0.86*h])
ds = {k: px["triad_" + k] - px["triad_o"] for k in ("ex", "ey", "ez")}
sc = 0.06*w/max(np.linalg.norm(d) for d in ds.values())
for k, s in (("ex", "x"), ("ey", "y"), ("ez", "z")):
    d = ds[k]*sc
    ax.annotate("", xy=o + d, xytext=o, arrowprops=dict(arrowstyle="-|>", lw=0.7, color="k",
                                                        shrinkA=0, shrinkB=0, mutation_scale=6))
    ax.text(*(o + d + 0.022*w*d/np.linalg.norm(d)), f"${s}$", fontsize=7, ha="center", va="center")
g0 = o + [0.09*w, -0.06*h]
dg = -ds["ey"]*sc
ax.annotate("", xy=g0 + dg, xytext=g0, arrowprops=dict(arrowstyle="-|>", lw=0.8, color="k",
                                                       shrinkA=0, shrinkB=0, mutation_scale=7))
ax.text(*(g0 + 0.5*dg + [0.012*w, 0]), r"$\mathbf{g}$", fontsize=8, ha="left", va="center")

cax = fig.add_axes([0.86, 0.25, 0.016, 0.55])
cb = fig.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(0, meta["umax"]), cmap="viridis"),
                  cax=cax)
cb.set_label(r"$|\mathbf{u}|$ (m/s)")
cb.ax.tick_params(direction="out", right=True)
fig.savefig(args.out, dpi=400)
if args.png:
    fig.savefig(args.png, dpi=200)
