"""Assemble the solution gallery (fig:solution_gallery of the verification
paper) from the panel scripts.

Each panel is drawn by the draw() function of its panel_*.py module on its
own axes, at true aspect ratio (the axes box has the aspect of the panel's
data window), with its colour bar to the right. The rows are:

    (a) ring added mass            (b) Mok lid-driven cavity
    (c) Womersley tube
    (d) collapsible channel
    (e) Turek-Hron FSI3
    (f) 3-D pressure-pulse tube    (g) flexible-bottom cavity
    (h) Hessenthaler flap          (i) blob in treacle

Usage:
  make_gallery.py OUT.pdf --womersley PROPS --collapsible CASE --fsi3 CASE
      --tube CASE --mok CASE --cavity CASE --blob CASE
      [--hessenthaler CASE] [--png OUT.png] [--width W]
See README.md for the data behind each argument."""
import argparse

import matplotlib.pyplot as plt

import common
import panel_blob
import panel_cavity
import panel_collapsible
import panel_fsi3
import panel_hessenthaler
import panel_mok
import panel_ring
import panel_tube3d
import panel_womersley

CB_GAP, CB_W, CB_SPACE = 0.05, 0.07, 0.50   # inches: gap, bar, bar + labels
TITLE_H, ROW_GAP = 0.18, 0.035             # inches
SQUARE_MAX = 1.4                          # inches, height cap of 1:1 panels


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out")
    ap.add_argument("--womersley", required=True,
                    help="womersleyTube constant/womersleyProperties")
    ap.add_argument("--collapsible", required=True)
    ap.add_argument("--fsi3", required=True)
    ap.add_argument("--tube", required=True)
    ap.add_argument("--mok", required=True)
    ap.add_argument("--cavity", required=True)
    ap.add_argument("--blob", required=True)
    ap.add_argument("--hessenthaler")
    ap.add_argument("--png")
    ap.add_argument("--width", type=float, default=372 / 72,
                    help="figure width in inches (default: text width, 372 pt)")
    a = ap.parse_args()

    print("loading ...")
    P = {
        "ring": (panel_ring, panel_ring.load()),
        "wom": (panel_womersley, panel_womersley.load(a.womersley)),
        "coll": (panel_collapsible, panel_collapsible.load(a.collapsible)),
        "fsi3": (panel_fsi3, panel_fsi3.load(a.fsi3)),
        "tube": (panel_tube3d, panel_tube3d.load(a.tube)),
        "mok": (panel_mok, panel_mok.load(a.mok)),
        "cav": (panel_cavity, panel_cavity.load(a.cavity)),
        "blob": (panel_blob, panel_blob.load(a.blob)),
    }
    if a.hessenthaler:
        P["hv"] = (panel_hessenthaler, panel_hessenthaler.load(a.hessenthaler))
    rows = [["ring", "mok"], ["wom"], ["coll"], ["fsi3"], ["tube", "cav"],
            ["hv", "blob"] if "hv" in P else ["blob"]]

    W = a.width
    # panel sizes (inches): each row is split equally between its panels
    layout, y = [], 0.0
    for row in rows:
        cell = W / len(row)
        sizes = []
        for key in row:
            mod = P[key][0]
            pw = cell - CB_SPACE - CB_GAP - CB_W
            ph = pw * mod.ASPECT
            if mod.ASPECT >= 0.9 and ph > SQUARE_MAX:
                ph = SQUARE_MAX
                pw = ph / mod.ASPECT
            sizes.append((pw, ph))
        rh = max(s[1] for s in sizes)
        layout.append((row, cell, sizes, y, rh))
        y += TITLE_H + rh + ROW_GAP
    H = y - ROW_GAP + 0.02

    common.use_style()
    fig = plt.figure(figsize=(W, H))
    letters = iter("abcdefghi")
    for row, cell, sizes, y0, rh in layout:
        for k, (key, (pw, ph)) in enumerate(zip(row, sizes)):
            mod, data = P[key]
            x0 = k * cell
            top = y0 + TITLE_H
            # panels top-aligned in their row, so paired titles line up
            yb = top + ph
            ax = fig.add_axes([x0 / W, 1 - yb / H, pw / W, ph / H])
            # colour bar slightly inset from the panel top, clear of the title
            ch = max(ph - 0.06, 0.3)
            cb_bottom = top + 0.06 + ch
            cax = fig.add_axes([(x0 + pw + CB_GAP) / W, 1 - cb_bottom / H,
                                CB_W / W, ch / H])
            mod.draw(fig, ax, cax, data, next(letters))
            print(f"  ({key}) {pw:.2f} x {ph:.2f} in")
    fig.savefig(a.out, dpi=400)
    print(f"wrote {a.out} ({W:.2f} x {H:.2f} in)")
    if a.png:
        fig.savefig(a.png, dpi=200)
        print("wrote", a.png)


if __name__ == "__main__":
    main()
