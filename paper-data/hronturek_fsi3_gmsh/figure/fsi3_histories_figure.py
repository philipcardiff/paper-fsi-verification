#!/usr/bin/env python3
"""Turek-Hron FSI3, Gmsh family: tip displacement and force histories (L1-L3).

Reads the CSV files written by extract_histories.py and draws, at the full
text width of the paper,

  (a) u_y(t), (b) u_x(t), (c) lift, (d) drag over the five-period QoI window,
  (e)-(g) zooms of one u_y maximum, one u_x minimum and one lift maximum,

with the references as horizontal lines at mean +/- amplitude: Failer and
Richter (2020), level 6, and the Featflow level-4 table. The references give
only the amplitude of u_y and of the lift, so their lines are centred on the
L3 (max+min)/2 of the window.

Each level is shifted in time so that its last upward u_y crossing before
t = 6 s (the end of its QoI window) is at t' = 0.

    python3 fsi3_histories_figure.py --data-dir data \
        --out ../../../figures/fsi3_gmsh_histories.pdf [--png preview.png]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib import font_manager  # noqa: E402

HERE = Path(__file__).resolve().parent
DEFAULT_STYLE = HERE.parents[1] / "figstyle" / "paper.mplstyle"
sys.path.insert(0, str(HERE.parents[1] / "figstyle"))
import figstyle  # noqa: E402  (paper widths and fonts)

# Failer & Richter (2020) level 6 and Featflow level 4: (mean, amplitude);
# None = mean not given (centred on the L3 mean). Units mm and N/m.
REF = {
    "uy": {"FR": (None, 35.03), "FF": (None, 34.99)},
    "ux": {"FR": (-2.890, 2.729), "FF": (-2.88, 2.72)},
    "lift": {"FR": (None, 156.6), "FF": (None, 153.91)},
    "drag": {"FR": (460.5, 27.91), "FF": (460.5, 27.74)},
}
LEVELS = ("L1", "L2", "L3")
STYLE = {"L1": dict(color="#0072B2", ls=(0, (5, 1.5)), lw=0.9),
         "L2": dict(color="#D55E00", ls=(0, (5, 1.2, 1, 1.2)), lw=0.9),
         "L3": dict(color="#000000", ls="-", lw=0.9)}
REF_STYLE = {"FR": dict(color="0.55", ls="-", lw=1.6, alpha=0.6, zorder=0),
             "FF": dict(color="0.35", ls=(0, (1, 1.2)), lw=0.8, zorder=0)}
LABEL = {"uy": r"$u_y$ (mm)", "ux": r"$u_x$ (mm)", "lift": r"Lift (N/m)", "drag": r"Drag (N/m)"}
SCALE = {"uy": 1e3, "ux": 1e3, "lift": 1.0, "drag": 1.0}
COL = {"ux": 1, "uy": 2, "drag": 3, "lift": 4}


def use_style(path: Path) -> None:
    figstyle.use()            # registers CMU Serif from the local TeX installation
    plt.style.use(str(path))
    names = {f.name for f in font_manager.fontManager.ttflist}
    if not names & {"Computer Modern Roman", "CMU Serif"}:
        # matplotlib ships Computer Modern as cmr10
        plt.rcParams["font.serif"] = ["cmr10"] + plt.rcParams["font.serif"]
        plt.rcParams["axes.formatter.use_mathtext"] = True


def crossings(t: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Upward crossings of the mean (as ht_coupled_analysis.py)."""
    mean = y.mean()
    band = 0.25 * (y.max() - y.min())
    out, armed = [], y[0] < mean - band
    for k in range(1, len(y)):
        if y[k] < mean - band:
            armed = True
        elif armed and y[k] >= mean:
            a, b = y[k - 1] - mean, y[k] - mean
            out.append(t[k - 1] + a / (a - b) * (t[k] - t[k - 1]))
            armed = False
    return np.array(out)


def load(data_dir: Path) -> dict:
    out = {}
    for lev in LEVELS:
        d = np.loadtxt(data_dir / f"fsi3_gmsh_{lev}_histories.csv", delimiter=",", skiprows=1)
        cr = crossings(d[:, 0], d[:, 2])
        bounds = cr[-6:]                      # five full periods: the QoI window
        out[lev] = {"t": d[:, 0] - bounds[-1], "bounds": bounds - bounds[-1], "t0": bounds[-1],
                    **{q: d[:, c] * SCALE[q] for q, c in COL.items()}}
    return out


def l3_mid(h: dict, q: str) -> float:
    """(max+min)/2 averaged over the five periods of the L3 window."""
    t, y, b = h["t"], h[q], h["bounds"]
    return float(np.mean([0.5 * (y[(t >= a) & (t <= c)].max() + y[(t >= a) & (t <= c)].min())
                          for a, c in zip(b[:-1], b[1:])]))


def ref_lines(ax, q: str, h3: dict, which=("max", "min")) -> None:
    for key in ("FR", "FF"):
        mean, amp = REF[q][key]
        if mean is None:
            mean = l3_mid(h3, q)
        if "max" in which:
            ax.axhline(mean + amp, **REF_STYLE[key])
        if "min" in which:
            ax.axhline(mean - amp, **REF_STYLE[key])


def panel_label(ax, s: str) -> None:
    ax.text(0.015, 0.97, s, transform=ax.transAxes, ha="left", va="top",
            bbox=dict(facecolor="white", edgecolor="none", pad=0.6))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--data-dir", type=Path, default=HERE / "data")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--png", type=Path, help="optional PNG preview")
    ap.add_argument("--style", type=Path, default=DEFAULT_STYLE)
    a = ap.parse_args()
    use_style(a.style)
    h = load(a.data_dir)
    h3 = h["L3"]
    for lev in LEVELS:
        print(f"{lev}: shift t0 = {h[lev]['t0']:.4f} s, window {h[lev]['bounds'][0]:.4f}-0 s")

    fig = plt.figure(figsize=(figstyle.TEXTWIDTH, 5.0))
    gs = fig.add_gridspec(3, 6, height_ratios=[1, 1, 1.05], hspace=0.42, wspace=1.6,
                          left=0.09, right=0.99, top=0.95, bottom=0.08)
    main_axes = {"uy": fig.add_subplot(gs[0, :3]), "ux": fig.add_subplot(gs[0, 3:]),
                 "lift": fig.add_subplot(gs[1, :3]), "drag": fig.add_subplot(gs[1, 3:])}
    t_lo = h3["bounds"][0]
    for (q, ax), lab in zip(main_axes.items(), "abcd"):
        ref_lines(ax, q, h3)
        for lev in LEVELS:
            m = h[lev]["t"] >= t_lo - 0.02
            ax.plot(h[lev]["t"][m], h[lev][q][m], label=lev, **STYLE[lev])
        ax.set_xlim(t_lo, 0.0)
        ax.set_ylabel(LABEL[q])
        panel_label(ax, f"({lab})")
        if q in ("lift", "drag"):
            ax.set_xlabel(r"$t - t_0$ (s)")
    yl = main_axes["drag"].get_ylim()
    main_axes["drag"].set_ylim(yl[0], yl[1] + 0.25 * (yl[1] - yl[0]))
    main_axes["ux"].set_ylim(-6.4, 1.3)
    main_axes["uy"].set_ylim(-45, 52)
    main_axes["lift"].set_ylim(-210, 260)

    # zooms around the extremum of the last period of the window
    zooms = [("uy", "max", "(e)"), ("ux", "min", "(f)"), ("lift", "max", "(g)")]
    zoom_axes = [fig.add_subplot(gs[2, 2 * k:2 * k + 2]) for k in range(3)]
    period = -h3["bounds"][-2]
    for (q, kind, lab), ax in zip(zooms, zoom_axes):
        half = {"uy": 0.012, "ux": 0.01, "lift": 0.015}[q]
        m3 = (h3["t"] >= -period + half) & (h3["t"] <= -half)
        k3 = np.argmax(h3[q][m3]) if kind == "max" else np.argmin(h3[q][m3])
        tc = h3["t"][m3][k3]
        ref_lines(ax, q, h3, which=(kind,))
        ys = []
        for lev in LEVELS:
            m = np.abs(h[lev]["t"] - tc) <= half
            ax.plot(h[lev]["t"][m], h[lev][q][m], **STYLE[lev])
            ys.append(h[lev][q][m])
        ys = np.concatenate(ys)
        refs = []
        for key in ("FR", "FF"):
            mean, amp = REF[q][key]
            mean = l3_mid(h3, q) if mean is None else mean
            refs.append(mean + amp if kind == "max" else mean - amp)
        ext = ys.max() if kind == "max" else ys.min()
        lo, hi = min(min(refs), ext), max(max(refs), ext)
        pad = 0.35 * (hi - lo)
        if kind == "max":
            ax.set_ylim(lo - 1.6 * pad, hi + 0.6 * pad)
        else:
            ax.set_ylim(lo - 0.6 * pad, hi + 1.6 * pad)
        ax.set_xlim(tc - half, tc + half)
        ax.set_xlabel(r"$t - t_0$ (s)")
        ax.set_ylabel(LABEL[q])
        ax.locator_params(axis="x", nbins=3)
        panel_label(ax, lab)
        # mark the zoom window in the main panel
        main_axes[q].axvspan(tc - half, tc + half, color="0.85", lw=0, zorder=-1)

    handles = [plt.Line2D([], [], label=lev, **STYLE[lev]) for lev in LEVELS]
    handles += [plt.Line2D([], [], label="Failer & Richter (2020), level 6", **REF_STYLE["FR"]),
                plt.Line2D([], [], label="Featflow, level 4", **REF_STYLE["FF"])]
    fig.legend(handles=handles, loc="upper center", ncol=5, bbox_to_anchor=(0.54, 1.005),
               handlelength=2.6, columnspacing=1.2)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(a.out)
    if a.png:
        fig.savefig(a.png, dpi=200)
    print(f"wrote {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
