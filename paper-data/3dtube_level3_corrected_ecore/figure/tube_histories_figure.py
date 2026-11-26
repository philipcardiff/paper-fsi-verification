#!/usr/bin/env python3
"""3-D pressure-pulse tube: radial and axial wall displacement at point A.

Draws, at the full text width of the paper,

  (a) u_r(A, t) and (b) u_z(A, t), 0-20 ms,

for the present levels 1-3 (corrected meshes, Robin-Neumann, BDF2, dt halved
with the mesh) against the published point-A curves: Tukovic et al. (2018),
BDF2, radial only, to 10 ms; Lozovskiy et al. (2019) and Eken (2016),
implicit Euler, dt = 1e-4 s. Panel (a) has an inset around the radial peak.

Inputs (CSV, SI units, one row per time step / sample):
  --level LABEL=<csv> (repeat): columns time_s, radial_displacement_m,
      axial_displacement_m[, axis_pressure_Pa], as written by the 3dTube
      verification driver (histories/robin_backward_mesh<N>.csv);
  --reference-dir: Tukovic2018_fig25a_pointA.csv, Lozovskiy2019_fig3_pointA.csv,
      Eken2016_fig5.6_pointA.csv (digitised; '#' comment headers).

    python3 tube_histories_figure.py \
        --level "Level 1=data/robin_backward_mesh1.csv" \
        --level "Level 2=data/robin_backward_mesh2.csv" \
        --level "Level 3=../histories/robin_backward_mesh3.csv" \
        --reference-dir data/reference --out ../../../figures/3dtube_pointA_histories.pdf
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

LEVEL_STYLE = [dict(color="#0072B2", ls=(0, (5, 1.5)), lw=0.9),
               dict(color="#D55E00", ls=(0, (5, 1.2, 1, 1.2)), lw=0.9),
               dict(color="#000000", ls="-", lw=0.9)]
REFS = [  # file, label, style, radial column, axial column
    ("Tukovic2018_fig25a_pointA.csv", r"Tukovi$\mathrm{\acute{c}}$ et al. (2018), BDF2",
     dict(color="#009E73", ls=(0, (1, 1)), lw=1.3), 1, None),
    ("Lozovskiy2019_fig3_pointA.csv", "Lozovskiy et al. (2019), Euler",
     dict(color="#CC79A7", ls="-", lw=0.8, marker="s", markersize=2.6, mfc="white",
          markevery=(0, 15)), 1, 2),
    ("Eken2016_fig5.6_pointA.csv", r"Eken and Sahin (2016), Euler",
     dict(color="#E69F00", ls="-", lw=0.8, marker="^", markersize=2.8, mfc="white",
          markevery=(7, 15)), 1, 2),
]


def use_style(path: Path) -> None:
    figstyle.use()            # registers CMU Serif from the local TeX installation
    plt.style.use(str(path))
    names = {f.name for f in font_manager.fontManager.ttflist}
    if not names & {"Computer Modern Roman", "CMU Serif"}:
        # matplotlib ships Computer Modern as cmr10
        plt.rcParams["font.serif"] = ["cmr10"] + plt.rcParams["font.serif"]
        plt.rcParams["axes.formatter.use_mathtext"] = True


def read_csv(path: Path) -> np.ndarray:
    lines = [ln for ln in path.read_text().splitlines() if ln.strip() and not ln.startswith("#")]
    return np.array([[float(v) for v in ln.split(",")] for ln in lines[1:]])


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--level", action="append", required=True, help="LABEL=history csv")
    ap.add_argument("--reference-dir", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--png", type=Path, help="optional PNG preview")
    ap.add_argument("--style", type=Path, default=DEFAULT_STYLE)
    a = ap.parse_args()
    use_style(a.style)

    fig, axes = plt.subplots(1, 2, figsize=(figstyle.TEXTWIDTH, 2.5))
    fig.subplots_adjust(left=0.085, right=0.99, bottom=0.17, top=0.79, wspace=0.28)
    ax_r, ax_z = axes
    inset = ax_r.inset_axes([0.62, 0.52, 0.36, 0.44])
    handles = []
    for (fname, label, st, cr, cz) in REFS:
        d = read_csv(a.reference_dir / fname)
        t = d[:, 0] * 1e3
        h, = ax_r.plot(t, d[:, cr] * 1e3, label=label, **st)
        inset.plot(t, d[:, cr] * 1e3, **({**st, "markevery": 3} if "marker" in st else st))
        if cz is not None:
            ax_z.plot(t, d[:, cz] * 1e3, **st)
        handles.append(h)
    for k, spec in enumerate(a.level):
        label, _, path = spec.partition("=")
        d = read_csv(Path(path))
        t = d[:, 0] * 1e3
        st = LEVEL_STYLE[k % len(LEVEL_STYLE)]
        h, = ax_r.plot(t, d[:, 1] * 1e3, label=label, **st)
        inset.plot(t, d[:, 1] * 1e3, **st)
        ax_z.plot(t, d[:, 2] * 1e3, **st)
        k_max = int(np.argmax(d[:, 1]))
        print(f"{label}: u_r,max = {d[k_max, 1] * 1e3:.5f} mm at {t[k_max]:.3f} ms; "
              f"u_z,min (t < 6.5 ms) = {d[t < 6.5, 2].min() * 1e3:.5f} mm")
        handles.insert(k, h)
    for ax, lab, yl in ((ax_r, "(a)", r"$u_r(A)$ (mm)"), (ax_z, "(b)", r"$u_z(A)$ (mm)")):
        ax.set_xlim(0, 20)
        ax.set_xlabel(r"$t$ (ms)")
        ax.set_ylabel(yl)
        ax.axhline(0, color="0.7", lw=0.5, zorder=-1)
        ax.text(0.025, 0.97, lab, transform=ax.transAxes, ha="left", va="top")
    ax_r.set_ylim(-0.11, 0.235)
    ax_z.set_ylim(-0.11, 0.15)
    inset.set_xlim(6.4, 8.0)
    inset.set_ylim(0.112, 0.164)
    inset.tick_params(labelsize=7, length=2)
    inset.set_xticks([6.5, 7.0, 7.5])
    ax_r.indicate_inset_zoom(inset, edgecolor="0.5", lw=0.5)
    # row 1: the present levels, row 2: the published curves (legend fills by column)
    n = len(a.level)
    if n == len(REFS):
        handles = [h for pair in zip(handles[:n], handles[n:]) for h in pair]
    fig.legend(handles=handles, loc="upper center", ncol=3, bbox_to_anchor=(0.5, 1.01),
               handlelength=2.2, columnspacing=0.6, handletextpad=0.5)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(a.out)
    if a.png:
        fig.savefig(a.png, dpi=200)
    print(f"wrote {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
