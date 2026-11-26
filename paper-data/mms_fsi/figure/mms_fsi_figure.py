"""Manufactured FSI solution, pressure-error figure of the verification paper
(fig:mms_fsi): the fluid pressure error p - p_exact at the end time on the
deformed fluid mesh, for the partitioned (Dirichlet-Neumann, Aitken) and the
quasi-monolithic (QM) runs on the same mesh, with one colour scale. The case
itself (domains, deformation, flow, exact interface of both variants) is
drawn by mms_fsi_case_figure.py (fig:mms_fsi_case).

The pressure error is the pError field written by the case's error function
object (p - p_exact at the cell centres, not mean-corrected).

Usage:
  mms_fsi_figure.py AITKEN_RUN QM_RUN OUT.pdf [TIME] [LABEL]
    AITKEN_RUN, QM_RUN: case directories with TIME/fluid/{pError,polyMesh/points}
      and constant/fluid/polyMesh/{faces,owner,boundary}
    TIME: time directory of the error fields (default 0.4)
    LABEL: mesh label for the panel titles (default "n = 128")
Prints the range and RMS of each error field."""
import re
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import Normalize
from matplotlib.ticker import FormatStrFormatter

plt.rcParams.update({"font.size": 8, "mathtext.fontset": "cm",
                     "font.family": "serif", "pdf.fonttype": 42})

L, H, T = 2.0, 1.0, 0.5


# ---- minimal ASCII OpenFOAM readers ----------------------------------------
def _body(path):
    txt = Path(path).read_text()
    txt = re.sub(r"/\*.*?\*/", "", txt, flags=re.S)
    return re.sub(r"//[^\n]*", "", txt)


def _list_after(txt, start):
    """The 'N ( ... )' list starting at or after index start: (string, end)."""
    m = re.compile(r"(\d+)\s*\(").search(txt, start)
    i, depth = m.end(), 1
    j = i
    while depth:
        depth += {"(": 1, ")": -1}.get(txt[j], 0)
        j += 1
    return int(m.group(1)), txt[i:j - 1]


def read_points(path):
    txt = _body(path)
    n, body = _list_after(txt, txt.index("}") + 1)
    p = np.array([[float(v) for v in s.split()]
                  for s in re.findall(r"\(([^()]*)\)", body)])
    assert len(p) == n
    return p


def read_faces(path):
    txt = _body(path)
    n, body = _list_after(txt, txt.index("}") + 1)
    f = [[int(v) for v in s.split()] for s in re.findall(r"\d+\s*\(([^()]*)\)", body)]
    assert len(f) == n
    return f


def read_labels(path):
    txt = _body(path)
    n, body = _list_after(txt, txt.index("}") + 1)
    v = np.array(body.split(), dtype=int)
    assert len(v) == n
    return v


def read_scalar_field(path):
    txt = _body(path)
    k = txt.index("internalField")
    n, body = _list_after(txt, k)
    v = np.array(body.split(), dtype=float)
    assert len(v) == n
    return v


def empty_patch(path):
    """(startFace, nFaces) of the empty (front/back) patch."""
    txt = _body(path)
    for m in re.finditer(r"(\w+)\s*\{([^{}]*)\}", txt):
        if re.search(r"\btype\s+empty\s*;", m.group(2)):
            s = int(re.search(r"startFace\s+(\d+)", m.group(2)).group(1))
            n = int(re.search(r"nFaces\s+(\d+)", m.group(2)).group(1))
            return s, n
    raise RuntimeError("no empty patch in " + str(path))


def fluid_cells(run, time):
    """Cell polygons (x, y) of the 2-D fluid mesh at time, and pError."""
    pm = Path(run) / "constant" / "fluid" / "polyMesh"
    pts = read_points(Path(run) / time / "fluid" / "polyMesh" / "points")
    faces = read_faces(pm / "faces")
    owner = read_labels(pm / "owner")
    s, n = empty_patch(pm / "boundary")
    err = read_scalar_field(Path(run) / time / "fluid" / "pError")
    zmin = pts[:, 2].min()
    polys = [None] * len(err)
    for fi in range(s, s + n):
        f = faces[fi]
        if np.allclose(pts[f, 2], zmin):
            polys[owner[fi]] = pts[f, :2]
    assert all(p is not None for p in polys)
    return polys, err


def panel_error(ax, polys, err, norm, title):
    pc = PolyCollection(polys, array=err, cmap="RdBu_r", norm=norm,
                        edgecolors="face", linewidths=0.05, rasterized=True)
    ax.add_collection(pc)
    xy = np.vstack(polys)
    ax.plot([0, L], [0, 0], color="black", lw=0.6)
    ax.set_xlim(-0.02, L + 0.02)
    ax.set_ylim(-0.02, xy[:, 1].max() + 0.03)
    ax.set_aspect("equal")
    ax.set_xlabel("$x$")
    ax.set_xticks([0, 0.5, 1, 1.5, 2])
    ax.set_yticks([0, 0.5, 1])
    ax.set_title(title, fontsize=8)
    ax.text(0.03, 0.06, "max $|p - p_\\mathrm{ex}| = " + sci(abs(err).max()) + "$",
            transform=ax.transAxes, fontsize=7, ha="left", va="bottom")
    return pc


def sci(v):
    """2-digit LaTeX scientific notation."""
    m, e = f"{v:.1e}".split("e")
    return f"{m}\\times10^{{{int(e)}}}"


def main():
    if len(sys.argv) < 4:
        sys.exit(__doc__)
    run_a, run_q, out = sys.argv[1:4]
    time = sys.argv[4] if len(sys.argv) > 4 else "0.4"
    label = sys.argv[5] if len(sys.argv) > 5 else "n = 128"

    pa, ea = fluid_cells(run_a, time)
    pq, eq = fluid_cells(run_q, time)
    for nm, e in (("Aitken", ea), ("QM", eq)):
        print(f"{nm}: pError min {e.min():.3e} max {e.max():.3e} "
              f"rms {np.sqrt(np.mean(e**2)):.3e} (unweighted)")
    vmax = max(abs(ea).max(), abs(eq).max())
    norm = Normalize(vmin=-vmax, vmax=vmax)

    fig = plt.figure(figsize=(372 / 72, 1.75))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.0, 0.035], left=0.075,
                          right=0.875, bottom=0.2, top=0.88, wspace=0.08)
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1], sharey=ax1)
    cax = fig.add_subplot(gs[2])

    pc = panel_error(ax1, pa, ea, norm,
                     f"(a) partitioned (Aitken), ${label}$")
    panel_error(ax2, pq, eq, norm,
                f"(b) quasi-monolithic, ${label}$")
    ax1.set_ylabel("$y$", labelpad=1)
    plt.setp(ax2.get_yticklabels(), visible=False)
    ax1.annotate("outlet–wall corner", xy=(1.98, 0.2), xytext=(1.3, 0.3),
                 fontsize=7, ha="center", va="center",
                 arrowprops=dict(arrowstyle="->", lw=0.6))
    cb = fig.colorbar(pc, cax=cax)
    for ax in (ax1, ax2):
        ax.xaxis.set_major_formatter(FormatStrFormatter("%g"))
        ax.yaxis.set_major_formatter(FormatStrFormatter("%g"))
    cb.ax.tick_params(labelsize=7)
    cb.set_label("$p - p_\\mathrm{ex}$", fontsize=8)

    fig.savefig(out, dpi=400)
    print("wrote", out)


if __name__ == "__main__":
    main()
