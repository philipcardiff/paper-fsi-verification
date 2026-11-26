"""Drawing helpers shared by the panel scripts of the solution gallery.

Every panel is drawn by a function draw(ax, cax, data) of its module, where
ax is the panel's axes and cax the axes of its colour bar; the modules can
also be run on their own to write a one-panel preview."""
import os

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import Normalize, TwoSlopeNorm  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
STYLE = os.path.join(HERE, "..", "figstyle", "paper.mplstyle")

SOLID = "0.78"          # solid regions: uniform light grey
EDGE = "0.1"            # domain outlines
UNDEFORMED = dict(color="black", lw=0.6, ls=(0, (3, 2)), zorder=8)
STREAM = dict(color="white", linewidth=0.45, arrowsize=0.5, density=0.9,
              zorder=4)


def use_style():
    plt.style.use(STYLE)
    # The style asks for Computer Modern; where CMU is not installed, use
    # matplotlib's own cmr10 (with mathtext minus signs) before DejaVu Serif
    serif = list(plt.rcParams["font.serif"])
    if "cmr10" not in serif:
        serif.insert(max(0, len(serif) - 1), "cmr10")
    plt.rcParams.update({"xtick.top": False, "ytick.right": False,
                         "font.serif": serif,
                         "axes.formatter.use_mathtext": True,
                         "contour.negative_linestyle": "solid"})


def bare(ax):
    """No ticks or frame: the panels carry a scale bar instead."""
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_aspect("equal")


def scalebar(ax, x0, y0, length, label, color="black", above=True,
             side=None):
    """Horizontal scale bar of `length` (data units) starting at (x0, y0)."""
    ax.plot([x0, x0 + length], [y0, y0], color=color, lw=1.0,
            solid_capstyle="butt", zorder=12)
    for x in (x0, x0 + length):
        tick = 0.012 * np.diff(ax.get_ylim())[0]
        ax.plot([x, x], [y0 - tick, y0 + tick], color=color, lw=0.6, zorder=12)
    if side == "right":
        ax.annotate(label, (x0 + length, y0), xytext=(3, 0),
                    textcoords="offset points", ha="left", va="center",
                    fontsize=7, color=color, zorder=12)
    elif side == "left":
        ax.annotate(label, (x0, y0), xytext=(-3, 0),
                    textcoords="offset points", ha="right", va="center",
                    fontsize=7, color=color, zorder=12)
    else:
        ax.annotate(label, (x0 + 0.5 * length, y0),
                    xytext=(0, 2 if above else -2), textcoords="offset points",
                    ha="center", va="bottom" if above else "top", fontsize=7,
                    color=color, zorder=12)


def note(ax, text, x=0.99, y=0.03, ha="right", va="bottom", color="black"):
    """Small in-panel note (e.g. a display magnification), axes coords."""
    ax.text(x, y, text, transform=ax.transAxes, ha=ha, va=va, fontsize=7,
            color=color, zorder=13)


def title(ax, letter, text):
    ax.set_title(f"({letter}) {text}", loc="left", fontsize=8, pad=2.5)


def colorbar(fig, mappable, cax, label, ticks=None, fmt=None, extend="neither"):
    cb = fig.colorbar(mappable, cax=cax, ticks=ticks, format=fmt, extend=extend,
                      extendfrac=0.06)
    h_in = cax.get_position().height * fig.get_figheight()
    if h_in < 0.45:
        # a short bar (thin panel): the label goes above it, left-aligned
        cb.ax.set_title(label, loc="left", fontsize=7, pad=6)
    else:
        cb.set_label(label, fontsize=7, labelpad=2)
    cb.ax.tick_params(labelsize=7, length=2, width=0.5, pad=1.5)
    cb.outline.set_linewidth(0.5)
    return cb


def sym_norm(vmax):
    return TwoSlopeNorm(vcenter=0.0, vmin=-vmax, vmax=vmax)


def nice_ceil(v, steps=(1, 1.2, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10)):
    """Round v up to a 'nice' number for colour-bar limits."""
    if v <= 0:
        return 1.0
    e = 10 ** np.floor(np.log10(v))
    for s in steps:
        if s * e >= v * (1 - 1e-9):
            return s * e
    return 10 * e


def streamlines(ax, tri, vec, xlim, ylim, n=220, density=0.9, **kw):
    """Streamlines of the in-plane vector `vec` (N, 2) given at the points of
    `tri`, interpolated on a regular grid (masked outside the domain)."""
    from foamio import grid_interp
    xr, yr = xlim[1] - xlim[0], ylim[1] - ylim[0]
    nx = n if xr >= yr else max(20, int(n * xr / yr))
    ny = n if yr >= xr else max(20, int(n * yr / xr))
    x = np.linspace(xlim[0], xlim[1], nx)
    y = np.linspace(ylim[0], ylim[1], ny)
    u = grid_interp(tri, vec[:, 0], x, y)
    v = grid_interp(tri, vec[:, 1], x, y)
    opts = dict(STREAM)
    opts["density"] = density
    opts.update(kw)
    min_sep = opts.pop("arrow_sep", 0.07)
    before = set(map(id, ax.patches))
    ax.streamplot(x, y, u, v, broken_streamlines=True, **opts)
    # streamplot puts one arrow at the middle of every streamline, which
    # bunches them; keep only arrows at least min_sep (fraction of the
    # window diagonal) apart
    new = [p for p in ax.patches if id(p) not in before]
    diag = np.hypot(xr, yr)
    kept = []
    for p in new:
        pos = np.asarray(p._posA_posB[1])
        if all(np.hypot(*(pos - q)) > min_sep * diag for q in kept):
            kept.append(pos)
        else:
            p.remove()


def outline(ax, tri, **kw):
    """Draw the boundary loops of a triangulation."""
    from foamio import boundary_loops
    opts = dict(color=EDGE, lw=0.6, zorder=7)
    opts.update(kw)
    for loop in boundary_loops(tri):
        ax.plot(tri.x[loop], tri.y[loop], **opts)


def solid_patch(ax, tri, **kw):
    """Fill a (deformed) solid triangulation in the solid colour."""
    opts = dict(color=SOLID, zorder=5, rasterized=True)
    opts.update(kw)
    z = np.zeros(len(tri.x))
    ax.tripcolor(tri, z, cmap=matplotlib.colors.ListedColormap([opts.pop("color")]),
                 shading="flat", edgecolors="face", linewidth=0.1, **opts)
    outline(ax, tri, zorder=opts.get("zorder", 5) + 1, lw=0.5)


def preview(module, draw_args, out, width=2.6, height=1.6, letter="a"):
    """Write a one-panel preview (used when a panel script runs on its own)."""
    use_style()
    fig = plt.figure(figsize=(width, height))
    ax = fig.add_axes([0.02, 0.05, 0.80, 0.82])
    cax = fig.add_axes([0.86, 0.15, 0.03, 0.62])
    module.draw(fig, ax, cax, draw_args, letter)
    fig.savefig(out, dpi=400)
    print("wrote", out)


__all__ = ["use_style", "bare", "scalebar", "title", "colorbar", "sym_norm",
           "nice_ceil", "streamlines", "outline", "solid_patch", "preview",
           "Normalize", "SOLID", "UNDEFORMED"]
