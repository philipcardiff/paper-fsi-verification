"""Loading and drawing of a planar snapshot of a solids4foam FSI case: the
fluid on its (moving) mesh coloured by a field, with streamlines, and the
solid on its deformed configuration (undeformed outline dashed).

Used by the panel scripts of the collapsible channel, Turek-Hron FSI3, Mok
cavity, flexible-bottom cavity, Hessenthaler flap and blob in treacle."""
import numpy as np
from matplotlib.tri import Triangulation

import common
import foamio


def load(case, time=None, origin=None, normal=(0, 0, 1), axes=(0, 1),
         fluid_fields=("U", "p"), vorticity=False):
    """Fluid and solid of `case` on the cutting plane. The solid is displaced
    by its point displacement pointD (total Lagrangian), the fluid mesh has
    already moved."""
    _, tf, df, t = foamio.plane_cut(case, "fluid", time, fluid_fields,
                                    origin, normal, axes)
    _, ts, ds, _ = foamio.plane_cut(case, "solid", t, ("pointD",), origin,
                                    normal, axes)
    D = ds["pointD"]
    tsd = Triangulation(ts.x + D[:, axes[0]], ts.y + D[:, axes[1]],
                        ts.triangles)
    out = dict(t=t, fluid=tf, solid=tsd, solid0=ts, D=D, axes=axes, **df)
    if vorticity:
        out["vort"] = vorticity_z(tf, df["U"][:, axes[0]], df["U"][:, axes[1]])
    return out


def vorticity_z(tri, u, v):
    """Out-of-plane vorticity at the points: area-weighted average of the
    piecewise-linear (per-triangle) gradients."""
    x, y, T = tri.x, tri.y, tri.triangles
    x0, x1, x2 = x[T[:, 0]], x[T[:, 1]], x[T[:, 2]]
    y0, y1, y2 = y[T[:, 0]], y[T[:, 1]], y[T[:, 2]]
    det = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0)

    def grad(f):
        f0, f1, f2 = f[T[:, 0]], f[T[:, 1]], f[T[:, 2]]
        gx = ((f1 - f0) * (y2 - y0) - (f2 - f0) * (y1 - y0)) / det
        gy = ((f2 - f0) * (x1 - x0) - (f1 - f0) * (x2 - x0)) / det
        return gx, gy

    _, duy = grad(u)
    dvx, _ = grad(v)
    w = dvx - duy
    area = 0.5 * np.abs(det)
    num = np.zeros(len(x))
    den = np.zeros(len(x))
    for k in range(3):
        np.add.at(num, T[:, k], w * area)
        np.add.at(den, T[:, k], area)
    return num / np.maximum(den, 1e-300)


def speed(d):
    a = d["axes"]
    return np.hypot(d["U"][:, a[0]], d["U"][:, a[1]])


def draw(fig, ax, cax, d, field, norm, cmap, cb_label, xlim, ylim,
         stream=True, stream_density=0.9, stream_n=260, cb_ticks=None,
         solid_outline0=True, scale=1.0, extend="neither"):
    """Colour the fluid by `field` (point values), add streamlines and the
    deformed solid. `scale` multiplies all coordinates (e.g. 1000 for mm)."""
    tf = d["fluid"]
    if scale != 1.0:
        tf = Triangulation(tf.x * scale, tf.y * scale, tf.triangles)
    pc = ax.tripcolor(tf, field, cmap=cmap, norm=norm, shading="gouraud",
                      rasterized=True, zorder=1)
    if stream:
        a = d["axes"]
        vec = np.c_[d["U"][:, a[0]], d["U"][:, a[1]]]
        common.streamlines(ax, tf, vec, xlim, ylim, n=stream_n,
                           density=stream_density)
    ts = d["solid"]
    if scale != 1.0:
        ts = Triangulation(ts.x * scale, ts.y * scale, ts.triangles)
    common.solid_patch(ax, ts)
    if solid_outline0:
        s0 = d["solid0"]
        from foamio import boundary_loops
        for loop in boundary_loops(s0):
            ax.plot(s0.x[loop] * scale, s0.y[loop] * scale, **common.UNDEFORMED)
    common.outline(ax, tf, lw=0.6)
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    common.bare(ax)
    common.colorbar(fig, pc, cax, cb_label, ticks=cb_ticks, extend=extend)
    return pc
