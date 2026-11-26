"""The manufactured FSI case (fig:mms_fsi_case of the verification paper):
snapshots of the exact solution of the full variant, and the space-time
evolution of the interface height.

Each snapshot shows, in the current configuration,
  - the deformed fluid channel coloured by the speed |v_F|, with the
    instantaneous streamlines (contours of the stream function psi, obtained
    by integrating v_F,x across the channel from psi = 0 on the wall) and
    arrows of the velocity direction;
  - the deformed solid wall coloured by the Frobenius norm of the
    Green-Lagrange strain E = (F^T F - I)/2 (plane strain, F_zz = 1), with
    F = I + (grad_0 D_S)^T from central differences of the exact displacement;
  - deformed material lines X = const and Y = const of the solid;
  - the undeformed interface Y = H as a dashed line.
Bottom left: the interface displacement eta = h(x, t) - H of the full
variant over [0, tmax], with the end of the runs and the snapshot times
marked. Bottom right: eta(x) of the full and flat variants at the --compare
times.

All fields are evaluated from the closed-form solution with the case's exact
evaluator (case/src/evalExact.C compiled against derive/generated or
derive/generated_flat), which reads lines
  fluid x y t      -> Ux Uy p fx fy dpdx dpdy gradU(4)
  solid X Y t      -> Dx Dy Ux Uy bx by
  interface X t    -> x y tx ty nx ny
from standard input.

Usage:
  mms_fsi_case_figure.py EVAL_FULL EVAL_FLAT OUT.pdf [--times T1 T2 T3 T4]
      [--compare T ...] [--tmax 1] [--trun 0.4] [--dpsi 0.06] [--width W]
      [--png OUT.png]
  EVAL_FULL, EVAL_FLAT: the exact evaluators of the full and flat variants.
Prints the peak speed, the strain, J and the interface height of each
snapshot, as a check against the values quoted in the paper."""
import argparse
import subprocess

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.ticker import FormatStrFormatter

plt.rcParams.update({"font.size": 8, "mathtext.fontset": "cm",
                     "font.family": "serif", "pdf.fonttype": 42,
                     "axes.linewidth": 0.6, "xtick.major.width": 0.6,
                     "ytick.major.width": 0.6, "xtick.major.size": 2.5,
                     "ytick.major.size": 2.5})

L, H, T = 2.0, 1.0, 0.5
FD_STEP = 1e-6          # central-difference step for F (relative to L = 2)


# ---- exact evaluator ---------------------------------------------------------
def evaluate(exe, kind, cols):
    """Run the evaluator on the points cols = (a, b[, t]) (equal-shape arrays);
    returns an array of shape (npts, nout)."""
    flat = [np.asarray(c, dtype=float).ravel() for c in cols]
    lines = "\n".join(kind + " " + " ".join(f"{v:.17g}" for v in row)
                      for row in zip(*flat)) + "\n"
    out = subprocess.run([exe], input=lines, capture_output=True, text=True,
                         check=True).stdout
    return np.array(out.split(), dtype=float).reshape(len(flat[0]), -1)


def interface(exe, t, nX=401):
    X = np.linspace(0.0, L, nX)
    v = evaluate(exe, "interface", (X, np.full_like(X, t)))
    return v[:, 0], v[:, 1]


def channel_height(exe, x, t):
    """h(x, t): the interface is the graph y = h(x); interpolate the exact
    interface curve (fine parametrisation) at the abscissae x."""
    xi, yi = interface(exe, t, nX=2001)
    assert np.all(np.diff(xi) > 0), "interface abscissa not monotonic"
    return np.interp(x, xi, yi)


def fluid_fields(exe, t, nx=241, ns=121):
    """Fluid on the deformed channel: a grid (x, s) with y = s h(x, t).
    Returns x, y, |v|, vx, vy, psi on the (ns, nx) grid."""
    x1 = np.linspace(0.0, L, nx)
    s1 = np.linspace(0.0, 1.0, ns)
    h = channel_height(exe, x1, t)
    assert np.all(h > 0)
    xx, ss = np.meshgrid(x1, s1)
    yy = ss * h[None, :]
    v = evaluate(exe, "fluid", (xx, yy, np.full_like(xx, t)))
    vx = v[:, 0].reshape(xx.shape)
    vy = v[:, 1].reshape(xx.shape)
    # psi(x, y) = int_0^y vx dy' (psi = 0 on the no-slip wall); dy = h ds
    dy = np.diff(yy, axis=0)
    psi = np.vstack([np.zeros((1, nx)),
                     np.cumsum(0.5 * (vx[1:] + vx[:-1]) * dy, axis=0)])
    return xx, yy, np.hypot(vx, vy), vx, vy, psi


def solid_fields(exe, t, nX=161, nY=41):
    """Solid on the reference grid (X, Y): current positions and strain
    measures from the exact displacement. Returns X, Y, x, y, |E|_F, J."""
    X1 = np.linspace(0.0, L, nX)
    Y1 = np.linspace(H, H + T, nY)
    XX, YY = np.meshgrid(X1, Y1)
    d = FD_STEP
    pts = [(XX, YY), (XX + d, YY), (XX - d, YY), (XX, YY + d), (XX, YY - d)]
    D = [evaluate(exe, "solid", (a, b, np.full_like(a, t)))[:, :2]
         .reshape(XX.shape + (2,)) for a, b in pts]
    dDdX = (D[1] - D[2]) / (2 * d)       # (.., 2): dDx/dX, dDy/dX
    dDdY = (D[3] - D[4]) / (2 * d)
    F11, F21 = 1 + dDdX[..., 0], dDdX[..., 1]
    F12, F22 = dDdY[..., 0], 1 + dDdY[..., 1]
    C11 = F11**2 + F21**2
    C22 = F12**2 + F22**2
    C12 = F11 * F12 + F21 * F22
    E11, E22, E12 = 0.5 * (C11 - 1), 0.5 * (C22 - 1), 0.5 * C12
    Enorm = np.sqrt(E11**2 + E22**2 + 2 * E12**2)
    J = F11 * F22 - F12 * F21
    return XX, YY, XX + D[0][..., 0], YY + D[0][..., 1], Enorm, J


def material_lines(exe, t, nX_lines=9, nY_lines=5, n=161):
    """Deformed images of the reference lines X = const and Y = const."""
    lines = []
    s = np.linspace(0, 1, n)
    for X0 in np.linspace(0, L, nX_lines):
        X, Y = np.full_like(s, X0), H + T * s
        D = evaluate(exe, "solid", (X, Y, np.full_like(s, t)))
        lines.append((X + D[:, 0], Y + D[:, 1]))
    for Y0 in np.linspace(H, H + T, nY_lines):
        X, Y = L * s, np.full_like(s, Y0)
        D = evaluate(exe, "solid", (X, Y, np.full_like(s, t)))
        lines.append((X + D[:, 0], Y + D[:, 1]))
    return lines


# ---- drawing -----------------------------------------------------------------
def streamline_arrows(ax, cs, xx, yy, vx, vy, min_len=0.15):
    """One arrowhead on each streamline (contour path of psi) longer than
    min_len, oriented by the local velocity; the position along the path is
    staggered from line to line."""
    pts = np.column_stack([xx.ravel(), yy.ravel()])
    vel = np.column_stack([vx.ravel(), vy.ravel()])
    k = 0
    for segs in cs.allsegs:
        for seg in segs:
            if len(seg) < 3:
                continue
            ds = np.hypot(*np.diff(seg, axis=0).T)
            arc = np.concatenate([[0.0], np.cumsum(ds)])
            if arc[-1] < min_len:
                continue
            frac = 0.3 + 0.4 * ((k * 0.618034) % 1.0)
            k += 1
            i = min(np.searchsorted(arc, frac * arc[-1]), len(seg) - 2)
            p0 = seg[i]
            tang = seg[i + 1] - seg[i]
            u = vel[np.argmin(np.sum((pts - p0) ** 2, axis=1))]
            if np.dot(tang, u) < 0:
                tang = -tang
            tang = tang / np.hypot(*tang)
            ln = min(0.07, 0.2 * arc[-1])
            ms = 10.0 * min(1.0, max(0.45, arc[-1] / 0.6))   # head size, small loops
            ax.annotate("", xy=(p0[0] + 0.5 * ln * tang[0], p0[1] + 0.5 * ln * tang[1]),
                        xytext=(p0[0] - 0.5 * ln * tang[0], p0[1] - 0.5 * ln * tang[1]),
                        arrowprops=dict(arrowstyle="-|>,head_length=0.35,head_width=0.18",
                                        color="white", lw=0.6, shrinkA=0, shrinkB=0,
                                        mutation_scale=ms),
                        zorder=6)


def vortex_level(f, ring=8, frac=0.7):
    """A stream-function level that draws a closed streamline around the
    recirculation centre: the interior point of least speed, away from the
    walls, if psi there is an extremum relative to a square ring of grid
    points around it. The level lies between the centre value and the
    ring's closest value, so the contour closes inside the ring."""
    xx, yy, speed, vx, vy, psi = f
    ns, nx = psi.shape
    j0, j1, i0, i1 = int(0.15 * ns), int(0.9 * ns), int(0.08 * nx), int(0.92 * nx)
    sub = speed[j0:j1, i0:i1]
    j, i = np.unravel_index(np.argmin(sub), sub.shape)
    j, i = j + j0, i + i0
    if speed[j, i] > 0.05 * speed.max():
        return []
    c, k = psi[j, i], ring
    r = np.concatenate([psi[j - k, i - k:i + k + 1], psi[j + k, i - k:i + k + 1],
                        psi[j - k:j + k + 1, i - k], psi[j - k:j + k + 1, i + k]])
    if np.all(r > c):
        return [c + frac * (r.min() - c)]
    if np.all(r < c):
        return [c + frac * (r.max() - c)]
    return []


def snapshot(ax, exe, f, s, mlines, norm_v, norm_e, psi_levels):
    xx, yy, speed, vx, vy, psi = f
    XX, YY, xs, ys, Enorm, J = s

    pv = ax.pcolormesh(xx, yy, speed, cmap="viridis", norm=norm_v,
                       shading="gouraud", rasterized=True)
    cs = ax.contour(xx, yy, psi, levels=np.union1d(psi_levels, vortex_level(f)),
                    colors="white", linewidths=0.5, alpha=0.9)
    streamline_arrows(ax, cs, xx, yy, vx, vy)

    pe = ax.pcolormesh(xs, ys, Enorm, cmap="YlOrRd", norm=norm_e,
                       shading="gouraud", rasterized=True)
    for lx, ly in mlines:
        ax.plot(lx, ly, color="0.15", lw=0.45)

    # undeformed interface, boundaries of the current domains
    ax.plot([0, L], [H, H], color="black", lw=0.6, ls=(0, (3, 2)), zorder=7)
    ax.plot(xx[-1], yy[-1], color="black", lw=0.9, zorder=7)
    ax.plot([0, L, L, 0, 0], [H + T, H + T, 0, 0, H + T], color="black",
            lw=0.6, zorder=7)
    ax.set_xlim(-0.03, L + 0.03)
    ax.set_ylim(-0.03, H + T + 0.03)
    ax.set_aspect("equal")
    ax.set_xticks([0, 0.5, 1, 1.5, 2])
    ax.set_yticks([0, 0.5, 1, 1.5])
    ax.xaxis.set_major_formatter(FormatStrFormatter("%g"))
    ax.yaxis.set_major_formatter(FormatStrFormatter("%g"))
    return pv, pe


def spacetime(exe, tmax, nt=201, nx=201):
    ts = np.linspace(0.0, tmax, nt)
    x = np.linspace(0.0, L, nx)
    eta = np.array([channel_height(exe, x, t) - H for t in ts])
    return x, ts, eta


def label_boundaries(ax):
    """Boundary names, drawn inside the first snapshot."""
    kw = dict(fontsize=7, ha="center", va="center", zorder=11,
              bbox=dict(boxstyle="round,pad=0.12", fc="#440154", ec="none",
                        alpha=0.75))
    ax.text(0.05, 0.45, "inlet", rotation=90, color="white", **kw)
    ax.text(L - 0.05, 0.45, "outlet", rotation=90, color="white", **kw)
    ax.text(1.0, 0.05, "no-slip wall", color="white", **kw)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("eval_full", help="exact evaluator, full variant")
    ap.add_argument("eval_flat", help="exact evaluator, flat variant")
    ap.add_argument("out", help="output PDF")
    ap.add_argument("--times", type=float, nargs=4, default=[0.1, 0.2, 0.3, 0.4],
                    help="the four snapshot times (default 0.1 0.2 0.3 0.4)")
    ap.add_argument("--compare", type=float, nargs="+", default=[0.2, 0.4],
                    help="times of the full/flat interface comparison (default 0.2 0.4)")
    ap.add_argument("--tmax", type=float, default=1.0,
                    help="end time of the space-time panel (default 1, one period)")
    ap.add_argument("--trun", type=float, default=0.4,
                    help="end time of the runs, marked in the space-time panel")
    ap.add_argument("--dpsi", type=float, default=0.06,
                    help="stream-function spacing of the streamlines (default 0.06)")
    ap.add_argument("--png", help="also write a PNG preview (300 dpi)")
    ap.add_argument("--width", type=float, default=372 / 72,
                    help="figure width in inches (default: the paper's text width, 372 pt)")
    a = ap.parse_args()
    cols = ["#0072B2", "#D55E00", "#009E73", "#CC79A7"]   # Okabe-Ito
    if len(a.compare) > len(cols):
        ap.error(f"--compare takes at most {len(cols)} times")

    # fields of the snapshots; common colour scales from the plotted grids
    flu = [fluid_fields(a.eval_full, t) for t in a.times]
    sol = [solid_fields(a.eval_full, t) for t in a.times]
    mls = [material_lines(a.eval_full, t) for t in a.times]
    vmax = max(f[2].max() for f in flu)
    psimin = min(f[5].min() for f in flu)
    psimax = max(f[5].max() for f in flu)
    emax = max(s_[4].max() for s_ in sol)
    norm_v = Normalize(0.0, np.ceil(vmax * 10) / 10)
    norm_e = Normalize(0.0, np.ceil(emax * 50) / 50)
    dpsi = a.dpsi
    psi_levels = np.arange(np.floor(psimin / dpsi) * dpsi + dpsi / 2,
                           psimax, dpsi)

    W = a.width
    fig = plt.figure(figsize=(W, 1.0 * W))
    left, right, cbx = 0.085, 0.855, 0.872
    gs = fig.add_gridspec(2, 2, left=left, right=right, top=0.99, bottom=0.40,
                          wspace=0.06, hspace=0.10)
    axs = [fig.add_subplot(gs[i // 2, i % 2]) for i in range(4)]
    letters = "abcd"
    for k, (ax, t) in enumerate(zip(axs, a.times)):
        pv, pe = snapshot(ax, a.eval_full, flu[k], sol[k], mls[k], norm_v, norm_e,
                          psi_levels)
        ax.text(0.02, 0.975, f"({letters[k]}) $t = {t:g}$", transform=ax.transAxes,
                ha="left", va="top", fontsize=8,
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85),
                zorder=12)
        if k % 2:
            plt.setp(ax.get_yticklabels(), visible=False)
        else:
            ax.set_ylabel("$y$", labelpad=1)
        if k < 2:
            plt.setp(ax.get_xticklabels(), visible=False)
        else:
            ax.set_xlabel("$x$", labelpad=1)
    label_boundaries(axs[0])

    # colour bars to the right of the snapshot grid
    b0, b1 = axs[3].get_position(), axs[1].get_position()
    cax_e = fig.add_axes([cbx, b1.y0 + 0.01, 0.017, b1.height - 0.02])
    cax_v = fig.add_axes([cbx, b0.y0 + 0.01, 0.017, b0.height - 0.02])
    cb = fig.colorbar(pe, cax=cax_e)
    cb.set_label(r"solid: $\|\mathbf{E}\|_F$", fontsize=8)
    cb = fig.colorbar(pv, cax=cax_v)
    cb.set_label(r"fluid: $|\mathbf{v}_F|$", fontsize=8)
    for c in (cax_e, cax_v):
        c.tick_params(labelsize=7)

    # (e) space-time interface displacement, full variant
    bot, top = 0.075, 0.305
    axe = fig.add_axes([left, bot, 0.40, top - bot])
    x, ts, eta = spacetime(a.eval_full, a.tmax)
    lim = np.ceil(eta.max() * 100) / 100
    lv = np.linspace(0, lim, 14)
    lv[0] = -1e-9
    cf = axe.contourf(x, ts, eta, levels=lv, cmap="cividis")
    cf.set_rasterized(True)
    axe.contour(x, ts, eta, levels=lv[1:], colors="white", linewidths=0.3,
                alpha=0.5)
    axe.plot([0, L - 0.32], [a.trun, a.trun], color="white", lw=0.8)
    for k, t in enumerate(a.times):
        axe.plot([L - 0.1, L], [t, t], color="white", lw=0.8)
        axe.text(L - 0.13, t, f"({letters[k]})", color="white", fontsize=7,
                 ha="right", va="center")
    axe.text(0.05, a.trun + 0.015, "end of runs", color="white", fontsize=7,
             ha="left", va="bottom")
    axe.set_xlabel("$x$", labelpad=1)
    axe.set_ylabel("$t$", labelpad=1)
    axe.set_xticks([0, 0.5, 1, 1.5, 2])
    axe.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    axe.xaxis.set_major_formatter(FormatStrFormatter("%g"))
    axe.yaxis.set_major_formatter(FormatStrFormatter("%g"))
    axe.set_title(r"(e) $\eta = h(x, t) - H$, full variant", fontsize=8, loc="left",
                  pad=3)
    pa = axe.get_position()
    cax_h = fig.add_axes([pa.x1 + 0.012, pa.y0, 0.017, pa.height])
    cb = fig.colorbar(cf, cax=cax_h, ticks=np.arange(0, lim + 1e-9, 0.04))
    cb.ax.yaxis.set_major_formatter(FormatStrFormatter("%g"))
    cb.ax.set_title(r"$\eta$", fontsize=8, pad=3)
    cax_h.tick_params(labelsize=7)

    # (f) full and flat interfaces, vertical scale enlarged
    axf = fig.add_axes([0.665, bot, 0.98 - 0.665, top - bot])
    for c, t in zip(cols, a.compare):
        for exe, ls, nm in ((a.eval_full, "-", "full"), (a.eval_flat, "--", "flat")):
            xi, yi = interface(exe, t)
            axf.plot(xi, yi - H, color=c, ls=ls, lw=1.0)
            if nm == "full":
                k = np.argmax(yi)
                axf.text(xi[k], yi[k] - H + 0.006, f"$t = {t:g}$", color=c,
                         fontsize=7, ha="center", va="bottom")
    axf.plot([], [], color="0.3", ls="-", lw=1.0, label="full")
    axf.plot([], [], color="0.3", ls="--", lw=1.0, label="flat")
    axf.legend(loc="upper left", fontsize=7, frameon=False, handlelength=2.0,
               borderaxespad=0.2, labelspacing=0.2)
    axf.set_xlim(0, L)
    axf.set_ylim(-0.005, 0.16)
    axf.set_xticks([0, 0.5, 1, 1.5, 2])
    axf.set_yticks([0, 0.05, 0.1, 0.15])
    axf.xaxis.set_major_formatter(FormatStrFormatter("%g"))
    axf.yaxis.set_major_formatter(FormatStrFormatter("%g"))
    axf.set_xlabel("$x$", labelpad=1)
    axf.set_ylabel(r"$\eta$", labelpad=1)
    axf.set_title("(f) full and flat variants", fontsize=8, loc="left",
                  pad=3)

    print("  t   max|v|  max|E|  minJ   maxJ   max eta  min eta")
    for t, f, s_ in zip(a.times, flu, sol):
        print("  %.2f %6.3f  %6.3f  %5.3f  %5.3f  %6.4f  %7.4f"
              % (t, f[2].max(), s_[4].max(), s_[5].min(), s_[5].max(),
                 f[1][-1].max() - H, f[1][-1].min() - H))
    print(f"  space-time (0 <= t <= {a.tmax:g}): max eta {eta.max():.4f}, "
          f"min eta {eta.min():.4f}")
    fig.savefig(a.out)
    if a.png:
        fig.savefig(a.png, dpi=300)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
