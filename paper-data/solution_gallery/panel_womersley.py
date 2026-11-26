"""Gallery panel: pulsatile compliant (Womersley) tube, exact solution.

The exact linear travelling-wave solution of the thick elastic tube coupled
to the linearised Navier-Stokes equations is read from the tutorial's
constant/womersleyProperties (written by solids4foam tutorials/
fluidSolidInteraction/womersleyTube/verification/scripts/womersley_exact.py),
which gives omega, the complex wave number k, s^2 = k^2 + i omega/nu, the
(kinematic) pressure coefficient A and the velocity coefficient B:

    p   = A I0(k r),
    u_x = k A/omega I0(k r) + B I0(s r),
    u_r = i k A/omega I1(k r) + i k B I1(s r)/s,

all times exp(i (omega t - k x)). The wall's radial displacement at r = Ri is
u_r(Ri)/(i omega). The panel shows a meridional section (r mirrored) of the
whole tube, 0 <= x <= L, coloured by u_x, with u_x(r) profiles at five
stations, and the wall displaced radially with the displacement exaggerated.

Usage: panel_womersley.py WOMERSLEY_PROPERTIES OUT.pdf [--phase DEG]"""
import argparse
import re

import numpy as np

import common

TITLE = "Womersley tube"
LENGTH = 15.0
RHO = 1000.0
EXAGGERATION = 400.0
ASPECT = 2.8 / 15.2


def read_properties(path):
    txt = open(path).read()

    def scalar(key):
        return float(re.search(rf"^{key}\s+([^;]+);", txt, re.M).group(1))

    def cplx(key):
        m = re.search(rf"^{key}\s+\(\s*(\S+)\s+(\S+)\s*\);", txt, re.M)
        return complex(float(m.group(1)), float(m.group(2)))

    return dict(omega=scalar("omega"), k=cplx("k"), s=cplx("s"), A=cplx("A"),
                B=cplx("B"), ri=scalar("innerRadius"), ro=scalar("outerRadius"))


def bessel_i(n, z):
    """Modified Bessel I_n (n = 0, 1) of complex argument, power series."""
    z = np.asarray(z, dtype=complex)
    term = (z / 2) ** n / (1 if n == 0 else 1)
    total = term.copy()
    q = (z / 2) ** 2
    for m in range(1, 200):
        term = term * q / (m * (m + n))
        total += term
    return total


def load(props, phase_deg=140.0):
    c = read_properties(props)
    w, k, s, A, B, ri, ro = (c[x] for x in ("omega", "k", "s", "A", "B", "ri", "ro"))
    x = np.linspace(0.0, LENGTH, 301)
    r = np.linspace(0.0, ri, 61)
    X, R = np.meshgrid(x, r, indexing="ij")
    wave = np.exp(1j * (np.deg2rad(phase_deg) - k * X))
    ux = ((k * A / w * bessel_i(0, k * R) + B * bessel_i(0, s * R)) * wave).real
    p = (RHO * A * bessel_i(0, k * R) * wave).real
    urw = 1j * k * A / w * bessel_i(1, k * ri) + 1j * k * B * bessel_i(1, s * ri) / s
    xi = (urw / (1j * w) * np.exp(1j * (np.deg2rad(phase_deg) - k * x))).real
    # deformed meridional section: the fluid radius stretched with the wall
    Rd = R * (1 + EXAGGERATION * xi[:, None] / ri)
    return dict(x=x, X=X, R=R, Rd=Rd, ux=ux, p=p, xi=xi, ri=ri, ro=ro,
                omega=w, k=k, phase=phase_deg)


def draw(fig, ax, cax, d, letter="b"):
    vmax = common.nice_ceil(1e3 * np.abs(d["ux"]).max())
    norm = common.sym_norm(vmax)
    for sgn in (1, -1):
        pc = ax.pcolormesh(d["X"], sgn * d["Rd"], 1e3 * d["ux"], cmap="RdBu_r",
                           norm=norm, shading="gouraud", rasterized=True)
    ex = EXAGGERATION
    ri, ro, x, xi = d["ri"], d["ro"], d["x"], d["xi"]
    for sgn in (1, -1):
        lo, hi = sgn * (ri + ex * xi), sgn * (ro + ex * xi)
        ax.fill_between(x, lo, hi, color=common.SOLID, lw=0, zorder=5)
        ax.plot(x, lo, color="0.1", lw=0.5, zorder=6)
        ax.plot(x, hi, color="0.1", lw=0.5, zorder=6)
        ax.plot([0, x[-1]], [sgn * ri, sgn * ri], **common.UNDEFORMED)
    ax.plot([0, x[-1]], [0, 0], color="0.4", lw=0.4, ls=":", zorder=3)
    # u_x(r) profiles at five stations, drawn as arrows from the station
    umax = np.abs(d["ux"]).max()
    scale = 1.1 / umax
    for xs in (1.5, 4.5, 7.5, 10.5, 13.5):
        i = int(np.argmin(np.abs(x - xs)))
        rr = d["Rd"][i]
        u = d["ux"][i]
        rr_full = np.r_[-rr[::-1], rr[1:]]
        u_full = np.r_[u[::-1], u[1:]]
        ax.plot([x[i], x[i]], [-rr[-1], rr[-1]], color="0.15", lw=0.4, zorder=6)
        ax.plot(x[i] + scale * u_full, rr_full, color="black", lw=0.8, zorder=6)
        for kk in range(4, len(rr_full) - 1, 15):
            ax.annotate("", (x[i] + scale * u_full[kk], rr_full[kk]),
                        (x[i], rr_full[kk]),
                        arrowprops=dict(arrowstyle="-|>,head_length=0.25,head_width=0.12",
                                        lw=0.4, color="black",
                                        shrinkA=0, shrinkB=0), zorder=6)
    ax.set_xlim(-0.1, LENGTH + 0.1)
    ax.set_ylim(-1.5, 1.3)
    common.bare(ax)
    common.scalebar(ax, 0.0, -1.38, 2.0, "2 m", side="right")
    common.title(ax, letter, rf"{TITLE}, $\omega t = {d['phase']:g}^\circ$, "
                 rf"wall displacement $\times{EXAGGERATION:g}$")
    common.colorbar(fig, pc, cax, r"$u_x$ (mm/s)")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("props")
    ap.add_argument("out")
    ap.add_argument("--phase", type=float, default=140.0,
                    help="omega t in degrees (default 140)")
    a = ap.parse_args()
    d = load(a.props, a.phase)
    print(f"max |u_x| = {np.abs(d['ux']).max():.3e} m/s, "
          f"max |xi| = {np.abs(d['xi']).max():.3e} m, wave speed "
          f"{d['omega'] / d['k'].real:.6f} m/s")
    import sys
    common.preview(sys.modules[__name__], d, a.out, width=5.0, height=1.1)


if __name__ == "__main__":
    main()
