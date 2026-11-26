"""Gallery panel: ring added-mass problem, exact n = 2 mode (strong level).

Everything is evaluated from the exact solution, so no simulation data are
used:
  - the solid mode shape (U(r) cos 2 theta, V(r) sin 2 theta) is the
    eigenvector of the plane-strain continuum ring loaded by the potential-flow
    added-mass traction, by Chebyshev collocation (the method and matrices of
    solids4foam tutorials/fluidSolidInteraction/ringAddedMass/verification/
    scripts/ring_exact_frequencies.py, repeated here so that the panel needs
    only numpy);
  - the fluid is the potential flow in the annulus Ro < r < b,
    phi = (A r^2 + B r^-2) cos 2 theta, with the wall's radial velocity at
    r = Ro and none at r = b; the pressure is p = -rho_F d(phi)/dt and the
    streamlines are contours of psi = (A r^2 - B r^-2) sin 2 theta, whose
    shape does not change during the oscillation.
The snapshot is the instant of maximum displacement (velocity zero, pressure
extreme). Displacements are exaggerated by EXAGGERATION; the fluid annulus is
drawn on the correspondingly deformed region (radial blend).

Usage: panel_ring.py OUT.pdf [--level strong]"""
import argparse

import numpy as np

import common

TITLE = "Ring added mass"
GEOM = dict(ri=0.95, ro=1.05, b=1.5, E=1.0e6, nu=0.3, rho_s=1000.0)
RHO_F = {"weak": 14.0, "moderate": 140.0, "strong": 1400.0}
N = 2
EXAGGERATION = 0.10      # maximum displacement drawn as 0.10 m (10% of a)
ASPECT = 1.0             # height / width of the drawn window


def chebyshev(n, r_min, r_max):
    x = np.cos(np.pi * np.arange(n + 1) / n)
    c = np.hstack([2.0, np.ones(n - 1), 2.0]) * (-1.0) ** np.arange(n + 1)
    dx = np.subtract.outer(x, x)
    d = np.outer(c, 1.0 / c) / (dx + np.eye(n + 1))
    d -= np.diag(d.sum(axis=1))
    r = r_min + 0.5 * (x + 1.0) * (r_max - r_min)
    return r, d * 2.0 / (r_max - r_min)


def added_mass(radius, outer, rho_f):
    return rho_f * radius / N * (outer ** (2 * N) + radius ** (2 * N)) \
        / (outer ** (2 * N) - radius ** (2 * N))


def mode(rho_f, n_solid=24):
    """Wet n = 2 eigenvalue and the mode (r, U, V), U(Ro) = 1."""
    g = GEOM
    lam = g["E"] * g["nu"] / ((1 + g["nu"]) * (1 - 2 * g["nu"]))
    mu = g["E"] / (2 * (1 + g["nu"]))
    r, d = chebyshev(n_solid, g["ri"], g["ro"])
    ns = n_solid + 1
    eye, rinv = np.eye(ns), np.diag(1.0 / r)
    del_u, del_v = d + rinv, N * rinv
    rot_u, rot_v = N * rinv, d + rinv
    lru = (lam + 2 * mu) * d @ del_u - mu * N * rinv @ rot_u
    lrv = (lam + 2 * mu) * d @ del_v - mu * N * rinv @ rot_v
    ltu = -(lam + 2 * mu) * N * rinv @ del_u + mu * d @ rot_u
    ltv = -(lam + 2 * mu) * N * rinv @ del_v + mu * d @ rot_v
    srr_u, srr_v = lam * del_u + 2 * mu * d, lam * del_v
    srt_u, srt_v = -mu * N * rinv, mu * (d - rinv)
    a = np.zeros((4 * ns, 4 * ns))
    bm = np.zeros((4 * ns, 4 * ns))
    iu, iv, iwu, iwv = 0, ns, 2 * ns, 3 * ns
    a[iu:iu + ns, iwu:iwu + ns] = eye
    bm[iu:iu + ns, iu:iu + ns] = eye
    a[iv:iv + ns, iwv:iwv + ns] = eye
    bm[iv:iv + ns, iv:iv + ns] = eye
    a[iwu:iwu + ns, iu:iu + ns] = lru
    a[iwu:iwu + ns, iv:iv + ns] = lrv
    bm[iwu:iwu + ns, iwu:iwu + ns] = g["rho_s"] * eye
    a[iwv:iwv + ns, iu:iu + ns] = ltu
    a[iwv:iwv + ns, iv:iv + ns] = ltv
    bm[iwv:iwv + ns, iwv:iwv + ns] = g["rho_s"] * eye
    for row, su, sv in ((iwu, srr_u, srr_v), (iwv, srt_u, srt_v)):
        for k in (n_solid, 0):
            a[row + k, :] = 0
            bm[row + k, :] = 0
            a[row + k, iu:iu + ns] = su[k]
            a[row + k, iv:iv + ns] = sv[k]
    # interface: sigma_rr(Ro) = -m_a s w_r(Ro)
    bm[iwu + 0, iwu + 0] = -added_mass(g["ro"], g["b"], rho_f)
    # thin-ring estimate as the shift
    h, am = g["ro"] - g["ri"], 0.5 * (g["ro"] + g["ri"])
    D = g["E"] * h ** 3 / (12 * (1 - g["nu"] ** 2))
    wdry = np.sqrt(D * N ** 2 * (N ** 2 - 1) ** 2
                   / (g["rho_s"] * h * am ** 4 * (N ** 2 + 1)))
    mstar = added_mass(g["ro"], g["b"], rho_f) * g["ro"] * N ** 2 \
        / (g["rho_s"] * h * am * (N ** 2 + 1))
    sigma = 1j * wdry / np.sqrt(1 + mstar)
    ev, vec = np.linalg.eig(np.linalg.solve(a - sigma * bm, bm))
    k = np.argmax(np.abs(ev))
    s = sigma + 1.0 / ev[k]
    x = vec[:, k]
    U, V = x[iu:iu + ns], x[iv:iv + ns]
    scale = U[0]
    U, V = (U / scale).real, (V / scale).real
    order = np.argsort(r)
    return s, r[order], U[order], V[order], mstar


def load(level="strong"):
    rho_f = RHO_F[level]
    s, r, U, V, mstar = mode(rho_f)
    g = GEOM
    th = np.linspace(0, 2 * np.pi, 361)
    amp = EXAGGERATION / max(np.abs(U).max(), np.abs(V).max())
    # solid, deformed
    rs = np.linspace(g["ri"], g["ro"], 9)
    Ur, Vr = np.interp(rs, r, U), np.interp(rs, r, V)
    RS, TH = np.meshgrid(rs, th, indexing="ij")
    ur = amp * Ur[:, None] * np.cos(N * TH)
    ut = amp * Vr[:, None] * np.sin(N * TH)
    xs = (RS + ur) * np.cos(TH) - ut * np.sin(TH)
    ys = (RS + ur) * np.sin(TH) + ut * np.cos(TH)
    # fluid: unit radial wall velocity W = 1 at Ro (shape only)
    ro, b = g["ro"], g["b"]
    A = 1.0 / (2 * (ro - b ** 4 / ro ** 3))
    B = A * b ** 4
    rf = np.linspace(ro, b, 60)
    RF, THF = np.meshgrid(rf, th, indexing="ij")
    phi = (A * RF ** 2 + B / RF ** 2) * np.cos(N * THF)
    psi = (A * RF ** 2 - B / RF ** 2) * np.sin(N * THF)
    # at maximum outward displacement of theta = 0 the acceleration there is
    # -omega^2 times the displacement: p = -rho dphi/dt  ~  +phi (normalised)
    p = phi / np.abs(phi).max()
    blend = (b - RF) / (b - ro)
    urf = amp * U[-1] * np.cos(N * THF) * blend
    xf = (RF + urf) * np.cos(THF)
    yf = (RF + urf) * np.sin(THF)
    return dict(xs=xs, ys=ys, xf=xf, yf=yf, p=p, psi=psi, omega=s.imag,
                mstar=mstar, level=level, amp=amp)


def draw(fig, ax, cax, d, letter="a"):
    g = GEOM
    norm = common.sym_norm(1.0)
    pc = ax.pcolormesh(d["xf"], d["yf"], d["p"], cmap="RdBu_r", norm=norm,
                       shading="gouraud", rasterized=True)
    lv = np.linspace(-0.9, 0.9, 10) * np.abs(d["psi"]).max()
    ax.contour(d["xf"], d["yf"], d["psi"], levels=lv, colors="0.15",
               linewidths=0.4, zorder=4)
    # solid ring (deformed) and the undeformed ring surfaces
    ax.fill(np.r_[d["xs"][-1], d["xs"][0][::-1]], np.r_[d["ys"][-1], d["ys"][0][::-1]],
            color=common.SOLID, zorder=5, lw=0)
    ax.plot(d["xs"][0], d["ys"][0], color="0.1", lw=0.5, zorder=6)
    ax.plot(d["xs"][-1], d["ys"][-1], color="0.1", lw=0.5, zorder=6)
    th = np.linspace(0, 2 * np.pi, 361)
    for r in (g["ri"], g["ro"]):
        ax.plot(r * np.cos(th), r * np.sin(th), **common.UNDEFORMED)
    ax.plot(g["b"] * np.cos(th), g["b"] * np.sin(th), color="0.1", lw=0.9,
            zorder=7)
    # modelled quadrant
    ax.plot([0, 0], [0, g["b"]], color="0.4", lw=0.4, ls=":", zorder=3)
    ax.plot([0, g["b"]], [0, 0], color="0.4", lw=0.4, ls=":", zorder=3)
    lim = g["b"] + 0.04
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    common.bare(ax)
    common.scalebar(ax, -0.5, -0.25, 1.0, "1 m")
    common.title(ax, letter, f"{TITLE}, $m^* = {d['mstar']:.2f}$")
    common.colorbar(fig, pc, cax, r"$p/|p|_{\max}$", ticks=[-1, 0, 1])


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out")
    ap.add_argument("--level", default="strong", choices=sorted(RHO_F))
    a = ap.parse_args()
    d = load(a.level)
    print(f"omega = {d['omega']:.6f} rad/s, m* = {d['mstar']:.4g}")
    import sys
    common.preview(sys.modules[__name__], d, a.out, width=2.0, height=1.8)


if __name__ == "__main__":
    main()
