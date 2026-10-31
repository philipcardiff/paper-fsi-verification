#!/usr/bin/env python3
"""radialTorsionalRing (torsion member) references for the FEBio pilot (RTR agent).  Light numpy/scipy only (login-node safe).

PROVENANCE.  The solid law and the coupled semi-analytical reference are the user's own draft (solids4foam PR #637, head 30ae00eb, tutorials/.../verification/scripts/rt_reference.py,
sha256 f464e35e...), vendored UNCHANGED in reference_data/rtr/rt_reference_s4f637_30ae00eb.py (do not edit; hash is checked below).  This file adds
  (1) coupled_reference(): the proposal's own torsion-member reference on an arbitrary time grid (same call as the proposal's make_reference.py: Params(theta0=0.1, Omega=3,
      P0=0, T_ramp=1, ramp="Cinf"), 20 GLL solid nodes, 32 Chebyshev fluid nodes, Radau rtol 1e-11);
  (2) static_twist(): static finite-strain torsion of the annulus A<=R<=B (inner rotated by alpha_A, outer rotated by alpha_B, both radially clamped), by Newton on the
      vendored discrete energy (solid_forces) = the structural reference for stage RT-A (same law and kinematics as the proposal);
  (3) fluid_exact(): INDEPENDENT fixed-geometry fluid reference for RT-B, written here (not taken from the proposal): L = r v_theta obeys L_t = nu (L_rr - L_r/r) on 1<=r<=2,
      L(1,t) = b^2 omega_w(t) (no slip), L_r(2) = 2 L(2)/2 (zero shear stress, the proposal's open condition), solved by Chebyshev collocation + scipy Radau (rtol 1e-12) and
      cross-checked against a second discretisation (N and tolerances); torque T_B = b^2 sigma_rt = mu (L_r - 2 L) at r = b (per radian, per metre depth), H_f = rho int L r dr,
      centripetal pressure p(r) = -rho int_r^c v^2/r' dr' (p(c)=0, Q = 0).
Conventions follow the proposal: torques and angular momenta per radian and per metre of depth (extensive FEBio value / (2 pi dz)); T_B = torque of the fluid ON the solid.
"""
import sys, os, json, math, hashlib, importlib.util
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
VEND = os.path.join(ROOT, "reference_data", "rtr", "rt_reference_s4f637_30ae00eb.py")
VEND_SHA = "f464e35ebd5a1da87fcb52966cfe7868811e7521d2a615a4b16fa7e9e3193ba3"


def load_vendored():
    h = hashlib.sha256(open(VEND, "rb").read()).hexdigest()
    assert h == VEND_SHA, "vendored reference modified: " + h
    sys.dont_write_bytecode = True
    spec = importlib.util.spec_from_file_location("rt_reference_vendored", VEND)
    m = importlib.util.module_from_spec(spec)
    sys.modules["rt_reference_vendored"] = m
    spec.loader.exec_module(m)
    return m


rt = load_vendored()
A, B, C = 0.5, 1.0, 2.0
RHO, MU = 1000.0, 300.0
NU = MU / RHO
THETA0, OMEGA, TRAMP = 0.1, 3.0, 1.0


def ramp_S(t):
    """C-infinity ramp S(t/T_ramp) and derivative (per second), vendored definition."""
    S, S1, S2 = rt.ramp_cinf(t / TRAMP)
    return S, S1 / TRAMP, S2 / TRAMP**2


def alpha_drive(t, theta0=THETA0):
    """alpha(t) = theta0 S(t) sin(3 t), omega = d alpha/dt (same function the proposal applies to the inner solid surface)."""
    S, S1, _ = ramp_S(t)
    s, c = math.sin(OMEGA * t), math.cos(OMEGA * t)
    return theta0 * S * s, theta0 * (S1 * s + S * OMEGA * c)


def coupled_reference(tt, Ns=20, Nf=32):
    p = rt.Params(theta0=THETA0, Omega=OMEGA, P0=0.0, T_ramp=TRAMP, ramp="Cinf")
    m, sol = rt.solve(p, tt[-1], Ns, Nf, rtol=1e-11, t_eval=tt)
    H = rt.history(m, sol)
    return {k: np.asarray(v) for k, v in H.items()}


# ------------------------------------------------------------------ RT-A: static finite torsion of the solid
def static_twist(alA, alB, Ns=24, guess=None, law="neoHookean"):
    """Static equilibrium (dV/dphi_i = 0, dV/dalpha_i = 0, interior nodes) with u_0 = u_N = 0, alpha_0 = alA, alpha_N = alB.
    Returns dict with nodal R, u, alpha, torque T_A = dV/dalpha_0 (per radian, per metre), sigma fields; homotopy in the load handles large rotations."""
    p = rt.Params(theta0=0.0, law=law)
    m = rt.Model(p, Ns, 8)
    N = Ns

    def unpack(z):
        u = np.zeros(N + 1); al = np.zeros(N + 1)
        u[1:N] = z[:N - 1]; al[1:N] = z[N - 1:]
        al[0], al[N] = alA_cur, alB_cur
        return u, al

    def F(z):
        u, al = unpack(z)
        phi, Gphi, Gal, V, PrR, PtR = m.solid_forces(u, al)
        return np.concatenate((Gphi[1:N], Gal[1:N]))

    xr = (m.R - A) / (B - A)
    z = np.zeros(2 * (N - 1)) if guess is None else guess
    nstep = 40
    for s in np.linspace(0.0, 1.0, nstep + 1)[1:]:
        alA_cur, alB_cur = s * alA, s * alB
        z0 = z.copy(); z0[N - 1:] = (alA_cur + (alB_cur - alA_cur) * xr)[1:N]     # linear-twist guess for alpha, previous u
        sol = root(F, z0, method="lm", tol=1e-14)
        z = sol.x
    u, al = unpack(z)
    phi, Gphi, Gal, V, PrR, PtR = m.solid_forces(u, al, want_energy=True)
    res = float(np.max(np.abs(F(z))))
    return dict(R=m.R.copy(), u=u, alpha=al, T_A=float(Gal[0]), T_B=float(-Gal[N]), V=float(V), residual=res, model=m, Ns=N,
                PrR=PrR, PtR=PtR)


def eval_field(m, vals, R):
    """Spectral (Lagrange/GLL) evaluation of a nodal field at radii R (polynomial fit of degree Ns through the GLL nodes)."""
    from numpy.polynomial import legendre as L
    x = 2 * (m.R - A) / (B - A) - 1
    c = L.legfit(x, vals, m.Ns)
    return L.legval(2 * (np.asarray(R) - A) / (B - A) - 1, c)


# ------------------------------------------------------------------ RT-B: independent exact fixed-geometry fluid reference
def cheb(N):
    k = np.arange(N + 1)
    x = np.cos(np.pi * k / N)
    c = np.ones(N + 1); c[0] = c[-1] = 2; c *= (-1.0) ** k
    X = np.tile(x, (N + 1, 1)).T
    dX = X - X.T
    D = np.outer(c, 1 / c) / (dX + np.eye(N + 1))
    D -= np.diag(D.sum(axis=1))
    # Clenshaw-Curtis weights
    th = np.pi * k / N
    w = np.zeros(N + 1); v = np.ones(N - 1)
    if N % 2 == 0:
        w[0] = w[N] = 1.0 / (N**2 - 1)
        for kk in range(1, N // 2):
            v -= 2 * np.cos(2 * kk * th[1:-1]) / (4 * kk**2 - 1)
        v -= np.cos(N * th[1:-1]) / (N**2 - 1)
    else:
        w[0] = w[N] = 1.0 / N**2
        for kk in range(1, (N - 1) // 2 + 1):
            v -= 2 * np.cos(2 * kk * th[1:-1]) / (4 * kk**2 - 1)
    w[1:-1] = 2 * v / N
    return x, D, w


def omega_table(path=os.path.join(ROOT, "reference_data", "rtr", "omega_B_table_dt1e-3.csv")):
    """Wall angular-velocity history omega_w(t) = d alpha_B/dt of the PROPOSAL'S coupled torsion reference (Radau 1e-11), tabulated every 1 ms; created on first use."""
    if not os.path.exists(path):
        tt = np.round(np.arange(0, 2001) * 1e-3, 10)
        ref = coupled_reference(tt)
        np.savetxt(path, np.c_[tt, ref["alpha_B"], ref["omega_B"]], header="t alpha_B omega_B (rt_reference_s4f637_30ae00eb, torsion member, Ns=20 Nf=32 Radau rtol 1e-11)", fmt="%.17e")
    d = np.loadtxt(path)
    return d[:, 0], d[:, 1], d[:, 2]


def fluid_exact(tt, N=40, rtol=1e-12, atol=1e-14, b=B, c=C):
    """Exact (spectral) fixed-geometry swirl for wall angle history alpha_w(t) = THETA0 S sin(3t): returns dict of arrays on tt.
    nodes r_k (k=0 at r=b), L, v = L/r, torque T = mu (L_r - 2L) at r=b (per radian per m), H_f, wall shear stress, centripetal pressure."""
    TT, _, OM = omega_table()
    x, D, w = cheb(N)                     # x from +1 (k=0) to -1: map r = b + (1 - x)/2 (c - b)  so k=0 -> r=b
    ell = c - b
    r = b + (1 - x) * 0.5 * ell
    Dr = D * (-2.0 / ell)                 # d/dr
    Drr = Dr @ Dr
    rowN = Dr[N]                          # r=c at k=N
    # interior unknowns k=1..N-1; boundaries L_0 = b^2 omega_w(t), L_N from Robin L_r(c) = 2 L(c)/c
    def Lfull(t, Li):
        L0 = b * b * np.interp(t, TT, OM)
        LN = -(rowN[1:N] @ Li + rowN[0] * L0) / (rowN[N] - 2.0 / c)
        return np.concatenate(([L0], Li, [LN]))

    def rhs(t, Li):
        L = Lfull(t, Li)
        Lr = Dr @ L; Lrr = Dr @ Lr
        return NU * (Lrr - Lr / r)[1:N]
    sol = solve_ivp(rhs, (0, tt[-1]), np.zeros(N - 1), method="Radau", rtol=rtol, atol=atol, t_eval=tt, max_step=1e-3)
    assert sol.success, sol.message
    Ls = np.array([Lfull(t, sol.y[:, i]) for i, t in enumerate(tt)])         # (nt, N+1)
    Lr0 = (Dr[0] @ Ls.T)
    T = MU * (b * Lr0 - 2 * Ls[:, 0])                                         # b^2 sigma_rt, sigma_rt = mu (L_r/b - 2 L/b^2)
    wr = 0.5 * ell * w                                                        # weights for int dr (x descending -> positive)
    H = RHO * (Ls * r * wr).sum(axis=1)
    # pressure at the wall: p(b) = -rho int_b^c L^2 / r^3 dr  (Q = 0)
    Iw = (Ls**2 / r**3 * wr).sum(axis=1)
    return dict(t=np.asarray(tt), r=r, L=Ls, v=Ls / r, T_B=T, sig_rt=T / b**2, H_f=H, p_b=-RHO * Iw, Dr=Dr, wr=wr)


def pressure_profile(ex, i):
    """p(r_k) at time index i: p = -rho int_r^c L^2/r'^3 dr' (cumulative CC quadrature through spectral interpolation)."""
    r = ex["r"]; L = ex["L"][i]
    f = L**2 / r**3
    # integrate with the Chebyshev antiderivative via the inverse of Dr on interior (use cumulative from r=c): solve Dr q = f with q(c)=0
    Dr = ex["Dr"]; N = len(r) - 1
    Mx = Dr.copy()
    Mx[N, :] = 0; Mx[N, N] = 1
    rhs = f.copy(); rhs[N] = 0
    q = np.linalg.solve(Mx, rhs)      # q(r) = int_c^r f  (negative of int_r^c)
    return RHO * q                    # q(r) = -int_r^c f  ->  p = rho q


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    tt = np.round(np.arange(0, 201) * 0.01, 10)
    omega_table()
    e40 = fluid_exact(tt, N=40); e56 = fluid_exact(tt, N=56, rtol=1e-13, atol=1e-15)
    d = np.max(np.abs(e40["T_B"] - e56["T_B"])) / np.max(np.abs(e56["T_B"]))
    print("fluid_exact N=40 vs N=56 torque rel max diff", d, "max|T_B| =", np.max(np.abs(e56["T_B"])), "max|H_f|=", np.max(np.abs(e56["H_f"])))
    ref = coupled_reference(tt)
    print("coupled ref: max|T_B|", np.max(np.abs(ref["T_B"])), " fixed-geometry (wall = alpha_A history) max|T_B|", np.max(np.abs(e56["T_B"])),
          " max|alpha_B|", np.max(np.abs(ref["alpha_B"])), "max|u_B|", np.max(np.abs(ref["u_B"])), "max|p_b|", np.max(np.abs(ref["p_b"])))
