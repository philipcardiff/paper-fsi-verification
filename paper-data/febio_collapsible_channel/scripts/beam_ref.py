#!/usr/bin/env python3
"""Independent re-derivation of the static clamped-wall reference (CC-S1): geometrically nonlinear plane-strain Kirchhoff-Love beam under a uniform follower
(constant-magnitude, normal) pressure, solved as the stationary point (Newton, complex-step gradient) of the total potential Pi = U + p * int(y dx) (a uniform normal pressure on a clamped-clamped beam is conservative) with cubic Hermite elements; U uses the
Green-Lagrange membrane strain gamma = (|r'|^2 - 1)/2 and the curvature kappa = (x' y'' - y' x'')/|r'| (a Kirchhoff-Love beam of the class used by the oomph-lib reference; the oomph-lib beam source was NOT read for this re-derivation), no pre-stress.
Constants from the tutorial: L = 10, h = 0.05, E = 4.55e8, nu = 0.3 (E_eff = E/(1-nu^2) = 5.0e8), p_ext = 200 Pa.  Light (one core, a few seconds).
usage: beam_ref.py OUT.json [10,20,40]  (element counts; the 40-element run takes about 1.5 min on one core)"""
import sys, json, numpy as np
from scipy.optimize import root
L, H, E, NU, P = 10.0, 0.05, 4.55e8, 0.3, 200.0
EEFF = E / (1 - NU ** 2); EA = EEFF * H; EI = EEFF * H ** 3 / 12
GP, GW = np.polynomial.legendre.leggauss(6)
def hermite(s, le):   # basis and derivatives on [0, le], s in [0, 1] local xi
    xi = s
    N = np.array([1 - 3 * xi ** 2 + 2 * xi ** 3, le * (xi - 2 * xi ** 2 + xi ** 3), 3 * xi ** 2 - 2 * xi ** 3, le * (-xi ** 2 + xi ** 3)])
    d1 = np.array([-6 * xi + 6 * xi ** 2, le * (1 - 4 * xi + 3 * xi ** 2), 6 * xi - 6 * xi ** 2, le * (-2 * xi + 3 * xi ** 2)]) / le
    d2 = np.array([-6 + 12 * xi, le * (-4 + 6 * xi), 6 - 12 * xi, le * (-2 + 6 * xi)]) / le ** 2
    return N, d1, d2
def newton(f, x, tol=1e-6, maxit=60):
    for it in range(maxit):
        r = f(x); nr = np.linalg.norm(r)
        if nr < tol * (1 + 0): return x
        n = len(x); Jm = np.zeros((n, n))
        for j in range(n):
            h = 1e-7 * max(1.0, abs(x[j])); xp = x.copy(); xp[j] += h; Jm[:, j] = (f(xp) - r) / h
        dx = np.linalg.solve(Jm, -r); a = 1.0
        while np.linalg.norm(f(x + a * dx)) > nr and a > 1e-4: a /= 2
        x = x + a * dx
    raise RuntimeError("no convergence, |r|=%g" % nr)
def solve(ne, nsteps=20):
    le = L / ne; xi = (GP + 1) / 2; w = GW / 2 * le
    B = [hermite(s, le) for s in xi]
    # dofs per node: x, x', y, y' (x' = dx/dS0).  q shape (ne+1, 4)
    X0 = np.linspace(0, L, ne + 1); q0 = np.zeros((ne + 1, 4)); q0[:, 0] = X0; q0[:, 1] = 1.0
    free = np.ones((ne + 1, 4), bool); free[0, :] = False; free[-1, :] = False
    def pot(qe_all, p):      # element potentials (complex-safe); qe_all (ne, 8): x0, x0', x1, x1', y0, y0', y1, y1'
        Pi = 0
        for (N, d1, d2), wk in zip(B, w):
            x1 = qe_all[:, 0] * N[0] + qe_all[:, 1] * N[1] + qe_all[:, 2] * N[2] + qe_all[:, 3] * N[3]
            xs = qe_all[:, 0] * d1[0] + qe_all[:, 1] * d1[1] + qe_all[:, 2] * d1[2] + qe_all[:, 3] * d1[3]
            xss = qe_all[:, 0] * d2[0] + qe_all[:, 1] * d2[1] + qe_all[:, 2] * d2[2] + qe_all[:, 3] * d2[3]
            y1 = qe_all[:, 4] * N[0] + qe_all[:, 5] * N[1] + qe_all[:, 6] * N[2] + qe_all[:, 7] * N[3]
            ys = qe_all[:, 4] * d1[0] + qe_all[:, 5] * d1[1] + qe_all[:, 6] * d1[2] + qe_all[:, 7] * d1[3]
            yss = qe_all[:, 4] * d2[0] + qe_all[:, 5] * d2[1] + qe_all[:, 6] * d2[2] + qe_all[:, 7] * d2[3]
            g = 0.5 * (xs * xs + ys * ys - 1); sl = np.sqrt(xs * xs + ys * ys); kap = (xs * yss - ys * xss) / sl
            Pi = Pi + np.sum(wk * (0.5 * EA * g * g + 0.5 * EI * kap * kap + p * ys * 0 + p * y1 * xs))   # + p * int y dx : pressure from above does work when y falls
        return Pi
    def gather(q):
        qe = np.zeros((ne, 8))
        for e in range(ne):
            a, b = q[e], q[e + 1]; qe[e] = [a[0], a[1], b[0], b[1], a[2], a[3], b[2], b[3]]
        return qe
    def grad(qf, p):   # complex-step gradient of the total potential wrt the free dofs, element-wise
        q = q0.copy(); q[free] = qf; qe = gather(q); g = np.zeros_like(q)
        dq = np.zeros((ne, 8)); 
        for c in range(8):
            qc = qe.astype(complex); qc[:, c] += 1e-30j
            Pi_e = pot_elem(qc, p); dq[:, c] = Pi_e.imag / 1e-30
        for e in range(ne):
            for loc, (nn, dd) in enumerate([(e, 0), (e, 1), (e + 1, 0), (e + 1, 1), (e, 2), (e, 3), (e + 1, 2), (e + 1, 3)]): g[nn, dd] += dq[e, loc]
        return g[free]
    def pot_elem(qe_all, p):
        out = np.zeros(ne, complex)
        for (N, d1, d2), wk in zip(B, w):
            f = lambda c0, ba: qe_all[:, c0] * ba[0] + qe_all[:, c0 + 1] * ba[1] + qe_all[:, c0 + 2] * ba[2] + qe_all[:, c0 + 3] * ba[3]
            xs, xss, ys, yss, y1 = f(0, d1), f(0, d2), f(4, d1), f(4, d2), f(4, N)
            g = 0.5 * (xs * xs + ys * ys - 1); sl = np.sqrt(xs * xs + ys * ys); kap = (xs * yss - ys * xss) / sl
            out = out + wk * (0.5 * EA * g * g + 0.5 * EI * kap * kap + p * y1 * xs)
        return out
    qf = q0[free].copy(); hist = []
    for k in range(1, nsteps + 1):
        p = P * k / nsteps
        qf = newton(lambda v: grad(v, p), qf)
        q = q0.copy(); q[free] = qf
        hist.append((p, q[ne // 2, 2] - 0.0))
    q = q0.copy(); q[free] = qf
    mid = q[ne // 2, 2]; qt = None
    # quarter points through Hermite interpolation
    def val(Xp):
        e = min(int(Xp / le), ne - 1); s = (Xp - e * le) / le; N = hermite(s, le)[0]; a, b = q[e], q[e + 1]
        return N @ np.array([a[2], a[3], b[2], b[3]])
    return dict(ne=ne, mid=float(mid), q1=float(val(2.5)), q3=float(val(7.5)), curve=[(float(a), float(b)) for a, b in hist])
if __name__ == "__main__":
    nes = [int(v) for v in (sys.argv[2] if len(sys.argv) > 2 else "10,20,40").split(",")]
    res = [solve(ne) for ne in nes]
    for r in res: print(r["ne"], "mid", r["mid"], "q1", r["q1"], "q3", r["q3"])
    ext = {}
    if len(res) >= 3:   # Richardson with the observed order from the last three (element count doubles)
        a, b, c = [r["mid"] for r in res[-3:]]; ratio = (b - a) / (c - b); p = np.log2(ratio); ext = dict(observed_ratio=float(ratio), observed_order=float(p), richardson_mid=float(c + (c - b) / (ratio - 1)), richardson_mid_p2=float(c + (c - b) / 3.0))
        print("ratio", ratio, "order", p, "Richardson (observed order)", ext["richardson_mid"], " (p=2)", ext["richardson_mid_p2"])
        ext["rel_to_oomph_static_mid"] = ext["richardson_mid"] / -0.14362832 - 1
    out = sys.argv[1] if len(sys.argv) > 1 else "beam_ref.json"
    json.dump(dict(refs=dict(mid=-0.14362832, q=-0.09645925), E_eff=EEFF, EA=EA, EI=EI, runs=res, extrapolation=ext), open(out, "w"), indent=1)
