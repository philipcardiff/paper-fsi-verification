import os, os as _os; HERE = _os.path.dirname(_os.path.abspath(__file__)); PKG = _os.path.dirname(HERE)   # package root; run directories of retained HPC outputs: env CAVITY_RUNS (default runs/cavityNext)
"""DERIVED from paper-data/cavity_plate_material_law/plate_laws.py (PR #20, unmodified original untouched): adds an optional tabulated follower pressure pfun(x_ref) (load factor ramps 0..1).
No smoothing of the tabulated pressure is applied by this module."""
import sys
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

E, NU = 500.0, 0.4
MU = E / (2 * (1 + NU))
LAM = E * NU / ((1 + NU) * (1 - 2 * NU))
L, H, X0 = 4.0, 0.4, 2.0


def first_pk(F2, law):
    F = np.eye(3)
    F[:2, :2] = F2
    if law == "stvk":
        Eg = 0.5 * (F.T @ F - np.eye(3))
        S = LAM * np.trace(Eg) * np.eye(3) + 2 * MU * Eg
        P = F @ S
    elif law == "neohookean":
        J = np.linalg.det(F)
        Fit = np.linalg.inv(F).T
        P = MU * (F - Fit) + LAM * np.log(J) * Fit
    elif law in ("nh_s4f", "nh_ansys"):
        # isochoric mu/2 (tr(bbar) - 3) + U(J), K = lambda + 2 mu/3
        K = LAM + 2 * MU / 3
        J = np.linalg.det(F)
        bbar = J ** (-2 / 3) * (F @ F.T)
        dev = bbar - np.trace(bbar) / 3 * np.eye(3)
        if law == "nh_s4f":
            Jp = 0.5 * K * (J * J - 1) / J      # U = K/4 (J^2 - 1 - 2 ln J)
        else:
            Jp = K * (J - 1)                     # U = K/2 (J - 1)^2 (Ansys, d = 2/K)
        sigma = MU / J * dev + Jp * np.eye(3)
        P = J * sigma @ np.linalg.inv(F).T
    elif law == "hencky":
        b = F @ F.T
        w, v = np.linalg.eigh(b)
        h = v @ np.diag(0.5 * np.log(w)) @ v.T
        tau = LAM * np.trace(h) * np.eye(3) + 2 * MU * h
        P = tau @ np.linalg.inv(F).T
    elif law == "linear":
        # small strain, for checking against beam theory
        eps = 0.5 * (F + F.T) - np.eye(3)
        P = LAM * np.trace(eps) * np.eye(3) + 2 * MU * eps
    return P[:2, :2]


# 1-D quadratic Lagrange basis on [-1, 1]
def b1(s):
    return np.array([0.5 * s * (s - 1), 1 - s * s, 0.5 * s * (s + 1)])


def db1(s):
    return np.array([s - 0.5, -2 * s, s + 0.5])


g = np.sqrt(3 / 5)
GP = np.array([-g, 0.0, g])
GW = np.array([5 / 9, 8 / 9, 5 / 9])


def run(nx, ny, law, p_total, nsteps=20, pfun=None):
    nnx, nny = 2 * nx + 1, 2 * ny + 1
    xs = np.linspace(X0, X0 + L, nnx)
    ys = np.linspace(-H, 0.0, nny)
    X = np.array([[x, y] for y in ys for x in xs])  # node (j, i) -> j*nnx+i
    nid = lambda i, j: j * nnx + i
    elems = []
    for ey in range(ny):
        for ex in range(nx):
            i0, j0 = 2 * ex, 2 * ey
            elems.append([nid(i0 + a, j0 + c) for c in range(3) for a in range(3)])
    elems = np.array(elems)
    ndof = 2 * len(X)

    # shape-function derivatives at Gauss points (reference element)
    dN = []
    for a, sa in enumerate(GP):
        for c, sc in enumerate(GP):
            Nx = np.outer(b1(sc), db1(sa)).ravel()  # d/dxi
            Ny = np.outer(db1(sc), b1(sa)).ravel()  # d/deta
            dN.append((np.vstack([Nx, Ny]), GW[a] * GW[c]))

    # top edge elements (follower pressure on y = 0)
    top = [nid(2 * ex + a, nny - 1) for ex in range(nx) for a in range(3)]
    top = np.array(top).reshape(nx, 3)

    fixed = [n for n in range(len(X))
             if abs(X[n, 0] - X0) < 1e-12 or abs(X[n, 0] - X0 - L) < 1e-12]
    fixed_dofs = np.array(sorted([2 * n + k for n in fixed for k in range(2)]))
    free = np.setdiff1d(np.arange(ndof), fixed_dofs)

    def elem_res(xe, Xe):
        r = np.zeros(18)
        for (dNr, w) in dN:
            J0 = dNr @ Xe  # 2x2, rows d/dxi,d/deta
            dNX = np.linalg.solve(J0, dNr)  # d/dX
            F = (dNX @ xe).T  # F_ij = dx_i/dX_j
            P = first_pk(F, law)
            r += (dNX.T @ P.T).ravel() * np.linalg.det(J0) * w
        return r

    def edge_res(xe, p, Xr=None):
        # force on the solid from pressure p on the current top face;
        # nodes ordered left to right, outward normal +y in the reference
        r = np.zeros(6)
        for s, w in zip(GP, GW):
            N, dNs = b1(s), db1(s)
            t = dNs @ xe  # current tangent dx/ds
            n = np.array([-t[1], t[0]])  # outward (up) normal times da/ds
            pl = p if pfun is None else p * pfun(float(N @ Xr[:, 0]))   # p = load factor; pfun(x_ref) = tabulated nodal pressure at the reference position of the Gauss point
            f = -pl * n
            r -= np.outer(N, f).ravel() * w  # residual = internal - external
        return r

    u = np.zeros(ndof)
    h = 1e-7
    for step in range(1, nsteps + 1):
        p = p_total * step / nsteps
        for it in range(30):
            x = X + u.reshape(-1, 2)
            R = np.zeros(ndof)
            rows, cols, vals = [], [], []
            for e in elems:
                dofs = np.array([[2 * n, 2 * n + 1] for n in e]).ravel()
                xe, Xe = x[e], X[e]
                r0 = elem_res(xe, Xe)
                R[dofs] += r0
                Ke = np.zeros((18, 18))
                xf = xe.ravel().copy()
                for k in range(18):
                    xp = xf.copy(); xp[k] += h
                    Ke[:, k] = (elem_res(xp.reshape(9, 2), Xe) - r0) / h
                rows.append(np.repeat(dofs, 18)); cols.append(np.tile(dofs, 18))
                vals.append(Ke.ravel())
            for ed in top:
                dofs = np.array([[2 * n, 2 * n + 1] for n in ed]).ravel()
                xe = x[ed]
                r0 = edge_res(xe, p, X[ed])
                R[dofs] += r0
                Ke = np.zeros((6, 6))
                xf = xe.ravel().copy()
                for k in range(6):
                    xp = xf.copy(); xp[k] += h
                    Ke[:, k] = (edge_res(xp.reshape(3, 2), p, X[ed]) - r0) / h
                rows.append(np.repeat(dofs, 6)); cols.append(np.tile(dofs, 6))
                vals.append(Ke.ravel())
            K = sp.csr_matrix((np.concatenate(vals),
                               (np.concatenate(rows), np.concatenate(cols))),
                              shape=(ndof, ndof))
            res = np.linalg.norm(R[free])
            if it == 0:
                res0 = max(res, 1e-30)
            if res < 1e-10 * max(res0, 1.0) and it > 0:
                break
            du = spla.spsolve(K[free][:, free].tocsc(), -R[free])
            u[free] += du
    # point A: top surface, x = 4
    iA = nid(nnx // 2, nny - 1)
    iM = nid(nnx // 2, (nny - 1) // 2)
    return u[2 * iA + 1], u[2 * iA], u[2 * iM + 1]


if __name__ == "__main__":
    nx, ny = int(sys.argv[1]), int(sys.argv[2])
    q = float(sys.argv[3]) if len(sys.argv) > 3 else 5.26 / L
    laws = sys.argv[4].split(",") if len(sys.argv) > 4 else \
        ["stvk"]
    for law in laws:
        uy, ux, uym = run(nx, ny, law, q, nsteps=1 if law == "linear" else 20)
        print(f"{law:11s} nx={nx:3d} ny={ny:2d} q={q:.4f}  "
              f"uy_A={uy:.5f}  ux_A={ux:.2e}  uy_mid={uym:.5f}", flush=True)
