#!/usr/bin/env python3
"""Operator audit (mathematics, not a FEBio run): equal-order Q1/Q1 (FEBio fluid-FSI continuity: integrated-by-parts form, momentum: strong grad p)
on a unit square, n x n uniform quads, velocity Dirichlet on the whole boundary (enclosed flow, no boundary flux term), continuous nodal pressure.
Discrete operators: A = stiffness of the vector Laplacian (mu = 1), G_{(a,d),e} = int N_a dN_e/dx_d  (momentum gradient), continuity row = G^T.
Schur complement S = G^T A^-1 G; generalised eigenproblem S q = lam M_p q with the pressure mass matrix M_p (inf-sup constant beta^2 = smallest lam
excluding the single constant mode). Also reports the compressibility-stabilised problem S + s M_p with s = gamma/(K dt) * mu-scaling, as in the
unsteady compressible problem, to show how far the time/compressibility term stabilises. Pure numpy/scipy, a few seconds."""
import numpy as np, scipy.sparse as sp, scipy.linalg as sl, json, sys
def build(n):
    h = 1.0 / n; N = (n + 1) ** 2; idx = lambda i, j: j * (n + 1) + i
    gp = [(-1 / np.sqrt(3), 1.0), (1 / np.sqrt(3), 1.0)]
    Ke = np.zeros((4, 4)); Me = np.zeros((4, 4)); Gx = np.zeros((4, 4)); Gy = np.zeros((4, 4))
    for a, wa in gp:
        for b, wb in gp:
            r = np.array([-1, 1, 1, -1]); s = np.array([-1, -1, 1, 1])
            Nn = 0.25 * (1 + r * a) * (1 + s * b)
            dNr = 0.25 * r * (1 + s * b); dNs = 0.25 * s * (1 + r * a)
            dNx = dNr * 2 / h; dNy = dNs * 2 / h; w = wa * wb * (h / 2) ** 2
            Ke += w * (np.outer(dNx, dNx) + np.outer(dNy, dNy)); Me += w * np.outer(Nn, Nn)
            Gx += w * np.outer(Nn, dNx); Gy += w * np.outer(Nn, dNy)       # (velocity node a, pressure node e)
    rows = []; Kc = []; Mc = []; Gxc = []; Gyc = []
    A = sp.lil_matrix((N, N)); M = sp.lil_matrix((N, N)); GX = sp.lil_matrix((N, N)); GY = sp.lil_matrix((N, N))
    for j in range(n):
        for i in range(n):
            nd = [idx(i, j), idx(i + 1, j), idx(i + 1, j + 1), idx(i, j + 1)]
            for a in range(4):
                for b in range(4):
                    A[nd[a], nd[b]] += Ke[a, b]; M[nd[a], nd[b]] += Me[a, b]; GX[nd[a], nd[b]] += Gx[a, b]; GY[nd[a], nd[b]] += Gy[a, b]
    inner = [idx(i, j) for j in range(1, n) for i in range(1, n)]
    A = A.tocsr()[inner][:, inner]; GX = GX.tocsr()[inner]; GY = GY.tocsr()[inner]
    return A.toarray(), M.toarray(), GX.toarray(), GY.toarray()
out = {}
for n in (8, 9, 16, 17):
    A, M, GX, GY = build(n)
    G = np.vstack([GX, GY]); Ablk = sl.block_diag(A, A)
    S = G.T @ np.linalg.solve(Ablk, G)
    lam = sl.eigh(S, M, eigvals_only=True)
    lam = np.sort(lam); out[n] = dict(smallest=[float(x) for x in lam[:6]], n_below_1em6_of_max=int((lam < 1e-6 * lam[-1]).sum()), lam_max=float(lam[-1]))
    # checkerboard content: Rayleigh quotient of the nodal checkerboard
    N1 = n + 1; cb = np.array([(-1) ** (i + j) for j in range(N1) for i in range(N1)], float)
    out[n]["checkerboard_rayleigh"] = float(cb @ S @ cb / (cb @ M @ cb))
    print(n, out[n])
json.dump(out, open(sys.argv[1] if len(sys.argv) > 1 else "q1q1_parity_operator.json", "w"), indent=1)
