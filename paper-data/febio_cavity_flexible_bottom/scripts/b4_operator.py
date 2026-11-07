#!/usr/bin/env python3
import os, os as _os; HERE = _os.path.dirname(_os.path.abspath(__file__)); PKG = _os.path.dirname(HERE)   # package root; run directories of retained HPC outputs: env CAVITY_RUNS (default runs/cavityNext)
"""B4: discrete operator of FEBio v4.13's fluid-FSI element (FEFluidFSIDomain3D, equal-order hex8 w / e_f), linearised about the rest state of a FIXED mesh (w0 = 0, ef0 = 0, Lf = 0),
assembled exactly from the source blocks (FEFluidFSIDomain3D.cpp ElementStiffness + ElementMassMatrix) for z-uniform plane-strain fields (hex8 one layer, wz = 0, 8-pt Gauss -> 2-D Q1, 2x2 Gauss).
 Kww = int gradN_i . cv . gradN_j  (+ rho c N_i N_j I, c = alpham/(alphaf gamma dt));  kwJ = dp N_i gradN_j (dp = -K);  kJw = gradN_i N_j;  kJJ = c N_i N_j  (source l.523-535, 740, 536).
 In the pressure variable q = -K ef:   [ A   G ] [w]   [ b_w ]       G_ij = int N_i grad N_j  ;  kJw = G^T (identity checked numerically);  kJJ -> -(c/K) M.
 (A)(B)(C)(D) analyses + backward-Euler transient emulation.  usage: b4_operator.py OUTDIR"""
import sys, os, json, numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla, scipy.linalg as sl
R = PKG; sys.path.insert(0, HERE)
import make_cavity_fluid_feb as GEO
MU, RHO, K = 0.01, 1.0, 1e6; DT = 2.5; RHOI = 0.0; ALPHAF = 1.0 / (1 + RHOI); ALPHAM = (2 - RHOI) / (1 + RHOI); GAMMA = 0.5 + ALPHAM - ALPHAF; BETA = (1 + ALPHAM - ALPHAF) ** 2 / 4
CC = ALPHAM / (ALPHAF * GAMMA * DT)     # source c (l.463)
GP = [(-1 / np.sqrt(3), -1 / np.sqrt(3)), (1 / np.sqrt(3), -1 / np.sqrt(3)), (1 / np.sqrt(3), 1 / np.sqrt(3)), (-1 / np.sqrt(3), 1 / np.sqrt(3))]
def N(r, s): return np.array([(1 - r) * (1 - s), (1 + r) * (1 - s), (1 + r) * (1 + s), (1 - r) * (1 + s)]) / 4
def dN(r, s): return np.array([-(1 - s), (1 - s), (1 + s), -(1 + s)]) / 4, np.array([-(1 - r), -(1 + r), (1 + r), (1 - r)]) / 4
def assemble(X, quads, mu=MU, rho=RHO, c=CC, bulkvisc=0.0):
    n = len(X); rows_A, cols_A, val_A = [], [], []; Gi, Gj, Gv = [], [], []; Mi, Mj, Mv = [], [], []; lam = bulkvisc - 2 * mu / 3
    for q in quads:
        x = X[list(q)]; ke = np.zeros((8, 8)); ge = np.zeros((8, 4)); me = np.zeros((4, 4))
        for (r, s) in GP:
            dr, ds = dN(r, s); J = np.array([[dr @ x[:, 0], dr @ x[:, 1]], [ds @ x[:, 0], ds @ x[:, 1]]]); dJ = np.linalg.det(J); Ji = np.linalg.inv(J)
            grad = (Ji @ np.vstack([dr, ds])).T      # 4x2: gradN_a
            H = N(r, s)
            for a in range(4):
                for b in range(4):
                    ga, gb = grad[a], grad[b]
                    Kab = lam * np.outer(ga, gb) + mu * ((ga @ gb) * np.eye(2) + np.outer(gb, ga))     # v.T.v of c = lam I(x)I + 2mu I4s  (index order i,k = a_i b_k ...)
                    ke[2 * a:2 * a + 2, 2 * b:2 * b + 2] += (Kab + rho * c * H[a] * H[b] * np.eye(2)) * dJ
                    ge[2 * a:2 * a + 2, b] += H[a] * gb * dJ       # momentum row a (vector), pressure col b:  int N_a grad N_b
                    me[a, b] += H[a] * H[b] * dJ
        idx = np.array([[2 * i, 2 * i + 1] for i in q]).ravel()
        for i in range(8):
            for j in range(8): rows_A.append(idx[i]); cols_A.append(idx[j]); val_A.append(ke[i, j])
            for j in range(4): Gi.append(idx[i]); Gj.append(q[j]); Gv.append(ge[i, j])
        for i in range(4):
            for j in range(4): Mi.append(q[i]); Mj.append(q[j]); Mv.append(me[i, j])
    A = sp.csr_matrix((val_A, (rows_A, cols_A)), shape=(2 * n, 2 * n)); G = sp.csr_matrix((Gv, (Gi, Gj)), shape=(2 * n, n)); M = sp.csr_matrix((Mv, (Mi, Mj)), shape=(n, n))
    return A, G, M
_ORIG = (GEO.V, GEO.BLOCKS, GEO.ARCS)
def geometry(kind, f=1):
    GEO.V, GEO.BLOCKS, GEO.ARCS = _ORIG
    if kind == "straight":
        GEO.V = [(0, 1), (14, 1), (14, 2), (0, 2)]; GEO.BLOCKS = [((0, 1, 2, 3), 140 * f, 10 * f)]; GEO.ARCS = {}
    nodes, quads = GEO.mesh(f); X = np.array(nodes); eps = 1e-6; cnt = {}
    for q in quads:
        for e in ((q[0], q[1]), (q[1], q[2]), (q[2], q[3]), (q[3], q[0])): cnt.setdefault(frozenset(e), []).append(e)
    bedges = [v[0] for v in cnt.values() if len(v) == 1]; bn = {i for e in bedges for i in e}
    inlet = sorted(i for i in bn if abs(X[i, 0]) < eps); outlet = sorted(i for i in bn if abs(X[i, 0] - 14) < eps)
    iface = sorted((i for i in bn if abs(X[i, 1]) < eps), key=lambda i: X[i, 0]) if kind == "cavity" else []
    walls = sorted(bn - set(inlet) - set(outlet) - set(iface)); in_corner = [i for i in inlet if abs(X[i, 1] - 1) < eps or abs(X[i, 1] - 2) < eps]; out_corner = [i for i in outlet if abs(X[i, 1] - 1) < eps or abs(X[i, 1] - 2) < eps]
    return dict(X=X, quads=quads, inlet=inlet, outlet=outlet, iface=iface, walls=walls, in_corner=in_corner, out_corner=out_corner)
def bcs(g, pins=None, rim_free=False):
    n = len(g["X"]); wfix = set(g["walls"]) | set(g["iface"]) | set(g["in_corner"]) | set(g["out_corner"]) | set(g["inlet"])   # inlet: all w prescribed (nodal normal velocity)
    pfix = set(g["outlet"])
    if not rim_free: pfix |= set(g["in_corner"])    # rim pressure prescribed (value = mean of the other inlet nodes: treated as a Dirichlet constraint)
    if pins is not None: pfix |= {i for i in g["iface"] if abs(g["X"][i, 0] - 2) < 1e-6 or abs(g["X"][i, 0] - 6) < 1e-6}
    wfree = np.array([2 * i + k for i in range(n) if i not in wfix for k in (0, 1)]); pfree = np.array([i for i in range(n) if i not in pfix])
    return wfree, pfree, wfix, pfix
def inlet_profile(g):
    X = g["X"]; ins = [i for i in g["inlet"] if i not in g["in_corner"]]; ys = np.array(sorted(X[g["inlet"], 1])); phi = {}
    yy = sorted(g["inlet"], key=lambda i: X[i, 1]); y = X[yy, 1]; ph = 6 * (y - 1) * (2 - y)           # nodal-exact parabola; normalise the discrete (trapezoid) mean to 1 as FEParabolicMap does
    mean = np.trapz(ph, y) / (y[-1] - y[0]); ph = ph / mean; return {i: p for i, p in zip(yy, ph)}
def inlet_load(g):
    """FEFluidNormalVelocity::LoadVector: dilatation-equation load  int N_i VN da (VN = -phi inward) -> contribution to the continuity row (rows: kJw w = load)."""
    X = g["X"]; prof = inlet_profile(g); yy = sorted(g["inlet"], key=lambda i: X[i, 1]); b = np.zeros(len(X))
    for a, c in zip(yy[:-1], yy[1:]):
        L = abs(X[c, 1] - X[a, 1]); pa, pc = prof[a], prof[c]
        b[a] += L * (2 * pa + pc) / 6; b[c] += L * (pa + 2 * pc) / 6     # int N_i phi dy (linear phi); sign: flux in = +phi.  kJw w = -(flux in)  ->  R side below
    return b, prof
def rows_above_floor(g):
    quads = g["quads"]; X = g["X"]; e2q = {}
    for qi, q in enumerate(quads):
        for e in ((q[0], q[1]), (q[1], q[2]), (q[2], q[3]), (q[3], q[0])): e2q.setdefault(frozenset(e), []).append(qi)
    rows = {}  # floor node -> [row0,row1,row2,...]
    fl = set(g["iface"])
    def up(a, b, qi):    # given floor-side edge (a,b) of quad qi return the opposite edge
        q = list(quads[qi]); i = q.index(a)
        for (u, v) in ((q[0], q[1]), (q[1], q[2]), (q[2], q[3]), (q[3], q[0])):
            pass
        return None
    # vertical chains: walk from each floor node, choose in the first quad the neighbour node that lies above
    nbr = {}
    for q in quads:
        for k in range(4): nbr.setdefault(q[k], set()).update([q[(k + 1) % 4], q[(k - 1) % 4]])
    chains = {}
    for a in g["iface"]:
        ch = [a]
        for _ in range(7):
            cur = ch[-1]; cand = [m for m in nbr[cur] if X[m, 1] > X[cur, 1] + 1e-9 and abs(X[m, 0] - X[cur, 0]) < 0.35 * max(1e-9, (X[m, 1] - X[cur, 1])) + 1e-6]
            if not cand: break
            ch.append(min(cand, key=lambda m: X[m, 1]))
        chains[a] = ch
    return chains
def run_all(outdir):
    res = {}; os.makedirs(outdir, exist_ok=True)
    for kind, f, pins in (("straight", 1, None), ("cavity", 1, None), ("cavity", 1, 1.15)):
        g = geometry(kind, f); X = g["X"]; n = len(X); A0, G, M = assemble(X, g["quads"], c=0.0); At, _, _ = assemble(X, g["quads"], c=CC)
        wfree, pfree, wfix, pfix = bcs(g, pins)
        key = f"{kind}_f{f}_" + ("pin" if pins else "nopin"); r = dict(n_nodes=n, n_wfree=len(wfree), n_pfree=len(pfree), n_pfixed=len(pfix))
        # identity check kJw = G^T (source: kJw = gradN_i H_j ; here assembled independently as int gradN_i N_j)
        CJ = sp.lil_matrix((n, 2 * n))
        for q in g["quads"]:
            xq = X[list(q)]
            for (rr, ss) in GP:
                dr, ds = dN(rr, ss); J = np.array([[dr @ xq[:, 0], dr @ xq[:, 1]], [ds @ xq[:, 0], ds @ xq[:, 1]]]); dJ = np.linalg.det(J); grad = (np.linalg.inv(J) @ np.vstack([dr, ds])).T; H = N(rr, ss)
                for a in range(4):
                    for b in range(4): CJ[q[a], 2 * q[b]] += grad[a, 0] * H[b] * dJ; CJ[q[a], 2 * q[b] + 1] += grad[a, 1] * H[b] * dJ
        r["max_abs(kJw - G^T)"] = float(abs(CJ.tocsr() - G.T).max())
        Gf = G[wfree][:, pfree]; Mp = M[pfree][:, pfree].tocsr(); res_local = {}
        for label, A in (("steady", A0), ("transient_dt2.5", At)):
            Aff = A[wfree][:, wfree].tocsc(); lu = sla.splu(Aff); X_ = lu.solve(Gf.toarray()); S = Gf.T @ X_; S = 0.5 * (S + S.T)
            ev, V = sl.eigh(S, Mp.toarray(), eigvals_only=False, subset_by_index=[0, 15])
            # parity / localisation metrics of the lowest vectors
            pos = X[pfree]; info = []
            for k in range(12):
                v = V[:, k]; w2 = v ** 2 * np.array(Mp.diagonal()); w2 /= w2.sum(); pr = 1.0 / np.sum(w2 ** 2)
                # y-parity indicator: correlation with the nodal parity signs (-1)^{round(y/h_local)} is geometry dependent -> use local sign-alternation fraction along y-edges
                info.append(dict(k=k, lam=float(ev[k]), participation_nodes=float(pr), centroid_xy=[float((w2 * pos[:, 0]).sum()), float((w2 * pos[:, 1]).sum())], frac_weight_in_cavity_box=float(w2[(pos[:, 0] > 2 - 1e-6) & (pos[:, 0] < 6 + 1e-6) & (pos[:, 1] < 0.75)].sum()),
                                 edge_sign_alt_frac=float(edge_alt(g, pfree, v))))
            res_local[label] = dict(lowest16_generalised_eig_S_q_eq_lam_M_q=[float(x) for x in ev], lowest12=info, eps_CK_M=float(CC / K))
            r["median_eig_S"] = None
        r["Schur"] = res_local; res[key] = r
        # parity test vectors (rows from the floor in the cavity; (-1)^j along y in the channel)
        if True:
            if kind == "straight":
                jidx = np.round((X[:, 1] - 1) / 0.1).astype(int); xi = np.round(X[:, 0] / 0.1).astype(int)
            else:
                chains = rows_above_floor(g); jidx = np.full(n, -1); xi = np.round(X[:, 0] / 0.1).astype(int)
                for a, ch in chains.items():
                    for j, m in enumerate(ch): jidx[m] = j
            def evalv(name, q):
                qf = q[pfree]; qfull = np.zeros(n); qfull[pfree] = qf; mom = Gf @ qf; nrm = np.sqrt(qf @ (Mp @ qf)); wl = np.linalg.norm(mom) / nrm
                # momentum residual rows: G q at free velocity rows; continuity residual with w = 0: -(c/K) M q
                return dict(norm_q_M=float(nrm), momentum_residual_norm_over_norm=float(wl), continuity_residual_norm_over_norm_w0=float(CC / K * np.sqrt((Mp @ qf) @ (Mp @ qf)) / nrm),
                            Rayleigh_steady=float(qf @ (Gf.T @ sla.splu(A0[wfree][:, wfree].tocsc()).solve(Gf @ qf)) / (qf @ (Mp @ qf))) if True else None)
            tests = {}
            if kind == "straight" or kind == "cavity":
                ok = jidx >= 0 if kind == "cavity" else np.ones(n, bool)
                q_y = np.where(ok, (-1.0) ** jidx, 0.0); q_x = (-1.0) ** xi * (np.abs(X[:, 0] - 14) > 0.2); q_c = q_y * (-1.0) ** xi
                q_ys = q_y * np.minimum(1, (14 - X[:, 0]) / 3.0) * (X[:, 1] < 2.5)   # y-parity with a smooth envelope that vanishes at the outlet
                q_smooth = (14 - X[:, 0]) / 14.0 * 1.4        # a typical smooth channel pressure drop (Pa)
                q_smooth2 = 1.4 * np.cos(np.pi * X[:, 1] / 3) * (14 - X[:, 0]) / 14
                for nm, q in (("y_parity_envelope", q_ys), ("y_parity_raw", q_y), ("x_parity", q_x), ("checker_xy", q_c), ("smooth_drop", q_smooth), ("smooth_drop_cos_y", q_smooth2)):
                    q = np.where(np.isin(np.arange(n), pfree), q, 0.0)
                    for label, A in (("steady", A0), ("dt2.5", At)):
                        Aff = A[wfree][:, wfree].tocsc(); lu = sla.splu(Aff); qf = q[pfree]; mom = Gf @ qf; nrm = np.sqrt(qf @ (Mp @ qf)); Rq = float(qf @ (Gf.T @ lu.solve(mom)) / (qf @ (Mp @ qf)))
                        tests[f"{nm}|{label}"] = dict(momentum_residual_norm_over_Mnorm=float(np.linalg.norm(mom) / nrm), Rayleigh_quotient_S=Rq, ratio_to_eps_CK=Rq / (CC / K))
            r["test_vectors"] = tests
        print(key, "done", {k: v for k, v in r.items() if k in ("n_wfree", "n_pfree")}); sys.stdout.flush()
    json.dump(res, open(outdir + "/b4_operator_spectra.json", "w"), indent=1)
def edge_alt(g, pfree, v):
    """fraction of mesh edges (both ends free) along which the vector changes sign: 1 for a checkerboard, ~0.0-0.1 for smooth"""
    full = np.zeros(len(g["X"])); full[pfree] = v; fr = set(pfree.tolist()); c = t = 0
    for q in g["quads"]:
        for a, b in ((q[0], q[1]), (q[1], q[2]), (q[2], q[3]), (q[3], q[0])):
            if a in fr and b in fr: t += 1; c += (full[a] * full[b] < 0)
    return c / max(t, 1)
if __name__ == "__main__": run_all(sys.argv[1])
