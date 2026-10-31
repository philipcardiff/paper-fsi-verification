#!/usr/bin/env python3
"""Independent re-derivation of the RT-A and RT-B numbers from the stored FEBio text outputs.

Written for the archive package: it does NOT import the campaign analysis code (rtr_check_rta.py, rtr_analyze.py).  Own file parsers, own algebra, own quadrature,
own fixed-geometry fluid reference (second-order finite differences + Crank-Nicolson, a different method from the Chebyshev/Radau reference of rtr_reference.py).
It reads only the FEBio inputs (.feb: node coordinates, element connectivity) and the text logs.  The proposal's finite-strain static-torsion reference (vendored
rt_reference, hash-checked, also used by the campaign) is the one piece of reference code that is reused for RT-A twist; everything else is re-implemented.

usage: independent_check.py RAW_DIR INPUT_DIR OUT.json
  RAW_DIR   directory with RTA_law1c, RTA_rotcheck, RTA_twist_M{1,2,3}, RTB_fluid_M{1,2}, RTB_fluid_M1_bcshift sub-directories (text logs, plain or .xz)
  INPUT_DIR directory with the RTB_fluid_M1/M2 .feb files (node coordinates and connectivity are parsed from them)
"""
import sys, os, io, re, json, math, lzma
import numpy as np
from scipy.linalg import solve_banded

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
MU_S, K_S = 2.0e4 / (2 * 1.3), 2.0e4 / (3 * (1 - 0.6))       # E = 2e4, nu = 0.3
MU_F, RHO_F, DZ = 300.0, 1000.0, 0.1
TWO_PI_DZ = 2 * math.pi * DZ


def opentext(path):
    if os.path.exists(path):
        return open(path)
    if os.path.exists(path + ".xz"):
        return io.TextIOWrapper(lzma.open(path + ".xz"))
    raise FileNotFoundError(path)


def parse_log(path):
    """FEBio node_data text log -> (times[nt], list of 2-D arrays [rows, cols] (first col = node id))."""
    times, blocks, cur = [], [], None
    with opentext(path) as f:
        for ln in f:
            if ln.startswith("*Time"):
                times.append(float(ln.split("=")[1])); cur = []; blocks.append(cur)
            elif ln[0] == "*":
                continue
            elif ln.strip():
                cur.append(ln.split())
    return np.array(times), [np.array(b, dtype=float) for b in blocks]


def ramp_S(t):
    """C-infinity ramp f/(f+g), f = exp(-1/tau), g = exp(-1/(1-tau)), tau = clamp(t, 0.005, 0.995) (the expression written in the FEBio input); T_ramp = 1 s."""
    tau = min(max(t, 0.005), 0.995)
    f, g = math.exp(-1 / tau), math.exp(-1 / (1 - tau))
    return f / (f + g)


# ----------------------------------------------------------------------------------------------------------------------------------------- RT-A
def strain_energy(F):
    """The proposal's neo-Hookean density W = mu/2 (Ibar1 - 3) + K/4 (J^2 - 1 - 2 ln J), J = det F (3-D F, F33 = 1)."""
    J = np.linalg.det(F)
    I1b = J ** (-2.0 / 3.0) * np.trace(F @ F.T)
    return 0.5 * MU_S * (I1b - 3) + 0.25 * K_S * (J * J - 1 - 2 * math.log(J))


def cauchy_numeric(F, h=1e-6):
    """sigma = (1/J) dW/dF F^T by central differences of strain_energy (no closed-form stress used)."""
    P = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            Fp, Fm = F.copy(), F.copy(); Fp[i, j] += h; Fm[i, j] -= h
            P[i, j] = (strain_energy(Fp) - strain_energy(Fm)) / (2 * h)
    return P @ F.T / np.linalg.det(F)


def cauchy_closed(F):
    J = np.linalg.det(F); b = F @ F.T; bb = J ** (-2.0 / 3.0) * b
    return MU_S / J * (bb - np.trace(bb) / 3 * np.eye(3)) + 0.5 * K_S * (J * J - 1) / J * np.eye(3)


def rta_law(raw):
    meta = json.load(open(os.path.join(raw, "RTA_law1c", "model_meta.json")))
    Fm = np.array(meta["F"])
    T, B = parse_log(os.path.join(raw, "RTA_law1c", "stress.txt"))
    worst_closed = worst_num = scale = 0.0; nrows = 0
    for t, blk in zip(T, B):
        if t == 0:
            continue
        F = np.eye(3) + t * (Fm - np.eye(3))                         # t = 0.2 k, F_k = I + 0.2 k (Fm - I)
        s_c, s_n = cauchy_closed(F), cauchy_numeric(F)
        pick = lambda s: np.array([s[0, 0], s[1, 1], s[2, 2], s[0, 1], s[1, 2], s[0, 2]])
        sc, sn = pick(s_c), pick(s_n)
        scale = max(scale, np.abs(sc).max())
        for row in blk:
            nrows += 1
            worst_closed = max(worst_closed, np.abs(row[1:7] - sc).max()); worst_num = max(worst_num, np.abs(row[1:7] - sn).max())
            assert abs(row[7] - np.linalg.det(F)) < 1e-12
    return dict(n_steps=len(T) - 1, n_element_rows=nrows, max_abs_stress=scale, max_abs_diff_closed_form=worst_closed, rel_closed_form=worst_closed / scale,
                max_abs_diff_numeric_dW_dF=worst_num, rel_numeric_dW_dF=worst_num / scale, gate_rel=1e-8, PASS=bool(worst_closed / scale <= 1e-8))


def rta_rot(raw):
    T, B = parse_log(os.path.join(raw, "RTA_rotcheck", "nodes.txt"))
    worst = 0.0
    for t, blk in zip(T, B):
        a = 0.5 * ramp_S(t) * math.sin(3 * t) if 0 < t else 0.0
        for nid, x, y, ux, uy in blk[:, :5]:
            X, Y = x - ux, y - uy
            worst = max(worst, abs(ux - ((math.cos(a) - 1) * X - math.sin(a) * Y)), abs(uy - (math.sin(a) * X + (math.cos(a) - 1) * Y)))
    return dict(n_times=len(T), max_abs_node_diff_m=worst, gate_m=1e-10, PASS=bool(worst <= 1e-10))


def rta_twist(raw, ref_cache):
    """Static twist torque.  Reference T_A(t) from the vendored proposal energy (rtr_reference.static_twist; cached in ref_cache) at the logged times."""
    sys.path.insert(0, HERE)
    out = {}
    res = {}
    for m in (1, 2, 3):
        d = os.path.join(raw, f"RTA_twist_M{m}")
        T, Bi = parse_log(os.path.join(d, "inner.txt")); _, Bo = parse_log(os.path.join(d, "outer.txt"))
        TA = np.array([np.sum(b[:, 1] * b[:, 6] - b[:, 2] * b[:, 5]) / TWO_PI_DZ for b in Bi])      # sum (x Ry - y Rx) / (2 pi dz): reaction torque on the inner circle
        TO = np.array([np.sum(b[:, 1] * b[:, 6] - b[:, 2] * b[:, 5]) / TWO_PI_DZ for b in Bo])
        res[m] = (T, TA, TO)
    T = res[1][0]
    amp = 0.5
    cache = {}
    if os.path.exists(ref_cache):
        cache = {float(k): v for k, v in json.load(open(ref_cache)).items()}
    need = [t for t in T if t > 0 and float(t) not in cache]
    if need:
        import rtr_reference as R
        for t in need:
            S = R.ramp_S(t)[0] if t < 1 else 1.0
            cache[float(t)] = R.static_twist(amp * S * math.sin(3 * t), 0.0, Ns=24)["T_A"]
        json.dump({repr(k): v for k, v in cache.items()}, open(ref_cache, "w"), indent=0)
    Tref = np.array([cache[float(t)] for t in T if t > 0])
    for m in (1, 2, 3):
        T_, TA, TO = res[m]
        TA, TO = TA[T_ > 0], TO[T_ > 0]
        sgn = -1.0 if np.sum(TA * Tref) < 0 else 1.0                                    # FEBio reaction force = minus the load applied by the BC
        E = np.max(np.abs(sgn * TA - Tref)) / np.max(np.abs(Tref))
        out[f"M{m}"] = dict(E_max_TA=float(E), sign_of_reaction_vs_reference=sgn, peak_TA_reference=float(np.max(np.abs(Tref))), peak_TA_febio=float(np.max(np.abs(TA))),
                            inner_plus_outer_over_max_inner=float(np.max(np.abs(TA + TO)) / np.max(np.abs(TA))), n_logged_steps=int(len(TA)))
    out["order_M1_M2"] = math.log2(out["M1"]["E_max_TA"] / out["M2"]["E_max_TA"]); out["order_M2_M3"] = math.log2(out["M2"]["E_max_TA"] / out["M3"]["E_max_TA"])
    out["gate"] = "E_max_TA(M3) <= 1e-3 and order(M2,M3) >= 1.8"
    out["PASS"] = bool(out["M3"]["E_max_TA"] <= 1e-3 and out["order_M2_M3"] >= 1.8)
    return out


# ----------------------------------------------------------------------------------------------------------------------------------------- RT-B
def omega_table(path):
    d = np.loadtxt(path); return d[:, 0], d[:, 2]


def fd_exact(tab_t, tab_om, N=3200, dt=1e-3, shift=0.0, nu=MU_F / RHO_F, tend=2.0):
    """Independent fixed-geometry swirl: L = r v_theta, L_t = nu (L_rr - L_r/r) on [1,2], L(1,t) = omega_w(t - shift), L_r(2) = L(2) (zero shear stress).
    Second-order central differences (N cells) in r, BDF2 in time (step dt; implicit-Euler first step: the solution is zero to 1e-45 there).  BDF2 is L-stable; Crank-Nicolson
    was tried first and rejected because it rings at the wall (the tabulated omega has slope jumps every 1 ms), which pollutes the wall torque at the 1e-4 level.
    Returns r, times, L[nt, N+1], T_B(t) = mu (L_r(1) - 2 L(1)), H_f(t)."""
    r = np.linspace(1.0, 2.0, N + 1); h = 1.0 / N
    nt = int(round(tend / dt))
    om = lambda t: np.interp(max(t - shift, 0.0), tab_t, tab_om)
    n = N                                                  # unknowns L_1..L_N (L_0 prescribed); ghost node at N: L_{N+1} = L_{N-1} + 2 h L_N
    lo = nu * (1 / h ** 2 + 1 / (2 * h * r[1:]))           # coefficient of L_{i-1}
    up = nu * (1 / h ** 2 - 1 / (2 * h * r[1:]))           # coefficient of L_{i+1}
    dg = -2 * nu / h ** 2 * np.ones(n)
    lo_last = lo[-1] + up[-1]; dg[-1] += 2 * h * up[-1]

    def banded(c0):                                        # matrix c0 I - A
        ab = np.zeros((3, n)); ab[1] = c0 - dg
        ab[0, 1:] = -up[:-1]; ab[2, :-1] = -lo[1:]; ab[2, n - 2] = -lo_last
        return ab
    ab_be, ab_b2 = banded(1.0 / dt), banded(1.5 / dt)
    Ls = np.zeros((nt + 1, N + 1))
    Lm, Lc = np.zeros(n), np.zeros(n)                      # L^{n-1}, L^n
    for k in range(nt):
        L0n = om((k + 1) * dt)
        if k == 0:
            rhs = Lc / dt
            rhs[0] += lo[0] * L0n
            Ln = solve_banded((1, 1), ab_be, rhs)
        else:
            rhs = (4 * Lc - Lm) / (2 * dt)
            rhs[0] += lo[0] * L0n
            Ln = solve_banded((1, 1), ab_b2, rhs)
        Lm, Lc = Lc, Ln
        Ls[k + 1, 0] = L0n; Ls[k + 1, 1:] = Ln
    tt = np.arange(nt + 1) * dt
    Lr0 = (-3 * Ls[:, 0] + 4 * Ls[:, 1] - Ls[:, 2]) / (2 * h)
    T = MU_F * (Lr0 - 2 * Ls[:, 0])
    H = RHO_F * np.trapz(Ls * r, r, axis=1)
    return dict(r=r, t=tt, L=Ls, T=T, H=H)


def pressure_exact_at(L_row, r, rn):
    """p(r) = -rho int_r^2 L^2 / r'^3 dr' (trapezoid on the fine grid, evaluated at the radii rn which are grid nodes)."""
    f = L_row ** 2 / r ** 3
    cum = np.concatenate(([0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(r))))      # int_1^r
    tot = cum[-1]
    p = -RHO_F * (tot - cum)
    idx = np.rint((rn - 1.0) / (r[1] - r[0])).astype(int)
    return p[idx]


def parse_feb_mesh(path):
    """z = 0 nodes (id, x, y) and the quad connectivity (first four hex8 nodes) of the single fluid domain."""
    nodes, quads = {}, []
    node_re = re.compile(r'<node id="(\d+)">([^<]+)</node>'); el_re = re.compile(r'<elem id="\d+">([^<]+)</elem>')
    with open(path, encoding="latin-1") as f:
        for ln in f:
            m = node_re.search(ln)
            if m:
                x, y, z = map(float, m.group(2).split(","))
                if abs(z) < 1e-14:
                    nodes[int(m.group(1))] = (x, y)
                continue
            m = el_re.search(ln)
            if m:
                quads.append([int(v) for v in m.group(1).split(",")[:4]])
    return nodes, np.array(quads)


class Mesh:
    def __init__(self, feb, ids, nrf):
        nodes, quads = parse_feb_mesh(feb)
        self.pos = {i: k for k, i in enumerate(ids)}
        self.X = np.array([nodes[int(i)] for i in ids])
        self.Q = np.vectorize(self.pos.get)(quads)
        self.r = np.hypot(self.X[:, 0], self.X[:, 1]); self.nrf = nrf
        self.ri = np.rint((self.r - 1.0) * nrf).astype(int)
        self.nn = len(ids); self.nq = len(self.Q)
        G = 1 / math.sqrt(3); self.gp = [(-G, -G), (G, -G), (G, G), (-G, G)]
        # per-Gauss-point shape-function gradients in physical space: B[q, g, a, :] = dN_a/dx, and |J|
        self.B = np.zeros((self.nq, 4, 4, 2)); self.detJ = np.zeros((self.nq, 4)); self.Nsh = np.zeros((4, 4))
        for g, (xi, eta) in enumerate(self.gp):
            dN = 0.25 * np.array([[-(1 - eta), -(1 - xi)], [(1 - eta), -(1 + xi)], [(1 + eta), (1 + xi)], [-(1 + eta), (1 - xi)]])
            self.Nsh[g] = 0.25 * np.array([(1 - xi) * (1 - eta), (1 + xi) * (1 - eta), (1 + xi) * (1 + eta), (1 - xi) * (1 + eta)])
            xe = self.X[self.Q]                                      # (nq, 4, 2)
            J = np.einsum("qak,al->qkl", xe, dN)                     # (nq, 2, 2): dx_k/dxi_l
            self.detJ[:, g] = np.linalg.det(J)
            self.B[:, g] = np.einsum("al,qlk->qak", dN, np.linalg.inv(J))
        # face (xi = -1) gradients at two Gauss points in eta
        self.Bface = np.zeros((self.nq, 2, 4, 2))
        for g, eta in enumerate((-G, G)):
            xi = -1.0
            dN = 0.25 * np.array([[-(1 - eta), -(1 - xi)], [(1 - eta), -(1 + xi)], [(1 + eta), (1 + xi)], [-(1 + eta), (1 - xi)]])
            J = np.einsum("qak,al->qkl", self.X[self.Q], dN)
            self.Bface[:, g] = np.einsum("al,qlk->qak", dN, np.linalg.inv(J))
        self.wall = np.where(self.ri[self.Q[:, 0]] == 0)[0]
        self.outer = np.where(self.ri[self.Q[:, 1]] == nrf)[0]
        # consistent lumped weights: integral of N_a over the element = sum_g detJ N_a(g)
        wnode = np.zeros(self.nn)
        for g in range(4):
            for a in range(4):
                np.add.at(wnode, self.Q[:, a], self.detJ[:, g] * self.Nsh[g, a])
        self.wnode = wnode


def stress2d(Lgrad):
    """Newtonian stress for 2-D velocity-gradient array [..., i, j] = dv_i/dx_j with w_z = 0, no z dependence, kappa = 0 (matches the FEBio input)."""
    D = 0.5 * (Lgrad + np.swapaxes(Lgrad, -1, -2)); tr = D[..., 0, 0] + D[..., 1, 1]
    s = 2 * MU_F * D
    s[..., 0, 0] -= 2 * MU_F / 3 * tr; s[..., 1, 1] -= 2 * MU_F / 3 * tr
    return s


def wall_torques(M, W):
    """T1 (stress from the gradient on the wall face), T2 (Gauss mean over the whole element, the FEBio fluid-FSI traction recipe).  per radian per metre."""
    T1 = T2 = 0.0
    q = M.Q[M.wall]; Wq = W[q]                                       # (nw, 4, 2)
    a, b = M.X[q[:, 0]], M.X[q[:, 3]]
    e = b - a; ell = np.hypot(e[:, 0], e[:, 1]); nf = np.stack([-e[:, 1], e[:, 0]], 1) / ell[:, None]
    mid = 0.5 * (a + b)
    flip = (nf * mid).sum(1) > 0; nf[flip] *= -1                      # normal out of the fluid: towards the origin
    Lg = np.einsum("qgai,qgaj->qgij", np.broadcast_to(Wq[:, None], (len(q), 4, 4, 2)), M.B[M.wall])      # dv_i/dx_j at the 4 Gauss points
    sv = stress2d(Lg).mean(axis=1)                                    # (nw, 2, 2)
    fa = -np.einsum("qij,qj->qi", sv, nf) * (ell * DZ / 4)[:, None]  # force per node (two nodes on the wall edge, two z-layers)
    for node, xy in ((0, a), (3, b)):
        T2 += 2 * np.sum(xy[:, 0] * fa[:, 1] - xy[:, 1] * fa[:, 0])
    Lf = np.einsum("qai,qgaj->qgij", Wq, M.Bface[M.wall])            # (nw, 2, 2, 2)
    sf = stress2d(Lf)
    G = 1 / math.sqrt(3)
    for g, eta in enumerate((-G, G)):
        pos = 0.5 * (1 - eta) * a + 0.5 * (1 + eta) * b
        t = -np.einsum("qij,qj->qi", sf[:, g], nf)
        T1 += np.sum((pos[:, 0] * t[:, 1] - pos[:, 1] * t[:, 0]) * ell / 2 * DZ)
    return T1 / TWO_PI_DZ, T2 / TWO_PI_DZ


def ang_mom(M, W):
    X = M.X[M.Q]; Wq = W[M.Q]
    H = 0.0
    for g in range(4):
        x = np.einsum("a,qak->qk", M.Nsh[g], X); v = np.einsum("a,qak->qk", M.Nsh[g], Wq)
        H += np.sum(RHO_F * (x[:, 0] * v[:, 1] - x[:, 1] * v[:, 0]) * np.abs(M.detJ[:, g]) * DZ)
    return H / TWO_PI_DZ


def outer_flux(M, W):
    q = M.Q[M.outer]; a, b = M.X[q[:, 1]], M.X[q[:, 2]]
    e = b - a; ell = np.hypot(e[:, 0], e[:, 1]); n = np.stack([e[:, 1], -e[:, 0]], 1) / ell[:, None]
    mid = 0.5 * (a + b); flip = (n * mid).sum(1) < 0; n[flip] *= -1
    vm = 0.5 * (W[q[:, 1]] + W[q[:, 2]]); vn = (vm * n).sum(1)
    Q = np.sum(vn * ell * DZ)
    Hf = np.sum(RHO_F * (mid[:, 0] * vm[:, 1] - mid[:, 1] * vm[:, 0]) * vn * ell * DZ) / TWO_PI_DZ
    return Q, Hf


def rtb(raw, inputs, run, nrf, ex_fine, tab_t, tab_om, shift_aligned=0.005):
    d = os.path.join(raw, run)
    T, B = parse_log(os.path.join(d, "Mon0.txt"))
    meta = json.load(open(os.path.join(d, "model_meta.json")))
    ids = B[0][:, 0].astype(int)
    feb = os.path.join(inputs, run + ".feb")
    M = Mesh(feb, ids, nrf)
    K = meta["k"]
    nt = len(T)
    rn = 1 + np.arange(nrf + 1) / nrf
    # exact-solution quantities at the logged times (fine grid has dt = 1e-3, logged every 10 steps)
    dtf = ex_fine["t"][1] - ex_fine["t"][0]
    kk = np.rint(T / dtf).astype(int)
    r_f = ex_fine["r"]; ridx = np.rint((rn - 1.0) / (r_f[1] - r_f[0])).astype(int)
    out = {}
    W_all = np.array([b[:, 1:3] for b in B]); EF = np.array([b[:, 3] for b in B])
    vth_all = (-M.X[:, 1] * W_all[:, :, 0] + M.X[:, 0] * W_all[:, :, 1]) / M.r
    vr_all = (M.X[:, 0] * W_all[:, :, 0] + M.X[:, 1] * W_all[:, :, 1]) / M.r
    wall_nodes = np.where(M.ri == 0)[0]
    wall_vel = vth_all[:, wall_nodes].mean(axis=1)                    # r = 1: v_theta = omega_wall
    T1 = np.zeros(nt); T2 = np.zeros(nt); H = np.zeros(nt); Qo = np.zeros(nt); Hfl = np.zeros(nt)
    for n in range(nt):
        T1[n], T2[n] = wall_torques(M, W_all[n]); H[n] = ang_mom(M, W_all[n]); Qo[n], Hfl[n] = outer_flux(M, W_all[n])
    dH = np.gradient(H, T, edge_order=2)
    T3 = -(dH + Hfl)
    res = dict(run=run, n_times=nt, dt=float(T[1] - T[0]))

    def metrics(shift, label):
        # exact solution of the problem whose wall history is delayed by `shift`: evaluate the shifted fine solution
        ex = ex_fine
        ksh = np.clip(np.rint((T - shift) / dtf).astype(int), 0, None)
        Lx = ex["L"][ksh][:, ridx]; Vex = Lx / rn                      # (nt, nrf+1)
        Tex = ex["T"][ksh]; Hex = ex["H"][ksh]
        e2 = (vth_all - Vex[:, M.ri]) ** 2 + vr_all ** 2
        num = np.array([np.sum(M.wnode * e2[n]) for n in range(nt)]); den = np.array([np.sum(M.wnode * Vex[n, M.ri] ** 2) for n in range(nt)])
        trap = lambda f: float(np.sum(0.5 * (f[1:] + f[:-1]) * np.diff(T)))
        Em = lambda q, ref: float(np.max(np.abs(q - ref)) / np.max(np.abs(ref)))
        El2 = lambda q, ref: math.sqrt(trap((q - ref) ** 2) / trap(ref ** 2))
        o = dict(shift_s=shift, E_L2_velocity_areaweighted=math.sqrt(trap(num) / trap(den)),
                 E_max_T1=Em(T1, Tex), E_max_T2=Em(T2, Tex), E_max_T3=Em(T3, Tex), E_L2_T2=El2(T2, Tex), E_L2_T3=El2(T3, Tex), E_max_H=Em(H, Hex),
                 E_max_T2_sign_flipped=Em(-T2, Tex), max_abs_T_exact=float(np.max(np.abs(Tex))), max_abs_H_exact=float(np.max(np.abs(Hex))),
                 T2_over_Texact_at_peak=float(T2[np.argmax(np.abs(Tex))] / Tex[np.argmax(np.abs(Tex))]))
        return o, Tex
    res["as_run"], Tex0 = metrics(0.0, "as_run")
    res["aligned_dt2"], Tex1 = metrics(shift_aligned, "aligned")
    # BC time convention: logged wall tangential velocity against the table at t and at t - dt/2
    om_t = np.interp(T, tab_t, tab_om); om_h = np.interp(np.clip(T - 0.005, 0, None), tab_t, tab_om)
    sel = (T >= 0.6)
    res["wall_velocity_vs_table"] = dict(max_abs_diff_to_omega_t=float(np.max(np.abs(wall_vel - om_t))), max_abs_diff_to_omega_t_minus_dt2=float(np.max(np.abs(wall_vel - om_h))),
                                         max_rel_diff_to_omega_t_minus_dt2_t_ge_0p6=float(np.max(np.abs(wall_vel[sel] - om_h[sel]) / np.maximum(np.abs(om_h[sel]), 0.02))),
                                         example_t1p0=dict(logged=float(wall_vel[np.argmin(abs(T - 1.0))]), table_t=float(om_t[np.argmin(abs(T - 1.0))]), table_t_minus_dt2=float(om_h[np.argmin(abs(T - 1.0))])))
    # pressure parity at the time of max wall centripetal pressure of the exact solution (as-run time labels)
    pw = np.array([pressure_exact_at(ex_fine["L"][k], r_f, rn)[0] for k in kk])
    ip = int(np.argmax(np.abs(pw)))
    P = -K * EF[ip]
    grid = np.zeros((nrf + 1, M.nn // (nrf + 1)))
    rj = np.rint((np.arctan2(M.X[:, 1], M.X[:, 0]) % (2 * math.pi)) / (2 * math.pi / grid.shape[1])).astype(int) % grid.shape[1]
    grid[M.ri, rj] = P
    pbar = grid.mean(axis=1)
    pex = pressure_exact_at(ex_fine["L"][kk[ip]], r_f, rn)
    e = pbar - pex; hf = e[1:-1] - 0.5 * (e[:-2] + e[2:])
    res["pressure_at_t_pmax"] = dict(t=float(T[ip]), p_exact_max_abs_Pa=float(np.max(np.abs(pex))), parity_theta_Pa=float(np.abs(grid - pbar[:, None]).max()),
                                     parity_radial_Pa=float(np.abs(hf).max()), parity_theta_rel=float(np.abs(grid - pbar[:, None]).max() / np.max(np.abs(pex))),
                                     parity_radial_rel=float(np.abs(hf).max() / np.max(np.abs(pex))), smooth_error_max_Pa=float(np.abs(e).max()))
    vws = float(np.max(np.abs(wall_vel)))
    res["continuity_Qout_rel_to_wall_flux_scale"] = float(np.abs(Qo).max() / (vws * 2 * math.pi * 1.0 * DZ))
    res["series"] = dict(t=T.tolist(), T_exact=Tex0.tolist(), T_exact_shifted=Tex1.tolist(), T1=T1.tolist(), T2=T2.tolist(), T3=T3.tolist(), H=H.tolist(), wall_velocity=wall_vel.tolist())
    return res


def main():
    raw, inputs, outp = sys.argv[1], sys.argv[2], sys.argv[3]
    tab_t, tab_om = omega_table(os.path.join(ROOT, "reference_data", "rtr", "omega_B_table_dt1e-3.csv"))
    out = {"note": "independent re-derivation (own parsers, own algebra, own finite-difference fluid reference); see docstring"}
    out["RT-A_law_identity"] = rta_law(raw)
    out["RT-A_rotation_expression"] = rta_rot(raw)
    out["RT-A_twist"] = rta_twist(raw, os.path.join(ROOT, "results", "static_twist_reference_TA.json"))
    ex = fd_exact(tab_t, tab_om, N=6400, dt=5e-4)
    ex_check = fd_exact(tab_t, tab_om, N=3200, dt=1e-3)
    out["fd_reference_self_check"] = dict(N=6400, dt=5e-4, N_coarse=3200, dt_coarse=1e-3,
                                          max_rel_diff_T_between_resolutions=float(np.max(np.abs(ex["T"][::2] - ex_check["T"])) / np.max(np.abs(ex["T"]))),
                                          max_abs_T=float(np.max(np.abs(ex["T"]))), max_abs_H=float(np.max(np.abs(ex["H"]))))
    # compare with the spectral reference used in the campaign (cached arrays, if present)
    cache = os.path.join(ROOT, "reference_data", "rtr", "fluid_exact_N40_dt0.01.npz")
    if os.path.exists(cache):
        c = np.load(cache, allow_pickle=True)
        k10 = np.arange(0, 201) * 20
        out["fd_reference_self_check"]["max_rel_diff_T_vs_spectral_reference"] = float(np.max(np.abs(ex["T"][k10] - c["T_B"])) / np.max(np.abs(c["T_B"])))
        out["fd_reference_self_check"]["max_rel_diff_H_vs_spectral_reference"] = float(np.max(np.abs(ex["H"][k10] - c["H_f"])) / np.max(np.abs(c["H_f"])))
    for run, nrf in (("RTB_fluid_M1", 16), ("RTB_fluid_M2", 32), ("RTB_fluid_M1_bcshift", 16)):
        out[run] = rtb(raw, inputs, run, nrf, ex, tab_t, tab_om)
    for nm in ("as_run", "aligned_dt2"):
        a, b = out["RTB_fluid_M1"][nm], out["RTB_fluid_M2"][nm]
        out[f"orders_{nm}"] = {k: math.log2(a[k] / b[k]) for k in ("E_L2_velocity_areaweighted", "E_max_T1", "E_max_T2", "E_max_T3", "E_max_H")}
        out[f"orders_{nm}"]["velocity_ratio_M1_over_M2"] = a["E_L2_velocity_areaweighted"] / b["E_L2_velocity_areaweighted"]
    json.dump(out, open(outp, "w"), indent=1, default=float)
    brief = {k: v for k, v in out.items() if not k.startswith("RTB_")}
    print(json.dumps(brief, indent=1, default=float))
    for run in ("RTB_fluid_M1", "RTB_fluid_M2", "RTB_fluid_M1_bcshift"):
        print(run, json.dumps({k: v for k, v in out[run].items() if k != "series"}, indent=1, default=float))


if __name__ == "__main__":
    main()
