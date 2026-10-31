#!/usr/bin/env python3
"""RT-B (fluid, fixed geometry) analysis for the radialTorsionalRing pilot.  usage: rtr_analyze.py fluid RUN_DIR [OUT.json]
Definitions fixed in RTR_PILOT_PLAN.md / RTR_EQUIVALENCE_AUDIT.md section 3 (T1 face-evaluated, T2 FEBio-applied Gauss-mean, T3 angular-momentum balance).
Torque = torque of the fluid ON the solid about +z, per radian and per metre of depth (extensive / (2 pi dz)).  FEBio text logs hold 12 significant digits."""
import sys, os, json, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.interpolate import BarycentricInterpolator as BI
import rtr_reference as R
from rtr_check_rta import read_log

MU, RHO, DZ = 300.0, 1000.0, 0.1
G = 1 / math.sqrt(3)
GP = [(-G, -G), (G, -G), (G, G), (-G, G)]


def dN(xi, eta):
    return 0.25 * np.array([[-(1 - eta), -(1 - xi)], [(1 - eta), -(1 + xi)], [(1 + eta), (1 + xi)], [-(1 + eta), (1 - xi)]])


def stress_from_L(L):
    D = 0.5 * (L + L.T)
    return 2 * MU * D + (0.0 - 2 * MU / 3) * np.trace(D) * np.eye(2)       # kappa = 0 (Stokes), as the FEBio Newtonian fluid input


class Fluid2D:
    """z = 0 layer of the slab as a 2-D quad mesh; nodes of the log in file order."""
    def __init__(self, run, ids):
        m = np.load(os.path.join(run, "mesh_rtr.npz"))
        self.nth, self.nrf = int(m["nth"]), int(m["nrf"])
        nodes = m["nodes"]; self.fel = m["fluid_el"]
        z0 = {i: k for k, i in enumerate(ids)}                       # file row -> node id (1-based ids in log)
        self.id2row = {int(i) - 1: k for k, i in enumerate(ids)}
        self.xy = np.array([nodes[int(i) - 1][:2] for i in ids])
        self.Q = np.array([[self.id2row[int(n)] for n in e[:4]] for e in self.fel])
        r = np.hypot(self.xy[:, 0], self.xy[:, 1]); th = np.arctan2(self.xy[:, 1], self.xy[:, 0]) % (2 * math.pi)
        self.ri = np.rint((r - 1.0) * self.nrf).astype(int); self.rj = np.rint(th / (2 * math.pi / self.nth)).astype(int) % self.nth
        self.r = r
        assert np.abs(r - (1 + self.ri / self.nrf)).max() < 1e-9
        # elements adjacent to the wall: first radial layer (ri of node 0 == 0)
        self.wall = [k for k, q in enumerate(self.Q) if self.ri[q[0]] == 0]
        self.outer = [k for k, q in enumerate(self.Q) if self.ri[q[1]] == self.nrf]

    def grad(self, W, X, qi, xi, eta):
        d = dN(xi, eta); J = X[self.Q[qi]].T @ d
        return W[self.Q[qi]].T @ (d @ np.linalg.inv(J)), np.linalg.det(J)


def torques(F, W, X=None):
    """T1, T2 (viscous part), p-part torque, H_f at one time for velocity W (n,2), nodal pressure P (n,)."""
    X = F.xy if X is None else X
    T1 = T2 = 0.0
    for qi in F.wall:
        q = F.Q[qi]; a, b = q[0], q[3]                 # wall edge: local nodes 1 (xi=-1, eta=-1) and 4 (xi=-1, eta=+1)
        e = X[b] - X[a]; ell = math.hypot(*e)
        nf = np.array([-e[1], e[0]]) / ell              # candidate normal; fluid outward normal at the inner wall points to smaller r
        if np.dot(nf, 0.5 * (X[a] + X[b])) > 0: nf = -nf
        # T2: Gauss-mean stress of the element, nodal forces f_a = -(sv_mean n_f) ell dz / 4 per layer-node, two layers
        sv = np.zeros((2, 2))
        for xi, eta in GP:
            L, _ = F.grad(W, X, qi, xi, eta); sv += stress_from_L(L)
        sv /= 4
        fa = -(sv @ nf) * ell * DZ / 4.0
        for node in (a, b):
            T2 += 2 * (X[node, 0] * fa[1] - X[node, 1] * fa[0])
        # T1: gradient evaluated on the wall face (xi = -1), 2 Gauss points in eta, true position on the chord
        for eta in (-G, G):
            L, _ = F.grad(W, X, qi, -1.0, eta); s = stress_from_L(L)
            pos = 0.5 * (1 - eta) * X[a] + 0.5 * (1 + eta) * X[b]
            t = -(s @ nf)
            T1 += (pos[0] * t[1] - pos[1] * t[0]) * (ell / 2.0) * DZ
    return T1, T2


def angular_momentum(F, W, X=None):
    X = F.xy if X is None else X
    H = 0.0
    for qi in range(len(F.Q)):
        for xi, eta in GP:
            d = dN(xi, eta); J = X[F.Q[qi]].T @ d; dj = abs(np.linalg.det(J))
            N = 0.25 * np.array([(1 - xi) * (1 - eta), (1 + xi) * (1 - eta), (1 + xi) * (1 + eta), (1 - xi) * (1 + eta)])
            x = N @ X[F.Q[qi]]; v = N @ W[F.Q[qi]]
            H += RHO * (x[0] * v[1] - x[1] * v[0]) * dj * DZ
    return H


def outer_flux(F, W):
    Q = 0.0; Hflux = 0.0
    for qi in F.outer:
        q = F.Q[qi]; a, b = q[1], q[2]
        e = F.xy[b] - F.xy[a]; ell = math.hypot(*e); n = np.array([e[1], -e[0]]) / ell
        if np.dot(n, 0.5 * (F.xy[a] + F.xy[b])) < 0: n = -n
        vm = 0.5 * (W[a] + W[b]); Q += (vm @ n) * ell * DZ
        xm = 0.5 * (F.xy[a] + F.xy[b]); Hflux += RHO * (xm[0] * vm[1] - xm[1] * vm[0]) * (vm @ n) * ell * DZ
    return Q, Hflux


def exact_fields(times, shift=0.0):
    """Independent exact solution at the logged times (cached): L(r), T_B, H_f, p(r).  shift > 0 evaluates the exact solution at t - shift (the exact solution of the problem
    whose wall-velocity history is delayed by `shift`; the PDE is time invariant and starts from rest, so this is exact)."""
    cache = os.path.join(R.ROOT, "reference_data", "rtr", "fluid_exact_N40_dt0.01%s.npz" % ("" if shift == 0 else "_shift%g" % shift))
    tt = np.round(np.arange(0, 201) * 0.01, 10)
    tq = np.clip(tt - shift, 0.0, None)
    if os.path.exists(cache):
        d = np.load(cache, allow_pickle=True); ex = {k: d[k] for k in d.files}
    else:
        e = R.fluid_exact(tq, N=40)
        P = np.array([R.pressure_profile(e, i) for i in range(len(tt))])
        ex = dict(t=tt, r=e["r"], L=e["L"], T_B=e["T_B"], H_f=e["H_f"], p=P)
        np.savez(cache, **ex)
    assert np.allclose(ex["t"][:len(times)], times, atol=1e-9)
    return ex


def analyze_fluid(run, shift=0.0):
    run = os.path.realpath(run); K = 1e7
    meta = json.load(open(run + "/model_meta.json")); K = meta.get("k", K); k2 = 2 * math.pi * DZ
    T, D = read_log(run + "/Mon0.txt")
    ids = D[0][:, 0].astype(int)
    F = Fluid2D(run, ids)
    nt = len(T); ex = exact_fields(T, shift)
    rr = ex["r"][::-1]; Lex = ex["L"][:, ::-1]
    # exact v_theta at the nodal radii (spectral interpolation of L, v = L/r)
    rn = 1 + np.arange(F.nrf + 1) / F.nrf
    Vex = np.array([BI(rr, Lex[i])(rn) / rn for i in range(nt)])      # (nt, nrf+1)
    Pex = np.array([BI(rr, ex["p"][i][::-1])(rn) for i in range(nt)])
    out = dict(t=T.tolist(), T1=[], T2=[], H=[], Qout=[], Fluxtorque=[], par_theta=[], par_rad=[], p_smooth_err=[], Verr2=[], Vex2=[], Verr2_eq=[], Vex2_eq=[], vr_max=[], psum=[], P_part=[])
    # radial weights for area-weighted L2: trapezoid in r of r dr (per ring)
    w_r = np.zeros(F.nrf + 1); dr = 1.0 / F.nrf
    w_r[:] = rn * dr; w_r[0] *= 0.5; w_r[-1] *= 0.5
    for n in range(nt):
        d = D[n]; W = d[:, 1:3]; P = -K * d[:, 3]
        r = F.r; vth = (-F.xy[:, 1] * W[:, 0] + F.xy[:, 0] * W[:, 1]) / r; vr = (F.xy[:, 0] * W[:, 0] + F.xy[:, 1] * W[:, 1]) / r
        # error fields on the (nrf+1, nth) grid
        Vth = np.zeros((F.nrf + 1, F.nth)); Vr = np.zeros_like(Vth); Pg = np.zeros_like(Vth)
        Vth[F.ri, F.rj] = vth; Vr[F.ri, F.rj] = vr; Pg[F.ri, F.rj] = P
        e2 = ((Vth - Vex[n][:, None]) ** 2 + Vr ** 2); out["Verr2"].append(float((w_r[:, None] * e2).sum() * 2 * math.pi / F.nth))
        out["Vex2"].append(float((w_r[:, None] * np.repeat(Vex[n][:, None] ** 2, F.nth, 1)).sum() * 2 * math.pi / F.nth))
        out["Verr2_eq"].append(float(e2.mean())); out["Vex2_eq"].append(float((Vex[n] ** 2).mean()))
        out["vr_max"].append(float(np.abs(Vr).max()))
        pbar = Pg.mean(axis=1); out["par_theta"].append(float(np.abs(Pg - pbar[:, None]).max()))
        e = pbar - Pex[n]; hf = np.zeros_like(e); hf[1:-1] = e[1:-1] - 0.5 * (e[:-2] + e[2:]); out["par_rad"].append(float(np.abs(hf).max()))
        out["p_smooth_err"].append(float(np.abs(e).max()))
        t1, t2 = torques(F, W)
        out["T1"].append(t1 / k2); out["T2"].append(t2 / k2)
        out["H"].append(angular_momentum(F, W) / k2)
        q, hfl = outer_flux(F, W); out["Qout"].append(q); out["Fluxtorque"].append(hfl / k2)
    o = {k: np.array(v) for k, v in out.items()}
    tt = o["t"]; Tex = ex["T_B"][:nt]; Hex = ex["H_f"][:nt]
    # T3: torque of the fluid on the solid = -dH/dt + outer flux torque (H per radian per m; central differences)
    dH = np.gradient(o["H"], tt, edge_order=2)
    T3 = -dH - o["Fluxtorque"] * 0   # outer-flux angular momentum leaves the fluid: dH/dt = T_wall_on_fluid - flux ; T_on_solid = -T_wall_on_fluid = -(dH/dt + flux)
    T3 = -(dH + o["Fluxtorque"])
    sc = np.max(np.abs(Tex))
    def Emax(q, ref): return float(np.max(np.abs(q - ref)) / np.max(np.abs(ref)))
    def trap(f): return float(np.sum(0.5 * (f[1:] + f[:-1]) * np.diff(tt)))
    res = dict(stage="RT-B", run=os.path.basename(run), time_shift_of_exact_s=shift, n_times=nt, dt=float(tt[1] - tt[0]), K=K, mesh=(F.nth, F.nrf))
    res["E_L2_velocity_areaweighted"] = math.sqrt(trap(o["Verr2"]) / trap(o["Vex2"]))
    res["E_L2_velocity_nodal_equal"] = math.sqrt(trap(o["Verr2_eq"]) / trap(o["Vex2_eq"]))
    for nm, q in (("T1", o["T1"]), ("T2", o["T2"]), ("T3", T3)):
        res["E_max_" + nm] = Emax(q, Tex); res["E_L2_" + nm] = math.sqrt(trap((q - Tex) ** 2) / trap(Tex ** 2))
    res["E_max_T2_sign_flipped"] = Emax(-o["T2"], Tex)
    res["E_max_H"] = Emax(o["H"], Hex); res["E_L2_H"] = math.sqrt(trap((o["H"] - Hex) ** 2) / trap(Hex ** 2))
    i_p = int(np.argmax(np.abs(ex["p"][:nt, 0]))); pmax = float(np.abs(Pex[i_p]).max())
    res["t_pmax"] = float(tt[i_p]); res["p_exact_max_abs_at_t"] = pmax
    res["parity_theta_Pa_at_tpmax"] = float(o["par_theta"][i_p]); res["parity_radial_Pa_at_tpmax"] = float(o["par_rad"][i_p])
    res["parity_theta_rel"] = res["parity_theta_Pa_at_tpmax"] / pmax; res["parity_radial_rel"] = res["parity_radial_Pa_at_tpmax"] / pmax
    res["parity_theta_max_over_time_Pa"] = float(o["par_theta"].max()); res["parity_radial_max_over_time_Pa"] = float(o["par_rad"].max())
    res["pressure_smooth_err_max_over_time_Pa"] = float(o["p_smooth_err"].max()); res["pressure_smooth_err_rel_at_tpmax"] = float(o["p_smooth_err"][i_p] / pmax)
    vws = float(np.max(np.abs(Vex[:, 0]))); res["continuity_Qout_rel"] = float(np.abs(o["Qout"]).max() / (vws * 2 * math.pi * 1.0 * DZ))
    res["vr_max_over_max_vtheta_wall"] = float(o["vr_max"].max() / vws)
    res["max_abs_T_exact"] = float(sc); res["max_abs_H_exact"] = float(np.max(np.abs(Hex)))
    res["T2_over_Texact_at_peak"] = float(o["T2"][np.argmax(np.abs(Tex))] / Tex[np.argmax(np.abs(Tex))])
    res["series"] = dict(t=tt.tolist(), T_exact=Tex.tolist(), T1=o["T1"].tolist(), T2=o["T2"].tolist(), T3=T3.tolist(), H=o["H"].tolist(), H_exact=Hex.tolist())
    return res


if __name__ == "__main__":
    kind, run = sys.argv[1], sys.argv[2]
    shift = float(sys.argv[4]) if len(sys.argv) > 4 else 0.0
    res = {"fluid": analyze_fluid}[kind](run, shift)
    if len(sys.argv) > 3 and sys.argv[3] != "-": json.dump(res, open(sys.argv[3], "w"), indent=1, default=float)
    print(json.dumps({k: v for k, v in res.items() if k != "series"}, indent=1, default=float))
