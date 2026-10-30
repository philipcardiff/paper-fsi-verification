#!/usr/bin/env python3
import os, os as _os; HERE = _os.path.dirname(_os.path.abspath(__file__)); PKG = _os.path.dirname(HERE)   # package root; run directories of retained HPC outputs: env CAVITY_RUNS (default runs/cavityNext)
"""B4 transient emulation: linear Stokes + inertia (no convection) with the FEBio fluid-FSI element blocks (b4_operator.assemble), the cavity M1 mesh and the model's BCs
(no-slip walls/floor, built-in parabolic inlet with the inlet flux load of FEFluidNormalVelocity::LoadVector, outlet p=0, rim pressure = lagged mean of the free inlet nodes, optional corner pins 1.15 Pa),
backward-Euler in time (dt, ramp 25 s). Purpose: does the operator + BCs by themselves excite a wall-normal row-parity pressure mode of the observed size (D(M1)=0.36 Pa)?
usage: b4_transient_stokes.py OUT.json [K dt tend pin(0/1) rim(1/0) f]"""
import sys, json, numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla
sys.path.insert(0, HERE); import b4_operator as B
def go(K=1e6, dt=2.5, tend=1000.0, pin=True, rim=True, f=1, ramp=25.0, kind="cavity", conv_steps=None, perturb_at=None, amp=0.1, perturb_kind="y"):
    g = B.geometry(kind, f); X = g["X"]; n = len(X); c = 1.0 / dt
    A0, G, M = B.assemble(X, g["quads"], c=0.0); At, _, _ = B.assemble(X, g["quads"], c=c); Mv = (At - A0) / (B.RHO * c)
    wfree, pfree, wfix, pfix = B.bcs(g, 1.15 if pin else None, rim_free=not rim); wfix_l = np.array(sorted(set(range(2 * n)) - set(wfree.tolist())))
    prof = B.inlet_profile(g); bin_, _ = B.inlet_load(g)
    wp = np.zeros(2 * n)
    for i, ph in prof.items(): wp[2 * i] = ph           # prescribed inlet nodal velocity (ramped)
    GT = G.T.tocsr(); Gw = G[wfree]; Gp = G[wfree][:, pfree]; GTw = GT[pfree][:, wfree]; Mpp = M[pfree][:, pfree]; Mpall = M
    pf_list = list(pfree); pidx = {nd: k for k, nd in enumerate(pf_list)}
    Kmat = sp.bmat([[At[wfree][:, wfree], Gp], [GTw, -(1.0 / (K * dt)) * Mpp]]).tocsc(); lu = sla.splu(Kmat)
    Awp = At[wfree][:, wfix_l]; GTwp = GT[pfree][:, wfix_l]; Mvff = Mv[wfree][:, wfree]
    # pressure-fixed values
    pin_nodes = [i for i in g["iface"] if abs(X[i, 0] - 2) < 1e-6 or abs(X[i, 0] - 6) < 1e-6] if pin else []
    pfix_l = sorted(pfix); Gpfix = G[wfree][:, pfix_l]; Mpfix = M[pfree][:, pfix_l]
    inl_free = [i for i in g["inlet"] if i not in g["in_corner"] and i in pidx]
    chains = B.rows_above_floor(g); ch = [chains[a] for a in g["iface"]]
    w = np.zeros(2 * n); p = np.zeros(n); nstep = int(round(tend / dt)); hist = []
    ifs = g["iface"]; xi = X[ifs, 0]
    def forces(pv):
        # F1: wall row p (bilinear on facets, flat floor): -sum trapezoid;  F2: pair mean (row0,row1)
        p0 = pv[ifs]; p1 = np.array([pv[c_[1]] for c_ in ch]); dx = np.diff(xi)
        return float(-np.sum(0.5 * (p0[1:] + p0[:-1]) * dx)), float(-np.sum(0.5 * ((p0[1:] + p0[:-1] + p1[1:] + p1[:-1]) / 2) * dx))
    for s in range(1, nstep + 1):
        t = s * dt
        if perturb_at is not None and abs(t - perturb_at) < 1e-9:
            for a_, ch_ in chains.items():
                if 2.5 < X[a_, 0] < 5.5:
                    for j_, m_ in enumerate(ch_):
                        if m_ in pidx: p[m_] += amp * ((-1.0) ** j_ if perturb_kind == 'y' else 1.0)
            base = float(np.mean([p[c_[1]] - p[c_[0]] for c_ in ch[len(ch)//2-5:len(ch)//2+5]]))
        lc = min(1.0, t / ramp); wpt = lc * wp
        pfix_val = np.zeros(len(pfix_l))
        for k, nd in enumerate(pfix_l):
            if nd in pin_nodes: pfix_val[k] = 1.15 * lc
            elif nd in g["in_corner"] and rim: pfix_val[k] = np.mean(p[inl_free]) if inl_free else 0.0
        # rhs
        rw = -Awp @ wpt[wfix_l] - Gpfix @ pfix_val + (1.0 / dt) * B.RHO * (Mvff @ w[wfree])
        # continuity (rows of free p): GT w_all - (1/(K dt)) M (p - p_old)  = bload  (M includes prescribed-p columns)
        rp = -bin_[pfree] * lc - GTwp @ wpt[wfix_l] + (1.0 / (K * dt)) * (Mpp @ p[pfree]) - (1.0 / (K * dt)) * (Mpfix @ (pfix_val - p[pfix_l]) * 0.0) - (-(1.0 / (K * dt)) * (Mpfix @ pfix_val)) * 0.0
        # exact treatment of the prescribed-p terms in the mass rows: -(1/(K dt)) M_{free,fix} (p_fix_new - p_fix_old)
        rp += (1.0 / (K * dt)) * (Mpfix @ (pfix_val - p[pfix_l])) * 1.0 * (-1.0) * (-1.0)
        x = lu.solve(np.concatenate([rw, rp])); w[wfree] = x[:len(wfree)]; w[wfix_l] = wpt[wfix_l]; p[pfree] = x[len(wfree):]; p[pfix_l] = pfix_val
        if s % max(1, int(round((25 if perturb_at is None else 12.5) / dt))) == 0 or s == nstep:
            p0 = p[ifs]; p1 = np.array([p[c_[1]] for c_ in ch]); p2 = np.array([p[c_[2]] if len(c_) > 2 else np.nan for c_ in ch])
            Fy1, Fy2 = forces(p); inn = slice(4, -4)
            hist.append(dict(t=t, Dc=float(np.mean([p[c_[1]] - p[c_[0]] for c_ in ch[len(ch)//2-5:len(ch)//2+5]])),  D_row1_minus_row0=float(np.mean((p1 - p0)[inn])), D_row2_minus_row1=float(np.mean((p2 - p1)[inn])), p_row0_mean=float(np.mean(p0[inn])), F1_wallrow=Fy1, F2_pairmean=Fy2, Qin=float(np.sum(bin_ * 0) ), umax=float(np.abs(w).max())))
    return dict(K=K, dt=dt, pin=pin, rim=rim, hist=hist)
if __name__ == "__main__":
    out = sys.argv[1]; res = {}
    for name, kw in (("base_K1e6_dt2.5_pin", {}), ("K1e7", dict(K=1e7)), ("dt1.25", dict(dt=1.25)), ("nopin", dict(pin=False)), ("norim", dict(rim=False)),
                     ("PERTURB_parity_y_amp0.1_t300", dict(tend=1300.0, perturb_at=300.0)), ("PERTURB_parity_y_K1e7", dict(tend=1300.0, perturb_at=300.0, K=1e7)), ("PERTURB_parity_y_dt1.25", dict(tend=1300.0, perturb_at=300.0, dt=1.25)), ("PERTURB_smooth_const_rows_amp0.1", dict(tend=1300.0, perturb_at=300.0, perturb_kind="s"))):
        r = go(**kw); res[name] = r; h = r["hist"]; print(name, [(round(x["t"]), round(x["Dc"], 5), round(x["p_row0_mean"], 4), round(x["F1_wallrow"], 4), round(x["F2_pairmean"], 4)) for x in (h[0], h[len(h) // 4], h[len(h) // 2], h[-1])]); sys.stdout.flush()
    json.dump(res, open(out, "w"), indent=1)
