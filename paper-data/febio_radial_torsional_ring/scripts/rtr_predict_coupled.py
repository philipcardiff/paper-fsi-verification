#!/usr/bin/env python3
"""AE ONLY (not a FEBio result): predicted effect of FEBio's element-mean interface shear traction on the COUPLED torsion QoIs, obtained from the proposal's own reduced model
(vendored rt_reference, unmodified) with ONE change in a subclass: the shear stress transmitted to the solid is evaluated at r = b + h/2 (element-centre value, h = radial fluid cell
size) instead of at r = b.  Everything else (law, fluid, convection, solid) is the reference.  usage: rtr_predict_coupled.py [h_inverse ...] (default 16 32)."""
import sys, os, json, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import rtr_reference as R
rt = R.rt


def bary_eval(x, w, f, xq):
    d = xq - x
    if np.any(d == 0): return f[np.argmin(np.abs(d))]
    t = w / d
    return (t @ f) / t.sum()


class EMModel(rt.Model):
    h = 0.0
    def interface(self, t, b, vB, L, ell, vBdot=None):
        F_r0, m_add, T_B, ex = super().interface(t, b, vB, L, ell)
        r = ex["r"]; Lr = ex["Lr"]; mu = self.p.mu_f
        sig = mu * (Lr / r - 2 * L / r ** 2)
        k = np.arange(len(r)); wts = (-1.0) ** k; wts[0] *= 0.5; wts[-1] *= 0.5
        rc = b + self.h
        sig_c = bary_eval(r, wts, sig, rc)
        ex = dict(ex); ex["sig_rt"] = sig_c
        return F_r0, m_add, b * b * sig_c, ex


def run(hinv, tt):
    p = rt.Params(theta0=R.THETA0, Omega=R.OMEGA, P0=0.0, T_ramp=R.TRAMP, ramp="Cinf")
    m = EMModel(p, 20, 32); m.h = 0.5 / hinv
    y0 = np.zeros(m.n); th = p.theta0; ua = max(th * 1e-2 * p.B, 1e-30); Ls = max(th * p.Omega * p.B ** 2, ua * 10)
    atol = np.concatenate((np.full(m.Ns, ua), np.full(m.Ns, max(th, ua)), np.full(m.Ns, ua * 10), np.full(m.Ns, max(th, ua) * 10), np.full(m.Nf - 1, Ls))) * 1e-12
    from scipy.integrate import solve_ivp
    sol = solve_ivp(m.rhs, (0, tt[-1]), y0, method="Radau", rtol=1e-9, atol=atol * 1e3, t_eval=tt)
    assert sol.success, sol.message
    H = rt.history(m, sol)
    return {k: np.asarray(v) for k, v in H.items()}


if __name__ == "__main__":
    hs = [int(x) for x in sys.argv[1:]] or [16, 32]
    tt = np.round(np.arange(0, 201) * 0.01, 10)
    ref = R.coupled_reference(tt)
    out = {}
    for hinv in hs:
        H = run(hinv, tt)
        def E(k, qk=None): return float(np.max(np.abs(H[k] - ref[k])) / np.max(np.abs(ref[k])))
        out[f"h=1/{hinv}"] = dict(E_max_alpha_B=E("alpha_B"), E_max_T_B_applied=E("T_B"), E_max_H_f=E("H_f"), E_max_sig_rt=E("sig_rt_b"))
        print(f"h = 1/{hinv}", json.dumps(out[f"h=1/{hinv}"]))
    json.dump(out, open(os.environ.get("RTR_PREDICT_OUT", os.path.join(R.ROOT, "results", "predicted_coupled_element_mean_AE.json")), "w"), indent=1)
