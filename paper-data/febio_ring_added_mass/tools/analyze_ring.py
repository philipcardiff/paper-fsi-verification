#!/usr/bin/env python3
"""Post-process a FEBio ring run: frequency (two independent diagnostics), amplitude decay.

QoI as in the solids4foam driver: D(t) = (u_r(theta=0) - u_r(theta=90deg))/2 on the outer surface.
probes.txt rows: node_id x y ux uy vx vy (x,y are current coordinates). Probe 1 is the outer surface
node on the y=0 plane (u_r = ux), probe 2 the outer node on x=0 (u_r = uy).
Only numpy is used. Nothing here is tuned to the reference value.
"""
from __future__ import annotations
import json, math, sys
import numpy as np


def read_blocks(path):
    out, cur = [], None
    for line in open(path):
        if line.startswith("*Step"):
            cur = {"rows": []}; out.append(cur)
        elif line.startswith("*Time"):
            cur["t"] = float(line.split("=")[1])
        elif line.startswith("*"):
            pass
        elif line.strip():
            cur["rows"].append([float(v) for v in line.split()])
    return out


def series(probe_file):
    b = read_blocks(probe_file)
    t = np.array([x["t"] for x in b])
    r0 = np.array([x["rows"][0] for x in b]); r1 = np.array([x["rows"][1] for x in b])
    d = 0.5 * (r0[:, 3] - r1[:, 4])                   # (ux at theta=0 - uy at theta=90)/2
    v = 0.5 * (r0[:, 5] - r1[:, 6])
    return t, d, v


def zero_crossing_omega(t, d, t0, t1):
    m = (t >= t0) & (t <= t1)
    tt, dd = t[m], d[m]
    cr = []
    for i in range(len(tt) - 1):
        if dd[i] == 0.0:
            cr.append(tt[i])
        elif dd[i] * dd[i + 1] < 0:
            cr.append(tt[i] - dd[i] * (tt[i + 1] - tt[i]) / (dd[i + 1] - dd[i]))
    cr = np.array(cr)
    if len(cr) < 3:
        return float("nan"), len(cr)
    # least-squares line through crossing index vs time: half period = slope
    k = np.arange(len(cr))
    slope = np.polyfit(k, cr, 1)[0]
    return math.pi / slope, len(cr)


def fit_damped(t, d, t0, t1, w_guess):
    """d = exp(-s t)(a sin wt + b cos wt) ; variable projection over (w, s)."""
    m = (t >= t0) & (t <= t1)
    tt, dd = t[m], d[m]

    def resid(w, s):
        e = np.exp(-s * tt)
        M = np.column_stack([e * np.sin(w * tt), e * np.cos(w * tt)])
        c, *_ = np.linalg.lstsq(M, dd, rcond=None)
        r = dd - M @ c
        return float(r @ r), c

    best = (1e300, w_guess, 0.0)
    w_lo, w_hi, s_lo, s_hi = 0.9 * w_guess, 1.1 * w_guess, -0.02, 0.02
    for _ in range(40):                                # nested grid refinement
        ws = np.linspace(w_lo, w_hi, 21); ss = np.linspace(s_lo, s_hi, 21)
        for w in ws:
            for s in ss:
                r, _ = resid(w, s)
                if r < best[0]:
                    best = (r, w, s)
        dw, ds = (w_hi - w_lo) / 10, (s_hi - s_lo) / 10
        w_lo, w_hi = best[1] - dw, best[1] + dw
        s_lo, s_hi = best[2] - ds, best[2] + ds
    r, c = resid(best[1], best[2])
    amp = math.hypot(*c)
    return dict(omega=best[1], decay_rate=best[2], amplitude0=amp, rms_resid=math.sqrt(r / len(tt)))


def analyse(probe_file, omega_ref, periods=4.0):
    t, d, v = series(probe_file)
    T = 2 * math.pi / omega_ref
    res = dict(n=len(t), t_end=float(t[-1]), omega_ref=omega_ref, windows={})
    for name, (a, b) in {"0-4T": (0, 4 * T), "T-4T": (T, 4 * T), "0-3T": (0, 3 * T), "0-2T": (0, 2 * T)}.items():
        if b > t[-1] + 1e-6:
            continue
        wz, nz = zero_crossing_omega(t, d, a, b)
        f = fit_damped(t, d, a, b, omega_ref)
        res["windows"][name] = dict(
            omega_zero_crossing=wz, n_crossings=nz, err_zc=(wz / omega_ref - 1) if wz == wz else None,
            omega_fit=f["omega"], err_fit=f["omega"] / omega_ref - 1,
            amp_loss_per_period=1 - math.exp(-f["decay_rate"] * T), amplitude0=f["amplitude0"], rms_resid=f["rms_resid"])
    res["max_abs_D"] = float(np.max(np.abs(d)))
    res["finite"] = bool(np.all(np.isfinite(d)))
    return res


if __name__ == "__main__":
    pf, wref = sys.argv[1], float(sys.argv[2])
    print(json.dumps(analyse(pf, wref), indent=1))
