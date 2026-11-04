#!/usr/bin/env python3
"""RT-A checks (RTR pilot).  usage: rtr_check_rta.py law1|rotcheck|twist RUN_DIR [OUT.json].  Criteria fixed in RTR_PILOT_PLAN.md before the runs."""
import sys, os, json, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import rtr_reference as R


def read_log(path):
    """FEBio text log -> (times, list of arrays [n_rows, ncols] including the id column)."""
    T, D, cur = [], [], None
    for ln in open(path):
        if ln.startswith("*Time"):
            T.append(float(ln.split("=")[1])); cur = []; D.append(cur)
        elif ln.startswith("*") or not ln.strip():
            continue
        else:
            cur.append([float(x) for x in ln.split()])
    return np.array(T), [np.array(d) for d in D]


def law1(run):
    meta = json.load(open(run + "/model_meta.json")); Fm = np.array(meta["F"]); mu, K = R.rt.Params().G, R.rt.Params().kappa
    T, D = read_log(run + "/stress.txt"); out = []
    worst = 0.0; scale = 1.0
    for t, d in zip(T, D):
        if t == 0: continue
        F = np.eye(3) + t * (Fm - np.eye(3)); J = np.linalg.det(F); b = F @ F.T; bb = J ** (-2 / 3) * b
        dev = bb - np.trace(bb) / 3 * np.eye(3)
        sig = mu / J * dev + K / 2 * (J * J - 1) / J * np.eye(3)
        ex = np.array([sig[0, 0], sig[1, 1], sig[2, 2], sig[0, 1], sig[1, 2], sig[0, 2]])
        scale = max(scale, np.abs(ex).max())
        for row in d:
            diff = np.abs(row[1:7] - ex).max(); worst = max(worst, diff)
            out.append(dict(t=t, elem=int(row[0]), maxabs_diff=float(diff), J_febio=row[7], J_exact=J))
    res = dict(stage="RT-A law1", mu=mu, K=K, worst_abs_diff=worst, max_abs_stress=scale, rel=worst / scale, criterion=1e-8, PASS=bool(worst / scale <= 1e-8), rows=out)
    return res


def rotcheck(run):
    T, D = read_log(run + "/nodes.txt")
    worst = 0.0
    for t, d in zip(T, D):
        S, _, _ = R.ramp_S(t) if 0 < t < 1 else ((0.0 if t <= 0 else 1.0), 0, 0)
        a = 0.5 * S * math.sin(3 * t)
        for row in d:
            nid, x, y, ux, uy = row[:5]
            X, Y = x - ux, y - uy
            ex = (math.cos(a) - 1) * X - math.sin(a) * Y; ey = math.sin(a) * X + (math.cos(a) - 1) * Y
            worst = max(worst, abs(ux - ex), abs(uy - ey))
    return dict(stage="RT-A rotcheck", worst_abs_diff_m=worst, criterion=1e-10, n_times=len(T), PASS=bool(worst <= 1e-10))


def twist(run, Ns=24):
    meta = json.load(open(run + "/model_meta.json")); dz = 0.1; k2 = 2 * math.pi * dz
    T, Din = read_log(run + "/inner.txt"); _, Dout = read_log(run + "/outer.txt"); _, D0 = read_log(run + "/ray0.txt"); _, D45 = read_log(run + "/ray45.txt")
    amp = meta["amp"]; nrs = meta.get("nrs"); out = dict(t=[], TA=[], TAref=[], TBout=[], ur=[], urref=[], al=[], alref=[], ur45=[])
    guess = None
    for i, t in enumerate(T):
        if t <= 0: continue
        S = R.ramp_S(t)[0] if t < 1 else 1.0
        aA = amp * S * math.sin(3 * t)
        ref = R.static_twist(aA, 0.0, Ns=Ns)
        di = Din[i]; x, y, ux, uy, rx, ry = di[:, 1], di[:, 2], di[:, 3], di[:, 4], di[:, 5], di[:, 6]
        Tin = np.sum(x * ry - y * rx) / k2
        do = Dout[i]; Tout = np.sum(do[:, 1] * do[:, 6] - do[:, 2] * do[:, 5]) / k2
        d0 = D0[i]; r = np.hypot(d0[:, 1], d0[:, 2]); X0, Y0 = d0[:, 1] - d0[:, 3], d0[:, 2] - d0[:, 4]; R0 = np.hypot(X0, Y0)
        imid = int(np.argmin(np.abs(R0 - 0.75)))
        ur = r[imid] - R0[imid]; al = math.atan2(d0[imid, 2], d0[imid, 1]) - math.atan2(Y0[imid], X0[imid])
        d4 = D45[i]; r4 = np.hypot(d4[:, 1], d4[:, 2]); R4 = np.hypot(d4[:, 1] - d4[:, 3], d4[:, 2] - d4[:, 4]); i4 = int(np.argmin(np.abs(R4 - 0.75)))
        out["t"].append(t); out["TA"].append(Tin); out["TAref"].append(ref["T_A"]); out["TBout"].append(Tout)
        out["ur"].append(ur); out["urref"].append(float(R.eval_field(ref["model"], ref["u"], R0[imid]))); out["al"].append(al); out["alref"].append(float(R.eval_field(ref["model"], ref["alpha"], R0[imid])))
        out["ur45"].append(r4[i4] - R4[i4])
    o = {k: np.array(v) for k, v in out.items()}
    sgn = np.sign(np.sum(o["TA"] * o["TAref"])) or 1.0
    eT = np.max(np.abs(sgn * o["TA"] - o["TAref"])) / np.max(np.abs(o["TAref"]))
    eEq = np.max(np.abs(o["TA"] + o["TBout"])) / np.max(np.abs(o["TA"]))
    eU = np.max(np.abs(o["ur"] - o["urref"])) / np.max(np.abs(o["urref"])); eA = np.max(np.abs(o["al"] - o["alref"])) / np.max(np.abs(o["alref"]))
    ax = np.max(np.abs(o["ur"] - o["ur45"])) / np.max(np.abs(o["ur"]))
    return dict(stage="RT-A twist", mesh=run.split("_")[-1], n_steps=len(o["t"]), reaction_sign_vs_ref=float(sgn), E_max_TA=float(eT), equilibrium_inner_plus_outer_rel=float(eEq),
                E_max_ur_075=float(eU), E_max_alpha_075=float(eA), axisymmetry_ur_ray0_vs_ray45_rel=float(ax), max_TA_ref=float(np.max(np.abs(o["TAref"]))), max_TA_febio=float(np.max(np.abs(o["TA"]))),
                t=o["t"].tolist(), TA=o["TA"].tolist(), TAref=o["TAref"].tolist())


if __name__ == "__main__":
    kind, run = sys.argv[1], os.path.realpath(sys.argv[2])
    res = {"law1": law1, "rotcheck": rotcheck, "twist": twist}[kind](run)
    if len(sys.argv) > 3: json.dump(res, open(sys.argv[3], "w"), indent=1, default=float)
    print(json.dumps({k: v for k, v in res.items() if k not in ("rows", "t", "TA", "TAref")}, indent=1, default=float))
