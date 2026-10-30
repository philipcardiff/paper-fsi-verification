#!/usr/bin/env python3
import os, os as _os; HERE = _os.path.dirname(_os.path.abspath(__file__)); PKG = _os.path.dirname(HERE)   # package root; run directories of retained HPC outputs: env CAVITY_RUNS (default runs/cavityNext)
"""B3: independent force check by momentum balance on control volumes, from FEBio xplt fields (z=0 layer, per unit depth).
usage: b3_momentum_balance.py RUNDIR TIME OUT.json   (reads, never writes, the run directory; TIME = xplt state time)
Fields: nodal p = -K*ef; element-average (centroid) velocity; nodal inlet/outlet velocity from Mon_inlet/Mon_outlet.
Sign convention: forces ON THE FLUID by each boundary piece; the interface piece is compared with minus the force of the fluid on the plate.
Control volumes: WHOLE fluid domain; CAVITY (box x in [2,6], reference y<0.75; mouth line y=0.75: the mesh has vertical cavity walls up to y=0.75 and flared lips above).
Nothing is tuned; first-order viscous terms on open boundaries are bounds, listed as such."""
import sys, os, json, importlib.util, numpy as np
spec = importlib.util.spec_from_file_location("aud", HERE + "/cavity_xplt_audit.py"); aud = importlib.util.module_from_spec(spec); spec.loader.exec_module(aud)
rd, tsel, outp = sys.argv[1], float(sys.argv[2]), sys.argv[3]
meta = json.load(open(os.path.join(rd, "model_meta.json"))); tag, K, mu, TZ = meta["tag"], meta["k"], meta["mu"], meta["thickness_z"]; rho = 1.0
o = aud.read_xplt(os.path.join(rd, tag + ".xplt")); ns, quads_if = aud.read_feb(os.path.join(rd, tag + ".feb"))
S = o["states"]; k = [i for i, s in enumerate(S) if abs(s["time"] - tsel) < 1e-6][0]; V = S[k]["vars"]
X0 = o["X0"]; Xd = X0 + V["displacement"]; ef = V["fluid dilatation"][:, 0]; P = -K * ef
conn = np.array([e[1:] for e in o["domains"][1]["elems"]]); W = V["relative fluid velocity"][1][:, :2]
# 2-D quads from the z=0 face (hex8: nodes 0-3 and 4-7 are the two cyclic faces)
Q = []
for c in conn:
    f = c[:4] if np.all(np.abs(X0[c[:4], 2]) < 1e-9) else c[4:]
    assert np.all(np.abs(X0[f, 2]) < 1e-9); Q.append(f)
Q = np.array(Q); nq = len(Q); cen = Xd[Q][:, :, :2].mean(1); cen0 = X0[Q][:, :, :2].mean(1)
edge_el = {}
for j, q in enumerate(Q):
    for i in range(4): edge_el.setdefault(tuple(sorted((int(q[i]), int(q[(i + 1) % 4])))), []).append(j)
ifn = set((ns["InterfaceZ0"] - 1).tolist())
def region(mask):
    """boundary edges of the quad set 'mask' with their quad, edge index"""
    out = []
    for j in np.where(mask)[0]:
        q = Q[j]
        for i in range(4):
            a, b = int(q[i]), int(q[(i + 1) % 4]); nb = [e for e in edge_el[tuple(sorted((a, b)))] if e != j]
            if (not nb) or (not mask[nb[0]]): out.append((j, i, a, b, nb[0] if nb else None))
    return out
def kind(a, b, cavity):
    xa, ya, xb, yb = X0[a, 0], X0[a, 1], X0[b, 0], X0[b, 1]
    if a in ifn and b in ifn: return "interface"
    if abs(xa) < 1e-9 and abs(xb) < 1e-9: return "inlet"
    if abs(xa - 14) < 1e-9 and abs(xb - 14) < 1e-9: return "outlet"
    if cavity and abs(ya - 0.75) < 1e-9 and abs(yb - 0.75) < 1e-9 and 2 - 1e-9 <= min(xa, xb) and max(xa, xb) <= 6 + 1e-9: return "mouth"
    if abs(ya - 2) < 1e-9 and abs(yb - 2) < 1e-9: return "wall_top"
    if abs(ya - 1) < 1e-9 and abs(yb - 1) < 1e-9: return "wall_botUp" if max(xa, xb) <= 2 + 1e-9 else ("wall_botDown" if min(xa, xb) >= 6 - 1e-9 else "wall_cavity")
    return "wall_cavity"
def inout_nodal(fn):
    m = aud.read_mon_block(os.path.join(rd, fn), tsel); return {i - 1: v for i, v in m.items()}
Min, Mout = inout_nodal("Mon_inlet.txt"), inout_nodal("Mon_outlet.txt")
def balance(mask, cavity):
    E = region(mask); acc = {}; ker = lambda: {"F_p": np.zeros(2), "F_v": np.zeros(2), "F_v1": np.zeros(2), "flux": np.zeros(2), "L": 0.0, "n": 0}
    nodeu = lambda nd: np.array((Min.get(nd) or Mout.get(nd))[:2]) if (nd in Min or nd in Mout) else None
    for (j, i, a, b, nbq) in E:
        kd = kind(a, b, cavity); d = Xd[b, :2] - Xd[a, :2]; L = float(np.hypot(*d)); t = d / L
        n_a = np.array([t[1], -t[0]])                      # one of the two normals
        mid = 0.5 * (Xd[a, :2] + Xd[b, :2]); nout = n_a if (mid - cen[j]) @ n_a > 0 else -n_a   # outward from the CV
        A = acc.setdefault(kd, ker()); A["L"] += L; A["n"] += 1
        pbar = 0.5 * (P[a] + P[b]); Fp = -pbar * nout * L; A["F_p"] += Fp   # force ON the fluid, pressure part
        if kd.startswith("wall") or kd == "interface":
            nin = -nout   # into the fluid; s measured into the fluid; wall velocity zero
            q = Q[j]; opp = [int(q[(i + 2) % 4]), int(q[(i + 3) % 4])]; e2 = [e for e in edge_el[tuple(sorted(opp))] if e != j]
            s1 = float((cen[j] - mid) @ nin); w1 = W[j]
            if e2:
                c2 = e2[0]; s2 = float((cen[c2] - mid) @ nin); w2 = W[c2]; a_ = (w1 * s2 ** 2 - w2 * s1 ** 2) / (s1 * s2 * (s2 - s1))
            else: a_ = w1 / s1
            tfw = mu * (a_ + nin * (a_ @ nin))         # traction of the fluid on the wall
            A["F_v"] += -tfw * L                        # force on the fluid = minus
            a1 = w1 / s1; A["F_v1"] += -mu * (a1 + nin * (a1 @ nin)) * L
        elif kd in ("inlet", "outlet"):
            Mx = Min if kd == "inlet" else Mout; ua, ub = np.array(Mx[a][:2]), np.array(Mx[b][:2]); ubar = 0.5 * (ua + ub)
            A["flux"] += rho * ubar * float(ubar @ nout) * L if False else rho * 0.5 * (ua * float(ua @ nout) + ub * float(ub @ nout)) * L
            # first-order viscous bound on the open boundary: d/dn from nodal boundary value to the first centroid; tangential derivative of the nodal values
            s1 = float((cen[j] - mid) @ -nout); dn = (W[j] - ubar) / s1 if s1 != 0 else np.zeros(2); dt = (ub - ua) / L
            # traction mu*(grad u + grad u^T) n_out ; grad u = dn (x) (-nout)... assemble with local frame (n_out, t)
            g = np.outer(-dn, nout) + np.outer(dt, t)       # g[i,j] = d u_i / d x_j
            Fv = mu * ((g + g.T) @ nout) * L; A["F_v"] += Fv
        elif kd == "mouth":
            # elements on the other side of the edge
            c2 = nbq; ub = 0.5 * (W[j] + W[c2]); d_ = cen[c2] - cen[j]; dist = float(d_ @ -nout)   # distance between the two centroids along the outward normal... (c2 lies outside the CV)
            dn = (W[c2] - W[j]) / (np.hypot(*(cen[c2] - cen[j])) + 1e-30) * np.sign(float((cen[c2] - cen[j]) @ nout) + 1e-30)
            g = np.outer(dn, nout)    # d u_i / d n with only the normal derivative; tangential derivative of v neglected (bound reported)
            A["F_v"] += mu * ((g + g.T) @ nout) * L; A["flux"] += rho * ub * float(ub @ nout) * L
    return acc
def summarise(acc, label):
    r = {"label": label, "pieces": {}}
    flux = sum((v["flux"] for v in acc.values()), np.zeros(2)); r["momentum_flux_out_xy"] = flux.tolist()
    tot = np.zeros(2)
    for kd, v in acc.items():
        r["pieces"][kd] = {"n_edges": v["n"], "length": v["L"], "F_pressure_on_fluid": v["F_p"].tolist(), "F_viscous_on_fluid": v["F_v"].tolist(), "F_viscous_1pt_on_fluid": v["F_v1"].tolist(), "F_total_on_fluid": (v["F_p"] + v["F_v"]).tolist(), "momentum_flux_out": v["flux"].tolist()}; tot += v["F_p"] + v["F_v"]
    r["sum_boundary_forces_on_fluid"] = tot.tolist(); r["residual_flux_minus_forces"] = (flux - tot).tolist()
    I = acc["interface"]; r["interface_direct_force_on_plate_pressure"] = (-I["F_p"]).tolist(); r["interface_direct_force_on_plate_total"] = (-(I["F_p"] + I["F_v"])).tolist()
    r["interface_inferred_force_on_plate_total"] = (-( flux - (tot - (I["F_p"] + I["F_v"])))).tolist()
    return r
res = {"tag": tag, "time": tsel, "state_flag": S[k]["flag"], "units": "per unit depth (z=0 layer line integrals), N/m", "K": K, "mu": mu}
res["whole_domain"] = summarise(balance(np.ones(nq, bool), False), "whole fluid domain")
res["cavity"] = summarise(balance(cen0[:, 1] < 0.75, True), "cavity box x in [2,6], reference y<0.75, mouth line y=0.75")
# mass balance
for nm, M in (("whole_domain", None),):
    pass
json.dump(res, open(outp, "w"), indent=1)
for key in ("whole_domain", "cavity"):
    r = res[key]; print(key, "flux", np.round(r["momentum_flux_out_xy"], 5), "resid", np.round(r["residual_flux_minus_forces"], 5))
    for kd, v in r["pieces"].items(): print("  ", kd, v["n_edges"], "Fp", np.round(v["F_pressure_on_fluid"], 5), "Fv", np.round(v["F_viscous_on_fluid"], 5), "Fv1", np.round(v["F_viscous_1pt_on_fluid"], 5), "flux", np.round(v["momentum_flux_out"], 5))
    print("   interface direct on plate (p only)", np.round(r["interface_direct_force_on_plate_pressure"], 5), "total", np.round(r["interface_direct_force_on_plate_total"], 5), "inferred", np.round(r["interface_inferred_force_on_plate_total"], 5))
