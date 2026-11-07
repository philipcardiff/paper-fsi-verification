#!/usr/bin/env python3
import os, os as _os; HERE = _os.path.dirname(_os.path.abspath(__file__)); PKG = _os.path.dirname(HERE)   # package root; run directories of retained HPC outputs: env CAVITY_RUNS (default runs/cavityNext)
"""Read-only audit of a finished cavity FSI run from its FEBio .xplt (format 0x35, uncompressed) fields.  usage: cavity_xplt_audit.py RUNDIR OUT.json
Minimal xplt reader (chunk = u32 tag, u32 size, payload; format from FEBio v4.13 source FEBioPlot/FEBioPlotFile.{h,cpp}), then
reaction-force balance on the plate clamps, nodal-pressure force, viscous traction, pin-node effect, pressure along x. Nothing is tuned to a target."""
import sys, os, json, struct, numpy as np
BR = {0x01000000, 0x01010000, 0x01020000, 0x01021000, 0x01023000, 0x01024000, 0x01025000, 0x01026000, 0x01020001, 0x01040000, 0x01041000, 0x01041100, 0x01042000, 0x01042100, 0x01042101, 0x01042200,
      0x01043000, 0x01043100, 0x01043200, 0x01044000, 0x01044100, 0x01045000, 0x01045100, 0x01046000, 0x01046100, 0x01047000, 0x01047100, 0x01048000, 0x01048100, 0x01048200, 0x01050000,
      0x02000000, 0x02010000, 0x02020000, 0x02020001, 0x02020003, 0x02020300, 0x02020400, 0x02020500, 0x02020600, 0x02030000, 0x02040000, 0x01051000, 0x01052000}
# tags that are branches in the writer; PLT_STATE_VAR_DATA (0x02020003) holds either one raw float block (node data: tag 0) or (domain id, floats) sub-chunks.
def chunks(b, i, end):
    while i + 8 <= end:
        tag, sz = struct.unpack_from("<II", b, i); yield tag, i + 8, sz; i += 8 + sz
def read_xplt(path):
    b = open(path, "rb").read(); assert b[:4] == b"BEF\0"
    out = {"dic": {"node": [], "domain": []}, "states": [], "domains": []}; trunc = False
    def dic(i, sz, key):
        for t, p, s in chunks(b, i, i + sz):
            it = {}
            for t2, p2, s2 in chunks(b, p, p + s):
                if t2 == 0x01020002: it["type"] = struct.unpack_from("<I", b, p2)[0]
                elif t2 == 0x01020003: it["fmt"] = struct.unpack_from("<I", b, p2)[0]
                elif t2 == 0x01020004: it["name"] = b[p2:p2 + s2].split(b"\0")[0].decode()
            out["dic"][key].append(it)
    NV = {0: 1, 1: 3, 2: 6, 3: 3, 4: 21, 5: 9}
    for t, p, s in chunks(b, 4, len(b)):
        if t == 0x01000000:
            for t2, p2, s2 in chunks(b, p, p + s):
                if t2 == 0x01020000:
                    for t3, p3, s3 in chunks(b, p2, p2 + s2):
                        if t3 == 0x01023000: dic(p3, s3, "node")
                        if t3 == 0x01024000: dic(p3, s3, "domain")
        elif t == 0x01040000:
            for t2, p2, s2 in chunks(b, p, p + s):
                if t2 == 0x01041000:
                    for t3, p3, s3 in chunks(b, p2, p2 + s2):
                        if t3 == 0x01041200: a = np.frombuffer(b, np.float32, s3 // 4, p3).reshape(-1, 4); out["ids"] = a[:, 0].copy().view(np.int32); out["X0"] = a[:, 1:].astype(float)
                if t2 == 0x01042000:
                    for t3, p3, s3 in chunks(b, p2, p2 + s2):
                        if t3 != 0x01042100: continue
                        d = {"elems": []}
                        for t4, p4, s4 in chunks(b, p3, p3 + s3):
                            if t4 == 0x01042101:
                                for t5, p5, s5 in chunks(b, p4, p4 + s4):
                                    if t5 == 0x01042102: d["etype"] = struct.unpack_from("<I", b, p5)[0]
                                    if t5 == 0x01042103: d["part"] = struct.unpack_from("<I", b, p5)[0]
                            if t4 == 0x01042200:
                                for t5, p5, s5 in chunks(b, p4, p4 + s4):
                                    d["elems"].append(np.frombuffer(b, np.int32, s5 // 4, p5).copy())
                        out["domains"].append(d)
        elif t == 0x02000000:
            st = {"vars": {}, "complete": True}
            if p + s > len(b): trunc = True; break
            for t2, p2, s2 in chunks(b, p, p + s):
                if t2 == 0x02010000:
                    for t3, p3, s3 in chunks(b, p2, p2 + s2):
                        if t3 == 0x02010002: st["time"] = struct.unpack_from("<f", b, p3)[0]
                        if t3 == 0x02010003: st["flag"] = struct.unpack_from("<I", b, p3)[0]
                if t2 == 0x02020000:
                    for t3, p3, s3 in chunks(b, p2, p2 + s2):
                        kind = {0x02020300: "node", 0x02020400: "domain"}.get(t3)
                        if kind is None: continue
                        for t4, p4, s4 in chunks(b, p3, p3 + s3):   # PLT_STATE_VARIABLE
                            vid = None
                            for t5, p5, s5 in chunks(b, p4, p4 + s4):
                                if t5 == 0x02020002: vid = struct.unpack_from("<I", b, p5)[0]
                                if t5 == 0x02020003:
                                    it = out["dic"][kind][vid - 1]; nv = NV[it["type"]]
                                    dd = {}; _dom(b, p5, s5, nv, dd)   # sub-chunks (region id, floats); node data uses id 0 -> key -1
                                    st["vars"][it["name"]] = dd[-1] if kind == "node" else dd
            out["states"].append(st)
    out["truncated"] = trunc; return out
def _dom(b, p, s, nv, dct):
    for d, p6, s6 in chunks(b, p, p + s): dct[d - 1] = np.frombuffer(b, np.float32, s6 // 4, p6).reshape(-1, nv).astype(float)

# ---------------------------------------------------------------- analysis
import re
HEX_EDGES = [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4), (0, 4), (1, 5), (2, 6), (3, 7)]
def read_feb(path):
    f = open(path, encoding="latin-1").read()
    ns = {m.group(1): np.array(m.group(2).split(","), int) for m in re.finditer(r'<NodeSet name="(\w+)">([^<]+)</NodeSet>', f)}
    sb = f[f.index('<Surface name="Interface">'):]; sb = sb[:sb.index("</Surface>")]
    quads = np.array([m.group(1).split(",") for m in re.finditer(r"<quad4 id=\"\d+\">([^<]+)</quad4>", sb)], int)
    return ns, quads
def read_mon_block(path, tsel):
    """one time block of a Mon_*.txt (id then values), streaming; returns {id: tuple}"""
    out, on = {}, False
    for ln in open(path):
        if ln.startswith("*Time"): on = abs(float(ln.split("=")[1]) - tsel) < 1e-9; continue
        if on and not ln.startswith("*") and ln.strip(): f = ln.split(); out[int(f[0])] = tuple(map(float, f[1:]))
        elif on and ln.startswith("*Step") and out: break
    return out
def poly_force(X, Y, P):   # per unit depth, pressure on polyline, force of fluid on plate (outward normal of plate = +y side): F = -int p n dl, n = (-dy, dx)/L
    dx, dy = np.diff(X), np.diff(Y); pb = 0.5 * (P[1:] + P[:-1]); return float(-np.sum(pb * dx)), float(np.sum(pb * dy))   # (Fy, Fx): Fx = -int p n_x dl = +int p dy
def audit(rd, tsel=None):
    meta = json.load(open(os.path.join(rd, "model_meta.json"))); tag, K, mu, TZ = meta["tag"], meta["k"], meta["mu"], meta["thickness_z"]
    o = read_xplt(os.path.join(rd, tag + ".xplt")); ns, quads = read_feb(os.path.join(rd, tag + ".feb"))
    S = o["states"]; res = {"tag": tag, "states": [{"time": s["time"], "flag": s["flag"]} for s in S], "n_nodes": len(o["ids"]), "TZ": TZ, "K": K, "mu": mu}
    assert (o["ids"] == np.arange(1, len(o["ids"]) + 1)).all()
    X0 = o["X0"]; plate_el, fl_el = o["domains"][0]["elems"], o["domains"][1]["elems"]
    conn = np.array([e[1:] for e in fl_el]); eid = np.array([e[0] for e in fl_el])   # element node lists are 0-based node indices
    fl_nodes = np.unique(conn); pl_nodes = np.unique(np.array([e[1:] for e in plate_el]))
    comp = [k for k, s in enumerate(S) if s["flag"] == 0]; kl = comp[-1] if tsel is None else [k for k in comp if S[k]["time"] == tsel][0]; res["state_used"] = {"index": kl, "time": S[kl]["time"], "flag": S[kl]["flag"], "note": "last state with flag 0 (complete, converged step); later states have flag 2 = written on solver failure"}
    def one(k):
        s = S[k]; V = s["vars"]; U = V["displacement"]; ef = V["fluid dilatation"][:, 0]; Xd = X0 + U; r = {"time": s["time"], "flag": s["flag"]}
        # reaction forces
        R = V["reaction forces"]; r["reaction_forces"] = {"max_abs_all_nodes": float(np.abs(R).max()), "sum_clamp_xyz": [float(R[np.r_[ns["ClampL"], ns["ClampR"]] - 1].sum(0)[j]) for j in range(3)], "n_nonzero_nodes": int((np.abs(R).sum(1) > 0).sum())}
        # interface z0 nodes ordered by x
        ids = ns["InterfaceZ0"] - 1; ids = ids[np.argsort(X0[ids, 0])]; Xi, Yi = Xd[ids, 0], Xd[ids, 1]; Pn = -K * ef[ids]
        Fy, Fx = poly_force(Xi, Yi, Pn); r["pressure_only"] = {"Fy": Fy, "Fx": Fx, "n_nodes": len(ids), "p_min": float(Pn.min()), "p_max": float(Pn.max()), "p_mean": float(Pn.mean()), "dtheta_max_slope": float(np.abs(np.diff(Yi) / np.diff(Xi)).max())}
        r["interface_p_along_x"] = {"x_ref": X0[ids, 0].tolist(), "x_def": Xi.tolist(), "y_def": Yi.tolist(), "p": Pn.tolist()}
        idx = np.arange(len(ids)); ev, od = Pn[idx % 2 == 0], Pn[idx % 2 == 1]; r["oddeven_amp_even_minus_odd"] = float(ev.mean() - od.mean())
        # element-wise pressure vs -K*mean nodal ef
        Pe = V["fluid pressure"][1][:, 0]; r["pressure_elem_vs_nodal"] = {"max_abs_diff_Pa": float(np.abs(Pe + K * ef[conn].mean(1)).max()), "max_abs_p_Pa": float(np.abs(Pe).max())}
        # pin
        pin = ns["CornerPin"] - 1; lo, hi = ids[0], ids[-1]; pn = {"pinned_ef_x1e6": [float(ef[j] * 1e6) for j in pin], "pinned_p": [float(-K * ef[j]) for j in pin]}
        pe = Pn.copy(); xs = Xi
        pe[0] = Pn[1] + (Pn[1] - Pn[2]) * (xs[0] - xs[1]) / (xs[1] - xs[2]); pe[-1] = Pn[-2] + (Pn[-2] - Pn[-3]) * (xs[-1] - xs[-2]) / (xs[-2] - xs[-3])
        Fy2, Fx2 = poly_force(Xi, Yi, pe); pn.update({"p_left_nodes_0to4": Pn[:5].tolist(), "p_right_nodes_last5": Pn[-5:].tolist(), "p_pin_left_right": [float(Pn[0]), float(Pn[-1])], "p_extrap_left_right": [float(pe[0]), float(pe[-1])],
              "Fy_extrapolated": Fy2, "Fx_extrapolated": Fx2, "dFy": Fy2 - Fy, "dFx": Fx2 - Fx, "note": "bound on pin effect only (linear extrapolation of the two nearest non-pinned interface nodes), NOT a corrected result"}); r["pin"] = pn
        # facets and adjacent elements
        z0 = set((ns["InterfaceZ0"]).tolist()); n2e = {}
        for j, c in enumerate(conn):
            for n in c: n2e.setdefault(int(n), []).append(j)
        rows = []
        W = V["relative fluid velocity"][1]
        for q in quads:
            za = [n for n in q if n in z0]; assert len(za) == 2; a, b = sorted(za, key=lambda n: X0[n - 1, 0]); a -= 1; b -= 1
            fq = set(int(n) - 1 for n in q); E1 = [j for j in n2e[a] if fq <= set(conn[j].tolist())]; assert len(E1) == 1; E1 = E1[0]; c1 = conn[E1]
            pos = {int(n): k for k, n in enumerate(c1)}; opp = []
            for n in fq:
                k = pos[n]; nb = [c1[m if l == k else l] for l, m in HEX_EDGES if k in (l, m)]; opp.append([int(x) for x in nb if int(x) not in fq][0])
            oppmap = {}
            for n in fq:
                k = pos[n]; nb = [int(c1[m if l == k else l]) for l, m in HEX_EDGES if k in (l, m)]; oppmap[n] = [x for x in nb if x not in fq][0]
            E2 = [j for j in set(n2e[opp[0]]) if j != E1 and set(opp) <= set(conn[j].tolist())]; E2 = E2[0] if E2 else None
            d = Xd[b, :2] - Xd[a, :2]; L = np.hypot(*d); n = np.array([-d[1], d[0]]) / L; mid = 0.25 * Xd[list(fq)].mean(0)[:2] * 4 if False else Xd[list(fq)][:, :2].mean(0)
            cen = lambda E: Xd[conn[E]][:, :2].mean(0)
            s1 = float((cen(E1) - mid) @ n); w1 = W[E1][:2]; row = {"L": float(L), "n": n, "s1": s1, "w1": w1, "E1": E1, "dxy_ratio": None}
            if E2 is not None: row["s2"] = float((cen(E2) - mid) @ n); row["w2"] = W[E2][:2]; row["E2"] = E2; row["a"] = a; row["b"] = b; row["opp"] = oppmap; row["fq"] = fq
            row["xm"] = float(mid[0]); rows.append(row)
        ratio = []
        def visc(mode):
            Fv = np.zeros(2); out = []; ratio.clear()
            for r_ in rows:
                if mode == 1 or "s2" not in r_: a_ = r_["w1"] / r_["s1"]
                else: s1_, s2_ = r_["s1"], r_["s2"]; a_ = (r_["w1"] * s2_ ** 2 - r_["w2"] * s1_ ** 2) / (s1_ * s2_ * (s2_ - s1_))
                t = mu * (a_ + r_["n"] * (a_ @ r_["n"])); Fv += t * r_["L"]; out.append(t); ratio.append(float(abs(a_ @ r_["n"]) / (np.linalg.norm(a_) + 1e-30)))
            return Fv, np.array(out)
        v1, t1 = visc(1); v2, t2 = visc(2)
        # pressure force facet-wise (per unit depth) with facet-mean nodal pressure, same as polyline
        r["viscous"] = {"method_1pt": "a=w_c1/s1 (w=0 at wall), first order", "method_2pt": "quadratic w(s)=a s+b s^2 through element-average values of the 1st and 2nd element row and w(wall)=0, second order in h on smooth fields",
                        "traction": "t = mu (grad w + grad w^T).n = mu (a + n (a.n)), grad w = a (x) n at a no-slip wall; Fvisc = sum t L (per unit depth)", "Fx_visc_1pt": float(v1[0]), "Fy_visc_1pt": float(v1[1]), "Fx_visc_2pt": float(v2[0]), "Fy_visc_2pt": float(v2[1]),
                        "s1_first_row_centroid_distance_range": [min(r_["s1"] for r_ in rows), max(r_["s1"] for r_ in rows)], "max_abs_a_dot_n_over_abs_a_2pt": max(ratio), "n_facets": len(rows),
                        "tx_2pt_vs_x": [[r_["xm"] for r_ in rows], t2[:, 0].tolist()], "ty_2pt_vs_x": t2[:, 1].tolist()}
        r["total"] = {"Fy_1pt": Fy + float(v1[1]), "Fx_1pt": Fx + float(v1[0]), "Fy_2pt": Fy + float(v2[1]), "Fx_2pt": Fx + float(v2[0])}
        r["_rows"] = rows; r["_conn"] = conn; r["_W"] = W; r["_Xd"] = Xd
        return r
    res["last_complete"] = one(kl); prev = [k for k in comp if k < kl]; res["previous_complete"] = one(prev[-1]) if prev else None
    res["failed_state"] = one(len(S) - 1) if S[-1]["flag"] != 0 else None
    for key in ("last_complete", "previous_complete", "failed_state"):
        if res[key] and key != "last_complete":
            for q in ("_rows", "_conn", "_W", "_Xd"): res[key].pop(q, None)
    lc = res["last_complete"]; rows = lc.pop("_rows"); cn = lc.pop("_conn"); Wl = lc.pop("_W"); Xdl = lc.pop("_Xd")
    mall = os.path.join(rd, "Mon_all.txt")
    if os.path.exists(mall) and os.path.exists(mall):
        M = read_mon_block(mall, S[kl]["time"]); w_n = lambda nd: np.array(M[nd + 1][2:4])   # Mon_all columns after id: ux uy wx wy ef
        d1, dF, anr = [], np.zeros(2), []
        for r_ in rows:
            E1 = r_["E1"]; z0n = [int(x) for x in cn[E1] if abs(X0[x, 2]) < 1e-9]; d1.append(float(np.abs(np.mean([w_n(x) for x in z0n], 0) - Wl[E1][:2]).max()))
        lc["viscous_nodal_check"] = {"time": S[kl]["time"], "elem_avg_vs_mean_of_4_z0_nodes_Mon_all_max_abs_diff": float(max(d1)), "max_abs_wx_wy_element": float(np.abs(Wl[:, :2]).max())}
        Fn = np.zeros(2)
        for r_ in rows:
            if "E2" not in r_: continue
            ac = []
            for nd in (r_["a"], r_["b"]):
                o1 = r_["opp"][nd]; E2 = r_["E2"]; c2 = cn[E2]; f1 = set(r_["opp"].values()); pos2 = {int(x): k for k, x in enumerate(c2)}
                k = pos2[o1]; nb = [int(c2[m if l == k else l]) for l, m in HEX_EDGES if k in (l, m)]; o2 = [x for x in nb if x not in f1][0]
                n = r_["n"]; s1_ = (Xdl[o1, :2] - Xdl[nd, :2]) @ n; s2_ = (Xdl[o2, :2] - Xdl[nd, :2]) @ n; w1_, w2_ = w_n(o1), w_n(o2)
                a_ = (w1_ * s2_ ** 2 - w2_ * s1_ ** 2) / (s1_ * s2_ * (s2_ - s1_)); ac.append(a_)
            a_ = 0.5 * (ac[0] + ac[1]); anr.append(float(abs(a_ @ r_["n"]) / (np.linalg.norm(a_) + 1e-30))); Fn += mu * (a_ + r_["n"] * (a_ @ r_["n"])) * r_["L"]
        lc["viscous_nodal_check"].update({"Fx_visc_nodal_3pt_quadratic": float(Fn[0]), "Fy_visc_nodal_3pt_quadratic": float(Fn[1]), "max_abs_a_dot_n_over_abs_a": float(max(anr)), "method": "w at wall = 0; quadratic through the 1st and 2nd nodes of the element edge normal to the wall (Mon_all nodal wx,wy, z=0 layer), trapezoid over the facet's two nodes; second order"})

    # ---- domain-wide pressure sections, inlet/outlet
    s = S[kl]; V = s["vars"]; ef = V["fluid dilatation"][:, 0]; Xd = X0 + V["displacement"]; K_ = K
    z0m = np.abs(X0[:, 2]) < 1e-9; fl0 = np.zeros(len(X0), bool); fl0[fl_nodes] = True; sel = z0m & fl0
    out = {"x_range_fluid_ref": [float(X0[sel, 0].min()), float(X0[sel, 0].max())], "y_range_fluid_ref": [float(X0[sel, 1].min()), float(X0[sel, 1].max())]}
    for nm, f in (("inlet", "Mon_inlet.txt"), ("outlet", "Mon_outlet.txt")):
        m = read_mon_block(os.path.join(rd, f), S[kl]["time"]); idn = np.array(sorted(m)) - 1; p = -K * ef[idn]; y = Xd[idn, 1]; o_ = np.argsort(y)
        out[nm] = {"n_nodes": len(idn), "x_ref": float(X0[idn[0], 0]), "p_mean_trapz_y": float(np.trapz(p[o_], y[o_]) / (y.max() - y.min())), "p_min": float(p.min()), "p_max": float(p.max()), "max_abs_p": float(np.abs(p).max()),
                   "max_abs_ef": float(np.abs(ef[idn]).max()), "p_nodes_by_y": p[o_].tolist(), "y": y[o_].tolist(), "wx_mean_flux": float(np.trapz([m[i + 1][0] for i in idn[o_]], y[o_]))}
    out["xplt_vs_Mon_ef_max_abs"] = {"iface": float(max(abs(ef[i - 1] - v[4]) for i, v in read_mon_block(os.path.join(rd, "Mon_iface.txt"), S[kl]["time"]).items()))}
    out["Mon_iface_w_max_abs"] = float(max(max(abs(v[2]), abs(v[3])) for v in read_mon_block(os.path.join(rd, "Mon_iface.txt"), S[kl]["time"]).values()))
    # channel band (reference y >= 1): one value per mesh column = trapezoid mean over y of nodal -K ef (z = 0 layer); cavity band reported at 2 stations
    ch = np.where(sel & (X0[:, 1] >= 1.0 - 1e-9))[0]; secs = []
    mf = meta["mesh_factor"]; hw = 0.1 / mf; dsl = 0.1 / mf; NS = int(round(1.0 / dsl)) + 1
    for x in np.arange(0.25, 13.76, 0.25):   # x-bin of width 2*0.1/mesh_factor around x (reference x); y-slab (0.1/mesh_factor) average first so non-uniform y node density does not bias; deformed mesh motion in the channel is O(1e-3)
        jj = ch[np.abs(X0[ch, 0] - x) <= hw]; slab = np.clip(((X0[jj, 1] - 1.0) / dsl).round().astype(int), 0, NS - 1); p = -K * ef[jj]
        sl = [p[slab == k].mean() for k in range(NS) if (slab == k).any()]
        secs.append({"x_ref": float(x), "n_nodes": int(len(jj)), "n_slabs": len(sl), "p_mean_slabavg": float(np.mean(sl)), "p_min": float(p.min()), "p_max": float(p.max())})
    def fit(a, b):
        r_ = [q for q in secs if a <= q["x_ref"] <= b and q["n_slabs"] >= 0.7 * NS]; xx = np.array([q["x_ref"] for q in r_]); pp = np.array([q["p_mean_slabavg"] for q in r_]); c = np.polyfit(xx, pp, 1)
        return {"x_range": [a, b], "n_cols": len(r_), "dp_dx": float(c[0]), "p_at_xmid": float(np.polyval(c, 0.5 * (a + b))), "max_resid": float(np.abs(pp - np.polyval(c, xx)).max())}
    out["channel_gradient_fits"] = {"upstream_0.3_1.7": fit(0.3, 1.7), "above_cavity_2.3_5.7": fit(2.3, 5.7), "downstream_6.3_13.7": fit(6.3, 13.7), "Poiseuille_plane_expected_dp_dx_H1_U1": -12 * mu * 1.0 / 1.0, "note": "expected value assumes channel height 1 and mean inlet speed 1 (as implied by the 0.12 Pa/m quoted)"}
    out["channel_sections_all"] = secs
    lc["domain"] = out
    out["outlet_all_42_nodes_max_abs_ef"] = float(np.abs(ef[ns["Outlet"] - 1]).max())
    A = {}; t_ = None; aid = meta["point_A_node"] - 1
    for ln in open(os.path.join(rd, "pointA.txt")):
        if ln.startswith("*Time"): t_ = float(ln.split("=")[1])
        elif not ln.startswith("*") and ln.strip(): A[t_] = list(map(float, ln.split()[1:]))
    res["validation"] = {"pointA_xplt_vs_pointA_txt": [{"time": q["time"], "uy_xplt": float(q["vars"]["displacement"][aid, 1]), "uy_txt": A[q["time"]][3], "abs_diff": float(abs(q["vars"]["displacement"][aid, 1] - A[q["time"]][3]))} for q in S if q["time"] in A],
        "pointA_last_logged_time": max(A), "uy_A_last_logged": A[max(A)][3], "p_eq_minusK_ef": "pressure_elem_vs_nodal: element 'fluid pressure' vs -K*mean of 8 nodal ef", "uz_max_abs": float(np.abs(V["displacement"][:, 2]).max())}
    return res, o, rows
def main():
    rd = sys.argv[1]; res, o, rows = audit(rd)
    json.dump(res, open(sys.argv[2], "w"), indent=1, default=lambda x: x.tolist() if hasattr(x, "tolist") else float(x))
    lc = res["last_complete"]; print(json.dumps({k: lc[k] for k in ("time", "reaction_forces", "pressure_only", "total", "pin")}, indent=1, default=float)); print(res["state_used"])
    print({k: v for k, v in lc["viscous"].items() if "vs_x" not in k}); d = lc["domain"]; print(json.dumps({k: d[k] for k in ("x_range_fluid_ref", "y_range_fluid_ref", "channel_gradient_fits", "xplt_vs_Mon_ef_max_abs", "Mon_iface_w_max_abs")}, default=float)); print({k: {a: b for a, b in d[k].items() if "nodes_by_y" not in a and a != "y"} for k in ("inlet", "outlet")})
if __name__ == "__main__": main()
