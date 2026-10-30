#!/usr/bin/env python3
"""B3: reconstruct the four interface forces F1..F4 from stored FEBio xplt fields (no simulation).  usage: CAVITY_RUNS=<dir with C2_fsi_pilot_f1/, C3_fsi_M2/, C3_fsi_M4/> b3_force_reconstruction.py OUT.json [FACET_DIR]   (needs the retained .xplt/.feb/model_meta.json: NOT reproducible from the package alone; FACET_DIR receives results/raw/facets/*.npz)
F1 = FEFluidFSITraction::LoadVector pressure part: p from SURFACE NODAL ef bilinearly interpolated at the 2x2 facet Gauss points (GetFluidDilatation), times gt = gr x gs (deformed), force on plate = -p gt.
     The viscous part (-sv.gt, sv = Gauss mean of the adjacent element viscous stress) needs element-Gauss velocity gradients: NOT stored; bounded by audit values.
F2 = FEPlotFluidSurfaceForce / pressure definition: adjacent-element Gauss-mean of m_sf; pressure part = element pressure (xplt 'fluid pressure') x area vector.
F3 = independent estimate from the xplt ELEMENT pressures of rows 1 and 2 (adjacent elements E1,E2) linearly extrapolated to the wall (centroid distances s1,s2), x area vector.
F4 = pair-row diagnostic: nodal (p_row0 + p_row1)/2, bilinear, Gauss integrated (non-physical until validated).
Depth normalisation: forces on the whole facet depth divided by TZ -> N/m.  Nothing is tuned."""
import sys, os, json, importlib.util, numpy as np
RUNS = os.environ.get("CAVITY_RUNS", "runs/cavityNext")   # directory holding the retained HPC run folders
spec = importlib.util.spec_from_file_location("aud", os.path.join(os.path.dirname(os.path.abspath(__file__)), "cavity_xplt_audit.py")); aud = importlib.util.module_from_spec(spec); spec.loader.exec_module(aud)
HEX_EDGES = aud.HEX_EDGES
G = 1 / np.sqrt(3); GP = [(-G, -G), (G, -G), (G, G), (-G, G)]
def shp(r, s): return np.array([(1 - r) * (1 - s), (1 + r) * (1 - s), (1 + r) * (1 + s), (1 - r) * (1 + s)]) / 4
def dshp(r, s): return (np.array([-(1 - s), (1 - s), (1 + s), -(1 + s)]) / 4, np.array([-(1 - r), -(1 + r), (1 + r), (1 - r)]) / 4)
def gt(x4, r, s):
    dr, ds = dshp(r, s); return np.cross(dr @ x4, ds @ x4)   # x4: 4x3 deformed coords; facet order (i,k0),(i,k1),(i+1,k1),(i+1,k0): r~z, s~x, gr x gs = +y (out of the plate)
def run(tag, rd, tsel=None):
    meta = json.load(open(f"{rd}/model_meta.json")); K, TZ = meta["k"], meta["thickness_z"]
    o = aud.read_xplt(f"{rd}/{tag}.xplt"); ns, quads = aud.read_feb(f"{rd}/{tag}.feb"); S = o["states"]
    comp = [k for k, s in enumerate(S) if s["flag"] == 0]; k = comp[-1] if tsel is None else [i for i in range(len(S)) if abs(S[i]["time"] - tsel) < 1e-6][0]
    V = S[k]["vars"]; X0 = o["X0"]; Xd = X0 + V["displacement"]; ef = V["fluid dilatation"][:, 0]; Pe = V["fluid pressure"][1][:, 0]
    conn = np.array([e[1:] for e in o["domains"][1]["elems"]]); n2e = {}
    for j, c in enumerate(conn):
        for n in c: n2e.setdefault(int(n), []).append(j)
    z0 = set(ns["InterfaceZ0"].tolist()); rows = []
    for q in quads:
        q0 = [int(n) - 1 for n in q]; fq = set(q0); za = [n for n in q if n in z0]; a, b = sorted(za, key=lambda n: X0[n - 1, 0]); a -= 1; b -= 1
        E1 = [j for j in n2e[a] if fq <= set(conn[j].tolist())][0]; c1 = conn[E1]; pos = {int(n): i for i, n in enumerate(c1)}
        opp = {}
        for n in fq:
            i = pos[n]; nb = [int(c1[m if l == i else l]) for l, m in HEX_EDGES if i in (l, m)]; opp[n] = [x for x in nb if x not in fq][0]
        E2 = [j for j in set(n2e[opp[q0[0]]]) if j != E1 and set(opp.values()) <= set(conn[j].tolist())]; E2 = E2[0] if E2 else None
        opp2 = {}
        if E2 is not None:
            c2 = conn[E2]; pos2 = {int(n): i for i, n in enumerate(c2)}
            for n in fq:
                i = pos2[opp[n]]; nb = [int(c2[m if l == i else l]) for l, m in HEX_EDGES if i in (l, m)]; opp2[n] = [x for x in nb if x not in set(opp.values())][0]
        rows.append(dict(q=q0, E1=E1, E2=E2, opp=opp, opp2=opp2, a=a, b=b))
    F1 = np.zeros(3); F2 = np.zeros(3); F3 = np.zeros(3); F4 = np.zeros(3); F1trap = np.zeros(3); A_tot = 0.0
    xs, p0, p1, p2 = [], [], [], []; pw3 = []
    for r in rows:
        x4 = Xd[r["q"]]; ef4 = ef[r["q"]]; p4 = -K * ef4; p14 = -K * ef[[r["opp"][n] for n in r["q"]]]; area_vec = np.zeros(3)
        pm4 = 0.5 * (p4 + p14)
        for (rr, ss) in GP:
            N = shp(rr, ss); g = gt(x4, rr, ss); area_vec += g
            F1 += -(N @ p4) * g; F4 += -(N @ pm4) * g
        n_c = gt(x4, 0, 0); n_c = n_c / np.linalg.norm(n_c); A = np.linalg.norm(area_vec); A_tot += A
        F2 += -Pe[r["E1"]] * n_c * A
        # F3: wall extrapolation from element-pressure rows 1,2 (centroid distances along the facet normal)
        c1 = Xd[conn[r["E1"]]].mean(0); mid = x4.mean(0); s1 = (c1 - mid) @ n_c
        if r["E2"] is not None:
            c2 = Xd[conn[r["E2"]]].mean(0); s2 = (c2 - mid) @ n_c; pw = Pe[r["E1"]] + (Pe[r["E1"]] - Pe[r["E2"]]) * s1 / (s2 - s1)
        else: pw = Pe[r["E1"]]
        F3 += -pw * n_c * A; pw3.append(pw)
        # trapezoid check of F1 (linear in x along facet, z-uniform)
        F1trap += -(p4.mean()) * area_vec
        # nodal rows for the parity decomposition (z0 nodes of the facet)
        z0n = [n for n in r["q"] if n in z0]; a_, b_ = r["a"], r["b"]
        r["p"] = {"x": [X0[a_, 0], X0[b_, 0]], "p0": [-K * ef[a_], -K * ef[b_]], "p1": [-K * ef[r["opp"][a_]], -K * ef[r["opp"][b_]]],
                  "p2": [-K * ef[r["opp2"][a_]] if r["opp2"] else np.nan, -K * ef[r["opp2"][b_]] if r["opp2"] else np.nan]}
    # nodal rows ordered by x
    node = {}
    for r in rows:
        for t, (n) in enumerate((r["a"], r["b"])): node[n] = (X0[n, 0], -K * ef[n], -K * ef[r["opp"][n if True else 0]] if False else None)
    # build row arrays by node
    nodes = sorted({r["a"] for r in rows} | {r["b"] for r in rows}, key=lambda n: X0[n, 0]); idx = {n: i for i, n in enumerate(nodes)}
    P0 = np.array([-K * ef[n] for n in nodes]); P1 = np.full(len(nodes), np.nan); P2 = P1.copy(); XN = np.array([X0[n, 0] for n in nodes]); YD = np.array([Xd[n, 1] for n in nodes])
    for r in rows:
        for n in (r["a"], r["b"]):
            zz = [m for m in r["q"] if abs(X0[m, 2]) < 1e-9 and abs(X0[m, 0] - X0[n, 0]) < 1e-9][0]
            P1[idx[n]] = -K * ef[r["opp"][zz]]
            if r["opp2"]: P2[idx[n]] = -K * ef[r["opp2"][zz]]
    # decomposition at interior nodes (rows 0,1,2): smooth s=(p0+2p1+p2)/4, parity a=(p0-2p1+p2)/4  (p0 = s + a exactly for p_yy=0)
    ok = np.isfinite(P2); ok[:4] = False; ok[-4:] = False
    S_ = (P0 + 2 * P1 + P2) / 4; A_ = (P0 - 2 * P1 + P2) / 4
    # x-parity along the wall row: a_x = (-1)^i (p_{i-1} - 2 p_i + p_{i+1})/4 at interior nodes
    d2 = np.zeros_like(P0); d2[1:-1] = (P0[:-2] - 2 * P0[1:-1] + P0[2:]) / 4; sgn = (-1.0) ** np.arange(len(P0)); axv = sgn * d2
    # forces of the smooth and parity parts (flat-floor approximation: F = -sum p_wall * (area vector) with the nodal weights of the F1 quadrature)
    # use the F1 operator linearly: F1[p0] for p0 -> replaced by (S_) / (A_) on interior nodes, others unchanged
    def F1_of(pn):
        F = np.zeros(3)
        for r in rows:
            x4 = Xd[r["q"]]; p4 = np.array([pn[idx[[r["a"], r["b"]][0 if abs(X0[m, 0] - X0[r["a"], 0]) < 1e-9 else 1]]] for m in r["q"]])
            for (rr, ss) in GP: F += -(shp(rr, ss) @ p4) * gt(x4, rr, ss)
        return F
    Fy_s = F1_of(np.where(ok | True, np.where(np.isfinite(S_), S_, P0), P0))   # smooth part where defined, raw at the 4 end nodes each side
    out = dict(tag=tag, time=S[k]["time"], flag=S[k]["flag"], TZ=TZ, K=K, n_facets=len(rows), area_per_depth=A_tot / TZ, previous_complete_time=S[comp[-2]]["time"] if len(comp) > 1 else None,
               F1_pressure_Fx_Fy_per_m=(F1[:2] / TZ).tolist(), F1_trapezoid_check=(F1trap[:2] / TZ).tolist(), F2_elementpressure_Fx_Fy_per_m=(F2[:2] / TZ).tolist(),
               F3_wall_extrapolated_element_pressure_Fx_Fy_per_m=(F3[:2] / TZ).tolist(), F4_pairrow_diagnostic_Fx_Fy_per_m=(F4[:2] / TZ).tolist(),
               F1_minus_F2_Fy=float(F1[1] / TZ - F2[1] / TZ), F1_minus_F4_Fy=float(F1[1] / TZ - F4[1] / TZ), F2_minus_F4_Fy=float(F2[1] / TZ - F4[1] / TZ),
               F_smooth_part_1_2_1_nodal_Fy=float(Fy_s[1] / TZ), F_parity_part_Fy=float((F1[1] - Fy_s[1]) / TZ),
               D_row1_minus_row0_mean_interior=float(np.mean((P1 - P0)[ok])), parity_amp_a_rms_interior=float(np.sqrt(np.mean(A_[ok] ** 2))), parity_amp_a_mean_interior=float(np.mean(A_[ok])),
               smooth_s_mean_interior=float(np.mean(S_[ok])), p0_mean_interior=float(np.mean(P0[ok])), x_parity_amp_rms_interior_nodes=float(np.sqrt(np.mean(axv[3:-3] ** 2))),
               x_parity_signed_mean=float(np.mean(axv[3:-3])), p_row0_x_smooth_rms_of_second_difference=float(np.sqrt(np.mean(d2[3:-3] ** 2))),
               ef_zlayer_consistency_max_abs_Pa=float(max(abs(K * (ef[r["q"][0]] - ef[r["q"][1]])) for r in rows)),
               nodal_viscous_not_reconstructable="F1 viscous term needs Gauss-point velocity gradients; xplt holds element-average velocity only (see audit 2-pt estimate for a bound)")
    dump = dict(tag=np.array(tag), time=np.array(S[k]["time"]), flag=np.array(S[k]["flag"]), K=np.array(K), TZ=np.array(TZ),
        XN=XN, P0=P0, P1=P1, P2=P2, YD=YD,
        q_cols=np.array([[idx[[r["a"], r["b"]][0 if abs(X0[m, 0] - X0[r["a"], 0]) < 1e-9 else 1]] for m in r["q"]] for r in rows]),
        X4=np.array([Xd[r["q"]] for r in rows]), p4=np.array([-K * ef[r["q"]] for r in rows]), p14=np.array([-K * ef[[r["opp"][n] for n in r["q"]]] for r in rows]),
        Pe1=np.array([Pe[r["E1"]] for r in rows]), Pe2=np.array([Pe[r["E2"]] if r["E2"] is not None else np.nan for r in rows]),
        C1=np.array([Xd[conn[r["E1"]]].mean(0) for r in rows]), C2=np.array([Xd[conn[r["E2"]]].mean(0) if r["E2"] is not None else [np.nan] * 3 for r in rows]))
    out["_dump"] = dump
    out["parity_force_estimate_from_D"] = float(-(out["D_row1_minus_row0_mean_interior"] / 2) * out["area_per_depth"] * -1)  # +D/2*L (sign: raising the wall-row pressure by D/2 adds -D/2*L to Fy)
    out["parity_force_estimate_from_D_note"] = "F2-F1 expected ~ -(D/2)*L (Fy more negative when the wall row is raised by D/2)"
    return out
if __name__ == "__main__":
    res = {}; fdir = sys.argv[2] if len(sys.argv) > 2 else None
    for key, (tag, tsel) in {"M1": ("C2_fsi_pilot_f1", 900.0), "M2": ("C3_fsi_M2", 750.0), "M3": ("C3_fsi_M4", None), "M3_t717p5_failure_exit_state": ("C3_fsi_M4", 717.5)}.items():
        res[key] = run(tag, f"{RUNS}/{tag}", tsel); dmp = res[key].pop("_dump")
        if fdir: os.makedirs(fdir, exist_ok=True); np.savez_compressed(f"{fdir}/facets_{key}.npz", **dmp)
        print(key, res[key]["time"], {k: res[key][k] for k in ("F1_pressure_Fx_Fy_per_m", "F2_elementpressure_Fx_Fy_per_m", "F3_wall_extrapolated_element_pressure_Fx_Fy_per_m", "F4_pairrow_diagnostic_Fx_Fy_per_m", "F_parity_part_Fy", "D_row1_minus_row0_mean_interior")})
    json.dump(res, open(sys.argv[1], "w"), indent=1)
