#!/usr/bin/env python3
"""Independent re-computation (own code; does NOT import or call the original analysis scripts) of the key numbers of the collapsibleChannel evidence from the stored raw outputs.
usage: reverify.py PACKAGE_ROOT GEN_DIR OUT.json
GEN_DIR holds the regenerated inputs (gen/<case>/<case>.feb, mesh2d.json) made by scripts/make_cc_feb.py.  Prints a table 'quantity | recomputed | stored/claimed | agrees'."""
import sys, os, re, gzip, json, numpy as np
ROOT, GEN, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
RAW = f"{ROOT}/raw"; REF = f"{ROOT}/reference"
out = {}; table = []
def chk(name, val, claim, tol, reported=None, inc=None):
    """claim = the value the recomputation is compared with; reported = the (different) number printed in the pilot report when a documented inconsistency exists (INCONSISTENCIES.md id inc)."""
    ok = None if claim is None else bool(abs(val - claim) <= tol)
    table.append((name + (f"  [report said {reported:g}; see {inc}]" if reported is not None else ""), val, claim, ok, tol)); out[name] = dict(value=val, claimed=claim, tol=tol, agrees=ok, reported_in_pilot_results=reported, inconsistency=inc)
def mon(path):
    d = {}; t = None
    for l in gzip.open(path, "rt"):
        if l.startswith("*Time"): t = float(l.split("=")[1]); d[t] = {}
        elif not l.startswith("*") and l.strip(): f = l.split(); d[t][int(f[0])] = np.array(f[1:], float)
    return d
def count(feb, name):
    txt = open(feb, encoding="latin-1").read(); m = re.search(r'<Elements type="hex8" name="%s">(.*?)</Elements>' % name, txt, re.S); return len(re.findall(r"<elem ", m.group(1)))
# ---- 1 mesh counts
chk("tutorial-config fluid cells (40+80+80 columns x 16 rows, generator defaults)", count(f"{GEN}/tut_fluid/tut_fluid.feb", "Fluid"), 3200, 0)
chk("tutorial-config solid cells (80 x 8)", count(f"{ROOT}/inputs/S1_80x8.feb", "Wall"), 640, 0)
chk("pilot fluid cells (40+320+80 columns x 16 rows)", count(f"{GEN}/D1_fsi_20steps/D1_fsi_20steps.feb", "Fluid"), 7040, 0)
chk("pilot wall cells (320 x 8)", count(f"{GEN}/D1_fsi_20steps/D1_fsi_20steps.feb", "Wall"), 2560, 0)
chk("pilot total cells", count(f"{GEN}/D1_fsi_20steps/D1_fsi_20steps.feb", "Fluid") + count(f"{GEN}/D1_fsi_20steps/D1_fsi_20steps.feb", "Wall"), 9600, 0)
chk("160x4 solid cells", count(f"{ROOT}/inputs/S1_160x4.feb", "Wall"), 640, 0)
# ---- 2 CC-S1
REFM, REFQ = -0.14362832, -0.09645925
S1 = {}
for tag, claim, cq in (("80x8", -0.11633, -0.19461), ("160x4", -0.030538, -0.062667), ("160x8", -0.029996, -0.061661), ("320x8", -0.0075661, -0.017086)):
    d = mon(f"{RAW}/S1_{tag}/Mon_bottom.txt.gz"); tl = max(d); A = np.array([d[tl][k] for k in sorted(d[tl])])
    x, uy = A[:, 0], A[:, 3]; mid = uy[np.argmin(abs(x - 10.0))]; q1 = uy[np.argmin(abs(x - 7.5))]; S1[tag] = mid
    chk(f"S1 {tag} load level reached (final step time)", tl, 1.0, 1e-12); chk(f"S1 {tag} midpoint error vs -0.14362832", mid / REFM - 1, claim, 2e-5); chk(f"S1 {tag} quarter-point error vs -0.09645925", q1 / REFQ - 1, cq, 2e-5)
e = {k: abs(v / REFM - 1) for k, v in S1.items()}
chk("S1 error ratio 80x8 -> 160x4/8 (160x8 used: same rows)", e["80x8"] / e["160x8"], 3.9, 0.06)
chk("S1 error ratio 160x8 -> 320x8", e["160x8"] / e["320x8"], 4.0, 0.06)
m160, m320 = S1["160x8"], S1["320x8"]; rich = m320 + (m320 - m160) / 3.0
chk("S1 Richardson (p = 2) limit, 160x8 + 320x8 (same thickness rows), relative to the beam value", rich / REFM - 1, -8.94e-5, 2e-6, reported=1e-4, inc="INC-3")
rich2 = S1["320x8"] + (S1["320x8"] - S1["160x4"]) / 3.0; chk("S1 Richardson (p = 2) limit, 160x4 + 320x8 (rows 4 vs 8, the pair used in the report), relative to the beam value", rich2 / REFM - 1, 9.11e-5, 1e-6)
chk("S1 160x4 gate margin: error minus 3 % (percentage points)", 100 * (e["160x4"] - 0.03), 0.054, 0.001)
# ---- 3 CC-F1 (own implementation)
M = json.load(open(f"{GEN}/F1_unpinned/mesh2d.json")); K = 1e7; nodes = np.array(M["nodes2"]); nr = M["nr"]; i0, i1 = M["i0w"], M["i1w"]; xs = np.array(M["xs"]); NX = len(xs) - 1
F = lambda i, j: i * (nr + 1) + j
rows = mon(f"{RAW}/F1_unpinned/Mon_rows.txt.gz"); tl = max(rows)
P = lambda i, j, d=rows, t=None: -K * d[tl if t is None else t][F(i, j) + 1][2]
p0 = np.array([P(i, nr) for i in range(i0, i1 + 1)]); p1 = np.array([P(i, nr - 1) for i in range(i0, i1 + 1)]); p2 = np.array([P(i, nr - 2) for i in range(i0, i1 + 1)]); p3 = None
Xw = nodes[[F(i, nr) for i in range(i0, i1 + 1)]]; dxf, dyf = np.diff(Xw[:, 0]), np.diff(Xw[:, 1]); xm = 0.5 * (Xw[:-1, 0] + Xw[1:, 0]); xref = 0.5 * (xs[i0:i1][:] + xs[i0 + 1:i1 + 1])
nodal = 0.5 * (p0[:-1] + p0[1:]); mean = 0.25 * (p0[:-1] + p0[1:] + p1[:-1] + p1[1:])
F1y, F2y = np.sum(nodal * dxf), np.sum(mean * dxf); phi = 0.5 * (1 - np.cos(2 * np.pi * (xref - 5.0) / 10.0)); G1, G2 = np.sum(nodal * phi * dxf), np.sum(mean * phi * dxf)
chk("F1 net F_y (nodal wall pressure) [N/m]", F1y, 31.6837, 5e-4); chk("F2 net F_y (element mean) [N/m]", F2y, 34.0665, 5e-4)
chk("F1 vs F2 relative, F_y", F1y / F2y - 1, -0.06995, 1e-4); chk("F1 vs F2 relative, first-mode load G", G1 / G2 - 1, -0.08629, 1e-4)
def flux(d, col):
    ids = [F(col, j) + 1 for j in range(nr + 1)]; return np.trapz([d[tl][i][0] for i in ids], [nodes[F(col, j)][1] for j in range(nr + 1)])
qi, qo = flux(rows, 0) if False else flux(mon(f"{RAW}/F1_unpinned/Mon_inlet.txt.gz"), 0), flux(mon(f"{RAW}/F1_unpinned/Mon_outlet.txt.gz"), NX)
chk("F1 mass balance Q_out/Q_in - 1", qo / qi - 1, 6.0e-7, 2e-8); chk("F1 Q_in [m^2/s per unit depth]", qi, 0.8810, 5e-4)
# steady drift of F1, F2 over the last 20 % of steps
ts = sorted(t for t in rows if t > 0); t80 = ts[int(round(0.8 * len(ts))) - 1]
def force_at(t):
    a = np.array([P(i, nr, t=t) for i in range(i0, i1 + 1)]); b = np.array([P(i, nr - 1, t=t) for i in range(i0, i1 + 1)])
    return np.sum(0.5 * (a[:-1] + a[1:]) * dxf), np.sum(0.25 * (a[:-1] + a[1:] + b[:-1] + b[1:]) * dxf)
f80, fl = force_at(t80), force_at(ts[-1]); chk("F1 steady drift last 20 % of steps (relative)", abs(fl[0] - f80[0]) / abs(fl[0]), 1.2e-11, 1e-9)
# streamwise parity, wall row: rms over interior nodes with 3 facets excluded
def at(p, skip=3): k = np.arange(1 + skip, len(p) - 1 - skip); return np.sqrt(np.mean(((p[k] - 0.5 * (p[k - 1] + p[k + 1])) / 2) ** 2))
inner = slice(4, len(p0) - 4); V = p0[inner].max() - p0[inner].min()
chk("F1 streamwise parity a_t (wall row) [Pa]", at(p0), 2.4e-4, 2e-5); chk("F1 streamwise parity a_t / V", at(p0) / V, 6.8e-5, 5e-6); chk("F1 smooth variation V along the wall [Pa]", V, 3.51, 0.01)
# wall-normal row parity
chk("F1 max |p(wall) - p(row 2)| (all wall nodes) [Pa]", np.max(np.abs(p0 - p2)), 4.09e-3, 1e-5, reported=1e-3, inc="INC-4"); chk("F1 rms |p(wall) - p(row 2)| [Pa]", float(np.sqrt(np.mean((p0 - p2) ** 2))), 2.38e-3, 1e-5); chk("F1 rms |p(wall) - p(row 1)| [Pa]", float(np.sqrt(np.mean((p0 - p1) ** 2))), None, 0)
d10 = p0 - p1; chk("F1 p(wall) - p(row1): minimum [Pa]", d10.min(), -1.338, 0.002); chk("F1 p(wall) - p(row1): maximum [Pa]", d10.max(), 0.551, 0.002, reported=0.47, inc="INC-4")
D = (p0 - 2 * p1 + p2)[inner]; chk("F1 D = p0 - 2 p1 + p2: minimum [Pa]", D.min(), -2.7, 0.05); chk("F1 D: maximum [Pa]", D.max(), 1.09, 0.02); chk("F1 D: rms [Pa]", np.sqrt(np.mean(D ** 2)), 1.75, 0.01)
# ---- 4 D1
M = json.load(open(f"{GEN}/D1_fsi_20steps/mesh2d.json")); nodes = np.array(M["nodes2"]); iface = M["iface"]; xs = np.array(M["xs"]); nw = i1 - i0
dr = mon(f"{RAW}/D1_fsi_20steps/Mon_rows.txt.gz"); di = mon(f"{RAW}/D1_fsi_20steps/Mon_iface.txt.gz"); T = np.array(sorted(t for t in di if t > 0))
U = np.array([[di[t][iface[k] + 1][1] for k in range(nw + 1)] for t in T]); UX = np.array([[di[t][iface[k] + 1][0] for k in range(nw + 1)] for t in T]); mid = U[:, nw // 2]
xi = xs[i0:i1 + 1]; xm0 = 0.5 * (xi[:-1] + xi[1:]); phi0 = 0.5 * (1 - np.cos(2 * np.pi * (xm0 - 5.0) / 10.0))
sp = {}
for j, t in enumerate(T):
    a = np.array([-K * dr[t][F(i, nr) + 1][2] for i in range(i0, i1 + 1)]); b = np.array([-K * dr[t][F(i, nr - 1) + 1][2] for i in range(i0, i1 + 1)]); c = np.array([-K * dr[t][F(i, nr - 2) + 1][2] for i in range(i0, i1 + 1)])
    X = xi + UX[j]; dX = np.diff(X); n_ = 0.5 * (a[:-1] + a[1:]); m_ = 0.25 * (a[:-1] + a[1:] + b[:-1] + b[1:])
    sp[round(t, 4)] = dict(F1=np.sum(n_ * dX), F2=np.sum(m_ * dX), G1=np.sum(n_ * phi0 * dX), G2=np.sum(m_ * phi0 * dX), D=np.mean((a - 2 * b + c)[4:-4]))
for t, f1, f2, sf, sg, dd in ((0.10, 791, 663, 0.194, 0.194, 51), (0.15, 1342, 1078, 0.245, 0.246, 106), (0.20, 1702, 1301, 0.308, 0.312, 161), (0.25, 1659, 1172, 0.416, 0.424, 196), (0.30, 1313, 808, 0.624, 0.651, 203), (0.40, -717, -865, -0.170, -0.186, 60), (0.45, -2093, -1843, 0.136, 0.144, -101), (0.50, -2842, -2234, 0.272, 0.293, -245)):
    r = sp[round(t, 4)]; chk(f"D1 t={t:.2f} F1_y [N/m]", r["F1"], f1, 1.0); chk(f"D1 t={t:.2f} F2_y [N/m]", r["F2"], f2, 1.0); chk(f"D1 t={t:.2f} (F1-F2)/F2", r["F1"] / r["F2"] - 1, sf, 1.5e-3)
    chk(f"D1 t={t:.2f} G spread", r["G1"] / r["G2"] - 1, sg, 1.5e-3); chk(f"D1 t={t:.2f} mean y-parity D [Pa]", r["D"], dd, 1.0)
t = 0.25; ids = [F(M["i0w"] + 0, 0)]; ix = int(np.argmin(abs(xs - 10.0)))
vals = [-K * dr[0.25][F(ix, nr - r) + 1][2] for r in range(3)] + [-K * dr[0.25][F(ix, 0) + 1][2]]
for nm, v, c in zip(("wall", "row 1", "row 2", "channel floor"), vals, (193.4, 77.0, 192.9, 191.7)): chk(f"D1 pressure at x=10, t=0.25, {nm} [Pa]", v, c, 0.06)
lg = gzip.open(f"{RAW}/D1_fsi_20steps/D1_fsi_20steps.log.gz", "rt").read(); its = [int(m.group(1)) for m in re.finditer(r"number of iterations\s*:\s*(\d+)", lg)]
chk("D1 Newton iterations per step: mean", float(np.mean(its)), 5.25, 0.005); chk("D1 Newton iterations per step: max", max(its), 6, 0); chk("D1 steps completed", len(its), 20, 0)
chk("D1 normal termination (1 = yes)", float("N O R M A L   T E R M I N A T I O N" in lg), 1.0, 0); chk("D1 'negative jacobian' occurrences in log (0 = none)", float(lg.lower().count("negative jacobian")), 0.0, 0)
# volume conservation
mi_, mo_ = mon(f"{RAW}/D1_fsi_20steps/Mon_inlet.txt.gz"), mon(f"{RAW}/D1_fsi_20steps/Mon_outlet.txt.gz")
def fx(d, col, t): return np.trapz([d[t][F(col, j) + 1][0] for j in range(nr + 1)], [nodes[F(col, j)][1] for j in range(nr + 1)])
dq = np.array([0.0] + [fx(mi_, 0, t) - fx(mo_, NX, t) for t in T]); tt = np.concatenate([[0.0], T]); cum = np.concatenate([[0.0], np.cumsum(0.5 * (dq[1:] + dq[:-1]) * np.diff(tt))])
vol = np.concatenate([[0.0], [np.trapz(U[j], xi) for j in range(len(T))]]); chk("D1 volume conservation residual / peak area change", np.max(np.abs(cum - vol)) / np.max(np.abs(vol)), 0.0024, 3e-4)
# ---- 5 reference CSVs (own extraction): 5-point parabolic fit around each extremum
def load(n):
    r = [l.split(",") for l in open(f"{REF}/collapsibleChannel_oomph_{n}.csv") if l[0] not in "#t"]; d = np.array(r, float); return d[:, 0], d[:, 2]
def ext(t, u, kind, lo, hi):
    m = (t >= lo) & (t <= hi); k = int(np.argmin(np.where(m, u, 9))) if kind == "min" else int(np.argmax(np.where(m, u, -9)))
    s = slice(max(k - 2, 0), k + 3); c = np.polyfit(t[s], u[s], 2); tt_ = -c[1] / (2 * c[0]); return tt_, np.polyval(c, tt_)
tr, ur = load("reference")
t1, u1 = ext(tr, ur, "min", 0.3, 0.75); tM, uM = ext(tr, ur, "max", 0.6, 1.0); t2, u2 = ext(tr, ur, "min", 0.9, 1.4); t3, u3 = ext(tr, ur, "min", 1.5, 2.0)
chk("oomph converged first trough value [m]", u1, -0.2160, 5e-4); chk("oomph converged first trough time [s]", t1, 0.4685, 1.5e-3)
chk("oomph first rebound maximum value [m]", uM, -0.0404, 5e-4); chk("oomph first rebound maximum time [s]", tM, 0.791, 3e-3)
chk("oomph second trough value [m]", u2, -0.2037, 5e-4); chk("oomph second trough time [s]", t2, 1.119, 2e-3); chk("oomph third trough time [s]", t3, 1.753, 3e-3)
chk("oomph trough-to-trough period T12 [s]", t2 - t1, 0.650, 3e-3); chk("oomph trough-to-trough period T23 [s]", t3 - t2, 0.634, 4e-3); chk("oomph trough to rebound-maximum spacing [s] (README '0.8 s')", tM - t1, None, 0)
t40, u40 = load("dt40"); a40, b40 = ext(t40, u40, "min", 0.3, 0.75); chk("oomph dt=0.025 first trough value [m]", b40, -0.2129, 5e-4); chk("oomph dt=0.025 first trough time [s]", a40, 0.4735, 3e-3)
# D1 lag vs converged oomph (interp), in % of 0.216
for tq, fe, oo in ((0.25, -0.0523, -0.0560), (0.30, -0.0874, -0.0979), (0.35, -0.1248, -0.1451), (0.425, -0.1764, -0.2041), (0.475, -0.1975, -0.2157), (0.50, -0.2016, -0.2101)):
    jf = int(np.argmin(abs(T - tq))); chk(f"D1 wall midpoint t={tq} [m]", mid[jf], fe, 1.5e-4); chk(f"oomph converged midpoint t={tq} [m]", float(np.interp(tq, tr, ur)), oo, 1.5e-4)
lag = np.array([mid[j] - np.interp(T[j], tr, ur) for j in range(len(T))]) / 0.216; chk("D1 max |history difference| vs converged oomph, t <= 0.5 s (fraction of 0.216)", float(np.max(np.abs(lag))), 0.128, 1.5e-3); chk("D1 time of max difference [s]", float(T[np.argmax(np.abs(lag))]), 0.425, 1e-9)
chk("D1 rms history difference (fraction of 0.216)", float(np.sqrt(np.mean(lag ** 2))), 0.065, 1.5e-3)
chk("D1 first trough bracketed? (min of last step is the minimum: 1 = still descending)", float(np.argmin(mid) == len(mid) - 1), 1.0, 0)
# ---- 6 free-body identification (valid while |u| <= |u_mid(S1_320x8, full load)|)
d = mon(f"{RAW}/S1_320x8/Mon_bottom.txt.gz"); ts_ = sorted(d); ids = sorted(d[ts_[0]]); um = np.array([d[t_][ids[len(ids) // 2]][3] for t_ in ts_]); pe = np.array([200.0 * t_ for t_ in ts_])
fin = np.array([d[ts_[-1]][i] for i in ids]); psi = fin[:, 3] / fin[len(ids) // 2, 3]; psim = 0.5 * (psi[:-1] + psi[1:]); dx0 = np.diff(xi); Ipsi = np.sum(psim * dx0); cin = np.sum(psim ** 2 * dx0) / Ipsi
pel = lambda u: np.interp(-u, -um, pe); pext = lambda t: 200.0 * (0.5 * (1 - np.cos(np.pi * t / 0.25)) if t < 0.25 else 1.0); acc = np.gradient(np.gradient(mid, T), T)
res = {}
for tq, c_id, c_f1, c_f2 in ((0.20, 173.7, 182.2, 139.0), (0.225, 183.3, 184.0, 134.7), (0.25, 179.8, 178.5, 125.5), (0.30, 141.3, 141.2, 85.8), (0.35, 58.9, 59.4, 14.0)):
    j = int(np.argmin(abs(T - tq))); pid = pel(mid[j]) + 0.05 * cin * acc[j] - pext(T[j]); t = T[j]
    a = np.array([-K * dr[t][F(i, nr) + 1][2] for i in range(i0, i1 + 1)]); b = np.array([-K * dr[t][F(i, nr - 1) + 1][2] for i in range(i0, i1 + 1)])
    f1e = np.sum(0.5 * (a[:-1] + a[1:]) * psim * dx0) / Ipsi; f2e = np.sum(0.25 * (a[:-1] + a[1:] + b[:-1] + b[1:]) * psim * dx0) / Ipsi
    chk(f"free-body t={tq} identified fluid pressure, stored sign (p_el + inertia - p_ext) [Pa]", pid, -c_id, 0.25, reported=c_id, inc="INC-5"); chk(f"free-body t={tq} F1 equivalent [Pa]", f1e, c_f1, 0.25); chk(f"free-body t={tq} F2 equivalent [Pa]", f2e, c_f2, 0.25)
    res[tq] = (pid, f1e, f2e)
chk("free-body validity limit |u_mid| of S1_320x8 at full 200 Pa [m]", -um[-1], 0.1425, 1e-4)
# ---- 7 sensitivity-law arithmetic (analytical estimate of BENCHMARK_REASSESSMENT_AFTER_FSI1.md section 4)
E_eff = 4.55e8 / (1 - 0.3 ** 2); h, L, u0 = 0.05, 10.0, 0.2160; EI = E_eff * h ** 3 / 12
chk("analytic: E_eff [Pa]", E_eff, 5.0e8, 1e5); chk("analytic: EI [N m]", EI, 5208.0, 1.0)
bend = 384 * EI / L ** 4 * u0; T_ = E_eff * h * (8 / 3) * (u0 / L) ** 2; mem = 8 * T_ * u0 / L ** 2
chk("analytic: bending restoring pressure at first trough [Pa]", bend, 43.0, 0.5); chk("analytic: membrane restoring pressure at first trough [Pa]", mem, 537.0, 3.0)
ub = 0.14363; chk("analytic: static check bending+membrane at 0.14363 m [Pa] (vs applied 200)", 384 * EI / L ** 4 * ub + 8 * E_eff * h * (8 / 3) * (ub / L) ** 2 * ub / L ** 2, 186.7, 0.5)
chk("analytic: linear bending-only static midpoint 200 L^4/(384 EI) [m]", 200 * L ** 4 / (384 * EI), 1.0, 0.01)
json.dump(out, open(OUT, "w"), indent=1, default=float)
w = max(len(r[0]) for r in table)
for name, v, c, ok, tol in table: print(f"{name:{w}s} | {v: .6g} | {'' if c is None else format(c, '.6g'):>10s} | {'-' if ok is None else ('yes' if ok else 'NO')}")
bad = [r for r in table if r[3] is False]; print(f"\n{len(table)} quantities, {len(bad)} outside tolerance"); sys.exit(1 if bad else 0)
