"""Extract the input CSVs of this package from the raw run output of the collapsibleChannel campaign.
usage:  python3 -B -I scripts/extract_from_raw.py EVIDENCE_DIR [OUT_DIR]
EVIDENCE_DIR is the directory holding P/, C1/, C2/, C3/, MAIN/, EXT/ (the extracted evidence tarball, see README.md); OUT_DIR defaults to ../data.
Needs numpy and scipy.  Not needed to regenerate the derived tables (scripts/analyze_collapsible.py) or the figures.
Definitions (fixed by the preregistration, PREREG_COLLAPSIBLE_PRODUCTION.md section 2): P = 0.216 m, '% P' = percent of P; x = 7.5, 10, 12.5 (25, 50, 75 % of the wall);
V1 = first trough, V2 = first maximum after V1, V3 = next minimum, found on a cubic Hermite interpolant of the stored wall displacement and interface velocity (output every step);
bottom face y = 1, midplane y = 1.025.  Fixed times 0.25, 0.46875, 0.79, 1.11875, 1.5 s."""
import sys, os, json, csv, glob
import numpy as np
from scipy.interpolate import CubicHermiteSpline, CubicSpline

EV = sys.argv[1]
OUT = sys.argv[2] if len(sys.argv) > 2 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
DATA = f"{EV}/P/benchmarks/collapsibleChannel/colP_data"
REF = f"{EV}/EXT/oomph_reference"
S4F = f"{EV}/EXT/s4f_mesh_study"
MS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "collapsible_mesh_study")   # oomph refinement files of PR #74 (in the repository)
P = 0.216
X = (7.5, 10.0, 12.5)
FIXED_T = (0.25, 0.46875, 0.79, 1.11875, 1.5)
NY = {"S1": 8, "S2": 16, "S3": 32}


def fl(v, n=10): return "" if v is None else ("%.*g" % (n, v))


def wr(name, header, rows):
    with open(os.path.join(OUT, name), "w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n"); w.writerow(header)
        for r in rows: w.writerow(r)


def load(name): return json.load(open(f"{DATA}/data_{name}.json"))


def hist(d, key, ukey):
    t = np.array([0.0] + [s["t"] for s in d["steps"]])
    D = np.vstack([np.zeros(3)] + [s[key] for s in d["steps"]]); U = np.vstack([np.zeros(3)] + [s[ukey] for s in d["steps"]])
    return t, D, U


def extrema(f, t0, t1, n=200001):
    tt = np.linspace(t0, t1, n); v = f(tt); dv = np.gradient(v, tt); out = []
    r = f.derivative()
    for k in np.where(np.sign(dv[:-1]) * np.sign(dv[1:]) < 0)[0]:
        a, b = tt[k], tt[k + 1]; ga = float(r(a))
        for _ in range(60):
            m = 0.5 * (a + b); gm = float(r(m))
            if np.sign(gm) == np.sign(ga): a, ga = m, gm
            else: b = m
        tm = 0.5 * (a + b); out.append((tm, float(f(tm)), "min" if dv[k] < 0 else "max"))
    return out


def V123(f, tmax=1.5):
    """V1: first local minimum deeper than 0.35 of the global minimum; V2: first maximum after V1 rising > 0.02 m; V3: next minimum falling > 0.02 m below V2."""
    ex = extrema(f, 0.0, tmax)
    gmin = min(e[1] for e in ex if e[2] == "min")
    v1 = [e for e in ex if e[2] == "min" and e[1] < 0.35 * gmin][0]
    maxs = [e for e in ex if e[2] == "max" and e[0] > v1[0] and e[1] - v1[1] > 0.02]
    v2 = maxs[0]
    mins3 = [e for e in ex if e[2] == "min" and e[0] > v2[0] and v2[1] - e[1] > 0.02]
    return v1, v2, (mins3[0] if mins3 else None)


def oomph(tag):
    f = {"ref": "oomph_reference", "dt40": "oomph_dt40", "dt80": "oomph_dt80", "dt160": "oomph_dt160", "dt320": "oomph_dt320"}[tag]
    A = np.genfromtxt(f"{REF}/collapsibleChannel_{f}.csv", delimiter=",", comments="#"); A = A[np.isfinite(A[:, 0])]
    return [CubicSpline(A[:, 0], A[:, 1 + i]) for i in range(3)]


def oomph_file(name):
    A = np.genfromtxt(os.path.join(MS, name), delimiter=",", comments="#", skip_header=0); A = A[np.isfinite(A[:, 0])]
    return [CubicSpline(A[:, 0], A[:, 1 + i]) for i in range(3)]


def s4f(F):
    A = np.genfromtxt(f"{S4F}/{F}.dat", comments="#"); A = A[np.isfinite(A[:, 0])]
    return CubicSpline(np.r_[0.0, A[:, 0]], np.r_[0.0, A[:, 2]])   # column 3 = D_y of wallMid


class Cell:
    def __init__(self, name):
        self.name = name; d = self.d = load(name)
        self.t, D, U = hist(d, "wall_d", "wall_u"); _, Dm, Um = hist(d, "midplane_d", "midplane_u")
        self.tend = float(self.t[-1])
        self.bot = [CubicHermiteSpline(self.t, D[:, i], U[:, i]) for i in range(3)]
        self.mid = [CubicHermiteSpline(self.t, Dm[:, i], Um[:, i]) for i in range(3)]
        tm = min(1.5, self.tend)
        try: self.v = {tag: [V123(sp[i], tm) for i in range(3)] for tag, sp in (("bottom", self.bot), ("midplane", self.mid))}
        except (ValueError, IndexError): self.v = None   # the 10-step default-path check never reaches a trough
    def V(self, face, i, k): return self.v[face][i][k][1] / P * 100


# ------------------------------------------------------------------ cell inventory
sacct_rows = list(csv.DictReader(open(f"{EV}/EXT/sacct_collapsible.psv"), delimiter="|"))
sacct = {}
for r in sacct_rows: sacct.setdefault(r["JobName"], r)    # job names of the production cells are unique; repeated names (smoke, debug) are only used in the cost table
names = sorted(os.path.basename(f)[5:-5] for f in glob.glob(f"{DATA}/data_S[123]_*.json"))
production = [n for n in names if n.split("_")[0] in NY]
cells = {n: Cell(n) for n in production}


def parse(n):
    p = n.split("_"); S = p[0]; dt = float(cells[n].d["args"]["dt"]) if n in cells else float(load(n)["args"]["dt"])
    var = "baseline"; extra = [q for q in p[1:] if not q.startswith("dt") and q != "r2"]
    if extra: var = extra[0]
    return S, var, dt


rows = []
for n, c in cells.items():
    d = c.d; S, var, dt = parse(n); st = d["steps"]
    sj = sacct.get("colP_" + n)
    grp = "baseline" if var == "baseline" else ("check" if var in ("eq10", "rerun") else "stage3")
    rows.append([n, grp, S, NY[S], 10 * NY[S], var, fl(dt), fl(c.tend, 6), len(st), d["ndof"], d["ranks"], d["status"], d["all_converged"],
                 max(s["newton"] for s in st), sum(1 for s in st if s["retried"]), fl(min(s["jmin"] for s in st), 6), "%.3e" % max(abs(s["mass_cons"]) for s in st),
                 "%.3e" % max(abs(s["mass_res"]) for s in st), fl(d["cpu_s_per_step_mean"], 5), sj["JobID"] if sj else "", sj["CPUTimeRAW"] if sj else ""])
wr("cells.csv", ["cell", "group", "mesh", "ny_fluid", "nx_wall", "variant", "dt_s", "t_end_s", "steps", "ndof", "ranks", "status", "all_converged", "max_newton_per_step",
                 "newton_retries", "min_jacobian", "max_abs_q_in_out_wall", "max_abs_bdf_gcl_residual", "driver_cpu_s_per_step", "slurm_job", "cpu_s_sacct"], rows)

# ------------------------------------------------------------------ troughs
rows = []
for n, c in cells.items():
    S, var, dt = parse(n)
    if c.v is None: continue
    for face in ("bottom", "midplane"):
        for i, x in enumerate(X):
            v = c.v[face][i]
            rows.append([n, S, var, fl(dt), face, x, fl(v[0][1] / P * 100, 8), fl(v[0][0], 8), fl(v[1][1] / P * 100, 8), fl(v[2][1] / P * 100, 8) if v[2] else ""])
wr("troughs_fe.csv", ["cell", "mesh", "variant", "dt_s", "face", "x_m", "V1_pctP", "t_V1_s", "V2_pctP", "V3_pctP"], rows)

rows = []
for tag in ("ref", "dt40", "dt80", "dt160", "dt320"):
    sp = oomph(tag)
    for i, x in enumerate(X):
        v = V123(sp[i], 1.5); rows.append(["oomph_" + tag, x, fl(v[0][1] / P * 100, 8), fl(v[0][0], 8), fl(v[1][1] / P * 100, 8), fl(v[2][1] / P * 100, 8)])
for tag, fn in (("extrap_dt160", "oomph_extrapolated_dt160.csv"), ("refine1_dt160", "oomph_refine1_dt160.csv"), ("refine2_dt160", "oomph_refine2_dt160.csv"), ("refine3_dt160", "oomph_refine3_dt160.csv"), ("refine4_dt160", "oomph_refine4_dt160.csv")):
    sp = oomph_file(fn)
    for i, x in enumerate(X):
        v = V123(sp[i], 1.5); rows.append(["oomph_" + tag, x, fl(v[0][1] / P * 100, 8), fl(v[0][0], 8), fl(v[1][1] / P * 100, 8), fl(v[2][1] / P * 100, 8)])
for F in ("F1", "F2", "F3", "F4"):
    v = V123(s4f(F), 1.5); rows.append(["s4f_" + F, 10.0, fl(v[0][1] / P * 100, 8), fl(v[0][0], 8), fl(v[1][1] / P * 100, 8), fl(v[2][1] / P * 100, 8)])
wr("troughs_reference.csv", ["series", "x_m", "V1_pctP", "t_V1_s", "V2_pctP", "V3_pctP"], rows)

# ------------------------------------------------------------------ fixed-time differences and history norms against oomph / s4f
oref = oomph("ref"); DTREF = {0.025: "dt40", 0.0125: "dt80", 0.00625: "dt160", 0.003125: "dt320"}
rows = []; hrows = []
for n, c in cells.items():
    S, var, dt = parse(n)
    if var not in ("baseline",): continue
    for i, x in enumerate(X):
        row = [n, S, fl(dt), x]
        for tt in FIXED_T:
            row.append(fl((float(c.bot[i](tt)) - float(oref[i](tt))) / P * 100, 6) if tt <= c.tend + 1e-9 else "")
        rows.append(row)
    refsets = [("oomph_reference", oref), ("oomph_same_dt", oomph(DTREF[dt]))]
    wins = [1.25] + ([1.5] if c.tend >= 1.5 - 1e-9 else [])
    for rn, rs in refsets:
        for w in wins:
            tt = np.linspace(0.0, w, 6001)
            for i, x in enumerate(X):
                e = (np.array([float(c.bot[i](q)) for q in tt]) - rs[i](tt)) / P * 100
                hrows.append([n, S, fl(dt), rn, w, x, fl(float(np.sqrt(np.mean(e ** 2))), 6), fl(float(np.abs(e).max()), 6)])
    if dt == 0.00625:
        sp = s4f("F4")
        for w in wins:
            tt = np.linspace(0.0, w, 6001)
            for i, x in enumerate(X[1:2], start=1):
                e = (np.array([float(c.bot[i](q)) for q in tt]) - sp(tt)) / P * 100
                hrows.append([n, S, fl(dt), "s4f_F4", w, x, fl(float(np.sqrt(np.mean(e ** 2))), 6), fl(float(np.abs(e).max()), 6)])
wr("fixed_time_vs_oomph_reference.csv", ["cell", "mesh", "dt_s", "x_m"] + ["FE_minus_oomph_t%g_pctP" % t for t in FIXED_T], rows)
wr("history_norms.csv", ["cell", "mesh", "dt_s", "reference", "window_t_max_s", "x_m", "L2_pctP", "Linf_pctP"], hrows)

# error-series correlation and s4f-oomph norms at x = 10 (dt 0.00625)
rows = []
sp4 = s4f("F4"); o160 = oomph("dt160")
for n in ("S2_dt0.00625", "S3_dt0.00625"):
    c = cells[n]
    for w in (1.25, 1.5):
        tt = np.linspace(0.0, w, 6001)
        fe = np.array([float(c.bot[1](q)) for q in tt]); o = o160[1](tt); s = sp4(tt)
        e1 = (fe - o) / P * 100; e2 = (s - o) / P * 100
        rows.append([n, w, fl(float(np.corrcoef(e1, e2)[0, 1]), 5), fl(float(np.abs(e1).max()), 5), fl(float(np.abs(e2).max()), 5), fl(float(np.abs((fe - s) / P * 100).max()), 5)])
wr("error_series_correlation.csv", ["fe_cell", "window_t_max_s", "corr_FEminusOomph_vs_s4fminusOomph", "Linf_FE_minus_oomph_dt160_pctP", "Linf_s4fF4_minus_oomph_dt160_pctP", "Linf_FE_minus_s4fF4_pctP"], rows)

# comparison against the spatially extrapolated oomph-lib solution (Delta t = 0.00625 only; PR #74), bottom face and midplane, plus s4f F4 and the oomph refinement levels
ex = oomph_file("oomph_extrapolated_dt160.csv"); rows = []
OM_V1 = [V123(ex[i], 1.5)[0][1] / P * 100 for i in range(3)]
for n in ("S1_dt0.00625_r2", "S2_dt0.00625", "S3_dt0.00625"):
    c = cells[n]
    for face, sp in (("bottom", c.bot), ("midplane", c.mid)):
        for w in (1.25, 1.5):
            tt = np.linspace(0.0, w, 6001)
            for i, x in enumerate(X):
                e = (np.array([float(sp[i](q)) for q in tt]) - ex[i](tt)) / P * 100
                rows.append([n, face, w, x, fl(float(np.sqrt(np.mean(e ** 2))), 6), fl(float(np.abs(e).max()), 6)])
sp4 = s4f("F4")
for w in (1.25, 1.5):
    tt = np.linspace(0.0, w, 6001); e = (sp4(tt) - ex[1](tt)) / P * 100
    rows.append(["s4f_F4", "wallMid", w, 10.0, fl(float(np.sqrt(np.mean(e ** 2))), 6), fl(float(np.abs(e).max()), 6)])
for k in (1, 2, 3, 4):
    rs = oomph_file("oomph_refine%d_dt160.csv" % k)
    for i, x in enumerate(X):
        tt = np.linspace(0.0, 1.5, 6001); e = (rs[i](tt) - ex[i](tt)) / P * 100
        rows.append(["oomph_refine%d" % k, "beam", 1.5, x, fl(float(np.sqrt(np.mean(e ** 2))), 6), fl(float(np.abs(e).max()), 6)])
wr("history_norms_vs_extrapolated_oomph.csv", ["series", "face", "window_t_max_s", "x_m", "L2_pctP", "Linf_pctP"], rows)
rows = []
for n in ("S2_dt0.00625", "S3_dt0.00625"):
    c = cells[n]
    for w in (1.25, 1.5):
        tt = np.linspace(0.0, w, 6001); fe = np.array([float(c.bot[1](q)) for q in tt]); e1 = (fe - ex[1](tt)) / P * 100; e2 = (sp4(tt) - ex[1](tt)) / P * 100
        rows.append([n, w, fl(float(np.corrcoef(e1, e2)[0, 1]), 5), fl(float(np.abs(e1).max()), 5), fl(float(np.abs(e2).max()), 5), fl(float(np.abs((fe - sp4(tt)) / P * 100).max()), 5)])
wr("error_series_correlation_extrapolated_oomph.csv", ["fe_cell", "window_t_max_s", "corr_FEminusExtrap_vs_s4fF4minusExtrap", "Linf_FE_minus_extrap_pctP", "Linf_s4fF4_minus_extrap_pctP", "Linf_FE_minus_s4fF4_pctP"], rows)

# D1: order sensitivity of the extrapolation (refine 3 and 4, h ~ 1/N, ratio 4/3) relative to the first-order file; D5: correlation against other references
r3, r4 = oomph_file("oomph_refine3_dt160.csv"), oomph_file("oomph_refine4_dt160.csv"); rows = []
tt = np.linspace(0.0, 1.5, 6001)
for pe in (0.8, 1.0, 1.1, 1.25, 1.5):
    f = 1.0 / ((4.0 / 3.0) ** pe - 1.0)
    for i, x in enumerate(X):
        sp = lambda t: r4[i](t) + (r4[i](t) - r3[i](t)) * f
        v = V123(CubicSpline(tt, sp(tt)), 1.5)[0][1] / P * 100
        rows.append([pe, x, fl(v - OM_V1[i], 6), fl(float(np.abs((sp(tt) - ex[i](tt)) / P * 100).max()), 6)])
wr("extrapolation_order_sensitivity.csv", ["assumed_order_p", "x_m", "first_trough_change_vs_first_order_file_pctP", "history_Linf_change_vs_first_order_file_pctP"], rows)
rows = []
refs = [("refine1", oomph_file("oomph_refine1_dt160.csv")), ("refine2", oomph_file("oomph_refine2_dt160.csv")), ("refine3", r3), ("refine4", r4), ("extrapolated", ex)]
for n in ("S2_dt0.00625", "S3_dt0.00625"):
    c = cells[n]
    for rn, rs in refs:
        for w0, w1 in ((0.0, 1.25), (0.0, 1.5), (0.5, 1.25), (0.8, 1.25)):
            tt2 = np.linspace(w0, w1, 6001); fe = np.array([float(c.bot[1](q)) for q in tt2]); e1 = (fe - rs[1](tt2)) / P * 100; e2 = (sp4(tt2) - rs[1](tt2)) / P * 100
            rows.append([n, rn, w0, w1, fl(float(np.corrcoef(e1, e2)[0, 1]), 4)])
wr("error_series_correlation_by_reference.csv", ["fe_cell", "reference", "window_start_s", "window_end_s", "corr_FEminusRef_vs_s4fF4minusRef_x10"], rows)

# history-level cross term and space effects (bottom face).  cross term X(t) = [S3-S2](dt a) - [S3-S2](dt b) over t <= 1.25 (S3 at the finest dt ends at 1.25 s)
def series(name, i, w):
    tt = np.linspace(0, w, 6001); return np.array([float(cells[name].bot[i](q)) for q in tt]) / P * 100
rows = []
for a, b in (("0.0125", "0.00625"), ("0.00625", "0.003125")):
    for i, x in enumerate(X):
        Xh = (series("S3_dt" + a, i, 1.25) - series("S2_dt" + a, i, 1.25)) - (series("S3_dt" + b, i, 1.25) - series("S2_dt" + b, i, 1.25))
        sp_ = series("S3_dt" + b, i, 1.25) - series("S2_dt" + b, i, 1.25); tm_ = series("S2_dt" + b, i, 1.25) - series("S2_dt" + a, i, 1.25)
        rows.append([a, b, x, fl(float(np.abs(Xh).max()), 5), fl(float(np.abs(sp_).max()), 5), fl(float(np.abs(tm_).max()), 5)])
wr("history_cross_term_S2_S3.csv", ["dt_a_s", "dt_b_s", "x_m", "max_abs_cross_term_pctP", "max_abs_space_effect_S3minusS2_at_dt_b_pctP", "max_abs_time_effect_S2_dt_b_minus_dt_a_pctP"], rows)
rows = []
for dt, s1 in (("0.0125", "S1_dt0.0125"), ("0.00625", "S1_dt0.00625_r2")):
    for i, x in enumerate(X):
        a_ = series(s1, i, 1.5); b_ = series("S2_dt" + dt, i, 1.5); c_ = series("S3_dt" + dt, i, 1.5)
        k = int(np.argmax(np.abs(c_ - b_)))
        rows.append([dt, x, fl(float(np.abs(b_ - a_).max()), 5), fl(float(np.abs(c_ - b_).max()), 5), fl(float(np.linspace(0, 1.5, 6001)[k]), 5), fl(float(np.abs(c_ - a_).max()), 5)])
wr("history_space_effects.csv", ["dt_s", "x_m", "Linf_S2_minus_S1_pctP", "Linf_S3_minus_S2_pctP", "t_of_max_S3_minus_S2_s", "Linf_S3_minus_S1_pctP"], rows)

# ------------------------------------------------------------------ stage-3 variants against the S2 baseline at the same dt
rows = []
for n, c in cells.items():
    S, var, dt = parse(n)
    if var in ("baseline", "eq10", "rerun") or S != "S2": continue
    b = cells[f"S2_dt{dt:g}"]
    dv = [c.V("bottom", i, 0) - b.V("bottom", i, 0) for i in range(3)]
    fd = max(abs(float(c.bot[i](tt)) - float(b.bot[i](tt))) for i in range(3) for tt in FIXED_T) / P * 100
    tt = np.linspace(0, 1.5, 6001); e = np.abs(np.array([float(c.bot[1](q)) - float(b.bot[1](q)) for q in tt])).max() / P * 100
    st = c.d["steps"]
    rows.append([n, var, fl(dt), c.d["status"], c.d["all_converged"], max(s["newton"] for s in st), fl(min(s["jmin"] for s in st), 5), "%.2e" % max(abs(s["mass_cons"]) for s in st),
                 fl(dv[0], 5), fl(dv[1], 5), fl(dv[2], 5), fl(c.V("bottom", 1, 1) - b.V("bottom", 1, 1), 5), fl(c.V("bottom", 1, 2) - b.V("bottom", 1, 2), 5), fl(fd, 5), fl(float(e), 5)])
wr("stage3_variants_vs_S2.csv", ["cell", "variant", "dt_s", "status", "all_converged", "max_newton", "min_jacobian", "max_abs_q_in_out_wall", "dV1_x7.5_pctP", "dV1_x10_pctP", "dV1_x12.5_pctP",
                                  "dV2_x10_pctP", "dV3_x10_pctP", "max_abs_fixed_time_change_pctP", "Linf_history_change_x10_pctP"], rows)

# ------------------------------------------------------------------ H8 checks (every 10th step carries wall pressure and interface-force diagnostics)
rows = []
for n in ("S1_dt0.0125", "S1_dt0.00625_r2", "S2_dt0.0125", "S2_dt0.00625", "S2_dt0.003125", "S3_dt0.0125", "S3_dt0.00625", "S3_dt0.003125"):
    d = load(n); Sx = [s for s in d["steps"] if "pw" in s]; par = []; F = []; bal = []
    for s in Sx:
        x = np.array(s["pw_x"]); p = np.array(s["pw"]); m = (x > 5.0 + 1e-9) & (x < 15.0 - 1e-9); pp = p[m]
        d2 = pp[:-2] - 2 * pp[1:-1] + pp[2:]; jj = np.arange(len(d2))
        par.append(abs(np.mean((-1.0) ** jj * d2)) / 4.0 / max(np.abs(pp).max(), 1e-9))
        F.append((s["F_if_weak_y"], s["F_if_dir_y"], s["F_if_dir_pressure_y"])); bal.append(s["balance_resid"])
    F = np.array(F); sc = np.abs(F[:, 1]).max()
    gcl = max(abs(s["mass_res"]) for s in d["steps"]) / max(abs(s["dVdt"]) for s in d["steps"])
    rows.append([n, len(Sx), "%.3e" % max(par), "%.4f" % (np.abs(F[:, 0] - F[:, 1]).max() / sc * 100), "%.4f" % (np.abs(F[:, 2] - F[:, 1]).max() / sc * 100), "%.2e" % max(bal), "%.1f" % sc,
                 "%.2e" % max(abs(s["mass_cons"]) for s in d["steps"]), "%.2e" % gcl])
wr("gates_H8.csv", ["cell", "diagnostic_steps", "max_pressure_parity_over_max_p", "weak_vs_direct_interface_force_pct_of_max_F", "pressure_part_vs_direct_pct_of_max_F", "max_balance_residual_N_per_m",
                    "max_abs_interface_force_N_per_m", "max_abs_q_in_out_wall", "bdf_gcl_residual_over_max_abs_dVdt"], rows)

# ------------------------------------------------------------------ MF1: open-loop beam against continuum replay
c = json.load(open(f"{DATA}/data_MF1_S2_dt0.00625.json")); s = load("S2_dt0.00625"); r = json.load(open(f"{DATA}/replay_MF1.json"))
b320 = json.load(open(f"{DATA}/replay_MF1_beam320.json"))      # 320 beam elements (PREREG_MF1_BEAM320.md); same driver code, separate job
rp160 = json.load(open(f"{DATA}/replay_MF1_beam_repro.json")); rp80 = json.load(open(f"{DATA}/replay_MF1_beam_repro80.json"))    # reproduction jobs of ne160 and ne80
r["beam"]["ne320_d2dt2"] = b320["beam"]["ne320_d2dt2"]
ctl160 = json.load(open(f"{DATA}/replay_MF1_beam160_amend1.json")); fail320 = json.load(open(f"{DATA}/replay_MF1_beam320_FAILED_tol1e-5_run18904.json")); diag320 = json.load(open(f"{DATA}/replay_MF1_diag320.json"))
dc = np.array([x["wall_d"] for x in c["steps"]]); mc = np.array([x["midplane_d"] for x in c["steps"]]); ds_ = np.array([x["wall_d"] for x in s["steps"]])
t = np.array([0.0] + [x["t"] for x in c["steps"]])
cb = np.array([r["continuum"]["d"][f"{x}_bot"] for x in X]).T; cm = np.array([r["continuum"]["d"][f"{x}_mid"] for x in X]).T
rows = [["C1_rerun_vs_stored_S2_max_abs_diff_over_max_abs_d", "%.3e" % (abs(dc - ds_).max() / abs(dc).max()), "tolerance 1e-8"],
        ["C2_replayed_continuum_vs_coupled_bottom_max_abs_diff_over_max_abs_d", "%.3e" % (abs(cb[1:] - dc).max() / abs(dc).max()), "tolerance 1e-6"],
        ["C2_replayed_continuum_vs_coupled_midplane_max_abs_diff_over_max_abs_d", "%.3e" % (abs(cm[1:] - mc).max() / abs(mc).max()), "tolerance 1e-6"],
        ["C2_abs_max_diff_bottom_m", "%.3e" % abs(cb[1:] - dc).max(), ""],
        ["rerun_status", c["status"] + " all_converged=" + str(c["all_converged"]), "%d steps, max Newton %d, min Jacobian %.3f, max |q-sum| %.2e" % (len(c["steps"]), max(x["newton"] for x in c["steps"]), min(x["jmin"] for x in c["steps"]), max(abs(x["mass_cons"]) for x in c["steps"]))],
        ["beam_reproduction_ne160_max_abs_diff_m", "%.3e" % abs(np.array(rp160["beam"]["ne160_d2dt2"]["Y"]) - np.array(r["beam"]["ne160_d2dt2"]["Y"])).max(), "rerun from the saved load file, job 18902"],
        ["beam_reproduction_ne80_max_abs_diff_m", "%.3e" % abs(np.array(rp80["beam"]["ne80_d2dt2"]["Y"]) - np.array(r["beam"]["ne80_d2dt2"]["Y"])).max(), "job 18903"],
        ["beam320_run", b320["status"], "ok=%s reached_end=%s steps=%d of %d max Newton %d wall %.1f s (job 18904)" % (b320["beam"]["ne320_d2dt2"]["ok"], b320["beam"]["ne320_d2dt2"]["reached_end"], b320["beam"]["ne320_d2dt2"]["steps_done"], b320["nsteps"], b320["beam"]["ne320_d2dt2"]["newton_max"], b320["beam"]["ne320_d2dt2"]["wall_s"])],
        ["amend1_control_ne160_max_abs_diff_m", "%.3e" % abs(np.array(ctl160["beam"]["ne160_d2dt2"]["Y"]) - np.array(r["beam"]["ne160_d2dt2"]["Y"])).max(), "ne160 rerun with the amended criterion (residual < 1e-4 and update < 1e-10 m), job 18911; AMEND1 asked <= 1e-9: exceeded by 25 %, see README"],
        ["amend1_ne320_vs_floor_run_max_abs_diff_m", "%.3e" % abs(np.array(b320["beam"]["ne320_d2dt2"]["Y"]) - np.array(diag320["beam"]["ne320_d2dt2"]["Y"])).max(), "amended ne320 (job 18912) against the 30-iteration state of diagnostic job 18908 (residual at its roundoff floor, update 4e-14 m)"],
        ["first_ne320_run_tol1e-5", fail320["status"], "job 18904: all 240 steps reached the 30-iteration cap (residual floor 1.0-1.6e-5 > 1e-5); result retained, not used"],
        ["amend1_ne320_final_residual_and_update_max", "%.2e / %.2e" % (max(q[-1] for q in b320["beam"]["ne320_diag"]["res"]), max(q[-1] for q in b320["beam"]["ne320_diag"]["upd"])), "residual 2-norm (floor 1.0-1.6e-5) / max-norm of the last Newton update in m, over all 240 steps"],
        ["replay_continuum_max_newton", str(max(r["continuum"]["newton"])), "ok=" + str(r["continuum"]["ok"])]]
wr("mf1_controls.csv", ["item", "value", "note"], rows)
spl = lambda y: CubicSpline(t, y)
v123 = lambda y: [e[1] / P * 100 for e in V123(spl(y), 1.5)]
rows = []; res = {}
for nm, b in r["beam"].items():
    Y = np.array(b["Y"]); res[nm] = Y
    for i, x in enumerate(X):
        e = (Y[:, i] - cm[:, i]) / P * 100; eb = (Y[:, i] - cb[:, i]) / P * 100; vb, vc = v123(Y[:, i]), v123(cm[:, i])
        rows.append([nm, b["ok"], b["newton_max"], x, fl(float(abs(e).max()), 5), fl(float(abs(eb).max()), 5), fl(float(np.sqrt(np.mean(e ** 2))), 5), fl(vb[0] - vc[0], 5), fl(vb[1] - vc[1], 5), fl(vb[2] - vc[2], 5)])
wr("mf1_gap.csv", ["beam_variant", "beam_converged", "beam_max_newton", "x_m", "gap_vs_continuum_midplane_max_abs_pctP", "gap_vs_continuum_bottom_max_abs_pctP", "history_L2_pctP",
                   "dV1_beam_minus_continuum_pctP", "dV2_pctP", "dV3_pctP"], rows)
rows = []
for i, x in enumerate(X):
    rows.append([x, fl(float(abs(res["ne160_d2dt2"][:, i] - res["ne80_d2dt2"][:, i]).max() / P * 100), 5), fl(float(abs(res["ne160_d2dt2"][:, i] - res["ne160_newmark"][:, i]).max() / P * 100), 5), fl(float(abs(cm[:, i] - cb[:, i]).max() / P * 100), 5)])
# refinement table 80 / 160 / 320: first troughs, signed difference, gap, changes, time alignment
rows = []; vt = lambda y: V123(spl(y), 1.5)[0]
for i, x in enumerate(X):
    c1 = vt(cm[:, i]); prev = None
    for ne in (80, 160, 320):
        Y = res["ne%d_d2dt2" % ne][:, i]; b1 = vt(Y); gmax = float(abs((Y - cm[:, i]) / P * 100).max())
        dv = (b1[1] - c1[1]) / P * 100
        rows.append([x, ne, fl(b1[1] / P * 100, 6), fl(b1[0], 5), fl(c1[1] / P * 100, 6), fl(c1[0], 5), fl(dv, 5), fl(-dv, 5), fl(gmax, 5),
                     "" if prev is None else fl(b1[1] / P * 100 - prev[0], 5), "" if prev is None else fl(-dv + prev[1], 5), "" if prev is None else fl(gmax - prev[2], 5),
                     "" if prev is None else fl(float(abs(Y - prev[3]).max() / P * 100), 5)])
        prev = (b1[1] / P * 100, -dv, gmax, Y)
wr("mf1_beam_refinement.csv", ["x_m", "beam_elements", "beam_V1_pctP", "t_V1_beam_s", "continuum_midplane_V1_pctP", "t_V1_continuum_s", "dV1_beam_minus_continuum_pctP", "continuum_minus_beam_pctP", "gap_max_abs_over_t_pctP",
                               "change_of_beam_V1_vs_previous_level", "change_of_continuum_minus_beam_vs_previous_level", "change_of_gap_vs_previous_level", "beam_vs_previous_beam_max_abs_pctP"], rows)
rows = []
for i, x in enumerate(X):
    rows.append([x, fl(float(abs(res["ne160_d2dt2"][:, i] - res["ne80_d2dt2"][:, i]).max() / P * 100), 5), fl(float(abs(res["ne160_d2dt2"][:, i] - res["ne160_newmark"][:, i]).max() / P * 100), 5),
                 fl(float(abs(cm[:, i] - cb[:, i]).max() / P * 100), 5), fl(float(abs(res["ne320_d2dt2"][:, i] - res["ne160_d2dt2"][:, i]).max() / P * 100), 5)])
wr("mf1_resolution.csv", ["x_m", "beam_ne160_vs_ne80_max_abs_pctP", "beam_d2dt2_vs_newmark_ne160_max_abs_pctP", "continuum_midplane_vs_bottom_max_abs_pctP", "beam_ne320_vs_ne160_max_abs_pctP"], rows)

# ------------------------------------------------------------------ P2-triangle wall alone (COL-3 S1 static gate), recomputed from the raw JSON
d = json.load(open(f"{EV}/C3/benchmarks/collapsibleChannel/data_S1_static.json")); rows = []
for q in d["results"]:
    if q["p"] != 200.0: continue
    _, nx, ny, _ = q["case"].split(":"); rows.append([q["case"], nx, ny, q["ndof"], fl(q["u_x10"], 10), fl((q["u_x10"] + 0.14362832) / 0.14362832 * 100, 5), fl(q["u_x75"], 10), fl((q["u_x75"] + 0.09645925) / 0.09645925 * 100, 5), q["converged"]])
wr("wall_static_p2.csv", ["case_deg_nx_ny_load", "nx", "ny", "ndof", "u_mid_m", "rel_mid_vs_beam_pct", "u_quarter_m", "rel_quarter_vs_beam_pct", "converged"], rows)

# ------------------------------------------------------------------ cost from sacct
rows = []
for rr in sorted(sacct_rows, key=lambda q: int(q["JobID"])):
    jn = rr["JobName"]
    camp = {"col1": "COL-1 preflight", "col2": "COL-2 preflight", "col3": "COL-3 preflight/solver", "colP": "COL-P production", "colA": "auditor"}[jn[:4]]
    if jn.startswith("colP_MF1"): camp = "COL-P MF1"
    rows.append([rr["JobID"], jn, camp, rr["AllocCPUS"], rr["Submit"], rr["Elapsed"], rr["CPUTimeRAW"], "%.3f" % (int(rr["CPUTimeRAW"]) / 3600.0), rr["State"]])
wr("cost_cpuh.csv", ["slurm_job", "job_name", "campaign", "alloc_cpus", "submit_time_node_local", "elapsed", "cpu_s_sacct", "cpu_h", "state"], rows)
print("ok", len(cells), "cells")
