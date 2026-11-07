"""Regenerate the derived CSVs of this package from the input CSVs in ../data (standard library only).
usage (from this package directory):  python3 -B -I scripts/analyze_cavity.py
Inputs : plate_only_nt32_control.csv, coupled_levels.csv, s4f_reference_variants_levels.csv, plate_only_ladder.csv, fluid_only_ladder.csv, context_references.csv, cost_cpuh.csv
Outputs: plate_only_q2_nt32_orders.csv, williams_exponent.csv, coupled_orders.csv, coupled_uy_limit_vs_order.csv, reaction_deficit.csv, reference_comparison.csv, context_comparison.csv,
         plate_only_orders.csv, fluid_only_orders.csv, cost_summary.csv
Acceptance rule: scripts/convergence.py (monotone, 0.5 <= p <= 3, Fs = 1.25); a triple that fails is reported as NOT ACCEPTED with no limit."""
import csv, os, sys, math
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import convergence as C
D = os.path.join(HERE, "..", "data")


def rd(n): return list(csv.DictReader(open(os.path.join(D, n))))


def wr(n, header, rows):
    with open(os.path.join(D, n), "w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n"); w.writerow(header)
        for r in rows: w.writerow([("%.10g" % v) if isinstance(v, float) else v for v in r])


def tri(f, h, label):
    a = C.analyse(f, h, label)
    t = a["triples"][-1]
    if a.get("accepted"):
        return [t["p"], t["status"], "accepted", a["f_extrap"], 100 * a["gci_fine"], 100 * a["rel_diff_fine_vs_extrap"], a["asymptotic_ratio"]]
    return [t["p"], t["status"], "NOT ACCEPTED (%s)" % a.get("reason"), "", "", "", ""]


lv = rd("coupled_levels.csv"); H = [float(r["dx_m"]) for r in lv]
col = lambda k: [float(r[k]) for r in lv]

# 1. orders for the coupled FEniCSx M1 sequence: both triples
rows = []
for q, k in (("uy_A", "uy_A_m"), ("ux_A", "ux_A_m"), ("Fy_total", "Fy_N_per_m"), ("Fy_pressure", "Fy_pressure")):
    v = col(k)
    for name, ix in (("L2-L4", [1, 2, 3]), ("L1-L3", [0, 1, 2])):
        rows.append([q, name] + tri([v[i] for i in ix], [H[i] for i in ix], q))
wr("coupled_orders.csv", ["quantity", "triple", "observed_order_p", "status", "acceptance", "richardson_limit", "GCI_fine_percent", "fine_minus_limit_percent", "asymptotic_ratio"], rows)

# 2. sensitivity of the u_y and Fy limits to the assumed order (L2-L4, r = 2)
rows = []
for q, k in (("uy_A", "uy_A_m"), ("Fy_total", "Fy_N_per_m")):
    v = col(k); f4, f3 = v[3], v[2]
    base = C.richardson(f4, f3, tri([v[1], v[2], v[3]], H[1:], q)[0], 2.0)
    for p in (1.0, 1.28, 1.5, 2.0, 2.18):
        lim = C.richardson(f4, f3, p, 2.0)
        rows.append([q, p, lim, 100 * (lim - base) / abs(base)])
wr("coupled_uy_limit_vs_order.csv", ["quantity", "assumed_order_p", "limit_(L3,L4_two_level_Richardson)", "change_vs_observed_order_limit_percent"], rows)

# 3. support reactions: artefact check (clamp-corner dofs excluded from the reported reactions)
rows = []
prev = None
for r in lv:
    ry = float(r["R_left_y"]) + float(r["R_right_y"]); fy = float(r["Fy_N_per_m"]); d = ry + fy
    rows.append([r["level"], float(r["dx_m"]), ry, fy, d, d / float(r["dx_m"]), (prev / d) if prev else ""])
    prev = d
wr("reaction_deficit.csv", ["level", "dx_m", "R_left_y_plus_R_right_y", "Fy_direct", "sum_(R+Fy)", "sum_divided_by_dx", "ratio_to_previous_coarser_level"], rows)

# 4. comparison with the solids4foam variants (limits of the finest triple, each from its own levels)
fx_uy = tri([float(lv[i]["uy_A_m"]) for i in (1, 2, 3)], H[1:], "uy")[3]
fx_fy = tri([float(lv[i]["Fy_N_per_m"]) for i in (1, 2, 3)], H[1:], "Fy")[3]
ref = rd("s4f_reference_variants_levels.csv"); vs = []
for r in ref:
    if r["variant"] not in vs: vs.append(r["variant"])
crit = {
    "partitioned_strict": "partitioned, requireAllResidualMeasures; dt 0.5/0.25/0.125/0.0625 s halves with dx; 200 s; steady = plateau spread over last 10 % <= 1e-4",
    "partitioned_strict_paper_rounded_levels": "same values as partitioned_strict with Fy rounded to 4 decimals (the paper's reference_levels.csv)",
    "partitioned_legacy_test": "partitioned, default residual test; level-4 coupling not converged (r1 ~ 1): Fy not accepted",
    "monolithic_pre597": "quasi-monolithic Newton before solids4foam #597; dt 0.5 s; inlet ramp 20 s; 200 s",
    "monolithic_597": "quasi-monolithic Newton including #597; dt 0.5 s; inlet ramp 20 s; 300 s",
    "monolithic_618": "quasi-monolithic on the #618 build (paper PR 23 file); dt not recorded in that file"}
rows = []
for v in vs:
    rr = [x for x in ref if x["variant"] == v]
    uy = [float(x["uy_top_m"]) for x in rr]; fy = [float(x["Fy_N_per_m"]) for x in rr]; h = [float(x["dx_m"]) for x in rr]
    a, b = tri(uy[1:], h[1:], "uy"), tri(fy[1:], h[1:], "Fy")
    du = 100 * (fx_uy - a[3]) / a[3] if a[3] != "" else ""
    df = 100 * (fx_fy - b[3]) / b[3] if b[3] != "" else ""
    ux = [x["ux_top_m"] for x in rr]
    if all(u != "" for u in ux):
        ax = tri([float(u) for u in ux[1:]], h[1:], "ux"); fx_ux = tri([float(lv[i]["ux_A_m"]) for i in (1, 2, 3)], H[1:], "ux")[3]
        dux = 100 * (fx_ux - ax[3]) / ax[3] if ax[3] != "" else ""; gux = ax[4]; dux4 = 100 * (float(lv[3]["ux_A_m"]) - float(ux[3])) / float(ux[3])
    else: dux = dux4 = gux = ""
    rows.append([v, crit[v], a[0], a[3], a[4], b[0] if b[3] != "" else "", b[3], b[4] if b[3] != "" else "", du, df,
                 100 * (float(lv[3]["uy_A_m"]) - uy[3]) / uy[3], 100 * (float(lv[3]["Fy_N_per_m"]) - fy[3]) / fy[3], dux, dux4, gux])
wr("reference_comparison.csv", ["reference_variant", "dt_and_steady_criterion", "ref_uy_order_p", "ref_uy_limit", "ref_uy_GCI_percent", "ref_Fy_order_p", "ref_Fy_limit", "ref_Fy_GCI_percent",
                                "FEniCSx_M1_limit_minus_ref_limit_uy_percent", "FEniCSx_M1_limit_minus_ref_limit_Fy_percent", "FEniCSx_L4_minus_ref_L4_uy_percent", "FEniCSx_L4_minus_ref_L4_Fy_percent", "FEniCSx_M1_limit_minus_ref_limit_ux_percent", "FEniCSx_L4_minus_ref_L4_ux_percent", "ref_ux_GCI_percent"], rows)

# 5. other references: context only, no pass/fail, not model-equivalent
rows = []
for r in rd("context_references.csv"):
    rows.append([r["reference"], float(r["uy_A_m"]), 100 * (float(r["uy_A_m"]) - fx_uy) / abs(fx_uy), float(r["Fy_N_per_m"]), 100 * (float(r["Fy_N_per_m"]) - fx_fy) / abs(fx_fy), "context only; no equivalence claim"])
wr("context_comparison.csv", ["reference", "uy_A_m", "difference_to_FEniCSx_M1_uy_limit_percent_(ref-fx)/|fx|", "Fy_N_per_m", "difference_to_FEniCSx_M1_Fy_limit_percent", "note"], rows)

# 6. plate-only (C1) and fluid-only (C2) ladders
pl = rd("plate_only_ladder.csv"); rows = []
ser = []
for r in pl:
    if r["series"] not in ser: ser.append(r["series"])
nt = [2, 4, 8, 16]
for s in ser:
    v = [float(r["uy_A_m"]) for r in pl if r["series"] == s]
    rows.append([s] + tri(v, [1.0 / n for n in nt], s))
up = {s: [float(r["uy_A_m"]) for r in pl if r["series"] == s][-1] for s in ser}
rows.append(["ratio dead_p_Q3 / follower_Q3 at nt = 16 (reference-pressure load over follower)", up["dead_p_Q3"] / up["follower_Q3"] - 1.0, "", "", "", "", "", ""])
rows.append(["ratio dead_vec_Q3 / follower_Q3 at nt = 16", up["dead_vec_Q3"] / up["follower_Q3"] - 1.0, "", "", "", "", "", ""])
wr("plate_only_orders.csv", ["series", "observed_order_p_(h=1/nt)", "status", "acceptance", "richardson_limit", "GCI_fine_percent", "fine_minus_limit_percent", "asymptotic_ratio"], rows)
fl = rd("fluid_only_ladder.csv"); rows = []
for fv in ("flat", "dump"):
    rr = [r for r in fl if r["floor"] == fv]; n = [float(r["cells"]) for r in rr]; h = [C.h_from_cells(x) for x in n]
    for q, k in (("Fy", "Fy_N_per_m"), ("p_floor_mean", "p_floor_mean")):
        rows.append([fv, q] + tri([float(r[k]) for r in rr], h, q))
fy_flat = float([r for r in fl if r["floor"] == "flat"][-1]["Fy_N_per_m"]); fy_dump = float([r for r in fl if r["floor"] == "dump"][-1]["Fy_N_per_m"])
rows.append(["flat vs dump finest", "Fy flat/dump - 1 (percent)", 100 * (fy_flat / fy_dump - 1), "", "", "", "", "", ""])
rows.append(["dump finest vs s4f strict level 2 (-5.238953, 80x8 plate)", "Fy relative difference (percent)", 100 * (fy_dump - (-5.238953)) / 5.238953 * -1, "", "", "", "", "", ""])
wr("fluid_only_orders.csv", ["floor", "quantity", "observed_order_p", "status", "acceptance", "richardson_limit", "GCI_fine_percent", "fine_minus_limit_percent", "asymptotic_ratio"], rows)

# 6b. uniform-mesh plate-only Q2 ladder extended to nt = 32 (h = 1/nt) and the Williams exponent of the linearised 90-degree clamped/free wedge
nt = rd("plate_only_nt32_control.csv"); rows = []
for ld in ("uniform", "follower"):
    v = [float(r["uy_A_m"]) for r in nt if r["load"] == ld and r["source"].endswith("(job 16716)") or (r["load"] == ld and r["nt"] in ("4", "8", "16"))]
    for name, ix in (("nt 4,8,16", [0, 1, 2]), ("nt 8,16,32", [1, 2, 3])):
        rows.append([ld, name] + tri([v[i] for i in ix], [1.0 / [4, 8, 16, 32][i] for i in ix], ld))
wr("plate_only_q2_nt32_orders.csv", ["load", "triple", "observed_order_p", "status", "acceptance", "richardson_limit", "GCI_fine_percent", "fine_minus_limit_percent", "asymptotic_ratio"], rows)
def _root(k):
    lo, hi = 0.3, 1.0
    g = lambda l: 1 + k * k + 2 * k * math.cos(math.pi * l) - 4 * l * l
    for _ in range(200):
        m = 0.5 * (lo + hi)
        if g(lo) * g(m) <= 0: hi = m
        else: lo = m
    return 0.5 * (lo + hi)
wr("williams_exponent.csv", ["poisson_ratio", "kappa_plane_strain_3-4nu", "lambda_smallest_root", "two_lambda"], [[nu, 3 - 4 * nu, _root(3 - 4 * nu), 2 * _root(3 - 4 * nu)] for nu in (0.3, 0.4, 0.45, 0.49)])

# 6c. sensitivity of the strict-partitioned limit to an unseen corner term (a least-squares fit, NOT a measurement): S = L + a h^2 + s c h^(2 lambda), c from the FEniCSx four-level fit
def _lsq(cols, y):
    n = len(cols); A = [[sum(ci * cj for ci, cj in zip(cols[i], cols[j])) for j in range(n)] for i in range(n)]; b = [sum(ci * yi for ci, yi in zip(cols[i], y)) for i in range(n)]
    for i in range(n):
        piv = A[i][i]
        for j in range(i + 1, n):
            f_ = A[j][i] / piv; A[j] = [x - f_ * w for x, w in zip(A[j], A[i])]; b[j] -= f_ * b[i]
    sol = [0.0] * n
    for i in reversed(range(n)): sol[i] = (b[i] - sum(A[i][j] * sol[j] for j in range(i + 1, n))) / A[i][i]
    return sol
lam2 = 1.2998200018987194; hh = H; fe_uy = col("uy_A_m")
c_fe = _lsq([[1.0] * 4, [h ** lam2 for h in hh]], fe_uy)[1]
sy = [float(x["uy_top_m"]) for x in ref if x["variant"] == "partitioned_strict"]
base = None; rows = []
for sc in (0.0, 0.5, 1.0, 2.0):
    L_ = _lsq([[1.0] * 4, [h ** 2 for h in hh]], [v - sc * c_fe * h ** lam2 for v, h in zip(sy, hh)])[0]
    if base is None: base = L_
    rows.append([sc, c_fe, L_, 100 * (L_ - base) / abs(base), 100 * (L_ - float(tri([float(lv[i]["uy_A_m"]) for i in (1, 2, 3)], H[1:], "uy")[3])) / abs(L_)])
wr("s4f_corner_term_sensitivity.csv", ["corner_term_scale_s", "FEniCSx_corner_amplitude_c", "s4f_strict_limit_from_fit", "shift_vs_s0_fit_percent_(negative=larger_|uy|)", "fit_limit_minus_FEniCSx_L2-L4_limit_percent_(over_|fit|)"], rows)

# 7. cost
cs = rd("cost_cpuh.csv")
prod = sum(float(r["sacct_cpu_h_elapsed_x_alloccpus"]) for r in cs if r["group"] == "production")
led = sum(float(r["ledger_cpu_h_wall_x_ranks"]) for r in cs if r["group"] in ("production", "preflight/geometry/scan/rigid"))
oth = sum(float(r["ledger_cpu_h_wall_x_ranks"]) for r in cs if r["group"] == "preflight/geometry/scan/rigid")
wr("cost_summary.csv", ["item", "cpu_h"], [["production runs P1-P4,T1,S1,M2 (sacct)", prod], ["coupled campaign incl. geometry/preflight/scan/rigid/scaling, all jobs incl. failed (ledger wall x ranks)", led],
                                           ["of which non-production (ledger)", oth], ["C1/C2 plate-only and fluid-only preflight (sacct)", 0.263]])
print("ok")
