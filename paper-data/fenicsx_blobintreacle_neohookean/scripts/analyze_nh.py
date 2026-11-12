#!/usr/bin/env python3
"""Derived tables of the neo-Hookean blobInTreacle supplement from the input CSVs in data/.

usage (from this directory): python3 -B -I scripts/analyze_nh.py
Inputs (data/): fe_steady_levels.csv, fe_levels_t1_t3.csv, fe_time_study_t1.csv, fe_run_validity.csv (regenerable from the raw runs with
extract_from_raw.py), s4f_rerun_steady.csv, s4f_rerun_time.csv (extracted from the solids4foam reruns in the solver repository, not regenerable
from this repository), solid_only_apex.csv, law_checks.csv (extract_phase0.py), cpu_ledger_sacct.csv, pr_text_values.csv (transcribed).
Outputs (data/): the files named in README.md under 'Derived'.  Order and limits use scripts/convergence.py (rules of the main article:
monotone, 0.5 <= p <= 3, Fs = 1.25, h = N^-1/2)."""
import csv, datetime as dt, hashlib, math, os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from convergence import analyse
D = os.path.join(HERE, "..", "data") + os.sep
LV = ["0.5", "1", "2", "4"]; CELLS = {"0.5": 150, "1": 600, "2": 2400, "4": 9600}; H = [CELLS[l] ** -0.5 for l in LV]
OBS = ["ux", "uy", "Fx", "Fy", "dp"]


def rd(name): return list(csv.DictReader(open(D + name)))


def num(x): return float(x)


def wr(name, header, rows):
    with open(D + name, "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(header)
        for r in rows: w.writerow(r)


def g(x, n=6): return ("%." + str(n) + "g") % x


fe = {r["level"]: r for r in rd("fe_steady_levels.csv")}
fe13 = {(r["level"], r["t"]): r for r in rd("fe_levels_t1_t3.csv")}
s4 = {r["level"]: r for r in rd("s4f_rerun_steady.csv") if r["level"] in LV}
pr = {(r["level"], r["obs"]): num(r["value"]) for r in rd("pr_text_values.csv") if r["level"] in LV}

# 1. matched-level differences, (FE - s4f)/s4f in percent. FE at t=6 (steady) and at t=3; s4f at t=3.
rows = []
for l in LV:
    for k in OBS:
        a, b, c = num(fe[l][k]), num(s4[l][k]), num(fe13[(l, "3")][k])
        rows.append([l, CELLS[l], k, g(a, 8), g(c, 8), g(b, 8), g(100 * (a / b - 1), 4), g(100 * (c / b - 1), 4), g(pr[(l, k)], 6) if (l, k) in pr else ""])
wr("steady_matched_differences.csv", ["level", "cells", "obs", "FE_t6", "FE_t3", "s4f_rerun_t3", "diff_FE_t6_vs_s4f_pct", "diff_FE_t3_vs_s4f_pct", "PR_text_value"], rows)
# rerun reproduces the PR-text numbers (max relative deviation over the printed digits)
dev = max(abs(num(s4[l][k]) / v - 1) for (l, k), v in pr.items() for l in [l] if (l, k) in pr and True) if False else max(abs(num(s4[l][k]) / pr[(l, k)] - 1) for (l, k) in pr)
# 2. orders, limits, GCI per observable and code, and difference of limits
conv, seqs = [], {}
for k in OBS:
    a_f = analyse([num(fe[l][k]) for l in LV], H, "FE " + k); a_s = analyse([num(s4[l][k]) for l in LV], H, "s4f " + k)
    seqs[k] = (a_f, a_s)
    row = [k]
    for a in (a_f, a_s):
        row += [g(a["triples"][0]["p"], 3), g(a["triples"][1]["p"], 3), g(a["f_extrap"], 7), g(100 * a["gci_fine"], 3)] if a.get("accepted") else ["", "", "", ""]
    d = abs(a_f["f_extrap"] - a_s["f_extrap"]); gc = math.hypot(a_f["gci_fine"] * abs(a_f["f_extrap"]), a_s["gci_fine"] * abs(a_s["f_extrap"]))
    sgn = 100 * (a_f["f_extrap"] / a_s["f_extrap"] - 1)
    row += [g(sgn, 3), g(d / gc, 3)]; conv.append(row)
wr("steady_convergence.csv", ["obs", "FE_p_coarse_triple", "FE_p_finest_triple", "FE_limit", "FE_GCI_pct", "s4f_p_coarse_triple", "s4f_p_finest_triple", "s4f_limit", "s4f_GCI_pct",
                              "limit_diff_FE_minus_s4f_pct", "limit_diff_over_combined_GCI"], conv)
# the PR-text u_x sequence with the same rule
a_pr = analyse([pr[(l, "ux")] for l in LV], H, "PR ux")
# 3. sensitivity of the limit difference to the assumed order: limit = x4 + (x4 - x2)/(2^p - 1); 4-level least squares fit f0 + c h^p (p scanned 0.5..3 in steps of 0.01)
def rich(f2, f4, p): return f4 + (f4 - f2) / (2 ** p - 1)


def lsq(v):
    best = None
    for i in range(0, 251):
        p = 0.5 + 0.01 * i; hp = [x ** p for x in H]; n = 4
        sh, sv, shh, shv = sum(hp), sum(v), sum(x * x for x in hp), sum(x * y for x, y in zip(hp, v)); den = n * shh - sh * sh
        c = (n * shv - sh * sv) / den; f0 = (sv - c * sh) / n; e = sum((f0 + c * x - y) ** 2 for x, y in zip(hp, v))
        if best is None or e < best[0]: best = (e, p, f0)
    return best[1], best[2]


rows = []
for k in OBS:
    a = [num(fe[l][k]) for l in LV]; b = [num(s4[l][k]) for l in LV]
    for p in (1.0, 1.5, 2.0):
        la, lb = rich(a[2], a[3], p), rich(b[2], b[3], p); rows.append([k, "fixed p=%g" % p, g(la, 7), g(lb, 7), g(100 * (la / lb - 1), 3)])
    pa, fa = lsq(a); pb, fb = lsq(b); rows.append([k, "least squares f0+c h^p (FE p=%.2f, s4f p=%.2f)" % (pa, pb), g(fa, 7), g(fb, 7), g(100 * (fa / fb - 1), 3)])
    rows.append([k, "data-driven p, finest triple (preregistered rule)", g(seqs[k][0]["f_extrap"], 7), g(seqs[k][1]["f_extrap"], 7), g(100 * (seqs[k][0]["f_extrap"] / seqs[k][1]["f_extrap"] - 1), 3)])
    rows.append([k, "no extrapolation: matched x4 mesh", g(a[3], 7), g(b[3], 7), g(100 * (a[3] / b[3] - 1), 3)])
wr("limit_sensitivity_to_order.csv", ["obs", "estimator", "FE", "s4f", "diff_FE_minus_s4f_pct"], rows)
sens = {}
for r in rows:
    if r[1].startswith("fixed p"): sens.setdefault(r[0], []).append(float(r[4]))
# 4. transient at t = 1 on the mesh-study runs (FE dt = 0.04 vs s4f dt = 0.01, x4 0.005)
rows = []
for l in LV:
    for k in OBS[:4]:
        a, b = num(fe13[(l, "1")][k]), num(s4[l][k + "_t1"]); rows.append([l, k, g(a, 8), g(b, 8), g(100 * (a / b - 1), 3)])
wr("transient_t1_matched_differences.csv", ["level", "obs", "FE_t1_dt0.04", "s4f_t1", "diff_pct"], rows)
# 5. time study at t = 1 (x1)
fet = rd("fe_time_study_t1.csv"); s4t = rd("s4f_rerun_time.csv")


def orders(v): return [math.log(a / b) / math.log(2) for a, b in zip(v, v[1:])]


fux = [num(r["dux"]) for r in fet[1:]]; sux = [num(r["dux"]) for r in s4t[1:]]
fif = [num(r["iface_l2_vertices"]) for r in fet[1:]]; sif = [num(r["iface_l2"]) for r in s4t[1:]]
rows = []
for i, r in enumerate(fet):
    s = [x for x in s4t if abs(num(x["dt"]) - num(r["dt"])) < 1e-12]
    rows.append([r["dt"], g(num(r["ux"]), 8), g(num(s[0]["ux"]), 8) if s else "", g(100 * (num(r["ux"]) / num(s[0]["ux"]) - 1), 3) if s else "",
                 g(100 * (num(r["uy"]) / num(s[0]["uy"]) - 1), 3) if s else "", g(100 * (num(r["Fx"]) / num(s[0]["Fx"]) - 1), 3) if s else "", g(num(r["dux"]), 4) if r["dux"] else "",
                 g(num(s[0]["dux"]), 4) if s and s[0]["dux"] else ""])
wr("time_study_t1.csv", ["dt", "FE_ux_t1", "s4f_ux_t1", "ux_FE_vs_s4f_pct", "uy_FE_vs_s4f_pct", "Fx_FE_vs_s4f_pct", "FE_abs_change_ux", "s4f_abs_change_ux"], rows)
rows = [["FE u_x change"] + [g(x, 3) for x in orders(fux)], ["FE interface change (vertices)"] + [g(x, 3) for x in orders(fif)], ["s4f u_x change"] + [g(x, 3) for x in orders(sux)], ["s4f interface change"] + [g(x, 3) for x in orders(sif)]]
wr("time_orders_t1.csv", ["series", "order_pair1", "order_pair2", "order_pair3", "order_pair4"], [r + [""] * (5 - len(r)) for r in rows])
# error against the finest run (reference), ratios of successive absolute errors
rows = []
for k in OBS:
    v = [num(r[k]) for r in fet]; ref = v[-1]; e = [x - ref for x in v[:-1]]
    rows.append([k, g(ref, 8)] + [g(x, 3) for x in e] + [g(abs(e[i] / e[i + 1]), 3) for i in range(4)])
wr("time_error_vs_finest.csv", ["obs", "ref_dt0.0078125", "err_dt0.25", "err_dt0.125", "err_dt0.0625", "err_dt0.03125", "err_dt0.015625", "ratio1", "ratio2", "ratio3", "ratio4"], rows)
# 5b. addenda: tight SNES tolerance and 16-substep start-up
ad = rd("fe_addenda_t1_values.csv"); A = {}
for r in ad: A.setdefault(r["series"], []).append(r)
rel = max(abs(num(a[k]) / num(b[k]) - 1) for a, b in zip(A["A1_tight_SNES"], A["main"]) for k in OBS)
a2 = [abs(num(b["ux"]) - num(a["ux"])) for a, b in zip(A["A2_16_substeps"], A["A2_16_substeps"][1:])]
a1 = [abs(num(b["ux"]) - num(a["ux"])) for a, b in zip(A["A1_tight_SNES"], A["A1_tight_SNES"][1:])]
wr("time_addenda.csv", ["check", "value"], [["A1: max relative difference of all five observables at t=1, tight SNES vs main, six dt", "%.2e" % rel], ["A1: u_x change orders, tight SNES"] + [", ".join(g(x, 3) for x in orders(a1))],
   ["A2: u_x change orders, 16 backward-Euler substeps in step 1", ", ".join(g(x, 3) for x in orders(a2))], ["A2: u_x change values", ", ".join(g(x, 4) for x in a2)]])
print("A1 max rel diff %.2e; A2 orders %s" % (rel, [round(x, 2) for x in orders(a2)]))
# 6. solid-only gate
sol = rd("solid_only_apex.csv"); rows = []; sg = {}
for ld in ("L1", "L2"):
    for k in ("ux", "uy"):
        f_ = [num(r[k]) for r in sol if r["load"] == ld and r["code"] == "FEniCSx"]; s_ = [num(r[k]) for r in sol if r["load"] == ld and r["code"] == "solids4foam"]
        d = [100 * (a / b - 1) for a, b in zip(f_, s_)]
        inc_f = [f_[i + 1] - f_[i] for i in range(3)]; inc_s = [s_[i + 1] - s_[i] for i in range(3)]
        of = [abs(math.log(abs(inc_f[i] / inc_f[i + 1])) / math.log(2)) for i in range(2)]
        os_ = [math.log(abs(inc_s[i] / inc_s[i + 1])) / math.log(2) if inc_s[i + 1] != 0 else float("nan") for i in range(2)]
        mono = "monotone" if all(x > 0 for x in inc_s) or all(x < 0 for x in inc_s) else "non-monotone"
        rows.append([ld, k] + [g(x, 3) for x in d] + [g(x, 3) for x in of] + [g(x, 3) for x in os_] + [mono]); sg[(ld, k)] = (d, of, os_, mono)
wr("solid_only_differences.csv", ["load", "obs", "diff_x0.5_pct", "diff_x1_pct", "diff_x2_pct", "diff_x4_pct", "FE_increment_order_1", "FE_increment_order_2", "s4f_increment_order_1", "s4f_increment_order_2", "s4f_increments"], rows)
# 7. preregistered predictions: scored values
mx = lambda k, ls: max(abs(100 * (num(fe[l][k]) / num(s4[l][k]) - 1)) for l in ls)
sc = []
sc.append(["P1", "u_x within 3 % of s4f at every matched level and shrinking", "FAILED at x0.5; met at x1, x2, x4", "x0.5 %.2f %%; x1 %.2f %%; x2 %.2f %%; x4 %.2f %%" % tuple(100 * (num(fe[l]["ux"]) / num(s4[l]["ux"]) - 1) for l in LV)])
sc.append(["P2", "u_y within 8 %", "met", "max %.2f %% (x0.5)" % mx("uy", LV)])
sc.append(["P3", "F_x, F_y within 3 % at every matched level", "FAILED at x0.5; met at x1 (F_x close to the limit), x2, x4", "x0.5 F_x %.2f %%, F_y %.2f %%; x1 F_x %.2f %%, F_y %.2f %%" % (100 * (num(fe["0.5"]["Fx"]) / num(s4["0.5"]["Fx"]) - 1), 100 * (num(fe["0.5"]["Fy"]) / num(s4["0.5"]["Fy"]) - 1), 100 * (num(fe["1"]["Fx"]) / num(s4["1"]["Fx"]) - 1), 100 * (num(fe["1"]["Fy"]) / num(s4["1"]["Fy"]) - 1))])
tt = [100 * (num(r["ux"]) / num(s["ux"]) - 1) for r in fet for s in s4t if abs(num(r["dt"]) - num(s["dt"])) < 1e-12]
sc.append(["P4", "(a) FE u_x(1) within 3 % of s4f at dt <= 0.0625; (b) larger difference at dt = 0.25", "(a) met; (b) FAILED", "differences at dt 0.25..0.03125: " + ", ".join("%.2f %%" % x for x in tt)])
fl = 100 * (seqs["ux"][0]["f_extrap"] / a_pr["f_extrap"] - 1)
sc.append(["P5", "FE u_x limit within 3 % of the PR extrapolated 0.14024", "met", "%.6f m, %+.3f %% from %.6f m" % (seqs["ux"][0]["f_extrap"], fl, a_pr["f_extrap"])])
sc.append(["V4", "u_x, F_x, dp order 1.5-2; u_y, F_y 0.7-1.5", "u_x FAILED; F_x, dp, F_y above the range; u_y met", "finest-triple p: " + ", ".join("%s %.2f" % (k, seqs[k][0]["p"]) for k in OBS)])
o5 = orders(fux)
sc.append(["V5", "first two time-order pairs in [1.6, 3.0] (u_x and interface change)", "FAILED for u_x (second pair); interface met", "u_x %.2f, %.2f; interface %.2f, %.2f" % (o5[0], o5[1], orders(fif)[0], orders(fif)[1])])
val = rd("fe_run_validity.csv"); byid = {r["run"]: r for r in val}
prod = [r for r in val if not r["run"].startswith(("A1", "A2"))]
cmax = max(max(num(r["closure_iface_max"]), num(r["closure_global_max"])) for r in val if not r["run"].startswith("A2"))
sc.append(["V1", "every run completes, every SNES step converged", "met", "%d runs, all reached_end, none failed, min SNES reason %s" % (len(val), min(int(r["snes_min_reason"]) for r in val))])
sc.append(["V2", "steady drift < 1e-5 over the last 0.5", "met", "max drift %.2e" % max(max(num(r["drift_" + k]) for k in OBS) for r in fe.values())])
sc.append(["V6", "closure <= 1e-6 N/m", "met (internal consistency of the weak form)", "max %.2e over all non-A2 runs; A2 (from step 2) max %.2e" % (cmax, max(max(num(r["closure_iface_max"]), num(r["closure_global_max"])) for r in val if r["run"].startswith("A2")))])
sc.append(["V7", "no inverted cell, solid J > 0.2", "met (fluid ALE warning at x4)", "solid J min %.3f (x4), fluid J min %.4f (x4)" % (num(byid["A1S4"]["Jmin_solid"]), num(byid["A1S4"]["Jmin_fluid"]))])
wr("preregistered_scorecard.csv", ["id", "prediction_or_criterion", "outcome", "values"], sc)
# 8. ALE and mesh quality per level (A1 repeats carry the extended logging)
wr("ale_quality.csv", ["level", "Jmin_solid", "Jmin_fluid", "scaled_jac_min_solid", "scaled_jac_min_fluid", "area_ratio_min_solid", "area_ratio_min_fluid"],
   [[l, byid["A1S" + l]["Jmin_solid"], byid["A1S" + l]["Jmin_fluid"], byid["A1S" + l]["sj_solid_min"], byid["A1S" + l]["sj_fluid_min"], byid["A1S" + l]["area_ratio_solid_min"], byid["A1S" + l]["area_ratio_fluid_min"]] for l in LV])
# 8b. share of the interface force carried by the two corner nodes (direct fluid-side integral vs weak-form force without corner dofs)
wr("interface_force_weak_vs_direct.csv", ["level", "Fx_direct", "Fx_weak", "corner_share_Fx_pct", "Fy_direct", "Fy_weak", "corner_share_Fy_pct"],
   [[l, g(num(fe[l]["Fx"]), 8), g(num(fe[l]["Fx_weak"]), 8), g(100 * (1 - num(fe[l]["Fx_weak"]) / num(fe[l]["Fx"])), 3), g(num(fe[l]["Fy"]), 8), g(num(fe[l]["Fy_weak"]), 8), g(100 * (1 - num(fe[l]["Fy_weak"]) / num(fe[l]["Fy"])), 3)] for l in LV])
# 9. cost
led = rd("cpu_ledger_sacct.csv"); cl = {}
for r in led: cl.setdefault(r["class"], [0, 0.0]); cl[r["class"]][0] += 1; cl[r["class"]][1] += num(r["cpu_h"])
ev = []
for r in led:
    s = dt.datetime.fromisoformat(r["start"]); h_, m_, s_ = map(int, r["elapsed"].split(":")); e = s + dt.timedelta(hours=h_, minutes=m_, seconds=s_)
    ev += [(s, int(r["cpus"])), (e, -int(r["cpus"]))]
ev.sort(key=lambda a: (a[0], a[1])); cur = pk = 0
for _, d_ in ev: cur += d_; pk = max(pk, cur)
canc = sum(num(r["cpu_h"]) for r in led if r["state"].startswith("CANCELLED"))
wr("cost_summary.csv", ["class", "jobs", "cpu_h"], [[k, v[0], g(v[1], 4)] for k, v in sorted(cl.items())] + [["TOTAL", len(led), g(sum(v[1] for v in cl.values()), 4)], ["of which cancelled jobs", sum(1 for r in led if r["state"].startswith("CANCELLED")), g(canc, 3)], ["peak concurrent cores (from start/elapsed)", "", str(pk)]])
# 10. preregistration file hashes
pdir = os.path.join(HERE, "..", "preregistration"); rows = []
for n in ("PREREG_NH_BLOB", "PREREG_NH_BLOB_ADDENDUM1", "PREREG_NH_BLOB_ADDENDUM2"):
    h = hashlib.sha256(open(os.path.join(pdir, n + ".md"), "rb").read()).hexdigest(); rec = open(os.path.join(pdir, n + ".sha256")).read().split()[0]
    rows.append([n + ".md", h, "equal" if h == rec else "DIFFERENT"])
wr("preregistration_hash_check.csv", ["file", "sha256_computed", "matches_recorded_hash"], rows)
print("rerun vs PR text, max relative deviation of the printed values: %.2e" % dev)
print("PR-text u_x sequence: p=%.3f limit=%.6f; FE limit vs it %+.3f %%" % (a_pr["p"], a_pr["f_extrap"], fl))
for k in OBS: print(k, "FE", seqs[k][0]["triples"][0]["p"], seqs[k][0]["p"], "s4f", seqs[k][1]["triples"][0]["p"], seqs[k][1]["p"], "fixed-p diffs (p=1,1.5,2):", sens[k])
print("orders FE ux", orders(fux), "iface", orders(fif), "s4f ux", orders(sux), "iface", orders(sif))
print("cost", {k: (v[0], round(v[1], 3)) for k, v in cl.items()}, "total", round(sum(v[1] for v in cl.values()), 3), "peak", pk, "cancelled", round(canc, 3))
for k, v in sg.items(): print("solid", k, [round(x, 2) for x in v[0]], [round(x, 2) for x in v[1]], [round(x, 2) for x in v[2]], v[3])
