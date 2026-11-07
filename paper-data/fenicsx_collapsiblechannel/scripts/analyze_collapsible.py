"""Regenerate the derived CSVs of this package from the input CSVs in ../data (standard library only).
usage (from this package directory):  python3 -B -I scripts/analyze_collapsible.py
Inputs : troughs_fe.csv, troughs_reference.csv, history_space_effects.csv, wall_static_p2.csv, cost_cpuh.csv, cells.csv
Outputs: factorial_V1.csv, cross_term_V1.csv, convergence_analyse.csv, time_limit_p_sensitivity.csv, H2_deficits.csv, successive_trough_differences.csv,
         same_dt_offset_to_oomph.csv, three_way_x10.csv, wall_static_p2_orders.csv, H1_prediction_bands.csv, cost_summary.csv
Acceptance rule: scripts/convergence.py (>= 3 levels, monotone with decreasing differences, 0.5 <= p <= 3, Fs = 1.25); a series that fails is reported as NOT ACCEPTED with no limit.
All values are % of P (P = 0.216 m) unless a column says otherwise."""
import csv, os, sys, math
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import convergence as C
D = os.path.join(HERE, "..", "data")
XS = ("7.5", "10", "12.5")


def rd(n): return list(csv.DictReader(open(os.path.join(D, n))))


def g(v, n=10): return ("%.*g" % (n, v)) if isinstance(v, float) else str(v)


def wr(n, header, rows):
    with open(os.path.join(D, n), "w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n"); w.writerow(header)
        for r in rows: w.writerow([g(v) for v in r])


tr = rd("troughs_fe.csv"); trr = rd("troughs_reference.csv")
V = {}
for r in tr:
    if r["variant"] == "baseline" and r["face"] == "bottom":
        V[(r["mesh"], r["dt_s"], g(float(r["x_m"])))] = (float(r["V1_pctP"]), float(r["V2_pctP"]), float(r["V3_pctP"]) if r["V3_pctP"] else None)
OM = {}
for r in trr:
    if r["series"].startswith("oomph_"): OM[(r["series"][6:], g(float(r["x_m"])))] = (float(r["V1_pctP"]), float(r["V2_pctP"]), float(r["V3_pctP"]))
    else: OM[(r["series"], g(float(r["x_m"])))] = (float(r["V1_pctP"]), float(r["V2_pctP"]), float(r["V3_pctP"]))
DTS = ("0.025", "0.0125", "0.00625", "0.003125"); DTREF = dict(zip(DTS, ("dt40", "dt80", "dt160", "dt320"))); SS = ("S1", "S2", "S3"); NY = {"S1": 8, "S2": 16, "S3": 32}

# 1. factorial table with same-dt oomph difference and difference to the oomph reference
rows = []
for x in XS:
    for S in SS:
        for dt in DTS:
            if (S, dt, x) in V:
                v = V[(S, dt, x)][0]; rows.append([x, S, dt, v, OM[(DTREF[dt], x)][0], v - OM[(DTREF[dt], x)][0], OM[("ref", x)][0], v - OM[("ref", x)][0]])
wr("factorial_V1.csv", ["x_m", "mesh", "dt_s", "FE_V1_pctP", "oomph_same_dt_V1_pctP", "FE_minus_oomph_same_dt_pctP", "oomph_reference_V1_pctP", "FE_minus_oomph_reference_pctP"], rows)

# 2. cross term (registered additivity test) at V1
rows = []
for x in XS:
    for Sa, Sb in (("S1", "S2"), ("S2", "S3")):
        for a, b in (("0.0125", "0.00625"), ("0.00625", "0.003125")):
            if all((S, d, x) in V for S in (Sa, Sb) for d in (a, b)):
                v = lambda S, d: V[(S, d, x)][0]
                sp1, sp2 = v(Sb, a) - v(Sa, a), v(Sb, b) - v(Sa, b); t1, t2 = v(Sa, b) - v(Sa, a), v(Sb, b) - v(Sb, a); Xc = sp1 - sp2
                thr = max(0.05, 0.25 * min(abs(sp2), abs(t2)))
                rows.append([x, Sa + "/" + Sb, a, b, sp1, sp2, t1, t2, Xc, thr, "additive" if abs(Xc) <= thr else "NOT additive", "yes" if thr > max(abs(sp1), abs(sp2)) else "no"])
wr("cross_term_V1.csv", ["x_m", "mesh_pair", "dt_a_s", "dt_b_s", "space_effect_at_dt_a", "space_effect_at_dt_b", "time_effect_lower_mesh", "time_effect_upper_mesh", "cross_term_X", "registered_threshold",
                         "verdict", "threshold_exceeds_both_space_effects"], rows)

# 3. convergence.analyse on V1, V2, V3 (time: h = dt; space: h = 1/ny)
rows = []


def an(label, f, h, kind, x, series, qty):
    a = C.analyse(f, h, label); t = a["triples"][-1]
    if a.get("accepted"):
        rows.append([kind, x, series, qty, "|".join(g(q, 8) for q in f), t["p"], t["status"], "accepted", a["f_extrap"], 100 * a["gci_fine"]])
    else:
        rows.append([kind, x, series, qty, "|".join(g(q, 8) for q in f), t["p"], t["status"], "NOT ACCEPTED (%s)" % a.get("reason"), "", ""])


for k, qty in ((0, "V1"), (1, "V2"), (2, "V3")):
    for x in XS:
        for S in SS:
            L = [dt for dt in DTS[1:] if (S, dt, x) in V]
            if len(L) >= 3: an("t", [V[(S, dt, x)][k] for dt in L], [float(dt) for dt in L], "time", x, S, qty)
        for dt in DTS[1:]:
            L = [S for S in SS if (S, dt, x) in V]
            if len(L) >= 3: an("s", [V[(S, dt, x)][k] for S in L], [1.0 / NY[S] for S in L], "space", x, "dt " + dt, qty)
wr("convergence_analyse.csv", ["direction", "x_m", "series", "quantity", "values_coarse_to_fine", "p_finest_triple", "monotonicity", "acceptance", "richardson_limit", "GCI_fine_percent"], rows)

# 4. p-sensitivity of the time limit (two finest dt, r = 2); band quoted with every limit
rows = []
for S in ("S2", "S3"):
    for x in ("10", "12.5"):
        f3, f2, f1 = [V[(S, dt, x)][0] for dt in DTS[1:]]
        pobs = C.observed_order([f3, f2, f1], [0.0125, 0.00625, 0.003125])[0]
        for lab, p in (("1.0", 1.0), ("observed", pobs), ("1.5", 1.5), ("2.0", 2.0)):
            rows.append([S, x, lab, p, C.richardson(f1, f2, p, 2.0)])
wr("time_limit_p_sensitivity.csv", ["mesh", "x_m", "assumed_order_label", "p", "limit_V1_pctP"], rows)

# 5. H2 deficits: registered tolerance, two sets of oomph deficits (as quoted in the preregistration / recomputed from the reference CSV), p-sensitivity of the FE limit
TOL = {"0.025": 0.2, "0.0125": 0.1, "0.00625": 0.05, "0.003125": 0.05}
PREREG_OOMPH = {"0.025": 1.40, "0.0125": 0.375, "0.00625": 0.125, "0.003125": 0.019}
rows = []
oref = OM[("ref", "10")][0]
for S in ("S2", "S3"):
    f3, f2, f1 = [V[(S, dt, "10")][0] for dt in DTS[1:]]
    pobs = C.observed_order([f3, f2, f1], [0.0125, 0.00625, 0.003125])[0]
    for lab, p in (("observed", pobs), ("1.0", 1.0), ("1.5", 1.5), ("2.0", 2.0)):
        lim = C.richardson(f1, f2, p, 2.0)
        for dt in DTS:
            if (S, dt, "10") not in V: continue
            fe_def = abs(V[(S, dt, "10")][0] - lim)
            for src, od in (("prereg_quoted", PREREG_OOMPH[dt]), ("recomputed_from_csv", abs(OM[(DTREF[dt], "10")][0] - oref))):
                diff = fe_def - od
                rows.append([S, lab, p, lim, dt, fe_def, src, od, diff, TOL[dt], "pass" if abs(diff) <= TOL[dt] else "FAIL"])
wr("H2_deficits.csv", ["mesh", "extrapolation_order_label", "p", "FE_time_limit_V1_pctP", "dt_s", "FE_deficit_pctP", "oomph_deficit_source", "oomph_deficit_pctP", "FE_minus_oomph_deficit",
                       "registered_tolerance", "outcome_as_written"], rows)

# 6. successive trough differences (model free)
rows = []
for S in SS:
    for x in XS:
        for a, b in zip(DTS[:-1], DTS[1:]):
            if (S, a, x) in V and (S, b, x) in V:
                fe = V[(S, b, x)][0] - V[(S, a, x)][0]; om = OM[(DTREF[b], x)][0] - OM[(DTREF[a], x)][0]
                rows.append([S, x, a, b, fe, om, fe - om])
wr("successive_trough_differences.csv", ["mesh", "x_m", "dt_from_s", "dt_to_s", "FE_V1_change_pctP", "oomph_V1_change_pctP", "FE_minus_oomph_pctP"], rows)

# 7. same-dt offset (FE minus oomph at the same dt) at x = 10 and all x: constant, limit difference
rows = []
for S in SS:
    for x in XS:
        for dt in DTS:
            if (S, dt, x) in V: rows.append([S, x, dt, V[(S, dt, x)][0] - OM[(DTREF[dt], x)][0]])
wr("same_dt_offset_to_oomph.csv", ["mesh", "x_m", "dt_s", "FE_minus_oomph_same_dt_V1_pctP"], rows)

# 8. three-way table at x = 10 (no composite limit)
rows = [["oomph reference (Richardson in time)", OM[("ref", "10")][0], "U_ref 0.30 % P; spatial bias 0.26 % P unremoved"],
        ["oomph dt160", OM[("dt160", "10")][0], "same-dt"], ["oomph dt320", OM[("dt320", "10")][0], "same-dt"]]
for S, dt in (("S2", "0.00625"), ("S3", "0.00625"), ("S3", "0.003125")): rows.append(["FE %s dt %s" % (S, dt), V[(S, dt, "10")][0], "finest S3 at dt 0.003125 covers t <= 1.25 s" if S == "S3" else ""])
for F in ("F1", "F2", "F3", "F4"): rows.append(["s4f %s (Robin, dt 0.00625; pre-#671 legacy run)" % F, OM[("s4f_" + F, "10")][0], "wallMid; solids4foam a421f063d, 2 Nov 2026"])
wr("three_way_x10.csv", ["item", "V1_pctP", "note"], rows)
# finest-cell minus oomph reference at all x
rows = []
for x in XS: rows.append([x, V[("S3", "0.003125", x)][0], OM[("ref", x)][0], V[("S3", "0.003125", x)][0] - OM[("ref", x)][0]])
wr("finest_cell_minus_oomph_reference.csv", ["x_m", "FE_S3_dt0.003125_V1_pctP", "oomph_reference_V1_pctP", "difference_pctP"], rows)

# 9. P2-triangle wall alone: order on 80/160/320 x 4
ws = {(r["nx"], r["ny"]): r for r in rd("wall_static_p2.csv")}
rows = []
for q, k in (("mid", "u_mid_m"), ("quarter", "u_quarter_m")):
    f = [float(ws[(n, "4")][k]) for n in ("80", "160", "320")]; a = C.analyse(f, [1 / 80, 1 / 160, 1 / 320], q)
    rows.append([q, "80x4|160x4|320x4", a["p"], a["triples"][-1]["status"], "accepted" if a.get("accepted") else "NOT ACCEPTED", a.get("f_extrap", ""), 100 * a["gci_fine"] if a.get("accepted") else ""])
wr("wall_static_p2_orders.csv", ["quantity", "levels", "p", "monotonicity", "acceptance", "richardson_limit_m", "GCI_fine_percent"], rows)

# 10. H1: registered ranges against the trough (V1) and against the history norm (L-inf, t <= 1.5, dt 0.00625)
he = {(r["dt_s"], g(float(r["x_m"]))): r for r in rd("history_space_effects.csv")}
rows = []
for x in XS:
    s12 = abs(V[("S2", "0.00625", x)][0] - V[("S1", "0.00625", x)][0]); s23 = abs(V[("S3", "0.00625", x)][0] - V[("S2", "0.00625", x)][0])
    h = he[("0.00625", x)]
    rows.append([x, s12, s23, float(h["Linf_S2_minus_S1_pctP"]), float(h["Linf_S3_minus_S2_pctP"]), "0.5-3 / 0.1-0.5",
                 "ratio of history L-inf changes %.1f" % (float(h["Linf_S2_minus_S1_pctP"]) / float(h["Linf_S3_minus_S2_pctP"]))])
wr("H1_prediction_bands.csv", ["x_m", "V1_change_S1_to_S2", "V1_change_S2_to_S3", "history_Linf_S1_to_S2", "history_Linf_S2_to_S3", "registered_ranges_pctP_(S1-S2 / S2-S3)", "note"], rows)

# 11. cost
cost = rd("cost_cpuh.csv"); cells = rd("cells.csv"); cs = {}
for r in cost: cs.setdefault(r["campaign"], [0.0, 0, 0]); cs[r["campaign"]][0] += int(r["cpu_s_sacct"]) / 3600.0; cs[r["campaign"]][1] += 1; cs[r["campaign"]][2] += (r["state"] != "COMPLETED")
rows = [[k, v[1], v[2], v[0]] for k, v in sorted(cs.items())]
prod = sum(int(r["cpu_s_sacct"]) for r in cost if r["campaign"] == "COL-P production") / 3600.0
mf1 = sum(int(r["cpu_s_sacct"]) for r in cost if r["campaign"] == "COL-P MF1") / 3600.0
aud = sum(int(r["cpu_s_sacct"]) for r in cost if r["campaign"] == "auditor") / 3600.0
rows += [["COL-P production + MF1", "", "", prod + mf1], ["COL-P production + MF1 + auditor", "", "", prod + mf1 + aud]]
wr("cost_summary.csv", ["campaign", "jobs", "jobs_not_COMPLETED", "cpu_h_sacct"], rows)
print("ok")
