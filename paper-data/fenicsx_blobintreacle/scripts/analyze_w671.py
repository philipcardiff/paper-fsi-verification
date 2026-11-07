"""Analysis of the solids4foam #671 A/B data (proposed branch head da1044a4f, unmerged; NOT an official release).
Usage: python3 -B scripts/analyze_w671.py     (standard library; needs data/s4f_w671_levels.csv, data/steady_convergence.csv)
Writes data/s4f_w671_convergence.csv, s4f_w671_effect.csv, s4f_w671_crosscode.csv, s4f_w671_ratios.csv
Conventions: effect = (CORRECTED - LEGACY)/|LEGACY| per level, 'mag' = (|C|-|L|)/|L| for negative quantities;
limit-effect resolved only if |dL| > 2u, u = sqrt(GCI_C^2 + GCI_L^2) (the campaign's preregistered rule B).
gap = s4f/FEniCSx - 1 (the campaign's convention in WALL671_RESULTS.md; the supplement's Table S2 uses the
opposite sign, FEniCSx/s4f - 1)."""
import csv, math, os, sys
sys.path.insert(0, os.path.dirname(__file__))
from convergence import analyse
here = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
R = list(csv.DictReader(open(os.path.join(here, "s4f_w671_levels.csv"))))
S = {}
for r in R:
    if r["state"] == "steady_t4":
        S.setdefault((r["model"], r["arm"], r["observable"]), []).append((r["level"], int(r["cells_total"]), float(r["value"])))
res = {}; out = []
for k, s in S.items():
    s.sort(key=lambda x: x[1]); f = [v for *_, v in s]; h = [c ** -0.5 for _, c, _ in s]
    a = analyse(f, h, "/".join(k))
    ok = a.get("accepted")
    if ok: res[k] = a
    out.append([k[0], k[1], k[2], "/".join(l for l, _, _ in s[-3:]), "%.3f" % a["p"] if "p" in a else "", "%.7g" % a["f_extrap"] if ok else "", "%.4f" % (100 * a["gci_fine"]) if ok else "", "accepted" if ok else "not accepted"])
with open(os.path.join(here, "s4f_w671_convergence.csv"), "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["model", "arm", "observable", "levels(finest3)", "p", "limit", "GCI_percent", "status"]); w.writerows(out)
def val(m, arm, o, lev): return next(v for l, c, v in S[(m, arm, o)] if l == lev)
eff = []
for m in ("orig", "earea"):
    for o in ("apex_ux", "apex_uy", "Fx_fluid", "Fy_fluid", "dp"):
        row = [m, o]
        for lev in ("x0.5", "x1", "x2", "x4"):
            try: row.append("%+.2f" % (100 * (abs(val(m, "CORRECTED", o, lev)) - abs(val(m, "LEGACY", o, lev))) / abs(val(m, "LEGACY", o, lev))))
            except StopIteration: row.append("")
        c, l = res.get((m, "CORRECTED", o)), res.get((m, "LEGACY", o))
        if c and l:
            d = (abs(c["f_extrap"]) - abs(l["f_extrap"])) / abs(l["f_extrap"]); u = math.hypot(c["gci_fine"], l["gci_fine"])
            row += ["%+.3f" % (100 * d), "%.3f" % (100 * u), "yes" if abs(d) > 2 * u else "no"]
        else: row += ["", "", ""]
        eff.append(row)
with open(os.path.join(here, "s4f_w671_effect.csv"), "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["model", "observable", "mag_change_x0.5_pct", "x1", "x2", "x4", "limit_change_pct", "u_pct", "resolved(|d|>2u)"]); w.writerows(eff)
# cross-code against FEniCSx M2 (steady_convergence.csv limits and GCI)
FE = {}
for r in csv.DictReader(open(os.path.join(here, "steady_convergence.csv"))):
    if r["model"] == "fenicsx_M2_BE" and r["triple"] == "finest3" and r["richardson_limit"]: FE[r["observable"]] = (float(r["richardson_limit"]), float(r["GCI_fine_percent"]) / 100)
FEv = {}
for r in csv.DictReader(open(os.path.join(here, "steady_levels.csv"))):
    if r["model"] == "fenicsx_M2_BE": FEv[(r["observable"], r["level"])] = float(r["value"])
rows = []
for o, fo in (("apex_ux", "apex_ux"), ("apex_uy", "apex_uy"), ("Fx_fluid", "Fx_fluid"), ("dp", "dp")):
    for arm in ("LEGACY", "CORRECTED"):
        a = res[("orig", arm, o)]; fe, gfe = FE[fo]
        g = a["f_extrap"] / fe - 1; u = math.hypot(a["gci_fine"], gfe)
        rows.append([o, arm, "limits", "%.7g" % a["f_extrap"], "%.7g" % fe, "%+.3f" % (100 * g), "%.3f" % (100 * u), "%.2f" % (abs(g) / u), "within combined GCI" if abs(g) <= u else "outside combined GCI (x1)" if abs(g) <= 2 * u else "outside 2u"])
        for lev, flev in (("x1", "1"), ("x2", "2"), ("x4", "4")):
            s4 = val("orig", arm, o, lev); f = FEv[(o, flev)]
            rows.append([o, arm, lev, "%.7g" % s4, "%.7g" % f, "%+.3f" % (100 * (s4 / f - 1)), "", "", ""])
with open(os.path.join(here, "s4f_w671_crosscode.csv"), "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["observable", "arm", "level_or_limits", "s4f", "fenicsx_M2", "gap_s4f_over_fenicsx_minus_1_pct", "u_combined_GCI_pct", "gap_over_u", "reading"]); w.writerows(rows)
# E-AREA/M2 like-for-like ratios, coarse triple x0.5,x1,x2 in both
rat = []
def lim3(m, arm, o):
    s = sorted(S[(m, arm, o)], key=lambda x: x[1])[:3]
    a = analyse([v for *_, v in s], [c ** -0.5 for _, c, _ in s], "")
    return a["f_extrap"] if a.get("accepted") else None
for arm in ("LEGACY", "CORRECTED"):
    for o in ("apex_ux", "apex_uy"):
        r3 = lim3("earea", arm, o) / lim3("orig", arm, o)
        r1 = val("earea", arm, o, "x1") / val("orig", arm, o, "x1"); r2 = val("earea", arm, o, "x2") / val("orig", arm, o, "x2")
        rat.append([arm, o, "%.5f" % r3, "%.5f" % r1, "%.5f" % r2])
with open(os.path.join(here, "s4f_w671_ratios.csv"), "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["arm", "observable", "E-AREA/M2 like-for-like limits (x0.5,x1,x2)", "level x1", "level x2"]); w.writerows(rat)
