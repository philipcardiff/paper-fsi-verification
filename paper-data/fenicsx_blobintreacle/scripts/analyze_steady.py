"""Steady-state spatial convergence, cross-code and model-form tables for the
blobInTreacle FEniCSx evidence.

Usage (from this directory):  python3 scripts/analyze_steady.py
Reads  data/steady_levels.csv   (hand-transcribed from the signed-off run records,
                                 one row per model, level, observable)
Writes data/steady_convergence.csv   observed order, Richardson limit, GCI
       data/steady_crosscode.csv     FEniCSx M2 against solids4foam, same level and limits
       data/modelform_M1_M2.csv      M1/M2 ratios at each level and for the limits
       data/fourway_M1_M2.csv        s4f M2, FEniCSx M2, s4f E-AREA, FEniCSx M1 (finest level and limits)
h = (total cells)^(-1/2).  The analysis is verification/convergence.py of the
FEniCSx campaign, copied unchanged to scripts/convergence.py (finest-three and
coarsest-three triples; accepted only if monotone and 0.5 <= p <= 3; Fs = 1.25).
Needs only the Python standard library."""
import csv, os, sys
sys.path.insert(0, os.path.dirname(__file__))
from convergence import analyse

here = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
rows = list(csv.DictReader(open(os.path.join(here, "steady_levels.csv"))))
series = {}
for r in rows:
    series.setdefault((r["model"], r["observable"]), []).append((r["level"], int(r["cells_total"]), float(r["value"])))

out, res = [], {}
for (model, obs), s in series.items():
    s.sort(key=lambda x: x[1])
    f = [v for _, _, v in s]; h = [c ** -0.5 for _, c, _ in s]
    for name, sl in (("finest3", slice(-3, None)), ("coarsest3", slice(0, 3))):
        if name == "coarsest3" and len(f) < 4:
            continue
        a = analyse(f[sl], h[sl], "%s %s %s" % (model, obs, name))
        levels = "/".join(l for l, _, _ in s[sl])
        if a.get("accepted"):
            res[(model, obs, name)] = a
            out.append([model, obs, name, levels, "%.3f" % a["p"], "%.7g" % a["f_extrap"], "%.4f" % (100 * a["gci_fine"]), "%.3f" % a["asymptotic_ratio"], "accepted"])
        else:
            out.append([model, obs, name, levels, "%.3f" % a["triples"][-1]["p"], "", "", "", "not accepted: " + a.get("reason", "")])
with open(os.path.join(here, "steady_convergence.csv"), "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["model", "observable", "triple", "levels", "p_observed", "richardson_limit", "GCI_fine_percent", "asymptotic_ratio", "status"])
    w.writerows(out)

def val(model, obs, level):
    for l, c, v in series[(model, obs)]:
        if l == level: return v

# cross-code: FEniCSx M2 (finest level F) against solids4foam (same level label)
obs_x = ["apex_ux", "apex_uy", "Fx_fluid", "Fy_fluid", "Rx_base", "Ry_base", "dp"]
with open(os.path.join(here, "steady_crosscode.csv"), "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["row", "observable", "level_or_triple", "s4f_M2", "fenicsx_M2", "fenicsx_over_s4f_minus_1_percent", "s4f_p", "s4f_GCI_percent", "fenicsx_p", "fenicsx_GCI_percent", "combined_GCI_percent"])
    for o in obs_x:
        for L in ["0.5", "1", "2", "4"]:
            a, b = val("s4f_M2_original", o, L), val("fenicsx_M2_BE", o, L)
            w.writerow(["same_level", o, "x%s/F%s" % (L, L), "%.7g" % a, "%.7g" % b, "%.3f" % (100 * (b / a - 1)), "", "", "", "", ""])
        a, b = res.get(("s4f_M2_original", o, "finest3")), res.get(("fenicsx_M2_BE", o, "finest3"))
        if a and b:
            w.writerow(["extrapolated", o, "finest3 (1,2,4)", "%.7g" % a["f_extrap"], "%.7g" % b["f_extrap"], "%.3f" % (100 * (b["f_extrap"] / a["f_extrap"] - 1)),
                        "%.2f" % a["p"], "%.3f" % (100 * a["gci_fine"]), "%.2f" % b["p"], "%.4f" % (100 * b["gci_fine"]), "%.3f" % (100 * (a["gci_fine"] ** 2 + b["gci_fine"] ** 2) ** 0.5)])

# model form: FEniCSx M1 / M2 at each level and for the limits
with open(os.path.join(here, "modelform_M1_M2.csv"), "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["observable", "level_or_limit", "fenicsx_M2", "fenicsx_M1", "M1_over_M2", "s4f_EAREA_over_s4f_M2"])
    for o in obs_x:
        for L in ["0.5", "1", "2", "4"]:
            a, b = val("fenicsx_M2_BE", o, L), val("fenicsx_M1_BE", o, L)
            e = ""
            if L in ("0.5", "1", "2") and ("s4f_EAREA_diagnostic", o) in series:
                e = "%.4f" % (val("s4f_EAREA_diagnostic", o, L) / val("s4f_M2_original", o, L))
            w.writerow([o, "F%s" % L, "%.7g" % a, "%.7g" % b, "%.4f" % (b / a), e])
        a, b = res.get(("fenicsx_M2_BE", o, "finest3")), res.get(("fenicsx_M1_BE", o, "finest3"))
        if a and b:
            w.writerow([o, "limit (finest3)", "%.7g" % a["f_extrap"], "%.7g" % b["f_extrap"], "%.4f" % (b["f_extrap"] / a["f_extrap"]), ""])
# four-way table: s4f M2 (original), FEniCSx M2, s4f E-AREA (modified coupling, same-lineage diagnostic), FEniCSx M1
cols = [("s4f_M2_original", "x4", "s4f M2 (original)"), ("fenicsx_M2_BE", "4", "FEniCSx M2"),
        ("s4f_EAREA_diagnostic", "2", "s4f E-AREA (diagnostic)"), ("fenicsx_M1_BE", "4", "FEniCSx M1")]
with open(os.path.join(here, "fourway_M1_M2.csv"), "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["observable"] + ["%s finest level (%s%s)" % (c[2], "x" if c[0].startswith("s4f") else "F", c[1].replace("x", "")) for c in cols] + ["%s limit" % c[2] for c in cols])
    for o in obs_x:
        row = [o]
        for m, L, _ in cols:
            v = val(m, o, L.replace("x", "")) if (m, o) in series else None
            row.append("" if v is None else "%.7g" % v)
        for m, L, _ in cols:
            a = res.get((m, o, "finest3"))
            row.append("" if a is None else "%.7g" % a["f_extrap"])
        w.writerow(row)
# ratio table: every M1/M2-type ratio with its exact numerator, denominator and levels
def lim(model, obs, trip):
    a = res.get((model, obs, trip))
    return None if a is None else a["f_extrap"]
tab = []
def add(label, num, den, levels, use, getter):
    r = []
    for o in ("apex_ux", "apex_uy", "Fx_fluid", "Fy_fluid"):
        a, b = getter(num, o), getter(den, o)
        r.append("" if a is None or b is None else "%.4f" % (a / b))
    tab.append([label, num[0], den[0], levels, use] + r)
L = lambda trip: (lambda m, o: lim(m[0] if isinstance(m, tuple) else m, o, trip))
add("FEniCSx M1 / M2, extrapolated limits", ("fenicsx_M1_BE",), ("fenicsx_M2_BE",), "F1,F2,F4 in both", "preferred FEniCSx ratio", L("finest3"))
add("FEniCSx M1 / M2, extrapolated limits", ("fenicsx_M1_BE",), ("fenicsx_M2_BE",), "F0.5,F1,F2 in both", "coarse-triple check", L("coarsest3"))
add("solids4foam E-AREA / M2, like-for-like", ("s4f_EAREA_diagnostic",), ("s4f_M2_original",), "x0.5,x1,x2 in both (E-AREA has no x4)", "preferred s4f ratio", lambda m, o: lim(m[0], o, "finest3") if m[0].startswith("s4f_EAREA") else lim(m[0], o, "coarsest3"))
add("solids4foam E-AREA / M2, MIXED level", ("s4f_EAREA_diagnostic",), ("s4f_M2_original",), "E-AREA x0.5,x1,x2 over M2 x1,x2,x4", "NOT to be used", lambda m, o: lim(m[0], o, "finest3"))
with open(os.path.join(here, "modelform_ratio_table.csv"), "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["ratio", "numerator", "denominator", "levels", "use", "apex_ux", "apex_uy", "Fx_fluid", "Fy_fluid"])
    w.writerows(tab)
# base-reaction ratios: why M1 Rx converges slowly while M2 Rx/Fx is constant
with open(os.path.join(here, "base_reaction_ratio.csv"), "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["level", "M2_abs_Rx_over_Fx_fluid", "M1_abs_Rx_over_Fx_applied", "M1_Rx_over_M2_Rx"])
    for L_ in ["0.5", "1", "2", "4"]:
        a2 = abs(val("fenicsx_M2_BE", "Rx_base", L_)) / val("fenicsx_M2_BE", "Fx_fluid", L_)
        a1 = abs(val("fenicsx_M1_BE", "Rx_base", L_)) / val("fenicsx_M1_BE", "Fx_fluid", L_)
        w.writerow(["F" + L_, "%.4f" % a2, "%.4f" % a1, "%.4f" % (val("fenicsx_M1_BE", "Rx_base", L_) / val("fenicsx_M2_BE", "Rx_base", L_))])
print("written", os.path.abspath(here))
