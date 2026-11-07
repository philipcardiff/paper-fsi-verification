"""Extract the small input tables of this package from the raw solver output (needs the solver repository worktrees; NOT needed to regenerate the derived CSVs).
usage: python3 -B extract_from_raw.py COUPLED_DIR C1_DIR OUT_DATA_DIR
  COUPLED_DIR: .../benchmarks/cavityFlexibleBottom/coupled   (prod_*/history.csv, summary.json, scan_*, rigid_*)
  C1_DIR     : .../benchmarks/cavityFlexibleBottom/structural_fluid  (c1_c2_gates.json)
Writes coupled_levels.csv, coupled_monitors.csv, pseudo_time_tail.csv, m1_vs_labelled_sensitivity.csv, time_solver_rigid_checks.csv,
plate_only_ladder.csv, fluid_only_ladder.csv.  Values are copied with full precision (repr)."""
import sys, os, csv, json, math

cdir, c1dir, out = sys.argv[1:4]


def hist(d):
    return list(csv.DictReader(open(os.path.join(cdir, d, "history.csv"))))


def summ(d):
    try:
        return json.load(open(os.path.join(cdir, d, "summary.json")))
    except Exception:
        return None


def f(x): return float(x)


def writecsv(name, header, rows):
    with open(os.path.join(out, name), "w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n"); w.writerow(header)
        for r in rows: w.writerow([repr(v) if isinstance(v, float) else v for v in r])


LV = [("L1", 0.1, "prod_P1_L1_dt10"), ("L2", 0.05, "prod_P2_L2_dt10"), ("L3", 0.025, "prod_P3_L3_dt10"), ("L4", 0.0125, "prod_P4_L4_dt10")]
rows, mon, tail = [], [], []
for lv, dx, d in LV:
    h, s = hist(d), summ(d); last = h[-1]
    rows.append([lv, dx, s["config"]["n_fluid_cells"], s["config"]["n_solid_cells"], s["config"]["ndofs"], s["cores"], s["steps"],
                 sum(int(float(x["snes_its"])) for x in h)] + [f(last[k]) for k in ("uy_A", "ux_A", "Fy", "Fy_p", "Fy_v", "Fx", "Fy_weak_all", "RLy", "RRy", "p_floor_mean", "V_fluid")])
    ipk = max(range(len(h)), key=lambda i: abs(f(h[i]["mass_defect"])))
    mon.append([lv, max(f(x["ifvel_err"]) for x in h), max(abs(f(x["P_f"]) - f(x["P_s"])) for x in h), max(f(x["balance_rel"]) for x in h),
                min(f(x["Jmin_f"]) for x in h), min(f(x["Jmin_s"]) for x in h),
                max(abs(f(x["mass_defect"]) - f(x["geom_closed_form"])) for x in h),
                max(abs(f(x["mass_defect"])) for x in h), h[ipk]["step"], abs(f(h[ipk]["Q_in"])),
                max(f(x["mass_norm"]) for x in h), abs(f(last["mass_norm"])), f(last["rate_uy"]), f(last["rate_Fy"])])
    if lv in ("L2", "L4"):
        for i in range(len(h) - 8, len(h)):
            tail.append([lv, int(h[i]["step"]), f(h[i]["t"]), f(h[i]["uy_A"]), f(h[i]["uy_A"]) - f(h[i - 1]["uy_A"]), f(h[i]["rate_uy"])])
writecsv("coupled_levels.csv", ["level", "dx_m", "fluid_tri", "plate_tri", "ndofs", "ranks", "steps", "newton_its_total", "uy_A_m", "ux_A_m", "Fy_N_per_m", "Fy_pressure", "Fy_viscous",
                                "Fx", "Fy_weak_all_dofs", "R_left_y", "R_right_y", "p_floor_mean", "V_fluid"], rows)
writecsv("coupled_monitors.csv", ["level", "max_ifvel_err", "max_abs_Pf_minus_Ps", "max_balance_rel", "min_Jmin_fluid", "min_Jmin_plate", "max_abs_mass_minus_closed_form",
                                  "peak_abs_mass_defect", "peak_abs_step", "abs_Qin_at_peak", "peak_mass_norm_instantaneous_Qin", "final_mass_norm", "final_rate_uy_per_s", "final_rate_Fy_per_s"], mon)
writecsv("pseudo_time_tail.csv", ["level", "step", "t_s", "uy_A_m", "increment_uy_m", "rate_uy_per_s"], tail)

# M1 vs labelled sensitivity (L2, same mesh and dt)
a, b = hist("prod_P2_L2_dt10"), hist("prod_M2_L2_dt10")
la, lb = a[-1], b[-1]
srow = []
for k in ("uy_A", "ux_A", "Fy", "Fy_weak_all", "RLy", "RRy", "balance_rel"):
    srow.append([k, f(la[k]), f(lb[k]), "" if k == "balance_rel" else (f(lb[k]) - f(la[k])) / abs(f(la[k])), max(f(x[k]) for x in a) if k == "balance_rel" else "", max(f(x[k]) for x in b) if k == "balance_rel" else ""])
srow.append(["max_abs_Pf_minus_Ps", max(abs(f(x["P_f"]) - f(x["P_s"])) for x in a), max(abs(f(x["P_f"]) - f(x["P_s"])) for x in b), "", "", ""])
writecsv("m1_vs_labelled_sensitivity.csv", ["quantity", "M1_L2", "M2_sensitivity_L2", "relative_change_(M2-M1)/|M1|", "M1_max_over_steps", "M2_max_over_steps"], srow)

# time-step, solver tolerance, rigid-wall and scan checks
chk = []
def last(d): return hist(d)[-1]
p1, p2 = last("prod_P1_L1_dt10"), last("prod_P2_L2_dt10")
t1, s1 = last("prod_T1_L2_dt5"), last("prod_S1_L1_tol")
for name, run, ref, lab in (("T1 L2 dt=5 vs dt=10", t1, p2, "prod_T1_L2_dt5"), ("S1 L1 tight SNES tol vs default", s1, p1, "prod_S1_L1_tol")):
    for k in ("uy_A", "Fy", "Fx", "ux_A"):
        chk.append([name, k, f(run[k]), f(ref[k]), abs(f(run[k]) - f(ref[k])) / abs(f(ref[k])), summ(lab)["steps"], "converged"])
for dtn in ("2.5", "5", "20"):
    d = "scan_L1_dt" + dtn; h = hist(d); s = summ(d)
    if dtn == "20":
        chk.append(["scan L1 dt=20 (diverged, kept)", "uy_A", "nan", f(p1["uy_A"]), "", len(h), "SNES diverged at step 1"]); continue
    for k in ("uy_A", "Fy"):
        chk.append(["scan L1 dt=" + dtn + " vs production dt=10", k, f(h[-1][k]), f(p1[k]), abs(f(h[-1][k]) - f(p1[k])) / abs(f(p1[k])), s["steps"], "converged"])
rf, r5 = hist("rigid_L1_fixed")[-1], hist("rigid_L1_stiff1e5")[-1]
chk.append(["rigid limit: E x1e5 vs frozen rigid floor", "Fy", f(r5["Fy"]), f(rf["Fy"]), abs(f(r5["Fy"]) - f(rf["Fy"])) / abs(f(rf["Fy"])), 40, "converged"])
chk.append(["rigid limit: E x1e5 vs frozen rigid floor", "p_floor_mean", f(r5["p_floor_mean"]), f(rf["p_floor_mean"]), abs(f(r5["p_floor_mean"]) - f(rf["p_floor_mean"])) / abs(f(rf["p_floor_mean"])), 40, "converged"])
hs, hf = hist("rigid_L1_stiff"), hist("rigid_L1_fixed")
ok = [i for i, x in enumerate(hs) if x["Fy"] not in ("nan", "")]
i14 = ok[-1]
chk.append(["rigid limit: E x1e7 at last converged step %d vs frozen at same step" % int(hs[i14]["step"]), "Fy", f(hs[i14]["Fy"]), f(hf[i14]["Fy"]), abs(f(hs[i14]["Fy"]) - f(hf[i14]["Fy"])) / abs(f(hf[i14]["Fy"])), int(hs[i14]["step"]), "then SNES failed (reason -5) at step %d; kept" % (int(hs[i14]["step"]) + 1)])
writecsv("time_solver_rigid_checks.csv", ["check", "quantity", "value_tested", "value_reference", "relative_difference", "steps", "status"], chk)

# plate-only and fluid-only ladders (preflight C1/C2)
g = json.load(open(os.path.join(c1dir, "c1_c2_gates.json")))
pr = []
for key, vals in g["ladder_uyA"].items():
    for nt, v, ux in zip((2, 4, 8, 16), vals, g["ladder_ux_A"][key]): pr.append([key, nt, v, ux])
writecsv("plate_only_ladder.csv", ["series", "elements_through_thickness", "uy_A_m", "ux_A_m"], pr)
fr = []
for fl, tab in g["fluid_table"].items():
    for row in tab:
        fr.append([fl] + row[:4])
writecsv("fluid_only_ladder.csv", ["floor", "cells", "Fy_N_per_m", "Fx_N_per_m", "p_floor_mean"], fr)
print("ok")
