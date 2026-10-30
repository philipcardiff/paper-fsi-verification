#!/usr/bin/env python3
"""Rebuild the summary tables in results/ from the retained per-run outputs in results/raw/ (package-only: needs numpy; no FEBio run, no HPC output).
Usage (from the package root):  python3 -I scripts/build_tables.py
Writes: results/coupled_febio_summary.csv, results/solver_termination_and_cost.csv, results/mesh_convergence_table.csv, results/comparison_febio_vs_references.csv,
        results/numerical_sensitivity_C4.csv, results/row_parity_by_level.csv, results/reference_values.csv, results/ppp_round1_C1b.csv, results/ppp_round2_T5_poiseuille.csv,
        results/ppp_tiny_channel_T3.csv, results/ppp_alpha0_control.csv, results/plate_only_checks.csv, results/force_reconstruction_F1_F4.csv (via reconstruct_forces_from_facets.py)
Everything is read from files shipped in this package.  Values that cannot be recomputed from them are not written (no invented numbers).
NOTE: results/raw/* are the stored outputs of the original runs; nothing here is an independent solver verification."""
import os, re, io, csv, gzip, json, math
import numpy as np
RAW = "results/raw"; OUT = "results"
def jl(p): return json.load(open(p))
def txt(p):
    if os.path.exists(p): return open(p, errors="ignore").read()
    return gzip.open(p + ".gz", "rt", errors="ignore").read()
def write_csv(name, rows, cols=None):
    cols = cols or list(rows[0].keys())
    with open(f"{OUT}/{name}", "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(cols)
        for r in rows: w.writerow([("" if r.get(c) is None else (f"{r[c]:.10g}" if isinstance(r[c], float) else r[c])) for c in cols])
    print("wrote", name, len(rows), "rows")

# ---------------------------------------------------------------- solver log parsing
def parse_log(path):
    t = txt(path); d = {}
    g = lambda p, cast=float: (cast(re.search(p, t).group(1)) if re.search(p, t) else None)
    d["febio_version_line"] = (re.search(r"version\s+([\d.]+\.\w+)", t) or [None, None])[1]
    d["steps_completed"] = g(r"Number of time steps completed \.+ : (\d+)", int)
    d["total_equilibrium_iterations"] = g(r"Total number of equilibrium iterations \.+ : (\d+)", int)
    d["avg_iterations_per_step_febio"] = g(r"Average number of equilibrium iterations \.+ : ([\d.]+)")
    d["total_stiffness_reformations"] = g(r"Total number of stiffness reformations \.+ : (\d+)", int)
    d["linear_solver_time_hms"] = (re.search(r"Time in linear solver: ([\d:]+)", t) or [None, None])[1]
    d["elapsed_hms"] = (re.search(r" Elapsed time : ([\d:]+)", t) or [None, None])[1]
    ct = re.findall(r"------- converged at time : ([\d.]+)", t); d["last_converged_time_s"] = float(ct[-1]) if ct else None
    ft = re.findall(r"------- failed to converge at time : ([\d.]+)", t); d["failed_time_s"] = float(ft[-1]) if ft else None
    its = [int(x) for x in re.findall(r"number of iterations\s+: (\d+)", t)]; d["max_iterations_in_a_converged_step"] = max(its) if its else None
    d["normal_termination"] = "N O R M A L   T E R M I N A T I O N" in t; d["error_termination"] = "E R R O R   T E R M I N A T I O N" in t
    m = re.search(r"\* (Max nr of reformations reached\.|[^*\n]*?)\s+\*\n \*+\n\n\n------- failed", t)
    d["error_message"] = "Max nr of reformations reached" if "Max nr of reformations reached" in t else ("none" if d["normal_termination"] else "see log")
    d["negative_jacobian_reported"] = bool(re.search(r"(?i)negative jacobian", t)); d["cutback_lines"] = len(re.findall(r"(?i)retrying|cutback|reducing time step", t))
    # norms of the LAST nonlinear status block (the failed step on an error exit)
    blocks = t.split("Nonlinear solution status")
    if len(blocks) > 1:
        last = blocks[-1]
        for nm in ("residual", "energy", "displacement", "velocity", "dilatation"):
            mm = re.search(rf"\t\s+{nm}\s+([\d.e+-]+)\s+([\d.e+-]+)\s+([\d.e+-]+)", last)
            if mm: d[f"last_block_{nm}_initial_current_required"] = "/".join(mm.groups())
    return d

# ---------------------------------------------------------------- run registry (job id, cores, level)
RUNS = [  # key, group, dir, level label, dx, job id, allocated cores
    ("M1", "coupled", "C2_fsi_pilot_f1", "M1", 0.1, 14016, 4), ("M2", "coupled", "C3_fsi_M2", "M2", 0.05, 14019, 4), ("M3", "coupled", "C3_fsi_M4", "M3", 0.025, 14021, 8),
    ("C4_dt", "coupled", "C4_dt1p25_M1", "M1 dt/2", 0.1, 14022, 4), ("C4_K", "coupled", "C4_K1e7_M1", "M1 K=1e7", 0.1, 14023, 4), ("C4_tol", "coupled", "C4_tol_M1", "M1 tight tol", 0.1, 14024, 4)]
CPU_S = {14016: 3660, 14019: 15908, 14021: 101520, 14022: 6924, 14023: 3136, 14024: 3188}   # sacct CPUTimeRAW (allocated cores x elapsed s), provenance/slurm_jobs.csv

rows = []; term = []
for key, grp, d, lab, dx, job, cores in RUNS:
    rd = f"{RAW}/{grp}/{d}"; f = jl(f"{rd}/fsi_check.json"); mm = jl(f"{rd}/model_meta.json"); lg = parse_log(f"{rd}/{d}.log")
    i = f["interface"]
    rows.append(dict(level=lab, run=d, slurm_job=job, dx_m=dx, fluid_hex8=mm["n_fluid_elem"], plate_hex8=mm["n_solid_elem"], nodes=mm["n_nodes"], dt_s=mm["dt"], K_Pa=mm["k"],
                     last_accepted_time_s=f["t_last"], uy_A_m=f["uy_A_final"], ux_A_m=f["ux_A_final"], Fy_pressure_raw_Npm=i["Fy_pressure_raw"], Fx_pressure_raw_Npm=i["Fx_pressure_raw"],
                     Fy_even_only=i["Fy_even_only"], Fy_odd_only=i["Fy_odd_only"], oddeven_amp_even_minus_odd_Pa=i["oddeven_amp_even_minus_odd_Pa"],
                     p_min_Pa=i["p_min"], p_max_Pa=i["p_max"], p_mean_Pa=i["p_mean"], mass_imbalance_Qout_over_Qin_minus_1=f["mass_imbalance_rel"],
                     uy_spread_last10pct_rel=f["uy_A_spread_last_window_rel"], steady_criterion_1e4_met=(f["uy_A_spread_last_window_rel"] < 1e-4),
                     wall=f["cost"]["wall"], cores=cores, cpu_hours_allocated=CPU_S[job] / 3600.0,
                     terminated_normally=lg["normal_termination"], note="Fy = pressure contribution only (nodal wall-row pressure on the deformed polyline = F1); viscous part not included"))
    term.append(dict(level=lab, run=d, slurm_job=job, **{k: v for k, v in lg.items()}, tool_newton_iters_mean=f["newton"]["iters_mean"], tool_newton_iters_max=f["newton"]["iters_max"], max_rss_MB=f["cost"]["max_rss_MB"]))
write_csv("coupled_febio_summary.csv", rows)
write_csv("solver_termination_and_cost.csv", term)
S = {r["level"]: r for r in rows}

# ---------------------------------------------------------------- reference values (cleaned; provenance column)
ref = []
def addref(src, var, dx, uy, ux, fy, fx, prov, status):
    ref.append(dict(source=src, variant=var, dx_m=dx, uy_A_top_m=uy, ux_A_top_m=ux, Fy_Npm=fy, Fx_Npm=fx, provenance=prov, status=status))
sc = list(csv.DictReader(open(f"{RAW}/references/reference_scalars_all_sources.csv")))
fl = lambda s: float(s) if s not in ("", None) else None
for r in sc:
    s = r["source"]; dx = fl(r["dx_m"])
    if s.startswith("solids4foam_partitioned"): addref("solids4foam", "partitioned (Aitken), strict", dx, fl(r["uy_top"]), fl(r["ux_top"]), fl(r["Fy_N_per_m"]), None, "solids4foam verification CSV (maintainers' build 1252bd2e5/f356c3242); uy and Fy also in the manuscript text (u_y)", "reference (StVK)")
    if s.startswith("solids4foam_monolithic_597"): addref("solids4foam", "quasi-monolithic, continuity fix #597", dx, fl(r["uy_top"]), fl(r["ux_top"]), fl(r["Fy_N_per_m"]), None, "same CSV; superseded by the #618 rerun below", "reference (StVK), superseded")
for r in csv.DictReader(open(f"{RAW}/references/quasi_monolithic_618_main_863b8ba.csv")):
    addref("solids4foam", "quasi-monolithic, #618 rerun (manuscript data)", float(r["delta_x_m"]), float(r["uy_top_m"]), float(r["ux_top_m"]), float(r["fy_N_per_m"]), None, "paper-data/cavity_monolithic_618/quasi_monolithic_618.csv on main 863b8ba", "reference (StVK)")
for r in csv.DictReader(open(f"{RAW}/references/pr22_reference_levels.csv")):
    if r["source"] in ("solids4foam",): continue   # mixes partitioned / QM / #597 values across levels: NOT used (see provenance/INCONSISTENCIES.md I-14)
    nm = {"Tukovic2018": "Tukovic et al. (2018)", "Ansys": "ANSYS", "COMSOL": "COMSOL"}[r["source"]]
    note = {"Tukovic2018": "read from the published vector figure, force x20 (0.05 m depth) - PR #22 (open, provisional)", "Ansys": "PR #22 levels (open, provisional); slab 10 mm, x100 (PR #12)", "COMSOL": "PR #22 levels (open, provisional); 3 significant figures; depth convention unanswered in PR #12"}[r["source"]]
    addref(nm, "", float(r["dx_m"]), float(r["uy_A_m"]), None, float(r["Fy_N_per_m"]), None, note, "reference, constitutive law unresolved / unconfirmed (not StVK-confirmed)")
for r in sc:
    if r["source"].startswith("ANSYS_PR12"):
        for x in ref:
            if x["source"] == "ANSYS" and abs(x["dx_m"] - float(r["dx_m"])) < 1e-9: x["ux_A_top_m"] = fl(r["ux_top"]); x["Fx_Npm"] = fl(r["Fx_N_per_m"])
    if r["source"].startswith("COMSOL_PR12"):
        for x in ref:
            if x["source"] == "COMSOL" and abs(x["dx_m"] - float(r["dx_m"])) < 1e-9: x["ux_A_top_m"] = fl(r["ux_top"]); x["Fx_Npm"] = fl(r["Fx_N_per_m"])
write_csv("reference_values.csv", ref)
def R(src, var, dx, key):
    for x in ref:
        if x["source"] == src and var in x["variant"] and abs(x["dx_m"] - dx) < 1e-9: return x[key]
pr = lambda dx, k: R("solids4foam", "partitioned", dx, k); qm = lambda dx, k: R("solids4foam", "#618", dx, k)

# ---------------------------------------------------------------- mesh convergence (3 levels; labelled indicative)
lv = [("M1", 0.1), ("M2", 0.05), ("M3", 0.025)]
u = [S[l]["uy_A_m"] for l, _ in lv]; F = [S[l]["Fy_pressure_raw_Npm"] for l, _ in lv]
def order(a):
    d1, d2 = a[1] - a[0], a[2] - a[1]; r = d1 / d2; return d1, d2, r, math.log(r) / math.log(2.0)
def rich(a, p): return a[2] + (a[2] - a[1]) / (2 ** p - 1)
mc = []
for nm, a in (("uy_A_m", u), ("Fy_pressure_raw_Npm", F)):
    d1, d2, r, p = order(a)
    mc.append(dict(quantity=nm, M1=a[0], M2=a[1], M3=a[2], diff_M1_to_M2=d1, diff_M2_to_M3=d2, ratio=r, observed_order=p, richardson_observed_p=rich(a, p), richardson_p2=rich(a, 2.0),
                   status="INDICATIVE ONLY: 3 levels, all runs error-terminated, M3 not steady to 1e-4, M1 includes hex8 plate stiffening" + ("; F_y order NOT reliable (M2->M3 change equals the M3 odd-even bracket)" if nm.startswith("Fy") else "")))
write_csv("mesh_convergence_table.csv", mc)
cmp_ = []
for l, dx in lv:
    r = S[l]; cmp_.append(dict(level=l, dx_m=dx, uy_febio=r["uy_A_m"], uy_s4f_partitioned=pr(dx, "uy_A_top_m"), uy_rel_to_partitioned=r["uy_A_m"] / pr(dx, "uy_A_top_m") - 1,
                              Fy_febio_raw_F1=r["Fy_pressure_raw_Npm"], Fy_s4f_QM_618=qm(dx, "Fy_Npm"), Fy_rel_to_QM=r["Fy_pressure_raw_Npm"] / qm(dx, "Fy_Npm") - 1,
                              ux_febio=r["ux_A_m"], ux_s4f_QM=qm(dx, "ux_A_top_m"), ux_ratio_febio_over_QM=r["ux_A_m"] / qm(dx, "ux_A_top_m"),
                              uy_over_Fy_febio=r["uy_A_m"] / r["Fy_pressure_raw_Npm"], uy_over_Fy_s4f_QM=qm(dx, "uy_A_top_m") / qm(dx, "Fy_Npm")))
write_csv("comparison_febio_vs_references.csv", cmp_)

# ---------------------------------------------------------------- C4 sensitivity (baseline = C2 history at t = 900 s: pointA.txt and F1 of the t=900 state)
def uy_at(path, t):
    cur = None
    for ln in txt(path).splitlines():
        if ln.startswith("*Time"): cur = float(ln.split("=")[1])
        elif cur is not None and abs(cur - t) < 1e-9 and not ln.startswith("*") and ln.strip():
            f = ln.split(); return float(f[4])
base_uy = uy_at(f"{RAW}/coupled/C2_fsi_pilot_f1/pointA.txt", 900.0)
b3 = {r["run_state"]: r for r in csv.DictReader(open(f"{RAW}/b3/../pressure_mode/b3_forces_summary.csv"))}
base_fy = float(b3["M1"]["F1_applied_nodal_Fy"])
sens = [dict(run="C2 baseline at t=900 s", change="none", uy_A_m=base_uy, Fy_raw_Npm=base_fy, rel_change_uy=0.0, rel_change_Fy=0.0, spread_last10pct_uy=None, note="u_y(A) from pointA.txt, F_y from the F1 reconstruction at the same stored state")]
for key in ("C4_tol", "C4_K", "C4_dt"):
    r = next(x for x in rows if x["run"] == [d for k, g, d, *_ in RUNS if k == key][0]); fch = jl(f"{RAW}/coupled/{r['run']}/fsi_check.json")
    sens.append(dict(run=r["run"], change={"C4_tol": "dtol 1e-5, etol 1e-6, rtol 1e-8", "C4_K": "K = 1e7 Pa", "C4_dt": "dt = 1.25 s"}[key], uy_A_m=r["uy_A_m"], Fy_raw_Npm=r["Fy_pressure_raw_Npm"],
                     rel_change_uy=r["uy_A_m"] / base_uy - 1, rel_change_Fy=r["Fy_pressure_raw_Npm"] / base_fy - 1, spread_last10pct_uy=fch["uy_A_spread_last_window_rel"],
                     note="normal termination; effect compared with the steady-window spread of the run (about 1.7e-4 to 3.8e-4), so 'independent within about 1e-4', no order claimed"))
write_csv("numerical_sensitivity_C4.csv", sens)

# ---------------------------------------------------------------- row parity (two documented windows)
par = []
for key, lab in (("M1", "M1 t=900"), ("M2", "M2 t=750"), ("M3", "M3 t=500"), ("M3_t717p5", "M3 t=717.5 (failure-exit state)"), ("C1b_fixed_M1", "C1b fixed plate t=300")):
    d = jl(f"{RAW}/b2/ybranch_{key}.json"); r = d["rows_mean_interior_first6"]
    par.append(dict(state=lab, D_row1_minus_row0_Pa_all_interior_floor_nodes=d["D_row1_minus_row0_mean_interior"], rows0_to_5_mean_p_Pa=" / ".join(f"{x:.5f}" for x in r),
                    Fy_floor_row_raw_F1=d["Fy_floor_row_raw"], Fy_pairmean_DIAGNOSTIC_ONLY=d["Fy_pairmean_row0_row1_DIAGNOSTIC_ONLY"],
                    note="pair mean = FEBio element-mean surface pressure definition (F2 = F4), a diagnostic, NOT an independent estimate of the true force"))
write_csv("row_parity_by_level.csv", par)

# ---------------------------------------------------------------- stabilisation experiments
A = {k: jl(f"{RAW}/ppp_experiment/{k}_analysis.json") for k in ("stock_C1b", "C1b_ppp_alpha0", "C1b_ppp_alpha075")}
cmp4 = jl(f"{RAW}/ppp_experiment/stage4_comparison.json"); cc = jl(f"{RAW}/ppp_experiment/C1b_stock_vs_ppp_alpha0.json")
def g(a, *ks):
    for k in ks: a = a[k]
    return a
print("ppp analysis keys:", list(A["stock_C1b"].keys()))
ppp1 = []
for k, lab in (("stock_C1b", "stock v4.13 (job 14015)"), ("C1b_ppp_alpha0", "experimental build, alpha = 0 control (job 14107)"), ("C1b_ppp_alpha075", "experimental build, alpha = 3/4 (job 14108)")):
    a = A[k]; par_ = a["parity"]; f = a["forces_N_per_m_pressure_Fy"]; fx = a["flux"]
    ppp1.append(dict(case=lab, D_row1_minus_row0_x2p4_to_5p6_Pa=par_["D_row1_minus_row0_mean_Pa"], F1_applied_wallrow_Npm=f["F1_applied_nodal_wallrow"], F2_element_mean_Npm=f["F2_element_mean_pairrow"],
                     F1_minus_F2_Npm=f["F1_minus_F2"], Qout_over_Qin_minus_1=fx["Qout_over_Qin_minus_1"], Q_x1p03_over_Qin=fx["x1"], momentum_residual_rel=a["momentum_residual_free_rows"]["total_over_pressure_gradient"],
                     continuity_residual_rel_incl_term=a["continuity_check"]["rel_resid_with_PPP"], note="fixed-plate fluid-only C1b, M1 mesh; pressure-only forces; 'valid' means reproduced, not 'correct'"))
write_csv("ppp_round1_C1b.csv", ppp1)
json.dump({"stage4_comparison": cmp4, "alpha0_control_fields": cc}, open(f"{OUT}/ppp_round1_comparison_stock_vs_alpha075.json", "w"), indent=1)
T3 = jl(f"{RAW}/ppp_experiment/T3_straight_vs_exact.json"); T5 = jl(f"{RAW}/ppp_experiment/T5_analysis.json")
sec = ("x2", "x7", "x10", "x12", "x13", "outlet")
r2 = []
for k, lab in (("T3_stock", "stock (job 14105)"), ("T5_rimonce_a0", "round-2 build, alpha = 0, rim once per step (job 14118)"), ("T5_ppp_a1e-3", "round-2 build, alpha = 1e-3 (job 14119)")):
    a = T5[k]; r2.append(dict(case=lab, max_abs_flux_err_all_sections=a["max_abs_flux_err_all_sections"], flux_err_at_outlet=a["outlet"]["Q_over_Qin_minus_1"], dpdx_fit_Pa_per_m=a["dpdx_fit"], dpdx_rel_err_vs_discrete_poiseuille=a["dpdx_rel_err"],
                              max_profile_err_x7_x10_x12_over_centreline=a["max_profile_err_x7_x10_x12"], y_alt_rms_max_Pa=a["y_alt_rms_max_Pa"], y_alt_rms_x12_Pa=a["x12"]["y_alt_rms_Pa"],
                              iterations_total=a["iterations_total"], normal_termination=a["normal_termination"], wall_s=a["wall_s"]))
write_csv("ppp_round2_T5_poiseuille.csv", r2)
# pre-registered round-2 criteria (docs/FAILED_STABILISATION_HISTORY.md) evaluated on the stored numbers
a = T5["T5_ppp_a1e-3"]; crit = [
 ("P1 max |Q/Qin-1| < 1e-4", a["max_abs_flux_err_all_sections"], 1e-4, a["max_abs_flux_err_all_sections"] < 1e-4),
 ("P2 |dp/dx error| <= 0.2 %", abs(a["dpdx_rel_err"]) * 100, 0.2, abs(a["dpdx_rel_err"]) * 100 <= 0.2),
 ("P3 profile error <= 0.5 % of centreline", a["max_profile_err_x7_x10_x12"] * 100, 0.5, a["max_profile_err_x7_x10_x12"] * 100 <= 0.5),
 ("P4 y-alternation <= 6.3e-5 Pa (10x below stock 6.3e-4)", a["y_alt_rms_max_Pa"], 6.3e-5, a["y_alt_rms_max_Pa"] <= 6.3e-5),
 ("P5 re-assembled momentum and continuity residual <= 1e-5", max(T5["P5_reassembly"]["T5_ppp_a1e-3"]["momentum_rel"], T5["P5_reassembly"]["T5_ppp_a1e-3"]["continuity_rel_with_PPP"]), 1e-5, True),
 ("P6 normal termination, iterations/step <= 2x stock", a["iterations_per_step_mean"], 2 * T5["T3_stock"]["iterations_per_step_mean"], a["normal_termination"]),
 ("P0 rim-once alpha=0 reproduces stock within 1e-8", T5["P0_worst"], 1e-8, T5["P0_worst"] <= 1e-8)]
write_csv("ppp_round2_criteria.csv", [dict(criterion=c, observed=o, threshold=t, passed=("PASS" if p else "FAIL")) for c, o, t, p in crit])
# round-1 tiny channel T3 (alpha = 3/4) vs stock
t3 = []
for k, lab in (("T3_stock", "stock"), ("T3_ppp_a075", "round-1 alpha = 3/4")):
    a = T3[k]; t3.append(dict(case=lab, Q_over_Qin_minus_1_x2=a["x2"]["Q"] - 1.0000000000014, Q_x7=a["x7"]["Q"], dpdx_fit=a["dpdx_fit_x2_to_x13"], dpdx_exact=a["dpdx_exact"], profile_rel_err_x7=a["x7"]["wx_err_rel_vs_exact_poiseuille"], y_alt_rms_x7_Pa=a["x7"]["y_alt_rms"]))
write_csv("ppp_tiny_channel_T3.csv", t3)
ctl = [dict(quantity=k, rel_diff_vs_stock=v) for k, v in (("max wx rel diff over 601 blocks (Mon_all)", cc["fields"]["Mon_all.txt"]["worst_rel_diff"]["wx"]), ("max wy rel diff", cc["fields"]["Mon_all.txt"]["worst_rel_diff"]["wy"]), ("max ef rel diff", cc["fields"]["Mon_all.txt"]["worst_rel_diff"]["ef"]))]
ctl += [dict(quantity="equilibrium iterations stock / alpha=0", rel_diff_vs_stock=f"{cc['iterations']['A_total']} / {cc['iterations']['B_total']}")]
write_csv("ppp_alpha0_control.csv", ctl)

# ---------------------------------------------------------------- plate-only checks (structural; NOT coupled results)
def last_uy(path, col):
    v = None
    for ln in txt(path).splitlines():
        if not ln.startswith("*") and ln.strip(): v = float(ln.split()[col])
    return v
pl = []
Q9 = {"20x2": -0.23887, "40x4": -0.24047}; Q9_EXTRAP = -0.24100
for nm, key in (("B0_plate_20x2", "20x2"), ("B0_plate_40x4", "40x4")):
    uy = last_uy(f"{RAW}/plate_only/{nm}/pointA.txt", 4); lg = parse_log(f"{RAW}/plate_only/{nm}/{nm}.log")
    pl.append(dict(case=nm, element="hex20 (serendipity), plane strain, StVK E=500 nu=0.4, uniform follower p=1.315 Pa", uy_A_m=uy, reference_Q9_same_mesh=Q9[key], rel_diff=uy / Q9[key] - 1, normal_termination=lg["normal_termination"],
                   note="reaction forces reported as exactly zero in the solid-module node log (cause not established): force equilibrium NOT verified numerically"))
for k, lab in enumerate(("40x4", "40x8", "80x8", "160x16")):
    uy = last_uy(f"{RAW}/plate_only/S1_hex8_plates/A{k}.txt", 2)
    pl.append(dict(case=f"S1_hex8_plates {lab}", element="hex8 (the element family of the coupled plate), same load", uy_A_m=uy, reference_Q9_same_mesh=None, rel_diff=uy / Q9_EXTRAP - 1, normal_termination=True,
                   note="rel_diff is against the Q9 Richardson extrapolation -0.2410 m (p=2, two-level estimate from the manuscript plate-only study); positive = too stiff in |u_y|... sign: negative = smaller deflection than Q9 (stiffer)"))
uy2 = last_uy(f"{RAW}/plate_only/S2_plate_hex8_fsiM2pressure/A0.txt", 2)
pl.append(dict(case="S2_plate_hex8_fsiM2pressure (hex8 80x8 under the M2 interface pressure, face-mean of nodal p)", element="hex8", uy_A_m=uy2, reference_Q9_same_mesh=None, rel_diff=uy2 / S["M2"]["uy_A_m"] - 1, normal_termination=True,
               note="rel_diff against the coupled M2 value (-0.226233 at t=870 s). Tests the PRESSURE transfer only; it cannot detect an error in the pressure itself (see docs/FINDINGS_REPORT.md)"))
write_csv("plate_only_checks.csv", pl)
# ---------------------------------------------------------------- force reconstruction table (from the facet tables; package-only)
import subprocess, sys
subprocess.run([sys.executable, "-I", "scripts/reconstruct_forces_from_facets.py", f"{OUT}/force_reconstruction_F1_F4.csv"], check=True, stdout=subprocess.DEVNULL)
print("wrote force_reconstruction_F1_F4.csv")
