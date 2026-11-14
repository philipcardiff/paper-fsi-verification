#!/usr/bin/env python3
"""QoIs of a Kratos 3dTube run, with the definitions of the solids4foam driver, and a comparison
table against solids4foam and the published implicit-Euler solutions.

The QoI helpers (peak, first_peak, arrival_time, front_arrival, fit_speed) and the reference data
(stations, wave-front rules, published values and digitised histories) are imported from the
solids4foam driver tube_3d_verification.py and its reference/ directory (commit a421f063d).

Usage: python tube_qoi.py <driver.py> <run dir> [<run dir> ...] [--out summary.md]
"""
import csv
import importlib.util
import json
import math
import os
import sys

# solids4foam level-1 values (paper-data/3dtube_level3_findings.md, Sections 6 and 12.6-12.7).
# The Euler u_z,min and u_r,min were recorded before the two solver fixes of Section 12, which
# changed u_z on the platform used; only u_r,max was re-checked after the fixes (0.12241 mm).
S4F = {
    "Euler L1, dt 1e-4 (after the fixes; only u_r,max re-checked)": {"ur_max_A_m": 0.12241e-3},
    "Euler L1, dt 1e-4 (before the fixes, Section 6)": {
        "ur_max_A_m": 0.12233e-3, "t_ur_max_A_s": 7.005e-3, "uz_min_A_m": -0.07601e-3,
        "ur_min_late_A_m": -0.07756e-3},
    "Euler L2, dt 1e-4 (before the fixes)": {
        "ur_max_A_m": 0.12107e-3, "t_ur_max_A_s": 7.037e-3, "uz_min_A_m": -0.07515e-3,
        "ur_min_late_A_m": -0.07611e-3},
    "Backward L1, dt 2.5e-5 (after the fixes)": {
        "ur_max_A_m": 0.15978e-3, "uz_min_A_m": -0.08819e-3, "t_arrival_A_s": 5.873e-3,
        "wave_speed_pressure_m_s": 4.654},
    "Backward L3, dt 6.25e-6 (after the fixes)": {
        "ur_max_A_m": 0.15851e-3, "t_ur_max_A_s": 7.252e-3, "uz_min_A_m": -0.08571e-3,
        "t_arrival_A_s": 5.825e-3, "wave_speed_pressure_m_s": 4.594},
}
QOIS = [("ur_max_A_m", "u_r,max(A)", 1e3, "mm"), ("t_ur_max_A_s", "t(u_r,max)", 1e3, "ms"),
        ("uz_min_A_m", "u_z,min(A)", 1e3, "mm"), ("ur_min_late_A_m", "u_r,min(A), t > 14 ms", 1e3, "mm"),
        ("t_arrival_A_s", "t_arr(A)", 1e3, "ms"), ("wave_speed_wall_m_s", "c_wall", 1.0, "m/s"),
        ("wave_speed_pressure_m_s", "c_p", 1.0, "m/s")]


def load_driver(path):
    spec = importlib.util.spec_from_file_location("tubedriver", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def read_run(run_dir):
    with open(os.path.join(run_dir, "history.csv")) as f:
        rows = list(csv.DictReader(f))
    case = json.load(open(os.path.join(run_dir, "case.json")))
    info = {}
    if os.path.exists(os.path.join(run_dir, "run_info.json")):
        info = json.load(open(os.path.join(run_dir, "run_info.json")))
    return rows, case, info


def extract(drv, ref, rows, end_time):
    """Mirror of the driver's extract(), on the Kratos history.csv."""
    stations = [float(z) for z in ref["stations_z_m"]]
    point_a = float(ref["point_A_z_m"])
    t = [float(r["time"]) for r in rows]

    def col(name):
        return [float(r[name]) for r in rows]

    def tag(z):
        return f"z{round(z * 1e4):03d}"

    radial = {z: list(zip(t, col(f"ur_{tag(z)}"))) for z in stations}
    axial_a = list(zip(t, col(f"uz_{tag(point_a)}")))
    pressures = {z: list(zip(t, col(f"p_{tag(z)}"))) for z in stations}
    out = {}
    out["t_ur_max_A_s"], out["ur_max_A_m"] = drv.peak(radial[point_a])
    incident_end = float(ref["incident_axial_trough_before_s"])
    out["t_uz_min_A_s"], out["uz_min_A_m"] = drv.peak([(a, b) for a, b in axial_a if a <= incident_end], -1.0)
    late = [(a, b) for a, b in radial[point_a] if a >= float(ref["late_minimum_after_s"])]
    if len(late) > 2 and end_time >= float(ref["end_time_s"]) - 1e-12:
        out["t_ur_min_late_A_s"], out["ur_min_late_A_m"] = drv.peak(late, -1.0)
    else:
        out["ur_min_late_A_m"] = out["t_ur_min_late_A_s"] = math.nan
    fit = ref["wave_front"]["fit_z_range_m"]
    in_fit = [fit[0] - 1e-9 <= z <= fit[1] + 1e-9 for z in stations]
    wall_arr, p_arr = [], []
    for z in stations:
        arrival, _, _ = drv.front_arrival(radial[z], ref)
        wall_arr.append(arrival)
        p_arr.append(drv.arrival_time(pressures[z], float(ref["wave_front"]["arrival_fraction"])
                                      * float(ref["inlet_pressure_Pa"])))
        if math.isclose(z, point_a):
            out["t_arrival_A_s"] = arrival
    fz = [z for z, k in zip(stations, in_fit) if k]
    out["wave_speed_wall_m_s"] = drv.fit_speed(fz, [a for a, k in zip(wall_arr, in_fit) if k])
    out["wave_speed_pressure_m_s"] = drv.fit_speed(fz, [a for a, k in zip(p_arr, in_fit) if k])
    out["_radial_A"], out["_axial_A"] = radial[point_a], axial_a
    return out


def fmt(v, scale):
    return "–" if v is None or (isinstance(v, float) and math.isnan(v)) else f"{v * scale:.5g}"


def main():
    argv = sys.argv[1:]
    out_path = None
    if "--out" in argv:
        i = argv.index("--out")
        out_path = argv[i + 1]
        del argv[i:i + 2]
    drv = load_driver(argv[0])
    here = os.path.dirname(os.path.abspath(argv[0]))
    # reference/ next to the driver (the copy used on smartsimet) or ../reference/ (solids4foam layout)
    ref_dir = next(d for d in (os.path.join(here, "reference"), os.path.join(here, "..", "reference"))
                   if os.path.exists(os.path.join(d, "3dTube_verification_references.json")))
    ref = json.load(open(os.path.join(ref_dir, "3dTube_verification_references.json")))
    lit = ref["published_checks"]["literature"]
    published = {"Lozovskiy et al. (2019), Euler 1e-4": {k: v["value"] for k, v in lit["lozovskiy2019"].items()},
                 "Eken (2016), Euler 1e-4": {k: v["value"] for k, v in lit["eken2016"].items()},
                 "Tuković et al. (2018), backward 2.5e-5": {
                     k: v["value"] for k, v in ref["published_checks"]["mesh"]["tukovic2018"].items()}}
    runs = []
    for run_dir in argv[1:]:
        rows, case, info = read_run(run_dir)
        q = extract(drv, ref, rows, float(rows[-1]["time"]))
        solves = [int(r["coupling_solves"]) for r in rows]
        q["_meta"] = {
            "name": os.path.basename(run_dir.rstrip("/")), "factor": case["factor"], "dt": case["dt"],
            "scheme": case["scheme"], "end": float(rows[-1]["time"]), "steps": len(rows),
            "solves_mean": sum(solves) / len(solves), "solves_max": max(solves),
            "not_converged": sum(1 for r in rows if r["converged"] != "1"),
            "max_residual": max(float(r["rms_residual_at_acceptance"]) for r in rows),
            "wall_s": info.get("wall_total_s"), "threads": info.get("kratos_threads"),
            "fluid_tets": case["fluid_tets"], "fluid_nodes": case["fluid_nodes"],
            "hex_equiv": case["fluid_hex_equivalent"], "solid_hexes": case["solid_hexes"]}
        runs.append(q)

    L = ["# Kratos 3dTube QoIs (definitions of the solids4foam driver)\n"]
    L.append("| run | level factor | scheme | dt (s) | to t (ms) | fluid tets (hex equiv.) | wall hexes | "
             "coupling solves mean/max | steps over tol | max residual (m) | wall (min) | threads |")
    L.append("|---|---|---|---|---|---|---|---|---|---|---|---|")
    for q in runs:
        m = q["_meta"]
        wall = f"{m['wall_s'] / 60:.1f}" if m["wall_s"] else "running"
        L.append(f"| {m['name']} | {m['factor']} | {m['scheme']} | {m['dt']:g} | {m['end'] * 1e3:.3g} | "
                 f"{m['fluid_tets']} ({m['hex_equiv']}) | {m['solid_hexes']} | {m['solves_mean']:.2f}/{m['solves_max']} | "
                 f"{m['not_converged']} | {m['max_residual']:.1e} | {wall} | {m['threads']} |")
    L.append("")
    header = "| quantity | " + " | ".join(q["_meta"]["name"] for q in runs) + " | " + \
             " | ".join(published) + " | " + " | ".join(f"solids4foam {k}" for k in S4F) + " |"
    L.append(header)
    L.append("|" + "---|" * (1 + len(runs) + len(published) + len(S4F)))
    for key, label, scale, unit in QOIS:
        cells = [fmt(q.get(key), scale) for q in runs]
        cells += [fmt(p.get(key), scale) for p in published.values()]
        cells += [fmt(s.get(key), scale) for s in S4F.values()]
        L.append(f"| {label} ({unit}) | " + " | ".join(cells) + " |")
    if runs:
        q = runs[0]
        L.append("\nDifferences of the first run from the references:\n")
        L.append("| quantity | " + " | ".join(list(published) + ["solids4foam Euler L1"]) + " |")
        L.append("|" + "---|" * (1 + len(published) + 1))
        s4f_euler = dict(S4F["Euler L1, dt 1e-4 (before the fixes, Section 6)"],
                         **S4F["Euler L1, dt 1e-4 (after the fixes; only u_r,max re-checked)"])
        for key, label, scale, unit in QOIS:
            v = q.get(key)
            cells = []
            for refv in list(published.values()) + [s4f_euler]:
                r = refv.get(key)
                cells.append("–" if r is None or v is None or math.isnan(v) else f"{100 * (v - r) / abs(r):+.1f}%")
            L.append(f"| {label} | " + " | ".join(cells) + " |")
    text = "\n".join(L) + "\n"
    print(text)
    if out_path:
        with open(out_path, "w") as f:
            f.write(text)
        with open(os.path.splitext(out_path)[0] + "_histories.csv", "w") as f:
            f.write("run,time_s,radial_A_m,axial_A_m\n")
            for q in runs:
                for (t, r), (_, a) in zip(q["_radial_A"], q["_axial_A"]):
                    f.write(f"{q['_meta']['name']},{t:.8g},{r:.8e},{a:.8e}\n")


if __name__ == "__main__":
    main()
