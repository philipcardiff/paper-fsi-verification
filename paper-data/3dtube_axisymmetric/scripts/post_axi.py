#!/usr/bin/env python3
"""Point-A QoIs of the axisymmetric 3dTube runs, with the 3-D driver's logic.

The peak, first-peak, arrival-time and wave-speed functions are imported from
the 3-D verification driver (solids4foam 206437983 plus the monitor-time
snapping patch), so the definitions are those of Table tab:3dtube. The only
difference is the orientation: the axis is x, and the wedge monitors sit on
a wedge face at half the wedge angle to the z = 0 plane, so u_r = sign(Dy) *
hypot(Dy, Dz) and u_z = Dx, instead of the 3-D driver's Dy and Dz.

Usage:
  post_axi.py <case> [<case> ...] [--csv out.csv] [--histories dir]
              [--end-time 0.02] [--allow-partial]
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location(
    "tube_3d_verification", HERE / "driver" / "tube_3d_verification.py")
drv = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(drv)
REFERENCE = json.loads(
    (HERE / "driver" / "3dTube_verification_references.json").read_text())


TRUNCATE = None  # set by --truncate: drop rows with t > TRUNCATE


def _cut(rows):
    if TRUNCATE is None:
        return rows
    return [r for r in rows if r[0] <= TRUNCATE + 1e-12]


def wall_history(case: Path, name: str, end_time: float, complete: bool):
    files = sorted(case.glob(
        f"postProcessing/**/solidPointDisplacement_{name}.dat"))
    if not files:
        drv.fail(f"monitor {name} not found in {case}")
    rows = drv.last_per_time(drv.numeric_rows(files[0]))
    if complete:
        drv.require_complete(files[0], rows, end_time, 5)
    rows = _cut(rows)
    # axis x, radial y, wedge across z
    return [(r[0], math.copysign(math.hypot(r[2], r[3]), r[2]), r[1])
            for r in rows]


def pressure_hist(case: Path, n: int, end_time: float, rho: float,
                  complete: bool):
    files = sorted(case.glob("postProcessing/**/verificationAxisPressure/**/p"))
    if not files:
        drv.fail(f"axis probes not found in {case}")
    rows = drv.last_per_time(drv.numeric_rows(files[0]))
    if complete:
        drv.require_complete(files[0], rows, end_time, n + 1)
    rows = _cut(rows)
    return [[(r[0], rho * r[i + 1]) for r in rows] for i in range(n)]


def extract_axi(case: Path, end_time: float, complete: bool = True) -> dict:
    """drv.extract() with the wedge radial displacement."""
    ref = REFERENCE
    stations = [float(z) for z in ref["stations_z_m"]]
    point_a = float(ref["point_A_z_m"])
    fr = ref["wave_front"]["fit_z_range_m"]
    in_fit = [fr[0] - 1e-9 <= z <= fr[1] + 1e-9 for z in stations]
    wall = {z: wall_history(case, drv.station_name(z), end_time, complete)
            for z in stations}
    radial = {z: [(t, r) for t, r, _ in wall[z]] for z in stations}
    axial_a = [(t, a) for t, _, a in wall[point_a]]
    pressures = pressure_hist(case, len(stations), end_time,
                              float(ref["fluid"]["rho_kg_m3"]), complete)
    row: dict = {}
    row["ur_max_A_m"], row["t_ur_max_A_s"] = reversed(drv.peak(radial[point_a]))
    inc = float(ref["incident_axial_trough_before_s"])
    row["uz_min_A_m"], row["t_uz_min_A_s"] = reversed(drv.peak(
        [(t, a) for t, a in axial_a if t <= inc], -1.0))
    late0 = float(ref["late_minimum_after_s"])
    late = [(t, r) for t, r in radial[point_a] if t >= late0]
    if len(late) > 2 and late[-1][0] >= float(ref["end_time_s"]) - 1e-9:
        row["ur_min_late_A_m"], row["t_ur_min_late_A_s"] = reversed(
            drv.peak(late, -1.0))
    else:
        row["ur_min_late_A_m"] = row["t_ur_min_late_A_s"] = math.nan
    wa, pa = [], []
    for i, z in enumerate(stations):
        arr, pt, pv = drv.front_arrival(radial[z], ref)
        wa.append(arr)
        _, _, ppk = drv.front_arrival(pressures[i], ref)
        parr = drv.arrival_time(
            pressures[i], float(ref["wave_front"]["arrival_fraction"])
            * float(ref["inlet_pressure_Pa"]))
        pa.append(parr)
        tag = f"z{round(z * 1e3):02d}mm"
        row[f"wall_arrival_{tag}_s"] = arr
        row[f"wall_first_peak_{tag}_m"] = pv
        row[f"pressure_arrival_{tag}_s"] = parr
        if math.isclose(z, point_a):
            row["t_arrival_A_s"] = arr
    fz = [z for z, k in zip(stations, in_fit) if k]
    row["wave_speed_wall_m_s"] = drv.fit_speed(
        fz, [t for t, k in zip(wa, in_fit) if k])
    row["wave_speed_pressure_m_s"] = drv.fit_speed(
        fz, [t for t, k in zip(pa, in_fit) if k])
    row["_histories"] = {"radial_A": radial[point_a], "axial_A": axial_a,
                         "pressure_A": pressures[stations.index(point_a)]}
    return row


def coupling_stats(case: Path) -> dict:
    out = {}
    try:
        files = sorted(case.glob("postProcessing/**/fsiConvergenceData.dat"))
        rows = drv.last_per_time(drv.numeric_rows(files[0]))
        its = [r[1] for r in rows]
        out.update(time_steps=len(its), mean_fsi_iterations=sum(its) / len(its),
                   maximum_fsi_iterations=max(its))
        path = case / "postProcessing/fsiResiduals.dat"
        header = path.read_text(errors="replace").splitlines()[0].split()
        res = drv.last_per_time(drv.numeric_rows(path))
        st = header.index("robinConvergenceState")
        out["unconverged_steps"] = sum(r[st] <= 0 for r in res)
        out["stalled_steps"] = sum(r[st] == 2 for r in res)
    except Exception as e:  # noqa: BLE001
        out["coupling_error"] = str(e)
    try:
        out["clock_time_s"] = drv.clock_time(case)
    except Exception:  # noqa: BLE001
        out["clock_time_s"] = math.nan
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("cases", nargs="+", type=Path)
    ap.add_argument("--csv", type=Path)
    ap.add_argument("--histories", type=Path)
    ap.add_argument("--end-time", type=float, default=0.02)
    ap.add_argument("--allow-partial", action="store_true")
    ap.add_argument("--truncate", type=float,
                    help="evaluate the histories only up to this time (s); "
                         "the point-A QoIs and c_p are fixed by 10.3 ms")
    args = ap.parse_args()
    global TRUNCATE
    TRUNCATE = args.truncate
    rows = []
    for case in args.cases:
        meta = json.loads((case / "axi_settings.json").read_text())
        row = extract_axi(case, args.end_time, not args.allow_partial)
        hist = row.pop("_histories")
        row = {"case": case.name, **{k: meta[k] for k in (
            "level", "law", "fluid_cells", "solid_cells", "delta_t_s",
            "wedge_deg", "cores")}, **row, **coupling_stats(case)}
        rows.append(row)
        if args.histories:
            args.histories.mkdir(parents=True, exist_ok=True)
            ax, pr = dict(hist["axial_A"]), dict(hist["pressure_A"])
            with (args.histories / f"{case.name}.csv").open("w") as f:
                f.write("time_s,radial_displacement_m,axial_displacement_m,"
                        "axis_pressure_Pa\n")
                for t, r in hist["radial_A"]:
                    f.write(f"{t:.10g},{r:.10g},{ax.get(t, math.nan):.10g},"
                            f"{pr.get(t, math.nan):.10g}\n")
        print(f"{case.name}: L{meta['level']} {meta['law']} "
              f"ur_max={row['ur_max_A_m']*1e3:.6f} mm "
              f"t(ur_max)={row['t_ur_max_A_s']*1e3:.4f} ms "
              f"uz_min={row['uz_min_A_m']*1e3:.6f} mm "
              f"t_arr={row['t_arrival_A_s']*1e3:.4f} ms "
              f"c_p={row['wave_speed_pressure_m_s']:.4f} "
              f"c_wall={row['wave_speed_wall_m_s']:.4f} m/s "
              f"its={row.get('mean_fsi_iterations', math.nan):.2f} "
              f"unconv={row.get('unconverged_steps')}")
    if args.csv:
        cols = []
        for r in rows:
            cols += [k for k in r if k not in cols]
        with args.csv.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=cols)
            w.writeheader()
            w.writerows(rows)
    return 0


if __name__ == "__main__":
    sys.exit(main())
