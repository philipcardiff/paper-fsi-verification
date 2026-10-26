"""cavityFlexibleBottom: Richardson extrapolation of the four references
(solids4foam, Tukovic et al. 2018, Ansys, COMSOL) at point A.

Writes reference_levels.csv (per-level values) and reference_richardson.csv
(observed order and extrapolated value from each triplet of consecutive
levels). Run: python reference_convergence.py

Sources (per unit depth, steady state):
- solids4foam: u_y from sectionBenchmarkSuite.tex; F_y from the
  quasi-monolithic rerun (solids4foam f356c3242,
  cavityFlexibleBottom/verification/results/four_level_partitioned_monolithic.csv).
- Ansys: PR #12 (E. Tandis), AnsysResults/cavityFlexibleBottom-Ansys-Dx*.csv
  at fe69944, last row (t = 1000 s), F_y x 100 for the 10 mm slab depth.
- COMSOL: PR #12, ComsolResults/cavityFlexibleBottom_comsol_steadyStateResults.csv
  at 551e91e (steady solver). u_y is given to 3 significant figures only.
- Tukovic et al. (2018), Fig. 12(a): read from the vector paths of the
  published PDF (page 16), not hand digitised. The page coordinates of the
  four markers and of the axis ticks are listed below. The published forces
  are per 0.05 m depth (confirmed by Z. Tukovic, Oct 2026), so F_y is
  divided by 0.05 here. The extraction reproduces the paper's stated orders
  (1.85 for u_y, 1.91 for F_y).
"""
import csv
import math
import os

DX = [0.1, 0.05, 0.025, 0.0125]

# Tukovic Fig. 12(a), page coordinates (pt, y downwards), levels Dx = 0.1 .. 0.0125
TUK_UY_Y = [475.16, 533.48, 551.43, 556.42]
TUK_FY_Y = [500.33, 553.11, 568.76, 572.91]
# left axis: y = 576.92 -> u_y = -0.255, y = 483.39 -> u_y = -0.235
# right axis: y = 576.92 -> F_y = -0.263 N, y = 469.36 -> F_y = -0.256 N
TUK_DEPTH = 0.05


def tuk_uy(y):
    return -0.255 + (576.92 - y) / (576.92 - 483.39) * 0.020


def tuk_fy(y):
    return (-0.263 + (576.92 - y) / (576.92 - 469.36) * 0.007) / TUK_DEPTH


SOURCES = {
    "solids4foam": ([-0.2036, -0.2314, -0.2385, -0.2401],
                    [-5.123, -5.228, -5.251, -5.256]),
    "Tukovic2018": ([round(tuk_uy(y), 5) for y in TUK_UY_Y],
                    [round(tuk_fy(y), 4) for y in TUK_FY_Y]),
    "Ansys": ([-0.24727, -0.25054, -0.25197, -0.25274],
              [-5.144, -5.218, -5.250, -5.263]),
    "COMSOL": ([-0.236, -0.241, -0.243, -0.243],
               [-5.1088, -5.2276, -5.2956, -5.2952]),
}


def richardson(f3, f2, f1, r=2.0):
    """Coarse f3, medium f2, fine f1. Returns (p, extrapolated) or None
    when the triplet is not monotone (Richardson not applicable)."""
    d32, d21 = f3 - f2, f2 - f1
    if d21 == 0 or d32 / d21 <= 0:
        return None
    p = math.log(d32 / d21) / math.log(r)
    return p, f1 + (f1 - f2) / (r**p - 1)


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, "reference_levels.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["source", "dx_m", "uy_A_m", "Fy_N_per_m"])
        for src, (uy, fy) in SOURCES.items():
            for dx, u, f in zip(DX, uy, fy):
                w.writerow([src, dx, u, f])
    with open(os.path.join(here, "reference_richardson.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["source", "quantity", "levels_dx", "observed_order",
                    "extrapolated", "note"])
        for src, (uy, fy) in SOURCES.items():
            for q, v in (("uy_A_m", uy), ("Fy_N_per_m", fy)):
                for i in range(len(v) - 2):
                    lev = f"{DX[i]}/{DX[i+1]}/{DX[i+2]}"
                    res = richardson(v[i], v[i + 1], v[i + 2])
                    if res is None:
                        w.writerow([src, q, lev, "", "",
                                    "not monotone: Richardson not applicable"])
                        print(f"{src:12s} {q:11s} {lev:17s} not monotone")
                    else:
                        p, ext = res
                        w.writerow([src, q, lev, f"{p:.2f}", f"{ext:.5f}", ""])
                        print(f"{src:12s} {q:11s} {lev:17s} p={p:5.2f} ext={ext:.5f}")


if __name__ == "__main__":
    main()
