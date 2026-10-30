#!/usr/bin/env python3
"""Recompute the interface force definitions F1..F4 and the row-parity decomposition from the per-facet tables stored in results/raw/facets/
(extracted once from the FEBio .xplt files by scripts/hpc/b3_force_reconstruction.py).  Package-only: needs numpy, no FEBio output files.
This is a RE-EVALUATION OF STORED SOLVER FIELDS with the quadrature of the FEBio source routines.  It is NOT an independent solver verification.
  F1 = pressure part of FEFluidFSITraction::LoadVector (nodal wall-row pressure, 2x2 Gauss, deformed facet)   = load that moves the plate
  F2 = element-mean pressure of the adjacent fluid element x area vector (FEBio's built-in 'fluid pressure' / surface-force definition)
  F3 = element pressures of element rows 1,2 linearly extrapolated to the wall (parity-immune by construction)
  F4 = nodal pair-row mean (row0+row1)/2, bilinear Gauss (F2 = F4 for planar facets, so F4 is NOT an independent estimator)
Viscous contribution: NOT reconstructed here (needs Gauss-point velocity gradients that are not stored).  Usage (from the package root):
  python3 -I scripts/reconstruct_forces_from_facets.py [OUT.csv]"""
import sys, os, glob, json, numpy as np
G = 1 / np.sqrt(3); GP = [(-G, -G), (G, -G), (G, G), (-G, G)]
def shp(r, s): return np.array([(1 - r) * (1 - s), (1 + r) * (1 - s), (1 + r) * (1 + s), (1 - r) * (1 + s)]) / 4
def dshp(r, s): return (np.array([-(1 - s), (1 - s), (1 + s), -(1 + s)]) / 4, np.array([-(1 - r), -(1 + r), (1 + r), (1 - r)]) / 4)
def gt(x4, r, s):
    dr, ds = dshp(r, s); return np.cross(dr @ x4, ds @ x4)
def recompute(d):
    TZ = float(d["TZ"]); X4 = d["X4"]; p4 = d["p4"]; p14 = d["p14"]; Pe1 = d["Pe1"]; Pe2 = d["Pe2"]; C1 = d["C1"]; C2 = d["C2"]
    F1 = np.zeros(3); F2 = np.zeros(3); F3 = np.zeros(3); F4 = np.zeros(3); A_tot = 0.0
    for i in range(len(X4)):
        x4 = X4[i]; pm4 = 0.5 * (p4[i] + p14[i]); area = np.zeros(3)
        for (rr, ss) in GP:
            N = shp(rr, ss); g = gt(x4, rr, ss); area += g; F1 += -(N @ p4[i]) * g; F4 += -(N @ pm4) * g
        n_c = gt(x4, 0, 0); n_c = n_c / np.linalg.norm(n_c); A = np.linalg.norm(area); A_tot += A
        F2 += -Pe1[i] * n_c * A
        mid = x4.mean(0); s1 = (C1[i] - mid) @ n_c
        if np.isfinite(Pe2[i]): s2 = (C2[i] - mid) @ n_c; pw = Pe1[i] + (Pe1[i] - Pe2[i]) * s1 / (s2 - s1)
        else: pw = Pe1[i]
        F3 += -pw * n_c * A
    P0, P1, P2 = d["P0"], d["P1"], d["P2"]; ok = np.isfinite(P2); ok[:4] = False; ok[-4:] = False
    S_ = (P0 + 2 * P1 + P2) / 4
    pn = np.where(np.isfinite(S_), S_, P0); Fs = np.zeros(3)
    for i in range(len(X4)):
        pp = pn[d["q_cols"][i]]
        for (rr, ss) in GP: Fs += -(shp(rr, ss) @ pp) * gt(X4[i], rr, ss)
    f = lambda v: v / TZ
    return dict(F1_Fy=f(F1[1]), F2_Fy=f(F2[1]), F3_Fy=f(F3[1]), F4_Fy=f(F4[1]), F1_Fx=f(F1[0]), F2_Fx=f(F2[0]),
                F1_minus_F2=f(F1[1] - F2[1]), parity_part_of_F1=f(F1[1] - Fs[1]), D_row1_minus_row0_Pa=float(np.mean((P1 - P0)[ok])), wetted_length_m=A_tot / TZ)
if __name__ == "__main__":
    out = sys.argv[1] if len(sys.argv) > 1 else "results/force_reconstruction_F1_F4.csv"
    keys = ["M1", "M2", "M3", "M3_t717p5_failure_exit_state"]; rows = []
    for k in keys:
        d = np.load(f"results/raw/facets/facets_{k}.npz"); r = recompute(d)
        rows.append(dict(run_state=k, time_s=float(d["time"]), solver_flag=int(d["flag"]), **r))
    cols = ["run_state", "time_s", "solver_flag", "F1_Fy", "F2_Fy", "F3_Fy", "F4_Fy", "F1_minus_F2", "parity_part_of_F1", "D_row1_minus_row0_Pa", "F1_Fx", "F2_Fx", "wetted_length_m"]
    with open(out, "w") as fh:
        fh.write(",".join(cols) + ",units\n")
        for r in rows: fh.write(",".join(f"{r[c]:.6f}" if isinstance(r[c], float) else str(r[c]) for c in cols) + ",N/m (force of the fluid on the plate per m depth; pressure contribution only)\n")
    print(open(out).read())
