#!/usr/bin/env python3
"""Stage B equilibrium evidence from STORED output only (no simulation).  FEBio's node_data 'Rx;Ry;Rz' is zero because the 'zero displacement' boundary condition makes the DOFs
DOF_FIXED (no equation number): FESolidSolver2::Residual sets node loads only for nodes with m_ID >= 0 or the prescribed code -ID-2 >= 0 (FEBioMech/FESolidSolver2.cpp l.1630-1642), so a fixed DOF
never receives a reaction.  The check done here instead: free-body equilibrium of the part of the flag beyond a vertical section at x = x_s from the stored element Cauchy stress (xplt 'stress',
one value per element): sum over the section of  0.5(sigma_xy[col k] + sigma_xy[col k+1]) dy TZ  must equal -(weight of the free part) = -rho |g| (L - s) h TZ (small-slope approximation:
the true section is rotated, so the ratio is expected to deviate by about 1-cos(theta) and the axial stress resultant by sin(theta) x shear), and the moment of sigma_xx about the section centre
against the beam value rho|g| h TZ (L-s)^2/2 (magnitude).  usage: equilibrium_from_xplt.py RUNS_DIR OUT.csv [case ...]   (cases default: all B_CSM*_h8/h20 runs with an .xplt)"""
import sys, os, json, glob, csv
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from xplt_reader import read_xplt
runs = os.path.abspath(sys.argv[1]); out = sys.argv[2]
cases = sys.argv[3:] or sorted(os.path.basename(p)[:-5] for p in glob.glob(runs + "/B_CSM*/*.xplt") if "root025" not in p)
rows = []
for c in cases:
    d = f"{runs}/{c}"; m = json.load(open(d + "/model_meta.json")); nt, nx, L = m["nt"], m["nx"], m["length"]; h, TZ, rho, g = 0.02, 0.015, 1000.0, abs(m["g"])
    if m["elem"] != "hex8" and m["elem"] != "hex20": continue
    S = read_xplt(f"{d}/{c}.xplt")["states"][-1]; st = S["vars"]["stress"][0]; sxx, sxy = st[:, 0], st[:, 3]
    dy, dx = h / nt, L / nx; yy = (np.arange(nt) + 0.5) * dy - h / 2
    for ex in (1, nx // 4, nx // 2, 3 * nx // 4, nx - 2):
        col = lambda k, v: v[k * nt:(k + 1) * nt]
        V = float(np.sum(0.5 * (col(ex, sxy) + col(ex + 1, sxy))) * dy * TZ); N = float(np.sum(0.5 * (col(ex, sxx) + col(ex + 1, sxx))) * dy * TZ)
        Mm = float(np.sum(0.5 * (col(ex, sxx) + col(ex + 1, sxx)) * yy) * dy * TZ)
        s = (ex + 1) * dx; W = rho * g * (L - s) * h * TZ; Mb = rho * g * h * TZ * (L - s) ** 2 / 2
        rows.append(dict(case=c, E=m["E"], elem=m["elem"], nt=nt, x_section_from_root=round(s, 5), shear_resultant=V, minus_weight_free_part=-W, shear_ratio=V / (-W),
                         axial_resultant=N, axial_over_weight=N / W, moment_sxx=Mm, moment_beam_magnitude=Mb, moment_ratio=abs(Mm) / Mb))
with open(out, "w", newline="") as fh:
    w = csv.DictWriter(fh, list(rows[0])); w.writeheader(); w.writerows(rows)
for r in rows:
    if r["nt"] == max(x["nt"] for x in rows if x["case"].split("_")[1] == r["case"].split("_")[1] and x["elem"] == "hex8") and r["elem"] == "hex8":
        print(r["case"], r["x_section_from_root"], "shear ratio %.5f moment ratio %.4f" % (r["shear_ratio"], r["moment_ratio"]))
