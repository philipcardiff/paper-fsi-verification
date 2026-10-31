#!/usr/bin/env python3
"""Phase-2 reconstruction of the FEBio v4.13 fluid-FSI interface traction (FEFluidFSITraction::LoadVector) applied to the stored RIGID-flag CFD1 fields (NO coupled run).
Source facts (upstream/febio-v4.13/FEBioFluid/FEFluidFSITraction.cpp, l.~100-230): f = gt*p - sv*gt at each facet Gauss point, gt = gr x gs (covariant surface basis on the CURRENT facet
geometry, oriented out of the adjacent fluid element via FacePointing), p = Pressure(ef) with ef interpolated from the facet's NODAL dilatation (H shape functions), sv = AVERAGE over all
Gauss points of the adjacent fluid-FSI element of the Newtonian viscous stress  s = (kappa - 2mu/3) tr(D) I + 2 mu D  (FENewtonianFluid::Stress l.68-80), nodal force a: fa = H_a f * w_gp.
Reconstruction here for the z-extruded slab (2-D): pressure part = integral of the linear nodal interpolant (exact for the straight facet) ; viscous part = element-mean over the 2x2 in-plane
Gauss points of sv (z-uniform).  Variants: with and without the -2mu/3 tr(D) term.  usage: hronturek_fsi_traction_recon.py RUN_DIR [OUT.json]"""
import sys, os, re, json, math, glob
import numpy as np
RUN = os.path.realpath(sys.argv[1]); OUT = sys.argv[2] if len(sys.argv) > 2 else RUN + "/fsi_traction_reconstruction.json"
feb = glob.glob(RUN + "/*.feb")[0]; txt = open(feb, encoding="latin-1").read()
K = float(re.search(r"<k>([^<]+)</k>", txt).group(1)); MU = float(re.search(r"<mu>([^<]+)</mu>", txt).group(1)); KAPPA = float(re.search(r"<kappa>([^<]+)</kappa>", txt).group(1))
nodes = np.array([[float(v) for v in m.group(2).split(",")] for m in re.finditer(r'<node id="(\d+)">([^<]+)</node>', txt)])
els = np.array([[int(v) - 1 for v in m.group(1).split(",")] for m in re.finditer(r'<elem id="\d+">([^<]+)</elem>', txt)])
z0 = np.where(np.abs(nodes[:, 2]) < 1e-12)[0]; remap = -np.ones(len(nodes), int); remap[z0] = np.arange(len(z0)); XY = nodes[z0, :2]; Q = remap[els[:, :4]]; n2 = len(z0)
rows = None
for ln in open(RUN + "/Mon_all.txt"):
    if ln.startswith("*Time"): T = float(ln.split("=")[1]); rows = []
    elif not ln.startswith("*") and ln.strip(): rows.append(ln.split())
A = np.array([[float(x) for x in r[1:4]] for r in rows]); W = A[:, :2]; P = -K * A[:, 2]
g = 1 / math.sqrt(3); GP = [(-g, -g), (g, -g), (g, g), (-g, g)]
def dN(xi, eta): return 0.25 * np.array([[-(1 - eta), -(1 - xi)], [(1 - eta), -(1 + xi)], [(1 + eta), (1 + xi)], [-(1 + eta), (1 - xi)]])
def Lgrad(qi, xi, eta):
    d = dN(xi, eta); J = XY[Q[qi]].T @ d; return W[Q[qi]].T @ (d @ np.linalg.inv(J))
def sv_mean(qi, trace):
    s = np.zeros((2, 2))
    for xi, eta in GP:
        L = Lgrad(qi, xi, eta); D = 0.5 * (L + L.T); s += 2 * MU * D + (np.eye(2) * np.trace(D) * (KAPPA - 2 * MU / 3) if trace else 0)
    return s / 4
cnt = {}
for qi, q in enumerate(Q):
    for k in range(4): a, b = int(q[k]), int(q[(k + 1) % 4]); cnt.setdefault((min(a, b), max(a, b)), []).append((a, b, qi))
bed = [v[0] for v in cnt.values() if len(v) == 1]
def obstacle(a, b):
    c = lambda i: (abs(XY[i, 0]) < 1e-9, abs(XY[i, 0] - 2.5) < 1e-9, abs(XY[i, 1]) < 1e-9, abs(XY[i, 1] - 0.41) < 1e-9)
    return not any(x and y for x, y in zip(c(a), c(b)))
Fp = np.zeros(2); Fv = {True: np.zeros(2), False: np.zeros(2)}; Fvc = np.zeros(2)
for a, b, qi in bed:
    if not obstacle(a, b): continue
    d = XY[b] - XY[a]; ln = math.hypot(*d); n = np.array([d[1], -d[0]]) / ln
    Fp += 0.5 * (P[a] + P[b]) * ln * n
    for tr in (True, False): Fv[tr] += -(sv_mean(qi, tr) @ n) * ln
    Lc = Lgrad(qi, 0.0, 0.0); Dc = 0.5 * (Lc + Lc.T); Fvc += -(2 * MU * Dc @ n) * ln
res = dict(run=os.path.basename(RUN), t=T, pressure_nodal_facet=Fp.tolist(), viscous_GPmean_FEBio_form=Fv[True].tolist(), viscous_GPmean_2muD_only=Fv[False].tolist(), viscous_centre_2muD=Fvc.tolist(),
           F_actual_recon=(Fp + Fv[True]).tolist(), F_actual_recon_no_trace_term=(Fp + Fv[False]).tolist())
# z-uniformity / data-path check: element-average velocity of the xplt vs mean of the logged z=0 nodal velocities
try:
    import importlib.util
    spec = importlib.util.spec_from_file_location("xa", os.path.join(os.path.dirname(os.path.abspath(__file__)), "xplt_reader.py")); xa = importlib.util.module_from_spec(spec); sys.argv = ["x"]
    try: spec.loader.exec_module(xa)
    except SystemExit: pass
    r = xa.read_xplt(glob.glob(RUN + "/*.xplt")[0]); S = r["states"][-1]; ve = S["vars"]["fluid velocity"][0][:, :2]
    res["xplt_element_velocity_vs_nodal_mean_max_abs_diff"] = float(np.abs(ve - W[Q].mean(axis=1)).max()); res["velocity_scale"] = float(np.abs(W).max())
    res["xplt_has_velocity_gradients"] = False
except Exception as exc: res["xplt_check"] = str(exc)
# magnitude of tr(D): max over element centres
trD = np.array([np.trace(0.5 * (Lgrad(i, 0, 0) + Lgrad(i, 0, 0).T)) for i in range(len(Q))]); res["trD_max_abs_element_centre"] = float(np.abs(trD).max()); res["trD_rms"] = float(np.sqrt(np.mean(trD**2)))
json.dump(res, open(OUT, "w"), indent=1); print(json.dumps(res, indent=1))
