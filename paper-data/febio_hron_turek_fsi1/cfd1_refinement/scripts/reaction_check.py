#!/usr/bin/env python3
"""Post-processing check of the reaction-transfer idea on STOCK output (no FEBio run): compute the fluid momentum residual (strong-pressure-gradient form, FEBio's
own w-row) at the interface nodes from logged nodal fields, R_a = int N_a (rho a + grad p) dV (viscous term negligible at nu = 1e-7, reported separately by the
traction audit), over the layer-0 fluid elements, and compare the work it would do on the ring, sum R_a . du, with the ring's energy change (which the stock traction
work reproduces to 0.02-0.09 %). usage: reaction_check.py <run_dir> <energy_npz> <out_json>"""
import sys, json, math
import numpy as np, xml.etree.ElementTree as ET
sys.path.insert(0, "<project-root>/tools"); import analyze_ring as AR
run = sys.argv[1].rstrip("/"); en = np.load(sys.argv[2]); out = sys.argv[3]
meta = json.load(open(run + "/model_meta.json")); tag = meta["tag"]; nrs, nt, nrf = meta["dims_nrs_nt_nrf"]; K = meta["k"]; RHO = meta["rho_f"]
root = ET.parse(f"{run}/{tag}.feb").getroot()
X = {int(n.get("id")): np.array([float(v) for v in n.text.split(",")]) for n in root.find("Mesh/Nodes")}
fl = []; first = None
for el in root.findall("Mesh/Elements"):
    if el.get("name") == "Fluid":
        for e in el: fl.append((int(e.get("id")), [int(v) for v in e.text.split(",")]))
first = fl[0][0]; layer0 = [fl[j * nrf] for j in range(nt)]
nodes = sorted({n for _, c in layer0 for n in c}); idx = {n: i for i, n in enumerate(nodes)}
fb = AR.read_blocks(run + "/fluid_nodes.txt"); eb = AR.read_blocks(run + "/fluid_ef.txt"); rb = AR.read_blocks(run + "/ring_nodes.txt")
def tab(blocks, ids, c):
    o = np.zeros((len(blocks), len(ids), c)); pos = None
    for k, b in enumerate(blocks):
        r = np.array(b["rows"]); m = {int(x): i for i, x in enumerate(r[:, 0])}
        if pos is None: pos = [m[n] for n in ids]
        o[k] = r[pos, 1:1 + c]
    return o
t = np.array([b["t"] for b in fb]); nT = len(t)
F = tab(fb, nodes, 6); EF = tab(eb, nodes, 1)[:, :, 0]                      # ux uy vx vy wx wy ; ef
vabs = F[:, :, 2:4] + F[:, :, 4:6]; U = F[:, :, 0:2]
gp = 1 / math.sqrt(3); SG = np.array([[-1,-1,-1],[1,-1,-1],[1,1,-1],[-1,1,-1],[-1,-1,1],[1,-1,1],[1,1,1],[-1,1,1]], float)
def shp(r, s, tt): return np.array([(1+SG[a,0]*r)*(1+SG[a,1]*s)*(1+SG[a,2]*tt)/8 for a in range(8)])
def dshp(r, s, tt): return np.array([[SG[a,0]*(1+SG[a,1]*s)*(1+SG[a,2]*tt)/8, SG[a,1]*(1+SG[a,0]*r)*(1+SG[a,2]*tt)/8, SG[a,2]*(1+SG[a,0]*r)*(1+SG[a,1]*s)/8] for a in range(8)])
GPS = [(a, b, c) for a in (-gp, gp) for b in (-gp, gp) for c in (-gp, gp)]
elem = []
for _, conn in layer0:
    P = np.array([X[n] for n in conn]); li = [idx[n] for n in conn]
    q = []
    for (a, b, c) in GPS:
        N = shp(a, b, c); dN = dshp(a, b, c); J = P.T @ dN; dj = np.linalg.det(J); q.append((N, dN @ np.linalg.inv(J), dj))
    elem.append((li, q))
iface = [n for n in nodes if abs(math.hypot(X[n][0], X[n][1]) - 1.05) < 1e-9]
isf = [idx[n] for n in iface]
R_in = np.zeros((nT - 1, len(nodes), 2)); R_pg = np.zeros_like(R_in)
for k in range(nT - 1):
    dt = t[k + 1] - t[k]; acc = (vabs[k + 1] - vabs[k]) / dt; p = -K * 0.5 * (EF[k] + EF[k + 1])
    for li, q in elem:
        for (N, B, dj) in q:
            a_gp = N @ acc[li]; gradp = B.T @ p[li]; w = dj
            for a, ia in enumerate(li):
                R_in[k, ia] += N[a] * RHO * a_gp * w; R_pg[k, ia] += N[a] * gradp[:2] * w
du = U[1:] - U[:-1]
W_in = np.array([(R_in[k, isf] * du[k, isf]).sum() for k in range(nT - 1)]); W_pg = np.array([(R_pg[k, isf] * du[k, isf]).sum() for k in range(nT - 1)])
Wc = np.concatenate([[0], np.cumsum(W_in + W_pg)]); Wi = np.concatenate([[0], np.cumsum(W_in)]); Wp = np.concatenate([[0], np.cumsum(W_pg)])
Er = en["KEr"] + en["PEr"]; dE = Er - Er[0]
# the reaction force acts on the ring with sign s; report both signs and the best-fit scale
for name, W in (("reaction_total", Wc), ("reaction_inertia", Wi), ("reaction_pressure_gradient", Wp)):
    pass
A = np.column_stack([Wc]); scale = float(np.linalg.lstsq(A, dE, rcond=None)[0][0])
res = dict(run=run, dE_ring_range=float(dE.max() - dE.min()), W_react_range=float(Wc.max() - Wc.min()), W_inertia_range=float(Wi.max() - Wi.min()), W_pgrad_range=float(Wp.max() - Wp.min()),
           best_scale_of_reaction_to_match_ring_energy=scale, max_abs_mismatch_plus=float(np.abs(dE - Wc).max()), max_abs_mismatch_minus=float(np.abs(dE + Wc).max()),
           corr_plus=float(np.corrcoef(dE, Wc)[0, 1]))
json.dump(res, open(out, "w"), indent=1); print(json.dumps(res, indent=1))
np.savez_compressed(out.replace(".json", ".npz"), t=t, dE=dE, Wc=Wc, Wi=Wi, Wp=Wp)
