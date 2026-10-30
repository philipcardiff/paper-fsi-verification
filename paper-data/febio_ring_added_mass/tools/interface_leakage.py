#!/usr/bin/env python3
"""Normal leakage through the ring-fluid interface from logged nodal data. usage: interface_leakage.py <run_dir>
Uses fluid_nodes.txt rows: id ux uy vx vy nfvx nfvy (v = mesh velocity, nfv = relative fluid velocity w).
Interface nodes are those at reference radius 1.05. The unit normal is the current radial direction.
w_n = w . n (relative flux through the moving ring surface); v_s,n = ring normal velocity at the node.
Reports: relative rms leakage, max |w_n|/V0, net (n=0) leaked flux and its ratio to the peak ring flux, n=2 content of w_n.
"""
import sys, math, json, numpy as np
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__))); import analyze_ring as AR
import xml.etree.ElementTree as ET
run = sys.argv[1].rstrip('/'); meta = json.load(open(run + '/model_meta.json')); tag = meta['tag']
nrs, nt, nrf = meta['dims_nrs_nt_nrf']; V0 = 1.5e-3; Ro = 1.05; h = 0.1
root = ET.parse(run + '/' + tag + '.feb').getroot()
N = {int(n.get('id')): np.array([float(v) for v in n.text.split(',')]) for n in root.find('Mesh/Nodes')}
fl = sorted({int(v) for e in root.find("Mesh/Elements[@name='Fluid']") for v in e.text.split(',')})
iface = [i for i in fl if abs(math.hypot(*N[i][:2]) - Ro) < 1e-9]
pos = {i: k for k, i in enumerate(fl)}
blocks = AR.read_blocks(run + '/fluid_nodes.txt'); t = np.array([b['t'] for b in blocks])
A = np.array([b['rows'] for b in blocks])                      # steps x nodes x cols
wn = np.zeros((len(t), len(iface))); vn = np.zeros_like(wn); th = np.zeros(len(iface)); z = np.zeros(len(iface))
for c, i in enumerate(iface):
    r = A[:, pos[i], :]; X = N[i][0] + r[:, 1]; Y = N[i][1] + r[:, 2]; rr = np.hypot(X, Y); nx, ny = X / rr, Y / rr
    wn[:, c] = r[:, 5] * nx + r[:, 6] * ny; vn[:, c] = r[:, 3] * nx + r[:, 4] * ny
    th[c] = math.atan2(N[i][1], N[i][0]); z[c] = N[i][2]
sl = slice(1, None)
rel_rms = np.sqrt((wn[sl] ** 2).mean()) / np.sqrt((vn[sl] ** 2).mean())
wmax = np.abs(wn[sl]).max() / V0
# area weights on the quarter: trapezoid in theta; both z layers carry equal weight (slab thickness h)
order = np.argsort(th); ths = np.unique(th); wts = {}
dth = np.diff(ths); wt = np.zeros(len(ths)); wt[:-1] += dth / 2; wt[1:] += dth / 2
wmap = {v: k for k, v in enumerate(ths)}
Wt = np.array([wt[wmap[v]] for v in th]) * Ro * (h / 2)       # per node (2 z-layers share h)
net_leak = (wn * Wt).sum(axis=1); net_ring = (vn * Wt).sum(axis=1); abs_ring = (np.abs(vn) * Wt).sum(axis=1)
net_rel = np.abs(net_leak[sl]).max() / abs_ring[sl].max()
c2 = np.cos(2 * th)
n2_leak = (wn * c2 * Wt).sum(axis=1); n2_ring = (vn * c2 * Wt).sum(axis=1)
n2_rel = np.abs(n2_leak[sl]).max() / np.abs(n2_ring[sl]).max()
out = dict(n_interface_nodes=len(iface), rms_leak_over_rms_ring_normal_velocity=float(rel_rms), max_abs_wn_over_V0=float(wmax),
           net_leaked_flux_peak_over_peak_abs_ring_flux=float(net_rel), n2_leaked_flux_peak_over_n2_ring_flux_peak=float(n2_rel),
           net_ring_flux_peak=float(np.abs(net_ring[sl]).max()), net_leak_flux_peak=float(np.abs(net_leak[sl]).max()),
           integrated_leaked_volume_peak=float(np.abs(np.cumsum(net_leak[sl]) * (t[1] - t[0])).max()))
json.dump(out, open(run + '/interface_leakage.json', 'w'), indent=1)
for k, v in out.items(): print(f'  {k}: {v:.4e}' if isinstance(v, float) else f'  {k}: {v}')
