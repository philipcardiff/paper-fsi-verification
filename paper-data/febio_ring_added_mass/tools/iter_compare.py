#!/usr/bin/env python3
"""Compare a tightened-tolerance run with its baseline over the common time window. usage: iter_compare.py <baseline_run> <tight_run>
Quantities: ring ovalling displacement D(t) = (u_r(0)-u_r(90))/2 (analyze_ring.series), max |dD|/max|D| over the window, the time of the
first zero crossing (if present), the n=2 mean-pressure probe, and the number of residual evaluations per step in the window."""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import analyze_ring as AR
b, t_ = sys.argv[1].rstrip('/'), sys.argv[2].rstrip('/')
tb, Db, vb = AR.series(b + '/probes.txt'); tt, Dt, vt = AR.series(t_ + '/probes.txt')
n = min(len(tb), len(tt)); assert np.allclose(tb[:n], tt[:n])
d = Dt[:n] - Db[:n]; sc = np.abs(Db[:n]).max()
print(f'common window: {n - 1} steps, t = 0 .. {tb[n-1]:.5f} s ({tb[n-1] / (2 * np.pi / 1.8159716):.3f} periods)')
print(f'max |D_tight - D_base| / max|D_base| = {np.abs(d).max() / sc:.3e}  (rms {np.sqrt((d**2).mean()) / sc:.3e})')
print(f'max |v_tight - v_base| / max|v_base| = {np.abs(vt[:n] - vb[:n]).max() / np.abs(vb[:n]).max():.3e}')
def zc(t, D):
    for i in range(len(t) - 1):
        if D[i] * D[i + 1] < 0: return t[i] - D[i] * (t[i + 1] - t[i]) / (D[i + 1] - D[i])
zb, zt = zc(tb[:n], Db[:n]), zc(tt[:n], Dt[:n])
if zb and zt: print(f'first zero crossing: base {zb:.8f} tight {zt:.8f}  relative diff {(zt - zb) / zb:+.3e}  (-> omega estimate change {-(zt - zb) / zb:+.3e})')
