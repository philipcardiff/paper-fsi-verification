#!/usr/bin/env python3
"""Angular-mode decomposition of the first fluid layer pressure: usage pressure_modes.py <run_dir> <omega> [nt] [nr]"""
import sys, math, numpy as np
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__))); import analyze_ring as AR
run, w = sys.argv[1], float(sys.argv[2]); NT = int(sys.argv[3]) if len(sys.argv) > 3 else 48; NR = int(sys.argv[4]) if len(sys.argv) > 4 else 12
dth = (math.pi / 2) / NT; th = (np.arange(NT) + 0.5) * dth
b = AR.read_blocks(run + '/fluid_elems.txt'); t = np.array([x['t'] for x in b]); dt = t[1] - t[0]
P = np.array([[r[1] for r in x['rows']] for x in b]).reshape(len(t), NT, NR); L = P[:, :, 0]
out = {}
for name, x in (('n0', L.mean(axis=1)), ('n2', (L * np.cos(2 * th)).sum(axis=1) * dth / (math.pi / 4)), ('n4', (L * np.cos(4 * th)).sum(axis=1) * dth / (math.pi / 4))):
    M = np.column_stack([np.sin(w * t), np.cos(w * t), np.sin(2 * w * t), np.cos(2 * w * t), np.ones_like(t)]); c, *_ = np.linalg.lstsq(M, x, rcond=None); r = x - M @ c
    F = np.abs(np.fft.rfft((r - r.mean()) * np.hanning(len(r)))) ** 2; fr = np.fft.rfftfreq(len(r), dt) * 2 * math.pi
    out[name] = dict(amp_w=math.hypot(c[0], c[1]), amp_2w=math.hypot(c[2], c[3]), resid_rms=float(r.std()), hf_fraction=float(F[fr > 20].sum() / F.sum()))
    print(f"{name}: amp@w {out[name]['amp_w']:.4f} Pa, amp@2w {out[name]['amp_2w']:.4f} Pa, residual rms {out[name]['resid_rms']:.4f} Pa ({100*out[name]['hf_fraction']:.0f}% of residual power > 20 rad/s)")
import json; json.dump(out, open(run + '/pressure_modes.json', 'w'), indent=1)
