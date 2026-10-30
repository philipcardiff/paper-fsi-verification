#!/usr/bin/env python3
"""Per-step convergence evidence from a verbose FEBio log: for each accepted step, the LAST nonlinear-iteration table;
margin = (achieved norm)/(tolerance norm) for residual, energy, displacement, velocity, dilatation (<1 means satisfied).
usage: convergence_margins.py <log> -> prints summary and writes <log>.margins.json"""
import re, sys, json, numpy as np
log = open(sys.argv[1], errors='ignore').read()
steps = log.split('Nonlinear solution status: time=')[1:]
rows = {}; last = {}
cur_t = None; per_step = {}
for blk in steps:
    t = float(blk.split()[0]); vals = {}
    for key in ('residual', 'energy', 'displacement', 'velocity', 'dilatation'):
        m = re.search(r'\t\s+%s\s+([-\d.eE+]+)\s+([-\d.eE+]+)\s+([-\d.eE+]+)' % key, blk)
        if m: vals[key] = (float(m.group(1)), float(m.group(2)), float(m.group(3)))
    nit = re.search(r'right hand side evaluations\s+=\s+(\d+)', blk)
    per_step[t] = (vals, int(nit.group(1)) if nit else None)      # later blocks of the same time overwrite -> last iteration
out = {}
for key in ('residual', 'energy', 'displacement', 'velocity', 'dilatation'):
    mg = np.array([v[0][key][1] / v[0][key][2] for v in per_step.values() if key in v[0] and v[0][key][2] > 0])
    out[key] = dict(n=int(len(mg)), margin_max=float(mg.max()), margin_median=float(np.median(mg)), margin_min=float(mg.min()))
its = [v[1] for v in per_step.values() if v[1]]
out['iterations_per_step'] = dict(min=min(its), median=float(np.median(its)), max=max(its))
out['n_steps_with_table'] = len(per_step)
json.dump(out, open(sys.argv[1] + '.margins.json', 'w'), indent=1)
print(f"steps {len(per_step)}; iterations/step min/median/max {min(its)}/{np.median(its):.0f}/{max(its)}")
for k in ('residual', 'energy', 'displacement', 'velocity', 'dilatation'):
    o = out[k]; print(f"  {k:13s} achieved/tolerance: worst {o['margin_max']:.2e}, median {o['margin_median']:.2e}, best {o['margin_min']:.2e}  (n={o['n']})")
