#!/usr/bin/env python3
"""Investigate the dilatation-increment noise floor from FEBio logs. usage: dilatation_floor.py <log> [<log> ...]
FEBio prints the SQUARED norms: CURRENT = |increment|^2, REQUIRED = tol^2 * |accumulated step increment|^2.
For each step the script prints the dilatation CURRENT per Newton iteration; the 'floor' is the plateau reached
once the quadratic decay stops (later iterations of the same step, or the last iteration of every step in a normal run).
Output: per-log summary (min/median/max of the final-iteration increment^2 and of accum^2 recovered as REQUIRED/tol^2).
usage with --tol <ftol> sets the tolerance that produced the log (default 1e-6).
"""
import re, sys, numpy as np
args = sys.argv[1:]; tol = 1e-6
if '--tol' in args: i = args.index('--tol'); tol = float(args[i + 1]); del args[i:i + 2]
for p in args:
    steps = []; cur = None
    for l in open(p, errors='ignore'):
        if l.startswith('===== beginning time step'): cur = []; steps.append(cur)
        m = re.match(r'\s+dilatation\s+(\S+)\s+(\S+)\s+(\S+)', l)
        if m and cur is not None: cur.append((float(m.group(2)), float(m.group(3))))
    steps = [s for s in steps if s]
    fin = np.array([s[-1][0] for s in steps]); acc = np.array([s[-1][1] / tol ** 2 for s in steps])
    nit = np.array([len(s) for s in steps])
    # plateau: in steps with >=4 iterations, the last two increments^2
    plat = np.array([s[-1][0] for s in steps if len(s) >= 5])
    print(p)
    print(f'  steps {len(steps)}; iterations/step min/median/max {nit.min()}/{int(np.median(nit))}/{nit.max()}')
    print(f'  final-iteration |increment|^2: min {fin.min():.2e} median {np.median(fin):.2e} max {fin.max():.2e}')
    print(f'  accumulated step |increment|^2 (REQUIRED/tol^2): min {acc.min():.2e} median {np.median(acc):.2e} max {acc.max():.2e}')
    rel = np.sqrt(fin / acc)
    print(f'  final relative increment |inc|/|accum|: min {rel.min():.2e} median {np.median(rel):.2e} max {rel.max():.2e}  (tolerance {tol:g})')
    if len(plat): print(f'  steps with >=5 iterations: {len(plat)}; plateau |inc|^2 median {np.median(plat):.2e}')
