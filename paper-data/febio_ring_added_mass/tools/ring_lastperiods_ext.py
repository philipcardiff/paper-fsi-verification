#!/usr/bin/env python3
"""Part 1c: sensitivity of omega and ratio error to extra nuisance terms (DC, linear drift, 3rd harmonic) in the fit model. Stored data only."""
import sys, os, json, math, csv
import numpy as np
from scipy.optimize import minimize
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import analyze_ring as AR
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
REF = json.load(open('reference/s4f_ringAddedMass/verification/reference/ringAddedMass_verification_references.json'))['exact']
LV = REF['levels']
DRY = 'runs/phaseB/D1_dry_f2'
WET = {'weak': 'runs/phaseD/D4_weak_M2', 'moderate': 'runs/phaseC/C2_Emesh1e-3_wet_f2', 'strong': 'runs/phaseD/D3_strong_M2'}
MODELS = {'base(2 par)': (), 'DC': ('dc',), 'DC+lin': ('dc', 'lin'), 'DC+lin+3w': ('dc', 'lin', 'h3'), '3w only': ('h3',)}
WIN = {'0-4': (0, 4), '1-4': (1, 4), '2-4': (2, 4), '0-3': (0, 3), '1-3': (1, 3), '0-2': (0, 2)}

def fit(t, d, a, b, w0, terms):
    m = (t >= a) & (t <= b); tt, dd = t[m], d[m]; tm = tt - tt.mean()
    def cols(w, s):
        e = np.exp(-s * tt); c = [e * np.sin(w * tt), e * np.cos(w * tt)]
        if 'dc' in terms: c.append(np.ones_like(tt))
        if 'lin' in terms: c.append(tm)
        if 'h3' in terms: c += [np.sin(3 * w * tt), np.cos(3 * w * tt)]
        return np.column_stack(c)
    def f(p):
        M = cols(p[0] * w0, p[1] * 1e-3); c, *_ = np.linalg.lstsq(M, dd, rcond=None); r = dd - M @ c; return float(r @ r) / float(dd @ dd)
    best = None
    for w in (0.995, 1.0, 1.005):
        r = minimize(f, [w, 0.0], method='Nelder-Mead', options=dict(xatol=1e-12, fatol=1e-30, maxiter=3000))
        if best is None or r.fun < best.fun: best = r
    return best.x[0] * w0

rows = []
dt_ = AR.series(DRY + '/probes.txt'); dm = json.load(open(DRY + '/model_meta.json')); Td = 2 * math.pi / dm['omega_ref']
for lev, run in WET.items():
    tw, dw, _ = AR.series(run + '/probes.txt'); wm = json.load(open(run + '/model_meta.json')); Tw = 2 * math.pi / wm['omega_ref']
    for mn, terms in MODELS.items():
        for wn, (a, b) in WIN.items():
            ow = fit(tw, dw, a * Tw, b * Tw, wm['omega_ref'], terms); od = fit(dt_[0], dt_[1], a * Td, b * Td, dm['omega_ref'], terms)
            rows.append(dict(regime=lev, model=mn, window=wn, omega_wet=ow, omega_dry=od, ratio_err_pct=100 * (ow / od / LV[lev]['ratio'] - 1)))
            print(rows[-1], flush=True)
with open('results/ring_lastperiods_extmodels.csv', 'w', newline='') as fh:
    wr = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); wr.writeheader()
    for r in rows: wr.writerow({k: (f'{x:.10g}' if isinstance(x, float) else x) for k, x in r.items()})
