#!/usr/bin/env python3
"""Part 1 analysis (post-processing only, stored data): window / estimator dependence of the ring frequency and wet/dry ratio.
Uses tools/analyze_ring.py's ORIGINAL estimator (fit_damped: nested-grid variable projection of d=exp(-s t)(a sin wt + b cos wt))
unchanged, plus a fast scipy implementation of the same model that is cross-checked against it. Writes results/ring_lastperiods_*.csv."""
import sys, os, json, math, csv
import numpy as np
from scipy.optimize import minimize
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import analyze_ring as AR
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); os.chdir(ROOT)
REF = json.load(open('reference/s4f_ringAddedMass/verification/reference/ringAddedMass_verification_references.json'))['exact']
OM_DRY = REF['dry']['omega']; LV = REF['levels']
RHO = {'weak': 14, 'moderate': 140, 'strong': 1400}
RUNS = {  # (wet, dry, steps/period, mesh)
 ('weak', 'M2', 100): ('runs/phaseD/D4_weak_M2', 'runs/phaseB/D1_dry_f2'),
 ('moderate', 'M2', 100): ('runs/phaseC/C2_Emesh1e-3_wet_f2', 'runs/phaseB/D1_dry_f2'),
 ('strong', 'M2', 100): ('runs/phaseD/D3_strong_M2', 'runs/phaseB/D1_dry_f2'),
 ('strong', 'M2', 50): ('runs/ringNext/R_TIME_strong_M2_s50', 'runs/ringNext/R_TIME_dry_M2_s50'),
 ('strong', 'M2', 25): ('runs/ringNext/R_TIME_strong_M2_s25', 'runs/ringNext/R_TIME_dry_M2_s25'),
 ('moderate', 'M2', 50): ('runs/ringNext/R_TIME_moderate_M2_s50', 'runs/ringNext/R_TIME_dry_M2_s50'),
 ('weak', 'M2', 50): ('runs/ringNext/R_TIME_weak_M2_s50', 'runs/ringNext/R_TIME_dry_M2_s50'),
 ('strong', 'M1', 100): ('runs/phaseD/D2_strong_M1', 'runs/phaseB/B1_dry_hex8_f1'),
 ('weak', 'M1', 100): ('runs/phaseD/D1_weak_M1', 'runs/phaseB/B1_dry_hex8_f1'),
}

def load(run):
    t, d, v = AR.series(run + '/probes.txt'); return t, d, v

def vp(tt, dd, w, s):
    e = np.exp(-s * tt); M = np.column_stack([e * np.sin(w * tt), e * np.cos(w * tt)])
    c, *_ = np.linalg.lstsq(M, dd, rcond=None); r = dd - M @ c
    return float(r @ r), c

def fit_fast(t, d, a, b, w0):
    m = (t >= a) & (t <= b); tt, dd = t[m], d[m]
    f = lambda p: vp(tt, dd, p[0] * w0, p[1] * 1e-3)[0] / float(dd @ dd)
    # multi-start (a coarse scan avoids local minima in short windows)
    best = None
    for w in np.linspace(0.97, 1.03, 13):
        r = minimize(f, [w, 0.0], method='Nelder-Mead', options=dict(xatol=1e-12, fatol=1e-30, maxiter=4000, maxfev=8000))
        if best is None or r.fun < best.fun: best = r
    w, s = best.x[0] * w0, best.x[1] * 1e-3
    _, c = vp(tt, dd, w, s)
    return dict(omega=w, decay=s, amp0=math.hypot(*c), rms=math.sqrt(best.fun * float(dd @ dd) / len(tt)), rel_rms=math.sqrt(best.fun), n=len(tt))

def fit_undamped(t, d, a, b, w0):
    m = (t >= a) & (t <= b); tt, dd = t[m], d[m]
    f = lambda p: vp(tt, dd, p[0] * w0, 0.0)[0] / float(dd @ dd)
    best = None
    for w in np.linspace(0.97, 1.03, 13):
        r = minimize(f, [w], method='Nelder-Mead', options=dict(xatol=1e-13, fatol=1e-30, maxiter=4000))
        if best is None or r.fun < best.fun: best = r
    return best.x[0] * w0

def peaks_omega(t, d, a, b):
    """times of extrema by parabolic interpolation of d(t) (3 pts); omega from slope of extremum index vs time (half period)."""
    m = np.where((t >= a) & (t <= b))[0]; ex = []
    for i in m:
        if 0 < i < len(t) - 1 and ((d[i] - d[i - 1]) * (d[i + 1] - d[i]) < 0):
            y0, y1, y2 = d[i - 1], d[i], d[i + 1]; den = y0 - 2 * y1 + y2
            ex.append(t[i] + 0.5 * (t[i + 1] - t[i - 1]) / 2 * (y0 - y2) / den * 1.0 if False else t[i] + (t[i + 1] - t[i]) * 0.5 * (y0 - y2) / den)
    if len(ex) < 3: return float('nan')
    return math.pi / np.polyfit(np.arange(len(ex)), ex, 1)[0]

def windows():
    W = {}
    for a in (0, 1, 2, 3): W[f'{a}-4'] = (a, 4)
    W['0-3'] = (0, 3); W['0-2'] = (0, 2); W['1-3'] = (1, 3); W['2-3'] = (2, 3); W['0-1'] = (0, 1); W['1-2'] = (1, 2); W['3-4'] = (3, 4)
    W['0.5-4'] = (0.5, 4); W['1.5-4'] = (1.5, 4); W['2.5-4'] = (2.5, 4)
    return W

def main():
    # 1) reproduction of the stored numbers with the ORIGINAL estimator
    rep_rows = []
    cache = {}
    for key, (wet, dry) in RUNS.items():
        lev, mesh, spp = key
        for run, om_ref in ((wet, LV[lev]['omega']), (dry, OM_DRY)):
            if run in cache: continue
            t, d, v = load(run); cache[run] = (t, d, v)
            meta = json.load(open(run + '/model_meta.json')); rep = json.load(open(run + '/report.json'))
            T = 2 * math.pi / meta['omega_ref']
            f_orig = AR.fit_damped(t, d, 0, 4 * T, meta['omega_ref'])
            ff = fit_fast(t, d, 0, 4 * T, meta['omega_ref'])
            rep_rows.append(dict(run=os.path.basename(run), omega_stored=rep['omega_fit'], omega_original_estimator=f_orig['omega'],
                                 rel_diff_orig_vs_stored=f_orig['omega'] / rep['omega_fit'] - 1, omega_fast=ff['omega'], rel_diff_fast_vs_stored=ff['omega'] / rep['omega_fit'] - 1,
                                 decay_orig=f_orig['decay_rate'], decay_fast=ff['decay'], T_ref=T, n_samples=len(t), t_end_over_T=float(t[-1]) / T))
    with open('results/ring_lastperiods_reproduction.csv', 'w', newline='') as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rep_rows[0].keys())); wr.writeheader(); wr.writerows(rep_rows)
    # 2) window study
    out = []
    for key, (wet, dry) in RUNS.items():
        lev, mesh, spp = key
        res = {}
        for tag, run, om_ref in (('wet', wet, LV[lev]['omega']), ('dry', dry, OM_DRY)):
            t, d, v = cache[run]; meta = json.load(open(run + '/model_meta.json')); T = 2 * math.pi / meta['omega_ref']
            for wn, (a, b) in windows().items():
                if b * T > t[-1] + 1e-6 * T: continue
                f = fit_fast(t, d, a * T, b * T, meta['omega_ref'])
                res[(tag, wn)] = dict(f, undamped=fit_undamped(t, d, a * T, b * T, meta['omega_ref']),
                                      zc=AR.zero_crossing_omega(t, d, a * T, b * T)[0], pk=peaks_omega(t, d, a * T, b * T))
        ex_r = LV[lev]['ratio']
        for wn in windows():
            if ('wet', wn) not in res: continue
            w, dr = res[('wet', wn)], res[('dry', wn)]
            row = dict(regime=lev, mesh=mesh, steps_per_period=spp, window_periods=wn,
                       omega_wet=w['omega'], omega_dry=dr['omega'], ratio=w['omega'] / dr['omega'],
                       ratio_err_pct=100 * (w['omega'] / dr['omega'] / ex_r - 1),
                       wet_err_pct=100 * (w['omega'] / LV[lev]['omega'] - 1), dry_err_pct=100 * (dr['omega'] / OM_DRY - 1),
                       decay_wet=w['decay'], decay_dry=dr['decay'], amp0_wet=w['amp0'], rel_rms_wet=w['rel_rms'], rel_rms_dry=dr['rel_rms'],
                       ratio_err_undamped_pct=100 * (w['undamped'] / dr['undamped'] / ex_r - 1),
                       ratio_err_zc_pct=100 * (w['zc'] / dr['zc'] / ex_r - 1) if w['zc'] == w['zc'] and dr['zc'] == dr['zc'] else float('nan'),
                       ratio_err_peaks_pct=100 * (w['pk'] / dr['pk'] / ex_r - 1) if w['pk'] == w['pk'] and dr['pk'] == dr['pk'] else float('nan'))
            out.append(row)
    keys = list(out[0].keys())
    with open('results/ring_lastperiods_windows.csv', 'w', newline='') as fh:
        wr = csv.writer(fh); wr.writerow(keys)
        for r in out: wr.writerow([f'{r[k]:.10g}' if isinstance(r[k], float) else r[k] for k in keys])
    print('wrote', len(out), 'rows')
    for r in rep_rows: print(r['run'], f"{r['rel_diff_orig_vs_stored']:.2e} {r['rel_diff_fast_vs_stored']:.2e}")

if __name__ == '__main__':
    main()
