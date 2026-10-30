#!/usr/bin/env python3
"""Compare the half-amplitude strong M2 run with the baseline (same estimators, same dry reference). Writes results/ring_amplitude_comparison.csv and figures/ring_lastperiods_amplitude.png."""
import sys, os, json, math, csv
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import analyze_ring as AR
from scipy.optimize import minimize
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import importlib.util
spec = importlib.util.spec_from_file_location('ext', 'tools/ring_lastperiods_ext.py')
src = open('tools/ring_lastperiods_ext.py').read().split('rows = []')[0]   # reuse fit() without running the study
ns = {'__file__': os.path.abspath('tools/ring_lastperiods_ext.py')}; exec(compile(src, 'ext', 'exec'), ns); fit = ns['fit']
REF = json.load(open('reference/s4f_ringAddedMass/verification/reference/ringAddedMass_verification_references.json'))['exact']
RATIO = REF['levels']['strong']['ratio']; OMW = REF['levels']['strong']['omega']
DRY = 'runs/phaseB/D1_dry_f2'; RUNS = {'baseline V0=1.5e-3': 'runs/phaseD/D3_strong_M2', 'half V0=0.75e-3': 'runs/ringNext/R_AMP_strong_M2_V0half'}
td, dd, _ = AR.series(DRY + '/probes.txt'); dm = json.load(open(DRY + '/model_meta.json')); Td = 2 * math.pi / dm['omega_ref']
rows = []; harm = {}
for name, run in RUNS.items():
    t, d, v = AR.series(run + '/probes.txt'); m = json.load(open(run + '/model_meta.json')); T = 2 * math.pi / m['omega_ref']
    for mn, terms in {'damped (original)': (), 'DC+lin+3w': ('dc', 'lin', 'h3')}.items():
        for wn, (a, b) in {'0-4': (0, 4), '1-4': (1, 4), '2-4': (2, 4), '0-3': (0, 3), '1-3': (1, 3), '0-2': (0, 2)}.items():
            ow = fit(t, d, a * T, b * T, m['omega_ref'], terms); od = fit(td, dd, a * Td, b * Td, dm['omega_ref'], terms)
            rows.append(dict(case=name, model=mn, window=wn, omega_wet=ow, omega_wet_sq=ow ** 2, omega_dry=od, ratio_err_pct=100 * (ow / od / RATIO - 1), wet_err_pct=100 * (ow / OMW - 1),
                             d_omega2_vs_exact=ow ** 2 - OMW ** 2))
    rep = json.load(open(run + '/report.json'))
    w = rep['omega_fit']; sdec = -math.log(1 - rep['amp_loss_per_period']) / (2 * math.pi / w); e = np.exp(-sdec * t)
    M = np.column_stack([e * np.sin(w * t), e * np.cos(w * t)] + [f(k * w * t) for k in (2, 3) for f in (np.sin, np.cos)] + [np.ones_like(t)])
    c, *_ = np.linalg.lstsq(M, d, rcond=None); A = math.hypot(c[0], c[1])
    harm[name] = dict(A=A, h2=math.hypot(c[2], c[3]) / A, h3=math.hypot(c[4], c[5]) / A, dc=c[6] / A)
with open('results/ring_amplitude_comparison.csv', 'w', newline='') as fh:
    wr = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); wr.writeheader()
    for r in rows: wr.writerow({k: (f'{x:.12g}' if isinstance(x, float) else x) for k, x in r.items()})
for mn in ('damped (original)', 'DC+lin+3w'):
    for wn in ('0-4', '2-4', '0-2'):
        b = [r for r in rows if r['case'].startswith('baseline') and r['model'] == mn and r['window'] == wn][0]; h = [r for r in rows if r['case'].startswith('half') and r['model'] == mn and r['window'] == wn][0]
        print(mn, wn, f"base {b['ratio_err_pct']:.5f} half {h['ratio_err_pct']:.5f} delta {h['ratio_err_pct']-b['ratio_err_pct']:.5f}  domega2 base {b['d_omega2_vs_exact']:.3e} half {h['d_omega2_vs_exact']:.3e}")
for k, v in harm.items(): print(k, v)
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 2, figsize=(10.5, 4.2))
for j, (mn, col) in enumerate((('damped (original)', '#d95f02'), ('DC+lin+3w', '#1b9e77'))):
    for k, (nm, am) in enumerate((('baseline V0=1.5e-3', 1.94e-3), ('half V0=0.75e-3', 0.968e-3))):
        ys = [r['ratio_err_pct'] for r in rows if r['case'] == nm and r['model'] == mn]
        ax[0].plot([am * 1e3] * len(ys), ys, 'o', color=col, alpha=.7, label=mn if k == 0 else None)
ax[0].set_xlabel('wet displacement amplitude A (mm)'); ax[0].set_ylabel('strong M2 ratio error (%)'); ax[0].set_xlim(0, 2.2); ax[0].set_ylim(-0.005, 0.045)
xx = np.linspace(0, 2.2, 50); ax[0].plot(xx, 0.03987 * (xx / 1.9356) ** 2, 'k--', lw=1, label='c A^2 through baseline'); ax[0].legend(); ax[0].grid(alpha=.3); ax[0].set_title('(a) six windows per case, two fit models')
ax[1].plot([0.611, 0.825, 1.9356], [0.00413, 0.00723, 0.03987], 'o', color='#444', label='baseline weak / moderate / strong M2')
ax[1].plot([0.968], [rows[[r['case'] for r in rows].index('half V0=0.75e-3')]['ratio_err_pct']], 's', color='#d95f02', label='strong M2, V0/2')
xx = np.linspace(0, 2.2, 50); ax[1].plot(xx, 0.03987 * (xx / 1.9356) ** 2, 'k--', lw=1, label='c A^2'); ax[1].set_xlabel('wet displacement amplitude A (mm)'); ax[1].set_ylabel('ratio error (%)')
ax[1].legend(fontsize=8); ax[1].grid(alpha=.3); ax[1].set_title('(b) all regimes against A (floor / A^2 = const)')
plt.tight_layout(); plt.savefig('figures/ring_lastperiods_amplitude.png', dpi=130)
