#!/usr/bin/env python3
"""Figures for the ring final comparison, from RING_*_STUDY.csv only."""
import csv, os, math
import matplotlib.ticker as mt
import numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
os.chdir(os.environ.get('RING_ROOT', os.getcwd())); os.makedirs('figures', exist_ok=True)
rd = lambda p: list(csv.DictReader(open(p)))
fl = lambda x: float(x) if x not in ('', None) else float('nan')
SP, TP, CP = rd('RING_SPATIAL_STUDY.csv'), rd('RING_TEMPORAL_STUDY.csv'), rd('RING_COMPRESSIBILITY_STUDY.csv')
COL = {'weak': '#1b9e77', 'moderate': '#d95f02', 'strong': '#7570b3', 'dry': '0.3'}
plt.rcParams.update({'font.size': 9, 'axes.grid': True, 'grid.alpha': .3})
H = {'M1': 0.45 / 12, 'M2': 0.45 / 24, 'M3': 0.45 / 48}   # fluid radial cell size (m), proxy for h (ring cell: 0.1/ (1,2,4))
hh = {'M1': 0.1, 'M2': 0.05, 'M3': 0.025}
# ---- Fig 1: spatial: absolute errors vs ratio error ----
fig, ax = plt.subplots(1, 2, figsize=(10, 3.9))
for lev in ('weak', 'moderate', 'strong'):
    rs = [r for r in SP if r['regime'] == lev]; x = [hh[r['mesh']] for r in rs]
    ax[0].plot(x, [fl(r['err_wet_pct']) for r in rs], 'o-', color=COL[lev], label=f'{lev} wet $\\omega$')
    ax[1].plot(x, [fl(r['ratio_err_pct']) for r in rs], 'o-', color=COL[lev], label=f'FEBio {lev}')
d = [r for r in SP if r['regime'] == 'dry']; ax[0].plot([hh[r['mesh']] for r in d], [fl(r['err_pct']) for r in d], 's--', color=COL['dry'], label='dry $\\omega$')
for lev, mk in (('moderate', '^'), ('strong', 'v')):
    for key, ls, nm in (('s4f_standard_iqnils_ratio_err_pct', 'x', 'solids4foam IQN-ILS (standard solid, PR #13)'), ('s4f_standard_monolithic_ratio_err_pct', '+', 'solids4foam quasi-monolithic (standard solid, PR #13)')):
        rs = [r for r in SP if r['regime'] == lev and r[key] != '']
        ax[1].plot([hh[r['mesh']] for r in rs], [fl(r[key]) for r in rs], ls, ms=7, color=COL[lev], alpha=.8, label=nm if lev == 'moderate' else None)
ax[0].set_ylabel('error of $\\omega$ vs exact (%)'); ax[1].set_ylabel('wet/dry ratio error vs exact (%)')
for a in ax: a.set_xscale('log'); a.set_xticks([0.1, 0.05, 0.025]); a.set_xticklabels(['0.1\nM1', '0.05\nM2', '0.025\nM3']); a.xaxis.set_minor_formatter(mt.NullFormatter()); a.set_xlabel('ring cell size $h$ (m)'); a.invert_xaxis()
ax[1].axhline(0.1, color='k', lw=.8, ls=':'); ax[1].text(0.099, 0.103, '0.1 % band', fontsize=7)
ax[0].legend(fontsize=7); ax[1].legend(fontsize=6.5, loc='upper right', ncol=1)
ax[0].set_title('(a) absolute errors (hex8 solid: spatial error dominates)'); ax[1].set_title('(b) ratio error: offset, not mesh-converging to 0')
plt.tight_layout(); plt.savefig('figures/ring_fig1_spatial.png', dpi=140); plt.close()
# ---- Fig 2: temporal ----
fig, ax = plt.subplots(1, 2, figsize=(10, 3.9))
dr = [r for r in TP if r['regime'] == 'dry']
ax[0].plot([int(r['steps_per_period']) for r in dr], [fl(r['err_pct']) for r in dr], 's--', color=COL['dry'], label='dry')
for lev in ('weak', 'moderate', 'strong'):
    rs = [r for r in TP if r['regime'] == lev]; ax[0].plot([int(r['steps_per_period']) for r in rs], [fl(r['err_wet_pct']) for r in rs], 'o-', color=COL[lev], label=lev)
    ax[1].plot([int(r['steps_per_period']) for r in rs], [fl(r['ratio_err_pct']) for r in rs], 'o-', color=COL[lev], label=f'FEBio {lev}')
    if lev in ('moderate', 'strong'):
        ax[1].plot([int(r['steps_per_period']) for r in rs], [fl(r['s4f_standard_iqnils_ratio_err_pct']) for r in rs], 'x--', color=COL[lev], alpha=.8, label='s4f IQN-ILS (std solid)' if lev == 'moderate' else None)
        ax[1].plot([int(r['steps_per_period']) for r in rs], [fl(r['s4f_standard_monolithic_ratio_err_pct']) for r in rs], '+:', color=COL[lev], ms=8, alpha=.8, label='s4f quasi-monolithic (std solid)' if lev == 'moderate' else None)
for a in ax: a.set_xscale('log'); a.set_xticks([25, 50, 100]); a.set_xticklabels(['25', '50', '100']); a.xaxis.set_minor_formatter(mt.NullFormatter()); a.set_xlabel('time steps per period (M2)')
ax[0].set_ylabel('$\\omega$ error vs exact (%)'); ax[1].set_ylabel('wet/dry ratio error vs exact (%)'); ax[0].legend(fontsize=7)
ax[1].legend(fontsize=6.5); ax[0].set_title('(a) $\\omega$: 2nd-order time error (order 1.99) on top of the spatial offset'); ax[1].set_title('(b) ratio error: no dt trend in FEBio')
ax[0].title.set_fontsize(8)
plt.tight_layout(); plt.savefig('figures/ring_fig2_temporal.png', dpi=140); plt.close()
# ---- Fig 3: plateau diagnostic: equivalent additive d(omega^2) ----
fig, ax = plt.subplots(1, 2, figsize=(10, 3.9))
for lev in ('weak', 'moderate', 'strong'):
    rs = [r for r in SP if r['regime'] == lev]
    ax[0].plot([hh[r['mesh']] for r in rs], [fl(r['equiv_additive_domega2_rad2_s2']) * 1e4 for r in rs], 'o-', color=COL[lev], label=lev)
ax[0].set_xscale('log'); ax[0].set_xticks([0.1, 0.05, 0.025]); ax[0].set_xticklabels(['0.1\nM1', '0.05\nM2', '0.025\nM3']); ax[0].xaxis.set_minor_formatter(mt.NullFormatter()); ax[0].invert_xaxis(); ax[0].set_xlabel('ring cell size $h$ (m)'); ax[0].set_ylabel('$2\\,\\omega_{exact}^2\\,\\epsilon_{ratio}$  ($10^{-4}$ rad$^2$/s$^2$)'); ax[0].set_ylim(0, 7); ax[0].legend(fontsize=8)
ax[0].set_title('(a) equivalent additive shift of $\\omega_{wet}^2$ (all regimes, all meshes)')
x = [1 / (2 * fl(r['omega_wet_exact']) ** 2) for r in SP if r['regime'] != 'dry' and r['mesh'] == 'M2']; y = [fl(r['ratio_err_pct']) / 100 for r in SP if r['regime'] != 'dry' and r['mesh'] == 'M2']
ax[1].plot(x, y, 'o', color='k'); xx = np.linspace(0, max(x) * 1.05, 5); k = np.mean([b / a for a, b in zip(x, y)]); ax[1].plot(xx, k * xx, '-', color='0.5', label=f'$\\epsilon = \\Delta/(2\\omega^2)$, $\\Delta$={k:.2e}')
for lev, xi, yi in zip(('weak', 'moderate', 'strong'), x, y): ax[1].annotate(lev, (xi, yi), textcoords='offset points', xytext=(5, -10), fontsize=8)
ax[1].set_xlabel('$1/(2\\omega_{exact}^2)$  (s$^2$/rad$^2$)'); ax[1].set_ylabel('ratio error $\\epsilon$ at M2 (fraction)'); ax[1].legend(fontsize=8); ax[1].set_title('(b) plateau scales as $1/\\omega_{wet}^2$ (descriptive, cause not identified)'); ax[1].title.set_fontsize(8)
plt.tight_layout(); plt.savefig('figures/ring_fig3_plateau_diagnostic.png', dpi=140); plt.close()
# ---- Fig 4: compressibility ----
fig, ax = plt.subplots(1, 2, figsize=(10, 3.9))
pairs = [(r['regime'] + ' M1\nK ' + ('1e8→1e9' if r['regime'] == 'strong' else '1e8→1e7'), fl(r['measured_domega_over_omega_vs_first_row']), fl(r['theory_domega_over_omega']), fl(r['measurement_floor_window_spread_rel'])) for r in CP if r.get('measured_domega_over_omega_vs_first_row', '') != '']
xs = np.arange(len(pairs)); ax[0].bar(xs - .18, [p[1] * 1e6 for p in pairs], .36, label='FEBio measured', color='C0'); ax[0].bar(xs + .18, [p[2] * 1e6 for p in pairs], .36, label='linear acoustics theory', color='C1')
ax[0].errorbar(xs - .18, [p[1] * 1e6 for p in pairs], yerr=[p[3] * 1e6 for p in pairs], fmt='none', color='k', capsize=3); ax[0].set_xticks(xs); ax[0].set_xticklabels([p[0] for p in pairs]); ax[0].set_ylabel('$\\Delta\\omega/\\omega$ for the K change ($10^{-6}$)'); ax[0].legend(fontsize=8)
ax[0].set_title('(a) K effect on $\\omega$: measured vs theory (error bar: fit-window spread)'); ax[0].title.set_fontsize(8)
s = [r for r in CP if r['regime'] == 'strong' and r['run'].endswith(('M1', 'K1e9')) and r['ratio_err_pct'] != '']
ax[1].plot([fl(r['K']) for r in s], [fl(r['ratio_err_pct']) for r in s], 'o-', color=COL['strong']); ax[1].set_xscale('log'); ax[1].set_ylim(0, 0.05); ax[1].axhline(0.0355, color='0.6', ls=':')
ax[1].set_xlabel('K (Pa)'); ax[1].set_ylabel('strong M1 ratio error (%)'); ax[1].set_title('(b) the strong-case plateau does not move with K'); ax[1].title.set_fontsize(8)
plt.tight_layout(); plt.savefig('figures/ring_fig4_compressibility.png', dpi=140); plt.close()
# ---- Fig 5: cost ----
cs = rd('RING_COMPUTATIONAL_COST.csv'); fig, ax = plt.subplots(1, 2, figsize=(10, 3.6))
wet = [r for r in cs if r['regime'] != 'dry' and r['steps_per_period'] == '100' and r['group'] in ('phaseD', 'phaseC', 'ringNext') and 'loose' not in r['run'] and 'K1e9' not in r['run']]
ax[0].scatter([fl(r['mesh_nodes']) for r in wet], [fl(r['allocated_cpu_h']) for r in wet], c=[COL[r['regime']] for r in wet]); ax[0].set_xscale('log'); ax[0].set_yscale('log'); ax[0].set_xlabel('mesh nodes'); ax[0].set_ylabel('allocated CPU-hours')
for lev in COL:
    if lev != 'dry': ax[0].scatter([], [], c=COL[lev], label=lev)
[ax[0].annotate(f"{r['cores']}c", (fl(r['mesh_nodes']), fl(r['allocated_cpu_h'])), textcoords='offset points', xytext=(6, -3), fontsize=7) for r in wet]
ax[0].legend(fontsize=8); ax[0].set_title('(a) wet runs, 100 steps/period; label = cores used'); ax[0].title.set_fontsize(8)
ax[1].scatter([fl(r['mesh_nodes']) for r in wet], [fl(r['rss_MB']) for r in wet], c=[COL[r['regime']] for r in wet]); ax[1].set_xscale('log'); ax[1].set_yscale('log'); ax[1].set_xlabel('mesh nodes'); ax[1].set_ylabel('peak RSS (MB)'); ax[1].set_title('(b) memory')
plt.tight_layout(); plt.savefig('figures/ring_fig5_cost.png', dpi=140); plt.close()
print('figures done')
