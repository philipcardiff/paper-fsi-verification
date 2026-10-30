#!/usr/bin/env python3
import os, csv, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
W = [r for r in csv.DictReader(open('results/ring_lastperiods_windows.csv')) if r['mesh'] == 'M2' and r['steps_per_period'] == '100']
X = [r for r in csv.DictReader(open('results/ring_lastperiods_extmodels.csv'))]
C = {'weak': '#1b9e77', 'moderate': '#7570b3', 'strong': '#d95f02'}
starts = {'0-4': 0, '0.5-4': .5, '1-4': 1, '1.5-4': 1.5, '2-4': 2, '2.5-4': 2.5, '3-4': 3}
fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
for lev in C:
    pts = [(starts[r['window_periods']], float(r['ratio_err_pct'])) for r in W if r['regime'] == lev and r['window_periods'] in starts]
    pts.sort(); xs, ys = zip(*pts); ax[0].plot(xs, np.array(ys) / ys[0], 'o-', color=C[lev], label=lev)
ax[0].set_xlabel('window start (periods), window ends at 4 periods'); ax[0].set_ylabel('ratio error / ratio error of the 4-period fit')
ax[0].set_title('(a) window-start dependence (damped-sinusoid fit)'); ax[0].legend(); ax[0].grid(alpha=.3)
mods = ['base(2 par)', 'DC', 'DC+lin', 'DC+lin+3w']; wins = ['0-4', '1-4', '2-4', '0-3', '1-3', '0-2']
for k, lev in enumerate(C):
    for j, m in enumerate(mods):
        ys = [float(r['ratio_err_pct']) for r in X if r['regime'] == lev and r['model'] == m]
        if lev == 'strong':
            ax[1].plot([j] * len(ys), ys, 'o', color=C[lev], alpha=.7)
ax[1].set_xticks(range(len(mods))); ax[1].set_xticklabels(mods, rotation=15); ax[1].set_ylabel('strong M2 ratio error (%)')
ax[1].set_title('(b) strong M2, six windows per fit model'); ax[1].grid(alpha=.3)
# (c) resolution vs discrimination (log scale, percentage points of the strong M2 ratio error)
strong = [float(r['ratio_err_pct']) for r in W if r['regime'] == 'strong' and r['window_periods'] in ('0-4', '1-4', '2-4', '0-3', '1-3', '0-2')]
ext = [float(r['ratio_err_pct']) for r in X if r['regime'] == 'strong']
last2 = float([r for r in W if r['regime'] == 'strong' and r['window_periods'] == '2-4'][0]['ratio_err_pct']) - float([r for r in W if r['regime'] == 'strong' and r['window_periods'] == '0-4'][0]['ratio_err_pct'])
bars = [('window spread\n(6 windows >= 2 T)', max(strong) - min(strong)), ('last 2 T minus\n4 T fit', abs(last2)), ('fit-model spread\n(base..+3w, 4 T)', 0.03987 - 0.03979),
        ('registered sigma', 0.0003), ('10 % A^2 share\n(delta)', 0.1 * 0.75 * 0.03987), ('full A^2 delta\n(0.75 x floor)', 0.75 * 0.03987)]
cols = ['#999999'] * 3 + ['#444444', '#d95f02', '#d95f02']
ax[2].barh(range(len(bars)), [b[1] for b in bars], color=cols); ax[2].set_xscale('log'); ax[2].set_yticks(range(len(bars))); ax[2].set_yticklabels([b[0] for b in bars], fontsize=8)
ax[2].set_xlabel('percentage points of the strong M2 ratio error'); ax[2].set_title('(c) fit uncertainty vs change'); ax[2].grid(alpha=.3, which='both', axis='x')
plt.tight_layout(); plt.savefig('figures/ring_lastperiods_windows.png', dpi=130)
# diagnostics figure
Z = np.load('results/ring_lastperiods_resid.npz')
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
for n, c in (('dry_M2', 'k'), ('weak_M2', C['weak']), ('moderate_M2', C['moderate']), ('strong_M2', C['strong'])):
    ax[0].plot(Z[n + '_t'] / Z[n + '_t'][-1] * 4, Z[n + '_r'], color=c, label=n, lw=1)
    ax[1].semilogy(Z[n + '_f'], Z[n + '_amp'], color=c, lw=1, label=n)
ax[0].set_xlabel('time (periods)'); ax[0].set_ylabel('residual of the damped-sinusoid fit / amplitude'); ax[0].legend(); ax[0].grid(alpha=.3)
ax[1].set_xlabel('frequency / omega_fit'); ax[1].set_ylabel('best single-sinusoid amplitude in residual / A'); ax[1].grid(alpha=.3)
ax[0].set_title('(a) fit residual'); ax[1].set_title('(b) residual spectrum (least squares scan)')
plt.tight_layout(); plt.savefig('figures/ring_lastperiods_residual.png', dpi=130)
