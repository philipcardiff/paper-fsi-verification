#!/usr/bin/env python3
"""Phase C results builder: CSVs + figures from runs/*/report.json (+ pressure_modes.json, interface_leakage.json). No hand-typed values."""
import json, os, csv, math, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OM = {'dry': 2.55504219581415, 'wet': 1.8045131724003833}; RATIO_REF = 0.7062557226478153   # exact values from ringAddedMass_verification_references.json (earlier version used 8-digit roundings, relative difference 1.5e-7)
def L(path):
    p = f'{path}/report.json'
    return json.load(open(p)) if os.path.exists(p) else None
def J(path, name):
    p = f'{path}/{name}'
    return json.load(open(p)) if os.path.exists(p) else None

# registry: (label, family, mesh factor, dry run, wet run)
slip = [('M1', 1, 'runs/phaseB/B1_dry_hex8_f1', 'runs/phaseC/C0_wet_f1_wallslip_ifslip'),
        ('M2', 2, 'runs/phaseB/D1_dry_f2', 'runs/phaseC/C1_M2_wet_f2'),
        ('M3', 4, 'runs/phaseC/C1_M3_dry_f4', 'runs/phaseC/C1_M3_wet_f4')]
noslip = [('M1', 1, 'runs/phaseB/B1_dry_hex8_f1', 'runs/phaseB/R1_wet_f1_K1e8_icv2'),
          ('M2', 2, 'runs/phaseB/D1_dry_f2', 'runs/phaseB/D2_wet_f2_K1e8_icv2')]
wallslip_only = [('M1', 1, 'runs/phaseB/B1_dry_hex8_f1', 'runs/phaseB/D5_wet_f1_wallslip')]
fam = {'slip-wall + slip-interface (inviscid-limit formulation)': slip, 'no-slip wall + no-slip interface (R1, D2)': noslip, 'slip wall + no-slip interface (D5)': wallslip_only}

rows = []
for fname, lst in fam.items():
    for lab, f, dry, wet in lst:
        d, w = L(dry), L(wet)
        if d is None and w is None: continue
        h_ring = 0.1 / (4 * f)
        r = dict(family=fname, level=lab, mesh_factor=f, h_ring_radial_m=h_ring, h_fluid_radial_m=0.45 / (12 * f), dry_run=os.path.basename(dry), wet_run=os.path.basename(wet))
        if d: r.update(omega_dry=d['omega_fit'], err_dry_pct=d['err_fit_pct'], dry_nodes=d['n_nodes'], dry_wall_s=d['wall_s'], dry_dt=d['dt'])
        if w:
            r.update(omega_wet=w['omega_fit'], err_wet_pct=w['err_fit_pct'], wet_nodes=w['n_nodes'], wet_wall_s=w['wall_s'], wet_rss_MB=w['max_rss_MB'], wet_dt=w['dt'],
                     mu_eff=w.get('mu_eff'), ma_l0=w.get('ma_eff_ratio_layer0'), ma_l1=w.get('ma_eff_ratio_layer1'), n2_force_err_pct=100 * (w['n2_force_N'] / w['n2_force_analytic_N'] - 1),
                     p_noise_rel=w.get('p_noise_over_analytic'), n0_amp_w=(J(wet, 'pressure_modes.json') or {}).get('n0', {}).get('amp_w'),
                     leak_rms=(J(wet, 'interface_leakage.json') or {}).get('rms_leak_over_rms_ring_normal_velocity'),
                     leak_max=(J(wet, 'interface_leakage.json') or {}).get('max_abs_wn_over_V0'), KE_integer_periods=w['KE_total_at_zero_crossings_over_initial'][2::2],
                     vol_meanJ=w['volume_meanJm1_max'], amp_loss=w['amp_loss_per_period'])
        if d and w:
            r['ratio'] = w['omega_fit'] / d['omega_fit']; r['ratio_err_pct'] = 100 * (r['ratio'] / RATIO_REF - 1)
        rows.append(r)

os.makedirs('results', exist_ok=True); os.makedirs('figures', exist_ok=True)
def write(name, rows_, keys):
    with open(f'results/{name}', 'w', newline='') as fh:
        w = csv.writer(fh); w.writerow(keys)
        for r in rows_: w.writerow([('' if r.get(k) is None else (f'{r[k]:.8g}' if isinstance(r[k], float) else (json.dumps(r[k]) if isinstance(r[k], list) else r[k]))) for k in keys])
write('phaseC_mesh_convergence.csv', rows, ['family', 'level', 'mesh_factor', 'h_ring_radial_m', 'h_fluid_radial_m', 'dry_run', 'wet_run', 'omega_dry', 'err_dry_pct', 'omega_wet', 'err_wet_pct', 'ratio', 'ratio_err_pct', 'mu_eff'])
# exact radial profile of the n=2 potential-flow pressure at the fluid element-centre radii (layer i centre = Ro + (i+1/2) dr)
A_, B_ = -0.00022568759047243063, -0.0011425434267666801
prof = lambda r: (A_ * r * r + B_ / (r * r)) / (A_ * 1.05 ** 2 + B_ / 1.05 ** 2)
for r in rows:
    if r.get('ma_l0') is None: continue
    dr = 0.45 / (12 * r['mesh_factor']); r['profile_l0'] = prof(1.05 + 0.5 * dr); r['profile_l1'] = prof(1.05 + 1.5 * dr)
    r['meas_over_profile_l0'] = r['ma_l0'] / r['profile_l0']; r['meas_over_profile_l1'] = r['ma_l1'] / r['profile_l1']
    # exact layer-0 profile factor implies this force error for a perfect solution
    r['exact_force_err_pct_layer0'] = 100 * (prof(1.05 + 0.5 * dr) - 1)
write('phaseC_pressure_traction.csv', rows, ['family', 'level', 'wet_run', 'ma_l0', 'ma_l1', 'profile_l0', 'profile_l1', 'meas_over_profile_l0', 'meas_over_profile_l1', 'n2_force_err_pct', 'exact_force_err_pct_layer0', 'n0_amp_w', 'p_noise_rel', 'leak_rms', 'leak_max', 'KE_integer_periods', 'vol_meanJ', 'amp_loss'])
write('phaseC_costs.csv', rows, ['family', 'level', 'dry_run', 'dry_nodes', 'dry_wall_s', 'wet_run', 'wet_nodes', 'wet_wall_s', 'wet_rss_MB'])

# sensitivity table (coarse mesh, slip family baseline C0)
base = L('runs/phaseC/C0_wet_f1_wallslip_ifslip'); dryb = L('runs/phaseB/B1_dry_hex8_f1'); sens = []
def ratio_of(w, dry): return w['omega_fit'] / dry['omega_fit']
if base:
    rb = ratio_of(base, dryb)
    for lab, wetp, dryp in (('time step halved (200 steps/period), paired dry+wet', 'runs/phaseC/C1_dt200_wet_f1', 'runs/phaseC/C1_dt200_dry_f1'), ('K 1e8 -> 1e7', 'runs/phaseC/C1_K1e7_wet_f1', None), ('mesh-regulariser E_mesh 1 -> 1e-3 Pa', 'runs/phaseC/C2_Emesh1e-3_wet_f1', None)):
        w = L(wetp); d = L(dryp) if dryp else dryb
        if w and d:
            r = ratio_of(w, d)
            sens.append(dict(effect=lab, d_omega_wet_pct=100 * (w['omega_fit'] / base['omega_fit'] - 1), d_omega_dry_pct=(100 * (d['omega_fit'] / dryb['omega_fit'] - 1)) if dryp else 0.0,
                             d_ratio_pct=100 * (r / rb - 1), ratio_err_pct=100 * (r / RATIO_REF - 1)))
    mm = [r for r in rows if r['family'].startswith('slip-wall') and r['level'] == 'M2']
    if mm and 'ratio' in mm[0]: sens.append(dict(effect='mesh M1 -> M2 (spatial)', d_omega_wet_pct=100 * (mm[0]['omega_wet'] / base['omega_fit'] - 1), d_omega_dry_pct=100 * (mm[0]['omega_dry'] / dryb['omega_fit'] - 1),
                                                   d_ratio_pct=100 * (mm[0]['ratio'] / rb - 1), ratio_err_pct=mm[0]['ratio_err_pct']))
write('phaseC_time_iter_sensitivity.csv', sens, ['effect', 'd_omega_wet_pct', 'd_omega_dry_pct', 'd_ratio_pct', 'ratio_err_pct'])

# ---------------- figures
col = {'slip-wall + slip-interface (inviscid-limit formulation)': 'C0', 'no-slip wall + no-slip interface (R1, D2)': 'C3', 'slip wall + no-slip interface (D5)': 'C2'}
mk = {'slip-wall + slip-interface (inviscid-limit formulation)': 'o', 'no-slip wall + no-slip interface (R1, D2)': 's', 'slip wall + no-slip interface (D5)': '^'}
def series(f, key):
    rr = [r for r in rows if r['family'] == f and r.get(key) is not None]; return [r['h_ring_radial_m'] for r in rr], [r[key] for r in rr]
plt.rcParams.update({'font.size': 9, 'figure.dpi': 130})
fig, ax = plt.subplots(1, 2, figsize=(9, 3.4))
for f in fam:
    h, y = series(f, 'omega_wet'); ax[1].plot(h, y, mk[f] + '-', color=col[f], label=f)
h, y = series(list(fam)[0], 'omega_dry'); ax[0].plot(h, y, 'ko-', label='dry hex8 ring (matched across families)')
ax[0].axhline(OM['dry'], color='gray', ls='--', label='exact 2.555042'); ax[1].axhline(OM['wet'], color='gray', ls='--', label='exact 1.804513')
for a, t in zip(ax, ('dry ring', 'wet ring (rho_f = 140)')): a.set_xscale('log'); a.set_xlabel('ring radial element size h (m)  (finer to the right)'); a.set_ylabel('omega (rad/s)'); a.set_title(t); a.invert_xaxis(); a.legend(fontsize=6)
fig.tight_layout(); fig.savefig('figures/fig_C1_frequency_vs_mesh.png'); plt.close(fig)

fig, ax = plt.subplots(figsize=(5.6, 3.6))
for f in fam:
    h, y = series(f, 'ratio'); ax.plot(h, y, mk[f] + '-', color=col[f], label=f)
ax.axhline(RATIO_REF, color='gray', ls='--', label='exact 0.706256'); ax.axhspan(RATIO_REF * 0.999, RATIO_REF * 1.001, color='gray', alpha=0.15, label='±0.1 % official ratio tolerance')
ax.set_xscale('log'); ax.invert_xaxis(); ax.set_xlabel('ring radial element size h (m)  (finer to the right)'); ax.set_ylabel('wet/dry frequency ratio (matched dry)'); ax.set_title('ratio (full range; see zoom in fig_C3 right)', fontsize=8); ax.legend(fontsize=6); fig.tight_layout(); fig.savefig('figures/fig_C2_ratio_vs_mesh.png'); plt.close(fig)

fig, ax = plt.subplots(1, 2, figsize=(9, 3.4))
for f in fam:
    h, y = series(f, 'ratio_err_pct'); ax[1].plot(h, y, mk[f] + '-', color=col[f], label=f)
h, y = series(list(fam)[0], 'err_dry_pct'); ax[0].plot(h, y, 'ko-', label='dry omega error'); h, y = series(list(fam)[0], 'err_wet_pct'); ax[0].plot(h, y, 'o-', color='C0', label='wet omega error (slip family)')
ax[0].axhline(0, color='gray', lw=0.8); ax[1].axhline(0, color='gray', lw=0.8); ax[1].axhspan(-0.1, 0.1, color='gray', alpha=0.15)
ax[0].set_ylabel('omega error vs exact (%)  [dry 2.555042, wet 1.804513]'); ax[0].set_title('frequency error', fontsize=9); ax[1].set_title('ratio error (shaded: +/-0.1 % official ratio tolerance)', fontsize=9); ax[1].set_ylabel('wet/dry ratio error vs 0.706256 (%)')
for a in ax: a.set_xscale('log'); a.invert_xaxis(); a.set_xlabel('ring radial element size h (m)  (finer to the right)'); a.legend(fontsize=6)
fig.tight_layout(); fig.savefig('figures/fig_C3_errors.png'); plt.close(fig)

SHORT = {'slip-wall + slip-interface (inviscid-limit formulation)': 'slip/slip\n(inviscid)', 'no-slip wall + no-slip interface (R1, D2)': 'no-slip/\nno-slip', 'slip wall + no-slip interface (D5)': 'D5: wall slip,\niface no-slip'}
fig, ax = plt.subplots(1, 3, figsize=(13, 3.8)); x = 0
labels, l0, l1, fe, lk = [], [], [], [], []
for r in rows:
    if r.get('ma_l0') is None: continue
    labels.append(SHORT[r['family']] + '\n' + r['level']); l0.append(r['ma_l0'] / 1.0); l1.append(r['ma_l1']); fe.append(r['n2_force_err_pct']); lk.append(r['leak_rms'] * 100 if r.get('leak_rms') else float('nan'))
xx = np.arange(len(labels)); ax[0].bar(xx - 0.2, l0, 0.4, label='layer 0'); ax[0].bar(xx + 0.2, l1, 0.4, label='layer 1')
pr = [r for r in rows if r.get('ma_l0') is not None]; ax[0].plot(xx - 0.2, [r['profile_l0'] for r in pr], 'k_', ms=22, mew=2); ax[0].plot(xx + 0.2, [r['profile_l1'] for r in pr], 'k_', ms=22, mew=2, label='exact profile at that layer'); ax[0].set_ylim(0.9, 1.1); ax[0].set_ylabel('n=2 pressure amplitude / (m_a w^2 D)'); ax[0].set_title('black ticks: exact n=2 profile at each layer position', fontsize=8); ax[0].legend(fontsize=6)
ax[1].bar(xx, fe, label='measured'); ax[1].plot(xx, [r['exact_force_err_pct_layer0'] for r in rows if r.get('ma_l0') is not None], 'k_', ms=22, mew=2, label='exact profile sampled at layer 0 (sampling effect only)'); ax[1].axhline(0, color='gray', lw=0.8); ax[1].set_ylim(-3, 7); ax[1].legend(fontsize=6, loc='upper right'); ax[1].set_ylabel('n=2 integrated radial force error (%)'); ax[1].set_title('n=2 interface force vs exact m_a*u_r', fontsize=8)
ax[2].bar(xx, np.nan_to_num(lk)); ax[2].set_ylabel('interface w_n rms (% of ring normal speed)'); ax[2].set_title('leakage (only where w_n is a natural BC; w=0 runs: n/a, enforced)', fontsize=7)
for i, v in enumerate(lk):
    if v != v: ax[2].text(i, 0.05, 'n/a', ha='center', fontsize=7)
for a in ax: a.set_xticks(xx); a.set_xticklabels(labels, fontsize=6)
fig.tight_layout(); fig.savefig('figures/fig_C4_pressure_traction.png'); plt.close(fig)

if sens:
    names = [s_['effect'].replace(', paired dry+wet', '') for s_ in sens]
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.4))
    ax[0].barh(names, [max(abs(s_['d_ratio_pct']), 1e-5) for s_ in sens]); ax[0].set_xscale('log'); ax[0].set_xlabel('|change in wet/dry ratio| (percentage points, log)'); ax[0].set_title('ratio sensitivity (M1 baseline C0)', fontsize=9)
    ax[1].barh(names, [max(abs(s_['d_omega_wet_pct']), 1e-5) for s_ in sens]); ax[1].set_xscale('log'); ax[1].set_yticklabels([]); ax[1].set_xlabel('|change in wet omega| (%, log)'); ax[1].set_title('wet frequency sensitivity', fontsize=9)
    fig.text(0.5, 0.005, 'Iterative (solver-tolerance) error NOT quantified: the tightened-tolerance run failed (dilatation noise floor). Inviscid-limit (slip) formulation.', ha='center', fontsize=7)
    fig.tight_layout(rect=(0, 0.04, 1, 1)); fig.savefig('figures/fig_C5_temporal_iterative_sensitivity.png'); plt.close(fig)

fig, ax = plt.subplots(figsize=(5, 3.4))
cr = [r for r in rows if r['family'].startswith('slip-wall') and r.get('wet_nodes')]
ax.loglog([r['wet_nodes'] for r in cr], [r['wet_wall_s'] for r in cr], 'o-', label='wet, 1 core'); cd = [r for r in rows if r['family'].startswith('slip-wall') and r.get('dry_nodes')]
ax.loglog([r['dry_nodes'] for r in cd], [r['dry_wall_s'] for r in cd], 's-', label='dry, 1 core'); [ax.annotate(r['level'], (r['wet_nodes'], r['wet_wall_s']), textcoords='offset points', xytext=(4, -10), fontsize=7) for r in cr]; ax.set_title('cost, slip family (M1, M2, M3)', fontsize=9); ax.set_xlabel('nodes'); ax.set_ylabel('wall time (s), 1 core'); ax.legend(fontsize=7); fig.tight_layout(); fig.savefig('figures/fig_C6_cost.png'); plt.close(fig)
print(open('results/phaseC_mesh_convergence.csv').read())
print(open('results/phaseC_time_iter_sensitivity.csv').read())
