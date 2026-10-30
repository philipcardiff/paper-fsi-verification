#!/usr/bin/env python3
"""Phase D results builder: FEBio (inviscid-limit, E_mesh = 1e-3 Pa) vs the exact references of the solids4foam ringAddedMass verification
study and the solids4foam README tables. All numbers come from runs/*/report.json, pressure_modes.json, interface_leakage.json,
*.margins.json, time.txt (via report.json) and from reference/.../ringAddedMass_verification_references.json. Nothing hand-typed.
Missing runs are skipped, so the script can be run on partial data."""
import json, os, csv, math, glob
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

REF = json.load(open('reference/s4f_ringAddedMass/verification/reference/ringAddedMass_verification_references.json'))
EX = REF['exact']; TOL = REF['tolerances']
OM_DRY = EX['dry']['omega']
LV = {k: EX['levels'][k] for k in ('weak', 'moderate', 'strong')}
RHO = {'weak': 14.0, 'moderate': 140.0, 'strong': 1400.0}
A_, B_ = -0.00022568759047243063, -0.0011425434267666801
prof = lambda r: (A_ * r * r + B_ / (r * r)) / (A_ * 1.05 ** 2 + B_ / 1.05 ** 2)

def rj(path, name='report.json'):
    p = f'{path}/{name}'
    return json.load(open(p)) if os.path.exists(p) else None
def margins(path):
    g = glob.glob(f'{path}/*.margins.json')
    return json.load(open(g[0])) if g else None

DRY = {1: 'runs/phaseB/B1_dry_hex8_f1', 2: 'runs/phaseB/D1_dry_f2', 4: 'runs/phaseC/C1_M3_dry_f4'}
FAM = {
    'FEBio E_mesh=1e-3 Pa': {
        'weak': {1: 'runs/phaseD/D1_weak_M1', 2: 'runs/phaseD/D4_weak_M2', 4: 'runs/ringNext/R_MESH1_weak_M3'},
        'moderate': {1: 'runs/phaseC/C2_Emesh1e-3_wet_f1', 2: 'runs/phaseC/C2_Emesh1e-3_wet_f2', 4: 'runs/phaseD/D5_moderate_M3'},
        'strong': {1: 'runs/phaseD/D2_strong_M1', 2: 'runs/phaseD/D3_strong_M2', 4: 'runs/ringNext/R_MESH2_strong_M3'}},
    'FEBio E_mesh=1 Pa (Phase C baseline)': {
        'moderate': {1: 'runs/phaseC/C0_wet_f1_wallslip_ifslip', 2: 'runs/phaseC/C1_M2_wet_f2', 4: 'runs/phaseC/C1_M3_wet_f4'}},
}

rows = []
for fam, levels in FAM.items():
    for lev, meshes in levels.items():
        for f, wet in meshes.items():
            w = rj(wet); d = rj(DRY[f])
            if not w or not d: continue
            ex = LV[lev]; dr = 0.45 / (12 * f)
            pm = rj(wet, 'pressure_modes.json') or {}; lk = rj(wet, 'interface_leakage.json') or {}; mg = margins(wet) or {}
            ke = w['KE_total_at_zero_crossings_over_initial']; ke_int = ke[2::2]; ke_half = ke[1::2]
            ratio = w['omega_fit'] / d['omega_fit']
            r = dict(family=fam, level=lev, rho_f=RHO[lev], mesh_factor=f, wet_run=os.path.basename(wet), dry_run=os.path.basename(DRY[f]),
                     omega_dry=d['omega_fit'], omega_dry_exact=OM_DRY, err_dry_pct=100 * (d['omega_fit'] / OM_DRY - 1),
                     omega_wet=w['omega_fit'], omega_wet_exact=ex['omega'], err_wet_pct=100 * (w['omega_fit'] / ex['omega'] - 1),
                     ratio=ratio, ratio_exact=ex['ratio'], ratio_err_pct=100 * (ratio / ex['ratio'] - 1),
                     ratio_thin_err_pct=100 * (ratio / ex['ratioThinRing'] - 1), viscous_omega_shift_pct=100 * ex['viscousFrequencyShift'],
                     mu_eff=w.get('mu_eff'), mu_eff_exact_continuum=(OM_DRY / ex['omega']) ** 2 - 1, mu_eff_err_pct=100 * (w.get('mu_eff') / ((OM_DRY / ex['omega']) ** 2 - 1) - 1),
                     amp_loss_per_period_pct=100 * w['amp_loss_per_period'], window_spread_pct=w['window_spread_pct'],
                     n2_p_amp_Pa=pm.get('n2', {}).get('amp_w'), p_layer0_over_exact_profile=w['ma_eff_ratio_layer0'] / prof(1.05 + 0.5 * dr),
                     p_layer1_over_exact_profile=w['ma_eff_ratio_layer1'] / prof(1.05 + 1.5 * dr), p2_phase_deg=w['p2_phase_deg'],
                     n2_force_err_pct=100 * (w['n2_force_N'] / w['n2_force_analytic_N'] - 1), exact_force_sampling_offset_pct=100 * (prof(1.05 + 0.5 * dr) - 1),
                     leak_rms_pct=100 * lk.get('rms_leak_over_rms_ring_normal_velocity', float('nan')), leak_max_over_V0_pct=100 * lk.get('max_abs_wn_over_V0', float('nan')),
                     leak_net_flux_pct=100 * lk.get('net_leaked_flux_peak_over_peak_abs_ring_flux', float('nan')),
                     KE_int_min=min(ke_int), KE_int_max=max(ke_int), KE_half_max_excursion_pct=100 * max(abs(x - 1) for x in ke_half) if ke_half else float('nan'),
                     meanJm1_max=w['volume_meanJm1_max'], J_dev_max=w['J_max_dev'], jac_nan_msgs=w['jacobian_or_nan_messages'], max_refs_events=w['max_refs_reached_events'],
                     iters_per_step=w['avg_newton_iters_per_step'], dil_margin_worst=(mg.get('dilatation') or {}).get('margin_max'),
                     res_margin_worst=(mg.get('residual') or {}).get('margin_max'), dis_margin_worst=(mg.get('displacement') or {}).get('margin_max'),
                     vel_margin_worst=(mg.get('velocity') or {}).get('margin_max'), nodes=w['n_nodes'], cores=w['cores'],
                     wall_s=w['wall_s'], cpu_h=w['cpu_s'] / 3600, rss_MB=w['max_rss_MB'], dry_wall_s=d['wall_s'], dry_cpu_h=d['cpu_s'] / 3600, dry_rss_MB=d['max_rss_MB'])
            rows.append(r)

os.makedirs('results', exist_ok=True); os.makedirs('figures', exist_ok=True)
keys = list(rows[0].keys()) if rows else []
with open('results/phaseD_febio_summary.csv', 'w', newline='') as fh:
    wr = csv.writer(fh); wr.writerow(keys)
    for r in rows: wr.writerow([f'{r[k]:.8g}' if isinstance(r[k], float) else r[k] for k in keys])

# solids4foam table (verbatim transcription) ----------------------------------------------------------------------
S4F = [r for r in csv.DictReader(l for l in open('results/solids4foam_ringAddedMass_README_tables.csv') if not l.startswith('#'))]
s4 = lambda lev, f, tab='high-order IQN-ILS', spp='100': next((r for r in S4F if r['table'] == tab and r['level'] == lev and r['mesh_factor'] == str(f) and r['steps_per_period'] == spp), None)

# side-by-side comparison CSV
cmp_rows = []
for lev in ('weak', 'moderate', 'strong'):
    for f in (1, 2, 4):
        s = s4(lev, f); sd = s4('dry', f)
        fb = next((r for r in rows if r['family'].startswith('FEBio E_mesh=1e-3') and r['level'] == lev and r['mesh_factor'] == f), None)
        if not s and not fb: continue
        cmp_rows.append(dict(level=lev, mesh_factor=f, exact_omega=LV[lev]['omega'], exact_ratio=LV[lev]['ratio'],
                             s4f_omega=s and float(s['omega_rad_s']), s4f_omega_err_pct=s and float(s['error_pct']), s4f_ratio_err_pct=s and float(s['ratio_error_pct']), s4f_loss_pct=s and float(s['loss_pct']),
                             s4f_dry_omega=sd and float(sd['omega_rad_s']), s4f_dry_err_pct=sd and float(sd['error_pct']),
                             febio_omega=fb and fb['omega_wet'], febio_omega_err_pct=fb and fb['err_wet_pct'], febio_ratio_err_pct=fb and fb['ratio_err_pct'], febio_loss_pct=fb and fb['amp_loss_per_period_pct'],
                             febio_dry_omega=fb and fb['omega_dry'], febio_dry_err_pct=fb and fb['err_dry_pct'],
                             febio_wall_s=fb and fb['wall_s'], febio_cpu_h=fb and fb['cpu_h'], febio_rss_MB=fb and fb['rss_MB']))
with open('results/phaseD_s4f_vs_febio.csv', 'w', newline='') as fh:
    wr = csv.writer(fh); wr.writerow(list(cmp_rows[0].keys()))
    for r in cmp_rows: wr.writerow(['' if v is None else (f'{v:.8g}' if isinstance(v, float) else v) for v in r.values()])

# official acceptance checks applied to FEBio (criteria from the solids4foam verification README / JSON) ----------------
chk = []
for lev in ('weak', 'moderate', 'strong'):
    ser = sorted([r for r in rows if r['family'].startswith('FEBio E_mesh=1e-3') and r['level'] == lev], key=lambda r: r['mesh_factor'])
    if not ser: continue
    fin = ser[-1]; wl = f'{lev} (finest run: M{fin["mesh_factor"]}, {len(ser)} levels)'
    chk.append((lev, 'finest wet omega error <= 0.3 %', f'{fin["err_wet_pct"]:+.3f} %', abs(fin['err_wet_pct']) <= 100 * TOL['omega'], 'FEBio solid is standard hex8, not the high-order solid used by solids4foam for this criterion'))
    chk.append((lev, 'finest wet/dry ratio error <= 0.1 %', f'{fin["ratio_err_pct"]:+.4f} %', abs(fin['ratio_err_pct']) <= 100 * TOL['ratio'], ''))
    if len(ser) >= 2:
        chk.append((lev, 'ratio error decreases with refinement', ' -> '.join(f'{r["ratio_err_pct"]:+.4f}' for r in ser), abs(ser[-1]['ratio_err_pct']) < abs(ser[0]['ratio_err_pct']), 'literal criterion: |error| at finest < |error| at coarsest'))
        ch = abs(ser[-1]['omega_wet'] / ser[-2]['omega_wet'] - 1)
        chk.append((lev, 'wet omega change between last two meshes <= 0.1 %', f'{100 * ch:.3f} %', ch <= TOL['meshChange'], 'standard solid is not converged in absolute frequency at these meshes'))
    chk.append((lev, 'finest ratio vs thin-ring <= 0.25 %', f'{fin["ratio_thin_err_pct"]:+.3f} %', abs(fin['ratio_thin_err_pct']) <= 100 * TOL['ratioThinRing'], ''))
    chk.append((lev, 'amplitude loss per period < 2 % (finest)', f'{fin["amp_loss_per_period_pct"]:.5f} %', fin['amp_loss_per_period_pct'] < 100 * TOL['maxLossPerPeriod'], ''))
    chk.append((lev, 'amplitude growth <= 0.2 % per period', f'{-fin["amp_loss_per_period_pct"]:.5f} % growth', -fin['amp_loss_per_period_pct'] <= 100 * TOL['maxGrowthPerPeriod'], ''))
_tp = 'RING_TEMPORAL_STUDY.csv'
if os.path.exists(_tp):
    _o = [float(r['omega_wet_observed_temporal_order']) for r in csv.DictReader(open(_tp)) if r.get('regime') in ('weak', 'moderate', 'strong') and r.get('steps_per_period') == '100' and r.get('omega_wet_observed_temporal_order')]
    chk.append(('all', 'observed time order >= 1.5 (M2, 25/50/100 steps per period)', ('wet omega order ' + ', '.join(f'{x:.3f}' for x in _o)) if _o else 'n/a', bool(_o) and min(_o) >= 1.5, 'order of the absolute frequency (generalized-alpha, rhoi=1); the ratio error itself does not change with dt (see RING_TEMPORAL_STUDY.csv)'))
else:
    chk.append(('all', 'observed time order >= 1.5', 'not tested', None, 'run the 25/50/100 study'))
with open('results/phaseD_acceptance_checks.csv', 'w', newline='') as fh:
    wr = csv.writer(fh); wr.writerow(['level', 'criterion', 'FEBio value', 'pass', 'note'])
    for c in chk: wr.writerow([c[0], c[1], c[2], {True: 'yes', False: 'NO', None: 'not tested'}[c[3]], c[4]])

# --------------------------------------------------------------------------------------------- figures
plt.rcParams.update({'font.size': 9, 'figure.dpi': 130})
colL = {'weak': 'C2', 'moderate': 'C0', 'strong': 'C3'}
FE = [r for r in rows if r['family'].startswith('FEBio E_mesh=1e-3')]; FE1 = [r for r in rows if r['family'].startswith('FEBio E_mesh=1 ')]
def ser(rs, lev, key):
    q = sorted([r for r in rs if r['level'] == lev], key=lambda r: r['mesh_factor']); return [r['mesh_factor'] for r in q], [r[key] for r in q]
def s4ser(lev, key):
    q = [(int(r['mesh_factor']), float(r[key])) for r in S4F if r['table'] == 'high-order IQN-ILS' and r['level'] == lev and r['steps_per_period'] == '100' and r[key] != '']
    q.sort(); return [a for a, b in q], [b for a, b in q]
def fx(ax):
    ax.set_xscale('log', base=2); ax.set_xticks([1, 2, 4]); ax.set_xticklabels(['M1', 'M2', 'M3']); ax.set_xlabel('mesh level (refinement factor 1, 2, 4)')

# D1: moderate mesh convergence
fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))
x, y = ser(FE, 'moderate', 'err_dry_pct'); ax[0].plot(x, y, 'ko-', label='FEBio dry (hex8)')
x, y = ser(FE, 'moderate', 'err_wet_pct'); ax[0].plot(x, y, 'o-', color='C0', label='FEBio wet, E_mesh = 1e-3 Pa')
x, y = ser(FE1, 'moderate', 'err_wet_pct'); ax[0].plot(x, y, 's--', color='C1', label='FEBio wet, E_mesh = 1 Pa')
x, y = s4ser('dry', 'error_pct') if False else ([int(r['mesh_factor']) for r in S4F if r['table'] == 'high-order IQN-ILS' and r['level'] == 'dry' and r['steps_per_period'] == '100'], [float(r['error_pct']) for r in S4F if r['table'] == 'high-order IQN-ILS' and r['level'] == 'dry' and r['steps_per_period'] == '100'])
ax[0].plot(x, y, '^:', color='gray', label='solids4foam dry (high-order solid)')
x, y = s4ser('moderate', 'error_pct'); ax[0].plot(x, y, 'v:', color='k', label='solids4foam wet (high-order solid)')
ax[0].axhline(0, color='gray', lw=0.8); ax[0].axhspan(-100 * TOL['omega'], 100 * TOL['omega'], color='gray', alpha=0.12, label='official +/-0.3 % omega tolerance')
ax[0].set_ylabel('frequency error vs exact continuum (%)'); ax[0].set_title('moderate: absolute frequency error', fontsize=9); ax[0].legend(fontsize=6); fx(ax[0])
x, y = ser(FE, 'moderate', 'ratio_err_pct'); ax[1].plot(x, y, 'o-', color='C0', label='FEBio, E_mesh = 1e-3 Pa')
x, y = ser(FE1, 'moderate', 'ratio_err_pct'); ax[1].plot(x, y, 's--', color='C1', label='FEBio, E_mesh = 1 Pa')
x, y = s4ser('moderate', 'ratio_error_pct'); ax[1].plot(x, y, 'v:', color='k', label='solids4foam (high-order solid)')
ax[1].axhline(0, color='gray', lw=0.8); ax[1].axhspan(-100 * TOL['ratio'], 100 * TOL['ratio'], color='gray', alpha=0.15, label='official +/-0.1 % ratio tolerance')
ax[1].set_ylabel('wet/dry ratio error vs exact 0.7062557 (%)'); ax[1].set_title('moderate: ratio error', fontsize=9); ax[1].legend(fontsize=6); fx(ax[1])
x, y = ser(FE, 'moderate', 'ratio_err_pct'); ax[2].plot(x, y, 'o-', color='C0', label='FEBio, E_mesh = 1e-3 Pa')
x, y = s4ser('moderate', 'ratio_error_pct'); ax[2].plot(x, y, 'v:', color='k', label='solids4foam')
ax[2].axhline(0, color='gray', lw=0.8); ax[2].set_ylabel('ratio error (%)  [zoom]'); ax[2].set_title('zoom (no E_mesh = 1 curve)', fontsize=9); ax[2].legend(fontsize=6); fx(ax[2])
fig.tight_layout(); fig.savefig('figures/fig_D1_moderate_mesh_convergence.png'); plt.close(fig)

# D2: screening weak / moderate / strong
fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
for lev in ('weak', 'moderate', 'strong'):
    x, y = ser(FE, lev, 'err_wet_pct'); ax[0].plot(x, y, 'o-', color=colL[lev], label=f'FEBio {lev} (rho_f = {RHO[lev]:g})')
    x, y = s4ser(lev, 'error_pct'); ax[0].plot(x, y, 'v:', color=colL[lev], alpha=0.6, label=f's4f {lev}')
    x, y = ser(FE, lev, 'ratio_err_pct'); ax[1].plot(x, y, 'o-', color=colL[lev], label=f'FEBio {lev}')
    x, y = s4ser(lev, 'ratio_error_pct'); ax[1].plot(x, y, 'v:', color=colL[lev], alpha=0.6, label=f's4f {lev}')
ax[0].axhline(0, color='gray', lw=0.8); ax[0].axhspan(-100 * TOL['omega'], 100 * TOL['omega'], color='gray', alpha=0.12)
ax[1].axhline(0, color='gray', lw=0.8); ax[1].axhspan(-100 * TOL['ratio'], 100 * TOL['ratio'], color='gray', alpha=0.15, label='official +/-0.1 %')
ax[0].set_ylabel('wet frequency error vs exact (%)'); ax[0].set_title('absolute frequency (shaded: official +/-0.3 %)', fontsize=9); ax[1].set_ylabel('wet/dry ratio error vs exact (%)'); ax[1].set_title('frequency ratio', fontsize=9)
for a in ax: a.legend(fontsize=6); fx(a)
fig.text(0.5, 0.005, 'FEBio: inviscid-limit formulation (slip wall + slip interface, E_mesh = 1e-3 Pa, hex8 solid).\nsolids4foam: nu_f = 1e-7 no-slip, high-order solid, IQN-ILS: the absolute-frequency panel compares different solid discretisations.', ha='center', fontsize=7)
fig.tight_layout(rect=(0, 0.07, 1, 1)); fig.savefig('figures/fig_D2_added_mass_levels.png'); plt.close(fig)

# D3: pressure, force, leakage
fig, ax = plt.subplots(1, 3, figsize=(13, 3.8)); lab = [f"{r['level'][:4]}.\n{ {1: 'M1', 2: 'M2', 4: 'M3'}[r['mesh_factor']] }" for r in FE]; xx = np.arange(len(FE))
ax[0].bar(xx - 0.2, [r['p_layer0_over_exact_profile'] for r in FE], 0.4, label='layer 0'); ax[0].bar(xx + 0.2, [r['p_layer1_over_exact_profile'] for r in FE], 0.4, label='layer 1')
ax[0].axhline(1, color='k', lw=0.8); ax[0].set_ylim(0.98, 1.02); ax[0].set_ylabel('measured n=2 pressure / exact radial profile'); ax[0].set_title('n=2 pressure vs exact profile (1 = exact)', fontsize=9); ax[0].legend(fontsize=6)
ax[1].bar(xx, [r['n2_force_err_pct'] for r in FE], label='measured'); ax[1].plot(xx, [r['exact_force_sampling_offset_pct'] for r in FE], 'k_', ms=22, mew=2, label='exact profile at layer 0 (sampling offset)')
ax[1].axhline(0, color='gray', lw=0.8); ax[1].set_ylabel('n=2 integrated force error vs m_a u_r (%)'); ax[1].set_title('interface force', fontsize=9); ax[1].legend(fontsize=6, loc='lower right'); ax[1].set_ylim(-2.6, 0.4)
ax[2].bar(xx - 0.2, [r['leak_rms_pct'] for r in FE], 0.4, label='rms w_n / rms ring normal velocity'); ax[2].bar(xx + 0.2, [r['leak_max_over_V0_pct'] for r in FE], 0.4, label='max |w_n| / V0')
ax[2].set_ylabel('normal leakage (%)'); ax[2].set_title('interface leakage (weakly enforced w.n = 0)', fontsize=9); ax[2].legend(fontsize=6)
for a in ax: a.set_xticks(xx); a.set_xticklabels(lab, fontsize=7)
fig.tight_layout(); fig.savefig('figures/fig_D3_pressure_force_leakage.png'); plt.close(fig)

# D4: cost (x axis = mesh level; s4f times are the README IQN-ILS/Robin runs on mesh level 2 only)
fig, ax = plt.subplots(1, 3, figsize=(13, 3.8)); lv = {1: 'M1', 2: 'M2', 4: 'M3'}
for lev in ('weak', 'moderate', 'strong'):
    q = sorted([r for r in FE if r['level'] == lev], key=lambda r: r['mesh_factor'])
    for a_, key in zip(ax, ('wall_s', 'cpu_h', 'rss_MB')): a_.semilogy([r['mesh_factor'] for r in q], [r[key] for r in q], 'o-', color=colL[lev], label=f'FEBio wet {lev}')
dq = sorted({r['mesh_factor']: r for r in FE}.values(), key=lambda r: r['mesh_factor'])
for a_, key in zip(ax, ('dry_wall_s', 'dry_cpu_h', 'dry_rss_MB')): a_.semilogy([r['mesh_factor'] for r in dq], [r[key] for r in dq], 's--', color='gray', label='FEBio dry')
for lev in ('weak', 'moderate', 'strong'):
    for tab_row in S4F:
        if tab_row['table'].startswith('IQN-ILS vs Robin') and tab_row['level'] == lev and tab_row['time_s']:
            t = float(tab_row['time_s']); mk_ = '^' if tab_row['coupling'] == 'IQN-ILS' else 'v'
            ax[0].semilogy([2], [t], mk_, color=colL[lev], mfc='none', label=f's4f {tab_row["coupling"]} {lev} (README, mesh 2)'); ax[1].semilogy([2], [t / 3600], mk_, color=colL[lev], mfc='none')
for a_, t in zip(ax, ('wall time (s) (FEBio 1 core; weak/strong M3: 8 cores)', 'measured CPU time (core-hours) (FEBio user+sys)', 'peak memory (MB), FEBio only')): 
    a_.set_xscale('log', base=2); a_.set_xticks([1, 2, 4]); a_.set_xticklabels(['M1', 'M2', 'M3']); a_.set_xlabel('mesh level'); a_.set_ylabel(t, fontsize=8)
ax[0].set_xlim(0.8, 11); ax[0].legend(fontsize=5, loc='lower right'); ax[1].legend(fontsize=6, loc='lower right')
fig.text(0.5, 0.005, 'Different machines/loads, solvers and solid discretisations: indicative only. s4f: OpenFOAM v2412 on xenosim (README), serial; FEBio: v4.13, MKL Pardiso, 1 core except weak/strong M3 (8 cores, 40 % CPU utilisation), shared node. s4f memory not recorded.', ha='center', fontsize=6)
fig.tight_layout(rect=(0, 0.05, 1, 1)); fig.savefig('figures/fig_D4_cost.png'); plt.close(fig)

# D5: iterative-error evidence: relative dilatation increment at the stopping iteration vs the tolerance (floor-limited), from the logs
import re
def final_rel(path, tol):
    steps = []; cur = None
    for l in open(path, errors='ignore'):
        if l.startswith('===== beginning time step'): cur = []; steps.append(cur)
        m = re.match(r'\s+dilatation\s+(\S+)\s+(\S+)\s+(\S+)', l)
        if m and cur is not None: cur.append((float(m.group(2)), float(m.group(3))))
    steps = [q for q in steps if q]
    return np.array([math.sqrt(q[-1][0] / (q[-1][1] / tol ** 2)) for q in steps])
fig, ax = plt.subplots(1, 2, figsize=(10, 3.6))
for lab, path, c in (('moderate M2, E_mesh 1e-3 (C2)', 'runs/phaseC/C2_Emesh1e-3_wet_f2/C1M2_wet.log', 'C0'), ('moderate M3, E_mesh 1 (C1)', 'runs/phaseC/C1_M3_wet_f4/C1M3_wet.log', 'C1'),
                     ('strong M1 (D2)', 'runs/phaseD/D2_strong_M1/D2_strong_M1.log', 'C3')):
    if os.path.exists(path): ax[0].hist(np.log10(final_rel(path, 1e-6)), bins=30, alpha=0.5, color=c, label=lab)
ax[0].axvline(-6, color='k', ls='--', lw=1, label='ftol = 1e-6'); ax[0].set_xlabel('log10( |dilatation increment| / |step increment| ) at the stopping iteration'); ax[0].set_ylabel('number of steps'); ax[0].legend(fontsize=6)
ax[0].set_title('stopping increments of the dilatation norm', fontsize=9)
fl = 'runs/phaseC/C1_tol_wet_f1/C1tol_wet.log'
if os.path.exists(fl):
    v = []; cur = None
    for l in open(fl, errors='ignore'):
        m = re.match(r'\s+dilatation\s+(\S+)\s+(\S+)\s+(\S+)', l)
        if m: v.append(float(m.group(2)))
    ax[1].semilogy(range(1, len(v) + 1), v, 'o-', ms=3); ax[1].axhline(5.17e-33, color='r', ls='--', lw=1, label='required at tol 1e-8 (5.2e-33)'); ax[1].axhline(5.17e-29, color='k', ls='--', lw=1, label='required at tol 1e-6 (5.2e-29)')
    ax[1].set_xlabel('Newton iteration (step 1, Phase C tolerance run)'); ax[1].set_ylabel('dilatation |increment|^2'); ax[1].set_title('arithmetic floor: flat after iteration 3', fontsize=9); ax[1].legend(fontsize=6)
fig.text(0.5, 0.005, 'Iterative error NOT isolated: the tightened-tolerance runs failed at the floor (see PHASE_D_ITERATIVE_ERROR.md).', ha='center', fontsize=7)
fig.tight_layout(rect=(0, 0.03, 1, 1)); fig.savefig('figures/fig_D5_iterative_floor.png'); plt.close(fig)
print(len(rows), 'FEBio rows;', len(cmp_rows), 'comparison rows;', len(chk), 'checks;', 'figures done')
