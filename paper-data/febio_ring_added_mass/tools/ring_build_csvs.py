#!/usr/bin/env python3
"""Builds RING_SPATIAL_STUDY.csv, RING_TEMPORAL_STUDY.csv, RING_COMPRESSIBILITY_STUDY.csv, RING_COMPUTATIONAL_COST.csv from run reports.
Sources: runs/*/report.json (+pressure_modes.json, interface_leakage.json), the official reference JSON (exact values),
solids4foam README tables (results/solids4foam_ringAddedMass_README_tables.csv, high-order solid, printed values) and the raw PR #13 CSVs
(reference/pr13_ringAddedMass/*.csv, solids4foam STANDARD solid; ratio_error there is a fraction, converted to %). Nothing typed by hand. Missing runs are skipped."""
import json, os, csv, math, sys, glob, subprocess
sys.path.insert(0, 'tools')
os.chdir(os.environ.get('RING_ROOT', os.getcwd()))
REF = json.load(open('reference/s4f_ringAddedMass/verification/reference/ringAddedMass_verification_references.json'))['exact']
OD = REF['dry']['omega']; LV = REF['levels']; RHO = {'weak': 14, 'moderate': 140, 'strong': 1400}
MLAB = {1: 'M1', 2: 'M2', 4: 'M3'}
rj = lambda d, n='report.json': json.load(open(f'{d}/{n}')) if os.path.exists(f'{d}/{n}') else None
DRY = {1: 'runs/phaseB/B1_dry_hex8_f1', 2: 'runs/phaseB/D1_dry_f2', 4: 'runs/phaseC/C1_M3_dry_f4'}
SP = {  # steps per period -> dry run (M2)
    25: 'runs/ringNext/R_TIME_dry_M2_s25', 50: 'runs/ringNext/R_TIME_dry_M2_s50', 100: DRY[2]}
WET = {('weak', 1): 'runs/phaseD/D1_weak_M1', ('weak', 2): 'runs/phaseD/D4_weak_M2', ('weak', 4): 'runs/ringNext/R_MESH1_weak_M3',
       ('moderate', 1): 'runs/phaseC/C2_Emesh1e-3_wet_f1', ('moderate', 2): 'runs/phaseC/C2_Emesh1e-3_wet_f2', ('moderate', 4): 'runs/phaseD/D5_moderate_M3',
       ('strong', 1): 'runs/phaseD/D2_strong_M1', ('strong', 2): 'runs/phaseD/D3_strong_M2', ('strong', 4): 'runs/ringNext/R_MESH2_strong_M3'}
TW = {(l, s): (WET[(l, 2)] if s == 100 else f'runs/ringNext/R_TIME_{l}_M2_s{s}') for l in RHO for s in (25, 50, 100)}
S4F = [r for r in csv.DictReader(l for l in open('results/solids4foam_ringAddedMass_README_tables.csv') if not l.startswith('#'))]
P13R = list(csv.DictReader(open('reference/pr13_ringAddedMass/ringAddedMass_monolithic_runs.csv'))); P13M = list(csv.DictReader(open('reference/pr13_ringAddedMass/ringAddedMass_monolithic_mesh.csv')))
def s4(lev, f, spp, tab='high-order IQN-ILS'):
    return next((r for r in S4F if r['table'] == tab and r['level'] == lev and r['mesh_factor'] == str(f) and r['steps_per_period'] == str(spp)), None)
def p13r(kind, lev, coup, f, n):
    for r in P13R:
        if r['kind'] == kind and r['factor'] == str(f) and r['steps'] == str(n) and (kind == 'dry' or (r['level'] == lev and r['coupling'] == coup)): return r
def p13m(lev, coup, f):
    key = f'wet_{lev}_{"monolithic" if coup == "monolithic" else "iqnils"}_standard_m{f}_n100_p4' + ('_fin' if coup == 'monolithic' else '')
    return next((r for r in P13M if r['name'] == key), None)
f_ = lambda x, n=8: ('' if x is None or (isinstance(x, float) and math.isnan(x)) else (f'{x:.{n}g}' if isinstance(x, float) else x))
def write(path, rows):
    keys = []; [keys.append(k) for r in rows for k in r if k not in keys]
    with open(path, 'w', newline='') as fh:
        w = csv.writer(fh); w.writerow(keys); [w.writerow([f_(r.get(k)) for k in keys]) for r in rows]
    print(path, len(rows), 'rows')

# ---------------- spatial ----------------
rows = []
for f in (1, 2, 4):
    d = rj(DRY[f])
    rows.append(dict(regime='dry', mesh=MLAB[f], febio_run=os.path.basename(DRY[f]), omega=d['omega_fit'], omega_exact=OD, err_pct=100 * (d['omega_fit'] / OD - 1), amp_loss_pct=100 * d['amp_loss_per_period'],
                     s4f_highorder_iqnils_err_pct=(s4('dry', f, 100) or {}).get('error_pct'), s4f_standard_iqnils_err_pct=(lambda r: r and 100 * (float(r['omega']) / OD - 1))(p13m('dry', None, f) if False else next((r for r in P13M if r['name'] == f'dry_standard_m{f}_n100_p4'), None)),
                     wall_s=d['wall_s'], cores=d['cores'], cpu_h=d['cpu_s'] / 3600, rss_MB=d['max_rss_MB']))
for lev in ('weak', 'moderate', 'strong'):
    L = LV[lev]; prev = None
    for f in (1, 2, 4):
        w = rj(WET[(lev, f)]); d = rj(DRY[f])
        if not w: continue
        ratio = w['omega_fit'] / d['omega_fit']; ew = 100 * (w['omega_fit'] / L['omega'] - 1); ed = 100 * (d['omega_fit'] / OD - 1); er = 100 * (ratio / L['ratio'] - 1)
        pm = rj(WET[(lev, f)], 'pressure_modes.json') or {}; lk = rj(WET[(lev, f)], 'interface_leakage.json') or {}
        a = s4(lev, f, 100); mono = p13m(lev, 'monolithic', f); iq = p13m(lev, 'iqnils', f)
        r = dict(regime=lev, rho_f=RHO[lev], mesh=MLAB[f], febio_run=os.path.basename(WET[(lev, f)]), dry_run=os.path.basename(DRY[f]), K=w['K'], E_mesh=w['e_mesh'], steps_per_period=100, cores=w['cores'],
                 omega_dry=d['omega_fit'], omega_wet=w['omega_fit'], omega_wet_exact=L['omega'], err_dry_pct=ed, err_wet_pct=ew, ratio=ratio, ratio_exact=L['ratio'], ratio_err_pct=er,
                 abs_err_difference_wet_minus_dry_pct=ew - ed, cancellation_factor_abs_err_over_ratio_err=(abs(ew) / abs(er) if er else None),
                 equiv_additive_domega2_rad2_s2=2 * L['omega'] ** 2 * er / 100, mu_eff=w.get('mu_eff'),
                 n2_force_err_pct=100 * (w['n2_force_N'] / w['n2_force_analytic_N'] - 1), n2_p_amp_Pa=pm.get('n2', {}).get('amp_w'), p2_phase_deg=w['p2_phase_deg'],
                 leak_rms_pct=100 * lk.get('rms_leak_over_rms_ring_normal_velocity', float('nan')), net_flux_leak_pct=100 * lk.get('net_leaked_flux_peak_over_peak_abs_ring_flux', float('nan')),
                 amp_loss_per_period_pct=100 * w['amp_loss_per_period'], window_spread_pct=w['window_spread_pct'], meanJm1_max=w['volume_meanJm1_max'], iters_per_step=w['avg_newton_iters_per_step'],
                 wall_s=w['wall_s'], cpu_h=w['cpu_s'] / 3600, rss_MB=w['max_rss_MB'], nodes=w['n_nodes'],
                 s4f_highorder_iqnils_omega_err_pct=a and a['error_pct'], s4f_highorder_iqnils_ratio_err_pct=a and a['ratio_error_pct'],
                 s4f_standard_iqnils_ratio_err_pct=iq and 100 * float(iq['ratio_error']), s4f_standard_monolithic_ratio_err_pct=mono and 100 * float(mono['ratio_error']))
        rows.append(r)
# observed omega order over M1..M3 (consistency check only: standard hex8 solid; wet omega converges at second order only if three levels exist)
for lev in ('weak', 'moderate', 'strong'):
    ser = [r for r in rows if r['regime'] == lev]
    if len(ser) == 3:
        w = [r['omega_wet'] for r in ser]; d1, d2 = w[0] - w[1], w[1] - w[2]; o = math.log(d1 / d2, 2) if d1 > 0 and d2 > 0 else None
        for r in ser: r['omega_wet_observed_order_M1M2M3'] = o
        ser[-1]['omega_wet_richardson_limit'] = w[2] - d2 / (2 ** o - 1) if o else None
for lev in ('weak', 'moderate', 'strong'):
    ser = [r for r in rows if r['regime'] == lev]
    if len(ser) == 3:
        e = [r['ratio_err_pct'] for r in ser]; d1, d2 = e[1] - e[0], e[2] - e[1]
        if d1 * d2 > 0:
            o = math.log(d1 / d2, 2); lim = e[2] + d2 / (2 ** o - 1)
            for r in ser: r['ratio_err_observed_order_M1M2M3'] = o
            ser[-1]['ratio_err_richardson_floor_pct'] = lim; ser[-1]['equiv_additive_domega2_at_floor'] = 2 * LV[lev]['omega'] ** 2 * lim / 100
write('RING_SPATIAL_STUDY.csv', rows)

# ---------------- temporal ----------------
rows = []
for s in (25, 50, 100):
    d = rj(SP[s]); 
    if not d: continue
    T = 2 * math.pi / OD
    rows.append(dict(regime='dry', steps_per_period=s, dt_s=2 * math.pi / OD / s, omega=d['omega_fit'], err_pct=100 * (d['omega_fit'] / OD - 1), amp_loss_per_period_pct=100 * d['amp_loss_per_period'],
                     s4f_highorder_iqnils_err_pct=(s4('dry', 2, s) or {}).get('error_pct'), s4f_standard_iqnils_omega_err_pct=(lambda r: r and float(r['error']))(p13r('dry', None, None, 2, s)),
                     wall_s=d['wall_s'], cores=d['cores'], cpu_h=d['cpu_s'] / 3600))
for lev in ('weak', 'moderate', 'strong'):
    ser = []
    for s in (25, 50, 100):
        w = rj(TW[(lev, s)]); d = rj(SP[s])
        if not (w and d): continue
        L = LV[lev]; ratio = w['omega_fit'] / d['omega_fit']; a = s4(lev, 2, s)
        pi_ = p13r('wet', lev, 'iqnils', 2, s); pm_ = p13r('wet', lev, 'monolithic', 2, s)
        r = dict(regime=lev, rho_f=RHO[lev], steps_per_period=s, dt_s=w['dt'], febio_run=os.path.basename(TW[(lev, s)]), omega_dry=d['omega_fit'], omega_wet=w['omega_fit'], err_dry_pct=100 * (d['omega_fit'] / OD - 1),
                 err_wet_pct=100 * (w['omega_fit'] / L['omega'] - 1), ratio=ratio, ratio_err_pct=100 * (ratio / L['ratio'] - 1), amp_loss_per_period_pct=100 * w['amp_loss_per_period'],
                 n2_force_err_pct=100 * (w['n2_force_N'] / w['n2_force_analytic_N'] - 1), iters_per_step=w['avg_newton_iters_per_step'], wall_s=w['wall_s'], cores=w['cores'], cpu_h=w['cpu_s'] / 3600, rss_MB=w['max_rss_MB'],
                 s4f_highorder_iqnils_omega_err_pct=a and a['error_pct'], s4f_highorder_iqnils_ratio_err_pct=a and a['ratio_error_pct'], s4f_highorder_iqnils_loss_pct=a and a['loss_pct'],
                 s4f_standard_iqnils_ratio_err_pct=pi_ and 100 * float(pi_['ratio_error']), s4f_standard_monolithic_ratio_err_pct=pm_ and 100 * float(pm_['ratio_error']))
        ser.append(r); rows.append(r)
    if len(ser) == 3:
        w = [r['omega_wet'] for r in ser]; d1, d2 = w[1] - w[0], w[2] - w[1]; o = math.log(d1 / d2, 2)
        for r in ser: r['omega_wet_observed_temporal_order'] = o
        ser[-1]['omega_wet_richardson_dt0'] = w[2] + d2 / (2 ** o - 1)
        e = [r['ratio_err_pct'] for r in ser]; ser[-1]['ratio_err_pct_richardson_dt0'] = e[2] + (e[2] - e[1]) / (2 ** o - 1); ser[-1]['ratio_err_change_25_to_100_pp'] = e[2] - e[0]
dr = [r for r in rows if r['regime'] == 'dry']
if len(dr) == 3:
    w = [r['omega'] for r in dr]; d1, d2 = w[1] - w[0], w[2] - w[1]; o = math.log(d1 / d2, 2)
    for r in dr: r['omega_observed_temporal_order'] = o
    dr[-1]['omega_richardson_dt0'] = w[2] + d2 / (2 ** o - 1); dr[-1]['err_richardson_dt0_pct'] = 100 * (dr[-1]['omega_richardson_dt0'] / OD - 1)
write('RING_TEMPORAL_STUDY.csv', rows)

# ---------------- compressibility ----------------
import numpy as np
from scipy.special import jv, yv, jvp, yvp
def fk(k, a=1.05, b=1.5, n=2):
    C = -jvp(n, k * b) / yvp(n, k * b); return (jv(n, k * a) + C * yv(n, k * a)) / (k * (jvp(n, k * a) + C * yvp(n, k * a)))
f0 = fk(1e-4)
def theory_omega_shift_rel(lev, K1, K2):
    """predicted d omega/omega when the bulk modulus goes K1 -> K2 (linear acoustics, rigid wall at 1.5 m, thin-ring frac = 1 - ratio^2)"""
    L = LV[lev]; rho = RHO[lev]; w = L['omega']; frac = 1 - L['ratio'] ** 2
    d = lambda K: fk(w / math.sqrt(K / rho)) / f0 - 1
    return -0.5 * frac * (d(K2) - d(K1)), d(K1), d(K2)
rows = []
pairs = [('strong', 1, 'runs/phaseD/D2_strong_M1', 'runs/ringNext/R1_strong_M1_K1e9', 1e8, 1e9, 'E_mesh=1e-3, slip/slip, hex8 (this campaign)'),
         ('moderate', 1, 'runs/phaseC/C0_wet_f1_wallslip_ifslip', 'runs/phaseC/C1_K1e7_wet_f1', 1e8, 1e7, 'E_mesh=1 Pa, slip/slip (Phase C); wet/dry ratio shift')]
for lev, f, a_, b_, K1, K2, note in pairs:
    ra, rb = rj(a_), rj(b_)
    if not (ra and rb): continue
    dry = rj(DRY[f]); L = LV[lev]
    pred, d1, d2 = theory_omega_shift_rel(lev, K1, K2)
    for r, K in ((ra, K1), (rb, K2)):
        ratio = r['omega_fit'] / dry['omega_fit']; pm = rj(a_ if r is ra else b_, 'pressure_modes.json') or {}; lk = rj(a_ if r is ra else b_, 'interface_leakage.json') or {}
        rows.append(dict(regime=lev, mesh=f'M{f}', run=r['run'], K=K, sound_speed_m_s=math.sqrt(K / RHO[lev]), omega_wet=r['omega_fit'], omega_err_pct=100 * (r['omega_fit'] / L['omega'] - 1), ratio=ratio, ratio_err_pct=100 * (ratio / L['ratio'] - 1),
                         n2_p_amp_Pa=pm.get('n2', {}).get('amp_w'), n2_force_N=r['n2_force_N'], leak_rms_pct=100 * lk.get('rms_leak_over_rms_ring_normal_velocity', float('nan')), amp_loss_per_period_pct=100 * r['amp_loss_per_period'],
                         window_spread_pct=r['window_spread_pct'], iters_per_step=r['avg_newton_iters_per_step'], wall_s=r['wall_s'], cores=r['cores'], note=note))
    rows[-1]['measured_domega_over_omega_vs_first_row'] = rb['omega_fit'] / ra['omega_fit'] - 1
    rows[-1]['theory_domega_over_omega'] = pred
    rows[-1]['theory_added_mass_rel_change_at_K1'] = d1; rows[-1]['theory_added_mass_rel_change_at_K2'] = d2
    rows[-1]['measurement_floor_window_spread_rel'] = max(ra['window_spread_pct'], rb['window_spread_pct']) / 100
# theory-only table (actual sound speeds) for all regimes and K
for lev in ('weak', 'moderate', 'strong'):
    for K in (1e7, 1e8, 1e9):
        L = LV[lev]; c = math.sqrt(K / RHO[lev]); k = L['omega'] / c; d = fk(k) / f0 - 1
        rows.append(dict(regime=lev, run='THEORY (linear acoustics, no FEBio)', K=K, sound_speed_m_s=c, omega_wet=L['omega'], note=f'(omega b/c)^2={(k*1.5)**2:.3e}; added-mass relative change vs incompressible {d:+.3e}; d(ratio)/ratio = {-0.5*(1-L["ratio"]**2)*d:+.3e}'))
write('RING_COMPRESSIBILITY_STUDY.csv', rows)

# ---------------- cost ----------------
rows = []
for d in sorted(set(glob.glob('runs/phaseD/*') + glob.glob('runs/ringNext/*') + [DRY[1], DRY[2], DRY[4], 'runs/phaseC/C2_Emesh1e-3_wet_f1', 'runs/phaseC/C2_Emesh1e-3_wet_f2'])):
    r = rj(d)
    if not r: continue
    rows.append(dict(run=os.path.basename(d), group=d.split('/')[1], regime=('dry' if not r['wet'] else {14.0: 'weak', 140.0: 'moderate', 1400.0: 'strong'}[r['rho_f']]), mesh_nodes=r['n_nodes'], ring_elems=r['n_ring_el'], fluid_elems=r['n_fluid_el'],
                     steps=r['steps_completed'], steps_per_period=r['steps_per_period'], K=r['K'], cores=r['cores'], wall_s=r['wall_s'], cpu_s=r['cpu_s'], allocated_cpu_h=r['cores'] * r['wall_s'] / 3600, cpu_utilisation=r['cpu_s'] / (r['cores'] * r['wall_s']),
                     rss_MB=r['max_rss_MB'], newton_iters_per_step=r['avg_newton_iters_per_step'], residual_evals=r['residual_evals'], stiffness_reforms=r['stiffness_reforms'], linear_solves=r['linear_solves'], slurm_job=r['slurm_job'], febio_sha=r['febio_sha256']))
write('RING_COMPUTATIONAL_COST.csv', rows)
