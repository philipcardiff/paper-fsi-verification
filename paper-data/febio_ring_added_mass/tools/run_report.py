#!/usr/bin/env python3
"""Uniform per-run report for ring runs. usage: run_report.py <run_dir> [matching_dry_run_dir]
Writes <run_dir>/report.json and prints a summary. Post-processing only."""
import sys, os, re, json, math
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import analyze_ring as AR
import xml.etree.ElementTree as ET

run = sys.argv[1].rstrip('/')
dry_dir = sys.argv[2].rstrip('/') if len(sys.argv) > 2 else None
meta = json.load(open(run + '/model_meta.json'))
tag = meta['tag']; wet = meta['rho_f'] is not None
nrs, nt, nrf = meta['dims_nrs_nt_nrf']
omega_ref = meta['omega_ref']
rep = dict(run=os.path.basename(run), tag=tag, wet=wet, mesh=f"ring {nrs}x{nt}" + (f", fluid {nrf}x{nt}" if wet else ""),
           K=meta['k'] if wet else None, rho_f=meta['rho_f'], rhoi=meta['rhoi'], ic_field=meta.get('ic_field', 'potential(v1)') if wet else None,
           wall=meta.get('wall', 'noslip') if wet else None, e_mesh=meta['e_mesh'] if wet else None, steps_per_period=meta['steps_per_period'],
           dt=meta['dt'], n_nodes=meta['n_nodes'], n_ring_el=meta['n_ring_el'], n_fluid_el=meta['n_fluid_el'], feb_sha256=meta['feb_sha256'])
# run metadata
s = json.load(open(run + '/run_meta_start.json')); e = json.load(open(run + '/run_meta_end.json'))
rep.update(slurm_job=s['slurm_job_id'], partition=s['partition'], cores=s['cpus_per_task'], febio_sha256=s['febio_sha256'][:12], exit_code=e['exit_code'])
t = open(run + '/time.txt').read()
w = re.search(r'Elapsed \(wall clock\) time.*: ([\d:.]+)', t).group(1); parts = [float(x) for x in w.split(':')]
rep['wall_s'] = sum(x * 60 ** i for i, x in enumerate(reversed(parts)))
rep['cpu_s'] = float(re.search(r'User time \(seconds\): ([\d.]+)', t).group(1)) + float(re.search(r'System time \(seconds\): ([\d.]+)', t).group(1))
rep['max_rss_MB'] = int(re.search(r'Maximum resident set size \(kbytes\): (\d+)', t).group(1)) / 1024
log = open(run + '/' + tag + '.log', errors='ignore').read()
rep['steps_completed'] = len(re.findall(r'converged at time', log)); rep['steps_expected'] = meta['n_steps']
rep['normal_termination'] = 'N O R M A L   T E R M I N A T I O N' in log
rep['jacobian_or_nan_messages'] = len(re.findall(r'(?i)negative jacobian|NaN|inverted', log.replace('NaN-', '')))
m = re.search(r'Total number of right hand evaluations \.* : (\d+)', log); rep['residual_evals'] = int(m.group(1)) if m else None
m = re.search(r'Total number of stiffness reformations \.* : (\d+)', log); rep['stiffness_reforms'] = int(m.group(1)) if m else None
m = re.search(r'Total calls to linear solver \.* : (\d+)', log); rep['linear_solves'] = int(m.group(1)) if m else None
rep['avg_newton_iters_per_step'] = (rep['residual_evals'] / rep['steps_completed']) if rep['residual_evals'] and rep['steps_completed'] else None
rep['max_refs_reached_events'] = log.count('Max nr of iterations reached')
# frequency
an = AR.analyse(run + '/probes.txt', omega_ref)
W = an['windows']['0-4T']
rep.update(omega_fit=W['omega_fit'], omega_zc=W['omega_zero_crossing'], omega_ref=omega_ref, err_fit_pct=100 * W['err_fit'], err_zc_pct=100 * W['err_zc'],
           window_spread_pct=100 * (max(v['omega_fit'] for v in an['windows'].values()) - min(v['omega_fit'] for v in an['windows'].values())) / omega_ref,
           amp_loss_per_period=W['amp_loss_per_period'], amplitude0=W['amplitude0'])
T = 2 * math.pi / W['omega_fit']
# velocity consistency
pb = AR.read_blocks(run + '/probes.txt'); tt = np.array([b['t'] for b in pb]); r0 = np.array([b['rows'][0] for b in pb]); u, v = r0[:, 3], r0[:, 5]
dtm = tt[1] - tt[0]; d = v[1:-1] - (u[2:] - u[:-2]) / (2 * dtm)
rep['vel_vs_dudt_rms_pct'] = 100 * d.std() / np.abs(v).max()
if wet:
    # matching dry
    if dry_dir:
        dr = json.load(open(dry_dir + '/report.json')); rep['dry_ref_run'] = os.path.basename(dry_dir)
        rep['wet_dry_ratio'] = W['omega_fit'] / dr['omega_fit']; rep['ratio_err_pct'] = 100 * (rep['wet_dry_ratio'] / (omega_ref / 2.55504257) - 1)
        rep['mu_eff'] = (dr['omega_fit'] / W['omega_fit']) ** 2 - 1
    # pressure
    NR = nrf
    def blocks(p): return AR.read_blocks(p)
    fb = blocks(run + '/fluid_elems.txt'); tf = np.array([b['t'] for b in fb])
    P = np.array([[r[1] for r in b['rows']] for b in fb]); J = np.array([[r[2] for r in b['rows']] for b in fb])
    D = np.array([0.5 * (b['rows'][0][3] - b['rows'][1][4]) for b in pb])
    dth = (math.pi / 2) / nt; c2 = np.cos(2 * (np.arange(nt) + 0.5) * dth)
    P2 = lambda i: (P[:, [j * NR + i for j in range(nt)]] @ c2) * dth / (math.pi / 4)
    def harm(x, wv):
        M = np.column_stack([np.sin(wv * tf), np.cos(wv * tf), np.ones_like(tf)]); c, *_ = np.linalg.lstsq(M, x, rcond=None); return c, (x - M @ c).std()
    wv = W['omega_fit']; cD, _ = harm(D, wv); aD = math.hypot(cD[0], cD[1])
    rho, Ro, bb, h = 140.0 if meta['rho_f'] is None else meta['rho_f'], 1.05, 1.5, 0.1
    ma = rho * Ro / 2 * (bb ** 4 + Ro ** 4) / (bb ** 4 - Ro ** 4)
    c0, n0 = harm(P2(0), wv); c1, n1 = harm(P2(1), wv)
    rep['p2_layer0_amp_Pa'] = math.hypot(c0[0], c0[1]); rep['p2_layer1_amp_Pa'] = math.hypot(c1[0], c1[1])
    rep['p2_pred_Pa'] = ma * wv ** 2 * aD
    rep['ma_eff_ratio_layer0'] = -((c0[0] * cD[0] + c0[1] * cD[1]) / aD ** 2) / wv ** 2 / ma
    rep['ma_eff_ratio_layer1'] = -((c1[0] * cD[0] + c1[1] * cD[1]) / aD ** 2) / wv ** 2 / ma
    rep['p2_phase_deg'] = math.degrees(math.atan2(c0[1], c0[0]) - math.atan2(cD[1], cD[0])) % 360
    rep['n2_force_N'] = math.hypot(c0[0], c0[1]) * Ro * h * math.pi / 4; rep['n2_force_analytic_N'] = ma * wv ** 2 * aD * Ro * h * math.pi / 4
    pp = P[:, 0]; c, r = harm(pp, wv)
    rep['p_probe_amp_Pa'] = math.hypot(c[0], c[1]); rep['p_probe_noise_rms_Pa'] = r; rep['p_max_any_Pa'] = float(np.abs(P).max())
    rep['p_max_first5_Pa'] = float(np.abs(P[:6]).max()); rep['p_max_after5_Pa'] = float(np.abs(P[6:]).max())
    rep['p_noise_over_analytic'] = r / (ma * wv * 1.5e-3)
    rep['volume_meanJm1_max'] = float(np.abs((J - 1).mean(axis=1)).max()); rep['J_max_dev'] = float(np.abs(J - 1).max())
    # kinetic energy at displacement zero crossings (lumped mass)
    root = ET.parse(run + '/' + tag + '.feb').getroot()
    N = {int(n.get('id')): np.array([float(x) for x in n.text.split(',')]) for n in root.find('Mesh/Nodes')}
    dom = {el.get('name'): [[int(x) for x in q.text.split(',')] for q in el] for el in root.findall('Mesh/Elements')}
    sx = [(-1,-1,-1),(1,-1,-1),(1,1,-1),(-1,1,-1),(-1,-1,1),(1,-1,1),(1,1,1),(-1,1,1)]; g = [-1/math.sqrt(3), 1/math.sqrt(3)]
    def lumped(els, rho_):
        mm = {}
        for el in els:
            p = np.array([N[i] for i in el]); vol = 0
            for a in g:
                for b2 in g:
                    for c3 in g:
                        dN = np.array([[s_[0]*(1+s_[1]*b2)*(1+s_[2]*c3), s_[1]*(1+s_[0]*a)*(1+s_[2]*c3), s_[2]*(1+s_[0]*a)*(1+s_[1]*b2)] for s_ in sx]) / 8; vol += np.linalg.det(p.T @ dN)
            for i in el: mm[i] = mm.get(i, 0) + rho_ * vol / 8
        return mm
    mr, mf = lumped(dom['Ring'], 1000.0), lumped(dom['Fluid'], meta['rho_f'])
    rb = AR.read_blocks(run + '/ring_nodes.txt'); fnb = AR.read_blocks(run + '/fluid_nodes.txt')
    ids_r, ids_f = sorted(mr), sorted(mf)
    def KE(k):
        rn = np.array(rb[k]['rows']); fn = np.array(fnb[k]['rows'])
        kr = 0.5 * sum(mr[i] * (r_[3]**2 + r_[4]**2) for i, r_ in zip(ids_r, rn)); kf = 0.5 * sum(mf[i] * ((r_[3]+r_[5])**2 + (r_[4]+r_[6])**2) for i, r_ in zip(ids_f, fn))
        return kr, kf
    idx = [int(round(k * T / 2 / dtm)) for k in range(0, int(4 * 2) + 1) if int(round(k * T / 2 / dtm)) < len(rb)]
    ke = [KE(k) for k in idx]; tot = [a + b2 for a, b2 in ke]
    rep['KE0_ring_J'] = ke[0][0]; rep['KE0_fluid_J'] = ke[0][1]
    rep['KE_total_at_zero_crossings_over_initial'] = [round(x / tot[0], 4) for x in tot]
    rep['KE_total_spread_after_T'] = (max(tot[2:]) - min(tot[2:])) / tot[0] if len(tot) > 3 else None
json.dump(rep, open(run + '/report.json', 'w'), indent=1)
for k in ('run', 'mesh', 'steps_completed', 'steps_expected', 'cores', 'wall_s', 'omega_fit', 'omega_zc', 'err_fit_pct', 'window_spread_pct', 'amp_loss_per_period', 'vel_vs_dudt_rms_pct',
          'wet_dry_ratio', 'ratio_err_pct', 'mu_eff', 'ma_eff_ratio_layer0', 'ma_eff_ratio_layer1', 'p2_phase_deg', 'p_noise_over_analytic', 'p_max_first5_Pa', 'p_max_after5_Pa', 'volume_meanJm1_max',
          'KE_total_at_zero_crossings_over_initial', 'jacobian_or_nan_messages', 'max_refs_reached_events'):
    if k in rep: print(f'  {k}: {rep[k]}')
