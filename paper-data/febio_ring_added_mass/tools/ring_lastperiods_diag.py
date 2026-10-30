#!/usr/bin/env python3
"""Part 1b: residual / harmonic / transient diagnostics of stored ring displacement histories (post-processing only)."""
import sys, os, json, math, csv
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import analyze_ring as AR
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RUNS = {'dry_M2': 'runs/phaseB/D1_dry_f2', 'weak_M2': 'runs/phaseD/D4_weak_M2', 'moderate_M2': 'runs/phaseC/C2_Emesh1e-3_wet_f2', 'strong_M2': 'runs/phaseD/D3_strong_M2'}
rows = []; res_store = {}
for name, run in RUNS.items():
    t, d, v = AR.series(run + '/probes.txt'); rep = json.load(open(run + '/report.json'))
    w = rep['omega_fit']; sdec = -math.log(1 - rep['amp_loss_per_period']) / (2 * math.pi / w)
    e = np.exp(-sdec * t)
    # fundamental-only damped fit, residual
    M1 = np.column_stack([e * np.sin(w * t), e * np.cos(w * t)]); c1, *_ = np.linalg.lstsq(M1, d, rcond=None); r1 = d - M1 @ c1
    A = math.hypot(*c1)
    # + harmonics 2,3 and DC (diagnostic of nonlinearity)
    cols = [e * np.sin(w * t), e * np.cos(w * t)]
    for k in (2, 3): cols += [np.sin(k * w * t), np.cos(k * w * t)]
    cols += [np.ones_like(t)]
    M = np.column_stack(cols); c, *_ = np.linalg.lstsq(M, d, rcond=None); r2 = d - M @ c
    h2 = math.hypot(c[2], c[3]) / A; h3 = math.hypot(c[4], c[5]) / A; dc = c[6] / A
    # best-sinusoid scan of fundamental residual
    fr = np.linspace(0.05, 30, 6000); amp = []
    for f in fr:
        Mx = np.column_stack([np.sin(f * w * t), np.cos(f * w * t)]); cc, *_ = np.linalg.lstsq(Mx, r1, rcond=None); amp.append(math.hypot(*cc))
    amp = np.array(amp); i = int(np.argmax(amp))
    # amplitude evolution: envelope = sqrt(d^2 + (v/w)^2) using velocity of probes
    env = np.sqrt(d ** 2 + (v / w) ** 2); n = len(t)
    spp = int(round(n - 1) / 4)
    envp = [float(np.mean(env[k * spp:(k + 1) * spp + 1])) for k in range(4)]
    rows.append(dict(case=name, omega=w, A_fund=A, rel_rms_fund_only=float(np.sqrt(np.mean(r1 ** 2))) / A, rel_rms_with_h2_h3_dc=float(np.sqrt(np.mean(r2 ** 2))) / A,
                     h2_over_A=h2, h3_over_A=h3, dc_over_A=dc, resid_peak_amp_over_A=float(amp[i]) / A, resid_peak_freq_over_omega=float(fr[i]),
                     env_mean_P1=envp[0], env_mean_P2=envp[1], env_mean_P3=envp[2], env_mean_P4=envp[3], decay_rate=sdec,
                     duffing_shift_bound_from_h3=12 * h3))
    res_store[name] = (t, r1 / A, fr, amp / A)
with open('results/ring_lastperiods_diagnostics.csv', 'w', newline='') as fh:
    wr = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); wr.writeheader()
    for r in rows: wr.writerow({k: (f'{x:.8g}' if isinstance(x, float) else x) for k, x in r.items()})
for r in rows: print({k: (f'{x:.3g}' if isinstance(x, float) else x) for k, x in r.items()})
np.savez('results/ring_lastperiods_resid.npz', **{f'{k}_{j}': a for k, v in res_store.items() for j, a in zip(('t', 'r', 'f', 'amp'), v)})
