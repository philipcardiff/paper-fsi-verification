#!/usr/bin/env python3
"""Audit item G: n=2 Fourier component of the interface-adjacent pressure vs the added-mass relation.

Linear potential-flow theory: p(R,theta,t) = m_a * d2/dt2 u_r(R,theta,t),  u_r = D(t) cos 2theta
   =>  P2(t) = m_a * Ddd(t)   with P2 the cos2theta Fourier coefficient of p on the interface.
We project the first fluid layer element pressures (j = 0..47, centred at theta=(j+1/2)dtheta) onto cos 2theta,
fit both P2(t) and D(t) with a harmonic at the measured frequency, and report the effective added mass
m_a,eff = |P2 amplitude| / (omega^2 |D amplitude|), plus the sign/phase relation and the integrated radial force.
Only post-processing of existing logs; nothing is run.
"""
import sys, math, json
import numpy as np

run, w_fit = sys.argv[1], float(sys.argv[2])
NT, NR = 48, 12
def blocks(p):
    out, cur = [], None
    for line in open(p):
        if line.startswith('*Step'): cur = {'rows': []}; out.append(cur)
        elif line.startswith('*Time'): cur['t'] = float(line.split('=')[1])
        elif line.startswith('*'): pass
        elif line.strip(): cur['rows'].append([float(x) for x in line.split()])
    return out
fb = blocks(run + '/fluid_elems.txt'); t = np.array([b['t'] for b in fb])
P = np.array([[r[1] for r in b['rows']] for b in fb]); J = np.array([[r[2] for r in b['rows']] for b in fb])
pb = blocks(run + '/probes.txt'); D = np.array([0.5 * (b['rows'][0][3] - b['rows'][1][4]) for b in pb])
dth = (math.pi / 2) / NT
layer = lambda i: P[:, [j * NR + i for j in range(NT)]]       # element k = j*NR + i
c2 = np.cos(2 * (np.arange(NT) + 0.5) * dth)
def P2_of(i): return (layer(i) @ c2) * dth / (math.pi / 4)     # (1/norm) * int p cos2theta dtheta, norm = pi/4 on the quarter
def harm(x, w):
    M = np.column_stack([np.sin(w * t), np.cos(w * t), np.ones_like(t)]); c, *_ = np.linalg.lstsq(M, x, rcond=None)
    r = x - M @ c; return c, r.std()
rho, Ro, b, h = 140.0, 1.05, 1.5, 0.1
ma = rho * Ro / 2 * (b**4 + Ro**4) / (b**4 - Ro**4)
cD, rD = harm(D, w_fit); ampD = math.hypot(cD[0], cD[1])
print(f"run {run}: omega_fit {w_fit}; D amplitude {ampD:.4e} m (fit rms resid {rD:.1e})")
print(f"analytic added mass m_a = {ma:.3f} kg/m2 ; predicted P2 amplitude at the measured omega = m_a*omega^2*D = {ma*w_fit**2*ampD:.4f} Pa")
res = {}
for i in (0, 1, 2):
    c, r = harm(P2_of(i), w_fit); amp = math.hypot(c[0], c[1])
    # phase relative to D: P2 should equal -m_a w^2 D  (anti-phase with D)
    ph = math.degrees(math.atan2(c[1], c[0]) - math.atan2(cD[1], cD[0])) % 360
    # sign: projection of P2 on -D
    proj = (c[0] * cD[0] + c[1] * cD[1]) / (ampD ** 2)          # P2 = proj * D  (expect -m_a w^2)
    # correct to the interface: dp/dr = -rho * a_r ; layer centre at dr*(i+0.5) from the interface (dr = 0.45/12)
    dr = (b - Ro) / NR * (i + 0.5)
    print(f"  fluid layer {i}: P2 amplitude {amp:.4f} Pa (noise rms of P2 after harmonic fit {r:.4f}); P2/D = {proj:+.1f} Pa/m -> m_a,eff = {-proj/w_fit**2:.2f} kg/m2 = {-proj/w_fit**2/ma:.3f} x analytic;  phase vs D {ph:.1f} deg;  centre-offset correction ~ rho*omega^2*D*dr/(P2) = {rho*w_fit**2*ampD*dr/amp if amp>0 else float('nan'):.3f}")
    res[i] = -proj / w_fit**2 / ma
# integrated radial force (quarter, per layer-0 pressure) F = int p * (outward normal force on ring is -p) ; cos2theta-weighted
c0, _ = harm(P2_of(0), w_fit)
Fr = c0 * (Ro * h * math.pi / 4)      # coefficient of the n=2 radial force resultant (N, quarter, cos2theta weighted)
print(f"  integrated n=2 radial loading amplitude on the ring (layer-0 pressure): {math.hypot(Fr[0],Fr[1]):.3e} N ; analytic m_a w^2 D R h pi/4 = {ma*w_fit**2*ampD*Ro*h*math.pi/4:.3e} N")
# mass / volume conservation
dV = (J - 1).mean(axis=1)
print(f"  fluid volume: mean(J-1) over elements, max |.| over time {np.abs(dV).max():.2e}; max |J-1| any element {np.abs(J-1).max():.2e}")
print(f"  mean pressure over all elements at each time: max |mean p| {np.abs(P.mean(axis=1)).max():.3f} Pa ; p at t=first nonzero block max|p| = {np.abs(P[1]).max():.3f} Pa (shock indicator)")
print(f"  max |p| in first 5 steps {np.abs(P[:6]).max():.2f} Pa vs max |p| thereafter {np.abs(P[6:]).max():.2f} Pa")
