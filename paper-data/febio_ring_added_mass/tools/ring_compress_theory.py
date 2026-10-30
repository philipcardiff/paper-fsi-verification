#!/usr/bin/env python3
"""Compressibility correction to the n=2 added mass of the ring annulus (inner radius a=Ro=1.05 m interface, rigid wall b=1.5 m).
Inviscid, linear acoustics: p'' + p'/r - n^2 p/r^2 + k^2 p = 0, k = w/c, c = sqrt(K/rho_f); rigid wall dp/dr(b)=0;
at the interface dp/dr(a) = rho w^2 u_r (ring radial displacement u_r).  The 'added pressure per unit displacement' at r=a:
  Z(k) = p(a)/u_r = rho w^2 * f(k),  f = [A J + B Y](a) / d/dr[...](a).   Incompressible limit: f0 = f(k->0).
Added-mass-like quantity m_a ~ Z/w^2 = rho f(k).  Relative change d = f(k)/f0 - 1 (> 0 ... sign printed).
Then the effect on the wet/dry frequency ratio via  w_wet^2 = k_s/(m_s + m_a):  d ln(ratio) = -1/2 * (m_a/(m_s+m_a)) * d ,
with m_a/(m_s+m_a) = 1 - (w_wet/w_dry)^2 (taken from the exact wet/dry ratio).  Uses scipy (Bessel J,Y).
"""
import numpy as np
from scipy.special import jv, yv, jvp, yvp
a, b, n = 1.05, 1.5, 2
refs = {"weak": (14, 2.4356148807, 0.9532581829), "moderate": (140, 1.8045131724, 0.7062557226), "strong": (1400, 0.7686544910, 0.3008382767)}
def f(k):
    # p = J(kr)+C Y(kr) with dp/dr(b)=0 -> C = -J'(kb)/Y'(kb)
    C = -jvp(n, k * b) / yvp(n, k * b)
    p = lambda r: jv(n, k * r) + C * yv(n, k * r); dp = lambda r: k * (jvp(n, k * r) + C * yvp(n, k * r))
    return p(a) / dp(a)
f0 = f(1e-4) * (1 + 0)   # k -> 0 limit (k small but nonzero for stability)
print("incompressible f0 =", f0)
for K in (1e7, 1e8, 1e9):
    for name, (rho, w, ratio) in refs.items():
        c = np.sqrt(K / rho); k = w / c; d = f(k) / f0 - 1
        frac = 1 - ratio ** 2; dr = -0.5 * frac * d
        print(f"K={K:.0e} {name:9s} c={c:8.1f} m/s  k*b={k*b:.3e}  (w b/c)^2={(k*b)**2:.3e}  d(added mass)/m_a={d:+.3e}  m_a/(m_s+m_a)={frac:.4f}  d(ratio)/ratio={dr:+.3e} ({100*dr:+.5f} %)")
