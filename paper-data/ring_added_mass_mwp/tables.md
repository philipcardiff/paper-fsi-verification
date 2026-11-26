# Ring added-mass problem: rerun with `movingWallPressure` on the IQN-ILS interface

Build: solids4foam development 7a2f0a63e + #721 fix (branch fix/movingWallPressure-721), OpenFOAM v2512, UCD Sonic (Linux, GCC), `./Allverify --cores 16`; tutorial `ringAddedMass` with `0/fluid/p.dirichletNeumann` interface `movingWallPressure` (viscousTerm default yes). Old values: zeroGradient interface, solids4foam README (v2412, xenosim), as in the paper.

## Table tab:ring (high-order solid, IQN-ILS)

| Level | Mesh | Steps | ω new | Error new | Ratio error new | Loss/period new | ω old | Error old | Ratio error old | Loss old |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| dry | 1 | 100 | 2.55254 | -0.098% | – | 0.039% | 2.55254 | -0.098% | – | 0.04% |
| dry | 2 | 100 | 2.55222 | -0.110% | – | 0.039% | – | – | – | – |
| dry | 4 | 100 | 2.55173 | -0.129% | – | 0.039% | 2.55173 | -0.129% | – | 0.04% |
| dry | 2 | 25 | 2.50576 | -1.929% | – | 2.176% | 2.50576 | -1.929% | – | 2.18% |
| dry | 2 | 50 | 2.54240 | -0.495% | – | 0.301% | – | – | – | – |
| weak | 1 | 100 | 2.43309 | -0.104% | -0.006% | 0.034% | 2.43341 | -0.091% | +0.007% | 0.25% |
| weak | 2 | 100 | 2.43289 | -0.112% | -0.002% | 0.039% | 2.43298 | -0.108% | +0.002% | 0.15% |
| weak | 4 | 100 | 2.43245 | -0.130% | -0.000% | 0.040% | 2.43248 | -0.129% | +0.001% | 0.10% |
| weak | 2 | 25 | 2.38861 | -1.930% | -0.001% | 2.176% | 2.38894 | -1.916% | +0.013% | 2.30% |
| weak | 2 | 50 | 2.42353 | -0.496% | -0.002% | 0.302% | 2.42378 | -0.486% | +0.009% | 0.41% |
| moderate | 1 | 100 | 1.80219 | -0.129% | -0.031% | 0.017% | 1.80344 | -0.059% | +0.038% | 1.18% |
| moderate | 2 | 100 | 1.80237 | -0.119% | -0.009% | 0.041% | 1.80273 | -0.099% | +0.011% | 0.64% |
| moderate | 4 | 100 | 1.80214 | -0.132% | -0.002% | 0.050% | 1.80224 | -0.126% | +0.003% | 0.35% |
| moderate | 2 | 25 | 1.76957 | -1.936% | -0.008% | 2.179% | 1.77093 | -1.861% | +0.069% | 2.85% |
| moderate | 2 | 50 | 1.79542 | -0.504% | -0.009% | 0.305% | 1.79645 | -0.447% | +0.048% | 0.92% |
| strong | 1 | 100 | 0.76743 | -0.159% | -0.061% | -0.001% | 0.76843 | -0.029% | +0.068% | 2.11% |
| strong | 2 | 100 | 0.76768 | -0.126% | -0.016% | 0.056% | 0.76796 | -0.090% | +0.021% | 1.14% |
| strong | 4 | 100 | 0.76763 | -0.134% | -0.004% | 0.087% | 0.76770 | -0.124% | +0.006% | 0.64% |
| strong | 2 | 25 | 0.75371 | -1.944% | -0.016% | 2.196% | 0.75477 | -1.806% | +0.125% | 3.41% |
| strong | 2 | 50 | 0.76472 | -0.512% | -0.017% | 0.322% | 0.76552 | -0.408% | +0.087% | 1.43% |

## Orders

| Series | new values | new order | old order |
|---|---|---:|---:|
| dry: frequency, time step (m2) | 2.505764, 2.542401, 2.552224 | 1.90 | 1.90 |
| dry: Richardson (order 2, n50/n100) error | 2.555498 | +0.018% | +0.018% max |
| dry: standard solid frequency error, mesh 1/2/4 | +20.66%, +6.17%, +1.58% | 1.66 | 1.66 |
| weak: frequency, time step (m2) | 2.388606, 2.423527, 2.432890 | 1.90 | 1.92 |
| weak: Richardson (order 2, n50/n100) error | 2.436011 | +0.016% | +0.018% max |
| weak: ratio error, mesh 1/2/4 | -0.006%, -0.002%, -0.000% | 1.90 | 1.77 |
| weak: standard solid frequency error, mesh 1/2/4 | +20.65%, +6.16%, +1.57% | 1.66 | 1.66 |
| weak: standard solid ratio error, mesh 1/2/4 | -0.0039%, -0.0007%, -0.0001% | 2.56 | – |
| moderate: frequency, time step (m2) | 1.769574, 1.795421, 1.802369 | 1.90 | 2.02 |
| moderate: Richardson (order 2, n50/n100) error | 1.804685 | +0.010% | +0.018% max |
| moderate: ratio error, mesh 1/2/4 | -0.031%, -0.009%, -0.002% | 1.84 | 1.73 |
| moderate: standard solid frequency error, mesh 1/2/4 | +20.63%, +6.16%, +1.57% | 1.66 | 1.66 |
| moderate: standard solid ratio error, mesh 1/2/4 | -0.0233%, -0.0044%, -0.0009% | 2.42 | – |
| strong: frequency, time step (m2) | 0.753712, 0.764722, 0.767683 | 1.89 | 2.14 |
| strong: Richardson (order 2, n50/n100) error | 0.768670 | +0.002% | +0.018% max |
| strong: ratio error, mesh 1/2/4 | -0.061%, -0.016%, -0.004% | 1.90 | 1.68 |
| strong: standard solid frequency error, mesh 1/2/4 | +20.60%, +6.16%, +1.57% | 1.66 | 1.66 |
| strong: standard solid ratio error, mesh 1/2/4 | -0.0450%, -0.0089%, -0.0018% | 2.36 | – |

## Coupling (mesh 2, 100 steps per period, high-order solid)

| Level | Coupling | ω | Error | Loss/period | Iter./step |
|---|---|---:|---:|---:|---:|
| weak | iqnils | 2.432890 | -0.112% | 0.039% | 2.11 |
| weak | robin | 2.432894 | -0.112% | 0.039% | 3.84 |
| moderate | iqnils | 1.802369 | -0.119% | 0.041% | 2.38 |
| moderate | robin | 1.802375 | -0.118% | 0.040% | 4.80 |
| strong | iqnils | 0.767683 | -0.126% | 0.056% | 2.75 |
| strong | robin | 0.767686 | -0.126% | 0.056% | 10.21 |

## Amplitude loss per period and decay rate (IQN-ILS, high-order solid)

Exact viscous decay rates α_V: weak 1.005e-4, moderate 4.748e-4, strong 5.621e-4 1/s. BDF2 (solids4foam backward d2dt2) loss per period: 2.174%, 0.301%, 0.0386% at 25, 50, 100 steps.

| Level | Mesh | Steps | Loss/period | Decay rate (1/s) | Decay/α_V |
|---|---:|---:|---:|---:|---:|
| weak | 1 | 100 | 0.034% | 1.316e-04 | 1.31 |
| weak | 2 | 100 | 0.039% | 1.498e-04 | 1.49 |
| weak | 4 | 100 | 0.040% | 1.552e-04 | 1.54 |
| weak | 2 | 25 | 2.176% | 8.364e-03 | 83.23 |
| weak | 2 | 50 | 0.302% | 1.165e-03 | 11.59 |
| moderate | 1 | 100 | 0.017% | 4.910e-05 | 0.10 |
| moderate | 2 | 100 | 0.041% | 1.166e-04 | 0.25 |
| moderate | 4 | 100 | 0.050% | 1.424e-04 | 0.30 |
| moderate | 2 | 25 | 2.179% | 6.204e-03 | 13.07 |
| moderate | 2 | 50 | 0.305% | 8.715e-04 | 1.84 |
| strong | 1 | 100 | -0.001% | -1.304e-06 | -0.00 |
| strong | 2 | 100 | 0.056% | 6.880e-05 | 0.12 |
| strong | 4 | 100 | 0.087% | 1.062e-04 | 0.19 |
| strong | 2 | 25 | 2.196% | 2.664e-03 | 4.74 |
| strong | 2 | 50 | 0.322% | 3.928e-04 | 0.70 |

## Coupling-tolerance sweep (tab:tol_sweeps; strong, mesh 2, 100 steps per period, standard solid, strict test, couplingReuse 0)

Dry ring at the same discretisation (Sonic run): ω = 2.7125759973. Robin–Neumann columns unchanged from the paper (elasticWallPressure is not affected).

| Tol. | IQN-ILS ratio error new | IQN-ILS loss/period new | Iter./step new | IQN-ILS ratio error old | Loss old | Iter. old | RN ratio error (paper) | RN iter. |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 1e-4 | -0.00884% | 0.0630% | 4.70 | +0.03316% | 1.148% | 4.68 | -0.00949% | 2.35 |
| 1e-5 | -0.00886% | 0.0630% | 5.03 | +0.03315% | 1.148% | 5.00 | -0.00853% | 3.78 |
| 1e-6 | -0.00886% | 0.0630% | 5.27 | +0.03315% | 1.148% | 5.24 | -0.00851% | 11.73 |
| 1e-7 | -0.00886% | 0.0630% | 6.01 | +0.03315% | 1.148% | 5.91 | -0.00851% | 17.32 |
| 1e-8 | -0.00886% | 0.0630% | 7.14 | +0.03315% | 1.148% | 6.96 | -0.00851% | 23.48 |
