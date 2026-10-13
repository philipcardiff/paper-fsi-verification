# beamInCrossFlow: Richter (2012) 3-D steady FSI benchmark — findings

Source of truth: solids4foam branch `verification/beam-richter` (also pushed
to PR #544's branch `verification/beam-crossflow-graded-mesh`), file
`tutorials/fluidSolidInteraction/beamInCrossFlow/verification/README.md`, and
the compact evidence in `.../verification/reference/campaign/`. Data copied
here: `paper-data/beam_richter/`. Status: L0-L3 complete (2026-10-13).

## Platform and build

MeluXina, OpenFOAM v2412 (foss-2024a, GCC 13.3), PETSc 3.22.0, frozen copy of
the solids4foam build `aee0c8e350267` (-O2). One build for all definitive runs.

## Definition

Richter (2012) section 7.2: duct (0,1.5)x(0,0.4)x(-0.4,0.4) m, z = 0
symmetry half; St Venant-Kirchhoff beam (0.4,0.5)x(0,0.2)x(-0.2,0.2),
E = 1.4 MPa, nu = 0.4; fluid rho = 1000, nu = 1e-3; printed bi-parabolic inlet
with peak 0.3 m/s (area mean 0.133 m/s; Richter's "mean 0.2, Re = 40" is
inconsistent with the printed profile). Point A = (0.45, 0.15, -0.15) on the
computed half; F = traction integral on the half interface. Transient-to-steady
partitioned continuation (IQN-ILS, PIMPLE, velocity-Laplacian ALE, Newton solid).
Caveat: an earlier 0.2 m/s peak sequence landed near Richter's Table 6; the
comparison with Richter is made only after our own sequence is resolved.

## Key findings

1. Coupling stall of the October attempt: caused by PIMPLE `nOuterCorrectors 1`
   (history-dependent fluid map), not by the coupling tolerance; 3 correctors
   meet the strict test in 8-15 iterations and are 1.8x faster.
2. Near the steady state the relative measure r1 is unattainable (r2 ~ 1e-14).
   Production: strict max(r1, r2) <= 1e-5, cap 30 (counted). Coupling
   tolerances 1e-4..1e-8, legacy test, PIMPLE 6 and SNES 1e-9 change u_x(A) and
   F_x by <= ~2e-6 (L0-L2).
3. With fluid solver tolerance 1e-6 a workstation L0 run froze on
   history-dependent states (3.6e-4 in F_x); p/U 1e-10 removes it. Not
   reproduced on MeluXina at L1. History independence (rest vs rigid-flow
   start) 1.3e-6 (L0), 1.6e-6 (L1).
4. Error separation: u_x(A) error is dominated by through-thickness solid
   resolution (static diagnostic: 26/12/4.5/1.6 % too stiff for N_t = 4f). The
   production family refines the solid 2x in thickness (N_t = 8f).
5. The steady state depends on dt through Rhie-Chow (FSI L1: ~1e-3 for dt x2
   or /2; rigid: 1e-4 (L1), 2-4e-5 (L2), 3e-6 (L3)); dt = 0.00625 fixed.

## Production sequence (L0-L3)

| level | u_x(A) (m) | u_y(A) (m) | u_z(A) (m) | F_x (N) | F_y (N) |
| --- | ---: | ---: | ---: | ---: | ---: |
| L0 | 8.13919e-5 | 1.95723e-6 | -2.03080e-6 | 2.15736 | 0.32694 |
| L1 | 9.62607e-5 | 1.94057e-6 | -2.15035e-6 | 2.30568 | 0.33773 |
| L2 | 1.024547e-4 | 1.92632e-6 | -2.20300e-6 | 2.36363 | 0.34428 |
| L3* | 1.044776e-4 | 1.91580e-6 | -2.21539e-6 | 2.38305 | 0.34577 |

*L3: from the converged rigid flow, two segments (0-3.69 s, 3-4.5 s), mean
over 4.0-4.5 s. The 7.9 Hz damped fluid-loaded bending oscillation (8.0 Hz at
L2, not shedding) decays 2.1x per 0.25 s: u_x half-amplitude 9.4e-5 (3.0-3.25
s) to 2.1e-6 (4.25-4.5 s). Remaining slow wake drift: ~2e-6 (u_x), ~3e-6 (F_x),
~1e-4 (F_y, F_z). The plain restart at t = 3 s did not re-excite the mode; the
L2 restart with re-decomposition did (1.4e-4, then decayed).

Increments L0-1/1-2/2-3: u_x 15.4/6.0/1.9 %, F_x 6.4/2.5/0.81 %. Observed
orders (L0-L2, L1-L3): u_x 1.26, 1.61; F_x 1.36, 1.58 (rigid-beam F_x 1.36,
1.55); F_y 0.72, 2.16; u_z 1.18, 2.09; u_y 0.23, 0.44. Wording: monotone,
orders increasing towards 2, asymptotic regime suggested but not
demonstrated; no GCI. Conditional Richardson values (L1-L3 order): u_x
1.0546e-4 m (+0.94 %), F_x 2.3929 N (+0.41 %); with p = 2: 1.0515e-4 (+0.65 %),
2.3895 (+0.27 %). F_z orders 0.73, 1.84.

Non-spatial errors: coupling <= 2e-6, history <= 1.6e-6, steadiness 2.5e-7
(L2), fixed-dt Rhie-Chow sensitivity (dt x2, /2) 2-4e-4 at L2, falling ~4x
per level.

## Richter comparison

Table 6: drag 1.327 N, u_x(A) 5.924e-5 m. Our printed-inlet (peak 0.3)
sequence approaches ~2.39 N and ~1.05e-4 m (1.80x, 1.78x), far beyond the
numerical corrections. Earlier peak-0.2 runs (coarser solid) gave 1.3311 N and
5.915e-5 m, close to Table 6: a possible inlet-definition mismatch in the
source, not established. No tuning to Table 6.

## Cost

73.0 node-hours = 9,350 physical core-h (node-hours x 128) in total; L3 44.6
node-hours (segment 1 36.0, segment 2 8.6).

## Astra

Design review and final interpretation in `paper-data/beam_richter/astra/`.
Open: coupled L3 dt/2 sensitivity (>= 60 node-hours) not run.
