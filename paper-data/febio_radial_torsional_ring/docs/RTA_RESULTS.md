# RT-A: structural law, rotation expression and static finite-strain torsion

Scope: the solid only, **static**. Inputs `../inputs/final/RTA_*.feb`; raw outputs `../results/raw/RTA_*`; figure `../figures/fig1_rta_static_twist.png`. Numbers below are direct FEBio output reduced by the campaign script (`scripts/rtr_check_rta.py`) and re-derived with independent code (`scripts/independent_check.py`); the two agree to all printed digits (the twist reference is the same vendored function).

## Law identity (job 14189, `RTA_law1c`)
2x2 hex8 block (2 x 1 x 1 m), plane strain, all boundary nodes prescribed to F_k = I + 0.2 k (Fm - I), Fm = [[1.1, 0.4, 0], [0.1, 0.9, 0], [0, 0, 1]], k = 1..5 (J = 0.998 at k = 1 down to 0.95 at k = 5), FEBio `Mooney-Rivlin` c1 = 3846.154, c2 = 0, k = 16666.67, NIKE3D, three-field hex8 domain. FEBio Cauchy stress (all six components, all elements, all steps) against (a) the hand formula sigma = mu/J dev(bbar) + K/2 (J^2 - 1)/J I and (b) a central finite-difference derivative sigma = (1/J) dW/dF F^T of the proposal's strain energy:

| Comparison | max abs difference | relative to max abs stress (3,938 Pa) | gate |
|---|---|---|---|
| hand formula | 4.8e-9 Pa | **1.2e-12** | 1e-8 |
| finite-difference dW/dF | 1.3e-6 Pa | 3.3e-10 (limited by the derivative's truncation) | not gated |

Superseded designs (not evidence): jobs 14186 and 14188 left the end faces free, so a homogeneous F is not an equilibrium state; the mismatch they showed was the test design, not the law.

## Rotation expression (job 14187, `RTA_rotcheck`)
One hex8 under the exact rotation expressions with the C-infinity ramp, amplitude 0.5 rad, 40 steps of 0.05 s: maximum node-displacement difference to the closed-form rotation 7.0e-13 m (gate 1e-10 m). It confirms the `math` boundary condition reads the end-of-step time with no offset.

## Static twist (jobs 14193 M1, 14195 M2, 14196 M3)
Annulus 0.5 <= R <= 1 m, hex8, nth x nrs = 64x8, 128x16, 256x32, inner circle rotated by alpha_A(t) = 0.5 S(t) sin(3t) (both in-plane components prescribed as a rigid rotation: the radially clamped condition of `EQUIVALENCE.md` section 2), outer circle prescribed to zero displacement (so that FEBio logs its reaction force; reaction forces exist only for prescribed-displacement BCs), static, 20 steps of 0.1 s. FEBio quantity: drive torque T_A = sum (x R_y - y R_x)/(2 pi dz) over the inner nodes (reaction sign reversed: FEBio's reaction is minus the load the constraint applies). Reference: the proposal's discrete finite-strain energy (vendored `rt_reference`, hash-checked) solved statically on 24 GLL nodes (residual 1e-11).

| Mesh | E_max(T_A) | observed order | E_max u_r(R = 0.75) | E_max alpha(0.75) | (T_inner + T_outer)/max T_inner | ray 0 vs 45 deg |
|---|---|---|---|---|---|---|
| M1 | 3.396e-3 | 1.95 (M1-M2) | 1.76e-3 | 3.11e-3 | 1.8e-12 | 7e-11 |
| M2 | 8.779e-4 | 1.99 (M2-M3) | 5.69e-4 | 8.08e-4 | 4.3e-13 | 5e-11 |
| M3 | 2.215e-4 | | 1.51e-4 | 2.04e-4 | 3.4e-13 | 4e-11 |

Peak |T_A|: 2,551.05 (reference), 2,559.7 / 2,553.3 / 2,551.6 (M1 / M2 / M3). Gate: E_max(T_A) at M3 at most 1e-3 and order (M2, M3) at least 1.8: **PASS** (all three T_A values are from the log of the final run; the first attempts 14190/14191 failed on a duplicate node-set name in the input, 14192 omitted the outer reaction log).

**Classification: VERIFIED AGREEMENT for the static solid (law, rotation expression, finite-strain torsion stiffness), time not applicable.** Evidence: law identity and rotation expression PROVEN to round-off; mesh-converged torsion stiffness SUPPORTED (three meshes at second order; the reference is the proposal's own, but the law identity is checked independently of it). Not covered: the dynamic solid, any coupled behaviour.

Note: the campaign file `check.json` of the final M1 run stored an equilibrium value of 1.0 (stale, computed from the earlier attempt that did not log the outer reaction); the recomputed value from the final data is 1.8e-12 (`../provenance/INCONSISTENCIES.md`).
