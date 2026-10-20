# Quasi-monolithic cross-checks: mokLidDrivenCavity, cavityFlexibleBottom, blobInTreacle (2026-10-16)

After the ringAddedMass study (`ringAddedMass_monolithic_findings.md`), the
corrected quasi-monolithic Newton coupling was run on three further cases as a
cross-check of the partitioned results. The partitioned results of the paper
are unchanged; the monolithic results are a same-code, different-algorithm
comparison, not an independent reference.

## Provenance

- solids4foam branch `feature-ring-monolithic-verification` (PR #582), on
  `feature-petsc-snes-quasi-monolithic-update` (PR #580):
  - `ba69872bd`: verification drivers run the quasi-monolithic coupling;
  - `dbfda65ba`: interface traction is the momentum equation's wall flux
    (third correction, below);
  - `c0c91c76b`, `22c123755`: mokLidDrivenCavity `verification/README.md`
    (monolithic mesh study, 128 x 128 level);
  - `291aaa968`: cavityFlexibleBottom `verification/README.md`;
  - earlier corrections: BDF2 solid `c008e8dd1`, interface pressure gradient
    `3e93daef6`.
- OpenFOAM v2512, macOS (Apple M-series), serial unless stated; execution
  times indicative only.
- Data: `mokLidDrivenCavity_monolithic.csv`,
  `cavityFlexibleBottom_monolithic.csv`, `blobInTreacle_monolithic.csv`.

## Third correction: interface viscous traction (found on blobInTreacle)

- The solver passed the viscous traction `n & devReff` to the solid, but the
  fluid momentum equation contains only `laplacian(nuEff, U)`, whose wall flux
  is `nuEff*snGrad(U)`: the traction contained a first-order force that the
  fluid's momentum balance does not. The traction is now the momentum
  equation's wall flux (`dbfda65ba`); blobInTreacle then converges at second
  order (below).
- Rejected alternative: a quadratic "second-order" wall-shear estimate through
  the cell values. Applied as an experiment to the partitioned `pimpleFluid`,
  it degraded the blobInTreacle mesh convergence to first order (errors
  against 0.1204 m on x1/x2/x4: 3.1e-3, 0.70e-3, 0.18e-3 m before;
  3.3e-3, 2.7e-3, 1.5e-3 m with it). In a converged conservative finite volume
  solution the one-sided wall flux is consistent with the momentum balance and
  converges at second order (supraconvergence); a fit through cell values is
  not conservative.
- Lesson: the traction passed to the solid must be built from exactly the
  operators and boundary values of the fluid's momentum equation.
- Recorded in solids4foam issue #583. The partitioned results of the paper
  are not affected (the experiment was not adopted).

## mokLidDrivenCavity

Same members as the partitioned mesh study; the monolithic velocity solid has
the standard residual (the high-order setting is ignored). Partitioned and
standard/high-order solids agree to 0.04%, so the solid difference is
immaterial.

| Mesh (Δt) | Peak mono | Trough mono | Mean mono | Time (s) | Peak part | Trough part | Mean part |
|---|---:|---:|---:|---:|---:|---:|---:|
| 16² (0.2) | 0.2950 | 0.2053 | 0.2495 | 16 | 0.2934 | 0.1942 | 0.2426 |
| 32² (0.1) | 0.2903 | 0.2027 | 0.2460 | 96 | 0.2890 | 0.1986 | 0.2427 |
| 64² (0.05) | 0.2877 | 0.2026 | 0.2447 | 606 | 0.2870 | 0.2022 | 0.2437 |
| 96² (0.025) | 0.2867 | 0.2025 | 0.2441 | 2499 | 0.2861 | 0.2044 | 0.2445 |
| 128² (0.0125) | 0.2860 | 0.2025 | 0.2438 | 9109 | diverged | | |

- Monolithic trough mesh-independent to 0.1% from 32² (0.2025 m); the
  partitioned trough is still rising and passes it on 96² (0.9% above). The
  slow trough convergence is a property of the partitioned solution, not of
  the benchmark. Why the partitioned trough moves past the monolithic value is
  not explained *(open)*.
- Peaks agree to 0.2% on 96². Monolithic peak changes 0.0047, 0.0026 m
  (16/32/64; order 0.85 *(derived here)*), then 0.0017 m from 64² to 128²:
  order below 1, no Richardson estimate.
- Means approach ≈0.244 m from either side.
- 128² (diverged with partitioned coupling) runs to completion with the
  monolithic coupling: the divergence is specific to the partitioned coupling.
- On 128² the monolithic trough and mean lie inside the group-B envelope; the
  peak is 3.3% above its maximum *(derived here)*.

## cavityFlexibleBottom

Steady u_y at (4, −1) m and F_y (N, 1 m depth). Monolithic at Δt = 0.5 s with
the inlet ramped over 20 s (steady state unaffected); partitioned Aitken local
rerun at Δt = 0.5 s (level 3 at 0.25 s, at which Aitken converges).

| Level | Δx | u_y mono | u_y part | F_y mono | F_y part |
|---:|---:|---:|---:|---:|---:|
| 1 | 0.1 | −0.201557 | −0.202495 | −5.1225 | −5.1527 |
| 2 | 0.05 | −0.228817 | −0.228799 | −5.2381 | −5.2378 |
| 3 | 0.025 | −0.235676 | −0.235628 | −5.2576 | −5.2564 |
| 4 | 0.0125 | −0.237336 | – | −5.2606 | – |

- The couplings agree to 0.02% in u_y and F_y on level 3.
- Monolithic u_y: observed orders 1.99, 2.05; extrapolated −0.23787 m.
- The old 24% monolithic/partitioned difference (`simulationData/cavityFlexibleBottom`,
  −0.29192 m at Δx = 0.025, pre-dating the solver corrections) is resolved:
  those data remain unused.
- Tuković et al. (2018) mesh-study values are about 5% larger in magnitude
  (−0.2496 m at 0.025); forces agree with them to 0.1–0.2%.
- Independent ANSYS calculation (E. Tandis, paper PR #12,
  `AnsysResults/cavityFlexibleBottom-Ansys-*.csv` on that PR's branch):
  u_y ≈ −0.250, −0.2535, −0.2525 m (finest not yet steady), but its horizontal
  displacement at the monitor point (+0.0325 m against −1.7e-4 m here)
  suggests a different monitor point or definition; questions open on PR #12.
  The remaining 5% is a case-definition question, not a discretisation or
  coupling error of either solids4foam coupling.

## blobInTreacle

Steady top-point u_x (m), levels x1/x2/x4, time steps of the mesh study.

| Coupling | x1 | x2 | x4 | Order | Extrapolated |
|---|---:|---:|---:|---:|---:|
| partitioned (macOS rerun) | 0.117301 | 0.119696 | 0.120224 | 2.2 | 0.12037 |
| partitioned (paper, Linux) | 0.117301 | 0.119698 | 0.120152 | 2.40 | (0.12026) |
| monolithic (final) | 0.117429 | 0.119524 | 0.120053 | 2.0 | 0.12023 |

- Same-platform limits agree to 0.12%; monolithic within 0.15% of the paper's
  partitioned values on every level.
- The partitioned x4 value differs by 0.06% between the Linux and macOS
  runs *(open: record which run the paper quotes and why they differ)*.
- This replaces the inventory note "Not cross-checked against the
  quasi-monolithic solver after the StVK change".
