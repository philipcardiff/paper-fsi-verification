# 3-D pressure-pulse tube: Hookean versus St Venant–Kirchhoff wall on L1

Effect of the wall law on the point-A QoIs of the 3-D pressure-pulse tube
(paper Section `sec:3dtube`, Problem paragraph and Lessons; Appendix E).

## Set-up

- L1 of the paper's mesh family (16 000 fluid, 6 400 solid cells), Δt = 2.5e-5 s,
  Robin–Neumann coupling, BDF2 in fluid and solid, end time 0.02 s.
  Both runs use the tutorial verification driver:
  `Allverify --study mesh --levels 1 --cores 4`.
- Build: solids4foam 206437983 (plus one unused outlet-BC commit), OpenFOAM v2512,
  on UCD Sonic. Slurm jobs 1192829 (Hookean) and 1192830 (StVK).
- **Hookean:** `linearElastic` with `linearGeometryTotalDisplacement` (the paper's case).
- **StVK:** `StVenantKirchhoffElastic` with
  `nonLinearGeometryTotalLagrangianTotalDisplacement` (finite strain), with the same
  E = 0.3 MPa, ν = 0.3 and ρ = 1200 kg/m³. Only `constant/solid` differs; the
  coefficient block (PETSc SNES, tolerances, predictor, momentum stabilisation) is unchanged.
- Control: the Hookean run reproduces the paper's L1 (0.15958 mm, −0.08817 mm,
  5.872 ms, 4.664 m/s) to all reported digits.
- Coupling: every step met the Robin criteria in both runs; mean 4.78 iterations per
  step in both (maximum 9).

## Results

`stvk_vs_hookean_L1.csv`. Changes are relative, (StVK − Hookean)/|Hookean|, as in the
paper's Table `tab:3dtube`. The mesh columns are the paper's L1→L2 and L2→L3 changes
(Hookean), and the ratios divide the magnitude of the law effect by the magnitude of
each mesh change.

- The wall law changes u_r,max, u_z,min, t_arr and c_p by 1.2–5 times the mesh
  changes; only c_wall is insensitive (0.04 %).
- With StVK the radial peak is lower and earlier, the axial trough deeper and the
  front earlier. Wall-front arrivals are earlier by 7–48 µs at z = 5–25 mm and by
  less at 30–40 mm, so the least-squares c_wall barely changes while t_arr does.
- The point-A histories differ by up to 5.2 % of the radial peak (8.05 ms, falling
  flank), 4.0 % of the axial peak and 2.9 % of the axis-pressure peak.
- Against Tukovic et al. (2018, StVK): the StVK L1 peak is 0.6 % above the published
  peak (Hookean L1: 1.7 %); peak time 7.14 ms against 7.23 ms. One level only.

## Caveat

Measured on L1 only. The law effect is a modelling difference and is expected to vary
little with the mesh, but this is not verified.

## Raw data

Not copied here: `sonic:/scratch/pcardiff/tube_stvk/{L1hooke,L1stvk}/`
(`log.Allverify`, `log.environment`, `postProcessing/` with `mesh_robin_backward.csv`
and the histories); driver script `sonic:/scratch/pcardiff/tube_stvk/run_driver.sh`.
