# Figures
Regenerated from retained data by `scripts/make_figures.py` (labels and units checked by viewing each file):
- `fig_mesh_convergence_vs_references.png` - u_y(A) and F_y vs dx, all reference classes, FEBio F1 and F2 (diagnostic), odd-even bracket; the header states that all FEBio runs error-terminated.
- `fig_force_definitions.png` - F1, F2 (= F4) and the solids4foam QM force at the same stored state per level.
- `fig_row_pressure_alternation.png` - mean nodal pressure by row above the floor; floor pressure along x (M3 failure-exit state).
- `fig_numerical_sensitivity.png` - K, tolerance, dt effects at M1.
- `fig_stabilisation_outcomes.png` - round 1 / round 2 outcomes.
Original figures produced from HPC-resident outputs by `scripts/orig_cavity_figures_*` (not regenerable from the package; viewed and checked): `cavity_geometry_and_flow.png` (mesh level 1 and fixed-plate C1b velocity magnitude), `cavity_interface_pressure_and_deflection.png` (C1b raw nodal floor pressure with even/odd nodes; M1 plate deflection; read panel (b) as a conditional estimate of an error-terminated run), `cavity_poiseuille_diagnostic.png`, `cavity_cost_summary.png` (allocated core-hours, shared loaded node). The earlier originals `cavity_mesh_convergence_cross_solver.png` and `cavity_numerical_sensitivity.png` are superseded by the regenerated figures and not included.
