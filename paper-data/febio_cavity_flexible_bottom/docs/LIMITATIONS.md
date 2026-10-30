# Limitations (read before using any number in this package)

1. **All three coupled runs terminated abnormally** (job exit status 1, 'Max nr of reformations reached' on the dilatation-increment criterion at its arithmetic floor). Values are the last accepted states, labelled conditional steady estimates. The cause is judged a stopping-rule artefact at steady state but it was not tested by changing the tolerance.
2. **M3 is not steady to the stated criterion** (relative u_y spread 1.05e-3 vs 1e-4); about +/-2e-4 m and +/-1e-3 relative. Observed orders and extrapolations use three levels, M1 includes hex8 plate stiffening (4.5 %), and the F_y order is not reliable. They are indicative only.
3. **No mesh level finer than M3**, no M2/M3 sensitivity runs, no variation of rhoi, ramp length, E_mesh, pin value or the dilatation tolerance.
4. **The strict fluid-only gate was not passed** (C1: drift; C1b: upstream gradient criterion failed as written, criterion defect diagnosed from the same data). The coupled runs proceeded as a documented diagnostic.
5. **No independent StVK reference exists**; solids4foam partitioned and quasi-monolithic share a code lineage. Tukovic/ANSYS/COMSOL displacements follow other (or unconfirmed) plate laws. The solids4foam pressure and traction fields were never available, so no local comparison was possible.
6. **The cause of the 5 % force gap is not established.** The row-parity mode is real and contaminates the applied wall load (F1); that it explains the gap is not shown (F2 matches at M3 only and is FEBio's own averaging, not an independent estimator).
7. **Viscous contribution to F1 not reconstructed independently**; a 2-point estimate gives <=0.0015 N/m (<=0.03 %).
8. **FEBio reaction forces are identically zero**, so no plate force balance exists; S2 tests only the pressure part of the transfer and cannot detect an error in the pressure.
9. **The four corner pins are an added numerical constraint** (required, bound 0.1-0.4 % of F_y); the inlet rim pressure is lagged.
10. **Binary provenance** is strongly supported, not proven bit-exactly (BUILD_PROVENANCE.md). Runs were made on a heavily loaded shared node; wall times are not performance data.
11. **Stabilisation experiments** are rejected prototypes: one alpha per round, fixed-plate / tiny-channel tests only, the FSI-domain term exercised for 5 steps; scaling statements are derivations, not tests.
12. **Reproducibility**: re-tabulating from stored outputs is bookkeeping, not verification; 0.51 GB of HPC outputs are not shipped (hashes provided).
13. No claim is made about the benchmark specification, the manuscript, solids4foam, or other FEBio versions.
