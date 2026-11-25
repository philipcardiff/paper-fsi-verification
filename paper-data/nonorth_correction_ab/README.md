# Boundary non-orthogonal correction: A/B on Turek–Hron FSI3 1x′

On OpenFOAM.com (v2512), the moving-wall and interface velocity conditions build their non-orthogonal correction vector from `patch().delta()`, which is purely normal for non-coupled patches. The correction is therefore zero (solids4foam issue #748). The fix is in solids4foam PR #760, commit `f5e355c56`.

## Runs

Both runs use the corrected FSI3 set-up at level 1x′. The case generator is `paper-data/hronturek_production_scripts/ht_coupled.py`, with `--level 1 --solid-level 2 --delta-t 5e-4 --end 6.2 --cores 8`. Both ran on the Mac Studio, with OpenFOAM v2512 and Apple Clang.

| Run | Build | Case |
|---|---|---|
| base | `206437983` (post-#671) | `/Users/philipc/claude-runs/fsi3ho/cases/c1x_std_np8_t6.2` |
| corrected | `f5e355c56` (PR #760) | `/Users/philipc/claude-runs/fsi3ho/cases/c1x_std_np8_t6.2_pr760` |

The set-ups are identical: `controlDict`, `fsiProperties` and `fvSchemes` were diffed.

## Results

The analysis is `ht_coupled_analysis.py --t-end 6`: five u_y periods before t = 6 s, with the benchmark definitions. Displacements are in mm, forces in N/m and the frequency in Hz.

| Quantity | base | corrected | change |
|---|---:|---:|---:|
| u_x time mean | −2.7238 | −2.7146 | +0.34 % |
| u_x (max+min)/2 | −2.7521 | −2.7419 | +0.37 % |
| u_x amplitude | 2.6506 | 2.6413 | −0.35 % |
| u_y amplitude | 33.9984 | 33.9508 | −0.14 % |
| frequency | 5.5777 | 5.5812 | +0.06 % |
| drag mean | 457.3173 | 457.3512 | +0.01 % |
| drag amplitude | 27.4289 | 27.4789 | +0.18 % |
| lift amplitude | 191.3100 | 191.2455 | −0.03 % |

Coupling iterations per step: 11.61 for base, 11.55 for corrected.

## Conclusion

The effect is at most 0.37 %, well below the 1x′→2x′ mesh change. The paper's results use the uncorrected conditions and say so (Section 3, Appendix "Lessons").
