# cavityFlexibleBottom: plate material law and the reference differences — findings

Data and script: `paper-data/cavity_plate_material_law/`
(`plate_laws.py`, `results.csv`). Status: plate-only check complete
(2026-10-20); the law behind the Tuković et al. (2018) values and the Ansys
material model are not yet confirmed (PR #12).

## The difference to explain

Steady state, point A on the wetted surface of the plate, four levels:

| Δx (m) | Ansys u_y (m) | solids4foam u_y (m) | Ansys F_y ×100 (N/m) | solids4foam F_y (N/m) |
|---|---:|---:|---:|---:|
| 0.1 | −0.24727 | −0.2036 | −5.144 | −5.123 |
| 0.05 | −0.25054 | −0.2314 | −5.218 | −5.228 |
| 0.025 | −0.25197 | −0.2385 | −5.250 | −5.251 |
| 0.0125 | −0.25274 | −0.2401 | −5.263 | −5.256 |

- **Sources:**
  - Ansys: PR #12 (E. Tandis), `AnsysResults/cavityFlexibleBottom-Ansys-Dx*.csv` at `fe69944`, last row (t = 1000 s), with F_y × 100 for the 10 mm slab depth.
  - solids4foam u_y: `sectionBenchmarkSuite.tex` (top surface, partitioned and quasi-monolithic).
  - solids4foam F_y: the quasi-monolithic rerun after #597 (solids4foam `f356c3242`, `cavityFlexibleBottom/verification/results/four_level_partitioned_monolithic.csv`).
- **Extrapolated u_y:**
  - solids4foam: −0.2407 m (observed orders 1.95–2.05).
  - Ansys: changes of 3.27e-3, 1.43e-3 and 7.7e-4 m between levels, i.e. about first order, giving roughly −0.2535 m.
  - Tuković et al. (2018), finest mesh: −0.2506 m.
- **Force:** agrees to 0.1–0.2% on every level, so the fluid load on the plate is the same.
- **Displacement:** Ansys is 5.3% above solids4foam and Tuković 4.1% above.
- **Horizontal displacement at A:** Ansys 4.1e-5 m on the finest level, solids4foam below 1.1e-4 m. The earlier large Ansys u_x (0.033 m) came from a mid-thickness monitor point and has been corrected.

## Hypothesis

The plate is 4 m × 0.4 m (L/h = 10), clamped at both ends, and deflects about 0.24 m. For a clamped–clamped beam under uniform load, the end curvature is about 32 w/L², so the bending strain at the clamped ends is about 32 × 0.24/16 × 0.2 ≈ 0.10.

At 10% strain, large-strain laws with the same small-strain E and ν differ in stiffness by several per cent:
- the case specifies St Venant–Kirchhoff (StVK);
- as far as we know, Ansys Mechanical has no built-in StVK;
- its "linear elastic + large deflection" setting behaves like a logarithmic-strain (Hencky-type) law.

## Plate-only check

`plate_laws.py`:
- **Geometry and supports:** plane strain; x ∈ [2, 6] m, y ∈ [−0.4, 0] m; both ends clamped over the full thickness; underside traction free.
- **Material:** E = 500 Pa, ν = 0.4, so μ = 178.57 Pa, λ = 714.29 Pa and K = 833.33 Pa.
- **Load:** a uniform follower pressure of q = 1.315 Pa on the top face, i.e. the FSI total force 5.26 N/m spread over 4 m.
- **Discretisation and solution:** total Lagrangian, 9-node quadratic quadrilaterals, 3 × 3 Gauss points, Newton with a finite-difference element tangent and 20 load steps (one for the linear law).

Run (numpy 2.5.4, scipy 1.18.1; about 25 s per law on 20 × 2 and 2 min per law on 40 × 4):

```text
python plate_laws.py 20 2
python plate_laws.py 40 4 1.315 stvk,neohookean,hencky
python plate_laws.py 20 2 1.315 neohookean,nh_s4f,nh_ansys
```

| Law | u_y at A, 20 × 2 (m) | u_y at A, 40 × 4 (m) | Ratio to StVK (40 × 4) |
|---|---:|---:|---:|
| Small strain (linear) | −0.28816 | – | – |
| St Venant–Kirchhoff | −0.23887 | −0.24047 | 1 |
| Neo-Hookean, P = μ(F − F⁻ᵀ) + λ ln J F⁻ᵀ | −0.25009 | −0.25191 | 1.048 |
| Hencky (logarithmic strain) | −0.25231 | −0.25420 | 1.057 |
| Neo-Hookean, split, U = K/4 (J² − 1 − 2 ln J) (solids4foam `neoHookeanElastic`) | −0.25119 | – | (1.052 on 20 × 2) |
| Neo-Hookean, split, U = K/2 (J − 1)² (Ansys Neo-Hookean, d = 2/K) | −0.25106 | – | (1.051 on 20 × 2) |

- **Mesh dependence:** the ratios between laws change by at most 0.1% from 20 × 2 to 40 × 4. The absolute values still change by about 0.7%; the 9-node elements are not converged in absolute terms on 2–4 elements through the thickness, but the comparison between laws is.
- **Not the FSI load:** the uniform pressure is not the real load distribution. The ratios are the evidence, not the absolute values (the 40 × 4 StVK value happens to be within 0.1% of the FSI value).
- **Neo-Hookean variants:** they differ from each other by at most 0.45%. The solids4foam and Ansys forms agree to 0.05%.
- **Linear check:** the small-strain deflection agrees with clamped–clamped beam theory with shear. Bending alone gives qL⁴/(384 E′ I) = 0.276 m, with E′ = E/(1 − ν²); shear adds about 0.01 m.

## Comparison

| Source | u_y at A (m) | Ratio to solids4foam (StVK) | Plate-only ratio of the matching law |
|---|---:|---:|---:|
| solids4foam, StVK, extrapolated | −0.2407 | 1 | StVK: 1 |
| Tuković et al. (2018), finest | −0.2506 | 1.041 | neo-Hookean: 1.047–1.052 |
| Ansys, extrapolated | ≈ −0.2535 | 1.053 | Hencky: 1.057 |

## Status and open points

1. **Which law did Tuković et al. use?** Their published displacement matches a neo-Hookean plate, although the case is stated as StVK. Philip is checking this with Ž. Tuković.
2. **Which material model did Ansys use?** Asked on PR #12.
3. **COMSOL with StVK.** COMSOL's "Linear Elastic Material" with geometric nonlinearity is StVK. Requested on PR #12; the expected value is about −0.241 m.
4. **Proposed redefinition if StVK is impractical in Ansys.** Use a compressible neo-Hookean plate with μ = 178.57 Pa and K = 833.33 Pa (Ansys Neo-Hookean with d = 0.0024 Pa⁻¹; solids4foam `neoHookeanElastic`). Ansys, COMSOL, solids4foam and, if confirmed, Tuković et al. (2018) would then solve the same problem. The case was proposed recently, by Tuković et al. (2018), so there is no long legacy to preserve. The solids4foam four-level study would then need rerunning with `neoHookeanElastic`.
5. **For the paper:** the case illustrates how an independent reference can agree in force and still solve a different constitutive problem. Matching the integrated load does not verify the structural response.
