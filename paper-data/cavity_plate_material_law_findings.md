# cavityFlexibleBottom: plate material law and the reference differences — findings

Data and scripts: `paper-data/cavity_plate_material_law/`
(`plate_laws.py`, `results.csv`; `reference_convergence.py`,
`reference_levels.csv`, `reference_richardson.csv`). Status: plate-only check
complete (2026-10-20); COMSOL results and the Tuković force correction added
(2026-10-26, see "Four references"); COMSOL to 5 decimal places and the
Ansys and COMSOL material models confirmed by E. Tandis on PR #12
(2 November 2026): both used linear elastic material with geometric
nonlinearity, Ansys as a 3-D slab with symmetry planes front and back,
COMSOL in 2-D.

**Update, 27 November 2026 (Tuković law, open point 1).** The released
foam-extend-4.0 FSI toolkit of Tuković et al. (Fsi_40, June 2016) was rebuilt
and run on the cavity: see `cavity_tukovic_fsi40/README.md`. Its finite-strain
"linearElastic" plate is exactly StVK (no neo-Hookean or hypoelastic law), and
on L1–L3 it converges to u_y −0.24061 m and F_y −5.2594 N/m, within 0.01 % and
0.03 % of our limits. The published displacement is 0.8 %, 3.5 % and 4.2 %
above it on L1–L3; small strain (+22 %), plane stress (+12 %) and the
`updatedLagrangian` keyword (−27 %) do not reproduce it either. The source of
the published 4.3 % remains unidentified.

**Update, 20 November 2026.** The solids4foam rows of `reference_levels.csv`
and `reference_richardson.csv` are now the partitioned sequence with a
skew-corrected continuity flux (`paper-data/cavity_partitioned_skew/`). The
tables below keep the earlier sequence. With the new sequence:

- **Limits (levels 2–4):** u_y −0.240583 m (order 1.98) and F_y −5.257947 N/m (order 1.86).
- **Forces:** Tuković agrees to 0.02% (was 0.07%). ANSYS is 0.27% and COMSOL 0.71% above (relative to solids4foam).
- **COMSOL displacement:** COMSOL is 1.0% above, and 0.55% above at its load.
- **Displacement ratios:** ANSYS 1.054 and Tuković 1.043, unchanged.

## The difference to explain

Steady state, point A on the wetted surface of the plate, four levels:

| Δx (m) | Ansys u_y (m) | solids4foam u_y (m) | Ansys F_y ×100 (N/m) | solids4foam F_y (N/m) |
|---|---:|---:|---:|---:|
| 0.1 | −0.24727 | −0.20457 | −5.144 | −5.1528 |
| 0.05 | −0.25054 | −0.23145 | −5.218 | −5.2388 |
| 0.025 | −0.25197 | −0.23840 | −5.250 | −5.2566 |
| 0.0125 | −0.25274 | −0.24011 | −5.263 | −5.2599 |

- **Sources:**
  - Ansys: PR #12 (E. Tandis), `AnsysResults/cavityFlexibleBottom-Ansys-Dx*.csv` at `fe69944`, last row (t = 1000 s), with F_y × 100 for the 10 mm slab depth.
  - solids4foam u_y and F_y: the partitioned strict-test sequence (solids4foam `f356c3242`, `cavityFlexibleBottom/verification/results/four_level_partitioned_monolithic.csv`, recorded at `5f6f7a6da`). Until 4 November 2026 the solids4foam column here mixed the pre-#597 quasi-monolithic u_y with the #597 quasi-monolithic F_y.
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
python plate_laws.py 20 2 1.3245 stvk   # load sensitivity, +0.72% (COMSOL force)
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

## Four references: Richardson extrapolation (2026-10-26)

`python reference_convergence.py` writes `reference_levels.csv` and
`reference_richardson.csv`. All forces are per unit depth.

- **Tuković et al. (2018) forces:** the published F_y values are per 0.05 m depth, not per unit depth (confirmed by Ž. Tuković, October 2026), so they are divided by 0.05. The Fig. 12(a) values are read from the vector paths of the published PDF (page 16), not hand digitised; they reproduce the paper's stated orders (1.85 for u_y, 1.91 for F_y).
- **COMSOL:** PR #12, `ComsolResults/cavityFlexibleBottom_comsol_steadyStateResults.csv` at `4f8a582` (steady solver), to 5 decimal places (updated 2 November 2026; the earlier `551e91e` gave 3 significant figures). Solid model: linear elastic material with geometric nonlinearity, i.e. St Venant–Kirchhoff, in 2-D (E. Tandis, PR #12). Element types and the coupling are not yet stated.

| Δx (m) | solids4foam u_y | Tuković u_y | Ansys u_y | COMSOL u_y | solids4foam F_y | Tuković F_y | Ansys F_y | COMSOL F_y |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.1 | −0.20457 | −0.23324 | −0.24727 | −0.23632 | −5.1528 | −5.1603 | −5.144 | −5.1088 |
| 0.05 | −0.23145 | −0.24571 | −0.25054 | −0.24076 | −5.2388 | −5.2290 | −5.218 | −5.2276 |
| 0.025 | −0.23840 | −0.24955 | −0.25197 | −0.24315 | −5.2566 | −5.2494 | −5.250 | −5.2956 |
| 0.0125 | −0.24011 | −0.25062 | −0.25274 | −0.24306 | −5.2599 | −5.2548 | −5.263 | −5.2952 |

Richardson extrapolation from the three finest levels (r = 2):

| Source | u_y extrapolated (m) | Order | F_y extrapolated (N/m) | Order |
|---|---:|---:|---:|---:|
| solids4foam | −0.24067 | 2.03 | −5.2607 | 2.43 |
| Tuković et al. (2018) | −0.25103 | 1.84 | −5.2567 | 1.92 |
| Ansys | −0.25364 | 0.89 | −5.2719 | 1.30 |
| COMSOL | −0.2431 (plateau, two finest within 0.04%) | – | −5.2954 (plateau) | – |

- **COMSOL extrapolation is not possible on the finest three levels**, now with 5 decimal places as well. Both u_y and F_y converge and then reverse slightly on the finest level: u_y changes by −4.44e-3 and −2.39e-3 m, then by +9e-5 m (0.04%); F_y by −0.119 and −0.068, then by +0.0004 N/m (0.008%). The three coarse levels give orders of 0.89 (u_y) and 0.80 (F_y), so the coarse-level extrapolations (−0.2459 m, −5.387 N/m) are not credible. The two finest levels agree to 0.04%, so the plateau −0.2431 m (F_y −5.2954 N/m) is used, with an uncertainty of about 1e-4 m from their spread. A floor of this kind is typical of a steady solver tolerance or of a non-nested mesh sequence; COMSOL's element types and tolerances are still to be confirmed.
- **Tuković force = solids4foam force.** With the depth correction, the extrapolated forces agree to 0.07% (−5.2567 vs −5.2607 N/m). Tuković and solids4foam apply the same fluid load, so the 4.3% displacement difference (−0.2510 vs −0.2407 m) is entirely a plate response, consistent with a neo-Hookean plate (plate-only ratio 1.047–1.052).
- **Ansys force** is 0.2% above solids4foam (0.3% above Tuković); **COMSOL force** is 0.7% above. The cause of the COMSOL offset (integration surface, viscous traction, depth convention) is asked on PR #12.
- **COMSOL displacement after allowing for its load.** In the plate-only StVK model, raising the load by the same 0.72% (q = 1.315 → 1.3245 Pa, 20 × 2) raises u_y by 0.50% (−0.23887 → −0.24007 m), less than proportional because the plate stiffens. A StVK solution under COMSOL's load (0.66% above solids4foam's) would therefore give about −0.2407 × 1.0046 ≈ −0.2418 m. COMSOL's −0.2431 m is 1.0% above solids4foam's extrapolated −0.24067 m and 0.55% above the load-corrected −0.2418 m, against 4–5% for a neo-Hookean or Hencky plate. With the solid model now confirmed as StVK, COMSOL converges to the StVK solution of solids4foam to within about 0.5% after allowing for its load, and within 1.0% without that allowance.

## Status and open points

1. **Which law did Tuković et al. use?** Their published displacement matches a neo-Hookean plate, although the case is stated as StVK. Philip is checking this with Ž. Tuković. (27 November 2026: the released Fsi_40 code reproduces StVK, not the published value; see the update at the top and `cavity_tukovic_fsi40/`.) Their forces (per 0.05 m depth, now corrected) match solids4foam to 0.01%.
2. **Ansys material model: answered.** Linear elastic material with geometric nonlinearity (large deflection), 3-D slab with symmetry planes. As noted above, Ansys Mechanical's large-deflection linear elastic setting is, as far as we know, a logarithmic-strain (Hencky-type) law, not StVK. Its displacement ratio to solids4foam, 1.053, matches the plate-only Hencky ratio, 1.057. The element type and formulation (e.g. SOLID186 or SOLID185 with enhanced strain) are still to be confirmed; if Ansys has a StVK option (for example a user material), a StVK Ansys run would close the comparison.
3. **COMSOL with StVK: confirmed.** The solid model is linear elastic with geometric nonlinearity (StVK), in 2-D. COMSOL's −0.2431 m is within 1.0% of solids4foam (0.55% at its 0.7% higher load). Still needed from PR #12: how F_y is integrated (the cause of the 0.7% force offset), element types and orders, and the coupling and its tolerance.
4. **Proposed redefinition if StVK is impractical in Ansys.** Use a compressible neo-Hookean plate with μ = 178.57 Pa and K = 833.33 Pa (Ansys Neo-Hookean with d = 0.0024 Pa⁻¹; solids4foam `neoHookeanElastic`). Ansys, COMSOL, solids4foam and, if confirmed, Tuković et al. (2018) would then solve the same problem. The case was proposed recently, by Tuković et al. (2018), so there is no long legacy to preserve. The solids4foam four-level study would then need rerunning with `neoHookeanElastic`.
5. **For the paper:** the case illustrates how an independent reference can agree in force and still solve a different constitutive problem. Matching the integrated load does not verify the structural response.
