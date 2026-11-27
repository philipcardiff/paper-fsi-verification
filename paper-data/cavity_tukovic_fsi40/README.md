# cavityFlexibleBottom with the released code of Tuković et al. (2018) (November 2026)

The flexible-bottom cavity was rebuilt and run with the fluid–solid interaction
toolkit for foam-extend-4.0 released by the authors of Tuković et al. (2018),
to find out which plate law produced the published displacement. The toolkit is
the same lineage as solids4foam (solids4foam descends from it), so this is a
same-lineage check, not an independent one.

## Code

- **Source:** `https://openfoamwiki.net/images/d/d6/Fsi_40.tar.gz` (tarball
  dated 24 June 2016; the newest FSI tarball on the openfoamwiki
  fluid–structure interaction toolkit page; the others are `Fsi_31` and an
  unrelated turbine tutorial). It is not
  redistributed here; download it from that URL.
- **Contents used:** `fsiFoam` and `solidFoam` with the solid solver
  `unsTotalLagrangianSolid` and the `linearElastic` rheology law. The tarball
  has no cavity tutorial; the case was rebuilt from the solids4foam tutorial.
- **Plate laws available:** `linearElastic` with `nonLinear off` is small
  strain; with `nonLinear totalLagrangian` it is total Lagrangian with the
  Green strain and S = λ tr(E) I + 2μE, i.e. St Venant–Kirchhoff
  (`unsTotalLagrangianSolid.C`, lines 1260–1275). Other `nonLinear` keywords
  (e.g. `updatedLagrangian`) fall through to the small-strain equation in the
  solver. There is no hyperelastic (neo-Hookean) and no hypoelastic law.
  `planeStress yes` replaces λ by νE/((1+ν)(1−ν)).

## Build (foam-extend-4.0, Gcc7, `linux64Gcc7DPOpt`)

Built against an existing foam-extend-4.0 installation into an isolated
`$HOME/tukovic-fsi40/platform/{lib,bin}` (`scripts/env-fe40.sh`), then
`./Allwmake` in the extracted `FluidSolidInteraction/`. Three porting fixes
were needed; none touches the solid or coupling physics:

1. Removed `-Werror` from three `Make/options` files (Eigen warnings under
   gcc 7 broke the build).
2. `RBFMeshMotionSolver.C`: added `#include <vector>` and replaced
   `shared_ptr(false)` by `nullptr` (RBF mesh motion is not used here).
3. `icoFluid.C` constructor: added
   `mesh.schemesDict().setFluxRequired(p_.name());`. The installed
   foam-extend-4.0 revision no longer reads `fluxRequired` from `fvSchemes`,
   and without it the first pressure-flux call aborts.

## Case

`scripts/cavity/mkcase.py OUTDIR R DT ENDTIME` writes an `fsiFoam` case:

- Fluid mesh: the solids4foam tutorial `blockMeshDict`
  (`scripts/cavity/s4f_fluid_blockMeshDict`) with all cell counts × R
  (R = 1, 2, 4 for L1–L3, Δx = 0.1, 0.05, 0.025 m). Solid mesh 40R × 4R, one
  cell deep (`empty`).
- Fluid: ρ = 1, ν = 0.01, `consistentIcoFluid`, Euler; parabolic inlet
  (peak 1.5 m/s, 10 s ramp), outlet p = 0, no-slip walls with
  `extrapolatedPressure`, `movingWallVelocity` on the interface; quadratic
  inverse-distance Laplacian mesh motion; skew-corrected schemes as in the
  toolkit's HronTurekFsi3 tutorial.
- Solid: `unsTotalLagrangianSolid`, `nonLinear totalLagrangian`,
  `linearElastic` E = 500 Pa, ν = 0.4, ρ = 1000 kg/m³, `planeStress no`;
  clamped ends, traction-free underside; damping 0.1 s⁻¹.
- Coupling: Aitken (initial factor 0.1), tolerance 1e-5, GGI; at most 8
  coupling iterations per step on every level.
- Pseudo-time Δt = 0.5, 0.25, 0.125 s (L1–L3) to 250 s; u_y steady to
  < 1e-6 m over the last 50 s on every level.
- QoIs: `scripts/cavity/extract.py CASE R` (u_y at the top mid-span vertex
  from `pointD`; F_y from the last-iteration "Total force (fluid)", per 1 m
  depth; its printed sign is opposite to ours).

Runs were serial on xenosim (4–45 min each).

## Results

| Level (Δx, m) | Fsi_40 StVK u_y (m) | Fsi_40 F_y (N/m) | Tuković (2018) u_y (m) | Tuković F_y (N/m) | solids4foam u_y (m) | solids4foam F_y (N/m) |
|---|---:|---:|---:|---:|---:|---:|
| L1 (0.1) | −0.231358 | −5.16296 | −0.23324 | −5.1603 | −0.20457 | −5.1528 |
| L2 (0.05) | −0.237310 | −5.23740 | −0.24571 | −5.2290 | −0.231268 | −5.233334 |
| L3 (0.025) | −0.239432 | −5.25437 | −0.24955 | −5.2494 | −0.238225 | −5.251155 |
| extrapolated | **−0.24061** (p = 1.49) | **−5.2594** (p = 2.13) | −0.25103 (p = 1.84) | −5.2567 | −0.240583 | −5.257947 |

- Fsi_40 limits are three-level Richardson estimates (L1–L3; L4 not run). They
  lie 0.01 % (u_y) and 0.03 % (F_y) from the solids4foam limits.
- Tuković / Fsi_40 u_y ratio: 1.008 (L1), 1.035 (L2), 1.042 (L3): the gap
  grows under refinement.
- Tuković and our L1 values: `../cavity_plate_material_law/reference_levels.csv`;
  our L2–L3 and limits: `../cavity_partitioned_skew/`.

Other plate options of the same code, L2 (F_y within 0.1 % of the StVK run):

| Option | u_y (m) | ratio to StVK |
|---|---:|---:|
| `planeStress yes` | −0.265107 | 1.117 |
| `nonLinear off` (small strain) | −0.290320 | 1.223 |
| `nonLinear updatedLagrangian` (small strain in the solver, finite-strain terms in the traction condition) | −0.174309 | 0.735 |

The incremental total Lagrangian solver (`unsIncrTotalLagrangianSolid`) is
StVK by construction; its L2 run stopped without an error message at
t ≈ 52 s. None of the options reproduces the published +4.3 %.

## Plate-only check (solidFoam, uniform follower pressure 1.315 Pa)

`scripts/plate/mkplate.sh DIR NX NY NONLINEAR PLANESTRESS`, then
`solidFoam`, then `scripts/plate/extract_A.py DIR NX NY`. Same problem as
`../cavity_plate_material_law/plate_laws.py`.

| Mesh | Fsi_40 StVK | Fsi_40 small strain | FE StVK (Q9) |
|---|---:|---:|---:|
| 40 × 4 | −0.235187 | −0.297011 | −0.24047 |
| 80 × 8 | −0.238765 | −0.305254 | −0.24101 |
| 160 × 16 | −0.240350 | −0.308842 | – |
| 320 × 32 | −0.240980 | – | – |
| extrapolated | −0.24140 (p = 1.33) | −0.3116 (p = 1.20) | ≈ −0.2412 |

The code's TL solid converges to the independent FE StVK plate to about 0.1 %.
Plane stress (160 × 16): −0.26743 m, 1.113 × StVK (FE check
`scripts/plate/planestress_check.py`: 1.117). FE reference outputs:
`scripts/plate/fe_*.txt`. The small-strain ratio (1.29 here, 1.21 in
`plate_laws.py`) depends on whether the pressure acts on the reference or the
deformed edge.

## Conclusion

The released code, contemporary with the publication, reproduces the
St Venant–Kirchhoff solution of solids4foam, SU2–CalculiX and FEniCSx, not the
published displacement, and no plate formulation available in it does. The
source of the published 4.3 % is not identified; it is not in this snapshot
(possibly a later version of the code or a different case setting).
