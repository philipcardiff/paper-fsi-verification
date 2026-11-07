# REJECTED RESEARCH PROTOTYPES - not approved, not production-ready FEBio patches

These diffs are the source changes of two **experimental pressure-projection stabilisations** (Dohrmann-Bochev polynomial pressure projection added to the continuity row of FEBio v4.13's fluid dilatation formulation). Both **failed their pre-registered acceptance criteria** (see `../docs/FAILED_STABILISATION_HISTORY.md`). They are kept so that the failed experiments can be inspected and reproduced, not so that they can be applied.

- They are **not proposed for FEBio**, are not a correction of anything, were never merged into or installed over the stock source or binary, and carry no claim that they improve the benchmark result.
- No compiled binaries and no copied source trees are included. Each diff is a `diff -ruN` against a pristine export of the tag `v4.13`; the round-1 and round-2 diffs were checked to apply cleanly to such an export (`patch -p2`, on a scratch copy; the maintainers' source tree was not modified).
- The stabilisation is switched on only by the environment variable `FEBIO_PPP_ALPHA` (0 = the stock formulation, bitwise at element level).

| file | what | lines | sha256 |
|---|---|---|---|
| `round1_alpha075/ppp_v413_round1_alpha075.diff` | round 1: new `FEFluidPPP.h/.cpp`, hooks in `FEFluidDomain3D`, `FEFluidFSIDomain3D`, rim pressure re-evaluated once per step when alpha>0 (`FEFluidNormalVelocity`), version string. Predeclared alpha = 3/4. (`ppp_v413_experimental.diff` of the working tree is byte-identical.) | 266 | `d43859c1fee504c7a430fb2326fd5e8cc6cab2143076dc31a38c2a92f780f68a` |
| `round2_flux_preserving/ppp_v413_round2_flux_preserving.diff` | round 2 (cumulative against v4.13): constrained-node redistribution of the stabilisation rows (flux-preserving form), optional `RimOnce()`; predeclared alpha = 1e-3 | 297 | `71bffd9b2290ef67871259f60b30b8916f8884ea24f8dee0e11faebcfb123c98` |
| `superseded_unapplied_sketches/experimental_BP_stabilisation.diff` | early sketch (Brezzi-Pitkaranta type), never compiled or run | 61 | `5e83573ab5f74a0fc0f39b9306a606c410adce8669275704c9b0b4644f13e653` |
| `superseded_unapplied_sketches/experimental_PPP_stabilisation.diff` | early sketch of the projection term (FSI domain only; its tangent was judged inconsistent by the derivation report), never compiled or run | 78 | `ade62475e375b4b7fb0c25fba49312a8ca2ab45145028b372a7cba431828dace` |

Binaries that were built from them (not shipped): `febio4` 85c344a7308ecf9e1bd27a287277ae85f868f8db39db43ab7cc2c3b2649ca8a8; `libfebiofluid.so` 1bff35b8...0dca27 (round 1) and 6627dadd...3118bb (round 2) - full hashes in `../provenance/BUILD_PROVENANCE.md`.

Outcome in one line each: round 1 suppressed the row-parity mode 4500-fold but broke mass conservation (Q_out/Q_in - 1 = -1.55 %) and shifted the smooth pressure by 3-7 %; round 2 conserved the global row balance to machine precision but failed three of its six pre-registered criteria on the exact Poiseuille channel (nodal flux 4.9e-4, dp/dx -0.29 %, y-alternation only about 2x lower near the outlet). Neither was run on a coupled steady problem beyond 5 steps.
