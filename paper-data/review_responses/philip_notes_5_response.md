# Response to Philip's notes 5 (final review: Sections 5.2.4 to 8)

Branch `text/notes5`. There is one entry per comment in `philip_notes_5.md`. Unless stated otherwise, the edits are in `sectionBenchmarkSuite.tex`. `latexmk -pdf article.tex` builds with no errors and no undefined or multiply defined references (104 pages).

## 5.2.4 Flexible-bottom cavity

**Outlet condition.** Done. The outlet now has `p = 0` and a zero normal gradient of the velocity.

## 5.3 Partly quantified references (category C)

**Move FSI2 into 5.3.** Done.
- The new Section 5.3.1 "Turek–Hron FSI2" (`sec:fsi2`) has the usual headings. It refers back to 5.2.2 for the geometry, conditions, ramp and QoI definitions. It gives the FSI2 parameters, the category-C reasoning (Featflow level 3→4 changes of 0.9–2.5 %, no finer or independent sequence), the numerical model, Table `tab:hronturek_fsi2` and its results paragraph, and the lessons.
- All the `\hl{[pending: FSI2 on the Gmsh family]}` markers are kept.
- 5.2.2 is retitled "Turek–Hron family: definition, FSI1 and FSI3". Its introduction says that FSI2 appears in 5.3 because of its category, and that the definitions in 5.2.2 apply to all three problems. The FSI2 parts of its Reference, Numerical model and Lessons now point to 5.3.1.
- The 5.3 preamble now introduces two problems.

**"External flow".** Replaced. The beam is now "the only steady three-dimensional case in the suite in which the flow passes around a solid obstacle, here a plate clamped in a duct".

**"The distributed case" and the viscous term.**
- The text now reads "a coarse mesh of the plate in the duct", with no reference to the tutorial.
- The sentence now states the point directly: the effect of omitting the viscous term is small and decreases with refinement (+1.7 % on the coarse mesh, +0.35 % on that mesh refined twice), consistent with an O(h) boundary error.
- It then compares this with the beam's uncertainty: on L3/L4 the effect is expected to be a few tenths of a per cent or less. That is below the stated 1 % accuracy of the reference and below the 0.5–0.8 % by which our conditional Richardson values exceed it. The text also says the effect was not computed on L1–L4.
- I could not find in paper-data how that mesh relates to L1–L4, so I added `\hl{[check: relation of these two meshes to L1–L4]}`.
- The matching sentence in Section 3 (`sectionNumericalModel.tex`) is changed in the same way.

**"The comparison was made after the sequence was characterised".** Explained. The text now says that the comparison with the reference values was made only after the convergence of the sequence had been assessed, and that no setting was adjusted towards them.

**What the Kratos result confirms.** Added: "This confirms that the reference values correspond to an inflow peak of 0.2 m/s, not the printed 0.3 m/s."

**Why C, and how it could become B.** Added. Neither sequence is shown to be asymptotic, and the reference is published to three digits without its sequence. The missing evidence is a finer level: L5, about 3.3×10⁷ fluid cells, estimated at 340–930 node-hours (Appendix A). Faster implementations (GPUs, parallel-in-time integration) are likely to make it affordable, and with it, or with an independent converged sequence, the case could become B.

## 5.4 Unquantified references (category D)

**Validation, not verification.** Added to the D preamble: "…a separate comparison with measurements, which, as noted, is evidence of validation, not of verification."

**Wall-front speed.** Defined at first use, in the Problem paragraph. It is the speed at which the radial bulge of the wall travels along the tube. At each of a few inner-wall stations, the arrival time is the first time at which u_r reaches half its first peak. c_wall comes from a least-squares fit of these arrival times against axial position over z = 1–4 cm. The pressure-front speed is defined in parallel.

**StVK vs Hookean wall.** As agreed, there is no change yet. The paper keeps the Hookean wall until StVK results (ours and Kratos) are available. The existing L1 StVK-vs-Hookean sentence stays.

**Quarter-domain / axisymmetry.** Added to the Numerical model: the problem is axisymmetric and an axisymmetric model would suffice. The three-dimensional model is used deliberately, to exercise the three-dimensional curved non-conformal interface.

**"Switch-off" and "meshes".**
- **Problem paragraph.** It now explains the switch-off: the pulse has a fixed duration, 1333.2 Pa for 3 ms, and is then switched off by a linear ramp to zero over 0.1 ms.
- **Reference paragraph.** It says the references state no ramp, so their switch-off is taken as instantaneous, and that none of their meshes has a published mesh study.
- **Literature-discrepancy item.** It now names the two unchecked differences between the references' set-ups and ours:
  - *The switch-off.* Our ramp adds 1.7 % to the impulse and smooths only the trailing edge, so it is likely to change the QoIs by about a per cent or less. That is small against the 20–25 % discrepancy, but not necessarily against what remains once the time integration is matched.
  - *The references' meshes.* Their spatial errors are unknown and may account for part of the remaining differences.
- The item states that neither difference concerns our own mesh study.

**References to the solids4foam tutorials (all cases checked).**
- **3-D tube.** The digitised curves now point to `paper-data/3dtube_level3_corrected_ecore/figure/data/reference/`.
- **Ring.** The Table caption no longer cites `ringAddedMass/verification/`.
- **Collapsible channel.** "driver source is stored with the case" and "reference stored with the case" now point to `paper-data/collapsible_mesh_study/`.
- **Category-A preamble.** "Collocation scripts are stored with the cases" is removed.
- **Manufactured solution.** "The case is distributed with solids4foam [pending: release]" now reads: the case set-up, derivation and generated code are in `paper-data/mms_fsi/` (Appendix on the manufactured solution).
- **Blob in treacle.** "meshes of the tutorial" is now "baseline meshes".
- **Beam.** "distributed case" is replaced (see above).
- **Appendix D.** "tutorial filtering" / "used in the tutorials" now reads "default filtering".
- **Kept:**
  - Mentions of the oomph-lib and Kratos tutorials. These are published, cited examples of other codes, not our code tree.
  - Commented-out text and `\todo` notes.
  - The Data Availability statement, which still names the `tutorials/fluidSolidInteraction` directory of the repository. A code-availability statement needs to say where the cases are. Please say if you would rather give only the release tag and Zenodo DOI.
  - Table `tab:driver_properties` (Appendix F), which describes the solids4foam drivers.

**Digitisation paragraph.** Shortened from four sentences to three:
- how the curves were extracted (five-column average, resampled at 0.1 ms) and where they are in paper-data;
- the uncertainty: about ±3×10⁻⁶ m, about 2.5 % of the peaks, with the smoothing adding less than a third of this, as bounded on our L3 history;
- the comparison: this is much smaller than the 20–25 % difference from our time-converged results.

"The raw traces were not kept…" is removed. The existing `\hl{[check: pixel width of the two published figures]}` is kept, because the smoothing bound still depends on it.

## 5.4.2 Hessenthaler flap

**Deformed centreline.** The caption now says it is "the deformed flap centreline measured by MRI in the Phase I experiment, not computed".

**3-D perspective view.** As agreed, another agent is making the figure. Below Fig. 11 there is a placeholder, `\hl{[pending: 3-D perspective figure]}`, with a sentence saying the domain is fully three-dimensional and that the views shown are sections of it.

**Solid mesh error.** Added to the Numerical model.
- The fixed 12×8×65 solid mesh is not part of the fluid-mesh sensitivity, but its size is known.
- The shear modulus is the continuum limit of the four-nested-mesh calibration (GCI 0.7 %), so it carries no error from this mesh.
- In the earlier coupled study (`paper-data/hessenthaler_phaseI_verification.md`, on a fluid mesh of a previous family), refining the solid to 24×16×130 raised the tip by 0.10 mm (0.6 %), at an observed order of 2.2 over three nested solid meshes. The extrapolated tip is 0.12 mm (0.75 %) above the 12×8×65 value.
- So the solid error of the tip is about 0.7 %, the same size as the last fluid-mesh change. No solid-mesh study exists on the current L1–L3 fluid meshes.

**Flow-rate sentence.** Rephrased so it cannot be misread: "There the measured data are themselves inconsistent: in the experiment, the flow rates through the two planes obtained by integrating the MRI velocity measurements exceed the flow rate set at the pump by 8–10 %, an inconsistency in the measured data, not in the model."

**"Affordable".** The Lessons now give the cost of the finest mesh and point to Table `tab:cost` (Appendix F, "Provenance, cost and reproducibility"), which already has a Hessenthaler row. The recorded L3 wall times are 37 h and 54 h on 48 ranks (about two days of one node per run), not many days. The text gives "about two days of a 48-core node per run" so that it matches the recorded value.

## 5.4.3 Blob in treacle

Added: "referred to here as 'blob in treacle', a name introduced with the case in solids4foam, not used by Liu et al."

## Section 6

**FEBio cost sentence.** Removed, along with its `\todo` on per-DOF costs. It is replaced by a neutral pointer: the recorded costs of all cases are in Table `tab:cost` (Appendix F). They were measured on different machines and core counts, where the same calculation can differ by a factor of two or more, so they indicate orders of magnitude only and are not a comparison of codes. The FEBio timing (2 033 s against 224/315 s) no longer appears in the paper; it is still in `paper-data/febio_ring_added_mass/` if you want it back in the supplement.

**Backbone sentence.** Updated. The category-A problems cover small deformation (ring, Womersley) and large deformation (manufactured solution). The category-B problems extend this to thin bending walls with strong added mass (collapsible channel) and to steady large-deformation coupling of a flag behind a cylinder (FSI1). The reference uncertainties are unchanged.

**Scope paragraph.** Rephrased.
- The problems are solved here with a body-fitted ALE fluid, a Lagrangian solid and a sharp interface.
- The problem definitions and references do not depend on that approach. The suite is equally suitable for immersed-boundary and fictitious-domain methods, Eulerian and meshless formulations, and reduced solid models (beams, plates), provided the QoIs are evaluated as defined.
- It is naturally extendable to compressible or turbulent flow, non-Newtonian fluids, contact and free surfaces, although each would need specific cases of its own.

**Merge of 6.7.** Done. The subsection is removed, and its content is merged into the former 6.4, now "Defects, verification claims and reproducibility". It follows the five items of a verification claim as one paragraph:
- the motivation, in one sentence;
- the four recommended properties, generalised: the reference files also record the reference category, and the code and library versions replace "the OpenFOAM version";
- a link to item (v).

The label `sec:infrastructure` moves with the paragraph, so the references from Sections 4 and Appendix F still resolve.

## 7 Conclusions

**Heading.** "An exact reference at finite deformation qualifies this." is now "At finite deformation, an exact reference exposed errors that the linear exact problems could not."

**Long sentence.** Split into three: what the solution is; that it needs no interface source and is the first such solution for a cell-centred FV code; and the defects found. The list of four errors now starts a sentence of its own.

**"Unreleased".** Removed, here and in the orders table of Section 6.

**Future work.**
- **Axisymmetric tube.** Added to the first item: an axisymmetric model, much cheaper than the 3-D one, would allow a full convergence study and a converged reference for the 3-D calculations, with `\hl{[pending: axisymmetric tube]}`. The other agent is building it now.
- **Wider coverage.** Immersed-boundary approaches added, "which the present cases can already verify".
- **Higher order.** Added: a high-order solid was used here, but the overall procedure, with the fluid and the coupling, remains second order, so what is proposed is higher order in the fluid, the solid and the coupling together.

## Acknowledgments

16/RC/3872 removed. 21/RC/10295_P2, 21/SPP/3756 and the ERC grant are kept.
