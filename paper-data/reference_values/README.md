# Reference values (Appendix `app:reference_values`)

`reference_values.csv` holds, in machine-readable form, the values of the
"Reference values" appendix of the paper (`appendixReferenceValues.tex`,
Tables `tab:refvals_*`), so that another code can be compared with them
without the solids4foam case files.

Each row is one quantity of interest (QoI) of one case. Columns:

| Column | Content |
|---|---|
| `case` | benchmark (tutorial name where one exists) |
| `quantity` | the QoI |
| `definition` | how the QoI is evaluated: point, time window, statistic |
| `units` | units of `finest_value`, `extrapolated_value` and `reference_value` (a scale such as `1e-5 m` multiplies the value) |
| `finest_level` | our finest mesh level and time step (levels are named L1, L2, ... from the coarsest; `n` for uniform square meshes; 1x-8x for the uniform Turek-Hron family) |
| `finest_value` | our value on that level |
| `extrapolated_value` | Richardson value, where the rules of Section 4.2 of the paper allow one |
| `uncertainty`, `uncertainty_type` | GCI (labelled conditional or provisional where so in the paper), finest-level error against an exact reference, a last change, or a difference from the reference |
| `reference_value`, `reference_source`, `reference_category` | the reference and its category (A exact, B independently converged, C partly quantified, D unquantified; Table 1 of the paper) |
| `status` | `final`, or `pending` where the result is being recomputed |
| `paper_table`, `source` | the appendix table and the table, section or paper-data file the value is taken from |

Every value is copied from the paper's tables and text, except the collapsible-channel
extrema, which are read from the stored wall-midpoint histories in
`../collapsible_mesh_study/` (`F4.dat`, the L4 run, and `oomph_extrapolated_dt160.csv`)
as the vertex of the parabola through the three stored values nearest each extremum,
the method of the tutorial driver (it reproduces the driver's first trough,
`troughValue` in `mesh_study.csv`, to all digits).

Pending (marked `pending` here and `[pending]` in the paper): the computed ring values
(IQN-ILS rerun)
and Turek-Hron FSI2 (rerun on the Gmsh family). The SU2-CalculiX cavity limits are
quoted as in the paper, pending the merge of the SU2 data.

The partitioned manufactured-solution errors (Aitken, Robin-Neumann) are the final
results of `../mms_fsi/results/partitioned_final/` (Table `tab:mms_fsi`).

The exact fields of the manufactured solution are not tabulated: they are evaluated
by `../mms_fsi/derive/generated/mmsFsiExact.H` (Appendix C of the paper).
