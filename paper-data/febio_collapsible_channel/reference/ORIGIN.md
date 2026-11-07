# Origin of the reference data (oomph-lib solutions)

These seven files are unmodified copies of `tutorials/fluidSolidInteraction/collapsibleChannel/verification/reference/` of the solids4foam repository (`development` branch; the audit records commit `b200718f4f8390c90f9c7e3ded8b74f3e13941cc` of 2026-10-30; the working copies available when this package was prepared are other solids4foam checkouts whose copies of these files have identical sha256, so the files themselves are the same). The values were computed by the solids4foam maintainers with the oomph-lib driver `reference/oomph-lib/collapsible_channel_ref.cc` (clamped Kirchhoff-Love beam with the regularising wall density, Crouzeix-Raviart fluid, BDF2 + Newmark; "refinement 1" = 100 x 16 fluid and 40 beam elements); the oomph-lib driver was not re-run for this package. The converged reference is the Richardson extrapolation in time of dt = 0.003125 and 0.0015625 s. The static beam deflection (-0.14362832 m midpoint, -0.09645925 m at the quarter points) is in the JSON file (computed by the maintainers with `--q 0 --steady 20`). Stated precision (maintainers' README): about 0.3 % of the 0.216 m peak.

| file | bytes | sha256 |
|---|---|---|
| collapsibleChannel_oomph_reference.csv | 49,555 | aee4ad507b55d14e388843206254e593ab45bb0bc8a4f183bffa0020c84c62c6 |
| collapsibleChannel_oomph_massless_reference.csv | 49,555 | 16b6cdbbd2f6b8299e0568f576123fa3407e3a347a6da871aa5045790a30fd7f |
| collapsibleChannel_oomph_dt40.csv | 2,664 | 3ad35fca25db17b08fff45d1ad9292918c5606608d4bd042daf37ab6d6bf6677 |
| collapsibleChannel_oomph_dt80.csv | 5,237 | 4e60ea020d8ae23269b73b40cb1b0ae6832b0f61a8fc810b792558237fbcf968 |
| collapsibleChannel_oomph_dt160.csv | 10,508 | 07eb49cd1a8432dcb8ae4a0d165d1131bcaace0c5bb666cbb205324bf867eb38 |
| collapsibleChannel_oomph_dt320.csv | 21,286 | e091ed7c1ac0d5cb973289089f4438acdf6cdb90662739f1d6e7e6592796cf1f |
| collapsibleChannel_verification_references.json | 2,234 | 31e82f01a05538089dea06c77afbec1c512d6dfa46fa172c7dfe93555cc677fb |

Columns of the CSVs: time, wallQuarter, wallMid, wallThreeQuarter (vertical displacement in m at x = 7.5, 10, 12.5). The `massless` file is not used in this package. All numbers quoted from these files (first trough -0.2160 m at 0.4685 s converged, -0.2129 m at 0.4735 s at dt 0.025; rebound maximum -0.0404 m at 0.7905 s; second trough -0.2037 m at 1.120 s) are recomputed in `results/oomph_reference_features.csv` by a 5-point parabolic fit around each extremum.
