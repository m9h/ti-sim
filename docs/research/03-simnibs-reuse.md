# What SimNIBS already does for TI, and what this repository adds

SimNIBS 4.6.0 (python3-simnibs 4.6.0-21 in the mhough/neurofedora COPR, GPL-3) is the
volume-conductor engine here. Nothing in it is copied; it is imported.

## Reused as is

| SimNIBS piece | what it is | used by |
|---|---|---|
| `charm` head models (`m2m_<subject>`) | segmentation and meshing from T1/T2 | `ti_sim.head` (ernie) |
| `sim_struct.SESSION` + `TDCSLIST` | one FEM solve per electrode pair, electrodes placed on the skin at 10-10 names or coordinates | `ti_sim.head.TwoPairSimulation` |
| `utils.TI_utils.get_maxTI`, `get_dirTI` | Grossman's equation 1 and its closed-form maximum | cross-checked against `ti_sim.envelope` in `tests/test_envelope.py` (same numbers to 1e-12) |
| `TDCSLEADFIELD` + `TI_utils.load_leadfield/get_field` | lead fields over a cap; any pair by subtraction | `ti_sim.montage` |
| `mesh_tools.meshing.create_mesh` | CGAL tetrahedralisation of any label volume, surface tags, `relabel_internal_air` | `ti_sim.labelvol` (spheres), `ti_sim.mida` |
| custom-mesh simulations (`S.fnamehead`) with per-tag conductivities | TI on spheres and MIDA | `TwoPairSimulation(fnamehead=...)` |
| `mni2subject_coords` | MNI targets in subject space | `examples/ernie_ti.py` |
| `opt_struct.TesFlexOptimization` with `e_postproc="max_TI"` / `"dir_TI"` | differential-evolution placement of two free pairs maximising mean/focality of the envelope in an ROI (Weise et al.) | not wrapped yet; it is the right tool for free electrode placement on a charm head, and the natural next example |
| `examples/simulations/TI.py`, `examples/optimization/TI_fields_leadfields.py`, `examples/tes_flex_optimization/tes_flex_ti_*.py` | reference usage | our examples follow them |

Also from the Fedora stack: `simnibs-data` (templates), `python3-nibabel`, `python3-h5py`.

## Added here

* `ti_sim.envelope`: the formula with a brute-force direction sweep as an independent
  check, the directional envelope, modulation depth, tACS sum and the multi-pair bound.
* `ti_sim.phantom`: 2-D scikit-fem disc with wire electrodes, the Grossman geometry, peak
  and HWHM extraction; validated against the analytic Kelvin-image field.
* `ti_sim.labelvol`: synthetic layered spheres and the label-volume-to-mesh route with
  caching and extra tags.
* `ti_sim.mida`: MIDA label parsing, the 115-structure to tissue-tag mapping, crop,
  downsample, mesh.
* `ti_sim.montage`: exhaustive two-pair search with current ratios, ROI median scoring,
  focality re-scoring, target-neighbourhood candidate restriction.
* `ti_sim.membrane`: the passive low-pass / rectification arithmetic (Mirzakhalili) and
  the Esmaeilpour thresholds.
* `ti_sim.tissues`: named conductivity tables (SimNIBS, Rampersad, Grossman sphere) so a
  result records which one it used.
* Machine-checked examples and a BATS pipeline, a Fedora container.

## Not reused, and why

* **TI-Toolbox** (GPL-3): a complete application for individual MRIs on Docker Desktop with
  Electron; it pins SimNIBS 4.1/4.5 in its own image and its Python package is meant to
  run inside that image. For Fedora-native, scriptable, phantom/sphere/MIDA work it is the
  wrong shape, but its flex/exhaustive search and mTI kernel are worth reading, and its
  results on charm heads are a cross-check for `ti_sim.montage`.
* **Sim4Life**: commercial.
* SimNIBS's own TI examples are GPL-3 and are not copied; the same API calls are made.

## Licence note

ti-sim is Apache-2.0 and imports SimNIBS (GPL-3) at run time, as ae-sim does. Nothing of
SimNIBS is vendored. Distributing ti-sim as a Fedora package would make the package
depend on python3-simnibs; that is fine for a free COPR.
