# Replication targets

Which published TI results report numbers a simulator can be scored against, and where
this repository stands on each.

## 1. Grossman 2017, figure 2: the saline phantom (done, with an amplitude caveat)

Non-conductive cylinder 50 mm × 10 mm filled with saline, 1 mm silver wires on the rim,
1 kHz / 1.02 kHz, 1 mA per pair (2 mA sum with ratios in D-E). Five layouts, A-E, with
the electrodes on the corners of an inscribed trapezoid (bases a, b in units of R):

| layout | a | b | I1:I2 | paper peak E_AM,y (sim / meas, V/m) | paper position (x, y)/R | paper HWHM (x, y)/R |
|---|---|---|---|---|---|---|
| A | 0.39 | 1.96 | 1:1 | 2.5 / 1.1 | (0, -0.77) | (0.18, 0.37) |
| B | 1.11 | 1.96 | 1:1 | 0.9 / 0.7 | (0, -0.37) | (0.40, 0.52) |
| C | 1.96 | 1.96 | 1:1 | 0.21 / 0.23 | (0, 0) | (0.42, 0.49) |
| D | 1.96 | 1.96 | 1:2.5 | 0.20 / 0.22 | (-0.2, 0) | (0.38, 0.40) |
| E | 1.96 | 1.96 | 1:4 | 0.19 / 0.20 | (-0.35, 0) | (0.37, 0.34) |

`examples/grossman_phantom.py` (2-D scikit-fem disc, `ti_sim.phantom`; the FEM is checked
against the analytic line-source field with Kelvin images to 3 %):

| layout | ours E_AM,y (V/m) | position | HWHM | ratio to paper sim |
|---|---|---|---|---|
| A | 11.66 | (0, -0.79) | (0.17, 0.38) | 4.7 |
| B | 5.04 | (0, -0.38) | (0.40, 0.52) | 5.6 |
| C | 1.51 | (0, 0.00) | (0.42, 0.50) | 7.2 |
| D | 1.43 | (-0.23, 0.01) | (0.42, 0.39) | 7.2 |
| E | 1.34 | (-0.34, -0.01) | (0.43, 0.34) | 7.1 |

and E_AM,x: C at (0, -0.59) R, HWHM (0.37, 0.38) vs the paper's (0, -0.59), (0.37, 0.42);
steering D/E to x = -0.23 / -0.33 R vs -0.2 / -0.35 R. **Every position and half-width
matches**; the envelope *geometry* of the paper is reproduced. **Amplitudes are a constant
factor ~7 above the paper's simulated values** in the far-from-electrode layouts (and
less, 4.7-5.6, next to the wires where their electrodes were 0.16 R patches rather than
1 mm wires). A uniform factor with identical geometry means a different σ·h·I product:
the paper states 0.333 S/m, 10 mm, 1 mA, and we use exactly those. Candidates: the
reported values being envelope *amplitude* (half of peak-to-trough) at 1 mA *peak*
rather than our peak-to-trough at 1 mA amplitude (×2 and ×2), or the physical saline
(0.9 % NaCl is ~1.5 S/m, 4.5× 0.333) having been used in the simulation too. The
three together would over-explain it; none is stated. The pipeline asserts the geometry
and records the ratio (5-9) rather than hiding it.

## 2. Grossman 2017, figure S2: the 4-layer sphere (built, cross-checked)

`examples/layered_sphere.py`: the sphere with the paper's proportions, scaled to R =
85 mm (scalp 85.000 mm 0.333 S/m, skull 80.750 0.0083, CSF 73.525 1.79, brain 71.570
0.333), built as a 1 mm label volume, meshed with SimNIBS `create_mesh` and solved with
SimNIBS, two pairs on the equator in layout C (pair 1 at 180 ± 11.5°, pair 2 at ∓11.5°,
10 mm electrodes) and the 1:1, 1:2.5, 1:4 ratios. The paper gives no table for the sphere;
what the pipeline checks is qualitative and theirs.

**Cross-check (2026-10-26)** against an independent DOLFINx solve (CG + BoomerAMG, gmsh CAD
mesh, equipotential scalp patches; `ci/conductivity-admm/sphere_ti_forward.py` in the
neurofedora repo) and the analytic n = 1 Legendre solution of the 4-shell sphere with
point electrodes, 1 mA per pair at 1:1:

| quantity | analytic (point electrodes) | DOLFINx, 2 mm 1 S/m pads, 634k tets | SimNIBS on the same gmsh mesh (262k tets), corrected current | SimNIBS, our 1 mm voxel mesh, corrected current (49k brain tets, ~6 mm) |
|---|---|---|---|---|
| \|E\| per pair at the centre | 0.0283 V/m | 0.0259 | 0.0262 | 0.0303 (r < 10 mm, volume-weighted) |
| E_AM at the centre | 0.0565 | 0.0518 | 0.049–0.051 | 0.043 (r < 10 mm); shell medians 0.043, 0.045, 0.040 … 0.027 at 60–70 mm |
| E_AM maximum in the brain | at the centre | 0.0518, centre | 0.0565 nominal → 0.053, centre (depth 77.8 mm) | 0.052; top-1 % centroid 82 mm deep (the 2026-10-25 "55 mm" was the argmax of a plateau) |
| \|E₁\| maximum, depth | | 0.299 V/m, 15.0 mm | 0.318, 15.8 mm | 0.287, 16.3 mm (top-1 % centroid) |
| pair resistance | | 280 Ω | 279 Ω (ΔV / I_h) | 278 / 278 Ω (SimNIBS's nominal currents were 10.7 and 3.6 % off, which made the uncorrected run read 307 / 338) |

Two lessons. The "55 mm deep" maximum reported on 2026-10-25 was noise: the envelope is a
plateau over the central 30 mm (shell medians 0.044 → 0.030 V/m from the centre outward,
shell maxima 0.045–0.056 everywhere) and our voxel-derived mesh has ~6 mm brain elements,
so the element argmax lands anywhere on it. The finite 10 mm pads lower the centre field by
8 % against point electrodes because the current concentrates at the pad rim and the pair's
effective dipole shrinks (framework's analytic check: 80 % of the patch current within 2 mm
of the rim).

The second lesson is about SimNIBS and affects every head number in this repository:
**SimNIBS's current normalisation is biased high.** `fem.tdcs` solves with 0/1 V on the
electrode tops and rescales by a current estimated as the flux of the element-wise J
through the skin–electrode interface, averaged over the two electrodes. The exact current
of the discrete solution follows from the energy identity ∫σ|∇v|² dV = ΔV·I_h; on the gmsh
sphere I_h = 1.061 mA for a nominal 1 mA (SimNIBS's own "calibration error", the
disagreement between its two electrode fluxes, read 2.7 %), on ernie +4 to +13 %
depending on electrode size and layering (F3–F4 2 mm pads: +13.5 %; 4 mm gel + 2 mm rubber:
+9.3 %; the 20 mm SimNIBS example electrodes: +4.3 / +9.4 %). `TwoPairSimulation` now saves
v and J, computes I_h per pair (`head.discrete_current`) and rescales E by I_nominal / I_h;
`calibration` in every result JSON records both currents and the pair resistance. The lead
field is **not** affected: `tdcs_leadfield` prescribes each electrode's current through the
right-hand side (`assemble_rhs`, a Neumann-type formulation) and only logs the flux estimate
as a mesh-quality warning; `montage.calibrate_leadfield` re-solved all 75 ernie electrodes
with Dirichlet conditions and the energy-identity current and found SimNIBS/exact =
0.9985–1.022 (median 1.001). `TesFlexOptimization` uses node currents and is unverified.
Upstream note drafted in `docs/upstream-drafts.md`.

The voxel sphere's resistance and |E₁| now agree with the CAD mesh; its centre field is 17 %
high and the envelope there 17 % low because the brain is still meshed at ~6 mm (the
`elem_sizes` range passed to `create_mesh` did not refine the interior; open item) and the
two carriers' fields at one coarse element are neither equal nor parallel. At 1:2.5 and 1:4
(2 mA sum) the envelope hotspot moves to the weaker pair's side (x-profile peak at −0.53 R)
and the centre envelope falls 0.043 → 0.035 → 0.024 V/m, Grossman's steering.

## 3. Rampersad 2019: human head numbers (ernie stands in)

| quantity | paper | ours (ernie low-res, SimNIBS conductivities, isotropic; 2026-10-26, fields normalised to the exact discrete current, see §2) |
|---|---|---|
| standard montage, 1 mA per pair: hippocampus / pallidum / M1 E_TIfree | 0.064 / 0.11 / 0.25 V/m | SimNIBS example montage F5-P5 / F6-P6, 1 mA, 20 mm: hippocampus 0.14, putamen 0.10, M1 hand 0.135 V/m (ROI medians); GM median 0.097, p99 0.18 (SimNIBS's nominal-current fields were 4 and 9 % higher for the two pairs) |
| optimised, 2 mA total: L hippocampus / R pallidum / L M1 | 0.24 / 0.37 / 0.56 | L hippocampus: `examples/ernie_montage_search.py`, lead field over the 10-10 cap (75 solves, 6 min), 189,750 two-pair montages from the 25 electrodes within 90 mm of the target at ratios 1:1, 1:2, 1:3 (29 s): best FC5-P9 + FC3-PO9 at 1 + 1 mA, ROI median **0.293 V/m**, p99 of the rest of the brain 0.42, but 198 cm³ of brain above 0.25 V/m: the intensity optimum is not focal, as Rampersad found. The focal ranking (montages within 70 % of the best ROI median, least brain above 0.25 V/m outside the ROI, 300 shortlisted) picks T7-P9 + C5-PO7 at 1 + 1 mA: ROI median **0.245 V/m** with 57 cm³ outside, 3.5× less spill for 16 % less intensity, and the same number as Rampersad's 0.24 |
| tACS same electrodes | slightly higher, 10-30× the volume > 0.25 V/m | example montage: brain > 0.2 V/m is 20 cm³ for TI vs 330 cm³ for tACS (17×); > 0.1 V/m 804 vs 1214 cm³ |

And Wessel 2023's montage on ernie (F3-F4 / TP7-TP8, 2 mA each, 9.8 mm electrodes):
putamen medians **0.35 / 0.37 V/m** (left / right), hippocampus 0.35 / 0.32, M1 hand
0.31 / 0.34, GM median 0.19, p99 0.43; brain > 0.2 V/m: 728 cm³ TI vs 1235 cm³ tACS. The
current correction was +13.4 % for F3-F4 and −5.5 % for TP7-TP8 (SimNIBS's estimate can err
either way; pair resistances 214 and 169 Ω), so the uncorrected 2026-10-25 numbers (0.39 /
0.40) also had the two carriers in the wrong ratio.
Wessel on MIDA: striatum 0.26 ± 0.04 (p99 0.35), hand knob 0.20. Ours is 1.4× theirs in the
striatum and 1.6× in M1, with the low-res ernie (SimNIBS skull 0.010 S/m vs the IT'IS LF
database, no anisotropy, 8 mm spherical ROIs instead of the striatal label). The deep-to-GM
median ratio is 1.9, the same qualitative picture as theirs. The MIDA run will show how much
of the factor is the head and how much the tables.

Different head (ernie vs their 25-year-old male), electrodes and conductivity tables, so
the comparison is order-of-magnitude and ratio (deep/cortex, TI/tACS), not equality.

### Placement sensitivity (Yamada 2026) and free placement (TesFlexOptimization)

`examples/ernie_placement_sensitivity.py`: each Wessel electrode moved ±5 mm along two
tangential scalp directions, 16 runs (23 min; corrected currents, 2026-10-26). Putamen
median changes: **F3 up to 7.0 %** (left) / 2.5 % (right), F4 2.1 / 6.0 %, TP7 and TP8
**0.0 %**: with the currents corrected the F3-F4 pair is the weaker field throughout both
putamina, so E_AM = 2|E_F3-F4| there and moving the temporo-parietal pair changes nothing
until the branch flips (the uncorrected run, with F3-F4 13 % too strong and TP7-TP8 5 % too
weak, had TP8 at 12.3 %: that number was an artefact of the mis-ratioed carriers). Yamada
report 0.4-14.9 % for ±10 mm grids. The deep hotspot (volume-weighted centroid of the top 1 %
of the envelope within 25 mm of the targets' midpoint; a single element's maximum is a mesh
artefact that never moves) shifts by up to 3.8 mm (F3), 4.5 (F4), 2.1 (TP7), 0.5 (TP8);
Yamada: 0.6-6.9 mm. Lesson for the montage papers: in a two-pair TI montage the envelope in
the target is set by the weaker pair, so placement tolerance is asymmetric between pairs.

`examples/ernie_flex_ti.py`: SimNIBS's differential-evolution placement of two free pairs
(10 mm, 1 mA each) maximising mean E_AM,max in the 8 mm left-hippocampus sphere. The smoke
run (2 iterations, population 4) took 98 s and ended at 0.26 V/m during the search (node-wise
sources) and 0.15 V/m in the final simulation with real electrodes: already Rampersad's
range at half his total current (the optimiser's own electrode model prescribes node
currents, so its fields are not subject to the `fem.tdcs` normalisation bias of §2, but this
has not been checked against `TwoPairSimulation` yet). The 50 × 13 run converged after 22 generations and
2392 FEM evaluations (11 min, 2.6 GB): best node-wise value 0.407 V/m; the final simulation
with real 10 mm electrodes gives a **ROI mean of 0.453 V/m (median 0.423) at 1 mA per
pair**, 1.5× the cap-restricted lead-field optimum (0.293 at 2 mA) and 1.9× Rampersad's
0.24: free placement finds positions the 10-10 cap does not offer. The price is spill:
511 cm³ of other brain above 0.25 V/m and a GM p99 of 0.465 V/m, i.e. the hippocampus gets
about what the most exposed cortex gets. Focality needs the "focality" goal with
thresholds (SimNIBS's tes_flex_ti_focality example), which is the next run.

## 4. Violante 2023 and Wessel 2023 on MIDA (needs the data)

| paper | montage | model numbers |
|---|---|---|
| Violante 2023 | e1-e3 (left, nasion plane, ~5 cm apart) and e2-e4 (right, above the eyebrow, ~16 cm apart), 1 mA per pair, 2 cm electrodes, DTI tensors | hippocampus 0.26 ± 0.04 V/m; cortex anterior 0.18, middle 0.12, posterior 0.10 |
| Wessel 2023 | F3-F4 and TP7-TP8, 2 mA each, 9.8 mm electrodes, IT'IS LF database v4.0 | striatum 0.26 ± 0.04 (median 0.27, p99 0.35); precentral hand knob 0.20; tDCS control 0.17 |

`examples/mida_ti.py` runs both once `TI_SIM_MIDA` points at a MIDA download and
`TI_SIM_MIDA_ELECTRODES` gives the electrode coordinates (MIDA has no cap file; Violante's
exact positions are only in their figure). Expect agreement within the Ivanov 2026 mesh
noise (2-4 %) plus the conductivity-table and isotropy differences (Yatsuda: < 12 % in
deep targets). Wessel's montage is the cleaner target because its electrodes are named.

## 5. Phantom and cadaver measurements (planned)

* Ghabel Damirchi 2026: paediatric phantom with six dipoles, HD-EEG reconstruction; the
  phantom's geometry would have to be requested.
* Fathi 2026: cadaver insula with six depth electrodes, steerable with the current ratio.
* Violante 2023 cadaver: ~0.1 V/m per mA at 44 mm depth, hippocampus 75 % above cortex.

## 6. Neural response (planned)

Esmaeilpour 2021's slice thresholds (`ti_sim.membrane.ESMAEILPOUR_2021`) and Gomez-Tames
2021's mouse thresholds are the targets for a neuron-model readout.
