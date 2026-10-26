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

## 2. Grossman 2017, figure S2: the 4-layer sphere (built)

`examples/layered_sphere.py`: the sphere with the paper's proportions, scaled to R =
85 mm, built as a 1 mm label volume, meshed with SimNIBS `create_mesh` and solved with
SimNIBS, two pairs on the equator in layout C and the 1:1, 1:2.5, 1:4 ratios. The paper
gives no table for the sphere; what the pipeline checks is qualitative and theirs: the
envelope peaks deeper than the unmodulated field and moves towards the weaker pair.

## 3. Rampersad 2019: human head numbers (ernie stands in)

| quantity | paper | ours (ernie low-res, SimNIBS conductivities, isotropic; 2026-10-25) |
|---|---|---|
| standard montage, 1 mA per pair: hippocampus / pallidum / M1 E_TIfree | 0.064 / 0.11 / 0.25 V/m | SimNIBS example montage F5-P5 / F6-P6, 1 mA, 20 mm: hippocampus 0.15, putamen 0.11, M1 hand 0.145 V/m (ROI medians); GM median 0.10, p99 0.19, max 0.52 |
| optimised, 2 mA total: L hippocampus / R pallidum / L M1 | 0.24 / 0.37 / 0.56 | `examples/ernie_montage_search.py` (TI_SIM_LEADFIELD=1), not run yet |
| tACS same electrodes | slightly higher, 10-30× the volume > 0.25 V/m | example montage: brain > 0.2 V/m is 31 cm³ for TI vs 430 cm³ for tACS (14×); tACS GM p99 0.32 vs TI 0.19 |

And Wessel 2023's montage on ernie (F3-F4 / TP7-TP8, 2 mA each, 9.8 mm electrodes):
putamen medians **0.39 / 0.40 V/m** (left / right), hippocampus 0.40 / 0.36, M1 hand
0.36 / 0.37, GM median 0.20, p99 0.44; brain > 0.2 V/m: 795 cm³ TI vs 1235 cm³ tACS.
Wessel on MIDA: striatum 0.26 ± 0.04 (p99 0.35), hand knob 0.20. Ours is 1.5× theirs in the
striatum and 1.8× in M1, with the low-res ernie (SimNIBS skull 0.010 S/m vs the IT'IS LF
database, no anisotropy, 8 mm spherical ROIs instead of the striatal label). The deep-to-GM
median ratio is 1.9, the same qualitative picture as theirs. The MIDA run will show how much
of the factor is the head and how much the tables.

Different head (ernie vs their 25-year-old male), electrodes and conductivity tables, so
the comparison is order-of-magnitude and ratio (deep/cortex, TI/tACS), not equality.

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
