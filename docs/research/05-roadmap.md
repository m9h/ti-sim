# Roadmap

Ordered by what unblocks what.

## Done (2026-10-25)

- [x] Envelope maths with three independent checks (closed form, direction sweep, SimNIBS).
- [x] Grossman 2017 phantom: geometry reproduced, amplitude factor recorded.
- [x] SimNIBS head route on ernie: two montages, ROI statistics, focality vs tACS, gmsh output.
- [x] Label volume → SimNIBS mesh → TI (layered sphere), the route MIDA takes.
- [x] MIDA loader and mapping (all 115 structures), licence handling, check script.
- [x] Lead-field montage search with current ratios and candidate restriction.
- [x] Membrane low-pass / rectification arithmetic and the Esmaeilpour thresholds.
- [x] BATS pipeline, container recipe.

## Next

* Current normalisation (note 02 §2): `TwoPairSimulation` corrected, lead field verified
  unbiased (`montage.calibrate_leadfield`); check whether `TesFlexOptimization`'s
  node-current electrode model needs it at all (compare one montage against `TwoPairSimulation`).
* Sphere: the example now reports the top-1 % centroid depth, centre means and shell
  medians instead of an argmax, but `create_mesh(elem_sizes={"standard": {"range": [1.5, 3.0]}})`
  still gives ~6 mm brain elements (49k tets): find the parameter that refines the interior
  (or mesh the sphere with `scripts/gmsh-sphere-mesh.py`, which the cross-check already uses).

1. ~~**MIDA with the data**~~: done 2026-10-26 on the DGX Spark (see [02 §4](02-replication-targets.md)):
   parser fixed, meshed at 2, 1 and 0.5 mm (converged at 1 mm, 3.5 M tets, 44 min, 5.5 GB),
   Wessel's montage 0.49 V/m putamen with the SimNIBS table, 0.45 with Rampersad's; half the
   peak-to-trough envelope is 0.23, inside Wessel's 0.26 ± 0.04. Still open: Violante's e1-e4
   positions, and a Slicer check of the fiducial-derived cap.
2. ~~Electrode positions on MIDA~~: done 2026-10-25, `ti_sim.electrodes.ten_ten` (geodesic
   10-10 from Nz/Iz/LPA/RPA on any skin surface; 6.4 mm mean against ernie's template cap).
   Still needed: the four MIDA fiducial coordinates.
3. **The amplitude factor** in the phantom: write to Grossman's group? (Morgan's call;
   draft in docs/upstream-drafts.md) or find the definition in the Cell supplement.
4. ~~Lead field on ernie~~: done 2026-10-25 (6 min; L hippocampus optimum 0.293 V/m at
   2 mA, see [02](02-replication-targets.md)); `TesFlexOptimization` wrapped as
   `examples/ernie_flex_ti.py` (smoke run); a focality-aware objective for `ti_sim.montage`
   (Rampersad's volume-above-threshold, Lee 2020's ratio) is the next step.
5. **Mesh-realisation and placement sensitivity** (Ivanov 2026, Yamada 2026): placement
   done as `examples/ernie_placement_sensitivity.py` (±5 mm, 16 runs); mesh realisations
   need charm reruns (not on the low-res ernie).
6. **Anisotropy**: `anisotropy_type='vn'` on a charm head with DTI (ernie lowres has none).
7. **Neuron readout**: started with `ti_sim.hh`, a Hodgkin-Huxley point neuron polarised
   by lambda E(t) (Radman's 0.2 mm), with a bisection for the firing threshold vs carrier
   frequency (rises with the carrier, as Gomez-Tames 2021). NEURON is not in Fedora but
   `uv pip install neuron` resolves on Python 3.14; an MRG axon and a layer-5 cell under
   E(t) from a montage are the next step.
8. **Multi-pair and phase-shifted TI** (Savvateev 2026, Huang 2026 two-array): done
   2026-10-26 for the phase-free case, `envelope.multi_pair_envelope` (exact: with
   a_i = E_i·n the envelope is min(Σ|a_i|, 2(Σ|a_i| − max|a_i|)), which is Grossman's
   equation 1 for two carriers and ≥ his 2 min|E_i| for more), `max_multi_pair_envelope`
   (direction sweep; no closed form for k > 2) and `multi_tone_envelope_bounds` (the
   phase-locked value from the waveform, needed when the carriers are equally spaced and
   the beats are commensurate: the phase-free formula is then only an upper bound). Still to
   do: cross-check against TI-Toolbox's mTI kernel on one montage, and use it in `montage`.
9. **mepkit / open-neuronav tie-in**: TI montages as NaviNIBS targets; stimulator
   control stays out of scope.
10. **Fedora package** `python-ti-sim` once the API settles (depends on python3-simnibs,
    scikit-fem is not packaged: either package it or keep ti-sim a uv project).
