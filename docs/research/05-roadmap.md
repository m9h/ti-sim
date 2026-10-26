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
* Sphere: size field to ~3 mm in the brain (the 6 mm interior elements made the envelope
  argmax meaningless) and report the centre value and shell medians, not an argmax depth.

1. **MIDA with the data** (needs Morgan's download, or the KCL colleague's): run
   `scripts/check-mida.sh`, fix the label-file parser to the real layout, mesh at 1 mm,
   Wessel's montage, compare to 0.26 V/m in striatum; then Violante's.
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
8. **Multi-pair and phase-shifted TI** (Savvateev 2026, Huang 2026 two-array): extend
   `envelope` to n carriers with the exact envelope (TI-Toolbox's mTI kernel is the
   reference).
9. **mepkit / open-neuronav tie-in**: TI montages as NaviNIBS targets; stimulator
   control stays out of scope.
10. **Fedora package** `python-ti-sim` once the API settles (depends on python3-simnibs,
    scikit-fem is not packaged: either package it or keep ti-sim a uv project).
