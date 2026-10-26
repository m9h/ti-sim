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

1. **MIDA with the data** (needs Morgan's download, or the KCL colleague's): run
   `scripts/check-mida.sh`, fix the label-file parser to the real layout, mesh at 1 mm,
   Wessel's montage, compare to 0.26 V/m in striatum; then Violante's.
2. **Electrode positions on MIDA**: 10-10 from fiducials for custom meshes.
3. **The amplitude factor** in the phantom: write to Grossman's group? (Morgan's call;
   draft in docs/upstream-drafts.md) or find the definition in the Cell supplement.
4. **Lead field on ernie** (`TI_SIM_LEADFIELD=1`, 20-40 min) and the Rampersad
   comparison; then wrap `TesFlexOptimization` for free electrode placement.
5. **Mesh-realisation and placement sensitivity** (Ivanov 2026, Yamada 2026): repeat the
   ernie montages over N meshes / ±5 mm displacements; report CoV.
6. **Anisotropy**: `anisotropy_type='vn'` on a charm head with DTI (ernie lowres has none).
7. **Neuron readout**: NEURON (python3-neuron is in Fedora) with an MRG axon and a
   layer-5 pyramidal cell under E(t) from a montage; score montages by response, and
   reproduce Gomez-Tames's carrier-frequency threshold trend.
8. **Multi-pair and phase-shifted TI** (Savvateev 2026, Huang 2026 two-array): extend
   `envelope` to n carriers with the exact envelope (TI-Toolbox's mTI kernel is the
   reference).
9. **mepkit / open-neuronav tie-in**: TI montages as NaviNIBS targets; stimulator
   control stays out of scope.
10. **Fedora package** `python-ti-sim` once the API settles (depends on python3-simnibs,
    scikit-fem is not packaged: either package it or keep ti-sim a uv project).
