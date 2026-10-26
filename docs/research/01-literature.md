# Temporal interference stimulation: the literature, 2017 to 2026

*Survey date: 2026-10-25. Question from Morgan: evaluate the current work including the
original paper and everything since, and how much can be reused from SimNIBS. Sources:
PubMed (162 hits for "temporal interference" + simulation/model terms, 2017-2026, all
titles read, abstracts of the modelling papers read in full), the open-access full texts
of Grossman 2017, Rampersad 2019, Violante 2023 and Wessel 2023, and the SimNIBS 4.6
source tree. DOIs are in [references.md](references.md).*

## The founding paper (Cell, 1 June 2017, not 2016)

Grossman N, Bono D, Okun MS, ... Boyden ES. *Noninvasive deep brain stimulation via
temporally interfering electric fields.* Cell 169:1029-1041, published 1 June 2017
(PubMed lists 2017 Jun 01). The idea: two kHz fields E1 cos 2πf1t and E2 cos 2πf2t with
Δf = |f1 - f2| in the neural range (10 Hz) sum to a carrier whose amplitude beats at Δf.
Neurons are assumed not to follow the carrier (membrane low-pass) but to follow the beat,
and the beat is largest where the two fields are equal and parallel, which can be deep.
Along a direction n the modulation amplitude is

    E_AM(n) = | |(E1 + E2) · n| - |(E1 - E2) · n| |

and its maximum over n has a closed form (SimNIBS implements it; `ti_sim.envelope`).

What the paper contains, and what of it is checkable:

| part | content | checkable? |
|---|---|---|
| Saline phantom (fig. 2) | 50 mm × 10 mm saline disc, 1 mm silver wires on the rim, 1 kHz / 1.02 kHz, 1 mA per pair; five trapezoid layouts; a table of simulated and measured peaks, positions and half-widths | **yes, replicated here**: positions and half-widths match to 0.03 R; amplitudes are a constant ~7× above their simulated numbers (open question, see [02](02-replication-targets.md)) |
| Sphere models (fig. S2) | homogeneous and 4-layer spheres (scalp 0.05 R / 0.333 S/m, skull 0.085 R / 0.0083, CSF 0.023 R / 1.79, brain 0.83 R / 0.333) | yes: `examples/layered_sphere.py` builds the 4-layer sphere through the MIDA route |
| Mouse | IT'IS mouse model at 42 × 42 × 700 µm, Sim4Life; 2 kHz / 2.01 kHz; hippocampal stimulation at 125 µA per channel without cortical c-fos; motor steering with current ratios 1:8 to 8:1 at ~0.78 mA total | partly: Rampersad 2019 remodelled it (383 V/m in the mouse brain at 0.776 mA) |
| Neurons | demonstrated in vivo (patch, c-fos), no biophysical model of *why* neurons demodulate | no; this is where the critiques went |

## What happened since, by theme

### 1. Human head modelling: fields are weak, focality is real, steering is real

* **Rampersad 2019** (NeuroImage): the first careful human study. SimNIBS/SCIRun head with
  DTI anisotropy, 88 electrode positions, 146 M current patterns. Standard montage at
  1 mA per pair: 0.064 V/m in hippocampus, 0.11 in pallidum, 0.25 in M1; optimised
  four-electrode montages at 2 mA total: 0.24 (hippocampus), 0.37 (pallidum), 0.56 (M1).
  tACS with the same electrodes is slightly stronger but covers 10-30× more brain above
  0.25 V/m. Suprathreshold stimulation in humans would need ~38 mA per pair. Their mouse
  model confirms Grossman's mouse was suprathreshold (383 V/m).
* **Huang & Parra 2019/2020**: the modulation-depth optimum has an analytic solution
  (it is the HD-tES solution with a modulated source); focality can beat HD-tES but
  intensity is weak and the locus of maximal intensity is not the locus of maximal
  modulation. 2023: visualisations on the MNI-152 head; 2026 preprint: two *arrays* of
  electrodes instead of two pairs, optimised in under 30 s, 3.0 cm focality at the
  hippocampus.
* **Lee et al. 2020** (Sci Rep): individualised exhaustive two-pair search on three heads
  for the right hippocampal head; the optimum differs per head, so montages must be
  personalised. **Lee 2022**: multipair TI improves focality. **Lee 2026**: target-specific
  restricted exhaustive search, 30-100× faster with the same result.
* **von Conta 2021** (Sci Rep): 25 subjects, large inter-individual variability of the
  field in hippocampus, M1 and thalamus; tTIS more focal than tACS with comparable ROI
  fields. **von Conta 2022** (Cortex): no α-power aftereffect in 34 humans (first null).
* **Karimi 2019**, **Stoupis 2022** (genetic algorithm), **Song 2021** (six-channel,
  mice), **Terasawa 2021** (phase modulation), **Cao 2019** (STIMULUS: spatio-temporal
  patterns beat plain TI in a HH model), **Savvateev 2026** (Cell Systems): a third,
  180°-phase-shifted pair cancels off-target modulation in mice, with fMRI.
* **Numerical hygiene, 2025-2026**: **Ivanov 2026a** finds MNI152-template montages
  underestimate deep-target fields (16th-21st percentile of 132 CamCAN heads) while
  overestimating coverage ratios; **Ivanov 2026b** shows independent remeshing of the
  same labels gives 2-4 % CoV in the hippocampal envelope (average 5-10 meshes);
  **Yatsuda 2026** puts WM anisotropy at up to 18 % in WM, < 12 % in deep targets, and
  50 % montage divergence unless tolerances are allowed; **Yamada 2026**: ±10 mm
  electrode displacement changes ROI intensity by 0.4-15 % and moves the hotspot
  0.6-6.9 mm; **Meng 2025**: skull thickness and CSF fraction predict the hippocampal
  modulation depth without a simulation; **Noetscher 2026** (ECT, not TI) on hippocampal
  field convergence.

### 2. Biophysics: why (and whether) neurons demodulate

* **Mirzakhalili 2020** (Cell Systems): passive membrane filtering cannot demodulate a
  sum of two sinusoids (it is still a sum of two sinusoids); demodulation needs an
  ion-channel rectification, which is the same mechanism that produces kHz conduction
  block off target. Proposes controls for every TI experiment. `ti_sim.membrane` shows
  the arithmetic.
* **Esmaeilpour 2021** (Brain Stimul): rat hippocampal slices. Modulating carbachol
  gamma needs ~5 V/m with a 100 Hz carrier but ~60 V/m at 1 kHz and ~80 V/m at 2 kHz
  (~160-220 mA at the scalp); selectivity comes from network adaptation faster than the
  beat. **Esmaeilpour 2020** (eNeuro): no synaptic effect of unmodulated 1-10 kHz fields
  below 60 V/m.
* **Howell & McIntyre 2020**: at the human scale, interferential tACS cannot activate
  axons at ≤ 10 mA but can modulate spiking indirectly, less potently than plain tACS at
  matched intensity; no gain in selectivity. **Gomez-Tames 2021**: multiscale mouse
  model; threshold rises with carrier frequency, is independent of the beat frequency;
  activation follows regions of high modulation depth *and* low envelope minimum.
  **Cerpa 2022**, **Plovie 2022** (point neurons), **Budde 2023** and **Jabban 2025**
  (peripheral nerves: thresholds do not depend on modulation depth, so TI is not an
  envelope extractor there).
* Variants: pulse-width-modulated TI (**Luff 2023**), intersectional short pulses
  (**Botzanowski 2022**, ILTI), EM-wave TI (**Ahsan 2022**), TI-TMS (**Sorkhabi 2020**,
  **Wang 2026**), epidural TI (**Lee 2022**, with a skull-phantom validation).

### 3. Human validation with models on MIDA (the KCL connection)

* **Violante et al. 2023** (Nat Neurosci), Grossman's group with Inês Violante (now at
  King's College London, School of Biomedical Engineering & Imaging Sciences): Sim4Life
  on **MIDA** (0.5 mm, >100 tissues, DTI anisotropy), 2 kHz / 2.005 kHz, 1 mA per pair in
  the model, 2 cm electrodes: hippocampal envelope **0.26 ± 0.04 V/m**, cortex 0.10-0.18
  (30-60 % lower); a cadaver confirmed ~75 % more envelope in hippocampus than overlying
  cortex, ~0.1 V/m per mA at 44 mm depth; then fMRI (n = 20, 2 mA) and behaviour (n = 21):
  focal hippocampal modulation and better episodic memory accuracy.
* **Wessel et al. 2023** (Nat Neurosci, Hummel lab): Sim4Life 5.0 on **MIDA** (117
  tissues), F3-F4 and TP7-TP8, 2 mA each, 9.8 mm electrodes: striatum **0.26 ± 0.04 V/m**
  (median 0.27, p99 0.35), precentral hand knob 0.20; theta-burst TI (2.00 / 2.10 kHz
  bursts) raised striatal fMRI activity and motor learning.
* Violante's KCL group since: **Kurtin 2025** (Brain Stimul, hippocampal TI changes
  post-task resting connectivity), the 2025 bioRxiv memory-encoding preprint with Brno
  (Missey, Williamson) which modelled four-pair, four-carrier TI on the ViP **Jeduk V4.0**
  model rather than MIDA, and two 2026 syntheses she co-authored: a meta-analysis of 15
  behavioural TI studies (Hedges' g = 0.30; **Aziziaram 2026**) and a safety review of 24
  studies, 853 participants, no serious adverse events at ≤ 2 mA per channel
  (**Demchenko 2026**). A KCL seminar in 2026 was titled "Modulating memory and brain
  dynamics with temporal interference stimulation". So "a KCL researcher simulating TI on
  MIDA" fits that group; MIDA is also the model of the 2022 IT'IS/FDA safety preprint.
* Other human work: **Zhu 2022** (M1 connectivity), **Kurtin 2025**, **Zhang 2026**
  (nucleus accumbens fMRI), **Wang 2026** (Cell Rep Med, high-current sgACC TI RCT in
  depression), **Chen 2026** (DLPFC/dACC), **Lyu 2026** (mPFC, chronic pain), clinical
  reviews (**Guo 2023**, **Zhu 2023**, **Albantakis 2026**, **Wang 2026**).

### 4. Physical validation outside the brain

* **Ghabel Damirchi 2026** (NeuroImage): 3D-printed paediatric head phantom with six
  implanted dipoles; TI envelopes reconstructed from HD-EEG; 7-9 kHz carriers gave the
  most confined envelopes (7-23 mm error vs 50-64 mm for a low-frequency control).
* **Fathi 2026**: fresh-frozen cadaver, six depth electrodes in the insula, 990/1000 Hz,
  simulation-guided montages from the "TI Planning Tool"; steering with the current ratio.
* Grossman 2017's phantom and Violante 2023's cadaver (above); **Lee 2022**'s skull phantom.

### 5. Software

* **SimNIBS 4.x** has TI built in: `simnibs.utils.TI_utils` (`get_maxTI`, `get_dirTI`,
  lead-field helpers), two examples (`examples/simulations/TI.py`,
  `examples/optimization/TI_fields_leadfields.py`) and the 2024-25 flexible electrode
  optimiser `TesFlexOptimization` with `e_postproc = "max_TI"` / `"dir_TI"` (two
  freely-placed pairs, differential evolution). See [03](03-simnibs-reuse.md).
* **TI-Toolbox** (Haber, Jackson, Thielscher, Hai, Tononi; Brain Stimul 2026; GPL-3,
  github.com/idossha/TI-Toolbox, v3.1.2, pushed 2026-10-24): a Docker + Electron
  application around dcm2niix, FreeSurfer, SimNIBS charm, flex-search and exhaustive
  montage search, ROI analysis and reports; also a multi-pair (mTI) envelope kernel.
  It is the end-user tool for individual MRIs; it does not do phantoms, spheres or MIDA.
* **Sim4Life** (ZMT, commercial): what Grossman, Violante and Wessel used, with MIDA.
* ROAST/Huang's MATLAB code for the two-array optimiser; **TI Planning Tool** (Fathi).

## Assessment

1. **The physics is settled and reproducible.** The envelope is a geometric function of
   two quasi-static fields; everything from Grossman's phantom table to the MIDA numbers
   is a Laplace solve per pair plus equation 1. This repository reproduces the phantom
   geometry exactly and runs the head route with SimNIBS.
2. **The numbers in humans are small.** 0.2-0.4 V/m envelopes at 2 mA, comparable to
   tACS in the target but much more focal. Direct activation (Grossman's mouse, hundreds
   of V/m) does not scale to humans at safe currents (Rampersad: ~38 mA per pair).
3. **The mechanism is the open question.** Mirzakhalili, Esmaeilpour, Howell, Budde and
   Jabban all show that a passive membrane does not demodulate and that the kHz carrier
   itself matters (threshold rises with carrier frequency; peripheral thresholds ignore
   modulation depth). The human effects (Violante, Wessel, the 2026 meta-analysis) are
   real but modest (g ≈ 0.3) and their cellular route is not established.
4. **Modelling practice is converging** on: individual heads (not MNI templates),
   multiple mesh realisations, reported electrode-placement sensitivity, anisotropy where
   DTI exists, and MIDA as the shared reference head. The Ivanov/Yatsuda/Yamada papers
   are a methods checklist this repository's examples should satisfy.
5. **What is missing from the open tools** and worth building here: a validated phantom
   and sphere path (done), MIDA through SimNIBS instead of Sim4Life (built, untested
   without the data), the membrane/rectification arithmetic beside the field numbers,
   and a neuron-model readout (MRG axon or Blue Brain cortical cells under a modulated
   kHz field) so a montage is scored by a response, not only by E_AM.
