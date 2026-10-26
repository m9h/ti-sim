# MIDA in SimNIBS

**MIDA** (Multimodal Imaging-Based Detailed Anatomical model; Iacono et al. 2015, PLoS ONE
10:e0124126; IT'IS Foundation and the US FDA; DOI 10.13099/ViP-MIDA-V1.0): head and neck of
a healthy 29-year-old woman at 500 µm isotropic, 115 structures plus background, from
T1/T2/MRA/DTI; distributed as voxel (MAT, RAW, NII, TXT) and STL files. It is the head the
2023 human TI papers were modelled on (Sim4Life), and the one a KCL colleague is using.

## Licence (LicenseAgreementMIDA, IT'IS, read 2026-10-25)

Free, non-exclusive, unlimited in time, obtained through the form on the "MIDA V1.0" page.
Binding clauses for us:

* 2.3.2 **the model data or works derived from it must not be distributed or
  redistributed**, original or modified: so no copy, crop, mesh or lead field of MIDA goes
  into this repository, the COPR or a container image. `scripts/check-mida.sh` validates a
  user's own download; `ti_sim.mida` caches meshes under the user's data directory.
* 2.3.3 images based on the data may be published only with the face disguised.
* 2.3.8 publications must credit the FDA CDRH and IT'IS as creators and cite Iacono et al.
* 2.3.5 no military use. 2.3.10 use ends if the subject withdraws consent.

The PLoS ONE *article* is CC0; the *model* is not.

## The 115 structures to SimNIBS tags

`ti_sim.mida.RULES` maps structure names (Table 1 of the paper, embedded as `TABLE1` and
tested) to SimNIBS tags, first match wins:

| SimNIBS tag | MIDA structures |
|---|---|
| 0 air | Air Internal (7 sinuses/cavities), Ear Auditory Canal, Ear Pharyngotympanic Tube |
| 1 WM (0.126 S/m) | Brain and Cerebellum White Matter, Commissura Anterior/Posterior, Optic Chiasm, Optic Tract, Cranial Nerves I-XII (no nerve class in SimNIBS; WM is the nearest) |
| 2 GM (0.275) | Brain and Cerebellum Gray Matter, Pineal Body, Hypophysis |
| 3 CSF (1.654) | CSF General, CSF Ventricles, Dura (charm's convention), Ear Cochlea, Ear Semicircular Canals |
| 4/7/8 bone | Skull → 4; Skull Outer/Inner Table, Mandible, Hyoid, Vertebrae C1-C5, Teeth → 7 compact (0.008); Skull Diploë → 8 spongy (0.025) |
| 5 scalp (0.465) | Epidermis/Dermis, Mucosa, Parotid/Sublingual/Submandibular Gland |
| 6 eyes (0.5) | Eye Aqueous, Cornea, Lens, Retina/Choroid/Sclera, Vitreous |
| 9 blood (0.6) | Blood Arteries, Blood Veins |
| 10 muscle (0.16) | all 30 Muscle entries, Tongue, Tendons |
| 11 cartilage (0.88) | Ear Auricular Cartilage, Nasal Septum, Intervertebral Disc |
| 12 fat (0.078) | Adipose Tissue, Subcutaneous Adipose Tissue |
| 21-29 (GM, 0.275) | Hippocampus 21, Amygdala 22, Thalamus 23, Caudate 24, Putamen 25, Globus Pallidus 26, Nucleus Accumbens 27, Substantia Nigra 28, Hypothalamus + Mammillary Body 29: kept separate as ROIs |
| 30 brainstem (0.201, Rampersad) | Brainstem Medulla/Midbrain/Pons, Cerebral Peduncles, Spinal Cord |

Choices to revisit with the data in hand: nerves as WM, glands as scalp, dura as CSF, the
thalamus's nuclei (Table 2 of the paper) merged. Conductivities are the SimNIBS table by
default (`TI_SIM_COND=rampersad` for theirs); Violante used a meta-analysis table plus DTI
tensors, Wessel the IT'IS LF database v4.0, so a few percent of difference is expected
before any geometry.

## Mesh

`ti_sim.mida.prepare`: load NII + TXT, map, optionally crop below `z_min_mm` (the neck is
not needed for TI) and downsample 0.5 → 1 mm (nearest neighbour), then SimNIBS
`create_mesh` with the deep nuclei in the surface hierarchy right after WM, and
`relabel_internal_air`. Cached by a hash of the inputs. Expect a few million tetrahedra
at 1 mm with elements of 1.5-6 mm; the full 0.5 mm volume is possible but slow.

## Electrodes

MIDA has no EEG cap file. `examples/mida_ti.py` takes a JSON of electrode coordinates in
MIDA world space (`TI_SIM_MIDA_ELECTRODES`); for 10-10 names, the quickest route is to
place the fiducials (nasion, inion, preauricular points) on the skin surface and compute
the positions with SimNIBS's `eeg_positions` utilities (MNI cap warped to the head) or a
10-10 geodesic construction. Wessel 2023's montage (F3, F4, TP7, TP8) is the first target.
