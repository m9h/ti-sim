"""The MIDA head model (Iacono et al. 2015, PLoS ONE 10:e0124126) as a SimNIBS mesh.

MIDA is the model behind the two human TI validation papers (Violante 2023, Wessel 2023,
both Sim4Life). It is a 0.5 mm label volume with 115 structures (plus background) from the
IT'IS Foundation and the US FDA. Its licence (LicenseAgreementMIDA, IT'IS) forbids
redistributing the data or derivatives, so nothing of it is in this repository or any
package: `scripts/check-mida.sh` tells you where to register, and TI_SIM_MIDA points at
your copy. Published images must disguise the face (clause 2.3.3) and papers must cite the
FDA/IT'IS and Iacono et al. (clause 2.3.8).

`load` reads the NIfTI label volume and the label-name text file, `classify` maps the
structure names onto SimNIBS tissue tags (plus extra tags for the deep nuclei that TI
papers target, kept separate so they can be ROIs), `prepare` crops, optionally downsamples,
relabels and meshes it with SimNIBS create_mesh (`ti_sim.labelvol`). The name table is the
published Table 1; the mapping is explicit and listed in docs/research/04-mida.md.
"""

import glob
import os
import re

import numpy as np

from . import labelvol

# extra tags (after the 12 SimNIBS tissues) for structures kept as their own regions
HIPPOCAMPUS, AMYGDALA, THALAMUS, CAUDATE, PUTAMEN, PALLIDUM, ACCUMBENS, NIGRA, HYPOTHALAMUS, BRAINSTEM = range(21, 31)
EXTRA_NAMES = {21: "hippocampus", 22: "amygdala", 23: "thalamus", 24: "caudate", 25: "putamen",
               26: "globus pallidus", 27: "nucleus accumbens", 28: "substantia nigra",
               29: "hypothalamus", 30: "brainstem"}
# conductivities for the extra tags: grey matter for the nuclei (SimNIBS 0.275), brainstem
# 0.201 (Rampersad 2019)
EXTRA_COND = {**{t: 0.275 for t in range(21, 30)}, 30: 0.201}

# (lower-case substring of the MIDA structure name, tag); first match wins
RULES = [
    ("background", 0), ("air internal", 0), ("auditory canal", 0), ("pharyngotympanic", 0),
    ("hippocampus", HIPPOCAMPUS), ("amygdala", AMYGDALA), ("thalamus", THALAMUS),
    ("caudate", CAUDATE), ("putamen", PUTAMEN), ("globus pallidus", PALLIDUM),
    ("accumbens", ACCUMBENS), ("substantia nigra", NIGRA), ("hypothalamus", HYPOTHALAMUS),
    ("mammillary", HYPOTHALAMUS),
    ("brainstem", BRAINSTEM), ("spinal cord", BRAINSTEM), ("cerebral peduncles", BRAINSTEM),
    ("white matter", 1), ("commissura", 1), ("optic chiasm", 1), ("optic tract", 1),
    ("cranial nerve", 1),
    ("gray matter", 2), ("grey matter", 2), ("pineal", 2), ("hypophysis", 2), ("pituitary", 2),
    ("csf", 3), ("dura", 3), ("cochlea", 3), ("semicircular", 3),
    ("eye", 6),
    ("skull diplo", 8), ("skull outer table", 7), ("skull inner table", 7), ("skull", 4),
    ("mandible", 7), ("hyoid", 7), ("vertebra", 7), ("teeth", 7),
    ("intervertebral", 11), ("cartilage", 11), ("nasal septum", 11),
    ("blood", 9), ("arter", 9), ("vein", 9),
    ("muscle", 10), ("tongue", 10), ("tendon", 10),
    ("adipose", 12),
    ("epidermis", 5), ("dermis", 5), ("skin", 5), ("mucosa", 5), ("gland", 5),
]

# The download's MIDA_v1.txt numbering (116 entries including 50 = Background), as
# recorded by brain-fwi's loader from a licensed copy: names only, no model data.
MIDA_V1_LABELS = {1: "Dura", 2: "Cerebellum Gray Matter", 3: "Pineal Body", 4: "Amygdala", 5: "Hippocampus",
    6: "CSF Ventricles", 7: "Caudate Nucleus", 8: "Putamen", 9: "Cerebellum White Matter", 10: "Brain Gray Matter",
    11: "Brainstem Midbrain", 12: "Brain White Matter", 13: "Spinal Cord", 14: "Brainstem Pons", 15: "Brainstem Medulla",
    16: "Nucleus Accumbens", 17: "Globus Pallidus", 18: "Optic Tract", 19: "Hypophysis (Pituitary Gland)",
    20: "Mammillary Body", 21: "Hypothalamus", 22: "Commissura (Anterior)", 23: "Commissura (Posterior)",
    24: "Blood Arteries", 25: "Blood Veins", 26: "Air Internal - Ethmoidal Sinus", 27: "Air Internal - Frontal Sinus",
    28: "Air Internal - Maxillary Sinus", 29: "Air Internal - Sphenoidal Sinus", 30: "Air Internal - Mastoid",
    31: "Air Internal - Nasal/Pharynx", 32: "CSF General", 33: "Ear Cochlea", 34: "Ear Semicircular Canals",
    35: "Ear Auricular Cartilage (Pinna)", 36: "Mandible", 37: "Mucosa", 38: "Muscle (General)",
    39: "Nasal Septum (Cartilage)", 40: "Skull", 41: "Teeth", 42: "Tongue", 43: "Adipose Tissue",
    44: "Vertebra - C1 (atlas)", 45: "Vertebra - C2 (axis)", 46: "Vertebra - C3", 47: "Vertebra - C4", 48: "Vertebra - C5",
    49: "Intervertebral Discs", 50: "Background", 51: "Epidermis/Dermis", 52: "Skull Diploe", 53: "Skull Inner Table",
    54: "Skull Outer Table", 55: "Eye Lens", 56: "Eye Retina/Choroid/Sclera", 57: "Eye Vitreous", 58: "Eye Cornea",
    59: "Eye Aqueous", 60: "Muscle - Platysma", 61: "Tendon - Galea Aponeurotica", 62: "Subcutaneous Adipose Tissue",
    63: "Muscle - Temporalis/Temporoparietalis", 64: "Muscle - Occipitiofrontalis - Frontal Belly",
    65: "Muscle - Lateral Pterygoid", 66: "Muscle - Masseter", 67: "Muscle - Splenius Capitis",
    68: "Muscle - Sternocleidomastoid", 69: "Muscle - Occipitiofrontalis - Occipital Belly", 70: "Muscle - Trapezius",
    71: "Muscle - Mentalis", 72: "Muscle - Depressor Anguli Oris", 73: "Muscle - Depressor Labii", 74: "Muscle - Nasalis",
    75: "Muscle - Orbicularis Oris", 76: "Muscles - Procerus", 77: "Muscle - Levator Labii Superioris",
    78: "Muscle - Zygomaticus Major", 79: "Muscle - Orbicularis Oculi", 80: "Muscle - Levator Scapulae",
    81: "Muscle - Medial Pterygoid", 82: "Muscle - Zygomaticus Minor", 83: "Muscles - Risorius", 84: "Muscle - Buccinator",
    85: "Ear Auditory Canal", 86: "Ear Pharyngotympanic Tube", 87: "Hyoid Bone", 88: "Submandibular Gland",
    89: "Parotid Gland", 90: "Sublingual Gland", 91: "Muscle - Superior Rectus", 92: "Muscle - Medial Rectus",
    93: "Muscle - Lateral Rectus", 94: "Muscle - Inferior Rectus", 95: "Muscle - Superior Oblique",
    96: "Muscle - Inferior Oblique", 97: "Air Internal - Oral Cavity", 98: "Tendon - Temporalis Tendon",
    99: "Substantia Nigra", 100: "Cerebral Peduncles", 101: "Optic Chiasm", 102: "Cranial Nerve I - Olfactory",
    103: "Cranial Nerve II - Optic", 104: "Cranial Nerve III - Oculomotor", 105: "Cranial Nerve IV - Trochlear",
    106: "Cranial Nerve V - Trigeminal", 107: "Cranial Nerve V2 - Maxillary Division",
    108: "Cranial Nerve V3 - Mandibular Division", 109: "Cranial Nerve VI - Abducens", 110: "Cranial Nerve VII - Facial",
    111: "Cranial Nerve VIII - Vestibulocochlear", 112: "Cranial Nerve IX - Glossopharyngeal", 113: "Cranial Nerve X - Vagus",
    114: "Cranial Nerve XI - Accessory", 115: "Cranial Nerve XII - Hypoglossal", 116: "Thalamus"}

# Table 1 of Iacono et al. 2015 (the 115 structures), for tests and for checking a download
TABLE1 = """Adipose Tissue|Air Internal-Ethmoidal Sinus|Air Internal-Frontal Sinus|Air Internal-Mastoid|
Air Internal-Maxillary Sinus|Air Internal-Nasal/Pharynx|Air Internal-Oral Cavity|Air Internal-Sphenoidal Sinus|
Amygdala|Blood Arteries|Blood Veins|Brain Gray Matter|Brain White Matter|Brainstem Medulla|Brainstem Midbrain|
Brainstem Pons|Caudate Nucleus|Cerebellum Gray Matter|Cerebellum White Matter|Cerebral Peduncles|
Commissura (Anterior)|Commissura (Posterior)|Cranial Nerve I-Olfactory|Cranial Nerve II-Optic|
Cranial Nerve III-Oculomotor|Cranial Nerve IV-Trochlear|Cranial Nerve V-Trigeminal Nerve|
Cranial Nerve V2-Maxillary Division|Cranial Nerve V3-Mandibular Division|Cranial Nerve VI-Abducens|
Cranial Nerve VII-Facial|Cranial Nerve VIII-Vestibulocochlear|Cranial Nerve IX-Glossopharyngeal|
Cranial Nerve X-Vagus|Cranial Nerve XI-Accessory|Cranial Nerve XII-Hypoglossal|CSF General|CSF Ventricles|Dura|
Ear Auditory Canal|Ear Auricular Cartilage (Pinna)|Ear Cochlea|Ear Pharyngotympanic Tube|Ear Semicircular Canals|
Epidermis/Dermis|Eye Aqueous|Eye Cornea|Eye Lens|Eye Retina/Choroid/Sclera|Eye Vitreous|Globus Pallidus|
Hippocampus|Hyoid Bone|Hypophysis Or Pituitary Gland|Hypothalamus|Intervertebral Disc|Mammillary Body|Mandible|
Mucosa|Muscle (General)|Muscle-Buccinator|Muscle-Depressor Anguli Oris|Muscle-Depressor Labii|
Muscle-Inferior Oblique|Muscle-Inferior Rectus|Muscle-Lateral Pterygoid|Muscle-Lateral Rectus|
Muscle-Levator Labii Superioris|Muscle-Levator Scapulae|Muscle-Masseter|Muscle-Medial Pterygoid|
Muscle-Medial Rectus|Muscle-Mentalis|Muscle-Nasalis|Muscle-Occipitiofrontalis-Frontal Belly|
Muscle-Occipitiofrontalis-Occipital Belly|Muscle-Orbicularis Oculi|Muscle-Orbicularis Oris|Muscle-Platysma|
Muscle-Procerus|Muscle-Risorius|Muscle-Splenius Capitis|Muscle-Sternocleidomastoid|Muscle-Superior Oblique|
Muscle-Superior Rectus|Muscle-Temporalis/Temporoparietalis|Muscle-Trapezius|Muscle-Zygomaticus Major|
Muscle-Zygomaticus Minor|Nasal Septum (Cartilage)|Nucleus Accumbens|Optic Chiasm|Optic Tract|Parotid Gland|
Pineal Body|Putamen|Skull|Skull Diploe|Skull Inner Table|Skull Outer Table|Spinal Cord|Subcutaneous Adipose Tissue|
Sublingual Gland|Submandibular Gland|Substantia Nigra|Teeth|Tendon-Galea Aponeurotica|Tendon-Temporalis|Thalamus|
Tongue|Vertebra-C1 (atlas)|Vertebra-C2 (axis)|Vertebra-C3|Vertebra-C4|Vertebra-C5""".replace("\n", "").split("|")


def classify(names):
    """{label: tag} for {label: structure name}; raises on a name no rule covers."""
    out, missing = {}, []
    for label, name in names.items():
        low = name.lower().replace("—", "-").replace("–", "-")
        for pat, tag in RULES:
            if pat in low:
                out[int(label)] = tag
                break
        else:
            missing.append(name)
    if missing:
        raise ValueError(f"no tissue rule for MIDA structures: {missing}")
    return out


def parse_label_file(fn):
    """{label: name} from MIDA's label text file: lines of 'index name' or 'name index'
    (tab or space separated; the download's exact layout is accepted either way)."""
    names = {}
    for line in open(fn, encoding="utf-8", errors="replace"):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        m = re.match(r"^(\d+)\s+(.+?)$", line) or re.match(r"^(.+?)\s+(\d+)$", line)
        if not m:
            continue
        a, b = m.groups()
        label, name = (int(a), b) if a.isdigit() else (int(b), a)
        names[label] = name.strip().strip('",')
    if len(names) < 50:
        raise ValueError(f"{fn}: only {len(names)} labels parsed; not a MIDA label file?")
    return names


def find_files(mida_dir):
    """(nifti, labels txt or None) in a MIDA download directory (searched recursively);
    the download has MIDA_v1_voxels/MIDA_v1.nii and MIDA_v1.txt."""
    nii = sorted(glob.glob(os.path.join(mida_dir, "**", "*.nii*"), recursive=True))
    txt = sorted(glob.glob(os.path.join(mida_dir, "**", "*.txt"), recursive=True))
    if not nii:
        raise FileNotFoundError(f"need a .nii(.gz) label volume under {mida_dir}")
    nii = [f for f in nii if "mida" in os.path.basename(f).lower()] or nii
    txt = [f for f in txt if "mida" in os.path.basename(f).lower()] or txt
    return nii[0], (txt[0] if txt else None)


def load(mida_dir):
    """(labels uint16 (nx, ny, nz), affine mm, {label: name}) from the download."""
    import nibabel as nib

    fn_nii, fn_txt = find_files(mida_dir)
    img = nib.load(fn_nii)
    labels = np.asarray(img.dataobj).astype(np.uint16)
    names = dict(MIDA_V1_LABELS)
    if fn_txt:
        try:
            names = parse_label_file(fn_txt)
        except ValueError:
            pass  # unknown layout: the embedded v1.0 table applies to MIDA_v1.nii
    return labels, img.affine.astype(float), names


def relabel(labels, mapping):
    lut = np.zeros(int(labels.max()) + 1, np.uint16)
    for k, v in mapping.items():
        if k < len(lut):
            lut[k] = v
    return lut[labels]


def downsample(labels, affine, factor):
    """Nearest-neighbour decimation by an integer factor (0.5 mm -> 1 mm with factor 2)."""
    if factor == 1:
        return labels, affine
    out = labels[::factor, ::factor, ::factor]
    a = affine.copy()
    a[:3, :3] *= factor
    return out, a


def crop(labels, affine, z_min_mm=None):
    """Drop slices below z_min_mm (world), e.g. the neck, to keep the mesh small."""
    if z_min_mm is None:
        return labels, affine
    k = np.arange(labels.shape[2])
    z = (affine @ np.stack([np.zeros_like(k), np.zeros_like(k), k, np.ones_like(k)]))[2]
    keep = z >= z_min_mm
    k0 = int(np.argmax(keep))
    out = labels[:, :, k0:]
    a = affine.copy()
    a[:3, 3] = (affine @ np.array([0, 0, k0, 1]))[:3]
    return out, a


def prepare(mida_dir, workdir, downsample_factor=2, z_min_mm=None, threads=4, elem_sizes=None):
    """MIDA -> SimNIBS .msh (cached). Returns (msh path, {tag: name}, mapping)."""
    labels, affine, names = load(mida_dir)
    mapping = classify(names)
    vol = relabel(labels, mapping)
    vol, affine = crop(vol, affine, z_min_mm)
    vol, affine = downsample(vol, affine, downsample_factor)
    fn = labelvol.mesh_from_labels(vol, affine, workdir, "mida", elem_sizes=elem_sizes,
                                   extra_tags=sorted(EXTRA_COND), threads=threads,
                                   params={"downsample": downsample_factor, "z_min": z_min_mm})
    return fn, EXTRA_NAMES, mapping
