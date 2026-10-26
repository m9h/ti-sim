"""Conductivities (S/m) used by the TI literature, so a run can say which table it used.

TI carriers sit at 1-10 kHz, where tissue conductivity differs little from DC for the
head tissues that matter (Gabriel et al. 1996: grey matter 0.10 at 10 Hz vs 0.11 at 2 kHz
in their 4-Cole-Cole fit; skull and CSF flat). Every published TI head model uses
quasi-static DC-like values and solves one Laplace problem per electrode pair, so the two
carriers see the same volume conductor and the envelope is purely geometric. We do the
same and keep the tables explicit.

Tags follow SimNIBS ElementTags (1 WM, 2 GM, 3 CSF, 4 bone, 5 scalp, 6 eyes, 7 compact
bone, 8 spongy bone, 9 blood, 10 muscle, 11 cartilage, 12 fat) so a table can be handed to
a SimNIBS session (`apply_to_tdcslist`) or to a scikit-fem per-element array.
"""

import os

# SimNIBS 4 defaults (Thielscher/Saturnino; the values charm heads are solved with).
SIMNIBS = {1: 0.126, 2: 0.275, 3: 1.654, 4: 0.010, 5: 0.465, 6: 0.500, 7: 0.008, 8: 0.025,
           9: 0.600, 10: 0.160, 11: 0.880, 12: 0.078}

# Rampersad et al. 2019, NeuroImage 202:116124, table 1 (isotropic values; they also used
# DTI tensors in GM/WM). Brainstem 0.201 and muscle 0.4 have no SimNIBS tag; muscle -> 10.
RAMPERSAD = {1: 0.126, 2: 0.276, 3: 1.65, 4: 0.007, 5: 0.465, 6: 1.5, 7: 0.007, 8: 0.025,
             9: 0.6, 10: 0.4, 11: 0.88, 12: 0.078}

# Grossman et al. 2017, Cell 169:1029, STAR Methods: saline phantom and 4-layer sphere
# (scalp 0.333, skull 0.0083, CSF 1.79, brain 0.333). Keys are layer names, not tags.
GROSSMAN_SPHERE = {"scalp": 0.333, "skull": 0.0083, "csf": 1.79, "brain": 0.333}
GROSSMAN_SALINE = 0.333

TABLES = {"simnibs": SIMNIBS, "rampersad": RAMPERSAD}

NAMES = {1: "WM", 2: "GM", 3: "CSF", 4: "bone", 5: "scalp", 6: "eyes", 7: "compact bone",
         8: "spongy bone", 9: "blood", 10: "muscle", 11: "cartilage", 12: "fat"}


def table_from_env(default="simnibs"):
    """(name, {tag: sigma}) from TI_SIM_COND: a table name, or 'tag=value,...' overrides on
    the default table (e.g. 'simnibs,4=0.0085')."""
    spec = os.environ.get("TI_SIM_COND", default)
    parts = [p.strip() for p in spec.split(",") if p.strip()]
    name = parts[0] if parts and "=" not in parts[0] else default
    table = dict(TABLES[name])
    for p in parts:
        if "=" in p:
            k, v = p.split("=")
            table[int(k)] = float(v)
    label = name if len(parts) <= 1 else spec
    return label, table


def apply_to_tdcslist(tdcs, table):
    """Set the conductivities of a SimNIBS TDCSLIST from a {tag: sigma} table."""
    for tag, sigma in table.items():
        tdcs.cond[tag - 1].value = float(sigma)
        tdcs.cond[tag - 1].name = NAMES.get(tag, f"tag{tag}")
    return tdcs
