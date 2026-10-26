"""TI on the MIDA head: the Violante 2023 and Wessel 2023 montages on the model they used.

Needs your own MIDA copy (TI_SIM_MIDA, see scripts/check-mida.sh). The model is cropped
below z = TI_SIM_MIDA_ZMIN (mm, default: none), downsampled by TI_SIM_MIDA_DS (default 2,
1 mm) and meshed once (cached). MIDA has no EEG cap: give the four fiducials in MIDA world
coordinates (TI_SIM_MIDA_FIDUCIALS, a JSON file {"Nz": [x, y, z], "Iz": ..., "LPA": ...,
"RPA": ...}) and `ti_sim.electrodes.ten_ten` builds the 10-10 positions geodesically on the
meshed skin (6 mm mean against a template cap on ernie); or give positions directly
(TI_SIM_MIDA_ELECTRODES, {name: [x, y, z]}), which also override the constructed ones
(Violante's e1-e4 are custom positions). Without either the script only builds the mesh
and reports the tissue volumes, which is already a check of the mapping.

Reference: Violante 2023, 1 mA per pair, 2 cm electrodes: hippocampus E_AM 0.26 +- 0.04 V/m,
cortex 0.10-0.18 V/m. Wessel 2023, 2 mA per pair, 9.8 mm: striatum 0.26 +- 0.04 V/m.
"""

import json
import os

import numpy as np

from ti_sim import mida
from ti_sim.head import GM, WM, TwoPairSimulation
from ti_sim.report import save_results
from ti_sim.tissues import table_from_env

MIDA_DIR = os.environ["TI_SIM_MIDA"]
WORK = os.environ.get("TI_SIM_WORK", "data/ti-sim-work")
fn, extra_names, mapping = mida.prepare(MIDA_DIR, WORK, downsample_factor=int(os.environ.get("TI_SIM_MIDA_DS", 2)),
                                        z_min_mm=float(os.environ["TI_SIM_MIDA_ZMIN"]) if "TI_SIM_MIDA_ZMIN" in os.environ else None)
from simnibs import mesh_io  # noqa: E402

m = mesh_io.read_msh(fn)
tets = m.elm.elm_type == 4
vol = m.elements_volumes_and_areas().value[tets]
tags = m.elm.tag1[tets]
results = {"mesh": fn, "n_tets": int(tets.sum()), "mapping_size": len(mapping),
           "volumes_cm3": {str(t): float(vol[tags == t].sum() / 1e3) for t in np.unique(tags)}}
print(f"MIDA mesh {fn}: {results['n_tets']} tets; hippocampus {results['volumes_cm3'].get(str(mida.HIPPOCAMPUS), 0):.1f} cm3")
pos = {}
fn_fid = os.environ.get("TI_SIM_MIDA_FIDUCIALS")
if fn_fid:
    from ti_sim.electrodes import SkinSurface, ten_ten

    fid = {k: list(map(float, v)) for k, v in json.load(open(fn_fid)).items()}
    skin = SkinSurface.from_simnibs(m)
    pos = {k: list(map(float, v)) for k, v in ten_ten(skin.p, skin.t, fid).items()}
    results["ten_ten_from_fiducials"] = {k: pos[k] for k in ("F3", "F4", "TP7", "TP8", "Cz") if k in pos}
    print(f"10-10 positions from fiducials: {len(pos)} names; F3 {np.round(pos['F3'], 1)} TP7 {np.round(pos['TP7'], 1)}")
fn_el = os.environ.get("TI_SIM_MIDA_ELECTRODES")
if fn_el:
    pos.update({k: list(map(float, v)) for k, v in json.load(open(fn_el)).items()})
if pos:
    cond_label, cond = table_from_env()
    montages = {"violante2023_hippocampus": {"pairs": [("e1", "e3"), ("e2", "e4")], "currents": (1e-3, 1e-3), "radius": 10.0},
                "wessel2023_striatum": {"pairs": [("F3", "F4"), ("TP7", "TP8")], "currents": (2e-3, 2e-3), "radius": 9.8}}
    results["montages"] = {}
    for name, mt in montages.items():
        if not all(e in pos for p in mt["pairs"] for e in p):
            print(f"{name}: electrode positions missing, skipped")
            continue
        pairs = [(pos[a], pos[b]) for a, b in mt["pairs"]]
        sim = TwoPairSimulation(pairs, mt["currents"], WORK, fnamehead=fn, radius=mt["radius"], cond=cond,
                                cond_label=cond_label, extra_cond=mida.EXTRA_COND)
        row = {"rois": {"GM": sim.stats(sim.in_tags(GM)), "brain": sim.stats(sim.in_tags((GM, WM)))}}
        for tag, nm in extra_names.items():
            if sim.in_tags(tag).any():
                row["rois"][nm] = sim.stats(sim.in_tags(tag))
        results["montages"][name] = row
        print(name, {k: round(v["am"]["median"], 3) for k, v in row["rois"].items()})
save_results(results)
