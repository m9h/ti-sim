"""Electrode-placement sensitivity of a TI montage (Yamada et al. 2026, Phys Med Biol):
each electrode of Wessel 2023's striatal montage (F3-F4 / TP7-TP8, 2 mA, 9.8 mm) is moved
+-5 mm along two tangential directions on the scalp, one electrode at a time (16 runs plus
the reference), and the ROI median envelope and hotspot position are compared. Yamada
report 0.4-15 % intensity sensitivity and 0.6-6.9 mm hotspot displacement for +-10 mm
grids; this is the +-5 mm ring of that grid on ernie.

    TI_SIM_M2M=data/m2m_ernie TI_SIM_OUT=results/ernie_placement_sensitivity.json uv run examples/ernie_placement_sensitivity.py
"""

import os
import time

import numpy as np
from simnibs import mni2subject_coords

from ti_sim.electrodes import SkinSurface
from ti_sim.head import GM, WM, TwoPairSimulation, eeg_positions
from ti_sim.report import save_results
from ti_sim.tissues import table_from_env

M2M = os.environ["TI_SIM_M2M"]
WORK = os.environ.get("TI_SIM_WORK", os.path.join(os.path.dirname(M2M), "ti-sim-work"))
SHIFT = float(os.environ.get("TI_SIM_SHIFT_MM", 5.0))
COND_LABEL, COND = table_from_env()
MONTAGE = {"pairs": [("F3", "F4"), ("TP7", "TP8")], "currents": (2e-3, 2e-3), "radius": 9.8}
cap = eeg_positions(M2M)
from simnibs import mesh_io  # noqa: E402

skin = SkinSurface.from_simnibs(mesh_io.read_msh(os.path.join(M2M, "ernie.msh")))
targets = {"L putamen": mni2subject_coords((-25, 2, 0), M2M), "R putamen": mni2subject_coords((25, 2, 0), M2M)}


def tangents(p):
    """Two unit vectors tangent to the scalp at p (normal from the local skin patch)."""
    near = np.linalg.norm(skin.p - p, axis=1) < 15.0
    q = skin.p[near] - skin.p[near].mean(axis=0)
    n = np.linalg.svd(q)[2][-1]
    t1 = np.cross(n, [0, 0, 1.0])
    if np.linalg.norm(t1) < 0.3:
        t1 = np.cross(n, [0, 1.0, 0])
    t1 /= np.linalg.norm(t1)
    return t1, np.cross(n, t1)


def run(positions, label):
    t0 = time.time()
    pairs = [(positions[a], positions[b]) for a, b in MONTAGE["pairs"]]
    sim = TwoPairSimulation(pairs, MONTAGE["currents"], WORK, m2m_dir=M2M, radius=MONTAGE["radius"],
                            cond=COND, cond_label=COND_LABEL)
    brain = sim.in_tags((GM, WM))
    row = {"seconds": time.time() - t0, "positions": {k: list(map(float, v)) for k, v in positions.items()}}
    for name, c in targets.items():
        s = sim.stats(sim.in_sphere(c, 8.0, tags=(GM, WM)))
        row[name] = {"median": s["am"]["median"], "p99": s["am"]["p99"]}
    i = np.argmax(np.where(brain, sim.am_max, -1))
    row["hotspot_mm"] = list(map(float, sim.centres[i]))
    row["brain_max"] = float(sim.am_max[i])
    row["gm_median"] = sim.stats(sim.in_tags(GM))["am"]["median"]
    print(f"{label:12s} L putamen {row['L putamen']['median']:.3f} R {row['R putamen']['median']:.3f} V/m, "
          f"hotspot {np.round(row['hotspot_mm'], 1)} ({row['seconds']:.0f} s)")
    return row


names = [e for pair in MONTAGE["pairs"] for e in pair]
base = {n: skin.nearest(cap[n]) for n in names}
results = {"montage": MONTAGE, "shift_mm": SHIFT, "conductivities": COND_LABEL, "reference": run(base, "reference"), "runs": {}}
ref = results["reference"]
for n in names:
    t1, t2 = tangents(base[n])
    for tag, d in (("+t1", t1), ("-t1", -t1), ("+t2", t2), ("-t2", -t2)):
        pos = dict(base)
        pos[n] = skin.nearest(base[n] + SHIFT * d)
        row = run(pos, f"{n} {tag}")
        for roi in targets:
            row[roi]["change_pct"] = 100 * (row[roi]["median"] / ref[roi]["median"] - 1)
        row["hotspot_shift_mm"] = float(np.linalg.norm(np.array(row["hotspot_mm"]) - ref["hotspot_mm"]))
        results["runs"][f"{n} {tag}"] = row
summary = {}
for n in names:
    rows = [results["runs"][f"{n} {t}"] for t in ("+t1", "-t1", "+t2", "-t2")]
    summary[n] = {"intensity_sensitivity_pct": {roi: float(np.max([abs(r[roi]["change_pct"]) for r in rows])) for roi in targets},
                  "hotspot_shift_max_mm": float(np.max([r["hotspot_shift_mm"] for r in rows]))}
    print(f"{n}: max |change| L {summary[n]['intensity_sensitivity_pct']['L putamen']:.1f} % R "
          f"{summary[n]['intensity_sensitivity_pct']['R putamen']:.1f} %, hotspot shift up to {summary[n]['hotspot_shift_max_mm']:.1f} mm")
results["summary"] = summary
save_results(results)
