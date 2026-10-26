"""TI on the SimNIBS example head (ernie, low-res charm model): two published montages.

1. SimNIBS's own TI example, F5-P5 and F6-P6 at 1 mA each (20 mm electrodes).
2. Wessel et al. 2023 (Nat Neurosci 26:2005), striatum: F3-F4 and TP7-TP8 at 2 mA each,
   9.8 mm electrodes. They report, on MIDA, E_AM 0.26 +- 0.04 V/m (median 0.27, p99 0.35)
   in the bilateral striatum and 0.20 V/m in the precentral hand knob.

ROIs are 8 mm spheres at MNI targets mapped into subject space (hippocampus +-28,-20,-14;
putamen +-25,2,0; M1 hand +-37,-21,58), restricted to GM+WM, plus whole GM. Reference
numbers: Rampersad 2019 hippocampus 0.064 V/m (standard montage, 1 mA per pair) to
0.24 V/m (optimised, 2 mA total); Violante 2023 hippocampus 0.26 V/m at 1 mA per pair on
MIDA with cortex 0.10-0.18 V/m. Different heads, electrodes and segmentation, so the
check is "same order of magnitude and the deep/cortex ratio", not equality.

    TI_SIM_M2M=data/m2m_ernie TI_SIM_OUT=results/ernie_ti.json uv run examples/ernie_ti.py
"""

import os
import time

import numpy as np

from ti_sim.head import GM, WM, TwoPairSimulation
from ti_sim.report import save_results
from ti_sim.tissues import table_from_env

M2M = os.environ["TI_SIM_M2M"]
WORK = os.environ.get("TI_SIM_WORK", os.path.join(os.path.dirname(M2M), "ti-sim-work"))
COND_LABEL, COND = table_from_env()
MNI_TARGETS = {"L hippocampus": (-28, -20, -14), "R hippocampus": (28, -20, -14),
               "L putamen": (-25, 2, 0), "R putamen": (25, 2, 0),
               "L M1 hand": (-37, -21, 58), "R M1 hand": (37, -21, 58)}
ROI_RADIUS = 8.0
MONTAGES = {
    "simnibs_example": {"pairs": [("F5", "P5"), ("F6", "P6")], "currents": (1e-3, 1e-3), "radius": 20.0,
                        "reference": "SimNIBS TI.py example"},
    "wessel2023_striatum": {"pairs": [("F3", "F4"), ("TP7", "TP8")], "currents": (2e-3, 2e-3), "radius": 9.8,
                            "reference": "Wessel 2023: striatum 0.26 V/m (median 0.27, p99 0.35), M1 hand 0.20 V/m, on MIDA"},
}

from simnibs import mni2subject_coords  # noqa: E402

targets = {k: mni2subject_coords(v, M2M) for k, v in MNI_TARGETS.items()}
results = {"m2m": os.path.basename(M2M), "conductivities": COND_LABEL, "roi_radius_mm": ROI_RADIUS,
           "targets_subject_mm": {k: list(map(float, v)) for k, v in targets.items()}, "montages": {}}
for name, m in MONTAGES.items():
    t0 = time.time()
    sim = TwoPairSimulation(m["pairs"], m["currents"], WORK, m2m_dir=M2M, radius=m["radius"],
                            cond=COND, cond_label=COND_LABEL)
    row = {"pairs": m["pairs"], "currents_mA": [1e3 * c for c in m["currents"]], "radius_mm": m["radius"],
           "reference": m["reference"], "seconds": time.time() - t0, "rois": {}}
    gm = sim.in_tags(GM)
    row["rois"]["GM"] = sim.stats(gm)
    row["rois"]["brain"] = sim.stats(sim.in_tags((GM, WM)))
    for k, c in targets.items():
        row["rois"][k] = sim.stats(sim.in_sphere(c, ROI_RADIUS, tags=(GM, WM)))
    for thr in (0.1, 0.2, 0.25):
        row[f"focality_{thr}"] = sim.focality(thr)
    deep = np.mean([row["rois"][k]["am"]["median"] for k in ("L putamen", "R putamen", "L hippocampus", "R hippocampus")])
    row["deep_to_gm_median_ratio"] = deep / row["rois"]["GM"]["am"]["median"]
    row["deep_to_gm_p99_ratio"] = deep / row["rois"]["GM"]["am"]["p99"]
    row["calibration"] = sim.calibration
    results["montages"][name] = row
    print(f"{name}: {row['seconds']:.0f} s; GM median {row['rois']['GM']['am']['median']:.3f} "
          f"p99 {row['rois']['GM']['am']['p99']:.3f} max {row['rois']['GM']['am']['max']:.3f} V/m; "
          f"tACS GM p99 {row['rois']['GM']['tacs']['p99']:.3f}")
    for k in targets:
        s = row["rois"][k]
        print(f"  {k:14s} E_AM median {s['am']['median']:.3f} p99 {s['am']['p99']:.3f} | "
              f"|E1| {s['e1']['median']:.3f} |E2| {s['e2']['median']:.3f} tACS {s['tacs']['median']:.3f} V/m "
              f"({s['n_elements']} tets)")
    print(f"  deep/GM median ratio {row['deep_to_gm_median_ratio']:.2f}; brain > 0.2 V/m: "
          f"TI {row['focality_0.2']['am_mm3'] / 1e3:.1f} cm3, tACS {row['focality_0.2']['tacs_mm3'] / 1e3:.1f} cm3")
    out = os.environ.get("TI_SIM_OUT")
    if out:
        sim.write_msh(os.path.join(os.path.dirname(out), f"ernie_{name}.msh"))
save_results(results)
