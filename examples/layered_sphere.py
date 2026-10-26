"""Replication target 2: Grossman 2017 figure S2, the 4-layer sphere, through the same
label-volume -> create_mesh -> SimNIBS route MIDA takes.

The sphere is scaled to a human head (R = 85 mm) with the paper's proportions: scalp
0.05 R (0.333 S/m), skull 0.085 R (0.0083), CSF 0.023 R (1.79), brain 0.83 R (0.333). Two
pairs sit on the equator as legs of a trapezoid like the phantom's layout C, and the
current ratio is swept 1:1, 1:2.5, 1:4 at a 2 mA sum. What the paper shows and we check:
the envelope peaks between the pairs, deeper than the unmodulated field under the
electrodes, and the peak moves towards the weaker pair as the ratio changes.

    TI_SIM_OUT=results/layered_sphere.json uv run examples/layered_sphere.py
"""

import os

import numpy as np

from ti_sim import labelvol
from ti_sim.head import GM, TwoPairSimulation
from ti_sim.report import save_results
from ti_sim.tissues import GROSSMAN_SPHERE

R = float(os.environ.get("TI_SIM_SPHERE_R", 85.0))
SPACING = float(os.environ.get("TI_SIM_SPHERE_VOXEL", 1.0))
WORK = os.environ.get("TI_SIM_WORK", "data/ti-sim-work")
radii = (R, R * (1 - 0.05), R * (1 - 0.135), R * (1 - 0.158))
labels, aff = labelvol.sphere_labels(radii, (5, 4, 3, 2), spacing_mm=SPACING)
fn = labelvol.mesh_from_labels(labels, aff, WORK, f"sphere{R:g}",
                               elem_sizes={"standard": {"range": [1.5, 6.0], "slope": 1.0}})
cond = {1: 0.333, 2: GROSSMAN_SPHERE["brain"], 3: GROSSMAN_SPHERE["csf"], 4: GROSSMAN_SPHERE["skull"],
        5: GROSSMAN_SPHERE["scalp"], 6: 0.5, 7: 0.0083, 8: 0.0083, 9: 0.6, 10: 0.16, 11: 0.88, 12: 0.078}


def on_sphere(deg):
    a = np.radians(deg)
    return [R * np.cos(a), R * np.sin(a), 0.0]


# layout C of the phantom on the equator: pair 1 = left leg, pair 2 = right leg
pairs = [(on_sphere(180 + 11.5), on_sphere(180 - 11.5)), (on_sphere(-11.5), on_sphere(11.5))]
results = {"R_mm": R, "voxel_mm": SPACING, "mesh": fn, "ratios": {}}
for label, (i1, i2) in {"1:1": (1e-3, 1e-3), "1:2.5": (2e-3 / 3.5, 2e-3 * 2.5 / 3.5), "1:4": (2e-3 / 5, 2e-3 * 4 / 5)}.items():
    sim = TwoPairSimulation(pairs, (i1, i2), WORK, fnamehead=fn, radius=10.0, cond=cond, cond_label="grossman-sphere")
    gm = sim.in_tags(GM)
    c = sim.centres[gm]
    am = sim.am_max[gm]
    # profile along the x axis (through both pairs) in a 6 mm tube
    tube = (np.abs(c[:, 1]) < 3) & (np.abs(c[:, 2]) < 3)
    order = np.argsort(c[tube, 0])
    x, prof = c[tube, 0][order], am[tube][order]
    peak_x = float(x[np.argmax(prof)])
    e1 = np.linalg.norm(sim.e1[gm], axis=1)
    row = {"I1_mA": 1e3 * i1, "I2_mA": 1e3 * i2, "am_max_brain": float(am.max()),
           "am_peak_x_over_R": peak_x / R, "am_at_centre": float(np.median(am[np.linalg.norm(c, axis=1) < 8])),
           "e1_max_brain": float(e1.max()), "e1_at_centre": float(np.median(e1[np.linalg.norm(c, axis=1) < 8])),
           "am_max_depth_mm": float(R - np.linalg.norm(c[np.argmax(am)])),
           "e1_max_depth_mm": float(R - np.linalg.norm(c[np.argmax(e1)])),
           "nodes": int(sim.mesh.nodes.nr), "tets": int(gm.sum())}
    results["ratios"][label] = row
    print(f"{label}: E_AM max {row['am_max_brain']:.3f} V/m at depth {row['am_max_depth_mm']:.0f} mm, "
          f"x-profile peak at {row['am_peak_x_over_R']:+.2f} R; centre {row['am_at_centre']:.3f} V/m; "
          f"|E1| max {row['e1_max_brain']:.3f} at depth {row['e1_max_depth_mm']:.0f} mm")
save_results(results)
