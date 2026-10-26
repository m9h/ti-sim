"""Replication target 1: Grossman et al. 2017, figure 2, the saline-disc phantom.

Five electrode layouts (trapezoid bases a, b in units of R, current ratio I1:I2):
A 0.39/1.96 1:1, B 1.11/1.96 1:1, C 1.96/1.96 1:1, D 1.96/1.96 1:2.5, E 1.96/1.96 1:4;
1 mA per pair in A-C, a 2 mA sum in D-E. The paper tabulates, per layout, the peak of the
envelope along x and along y (simulated / measured, V/m), its position (x, y)/R and the
half-widths at half maximum (x, y)/R. We compute the same from the 2-D FEM and save
them beside the paper's simulated numbers (PAPER), with ratios for the pipeline.

    TI_SIM_OUT=results/grossman_phantom.json uv run examples/grossman_phantom.py
"""

import numpy as np

from ti_sim.phantom import DiscTI, R_DISC, trapezoid_electrodes
from ti_sim.report import figure_path, save_results

LAYOUTS = {  # name: (a, b, currents (A))
    "A": (0.39, 1.96, (1e-3, 1e-3)),
    "B": (1.11, 1.96, (1e-3, 1e-3)),
    "C": (1.96, 1.96, (1e-3, 1e-3)),
    "D": (1.96, 1.96, (2e-3 / 3.5, 2e-3 * 2.5 / 3.5)),
    "E": (1.96, 1.96, (2e-3 / 5, 2e-3 * 4 / 5)),
}
# Grossman 2017 figure 2 table: simulated (measured) peak |E_AM,x| and |E_AM,y| in V/m,
# peak position (x, y)/R and HWHM (x, y)/R.
PAPER = {
    "A": {"x": (3.7, 1.0, (0, -1), (0.18, 0.1)), "y": (2.5, 1.1, (0, -0.77), (0.18, 0.37))},
    "B": {"x": (0.9, 0.75, (0, -0.87), (0.33, 0.25)), "y": (0.9, 0.7, (0, -0.37), (0.4, 0.52))},
    "C": {"x": (0.13, 0.15, (0, -0.59), (0.37, 0.42)), "y": (0.21, 0.23, (0, 0), (0.42, 0.49))},
    "D": {"x": (0.13, 0.14, (-0.2, -0.59), (0.37, 0.42)), "y": (0.2, 0.22, (-0.2, 0), (0.38, 0.4))},
    "E": {"x": (0.12, 0.13, (-0.35, 0), (0.36, 0.42)), "y": (0.19, 0.2, (-0.35, 0), (0.37, 0.34))},
}
EXCLUDE = 1.5e-3  # m from the rim: the wires' own 1/r field

results = {"radius_m": R_DISC, "exclude_m": EXCLUDE, "layouts": {}}
maps = {}
for name, (a, b, currents) in LAYOUTS.items():
    ti = DiscTI(trapezoid_electrodes(a, b), currents, nrefs=7, refine_rounds=2)
    em = ti.envelope_map(n=201)
    maps[name] = em
    row = {"a": a, "b": b, "currents_mA": [1e3 * c for c in currents], "nodes": ti.mesh.p.shape[1]}
    for axis in ("x", "y"):
        v, pos, hw = em.peak(f"am_{axis}", exclude=EXCLUDE)
        sim, meas, ppos, phw = PAPER[name][axis]
        row[axis] = {"peak_V_per_m": v, "pos_R": pos, "hwhm_R": hw, "paper_sim": sim,
                     "paper_meas": meas, "paper_pos_R": ppos, "paper_hwhm_R": phw,
                     "ratio_to_sim": v / sim, "ratio_to_meas": v / meas}
        print(f"{name} E_AM,{axis}: {v:6.3f} V/m at ({pos[0]:+.2f}, {pos[1]:+.2f})R "
              f"hwhm ({hw[0]:.2f}, {hw[1]:.2f})R | paper sim {sim} meas {meas} "
              f"at {ppos} hwhm {phw}")
    v, pos, hw = em.peak("am_max", exclude=EXCLUDE)
    row["free"] = {"peak_V_per_m": v, "pos_R": pos, "hwhm_R": hw}
    results["layouts"][name] = row

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 5, figsize=(16, 6.5), constrained_layout=True)
    for j, name in enumerate(LAYOUTS):
        em = maps[name]
        for i, axis in enumerate(("x", "y")):
            ax = axes[i, j]
            f = getattr(em, f"am_{axis}")
            vmax = PAPER[name][axis][0] * 1.5
            ax.imshow(f, origin="lower", extent=[-1, 1, -1, 1], vmin=0, vmax=vmax, cmap="magma")
            for pair in trapezoid_electrodes(LAYOUTS[name][0], LAYOUTS[name][1], 1.0):
                ax.plot([p[0] for p in pair], [p[1] for p in pair], "wo-", ms=4, lw=0.8)
            ax.set_title(f"{name}: |E_AM,{axis}|  (clip {vmax:.2g} V/m)", fontsize=9)
            ax.set_xticks([]), ax.set_yticks([])
    fig.savefig(figure_path("grossman_phantom.png"), dpi=110)
except ImportError:
    pass

save_results(results)
