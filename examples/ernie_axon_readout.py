"""A NEURON cable axon at a montage's target, driven by the montage's own E1(t), E2(t):
how far the published montages sit below the axon's firing threshold.

Wessel 2023's striatal montage on ernie (F3-F4 / TP7-TP8, 2 mA each) is simulated once
(`TwoPairSimulation`, cached), E1 and E2 are read at the left putamen target, and a 10 mm
HH cable (`ti_sim.axon`) is aligned with the envelope's best direction there. The field
E1 cos(2 pi f t) + E2 cos(2 pi (f + 10) t) is then scaled until the axon fires: the gain
needed is the margin below direct activation (Grossman's point that kHz TI at human currents
does not fire neurons by itself; Esmaeilpour 2021: ~80 V/m at 2 kHz in slices).

    TI_SIM_M2M=data/m2m_ernie TI_SIM_OUT=results/ernie_axon_readout.json uv run examples/ernie_axon_readout.py
"""

import os
import time

import numpy as np
from simnibs import mni2subject_coords

from ti_sim import axon, envelope
from ti_sim.head import TwoPairSimulation
from ti_sim.report import save_results
from ti_sim.tissues import table_from_env

M2M = os.environ["TI_SIM_M2M"]
WORK = os.environ.get("TI_SIM_WORK", os.path.join(os.path.dirname(M2M), "ti-sim-work"))
CARRIERS = [float(c) for c in os.environ.get("TI_SIM_AXON_CARRIERS", "0,1000,2000").split(",")]
BEAT = 10.0
cond_label, cond = table_from_env()
target = mni2subject_coords((-25, 2, 0), M2M)  # left putamen

sim = TwoPairSimulation([("F3", "F4"), ("TP7", "TP8")], (2e-3, 2e-3), WORK, m2m_dir=M2M, radius=9.8,
                        cond=cond, cond_label=cond_label)
e1, e2, am = (x[0] for x in sim.fields_at([target]))
direction = envelope.envelope_direction(e1[None], e2[None])[0]
a1, a2 = float(e1 @ direction), float(e2 @ direction)
print(f"putamen target {np.round(target, 1)}: |E1| {np.linalg.norm(e1):.3f}, |E2| {np.linalg.norm(e2):.3f}, "
      f"E_AM,max {am:.3f} V/m along {np.round(direction, 2)}; projections {a1:.3f}, {a2:.3f}")

ax = axon.Axon(length_mm=10.0, diam_um=2.0, nseg=101, axis=direction)
t = np.arange(0.0, 120.0 + 0.0125, 0.0125)
results = {"montage": "wessel2023_striatum", "target_subject_mm": list(map(float, target)),
           "E1": e1.tolist(), "E2": e2.tolist(), "E_AM_max_V_per_m": float(am), "direction": direction.tolist(),
           "axon": {"length_mm": 10.0, "diam_um": 2.0, "celsius": 6.3}, "beat_hz": BEAT, "carriers": {}}
t0 = time.time()
for f in CARRIERS:
    def field(gain, f=f):
        ts = t * 1e-3
        if f == 0:
            return gain * (a1 + a2) * np.cos(2 * np.pi * BEAT * ts)
        return gain * (a1 * np.cos(2 * np.pi * f * ts) + a2 * np.cos(2 * np.pi * (f + BEAT) * ts))

    lo, hi = 1.0, 1e4
    if ax.spikes(t, field(hi)) < 1:
        gain = float("nan")
    elif ax.spikes(t, field(lo)) >= 1:
        gain = lo
    else:
        while hi / lo > 1.05:
            mid = float(np.sqrt(lo * hi))
            if ax.spikes(t, field(mid)) >= 1:
                hi = mid
            else:
                lo = mid
        gain = hi
    row = {"gain_to_fire": gain, "fires_at_montage_current": bool(ax.spikes(t, field(1.0)) >= 1)}
    if np.isfinite(gain):
        row["E_AM_max_at_threshold_V_per_m"] = float(am * gain)
        row["montage_fraction_of_threshold_percent"] = 100.0 / gain
    results["carriers"][str(int(f))] = row
    print(f"carrier {f:5.0f} Hz: the montage field must be scaled x{gain:.1f} to fire the axon "
          f"(E_AM,max {am * gain:.1f} V/m); at 2 mA it is {100 / gain:.2f} % of threshold" if np.isfinite(gain)
          else f"carrier {f:5.0f} Hz: no firing up to x{hi:.0f}")
results["seconds"] = time.time() - t0
save_results(results)
