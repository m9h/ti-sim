"""Firing threshold of a Hodgkin-Huxley point neuron vs carrier frequency under TI (ti_sim.hh),
10 Hz beat, equal pairs (E_AM,max = 2 A), warm kinetics (Q10 = 3 from 6.3 C, x3), against a
plain 10 Hz field and against Esmaeilpour 2021's slice thresholds (5 / 60 / 80 V/m at
100 Hz / 1 kHz / 2 kHz carriers, amplitude-modulated at 5 Hz). Also how far the montages
simulated on ernie are below that threshold (Wessel montage at 2 mA: ~0.4 V/m in the
striatum): the point TI makes about direct activation at human currents.

    TI_SIM_OUT=results/hh_thresholds.json uv run examples/hh_thresholds.py
"""

import os
import time

import numpy as np

from ti_sim import hh
from ti_sim.membrane import ESMAEILPOUR_2021
from ti_sim.report import save_results

CARRIERS = [float(c) for c in os.environ.get("TI_SIM_HH_CARRIERS", "0,100,500,1000,2000,5000").split(",")]
WARM = float(os.environ.get("TI_SIM_HH_WARM", 3.0))
t0 = time.time()
thr = {}
for f in CARRIERS:
    a = hh.threshold(f, 10.0, lo=0.5, hi=3000.0, duration_ms=150, warm=WARM)
    thr[f] = {"amplitude_per_carrier_V_per_m": a, "E_AM_max_V_per_m": 2 * a if np.isfinite(a) else None}
    print(f"carrier {f:6.0f} Hz: threshold A = {a:8.2f} V/m per carrier (E_AM,max {2*a:8.2f})" if np.isfinite(a) else f"carrier {f:6.0f} Hz: no firing below 3000 V/m per carrier")
ref = thr[0.0]["E_AM_max_V_per_m"] if 0.0 in thr else None
results = {"warm": WARM, "beat_hz": 10.0, "polarisation_mm": 0.2, "thresholds": {str(int(k)): v for k, v in thr.items()},
           "esmaeilpour_2021_V_per_m": {str(k): v[0] for k, v in ESMAEILPOUR_2021.items()},
           "ernie_wessel_striatum_V_per_m": 0.39, "seconds": time.time() - t0}
if ref:
    for f, row in thr.items():
        if row["E_AM_max_V_per_m"]:
            row["ratio_to_plain_field"] = row["E_AM_max_V_per_m"] / ref
            row["ernie_striatum_fraction_of_threshold"] = 0.39 / row["E_AM_max_V_per_m"]
    print(f"plain 10 Hz field threshold {ref:.2f} V/m; the ernie Wessel montage (0.39 V/m) is "
          f"{100*0.39/ref:.1f} % of it, and {100*0.39/thr[2000.0]['E_AM_max_V_per_m']:.2f} % of the 2 kHz TI threshold" if thr.get(2000.0, {}).get("E_AM_max_V_per_m") else "")
save_results(results)
