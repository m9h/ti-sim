"""Montage search for the left hippocampus on ernie: lead field over the 10-10 cap, then
every two-pair montage from the electrodes within TI_SIM_SEARCH_RADIUS mm (default 90) of
the target, at current ratios 1:1, 1:2, 1:3 with a 2 mA sum (Rampersad 2019: 2 mA total;
Lee 2020/2026: exhaustive search restricted to target-specific candidates).

The lead field is one SimNIBS solve per electrode (75 on the 10-10 cap; 6 min on the
low-res ernie with 6 threads) and is cached. Rampersad's optimum for the left hippocampus was 0.24 V/m
(free direction) with 2 mA total.

    TI_SIM_M2M=data/m2m_ernie TI_SIM_OUT=results/ernie_montage_search.json uv run examples/ernie_montage_search.py
"""

import os
import time

from simnibs import mni2subject_coords

from ti_sim.head import eeg_positions
from ti_sim.montage import Leadfield, compute_leadfield, near_target, search_two_pairs
from ti_sim.report import save_results

M2M = os.environ["TI_SIM_M2M"]
WORK = os.environ.get("TI_SIM_WORK", os.path.join(os.path.dirname(M2M), "ti-sim-work"))
RADIUS = float(os.environ.get("TI_SIM_SEARCH_RADIUS", 90.0))
target = mni2subject_coords((-28, -20, -14), M2M)
t0 = time.time()
fn = compute_leadfield(M2M, WORK, radius=10.0)
lf = Leadfield(fn)
t_lf = time.time() - t0
cands = near_target(eeg_positions(M2M), target, RADIUS)
cands = [c for c in cands if c in lf.idx and lf.idx[c] is not None]
roi = lf.sphere(target, 8.0)
t0 = time.time()
res = search_two_pairs(lf, roi, electrodes=cands, ratios=(1.0, 2.0, 3.0), total_current=2e-3, top=10, threshold=0.25,
                       shortlist=int(os.environ.get("TI_SIM_SHORTLIST", 300)), min_roi_fraction=0.7)
res.update({"leadfield": fn, "leadfield_seconds": t_lf, "search_seconds": time.time() - t0,
            "target_subject_mm": list(map(float, target)), "candidates": cands,
            "reference": "Rampersad 2019 L hippocampus optimum 0.24 V/m at 2 mA total"})
print(f"lead field {t_lf:.0f} s; {res['n_montages']} montages from {res['n_electrodes']} electrodes in {res['search_seconds']:.0f} s")
for label, key in (("by ROI intensity", "best"), ("focal (>= 70 % of the best ROI median, least brain > 0.25 V/m outside)", "focal")):
    print(label)
    for r in res.get(key, [])[:5]:
        print(f"  {r['pair1']} {r['I1_mA']:.2f} mA + {r['pair2']} {r['I2_mA']:.2f} mA: ROI median {r['roi_median']:.3f} V/m, "
              f"outside p99 {r['brain_p99_outside_roi']:.3f}, >0.25 V/m: roi {r['above_threshold_mm3']['roi']:.0f} / outside {r['above_threshold_mm3']['outside']:.0f} mm3")
save_results(res)
