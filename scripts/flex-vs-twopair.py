"""Does TesFlexOptimization's final real-electrode simulation need the discrete-current
correction? Read its two tDCS meshes, extract the electrode centres, and rerun the same
montage through TwoPairSimulation on the same head (which measures the exact currents).

    TI_SIM_M2M=data/m2m_ernie TI_SIM_FLEX_JSON=results/ernie_flex_ti.json TI_SIM_OUT=results/flex_vs_twopair.json \\
        uv run scripts/flex-vs-twopair.py
"""
import json, os, sys

import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
from simnibs import mesh_io
from ti_sim import envelope
from ti_sim.head import GM, WM, TwoPairSimulation, discrete_current

d = json.load(open(os.environ.get("TI_SIM_FLEX_JSON", "results/ernie_flex_ti.json")))
folder = d["output_folder"]; target = np.array(d["target_subject_mm"]); I = d["current_mA"] * 1e-3
out = {"flex_final_roi_mean": d["final_roi_mean_V_per_m"], "flex_final_roi_median": d["final_roi_median_V_per_m"]}
centres, factors, fields = [], [], []
for k in (0, 1):
    m = mesh_io.read_msh(os.path.join(folder, f"final_sim_{k}", "ernie_TDCS_1_scalar.msh"))
    print("fields in final sim", k, list(m.field.keys()))
    tri = m.elm.elm_type == 2
    tops = [t for t in np.unique(m.elm.tag1[tri]) if 2100 <= t < 2200]
    cs = []
    for t in tops:
        nodes = np.unique(m.elm.node_number_list[tri & (m.elm.tag1 == t)][:, :3]) - 1
        cs.append(m.nodes.node_coord[nodes].mean(0))
    centres.append(cs)
    if "J" in m.field and "v" in m.field:
        c = discrete_current(m, I); factors.append(c["factor"]); out[f"pair{k+1}_calibration"] = c
        print(f"final sim {k}: electrode tops {tops}, centres {np.round(cs, 1).tolist()}, SimNIBS bias {c['simnibs_bias_percent']:+.1f} %")
    else:
        factors.append(None); print(f"final sim {k}: no J/v field, cannot measure the current; centres {np.round(cs, 1).tolist()}")
# the optimiser's own head mesh: max_TI in the ROI, uncorrected
hm = mesh_io.read_msh(os.path.join(folder, "ernie_tes_flex_opt_head_mesh.msh"))
tets = hm.elm.elm_type == 4
am, roi = hm.field["max_TI"].value[tets], hm.field["ROI"].value[tets] > 0
vol = hm.elements_volumes_and_areas().value[tets]
print(f"flex head mesh: ROI mean {np.average(am[roi], weights=vol[roi]):.4f}, n ROI tets {roi.sum()}")
# same montage through TwoPairSimulation (its own SimNIBS runs, corrected to the discrete current)
pairs = [(centres[0][0].tolist(), centres[0][1].tolist()), (centres[1][0].tolist(), centres[1][1].tolist())]
sim = TwoPairSimulation(pairs, (I, I), os.environ.get("TI_SIM_WORK", "data/ti-sim-work"), m2m_dir=os.environ["TI_SIM_M2M"], radius=10.0, threads=8)
mask = sim.in_sphere(target, 8.0, tags=(GM, WM))
st = sim.stats(mask)
raw = envelope.max_envelope(sim.e1 / sim.calibration[0]["factor"], sim.e2 / sim.calibration[1]["factor"])
v = sim.volumes[mask]
out.update({"twopair_calibration": sim.calibration, "twopair_roi_mean_corrected": st["am"]["mean"],
            "twopair_roi_median_corrected": st["am"]["median"],
            "twopair_roi_mean_uncorrected": float(np.average(raw[mask], weights=v)), "pairs_mm": pairs})
print(json.dumps({k: v for k, v in out.items() if "calibration" not in k}, indent=1))
print("twopair calibration", [(round(c["factor"], 3), round(c["simnibs_bias_percent"], 1)) for c in sim.calibration])
json.dump(out, open(os.environ.get("TI_SIM_OUT", "results/flex_vs_twopair.json"), "w"), indent=1)
