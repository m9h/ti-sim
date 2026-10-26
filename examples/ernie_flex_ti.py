"""Free electrode placement for TI with SimNIBS's TesFlexOptimization (Weise et al.): two
freely movable pairs, 10 mm electrodes, 1 mA each (Violante 2023's currents), maximising
the mean E_AM,max in an 8 mm sphere at the left hippocampus (GM+WM volume ROI), by
differential evolution with a FEM solve per candidate.

TI_SIM_FLEX_MAXITER (default 2) and TI_SIM_FLEX_POPSIZE (default 4) make a 2 min smoke
run that still ends with the final electrode simulation and its ROI value; 50 and 13
converge in ~11 min on the low-res ernie (2392 FEM evaluations), which `ti_sim.montage`'s
exhaustive search and Rampersad's 0.24 V/m can be compared with.

    TI_SIM_M2M=data/m2m_ernie TI_SIM_OUT=results/ernie_flex_ti.json uv run examples/ernie_flex_ti.py
"""

import glob
import json
import os
import time

import numpy as np
from simnibs import mni2subject_coords, opt_struct

from ti_sim.report import save_results

M2M = os.environ["TI_SIM_M2M"]
WORK = os.environ.get("TI_SIM_WORK", os.path.join(os.path.dirname(M2M), "ti-sim-work"))
MAXITER = int(os.environ.get("TI_SIM_FLEX_MAXITER", 2))
POPSIZE = int(os.environ.get("TI_SIM_FLEX_POPSIZE", 4))
CURRENT = float(os.environ.get("TI_SIM_FLEX_MA", 1.0)) * 1e-3
target = mni2subject_coords((-28, -20, -14), M2M)

opt = opt_struct.TesFlexOptimization()
opt.subpath = M2M
opt.output_folder = os.path.join(WORK, f"flex_ti_lhipp_it{MAXITER}_pop{POPSIZE}_{CURRENT*1e3:g}mA")
opt.goal = "mean"
opt.e_postproc = "max_TI"
opt.open_in_gmsh = False
opt.run_final_electrode_simulation = True
opt.optimizer_options = {"maxiter": MAXITER, "popsize": POPSIZE, "disp": True, "tol": 0.1}
for _ in range(2):
    layout = opt.add_electrode_layout("ElectrodeArrayPair")
    layout.radius = [10]
    layout.current = [CURRENT, -CURRENT]
roi = opt.add_roi()
roi.method = "volume"
roi.tissues = [2, 1]
roi.roi_sphere_center_space = "subject"
roi.roi_sphere_center = list(map(float, target))
roi.roi_sphere_radius = 8
t0 = time.time()
opt.run(cpus=int(os.environ.get("OMP_NUM_THREADS", 4)))
seconds = time.time() - t0
results = {"target_subject_mm": list(map(float, target)), "maxiter": MAXITER, "popsize": POPSIZE,
           "current_mA": CURRENT * 1e3, "seconds": seconds, "output_folder": opt.output_folder,
           "goal_fun_value": float(np.ravel(opt.goal_fun_value)[0]) if opt.goal_fun_value is not None else None,
           "optim_funvalue": float(opt.optim_funvalue) if opt.optim_funvalue is not None else None,
           "electrode_pos_opt": np.asarray(opt.electrode_pos_opt, float).tolist() if opt.electrode_pos_opt is not None else None,
           "final_sim_files": [os.path.basename(f) for f in (opt.fn_final_sim or [])],
           "reference": "Rampersad 2019 L hippocampus optimum 0.24 V/m at 2 mA total; Violante 2023 0.26 V/m at 1 mA per pair on MIDA"}
summaries = glob.glob(os.path.join(opt.output_folder, "*summary*")) + glob.glob(os.path.join(opt.output_folder, "*.json"))
results["files"] = [os.path.basename(f) for f in summaries]
results["note"] = "goal_fun_value is the first evaluation; optim_funvalue the best (node-wise current sources); final_* is the real-electrode simulation"
# the optimiser's head mesh carries max_TI and the ROI of the final real-electrode simulation
fn_head = os.path.join(opt.output_folder, f"{os.path.basename(os.path.normpath(M2M)).removeprefix('m2m_')}_tes_flex_opt_head_mesh.msh")
if os.path.exists(fn_head):
    from simnibs import mesh_io

    hm = mesh_io.read_msh(fn_head)
    tets = hm.elm.elm_type == 4
    am, roi_f = hm.field["max_TI"].value[tets], hm.field["ROI"].value[tets]
    vol, tags = hm.elements_volumes_and_areas().value[tets], hm.elm.tag1[tets]
    roi, brain = roi_f > 0, np.isin(tags, (1, 2))
    results["final_roi_mean_V_per_m"] = float(np.average(am[roi], weights=vol[roi]))
    results["final_roi_median_V_per_m"] = float(np.median(am[roi]))
    results["final_gm_p99_V_per_m"] = float(np.percentile(am[tags == 2], 99))
    results["final_brain_above_0.25_outside_roi_cm3"] = float(vol[brain & ~roi & (am > 0.25)].sum() / 1e3)
    print(f"final real-electrode simulation: ROI mean {results['final_roi_mean_V_per_m']:.3f} V/m, median "
          f"{results['final_roi_median_V_per_m']:.3f}, GM p99 {results['final_gm_p99_V_per_m']:.3f}, "
          f"brain > 0.25 V/m outside the ROI {results['final_brain_above_0.25_outside_roi_cm3']:.0f} cm3")
print(f"flex TI optimisation: {seconds:.0f} s, first evaluation {results['goal_fun_value']}, best {results['optim_funvalue']}")
print(f"electrode positions: {results['electrode_pos_opt']}")
save_results(results)
