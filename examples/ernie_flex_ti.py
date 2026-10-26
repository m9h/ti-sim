"""Free electrode placement for TI with SimNIBS's TesFlexOptimization (Weise et al.): two
freely movable pairs, 10 mm electrodes, 1 mA each (Violante 2023's currents), maximising
the mean E_AM,max in an 8 mm sphere at the left hippocampus (GM+WM volume ROI), by
differential evolution with a FEM solve per candidate.

The full optimisation is hours on the low-res ernie; TI_SIM_FLEX_MAXITER (default 2) and
TI_SIM_FLEX_POPSIZE (default 4) make a smoke run of ~15-30 min that still ends with the
final electrode simulation and its ROI value, which `ti_sim.montage`'s exhaustive search
and Rampersad's 0.24 V/m can be compared with. Set them to 1000 and 13 (SimNIBS defaults)
for a real run on a bigger machine.

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
print(f"flex TI optimisation: {seconds:.0f} s, goal value {results['goal_fun_value']}, optimiser value {results['optim_funvalue']}")
print(f"electrode positions: {results['electrode_pos_opt']}")
save_results(results)
