"""TI on the SimNIBS ernie example head (slow: two FEM solves, ~1-2 min on the lowres head).

Checks the numbers every TI paper reports are in range: with 1 mA per pair at F5-P5 and
F6-P6 (SimNIBS's TI example), the envelope in grey matter peaks below 1 V/m, the median
deep-brain envelope is 0.05-0.4 V/m (Rampersad 2019: 0.06-0.24 V/m in hippocampus and
pallidum; Violante 2023: 0.26 V/m hippocampus on MIDA), and tACS with the same electrodes
has the larger cortical maximum (Rampersad: tTIS is weaker but more focal).
"""

import os

import numpy as np
import pytest

from ti_sim.head import GM, WM, TwoPairSimulation

ERNIE = os.environ.get("TI_SIM_M2M", "")


@pytest.mark.skipif(not os.path.isdir(ERNIE), reason="set TI_SIM_M2M=/path/to/m2m_ernie")
def test_two_pair_ti_on_ernie(tmp_path):
    work = os.environ.get("TI_SIM_WORK", str(tmp_path))
    sim = TwoPairSimulation([("F5", "P5"), ("F6", "P6")], (1e-3, 1e-3), work, m2m_dir=ERNIE, radius=20)
    gm = sim.in_tags(GM)
    s = sim.stats(gm)
    assert 0.02 < s["am"]["median"] < 0.5 and s["am"]["max"] < 2.0
    assert s["tacs"]["max"] > s["am"]["max"]  # unmodulated in-phase field is the stronger one
    # a deep sphere between the pairs: thalamus-ish, 20 mm radius at the head centre (GM+WM)
    centre = sim.centres[sim.in_tags((GM, WM))].mean(axis=0)
    deep = sim.in_sphere(centre, 20.0, tags=(GM, WM))
    d = sim.stats(deep)
    assert d["n_elements"] > 100 and 0.02 < d["am"]["median"] < 0.5
    foc = sim.focality(0.2)
    assert foc["am_mm3"] <= foc["tacs_mm3"]  # TI never covers more brain than tACS
    assert np.isfinite(sim.am_max).all()
    fn = sim.write_msh(os.path.join(str(tmp_path), "ti.msh"))
    assert os.path.getsize(fn) > 1000
