"""Label volume -> SimNIBS mesh -> TI, on a synthetic 4-layer sphere (Grossman 2017 figure
S2 proportions at a 40 mm radius): the whole custom-mesh route MIDA will take. Slow
(CGAL meshing plus two FEM solves, ~1-3 min); skipped without SimNIBS."""

import os

import numpy as np
import pytest

from ti_sim import labelvol
from ti_sim.head import GM, TwoPairSimulation
from ti_sim.tissues import GROSSMAN_SPHERE

simnibs = pytest.importorskip("simnibs")


def test_sphere_labels_geometry():
    labels, aff = labelvol.sphere_labels((10, 8, 6), (5, 4, 2), spacing_mm=1.0, pad_mm=2,
                                         blobs=[((0, 0, 0), 2, 21)])
    assert labels.shape == (25, 25, 25) and labels[12, 12, 12] == 21
    assert labels[12, 12, 12 + 7] == 4 and labels[12, 12, 12 + 9] == 5 and labels[0, 0, 0] == 0
    p = labelvol.surface_point(labels, aff, (0, 0, 1))
    assert np.isclose(p[2], 10.0, atol=1.0) and np.allclose(p[:2], 0, atol=1e-9)


@pytest.mark.slow
def test_layered_sphere_ti(tmp_path):
    R = 40.0
    radii = (R, R * (1 - 0.05), R * (1 - 0.05 - 0.085), R * (1 - 0.05 - 0.085 - 0.023))
    labels, aff = labelvol.sphere_labels(radii, (5, 4, 3, 2), spacing_mm=1.0,
                                         blobs=[((0, 0, 0), 5, 21)])
    work = os.environ.get("TI_SIM_WORK", str(tmp_path))
    fn = labelvol.mesh_from_labels(labels, aff, work, "sphere40", extra_tags=(21,),
                                   elem_sizes={"standard": {"range": [1.5, 5.0], "slope": 1.0}})
    m = simnibs.mesh_io.read_msh(fn)
    tags = set(np.unique(m.elm.tag1).tolist())
    assert {2, 3, 4, 5, 21, 1005} <= tags
    cond = {1: 0.333, 2: GROSSMAN_SPHERE["brain"], 3: GROSSMAN_SPHERE["csf"],
            4: GROSSMAN_SPHERE["skull"], 5: GROSSMAN_SPHERE["scalp"], 6: 0.5, 7: 0.0083,
            8: 0.0083, 9: 0.6, 10: 0.16, 11: 0.88, 12: 0.078}
    # two pairs on the equator, each spanning 90 degrees, 45 degrees apart (one leg each)
    def on_sphere(deg):
        a = np.radians(deg)
        return [R * np.cos(a), R * np.sin(a), 0.0]
    pairs = [(on_sphere(-20), on_sphere(200)), (on_sphere(20), on_sphere(160))]
    sim = TwoPairSimulation(pairs, (1e-3, 1e-3), work, fnamehead=fn, radius=4.0, cond=cond,
                            cond_label="grossman-sphere", extra_cond={21: GROSSMAN_SPHERE["brain"]})
    centre = sim.stats(sim.in_tags(21))
    brain = sim.stats(sim.in_tags(GM))
    assert centre["n_elements"] > 20 and np.isfinite(centre["am"]["median"])
    assert 0.01 < centre["am"]["median"] < 20.0  # ~1 mA across a 4 cm sphere: V/m scale
    assert brain["am"]["max"] >= centre["am"]["median"]
