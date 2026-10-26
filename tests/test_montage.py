"""The montage search on a synthetic lead field: it finds the montage whose pairs produce
equal, parallel fields at the target, and respects disjoint pairs and the current split."""

import numpy as np

from ti_sim import montage


class FakeLF:
    """Four electrodes; the field at 3 ROI points is a fixed vector per electrode."""

    def __init__(self):
        self.vec = {"A": np.array([1.0, 0, 0]), "B": np.array([-1.0, 0, 0]),
                    "C": np.array([0.5, 0, 0]), "D": np.array([-0.5, 0, 0]), "E": np.array([0, 0.1, 0])}
        self.idx = {k: i for i, k in enumerate(self.vec)}
        self.idx["ref"] = None
        self.names = list(self.idx)
        self.centres = np.zeros((3, 3))
        self.volumes = np.ones(3)
        self.tags = np.array([2, 2, 1])

    def field(self, a, b, current):
        v = (self.vec.get(a, 0) - self.vec.get(b, 0)) * current
        return np.tile(v, (3, 1))


def test_search_prefers_equal_parallel_fields():
    lf = FakeLF()
    res = montage.search_two_pairs(lf, np.array([True, True, False]), ratios=(1.0, 2.0),
                                   total_current=2e-3, top=3, threshold=1e-3)
    best = res["best"][0]
    # A-D and B-C both give 1.5 V/m per A along x (opposite signs, which the envelope
    # ignores): equal and parallel at 1:1, so 2 * 1.5 mA-fields = 3e-3 V/m
    assert {tuple(sorted(best["pair1"])), tuple(sorted(best["pair2"]))} == {("A", "D"), ("B", "C")}
    assert np.isclose(best["I1_mA"], 1.0) and np.isclose(best["I2_mA"], 1.0)
    assert np.isclose(best["roi_median"], 3e-3)
    # A-B (2 V/m per A) with C-D (1 V/m per A) is next, best at the 1:2 split (2/3 and 4/3 mA)
    second = res["best"][1]
    assert {tuple(sorted(second["pair1"])), tuple(sorted(second["pair2"]))} == {("A", "B"), ("C", "D")}
    assert np.isclose(second["roi_median"], 2 * 2.0 * (2.0 / 3) * 1e-3)  # 2 * (2 V/m/A * 2/3 mA)
    assert res["n_electrodes"] == 5 and res["n_pairs"] == 10
    assert "above_threshold_mm3" in best and best["above_threshold_mm3"]["roi"] == 2.0


def test_focal_ranking_prefers_less_outside_volume():
    lf = FakeLF()
    lf.tags = np.array([2, 2, 2])  # all three points are brain; ROI = first point only
    res = montage.search_two_pairs(lf, np.array([True, False, False]), ratios=(1.0, 2.0), threshold=1e-3,
                                   top=2, shortlist=10, min_roi_fraction=0.5)
    assert res["shortlist"] == 10 and len(res["focal"]) == 2  # 15 disjoint pair combos x 3 ratios, 10 shortlisted
    # uniform fake fields: everything above threshold outside too, so the focal ranking falls
    # back to the ROI median order among the admissible montages
    assert res["focal"][0]["roi_median"] >= res["focal"][1]["roi_median"]
    assert all(r["roi_median"] >= 0.5 * res["best"][0]["roi_median"] for r in res["focal"])


def test_electrode_subset_and_near_target():
    pos = {"A": (0, 0, 0), "B": (10, 0, 0), "C": (50, 0, 0)}
    assert montage.near_target(pos, (0, 0, 0), 20) == ["A", "B"]
    lf = FakeLF()
    res = montage.search_two_pairs(lf, np.array([True, True, True]), electrodes=["A", "B", "C", "E"])
    assert res["n_pairs"] == 6 and res["n_montages"] == 3  # 3 ways to split 4 into 2 disjoint pairs
