"""NEURON cable axon under a TI field (ti_sim.axon); skipped without NEURON."""

import numpy as np
import pytest

neuron = pytest.importorskip("neuron")

from ti_sim import axon  # noqa: E402


def test_axon_geometry_and_rest():
    ax = axon.Axon(length_mm=4.0, diam_um=2.0, nseg=41)
    assert len(ax.x_m) == 41 and np.isclose(ax.x_m[-1] - ax.x_m[0], 4e-3 * 40 / 41, rtol=1e-6)
    t = np.arange(0.0, 20.0, 0.025)
    _, v = ax.run(t, np.zeros_like(t))
    assert abs(v[-1] - v[0]) < 1.0  # stays at rest without a field


@pytest.mark.slow
def test_field_along_the_axis_fires_and_orthogonal_does_not():
    ax = axon.Axon(length_mm=10.0, diam_um=2.0, nseg=101)
    t = np.arange(0.0, 60.0, 0.025)
    strong = axon.ti_field(t, 100.0, 0.0, 10.0)  # plain 10 Hz field, 200 V/m: ~5x the cable's threshold
    assert ax.spikes(t, strong) >= 1
    orthogonal = np.stack([np.zeros_like(strong), strong, np.zeros_like(strong)], 1)
    assert ax.spikes(t, orthogonal) == 0


@pytest.mark.slow
def test_threshold_rises_with_the_carrier():
    ax = axon.Axon(length_mm=10.0, diam_um=2.0, nseg=101)
    thr = {f: axon.threshold(ax, f, 10.0, lo=2.0, hi=1500.0, duration_ms=120.0) for f in (0.0, 1000.0, 2000.0)}
    # 10 mm HH cable at 6.3 C: ~20 V/m per carrier for a plain 10 Hz field, ~200 at 1 kHz, ~500 at 2 kHz
    assert 10 < thr[0.0] < 40
    assert thr[1000.0] > 4 * thr[0.0] and thr[2000.0] > thr[1000.0]
