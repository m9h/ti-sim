"""HH point neuron under TI: rests quietly, fires for a strong DC-like drive, needs far
more field at a kHz carrier than at the beat frequency, and the threshold rises with the
carrier (Gomez-Tames 2021, Esmaeilpour 2021)."""

import numpy as np
import pytest

from ti_sim import hh


def test_rest_and_strong_drive():
    t = np.arange(0, 50, 0.01)
    v = hh.simulate(t, np.zeros_like(t))
    assert abs(v[-1] - hh.V_REST) < 1.0 and hh.spikes(v, t) == 0
    v = hh.simulate(t, np.full_like(t, 200.0))  # 200 V/m * 0.2 mm = 40 mV depolarisation
    assert hh.spikes(v, t) >= 1


def test_kilohertz_carrier_needs_much_more_field():
    thr_beat = hh.threshold(0.0, 10.0, hi=500.0, duration_ms=150)  # plain 10 Hz field
    thr_1k = hh.threshold(1000.0, 10.0, hi=4000.0, duration_ms=150)
    assert np.isfinite(thr_beat) and thr_beat < 150
    assert np.isnan(thr_1k) or thr_1k > 5 * thr_beat


@pytest.mark.slow
def test_threshold_rises_with_carrier():
    thr = {f: hh.threshold(f, 10.0, hi=4000.0, duration_ms=150, warm=3.0) for f in (100.0, 500.0, 1000.0, 2000.0)}
    vals = [thr[f] for f in sorted(thr)]
    finite = [v for v in vals if np.isfinite(v)]
    assert len(finite) >= 2 and all(np.diff(finite) >= 0)
