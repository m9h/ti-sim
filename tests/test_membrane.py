"""A passive membrane does not demodulate a TI signal; a rectifier does."""

import numpy as np

from ti_sim import membrane as mb


def _signal():
    t = np.arange(0, 0.5, 1 / 40000.0)  # 0.5 s at 40 kHz
    return t, mb.ti_signal(t, 1.0, 1.0, 2000.0, 2010.0)  # 10 Hz beat


def test_lowpass_attenuates_but_does_not_demodulate():
    t, x = _signal()
    v = mb.lowpass(t, x, tau=5e-3)  # somatic time constant
    assert v.std() < 0.05 * x.std()  # 2 kHz through a 5 ms RC: ~1/(2 pi f tau) = 1.6 %
    assert mb.beat_power_fraction(t, x, 10.0) < 1e-3  # the raw signal has no 10 Hz power
    assert mb.beat_power_fraction(t, v, 10.0) < 1e-3  # and neither does its low-pass


def test_rectification_recovers_the_beat():
    t, x = _signal()
    v = mb.lowpass(t, mb.rectify(x), tau=5e-3)
    assert mb.beat_power_fraction(t, v, 10.0) > 0.5


def test_lowpass_is_exact_for_a_step():
    t = np.linspace(0, 0.05, 5001)
    v = mb.lowpass(t, np.ones_like(t), tau=0.01)
    assert np.allclose(v, 1 - np.exp(-t / 0.01), atol=1e-12)  # input held from the previous sample


def test_esmaeilpour_table():
    assert mb.esmaeilpour_threshold(2000) == (80.0, 220.0)
    assert mb.esmaeilpour_threshold(120) == (5.0, 13.0)
