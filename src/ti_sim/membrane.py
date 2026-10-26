"""What a neuron sees: the kHz carrier through a passive membrane, and why that is not
enough (Mirzakhalili et al. 2020, Cell Systems 11:557).

A passive membrane is a first-order low-pass with tau ~ 1-10 ms for somata, faster for
axons. Driven by E(t) = E1 cos(2 pi f1 t) + E2 cos(2 pi f2 t) it attenuates both carriers
by ~1/(2 pi f tau) and produces *no* component at the beat frequency: the sum of two
sinusoids is still the sum of two sinusoids, however it is filtered. Demodulation needs a
nonlinearity (rectification by voltage-gated channels), which is also what causes kHz
conduction block off target. `lowpass` and `beat_power_fraction` make that quantitative,
so an example can show the arithmetic rather than assert it.

Esmaeilpour et al. 2021 (Brain Stimul 14:55) measured the consequence in slices: with a
2 kHz carrier, modulating gamma oscillations took ~80 V/m amplitude-modulated field, versus
~5 V/m at a 100 Hz carrier. `esmaeilpour_threshold` tabulates those numbers.
"""

import numpy as np
from scipy import signal

# Esmaeilpour et al. 2021, table from the abstract: carrier frequency (Hz) -> field
# amplitude (V/m) that significantly modulated carbachol gamma in rat hippocampal slices
# with 5 Hz amplitude modulation, and the scalp current they estimate it corresponds to (mA).
ESMAEILPOUR_2021 = {100: (5.0, 13.0), 1000: (60.0, 160.0), 2000: (80.0, 220.0)}


def esmaeilpour_threshold(carrier_hz):
    """(V/m, mA at the scalp) for the nearest tabulated carrier."""
    f = min(ESMAEILPOUR_2021, key=lambda k: abs(k - carrier_hz))
    return ESMAEILPOUR_2021[f]


def lowpass(t, x, tau):
    """First-order RC response v' = (x - v)/tau to x(t) sampled on t (exact for piecewise
    constant input); x may be (..., len(t))."""
    t, x = np.asarray(t, float), np.asarray(x, float)
    dt = np.diff(t)
    a = np.exp(-dt / tau)
    v = np.zeros_like(x)
    for i in range(1, len(t)):
        v[..., i] = a[i - 1] * v[..., i - 1] + (1 - a[i - 1]) * x[..., i - 1]
    return v


def rectify(x):
    """Half-wave rectification, the simplest stand-in for the channel nonlinearity."""
    return np.maximum(np.asarray(x, float), 0.0)


def band_power(t, x, f_lo, f_hi):
    """Power of x in [f_lo, f_hi] Hz (periodogram, Hann window)."""
    t, x = np.asarray(t, float), np.asarray(x, float)
    fs = 1 / np.mean(np.diff(t))
    f, p = signal.periodogram(x - x.mean(), fs, window="hann")
    sel = (f >= f_lo) & (f <= f_hi)
    return float(p[sel].sum())


def beat_power_fraction(t, x, f_beat, half_width=None):
    """Fraction of the (mean-removed) power of x within f_beat +- half_width."""
    hw = half_width if half_width is not None else max(0.5 * f_beat, 1.0)
    total = band_power(t, x, 0.0, np.inf)
    return band_power(t, x, f_beat - hw, f_beat + hw) / total if total > 0 else 0.0


def ti_signal(t, a1, a2, f1, f2):
    """Scalar two-carrier signal a1 cos(2 pi f1 t) + a2 cos(2 pi f2 t)."""
    t = np.asarray(t, float)
    return a1 * np.cos(2 * np.pi * f1 * t) + a2 * np.cos(2 * np.pi * f2 * t)
