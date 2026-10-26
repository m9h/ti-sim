"""A Hodgkin-Huxley point neuron under a temporal-interference field: the cheapest readout
that has the physics the envelope lacks (Mirzakhalili 2020: demodulation needs the
channels' nonlinearity; Gomez-Tames 2021, Esmaeilpour 2021: the threshold rises with the
carrier frequency).

The extracellular field polarises the membrane by lambda * E(t) (Radman et al. 2009:
~0.2 mV per V/m for cortical somata at DC; `polarisation_mm` is that length, 0.2 mm
default). The polarisation is applied as an injected current g_L * lambda * E(t) (the
membrane's leak conductance times the drive), so at DC the steady depolarisation is
lambda * E. Squid HH kinetics at 6.3 C (Hodgkin & Huxley 1952) with the standard
parameters; `warm` scales the rates (Q10 = 3) towards 37 C where the membrane is faster and
follows higher carriers (Esmaeilpour: faster time constants raise the sensitivity).

`threshold(carrier_hz, beat_hz, ...)` bisects the field amplitude (equal E1 = E2, so
E_AM,max = 2 E) at which the neuron fires at least `min_spikes` spikes in `duration`.
Units: mV, ms, uA/cm^2, V/m.
"""

import numpy as np

C_M, G_NA, G_K, G_L = 1.0, 120.0, 36.0, 0.3  # uF/cm^2, mS/cm^2
E_NA, E_K, E_L = 50.0, -77.0, -54.4  # mV
V_REST = -65.0


V_CLIP = 150.0  # mV: a diverging Euler step is clipped instead of overflowing the rate exponentials


def _rates(v):
    v = min(max(v, -V_CLIP), V_CLIP)
    an = 0.01 * (v + 55) / (1 - np.exp(-(v + 55) / 10)) if abs(v + 55) > 1e-6 else 0.1
    bn = 0.125 * np.exp(-(v + 65) / 80)
    am = 0.1 * (v + 40) / (1 - np.exp(-(v + 40) / 10)) if abs(v + 40) > 1e-6 else 1.0
    bm = 4.0 * np.exp(-(v + 65) / 18)
    ah = 0.07 * np.exp(-(v + 65) / 20)
    bh = 1 / (1 + np.exp(-(v + 35) / 10))
    return an, bn, am, bm, ah, bh


def simulate(t, e_field, polarisation_mm=0.2, warm=1.0, v0=V_REST):
    """Membrane potential (mV) over time t (ms) for a field e_field (V/m) sampled on t;
    forward Euler on the HH equations (dt <= 0.01 ms recommended for kHz carriers)."""
    t = np.asarray(t, float)
    e = np.asarray(e_field, float)
    lam = polarisation_mm * 1e-3  # m: mV per V/m is lam * 1e3 ... (V/m * m = V = 1e3 mV) -> mV = E * lam * 1e3
    v = v0
    an, bn, am, bm, ah, bh = _rates(v)
    n, m, h = an / (an + bn), am / (am + bm), ah / (ah + bh)
    out = np.empty(len(t))
    out[0] = v
    for i in range(1, len(t)):
        dt = t[i] - t[i - 1]
        drive = G_L * lam * 1e3 * e[i - 1]  # uA/cm^2: g_L * (mV of polarisation)
        i_ion = G_NA * m**3 * h * (v - E_NA) + G_K * n**4 * (v - E_K) + G_L * (v - E_L)
        v = min(max(v + dt * (drive - i_ion) / C_M, -V_CLIP), V_CLIP)
        an, bn, am, bm, ah, bh = (r * warm for r in _rates(v))
        n += dt * (an * (1 - n) - bn * n)
        m += dt * (am * (1 - m) - bm * m)
        h += dt * (ah * (1 - h) - bh * h)
        out[i] = v
    return out


def spikes(v, t, level=0.0):
    """Number of upward crossings of `level` (mV)."""
    above = v > level
    return int(np.count_nonzero(above[1:] & ~above[:-1]))


def ti_field(t_ms, amplitude, carrier_hz, beat_hz):
    """E(t) = A cos(2 pi f1 t) + A cos(2 pi f2 t) with f2 = f1 + beat (V/m); beat 0 gives an
    unmodulated carrier of amplitude 2A; carrier 0 gives a plain beat-frequency field."""
    t = np.asarray(t_ms, float) * 1e-3
    if carrier_hz == 0:
        return 2 * amplitude * np.cos(2 * np.pi * beat_hz * t)
    return amplitude * (np.cos(2 * np.pi * carrier_hz * t) + np.cos(2 * np.pi * (carrier_hz + beat_hz) * t))


def fires(amplitude, carrier_hz, beat_hz, duration_ms=200.0, dt_ms=0.005, min_spikes=2, **kw):
    t = np.arange(0.0, duration_ms, dt_ms)
    v = simulate(t, ti_field(t, amplitude, carrier_hz, beat_hz), **kw)
    return spikes(v, t) >= min_spikes


def threshold(carrier_hz, beat_hz=10.0, lo=1.0, hi=2000.0, tol=0.05, **kw):
    """Field amplitude A (V/m per carrier; E_AM,max = 2 A) at which the neuron fires, by
    bisection between lo and hi; NaN when hi does not fire."""
    if not fires(hi, carrier_hz, beat_hz, **kw):
        return float("nan")
    if fires(lo, carrier_hz, beat_hz, **kw):
        return lo
    while hi / lo > 1 + tol:
        mid = np.sqrt(lo * hi)
        if fires(mid, carrier_hz, beat_hz, **kw):
            hi = mid
        else:
            lo = mid
    return float(hi)
