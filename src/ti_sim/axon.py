"""A straight axon in NEURON under the field of a TI montage: the readout between the point
neuron (`ti_sim.hh`) and a full MRG / layer-5 model.

The extracellular potential along a straight axon in a locally uniform field E(t) is
phi(x, t) = -E(t) . x, applied through NEURON's `extracellular` mechanism (e_extracellular
at each segment), so the cable sees the activating function of the field's projection on
its axis (Rattay 1986). Hodgkin-Huxley membrane throughout (NEURON's `hh` mechanism, with
its own Q10 scaling from h.celsius), sealed ends; the field is any (len(t),) or (len(t), 3)
array, e.g. `ti_field` or a montage's E1, E2 at a point put through `envelope.waveform`.

Nothing here needs compiled mod files, so it runs wherever `pip install neuron` works
(aarch64 included). `threshold` bisects the field amplitude for one propagating spike.
"""

import numpy as np


def _neuron():
    try:
        from neuron import h
    except ImportError as e:  # pragma: no cover
        raise ImportError("ti_sim.axon needs NEURON: uv pip install neuron") from e
    h.load_file("stdrun.hoc")
    return h


class Axon:
    """Unmyelinated HH cable of `length_mm`, `diam_um`, `nseg` segments, along `axis`
    (unit vector), centred on the origin; `celsius` sets the kinetics (squid HH stops
    firing above ~30 C, so the default is Hodgkin-Huxley's 6.3 C; `ti_sim.hh` warms the
    point neuron the same way)."""

    def __init__(self, length_mm=10.0, diam_um=2.0, nseg=101, axis=(1.0, 0.0, 0.0), celsius=6.3):
        h = self.h = _neuron()
        self.axis = np.asarray(axis, float) / np.linalg.norm(axis)
        self.sec = h.Section(name="axon")
        self.sec.L, self.sec.diam, self.sec.nseg = length_mm * 1e3, diam_um, int(nseg)
        self.sec.Ra, self.sec.cm = 100.0, 1.0
        self.sec.insert("hh")
        self.sec.insert("extracellular")
        h.celsius = celsius
        # segment centres along the axis, metres from the cable midpoint
        self.x_m = np.array([(seg.x - 0.5) * self.sec.L * 1e-6 for seg in self.sec])
        self._segs = list(self.sec)

    def run(self, t_ms, e_field, record_at=0.5):
        """Membrane potential (mV) at `record_at` for a field e_field (V/m, shape (len(t),)
        along the axis or (len(t), 3) in space) sampled on t_ms; returns (t_rec, v_rec)."""
        h = self.h
        e = np.asarray(e_field, float)
        e_axis = e @ self.axis if e.ndim == 2 else e
        t_ms = np.asarray(t_ms, float)
        phi = -np.outer(e_axis, self.x_m) * 1e3  # mV: V/m * m = V -> mV
        tvec = h.Vector(t_ms)
        self._played = [tvec]  # keep the played vectors alive for the run
        for j, seg in enumerate(self._segs):
            v = h.Vector(phi[:, j])
            v.play(seg._ref_e_extracellular, tvec, True)
            self._played.append(v)
        rec_v, rec_t = h.Vector().record(self.sec(record_at)._ref_v), h.Vector().record(h._ref_t)
        h.dt = float(min(np.diff(t_ms).min(), 0.0125))
        h.finitialize(-65.0)
        h.continuerun(float(t_ms[-1]))
        return np.array(rec_t), np.array(rec_v)

    def spikes(self, t_ms, e_field, record_at=0.5, level=0.0):
        _, v = self.run(t_ms, e_field, record_at)
        above = v > level
        return int(np.count_nonzero(above[1:] & ~above[:-1]))


def ti_field(t_ms, amplitude, carrier_hz, beat_hz):
    """E(t) along the axon: A cos(2 pi f1 t) + A cos(2 pi (f1 + beat) t) (V/m); carrier 0
    gives a plain beat-frequency field of amplitude 2 A."""
    t = np.asarray(t_ms, float) * 1e-3
    if carrier_hz == 0:
        return 2 * amplitude * np.cos(2 * np.pi * beat_hz * t)
    return amplitude * (np.cos(2 * np.pi * carrier_hz * t) + np.cos(2 * np.pi * (carrier_hz + beat_hz) * t))


def threshold(axon, carrier_hz, beat_hz=10.0, lo=1.0, hi=5000.0, tol=0.05, duration_ms=150.0,
              dt_ms=0.0125, min_spikes=1):
    """Field amplitude A (V/m per carrier, E_AM,max = 2 A) at which `axon` fires at least
    `min_spikes` spikes at its midpoint, by bisection; NaN when `hi` does not fire."""
    t = np.arange(0.0, duration_ms + dt_ms, dt_ms)

    def fires(a):
        return axon.spikes(t, ti_field(t, a, carrier_hz, beat_hz)) >= min_spikes

    if not fires(hi):
        return float("nan")
    if fires(lo):
        return lo
    while hi / lo > 1 + tol:
        mid = float(np.sqrt(lo * hi))
        if fires(mid):
            hi = mid
        else:
            lo = mid
    return float(hi)
