"""The temporal interference envelope: Grossman et al. 2017, Cell 169:1029, equation 1.

Two fields E1(r) cos(2 pi f1 t) and E2(r) cos(2 pi f2 t) with |f1 - f2| = Delta f add to a
kHz carrier whose amplitude is modulated at Delta f. Along a unit vector n the modulation
amplitude (peak-to-trough of the envelope, in V/m) is

    E_AM(n) = | |(E1 + E2) . n| - |(E1 - E2) . n| |

and its maximum over n has the closed form used by SimNIBS (TI_utils.get_maxTI): with
|E1| >= |E2| and the angle alpha between them made acute by flipping E2 if needed,

    E_AM,max = 2 |E2|                                  if |E2| <= |E1| cos(alpha)
             = 2 |E2 x (E1 - E2)| / |E1 - E2|          otherwise.

`max_envelope` reimplements it, `max_envelope_sweep` checks it by brute force over
directions, `tacs_sum` is the same montage driven in phase (the comparison every TI
paper makes), and `multi_pair_lower_bound` is the n-field approximation 2 min|E_i| that
Grossman uses for n > 2 fields. Everything takes (..., 3) arrays in V/m.
"""

import numpy as np


def _unit(v, axis=-1):
    n = np.linalg.norm(v, axis=axis, keepdims=True)
    return np.divide(v, n, out=np.zeros_like(v, dtype=float), where=n > 0)


def directional_envelope(e1, e2, n):
    """E_AM along n (unit vector(s), broadcastable to e1): | |(E1+E2).n| - |(E1-E2).n| |."""
    e1, e2 = np.asarray(e1, float), np.asarray(e2, float)
    n = _unit(np.broadcast_to(np.asarray(n, float), e1.shape))
    return np.abs(np.abs(((e1 + e2) * n).sum(-1)) - np.abs(((e1 - e2) * n).sum(-1)))


def max_envelope(e1, e2):
    """Maximum of E_AM over all directions (closed form; SimNIBS get_maxTI)."""
    e1, e2 = np.asarray(e1, float).copy(), np.asarray(e2, float).copy()
    n1, n2 = np.linalg.norm(e1, axis=-1), np.linalg.norm(e2, axis=-1)
    swap = n2 > n1  # ensure |E1| >= |E2|
    e1[swap], e2[swap] = e2[swap], e1[swap]
    n1, n2 = np.where(swap, n2, n1), np.where(swap, n1, n2)
    flip = (e1 * e2).sum(-1) < 0  # ensure alpha <= pi/2
    e2[flip] = -e2[flip]
    with np.errstate(invalid="ignore", divide="ignore"):
        cos_a = np.where(n1 * n2 > 0, (e1 * e2).sum(-1) / (n1 * n2), 1.0)
        d = e1 - e2
        general = 2 * np.linalg.norm(np.cross(e2, d), axis=-1) / np.linalg.norm(d, axis=-1)
    out = np.where(n2 <= n1 * cos_a, 2 * n2, general)
    return np.where(np.isfinite(out), out, 0.0)


def fibonacci_directions(n=2000):
    """n roughly uniform unit vectors on the sphere (brute-force direction sweeps)."""
    i = np.arange(n) + 0.5
    z = 1 - 2 * i / n
    r = np.sqrt(1 - z * z)
    phi = np.pi * (1 + 5**0.5) * i
    return np.stack([r * np.cos(phi), r * np.sin(phi), z], -1)


def max_envelope_sweep(e1, e2, n_dirs=4000):
    """Maximum of E_AM by sampling directions; an independent check of max_envelope."""
    e1, e2 = np.asarray(e1, float), np.asarray(e2, float)
    dirs = fibonacci_directions(n_dirs)
    best = np.zeros(e1.shape[:-1])
    for n in dirs:
        best = np.maximum(best, directional_envelope(e1, e2, n))
    return best


def envelope_direction(e1, e2, n_dirs=4000):
    """Unit vector along which E_AM is largest (sampled), one per point."""
    e1, e2 = np.asarray(e1, float), np.asarray(e2, float)
    dirs = fibonacci_directions(n_dirs)
    vals = np.stack([directional_envelope(e1, e2, n) for n in dirs], -1)
    return dirs[np.argmax(vals, axis=-1)]


def tacs_sum(e1, e2):
    """|E1 + E2| or |E1 - E2|, whichever is larger: the same electrodes driven in phase at a
    single frequency, with the polarity of one pair chosen for the larger field
    (Rampersad et al. 2019, eq. 4). The unmodulated field TI is compared against."""
    e1, e2 = np.asarray(e1, float), np.asarray(e2, float)
    return np.maximum(np.linalg.norm(e1 + e2, axis=-1), np.linalg.norm(e1 - e2, axis=-1))


def multi_pair_lower_bound(*fields):
    """2 min_i |E_i|: Grossman's approximation of the envelope for n > 2 carriers."""
    return 2 * np.min([np.linalg.norm(np.asarray(f, float), axis=-1) for f in fields], axis=0)


def waveform(e1, e2, f1, f2, t):
    """The instantaneous field (..., len(t), 3) of two carriers, for demonstrations and
    for feeding a membrane model: E1 cos(2 pi f1 t) + E2 cos(2 pi f2 t)."""
    e1, e2, t = np.asarray(e1, float), np.asarray(e2, float), np.asarray(t, float)
    return (e1[..., None, :] * np.cos(2 * np.pi * f1 * t)[:, None]
            + e2[..., None, :] * np.cos(2 * np.pi * f2 * t)[:, None])


def modulation_depth(e1, e2):
    """Envelope maximum minus minimum along the direction of the larger total field,
    divided by the maximum: 1 means the carrier is fully modulated at the beat frequency,
    0 means an unmodulated kHz field (Gomez-Tames 2021 uses max - min and min itself)."""
    e1, e2 = np.asarray(e1, float), np.asarray(e2, float)
    n = _unit(e1 + e2)
    a, b = np.abs(((e1 + e2) * n).sum(-1)), np.abs(((e1 - e2) * n).sum(-1))
    hi, lo = np.maximum(a, b), np.minimum(a, b)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(hi > 0, (hi - lo) / hi, 0.0)


# --- n carriers (Savvateev 2026 multi-pair, Huang 2026 two-array, TI-Toolbox mTI) ---

def multi_pair_envelope(fields, n):
    """Exact peak-to-trough envelope along unit vector(s) n for any number of carriers at
    distinct frequencies, phase free: with a_i = E_i . n the baseband amplitude is the
    modulus of a sum of phasors of fixed lengths |a_i| whose relative phases drift through
    every value (incommensurate beat frequencies), so the envelope ranges from
    max(0, 2 max|a_i| - sum|a_i|) to sum|a_i| and the modulation amplitude is

        E_AM(n) = min( sum_i |a_i|,  2 (sum_i |a_i| - max_i |a_i|) ).

    For two carriers this is 2 min(|a_1|, |a_2|), Grossman's equation 1; for more it is at
    least Grossman's 2 min|E_i| approximation (`multi_pair_lower_bound`). When the beat
    frequencies are commensurate (equally spaced carriers) the phases are locked and this is an
    upper bound; `multi_tone_envelope_bounds` gives the phase-locked value from the waveform."""
    fields = [np.asarray(f, float) for f in fields]
    n = _unit(np.broadcast_to(np.asarray(n, float), fields[0].shape))
    a = np.stack([np.abs((f * n).sum(-1)) for f in fields])  # (k, ...)
    total = a.sum(0)
    return np.minimum(total, 2 * (total - a.max(0)))


def max_multi_pair_envelope(fields, n_dirs=4000):
    """Maximum of multi_pair_envelope over sampled directions (no closed form for k > 2)."""
    fields = [np.asarray(f, float) for f in fields]
    best = np.zeros(fields[0].shape[:-1])
    for d in fibonacci_directions(n_dirs):
        best = np.maximum(best, multi_pair_envelope(fields, d))
    return best


def multi_tone_envelope_bounds(fields, freqs, t, n, phases=None):
    """(max, min) of the instantaneous amplitude envelope along n over the times t for
    carriers E_i cos(2 pi f_i t + phi_i): the baseband |sum_i a_i exp(i(2 pi (f_i - f_1) t
    + phi_i))| with a_i = E_i . n, evaluated at one point (fields of shape (3,) or (1, 3)).
    With incommensurate beats max - min tends to multi_pair_envelope; with locked phases
    (commensurate carriers) it is what the tissue actually sees."""
    fields = [np.asarray(f, float).reshape(-1, 3)[0] for f in fields]
    n = np.asarray(n, float); n = n / np.linalg.norm(n)
    freqs = np.asarray(freqs, float); t = np.asarray(t, float)
    phases = np.zeros(len(fields)) if phases is None else np.asarray(phases, float)
    a = np.array([f @ n for f in fields])
    baseband = sum(a_i * np.exp(1j * (2 * np.pi * (f_i - freqs[0]) * t + p_i))
                   for a_i, f_i, p_i in zip(a, freqs, phases))
    amp = np.abs(baseband)
    return float(amp.max()), float(amp.min())


def mti_envelope(pairs, n, envelope_phases=None):
    """Multipolar TI (Botzanowski et al. 2023; TI-Toolbox's mTI): k electrode pairs in
    separate carrier bands, each pair beating at the same Delta f (e.g. 1975/2025 Hz and
    2975/3025 Hz). The inter-band beats (~1 kHz) are above what the membrane follows, so the
    tissue sees the sum of the k low-frequency envelopes. Along n each pair contributes its
    two-field envelope E_AM,k(n) (`directional_envelope`); with the envelopes in phase the
    total is sum_k E_AM,k(n), and with phase offsets phi_k between the Delta f envelopes it is
    |sum_k E_AM,k(n) exp(i phi_k)|. `pairs` = [(E1, E2), ...] with (..., 3) arrays."""
    amps = [directional_envelope(e1, e2, n) for e1, e2 in pairs]
    if envelope_phases is None:
        return np.sum(amps, axis=0)
    ph = np.asarray(envelope_phases, float)
    return np.abs(sum(a * np.exp(1j * p) for a, p in zip(amps, ph)))


def max_mti_envelope(pairs, n_dirs=4000, envelope_phases=None):
    """Maximum of mti_envelope over sampled directions."""
    first = np.asarray(pairs[0][0], float)
    best = np.zeros(first.shape[:-1])
    for d in fibonacci_directions(n_dirs):
        best = np.maximum(best, mti_envelope(pairs, d, envelope_phases))
    return best
