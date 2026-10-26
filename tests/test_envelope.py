"""Grossman's envelope formula: closed form vs brute force vs SimNIBS, and the limits."""

import numpy as np
import pytest

from ti_sim import envelope as env


def test_parallel_equal_fields_give_twice_the_field():
    e1 = np.array([[1.0, 0, 0]])
    assert np.isclose(env.directional_envelope(e1, e1, [1, 0, 0]), 2.0)
    assert np.isclose(env.max_envelope(e1, e1), 2.0)
    assert np.isclose(env.modulation_depth(e1, e1), 1.0)


def test_antiparallel_and_orthogonal_cases():
    e1 = np.array([[1.0, 0, 0]])
    assert np.isclose(env.max_envelope(e1, -e1), 2.0)  # polarity does not matter
    e2 = np.array([[0, 1.0, 0]])
    # orthogonal equal fields: the general branch, 2|E2 x (E1-E2)| / |E1-E2| = 2/sqrt(2)
    assert np.isclose(env.max_envelope(e1, e2), np.sqrt(2))
    assert np.isclose(env.directional_envelope(e1, e2, [1, 0, 0]), 0.0)  # along x no beat


def test_unequal_parallel_fields_modulate_at_twice_the_weaker():
    e1, e2 = np.array([[3.0, 0, 0]]), np.array([[1.0, 0, 0]])
    assert np.isclose(env.max_envelope(e1, e2), 2.0)
    assert np.isclose(env.modulation_depth(e1, e2), 0.5)  # (4 - 2) / 4


def test_closed_form_matches_direction_sweep():
    rng = np.random.default_rng(0)
    e1, e2 = rng.normal(size=(200, 3)), rng.normal(size=(200, 3)) * rng.uniform(0.1, 3, (200, 1))
    a, b = env.max_envelope(e1, e2), env.max_envelope_sweep(e1, e2, n_dirs=20000)
    assert np.all(b <= a * (1 + 1e-9))  # the sweep can only undershoot the true maximum
    assert np.allclose(a, b, rtol=1e-2)


def test_matches_simnibs_TI_utils():
    simnibs = pytest.importorskip("simnibs")
    from simnibs.utils import TI_utils

    rng = np.random.default_rng(1)
    e1, e2 = rng.normal(size=(500, 3)), rng.normal(size=(500, 3))
    assert np.allclose(env.max_envelope(e1, e2), TI_utils.get_maxTI(e1, e2))
    n = rng.normal(size=(500, 3))
    assert np.allclose(env.directional_envelope(e1, e2, n), TI_utils.get_dirTI(e1, e2, n))
    del simnibs


def test_zero_field_is_safe():
    z = np.zeros((2, 3))
    assert np.all(env.max_envelope(z, z) == 0)
    assert np.all(env.modulation_depth(z, z) == 0)


def test_tacs_and_multi_pair_helpers():
    e1, e2 = np.array([[1.0, 0, 0]]), np.array([[-0.5, 0, 0]])
    assert np.isclose(env.tacs_sum(e1, e2), 1.5)  # polarity flipped for the larger sum
    assert np.isclose(env.multi_pair_lower_bound(e1, e2, 2 * e1), 1.0)


def test_waveform_envelope_is_the_formula():
    # sample a 2 kHz / 2.01 kHz pair along x and measure the envelope peak-to-trough
    e1, e2 = np.array([[1.0, 0.2, 0]]), np.array([[0.6, -0.1, 0]])
    t = np.linspace(0, 0.1, 200001)
    ex = env.waveform(e1, e2, 2000.0, 2010.0, t)[0, :, 0]
    # envelope of the beat: max over the fast carrier within each ~0.5 ms window
    w = 2000  # 1 ms windows: two carrier cycles
    peaks = np.abs(ex[: len(ex) // w * w]).reshape(-1, w).max(axis=1)
    assert np.isclose(peaks.max() - peaks.min(), env.directional_envelope(e1, e2, [1, 0, 0])[0], rtol=0.02)


# --- n carriers (roadmap item 8): the exact phase-free envelope and phase-locked waveforms ---

def test_multi_pair_envelope_reduces_to_two_pair_formula():
    rng = np.random.default_rng(1)
    e1, e2 = rng.normal(size=(50, 3)), rng.normal(size=(50, 3))
    n = env._unit(rng.normal(size=(50, 3)))
    assert np.allclose(env.multi_pair_envelope([e1, e2], n), env.directional_envelope(e1, e2, n))
    assert np.allclose(env.max_multi_pair_envelope([e1, e2], n_dirs=6000), env.max_envelope(e1, e2), rtol=0.02)


def test_multi_pair_envelope_three_equal_parallel_carriers():
    e = np.array([[1.0, 0, 0]])
    # three equal phasors: |sum| ranges from 0 to 3 as the beats drift, so the envelope is 3
    assert np.isclose(env.multi_pair_envelope([e, e, e], [1, 0, 0]), 3.0)
    # one dominant carrier: 4 + 1 + 1 ranges from 2 to 6, envelope 4 = 2 x (sum of the others)
    assert np.isclose(env.multi_pair_envelope([4 * e, e, e], [1, 0, 0]), 4.0)
    # Grossman's 2 min|E_i| is a lower bound, equal only for two carriers
    assert env.multi_pair_envelope([4 * e, e, e], [1, 0, 0]) >= env.multi_pair_lower_bound(4 * e, e, e)[0]
    assert np.isclose(env.multi_pair_envelope([e, e, e], [1, 0, 0]), 3.0) and env.multi_pair_lower_bound(e, e, e)[0] == 2.0


def test_multi_tone_waveform_envelope_matches_phase_free_formula_for_incommensurate_beats():
    e = np.array([[1.0, 0, 0]])
    fields = [2 * e, e, 0.5 * e]
    f = [1000.0, 1010.0, 1010.0 + 10 * np.sqrt(2)]  # incommensurate beats: every relative phase occurs
    t = np.linspace(0, 20.0, 400001)  # 20 s covers the slow beat
    hi, lo = env.multi_tone_envelope_bounds(fields, f, t, n=[1, 0, 0])
    assert np.isclose(hi, 3.5, atol=0.02) and np.isclose(lo, 0.5, atol=0.02)
    assert np.isclose(hi - lo, env.multi_pair_envelope(fields, [1, 0, 0])[0], atol=0.03)


def test_phase_locked_commensurate_carriers_can_have_smaller_envelope():
    e = np.array([[1.0, 0, 0]])
    fields = [e, e, e]
    f = [1000.0, 1010.0, 1020.0]  # equally spaced: the beat is periodic and phases are locked
    t = np.linspace(0, 0.2, 200001)
    hi, lo = env.multi_tone_envelope_bounds(fields, f, t, n=[1, 0, 0], phases=[0, 0, 0])
    assert np.isclose(hi, 3.0, atol=0.01)
    # with a half-turn on the middle carrier the sum never vanishes: the phase-free formula is an upper bound
    hi2, lo2 = env.multi_tone_envelope_bounds(fields, f, t, n=[1, 0, 0], phases=[0, np.pi / 2, 0])
    assert hi2 - lo2 <= env.multi_pair_envelope(fields, [1, 0, 0])[0] + 1e-9
    assert lo2 > lo + 0.05
