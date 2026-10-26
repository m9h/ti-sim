"""The saline disc: FEM against the analytic line-source field, and Grossman's geometry."""

import numpy as np

from ti_sim import phantom


def _analytic_disc_field(pts, a, b, i_per_h, sigma, radius):
    """E of line sources +-I/h at a, b inside an insulating disc of radius R: free-space
    terms plus their Kelvin images at R^2 a/|a|^2 (Neumann image for a disc; the uniform
    term from the images cancels for a zero-net source pair)."""
    def one(p, c, sign):
        d = p - c
        e = sign * i_per_h / (2 * np.pi * sigma) * d / (d * d).sum(-1, keepdims=True)
        ci = c * radius**2 / (c @ c)
        di = p - ci
        e += sign * i_per_h / (2 * np.pi * sigma) * di / (di * di).sum(-1, keepdims=True)
        return e

    return one(pts, a, +1) + one(pts, b, -1)


def test_fem_matches_analytic_line_sources():
    R, sigma, iph = 25e-3, 0.333, 0.1
    a, b = np.array([-0.6, 0.3]) * R, np.array([0.5, -0.4]) * R  # interior sources
    m = phantom.disc_mesh(R, nrefs=6, refine_near=[a, b], refine_rounds=2)
    f = phantom.DiscField(m, sigma, a, b, iph, wire_sigma=0.3e-3)
    rng = np.random.default_rng(0)
    pts = rng.uniform(-0.6 * R, 0.6 * R, (300, 2))
    far = (np.linalg.norm(pts - a, axis=1) > 4e-3) & (np.linalg.norm(pts - b, axis=1) > 4e-3)
    pts = pts[far]
    ef, ea = f.field(pts), _analytic_disc_field(pts, a, b, iph, sigma, R)
    err = np.linalg.norm(ef - ea, axis=1) / np.linalg.norm(ea, axis=1)
    assert np.median(err) < 0.03 and np.percentile(err, 90) < 0.08


def test_trapezoid_geometry():
    (al, bl), (ar, br) = phantom.trapezoid_electrodes(0.39, 1.96, radius=1.0)
    for p in (al, bl, ar, br):
        assert np.isclose(np.linalg.norm(p), 1.0)  # all on the rim
    assert np.isclose(ar[0] - al[0], 0.39) and np.isclose(br[0] - bl[0], 1.96)
    assert al[1] < 0 < bl[1]  # the short base (close wires) is at the bottom, as in the paper


def test_envelope_map_symmetry_and_peak(tmp_path):
    pairs = phantom.trapezoid_electrodes(1.96, 1.96)
    ti = phantom.DiscTI(pairs, nrefs=5, refine_rounds=1)
    em = ti.envelope_map(n=81)
    # mirror symmetry in x for equal currents: E_AM,y(x) == E_AM,y(-x)
    ok = np.isfinite(em.am_y) & np.isfinite(em.am_y[:, ::-1])
    assert np.allclose(em.am_y[ok], em.am_y[:, ::-1][ok], rtol=0.05, atol=1e-3)
    v, (px, py), (hx, hy) = em.peak("am_y", exclude=2e-3)
    assert 0.05 < v < 5.0  # Grossman: 0.21 V/m simulated, 0.23 measured; see the example
    assert abs(px) < 0.1 and abs(py) < 0.15  # centred, as in the paper
    assert 0.3 < hx < 0.6 and 0.35 < hy < 0.65  # paper HWHM (0.42, 0.49) R
    vx, (px, py), _ = em.peak("am_x", exclude=2e-3)
    assert abs(px) < 0.1 and -0.75 < py < -0.45  # paper: (0, -0.59) R
