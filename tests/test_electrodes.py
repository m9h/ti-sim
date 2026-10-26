"""10-10 construction from fiducials: exact on a sphere, and within the known geodesic-vs-cap
difference on ernie's skin against the MNI-warped Jurak cap."""

import os

import numpy as np
import pytest

from ti_sim import electrodes as el

ERNIE = os.environ.get("TI_SIM_M2M", "")


def _sphere(n=40, r=1.0):
    """Unit sphere triangulation from a lat-long grid."""
    th = np.linspace(0, np.pi, n)
    ph = np.linspace(0, 2 * np.pi, 2 * n, endpoint=False)
    T, P = np.meshgrid(th, ph, indexing="ij")
    pts = np.stack([r * np.sin(T) * np.cos(P), r * np.sin(T) * np.sin(P), r * np.cos(T)], -1).reshape(-1, 3)
    tri = []
    for i in range(n - 1):
        for j in range(2 * n):
            a, b = i * 2 * n + j, i * 2 * n + (j + 1) % (2 * n)
            c, d = a + 2 * n, b + 2 * n
            tri += [(a, b, c), (b, d, c)]
    return pts, np.array(tri)


def test_sphere_positions_are_the_geodesic_angles():
    p, t = _sphere()
    fid = {"Nz": [0, 1, 0], "Iz": [0, -1, 0], "LPA": [-1, 0, 0], "RPA": [1, 0, 0]}
    pos = el.ten_ten(p, t, fid)
    assert np.allclose(pos["Cz"], [0, 0, 1], atol=0.05)
    # Fz at 30 % of the nasion-inion half circle: polar angle 90 - 180*0.3 = 36 deg from nasion
    a = np.radians(180 * 0.3)
    assert np.allclose(pos["Fz"], [0, np.cos(a), np.sin(a)], atol=0.05)
    assert np.allclose(pos["C3"], [-np.cos(np.radians(54)), 0, np.sin(np.radians(54))], atol=0.05)
    assert np.allclose(pos["T7"], [-np.cos(np.radians(18)), 0, np.sin(np.radians(18))], atol=0.05)
    # ring quarter Fpz -> T7 lies in the ring plane z = sin(18 deg); F7 at 60 % of the quarter
    assert np.allclose(pos["F7"][2], np.sin(np.radians(18)), atol=0.05)
    ang = np.degrees(np.arctan2(pos["F7"][1], -pos["F7"][0]))  # from the left (-x) towards the front (+y)
    assert abs(ang - 0.4 * 90) < 3
    # mirror symmetry and the full set
    for l, r in (("F3", "F4"), ("P7", "P8"), ("TP7", "TP8"), ("FC5", "FC6")):
        assert np.allclose(pos[l] * [-1, 1, 1], pos[r], atol=0.05)
    assert len(pos) == 4 + 9 + 8 + 16 + 28  # fiducials, midline, coronal, ring, rows


@pytest.mark.skipif(not os.path.isdir(ERNIE), reason="set TI_SIM_M2M=/path/to/m2m_ernie")
def test_ernie_against_the_cap_file():
    from simnibs import mesh_io

    m = mesh_io.read_msh(os.path.join(ERNIE, "ernie.msh"))
    skin = el.SkinSurface.from_simnibs(m)
    fid = el.read_csv_positions(os.path.join(ERNIE, "eeg_positions", "Fiducials.csv"))
    cap = el.read_csv_positions(os.path.join(ERNIE, "eeg_positions", "EEG10-10_UI_Jurak_2007.csv"))
    pos = el.ten_ten(skin.p, skin.t, {k: fid[k] for k in ("Nz", "Iz", "LPA", "RPA")})
    d = {k: np.linalg.norm(pos[k] - cap[k]) for k in pos if k in cap}
    assert len(d) >= 60
    worst = max(d, key=d.get)
    print({k: round(v, 1) for k, v in sorted(d.items(), key=lambda kv: -kv[1])[:8]})
    assert np.mean(list(d.values())) < 10.0, (np.mean(list(d.values())), worst, d[worst])
    for k in ("Cz", "Fz", "Pz", "C3", "C4", "T7", "T8", "F3", "F4", "TP7", "TP8"):
        assert d[k] < 15.0, (k, d[k])
