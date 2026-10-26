"""The MIDA name -> tissue mapping covers the published Table 1, and the label-file parser
reads both layouts. The model itself is not redistributable, so no data test here; the
mesh route is exercised on a synthetic sphere in test_labelvol.py."""

import numpy as np
import pytest

from ti_sim import mida


def test_every_table1_structure_is_classified():
    names = {i + 1: n for i, n in enumerate(mida.TABLE1)}
    assert len(names) == 115
    m = mida.classify(names)
    assert set(m) == set(names)
    inv = {}
    for k, v in m.items():
        inv.setdefault(v, []).append(names[k])
    assert inv[mida.HIPPOCAMPUS] == ["Hippocampus"]
    assert set(inv[0]) == {n for n in mida.TABLE1 if n.startswith("Air Internal") or "Auditory Canal" in n or "Pharyngotympanic" in n}
    assert "Brain White Matter" in inv[1] and "Brain Gray Matter" in inv[2]
    assert all(n.startswith("Muscle") or n in ("Tongue",) or n.startswith("Tendon") for n in inv[10])
    assert "Skull Diploe" in inv[8] and "Skull" in inv[4] and "Teeth" in inv[7]
    assert set(m.values()) <= set(range(0, 13)) | set(mida.EXTRA_COND)


def test_download_numbering_is_classified():
    m = mida.classify(mida.MIDA_V1_LABELS)
    assert len(m) == 116 and m[50] == 0 and m[5] == mida.HIPPOCAMPUS and m[116] == mida.THALAMUS
    assert m[12] == 1 and m[10] == 2 and m[52] == 8 and m[40] == 4 and m[1] == 3
    assert m[19] == 2 and m[88] == 5  # hypophysis is brain tissue, the salivary glands are scalp-like
    assert sum(v == 0 for v in m.values()) == 10  # background + 7 sinuses/cavities + 2 ear canals


def test_unknown_name_raises():
    with pytest.raises(ValueError, match="no tissue rule"):
        mida.classify({1: "Hippocampus", 2: "Flux Capacitor"})


def test_parse_label_file_both_layouts(tmp_path):
    names = {i + 1: n for i, n in enumerate(mida.TABLE1)}
    f1 = tmp_path / "a.txt"
    f1.write_text("# MIDA labels\n" + "\n".join(f"{k}\t{v}" for k, v in names.items()) + "\n")
    f2 = tmp_path / "b.txt"
    f2.write_text("\n".join(f"{v} {k}" for k, v in names.items()) + "\n")
    assert mida.parse_label_file(f1) == names
    assert mida.parse_label_file(f2) == names
    f3 = tmp_path / "c.txt"
    f3.write_text("1 one\n2 two\n")
    with pytest.raises(ValueError, match="not a MIDA label file"):
        mida.parse_label_file(f3)


def test_relabel_crop_downsample():
    labels = np.zeros((8, 8, 8), np.uint16)
    labels[2:6, 2:6, 2:6] = 7
    labels[3:5, 3:5, 3:5] = 50
    aff = np.diag([0.5, 0.5, 0.5, 1.0])
    aff[:3, 3] = -2.0
    v = mida.relabel(labels, {7: 5, 50: mida.HIPPOCAMPUS})
    assert v[3, 3, 3] == mida.HIPPOCAMPUS and v[2, 2, 2] == 5 and v[0, 0, 0] == 0
    c, a = mida.crop(v, aff, z_min_mm=-1.0)  # z of slice k is -2 + 0.5 k -> keep k >= 2
    assert c.shape[2] == 6 and np.isclose(a[2, 3], -1.0)
    d, a2 = mida.downsample(c, a, 2)
    assert d.shape == (4, 4, 3) and np.isclose(a2[0, 0], 1.0)
