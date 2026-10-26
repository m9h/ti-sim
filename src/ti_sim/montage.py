"""Montage search on a SimNIBS lead field: Lee et al. 2020 (Sci Rep 10:11730) style
exhaustive search over two electrode pairs and their current ratio, scored by the envelope
in a target region.

A lead field (one FEM solve per electrode against a reference, SimNIBS TDCSLEADFIELD) gives
the field of any pair by subtraction, so thousands of montages cost one linear algebra
pass each. `compute_leadfield` makes a tetrahedral (not surface) lead field in GM+WM so deep
targets can be scored; `Leadfield` wraps SimNIBS's loader; `search_two_pairs` scores every
combination of two disjoint pairs from an electrode list at the given current ratios by the
median E_AM,max over the ROI elements, then re-scores the best on the whole brain for
focality. Huang & Parra 2020 showed the analytic optimum of modulation depth is the HD-tES
solution; this search is the plain enumeration the TI papers actually use.
"""

import itertools
import os

import numpy as np

from . import envelope
from .head import GM, WM


def compute_leadfield(m2m_dir, workdir, cap="EEG10-10_UI_Jurak_2007.csv", radius=10.0,
                      tissues=(GM, WM), threads=4):
    """Run (or reuse) a SimNIBS tDCS lead field with circular electrodes of `radius` mm on
    every cap position; tetrahedral, in `tissues`. Returns the hdf5 path."""
    from simnibs import ElementTags, run_simnibs, sim_struct

    os.environ.setdefault("OMP_NUM_THREADS", str(threads))
    subj = os.path.basename(os.path.normpath(m2m_dir)).removeprefix("m2m_")
    out = os.path.join(str(workdir), f"leadfield_{subj}_r{radius:g}_{'-'.join(map(str, tissues))}")
    fn = os.path.join(out, f"{subj}_leadfield_{os.path.splitext(cap)[0]}.hdf5")
    if os.path.exists(fn):
        return fn
    lf = sim_struct.TDCSLEADFIELD()
    lf.subpath, lf.pathfem, lf.eeg_cap = m2m_dir, out, cap
    lf.interpolation = None  # tetrahedra, not the middle-GM surface
    lf.tissues = [ElementTags(t) for t in tissues]
    lf.electrode.shape, lf.electrode.dimensions, lf.electrode.thickness = "ellipse", [2 * radius] * 2, [2]
    lf.field = "E"
    run_simnibs(lf)
    if not os.path.exists(fn):
        cands = [f for f in os.listdir(out) if f.endswith(".hdf5")]
        if not cands:
            raise FileNotFoundError(f"no lead field written under {out}")
        fn = os.path.join(out, cands[0])
    return fn


class Leadfield:
    def __init__(self, fn):
        from simnibs.utils import TI_utils

        self.fn = fn
        self.lf, self.mesh, self.idx = TI_utils.load_leadfield(fn)
        self.names = [n for n in self.idx]
        self.tets = self.mesh.elm.elm_type == 4
        if self.lf.shape[1] == self.tets.sum():  # element-wise lead field
            self.centres = self.mesh.elements_baricenters().value[self.tets]
            self.volumes = self.mesh.elements_volumes_and_areas().value[self.tets]
            self.tags = self.mesh.elm.tag1[self.tets]
        else:  # node-wise (surface) lead field
            self.centres = self.mesh.nodes.node_coord
            self.volumes = np.ones(self.lf.shape[1])
            self.tags = np.full(self.lf.shape[1], GM)

    def field(self, a, b, current):
        """E (M, 3) in V/m for `current` A into a and out of b."""
        from simnibs.utils import TI_utils

        return TI_utils.get_field([a, b, current], self.lf, self.idx)

    def sphere(self, centre, radius):
        return np.linalg.norm(self.centres - np.asarray(centre, float), axis=1) <= radius


def search_two_pairs(lf, roi, electrodes=None, ratios=(1.0,), total_current=2e-3, top=10,
                     threshold=None, min_separation_mm=0.0):
    """Enumerate montages (pair1, pair2, I1, I2) with pair1 and pair2 disjoint, from
    `electrodes` (default: all), currents I1 = total / (1 + r), I2 = r I1 for r in `ratios`
    and 1/r. Score = median E_AM,max over `roi` (bool mask on lead-field points). The `top`
    montages are re-scored over GM+WM: p99 outside the ROI and, with `threshold` (V/m), the
    brain volume above it inside and outside the ROI (mm^3)."""
    names = list(electrodes) if electrodes is not None else [n for n in lf.idx if lf.idx[n] is not None]
    names = [n for n in names if n in lf.idx]  # the reference electrode may be named explicitly
    pairs = list(itertools.combinations(names, 2))
    roi = np.asarray(roi, bool)
    f_roi = np.stack([lf.field(a, b, 1.0)[roi] for a, b in pairs])  # (P, M, 3) at 1 A
    rset = sorted({float(r) for r in ratios} | {1 / float(r) for r in ratios})
    results = []
    for p, (a1, b1) in enumerate(pairs):
        others = [q for q in range(p + 1, len(pairs)) if not ({a1, b1} & set(pairs[q]))]
        if not others:
            continue
        f2 = f_roi[others]  # (Q, M, 3)
        for r in rset:
            i1 = total_current / (1 + r)
            i2 = r * i1
            am = envelope.max_envelope(np.broadcast_to(f_roi[p] * i1, f2.shape), f2 * i2)
            med = np.median(am, axis=1)
            for q, m in zip(others, med):
                results.append((float(m), p, q, i1, i2))
    results.sort(key=lambda x: -x[0])
    out = []
    brain = np.isin(lf.tags, (GM, WM))
    for med, p, q, i1, i2 in results[:top]:
        e1, e2 = lf.field(*pairs[p], i1), lf.field(*pairs[q], i2)
        am = envelope.max_envelope(e1, e2)
        row = {"pair1": pairs[p], "pair2": pairs[q], "I1_mA": 1e3 * i1, "I2_mA": 1e3 * i2,
               "roi_median": med, "roi_p99": float(np.percentile(am[roi], 99)),
               "brain_p99_outside_roi": float(np.percentile(am[brain & ~roi], 99)) if (brain & ~roi).any() else float("nan"),
               "brain_max": float(am[brain].max())}
        if threshold is not None:
            row["above_threshold_mm3"] = {"roi": float(lf.volumes[roi & (am > threshold)].sum()),
                                          "outside": float(lf.volumes[brain & ~roi & (am > threshold)].sum())}
        out.append(row)
    return {"n_electrodes": len(names), "n_pairs": len(pairs), "n_montages": len(results),
            "ratios": rset, "total_current_mA": 1e3 * total_current, "best": out}


def near_target(positions, target, max_distance_mm):
    """Electrode names within `max_distance_mm` of `target` (the Lee 2026 idea: restrict the
    candidates to the target's neighbourhood)."""
    t = np.asarray(target, float)
    return [n for n, p in positions.items() if np.linalg.norm(np.asarray(p) - t) <= max_distance_mm]
