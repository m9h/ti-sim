"""TI on a SimNIBS head model (charm `m2m_<subject>` folder, or any custom .msh).

`TwoPairSimulation` runs SimNIBS once per electrode pair (one Laplace solve each, the two
carriers see the same volume conductor) and caches the result meshes in `workdir` under a
key made from the montage, the currents, the electrode size and the conductivity table.
The envelope is then evaluated per tetrahedron with `ti_sim.envelope` (same numbers as
SimNIBS's TI_utils; see tests). ROI statistics follow what the papers report: median,
mean, 99th percentile of E_AM,max in a tissue, in a sphere, or in a label volume (MIDA
hippocampus etc.), plus the unmodulated |E1|, |E2| and the tACS sum for the same montage.

Electrodes are circles of `radius` on the scalp at 10-10 names or subject coordinates (mm),
as in SimNIBS's TI.py example. Everything stays in the SimNIBS mm frame; results are V/m
at the currents given (amperes).
"""

import hashlib
import json
import os

import numpy as np

from . import envelope, tissues

WM, GM, CSF, BONE, SCALP = 1, 2, 3, 4, 5


def _key(parts):
    return hashlib.sha1(json.dumps(parts, sort_keys=True, default=str).encode()).hexdigest()[:12]


class TwoPairSimulation:
    """SimNIBS tDCS solves for pairs ((A1, B1), (A2, B2)) with currents (I1, I2) in A."""

    def __init__(self, pairs, currents, workdir, m2m_dir=None, fnamehead=None, radius=10.0,
                 thickness=2.0, cond=None, cond_label=None, extra_cond=None, threads=4):
        if (m2m_dir is None) == (fnamehead is None):
            raise ValueError("give exactly one of m2m_dir (charm folder) or fnamehead (custom .msh)")
        os.environ.setdefault("OMP_NUM_THREADS", str(threads))
        self.pairs, self.currents = [tuple(p) for p in pairs], tuple(float(c) for c in currents)
        self.radius, self.thickness = float(radius), float(thickness)
        self.cond = dict(cond or tissues.SIMNIBS)
        self.cond_label = cond_label or ("simnibs" if cond is None else "custom")
        self.extra_cond = dict(extra_cond or {})  # {tag: sigma} beyond the 12 SimNIBS tissues
        self.m2m_dir, self.fnamehead = m2m_dir, fnamehead
        subj = (os.path.basename(os.path.normpath(m2m_dir)).removeprefix("m2m_") if m2m_dir
                else os.path.splitext(os.path.basename(fnamehead))[0])
        self.subject = subj
        key = _key([self.pairs, self.currents, self.radius, self.thickness, self.cond,
                    self.extra_cond, subj])
        self.workdir = os.path.join(str(workdir), f"ti_{subj}_{key}")
        self.files = [os.path.join(self.workdir, f"{subj}_TDCS_{k}_scalar.msh") for k in (1, 2)]
        if not all(os.path.exists(f) for f in self.files):
            self._run()
        self._load()

    def _run(self):
        from simnibs import run_simnibs, sim_struct

        os.makedirs(os.path.dirname(self.workdir), exist_ok=True)
        s = sim_struct.SESSION()
        if self.m2m_dir:
            s.subpath = self.m2m_dir
        else:
            s.fnamehead = self.fnamehead
        s.pathfem, s.open_in_gmsh, s.fields = self.workdir, False, "eE"
        for (a, b), i in zip(self.pairs, self.currents):
            t = s.add_tdcslist()
            t.currents = [i, -i]
            tissues.apply_to_tdcslist(t, self.cond)
            for tag, sigma in self.extra_cond.items():
                t.cond[tag - 1].value, t.cond[tag - 1].name = float(sigma), f"tag{tag}"
            for nr, pos in enumerate((a, b), 1):
                e = t.add_electrode()
                e.channelnr, e.shape = nr, "ellipse"
                e.centre = list(map(float, pos)) if not isinstance(pos, str) else pos
                e.dimensions, e.thickness = [2 * self.radius] * 2, [self.thickness]
        run_simnibs(s)

    def _load(self):
        from simnibs import ElementTags, mesh_io

        keep = np.hstack((np.arange(ElementTags.TH_START, ElementTags.SALINE_START - 1),
                          np.arange(ElementTags.TH_SURFACE_START, ElementTags.SALINE_TH_SURFACE_START - 1)))
        ms = [mesh_io.read_msh(f).crop_mesh(tags=keep) for f in self.files]
        tets = ms[0].elm.elm_type == 4
        if not np.array_equal(ms[0].elm.tag1, ms[1].elm.tag1):
            raise RuntimeError("the two SimNIBS runs produced different meshes")
        self.mesh = ms[0]
        self.tags = ms[0].elm.tag1[tets]
        self.e1 = ms[0].field["E"].value[tets]
        self.e2 = ms[1].field["E"].value[tets]
        self.centres = ms[0].elements_baricenters().value[tets]
        self.volumes = ms[0].elements_volumes_and_areas().value[tets]
        self.am_max = envelope.max_envelope(self.e1, self.e2)
        self.tacs = envelope.tacs_sum(self.e1, self.e2)

    # --- regions -------------------------------------------------------------------------

    def in_tags(self, tags):
        return np.isin(self.tags, np.atleast_1d(tags))

    def in_sphere(self, centre, radius, tags=None):
        sel = np.linalg.norm(self.centres - np.asarray(centre, float), axis=1) <= radius
        return sel & self.in_tags(tags) if tags is not None else sel

    def stats(self, mask, field=None):
        """Volume-weighted median, mean, 99th percentile and max of a per-element field
        (default E_AM,max) over `mask`, with |E1|, |E2| and tACS for comparison."""
        f = self.am_max if field is None else field
        v = self.volumes[mask]
        out = {"n_elements": int(mask.sum()), "volume_mm3": float(v.sum())}
        for name, x in (("am", f[mask]), ("e1", np.linalg.norm(self.e1[mask], axis=1)),
                        ("e2", np.linalg.norm(self.e2[mask], axis=1)), ("tacs", self.tacs[mask])):
            if len(x) == 0:
                continue
            order = np.argsort(x)
            cw = np.cumsum(v[order]) / v.sum()
            out[name] = {"median": float(x[order][np.searchsorted(cw, 0.5)]),
                         "mean": float(np.average(x, weights=v)),
                         "p99": float(x[order][min(np.searchsorted(cw, 0.99), len(x) - 1)]),
                         "max": float(x.max())}
        return out

    def focality(self, threshold, tags=(GM, WM)):
        """Brain volume (mm^3) with E_AM,max above `threshold`, and the same for tACS."""
        m = self.in_tags(tags)
        return {"threshold": threshold,
                "am_mm3": float(self.volumes[m & (self.am_max > threshold)].sum()),
                "tacs_mm3": float(self.volumes[m & (self.tacs > threshold)].sum())}

    def write_msh(self, fn):
        """Mesh with |E1|, |E2|, E_AM,max and tACS element fields for gmsh."""
        from simnibs import mesh_io

        m = self.mesh
        m.elmdata = []
        tets = m.elm.elm_type == 4
        for name, x in (("magnE1", np.linalg.norm(self.e1, axis=1)),
                        ("magnE2", np.linalg.norm(self.e2, axis=1)), ("TImax", self.am_max),
                        ("tACS", self.tacs)):
            full = np.zeros(m.elm.nr)
            full[tets] = x
            m.add_element_field(full, name)
        mesh_io.write_msh(m, fn)
        return fn


def eeg_positions(m2m_dir, cap="EEG10-10_UI_Jurak_2007.csv"):
    """{name: xyz mm} from the charm folder's EEG cap file."""
    import csv

    with open(os.path.join(m2m_dir, "eeg_positions", cap)) as f:
        return {r[-1]: np.array(r[1:4], float) for r in csv.reader(f) if len(r) >= 5}
