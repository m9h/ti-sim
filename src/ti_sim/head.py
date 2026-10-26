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

Current normalisation. SimNIBS solves the tDCS problem with 0/1 V on the two electrode tops
and rescales to the requested current using an estimate of the flux of the element-wise J
through the skin-electrode interface (fem.tdcs, _calc_flux_electrodes). That estimate is
biased on coarse electrode meshes: +6 % on a 262k-tet gmsh sphere, +9 to +13 % on ernie
with 2 mm pads (found 2026-10-26 against a DOLFINx solve and the analytic 4-shell sphere).
The exact current of the discrete solution follows from the energy identity
int sigma |grad v|^2 dV = dV * I_h, with dV the potential difference between the electrode
tops, so `TwoPairSimulation` saves v and J, computes I_h per pair and rescales E by
I_nominal / I_h. `calibration` records both currents and the pair resistance dV / I_h.
"""

import hashlib
import json
import os

import numpy as np

from . import envelope, tissues

WM, GM, CSF, BONE, SCALP = 1, 2, 3, 4, 5


def _key(parts):
    return hashlib.sha1(json.dumps(parts, sort_keys=True, default=str).encode()).hexdigest()[:12]


def discrete_current(m, i_nominal):
    """Exact current of a SimNIBS tDCS solution (fields v, E, J on the full mesh with its
    electrodes) from the energy identity P = int J . E dV = dV * I_h, with dV the potential
    difference between the two electrode tops (SimNIBS Dirichlet surfaces 2100+). Returns
    the nominal and discrete currents (A), the factor i_nominal / I_h that rescales SimNIBS's
    fields to the requested current, and the pair resistance dV / I_h (ohm)."""
    v = m.field["v"].value
    tri = m.elm.elm_type == 2
    tops = [t for t in np.unique(m.elm.tag1[tri]) if 2100 <= t < 2200]
    if len(tops) != 2:
        raise RuntimeError(f"expected two electrode top surfaces (tags 2100+), found {tops}")
    means = [v[np.unique(m.elm.node_number_list[tri & (m.elm.tag1 == t)][:, :3]) - 1].mean() for t in tops]
    d_v = abs(means[0] - means[1])
    tets = m.elm.elm_type == 4
    vol = m.elements_volumes_and_areas().value[tets]
    power = float((np.einsum("ij,ij->i", m.field["J"].value[tets], m.field["E"].value[tets]) * vol).sum() * 1e-9)
    i_h = power / d_v
    return {"i_nominal_A": float(i_nominal), "i_discrete_A": float(i_h), "factor": float(abs(i_nominal) / i_h),
            "resistance_ohm": float(d_v / i_h), "simnibs_bias_percent": float((i_h / abs(i_nominal) - 1) * 100)}


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
                    self.extra_cond, subj, "calibration=discrete"])
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
        s.pathfem, s.open_in_gmsh, s.fields = self.workdir, False, "veEJ"
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
        full = [mesh_io.read_msh(f) for f in self.files]
        self.calibration = [discrete_current(m, i) for m, i in zip(full, self.currents)]
        ms = [m.crop_mesh(tags=keep) for m in full]
        del full
        tets = ms[0].elm.elm_type == 4
        if not np.array_equal(ms[0].elm.tag1, ms[1].elm.tag1):
            raise RuntimeError("the two SimNIBS runs produced different meshes")
        self.mesh = ms[0]
        self.tags = ms[0].elm.tag1[tets]
        self.e1 = ms[0].field["E"].value[tets] * self.calibration[0]["factor"]
        self.e2 = ms[1].field["E"].value[tets] * self.calibration[1]["factor"]
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

    def fields_at(self, points):
        """(E1, E2, E_AM,max) at points (n, 3) mm: the containing tetrahedron's values, zero
        outside the mesh (SimNIBS's locator; fine for ROI probes, not for sub-element detail)."""
        p = np.asarray(points, float).reshape(-1, 3)
        idx = self.mesh.find_tetrahedron_with_points(p, compute_baricentric=False)
        tet_index = np.flatnonzero(self.mesh.elm.elm_type == 4)
        pos = np.full(len(p), -1)
        ok = idx > 0
        lookup = np.full(self.mesh.elm.nr + 1, -1)
        lookup[tet_index + 1] = np.arange(len(tet_index))
        pos[ok] = lookup[idx[ok]]
        good = pos >= 0
        e1, e2, am = np.zeros((len(p), 3)), np.zeros((len(p), 3)), np.zeros(len(p))
        e1[good], e2[good], am[good] = self.e1[pos[good]], self.e2[pos[good]], self.am_max[pos[good]]
        return e1, e2, am

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
