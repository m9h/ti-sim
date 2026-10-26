"""10-10 electrode positions constructed on a skin surface from the four fiducials.

Charm heads carry a cap file; custom meshes (MIDA, spheres) do not. The 10-10 system is
defined geodesically (Oostenveld & Praamstra 2001): the midline arc nasion -> inion carries
Fpz/Fz/Cz/Pz/Oz at 10/30/50/70/90 %, the coronal arc LPA -> RPA through Cz carries
T7/C3/Cz/C4/T8 at 10/30/50/70/90 %, the circumferential ring through Fpz, T7, Oz, T8 carries
Fp1/AF7/F7/FT7 at 20/40/60/80 % of each quarter, and the intermediate coronal rows (F7-Fz-F8
etc.) carry the 5/3/1 positions at 12.5/25/37.5 %. Arcs are plane-surface intersections of
the triangulated skin, parameterised by chord length; the midline and coronal planes
contain the head centre so the arcs run over the vertex.

`ten_ten(nodes, triangles, fiducials)` -> {name: xyz} in the mesh units. Accuracy against
MNI-warped cap files is the few-millimetre level the literature reports for geodesic vs
template caps; against a real digitiser it is the usual 10-10 placement uncertainty.
"""

import numpy as np

MIDLINE = {"Fpz": 0.1, "AFz": 0.2, "Fz": 0.3, "FCz": 0.4, "Cz": 0.5, "CPz": 0.6, "Pz": 0.7, "POz": 0.8, "Oz": 0.9}
CORONAL = {"T7": 0.1, "C5": 0.2, "C3": 0.3, "C1": 0.4, "C2": 0.6, "C4": 0.7, "C6": 0.8, "T8": 0.9}
RING = {  # quarter: names at 20/40/60/80 %
    ("Fpz", "T7"): ("Fp1", "AF7", "F7", "FT7"), ("T7", "Oz"): ("TP7", "P7", "PO7", "O1"),
    ("Fpz", "T8"): ("Fp2", "AF8", "F8", "FT8"), ("T8", "Oz"): ("TP8", "P8", "PO8", "O2"),
}
ROWS = {  # (left end, midline, right end): names at 12.5/25/37.5 % and mirrored
    ("AF7", "AFz", "AF8"): (None, "AF3", None, None, "AF4", None),
    ("F7", "Fz", "F8"): ("F5", "F3", "F1", "F2", "F4", "F6"),
    ("FT7", "FCz", "FT8"): ("FC5", "FC3", "FC1", "FC2", "FC4", "FC6"),
    ("TP7", "CPz", "TP8"): ("CP5", "CP3", "CP1", "CP2", "CP4", "CP6"),
    ("P7", "Pz", "P8"): ("P5", "P3", "P1", "P2", "P4", "P6"),
    ("PO7", "POz", "PO8"): (None, "PO3", None, None, "PO4", None),
}


class SkinSurface:
    def __init__(self, nodes, triangles):
        self.p = np.asarray(nodes, float)
        self.t = np.asarray(triangles, int)
        if self.t.min() == 1 and self.t.max() == len(self.p):
            self.t = self.t - 1
        self.centre = self.p.mean(axis=0)

    @classmethod
    def from_simnibs(cls, msh, tag=1005):
        s = msh.crop_mesh(tags=[tag])
        tri = s.elm.node_number_list[s.elm.elm_type == 2][:, :3] - 1
        return cls(s.nodes.node_coord, tri)

    def nearest(self, x):
        return self.p[np.argmin(np.linalg.norm(self.p - np.asarray(x, float), axis=1))]

    def plane_loops(self, point, normal):
        """[(polyline (n, 3), closed)] of the surface cut by a plane; chains are walked
        both ways from their first segment, so open cuts (a mesh with a hole) come out whole."""
        n = np.asarray(normal, float) / np.linalg.norm(normal)
        d = (self.p - np.asarray(point, float)) @ n
        d[d == 0] = 1e-9 * (np.abs(d).max() or 1.0)  # vertices on the plane count as one side
        tri = self.t
        sd = d[tri]
        cut = (sd.max(axis=1) > 0) & (sd.min(axis=1) < 0)
        segs = {}  # edge key -> crossing point; segments as pairs of edge keys
        pairs = []
        for a, b, c in tri[cut]:
            pts = []
            for i, j in ((a, b), (b, c), (c, a)):
                if (d[i] > 0) != (d[j] > 0):
                    key = (min(i, j), max(i, j))
                    if key not in segs:
                        w = d[i] / (d[i] - d[j])
                        segs[key] = self.p[i] + w * (self.p[j] - self.p[i])
                    pts.append(key)
            if len(pts) == 2:
                pairs.append(tuple(pts))
        adj = {}
        for k1, k2 in pairs:
            adj.setdefault(k1, []).append(k2)
            adj.setdefault(k2, []).append(k1)
        seen, loops = set(), []
        for start in adj:
            if start in seen:
                continue
            path, prev, cur, closed = [start], None, start, False
            seen.add(start)
            while True:  # walk forward until the chain closes or ends
                nxt = [k for k in adj[cur] if k != prev and k not in seen]
                if not nxt:
                    closed = start in adj[cur] and len(path) > 2
                    break
                prev, cur = cur, nxt[0]
                seen.add(cur)
                path.append(cur)
            if not closed:  # open chain: walk the other way from the start as well
                prev, cur, back = (path[1] if len(path) > 1 else None), start, []
                while True:
                    nxt = [k for k in adj[cur] if k != prev and k not in seen]
                    if not nxt:
                        break
                    prev, cur = cur, nxt[0]
                    seen.add(cur)
                    back.append(cur)
                path = back[::-1] + path
            if len(path) >= 3:
                loops.append((np.array([segs[k] for k in path]), closed))
        return loops

    def arc(self, a, b, third, prefer="up", up=None):
        """Polyline from a to b along the plane through a, b and `third`; of the two
        arcs of the loop, prefer "up" (mean position along `up`) or "short"."""
        a, b = self.nearest(a), self.nearest(b)
        normal = np.cross(b - a, np.asarray(third, float) - a)
        loops = self.plane_loops(a, normal)
        best = None
        for loop, closed in loops:
            ia, ib = np.argmin(np.linalg.norm(loop - a, axis=1)), np.argmin(np.linalg.norm(loop - b, axis=1))
            da, db = np.linalg.norm(loop[ia] - a), np.linalg.norm(loop[ib] - b)
            if best is None or da + db < best[0]:
                best = (da + db, loop, closed, ia, ib)
        _, loop, closed, ia, ib = best
        if ia == ib:
            raise ValueError("arc end points coincide on the cut")
        if not closed:
            lo, hi = min(ia, ib), max(ia, ib)
            seg = loop[lo:hi + 1]
            return seg if ia < ib else seg[::-1]
        fwd = np.vstack([loop[ia:], loop[:ib + 1]]) if ib < ia else loop[ia:ib + 1]
        back = np.vstack([loop[ib:], loop[:ia + 1]]) if ia < ib else loop[ib:ia + 1]
        back = back[::-1]
        if prefer == "short":
            return fwd if _length(fwd) <= _length(back) else back
        u = np.asarray(up, float) if up is not None else np.array([0, 0, 1.0])
        return fwd if (fwd.mean(axis=0) - self.centre) @ u >= (back.mean(axis=0) - self.centre) @ u else back


def _length(poly):
    return np.linalg.norm(np.diff(poly, axis=0), axis=1).sum()


def along(poly, fraction):
    """Point at `fraction` of the chord length of a polyline."""
    seg = np.linalg.norm(np.diff(poly, axis=0), axis=1)
    cum = np.concatenate([[0.0], np.cumsum(seg)])
    s = fraction * cum[-1]
    i = min(np.searchsorted(cum, s, side="right") - 1, len(seg) - 1)
    w = (s - cum[i]) / seg[i] if seg[i] > 0 else 0.0
    return poly[i] + w * (poly[i + 1] - poly[i])


def head_frame(fid):
    """(origin, x right, y anterior, z up) from Nz, Iz, LPA, RPA."""
    lpa, rpa, nz = (np.asarray(fid[k], float) for k in ("LPA", "RPA", "Nz"))
    origin = (lpa + rpa) / 2
    x = rpa - lpa
    x /= np.linalg.norm(x)
    y = nz - origin
    y -= (y @ x) * x
    y /= np.linalg.norm(y)
    z = np.cross(x, y)
    return origin, x, y, z


def ten_ten(nodes, triangles, fiducials):
    """{name: xyz} for the 10-10 positions reachable from the fiducials (Nz, Iz, LPA, RPA)."""
    skin = SkinSurface(nodes, triangles)
    fid = {k: skin.nearest(v) for k, v in fiducials.items()}
    origin, x, y, z = head_frame(fid)
    top = origin + 200.0 * z  # a point above the head: midline/coronal planes contain the up axis
    pos = {"Nz": fid["Nz"], "Iz": fid["Iz"], "LPA": fid["LPA"], "RPA": fid["RPA"]}
    mid = skin.arc(fid["Nz"], fid["Iz"], top, up=z)
    for name, f in MIDLINE.items():
        pos[name] = along(mid, f)
    cor = skin.arc(fid["LPA"], fid["RPA"], pos["Cz"], up=z)
    for name, f in CORONAL.items():
        pos[name] = along(cor, f)
    # ring: least-squares plane through Fpz, T7, Oz, T8; each quarter is the short arc
    ring_pts = np.array([pos["Fpz"], pos["T7"], pos["Oz"], pos["T8"]])
    c = ring_pts.mean(axis=0)  # the ring plane passes through the two anchors and the centroid
    for (a, b), names in RING.items():
        arc = skin.arc(pos[a], pos[b], c, prefer="short")
        for name, f in zip(names, (0.2, 0.4, 0.6, 0.8)):
            pos[name] = along(arc, f)
    for (left, midname, right), names in ROWS.items():
        arc = skin.arc(pos[left], pos[right], pos[midname], up=z)
        for name, f in zip(names, (0.125, 0.25, 0.375, 0.625, 0.75, 0.875)):
            if name:
                pos[name] = along(arc, f)
    return {k: np.asarray(v, float) for k, v in pos.items()}


def read_csv_positions(fn):
    """{name: xyz} from a SimNIBS positions csv (Electrode/ReferenceElectrode/Fiducial rows)."""
    import csv

    out = {}
    with open(fn) as f:
        for r in csv.reader(f):
            if len(r) >= 5 and r[0] in ("Electrode", "ReferenceElectrode", "Fiducial"):
                out[r[4]] = np.array(r[1:4], float)
    return out
