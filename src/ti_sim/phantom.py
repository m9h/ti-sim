"""Grossman et al. 2017, figure 2: TI in a saline disc, as a 2-D scikit-fem problem.

The phantom is a non-conductive cylinder, 50 mm in diameter and 10 mm high, filled with
saline (0.333 S/m), with 1 mm silver wire electrodes standing on the perimeter through the
full height. The wires make the field in-plane, so the cylinder is a disc of radius R with
line sources of strength I/h (A/m) and an insulating rim: -div(sigma grad phi) = s in 2-D.

The paper places the four electrodes on the corners of an isosceles trapezoid inscribed in
the disc, with bases a (bottom) and b (top) given in units of R; each electrode pair is one
left-hand or one right-hand leg of the trapezoid (figure 2A-E: a = 0.39, b = 1.96; a = 1.11,
b = 1.96; a = b = 1.96, the last with current ratios 1:1, 1:2.5 and 1:4 at a 2 mA sum).
`trapezoid_electrodes` gives the positions, `DiscTI` solves the two pairs and
`envelope_map` evaluates E_AM along x, y and the free direction on a grid, with the peak
location and half-widths the paper tabulates (in units of R).

Units: metres, amperes, V/m. The disc is meshed with scikit-fem's init_circle plus
refinement around the wires; the wires are Gaussian line sources of width `wire_sigma`.
"""

from dataclasses import dataclass

import numpy as np
from scipy.sparse.linalg import spsolve
from scipy.spatial import cKDTree
from skfem import Basis, ElementTriP1, LinearForm, MeshTri, asm, condense
from skfem.models.poisson import laplace

from .envelope import directional_envelope, max_envelope

R_DISC, H_DISC, SIGMA_SALINE, WIRE_DIAMETER = 25e-3, 10e-3, 0.333, 1e-3


def trapezoid_electrodes(a, b, radius=R_DISC):
    """Four electrode positions (m) on the rim: base a*R at the bottom (y < 0), base b*R at
    the top, the paper's orientation (figure 2A's close pair of wires is at the bottom and
    its E_AM,x peak at (0, -1)). Returns ((bottom_left, top_left), (bottom_right, top_right))
    = (pair 1, pair 2): each pair is one leg of the trapezoid."""
    ya = -np.sqrt(max(1 - (a / 2) ** 2, 0.0))
    yb = np.sqrt(max(1 - (b / 2) ** 2, 0.0))
    al, ar = np.array([-a / 2, ya]), np.array([a / 2, ya])
    bl, br = np.array([-b / 2, yb]), np.array([b / 2, yb])
    return ((al * radius, bl * radius), (ar * radius, br * radius))


def disc_mesh(radius=R_DISC, nrefs=6, refine_near=(), refine_radius=None, refine_rounds=2):
    """Unit-circle mesh scaled to `radius`, refined `refine_rounds` times within
    `refine_radius` of each point in `refine_near` (the wires)."""
    m = MeshTri.init_circle(nrefs).scaled(radius)
    rr = refine_radius if refine_radius is not None else 0.12 * radius
    for _ in range(refine_rounds):
        c = m.p[:, m.t].mean(axis=1)  # element centroids (2, nelems)
        near = np.zeros(m.t.shape[1], bool)
        for p in refine_near:
            near |= np.linalg.norm(c - np.asarray(p, float)[:, None], axis=0) < rr
        if not near.any():
            break
        m = m.refined(np.flatnonzero(near))
        rr *= 0.6
    # the refinement keeps the polygonal rim: push boundary nodes back onto the circle
    p = m.p.copy()
    bn = m.boundary_nodes()
    p[:, bn] *= radius / np.linalg.norm(p[:, bn], axis=0)
    return MeshTri(p, m.t)


def gauss2(x, c, s):
    """Unit-mass 2-D Gaussian at quadrature points x (2, nelems, nqp)."""
    r2 = ((x - np.asarray(c, float)[:, None, None]) ** 2).sum(axis=0)
    return np.exp(-r2 / (2 * s * s)) / (2 * np.pi * s * s)


class DiscField:
    """phi and E = -grad phi in the disc for a pair of line sources +-I/h at a, b."""

    def __init__(self, mesh, sigma, a, b, current_per_height, wire_sigma=0.5 * WIRE_DIAMETER):
        self.mesh, self.sigma = mesh, sigma
        self.basis = Basis(mesh, ElementTriP1())
        A = asm(laplace, self.basis) * sigma

        @LinearForm
        def rhs(v, w):
            return current_per_height * (gauss2(w.x, a, wire_sigma) - gauss2(w.x, b, wire_sigma)) * v

        f = asm(rhs, self.basis)
        # the rim is insulating (natural condition); the net source is zero, so pinning one
        # node only fixes the additive constant
        pin = np.array([int(np.argmin(np.linalg.norm(mesh.p, axis=0)))])
        self.phi = np.zeros(self.basis.N)
        free = np.setdiff1d(np.arange(self.basis.N), pin)
        Ac, fc = condense(A, f, D=pin, expand=False)
        self.phi[free] = spsolve(Ac.tocsc(), fc)
        # E per element, then nodal (area-weighted) values for probing
        grad = self.basis.interpolate(self.phi).grad  # (2, nelems, nqp)
        self.e_elem = -grad.mean(axis=2).T  # (nelems, 2)
        v1 = mesh.p[:, mesh.t[1]] - mesh.p[:, mesh.t[0]]
        v2 = mesh.p[:, mesh.t[2]] - mesh.p[:, mesh.t[0]]
        area = np.abs(v1[0] * v2[1] - v1[1] * v2[0]) / 2
        acc, w = np.zeros((self.basis.N, 2)), np.zeros(self.basis.N)
        for k in range(3):
            np.add.at(acc, mesh.t[k], self.e_elem * area[:, None])
            np.add.at(w, mesh.t[k], area)
        self.e_node = acc / w[:, None]

    def field(self, points, k=12):
        """E (n, 2) at points (n, 2) in metres, linear interpolation; NaN outside.

        scikit-fem's element_finder evaluates every point against every candidate element
        (memory ~ k n^2), so this locates points itself: the k nearest element centroids
        (KD-tree) are tested with barycentric coordinates, point by point, vectorised."""
        pts = np.asarray(points, float)
        if not hasattr(self, "_tree"):
            self._tree = cKDTree(self.mesh.p[:, self.mesh.t].mean(axis=1).T)
        out = np.full((len(pts), 2), np.nan)
        if len(pts) == 0:
            return out
        p, t = self.mesh.p, self.mesh.t
        cand = self._tree.query(pts, k=min(k, t.shape[1]))[1]  # (n, k)
        p0, p1, p2 = (p[:, t[j][cand]] for j in range(3))  # each (2, n, k)
        d = (p1[0] - p0[0]) * (p2[1] - p0[1]) - (p2[0] - p0[0]) * (p1[1] - p0[1])
        qx, qy = pts[:, 0, None] - p0[0], pts[:, 1, None] - p0[1]
        l1 = (qx * (p2[1] - p0[1]) - (p2[0] - p0[0]) * qy) / d
        l2 = ((p1[0] - p0[0]) * qy - qx * (p1[1] - p0[1])) / d
        l0 = 1 - l1 - l2
        tol = -1e-9
        ok = (l0 >= tol) & (l1 >= tol) & (l2 >= tol)
        hit = ok.any(axis=1)
        j = np.argmax(ok, axis=1)
        rows = np.flatnonzero(hit)
        tri = cand[rows, j[rows]]
        lam = np.stack([l0[rows, j[rows]], l1[rows, j[rows]], l2[rows, j[rows]]], 1)
        out[rows] = sum(lam[:, [m]] * self.e_node[t[m][tri]] for m in range(3))
        return out


@dataclass
class EnvelopeMap:
    x: np.ndarray  # grid (m)
    y: np.ndarray
    e1: np.ndarray  # (ny, nx, 2)
    e2: np.ndarray
    am_x: np.ndarray  # (ny, nx) envelope along x
    am_y: np.ndarray
    am_max: np.ndarray

    def peak(self, which="am_max", radius=R_DISC, exclude=0.0):
        """(value, (x, y)/R, (hwhm_x, hwhm_y)/R) of the field `which`, ignoring points
        within `exclude` (m) of the rim where the wires dominate."""
        f = np.array(getattr(self, which), float)
        xx, yy = np.meshgrid(self.x, self.y)
        f[~np.isfinite(f)] = -np.inf
        f[np.hypot(xx, yy) > radius - exclude] = -np.inf
        iy, ix = np.unravel_index(np.argmax(f), f.shape)
        v = f[iy, ix]
        return float(v), (float(xx[iy, ix] / radius), float(yy[iy, ix] / radius)), \
            (self._hwhm(f[iy, :], self.x, ix, v) / radius, self._hwhm(f[:, ix], self.y, iy, v) / radius)

    @staticmethod
    def _hwhm(line, coord, i0, v):
        """Half width at half maximum of a profile around index i0 (mean of both sides)."""
        half = v / 2
        widths = []
        for step in (1, -1):
            i = i0
            while 0 <= i + step < len(line) and line[i + step] >= half:
                i += step
            if 0 <= i + step < len(line) and np.isfinite(line[i + step]):  # interpolate
                t = (line[i] - half) / (line[i] - line[i + step])
                widths.append(abs(coord[i] + t * (coord[i + step] - coord[i]) - coord[i0]))
            else:
                widths.append(abs(coord[i] - coord[i0]))
        return float(np.mean(widths))


class DiscTI:
    """Two electrode pairs in the saline disc; currents in A (per pair), heights in m."""

    def __init__(self, pairs, currents=(1e-3, 1e-3), radius=R_DISC, height=H_DISC,
                 sigma=SIGMA_SALINE, nrefs=6, refine_rounds=2):
        self.pairs, self.currents, self.radius = pairs, tuple(currents), radius
        wires = [p for pair in pairs for p in pair]
        self.mesh = disc_mesh(radius, nrefs, refine_near=wires, refine_rounds=refine_rounds)
        self.fields = [DiscField(self.mesh, sigma, a, b, i / height)
                       for (a, b), i in zip(pairs, self.currents)]

    def envelope_map(self, n=201, margin=0.0):
        g = np.linspace(-self.radius, self.radius, n)
        xx, yy = np.meshgrid(g, g)
        pts = np.stack([xx.ravel(), yy.ravel()], 1)
        # the rim is a polygon: its inscribed circle is where every grid point is in a triangle
        r_in = self.radius * np.cos(np.pi / len(self.mesh.boundary_facets())) * (1 - 1e-9)
        inside = np.hypot(pts[:, 0], pts[:, 1]) <= min(r_in, self.radius - margin)
        e = [np.full((len(pts), 3), np.nan) for _ in self.fields]
        for k, f in enumerate(self.fields):
            e[k][inside, :2] = f.field(pts[inside])
            e[k][inside, 2] = 0.0
        e1, e2 = e
        am_x = directional_envelope(e1, e2, [1, 0, 0]).reshape(n, n)
        am_y = directional_envelope(e1, e2, [0, 1, 0]).reshape(n, n)
        am_max = max_envelope(np.nan_to_num(e1), np.nan_to_num(e2)).reshape(n, n)
        am_max[~inside.reshape(n, n)] = np.nan
        return EnvelopeMap(g, g, e1[:, :2].reshape(n, n, 2), e2[:, :2].reshape(n, n, 2),
                           am_x, am_y, am_max)
