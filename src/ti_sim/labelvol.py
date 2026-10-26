"""Label volumes to SimNIBS meshes: synthetic layered spheres and any labelled image.

SimNIBS 4.6's `create_mesh` (CGAL) turns a uint8/uint16 label image with a voxel-to-world
affine into a tetrahedral head mesh with surface tags, which SimNIBS then simulates like a
charm head (electrodes placed on the 1005 skin surface, conductivities by tag). That is the
route for the MIDA model (`ti_sim.mida`) and for the layered spheres the early TI papers
used (Grossman 2017 figure S2: scalp 0.05 R, skull 0.085 R, CSF 0.023 R, brain 0.83 R).

Label 0 is air; internal air cavities become holes whose boundary SimNIBS would tag as
skin, so `mesh_from_labels` runs `relabel_internal_air` (as charm does) before writing.
Tags follow SimNIBS ElementTags; extra tags (deep nuclei) go after GM in the surface
hierarchy. Meshes are cached in `workdir` under a key of the inputs.
"""

import hashlib
import json
import os

import numpy as np

# CGAL sizing: element size = slope * local tissue thickness, clipped to range (mm)
MESH_SIZES = {"standard": {"range": [1.5, 6.0], "slope": 1.0}}
FACET_DISTANCES = {"standard": {"range": [0.2, 3.0], "slope": 0.5}}
DEFAULT_HIERARCHY = (1, 2, 9, 3, 4, 8, 7, 6, 11, 10, 12, 5)


def sphere_labels(radii_mm, tags, spacing_mm=1.0, pad_mm=4.0, blobs=()):
    """Nested spheres: radii_mm outer to inner, one tag each (e.g. (85, 80.75, 73.5, 71.5)
    with (5, 4, 3, 2)); `blobs` = [(centre_mm, radius_mm, tag)] painted last (deep ROIs).
    Returns (labels uint16 (nx, ny, nz), affine mm) centred on the origin."""
    radii, tags = [float(r) for r in radii_mm], [int(t) for t in tags]
    assert len(radii) == len(tags) and radii == sorted(radii, reverse=True)
    half = radii[0] + pad_mm
    n = int(np.ceil(2 * half / spacing_mm)) + 1
    ax = -half + spacing_mm * np.arange(n)
    x, y, z = np.meshgrid(ax, ax, ax, indexing="ij")
    r = np.sqrt(x * x + y * y + z * z)
    labels = np.zeros(r.shape, np.uint16)
    for rad, tag in zip(radii, tags):
        labels[r <= rad] = tag
    for centre, rad, tag in blobs:
        c = np.asarray(centre, float)
        labels[np.sqrt((x - c[0]) ** 2 + (y - c[1]) ** 2 + (z - c[2]) ** 2) <= rad] = tag
    affine = np.diag([spacing_mm, spacing_mm, spacing_mm, 1.0])
    affine[:3, 3] = -half
    return labels, affine


def mesh_from_labels(labels, affine, workdir, name, elem_sizes=None, facet_distances=None,
                     extra_tags=(), threads=4, params=None):
    """SimNIBS create_mesh on a label volume -> path of the cached .msh in workdir."""
    from simnibs import mesh_io
    from simnibs.mesh_tools import meshing

    os.makedirs(str(workdir), exist_ok=True)
    key = hashlib.sha1(json.dumps([labels.shape, affine.tolist(), elem_sizes, facet_distances,
                                   list(extra_tags), params, hashlib.sha1(labels.tobytes()).hexdigest()],
                                  sort_keys=True, default=str).encode()).hexdigest()[:10]
    fn = os.path.join(str(workdir), f"{name}_{key}.msh")
    if os.path.exists(fn):
        return fn
    extra = [int(t) for t in extra_tags if int(t) not in DEFAULT_HIERARCHY]
    hierarchy = (1,) + tuple(extra) + DEFAULT_HIERARCHY[1:]  # deep nuclei right after WM/GM
    msh = meshing.create_mesh(labels.astype(np.uint16), affine, elem_sizes=elem_sizes or MESH_SIZES,
                              facet_distances=facet_distances or FACET_DISTANCES,
                              hierarchy=hierarchy, num_threads=threads)
    msh = msh.relabel_internal_air()
    mesh_io.write_msh(msh, fn)
    return fn


def surface_point(labels, affine, direction):
    """World coordinate (mm) of the outermost labelled voxel along `direction` from the
    volume centre: where to put an electrode on a synthetic head."""
    d = np.asarray(direction, float)
    d /= np.linalg.norm(d)
    centre = (np.array(labels.shape) - 1) / 2
    vox = np.array(labels.shape).max()
    for s in np.linspace(vox / 2, 0, int(vox / 2) + 1):
        ijk = np.round(centre + s * d / np.abs(affine[:3, :3]).max() * np.abs(affine[:3, :3]).max()).astype(int)
        if np.all(ijk >= 0) and np.all(ijk < labels.shape) and labels[tuple(ijk)] > 0:
            return (affine @ np.append(ijk, 1))[:3]
    raise ValueError("no labelled voxel along that direction")
