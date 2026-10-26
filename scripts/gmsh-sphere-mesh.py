"""The 4-layer sphere of examples/layered_sphere.py as a CAD mesh (gmsh OCC), written as a
SimNIBS head mesh: tets tagged 2-5, boundary triangles 1002-1005. Geometry and size field
follow ci/conductivity-admm/sphere_ti_forward.py in gitlab.com/morgan.hough/neurofedora
(the DOLFINx reference this was cross-checked against on 2026-10-26), so the two solvers
see the identical mesh: four cylinders of radius 10 mm cut the scalp where the electrodes
go, hmin at every interface growing to hmax.

    uv run python scripts/gmsh-sphere-mesh.py data/work/gmsh_sphere85.msh 1.5 6.0

Needs python3-gmsh (Fedora) and SimNIBS. Use the file as `fnamehead` of TwoPairSimulation.
"""
import sys, numpy as np, gmsh
out = sys.argv[1]; hmin, hmax = float(sys.argv[2]), float(sys.argv[3])
R = 85.0
LAYERS = [(5, 85.000), (4, 80.750), (3, 73.525), (2, 71.570)]
ANGLES = {"1a": 180 + 11.5, "1b": 180 - 11.5, "2a": -11.5, "2b": 11.5}
gmsh.initialize(); gmsh.option.setNumber("General.Terminal", 0)
occ = gmsh.model.occ
spheres = [occ.addSphere(0, 0, 0, r) for _, r in LAYERS]
cyls = {}
for name, deg in ANGLES.items():
    n = np.array([np.cos(np.radians(deg)), np.sin(np.radians(deg)), 0.0])
    cyls[name] = occ.addCylinder(*((R - 4.0) * n), *(8.0 * n), 10.0)
res, _ = occ.fragment([(3, s) for s in spheres], [(3, c) for c in cyls.values()])
occ.synchronize()
keep = []
for d, v in res:
    if d != 3: continue
    if np.linalg.norm(occ.getCenterOfMass(3, v)) > R: occ.remove([(3, v)], recursive=True)
    else: keep.append(v)
occ.synchronize()
by_tag = {t: [] for t, _ in LAYERS}
for v in keep:
    bb = occ.getBoundingBox(3, v); extent = max(abs(b) for b in bb)
    if np.linalg.norm(occ.getCenterOfMass(3, v)) > 50: by_tag[5].append(v); continue
    by_tag[min(LAYERS, key=lambda L: abs(L[1] - extent))[0]].append(v)
for t, vs in by_tag.items(): gmsh.model.addPhysicalGroup(3, vs, t)
# surfaces: boundary of each tissue group, inner-first so interfaces get the inner tag
assigned = set()
for t in (2, 3, 4, 5):
    b = gmsh.model.getBoundary([(3, v) for v in by_tag[t]], combined=True, oriented=False)
    ss = [s for d, s in b if d == 2 and s not in assigned]
    gmsh.model.addPhysicalGroup(2, ss, 1000 + t); assigned.update(ss)
interfaces = [s for d, s in gmsh.model.getEntities(2)]
f = gmsh.model.mesh.field
dist = f.add("Distance"); f.setNumbers(dist, "SurfacesList", interfaces)
thr = f.add("Threshold")
for k, v in (("InField", dist), ("SizeMin", hmin), ("SizeMax", hmax), ("DistMin", 0.0), ("DistMax", hmax - hmin)):
    f.setNumber(thr, k, v)
f.setAsBackgroundMesh(thr)
for opt in ("MeshSizeExtendFromBoundary", "MeshSizeFromPoints", "MeshSizeFromCurvature"):
    gmsh.option.setNumber(f"Mesh.{opt}", 0)
gmsh.option.setNumber("General.NumThreads", 1)
gmsh.model.mesh.generate(3)
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
gmsh.option.setNumber("Mesh.Binary", 0)
gmsh.write(out); gmsh.finalize()
from simnibs import mesh_io
m = mesh_io.read_msh(out)
m = m.crop_mesh(elm_type=[2, 4])
m.elm.tag2[:] = m.elm.tag1
tets = m.elm.elm_type == 4
print("tets", tets.sum(), "tris", (~tets).sum(), "tet tags", sorted(set(m.elm.tag1[tets].tolist())), "tri tags", sorted(set(m.elm.tag1[~tets].tolist())))
vol = m.elements_volumes_and_areas().value
for t in (2, 3, 4, 5): print(f" tag {t}: {vol[tets & (m.elm.tag1 == t)].sum()/1e3:.1f} cm3")
m.write(out)
