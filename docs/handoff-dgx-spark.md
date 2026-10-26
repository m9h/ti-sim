# Handoff: the MIDA stage on the DGX Spark

Written 2026-10-26 on t2.macbook.pro.lan (Fedora 44, 16 cores, 31 GB shared) for the agent
that will run ti-sim's MIDA stage on the DGX Spark (aarch64, 128 GB unified memory). MIDA is
on the lab NAS, which this laptop cannot reach, and the full-resolution mesh needs more
memory than this host has. Everything else in the repo has run here and is green.

## What is done

* 14 commits on `main`; `uv run pytest -q` and the 13-stage `tests/pipeline.bats` pass with
  the optional stages skipped (41 search, 42 sensitivity run on request; 50 MIDA never ran).
* `ti_sim.mida` is tested only against a synthetic label volume (`tests/test_mida.py`,
  `tests/test_labelvol.py`): the name→tag mapping for all 115 structures, label-file
  parsing, `prepare()` (downsample, neck crop, SimNIBS `create_mesh`, internal-air relabel).
* `ti_sim.electrodes.ten_ten` builds the 10-10 cap geodesically from Nz/Iz/LPA/RPA; exact on
  a sphere, 6.4 mm mean / 10.7 mm max against ernie's Jurak cap. It is the only placement
  input `examples/mida_ti.py` needs.
* Reference numbers to compare against are in `docs/research/02-replication-targets.md`
  (ernie: Wessel putamen 0.39 V/m at 2 mA; Violante hippocampus; Wessel on MIDA 0.26 V/m
  in the paper; Rampersad 0.24).

## What the Spark agent has to do

1. **Get the repo.** It exists only on this laptop (`~/Workspace/ti-sim`, not on GitHub).
   Morgan publishes it (`gh repo create m9h/ti-sim --public --source . --push`) or rsyncs
   it; `data/` and `results/` are gitignored and not needed (ernie is fetched by
   `scripts/fetch-ernie.sh`, sha256 in the script).
2. **Environment.** Either the container (`scripts/container-run.sh build`; base fedora:44,
   the COPR has aarch64 builds of `python3-simnibs 4.6.0-21` and `simnibs-data`, the repo
   file uses `$basearch`) or natively: `dnf install python3-simnibs simnibs-data bats` from
   the `mhough/neurofedora` COPR, `pip install --user uv`, then `uv venv
   --system-site-packages --python /usr/bin/python3 && uv sync`. Check:
   `uv run pytest -q -m "not slow"` (~10 s). Nothing uses the GPU; SimNIBS's FEM is CPU
   (pardiso is x86-only, so SimNIBS will use its PETSc/scipy path on aarch64: expect
   slower solves than the ernie timings in `docs/running.md`, which are x86 pardiso).
3. **MIDA.** Morgan's download is on the lab TrueNAS export `192.168.108.117:/mnt/tank/shared`
   (legion mounts it at `/data`): `datasets/MIDAv1-0/MIDA_v1.0/MIDA_v1_voxels/` holds
   `MIDA_v1.nii`, `.txt`, `.mat`, `.raw`, with `MIDA_v1_surfaces/` beside it; the `__MACOSX/`
   sibling is AppleDouble stubs, ignore it (per the 3dslicer session, 2026-10-26). Read from
   the NAS export, not through legion, and copy the voxel files to a local directory,
   e.g. `data/mida/` (gitignored). Needed: `MIDA_v1.nii` (or `.nii.gz`) and optionally `MIDA_v1.txt`, found
   recursively. Run `scripts/check-mida.sh data/mida`, then `export TI_SIM_MIDA=$PWD/data/mida`.
   Licence (`docs/research/04-mida.md`): never commit, package, containerise or upload the
   volume, any crop or mesh of it, or an undisguised face image. Caches go to
   `$TI_SIM_WORK` (default `data/ti-sim-work`), which stays on that machine.
4. **Fiducials.** `TI_SIM_MIDA_FIDUCIALS=fiducials.json` with
   `{"Nz": [x,y,z], "Iz": [x,y,z], "LPA": [x,y,z], "RPA": [x,y,z]}` in the NIfTI's world mm
   (`nib.load(...).affine` applied to voxel indices). The 3dslicer session may have a Slicer
   markups file; otherwise pick them on the skin: Nz nasion (frontonasal suture), Iz inion,
   LPA/RPA preauricular points. `ti_sim.electrodes.ten_ten` snaps them to the skin surface.
   `TI_SIM_MIDA_ELECTRODES` (JSON `{name: [x,y,z]}`) overrides individual positions.
5. **Runs**, cheapest first, each under a memory cap (`systemd-run --user --scope -p
   MemoryMax=…`, or `--memory` on the container):
   * `TI_SIM_MIDA_DS=4 TI_SIM_MIDA_ZMIN=<mm> TI_SIM_OUT=results/mida_ds4.json uv run
     examples/mida_ti.py` — 2 mm voxels, neck cropped: a smoke test of the mapping and
     mesh, minutes.
   * `TI_SIM_MIDA_DS=2` (1 mm, the default): the comparison run. On this laptop the
     ernie-sized meshes took ~2.5 min per montage; a 1 mm MIDA mesh has roughly 4-5× the
     elements. Budget 20-30 GB.
   * `TI_SIM_MIDA_DS=1` (500 µm native): only on the Spark. SimNIBS `create_mesh` at this
     resolution is the step that exceeded 10 GB here on synthetic data at 1 mm already with
     the laptop's other load; expect 40-80 GB and an hour or more. Set
     `OMP_NUM_THREADS` to the core count.
   * Then the whole pipeline: `TI_SIM_MIDA=… TI_SIM_MIDA_FIDUCIALS=… scripts/run-pipeline.sh`
     (stage 50 runs, 41/42 with `TI_SIM_LEADFIELD=1`, `TI_SIM_SENSITIVITY=1`).
6. **Report back** (into `docs/research/02-replication-targets.md`, MIDA section, and the
   JSON in `results/`): tissue volumes per SimNIBS tag against Iacono's Table 1 (the
   mapping's first sanity check), node/element counts and wall time per DS, the Violante and
   Wessel montage ROI statistics (hippocampus, putamen; median, p95, volume above 0.2 V/m
   in and outside the ROI) against the ernie numbers and the paper's 0.26 V/m, and the
   constructed electrode positions for a visual check in Slicer (the JSON holds F3, F4,
   TP7, TP8 and Cz under `ten_ten_from_fiducials`; dump the full dict from `ten_ten` if
   more are wanted) with the face disguised before any figure leaves the machine.

## Status on the Spark (2026-10-26)

Done by the neurospark session: repo cloned from GitHub; Docker image `ti-sim` built from the
Containerfile on aarch64 (plus `mmg3d meshfix cat-surface time`, which SimNIBS's external-bin
symlinks need when weak deps are off; simnibs 4.6.0-22 fixes the first two); MIDA copied from
the NAS export to the T9 (`/mnt/t9/datasets/mida`, never in git); label-file parser fixed for
the download layout (commit 4fda967); fiducials estimated from the label volume
(`fiducials.json` beside the MIDA copy); 2 mm and 1 mm runs with Wessel's montage completed
(`mida_ds4.json`, `mida_ds2.json` under `/mnt/t9/ti-sim-host/outputs`); numbers in
`docs/research/02-replication-targets.md` §4. The host dir is `/mnt/t9/ti-sim-host`
(`inputs/` → `/data`, `outputs/` → `/results`); runs use `TI_SIM_WORK=/data/ti-sim-work`
because `examples/mida_ti.py` defaults to a relative path inside the read-only package.
Open: the 0.5 mm run; Violante's e1-e4; a Slicer check of the constructed electrodes; the
phantom test `test_envelope_map_symmetry_and_peak` picks the mirror peak (y = +0.59) on
aarch64 for the symmetric layout, so its sign assertion should accept either.

## Known soft spots to watch

* `mida.RULES` maps cranial nerves to WM and dura to CSF (charm's convention); if a volume
  comparison looks off, that table in `src/ti_sim/mida.py` is the place to look.
* `relabel_internal_air` turns enclosed air (sinuses) into a SimNIBS air tag; if
  `create_mesh` complains about unconnected regions, raise `TI_SIM_MIDA_ZMIN` to cut the
  neck or check for label islands after downsampling (nearest-neighbour, no smoothing).
* The example's `montages` dict holds the two montages: Wessel (F3-F4 + TP7-TP8, 2 mA,
  9.8 mm radius) runs from the 10-10 construction alone; Violante's hippocampal montage
  expects four positions named `e1`..`e4` in `TI_SIM_MIDA_ELECTRODES` (their paper gives
  individualised positions, not 10-10 names) and is skipped with a message otherwise.
* The JSON written to `TI_SIM_OUT` has `volumes_cm3` per SimNIBS tag (21-30 are the extra
  deep structures, see `mida.EXTRA_COND`), `n_tets`, and per montage the `am` statistics
  for GM, brain and each extra structure.
* Nothing here is multi-pair exact: the envelope is Grossman's two-pair formula
  (`ti_sim.envelope`, equal to SimNIBS `TI_utils.get_maxTI`).

## Also for the Spark (separate repo, `~/Workspace/fem-trials`, also local only)

The complex-valued dolfinx + python-mumps solver scaled 68.9k → 275k dofs at 1.2 → 7.2 GB
on this laptop; the 1e6-dof Helmholtz-with-PML case for ae-sim is estimated at ~40 GB and
was deferred to the Spark. `dolfinx-complex/README.md` has the container recipe and the
timings; the Fedora `python3-dolfinx-openmpi` and `python-mumps` packages used are in the
same COPR for aarch64 too.
