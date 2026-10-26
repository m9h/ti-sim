# Running ti-sim

A `uv` project whose venv is created with `--system-site-packages` so it sees the system
SimNIBS (`python3-simnibs`, `simnibs-data` from the mhough/neurofedora COPR). scikit-fem,
numpy (< 2.5 for SimNIBS's numba), scipy, matplotlib and nibabel come from uv.

```sh
uv venv --system-site-packages --python /usr/bin/python3 && uv sync    # once
uv run pytest -q -m "not slow"                                          # ~10 s
uv run pytest -q                                                        # + sphere mesh route, ~1 min
TI_SIM_M2M=data/m2m_ernie uv run pytest -q tests/test_head.py           # + ernie, ~1 min
```

Heavy runs on a shared machine go under a capped scope:
`systemd-run --user --scope -p MemoryMax=10G nice -n 19 scripts/run-pipeline.sh`.
The whole pipeline needs about 2 GB; the full-resolution MIDA mesh will need more.

## Examples and their outputs

Each example writes JSON to `$TI_SIM_OUT` (`ti_sim.report.save_results`) and figures or
meshes beside it, and prints a human summary.

| example | needs | time | writes |
|---|---|---|---|
| `grossman_phantom.py` | nothing | 10 s | peaks/positions/HWHM for layouts A-E vs the paper, `grossman_phantom.png` |
| `layered_sphere.py` | SimNIBS | ~2 min (mesh cached) | focus depth and steering on the 4-layer sphere |
| `ernie_ti.py` | `TI_SIM_M2M` | ~2.5 min (cached) | ROI statistics for two montages, `ernie_*.msh` for gmsh |
| `ernie_montage_search.py` | `TI_SIM_M2M`, 20-40 min once for the lead field | | best two-pair montages for the left hippocampus |
| `mida_ti.py` | `TI_SIM_MIDA` (+ `TI_SIM_MIDA_ELECTRODES`) | mesh: minutes | tissue volumes; Violante/Wessel montages |

## Environment

| variable | meaning | default |
|---|---|---|
| `TI_SIM_DATA` | data directory (ernie, MIDA, caches) | `./data` |
| `TI_SIM_M2M` | SimNIBS `m2m_ernie` folder | `$TI_SIM_DATA/m2m_ernie` |
| `TI_SIM_WORK` | SimNIBS runs and meshes (cached by parameter hash) | `$TI_SIM_DATA/ti-sim-work` |
| `TI_SIM_RESULTS` | JSON/PNG/msh output (pipeline) | `./results` |
| `TI_SIM_OUT` | JSON file for one example | unset: no file |
| `TI_SIM_FETCH=1` | allow `scripts/fetch-ernie.sh` in the pipeline | off |
| `TI_SIM_COND` | conductivity table: `simnibs`, `rampersad`, or `simnibs,4=0.0085` overrides | `simnibs` |
| `TI_SIM_LEADFIELD=1` | run the lead-field montage search stage | off |
| `TI_SIM_MIDA`, `TI_SIM_MIDA_DS`, `TI_SIM_MIDA_ZMIN`, `TI_SIM_MIDA_ELECTRODES` | MIDA directory, downsample factor (2), neck cut (mm), electrode JSON | unset |
| `OMP_NUM_THREADS` | SimNIBS/FEM threads | 4 |

## Pipeline (`tests/pipeline.bats`)

Stages: 00 toolchain, 01 `uv sync`, 02 SimNIBS visible, 10 unit tests, 20 ernie present
(fetch with `TI_SIM_FETCH=1`), 30 Grossman phantom (geometry asserted, amplitude ratio
recorded), 31 layered sphere (focus between the pairs, deeper than |E₁|, steerable),
40 ernie montages (striatum 0.08-0.6 V/m at Wessel's 2 mA, TI never covers more brain
than tACS), 41 montage search (optional), 50 MIDA (optional), 90 summary.

## Container

```sh
scripts/container-run.sh build      # fedora:44 + COPR SimNIBS + uv venv -> ti-sim:latest
scripts/container-run.sh smoke      # unit tests + the phantom
scripts/container-run.sh run        # the pipeline; ./inputs -> /data, ./outputs -> /results
```

Put `m2m_ernie` (or let `TI_SIM_FETCH=1` fetch it) and a MIDA copy under `inputs/`.
MIDA must never be baked into an image (licence).

## gmsh

`TwoPairSimulation.write_msh` writes |E₁|, |E₂|, `TImax` and `tACS` element fields;
`gmsh results/ernie_wessel2023_striatum.msh` shows them (gmsh is in Fedora).
