#!/bin/bash
# Build / run the ti-sim container (CPU only: SimNIBS's FEM is CPU; nothing here needs a GPU).
#
#   scripts/container-run.sh build                 # image ti-sim:latest (TI_SIM_IMAGE)
#   scripts/container-run.sh run [bats args]       # whole pipeline; ./data and ./results mounted
#   scripts/container-run.sh shell                 # interactive shell in the image
#   scripts/container-run.sh smoke                 # unit tests + the phantom example
#
# TI_SIM_HOST_DIR: host directory holding inputs/ (ernie, MIDA) and outputs/.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
IMAGE="${TI_SIM_IMAGE:-ti-sim:latest}"
RT="${TI_SIM_RUNTIME:-}"
[ -n "$RT" ] || { command -v podman >/dev/null && RT=podman; } || { command -v docker >/dev/null && RT=docker; } || true
[ -n "$RT" ] || { echo "neither podman nor docker found" >&2; exit 1; }
HOST="${TI_SIM_HOST_DIR:-$ROOT}"
mkdir -p "$HOST/inputs" "$HOST/outputs"
MOUNTS=(-v "$HOST/inputs:/data:Z" -v "$HOST/outputs:/results:Z")
COMMON=(--rm --memory="${TI_SIM_MEMORY:-10g}" --shm-size="${TI_SIM_SHM:-4g}")
if [ "$RT" = podman ]; then COMMON+=(--userns=keep-id); else COMMON+=(--user "$(id -u):$(id -g)"); fi
ENVS=(-e "TI_SIM_FETCH=${TI_SIM_FETCH:-1}" -e "OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}"
      ${TI_SIM_LEADFIELD:+-e TI_SIM_LEADFIELD=$TI_SIM_LEADFIELD} ${TI_SIM_MIDA:+-e TI_SIM_MIDA=/data/mida})
case "${1:-run}" in
    build) $RT build -t "$IMAGE" ${TI_SIM_BASE:+--build-arg "BASE=$TI_SIM_BASE"} -f "$ROOT/Containerfile" "$ROOT" ;;
    run)   shift; $RT run "${COMMON[@]}" "${MOUNTS[@]}" "${ENVS[@]}" "$IMAGE" "$@" ;;
    shell) $RT run -it "${COMMON[@]}" "${MOUNTS[@]}" "${ENVS[@]}" --entrypoint /bin/bash "$IMAGE" ;;
    smoke) $RT run "${COMMON[@]}" "${MOUNTS[@]}" "${ENVS[@]}" --entrypoint /bin/bash "$IMAGE" -c \
             'cd /opt/ti-sim && uv run --no-sync pytest -q -m "not slow" --deselect tests/test_head.py && TI_SIM_OUT=/results/grossman_phantom.json uv run --no-sync python examples/grossman_phantom.py' ;;
    *) echo "usage: $0 build|run|shell|smoke" >&2; exit 2 ;;
esac
