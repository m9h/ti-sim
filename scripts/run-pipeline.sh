#!/bin/bash
# Run the whole ti-sim pipeline under BATS with a clean environment (no conda, no stray
# PYTHON* variables). Pass-through: TI_SIM_* and OMP_NUM_THREADS. Extra args go to bats.
#
#   scripts/run-pipeline.sh                      # phantom, sphere, ernie (if present)
#   TI_SIM_FETCH=1 scripts/run-pipeline.sh       # allow the ernie download (413 MB)
#   TI_SIM_LEADFIELD=1 scripts/run-pipeline.sh   # also the montage search (20-40 min)
#   TI_SIM_MIDA=/path scripts/run-pipeline.sh    # also the MIDA stage
set -u
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
UV_BIN="$(command -v uv || echo "$HOME/.local/bin/uv")"
env -i PATH="$(dirname "$UV_BIN"):/usr/bin:/bin" HOME="$HOME" TERM="${TERM:-dumb}" \
    $(env | grep -E '^(TI_SIM_[A-Z0-9_]+|OMP_NUM_THREADS)=') \
    bats --timing --print-output-on-failure "$@" "$ROOT/tests/pipeline.bats" < /dev/null
