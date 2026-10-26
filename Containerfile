# ti-sim container: SimNIBS 4.6 from the mhough/neurofedora COPR, scikit-fem and the
# project in a uv venv, BATS pipeline as the entrypoint. CPU only (SimNIBS's FEM is CPU).
#
#   scripts/container-run.sh build          # image ti-sim:latest
#   scripts/container-run.sh smoke          # unit tests + the Grossman phantom
#   scripts/container-run.sh run            # the whole pipeline (tests/pipeline.bats)
#
# Base: fedora:44 (python3-simnibs 4.6.0-21 has x86_64 and aarch64 builds in the COPR).
# The container runs as the invoking user (read-only venv: TI_SIM_NO_SYNC); caches in /tmp.
ARG BASE=registry.fedoraproject.org/fedora:44
FROM ${BASE}
COPY container/mhough-neurofedora.repo /etc/yum.repos.d/
RUN dnf -y --setopt=install_weak_deps=False install \
        python3-simnibs simnibs-data python3 bats git-core curl unzip \
    && dnf clean all \
    && python3 -c 'import simnibs; print("simnibs", simnibs.__version__)'
COPY --from=ghcr.io/astral-sh/uv:latest /uv /usr/local/bin/uv

WORKDIR /opt/ti-sim
COPY pyproject.toml uv.lock README.md ./
COPY src ./src
COPY tests ./tests
COPY examples ./examples
COPY scripts ./scripts
RUN uv venv --system-site-packages --python /usr/bin/python3 \
    && uv sync \
    && uv run --no-sync python -c "import simnibs, skfem; print('solvers ok')" \
    && chmod -R a+rX /opt/ti-sim

ENV TI_SIM_DATA=/data TI_SIM_RESULTS=/results TI_SIM_FETCH=1 TI_SIM_NO_SYNC=1 \
    OMP_NUM_THREADS=4 MPLBACKEND=Agg HOME=/tmp/ti-home MPLCONFIGDIR=/tmp/ti-home/mpl \
    UV_CACHE_DIR=/tmp/ti-home/uv NUMBA_CACHE_DIR=/tmp/ti-home/numba
RUN mkdir -p /tmp/ti-home && chmod 1777 /tmp/ti-home
VOLUME ["/data", "/results"]
ENTRYPOINT ["/opt/ti-sim/scripts/run-pipeline.sh"]
