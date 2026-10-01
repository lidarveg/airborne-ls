# syntax=docker/dockerfile:1

# ---- Versions (override with --build-arg) ----------------------------------
ARG GDAL_VERSION=3.12.4
# Pin this to a specific uv release (e.g. 0.9.x) for reproducible builds
ARG UV_VERSION=latest

FROM ghcr.io/astral-sh/uv:${UV_VERSION} AS uv
FROM ghcr.io/osgeo/gdal:ubuntu-full-${GDAL_VERSION}

ARG PYTHON_VERSION=3.12
ARG AIRBORNE_LS_REF=main

# ---- System build tools ----------------------------------------------------
# Drop the stale Apache Arrow apt source (if present) so apt-get update works.
# No python3-dev: uv-managed Python ships its own headers.
RUN  apt-get update \
 && apt-get install -y --no-install-recommends build-essential git ca-certificates \
 && rm -rf /var/lib/apt/lists/*

# ---- uv --------------------------------------------------------------------
COPY --from=uv /uv /uvx /usr/local/bin/

ENV UV_PYTHON_PREFERENCE=only-managed \
    UV_PYTHON_INSTALL_DIR=/opt/uv-python \
    UV_LINK_MODE=copy

# ---- Project environment ---------------------------------------------------
WORKDIR /opt/alsenv

RUN --mount=type=cache,target=/root/.cache/uv \
    uv python install "${PYTHON_VERSION}" \
 && uv init --bare --name alsenv --python "${PYTHON_VERSION}" \
 && printf '\n[tool.uv]\npython-preference = "only-managed"\nconstraint-dependencies = ["gdal==%s"]\n' \
        "$(gdal-config --version)" >> pyproject.toml \
 && cat pyproject.toml \
 && uv add "git+https://github.com/lidarveg/airborne-ls@${AIRBORNE_LS_REF}" \
 && uv run python -c "from osgeo import gdal, gdal_array; print('GDAL bindings', gdal.__version__)" \
 && uv run alsv_tile_products --help > /dev/null

# Make the venv's commands available without activation
ENV VIRTUAL_ENV=/opt/alsenv/.venv \
    PATH="/opt/alsenv/.venv/bin:${PATH}"

WORKDIR /work
CMD ["bash"]
