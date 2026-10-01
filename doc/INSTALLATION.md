# Installing airborne-ls

## Dependencies
The following packages are required by this package, and should be installed using your favourite
installation method. See below for some common approaches.

  * [gdal](https://gdal.org) (with Python/numpy bindings)
  * [scipy](https://scipy.org/)
  * [laspy](https://laspy.readthedocs.io/en/latest/)
  * [lazrs](https://github.com/laz-rs/laz-rs) with [Python bindings](https://github.com/laz-rs/laz-rs-python)
  * [pyproj](https://github.com/pyproj4/pyproj)
  * [numba](https://numba.pydata.org/)
  * [pynninterp](https://gitlab.com/jrsrp/sys/lidar/pynninterp) (installing this will require a C compiler)

The `airborne-ls` package itself is only hosted on its Github repository, and can be installed directly from there.

### Optional Dependencies
  * [richdem](https://richdem.readthedocs.io/). Used for hydrological flow accumulation from
    the computed DEM. Can be installed with conda or pip.

## Installation with conda
Most of the required packages can be installed using [conda](https://continuumio-docs.readthedocs-hosted.com/miniconda/), from the `conda-forge` channel. It is strongly recommended to install in a conda environment specific for this work. For example

```bash
conda create --channel conda-forge -n airborne-ls pip gdal scipy laspy lazrs-python pyproj numba
conda activate airborne-ls
pip install git+https://github.com/lidarveg/airborne-ls
```

The pip command installs both `pynninterp` and `airborne-ls`, but uses the other packages from those
installed by conda.

## Installation with pip
If installing everything via `pip`, you should work in some sort of virtual environment
(e.g. `venv` or `virtualenv`) to avoid clashes with other installations.

You will need to have the GDAL libraries already installed on your system by some other means,
as this cannot be installed via `pip`. This should include the Python bindings, with `numpy` support.

Once you have a full installation of GDAL, the remaining requirements can be installed using `pip`, as follows

```bash
pip install git+https://github.com/lidarveg/airborne-ls
```

## Installation with `uv`

Since this project depends on the gdal python bindings, which depend on your
system's gdal version, we need to be careful when adding the gdal package. The best
way of doing this is by adding a constraint dependency to your project.

Installation in a uv project might proceed like:

```sh
uv init airbproj
cd airbproj
# add a constraint matching your system's gdal version
cat >> pyproject.toml <<EOF
[tool.uv]
constraint-dependencies = ["gdal==$(gdal-config --version)"]
EOF

uv add git+https://github.com/lidarveg/airborne-ls
```

If your gdal version is < 3.9, then you will also need to restrict numpy to be < 2.
The constraint dependency then becomes:

```sh
cat >> pyproject.toml <<EOF
[tool.uv]
constraint-dependencies = ["gdal==$(gdal-config --version)", "numpy<2"]
EOF
```

The gdal python bindings on pypi are source distribution only, and `pynninterp` will
also need to be built from source, so you will need compilers and git to install the package.
On ubuntu like systems, this would mean `apt-get install build-essential python3-dev git`.
The python headers `python3-dev` are needed if you use the system python distribution. Python
distributed by uv includes headers, so if you prefer to use uv managed python the package
`python3-dev` is not needed. You would need to set up your uv project using something like;

```sh
uv init --managed-python --python 3.12 airbproj
```

## Docker

See the [`Dockerfile`](../examples/Dockerfile) for an example of how the installation can be done in a container.

```sh
docker build --tag airbornels:latest .
```

Then the commands can be run like:

```sh
docker run --rm \
  -u "$(id -u):$(id -g)" \
  -v "$PWD":/work \
  airbornels:latest \
  alsv_lidar_standardisation --indir alsProjectDir --outdir alsProjectDir-std \
    --intilesize 1000 --epsg 28356 --ii mp --project myproj \
    --year 2025 --binsize 50

```
