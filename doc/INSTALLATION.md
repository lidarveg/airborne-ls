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
