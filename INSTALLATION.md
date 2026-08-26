# Installing airborne-ls

## Dependencies
The following packages are required by this package, and should be installed using your favourite installation method.

  * [gdal](https://gdal.org) (with Python/numpy bindings)
  * [scipy](https://scipy.org/)
  * [laspy](https://laspy.readthedocs.io/en/latest/)
  * [lazrs](https://github.com/laz-rs/laz-rs) with [Python bindings](https://github.com/laz-rs/laz-rs-python)
  * [numba](https://numba.pydata.org/)
  * [pynninterp](https://gitlab.com/jrsrp/sys/lidar/pynninterp) (installing this will require a C compiler)

The `airborne-ls` package itself is only hosted on its Github repository, and can be installed directly from there.

## Installation with conda
Most of the required packages can be installed using [conda](https://continuumio-docs.readthedocs-hosted.com/miniconda/), from the `conda-forge` channel. It is strongly recommended to install in a conda environment specific for this work. For example

```bash
conda create --channel conda-forge -n airborne-ls pip gdal scipy laspy lazrs-python numba
conda activate airborne-ls
pip install git+https://gitlab.com/jrsrp/sys/lidar/pynninterp@1.0.1
pip install --no-deps git+https://github.com/lidarveg/airborne-ls
```

## Installation with pip
If installing everything via `pip`, you should work in some sort of virtual environment (e.g. `venv` or `virtualenv`) to avoid clashes with other installations.

You will need to have the GDAL libraries already installed on your system by some other means, as this cannot be installed via `pip`. This should include the Python bindings, with `numpy` support.

Once you have a full installation of GDAL, the remaining requirements can be installed using `pip`, as follows

```bash
pip install scipy
pip install laspy[lazrs]
pip install numba
pip install git+https://gitlab.com/jrsrp/sys/lidar/pynninterp@1.0.1
pip install --no-deps git+https://github.com/lidarveg/airborne-ls
```

The `pyproject.toml` specifies exact versions for all dependencies, because it was constructed
for use with `uv`. For this reason, it should not be used by itself to manage dependencies.

## Installation with uv
The provided `uv.lock` file will not work for anyone outside of JRSRP, as it requires credentials to access their registry. Good luck with that.

The `uv` installation instructions provided by the JRSRP people are reproduced below.

---

```sh
uv add gdal[numpy]==$(gdal-config --version) \
    "git+https://gitlab.com/jrsrp/sys/lidar/pynninterp@1.0.1"

uv add airborne-ls
```


For development, it might be more convenient to install these using `uv pip install`, eg

```sh
uv add gdal[numpy]==$(gdal-config --version) # make sure your version matches your system
uv pip install git+https://gitlab.com/jrsrp/sys/lidar/pynninterp@1.0.1
```
