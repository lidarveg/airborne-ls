# airborne-ls

Currently under development

Airborne-ls is designed to be run on individual laz tiles with final products mosaiced together

Scripts to be run and in order:

1.) run_lidar_stardardisation.py
2.) run_tile_products.py
3.) run_dem_correction.py
4.) run_product_mosaic.py

And is best managed by a workflow. An example will be provided.

## Installation

Depends on packages [rios](https://github.com/ubarsc/rios/tree/rios-2.0.9) and
[pynninterp](https://gitlab.com/jrsrp/sys/lidar/pynninterp) which are not available on pypi. Installing `gdal` also requires
having the system libraries available, and your python version should match.

Users will need to manually install these along with this packge. eg

```sh
uv add gdal[numpy]==$(gdal-config --version) \
    "git+https://github.com/ubarsc/rios@rios-2.0.9" \
    "git+https://gitlab.com/jrsrp/sys/lidar/pynninterp@1.0.1" \
    "git+https://github.com/ubarsc/pynninterp"

uv add airborne-ls
```


For development, it might be more convenient to install these using `uv pip install`, eg

```sh
uv add gdal[numpy]==$(gdal-config --version) # make sure your version matches your system
uv pip install "git+https://github.com/ubarsc/rios@rios-2.0.9"
uv pip install git+https://gitlab.com/jrsrp/sys/lidar/pynninterp@1.0.1
```
