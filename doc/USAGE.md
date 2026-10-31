# Usage
The package is designed to process a contiguous capture of airborne LiDAR data, generally
known as a single "project". Such data is typically supplied as a collection of LAS/LAZ
files, each file being the data for a single square tile, and the collection of tiles
making up the capture for the whole project. Such data will typically have been acquired
in multiple flight lines over the region, and then divided up into tiles by the 
data providers.

It is also expected that point return classes will have been assigned by the data providers,
in accordance with the LAS standards. At a minimum, this should include classification
of ground returns, distinct from non-ground returns.

These are processed using the commands described below. All of these commands will
accept a `-h` or `--help` option to print brief help on all options. Many command line
options will have sensible default values, but these should be checked for suitability.

## Re-tiling LAS/LAZ files
The supplied LAS/LAZ files may occasionally be provided in the form of un-tiled flightlines.
In this case, it is necessary to re-tile as an independent first step, before attempting to
ingest and standardize the tiles. This step is not generally necessary if the data is supplied
in tiles already.

Example:
```bash
alsv_retilelas --indir alsProjectDir --outdir alsProjectDir-tiled --tilesize 1000
```

The resulting tiles are uncompressed LAS files, and are intended to be a temporary step between
the raw supplied files and the standardized, indexed LAZ files produced in the next step.

## Standardise the LAS/LAZ files
The supplied LAS/LAZ files are expected to be in a single directory. The first command
to run will ingest these, and re-write them into a different directory, with a spatial
index included.

Example:
```bash
alsv_lidar_standardisation --indir alsProjectDir --outdir alsProjectDir-std \
    --intilesize 1000 --epsg 28356 --ii mp --project myproj \
    --year 2025 --binsize 50
```

The input data coordinates are expected to be in UTM projection, in a single zone,
as specified by the EPSG number given on the command line. If no `--epsg` is given,
the first input file will be checked for projection information, but if given, it will
over-ride anything found in the LAS files.

The `--ii` option (and also `--ss` and `--pp`) relate to the naming of output files, which
follows a strict structure. For more details on this, including known values for these
options, the user should read [QVF Lidar filenames](https://jrsrp.gitlab.io/sys/meta_info/lidar_filename_codes/).

The output files are always compressed (i.e. with `.laz` extension). They will be written
with a minimum LAS format (currently defaults to LAS 1.4), and the data will be
indexed into square bins of the given size. The file names conform to an internal 
structure, which includes the coordinates of the top-left tile corner, the year of capture,
and the project name. 

All subsequent processing works with these standardised LAZ files.

For further details see [alsv_lidar_standardisation](alsv_lidar_standardisation.md).

## Add Height Above Ground
For every point return which is not classified as "Ground", calculate its height above
the interpolated ground surface, and store this in the LAZ file. This is then available to
subsequent steps which calculate vegetation metrics. This is known as Height Above Ground (HAG).

```bash
alsv_addhag --indir alsProjectDir-std --tilesize 1000
```

The LAZ files are updated in place.

For further detail see [alsv_addhag](alsv_addhag.md).

## Producing per-tile raster products
The next step is to produce a range of raster images from these standardised LAZ files.

These products include a raster of ground elevations (i.e. a Digital Elevation Model, or DEM),
and a number of vegetation characteristics in raster form, such as maximum height above ground
per pixel, interpolated maximum height (canopy surface), Foliage Projective Cover (FPC), and
several others.

These are produced for each LAZ tile, in a subdirectory under the directory containing the
standardised LAZ files, each subdirectory being named to match its corresponding LAZ file
(without the `.laz` suffix). 

The various interpolations from irregular point data are done using Natural Neighbour interpolation,
at the level of individual index bins, and will use data from surrounding bins (including from
neighbouring tiles) as required. 

```bash
alsv_tile_products --indir alsProjectDir-std --tilesize 1000 --pixsize 0.5
```

For further detailes see [alsv_tile_products](alsv_tile_products.md).

## Canopy Height Model
Some of the products created above represent the top of the vegetation canopy, but have holes
and defects. Also, they operate on the Z coordinate, rather than height above ground (HAG).

The most sophisticated algorithm we have implemented for canopy height is the pit-free
canopy height model of Khosravipour et. al (2014). This is fairly slow, taking somewhat more time
than all the other tile products produced above, and not all users need it, so this is implemented
in a separate command

```bash
alsv_canopyheightmodel --indir alsProjectDir-std --tilesize 1000
```
This writes output raster tiles of just this Canopy Height Model (CHM) product.

A Khosravipour, AK Skidmore, M Isenburg, T Wang & YA Hussin (2014). Generating Pit-free Canopy Height
Models from Airborne Lidar. Photogrammetric Engineering & Remote Sensing 80(9), pp 863-872.
doi: 10.14358/PERS.80.9.863

## Fill any DEM gaps
For a range of reasons, the resulting DEM can sometimes contain holes which could not be
interpolated from the point data. One common reason for this is a water body which extends
across multiple index bins. A separate program is provided which will interpolate from the
surrounding raster elevation values, in order to fill any remaining gaps.

```bash
alsv_dem_gapfill --indir alsProjectDir-std --pixsize 0.5 --tilesize 1000
```

Any tiles which do have gaps will be filled in and written as a separate file, with a
`gapfilleddem` tag in the name. If present, these will be used by the next step, in preference
to the unfilled version.

## Produce project mosaics
The raster tiles can be mosaiced together to produce raster products for the whole project area

```bash
alsv_product_mosaic --indir alsProjectDir-std
```
This will produce a number of output rasters, writing into the `indir` directory. By default
the output files are GeoTiff files with CloudOptimizedTiff (COG) layout (different raster
formats can be used).

Currently there are 18 output products. The file names include a stage code indicating which
product it is. The various fields are separated by underscores. For example, the DEM
is named something like
```
apmpdr_rbrisba_2014_bb0m6_r50cm.tif
```
The stage code for the DEM is `bb0`. The 5th field, tagged with `r`, is for the resolution
(i.e. pixel size), in this case 50cm. The project name is `brisba` and the year of acquisition
is 2014.

See [stageCodes.md](stageCodes.md) for all the products and stage codes, and a reference to the
file naming conventions.

## Hydrological flow accumulation
A simple wrapper script is provided which uses the RichDEM package to compute hydrological flow
over the computed DEM, allowing the mapping of stream lines.

Most importantly, this wrapper includes a `--degradefactor` argument, allowing the flow accumulation
to be performed on a reduced resolution version of the DEM, which gives substantial improvements
in speed and memory use. Working with a very high resolution DEM for this can require very large
amounts of memory and compute time.

This command would normally be run on the mosaic of the whole project, to avoid mis-matched stream lines
at tile boundaries.

Example:
```bash
alsv_flowacc apmpdr_rbrisba_2014_bb0m6_r50cm.tif
```
will produce a flow accumulation raster `apmpdr_rbrisba_2014_bbqm6_r50cm.tif`

