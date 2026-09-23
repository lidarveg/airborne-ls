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

The output files are always compressed (i.e. with `.laz` extension). They will be written
with a minimum LAS format (currently defaults to LAS 1.4), and the data will be
indexed into square bins of the given size. The file names conform to an internal 
structure, which includes the coordinates of the top-left tile corner, the year of capture,
and the project name. 

All subsequent processing works with these standardised LAZ files.

## Producing per-tile raster products
The next step is to produce a range of raster images from these standardised LAZ files.

These products include a raster of ground elevations (i.e. a Digital Elevation Model, or DEM),
and a number of vegetation characteristics in raster form, such as a Canopy Height Model (CHM),
Foliage Projective Cover (FPC), and several others. These are produced for each LAZ tile,
in a subdirectory under the directory containing the standardised LAZ files, each subdirectory
being named to match its corresponding LAZ file (without the `.laz` suffix). 

The various interpolations from irregular point data are done using Natural Neighbour interpolation,
at the level of individual index bins, and will use data from surrounding bins (including from
neighbouring tiles) as required. 

```bash
alsv_tile_products --indir alsProjectDir-std --tilesize 1000 --pixsize 0.5
```

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

Currently there are 17 output products. The file names include a stage code indicating which
product it is. The various fields are separated by underscores. For example, the DEM
is named something like
```
apmpdr_rbrisba_2014_bb0m6_r50cm.tif
```
The stage code for the DEM is `bb0`. The 5th field, tagged with `r`, is for the resolution
(i.e. pixel size), in this case 50cm. The project name is `brisba` and the year of acquisition
is 2014.

See [stageCodes.md](stageCodes.md) for all the products and stage codes.

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

