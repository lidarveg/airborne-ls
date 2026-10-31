# alsv_lidar_standardisation
Re-writes the input LAS/LAZ files as compressed LAZ files, with an up-to-date file version,
and an internal spatial index.

## Filtering
A number of filters are applied to the data.

  * If `--outtilesize` is different to `--intilesize`, then each tile is split into the
    required sub-tiles. Only points whose (x, y) coordinates are within the tile are included
    in the output tile. This option only allows splitting input tiles into smaller tiles, and is not
    a generic re-tiling. The output tile size must be an exact divisor of the input tile size. For
    a more general re-tiling, see `alsv_retilelaz`. 
  * The `--excludeclasses` option gives point classification values which are excluded from
    the output. This has some defaults, but can be over-ridden from the command line. 
  * The `--minz` & `--maxz` options give a range of acceptable elevation values. Any points
    which have a Z value outside this range are assumed to be errors and excluded from the output.
  * If `--epsg` is given, this will be written into the header. If not given, it will copy any
    EPSG from the input, and if nothing is found an exception is raised.

## File naming
The output LAZ files are named according to the conventions of
[QVF Lidar filenames](https://jrsrp.gitlab.io/sys/meta_info/lidar_filename_codes/). This file naming
structure forms the basis of all file names for this and other commands, and is particularly
important in allowing some commands to identify which tile is a direct neighbour of another.

## Spatial Index
The spatial index is based on dividing the whole tile area into small square bins. The bins are ordered
in row-major order, meaning that the top-left bin is bin 0, the bin to its right is bin 1,
across the top row of bins. The next bin is the left-most bin in the next row down, and so on
to the bottom-right bin.

The points are then sorted so that all points for bin 0 are first, then all points in bin 1, etc.
Within the bin, order is preserved as it was supplied. The index is constructed simply as
the point number of the first point in each bin, and finished with the total number of points.
This is stored in a pair of VLR records in the output LAZ file, with the user_id &
record_id values given in `airborne_ls.const`. Selecting points for a given bin requires just
the range between two consecutive index values.


