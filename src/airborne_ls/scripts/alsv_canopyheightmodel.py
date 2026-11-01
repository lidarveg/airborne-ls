#!/usr/bin/env python
"""
Calculate Pit-free Canopy Height Model for one or more LAS/LAZ files. Output raster
of the resulting smooth, continuous canopy height.

"""
import os
import argparse
import glob

import numpy as np
import laspy
from scipy import ndimage

from airborne_ls import (
    const, filenaming_methods, gridding_methods, lazfile_rw, qvf, rw_image_methods
)


# We will count return classes up to this class
CLASSCOUNTS_MAXCLASS = const.PTCLASS_TEMPORALEXCLUSION
# Product name string for canopy height model
CHM_PRODUCT_NAME = "chm"


def getCmdargs():
    """
    Get command line arguments
    """
    p = argparse.ArgumentParser(description="""
        Calculate pit-free canopy height model, based on Khosravipour, A, et. al. (2014),
        Generating Pit-free Canopy Height Models from Airborne Lidar.
        Photogrammetric Engineering & Remote Sensing 80(9), pp 863-872.
        doi: 10.14358/PERS.80.9.863

        One or more input LAZ files are processed to corresponding raster output images. The
        LAZ files must be as output by alsv_lidar_standardisation, with the Height Above Ground
        added using the alsv_addhag command.
    """)
    p.add_argument("--indir", help="Directory of LAZ files to process")
    p.add_argument("--infile", help="Individual LAZ file to process")
    p.add_argument("--groupMofN", nargs=2, metavar=('M', 'N'), type=int,
        help=("Use only with --indir. Divide the set of tiles to be processed into N groups, " +
              "and then process only the files in the M-th group (group numbering starts at 1). " +
              "This helps support efficient batch processing, where the user decides how " +
              "many batch jobs will run, and the command for each group processes only " +
              "the tiles for that group. For example, with 5 groups, the first group " +
              "would be specified as '--groupMofN 1 5'"))
    p.add_argument("--tilesize", default=1000, type=int,
        help=("Tile size (metres). We need this so we can find neighbouring tiles " +
              "(default=%(default)s)"))
    p.add_argument("--pixsize", default=0.2, type=float,
        help="Pixel size of the CHM output raster (metres) (default=%(default)s)")
    p.add_argument("--binmargin", type=int, default=33,
        help=("Percentage of points from neighbouring bins to keep for per-bin " +
            "DEM interpolation (default=%(default)s). Smaller values will run faster, " +
            "but too small can leave extra holes in DEM"))
    p.add_argument("--minhag", type=float, default=0,
        help=("Minimum acceptable value of HAG (metres). Points below this are excluded. " +
              "(default=%(default)s)"))
    p.add_argument("--maxhag", type=float, default=120,
        help=("Maximum acceptable value of HAG (metres). Points above this are excluded. " +
              "(default=%(default)s)"))
    p.add_argument("--noparallelread", default=False, action="store_true",
        help=("Disable laspy's parallel reading of LAS/LAZ files. By default, laspy will " +
              "try to use as many CPU cores as it can"))
    p.add_argument("--driver", default='GTiff',
        help="GDAL driver for output image format (default=%(default)s)")
    p.add_argument("--skipexisting", default=False, action="store_true",
        help=("Skip input file if the output file already exists. Default will re-create " +
              "the output file, regardless of existence"))

    cmdargs = p.parse_args()

    if (round(cmdargs.tilesize / cmdargs.pixsize) * cmdargs.pixsize) != cmdargs.tilesize:
        msg = f"Tile size {cmdargs.tilesize} is not divisible by pixel size {cmdargs.pixsize}"
        raise ValueError(msg)

    return cmdargs


def main():
    """
    Main routine
    """
    cmdargs = getCmdargs()

    if cmdargs.infile is not None:
        filelist = [cmdargs.infile]
    elif cmdargs.indir is not None:
        pattern = f"{cmdargs.indir}/*.laz"
        filelist = sorted(glob.glob(pattern))
        if cmdargs.groupMofN is not None:
            (group, numGroups) = tuple(cmdargs.groupMofN)
            filelist = filenaming_methods.filelistGroupSubset(filelist, group, numGroups)

    for filename in filelist:
        outputDir = qvf.setsuffix(filename, "")
        if not os.path.isdir(outputDir):
            os.mkdir(outputDir)

        tileBasename = os.path.basename(outputDir)
        outfile = os.path.join(outputDir, tileBasename)
        stageCode = filenaming_methods.stageByProductName[CHM_PRODUCT_NAME]
        outfile = qvf.setstagecode(outfile, stageCode)
        outfile = qvf.setoptionfield(outfile, 'l', CHM_PRODUCT_NAME)
        resStr = filenaming_methods.resolutionStrFromMetres(cmdargs.pixsize)
        outfile = qvf.setoptionfield(outfile, 'r', resStr)
        suffix = filenaming_methods.getSuffixFromDriverName(cmdargs.driver)
        outfile = qvf.setsuffix(outfile, suffix)

        if cmdargs.skipexisting and os.path.exists(outfile):
            print(f"{filename} - skipping, outfile already exists")
        else:
            print(filename)

            processOneTile(filename, cmdargs, outfile)


def processOneTile(filename, cmdargs, outfile):
    """
    Process canopy height model for a single tile

    Parameters:
      filename (str): Name of input LAZ file
      cmdargs (argparse.Namespace): Command arguments object
      outfile (str): Name of output raster file
    """
    nullVal = -999.0
    pixSize = cmdargs.pixsize

    # Get the projection, preserved in output raster
    crs = laspy.open(filename).header.parse_crs()
    if crs is None:
        raise ValueError(f"File {filename} has no CRS (i.e. map projection)")
    epsg = crs.to_epsg()

    # Output raster array shape
    tileSizePix = int(cmdargs.tilesize / pixSize)
    tileShape = (tileSizePix, tileSizePix)
    chmTile = np.zeros(tileShape, dtype=np.float32)
    nonGroundClassesTile = np.full(tileShape, 254, dtype=np.uint8)

    binnedData = lazfile_rw.readBinnedData(filename,
        tileSize=cmdargs.tilesize, withNeighbours=True,
        noParallelRead=cmdargs.noparallelread)
    where = qvf.getwhere(filename)
    (easting, northing, utmZone) = filenaming_methods.decomposeWhereField(where)

    # The number of pixel rows/cols in each bin
    binSize = binnedData.ndx.binSize
    nRowsPerBin = int(np.ceil(binSize / pixSize))
    nColsPerBin = int(np.ceil(binSize / pixSize))

    binList = binnedData.getBinList()
    for (binRow, binCol) in binList:
        binTopLeftX = easting + binSize * binCol
        binTopLeftY = northing - binSize * binRow

        # Pixel row/col of this bin (top-left pixel)
        binTopRow = binRow * nRowsPerBin
        binLeftCol = binCol * nColsPerBin
        # A slice which locates the bin within the tile
        binSlice = (slice(binTopRow, (binTopRow + nRowsPerBin)),
                    slice(binLeftCol, (binLeftCol + nColsPerBin)))

        # Point data for this bin, with and without neighbouring margin
        binPts = binnedData.getData(binRow, binCol)
        binPtsWithMarg = gridding_methods.getBinPtsWithMargin(binnedData,
            binRow, binCol, cmdargs.binmargin, binTopLeftX, binTopLeftY)

        validHag = ((binPts["HAG"] >= cmdargs.minhag) & (binPts["HAG"] <= cmdargs.maxhag))
        validHagWithMarg = ((binPtsWithMarg["HAG"] >= cmdargs.minhag) &
                            (binPtsWithMarg["HAG"] <= cmdargs.maxhag))

        maxH_grid = np.full((nRowsPerBin, nColsPerBin), nullVal,
                           dtype=np.float32)
        (pixRow, pixCol) = gridding_methods.xyToRowCol(
            binPts["X"][validHag], binPts["Y"][validHag],
            binTopLeftX, binTopLeftY, pixSize)

        gridding_methods.maxH_array(pixRow, pixCol, binPts["HAG"][validHag], maxH_grid)

        vegPtsMask = np.logical_and(
            binPtsWithMarg["CLASSIFICATION"] >= const.PTCLASS_LOWVEGETATION,
            binPtsWithMarg["CLASSIFICATION"] <= const.PTCLASS_HIGHVEGETATION,
            validHagWithMarg
        )

        hag = binPtsWithMarg[vegPtsMask]['HAG']
        if len(hag) > 3:
            chmVeg = gridding_methods.chm_alg(
                binPtsWithMarg[vegPtsMask], hag, maxH_grid, pixSize,
                binTopLeftX, binTopLeftY, binSize, nRowsPerBin, nColsPerBin, nullVal)
            chmTile[binSlice] = chmVeg

        # Count return classes per-pixel
        classCountsShape = (CLASSCOUNTS_MAXCLASS + 1, nRowsPerBin, nColsPerBin)
        classCounts = np.zeros(classCountsShape, dtype=np.uint16)
        gridding_methods.makeClassCounts(pixRow, pixCol, binPts["X"],
            binPts["Y"][validHag], binPts["Z"][validHag],
            binPts["CLASSIFICATION"][validHag], classCounts)
        nonGroundClasses = np.full((nRowsPerBin, nColsPerBin), 254, dtype=np.uint8)
        # nonGroundClasses is most common non-ground class in each pixel
        classCountsNonGround = np.copy(classCounts)
        classCountsNonGround[const.PTCLASS_GROUND] = 0
        nonGroundClasses = classCountsNonGround.argmax(axis=0)
        del classCountsNonGround
        nonGroundClassesTile[binSlice] = nonGroundClasses

    # Interpolating over buildings & water can be a problem - mask out affected pixels here
    veg_msk = np.logical_or(nonGroundClassesTile == const.PTCLASS_BUILDING,
                            nonGroundClassesTile == const.PTCLASS_WATER)
    struct2 = ndimage.generate_binary_structure(2, 2)
    veg_msk = ndimage.binary_dilation(veg_msk, structure=struct2)
    veg_msk = ndimage.binary_dilation(veg_msk, structure=struct2)
    chmTile[veg_msk] = 0

    # Smooth out with a median filter, ignoring zeros
    masked_array = np.ma.masked_equal(chmTile, 0)
    chmTile = ndimage.median_filter(masked_array, size=3)

    # NOTE: The two steps above, i.e. binary dilation of the building/water mask and the
    # median filter of the CHM values, both have the potential to create edge-of-tile effects.
    # These seem to be very minor, and are probably not a serious problem, but the only fix
    # would be to run these steps later, after the tiles have been mosaiced together. This seems
    # like a lot of unnecessary complexity, so for the moment we leave them here.
    # See also the notes in alsv_testall.checkCHM() on coping with the difference.

    # Anything which is not canopy is set to null
    chmTile[chmTile < 0.5] = nullVal

    # Now write the image to file
    rw_image_methods.writeImage(np.round(chmTile, 3), outfile,
        tlx=easting, tly=northing, binsize=pixSize, epsg=epsg,
        nullVal=nullVal, parent_file=filename, driverName=cmdargs.driver)


if __name__ == "__main__":
    main()
