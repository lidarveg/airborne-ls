#!/usr/bin/env python
"""
Fill in gaps in a DEM tile. Uses any surrounding tiles to tie the edges. The input
tiles are the DEM tiles created by run_tile_products.py, and the outputs are
corresponding tiles with gaps filled.

Gaps are areas filled with the null value, and typically result from interpolations
across areas with a low density of ground returns. One significant reason for this can
be water bodies noticeably larger than the bin size, another can be small areas on
bin edges where the neighbouring bin was lacking in sufficient ground returns.

"""
import sys
import os
import argparse
import glob

import numpy
from osgeo import gdal
from scipy import ndimage
import pynninterp

from airborne_ls import filenaming_methods, qvf, gridding_methods, rw_image_methods


def getCmdargs():
    """
    Get command line arguments
    """
    p = argparse.ArgumentParser(description="""
        Interpolate to fill in gaps (i.e. regions filled with null values) in DEM tiles.
        Outputs files are same names as input DEM files, with modified 'l' field.
    """)
    p.add_argument("--tilesize", type=int, default=1000,
        help=("Size of tiles (metres). Used to identify surrounding tile filenames " +
              "(default=%(default)s)"))
    p.add_argument("--skipexisting", default=False, action="store_true",
        help=("Skip existing gap-filled DEM tiles. Default will re-compute any " +
              "output files which already exist"))
    p.add_argument("--tilemargin", type=int, default=100,
        help=("Margin (in pixels) of data to bring in from surrounding tiles " +
              "(default=%(default)s)"))
    p.add_argument("--clumpborder", type=int, default=10,
        help=("Width (in pixels) of border around a clump of nulls. The pixels in this " +
              "border region are used to interpolate across the clump. " +
              "(default=%(default)s)"))

    singleFileGroup = p.add_argument_group("Processing a single DEM tile")
    singleFileGroup.add_argument("--demfile", help="A single specific DEM tile to process")

    indirGroup = p.add_argument_group("Processing a directory of DEM tiles")
    indirGroup.add_argument("--indir",
        help="A directory of standardized LAZ files, whose DEM tiles will be processed")
    indirGroup.add_argument("--pixsize", default=0.5, type=float,
        help=("Pixel size (metres) of DEM files. Used to identify DEM file names " +
              "when using --indir. (default=%(default)s)"))
    indirGroup.add_argument("--driver", default="GTiff",
        help=("Name of GDAL driver. Used to identify file names as well as to format " +
              "the output (default=%(default)s)"))
    cmdargs = p.parse_args()

    if cmdargs.demfile is not None and cmdargs.indir is not None:
        print("Use either --demfile or --indir, but not both", file=sys.stderr)
        sys.exit(1)

    return cmdargs


def main():
    """
    Main routine
    """
    cmdargs = getCmdargs()
    demfileList = getDemImageFiles(cmdargs)

    filledLabel = "gapfilleddem"
    hillshadeProduct = "demHS"
    hillshadeStage = filenaming_methods.stageByProductName[hillshadeProduct]

    for demfile in demfileList:
        demGapFilledFile = qvf.setoptionfield(demfile, 'l', filledLabel)
        # If the output for this file already exists, skip to the next file
        if os.path.exists(demGapFilledFile) and cmdargs.skipexisting:
            print("Skipping", demfile)
        else:
            inHdr = rw_image_methods.imgH(demfile)
            epsgNum = int(inHdr['sr'].GetAuthorityCode())
            (dem, tileSlice, nullVal) = readWithMargin(demfile, cmdargs.tilemargin)
            dem = fillGaps(dem, nullVal, cmdargs.clumpborder)
            # Strip off the margins, round, and back to float32
            dem = numpy.round(dem[tileSlice], 3).astype(numpy.float32)

            outDemHSfile = qvf.setstagecode(demfile, hillshadeStage)
            outDemHSfile = qvf.setoptionfield(outDemHSfile, 'l', hillshadeProduct)

            rw_image_methods.writeImage(dem, demGapFilledFile, driver=cmdargs.driver,
                tlx=inHdr['tlx'], tly=inHdr['tly'], binsize=inHdr['pixel_s'],
                epsg=epsgNum, nullVal=nullVal)

            # Re-create the hillshade image
            creationoptions = rw_image_methods.creationOptionsByDriver.get(cmdargs.driver, [])
            demOptions = gdal.DEMProcessingOptions(computeEdges=True,
                format=cmdargs.driver, creationOptions=creationoptions)
            gdal.DEMProcessing(outDemHSfile, demGapFilledFile, "hillshade", options=demOptions)


def fillGaps(dem, nullVal, border):
    """
    Fill in gaps (areas filled with null value) in the given dem array, by
    interpolating from the surrounding pixels. Return a filled copy of the array.

    Parameters:
      dem (array): 2D float array of DEM values
      nullVal (float): Null value

    Returns:
      (array): Filled copy of dem array
    """
    demCopy = dem.copy()
    (nRows, nCols) = dem.shape
    nullmask = (dem == nullVal)
    nullClumps = clump(nullmask, 0)
    valIndices = ndimage.value_indices(nullClumps, ignore_value=0)

    for clumpId in valIndices:
        (clumpRow, clumpCol) = valIndices[clumpId]
        (rowMin, rowMax) = (clumpRow.min(), clumpRow.max())
        (colMin, colMax) = (clumpCol.min(), clumpCol.max())
        left = max(0, colMin - border)
        right = min(nCols, colMax + border)
        top = max(0, rowMin - border)
        bottom = min(nRows, rowMax + border)
        demSubset = dem[top:bottom, left:right]
        # Row and column numbers for the rectangular subset region, within the full dem array
        (r, c) = numpy.mgrid[top:bottom, left:right]

        # Make a data mask which is just a border around the clump, of width
        # <margin> pixels
        clumpRowSubset = clumpRow - top
        clumpColSubset = clumpCol - left
        clumpSubset = numpy.zeros(demSubset.shape, dtype=numpy.uint8)
        clumpSubset[(clumpRowSubset, clumpColSubset)] = 1
        dilatemask = gridding_methods.circleLocs(border)
        clumpSubset = ndimage.binary_dilation(clumpSubset, dilatemask)
        clumpSubset[(clumpRowSubset, clumpColSubset)] = 0
        dataMask = ((clumpSubset == 1) & (demSubset != nullVal))

        cData = c[dataMask].astype(numpy.float64)
        rData = r[dataMask].astype(numpy.float64)
        zData = demSubset[dataMask].astype(numpy.float64)

        colRowNull = (numpy.vstack([clumpCol, clumpRow]).T).astype(numpy.float64)

        demInterp = pynninterp.NaturalNeighbourPts(cData, rData, zData, colRowNull)
        # Insert these values into the copy of the original dem array
        demCopy[(clumpRow, clumpCol)] = demInterp

    demCopy[numpy.isnan(demCopy)] = nullVal
    return demCopy


def clump(img, nullVal):
    """
    Returns an array of the same shape as img, with contiguous clumps
    of equal values each given a unique clump id. Contiguousness is tested
    against all eight surrounding neighbours.

    Only works on an integer-like input array.

    """
    shape = img.shape
    clumps = numpy.zeros(shape, dtype=numpy.uint32)
    # 8-way connectedness
    connect = numpy.ones((3, 3), dtype=numpy.uint8)

    clumpid = numpy.uint32(0)
    imgvals = numpy.unique(img)
    for val in imgvals:
        if val != nullVal:
            mask = (img == val)
            labelledmask = numpy.zeros(shape, dtype=numpy.int32)
            numObj = ndimage.label(mask, structure=connect, output=labelledmask)
            clumps = numpy.where(mask, labelledmask + clumpid, clumps)
            clumpid += numpy.uint32(numObj)

    return clumps.astype(numpy.uint32)


def getDemImageFiles(cmdargs):
    """
    Use the given command line argument to work out the list of DEM tile files to process

    Parameters:
      cmdargs (argparse.Namespace): Command line arguments object

    Returns:
      demfileList (list[str]): List of DEM image files to be processed
    """
    if cmdargs.demfile is not None:
        demfileList = [cmdargs.demfile]
    else:
        lazdir = cmdargs.indir
        lazfileList = sorted(glob.glob(f"{lazdir}/*.laz"))

        productName = "dem"
        suffix = filenaming_methods.getSuffixFromDriverName(cmdargs.driver)
        demStage = filenaming_methods.stageByProductName[productName]
        resStr = filenaming_methods.resolutionStrFromMetres(cmdargs.pixsize)
        demfileList = []
        for lazfile in lazfileList:
            subdir = qvf.setsuffix(lazfile, '')
            demfile = qvf.setstagecode(os.path.basename(lazfile), demStage)
            demfile = qvf.setoptionfield(demfile, 'l', productName)
            demfile = qvf.setoptionfield(demfile, 'p', qvf.getoptionfield(lazfile, 'p'))
            demfile = qvf.setoptionfield(demfile, 'r', resStr)
            demfile = qvf.setsuffix(demfile, suffix)
            demfile = os.path.join(subdir, demfile)
            if not os.path.exists(demfile):
                msg = f"DEM file {demfile} not found"
                raise FileNotFoundError(msg)

            demfileList.append(demfile)
    return demfileList


def neighbourTileFilename(tilefile, xOffset, yOffset):
    """
    Given the file name of a single tile of data, return the file name for
    a neighbouring tile, based on the xOffset & yOffset parameters.

    The X & Y offsets are taken to be in metres, and are exactly one tile in
    some direction. So, for example, if the tile size is 1000m, then xOffset of -1000
    would indicate the tile to the west, and a yOffset of +1000 would indicate the
    tile to the north.

    Parameters:
      tilefile (str): File name of a single tile (as created by makeFilenameForTileProduct)
      xOffset (int): Offset (metres) in X direction to top-left of neighbour tile
      yOffset (int): Offset (metres) in Y direction to top-left of neighbour tile

    Returns:
      nbrfile (str): File name of requested neighbouring tile
    """
    where = qvf.getwhere(tilefile)
    newWhere = filenaming_methods.neighbourTileWhere(where, xOffset, yOffset)
    nbrfileFull = qvf.setwhere(tilefile, newWhere)

    # We probably also need to change the where field in the directory name
    (nbrdir, nbrfile) = os.path.split(nbrfileFull)
    if qvf.isQvf(nbrdir):
        if where == qvf.getwhere(nbrdir):
            nbrdir = qvf.setwhere(nbrdir, newWhere)
            nbrfileFull = os.path.join(nbrdir, nbrfile)

    return nbrfileFull


if __name__ == "__main__":
    main()
