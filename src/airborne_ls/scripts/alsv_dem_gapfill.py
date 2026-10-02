#!/usr/bin/env python
"""
Fill in gaps in a DEM tile. Uses any surrounding tiles to tie the edges. The input
tiles are the DEM tiles created by alsv_tile_products.py, and the outputs are
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

import numpy as np
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
    p.add_argument("--tilemargin", type=int, default=200,
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
    indirGroup.add_argument("--groupMofN", nargs=2, metavar=('M', 'N'), type=int,
        help=("Use only with --indir. Divide the set of tiles to be processed into N groups, " +
              "and then process only the files in the M-th group (group numbering starts at 1)." +
              "This helps support efficient batch processing, where the user decides how " +
              "many batch jobs will run, and the command for each group processes only " +
              "the tiles for that group. For example, with 5 groups, the first group " +
              "would be specified as '--groupMofN 1 5'"))
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

    if cmdargs.indir is None and cmdargs.groupMofN is not None:
        print("Using --groupMofN requires --indir", file=sys.stderr)
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
            epsgNum = int(inHdr['sr'].GetAuthorityCode('PROJCS'))
            (dem, tileSlice, nullVal) = readWithMargin(demfile, cmdargs)
            demFilled = fillGaps(dem, nullVal, cmdargs.clumpborder, tileSlice)
            # Strip off the margins, back to float32, and round to 3 places
            demFilled = demFilled[tileSlice]
            demFilled = demFilled.astype(np.float32)
            demFilled = np.round(demFilled, 3)
            tileChanged = (demFilled != dem[tileSlice]).any()

            if tileChanged:
                rw_image_methods.writeImage(demFilled, demGapFilledFile, driverName=cmdargs.driver,
                    tlx=inHdr['tlx'], tly=inHdr['tly'], binsize=inHdr['pixel_s'],
                    epsg=epsgNum, nullVal=nullVal)

                # Re-create the hillshade image
                outDemHSfile = qvf.setstagecode(demfile, hillshadeStage)
                outDemHSfile = qvf.setoptionfield(outDemHSfile, 'l', hillshadeProduct)
                creationoptions = rw_image_methods.creationOptionsByDriver.get(cmdargs.driver, [])
                demOptions = gdal.DEMProcessingOptions(computeEdges=True,
                    format=cmdargs.driver, creationOptions=creationoptions)
                gdal.DEMProcessing(outDemHSfile, demGapFilledFile, "hillshade", options=demOptions)


def fillGaps(dem, nullVal, border, tileSlice):
    """
    Fill in gaps (areas filled with null value) in the given dem array, by
    interpolating from the surrounding pixels. Return a filled copy of the array.

    Parameters:
      dem (array): 2D float array of DEM values
      nullVal (float): Null value
      border (int): Size (pixels) of border around a clump to be filled
      tileSlice (slice tuple): Index of sub-array of dem corresponding to
                               original tile

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

        # Check whether this clump is one we should try to fill
        outsideTile = clumpBBoxOutsideTile(tileSlice, rowMin, rowMax, colMin, colMax)
        intersectsTile = False
        if not outsideTile:
            intersectsTile = clumpIntersects(clumpRow, clumpCol, dem.shape, tileSlice)
        reachesEdge = ((rowMin == 0) or (rowMax == dem.shape[0]) or
                       (colMin == 0) or (colMax == dem.shape[1]))
        toBeFilled = (not outsideTile) and intersectsTile and (not reachesEdge)

        if toBeFilled:
            left = max(0, colMin - border)
            right = min(nCols, colMax + border)
            top = max(0, rowMin - border)
            bottom = min(nRows, rowMax + border)
            demSubset = dem[top:bottom, left:right]
            # Row and column numbers for the rectangular subset region, within the full dem array
            (r, c) = np.mgrid[top:bottom, left:right]

            # Make a data mask which is just a border around the clump, of width
            # <margin> pixels
            clumpRowSubset = clumpRow - top
            clumpColSubset = clumpCol - left
            clumpSubset = np.zeros(demSubset.shape, dtype=np.uint8)
            clumpSubset[(clumpRowSubset, clumpColSubset)] = 1
            dilatemask = gridding_methods.circleLocs(border)
            clumpSubset = ndimage.binary_dilation(clumpSubset, dilatemask)
            clumpSubset[(clumpRowSubset, clumpColSubset)] = 0
            dataMask = ((clumpSubset == 1) & (demSubset != nullVal))

            cData = c[dataMask].astype(np.float64)
            rData = r[dataMask].astype(np.float64)
            zData = demSubset[dataMask].astype(np.float64)

            colRowNull = (np.vstack([clumpCol, clumpRow]).T).astype(np.float64)

            demInterp = pynninterp.NaturalNeighbourPts(cData, rData, zData, colRowNull)
            # Insert these values into the copy of the original dem array
            demCopy[(clumpRow, clumpCol)] = demInterp

    demCopy[np.isnan(demCopy)] = nullVal
    return demCopy


def clumpBBoxOutsideTile(tileSlice, rowMin, rowMax, colMin, colMax):
    """
    Check if the clump bounding box is outside the tile bounds

    Parameters:
      tileSlice (slice tuple): Index of sub-array of dem corresponding to
                               original tile
      rowMin, rowMax, colMin, colMax (int): Bounding box within array of the
                                            clump

    Returns:
      (bool): True if clump bounding box is outside the tile bounds
    """
    (tileRowMin, tileRowMax) = (tileSlice[0].start, tileSlice[0].stop)
    (tileColMin, tileColMax) = (tileSlice[1].start, tileSlice[1].stop)

    outside = ((rowMax < tileRowMin) or (rowMin > tileRowMax) or
               (colMax < tileColMin) or (colMin > tileColMax))
    return outside


def clumpIntersects(clumpRow, clumpCol, demShape, tileSlice):
    """
    If the clump may intersect, this checks every pixel to see if it really does.

    Parameter:
      clumpRow, clumpCol (index tuple): Indices of clump in full dem array
      demShape (tuple): Shape of full dem array
      tileSlice (slice tuple): Slice for tile within full dem array

    Returns:
      (bool): True if clump intersects original tile
    """
    # Empty version of full array
    arr = np.zeros(demShape, dtype=np.uint8)
    # Set clump pixels to 1
    arr[(clumpRow, clumpCol)] = 1
    # Count how may are in original tile
    count = np.count_nonzero(arr[tileSlice])

    intersects = (count > 0)
    return intersects


def clump(img, nullVal):
    """
    Returns an array of the same shape as img, with contiguous clumps
    of equal values each given a unique clump id. Contiguousness is tested
    against all eight surrounding neighbours.

    Only works on an integer-like input array.

    """
    shape = img.shape
    clumps = np.zeros(shape, dtype=np.uint32)
    # 8-way connectedness
    connect = np.ones((3, 3), dtype=np.uint8)

    clumpid = np.uint32(0)
    imgvals = np.unique(img)
    for val in imgvals:
        if val != nullVal:
            mask = (img == val)
            labelledmask = np.zeros(shape, dtype=np.int32)
            numObj = ndimage.label(mask, structure=connect, output=labelledmask)
            clumps = np.where(mask, labelledmask + clumpid, clumps)
            clumpid += np.uint32(numObj)

    return clumps.astype(np.uint32)


def getDemImageFiles(cmdargs):
    """
    Use the given command line arguments to work out the list of DEM tile files to process

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

        if cmdargs.groupMofN is not None:
            (group, numGroups) = tuple(cmdargs.groupMofN)
            demfileList = filenaming_methods.filelistGroupSubset(demfileList, group, numGroups)

    return demfileList


def readWithMargin(demfile, cmdargs):
    """
    Read the DEM data for the given tile, augmented with a margin of data from
    surrounding tiles, where available.
    """
    tilesize = cmdargs.tilesize
    margin = cmdargs.tilemargin
    nbrOffset = {
        'N': (0, tilesize), 'E': (tilesize, 0),
        'S': (0, -tilesize), 'W': (-tilesize, 0),
        'NE': (tilesize, tilesize), 'SE': (tilesize, -tilesize),
        'SW': (-tilesize, -tilesize), 'NW': (-tilesize, tilesize)
    }

    # Find the existing neighbour files for each direction
    nbrTile = {}
    for key in nbrOffset:
        (xOffset, yOffset) = nbrOffset[key]
        filename = neighbourTileFilename(demfile, xOffset, yOffset)
        if os.path.exists(filename):
            nbrTile[key] = filename
    # Only keep diagonal neighbours if both its adjacent direct neighbours also exist
    for key in ['NE', 'SE', 'SW', 'NW']:
        if key in nbrTile and (key[0] not in nbrTile or key[1] not in nbrTile):
            # Missing at least one direct neighbour, so remove this diagonal
            nbrTile.pop(key)

    # Header of central tile
    tileInfo = rw_image_methods.imgH(demfile)

    # The tile shape, assumed to be the same for all tiles
    (nRowsPerTile, nColsPerTile) = (tileInfo['ydim'], tileInfo['xdim'])
    # The shape of the output array, and the position of the central tile within that array
    (nRows, nCols) = (nRowsPerTile, nColsPerTile)
    (ctrTop, ctrBottom, ctrLeft, ctrRight) = (0, nRows, 0, nCols)
    if 'N' in nbrTile:
        nRows += margin
        (ctrTop, ctrBottom) = (ctrTop + margin, ctrBottom + margin)
    if 'S' in nbrTile:
        nRows += margin
    if 'E' in nbrTile:
        nCols += margin
    if 'W' in nbrTile:
        nCols += margin
        (ctrLeft, ctrRight) = (ctrLeft + margin, ctrRight + margin)

    # For each possible neighbour, the row/cols to read, and where to place it in
    # the output array. Each value is a tuple (fpos, apos).
    # The fpos tuple is for reading from the file (xoff, yoff, win_xsize, win_ysize)
    # The apos tuple is location within the output array (startRow, endRow, startCol, endCol)
    nbrPos = {
        'N': ((0, (nRowsPerTile - margin), nColsPerTile, margin), (0, margin, ctrLeft, ctrRight)),
        'E': ((0, 0, margin, nRowsPerTile), (ctrTop, ctrBottom, ctrRight, nCols)),
        'S': ((0, 0, nColsPerTile, margin), (ctrBottom, nRows, ctrLeft, ctrRight)),
        'W': ((nColsPerTile - margin, 0, margin, nRowsPerTile), (ctrTop, ctrBottom, 0, ctrLeft)),
        'NE': ((0, nRowsPerTile - margin, margin, margin), (0, margin, ctrRight, nCols)),
        'SE': ((0, 0, margin, margin), (ctrBottom, nRows, ctrRight, nCols)),
        'SW': ((nColsPerTile - margin, 0, margin, margin), (ctrBottom, nRows, 0, margin)),
        'NW': ((nColsPerTile - margin, nRowsPerTile - margin, margin, margin),
               (0, ctrTop, 0, ctrLeft))
    }

    # Read in the central tile
    ds = gdal.Open(demfile)
    band = ds.GetRasterBand(1)
    nullVal = band.GetNoDataValue()
    dem = np.full((nRows, nCols), nullVal, dtype=np.float32)
    dem[ctrTop:ctrBottom, ctrLeft:ctrRight] = band.ReadAsArray().astype(np.float32)
    del band, ds

    # Now read in all margins
    for key in nbrTile:
        (fpos, apos) = nbrPos[key]
        nbrDs = gdal.Open(nbrTile[key])
        nbrBand = nbrDs.GetRasterBand(1)
        arr = nbrBand.ReadAsArray(fpos[0], fpos[1], fpos[2], fpos[3])
        dem[apos[0]:apos[1], apos[2]:apos[3]] = arr

    # A slice object to select the central tile from the augmented array
    tileSlice = (slice(ctrTop, ctrBottom), slice(ctrLeft, ctrRight))

    return (dem, tileSlice, nullVal)


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

    # We also need to change the where field in the directory name
    (nbrdir, nbrfile) = os.path.split(nbrfileFull)
    if qvf.isQvf(nbrdir):
        if where == qvf.getwhere(nbrdir):
            nbrdir = qvf.setwhere(nbrdir, newWhere)
            nbrfileFull = os.path.join(nbrdir, nbrfile)

    return nbrfileFull


if __name__ == "__main__":
    main()
