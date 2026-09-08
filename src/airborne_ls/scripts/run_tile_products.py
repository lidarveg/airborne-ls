#!/usr/bin/env python

"""
###################################################################################################
#
# Purpose: Operational processing of .las/.laz files for standardised productions.
#          This incorporates surrounding tiles to eliminate edge effects.
#
# Approach: Process all tiles individually and merge at the project scale,
#           ultimately leaving no intermediate products.
#
# Assumptions:
# - Northing and easting tile locations are included in the filename.
# - New sub-tile indexing has been applied to optimise memory usage.
#
# Note: The resolution and point density of the datasets will significantly influence
#       processing times. The computation of the CHM is slow, will seek to optimise.
#
####################################################################################################

example:

uv run python scripts/run_tile_products.py \
    --indir lidarveg_testing_data/Brisbane_2014_LGA_sub/indexed_tiles/ --epsg 28356 \
    --laz_flist laz_flist --tile_s 1000. --psize 0.5 --chm_psize 0.2

"""

import sys
import os
import argparse
import glob
import logging
from pathlib import Path

import laspy
import numpy as np
from scipy import ndimage
from osgeo import gdal

from airborne_ls import (
    const,
    filenaming_methods,
    fpc_method,
    gridding_methods,
    lazfile_rw,
    rw_image_methods,
)

logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.ERROR, format="%(asctime)s: %(name)20s: %(levelname)10s: %(message)s"
)


# Default constants
percentiles = [1, 5, 25, 50, 75, 95, 99]
rtnClassNull = 255
# We will count return classes up to this class
CLASSCOUNTS_MAXCLASS = const.PTCLASS_TEMPORALEXCLUSION


def getCmdargs(inputargs):
    """
    Parse command-line arguments for running the script.

    Parameters:
        inputargs (list): List of command-line arguments.

    Returns:
        Namespace: Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Process .las/.laz files for standardised productions."
    )

    # Input and output directories
    parser.add_argument("--indir",
        help="Directory containing standardised, fully indexed LAZ files to process")
    parser.add_argument("--infile",
        help="Name of a single input standardised LAZ file to process")
    # Will add a --skipexisting later ......
    parser.add_argument("--tile_s", type=float, help="XY dimensions of LAS tile in metres.")
    parser.add_argument("--epsg", type=int,
        help="EPSG number of map projection used in the input LAZ files")

    # Output resolutions
    parser.add_argument("--psize", default=0.5, type=float,
        help=("Pixel size of gridded DEM, Intensity, and maxH output layers (metres). " +
              "Default: %(default)s."))
    parser.add_argument("--ptile_s", default=5, type=float,
        help="Pixel size of percentile output layers (metres). Default: %(default)s.")
    parser.add_argument("--fpc_psize", default=10.0, type=float,
        help="Pixel size of the FPC layer (metres). Default: %(default)s.")
    parser.add_argument("--chm_psize", default=0.2, type=float,
        help="Pixel size of the CHM layer (metres). Default: %(default)s.")

    # File processing options
    parser.add_argument("--split_fpc", default=False, action=argparse.BooleanOptionalAction,
        help="Split flight lines for FPC calculations.")
    parser.add_argument("--binmargin", type=int, default=33,
        help=("Percentage of points from neighbouring bins to keep for per-bin " +
            "DEM interpolation (default=%(default)s). Smaller values will run faster, " +
            "but too small can leave extra holes in DEM"))
    parser.add_argument("--driver", default='GTiff',
        help="GDAL driver for output image format (default=%(default)s)")

    cmdargs = parser.parse_args(inputargs)

    if cmdargs.indir is not None and cmdargs.infile is not None:
        msg = "Use either --indir or --infile, but not both"
        raise ValueError(msg)

    if cmdargs.indir is None and cmdargs.infile is None:
        msg = "Must supply one of --indir or --infile"
        raise ValueError(msg)

    if cmdargs.epsg is None:
        raise ValueError("EPSG code must be supplied.")

    if cmdargs.tile_s is None:
        raise ValueError("No cmdargs.tile_s supplied.")

    return cmdargs


def run_tile_products(cmdargs):
    """
    Process LAS/LAZ tiles individually and write output layers for each tile to a new folder.

    Parameters:
        cmdargs (Namespace): Parsed command-line arguments containing processing parameters.
    """
    # Check if the supplied pixel sizes are divisible
    check_divisible([cmdargs.psize, cmdargs.ptile_s, cmdargs.fpc_psize])

    if cmdargs.infile is not None:
        infilelist = [cmdargs.infile]
    elif cmdargs.indir is not None:
        pattern = f"{cmdargs.indir}/*.laz"
        infilelist = sorted(glob.glob(pattern))

    nullVal = -999.0
    neigh8 = np.array(
        [[-1, 1], [0, 1], [1, 1], [-1, 0], [1, 0], [-1, -1], [0, -1], [1, -1]]
    )

    # Read the first file to get bin size and number of bins
    (binSize, nbins, _, _) = lazfile_rw.read_laz_index(infilelist[0], cmdargs.tile_s)

    logger.info(f"binSize = {binSize} and nbins = {nbins} for {infilelist[0]}")

    # Process each file in the batch
    for infileFull in infilelist:
        infile = os.path.basename(infileFull)
        print(infile)

        # Set up output directories and filenames
        outputDir = Path(infileFull).with_suffix("")
        tileBasename = Path(infile).with_suffix("")
        outputBasename = Path(outputDir).joinpath(tileBasename)
        if not Path.is_dir(outputDir):
            Path.mkdir(outputDir)
        outfnames = filenaming_methods.get_outfnames(
            outputBasename,
            psize=cmdargs.psize,
            ptile_s=cmdargs.ptile_s,
            fpc_psize=cmdargs.fpc_psize,
            chm_psize=cmdargs.chm_psize,
            pptiles=percentiles,
        )
        logger.info(f"LAZ file and header read in for {infileFull}")

        # arrays to store processing segments of tiles
        tileSizePix = int(cmdargs.tile_s / cmdargs.psize)
        tileShape = (tileSizePix, tileSizePix)
        tileSizeChmPix = int(cmdargs.tile_s / cmdargs.chm_psize)
        tileChmShape = (tileSizeChmPix, tileSizeChmPix)
        tileSizeFpcPix = int(cmdargs.tile_s / cmdargs.fpc_psize)
        tileFpcShape = (tileSizeFpcPix, tileSizeFpcPix)
        tileSizePctPix = int(cmdargs.tile_s / cmdargs.ptile_s)
        tilePctShape = (len(percentiles), tileSizePctPix, tileSizePctPix)

        demTile = np.full(tileShape, nullVal, dtype=np.float32)
        csmTile = np.full(tileShape, nullVal, dtype=np.float32)
        chmTile = np.zeros(tileChmShape, dtype=np.float32)
        maxhTile = np.full(tileShape, nullVal, dtype=np.float32)
        intensTile = np.full(tileShape, nullVal, dtype=np.float32)
        ptDenTile = np.zeros(tileShape, dtype=np.uint16)
        # rtnClassNull # testing addition of new code to infill holes
        grTile = np.full(tileShape, 254, dtype=np.uint8)
        non_grTile = np.full(tileShape, 254, dtype=np.uint8)
        pctTile = np.full(tilePctShape, nullVal, dtype=np.float32)
        fpcTile = np.full(tileFpcShape, rtnClassNull, dtype=np.float32)

        # read in data
        # bData stores the data in chunks/bins for rapid access
        bData, northing, easting, rowS = lazfile_rw.bin_data(
            infileFull, cmdargs.tile_s, nbins
        )

        logger.info(f"Indexed LAZ FILE read in for {infileFull}")

        nRows = int(np.ceil(binSize / cmdargs.psize))
        nCols = int(np.ceil(binSize / cmdargs.psize))

        for rowB in range(1, rowS + 1, 1):
            for colB in range(1, rowS + 1, 1):
                fnc = f"row_{rowB}_col_{colB}"
                binTopLeftX = easting + binSize * (colB - 1)
                binTopLeftY = northing - cmdargs.tile_s + binSize * rowB
                binChunk = np.copy(bData[fnc])  # reset as we don't need buffer
                chunk = np.copy(bData[fnc])
                if len(chunk) > 10:
                    for n8 in neigh8:
                        fn = f"row_{int(rowB + n8[0])}_col_{int(colB + n8[1])}"
                        nbrBinData = np.copy(bData[fn])
                        if (nbrBinData.shape)[0] > 1:
                            nbrBinData = trimNeighbourBin(nbrBinData, n8[0], n8[1],
                                cmdargs.binmargin, binTopLeftX, binTopLeftY, binSize)
                            chunk = np.concatenate((chunk, nbrBinData))

                    # run dem
                    grdhits = (chunk["CLASSIFICATION"] == const.PTCLASS_GROUND)
                    if np.sum(grdhits) > 10:
                        # Set up a slice object for the part of the main tile arrays
                        # covering the current bin. Applies only to arrays of tileShape
                        xst = int(int((colB - 1) * binSize) / cmdargs.psize)
                        yst = int(int((rowS - rowB) * binSize) / cmdargs.psize)
                        binSlice = (slice(yst, (yst + nRows)), slice(xst, (xst + nRows)))

                        xst_bin = easting + int((colB - 1) * binSize)
                        yst_bin = northing - int((rowS - rowB) * binSize)
                        if np.sum(grdhits) > 0:
                            dem = gridding_methods.makeDemTile(
                                chunk["X"][grdhits], chunk["Y"][grdhits],
                                chunk["Z"][grdhits], xst_bin, yst_bin, binSize,
                                cmdargs.psize, nullVal=nullVal)
                            demTile[binSlice] = dem

                        # Get pixel coords of each point return, within the array for the bin
                        (row, col) = gridding_methods.xyToRowCol(
                            binChunk["X"], binChunk["Y"], xst_bin, yst_bin,
                            cmdargs.psize)

                        # Count return classes per-pixel
                        classCountsShape = (CLASSCOUNTS_MAXCLASS + 1, nRows, nCols)
                        classCounts = np.zeros(classCountsShape, dtype=np.uint16)
                        gridding_methods.makeClassCounts(row, col, binChunk["X"],
                            binChunk["Y"], binChunk["Z"], binChunk["CLASSIFICATION"],
                            classCounts)

                        #######################################################################
                        # run csm
                        xValsA, yValsA, zValsA = gridding_methods.maxH_xyzLocs(
                            chunk["X"], chunk["Y"], chunk["Z"], cmdargs.psize
                        )
                        csm = gridding_methods.makeDemTile(
                            xValsA, yValsA, zValsA, xst_bin, yst_bin,
                            binSize, cmdargs.psize, nullVal=nullVal)
                        csm = csm - dem
                        csm[dem == nullVal] = nullVal
                        csmTile[binSlice] = csm
                        csmTile[csmTile < -5] = nullVal
                        #######################################################################
                        # interp to irregular grid
                        nonGround = (binChunk["CLASSIFICATION"] != const.PTCLASS_GROUND)

                        heightAboveGround = (
                            gridding_methods.createHeightAboveGround(nonGround,
                                binChunk["X"], binChunk["Y"], binChunk["Z"],
                                chunk["X"][grdhits], chunk["Y"][grdhits], chunk["Z"][grdhits])
                        )

                        pntIntensity = binChunk["INTENSITY"]
                        # pntClass = binChunk["CLASSIFICATION"]
                        zArr = np.full((nRows, nCols), nullVal, dtype=np.float32)
                        intensityAtMaxH = np.full((nRows, nCols), nullVal, dtype=np.int16)
                        # Everywhere that we actually have data, initialize to zero
                        zArr[row, col] = 0
                        intensityAtMaxH[row, col] = 0
                        xArr = np.full((nRows, nCols), nullVal, dtype=np.float32)
                        yArr = np.full((nRows, nCols), nullVal, dtype=np.float32)
                        haveGroundReturn = np.full((nRows, nCols), 254, dtype=np.uint8)
                        nonGroundClasses = np.full((nRows, nCols), 254, dtype=np.uint8)

                        gridding_methods.maxH_workflow_layers(
                            row, col, binChunk["X"], binChunk["Y"],
                            heightAboveGround, pntIntensity, binChunk["CLASSIFICATION"],
                            xArr, yArr, zArr, intensityAtMaxH, haveGroundReturn)

                        # nonGroundClasses is most common non-ground class in each pixel
                        classCountsNonGround = np.copy(classCounts)
                        classCountsNonGround[const.PTCLASS_GROUND] = 0
                        nonGroundClasses = classCountsNonGround.argmax(axis=0)
                        del classCountsNonGround

                        maxhTile[binSlice] = zArr
                        intensTile[binSlice] = intensityAtMaxH
                        grTile[binSlice] = haveGroundReturn
                        non_grTile[binSlice] = nonGroundClasses

                        #######################################################################
                        # CREATE Canopy Height Model
                        # refer Khosravipour_2014 pit-free
                        number_veg_rets = np.sum(
                            binChunk["CLASSIFICATION"] == const.PTCLASS_MEDIUMVEGETATION
                        ) + np.sum(
                            binChunk["CLASSIFICATION"] == const.PTCLASS_HIGHVEGETATION
                        )  # could drop / add classification value of 3 (i.e. low veg)
                        if number_veg_rets > 10:
                            (row_chm, col_chm) = gridding_methods.xyToRowCol(
                                binChunk["X"], binChunk["Y"], xst_bin, yst_bin,
                                cmdargs.chm_psize)
                            nRows_chm = int(np.ceil(binSize / cmdargs.chm_psize))
                            nCols_chm = int(np.ceil(binSize / cmdargs.chm_psize))
                            xst_chm = int(
                                int((colB - 1) * binSize) / cmdargs.chm_psize
                            )
                            yst_chm = int(
                                int((rowS - rowB) * binSize) / cmdargs.chm_psize
                            )

                            maxH_hag = np.full((nRows_chm, nCols_chm), nullVal,
                                               dtype=np.float32)
                            gridding_methods.maxH_array(
                                row_chm, col_chm, heightAboveGround, maxH_hag
                            )
                            maxH_vals = (
                                chunk["CLASSIFICATION"] <= const.PTCLASS_HIGHVEGETATION
                            )  # this includes unclassified returns.. not sure of zero?
                            chunk_hag = gridding_methods.createHeightAboveGround(
                                maxH_vals, chunk["X"], chunk["Y"], chunk["Z"],
                                chunk["X"][grdhits], chunk["Y"][grdhits], chunk["Z"][grdhits])
                            vals = np.logical_and(
                                chunk["CLASSIFICATION"] >= const.PTCLASS_LOWVEGETATION,
                                chunk["CLASSIFICATION"] <= const.PTCLASS_HIGHVEGETATION,
                            )
                            if np.sum(vals) > 5:
                                chmVeg = gridding_methods.chm_alg(
                                    chunk[vals], chunk_hag[vals], maxH_hag,
                                    cmdargs.chm_psize, xst_bin, yst_bin,
                                    binSize, nRows_chm, nCols_chm, nullVal)
                                chmTile[
                                    yst_chm : (yst_chm + nRows_chm),
                                    xst_chm : (xst_chm + nCols_chm),
                                ] = chmVeg

                        ##############################
                        # run pulse density density
                        density = np.zeros((nRows, nCols), dtype=np.uint32)
                        fstR = binChunk["RETURN_NUMBER"] == 1
                        if np.sum(fstR) > 0:
                            (row_1st, col_1st) = gridding_methods.xyToRowCol(
                                binChunk["X"][fstR], binChunk["Y"][fstR],
                                xst_bin, yst_bin, cmdargs.psize)
                            gridding_methods.fstR_density(
                                row_1st, col_1st,
                                binChunk["X"][fstR], binChunk["Y"][fstR],
                                density)
                            ptDenTile[binSlice] = density
                        ############################
                        # run percentiles
                        (percentile_arr, nRows_pct) = (
                            gridding_methods.doHeightPercentileOutputs_idx(
                                binChunk["X"], binChunk["Y"], xst_bin, yst_bin,
                                binSize, heightAboveGround, cmdargs.ptile_s,
                                pptiles=percentiles, nullVal=nullVal))
                        xst_pct = int(int((colB - 1) * binSize) / cmdargs.ptile_s)
                        yst_pct = int(int((rowS - rowB) * binSize) / cmdargs.ptile_s)
                        pctTile[
                            :,
                            yst_pct : (yst_pct + nRows_pct),
                            xst_pct : (xst_pct + nRows_pct),
                        ] = percentile_arr

                        ############################
                        # run FPC
                        xst_fpc = int(int((colB - 1) * binSize) / cmdargs.fpc_psize)
                        yst_fpc = int(int((rowS - rowB) * binSize) / cmdargs.fpc_psize)
                        nRows_fpc = int(np.ceil(binSize / cmdargs.fpc_psize))
                        # nCols_fpc = int(np.ceil(binSize / cmdargs.fpc_psize))

                        flightlines = fpc_method.check_pts_pulses(binChunk)
                        canopyThreshold = 1.7
                        hag = heightAboveGround

                        fpc = fpc_method.doFPC(xst_bin, yst_bin, binChunk,
                            flightlines, hag, cmdargs.fpc_psize,
                            binSize, cmdargs.split_fpc, canopyThreshold)
                        fpcTile[
                            yst_fpc : (yst_fpc + nRows_fpc),
                            xst_fpc : (xst_fpc + nRows_fpc),
                        ] = fpc
                        ############################
                        del chunk, binChunk, flightlines, hag
        del bData

        rw_image_methods.writeImage(
            np.round(demTile, 3),
            outfnames["dem"], tlx=easting, tly=northing,
            binsize=cmdargs.psize, epsg=cmdargs.epsg, nullVal=nullVal,
            parent_file=infileFull, driverName=cmdargs.driver)
        rw_image_methods.writeImage(
            np.round(csmTile, 3),
            outfnames["csm"], tlx=easting, tly=northing,
            binsize=cmdargs.psize, epsg=cmdargs.epsg, nullVal=nullVal,
            parent_file=infileFull, driverName=cmdargs.driver)
        #######################
        # interpolating over buildings can be a problem - msk out affected pixels here
        multi = cmdargs.psize / cmdargs.chm_psize
        veg_msk = (non_grTile == const.PTCLASS_BUILDING) * 1 + (non_grTile == const.PTCLASS_WATER) * 1
        veg_msk = ndimage.zoom(veg_msk, multi, order=0)
        struct2 = ndimage.generate_binary_structure(2, 2)
        veg_msk = ndimage.binary_dilation(veg_msk, structure=struct2)
        veg_msk = ndimage.binary_dilation(veg_msk, structure=struct2)
        chmTile[veg_msk > 0] = 0

        masked_array = np.ma.masked_equal(chmTile, 0)
        chmTile = ndimage.median_filter(
            masked_array, size=3
        )  # median filter ignoring zeros
        chmTile[chmTile < 0.5] = nullVal
        rw_image_methods.writeImage(
            np.round(chmTile, 3),
            outfnames["chm"], tlx=easting, tly=northing,
            binsize=cmdargs.chm_psize, epsg=cmdargs.epsg, nullVal=nullVal,
            parent_file=infileFull, driverName=cmdargs.driver)
        #########################
        rw_image_methods.writeImage(
            np.round(maxhTile, 3),
            outfnames["maxH"], tlx=easting, tly=northing,
            binsize=cmdargs.psize, epsg=cmdargs.epsg, nullVal=nullVal,
            parent_file=infileFull, driverName=cmdargs.driver)
        rw_image_methods.writeImage(
            np.round(intensTile, 4),
            outfnames["intens"], tlx=easting, tly=northing,
            binsize=cmdargs.psize, epsg=cmdargs.epsg, nullVal=nullVal,
            parent_file=infileFull, driverName=cmdargs.driver)
        rw_image_methods.writeImage(
            grTile, outfnames["grdR"], tlx=easting, tly=northing,
            binsize=cmdargs.psize, epsg=cmdargs.epsg, nullVal=rtnClassNull,
            parent_file=infileFull, overviewResampling="MODE",
            driverName=cmdargs.driver)
        rw_image_methods.writeImage(
            non_grTile, outfnames["NonGrdCodes"],
            tlx=easting, tly=northing, binsize=cmdargs.psize, epsg=cmdargs.epsg,
            nullVal=rtnClassNull, parent_file=infileFull, overviewResampling="MODE",
            driverName=cmdargs.driver)
        rw_image_methods.writeImage(
            ptDenTile, outfnames["fstDens"],
            tlx=easting, tly=northing, binsize=cmdargs.psize, epsg=cmdargs.epsg,
            nullVal=0, parent_file=infileFull, driverName=cmdargs.driver)
        for idx, pp in enumerate(percentiles):
            layer = pctTile[idx, :, :]
            productName = f"percentile{pp}"
            rw_image_methods.writeImage(
                np.round(layer, 3), outfnames[productName], tlx=easting, tly=northing,
                binsize=cmdargs.ptile_s, epsg=cmdargs.epsg, nullVal=nullVal,
                parent_file=infileFull, driverName=cmdargs.driver)
        rw_image_methods.writeImage(
            np.rint(fpcTile).astype(np.uint8), outfnames["fpc"],
            tlx=easting, tly=northing, binsize=cmdargs.fpc_psize, epsg=cmdargs.epsg,
            nullVal=rtnClassNull, parent_file=infileFull, driverName=cmdargs.driver)

        # Apply colour tables
        fpcClrTbl = fpc_method.makeFPCcolorTable()
        rw_image_methods.setColorTable(outfnames["fpc"], fpcClrTbl)
        nonGrdClrTbl = makeNonGroundClrTbl()
        rw_image_methods.setColorTable(outfnames["NonGrdCodes"], nonGrdClrTbl)

        logger.info(f"Tiles written to file {infileFull}")

        demOptions = gdal.DEMProcessingOptions(computeEdges=True)
        gdal.DEMProcessing(outfnames["demHS"], outfnames["dem"], "hillshade",
            options=demOptions)

        # calculate and write out to a temporary file return,pulse,area stats while
        # data held in memory
        with laspy.open(infileFull) as f:
            nReturns = f.header.point_count
        dem_area = np.count_nonzero(demTile != nullVal)
        nPulses = np.sum(ptDenTile)

        outfile = Path(
            str(outfnames["dem"]).replace(
                Path(outfnames["dem"]).suffix, "_tempStats.txt"
            )
        )

        with open(outfile, "w") as fout:
            fout.write(f"fname = {infile}\n")
            fout.write(f"nPulses = {nPulses}\n")
            fout.write(f"nReturns = {nReturns}\n")
            fout.write(f"DEM area (pixels) = {dem_area}\n")


###################################################################################################
def reorder_flist(lazlistfull):
    """
    Reorder a list of LAS/LAZ files by file size in descending order.

    Parameters:
        lazlistfull (list): List of file paths to LAS/LAZ files.

    Returns:
        np.ndarray: Array of file paths sorted by file size in descending order.
    """
    file_sizes = [Path(ff).stat().st_size for ff in lazlistfull]
    sorted_indices = np.argsort(-np.array(file_sizes))  # Sort in descending order
    return np.array(lazlistfull)[sorted_indices]


def check_divisible(psizes):
    """
    Check if the provided pixel sizes are divisible by common factors.

    Parameters:
        psizes (list): List of pixel sizes to check.

    Raises:
        SystemExit: If any of the pixel sizes are not divisible by the expected factors.
    """
    pName = ["psize", "ptile_s", "fpc_psize"]

    for loc, psize in enumerate(psizes):
        if psize <= 1.0:
            if (1.0 % psize) != 0:
                msg = f"error: cmdargs.{pName[loc]} {psize} not divisible"
                sys.exit(msg)
        elif psize <= 10.0:
            if (10.0 % psize) != 0:
                msg = f"error: cmdargs.{pName[loc]} {psize} not divisible"
                sys.exit(msg)
        else:
            if (100.0 % psize) != 0:
                msg = f"error: cmdargs.{pName[loc]} {psize} not divisible"
                sys.exit(msg)


def trimNeighbourBin(data, binRowOff, binColOff, binMargin, topLeftX, topLeftY, binSize):
    """
    Trim off the points in a neighbouring bin, so that we are left with only
    those points close to the central bin.

    Bin row/col offsets define which neighbour direction this bin lies from the
    central bin. They were added to the bin row/col number to get the neighour bin
    row/col. The LAZ file point index as presented with row/col numbering starting
    at 1 in the bottom-left bin, and increasing eastwards and northwards.

    Parameters:
      data: Point data for the whole of the bin to be trimmed
      binRowOff, binColOff: Bin row & col offsets, relative to central bin
      binMargin: Percentage of the data to keep in the trimmed data
      topLeftX, topLeftY: (X, y) coordinates of the top-left corner of the
                          central bin, in metres
      binSize: Size of the bin (along one edge) in metres

    Returns:
      trimmedData: The subset of the given point data which lies closest
                   to the central bin
    """
    # Point (X, Y) coords
    x = data['X']
    y = data['Y']

    # There ar two main cases. In one case, we want the whole of one edge of the bin,
    # on the other we want one corner of the bin. The first case will have one
    # of the offset values equal to zero, the second case will have both non-zero
    if 0 in (binRowOff, binColOff):
        # Choose the coordinate (either X or Y) on which to select points. We calculate
        # the perpendicular distance of each point from the relevant edge of the central bin
        if binColOff == -1:
            dist = topLeftX - x
        elif binColOff == 1:
            dist = x - (topLeftX + binSize)
        elif binRowOff == -1:
            dist = topLeftY - binSize - y
        elif binRowOff == 1:
            dist = y - topLeftY
    else:
        # The coordinate to select on is distance from the corner point, rather than a single
        # coordinate. Use the offset values to work out the corner (X, Y) coords, then
        # calculate Euclidean distance from that
        cnrX = topLeftX
        if binColOff == 1:
            cnrX = topLeftX + binSize
        cnrY = topLeftY
        if binRowOff == -1:
            cnrY = topLeftY - binSize

        # The Euclidean distance is the "coordinate" on which we will select points
        dist = np.sqrt((x - cnrX) ** 2 + (y - cnrY) ** 2)

        # Adjust the binMargin because we only want the corner. In effect we square the proportion
        binMargin = int(100 * (binMargin / 100) ** 2)

    # A mask for the points with distance smaller than the given percentile
    threshold = np.percentile(dist, binMargin)
    keepMask = (dist < threshold)

    trimmedData = data[keepMask]
    return trimmedData


def makeNonGroundClrTbl():
    """
    Create a colour table array for the nonGrdCodes image

    The returned array has shape (256, 3), and type uint8. The columns are
    red/green/blue values, in the range [0, 255]. This is suitable for use
    with the rw_image_methods.setColorTable function.

    Returns:
      clrTblArr: Array of RGB values
    """
    clrTblArr = np.full((256, 3), 0, dtype=np.uint8)
    clrTblArr[0] = (254, 254, 254)  # never classified:
    clrTblArr[1] = (200, 200, 200)  # unclassified: light gray
    clrTblArr[2] = (0, 0, 0)        # ground classification
    clrTblArr[3] = (0, 240, 0)      # low veg: green1
    clrTblArr[4] = (0, 160, 0)      # medium veg: green2
    clrTblArr[5] = (0, 80, 0)       # high veg: green3
    clrTblArr[6] = (255, 0, 0)      # building: red
    clrTblArr[7] = (255, 255, 0)
    clrTblArr[8] = (255, 255, 0)
    clrTblArr[9] = (0, 0, 255)      # blue for water
    clrTblArr[10] = (255, 0, 255)
    clrTblArr[11] = (255, 20, 255)
    clrTblArr[12] = (255, 30, 255)
    clrTblArr[13] = (255, 40, 255)
    clrTblArr[14] = (255, 50, 255)
    clrTblArr[15] = (255, 60, 255)
    clrTblArr[16] = (255, 70, 255)
    clrTblArr[17] = (255, 80, 255)
    clrTblArr[18] = (255, 90, 255)
    clrTblArr[19] = (255, 100, 255)
    clrTblArr[254] = (101, 67, 33)  # brown for background

    return clrTblArr


###############################################################################################
###############################################################################################
###############################################################################################


def main(args=None):
    """
    Main entry point for the script, allowing external calls.

    Parameters:
        args (list[str], optional): Command-line parameter list. Defaults to None.
    """
    if args is None:
        args = sys.argv[1:]

    logger.debug("Getting command-line arguments.")
    cmdargs = getCmdargs(args)
    logger.debug(f"Command-line arguments: {args}")
    logger.debug(f"Input directory: {cmdargs.indir}")

    # Run the ALS tile processing
    run_tile_products(cmdargs)


if __name__ == "__main__":
    main()
