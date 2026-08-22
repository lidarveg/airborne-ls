#!/usr/bin/env python

"""
This script performs DEM infill for missing areas in LAS files and converts the results to
Cloud Optimized GeoTIFF (COG) format.

Needs refinement / further development

example:  uv run python scripts/run_dem_correction.py \
            --indir lidarveg_testing_data/Brisbane_2014_LGA_sub/indexed_tiles/ \
            --laz_flist  laz_flist2 --tile_s 1000. --psize 0.5 --epsg 28356

"""
import os
import argparse
import logging
import subprocess
import sys
import glob
from pathlib import Path

import numpy as np
from osgeo import gdal

from airborne_ls import gridding_methods, rw_image_methods, filenaming_methods, qvf

# Configure logging
logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.ERROR, format="%(asctime)s: %(name)20s: %(levelname)10s: %(message)s"
)


def getCmdargs(inputargs):
    """
    Parse command-line arguments for running as a standalone script.

    Parameters:
        inputargs (list): List of command-line arguments.

    Returns:
        argparse.Namespace: Parsed arguments as a namespace object.
    """
    parser = argparse.ArgumentParser(description="Batch processing for LAS files.")

    parser.add_argument("--indir", required=True, help="Directory containing LAS files.")
    parser.add_argument("--laz_flist", required=True, help="List of LAS files to be processed.")
    parser.add_argument("--epsg", type=int, required=True,
        help="EPSG code for map information.")
    parser.add_argument("--psize", default=0.5, type=float,
        help=("Pixel size of gridded DEM, Intensity, and maxH output layers (metres). " +
              "Default: %(default)s."))
    parser.add_argument("--tile_s", type=float, required=True,
        help="XY dimensions of LAS tile in metres.")
    parser.add_argument("--startfilenum", default=0, type=int,
        help="Position within laz_flist to start batch processing. Default: %(default)s.")
    parser.add_argument("--stopfilenum", type=int,
        help="Position within laz_flist to stop batch processing.")

    return parser.parse_args(inputargs)


def run_dem_correction(cmdargs):
    """
    Perform DEM infill for missing areas in LAS files.

    Parameters:
        cmdargs (argparse.Namespace): Parsed command-line arguments containing input directory,
                                      file list, and processing parameters.
    """
    # Read input file list
    fn = Path(cmdargs.indir).joinpath(cmdargs.laz_flist)
    with open(fn) as f:
        infiles = [line.strip() for line in f]
    infilesFull = [os.path.join(cmdargs.indir, fn) for fn in infiles]

    # Determine the range of files to process
    if not cmdargs.stopfilenum:
        cmdargs.stopfilenum = len(infiles)
    cmdargs.stopfilenum = min(cmdargs.stopfilenum, len(infiles))

    infiles = infiles[cmdargs.startfilenum : cmdargs.stopfilenum]

    # Set constants
    nullVal = -999.0
    img_size_metres = int(cmdargs.tile_s)
    img_size_pixels = int(img_size_metres / cmdargs.psize)
    realD_thres = -10.0  # Minimum valid DEM value

    demProduct = "dem"
    demfileList = getDemImageFiles(infilesFull, demProduct)
    codesProduct = "NonGrdCodes"
    codesStage = filenaming_methods.stageByProductName[codesProduct]

    # Process each file
    for demfile in demfileList:
        img = rw_image_methods.imgRead(demfile)

        # Check if the DEM contains values below the threshold (excluding edges)
        if np.min(img[2:-3, 2:-3]) < realD_thres:
            # Define offsets for neighbouring tiles
            yidx = [img_size_metres, 0, -img_size_metres]
            xidx = [-img_size_metres, 0, img_size_metres]
            multi = len(yidx)

            # Initialise temporary arrays for DEM and codes
            temp_dem = np.zeros((img_size_pixels * multi, img_size_pixels * multi))
            temp_bd4 = np.zeros(
                (img_size_pixels * multi, img_size_pixels * multi), dtype=np.uint8
            )
            sImgCount = 0

            # Process neighbouring tiles
            for ct1, yoffset in enumerate(yidx):
                for ct2, xoffset in enumerate(xidx):
                    offsetDemFile = neighbourTileFilename(demfile, xoffset, yoffset)
                    offsetCodesFile = qvf.setstagecode(offsetDemFile, codesStage)
                    offsetCodesFile = qvf.setoptionfield(offsetCodesFile, 'l', codesProduct)

                    if Path.is_file(Path(offsetDemFile)):
                        img_bb0 = rw_image_methods.imgRead(offsetDemFile)
                        img_bb4 = rw_image_methods.imgRead(offsetCodesFile)

                        ys = int(ct1 * img_size_pixels)
                        xs = int(ct2 * img_size_pixels)
                        temp_dem[
                            ys : ys + img_size_pixels, xs : xs + img_size_pixels
                        ] = img_bb0
                        temp_bd4[
                            ys : ys + img_size_pixels, xs : xs + img_size_pixels
                        ] = img_bb4
                        sImgCount += 1

            # Perform DEM infill if enough neighbouring tiles are available
            if sImgCount > 5:
                res = gridding_methods.dem_infill(
                    temp_dem, temp_bd4, minElev=realD_thres
                )
                res = res[
                    img_size_pixels : 2 * img_size_pixels,
                    img_size_pixels : 2 * img_size_pixels,
                ]

                # Replace NaN values with the null value
                res[np.isnan(res)] = nullVal
            else:
                res = np.full((img_size_pixels, img_size_pixels), nullVal)

            # Save the infilled DEM if it contains valid data
            if np.max(res) > realD_thres:
                outDemfile = qvf.setoptionfield(demfile, 'l', "demInfilled")
                hillshadeProduct = "demHS"
                hillshadeStage = filenaming_methods.stageByProductName[hillshadeProduct]
                outDemHSfile = qvf.setstagecode(demfile, hillshadeStage)
                outDemHSfile = qvf.setoptionfield(outDemHSfile, 'l', hillshadeProduct)
                h = rw_image_methods.imgH(demfile)
                res[res < realD_thres] = nullVal

                rw_image_methods.writeImage(
                    np.round(res.astype(np.float32), 3),
                    outDemfile, cmdargs, tlx=h["tlx"], tly=h["tly"], binsize=h["pixel_s"],
                    epsg=cmdargs.epsg, nullVal=nullVal)

                args = ["gdaldem", "hillshade", outDemfile, outDemHSfile, "-compute_edges"]
                proc = subprocess.Popen(
                    args, stdout=subprocess.PIPE, stderr=subprocess.PIPE
                )
                stdout, stderr = proc.communicate()
                if proc.returncode != 0:
                    msg = f"{stderr.strip().decode('utf-8')}. Code: {proc.returncode}"
                    raise ValueError(msg)
                else:
                    print(stdout.decode("utf-8"))
                    print(stderr.decode("utf-8"))


def getDemImageFiles(infilesFull, productName):
    """
    Use the given list of laz files to deduce the list of corresponding DEM image
    files.

    Parameters:
      infilesFull (list[str]): List of full paths for indexed laz files
      demStage (str): QVF stage code for dem files

    Returns:
      demfileList (list[str]): List of DEM image files corresponding to the
                               given list of laz files
    """
    demStage = filenaming_methods.stageByProductName[productName]
    demfileList = []
    for lazfile in infilesFull:
        subdir = qvf.setsuffix(lazfile, '')
        demfilePattern = qvf.setstagecode(os.path.basename(lazfile), demStage)
        demfilePattern = qvf.setoptionfield(demfilePattern, 'l', productName)
        demfilePattern = qvf.setoptionfield(demfilePattern, 'p', qvf.getoptionfield(lazfile, 'p'))
        demfilePattern = qvf.setoptionfield(demfilePattern, 'r', '*')
        demfilePattern = qvf.setsuffix(demfilePattern, '*')
        demfilePattern = os.path.join(subdir, demfilePattern)
        candidateList = glob.glob(demfilePattern)
        candidateList = [fn for fn in candidateList if gdal.IdentifyDriver(fn) is not None]
        if len(candidateList) == 0:
            print(demfilePattern)
            msg = f"No matching DEM image for laz file {lazfile}"
            raise FileNotFoundError(msg)
        elif len(candidateList) > 1:
            msg = f"Unable to identify matching DEM, found {candidateList}"
            raise ValueError(msg)

        demfile = candidateList[0]
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
    xNdx = where.find('x')
    yNdx = where.find('y')
    zNdx = where.find('z')
    xCoord = int(where[xNdx + 1:yNdx])
    yCoord = int(where[yNdx + 2:zNdx])
    newX = xCoord + int(xOffset)
    newY = yCoord + int(yOffset)
    utmZone = int(where[zNdx + 1:])
    if where[yNdx + 1] == 's':
        utmZone = -utmZone
    newWhere = qvf.makeTileWhere(newX, newY, utmZone)
    nbrfileFull = qvf.setwhere(tilefile, newWhere)

    # We probably also need to change the where field in the directory name
    (nbrdir, nbrfile) = os.path.split(nbrfileFull)
    if where == qvf.getwhere(nbrdir):
        nbrdir = qvf.setwhere(nbrdir, newWhere)
        nbrfileFull = os.path.join(nbrdir, nbrfile)

    return nbrfileFull


def tif2cog(infile, outfile):
    """
    Convert a GeoTIFF file to a Cloud Optimized GeoTIFF (COG) format.

    Parameters:
        infile (str): Path to the input GeoTIFF file.
        outfile (str): Path to the output COG file.

    Returns:
        str: Confirmation message indicating the COG was created.

    Raises:
        Exception: If the gdal_translate command fails.
    """
    args = [
        "gdal_translate",
        infile,
        outfile,
        "-of",
        "COG",
        "-co",
        "BLOCKSIZE=256",
        "-co",
        "RESAMPLING=BILINEAR",
        "-co",
        "COMPRESS=DEFLATE",
    ]

    proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()

    if proc.returncode != 0:
        msg = f"Error creating COG: {stderr.strip().decode('utf-8')}. Exit code: {proc.returncode}"
        raise ValueError(msg)
    else:
        print(stdout.decode("utf-8"))
        print(stderr.decode("utf-8"))

    return "COG created successfully."


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

    # Run the DEM infill process
    run_dem_correction(cmdargs)


if __name__ == "__main__":
    main()
