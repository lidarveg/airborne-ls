#!/usr/bin/env python

"""
Purpose: Generate the product mosaics from individually processed LiDAR tiles.

example: uv run python scripts/run_product_mosaic.py \
    --indir lidarveg_testing_data/Brisbane_2014_LGA_sub/indexed_tiles/ \
    --laz_flist  laz_flist --tile_s 1000. --psize 0.5 --chm_psize 0.2

"""

import argparse
import logging
import os
import sys
from pathlib import Path

from osgeo import gdal

from airborne_ls import filenaming_methods, qvf, rw_image_methods


gdal.UseExceptions()

# Configure logging
logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.ERROR, format="%(asctime)s: %(name)20s: %(levelname)10s: %(message)s"
)


def getCmdargs(inputargs):
    """
    Parse command-line arguments.

    Parameters:
        inputargs (list): List of command-line arguments.

    Returns:
        argparse.Namespace: Parsed arguments as a namespace object.
    """
    parser = argparse.ArgumentParser(description="Generate mosaics from LiDAR tiles.")

    parser.add_argument("--indir", required=True,
        help="Top-level directory containing input tiles.")
    parser.add_argument("--outdr",
        help="Directory to write mosaics. Default is the current directory.")
    parser.add_argument("--laz_flist", required=True,
        help="File containing the list of LAS/LAZ files.")
    parser.add_argument("--tile_s", type=float, required=True,
        help="XY dimensions of LAS tile (metres).")
    parser.add_argument("--psize", default=0.5, type=float,
        help=("Pixel size of gridded DEM, Intensity, and maxH output layers (metres). " +
              "Default: %(default)s."))
    parser.add_argument("--ptile_s", default=5, type=float,
        help="Pixel size of percentile output layers (metres). Default: %(default)s.")
    parser.add_argument("--fpc_psize", default=10.0, type=float,
        help="Pixel size of the FPC layer (metres). Default: %(default)s.")
    parser.add_argument("--chm_psize", default=None, type=float,
        help="Optional: value estimated using pulse density (metres).")
    parser.add_argument("--stagecode",
        help="Three-letter stage code to run. If blank, run all stage codes.")
    parser.add_argument("--driver", default='GTiff',
        help=("GDAL driver for image format (default=%(default)s). If 'GTiff', then " +
              "use the COG variant for output, otherwise use as given"))

    cmdargs = parser.parse_args(inputargs)

    # Input checks
    if not cmdargs.outdr:
        cmdargs.outdr = Path(cmdargs.indir)

    if not Path(cmdargs.indir).exists():
        logger.error(f"Input directory {cmdargs.indir} not found.")
        sys.exit(1)

    cmdargs.laz_flist = Path(cmdargs.indir).joinpath(cmdargs.laz_flist)
    if not cmdargs.laz_flist.is_file():
        logger.error(f"File list {cmdargs.laz_flist} not found.")
        sys.exit(1)

    return cmdargs


def runMerge(cmdargs):
    """
    Main routine
    """
    if cmdargs.stagecode is not None:
        outStageList = [cmdargs.stagecode]
    else:
        outStageList = list(filenaming_methods.stageByProductName.values())
    productNameList = [filenaming_methods.productNameByStage[stage] for stage in outStageList]

    # read in list of las files and extract batch subset
    infiles = [line.strip() for line in open(cmdargs.laz_flist)]

    for productName in productNameList:
        missing_tiles = []
        tileFileList = []
        for lazfile in infiles:
            subdir = os.path.join(cmdargs.indir, lazfile.replace('.laz', ''))
            stage = filenaming_methods.stageByProductName[productName]
            projectName = qvf.getoptionfield(lazfile, 'p')
            productFile = qvf.setoptionfield(lazfile, 'l', productName)
            productFile = qvf.setstagecode(productFile, stage)
            suffix = filenaming_methods.getSuffixFromDriverName(cmdargs.driver)
            productFile = qvf.setsuffix(productFile, suffix)
            productFile = os.path.join(subdir, productFile)

            # Choose resolution based on productName
            res = cmdargs.psize
            if productName == "fpc":
                res = cmdargs.fpc_psize
            elif productName == "chm":
                res = cmdargs.chm_psize
            elif productName.startswith('percentile'):
                res = cmdargs.ptile_s
            resStr = filenaming_methods.resolutionStrFromMetres(res)
            productFile = qvf.setoptionfield(productFile, 'r', resStr)

            # For the "dem" product, check for infilled version first
            if productName == "dem":
                productFile = qvf.setoptionfield(productFile, 'l', "demInfilled")
                if not os.path.exists(productFile):
                    productFile = qvf.setoptionfield(productFile, 'l', "dem")

            if not os.path.exists(productFile):
                missing_tiles.append(productFile)
            else:
                tileFileList.append(productFile)
        if len(missing_tiles) > 0:
            print("Missing", missing_tiles)
        if len(tileFileList) == 0:
            print("Missing everything")

        # Create the mosaic output file, starting with the first input file name
        outFile = os.path.basename(tileFileList[0])
        outFile = qvf.setoptionfield(outFile, 'l', None)
        outFile = qvf.setoptionfield(outFile, 'p', None)
        outFile = qvf.setwhere(outFile, f"r{projectName}")
        outFile = os.path.join(cmdargs.outdr, outFile)

        vrtFilename = qvf.setsuffix(outFile, 'vrt')
        gdal.BuildVRT(vrtFilename, tileFileList)

        logger.info(f"Output mosaic name: {outFile}")
        logger.debug(f"Missing tiles in {outFile} include {missing_tiles}")

        driverName = cmdargs.driver
        if driverName == "GTiff":
            driverName = "COG"
        creationOptions = rw_image_methods.creationOptionsByDriver.get(driverName, [])
        if qvf.getstagecode(outFile) in ('bb3', 'bb4'):
            # We do NOT want BILINEAR overview resampling for these two stages
            BILINEAR_RESAMPLING = "RESAMPLING=BILINEAR"
            if BILINEAR_RESAMPLING in creationOptions:
                creationOptions = [co for co in creationOptions if co != BILINEAR_RESAMPLING]
                creationOptions.append('RESAMPLING=MODE')
        translateOptions = gdal.TranslateOptions(format=driverName,
            creationOptions=creationOptions)
        gdal.Translate(outFile, vrtFilename, options=translateOptions)

        os.remove(vrtFilename)


def main(args=None):
    """
    Main entry point for the script.

    Parameters:
        args (list[str], optional): Command-line parameter list. Defaults to None.
    """
    if args is None:
        args = sys.argv[1:]
    cmdargs = getCmdargs(args)
    runMerge(cmdargs)


if __name__ == "__main__":
    main()
