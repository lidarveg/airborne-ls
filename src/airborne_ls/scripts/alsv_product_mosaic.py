#!/usr/bin/env python

"""
Purpose: Generate the product mosaics from individually processed LiDAR tiles.

example: alsv_product_mosaic \
            --indir lidarveg_testing_data/Brisbane_2014_LGA_sub/indexed_tiles/ \
            --pixsize 0.5 --chm_pixsize 0.2

"""

import sys
import os
import argparse
import glob
import logging

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
    parser.add_argument("--outdir",
        help="Directory to write mosaics. Default is same as indir")
    parser.add_argument("--skipexisting", default=False, action="store_true",
        help="Skip any output files which already exist (default will over-write)")
    pixsizeGrp = parser.add_argument_group("Pixel sizes, used for names of input " +
        "and output files")
    pixsizeGrp.add_argument("--pixsize", default=0.5, type=float,
        help=("Pixel size of gridded DEM, Intensity, and maxH output layers (metres) " +
              "(default=%(default)s)"))
    pixsizeGrp.add_argument("--pcntile_pixsize", default=5.0, type=float,
        help="Pixel size of percentile output layers (metres) (default=%(default)s)")
    pixsizeGrp.add_argument("--fpc_pixsize", default=10.0, type=float,
        help="Pixel size of the FPC layer (metres) (default=%(default)s)")
    pixsizeGrp.add_argument("--chm_pixsize", default=0.2, type=float,
        help="Pixel size of the CHM layer (metres) (default=%(default)s)")
    parser.add_argument("--stagecode",
        help="Three-letter stage code to run. Default will run all stage codes.")
    parser.add_argument("--driver", default='GTiff',
        help=("GDAL driver for image format (default=%(default)s). Used to identify " +
              "suffix for input tile files, and for the output mosaic files. " +
              "If 'GTiff', then it will use the COG variant for output, " +
              "otherwise used as given"))

    cmdargs = parser.parse_args(inputargs)

    # Input checks
    if cmdargs.outdir is None:
        cmdargs.outdir = cmdargs.indir

    if not os.path.isdir(cmdargs.indir):
        logger.error(f"Input directory {cmdargs.indir} not found.")
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
    # Remove product names which are not supposed to be as tiles in the first place.
    productNameList = [prodName for prodName in productNameList
                       if prodName not in filenaming_methods.nonTileProducts]

    # Start with the full list of LAZ files
    infiles = sorted(glob.glob(f"{cmdargs.indir}/*.laz"))

    for productName in productNameList:
        missing_tiles = []
        tileFileList = []
        for lazfile in infiles:
            subdir = lazfile.replace('.laz', '')
            stage = filenaming_methods.stageByProductName[productName]
            projectName = qvf.getoptionfield(lazfile, 'p')
            productFile = os.path.basename(lazfile)
            productFile = qvf.setoptionfield(productFile, 'l', productName)
            productFile = qvf.setstagecode(productFile, stage)
            suffix = filenaming_methods.getSuffixFromDriverName(cmdargs.driver)
            productFile = qvf.setsuffix(productFile, suffix)
            productFile = os.path.join(subdir, productFile)

            # Choose resolution based on productName
            res = cmdargs.pixsize
            if productName == "fpc":
                res = cmdargs.fpc_pixsize
            elif productName == "chm":
                res = cmdargs.chm_pixsize
            elif productName.startswith('percentile'):
                res = cmdargs.pcntile_pixsize
            resStr = filenaming_methods.resolutionStrFromMetres(res)
            productFile = qvf.setoptionfield(productFile, 'r', resStr)

            # For the "dem" product, check for infilled version first
            if productName == "dem":
                productFile = qvf.setoptionfield(productFile, 'l', "gapfilleddem")
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
            sys.exit(1)

        # Create the mosaic output file, starting with the first input file name
        outFile = os.path.basename(tileFileList[0])
        outFile = qvf.setoptionfield(outFile, 'l', None)
        outFile = qvf.setoptionfield(outFile, 'p', None)
        outFile = qvf.setwhere(outFile, f"r{projectName}")
        outFile = os.path.join(cmdargs.outdir, outFile)

        if not (cmdargs.skipexisting and os.path.exists(outFile)):
            vrtFilename = qvf.setsuffix(outFile, 'vrt')
            gdal.BuildVRT(vrtFilename, tileFileList)

            logger.info(f"Output mosaic name: {outFile}")
            logger.debug(f"Missing tiles in {outFile} include {missing_tiles}")

            driverName = cmdargs.driver
            if driverName == "GTiff":
                driverName = "COG"
            creationOptions = rw_image_methods.creationOptionsByDriver.get(driverName, [])
            if qvf.getstagecode(outFile) in ('bb3', 'bb4'):
                # We do NOT want BILINEAR overview resampling for these stages
                BILINEAR_RESAMPLING = "RESAMPLING=BILINEAR"
                if BILINEAR_RESAMPLING in creationOptions:
                    creationOptions = [co for co in creationOptions if co != BILINEAR_RESAMPLING]
                    creationOptions.append('RESAMPLING=MODE')
            translateOptions = gdal.TranslateOptions(format=driverName,
                creationOptions=creationOptions)
            gdal.Translate(outFile, vrtFilename, options=translateOptions)

            os.remove(vrtFilename)
        else:
            print("Skipping", outFile)


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
