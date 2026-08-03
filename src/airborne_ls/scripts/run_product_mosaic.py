#!/usr/bin/env python

"""
Purpose: Generate the product mosaics from individually processed LiDAR tiles.
"""

import argparse
import logging
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
from rios import rat

from airborne_ls import filenaming_methods

# Configure logging
logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.ERROR, format="%(asctime)s: %(name)20s: %(levelname)10s: %(message)s"
)

# Constants
PERCENTILES = [1, 5, 25, 50, 75, 95, 99]
PRODUCT_DICT = {
    "bb0": "dem",
    "bb1": "maxH",
    "bb2": "intens",
    "bb3": "grdR",
    "bb4": "NonGrd_returns",
    "bb5": "fst_dens",
    "bb8": "1_percentile",
    "bb9": "5_percentile",
    "bba": "25_percentile",
    "bbb": "50_percentile",
    "bbc": "75_percentile",
    "bbd": "95_percentile",
    "bbe": "99_percentile",
    "bbh": "fpc",
    "bbi": "demHS",
    "bbm": "csm",
    "bbn": "chm",
}


def getCmdargs(inputargs):
    """
    Parse command-line arguments.

    Parameters:
        inputargs (list): List of command-line arguments.

    Returns:
        argparse.Namespace: Parsed arguments as a namespace object.
    """
    parser = argparse.ArgumentParser(description="Generate mosaics from LiDAR tiles.")

    parser.add_argument(
        "--indir", required=True, help="Top-level directory containing input tiles."
    )
    parser.add_argument(
        "--outdr", help="Directory to write mosaics. Default is the current directory."
    )
    parser.add_argument(
        "--laz_flist", required=True, help="File containing the list of LAS/LAZ files."
    )
    parser.add_argument(
        "--tile_s",
        type=float,
        required=True,
        help="XY dimensions of LAS tile (metres).",
    )
    parser.add_argument(
        "--psize",
        default=0.5,
        type=float,
        help="Pixel size of gridded DEM, Intensity, and maxH output layers (metres). Default: %(default)s.",
    )
    parser.add_argument(
        "--ptile_s",
        default=5,
        type=float,
        help="Pixel size of percentile output layers (metres). Default: %(default)s.",
    )
    parser.add_argument(
        "--fpc_psize",
        default=10.0,
        type=float,
        help="Pixel size of the FPC layer (metres). Default: %(default)s.",
    )
    parser.add_argument(
        "--chm_psize",
        default=None,
        type=float,
        help="Optional: value estimated using pulse density (metres).",
    )
    parser.add_argument(
        "--outStageList", help="Three-letter stage code. If blank, run all stage codes."
    )

    cmdargs = parser.parse_args(inputargs)

    # Input checks
    if not cmdargs.outdr:
        cmdargs.outdr = Path(cmdargs.indir)

    if not cmdargs.outStageList:
        cmdargs.outStageList = list(PRODUCT_DICT.keys())
    else:
        cmdargs.outStageList = [cmdargs.outStageList]

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

    # read in list of las files and extract batch subset
    with open(cmdargs.laz_flist) as f:
        infiles = [line.strip() for line in f]

    fn_dict = filenaming_methods.createTileDict(infiles[0], cmdargs.tile_s)
    project, year, zone, zone_prefix, sensor_code = (
        fn_dict["project"],
        fn_dict["date"],
        fn_dict["zoneCode"],
        fn_dict["zone_prefix"],
        fn_dict["instrument"],
    )
    psizes = filenaming_methods.get_psizeDict(
        psize=cmdargs.psize,
        ptile_s=cmdargs.ptile_s,
        fpc_psize=cmdargs.fpc_psize,
        chm_psize=cmdargs.chm_psize,
    )
    missing_tiles = []
    temp_output = Path(cmdargs.indir).joinpath("temp_output")
    os.makedirs(temp_output, exist_ok=True)
    for loc, outStage in enumerate(cmdargs.outStageList):
        input_layers = []
        stagec_def = filenaming_methods.get_stageDict()
        tempList = Path(temp_output) / f"temp_list{outStage}"

        with open(tempList, "w") as fout:
            for fn in infiles:
                infileFull = os.path.join(cmdargs.indir, fn)
                outputDir = Path((infileFull).split(".")[0])
                tileBasename = (fn).split(".")[0]
                outputBasename = os.path.join(outputDir, tileBasename)
                fnames = filenaming_methods.get_outfnames(
                    outputBasename,
                    psize=cmdargs.psize,
                    ptile_s=cmdargs.ptile_s,
                    fpc_psize=cmdargs.fpc_psize,
                    chm_psize=cmdargs.chm_psize,
                    pptiles=PERCENTILES,
                )

                if (stagec_def[outStage]).endswith("percentile"):
                    num = int(((stagec_def[outStage]).split("_"))[0])
                    layer = fnames["ptiles"][
                        np.argwhere(np.array(PERCENTILES) == num)[0][0]
                    ]
                else:
                    layer = fnames[stagec_def[outStage]]

                if Path(layer).is_file():
                    ##### testing
                    if outStage == "bb0":
                        infilled_demf = layer.replace("dem", "dem_infilled")
                        if Path(infilled_demf).is_file():
                            layer = infilled_demf
                    #####
                    input_layers.append(layer)
                    fout.write(f"{layer}\n")
                    logger.info(f"File added to mosaic: {layer}")
                else:
                    missing_tiles.append(layer)
        fout.close()

        if len(str(year)) > 4:
            makeMosaicFilename = f"ap{sensor_code}dr_r{project}_y{year}_{outStage}{zone_prefix}{zone}_{psizes[outStage]}.tif"
        else:
            makeMosaicFilename = f"ap{sensor_code}dr_r{project}_{year}_{outStage}{zone_prefix}{zone}_{psizes[outStage]}.tif"

        logger.info(f"Output mosaic name: {makeMosaicFilename}")
        logger.debug(f"Missing tiles in {makeMosaicFilename} include {missing_tiles}")

        makeMosaicPathname = Path(os.path.join(cmdargs.outdr, makeMosaicFilename))
        vrt_filename = makeMosaicPathname.with_suffix(".vrt")

        args = ["gdalbuildvrt", "-input_file_list", tempList, vrt_filename]
        proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = proc.communicate()
        if proc.returncode != 0:
            # an error happened!
            err_msg = f"{stderr.strip()}. Code: {proc.returncode}"
            raise ValueError(err_msg)
        else:
            print(stdout)
            print(stderr)

        ########## note -- to create COG with history required two files
        outf_cog = makeMosaicPathname  # correct name will be used for COG
        outf_tif = str(makeMosaicPathname).replace(
            Path(makeMosaicPathname).suffix, "_temp.tif"
        )
        #########

        command = [
            "gdal_translate",
            "-of",
            "GTiff",
            "-co",
            "COMPRESS=LZW",
            "-co",
            "BIGTIFF=YES",
            "-co",
            "NUM_THREADS=4",
            vrt_filename,
            str(outf_tif),
        ]

        result = subprocess.run(command, capture_output=True, text=True, check=False)
        if result.returncode != 0:
            msg = "failed"
            raise ValueError(msg)

        if outStage == "bbh":
            applyFPCcolor(str(outf_tif))

        # Move vrt_filename (the VRT file) to temp folder
        vrt_filename_basename = vrt_filename.name  # Get the filename
        vrt_temp_path = (
            temp_output / vrt_filename_basename
        )  # Construct the destination path
        vrt_filename.rename(vrt_temp_path)  # Move the file
        vrt_filename = vrt_temp_path

        ## convert to COG
        args = [
            "gdal_translate",
            outf_tif,
            outf_cog,
            "-of",
            "COG",
            "-co",
            "BLOCKSIZE=256",
            "-co",
            "RESAMPLING=BILINEAR",
            "-co",
            "COMPRESS=DEFLATE",
            "-co",
            "BIGTIFF=YES",
        ]
        proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = proc.communicate()
        if proc.returncode != 0:
            # an error happened!
            err_msg = f"{stderr.strip()}. Code: {proc.returncode}"
            raise ValueError(err_msg)
        else:
            print(stdout)
            print(stderr)
        Path(outf_tif).unlink()


def applyFPCcolor(fout):
    """
    Apply a colour map to the FPC (Foliage Profile Curve) output.

    Parameters:
        fout (str): Path to the FPC output file.
    """
    clrTbl = np.zeros((256, 4), dtype=np.uint8)
    clrTbl.fill(255)
    clrTbl[10:90, 0] = np.mgrid[255:0:-80j].round().astype(np.uint8)
    clrTbl[90:101, 0] = 0
    clrTbl[:, 2] = clrTbl[:, 0]
    clrTbl[10:90, 1] = np.mgrid[255:100:-80j].round().astype(np.uint8)
    clrTbl[90:101, 1] = 100
    clrTbl[0, :] = [210, 180, 140, 255]  # Add brown for zero FPC

    rat.setColorTable(fout, clrTbl)


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
