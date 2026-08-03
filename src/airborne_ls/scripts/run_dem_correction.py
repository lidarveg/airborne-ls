#!/usr/bin/env python

"""
This script performs DEM infill for missing areas in LAS files and converts the results to Cloud Optimized GeoTIFF (COG) format.

Needs refinement / further development

"""

import argparse
import logging
import subprocess
import sys
from pathlib import Path

import numpy as np

from airborne_ls import gridding_methods, rw_image_methods

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

    parser.add_argument(
        "--indir", required=True, help="Directory containing LAS files."
    )
    parser.add_argument(
        "--laz_flist", required=True, help="List of LAS files to be processed."
    )
    parser.add_argument(
        "--epsg", type=int, required=True, help="EPSG code for map information."
    )
    parser.add_argument(
        "--psize",
        default=0.5,
        type=float,
        help="Pixel size of gridded DEM, Intensity, and maxH output layers (metres). Default: %(default)s.",
    )
    parser.add_argument(
        "--tile_s",
        type=float,
        required=True,
        help="XY dimensions of LAS tile in metres.",
    )
    parser.add_argument(
        "--startfilenum",
        default=0,
        type=int,
        help="Position within laz_flist to start batch processing. Default: %(default)s.",
    )
    parser.add_argument(
        "--stopfilenum",
        type=int,
        help="Position within laz_flist to stop batch processing.",
    )

    return parser.parse_args(inputargs)


def run_dem_correction(cmdargs):
    """
    Perform DEM infill for missing areas in LAS files.

    Parameters:
        cmdargs (argparse.Namespace): Parsed command-line arguments containing input directory, file list, and processing parameters.
    """
    # Read input file list
    fn = Path(cmdargs.indir).joinpath(cmdargs.laz_flist)
    with open(fn) as f:
        infiles = (line.strip() for line in f)

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
    resolution = f"r{int(cmdargs.psize * 100)}cm"

    # Process each file
    for fn in infiles:
        infileFull = Path(cmdargs.indir).joinpath(fn)
        base = fn.replace(".laz", "")
        demf = f"{cmdargs.indir}/{base}/{base}_bb0_dem_{resolution}.tif"
        img = rw_image_methods.imgRead(demf)

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
            where = (fn.split("_"))[1]
            sImgCount = 0

            # Process neighbouring tiles
            for ct1, yoffset in enumerate(yidx):
                for ct2, xoffset in enumerate(xidx):
                    xloc = int(where[1:7]) + xoffset
                    yloc = int(where[9:]) + yoffset
                    whereOffset = f"x{xloc}ys{yloc}"
                    baseOffset = base.replace(where, whereOffset)
                    offsetF = f"{cmdargs.indir}/{baseOffset}/{baseOffset}_bb0_dem_{resolution}.tif"

                    if Path.is_file(Path(offsetF)):
                        img_bb0 = rw_image_methods.imgRead(offsetF)
                        img_bb4 = rw_image_methods.imgRead(
                            offsetF.replace("dem", "NonGrd_codes").replace("bb0", "bb4")
                        )
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
                outf_fn = demf.replace("dem", "dem_infilled")
                outf_cog = outf_fn
                outfile = str(outf_fn).replace(Path(outf_fn).suffix, "_temp.tif")
                h = rw_image_methods.imgH(demf)
                res[res < realD_thres] = nullVal

                rw_image_methods.writeImage(
                    np.round(res.astype(np.float32), 3),
                    outfile,
                    cmdargs,
                    tlx=h["tlx"],
                    tly=h["tly"],
                    binsize=h["pixel_s"],
                    epsg=cmdargs.epsg,
                    nullVal=nullVal,
                    parent_file=infileFull,
                )
                tif2cog(outfile, outf_cog)

                # Recreate hillshade
                demHS_cog = f"{cmdargs.indir}/{base}/{base}_bbi_demHS_{resolution}.tif"
                outfile_HS = str(demHS_cog).replace(Path(demHS_cog).suffix, "_temp.tif")

                args = ["gdaldem", "hillshade", outfile, outfile_HS, "-compute_edges"]
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

                tif2cog(outfile_HS, demHS_cog)

                # Clean up temporary files
                Path(outfile).unlink()
                Path(outfile_HS).unlink()


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
