#!/usr/bin/env python

"""

Purpose: run importation of supplied las/laz files for lidarveg

bin_size: in my experience 100m works well for lower density datasets (<10 pts/ mt2)
          and 50m works well for lower density datasets (> 10 pts/ mt2).
          Not sure whether a finer bin size will be required (e.g. 25m)?
challenge: most project information (i.e. header info) cannot be relied upon to after
           this processing stage
options: educated guess which one to use by user (current approach) but a better approach
         could be to look at file size and
point format from header (this should be correct) using the largest few files?



example:
    uv run python scripts/run_lidar_standardisation.py\
        --indir lidarveg_testing_data/Brisbane_2014_LGA_sub/ \
        --outdr lidarveg_testing_data/Brisbane_2014_LGA_sub/indexed_tiles/\
        --epsg 28356 --laz_flist laz_flist \
        --tile_s 1000. --out_tile_s 1000. --ii 'mp' --proj brisba --year 2014 --binSize 50.

"""

import argparse

# Configure logging
import logging
import sys
from pathlib import Path

import laspy
import numpy as np
from osgeo import osr

from airborne_ls import lazfile_rw

logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.ERROR, format="%(asctime)s: %(name)20s: %(levelname)10s: %(message)s"
)
osr.UseExceptions()


def getCmdargs(inputargs):
    """
    Parse command-line arguments for running as a standalone script.

    Parameters:
        inputargs (list): List of command-line arguments.

    Returns:
        argparse.Namespace: Parsed arguments as a namespace object.
    """
    parser = argparse.ArgumentParser(description="Standardise LAS/LAZ files for lidarveg.")

    # Input and output directories
    parser.add_argument("--indir", required=True,
        help="Full path to directory containing LAS/LAZ files.")
    parser.add_argument("--outdr", required=True,
        help="Directory for newly named and indexed LAS/LAZ files.")
    parser.add_argument("--laz_flist", required=True,
        help="Text file containing LAS/LAZ files to be processed; one file per row.")

    # Tile dimensions
    parser.add_argument("--tile_s", required=True, type=float,
        help="Maximum XY dimension of LAS/LAZ file (metres).")
    parser.add_argument("--out_tile_s", default=1000, type=float,
        help="Equal or smaller maximum XY dimension for output LAS/LAZ files (metres).")

    # EPSG and spatial database options
    parser.add_argument("--epsg", type=int, required=True,
        help="EPSG code for map information.")

    # Metadata options
    parser.add_argument("--ss", type=str, default="ap",
        help="Platform type (e.g., 'ap' for airborne platform).")
    parser.add_argument("--ii", type=str, required=True,
        help="Predefined sensor code; can use uk if unknown")
    parser.add_argument("--pp", type=str, default="dr",
        help="Product type (e.g., 'dr' for discrete return).")
    parser.add_argument("--proj", required=True,
        help="Six-character project name (e.g., 'brisba').")
    parser.add_argument("--year", type=int, required=True,
        help="Year of data capture (e.g., 2022).")

    # Tile indexing options
    parser.add_argument("--binSize", default=50.0, type=float,
        help="XY bin size for data indexing (metres).")
    parser.add_argument("--startfilenum", default=0, type=int,
        help="Start position in file list for batch processing.")
    parser.add_argument("--stopfilenum", type=int,
        help="Stop position in file list for batch processing.")

    # Metadata flags

    cmdargs = parser.parse_args(inputargs)

    # Validate input file list
    infilelist = Path(cmdargs.indir).joinpath(cmdargs.laz_flist)
    if not infilelist.is_file():
        laslist = list(Path(cmdargs.indir).glob("*.las"))
        lazlist = list(Path(cmdargs.indir).glob("*.laz"))

        if laslist:
            laslist = reorder_flist(laslist)
            cmdargs.laz_flist = "lazlist_auto"
            outfile = Path(cmdargs.indir).joinpath(cmdargs.laz_flist)
            with open(outfile, "w") as fout:
                fout.writelines(f"{Path(fn).name}\n" for fn in laslist)

        if lazlist:
            lazlist = reorder_flist(lazlist)
            cmdargs.laz_flist = "lazlist_auto"
            outfile = Path(cmdargs.indir).joinpath(cmdargs.laz_flist)
            print(f"outfile: {outfile}")
            with open(outfile, "w") as fout:
                fout.writelines(f"{Path(fn).name}\n" for fn in lazlist)

    infilelist = Path(cmdargs.indir).joinpath(cmdargs.laz_flist)
    if not infilelist.is_file():
        raise AssertionError(f" laz_flist is invalid: {infilelist} ")

    # Validate project name
    if len(cmdargs.proj) != 6:
        raise AssertionError("Project name must be exactly 6 characters.")

    # Validate output directory
    if not Path(cmdargs.outdr).exists():
        raise AssertionError("Output directory does not exist.")

    if (cmdargs.tile_s % cmdargs.out_tile_s) != 0:
        msg = (f"Input tile size {cmdargs.tile_s} not divisible by " +
               f"output tile size {cmdargs.out_tile_s}")
        raise ValueError(msg)

    if (cmdargs.out_tile_s % cmdargs.binSize) != 0:
        msg = f"Output tile size {cmdargs.out_tile_s} not divisible by bin size {cmdargs.binSize}"
        raise ValueError(msg)

    return cmdargs


def run_las_standardisation(cmdargs):
    """
    Prepare a new ALS project by processing LAS/LAZ files and generating indexed tiles.

    Parameters:
        cmdargs (argparse.Namespace): Parsed command-line arguments.
    """
    # Check input files
    lazlistfull, _ = check_input_fns(
        cmdargs.indir, Path(cmdargs.indir).joinpath(cmdargs.laz_flist)
    )

    # Determine the range of files to process
    if not cmdargs.stopfilenum:
        cmdargs.stopfilenum = len(lazlistfull)
    cmdargs.stopfilenum = min(cmdargs.stopfilenum, len(lazlistfull))
    lazlistfull = lazlistfull[cmdargs.startfilenum : cmdargs.stopfilenum]

    input_tileS = cmdargs.tile_s

    what = f"{cmdargs.ss}{cmdargs.ii}{cmdargs.pp}"
    when = f"{cmdargs.year}"
    srs = osr.SpatialReference(epsg=cmdargs.epsg)
    utmZone = srs.GetUTMZone()
    if utmZone == 0:
        msg = f"Unknown EPSG {cmdargs.epsg}. Cannot translate to filename zone code"
        raise ValueError(msg)
    stageCode = getStageCode(cmdargs)

    for inLazfile in lazlistfull:
        print(f"Processing LAS/LAZ file: {inLazfile}")
        data = laspy.read(inLazfile)

        # Calculate northing and easting
        northing = int(
            np.ceil(np.median(data.y[data.y > 0]) / input_tileS) * input_tileS)
        easting = int(
            np.floor(np.median(data.x[data.x > 0]) / input_tileS) * input_tileS)

        # Run chunked LAS filtering
        _ = lazfile_rw.standardise_lasf(what, when, utmZone, stageCode,
                cmdargs.proj, cmdargs.outdr, data, easting, northing, input_tileS,
                cmdargs.out_tile_s, cmdargs.binSize, inLazfile)

        del data


def reorder_flist(lazlistfull):
    """
    Reorder a list of LAS/LAZ files by file size in descending order.

    Parameters:
        lazlistfull (list): List of LAS/LAZ file paths.

    Returns:
        np.ndarray: Array of file paths sorted by file size in descending order.
    """
    file_sizes = [Path(ff).stat().st_size for ff in lazlistfull]
    sorted_indices = np.argsort(-np.array(file_sizes))  # Sort in descending order
    return np.array(lazlistfull)[sorted_indices]


def check_input_fns(indir, infilelist):
    """
    Check whether input files are valid and are .laz or .las files.

    Parameters:
        indir (str): Input directory containing LAS/LAZ files.
        infilelist (Path): Path to the file containing the list of LAS/LAZ files.

    Returns:
        tuple: A tuple containing the list of valid files and a boolean indicating if
               all files are valid.
    """
    with open(infilelist) as f:
        filelist = [line.strip() for line in f]

    flist = []
    laz_count = 0

    for fn in filelist:
        if fn.endswith((".laz", ".las")):
            laz_count += 1

        if not Path(fn).is_file():
            fn = Path(indir, fn)
            if not fn.is_file():
                logger.info(f"File does not exist: {fn}")
                raise FileNotFoundError(f"ALS file {fn} is invalid.")
            else:
                flist.append(fn)
        else:
            flist.append(fn)

    return flist, len(flist) == laz_count


def getStageCode(cmdargs):
    """
    Get the output stage code from the given command line arguments
    """
    # At the moment, it is always ba1, but in principle this depends on the
    # height datum in use, which will eventually be given on the command line
    # Stages defined as:  ??? Is this still correct ?????
    #     - ba0: Ellipsoid heights
    #     - ba1: Geoid heights (AHD)
    #     - ba3: New indexed LAS
    return "ba1"


def main(args=None):
    """
    Main entry point for the script, allowing external calls.

    Parameters:
        args (list[str], optional): Command-line parameter list. Defaults to None.
    """
    if args is None:
        args = sys.argv[1:]

    # Log the start of the process
    logger.debug("Parsing command-line arguments.")

    # Parse command-line arguments
    cmdargs = getCmdargs(args)

    # Log the parsed arguments
    logger.debug(f"Command-line arguments: {args}")
    logger.debug(f"Input directory: {cmdargs.indir}")

    # Start the ALS project preparation
    run_las_standardisation(cmdargs)


if __name__ == "__main__":
    main()
