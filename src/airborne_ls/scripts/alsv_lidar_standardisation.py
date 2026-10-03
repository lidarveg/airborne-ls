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
    alsv_lidar_standardisation \
        --indir lidarveg_testing_data/Brisbane_2014_LGA_sub/ \
        --outdir lidarveg_testing_data/Brisbane_2014_LGA_sub/indexed_tiles/\
        --epsg 28356 --intilesize 1000 --outtilesize 1000 --ii 'mp' \
        --project brisba --year 2014 --binsize 50.

"""

import argparse
import glob

# Configure logging
import logging
import sys
from pathlib import Path

import laspy
import numpy as np
from osgeo import osr
import pyproj

from airborne_ls import lazfile_rw, const, filenaming_methods

logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.ERROR, format="%(asctime)s: %(name)20s: %(levelname)10s: %(message)s"
)
osr.UseExceptions()


DFLT_CLASSESTOEXCLUDE = ",".join([
    str(const.PTCLASS_NOISE_LOWPOINT),
    str(const.PTCLASS_NOISE_HIGHPOINT),
    str(const.PTCLASS_NOISE_PROVIDERDEFINED)
])
# Default min/max acceptable values (metres) for point height. These are good for
# Australian continent, although the Antarctic Territories have some higher points
# (e.g. Mt McClintock, 3490m).
DFLT_MINZ = -20         # Lower than Lake Eyre/Kati Thanda
DFLT_MAXZ = 2300        # Higher than Mt Kosciuszko


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
    parser.add_argument("--indir", help=("Directory containing LAS/LAZ files " +
        "to process."))
    parser.add_argument("--infile", help="Name of a single LAS/LAZ file to process")
    parser.add_argument("--groupMofN", nargs=2, metavar=('M', 'N'), type=int,
        help=("Use only with --indir. Divide the set of tiles to be processed into N groups, " +
              "and then process only the files in the M-th group (group numbering starts at 1)." +
              "This helps support efficient batch processing, where the user decides how " +
              "many batch jobs will run, and the command for each group processes only " +
              "the tiles for that group. For example, with 5 groups, the first group " +
              "would be specified as '--groupMofN 1 5'"))
    parser.add_argument("--skipexisting", default=False, action="store_true",
        help=("Skip existing output LAZ files. Default will re-create any " +
              "output files which already exist"))
    parser.add_argument("--outdir", required=True,
        help="Directory for newly named and indexed LAS/LAZ files.")
    parser.add_argument("--minlasversion", default="1.4",
        help=("Minimum LAS format version for output files (default=%(default)s). " +
              "Newer format input files will preserve their newer version"))

    # Tile dimensions
    parser.add_argument("--intilesize", required=True, type=int,
        help="Maximum XY dimension of LAS/LAZ file (metres).")
    parser.add_argument("--outtilesize", type=int,
        help=("XY dimension for output LAZ files (metres). Must be <= --intilesize, and " +
              "--intilesize must be whole number multiple of --outtilesize. " +
              "Default is same as --intilesize"))

    # EPSG and spatial database options
    parser.add_argument("--epsg", type=int,
        help=("EPSG code for map projection of input files. Default will check in " +
              "input LAS files, but this over-rides"))

    # Metadata options
    parser.add_argument("--ss", type=str, default="ap",
        help="Platform type (e.g., 'ap' for airborne platform) (default=%(default)s)")
    parser.add_argument("--ii", type=str, required=True,
        help="Predefined sensor code; can use uk if unknown")
    parser.add_argument("--pp", type=str, default="dr",
        help="Product type (e.g., 'dr' for discrete return) (default=%(default)s)")
    parser.add_argument("--project", required=True,
        help="Six-character project name (e.g., 'brisba').")
    parser.add_argument("--year", type=int, required=True,
        help="Year of data capture (e.g., 2022).")
    parser.add_argument("--excludeclasses", default=DFLT_CLASSESTOEXCLUDE,
        help=("List of point class values to exclude from data. Comma-separated, " +
              "no spaces. (default=%(default)s)"))
    parser.add_argument("--minz", type=float, default=DFLT_MINZ,
        help="Minimum acceptable point height value (metres) (default=%(default)s)")
    parser.add_argument("--maxz", type=float, default=DFLT_MAXZ,
        help="Maximum acceptable point height value (metres) (default=%(default)s)")

    # Tile indexing options
    parser.add_argument("--binsize", default=50.0, type=float,
        help="XY bin size for data indexing (metres).")

    cmdargs = parser.parse_args(inputargs)

    if cmdargs.indir is not None and cmdargs.infile is not None:
        msg = "Use either --indir or --infile, but not both"
        raise ValueError(msg)

    if cmdargs.indir is None and cmdargs.infile is None:
        msg = "Must supply one of --indir or --infile"
        raise ValueError(msg)

    if cmdargs.indir is None and cmdargs.groupMofN is not None:
        msg = "Using --groupMofN requires --indir"
        raise ValueError(msg)

    # Validate project name
    if len(cmdargs.project) != 6:
        raise ValueError("Project name must be exactly 6 characters.")

    if cmdargs.indir == cmdargs.outdir:
        raise ValueError("--outdir cannot be the same as --indir")

    # Validate output directory
    if not Path(cmdargs.outdir).exists():
        raise ValueError(f"Output directory '{cmdargs.outdir}' does not exist")

    if cmdargs.outtilesize is None:
        cmdargs.outtilesize = cmdargs.intilesize

    if (cmdargs.intilesize % cmdargs.outtilesize) != 0:
        msg = (f"Input tile size {cmdargs.intilesize} not divisible by " +
               f"output tile size {cmdargs.outtilesize}")
        raise ValueError(msg)

    if (cmdargs.outtilesize % cmdargs.binsize) != 0:
        msg = f"Output tile size {cmdargs.outtilesize} not divisible by bin size {cmdargs.binsize}"
        raise ValueError(msg)

    lasVersNum = tuple([int(i) for i in cmdargs.minlasversion.split('.')])
    if len(lasVersNum) > 2:
        msg = f"LAS version {cmdargs.minlasversion} has too many components ({len(lasVersNum)})"
        raise ValueError(msg)
    cmdargs.minlasversion = laspy.header.Version(*lasVersNum)

    return cmdargs


def run_las_standardisation(cmdargs):
    """
    Prepare a new ALS project by processing LAS/LAZ files and generating indexed tiles.

    Parameters:
        cmdargs (argparse.Namespace): Parsed command-line arguments.
    """
    if cmdargs.infile is not None:
        infileList = [cmdargs.infile]
    elif cmdargs.indir is not None:
        pattern = f"{cmdargs.indir}/*.la[sz]"
        infileList = sorted(glob.glob(pattern))
        if cmdargs.groupMofN is not None:
            (group, numGroups) = tuple(cmdargs.groupMofN)
            infileList = filenaming_methods.filelistGroupSubset(infileList, group, numGroups)

    # Sort out input and output map projections
    if cmdargs.epsg is None:
        cmdargs.epsg = getEPSGfromLAS(infileList[0])
        if cmdargs.epsg is not None:
            print(f"Found EPSG {cmdargs.epsg} in input las file")

    if cmdargs.epsg is None:
        msg = "No EPSG found in first infile. Please supply --epsg"
        raise ValueError(msg)

    srs = osr.SpatialReference()
    srs.ImportFromEPSG(cmdargs.epsg)
    utmZone = srs.GetUTMZone()
    if utmZone == 0:
        msg = f"Unknown EPSG {cmdargs.epsg}. Cannot translate to filename zone code"
        raise ValueError(msg)

    inTilesize = cmdargs.intilesize
    what = f"{cmdargs.ss}{cmdargs.ii}{cmdargs.pp}"
    when = f"{cmdargs.year}"
    stageCode = getStageCode(cmdargs)
    classesToExclude = [int(i) for i in cmdargs.excludeclasses.split(',')]

    for inLazfile in infileList:
        print(f"Processing LAS/LAZ file: {inLazfile}")
        data = laspy.read(inLazfile)
        # Convert to older formats to the requested LAS version.
        if data.header.version < cmdargs.minlasversion:
            data = laspy.convert(data, file_version=cmdargs.minlasversion)
        if data.header.parse_crs() is None:
            crsObj = pyproj.CRS.from_epsg(cmdargs.epsg)
            data.header.add_crs(crsObj)

        checkElevationRange(data, cmdargs.minz, cmdargs.maxz)

        # Calculate northing and easting of top-left corner of input tile
        northing = int(
            np.ceil(np.median(data.y[data.y > 0]) / inTilesize) * inTilesize)
        easting = int(
            np.floor(np.median(data.x[data.x > 0]) / inTilesize) * inTilesize)

        # Run chunked LAS filtering
        lazfile_rw.standardise_lasf(what, when, utmZone, stageCode, cmdargs.project,
            cmdargs.outdir, data, easting, northing, inTilesize, cmdargs.outtilesize,
            cmdargs.binsize, inLazfile, classesToExclude, cmdargs.minz, cmdargs.maxz,
            cmdargs.skipexisting)

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


def checkElevationRange(data, minz, maxz):
    """
    Check the elevation data against the given acceptable range.

    Print warning messages if points are outside range.

    Parameters:
      data (laspy.lasdata.LasData): Data from LAS file
      minz, maxz (float): Acceptable range of elevation data
    """
    hdr = data.header
    zVals = data['Z'] * hdr.z_scale + hdr.z_offset
    if hdr.z_min < minz:
        belowMinPts = zVals[zVals < minz]
        nPts = len(belowMinPts)
        if nPts > 0:
            (lower, upper) = (belowMinPts.min(), belowMinPts.max())
            print(f"  {nPts} points are less than {minz}m (range {lower:.2f} to {upper:.2f})")
    if hdr.z_max > maxz:
        aboveMaxPts = zVals[zVals > maxz]
        nPts = len(aboveMaxPts)
        if nPts > 0:
            (lower, upper) = (aboveMaxPts.min(), aboveMaxPts.max())
            print(f"  {nPts} points are greater than {maxz}m (range {lower:.2f} to {upper:.2f})")


def getEPSGfromLAS(lasfile):
    """
    Check in the given las/laz file for a CRS, and get the EPSG number

    Parameters:
      lasfile (str): Name of LAS/LAZ file

    Returns:
      epsg (int): EPSG number of projection, or None
    """
    f = laspy.open(lasfile)
    crs = f.header.parse_crs()
    if crs is not None:
        epsg = crs.to_epsg()
    else:
        epsg = None
    return epsg


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
