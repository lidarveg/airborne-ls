#!/usr/bin/env python

"""

Purpose: run importation of supplied las/laz files for lidarveg

bin_size: in my experience 100m works well for lower density datasets (<10 pts/ mt2)
          and 50m works well for lower density datasets (> 10 pts/ mt2). 
          Not sure whether a finer bin size will be required (e.g. 25m)?
challenge: most project information (i.e. header info) cannot be relied upon to after this processing stage
options: educated guess which one to use by user (current approach) but a better approach could be to look at file size and 
point format from header (this should be correct) using the largest few files? 
        

"""

import argparse

# Configure logging
import logging
import sys
from pathlib import Path
import laspy
import numpy as np

from airborne_ls import lazfile_rw

logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.ERROR,
    format="%(asctime)s: %(name)20s: %(levelname)10s: %(message)s"
)

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
    parser.add_argument("--indir", required=True, help="Full path to directory containing LAS/LAZ files.")
    parser.add_argument("--outdr", required=True, help="Directory for newly named and indexed LAS/LAZ files.")
    parser.add_argument("--laz_flist", required=True, help="Text file containing LAS/LAZ files to be processed; one file per row.")

    # Tile dimensions
    parser.add_argument("--tile_s", required=True, type=float, help="Maximum XY dimension of LAS/LAZ file (metres).")
    parser.add_argument("--out_tile_s", default=1000, type=float, help="Equal or smaller maximum XY dimension for output LAS/LAZ files (metres).")

    # EPSG and spatial database options
    parser.add_argument("--epsg", type=int, required=True, help="EPSG code for map information.")
    
    # Metadata options
    parser.add_argument("--ss", type=str, default='ap', help="Platform type (e.g., 'ap' for airborne platform).")
    parser.add_argument("--ii", type=str, required=True, help="Predefined sensor code; can use uk if unknown")
    parser.add_argument("--pp", type=str, default='dr', help="Product type (e.g., 'dr' for discrete return).")
    parser.add_argument("--proj", required=True, help="Six-character project name (e.g., 'brisba').")
    parser.add_argument("--year", type=int, required=True, help="Year of data capture (e.g., 2022).")

    # Tile indexing options
    parser.add_argument("--binSize", default=50.0, type=float, help="XY bin size for data indexing (metres).")
    parser.add_argument("--startfilenum", default=0, type=int, help="Start position in file list for batch processing.")
    parser.add_argument("--stopfilenum", type=int, help="Stop position in file list for batch processing.")
    parser.add_argument("--memperjob", type=int, default=6, help="Memory limit in GB for each batch job.")
    parser.add_argument("--timeperjob", type=int, default=12, help="Time limit in hours for each batch job.")

    # Metadata flags 
        
    cmdargs = parser.parse_args(inputargs)

    # Validate input file list
    infilelist = Path(cmdargs.indir).joinpath(cmdargs.laz_flist)
    if not infilelist.is_file():
        laslist = list(Path(cmdargs.indir).glob('*.las'))
        lazlist = list(Path(cmdargs.indir).glob('*.laz'))

        if laslist:
            laslist = reorder_flist(laslist)
            cmdargs.laz_flist = 'lazlist_auto'
            outfile = Path(cmdargs.indir).joinpath(cmdargs.laz_flist)
            with open(outfile, 'w') as fout:
                fout.writelines(f"{Path(fn).name}\n" for fn in laslist)

        if lazlist:
            lazlist = reorder_flist(lazlist)
            cmdargs.laz_flist = 'lazlist_auto'
            outfile = Path(cmdargs.indir).joinpath(cmdargs.laz_flist)
            print(f"outfile: {outfile}")
            with open(outfile, 'w') as fout:
                fout.writelines(f"{Path(fn).name}\n" for fn in lazlist)

    infilelist = Path(cmdargs.indir).joinpath(cmdargs.laz_flist)
    print(f"infilelist: {infilelist}")
    if not infilelist.is_file():
        raise AssertionError("laz_flist is invalid")

    # Validate project name
    if len(cmdargs.proj) != 6:
        raise AssertionError("Project name must be exactly 6 characters.")

    # Validate output directory
    if not Path(cmdargs.outdr).exists():
        raise AssertionError("Output directory does not exist.")

    # required_args = ['year', 'epsg','tile_s',]
    # for argname in required_args:
    #     arg = cmdargs.__dict__.get(argname)
    #     if arg is None:
    #         raise argparse.ArgumentError(None, f"argument '--{argname}' is mandatory")

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
    lazlistfull = lazlistfull[cmdargs.startfilenum:cmdargs.stopfilenum]

    input_tileS = cmdargs.tile_s

    # Generate the filename base
    fn_base, ba3, zone_code = get_fn_base(cmdargs) # noqa

    for fn in lazlistfull:
        print(f"Processing LAS/LAZ file: {fn}")
        data = laspy.read(fn)
        
        # Calculate northing and easting
        northing = int(np.ceil(np.median(data.y[data.y > 0]) / input_tileS) * input_tileS)
        easting = int(np.floor(np.median(data.x[data.x > 0]) / input_tileS) * input_tileS)

        # Update filename base with northing and easting
        fn_base = fn_base.replace("EASTING_UL", str(int(easting))).replace(
            "NORTHING_UL", str(int(northing))
        )

        # Update filename with zone code
        pts = fn_base.split('_')
        fn_where = f"x{easting}ys{northing}z{zone_code}"
        fn_base = "_".join([pts[0], fn_where, pts[2], pts[3], pts[4]])
        #outfn = str(Path(cmdargs.outdr).joinpath(fn_base))

        # Check if tile size and bin size are divisible
        test = abs(
            ((input_tileS // cmdargs.binSize) * cmdargs.binSize)
            - ((input_tileS / cmdargs.binSize) * cmdargs.binSize)
        )
        if test >= 1.0:
            raise ValueError(f"Laz tile size and binSize are not divisible: {input_tileS} and {cmdargs.binSize}")

        # Run chunked LAS filtering
        _ = lazfile_rw.standardise_lasf(
            fn_base, cmdargs.outdr, data, easting, northing, input_tileS,
            cmdargs.out_tile_s, cmdargs.binSize, fn
        ) 

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


def get_fn_base(cmdargs):
    """
    Build the filename base, which will have the northing and easting updated.

    Stages defined as:
        - ba0: Ellipsoid heights
        - ba1: Geoid heights (AHD)
        - ba3: New indexed LAS

    Parameters:
        cmdargs (argparse.Namespace): Parsed command-line arguments.

    Returns:
        tuple: Filename base, ba3 filename, and zone code.
    """
    zone = str(cmdargs.epsg)[-1:]
    zone_code = str(cmdargs.epsg)[-2:]

    easting_ul = "EASTING_UL"
    northing_ul = "NORTHING_UL"

    # Determine zone prefix
    if 28350 < cmdargs.epsg < 28360:
        zone_prefix = 'm'
    elif 7850 < cmdargs.epsg < 7860:
        zone_prefix = 'd'
    else:
        raise ValueError(f"Unknown EPSG code: {cmdargs.epsg}")

    # Build filename base
    fn_what = f"{cmdargs.ss}{cmdargs.ii}{cmdargs.pp}"
    fn_where = f"x{easting_ul}ys{northing_ul}"
    fn_when = f"{cmdargs.year}_ba1{zone_prefix}{zone}_p{cmdargs.proj}.laz"    
    fn_base = f"{fn_what}_{fn_where}_{fn_when}"
    ba3 = f"{fn_what}_r{cmdargs.proj}_{cmdargs.year}_ba3{zone_prefix}{zone}.zip"

    return fn_base, ba3, zone_code


def check_input_fns(indir, infilelist):
    """
    Check whether input files are valid and are .laz or .las files.

    Parameters:
        indir (str): Input directory containing LAS/LAZ files.
        infilelist (Path): Path to the file containing the list of LAS/LAZ files.

    Returns:
        tuple: A tuple containing the list of valid files and a boolean indicating if all files are valid.
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
