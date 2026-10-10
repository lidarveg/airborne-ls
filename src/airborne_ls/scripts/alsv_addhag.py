#!/usr/bin/env python
"""
Add Height Above Ground (HAG) to one or more indexed LAZ files.

Adds the data as an EVLR record.

"""
import argparse
import glob

from airborne_ls import filenaming_methods, gridding_methods, lazfile_rw, qvf


def getCmdargs():
    """
    Get command line arguments
    """
    p = argparse.ArgumentParser()
    p.add_argument("--indir", help="Directory of LAZ files to process")
    p.add_argument("--infile", help="Individual LAZ file to process")
    p.add_argument("--groupMofN", nargs=2, metavar=('M', 'N'), type=int,
        help=("Use only with --indir. Divide the set of tiles to be processed into N groups, " +
              "and then process only the files in the M-th group (group numbering starts at 1)." +
              "This helps support efficient batch processing, where the user decides how " +
              "many batch jobs will run, and the command for each group processes only " +
              "the tiles for that group. For example, with 5 groups, the first group " +
              "would be specified as '--groupMofN 1 5'"))
    p.add_argument("--tilesize", default=1000, type=int,
        help=("Tile size (metres). We need this so we can find neighbouring tiles " +
              "(default=%(default)s"))
    p.add_argument("--binmargin", type=int, default=33,
        help=("Percentage of points from neighbouring bins to keep for per-bin " +
            "DEM interpolation (default=%(default)s). Smaller values will run faster, " +
            "but too small can leave extra holes in DEM"))

    cmdargs = p.parse_args()
    return cmdargs


def main():
    """
    Main routine
    """
    cmdargs = getCmdargs()

    if cmdargs.infile is not None:
        filelist = [cmdargs.infile]
    elif cmdargs.indir is not None:
        pattern = f"{cmdargs.indir}/*.laz"
        filelist = sorted(glob.glob(pattern))
        if cmdargs.groupMofN is not None:
            (group, numGroups) = tuple(cmdargs.groupMofN)
            filelist = filenaming_methods.filelistGroupSubset(filelist, group, numGroups)

    for filename in filelist:
        print(filename)
        tileName = qvf.getwhere(filename)
        (easting, northing, _) = filenaming_methods.decomposeWhereField(tileName)
        binnedData = lazfile_rw.readBinnedData(filename, tileSize=cmdargs.tilesize,
            withNeighbours=True, includeHAG=False)
        hag = gridding_methods.calcHeightAboveGroundAllBins(binnedData, cmdargs.binmargin,
                easting, northing)
        del binnedData

        lazfile_rw.addHagEvlr(filename, hag)


if __name__ == "__main__":
    main()
